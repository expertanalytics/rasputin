#pragma once

// The HBV snow routine (Lindström et al. 1997, HBV-96) on many independent
// points, here the triangles of a mesh, with the snow inside each spread over N
// equal-area bins of fixed relative depth and the melt optionally corrected by
// a radiation anomaly per triangle. Within a bin the solid depth is a straight
// line, so the cover of a bin whose melt front lies inside it is a fraction
// (Shyft's continuous snow-covered area); a flat bin is plain HBV-96. Parameter
// names are Shyft's (tx, ts, cx, lw, cfr). The derivation, the invariants and
// the fixtures are in docs/increments/46-hbv-snow.md, sections 4.2 and 12.3.
//
// No file, no path, no clock: the step sees arrays and a length in seconds.
// Each triangle's update writes only its own slots, so the step runs in
// parallel over triangles and its result does not depend on the thread count.

#include <terrain/parallel_util/chunks.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace terrain::snow {

struct HbvSnowParameters {
    double tx;   // degC: precipitation falls as snow below it
    double ts;   // degC: melt threshold
    double cx;   // mm per degC per day: degree-day factor, >= 0
    double lw;   // liquid water the pack holds, fraction of its solid part, >= 0
    double cfr;  // refreeze coefficient: refreeze = cfr*cx*dt*(ts - T), >= 0
};

enum class SunCorrection : std::uint8_t { none, melt_factor, threshold };

struct SunCorrectionParameters {
    SunCorrection mode;
    double r;  // melt_factor: mm per degC per day per W/m2
    double d;  // threshold: degC per W/m2
};

// N equal-area bins, bin k holding f_k times the triangle's mean snowfall, and
// deep_k times it at the bin's deep edge (its depth runs linearly from
// 2 f_k - deep_k to deep_k). The factors are refused unless 1 <= N <= 64,
// finite, non-negative, ascending and of mean within 1e-9 of 1; deep is refused
// unless empty (flat: deep = factors) or N finite values with
// f_k <= deep_k <= 2 f_k. Both are divided by the factors' mean, so the
// snowfall spread over the bins sums to the snowfall to rounding.
class SnowDistribution {
  public:
    explicit SnowDistribution(std::span<const double> factors, std::span<const double> deep = {})
        : f_(factors.begin(), factors.end()), d_(deep.begin(), deep.end()) {
        if (f_.empty() || f_.size() > 64)
            throw std::invalid_argument("the number of snow bins must be between 1 and 64");
        if (!std::ranges::all_of(f_, [](double v) { return std::isfinite(v); }))
            throw std::invalid_argument("every snow bin factor must be finite");
        if (std::ranges::any_of(f_, [](double v) { return v < 0.0; }))
            throw std::invalid_argument("a snow bin factor is negative");
        if (!std::ranges::is_sorted(f_))
            throw std::invalid_argument("the snow bin factors must be ascending");
        double sum = 0.0;
        for (double v : f_)
            sum += v;
        const double mean = sum / static_cast<double>(f_.size());
        if (!(std::abs(mean - 1.0) <= 1e-9))
            throw std::invalid_argument("the snow bin factors' mean must be 1 within 1e-9");
        if (d_.empty())
            d_ = f_;
        if (d_.size() != f_.size())
            throw std::invalid_argument("the snow bins' deep edges must be empty or one per bin");
        if (!std::ranges::all_of(d_, [](double v) { return std::isfinite(v); }))
            throw std::invalid_argument("every snow bin's deep edge must be finite");
        for (std::size_t k = 0; k < f_.size(); ++k)
            if (!(f_[k] <= d_[k] && d_[k] <= 2.0 * f_[k]))
                throw std::invalid_argument("a bin's deep edge must lie between its mean and twice its mean");
        for (std::size_t k = 0; k < f_.size(); ++k) {
            f_[k] /= mean;
            d_[k] /= mean;
        }
    }
    [[nodiscard]] std::span<const double> factors() const noexcept { return f_; }
    [[nodiscard]] std::span<const double> deep() const noexcept { return d_; }

  private:
    std::vector<double> f_, d_;
};

// The covered fraction of a bin of mean solid sp and deep-edge solid deep:
// 1 while the line's shallow edge 2 sp - deep is not below zero, else the
// share of the bin where the line is above zero.
[[nodiscard]] inline double cover(double sp, double deep) noexcept {
    if (sp == 0.0)
        return 0.0;
    return 2.0 * sp >= deep ? 1.0 : 2.0 * sp / deep;
}

struct StepForcing {
    std::span<const double> temperature;        // degC, n
    std::span<const double> precipitation;      // mm over the step, n, >= 0
    std::span<const double> radiation_anomaly;  // W/m2, n; empty iff mode none
    double dt_seconds;                          // > 0
};

enum class StepStatus : std::uint8_t {
    ok,
    invalid_parameters,
    size_mismatch,
    non_finite_forcing,
    negative_precipitation,
    bad_step_length
};

class SnowState;
StepStatus step(const HbvSnowParameters&, const SunCorrectionParameters&, SnowState&,
                const StepForcing&, unsigned threads);

// Per triangle i and bin k (triangle-major, slot i*N + k): the mean solid (sp),
// the solid at the deep edge (deep) and the liquid (sw) water of the pack in mm,
// with 0 <= sp <= deep and sp = 0 exactly when deep = 0; per triangle the bins'
// mean swe, their mean cover (sca) and the mean outflow of the last step.
class SnowState {
  public:
    // Every bin of every triangle holds f_k*initial_swe_mm of ice and no water.
    [[nodiscard]] static SnowState uniform(std::size_t triangles, const SnowDistribution& dist,
                                           double initial_swe_mm) {
        if (!(std::isfinite(initial_swe_mm) && initial_swe_mm >= 0.0))
            throw std::invalid_argument("the initial SWE must be finite and not negative");
        SnowState s;
        s.f_.assign(dist.factors().begin(), dist.factors().end());
        s.d_.assign(dist.deep().begin(), dist.deep().end());
        const std::size_t bins = s.f_.size();
        s.sp_.resize(triangles * bins);
        s.deep_.resize(triangles * bins);
        s.sw_.assign(triangles * bins, 0.0);
        s.swe_.resize(triangles);
        s.sca_.resize(triangles);
        s.outflow_.assign(triangles, 0.0);
        for (std::size_t i = 0; i < triangles; ++i) {
            for (std::size_t k = 0; k < bins; ++k) {
                s.sp_[i * bins + k] = s.f_[k] * initial_swe_mm;
                s.deep_[i * bins + k] = s.d_[k] * initial_swe_mm;
            }
            s.summarise(i);
        }
        return s;
    }

    [[nodiscard]] std::size_t triangles() const noexcept { return swe_.size(); }
    [[nodiscard]] std::size_t bins() const noexcept { return f_.size(); }
    [[nodiscard]] std::span<const double> sp() const noexcept { return sp_; }
    [[nodiscard]] std::span<const double> deep() const noexcept { return deep_; }
    [[nodiscard]] std::span<const double> sw() const noexcept { return sw_; }
    [[nodiscard]] std::span<const double> swe() const noexcept { return swe_; }
    [[nodiscard]] std::span<const double> sca() const noexcept { return sca_; }
    [[nodiscard]] std::span<const double> outflow() const noexcept { return outflow_; }

  private:
    friend StepStatus step(const HbvSnowParameters&, const SunCorrectionParameters&, SnowState&,
                           const StepForcing&, unsigned threads);

    SnowState() = default;

    // swe and sca of triangle i from its bins.
    void summarise(std::size_t i) noexcept {
        const std::size_t bins = f_.size();
        double water = 0.0, covered = 0.0;
        for (std::size_t k = i * bins; k < (i + 1) * bins; ++k) {
            water += sp_[k] + sw_[k];
            covered += cover(sp_[k], deep_[k]);
        }
        swe_[i] = water / static_cast<double>(bins);
        sca_[i] = covered / static_cast<double>(bins);
    }

    std::vector<double> f_, d_, sp_, deep_, sw_, swe_, sca_, outflow_;
};

// One step of length dt for every triangle. Anything other than ok is returned
// before a slot is written, so a refused step leaves the state as it was.
inline StepStatus step(const HbvSnowParameters& p, const SunCorrectionParameters& c,
                       SnowState& s, const StepForcing& f, unsigned threads) {
    const auto finite = [](double v) { return std::isfinite(v); };
    const auto finite_all = [&](std::span<const double> v) { return std::ranges::all_of(v, finite); };
    const bool known_mode = c.mode == SunCorrection::none || c.mode == SunCorrection::melt_factor ||
                            c.mode == SunCorrection::threshold;
    const double params[]{p.tx, p.ts, p.cx, p.lw, p.cfr, c.r, c.d};
    if (!known_mode || !finite_all(params) || p.cx < 0.0 || p.lw < 0.0 || p.cfr < 0.0)
        return StepStatus::invalid_parameters;
    const std::size_t n = s.triangles();
    const std::size_t want_anomaly = c.mode == SunCorrection::none ? 0 : n;
    if (f.temperature.size() != n || f.precipitation.size() != n ||
        f.radiation_anomaly.size() != want_anomaly)
        return StepStatus::size_mismatch;
    if (!(std::isfinite(f.dt_seconds) && f.dt_seconds > 0.0))
        return StepStatus::bad_step_length;
    if (!finite_all(f.temperature) || !finite_all(f.precipitation) || !finite_all(f.radiation_anomaly))
        return StepStatus::non_finite_forcing;
    if (std::ranges::any_of(f.precipitation, [](double v) { return v < 0.0; }))
        return StepStatus::negative_precipitation;

    const double days = f.dt_seconds / 86400.0;
    const std::size_t bins = s.bins();
    parallel_util::for_each_block(n, threads, {}, [&](std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            const double t = f.temperature[i], precip = f.precipitation[i];
            const double snow = t < p.tx ? precip : 0.0;
            const double rain = precip - snow;
            double m = p.cx, theta = p.ts;
            if (c.mode == SunCorrection::melt_factor)
                m = std::max(0.0, p.cx + c.r * f.radiation_anomaly[i]);
            else if (c.mode == SunCorrection::threshold)
                theta = p.ts - c.d * f.radiation_anomaly[i];
            double out = 0.0;
            for (std::size_t k = 0; k < bins; ++k) {
                double& sp = s.sp_[i * bins + k];
                double& sw = s.sw_[i * bins + k];
                double& deep = s.deep_[i * bins + k];
                // Snow on a partly bare bin first spreads the old snow over the
                // whole bin as a wedge from zero, keeping its mass (Shyft's rule).
                if (snow > 0.0 && 2.0 * sp < deep)
                    deep = 2.0 * sp;
                sp += s.f_[k] * snow;
                deep += s.d_[k] * snow;
                const double covered = cover(sp, deep);
                sw += covered * rain;
                out += (1.0 - covered) * rain;  // the bare part keeps no rain
                if (t > theta) {
                    // M per unit of covered area; each branch is the area under
                    // the line lowered by M, cut at zero.
                    const double melt_potential = m * days * (t - theta);
                    double melt = sp;
                    if (deep <= melt_potential) {  // a bare bin too: melts its 0
                        sp = deep = 0.0;
                    } else if (2.0 * sp >= deep && melt_potential <= 2.0 * sp - deep) {
                        melt = melt_potential;
                        sp -= melt;
                        deep -= melt;
                    } else {
                        const double front = deep - melt_potential, ratio = front / deep;
                        const double left = 2.0 * sp >= deep ? front * front / (4.0 * (deep - sp))
                                                             : sp * ratio * ratio;
                        melt = sp - left;
                        sp = left;
                        deep = left == 0.0 ? 0.0 : front;
                    }
                    sw += melt;
                } else if (t < p.ts) {
                    const double refreeze = std::min(covered * p.cfr * p.cx * days * (p.ts - t), sw);
                    if (2.0 * sp >= deep)
                        deep += refreeze;
                    else
                        deep = deep * (sp + refreeze) / sp;  // keeps the cover
                    sp += refreeze;
                    sw -= refreeze;
                }
                out += std::max(0.0, sw - p.lw * sp);
                sw = std::min(sw, p.lw * sp);
            }
            s.outflow_[i] = out / static_cast<double>(bins);
            s.summarise(i);
        }
    });
    return StepStatus::ok;
}

}  // namespace terrain::snow
