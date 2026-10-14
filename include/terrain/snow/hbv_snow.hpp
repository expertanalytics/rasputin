#pragma once

// The HBV snow routine (Lindström et al. 1997, HBV-96) on many independent
// points, here the triangles of a mesh, with the snow inside each spread over N
// equal-area bins of fixed relative depth and the melt optionally corrected by
// a radiation anomaly per triangle. Parameter names are Shyft's (tx, ts, cx,
// lw, cfr). The derivation, the invariants and the fixtures are in
// docs/increments/46-hbv-snow.md, section 4.2.
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

// N equal-area bins, bin k holding f_k times the triangle's mean snowfall.
// The factors are refused unless 1 <= N <= 64, finite, non-negative,
// ascending and of mean within 1e-9 of 1; accepted ones are divided by their
// mean, so the snowfall spread over the bins sums to the snowfall to rounding.
class SnowDistribution {
  public:
    explicit SnowDistribution(std::span<const double> factors) : f_(factors.begin(), factors.end()) {
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
        for (double& v : f_)
            v /= mean;
    }
    [[nodiscard]] std::span<const double> factors() const noexcept { return f_; }

  private:
    std::vector<double> f_;
};

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

// Per triangle i and bin k (triangle-major, slot i*N + k): the solid (sp) and
// liquid (sw) water of the pack in mm; per triangle the bins' mean swe, the
// fraction of bins with snow (sca) and the mean outflow of the last step.
class SnowState {
  public:
    // Every bin of every triangle holds f_k*initial_swe_mm of ice and no water.
    [[nodiscard]] static SnowState uniform(std::size_t triangles, const SnowDistribution& dist,
                                           double initial_swe_mm) {
        if (!(std::isfinite(initial_swe_mm) && initial_swe_mm >= 0.0))
            throw std::invalid_argument("the initial SWE must be finite and not negative");
        SnowState s;
        s.f_.assign(dist.factors().begin(), dist.factors().end());
        const std::size_t bins = s.f_.size();
        s.sp_.resize(triangles * bins);
        s.sw_.assign(triangles * bins, 0.0);
        s.swe_.resize(triangles);
        s.sca_.resize(triangles);
        s.outflow_.assign(triangles, 0.0);
        for (std::size_t i = 0; i < triangles; ++i) {
            for (std::size_t k = 0; k < bins; ++k)
                s.sp_[i * bins + k] = s.f_[k] * initial_swe_mm;
            s.summarise(i);
        }
        return s;
    }

    [[nodiscard]] std::size_t triangles() const noexcept { return swe_.size(); }
    [[nodiscard]] std::size_t bins() const noexcept { return f_.size(); }
    [[nodiscard]] std::span<const double> sp() const noexcept { return sp_; }
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
            covered += sp_[k] > 0.0 ? 1.0 : 0.0;
        }
        swe_[i] = water / static_cast<double>(bins);
        sca_[i] = covered / static_cast<double>(bins);
    }

    std::vector<double> f_, sp_, sw_, swe_, sca_, outflow_;
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
                sp += s.f_[k] * snow;
                if (sp > 0.0)
                    sw += rain;
                else
                    out += rain;  // a bare bin keeps no rain
                if (t > theta) {
                    const double melt = std::min(m * days * (t - theta), sp);
                    sp -= melt;
                    sw += melt;
                } else if (t < p.ts) {
                    const double refreeze = std::min(p.cfr * p.cx * days * (p.ts - t), sw);
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
