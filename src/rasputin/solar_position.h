//
// Created by Ola Skavhaug on 26/08/2019.
//

#pragma once

#include <cmath>
#include <chrono>
#ifndef __clang__
#include "date/date.h"
#endif
#include <ctime>

namespace rasputin::solar_position::collectors {

auto azimuth_and_elevation() {
    return [] (double e0, double Gamma, double Phi, const double alpha_mark, const double delta_mark, const double H_mark) {
        return std::make_tuple(Phi, e0);
    };
}

}

namespace rasputin::solar_position::delta_t_calculator {
const std::vector<double> _dts{31.1, 33.2, 35.7, 40.2, 45.5, 50.5, 54.3, 56.9, 60.8, 63.8, 64.7};

auto parabolic_dt_calc() {
    return [] (const double year) {
         return -20 + 32*pow((year - 1820)/20.0, 2);
    };
}

auto coarse_date_calc() {
    return [] (const unsigned int year,
               const unsigned int month,
               const unsigned int day) {
        if (year < 1955 || year > 2009)
            return parabolic_dt_calc()(year);
        const auto idx = (year - 1955) / 5;
        return _dts[idx];
    };
}

auto coarse_timestamp_calc() {
    using namespace std::chrono;
#ifndef __clang__
    using namespace date;
#endif

    auto ep = sys_days{January/1/1970};

    return [ep] (system_clock::time_point tp) {
        const auto year = round(duration_cast<seconds>(tp - ep).count()*1.0/duration_cast<seconds>(years(1)).count()) + 1970;
        if (year < 1955 || year > 2009)
            return parabolic_dt_calc()(year);
        const auto idx = (year - 1955) / 5;
        return _dts[idx];
    };
}

}

namespace rasputin::solar_position {

auto limit_degrees(const double degrees) {
    return fmod(fmod(degrees, 360) + 360, 360);
}

auto d2r(const double d) {
    return d*M_PI/180.0;
}

auto r2d(const double r) {
    return r*180.0/M_PI;
}

auto jd_from_cal(unsigned int year, unsigned int month, const double day){

    if (month < 3) {
        month += 12;
        --year;
    }

    // Equation (4)
    const auto jd_temp = trunc(365.25*(year + 4716)) + trunc(30.6001*(month + 1)) + day - 1524.5;
    if (jd_temp < 2299160)
        return jd_temp;
    const auto A = trunc(year/100);
    return jd_temp + 2 - A + trunc(A/4);
}

double jd_from_clock(const std::chrono::system_clock::time_point tp) {
    using namespace std::chrono;
#ifndef __clang__
    using namespace date;
#endif
    const auto dt = tp - sys_days{January/1/1970};
    const auto t = duration_cast<seconds>(dt).count();
    return t/86400.0 + 2440587.5;
}

double jde(const double julian_day, const double dt){
    return julian_day + dt/86400.0;
}

double jc(const double julian_day) {
    return (julian_day - 2451545.0)/36525.0;
}

double jce(const double julian_ephemeris_day) {
    return (julian_ephemeris_day - 2451545.0)/36525.0;
}

double jme(const double julian_ephemeris_centry) {
    return julian_ephemeris_centry/10.0;
}

double heliocentric_longitude(const double julian_ephemeris_millennium) {
    const double jme = julian_ephemeris_millennium;
    const double L00  = 175347046.0;
    const double L01  = 3341656*cos(4.6692568 +  6283.07585 *jme);
    const double L02  =   34894*cos(4.6261    + 12566.1517  *jme);
    const double L03  =    3497*cos(2.7441    +  5753.3849  *jme);
    const double L04  =    3418*cos(2.8289    +     3.5231  *jme);
    const double L05  =    3136*cos(3.6277    + 77713.7715  *jme);
    const double L06  =    2676*cos(4.4181    +  7860.4194  *jme);
    const double L07  =    2343*cos(6.1352    +  3930.2097  *jme);
    const double L08  =    1324*cos(0.7425    + 11506.7698  *jme);
    const double L09  =    1273*cos(2.0371    +   529.691   *jme);
    const double L010 =    1199*cos(1.1096    +  1577.3435  *jme);
    const double L011 =     990*cos(5.233     +  5884.927   *jme);
    const double L012 =     902*cos(2.045     +    26.298   *jme);
    const double L013 =     857*cos(3.508     +   398.149   *jme);
    const double L014 =     780*cos(1.179     +  5223.694   *jme);
    const double L015 =     753*cos(2.533     +  5507.553   *jme);
    const double L016 =     505*cos(4.583     + 18849.228   *jme);
    const double L017 =     492*cos(4.205     +   775.523   *jme);
    const double L018 =     357*cos(2.92      +     0.067   *jme);
    const double L019 =     317*cos(5.849     + 11790.629   *jme);
    const double L020 =     284*cos(1.899     +   796.298   *jme);
    const double L021 =     271*cos(0.315     + 10977.079   *jme);
    const double L022 =     243*cos(0.345     +  5486.778   *jme);
    const double L023 =     206*cos(4.806     +  2544.314   *jme);
    const double L024 =     205*cos(1.869     +  5573.143   *jme);
    const double L025 =     202*cos(2.458     +  6069.777   *jme);
    const double L026 =     156*cos(0.833     +   213.299   *jme);
    const double L027 =     132*cos(3.411     +  2942.463   *jme);
    const double L028 =     126*cos(1.083     +    20.775   *jme);
    const double L029 =     115*cos(0.645     +     0.98    *jme);
    const double L030 =     103*cos(0.636     +  4694.003   *jme);
    const double L031 =     102*cos(0.976     + 15720.839   *jme);
    const double L032 =     102*cos(4.267     +     7.114   *jme);
    const double L033 =      99*cos(6.21      +  2146.17    *jme);
    const double L034 =      98*cos(0.68      +   155.42    *jme);
    const double L035 =      86*cos(5.98      +161000.69    *jme);
    const double L036 =      85*cos(1.3       +  6275.96    *jme);
    const double L037 =      85*cos(3.67      + 71430.7     *jme);
    const double L038 =      80*cos(1.81      + 17260.15    *jme);
    const double L039 =      79*cos(3.04      + 12036.46    *jme);
    const double L040 =      75*cos(1.76      +  5088.63    *jme);
    const double L041 =      74*cos(3.5       +  3154.69    *jme);
    const double L042 =      74*cos(4.68      +   801.82    *jme);
    const double L043 =      70*cos(0.83      +  9437.76    *jme);
    const double L044 =      62*cos(3.98      +  8827.39    *jme);
    const double L045 =      61*cos(1.82      +  7084.9     *jme);
    const double L046 =      57*cos(2.78      +  6286.6     *jme);
    const double L047 =      56*cos(4.39      + 14143.5     *jme);
    const double L048 =      56*cos(3.47      +  6279.55    *jme);
    const double L049 =      52*cos(0.19      + 12139.55    *jme);
    const double L050 =      52*cos(1.33      +  1748.02    *jme);
    const double L051 =      51*cos(0.28      +  5856.48    *jme);
    const double L052 =      49*cos(0.49      +  1194.45    *jme);
    const double L053 =      41*cos(5.37      +  8429.24    *jme);
    const double L054 =      41*cos(2.4       + 19651.05    *jme);
    const double L055 =      39*cos(6.17      + 10447.39    *jme);
    const double L056 =      37*cos(6.04      + 10213.29    *jme);
    const double L057 =      37*cos(2.57      +  1059.38    *jme);
    const double L058 =      36*cos(1.71      +  2352.87    *jme);
    const double L059 =      36*cos(1.78      +  6812.77    *jme);
    const double L060 =      33*cos(0.59      + 17789.85    *jme);
    const double L061 =      30*cos(0.44      + 83996.85    *jme);
    const double L062 =      30*cos(2.74      +  1349.87    *jme);
    const double L063 =      25*cos(3.16      +  4690.48    *jme);


    const double L0 = L00 + L01 + L02 + L03 + L04 + L05 + L06 + L07 + L08 + L09 +
            L010 + L011 + L012 + L013 + L014 + L015 + L016 + L017 + L018 + L019 +
            L020 + L021 + L022 + L023 + L024 + L025 + L026 + L027 + L028 + L029 +
            L030 + L031 + L032 + L033 + L034 + L035 + L036 + L037 + L038 + L039 +
            L040 + L041 + L042 + L043 + L044 + L045 + L046 + L047 + L048 + L049 +
            L050 + L051 + L052 + L053 + L054 + L055 + L056 + L057 + L058 + L059 +
            L060 + L061 + L062 + L063;


    const double  L10 = 628331966747;
    const double  L11 = 206059*cos(2.678235 +  6283.07585*jme);
    const double  L12 =   4303*cos(2.6351   + 12566.1517 *jme);
    const double  L13 =    425*cos(1.59     +     3.523  *jme);
    const double  L14 =    119*cos(5.796    +    26.298  *jme);
    const double  L15 =    109*cos(2.966    +  1577.344  *jme);
    const double  L16 =     93*cos(2.59     + 18849.23   *jme);
    const double  L17 =     72*cos(1.14     +   529.69   *jme);
    const double  L18 =     68*cos(1.87     +   398.15   *jme);
    const double  L19 =     67*cos(4.41     +  5507.55   *jme);
    const double L110 =     59*cos(2.89     +  5223.69   *jme);
    const double L111 =     56*cos(2.17     +   155.42   *jme);
    const double L112 =     45*cos(0.4      +   796.3    *jme);
    const double L113 =     36*cos(0.47     +   775.52   *jme);
    const double L114 =     29*cos(2.65     +     7.11   *jme);
    const double L115 =     21*cos(5.34     +     0.98   *jme);
    const double L116 =     19*cos(1.85     +  5486.78   *jme);
    const double L117 =     19*cos(4.97     +   213.3    *jme);
    const double L118 =     17*cos(2.99     +  6275.96   *jme);
    const double L119 =     16*cos(0.03     +  2544.31   *jme);
    const double L120 =     16*cos(1.43     +  2146.17   *jme);
    const double L121 =     15*cos(1.21     + 10977.08   *jme);
    const double L122 =     12*cos(2.83     +  1748.02   *jme);
    const double L123 =     12*cos(3.26     +  5088.63   *jme);
    const double L124 =     12*cos(5.27     +  1194.45   *jme);
    const double L125 =     12*cos(2.08     +  4694      *jme);
    const double L126 =     11*cos(0.77     +   553.57   *jme);
    const double L127 =     10*cos(1.3      +  6286.6    *jme);
    const double L128 =     10*cos(4.24     +  1349.87   *jme);
    const double L129 =      9*cos(2.7      +   242.73   *jme);
    const double L130 =      9*cos(5.64     +   951.72   *jme);
    const double L131 =      8*cos(5.3      +  2352.87   *jme);
    const double L132 =      6*cos(2.65     +  9437.76   *jme);
    const double L133 =      6*cos(4.67     +  4690.48   *jme);

    const double L1 = L10 + L11 + L12 + L13 + L14 + L15 + L16 + L17 + L18 + L19 +
            L110 + L111 + L112 + L113 + L114 + L115 + L116 + L117 + L118 + L119 +
            L120 + L121 + L122 + L123 + L124 + L125 + L126 + L127 + L128 + L129 +
            L130 + L131 + L132 + L133;

     const double L20  = 52919;
     const double L21  = 8720*cos(1.0721 +  6283.0758*jme);
     const double L22  =  309*cos(0.867  + 12566.152 *jme);
     const double L23  =   27*cos(0.05   +     3.52  *jme);
     const double L24  =   16*cos(5.19   +    26.3   *jme);
     const double L25  =   16*cos(3.68   +   155.42  *jme);
     const double L26  =   10*cos(0.76   + 18849.23  *jme);
     const double L27  =    9*cos(2.06   + 77713.77  *jme);
     const double L28  =    7*cos(0.83   +   775.52  *jme);
     const double L29  =    5*cos(4.66   +  1577.34  *jme);
     const double L210 =    4*cos(1.03   +     7.11  *jme);
     const double L211 =    4*cos(3.44   +  5573.14  *jme);
     const double L212 =    3*cos(5.14   +   796.3   *jme);
     const double L213 =    3*cos(6.05   +  5507.55  *jme);
     const double L214 =    3*cos(1.19   +   242.73  *jme);
     const double L215 =    3*cos(6.12   +   529.69  *jme);
     const double L216 =    3*cos(0.31   +   398.15  *jme);
     const double L217 =    3*cos(2.28   +   553.57  *jme);
     const double L218 =    2*cos(4.38   +  5223.69  *jme);
     const double L219 =    2*cos(3.75   +     0.98  *jme);

    const double L2 = L20 + L21 + L22 + L23 + L24 + L25 + L26 + L27 + L28 + L29 +
                      L210 + L211 + L212 + L213 + L214 + L215 + L216 + L217 + L218 + L219;

    const double L30 = 289*cos(5.844 +  6283.076*jme);
    const double L31 =  35;
    const double L32 =  17*cos(5.49  + 12566.15 *jme);
    const double L33 =   3*cos(5.2   +   155.42 *jme);
    const double L34 =     cos(4.72  +     3.52 *jme);
    const double L35 =     cos(5.3   + 18849.23 *jme);
    const double L36 =     cos(5.97  +   242.73 *jme);

    const double L3 = L30 + L31 + L32 + L33 + L34 + L35 + L36;

    const double L40 = 114*cos(3.142);
    const double L41 =   8*cos(4.13 +  6283.08*jme);
    const double L42 =     cos(3.84 + 12566.15*jme);
    const double L4 = L40 + L41 + L42;

    const double L5 = cos(3.14);

    const double L_rad = (L0 + L1*jme + L2*pow(jme, 2) + L3*pow(jme, 3) + L4*pow(jme, 4) + L5*pow(jme, 5))*1e-8;
    const double L_deg = (L_rad*180)/M_PI;

    // Limit answer to [0, 360)
    return limit_degrees(L_deg);
}

double heliocentric_latitude(double julian_ephemeris_millennium) {
    const double jme = julian_ephemeris_millennium;

    const double B00 = 280*cos(3.199 + 84334.662*jme);
    const double B01 = 102*cos(5.422 +  5507.553*jme);
    const double B02 =  80*cos(3.88  +  5223.69 *jme);
    const double B03 =  44*cos(3.7   +  2352.87 *jme);
    const double B04 =  32*cos(4     +  1577.34 *jme);
    const double B0 = B00 + B01 + B02 + B03 + B04;
    
    const double B10 = 9*cos(3.9  + 5507.55*jme);
    const double B11 = 6*cos(1.73 + 5223.69*jme);
    const double B1 = B10 + B11;

    const double B_rad = (B0 + B1*jme)*1e-8;
    const double B_deg = (B_rad*180)/M_PI;
    return B_deg;
}

double heliocentric_radius_vector(double julian_ephemeris_millennium) {
    const double jme = julian_ephemeris_millennium;

    const double R00  = 100013989;
    const double R01  =   1670700*cos(3.0984635 +   6283.07585*jme);
    const double R02  =     13956*cos(3.05525   +  12566.1517 *jme);
    const double R03  =      3084*cos(5.1985    +  77713.7715 *jme);
    const double R04  =      1628*cos(1.1739    +   5753.3849 *jme);
    const double R05  =      1576*cos(2.8469    +   7860.4194 *jme);
    const double R06  =       925*cos(5.453     +  11506.77   *jme);
    const double R07  =       542*cos(4.564     +   3930.21   *jme);
    const double R08  =       472*cos(3.661     +   5884.927  *jme);
    const double R09  =       346*cos(0.964     +   5507.553  *jme);
    const double R010 =       329*cos(5.9       +   5223.694  *jme);
    const double R011 =       307*cos(0.299     +   5573.143  *jme);
    const double R012 =       243*cos(4.273     +  11790.629  *jme);
    const double R013 =       212*cos(5.847     +   1577.344  *jme);
    const double R014 =       186*cos(5.022     +  10977.079  *jme);
    const double R015 =       175*cos(3.012     +  18849.228  *jme);
    const double R016 =       110*cos(5.055     +   5486.778  *jme);
    const double R017 =        98*cos(0.89      +   6069.78   *jme);
    const double R018 =        86*cos(5.69      +  15720.84   *jme);
    const double R019 =        86*cos(1.27      + 161000.69   *jme);
    const double R020 =        65*cos(0.27      +  17260.15   *jme);
    const double R021 =        63*cos(0.92      +    529.69   *jme);
    const double R022 =        57*cos(2.01      +  83996.85   *jme);
    const double R023 =        56*cos(5.24      +  71430.7    *jme);
    const double R024 =        49*cos(3.25      +   2544.31   *jme);
    const double R025 =        47*cos(2.58      +    775.52   *jme);
    const double R026 =        45*cos(5.54      +   9437.76   *jme);
    const double R027 =        43*cos(6.01      +   6275.96   *jme);
    const double R028 =        39*cos(5.36      +   4694      *jme);
    const double R029 =        38*cos(2.39      +   8827.39   *jme);
    const double R030 =        37*cos(0.83      +  19651.05   *jme);
    const double R031 =        37*cos(4.9       +  12139.55   *jme);
    const double R032 =        36*cos(1.67      +  12036.46   *jme);
    const double R033 =        35*cos(1.84      +   2942.46   *jme);
    const double R034 =        33*cos(0.24      +   7084.9    *jme);
    const double R035 =        32*cos(0.18      +   5088.63   *jme);
    const double R036 =        32*cos(1.78      +    398.15   *jme);
    const double R037 =        28*cos(1.21      +   6286.6    *jme);
    const double R038 =        28*cos(1.9       +   6279.55   *jme);
    const double R039 =        26*cos(4.59      +  10447.39   *jme);
    const double R0 = R00 + R01 + R02 + R03 + R04 + R05 + R06 + R07 + R08 + R09 +
            R010 + R011 + R012 + R013 + R014 + R015 + R016 + R017 + R018 + R019 +
            R020 + R021 + R022 + R023 + R024 + R025 + R026 + R027 + R028 + R029 +
            R030 + R031 + R032 + R033 + R034 + R035 + R036 + R037 + R038 + R039;

    const double R10 = 103019*cos(1.10749 +  6283.07585*jme);
    const double R11 =   1721*cos(1.0644  + 12566.1517 *jme);
    const double R12 =    702*cos(3.142);
    const double R13 =     32*cos(1.02    + 18849.23   *jme);
    const double R14 =     31*cos(2.84    +  5507.55   *jme);
    const double R15 =     25*cos(1.32    +  5223.69   *jme);
    const double R16 =     18*cos(1.42    +  1577.34   *jme);
    const double R17 =     10*cos(5.91    + 10977.08   *jme);
    const double R18 =      9*cos(1.42    +  6275.96   *jme);
    const double R19 =      9*cos(0.27    +  5486.78   *jme);
    const double R1 = R10 + R11 + R12 + R13 + R14 + R15 + R16 + R17 + R18 + R19;

    const double R20 = 4359*cos(5.7846 +  6283.0758*jme);
    const double R21 =  124*cos(5.579  + 12566.152 *jme);
    const double R22 =   12*cos(3.14);
    const double R23 =    9*cos(3.63   + 77713.77  *jme);
    const double R24 =    6*cos(1.87    + 5573.14  *jme);
    const double R25 =    3*cos(5.47   + 18849.23  *jme);
    const double R2 = R20 + R21 + R22 + R23 + R24 + R25;

    const double R30 = 145*cos(4.273 +  6283.076*jme);
    const double R31 =   7*cos(3.92  + 12566.15 *jme);
    const double R3 = R30 + R31;
    const double R4 = 4*cos(2.56 + 6283.08*jme);

    return (R0 + R1*jme + R2*pow(jme, 2) + R3*pow(jme, 3) + R4*pow(jme, 4))*1e-8;
}

double geocentric_longitude(const double heliocentric_longitude) {
    return heliocentric_longitude + 180;
}

double geocentric_latitude(const double heliocentric_latitude) {
    return -heliocentric_latitude;
}


auto nutation(const double julian_ephemeris_centry) {
    const auto jce = julian_ephemeris_centry;
    const auto jce2 = jce*jce;
    const auto jce3 = jce2*jce;

    const auto X0 = 297.85036 + 445267.111480*jce - 0.0019142*jce2 + jce3/189474;
    const auto X1 = 357.52772 +  35999.050340*jce - 0.0001603*jce2 - jce3/300000;
    const auto X2 = 134.96298 + 477198.867398*jce + 0.0086972*jce2 + jce3/56250;
    const auto X3 = 93.27191  + 483202.017538*jce - 0.0036825*jce2 + jce3/327270;
    const auto X4 = 125.04452 -   1934.136261*jce + 0.0020708*jce2 + jce3/450000;

    const auto dksi_00 = (-171996 - 174.2*jce)*sin(d2r(                        X4  ));
    const auto dksi_01 = ( -13187 -   1.6*jce)*sin(d2r(-X0*2            +X3*2 +X4*2));
    const auto dksi_02 = (  -2274 -   0.2*jce)*sin(d2r(                  X3*2 +X4*2));
    const auto dksi_03 = (   2062 +   0.2*jce)*sin(d2r(                        X4*2));
    const auto dksi_04 = (   1426 -   3.4*jce)*sin(d2r(       X1                   ));
    const auto dksi_05 = (    712 +   0.1*jce)*sin(d2r(            X2              ));
    const auto dksi_06 = (   -517 +   1.2*jce)*sin(d2r(-X0*2 +X1        +X3*2 +X4*2));
    const auto dksi_07 = (   -386 -   0.4*jce)*sin(d2r(                  X3*2 +X4  ));
    const auto dksi_08 =     -301             *sin(d2r(            X2   +X3*2 +X4*2));
    const auto dksi_09 = (    217 -   0.5*jce)*sin(d2r(-X0*2 -X1        +X3*2 +X4*2));
    const auto dksi_10 =     -158             *sin(d2r(-X0*2      +X2              ));
    const auto dksi_11 = (    129 +   0.1*jce)*sin(d2r(-X0*2            +X3*2 +X4  ));
    const auto dksi_12 =      123             *sin(d2r(           -X2   +X3*2 +X4*2));
    const auto dksi_13 =       63             *sin(d2r( X0*2                       ));
    const auto dksi_14 = (     63 +   0.1*jce)*sin(d2r(            X2         +X4  ));
    const auto dksi_15 =      -59             *sin(d2r( X0*2      -X2   +X3*2 +X4*2));
    const auto dksi_16 = (    -58 -   0.1*jce)*sin(d2r(           -X2         +X4  ));
    const auto dksi_17 =      -51             *sin(d2r(            X2   +X3*2 +X4  ));
    const auto dksi_18 =       48             *sin(d2r(-X0*2      +X2*2            ));
    const auto dksi_19 =       46             *sin(d2r(           -X2*2 +X3*2 +X4  ));
    const auto dksi_20 =      -38             *sin(d2r( X0*2            +X3*2 +X4*2));
    const auto dksi_21 =      -31             *sin(d2r(            X2*2 +X3*2 +X4*2));
    const auto dksi_22 =       29             *sin(d2r(            X2*2            ));
    const auto dksi_23 =       29             *sin(d2r(-X0*2      +X2   +X3*2 +X4*2));
    const auto dksi_24 =       26             *sin(d2r(                  X3*2      ));
    const auto dksi_25 =      -22             *sin(d2r(-X0*2            +X3*2      ));
    const auto dksi_26 =       21             *sin(d2r(           -X2   +X3*2 +X4  ));
    const auto dksi_27 = (     17 -   0.1*jce)*sin(d2r(       X1*2                 ));
    const auto dksi_28 =       16             *sin(d2r( X0*2      -X2         +X4  ));
    const auto dksi_29 = (    -16 +   0.1*jce)*sin(d2r(-X0*2 +X1*2      +X3*2 +X4*2));
    const auto dksi_30 =      -15             *sin(d2r(       X1              +X4  ));
    const auto dksi_31 =      -13             *sin(d2r(-X0*2      +X2         +X4  ));
    const auto dksi_32 =      -12             *sin(d2r(      -X1              +X4  ));
    const auto dksi_33 =       11             *sin(d2r(            X2*2 -X3*2      ));
    const auto dksi_34 =      -10             *sin(d2r( X0*2      -X2   +X3*2 +X4  ));
    const auto dksi_35 =       -8             *sin(d2r( X0*2      +X2   +X3*2 +X4*2));
    const auto dksi_36 =        7             *sin(d2r(       X1        +X3*2 +X4*2));
    const auto dksi_37 =       -7             *sin(d2r(-X0*2 +X1  +X2              ));
    const auto dksi_38 =       -7             *sin(d2r(      -X1        +X3*2 +X4*2));
    const auto dksi_39 =       -7             *sin(d2r( X0*2            +X3*2 +X4  ));
    const auto dksi_40 =        6             *sin(d2r( X0*2      +X2              ));
    const auto dksi_41 =        6             *sin(d2r(-X0*2      +X2*2 +X3*2 +X4*2));
    const auto dksi_42 =        6             *sin(d2r(-X0*2      +X2   +X3*2 +X4  ));
    const auto dksi_43 =       -6             *sin(d2r( X0*2      -X2*2       +X4  ));
    const auto dksi_44 =       -6             *sin(d2r( X0*2                  +X4  ));
    const auto dksi_45 =        5             *sin(d2r(      -X1  +X2              ));
    const auto dksi_46 =       -5             *sin(d2r(-X0*2 -X1        +X3*2 +X4  ));
    const auto dksi_47 =       -5             *sin(d2r(-X0*2                  +X4  ));
    const auto dksi_48 =       -5             *sin(d2r(            X2*2 +X3*2 +X4  ));
    const auto dksi_49 =        4             *sin(d2r(-X0*2      +X2*2       +X4  ));
    const auto dksi_50 =        4             *sin(d2r(-X0*2 +X1        +X3*2 +X4  ));
    const auto dksi_51 =        4             *sin(d2r(            X2   -X3*2      ));
    const auto dksi_52 =       -4             *sin(d2r(-X0        +X2              ));
    const auto dksi_53 =       -4             *sin(d2r(-X0*2 +X1                   ));
    const auto dksi_54 =       -4             *sin(d2r( X0                         ));
    const auto dksi_55 =        3             *sin(d2r(            X2   +X3*2      ));
    const auto dksi_56 =       -3             *sin(d2r(           -X2*2 +X3*2 +X4*2));
    const auto dksi_57 =       -3             *sin(d2r(-X0   -X1  +X2              ));
    const auto dksi_58 =       -3             *sin(d2r(       X1  +X2              ));
    const auto dksi_59 =       -3             *sin(d2r(      -X1  +X2   +X3*2 +X4*2));
    const auto dksi_60 =       -3             *sin(d2r( X0*2 -X1  -X2   +X3*2 +X4*2));
    const auto dksi_61 =       -3             *sin(d2r(            X2*3 +X3*2 +X4*2));
    const auto dksi_62 =       -3             *sin(d2r( X0*2 -X1        +X3*2 +X4*2));

    const auto deps_00 = (92025 + 8.9*jce)*cos(d2r(                          X4  ));
    const auto deps_01 = ( 5736 - 3.1*jce)*cos(d2r(-X0*2              +X3*2 +X4*2));
    const auto deps_02 = (  977 - 0.5*jce)*cos(d2r(                    X3*2 +X4*2));
    const auto deps_03 = ( -895 + 0.5*jce)*cos(d2r(                          X4*2));
    const auto deps_04 = (   54 - 0.1*jce)*cos(d2r(       X1                     ));
    const auto deps_05 =     -7           *cos(d2r(              X2              ));
    const auto deps_06 = (  224 - 0.6*jce)*cos(d2r(-X0*2 +X1          +X3*2 +X4*2));
    const auto deps_07 =    200           *cos(d2r(                    X3*2 +X4  ));
    const auto deps_08 = (  129 - 0.1*jce)*cos(d2r(              X2   +X3*2 +X4*2));
    const auto deps_09 = (  -95 + 0.3*jce)*cos(d2r(-X0*2 -X1          +X3*2 +X4*2));
    const auto deps_10 =    -70           *cos(d2r(-X0*2              +X3*2 +X4  ));
    const auto deps_11 =    -53           *cos(d2r(             -X2   +X3*2 +X4*2));
    const auto deps_12 =    -33           *cos(d2r(              X2         +X4  ));
    const auto deps_13 =     26           *cos(d2r( X0*2        -X2   +X3*2 +X4*2));
    const auto deps_14 =     32           *cos(d2r(             -X2         +X4  ));
    const auto deps_15 =     27           *cos(d2r(              X2   +X3*2 +X4  ));
    const auto deps_16 =    -24           *cos(d2r(             -X2*2 +X3*2 +X4  ));
    const auto deps_17 =     16           *cos(d2r( X0*2              +X3*2 +X4*2));
    const auto deps_18 =     13           *cos(d2r(              X2*2 +X3*2 +X4*2));
    const auto deps_19 =    -12           *cos(d2r(-X0*2        +X2   +X3*2 +X4*2));
    const auto deps_20 =    -10           *cos(d2r(             -X2   +X3*2 +X4  ));
    const auto deps_21 =     -8           *cos(d2r( X0*2        -X2         +X4  ));
    const auto deps_22 =      7           *cos(d2r(-X0*2 +X1*2        +X3*2 +X4*2));
    const auto deps_23 =      9           *cos(d2r(       X1                +X4  ));
    const auto deps_24 =      7           *cos(d2r(-X0*2        +X2         +X4  ));
    const auto deps_25 =      6           *cos(d2r(      -X1                +X4  ));
    const auto deps_26 =      5           *cos(d2r( X0*2        -X2   +X3*2 +X4  ));
    const auto deps_27 =      3           *cos(d2r( X0*2        +X2   +X3*2 +X4*2));
    const auto deps_28 =     -3           *cos(d2r(       X1          +X3*2 +X4*2));
    const auto deps_29 =      3           *cos(d2r(      -X1          +X3*2 +X4*2));
    const auto deps_30 =      3           *cos(d2r( X0*2              +X3*2 +X4  ));
    const auto deps_31 =     -3           *cos(d2r(-X0*2        +X2*2 +X3*2 +X4*2));
    const auto deps_32 =     -3           *cos(d2r(-X0*2        +X2   +X3*2 +X4  ));
    const auto deps_33 =      3           *cos(d2r( X0*2        -X2*2       +X4  ));
    const auto deps_34 =      3           *cos(d2r( X0*2                    +X4  ));
    const auto deps_35 =      3           *cos(d2r(-X0*2 -X1          +X3*2 +X4  ));
    const auto deps_36 =      3           *cos(d2r(-X0*2                    +X4  ));
    const auto deps_37 =      3           *cos(d2r(              X2*2 +X3*2 +X4  ));

    const auto dksi = (dksi_00 + dksi_01 + dksi_02 + dksi_03 + dksi_04 + dksi_05 + dksi_06 + dksi_07 + dksi_08 + dksi_09
                     + dksi_10 + dksi_11 + dksi_12 + dksi_13 + dksi_14 + dksi_15 + dksi_16 + dksi_17 + dksi_18 + dksi_19
                     + dksi_20 + dksi_21 + dksi_22 + dksi_23 + dksi_24 + dksi_25 + dksi_26 + dksi_27 + dksi_28 + dksi_29
                     + dksi_30 + dksi_31 + dksi_32 + dksi_33 + dksi_34 + dksi_35 + dksi_36 + dksi_37 + dksi_38 + dksi_39
                     + dksi_40 + dksi_41 + dksi_42 + dksi_43 + dksi_44 + dksi_45 + dksi_46 + dksi_47 + dksi_48 + dksi_49
                     + dksi_50 + dksi_51 + dksi_52 + dksi_53 + dksi_54 + dksi_55 + dksi_56 + dksi_57 + dksi_58 + dksi_59
                     + dksi_60 + dksi_61 + dksi_62)/36000000.0;

    const auto deps = (deps_00 + deps_01 + deps_02 + deps_03 + deps_04 + deps_05 + deps_06 + deps_07 + deps_08 + deps_09
                     + deps_10 + deps_11 + deps_12 + deps_13 + deps_14 + deps_15 + deps_16 + deps_17 + deps_18 + deps_19
                     + deps_20 + deps_21 + deps_22 + deps_23 + deps_24 + deps_25 + deps_26 + deps_27 + deps_28 + deps_29
                     + deps_30 + deps_31 + deps_32 + deps_33 + deps_34 + deps_35 + deps_36 + deps_37)/36000000.0;

    return std::make_pair(dksi, deps);
}

double true_ecliptic_obliquity(const double julian_ephemeris_millennium,
                               const double obliquity_nutation) {
    const double U = julian_ephemeris_millennium/10;
    const double U2 = U*U;
    const double U3 = U2*U;
    const double U4 = U3*U;
    const double U5 = U4*U;
    const double U6 = U5*U;
    const double U7 = U6*U;
    const double U8 = U7*U;
    const double U9 = U8*U;
    const double U10 = U9*U;
    // Equation (24)
    const auto eps_0 = 84381.448 - 4680.93*U   - 1.55*U2 + 1999.25*U3
                        - 51.38*U4 -249.67*U5 - 39.05*U6 + 7.12*U7
                        + 27.87*U8 +  5.79*U9 + 2.45*U10;
    // Equation (25)
    return eps_0/3600.0 + obliquity_nutation;
}

double aberration_correction(const double R) {
    // Equation (26)
    return -20.4898/(3600.0*R);
}

double apparent_sun_longitude(const double geoc_long, const double nut_long,  const double ab_corr) {
    // Equation (27)
    return geoc_long + nut_long + ab_corr;
}
double apparent_Greenwich_sidereal_time(const double julian_day,
                                              const double julian_ephemeris_centry,
                                              const double nut_long,
                                              const double true_ecliptic_obliquity) {
    // Equation (28)
    const double nu0 = 280.46061837 + 360.98564736629*(julian_day - 2451545) + 0.000387933*pow(julian_ephemeris_centry, 2) -
                     pow(julian_ephemeris_centry, 3)/38710000;

    // Equation (29) with [0, 360) range limit
    return limit_degrees(nu0) + nut_long*cos(d2r(true_ecliptic_obliquity));
}

double geocentric_sun_right_ascension(const double app_sun_longitude,
                           const double true_ecliptic_obliq,
                           const double geocentric_lat) {
    // Equation (30)
    const auto lambda_rad = d2r(app_sun_longitude);
    const auto eps_rad = d2r(true_ecliptic_obliq);
    const auto beta_rad = d2r(geocentric_lat);
    const auto alpha = atan2(sin(lambda_rad)*cos(eps_rad) - tan(beta_rad)*sin(eps_rad), cos(lambda_rad));
    return limit_degrees(r2d(alpha));
}

double geocentric_sun_declination(const double app_sun_longitude,
                                  const double true_ecliptic_obliq,
                                  const double geocentric_lat) {
    // Equation (31)
    const auto lambda_rad = d2r(app_sun_longitude);
    const auto eps_rad = d2r(true_ecliptic_obliq);
    const auto beta_rad = d2r(geocentric_lat);
    const auto delta = asin(sin(beta_rad)*cos(eps_rad) + cos(beta_rad)*sin(eps_rad)*sin(lambda_rad));
    return r2d(delta);
}

double observer_local_hour_angle(const double app_Greenwich_sidereal_time,
                        const double geocentric_sun_right_asc,
                        const double geographical_lon) {
    // Equation (32)
    const double H = app_Greenwich_sidereal_time + geographical_lon - geocentric_sun_right_asc;
    return limit_degrees(H);
}

auto topocentric_values(const double R,
                        const double geographical_latitude,
                        const double masl,
                        const double geocentric_sun_right_asc,
                        const double geocentric_sun_decl,
                        const double obs_local_hour_angle) {

    // Convert to rad
    const double lat_rad = d2r(geographical_latitude);
    const double hour_angle_rad = d2r(obs_local_hour_angle);
    const double sun_decl_rad = d2r(geocentric_sun_decl);

    // Equation (33)
    const double xi = 8.794/(3600.0*R);
    const double xi_rad = d2r(xi);

    // Equation (34)
    const double u = atan(0.99664719*tan(lat_rad));
    // Equation (35)
    const double x = cos(u) + masl*cos(lat_rad)/6378140.0;
    // Equation (36)
    const double y = 0.99664719*sin(u) + masl*sin(lat_rad)/6378140.0;
    // Equation (37)
    const double dalpha_rad = atan2(                     - x*sin(xi_rad)*sin(hour_angle_rad),
                                       cos(sun_decl_rad) - x*sin(xi_rad)*cos(hour_angle_rad));
    const double dalpha = r2d(dalpha_rad);
    // Equation (38)
    const double topo_sun_right_asc = geocentric_sun_right_asc + dalpha;
    // Equation (39)
    const double topo_sun_declination = r2d(atan2(sin(sun_decl_rad) - y*sin(xi_rad)*cos(dalpha_rad),
                                              cos(sun_decl_rad) - x*sin(xi_rad)*cos(hour_angle_rad)));
    // Equation (40)
    const double topo_local_hour_angle = obs_local_hour_angle - dalpha;

    return std::make_tuple(topo_sun_right_asc , topo_sun_declination, topo_local_hour_angle);
}

auto uncorrected_topocentric_elevation_angle(const double geographical_latitude,
                                             const double topo_sun_declination,
                                             const double topo_local_hour_angle) {
    // Equation (41)
    const double lat_rad = d2r(geographical_latitude);
    const double sun_decl_rad = d2r(topo_sun_declination);
    const double hour_angle_rad = d2r(topo_local_hour_angle);
    const double eps_0 = r2d(asin(sin(lat_rad)*sin(sun_decl_rad) +
                                  cos(lat_rad)*cos(sun_decl_rad)*cos(hour_angle_rad)));
    return eps_0;
}

auto topocentric_zenith_angle(const double eps_0,
                              const double P,
                              const double T) {
    // Equation (42)
    const double deps = P/1010*283/(273 + T)*1.02/(60.0*tan(d2r(eps_0 + 10.3/(eps_0 + 5.11))));
    // Equation (43)
    const double topocentric_elevation_angle = eps_0 + deps;
    // Equation (44)
    const double topocentric_zenith_angle = 90 - topocentric_elevation_angle;
    return std::make_pair(topocentric_elevation_angle, topocentric_zenith_angle);
}

auto topocentric_astronomers_and_azimuth_angle(const double topo_local_hour_angle,
                                               const double geographical_latitude,
                                               const double topo_sun_declination) {
    // Equation (45)
    const double h_rad = d2r(topo_local_hour_angle);
    const double lat_rad = d2r(geographical_latitude);
    const double sun_decl_rad = d2r(topo_sun_declination);
    const double Gamma = limit_degrees(r2d(atan2(sin(h_rad),
                                    cos(h_rad)*sin(lat_rad) - tan(sun_decl_rad)*cos(lat_rad))));
    // Equation (46)
    const double Phi = limit_degrees(Gamma + 180);
    return std::make_tuple(Gamma, Phi);
}

double indicence_angle(const double topo_zenith_angle, const double surface_slope, const double surface_azimuth, const double astronomers_azimuth_angle) {
    // Equation (47)
    const double zenit_rad = d2r(topo_zenith_angle);
    const double slope_rad = d2r(surface_slope);
    const double astro_azimuth_rad = d2r(astronomers_azimuth_angle);
    const double surf_azimuth_rad = d2r(surface_azimuth);
    const double I = r2d(acos(cos(zenit_rad)*cos(slope_rad)
                            + sin(slope_rad)*sin(zenit_rad)*cos(astro_azimuth_rad - surf_azimuth_rad)));
    return I;
}

template<typename collector_t>
auto solar_position(const double julian_day,
                    const double DT,
                    const double geographic_latitude,
                    const double geographic_longitude,
                    const double masl,
                    collector_t collector) {
    // Note that the args should be in UT, and that |UT - UTC| < 1.0

    //const auto DT = Delta_T(year);
    const auto julian_ephemeris_day = jde(julian_day, DT);
    const auto julian_ephemeris_centry = jce(julian_ephemeris_day);
    const auto julian_ephemeris_millennium = jme(julian_ephemeris_centry);
    const auto L = heliocentric_longitude(julian_ephemeris_millennium);
    const auto B = heliocentric_latitude(julian_ephemeris_millennium);
    const auto R = heliocentric_radius_vector(julian_ephemeris_millennium);
    const auto beta = geocentric_latitude(B);
    const auto Theta = geocentric_longitude(L);


    const auto [longitude_nutation, obliquity_nutation] = nutation(julian_ephemeris_centry);
    const auto epsilon = true_ecliptic_obliquity(julian_ephemeris_millennium, obliquity_nutation);
    const auto delta_tau = aberration_correction(R);
    const auto lambda = apparent_sun_longitude(Theta, longitude_nutation, delta_tau);
    const auto nu = apparent_Greenwich_sidereal_time(julian_day,
                                                     julian_ephemeris_centry,
                                                     longitude_nutation,
                                                     epsilon);
    const auto alpha = geocentric_sun_right_ascension(lambda, epsilon, beta);
    const auto delta = geocentric_sun_declination(lambda, epsilon, beta);
    const auto H = observer_local_hour_angle(nu, alpha, geographic_longitude);
    const auto [alpha_mark, delta_mark, H_mark]  = topocentric_values(R,
                                                                      geographic_latitude,
                                                                      masl,
                                                                      alpha,
                                                                      delta,
                                                                      H);
    const double e0 = uncorrected_topocentric_elevation_angle(geographic_latitude,
                                                              delta_mark,
                                                              H_mark);
    const auto [Gamma, Phi] = topocentric_astronomers_and_azimuth_angle(H_mark,
                                                                        geographic_latitude,
                                                                        delta_mark);
    return collector(e0, Gamma, Phi, alpha_mark, delta_mark, H_mark);
}

auto corrected_solar_elevation(const double e0, const double P, const double T) {
    return topocentric_zenith_angle(e0, P, T);
}

template<typename collector_t, typename dt_calc_t>
auto calendar_solar_position(unsigned int year,
                             unsigned int month,
                             double day,
                             const double geographic_latitude,
                             const double geographic_longitude,
                             const double masl,
                             collector_t collector,
                             dt_calc_t cal_calc) {
    const double DT = cal_calc(year, month, day);
    const double julian_day = jd_from_cal(year, month, day);
    return solar_position(julian_day,
                          DT,
                          geographic_latitude,
                          geographic_longitude,
                          masl, collector);
}

template<typename collector_t, typename dt_calc_t>
auto time_point_solar_position(const std::chrono::system_clock::time_point time_point,
                                   const double geographic_latitude,
                                   const double geographic_longitude,
                                   const double masl,
                                   collector_t collector,
                                   dt_calc_t time_point_calc) {

    const auto DT = time_point_calc(time_point);
    const auto julian_day = jd_from_clock(time_point);
    return solar_position(julian_day,
                          DT,
                          geographic_latitude,
                          geographic_longitude,
                          masl, collector);
}

}
