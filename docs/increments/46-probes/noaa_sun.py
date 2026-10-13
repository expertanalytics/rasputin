# ruff: noqa: E501 -- a one-off design probe; its formulas and report lines are long
"""NOAA's solar calculator (the spreadsheet's Meeus-based formulas, about 0.01 deg), standing in for TR1.

Elevation without refraction; azimuth true, from north, clockwise; R in AU.
Shared by probe_spring.py and probe_sun.py.
"""

import datetime as dt
import math


def sun(tt: dt.datetime, lat_deg: float, lon_deg: float) -> tuple[float, float, float]:
    """(elevation deg, azimuth deg, distance AU) at the naive-UTC time `tt`."""
    jd = 2440587.5 + tt.replace(tzinfo=dt.UTC).timestamp() / 86400
    jc = (jd - 2451545) / 36525
    l0 = (280.46646 + jc * (36000.76983 + jc * 0.0003032)) % 360
    m = 357.52911 + jc * (35999.05029 - 0.0001537 * jc)
    e = 0.016708634 - jc * (0.000042037 + 0.0000001267 * jc)
    mr = math.radians(m)
    c = (
        math.sin(mr) * (1.914602 - jc * (0.004817 + 0.000014 * jc))
        + math.sin(2 * mr) * (0.019993 - 0.000101 * jc)
        + math.sin(3 * mr) * 0.000289
    )
    v = math.radians(m + c)
    r = 1.000001018 * (1 - e * e) / (1 + e * math.cos(v))
    om = math.radians(125.04 - 1934.136 * jc)
    lam = math.radians(l0 + c - 0.00569 - 0.00478 * math.sin(om))
    eps0 = 23 + (26 + (21.448 - jc * (46.815 + jc * (0.00059 - jc * 0.001813))) / 60) / 60
    eps = math.radians(eps0 + 0.00256 * math.cos(om))
    d = math.asin(math.sin(eps) * math.sin(lam))
    y = math.tan(eps / 2) ** 2
    l0r = math.radians(l0)
    eot = 4 * math.degrees(
        y * math.sin(2 * l0r)
        - 2 * e * math.sin(mr)
        + 4 * e * y * math.sin(mr) * math.cos(2 * l0r)
        - 0.5 * y * y * math.sin(4 * l0r)
        - 1.25 * e * e * math.sin(2 * mr)
    )
    minutes = tt.hour * 60 + tt.minute + tt.second / 60 + tt.microsecond / 6e7
    ha = math.radians(((minutes + eot + 4 * lon_deg) % 1440) / 4 - 180)
    lat = math.radians(lat_deg)
    el = math.asin(math.sin(lat) * math.sin(d) + math.cos(lat) * math.cos(d) * math.cos(ha))
    az = (
        math.atan2(math.sin(ha), math.cos(ha) * math.sin(lat) - math.tan(d) * math.cos(lat))
        + math.pi
    )
    return math.degrees(el), math.degrees(az) % 360, r
