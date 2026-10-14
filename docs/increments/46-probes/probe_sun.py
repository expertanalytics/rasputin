# ruff: noqa: E501 -- a one-off design probe; its formulas and report lines are long
"""Increment 46 design probe: NOAA's solar calculator at the synthetic domain's centre.

    .venv/bin/python docs/increments/46-probes/probe_sun.py > docs/increments/46-probes/probe_sun.out
noaa_sun.py, as probe_spring.py uses it; one-minute steps over each UTC day.
"""

import datetime as dt

from noaa_sun import sun

for day in (dt.datetime(2027, 3, 20), dt.datetime(2027, 4, 15), dt.datetime(2027, 6, 21)):
    best = max(
        ((*sun(day + dt.timedelta(minutes=m), 61.334, 15.0)[:2], m) for m in range(1440)),
        key=lambda e: e[0],
    )
    print(
        f"{day.date()}: max elevation {best[0]:.2f} deg at {best[2] // 60:02d}:{best[2] % 60:02d} UTC, azimuth there {best[1]:.2f} deg"
    )
