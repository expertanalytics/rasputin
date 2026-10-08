"""Design probe for increment 33: the field's query cost on the corridor, estimated.

Arithmetic only, no data read. Inputs are the probes' own figures (README items
1, 3 and 4) and one benchmark record; every cost constant is an assumption
stated beside it. Prints the share of final triangles per distance band, the
mean query cost per scanned triangle, and the corridor's CPU and wall seconds.
"""

import math

# README item 3: uniform meshes of the 267 km2 Geilo-Ål section.
SECTION_KM2 = 267.0
TOLS = (20.0, 10.0, 5.0, 2.0, 1.0)
COUNTS = (5967, 23373, 82286, 369284, 970037)

# README item 1: band areas (km2, cumulative) of the three lines' 5 km corridor.
BAND_EDGES_M = (0, 100, 500, 1000, 2000, 3000)
BAND_CUM_KM2 = (0, 95, 478, 955, 1897, 2828)
CORRIDOR_KM2 = 4689.0

# The ramp of section 3, Ola's case.
NEAR, FAR, END = 1.0, 20.0, 3000.0

# Assumed cost per scanned triangle, microseconds, (low, high), from the
# doubling search over BroadPhase's buckets: each occupied bucket met costs
# about 45 segments x 15 ns = 0.7 us; a box 2g across meets 1 to 4 occupied
# buckets up to g = 750 m, 2 to 4 at 1.5 km, 4 to 6 at 3 km; the doubling
# sums the queries (g = E/16, E/8, E/4, E/2, E, E + margin).
COST_US = {
    (0, 100): (0.7, 2.7),
    (100, 500): (1.4, 5.4),
    (500, 1000): (2.1, 8.1),
    (1000, 2000): (3.5, 10.8),
    (2000, 3000): (5.0, 15.0),
    "beyond E": (5.0, 15.0),
}

# Scans per final triangle: triangles created (3 per insert, 2 per flip) over
# final triangles (2 per insert); the bench's 1 m tile, threads 0 and 1,
# docs/benchmarks/2026-10-04/23b-fix-base-r3/run.json: 219 837 inserted,
# 445 675 flips.
INSERTED, FLIPS = 219837, 445675
SCANS = (3 * INSERTED + 2 * FLIPS) / (2 * INSERTED)

FINAL_TRIANGLES = 1_372_967  # README item 4, corr_bands.py's linear ramp
THREADS = 10


def density(tol: float) -> float:
    """Triangles per km2 at a uniform tolerance, log-log between the section's runs."""
    lt = math.log(tol)
    for hi, lo, c_hi, c_lo in zip(TOLS, TOLS[1:], COUNTS, COUNTS[1:], strict=False):
        if lo <= tol <= hi:
            f = (math.log(hi) - lt) / (math.log(hi) - math.log(lo))
            return math.exp(
                math.log(c_hi / SECTION_KM2) * (1 - f) + math.log(c_lo / SECTION_KM2) * f
            )
    return COUNTS[-1] / SECTION_KM2 if tol < 1 else COUNTS[0] / SECTION_KM2


def main() -> None:
    triangles: dict[object, float] = {}
    for i in range(len(BAND_EDGES_M) - 1):
        area = BAND_CUM_KM2[i + 1] - BAND_CUM_KM2[i]
        mid = (BAND_EDGES_M[i] + BAND_EDGES_M[i + 1]) / 2
        triangles[(BAND_EDGES_M[i], BAND_EDGES_M[i + 1])] = area * density(
            NEAR + (FAR - NEAR) * mid / END
        )
    triangles["beyond E"] = (CORRIDOR_KM2 - BAND_CUM_KM2[-1]) * density(FAR)
    total = sum(triangles.values())
    low = high = 0.0
    for band, count in triangles.items():
        share = count / total
        low += share * COST_US[band][0]
        high += share * COST_US[band][1]
        print(f"band {band}: share {share:.3f}, cost us {COST_US[band]}")
    print(f"scans per final triangle {SCANS:.2f}")
    print(f"mean cost per scanned triangle us {low:.1f} to {high:.1f}")
    cpu = [FINAL_TRIANGLES * SCANS * c / 1e6 for c in (low, high)]
    wall = [c / THREADS for c in cpu]
    print(f"corridor CPU s {cpu[0]:.1f} to {cpu[1]:.1f}")
    print(f"corridor wall s on {THREADS} threads {wall[0]:.1f} to {wall[1]:.1f}")


if __name__ == "__main__":
    main()
