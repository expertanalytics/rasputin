"""Increment 34 design probe: the estimates of sections 6 and 10, computed
from the measured figures (README items 2, 3 and 5) so that none is done by hand.

1. rasputin's triangle count for a rule, from the simulation: the rule's
   position between the simulation's uniform-F and uniform-N counts, applied
   between rasputin's own uniform counts on the same window (the simulation
   overcounts by 28 to 57 %, so only its relative position is used).
2. The whole Romsdalen tile (6901_3): its steep shares, and NumPy Horn's cost
   per node, the bound section 10 uses for the C++ steepness pass.
"""

import time

import numpy as np
import tifffile
from slope_stats import DTM, horn

SIM = {  # README item 2, word for word from greedy_sim.py's output
    ("romsdalen", "step 30"): {
        "uniform-F": 56193,
        "uniform-N": 412641,
        "node": 278199,
        "tri": 321685,
        "tri+halo": 344821,
        "normal": 271457,
    },
    ("romsdalen", "ramp 25-35"): {"node": 253808, "tri": 302249, "tri+halo": 325719},
    ("geilo-al", "step 30"): {
        "uniform-F": 16460,
        "uniform-N": 197623,
        "node": 42174,
        "tri": 59364,
        "tri+halo": 73792,
        "normal": 39019,
    },
}
RASPUTIN = {
    "romsdalen": (38935, 321985),
    "geilo-al": (10502, 143783),
}  # README item 3: F = 10, N = 2


def main() -> None:
    for (case, ramp), rules in SIM.items():
        f_sim, n_sim = SIM[(case, "step 30")]["uniform-F"], SIM[(case, "step 30")]["uniform-N"]
        f_r, n_r = RASPUTIN[case]
        for rule, count in rules.items():
            if rule.startswith("uniform"):
                continue
            pos = (count - f_sim) / (n_sim - f_sim)
            print(
                f"{case:9s} {ramp:10s} {rule:8s} position {pos:.3f}  "
                f"rasputin estimate {f_r + pos * (n_r - f_r):9.0f}"
            )
    z = tifffile.imread(f"{DTM}/6901_3_10m_z33.tif").astype(np.float64)
    t0 = time.perf_counter()
    s = horn(z)
    ns = (time.perf_counter() - t0) * 1e9 / z.size
    print(
        f"tile 6901_3: {z.shape[0]} x {z.shape[1]} nodes; NumPy Horn {ns:.1f} ns/node, one thread"
    )
    print("  shares: " + " ".join(f">={a}: {100 * np.mean(s >= a):.1f}%" for a in (25, 30, 35)))
    # README item 3: rasputin on the whole tile, F = 10 and N = 2. Its steep share
    # lies between the two windows', so its position is taken between theirs:
    # Geilo-Al's step-30 node rule and Romsdalen's ramp-25-35 node rule.
    f_t, n_t = 615459, 5916344
    lo = (SIM[("geilo-al", "step 30")]["node"] - 16460) / (197623 - 16460)
    hi = (SIM[("romsdalen", "ramp 25-35")]["node"] - 56193) / (412641 - 56193)
    print(
        f"  node rule, ramp 25-35: {f_t + lo * (n_t - f_t):.0f} to "
        f"{f_t + hi * (n_t - f_t):.0f} triangles"
    )


if __name__ == "__main__":
    main()
