"""pstats text dump: top cumulative, top tottime, and the callees of the clean-up functions."""
import pstats, sys
p = pstats.Stats(sys.argv[1], stream=sys.stdout)
p.strip_dirs()
print("### top 45 by cumulative time"); p.sort_stats("cumulative").print_stats(45)
print("### top 30 by own (tottime)"); p.sort_stats("tottime").print_stats(30)
for fn in ("_take", "_clean", "snap_to_outline", r"\bring\b", "simplify_borders", "_outline_lines", "_add", "_loops", "_place", "_join", "_cut", "pre_clip", "_runs", "_chains"):
    print(f"### callees of {fn}"); p.sort_stats("cumulative").print_callees(fn)
