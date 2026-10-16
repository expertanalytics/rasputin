import pstats, sys
p = pstats.Stats(sys.argv[1], stream=sys.stdout); p.strip_dirs(); p.sort_stats("cumulative")
p.print_stats(r"transform|pyproj|_transformer|prepared\.py|prepare\)|base\.py:815|contains_properly|is_simple")
for fn in (r"_take\)", r"feature_input\.py:389", r"feature_input\.py:402", r"_outline_lines\)", r"base\.py:815"):
    print("### callees of", fn); p.print_callees(fn)
print("### callers of covers"); p.print_callers(r"base\.py:815")
