import pstats, sys, io
p = pstats.Stats(sys.argv[1], stream=sys.stdout); p.strip_dirs(); p.sort_stats("cumulative")
pat = r"cli\.py:592\(mesh\)|_open_features|_take\)|_clean\)|_dem_mesh|open_dem\)|set_operations\.py.*\((intersection|union_all|coverage_union_all|symmetric_difference)\)|coverage_clean|make_valid|is_valid\)|is_simple|contains_properly|covers\)|_geometry\.py.*\(transform\)|strtree\.py:75|query\)|query_nearest|linestring\.py:45|creation\.py.*linestrings|measurement\.py.*distance\)|_core\.simplify_borders|get_rings|get_coordinates\)|get_parts\)|intersects\)|_untouched|_finite"
p.print_stats(pat)
for fn in (r"_take\)", r"feature_input\.py:389", r"feature_input\.py:402", r"feature_input\.py:448", r"_Outline|feature_input\.py:761"):
    print("### callees of", fn); p.print_callees(fn)
