# Design probes for increment 33 (`@architect`, 2026-10-08)

Run by hand from one working directory holding Bane NOR's national
Banenettverk file, `Samferdsel_0000_Norge_25833_Banenettverk_GML.gml` (from
`https://nedlasting.geonorge.no/geonorge/Samferdsel/Banenettverk/GML/Samferdsel_0000_Norge_25833_Banenettverk_GML.zip`,
NLOD 1.0, dated 2025-03-06 in Geonorge's metadata; whoever fetches it again
sends a generic User-Agent, never Ola's email address), with the DTM10 tiles at
`../rasputin_data/DTM10_UTM33_20260925` (the scripts name the absolute path).
Python is the repository's `.venv`; `rasputin` is master at `4cd7e050`. The Mac
was on AC power (`pmset -g batt`: "Now drawing from 'AC Power'").

1. `probe_line.py`: the line's length, tunnels, corridor areas and DTM10
   coverage. Output, word for word (the first line of the Counter is the
   `(banestatus, medium)` count of Bergensbanen's links; `I` is in operation,
   `U` a tunnel):

   ```
   parse s 0.8
   Bergensbanen (status, medium) link counts Counter({('I', 'U'): 484, ('I', None): 482, ('I', 'B'): 246, ('I', 'L'): 68, ('I', 'T'): 60})
   ['Bergensbanen'] in operation: links 1340 km 743.2 tunnel km 168.2 vertices 51130 mean seg m 14.9
     simplify 1.0 m: vertices 10996 s 0.0
     5 km corridor km2 3591 s 10.6
     3 km band km2 2171
       band 100 m km2 75
       band 500 m km2 371
       band 1000 m km2 737
       band 2000 m km2 1458
     DTM10 tiles touching corridor 11 uncovered corridor km2 0.0
   ['Bergensbanen', 'Randsfjordbanen', 'Drammenbanen'] in operation: links 1499 km 948.7 tunnel km 193.8 vertices 69479 mean seg m 14.0
     simplify 1.0 m: vertices 13628 s 0.0
     5 km corridor km2 4689 s 10.0
     3 km band km2 2828
       band 100 m km2 95
       band 500 m km2 478
       band 1000 m km2 955
       band 2000 m km2 1897
     DTM10 tiles touching corridor 14 uncovered corridor km2 0.0
   ```

   "km" sums every link and is about twice the route length: the area of a
   50 m buffer of Bergensbanen's links (simplified by 1 m), divided by 100 m,
   is 371.7 km, so most of the route is drawn twice in the file (why, the
   probe did not look into). The band areas are of the union and count
   nothing twice.

2. `section.py`: the first case, Geilo to Ål: the links inside
   x 127 000-150 000, y 6 720 000-6 745 000 (EPSG:25833), buffered 5 km, cut
   to that x range and simplified by 20 m: 266.9 km², 53 vertices
   (`section_domain.geojson`, with a `crs` member naming EPSG:25833). It also
   writes the corridor of item 4 (`corridor_domain.geojson`).

3. Uniform meshes of the section, `rasputin mesh --dem <DTM10> --domain
   section_domain.geojson --tolerance T`, then `bands.py` (triangles by the
   distance from their centroid to the section's line):

   | T (m) | triangles | wall s |
   |---|---|---|
   | 20 | 5 967 | 0.83 |
   | 10 | 23 373 | 0.50 |
   | 5 | 82 286 | 0.55 |
   | 2 | 369 284 | 0.86 |
   | 1 | 970 037 | 1.53 |

   `bands.py`'s estimate of a mesh under the ramp (each band's density
   interpolated in log-log between the measured tolerances at the band's
   middle distance): linear from 1 m on the line to 20 m at 3 km, **53 154**;
   1 m held to 100 m, then linear, 74 228; 1 m out to 3 km, then 20 m,
   630 967.

4. The corridor, Hokksund to Bergen (the 4 196 km² largest part of the
   5 km buffer of the three lines above, which is not connected at Drammen to
   Hokksund because that stretch is Sørlandsbanen; 1 154 vertices, 3 holes),
   uniform:

   | T (m) | triangles | wall s | peak memory |
   |---|---|---|---|
   | 20 | 229 611 | 2.53 | 2.43 GB |
   | 1 | 21 669 768 | 36.95 | 11.75 GB |

   At 20 m, decode took 1.63 s and refine 0.47 s (`--stats`); at 1 m, refine
   17.89 s of which 12.94 s the serial split and flip phase. `corr_bands.py`
   (400 000 triangles sampled per mesh, seed 1, at 20, 10, 5, 2 and 1 m; line
   simplified by 1 m, 12 129 segments) estimates **1 372 967** triangles for
   the linear ramp from 0 to 3 km, 1 747 968 holding 1 m to 100 m, and
   3 115 974 holding it to 500 m.

5. `query_cost.py`: arithmetic only, no data read. The field's query cost
   on the corridor, from items 1, 3 and 4, the scans per final triangle of
   the bench's 1 m tile
   (`docs/benchmarks/2026-10-04/23b-fix-base-r3/run.json`: 219 837 inserted,
   445 675 flips; created triangles, 3 per insert and 2 per flip, over final
   triangles, 2 per insert), and a cost per scanned triangle per distance
   band computed from the rule the script states as an assumption (0.675 µs
   per occupied bucket; the buckets a query meets by its `g`; each band
   charged the queries its far edge needs). Output, word for word:

   ```
   band (0, 100): share 0.300, cost us 0.7 to 2.7
   band (100, 500): share 0.368, cost us 2.0 to 8.1
   band (500, 1000): share 0.146, cost us 3.4 to 10.8
   band (1000, 2000): share 0.096, cost us 6.1 to 14.8
   band (2000, 3000): share 0.037, cost us 6.1 to 14.8
   band beyond E: share 0.053, cost us 8.8 to 18.9
   scans per final triangle 3.53
   mean cost per scanned triangle us 2.7 to 8.3
   corridor CPU s 13.1 to 40.4
   corridor wall s on 10 threads 1.3 to 4.0
   uniform 1 m refine wall us per output triangle 0.83
   queries' wall us per output triangle 0.96 to 2.94
   ratio to uniform 1 m 2.2 to 4.6
   ```

   The shares come from a cruder density model than item 4's (the section's
   uniform densities applied to the corridor's band areas; it totals about
   0.78 million triangles, not 1.37 million), and are used only as shares.
   The laziness (section 4.4 of the design) skips the query for a triangle
   whose error is above `F`; the script does not, so it errs high there.

The estimates ignore the extra refinement that the per-triangle rule costs
where a large triangle reaches towards the line (section 3 of the design), so
they understate the count; by how much was not measured.
