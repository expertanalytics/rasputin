# Sula: dense triangles along lake shores in a direct mesh run

Diagnosis probe, 2026-10-28, on master `aabac548`. No fix designed.

## Cause

A direct run takes the constraint lines from both sources as they are: the
CORINE polygon outlines (whose lake class 512 borders run along the lakes) and
the NVE shorelines. `--prefer lake=nve-lakes` replaces the lake *classes* in
the land cover, but not CORINE's lake *outlines* among the lines, so each lake
has two shorelines a few metres apart, crossing each other. The thin strips
between them cannot reach the start mesh's smallest angle, so the start-quality
step keeps splitting the two lines and adding points between them.

The two-stage run took its lines from a mesh whose classes already had the NVE
lakes put in, so there is one shoreline per lake.

## Measured

Inside the Sula domain: 18.6 km of CORINE outline runs within 80 m of an NVE
shoreline (6.5 km within 5 m, 9.3 km at 5 to 20 m, 2.8 km at 20 to 80 m), and
the two cross 269 times. The *band* below is where both run within 40 m of a
point (1.25 km²).

All runs at 20 m terrain tolerance and `--features-tolerance 20`; one change
from the direct command per row. Small = under 300 m².

| Run | Start angle | Triangles | Line splits in start quality | Small, in band | Small, elsewhere |
|---|---|---|---|---|---|
| direct, as recorded | 15 (default) | 6,429 | 923 | 1,898 | 125 |
| direct, `--start-min-angle 25` | 25 | 9,571 | 1,311 | 3,128 | 231 |
| direct, `--start-min-angle 0` | 0 | 2,706 | 0 | 220 | 57 |
| direct, `--start-quality-gain -1` (no line splits) | 15 | 3,867 | 0 | | |
| direct, `--no-constraint-feet` | 15 | 3,578 | 0 | | |
| direct, `--features-outline-snap 1` | 15 | 6,452 | 924 | | |
| NVE lakes only | 15 | 2,466 | 56 | 7 | 58 |
| NVE lakes only, angle 25 | 25 | 3,397 | 128 | | |
| CORINE only | 15 | 2,025 | 51 | | |
| two-stage, as recorded | 25 (default) | 4,943 | 269 | 289 | 336 |
| two-stage, angle 15 | 15 | 3,564 | 138 | 213 | 248 |
| two-stage, angle 0 | 0 | 2,372 | 0 | 77 | 149 |

Each source alone splits about 50 lines; together they split 923. At the same
15 degrees, direct makes 6,429 triangles against two-stage's 3,564. The angle
and the line splits multiply the effect but do not cause it: with both off
(angle 0) the band already holds 282 triangles under 15 degrees, against 94 in
the two-stage start mesh.

## What a fix would change

The direct run would build its constraint lines from the land cover after
`--prefer` has been applied, as the two-stage run does, so CORINE's lake
outline is dropped where an NVE shoreline replaces it.
