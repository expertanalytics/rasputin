# Rehearsal: Swedish catchments as Shyft cells through the MCP server

Status: **rehearsal note; its fixes are designed in
`docs/increments/106-wednesday-cells.md`.** `@architect`, 2026-10-26, master
`098e1c80`, for Ola's Wednesday-morning run. Checked: the MCP tool list
(`tin_engine.mcp.schema.derive` over the CLI app), and the same argv run by
the CLI on Kultsjön (an SVAR unit in the Opal folder,
`~/rasputin-output/opal-se2-2026-10-21/`, read only). The server builds each
tool from the CLI's own options, so a CLI run is the MCP call's run. GeoPackage
compatibility with the client's project is out of scope.

**What is there.** The `mesh` tool has `triangles`, `max_triangles`,
`land_cover`, `land_cover_legend`, `reservoirs`, `domain`, `out_crs` and the
`features_*` options; `network_catchment` joins SVAR units. With
`--land-cover clc2018` and a `.gpkg` output, the file has the `cells` layer:
the Opal Kultsjön outline (867 vertices) gave 2,438 cells in 11 s, and
Ljungan (12,839 km²) gave 9,985 cells in 20 s. On both, the five Shyft
fractions summed to 1 within 3.6e-8.

**What would block or surprise Ola**

1. **`network_catchment` does not simplify the boundary.** Kultsjön's joined
   outline has 19,556 vertices, and `mesh` then refuses any budget under
   **34,451 triangles** (0.05 km² cells). The Opal outlines were simplified
   by the delivery's own assembly script, outside rasputin (a shapely
   `simplify(tolerance, preserve_topology=True)`); the same unit simplified that way allows 2,254. rasputin's
   area-keeping reduction (`_core.reduce_ring`) is reached only from
   `catchment` and `station_catchments` (`--outline-tolerance`, terrain
   route). Fix: `--outline-tolerance` on `network-catchment`, calling the
   same reduction. **About 40 lines, no C++, short path.** Workaround on the
   day: give `mesh` the Opal outlines.
2. **One `catchment_id` per run.** Every cell gets 1. Sub-catchments are
   one run, and one file, each. Tagging cells after the fact is no fix: on
   Kultsjön's 75 SVAR units, 18.9 % of the area lands in the wrong unit.
   With the unit borders as constraints, no cell straddles. Fix: a
   `--catchments` option. **About 120 lines**; not by Wednesday.
3. **Reservoir share is 0 in Sweden.** No Swedish reservoir polygons are
   registered. SMHI's dams join SMHI's lakes on the water-body id, which
   finds 82.73 km² of reservoir in Kultsjön, but rasputin cannot do that
   join today. **About 50 lines, short path.**
4. **`network_catchment --network` takes a file path, not a source name.**
   `smhi-svar2022-subcatchments` is refused ("use .gpkg, .geojson or
   .json"). The full path to `SVAR2022_delavrinningsomraden.gpkg` works.
   Under 10 lines to fix; the workaround is the path.

Question: should the Wednesday run start from the Opal outlines (already
simplified, seven catchments), or from `network_catchment` per forecast
point? Default: the Opal outlines, because of item 1.
