# Required reading

Invoke the Skill tool for `modern-cxx`, `computational-geometry`,
`python-development` and `geospatial-data-formats` — whichever touch the task —
before writing, judging or planning any code.

This is stated as a step rather than left to the `skills:` frontmatter key
because a dispatch test showed the key does not reliably preload them, and
subagents do not inherit skills from the caller. A session that skips it
re-derives decisions the project has already written down: the GeoTIFF decode
siting was escalated to @architect as an open question while the answer was
already in `geospatial-data-formats/SKILL.md`.

Also read `docs/increments/README.md` — the increment protocol, and the cost
constraints on a round — plus the increment file for whatever you are working
on. Designs live on disk precisely so they are not re-derived from a prompt.
