# Required reading

Invoke the Skill tool for `modern-cxx`, `computational-geometry`,
`python-development` and `geospatial-data-formats` — whichever touch the task —
before writing, judging or planning any code.

This is stated as a step rather than left to the `skills:` frontmatter key
because that key does not reliably preload them and subagents do not inherit
skills from the caller. Both are re-testable in about a minute, which is why the
instruction is an imperative rather than a declaration.

A session that skips it re-derives decisions the project has already written
down. Each `SKILL.md` is the settled answer to a question that otherwise comes
back as an open one — `geospatial-data-formats/SKILL.md` rules on where GeoTIFF
decoding sits, `computational-geometry/SKILL.md` on degeneracy policy.

When you write a justification into this project's governance, name the command,
file or commit a sceptical reader would use to check it. If you cannot, state
the mechanism rather than the incident: a mechanism someone can re-run is
stronger evidence than an event they can only be told about.

Also read `docs/increments/README.md` — the increment protocol, and the cost
constraints on a round — plus the increment file for whatever you are working
on. Designs live on disk precisely so they are not re-derived from a prompt.
