# Harness step 2: citation pinning

Status: built on branch `worktree-citation-pinning` (red 7618d05, green
5e3d8a6); see *Review* at the end. Source: `docs/research/generic-harness.md`
§3.2, §6.3 (defect 5) and §7 (step 2).

**Goal.** After this step no citation that points by line into a rule file is
unpinned, in prose or in a code comment, and `tools/check_citations.py` fails
the build if one appears.

Examples in this file use placeholders (`<rev>`, `<n>`) wherever a literal
would be parsed as a citation. Every real citation here is pinned or resolves.

---

## 1. Prior art: legacy and literature

*Literature.* This is tooling; no novelty is claimed. Three established
conventions are combined:

- **Commit-pinned permalinks.** GitHub's "permanent link" (`blob/<sha>/path#L<a>-L<b>`)
  fixes a line reference to a commit so the text it names cannot move. The
  `path@<rev>:<line>` form below is the same idea in plain text, resolved with
  `git cat-file -p <rev>:<path>` instead of a URL.
- **Stable labels over positions.** LaTeX `\label`/`\ref`, Sphinx
  cross-reference labels and RFC section numbers cite a name, not a position.
  The ID form below does the same with the section labels the rule files
  already carry.
- **Link checkers** (lychee, markdown-link-check) check that a target and an
  anchor exist. They do not check that the quoted text still holds; neither
  does this tool. That limit stays, and the at-risk worklist stays with it.

*Legacy.* Nothing to carry across. The legacy tree has no citation tooling:

```
$ git grep -il citation legacy-archive -- legacy
legacy-archive:legacy/rasputin/reader.py
$ git grep -lE '\.md:[0-9]' legacy-archive -- legacy
(no output)
```

The one hit is the GeoTIFF key names (`GTCitationGeoKey` and siblings).

---

## 2. Citation forms

The tool knows five forms. Each is found on a single line of scanned text
(see §4). A citation that wraps across two lines is not seen.

Common pieces:

- `PATH` = `[\w./-]+\.(md|py|hpp|cpp|h|yaml|yml|toml|txt|cmake)`, preceded
  by start of text or a character outside `[\w/.-]`, as today. **A `PATH`
  that starts with `/` is not a citation** (absolute paths, and URLs, where
  the match would begin at `//`).
- `LINES` = `<a>` or `<a>-<b>`, decimal. The last line is `max(a, b)`, as
  today. A line of `0` is broken for every form.

### 2.1 Unpinned line citation, `PATH:LINES` (existing)

Resolution unchanged except for one refinement: a `PATH` that is not a file
at `REPO/PATH` is matched against tracked-or-untracked files whose path *ends
with the cited path components* (so `geospatial-data-formats/SKILL.md` names
one file, not four). One candidate resolves; several are *ambiguous* (broken);
none is *no such file* (broken). For a bare basename this is the same as
today. The search roots stay `SEARCH_ROOTS`.

**The `legacy/` case is kept, and becomes sugar**: an unpinned `PATH` starting
with `legacy/` is resolved exactly as `PATH@legacy-archive:LINES` (§2.2), never
against the working tree. Its messages keep naming the tag.

### 2.2 Pinned line citation, `PATH@REV:LINES` (new)

`REV` = `[\w./-]+`. It is one of two kinds, decided by its text:

| `REV` text | Kind | Accepted when |
|---|---|---|
| matches `[0-9a-f]{7,40}` | commit sha | `git rev-parse --verify --quiet <REV>^{commit}` succeeds (unique, present) |
| anything else | tag | `git rev-parse --verify --quiet refs/tags/<REV>^{commit}` succeeds |

Branch names, `HEAD`, and anything else that moves are **not pins**. A `REV`
that is not a sha and not a tag, but is `HEAD` or names a branch
(`refs/heads/<REV>` or `refs/remotes/*/<REV>` resolves), is broken with a
message saying it is not a pin. Tags with `/` in the name
(`archive/…`) are accepted.

`PATH` is taken literally as a repo-relative path **at `REV`**: no basename
search, no suffix match, no working-tree lookup. Content comes from
`git cat-file -p <REV>:<PATH>` (which peels an annotated tag); its line count
is `len(stdout.splitlines())`.

The pinned pattern is matched first; the unpinned pattern in §2.1 must not
also match inside a pinned citation's span.

### 2.3 ID citation, `PATH §ID` (new)

Syntax: `PATH` ending in `.md`, optionally closed by a backtick, then optional
spaces or tabs, then `§` **immediately** followed by

`ID` = `\d+(\.\d+)*[A-Z]?` | `[A-Z]{1,2}\d+` | `[A-Z]`, not followed by `\w`.

So `` `CLAUDE.md` §2 ``, `` `tester.md` §3C ``, `` `generic-harness.md` §3.2. ``
(the final period is punctuation) and `PRINCIPLES.md §E1` are ID citations.
`§ Constraint` (space after `§`) and "`CLAUDE.md` section 2" are not citations
and are ignored.

`PATH` resolves as in §2.1 (working tree). **Declared IDs.** Parse the target
file's ATX headings (`#` to `######`, then whitespace), skipping lines inside
fenced code blocks (between lines opening with three or more backticks or
tildes). Strip one optional leading `§`. The heading declares:

- a numeric label `\d+(\.\d+)*[A-Z]?` followed by `.`, whitespace or end of
  line (`## 1. Multi-Agent Ecosystem` → `1`; `### 3.2 Citation stability`
  → `3.2`; `### QW1.` is not numeric);
- an alphanumeric label `[A-Z]{1,2}\d+` followed by `.`, `:`, whitespace or
  end (`### E1 — …` → `E1`; `### QW1. …` → `QW1`);
- a letter label `[A-Z]` followed by `.` then whitespace or end
  (`### C. Data Source …` → `C`; `### A tie point` declares nothing). A letter
  label also declares `<N><L>`, where `N` is the numeric label of the letter
  heading's direct parent (the nearest preceding heading with fewer `#`), if
  that parent has one (`### C.` under `## 3.` → `3C`).

An ID citation **resolves** when `ID` is in the target's declared set.

### 2.4 Heading-name citation, `` `PATH`, *NAME* `` (new)

Syntax, exactly: a backticked `PATH` ending in `.md`, a comma, one space, and
`*NAME*` with single asterisks and no `*` inside `NAME`. Resolves when some
heading of the target (fence-aware, as in §2.3), with the `#`s and following
whitespace removed, equals `NAME`, or starts with `NAME` followed by `:`,
` [`, ` —` or ` (`. Case-sensitive. So `` `testing.md`, *Frameworks* ``
resolves against `## Frameworks [partly live]`, and
`` `docs/increments/README.md`, *Acceptance* `` against `## Acceptance: …`.

### 2.5 What is not a citation

Line references without a `PATH` (a bare `:14`, "line 220 of"), and any form
with a placeholder instead of digits (`<n>`, `N`). Examples in docs use
placeholders; there is no escape syntax.

---

## 3. Outcomes

One entry per citation, in the first category that applies, in this order:
**broken**, **unpinned**, **at risk**. Exit status is **1 if broken or
unpinned is non-empty**, else 0.

| Case | Category |
|---|---|
| Unpinned: no such file / ambiguous / last line past the end / line 0 | broken (unchanged) |
| Unpinned `legacy/`: tag or file missing at `legacy-archive`, or past the end | broken (unchanged messages) |
| Pinned: `REV` is neither a present sha, a tag, `HEAD` nor a branch, and the clone is not shallow | broken: "no commit or tag `<REV>`" |
| Pinned: `REV` is `HEAD`, or names a branch and is not a tag | broken: "not a pin" |
| Pinned: `REV` not found **and** `git rev-parse --is-shallow-repository` prints `true` | broken, message names the rev and says the clone is shallow and needs full history (`fetch-depth: 0`) |
| Pinned: sha prefix ambiguous | broken |
| Pinned: `PATH` absent at `REV` | broken: "no `<PATH>` at `<REV>`" |
| Pinned: last line past the end at `REV`, or line 0 | broken: "`<PATH>` has `<n>` lines at `<REV>`" |
| ID: `PATH` does not resolve | broken (as for unpinned) |
| ID: `ID` not declared in the target | broken: "no heading declares §`<ID>`" |
| Heading-name: `PATH` does not resolve, or no heading matches | broken |
| Unpinned line citation that resolves into a **governed file** (§3.1) | **unpinned** |
| Unpinned line citation that resolves into a file the branch edits | at risk (unchanged; exit unaffected) |
| Pinned, ID and heading-name citations that resolve | nothing: never at risk, never unpinned |
| Base not diffable | warning (unchanged): at-risk disabled |
| A scanned `.py` that `ast.parse` or `tokenize` rejects | warning: "`<file>`: not scanned, does not parse" |

CI is covered: the `governance` job in `.github/workflows/main.yaml` checks
out with `fetch-depth: 0`, so every commit on master's history and every tag is
present. The shallow-clone row exists for local and other runners.

Every entry names the citing location as `<source>:<line>` and quotes the
cited token. Section headers, in this order when non-empty:
`== at risk (<n>)`, `== unpinned (<n>)`, `== broken (<n>)`, `== warnings (<n>)`.
Tests assert on the exit status, the header and the presence of the citing
location and cited token in that section, not on the rest of the wording.

### 3.1 Governed files, and why unpinned is a failure

A resolved repo-relative path is **governed** when it is one of
`CLAUDE.md`, `docs/PRINCIPLES.md`, `.claude/REQUIRED-READING.md`,
`docs/increments/README.md`, `testing.md`, or starts with `.claude/agents/`,
`.claude/skills/`, `.claude/hooks/` or `tools/check_`.

This is `guard_governance.py`'s `GOVERNED` and `GOVERNED_PREFIXES` plus
`testing.md` and `.claude/skills/` (§6.3 defect 3 of the design). The list
lives as a constant in `check_citations.py` for this step; step 4 replaces
both copies with `profile.rule_files`. Settings globs are not citable
(`.json` is not a `PATH` suffix) and are omitted.

A failure, not the warning §3.2 of the design proposed: a warning is a
worklist, and the at-risk list already shows that worklists are not settled.
After the migration in §5 the tree has no such citation, so the failure costs
nothing until one is written, and then it names the fix.

---

## 4. Scan set and comment extraction

**Default sources** (no `--paths`): `git ls-files --cached --others
--exclude-standard`, keeping files that exist, whose suffix is in
`SCAN_SUFFIXES = (.md, .py, .h, .hpp, .cpp, .cmake)` or whose name is
`CMakeLists.txt`, and whose path does not start with `lib/` (vendored),
`legacy/` or `.claude/worktrees/`. This takes in the root `*.md` files,
`tests/`, `tools/`, `include/`, `src/`, `src_python/` and `bindings/`, which
today's default (`docs .claude CLAUDE.md`) misses.

**`--paths`** keeps its meaning: named files are scanned whatever their
suffix; named directories are walked (`rglob`), filtered by the rule above
and by `in_repo_proper`.

**What text is scanned, by type.** A citation is sought only in this text,
reported at its own line number:

| Type | Scanned text |
|---|---|
| `.md`, and any explicitly named file of another type | every line |
| `.py` | `tokenize` `COMMENT` tokens; docstrings (the first statement of a module, class, `def` or `async def`, when it is a string constant), every line from `lineno` to `end_lineno`. Other string literals are **not** scanned |
| `.h`, `.hpp`, `.cpp` | text after `//` to end of line, and inside `/* … */` (which may span lines). Outside comments, skip `"…"` (with `\` escapes), `'…'` and raw strings `R"d(…)d"` (which may span lines). A `'` between two alphanumeric characters is a digit separator (`0x9513'0000`), not a delimiter. An unterminated `'` or `"` ends at end of line |
| `.cmake`, `CMakeLists.txt` | text after a `#` that is outside a `"…"` argument (with `\` escapes), to end of line. Bracket comments `#[[…]]` are not recognised (none in the tree) |

This is the false-positive control: a citation inside a string literal, a
raw-string docstring in `bindings/`, or test-fixture data is not seen. The
tool's own docstring and its test module follow the same rules as every other
file; there is no exclusion list.

---

## 5. Migration of existing citations

`git grep -nE '(PRINCIPLES|testing|CLAUDE|README|REQUIRED-READING)\.md:[0-9]|guard_[a-z_]+\.py:[0-9]'`
(all suffixes) at master `7426b69` finds 13 outside `docs/research/`, plus
the design doc's own. None was converted by step 1. Pins below were read at
the rev given, each against the quotation in the citing sentence.

### 5.1 @tester, in the red commit (living comments, cited by heading)

| Citing line | Old target | New text |
|---|---|---|
| `tests/cpp/property/noding_generators.h@7426b69:5` | `testing.md@7426b69:220` | `` `testing.md`, *Frameworks* `` |
| `tests/cpp/property/prop_noding_no_crossings.cpp@7426b69:98` | `testing.md@7426b69:220` | `` `testing.md`, *Frameworks* `` |

Comment edits only. They stay correct under today's checker (it does not scan
`.cpp`) and resolve under the new one.

### 5.2 Doc conversion, its own commit after red and before green

Historical records, pinned. The old checker does not parse `@`, so this
commit keeps today's gate green; the green commit then checks every pin.

| Citing line | Old | New |
|---|---|---|
| `docs/increments/05-noder.md@7426b69:992` | `testing.md` line 301 | `testing.md@3245958:301` |
| `docs/increments/05b-noder-driver.md@7426b69:411` | `docs/increments/README.md` lines 69-71 | `docs/increments/README.md@3b78594:69-71` |
| `docs/retrospectives/2026-09-29-orchestrator-and-hooks-audit.md@7426b69:142` | `guard_governance.py` lines 24-25 | `.claude/hooks/guard_governance.py@180b9ac:24-25` |
| same file, `:146` | `.hooks/guard_governance.py` lines 56-59 (path wrong) | `.claude/hooks/guard_governance.py@180b9ac:56-59` |
| same file, `:164` | `guard_push.py` lines 20-30 | `.claude/hooks/guard_push.py@180b9ac:20-30` |
| same file, `:194` | `REQUIRED-READING.md` lines 49-54 | `.claude/REQUIRED-READING.md@180b9ac:49-54` |
| same file, `:208` | lines 82-98 | `.claude/REQUIRED-READING.md@180b9ac:82-98` |
| same file, `:261` | lines 118-126 | `.claude/REQUIRED-READING.md@180b9ac:118-126` |
| same file, `:271` | lines 100-126 | `.claude/REQUIRED-READING.md@180b9ac:100-126` |
| same file, `:302` | lines 68-70 | `.claude/REQUIRED-READING.md@180b9ac:68-70` |
| same file, `:320` | lines 110-112 | `.claude/REQUIRED-READING.md@180b9ac:110-112` |

`180b9ac` is the audit's last commit; every one of its nine targets reads
correctly there (checked with `git show 180b9ac:<path> | sed -n '<a>,<b>p'`).
The retrospective's line numbers themselves are unchanged.

**`docs/research/generic-harness.md`.** It is a snapshot, so it is pinned,
not reworded:

- every `PATH:LINES` into a governed file in §1.5 and §3.2 and in §6.3
  defects 1-3 becomes the full repo-relative path `@3245958`
  (`python3 tools/check_citations.py --paths docs/research` lists them once
  green; 76 parse today);
- the four "`file.md` line N" forms in §1.5 and every path-less continuation
  (`:14`, `:3`, `:8`, `:20-23`, `:283-288` and similar) are spelled out as
  full pinned citations at `3245958`;
- in §6.3 defect 5 and §7, a citation of the audit's hook lines is pinned at
  `180b9ac`, `testing.md` and `README.md` line citations at `3245958`, and the
  planted example becomes a placeholder (an unpinned `testing.md:<n>`);
- §7 step 5a cites section 5 of `CLAUDE.md` by ID, and that section does not
  exist yet, so the ID citation would be broken: it becomes "a new section 5
  of `CLAUDE.md`".

The converter reads each pinned line at its rev before writing the pin. A
citation whose quotation does not hold at `3245958` is pinned to the rev where
it does, and the commit message lists it.

**Acceptance of the migration.** On the branch after green:
`python3 tools/check_citations.py` exits 0 with no unpinned and no broken
section; and the §6.3 recipe —
`git grep -nE '(PRINCIPLES|testing|CLAUDE|README|REQUIRED-READING)\.md:[0-9]'`
over all suffixes — returns only pinned forms or placeholders.

---

## 6. Test cases for @tester

Black-box, as the existing suite: copy the tool into a scratch repository,
commit, run it as a subprocess, assert on exit status and report section. All
existing tests keep passing unchanged (they pass `--paths docs`).

**Pinned (§2.2)**
1. `PATH@<full sha>:<a>` and `@<7-char sha>:<a>-<b>` on a file that later
   changes or is deleted in the working tree: resolve, exit 0, not at risk
   even though the branch edits the file.
2. Pinned to a lightweight tag and to an annotated tag, including a tag name
   with `/`: resolve.
3. Last line past the end at `REV`; line 0: broken.
4. Path absent at `REV` (present in the working tree): broken, no fallback to
   the working tree.
5. A basename that exists at `REV` only under a directory: broken (no search).
6. Sha that exists in no object; a 7-40 hex string that is not a commit:
   broken.
7. `REV` = a branch name (`master`), and `HEAD`: broken, "not a pin".
8. Shallow clone (`git clone --depth 1 file://…` of a repo whose pin is older):
   broken, message says shallow.
9. A pinned citation must not also be reported by the unpinned pattern (one
   entry, not two).

**Legacy (§2.1)**
10. Existing legacy tests pass unchanged; `legacy/x.py@legacy-archive:<a>`
    resolves the same as the unpinned `legacy/` form.

**ID and heading-name (§2.3, §2.4)**
11. `§2`, `§3.2`, `§E1`, `§3C` (letter heading under a numeric one), `§C`:
    resolve; each with the ID removed from the target: broken.
12. A matching `## 7.` heading inside a fenced block does not declare `7`.
13. `§ Name` (space after `§`) and "section 2", naming a heading the target
    lacks: not citations, exit 0.
14. `` `<file>.md`, *<Name>* `` against `## Name [tag]` and `## Name: rest`:
    resolve; against `## Named`: broken.
15. Partial path `dir/SKILL.md` with several `SKILL.md` in the tree: resolves
    by suffix; a bare `SKILL.md` stays ambiguous.

**Governed (§3.1)**
16. Unpinned line citation into each governed kind (`CLAUDE.md`, `testing.md`,
    `.claude/agents/x.md`, `.claude/skills/y/SKILL.md`, `.claude/hooks/z.py`,
    `tools/check_w.py`): `== unpinned`, exit 1. The same file cited pinned, by
    ID or by heading name: exit 0.
17. Unpinned line citation into an ungoverned file the branch edits: at risk,
    exit 0 (regression).

**Scan set (§4)**
18. Default run (no `--paths`) finds a planted unpinned `testing.md:<n>` in:
    a `.cpp` `//` comment, a `.hpp` `/* */` block spanning lines, a `.py`
    `#` comment, a `.py` function docstring, a `CMakeLists.txt` `#` comment, a
    root-level `*.md`. Each reported at its own line; exit 1. (This is the
    plant §7 step 2 of the design names.)
19. Not reported: the same token in a C++ string literal, in a raw string
    spanning lines, in a string literal that follows a digit separator on the
    same line (`0x9513'0000`), in a Python string that is not a docstring, in
    a CMake quoted argument, and in a URL (`https://host/testing.md:<n>`).
    And conversely: a `//` comment that follows a string containing `//` on
    the same line is still scanned.
20. A missing target in a `.cpp` comment is broken (the scan feeds every
    form, not just the governed check).
21. Files under `lib/`, `legacy/` and `.claude/worktrees/`, and gitignored
    files: not scanned; an untracked, unignored file: scanned.
22. A `.py` that does not parse: `== warnings`, exit unaffected.

Mutation testing is not required (the suite is not invariant-critical).

---

## 7. LOC estimate

Production lines as `CLAUDE.md` §2 counts them, all in
`tools/check_citations.py`; tests excluded. The table is the estimate. Measured
at green by `@reviewer`: 139 → 332 counted lines, net **+193**, inside the
estimate.

| Part | Lines |
|---|---|
| Pinned pattern, rev classification (sha, tag, branch, shallow), `cat-file` | 40-50 |
| ID and heading-name patterns, fence-aware heading parser, nested IDs | 35-45 |
| Governed set and the unpinned category | 10-15 |
| Comment extraction: Python (tokenize, ast), C-family lexer, hash | 50-65 |
| Default sources via `git ls-files`, suffix-path resolution | 15-20 |
| Report sections and exit status | 10-15 |
| **Total** | **~160-210** |

Above the design's 80-120, because the C-family lexer and the ID parser were
not counted there. Well under the ceiling in `CLAUDE.md` §2. One file, kept as
one: the suite copies the tool into a scratch repository, and a split would
make every fixture copy two files for no testability gain (the suite is
black-box).

Order of commits: red (@tester: suite plus §5.1), doc conversion (§5.2),
green (@developer: tool only), then `@reviewer` with the planted cases of
test 18 re-run in the real tree.

## Review

**Round 1, `@reviewer`, range `7426b69..45d1bd4`: CHANGES REQUESTED.**
- **Code and gates pass.** LOC: `tools/check_citations.py` goes from 139 to 332 counted lines, net +193, inside the estimate. The only production file changed. Red (7618d05) precedes green (5e3d8a6); green touches no test file; 1ca498c sits between them. Red suite against the red-time tool: 77 failed and 70 passed. At green: 147 pass (138 + 9).
- **Plants.** An unpinned `.cpp` comment citation, a missing sha, `@HEAD` and `@master`, a citation inside a C++ string, a line past the end, a missing §ID and a missing heading were each caught or ignored as specified. A spot-check of 10 of the 92 distinct pins in 1ca498c: all quotations hold. ruff, `mypy --strict` on the tool, and the governance gates are all green.
- **Blocking, all prose:**
  - the tool's docstring said "five forms" for four;
  - this file's status line and LOC section had not been updated;
  - §2.3's letter-label rule differed from the code, which uses the direct parent only; the spec was aligned to the code;
  - two design-doc citations quoted the pre-red C++ comments; they are now pinned at 7426b69.
