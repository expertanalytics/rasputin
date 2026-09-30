# A generic harness layer: design

Status: design only (@architect, 2026-09-30). No code and no rule file is
changed by this document. Estimates are labelled as estimates; everything
labelled "checked" was run or read on this branch.

**Decision in one paragraph.** Split the harness into a **generic layer** (the
loop, the approval boundary, the in-flight protocol, claims discipline, the
hooks, and the tools that carry no domain knowledge) and a **project profile**
(stack, constraints, gates, paths, owner, domain skills, domain oracles). Keep
the generic layer in its own repository as a **copier template** and vendor it
into each project with `copier update` (a 3-way merge). Do not use a Claude Code
plugin as the carrier yet: a plugin cannot deliver `CLAUDE.md`, stores its
parameters per user rather than per project, renames every persona to
`<plugin>:<name>`, and takes the rule text out of the project's tree, where the
governance guard, `@reviewer` and `check_citations.py` would otherwise see it.
Parameters used at runtime go in one committed file, `.claude/profile.toml`.
Each rule is stated once, next to its check, and everything else points at
it (§1.2; the survey in §1.5 lists 22 rules that are stated more than once
today). Rules are cited by stable ID, never by line number or ordinal.

---

## 0. Prior art: what this builds on

No novelty is claimed here. The two mechanisms this builds on:

- **Claude Code's own extension points** (code.claude.com/docs/en, fetched
  2026-09-30). The facts the design depends on are listed below, and each is
  cited at the point where it decides something.
  - **Plugins.** A plugin can ship skills, agents, hooks, MCP servers and
    output styles. It does not load a `CLAUDE.md` at its root: "To include
    instructions in a plugin, write them as a skill"
    (`plugins/components`, Skills). Plugin agents are named `<plugin>:<name>`
    and ignore `permissionMode`, `hooks`, `mcpServers` and `initialPrompt`
    (`plugins/components`, "Frontmatter fields in plugin agents").
    `userConfig` values go to `pluginConfigs`, whose scope is "User or
    managed", so they cannot be set per project (`settings-reference`,
    `pluginConfigs`). A project-scope install commits `enabledPlugins`, but
    "each collaborator also runs `claude plugin install … --scope project`
    once" (`plugins/install`, "Choose an install scope"). Auto-update is off by
    default for each marketplace (`plugins/publish`). A marketplace plugin is
    copied into a versioned cache and cannot read paths above its own root
    (`plugins/loading`). Cloud sessions do not add the marketplaces a
    repository declares (`plugins/loading`).
  - **Subagents.** A custom subagent loads the CLAUDE.md hierarchy and a
    git-status snapshot, but not the conversation (`sub-agents`, "What loads
    at startup"). The `skills:` frontmatter "injects the full content"
    (`sub-agents`, supported frontmatter). When names collide, precedence is
    managed, then `--agents`, then `.claude/agents/`, then `~/.claude/agents/`,
    then plugin (`sub-agents`, "Choose the subagent scope").
  - **Hooks.** Hooks from settings files and from plugins also fire inside
    subagents, and the payload carries `agent_type` and `agent_id` (`hooks`,
    common input fields). `${CLAUDE_PROJECT_DIR}` stays at the main checkout
    in a worktree session (`worktrees`).
  - **Memory.** `CLAUDE.md` can `@import` files. `.claude/rules/*.md` loads
    unconditionally unless it is path-scoped. A symlink in rules that points
    outside the working directory needs the external-import approval
    (`memory`, "Import additional files" and "Share rules across projects with
    symlinks").
- **Copy-and-update templating: copier** (MIT; version 9.18.2 checked on PyPI).
  `copier update` regenerates the project from the template at the old
  version, diffs that against the project, and applies the diff onto the new
  version. Local edits therefore survive an upgrade, and conflicts come out as
  `.rej` files or inline markers (copier.readthedocs.io/en/stable/updating,
  "How the update works"). The answers are committed in
  `.copier-answers.yml`. cookiecutter + cruft does the same job; copier does
  it natively.

What differs from both: plugins and copier do not settle *what* is generic.
The design adds a classification rule (§1), a runtime profile file (§3) and a
drift gate that keeps generic regions upstream-only (§2).

---

## 1. The split

### 1.1 The classification rule

An item is **generic** if it would be true and useful in a project with a
different stack, domain and owner, and can be stated without naming any of
them. Otherwise it belongs to the **profile**. A generic statement that needs
a value (a name, a number, a command or a path) takes it from the profile and
does not state it.

Where one file mixes the two, it is divided into marked regions:
`<!-- harness:begin <id> -->` … `<!-- harness:end <id> -->` around the generic
text, and everything else belongs to the project. Python hooks and tools are
not divided: a generic script is generic as a whole, and it reads whatever it
needs from `.claude/profile.toml` (§3).

### 1.2 One statement per rule (generic principle)

> **A rule is stated once.** A constraint's statement, including its value
> (a number or a list), lives in one place. For a project constraint that
> place is the profile, beside the machine check that enforces it. Everything
> else points at it by ID, and a pointer never repeats the value, the list, or
> any part of either. A pointer may add *how* for its own role (for example,
> "read GeoTIFF with `tifffile`"), but not *what* is required. Where the
> statement and its gate both hold the same list, the gate checks that the
> two agree.

The last sentence already has a working instance in the tree:
`tools/check_prohibited_deps.py` parses the list in `CLAUDE.md` §2 and fails
when that list and the gate's table disagree (checked: `grep -n CLAUDE
tools/check_prohibited_deps.py`, lines 222-291). The generic layer makes this
the pattern for every listed value, and a cheap gate can enforce it (§7,
step 1b): each constraint value in `profile.toml` (the ceiling, the coverage
floor, each prohibited name) may appear in the rule files only in its home
file.

The same rule applies inside the generic layer, and every generic rule has
one home:

| Generic rule family | Home (the one statement) |
|---|---|
| Cold start, in-flight state, claims, approval boundary, pre-push review, harness, data-folder channel | `.claude/REQUIRED-READING.md` |
| The loop, red/green trace, merge policy, roadmap row, cost constraints | `docs/increments/README.md` |
| Role-specific duties | the persona's generic region |
| Principle IDs, index and retirement | `docs/PRINCIPLES.md`, as an **index only** (§4) |

### 1.3 Artifact classification

| Artifact | Class | What is profile inside it |
|---|---|---|
| `CLAUDE.md` | **Profile file** with one generic region | §1's roster wording ("geometry fuzzing", "C++20"), §2 in full, §4 in full. §3 (the loop pointer) is generic. The section numbers are frozen as IDs (§3.2) |
| `.claude/REQUIRED-READING.md` | **Generic**, with three profile sections | see 1.4 |
| `.claude/agents/orchestrator.md` | Generic | step 2's geometry list and its §3C/§3D references; step 6's refine/mesh trigger |
| `.claude/agents/architect.md` | Mixed | §1-§3 (SoC across C++/Python, the pybind11 firewall, Concepts/Protocol, No-GDAL) are profile. §4 (literature) and §5 (blueprint before code) are generic |
| `.claude/agents/tester.md` | Mixed | The floor's value, §2 tiers (Catch2/pytest), §3A geometry, §3C ingestion, and §3D's two oracles as instances are profile. The rest of §1, the §3D *pattern* (§5a) and §4 are generic |
| `.claude/agents/developer.md` | Mixed | The skills list and the verification checklist (asan/ubsan, `-Werror`, `-ffp-contract=off`) are profile; the green-step duties are generic |
| `.claude/agents/reviewer.md` | Generic | the gate names (ruff, mypy, `-Werror`) move to the profile |
| `.claude/agents/perf.md` | Optional generic skeleton | §1-§2 (bench.py, refine, `pmset`, the ASan/DYLD note) are profile; §3 (evidence) and §4 (reporting) are generic |
| `.claude/skills/*` | **Profile**, all four | the generic layer ships no skills, only the rule for writing one (§5d) |
| `.claude/settings.json` | Generic | none |
| `.claude/hooks/guard_push.py` | Generic | none; the forge commands assume git + `gh`, which becomes a profile `forge` value |
| `.claude/hooks/guard_governance.py` | Generic engine | the governed set = the generic manifest ∪ `profile.rule_files` |
| `.claude/hooks/gates_after_commit.py` | Generic runner | the gate list and tool discovery (ruff in `.venv`) |
| `tools/session_state.py` | Generic | owner, `ASK <OWNER>` marker, base branch, roadmap path and its "done" words |
| `tools/check_citations.py` | Generic | search roots, the legacy prefix and tag, the base branch |
| `tools/check_prohibited_deps.py` | Generic engine (the import/include/pyproject/CMake scanners and the RULED/SPELLING/PENDING authority model) | the `PROHIBITED` table and the name of its home section |
| `tools/check_detria_boundary.py` | Profile | a later "confinement" gate could generalise it, but not from one instance (§5b) |
| `tools/bench.py` | Profile | — |
| `docs/increments/README.md` | Generic, with profile paragraphs | see 1.4 |
| `docs/PRINCIPLES.md` | Generic index | A5 (stale artifacts, C++-specific), E3 (I/O boundary) and E4 (ceiling) are profile pointers |
| `testing.md` | **Profile** | its conventions are generic and go into a template skeleton: `[live]`/`[planned]` markers, "never lose a bug twice", and the characterisation-test exception to "we do not test third-party code" |
| `.claude/current-task/` protocol | Generic | the owner marker only |
| `.gitignore` lines (`.claude/current-task/`, `.claude/worktrees/`) | Generic | — |

### 1.4 Mixed files, by section

**`.claude/REQUIRED-READING.md`**

| Section | Class | Profile content to move |
|---|---|---|
| Before you act in a cold or resumed session | Generic | the owner's name |
| While you act: in-flight state goes on disk | Generic | the `ASK <OWNER>:` marker |
| Before you write, judge or plan code | Generic rule | the four skill names go to the profile |
| Claims | Generic | the worked example (plant into `CMakeLists.txt`, run `check_prohibited_deps.py`) becomes "plant what the gate forbids" |
| Stale artifacts | **Profile** | all three items. Generic rule: "the profile lists every artifact that goes stale silently, with its rebuild command; rebuild before you measure" |
| Before you publish | Generic | — |
| The harness | Generic | the list of proposed hooks moves to the generic layer's open items |
| Data, scratch and temp folders | Generic rule | the two paths go to the profile |

**`docs/increments/README.md`**

| Paragraph | Class |
|---|---|
| Intro, "Why on disk", the loop's steps 1-4, amendments, ROADMAP row, merge commits | Generic. "post-CGAL", `legacy-archive` and the refine/mesh trigger are profile values |
| Step 1, *Legacy* | Generic shape ("if a predecessor codebase exists, paste the grep and the file list it returned"); where that codebase lives is profile |
| Acceptance: refine or mesh code | **Profile**. Generic shape: the profile may declare an acceptance run with trigger paths, a tool, an environment probe and an evidence directory |
| Cost: reference, do not restate; mutation testing; model matching; documentation defects | Generic; the examples (predicates, kernels) are profile |
| Cost: template cross products | Profile (C++) |
| Cost: "Independent suites run as parallel agents" | **Contradicts PRINCIPLES D1**; see §6, defect 1 |

### 1.5 Survey: rules stated more than once (checked on this branch)

Found with `grep -nE` over the rule files (`CLAUDE.md`,
`.claude/REQUIRED-READING.md`, `.claude/agents/*.md`,
`.claude/skills/*/SKILL.md`, `docs/increments/README.md`,
`docs/PRINCIPLES.md`, `testing.md`). Each line citation is pinned to
`3245958`, the tree the survey read (`git show 3245958:<path>`), so later
edits to these files do not move it. "Home" is where the rule should be stated
once. Each line under "Restated at" repeats the rule; lines marked *ok* only
point at it.

| Rule | Home | Restated at | Fix |
|---|---|---|---|
| Prohibited dependencies | `CLAUDE.md@3245958:23-26`, checked by `tools/check_prohibited_deps.py` | `.claude/agents/architect.md@3245958:30` ("No-GDAL and No-CGAL"); `.claude/skills/computational-geometry/SKILL.md@3245958:24` ("Do not use CGAL"); `.claude/skills/geospatial-data-formats/SKILL.md@3245958:10-11` (GDAL, OGR, Fiona) and `.claude/skills/geospatial-data-formats/SKILL.md@3245958:14` ("Never `rasterio`"), plus "Zero-GDAL" in its description and summary (`.claude/skills/geospatial-data-formats/SKILL.md@3245958:3`, `.claude/skills/geospatial-data-formats/SKILL.md@3245958:8`). Pointers, *ok*: `.claude/agents/orchestrator.md@3245958:56`, `testing.md@3245958:351` | Each becomes "see `CLAUDE.md` §2". The geospatial skill keeps its how-to (tifffile, shapely, pyproj) |
| Size ceiling | `CLAUDE.md@3245958:15-22` | `docs/PRINCIPLES.md@3245958:220` (its heading states "700"). Pointers, *ok*: `.claude/agents/developer.md@3245958:30`, `.claude/agents/reviewer.md@3245958:16`, `.claude/agents/tester.md@3245958:68`, `.claude/agents/orchestrator.md@3245958:25`, and the three skills' "Change Limit" lines (redundant but harmless) | Rename E4 without the number |
| Reconcile LOC against the estimate | `.claude/agents/reviewer.md@3245958:43-45` (check 3) | `.claude/agents/reviewer.md@3245958:16-20` (the same file, again); `.claude/agents/developer.md@3245958:30-33`; `.claude/agents/orchestrator.md@3245958:29-30` | Keep check 3; the others point at it |
| Coverage floor | `testing.md@3245958:279`, checked by `pyproject.toml:102` (`--cov-fail-under=85`) | `CLAUDE.md@3245958:9` ("85%"); `.claude/agents/tester.md@3245958:13-14` and `.claude/agents/tester.md@3245958:20-23` (twice in one file); `testing.md@3245958:18` and `testing.md@3245958:283-288` | State it once at `testing.md@3245958:279`; drop the number everywhere else |
| Approval boundary | `.claude/REQUIRED-READING.md@3245958:100-116` | `docs/PRINCIPLES.md@3245958:200-206` (E1 restates the list of acts); the message in `guard_push.py`. Pointers, *ok*: `CLAUDE.md@3245958:43-45`, `.claude/agents/orchestrator.md@3245958:38-43` | E1 becomes an index entry (§4) |
| `@reviewer` before a push | `.claude/REQUIRED-READING.md@3245958:118-126` | `docs/PRINCIPLES.md@3245958:208-211` (E2) | Index entry |
| Claims discipline | `.claude/REQUIRED-READING.md@3245958:59-80` | `docs/PRINCIPLES.md@3245958:24-57` and `docs/PRINCIPLES.md@3245958:91-108` (A1-A4, B1, B2 restate it in different words) | Index entries |
| Stale artifacts | `.claude/REQUIRED-READING.md@3245958:81-98` | `docs/PRINCIPLES.md@3245958:59-72` (A5, which adds a fourth artifact, the stale `.pyc`, that is stated nowhere else) | Move the fourth into the home, then index |
| Recap contents | `.claude/REQUIRED-READING.md@3245958:12-14` | `.claude/agents/orchestrator.md@3245958:45-49` (the same four items, plus "at every round") | The orchestrator keeps the trigger and points for the contents |
| Red/green trace | `docs/increments/README.md@3245958:49-52` | `.claude/agents/developer.md@3245958:23-26`; `docs/PRINCIPLES.md@3245958:138-144` (C2) | The developer keeps its duty in one line and points |
| Merge commits, never squash | `docs/increments/README.md@3245958:75-78` | `docs/PRINCIPLES.md@3245958:141-142` | Index entry |
| CI is authoritative | `CLAUDE.md@3245958:75-81` | `.claude/agents/reviewer.md@3245958:24-31` (with the command); `docs/increments/README.md@3245958:53` | The reviewer keeps its precondition as a pointer |
| Literature before design | `docs/increments/README.md@3245958:28-34` | `.claude/agents/architect.md@3245958:33-44` (§4, which is more detailed than its home); `.claude/skills/computational-geometry/SKILL.md@3245958:46-52` (with a pointer) | Merge §4's three points into the home. §4 keeps the duty ("you own it") and a pointer; the skill keeps only its domain examples |
| Power state with each run | `docs/increments/README.md@3245958:82` (the Acceptance section, which calls itself "the one statement"; the power-state bullet is at 88-90) | `.claude/agents/perf.md@3245958:39-42` | `perf.md` points |
| Adversarial geometry list | `.claude/agents/tester.md@3245958:37-42` (§3A) | `.claude/agents/orchestrator.md@3245958:19-20`; `.claude/skills/computational-geometry/SKILL.md@3245958:29`; `.claude/agents/tester.md@3245958:3` (the description) | The orchestrator says "per tester §3A" |
| Pybind11 isolation | `.claude/skills/modern-cxx/SKILL.md@3245958:27` | `.claude/agents/architect.md@3245958:22` | The architect points |
| Concepts / `typing.Protocol` | `.claude/skills/modern-cxx/SKILL.md@3245958:11`, `.claude/skills/python-development/SKILL.md@3245958:34` | `.claude/agents/architect.md@3245958:24` | The architect points |
| `-Werror` | `CMakeLists.txt:52,70` (the check) | `.claude/skills/modern-cxx/SKILL.md@3245958:32`; named in `.claude/agents/developer.md@3245958:35`, `.claude/agents/reviewer.md@3245958:13` and `.claude/agents/reviewer.md@3245958:35` and `.claude/agents/orchestrator.md@3245958:31` | State it once in `CLAUDE.md` §4, beside the build; the others name "the compiler gate" |
| I/O boundary | `CLAUDE.md@3245958:28` | `docs/PRINCIPLES.md@3245958:215-216` (E3 calls itself a pointer but restates the rule) | Index entry |
| Documentation defects | `docs/increments/README.md@3245958:130-133` | `docs/PRINCIPLES.md@3245958:146-150` (C3) | Index entry |
| Roadmap row | `docs/increments/README.md@3245958:69-73` | `docs/PRINCIPLES.md@3245958:152-156` (C4) | Index entry |
| Skills to invoke | `.claude/REQUIRED-READING.md@3245958:48-53` | `.claude/agents/developer.md@3245958:13-17` (the same four, as a list) | `developer.md` points |

Two patterns stand out. `docs/PRINCIPLES.md` restates nearly every rule it
indexes, which is why §4 turns it into an index. And the personas restate
their home rule where they could name it and add only their own duty.

---

## 2. Mechanism for reuse

### 2.1 Comparison

| | A. GitHub template repo | **B. copier template** | C. Claude Code plugin | D. git submodule | E. `~/.claude` (user level) |
|---|---|---|---|---|---|
| New project adopts | "Use this template"; greenfield only | `copier copy <harness> .`, answering 4-5 questions; works on an existing repo | add the marketplace, `install --scope project`, commit `enabledPlugins`; every collaborator installs once (`plugins/install`) | `git submodule add`, plus symlinks or `@import`s into `.claude/` | copy files into `~/.claude` |
| Generic improvement reaches projects | Never; the copy is one-shot | `copier update`: a 3-way merge, landing as an ordinary PR with the full text diff | a version bump; auto-update is off by default (`plugins/publish`); lands outside the project's git history | a pointer bump; the PR shows a SHA, not the rule text | at once, on one machine, for every project |
| Project improvement flows back | by hand | a PR to the harness repo; the drift gate (2.3) makes a local edit to a generic region fail, so the edit must go upstream | a PR to the plugin repo | a commit inside the submodule | — |
| Rule text in the project tree, where `guard_governance`, `@reviewer` and `check_citations` see it | yes | **yes** | **no** | yes, but a change arrives as an opaque bump | no |
| The rules that governed a commit can be recovered from git | yes | **yes** | no | yes, via the pointer | no |
| Per-project parameters | find and replace | copier answers (committed) + `profile.toml` | `userConfig` is **user or managed scope only** (`settings-reference`, `pluginConfigs`); project values need a repo file anyway | a repo file | none |
| Delivers `CLAUDE.md` / always-loaded rules | yes | yes | **no**: "doesn't load a CLAUDE.md at the plugin root" (`plugins/components`) | via `@import` | `~/.claude/CLAUDE.md` and `~/.claude/rules/`, personal |
| Persona names | `tester` | `tester` | `harness:tester` (`plugins/components`, Agents): every `@tester` in the rule text and every `agent_type` match in R-B changes | `tester` | `tester` |
| Worktrees and CI | tracked files, fine | tracked files, fine | per-user install; cloud sessions skip repo-declared marketplaces (`plugins/loading`) | `git worktree add` does not initialise submodules, so Claude's isolated worktrees start without the harness; CI needs `submodules: true` | not in CI; one machine only |
| Extra dependency | none | copier (MIT), a dev-time CLI; hooks and gates stay stdlib-only | Claude Code only | none | none |

### 2.2 Recommendation: B, a copier template, vendored

Option B is the only one that keeps the rule text in the project's tree and
still lets improvements flow in both directions. Keeping the text in the tree
is not a convenience. The approval boundary, `guard_governance.py`, the
pre-push `@reviewer` pass and `check_citations.py` all work because a rule
change is a diff in the project's own repository. A plugin moves that diff to
another repository and applies it silently at session start, so each of those
controls stops seeing rule changes. Option C should be reconsidered for the
hooks alone if Claude Code gains project-scoped plugin configuration. The
rules themselves should not move to a plugin while the controls above depend
on seeing them in the tree.

### 2.3 Blueprint

```
harness repo (copier template, tagged vMAJOR.MINOR)
  copier.yml                  questions: owner, project_name, forge, base_branch
  template/
    CLAUDE.md.jinja                   skeleton; _skip_if_exists (project-owned after copy)
    .claude/REQUIRED-READING.md.jinja generic, owner rendered
    .claude/agents/*.md.jinja         generic region + empty profile region
    .claude/hooks/*.py                generic, read profile.toml
    .claude/settings.json             hook wiring
    .claude/profile.toml.jinja        skeleton; _skip_if_exists
    tools/session_state.py, check_citations.py, check_prohibited_deps.py (engine),
    tools/check_harness.py            drift gate, stdlib only
    docs/PRINCIPLES.md                index (§4)
    docs/increments/README.md.jinja   the loop, with a marked profile tail
    testing.md.jinja                  skeleton; _skip_if_exists
        │  copier copy / copier update (3-way)
        ▼
project repo
  .copier-answers.yml         committed; copier owns it
  .claude/harness.lock        sha256 of every generic file and region, written by the template
  generic files and regions   upstream-only: check_harness.py fails on a local edit
  profile regions and files   project-owned
  .claude/profile.toml        read at runtime by the hooks and tools
```

- **Drift gate.** `tools/check_harness.py` hashes each generic file and each
  `harness:begin/end` region, and compares the hashes with
  `.claude/harness.lock`. It runs among the fast gates. It needs no network and
  no copier: the lock is the template's own output.
- **Regions, not includes, for personas.** A persona's body is its system
  prompt, so the generic and profile text must end up in one file. The
  `skills:` preload could deliver the profile half as a skill (§6, defect 2),
  but the local rules say that preload is unreliable. Regions work whichever
  of the two is true.
- **An update is a PR like any other.** `copier update` on a branch, then the
  gates, then `@reviewer` with the rule-file scope from `REQUIRED-READING`
  ("Before you publish"), then Ola's yes. `guard_governance` asks on every
  file it touches, as intended.

---

## 3. Parameters and citation stability

### 3.1 Parameters

Two carriers, with no overlap between them:

- **copier answers.** Only the values that must appear *in rule text*:
  `owner`, `project_name`, `forge`, `base_branch`. They are rendered once, at
  copy or update time.
- **`.claude/profile.toml`.** Everything the hooks and tools read at runtime.
  It is TOML because `tomllib` is in the standard library, so the hooks stay
  dependency-free. The file is the *machine* statement of each value; where a
  value also needs a prose statement, that statement sits beside its check
  (1.2).

| Parameter | Where it is today | Generic handling | rasputin value |
|---|---|---|---|
| Owner's name | `REQUIRED-READING.md` (3×), `session_state.py` ("ASK OLA", "Waiting on Ola"), `tests/python/test_session_state.py` | answer `owner`, rendered into the text; `profile.owner` for the tools; marker `ASK <OWNER>:` | Ola |
| Size ceiling | `CLAUDE.md` §2 | profile prose only (`CLAUDE.md` §2); the generic text says "the size ceiling (`CLAUDE.md` §2)" and never the number | 700, with its counting rule |
| Coverage floor | `testing.md`, `pyproject.toml` | stated once in `testing.md`, beside the `--cov-fail-under` that enforces it; generic `@tester` says "the coverage floor (`testing.md`)" | 85 |
| Gate commands | `CLAUDE.md` §4; `gates_after_commit.py` hard-codes five | `[[gates]]` in `profile.toml` (`name`, `cmd`, `fast`); `gates_after_commit.py` runs the `fast` ones; `CLAUDE.md` §4 points at the table and keeps only the build and test commands that are not gates | today's five fast gates |
| Environment probe for timings | `README` acceptance, `perf.md` (`pmset -g batt`, macOS only) | generic rule: "record every condition that changes timings with each run, and compare like with like"; `profile.perf.env_probe` gives the command | `pmset -g batt` |
| Data and scratch dirs | `REQUIRED-READING.md` | `profile.channels.forbidden_dirs`; the generic rule names them through the profile | `../rasputin_data`, `../rasputin_scratch`, `/tmp` |
| Base branch, forge | `session_state.py` (`origin/master`), `check_citations.py` (`--base master`), `CLAUDE.md` §4 (`gh`) | answers + profile | master, github |
| Roadmap | `session_state.py` (`ROADMAP.md`; "shipped", "landed", "unscheduled") | `profile.roadmap` | as today |
| Governed rule files | `guard_governance.py` `GOVERNED` | the generic manifest ∪ `profile.rule_files` | adds `testing.md` and `.claude/skills/` (§6, defect 3) |
| Skills to invoke | `REQUIRED-READING.md`, `developer.md` | `profile.skills` and one line in `CLAUDE.md` | the four |
| Stale artifacts | `REQUIRED-READING.md`, `PRINCIPLES` A5 | a profile section of `CLAUDE.md` (§7, step 5) | four items |
| Legacy archive | `check_citations.py`, `README` step 1 | `profile.legacy = {prefix, tag}` | `legacy/`, `legacy-archive` |
| Persona path classes (R-B) | the audit's §7 table | `profile.paths.{production,tests,design,evidence}` (§6) | `src/ include/ bindings/ src_python/`, `tests/`, `docs/increments/`, `docs/benchmarks/` |

### 3.2 Citation stability

Rule text changes, and a citation into it has to survive the change. Today
the tree answers this in three different ways: line numbers (checked by
`check_citations.py`), section ordinals kept with gaps ("There is no section
B; letters are kept stable for citations", `.claude/agents/tester.md@3245958:45-46`; "There are no
sections 2-4", `.claude/agents/reviewer.md@3245958:13`), and principle IDs (A1…E4). The generic rule
keeps the third and drops the other two:

1. **Living rule text is cited by ID.** An ID is a short label in a heading
   (`E1`, `§3D`, `§2`). It is assigned once and never renumbered or reused;
   retiring a rule leaves its ID dead. An ID is a name, not a position, so a
   gap needs no explanation and the "kept for citations" notes go.
2. **Line citations into a rule file are for history only, and are pinned.**
   A retrospective or review that quotes a rule as it stood writes
   `path@<sha>:<line>`, and `check_citations.py` resolves that through
   `git show`. The tool already does this for the `legacy-archive` tag; the
   design generalises it to any revision.
3. **The checker enforces both.** It fails when an ID citation names no
   heading, or a pinned citation does not resolve. It warns on an unpinned
   `file:line` that points into a governed rule file.

For rasputin, `CLAUDE.md` §1-§4 become IDs as they stand. The numbers must not
move; `git grep -hoE 'CLAUDE\.md`?,? *(§ *|section )[0-9]+' -- . ':!legacy' | wc -l`
counts the references that depend on them. New
profile sections are appended as §5 and later, never inserted. The existing
labels stay as IDs: tester `3A`/`3C`/`3D`, reviewer `1`/`5`, architect `4`,
PRINCIPLES `A1`-`E4`.

---

## 4. `docs/PRINCIPLES.md`: replace the fields, keep the IDs

**Verdict: replace.** The Origin and Last-exercised mechanism belongs in
neither layer, because it is not operating. Checked: `grep -o 'Last exercised:
[^.]*' docs/PRINCIPLES.md` returns 5b or 5c for 20 of the 21 entries that have
the field, and "not yet" for the other. Meanwhile the increment files run to
`22-auto-catchment.md`. By the file's own rule ("not exercised in three
increments is reviewed for retirement"), every principle is overdue, and no
retrospective has run the review. The field is a resolved value, which B2
itself forbids ("write the rule and the command, not the resolved value"),
and it went stale in exactly the way B2 predicts. The file also restates most
of the rules it indexes (1.5).

**Generic replacement:**

| Today | Generic layer |
|---|---|
| Rule + Apply text | **removed**. The rule is stated in its home (1.2), and the entry links to it |
| `Origin:` commits and retrospectives | **removed** from the generic file. History stays in the project's retrospectives. A project that wants the pointers can keep a profile file, `docs/principles-origin.md`, mapping ID → retrospective |
| `Last exercised:` a stored value | **a command**: `git grep -nw <ID> docs/retrospectives docs/increments`. A retrospective or increment that applies a principle names its ID, so the grep is the record |
| Retirement after 3 increments unexercised | **event-driven, at a retrospective**. Retire a principle when (i) a gate now enforces it, in which case the entry's check becomes that gate; (ii) the owner rules on a conflict between it and another principle; or (iii) the grep finds no citation across the last *N* retrospectives, with *N* in the profile and 3 by default. A retired ID stays dead |
| "Only a retrospective promotes a log entry to a principle" | kept, generic |

Each entry becomes one table row: `ID | title | home (file § heading) | check
(a command, or "object-identity question")`. The IDs and the A-E families
stay, so `guard_push.py`'s "PRINCIPLES.md E1" and every other ID citation keep
resolving. A5, E3 and E4 are profile rows: rasputin's index keeps them, and
the generic template does not ship them.

---

## 5. External input (Gemini): adopt, adapt or reject

| Idea | Verdict | Reason |
|---|---|---|
| (a) Dual-oracle verification pattern | **Adapt** | The pattern is generic, but the number two is not. See below |
| (b) Map domain rules onto generic patterns | **Adopt two, adapt one, reject one name** | See below |
| (c) Drop the persona "Required reading" header | **Reject** | See below; the evidence is new |
| (d) Trim generic language advice from the skills | **Adopt, with a test** | See below |

**(a) Oracle completeness, not "dual oracle".** The generic rule for the
tester's generic region:

> The design names the oracles for each property it claims. A property test
> checks every named oracle on every path it exercises, meaning every
> configuration switch, thread count and tolerance it varies. An oracle is
> recomputed from the output using the producer's relation, and never from
> the producer's own records.

The count is a fact about the domain. An approximation algorithm has a
structural invariant and a bound; a parser has a round-trip oracle and a
rejection oracle; a CRUD layer may have one. Fixing it at two would make
projects invent a second oracle or skip the rule. "SLA" is rejected as a
name: a service-level bound is about timing or availability, and it belongs
in the acceptance run (`@perf`), not in a property test. A property test that
asserts timing becomes flaky, which breaks the tester's own determinism
mandate. Rasputin's two oracles stay in tester §3D as the profile's
instances.

**(b) The mapping.**
- Delaunay → *invariant oracle*: adopt, as part of (a).
- Tolerance → *bound oracle*, not "SLA oracle": adopt, renamed as in (a).
- No-GDAL/No-CGAL → *dependency-boundary invariant, machine-checked*: adopt.
  The generic layer ships `check_prohibited_deps.py`'s engine and its
  authority model (each prohibited key must quote the human who ruled it), and
  the profile supplies the table. The detria confinement gate is the same
  shape ("this header in these translation units only"), but it is one
  instance, so it stays in the profile until a second one appears.
- `-ffp-contract=off` → *deterministic execution*: **adapt with a new name**.
  FMA contraction is deterministic for a given build. The hazard is that
  results differ across toolchains and platforms. The generic rule is
  "toolchain-independent results: a numeric test must pass under the build
  configuration that removes the compiler's freedom (the profile names it),
  as well as under the default". The flag itself stays in the profile.

**(c) Keep the pointer, one line, in the generic region.** The main session's
argument ("subagents do not inherit context") is imprecise. Custom subagents
do load the CLAUDE.md hierarchy (`sub-agents`, "What loads at startup"), and
`CLAUDE.md` §3 already sends every persona to `REQUIRED-READING.md`. So the
header is a duplicate pointer, not the only route. It should still stay,
because of *which* `CLAUDE.md` a subagent receives. The docs list "every level
of the CLAUDE.md hierarchy the main conversation loads" (`sub-agents`, "What
loads at startup"). That is the copy the **main session loaded when it
started**, not a fresh read of any checkout. So it can be stale against the
main checkout and the worktree alike, and nothing in the subagent's context
shows that it is.

Reproduced twice. This spawn's injected `CLAUDE.md` carried
`@migration-expert` and dated history, which the worktree's file no longer
had. `@reviewer` later found the main checkout's file byte-identical to the
worktree's, yet its own injected copy still carried `@migration-expert`, "20b
counted blank lines" and "`legacy/` is exempt". To reproduce: change
`CLAUDE.md` on disk after the main session has started, spawn any persona,
and have it quote a line that the change touched. The spawned persona quotes
the old text.

The persona's pointer is a relative path that the subagent reads with its own
tools, in its own working tree, when it acts. It is the only route that
yields the rules as they are on the branch now. On any branch that changes
the rules, or in any session older than the last rule change, the two routes
disagree, and only the pointer is right. One line in a file that is identical
in every project costs nothing to maintain. For the same reason, do not
replace the pointer with a `CLAUDE.md` `@import`, which is expanded at the
main session's start along with `CLAUDE.md` itself.

**(d) The skill-writing rule (generic).** A line in a skill stays only if it
(i) records a project decision that a competent engineer could reasonably
make the other way (for example "Concepts, not vtables", or "release the GIL
in kernels"); (ii) is domain how-to the model would not otherwise apply here
(for example `tifffile` instead of `rasterio`, or `always_xy=True`); or (iii)
points at a statement. Everything else goes: advice the model follows anyway
(RAII, "use `auto` where it reads better", "never `eval`"), constraints that
are restated rather than pointed at (1.5), and aspirations the code does not
bear out. Each line is checked against the tree before it is kept, per A1.
Two examples, from `git grep` on this branch: `modern-cxx/SKILL.md` prescribes
C++20 coroutines (`std::generator`, `co_yield`) for sampling, and no file
under `include/`, `src/` or `bindings/` uses one. `geospatial-data-formats`
and `python-development` prescribe `ijson`/`ujson`, `iterparse` and
`pydantic-settings`, and nothing in `src_python/` imports them. These may be
intentions rather than errors, which is for the trim to establish.

---

## 6. Open harness items, and defects found on the way

### 6.1 R-B, the per-persona path guard: generic, driven by the profile

Place it in the generic layer as `.claude/hooks/guard_persona_paths.py`, one
`PreToolUse` hook on `Edit|Write|NotebookEdit|Bash`. The split between the
layers:

- **Generic: the role table**, stated by *path class*, never by path:

  | `agent_type` | may write |
  |---|---|
  | `tester` | `tests`, plus the amendment sections of `design` |
  | `developer` | `production` |
  | `perf` | `evidence` and its own tool |
  | `architect` | `design` and `docs/research` |
  | `reviewer` | nothing; its tools are already read-only |
  | `orchestrator` | the ledger only (6.2) |
  | no `agent_type` (main session) | anything, but `ask` on `production` and `tests` |
  | any subagent | never `.md`/`.txt` under `profile.channels.forbidden_dirs` |

- **Profile: the classes**, as globs in `profile.paths`.
- **Decision:** `deny` inside a subagent, `ask` in the main session, as the
  audit (§7 of `docs/retrospectives/2026-09-29-orchestrator-and-hooks-audit.md`)
  argues. A denied subagent hands the work back; an `ask` inside a subagent
  would wait for an owner who may be away.
- **Precondition:** the hooks reference documents `agent_type` inside
  subagents (`hooks`, common input fields), but nobody has seen it live here.
  Step 1 of the build is a logging hook in one spawned `@tester`. That changes
  `.claude/settings.json`, so it needs Ola's fresh yes.
- **This depends on the carrier.** Under a plugin, `agent_type` would be
  `harness:tester`. That is one more reason for option B.

### 6.2 The orchestrator ledger (Ola's idea)

Generic, optional, and it depends on R-B, because only the guard can make
"Write, but only to one file" true.

- **Who writes it.** The dispatcher. Both the local rules and the user memory
  say the main session is not the orchestrator persona, but the main session
  is what dispatches. "Orchestrator-only" therefore has to become
  "dispatcher-only": the guard allows the main session (no `agent_type`) and
  `orchestrator`, and denies every other `agent_type` both Read and Write. An
  alternative is to run the main session as `claude --agent orchestrator`,
  which makes the hook see `agent_type = orchestrator` (`hooks`, common input
  fields); that is Ola's call (§8).
- **What it holds.** One line per dispatch, appended: time, persona, ask,
  expected file, outcome. It does not replace `.claude/current-task/`. That
  directory holds what is in flight and is deleted when the work lands; the
  ledger is the history, and it is read for self-inspection.
- **Where it lives.** It is tracked, one file per branch
  (`.claude/ledger/<branch>.md`), and committed with the branch. It survives a
  change of machine, can be audited in git, and does not conflict at merge
  because each branch writes its own file.

### 6.3 Defects found while classifying

These are recorded for the migration steps, not fixed here.

1. **Contradiction.** `docs/increments/README.md@3245958:126` says "Independent suites
   run as parallel agents". `docs/PRINCIPLES.md@3245958:162-168` (D1) says "Dispatch
   serially. Disjoint files do not make parallel dispatch safe." Only one of
   them can be the rule.
2. **A claim the docs now contradict.** `.claude/REQUIRED-READING.md@3245958:48-53` says the
   `skills:` frontmatter "does not reliably preload them". The current docs
   say preloading injects the skill's full content (`sub-agents`, supported
   frontmatter, `skills`). That needs a one-minute live re-test, per the
   rule's own history. If preload works, a generic persona could preload its
   profile skill instead of carrying a profile region.
3. **The governed set misses rule files.** `guard_governance.py`'s `GOVERNED`
   and `GOVERNED_PREFIXES` do not include `testing.md` or `.claude/skills/`,
   although the 2026-09-30 sweep treats both as rule files.
4. **History remains in code docstrings.** `guard_governance.py` (d1db597),
   `gates_after_commit.py` (dated activation findings), `check_citations.py`
   (the 2026-09-17 sweep), `session_state.py` ("retrospective rule 5, Ola")
   and `check_prohibited_deps.py` ("the reason this project exists") all carry
   incidents in their docstrings. The authority quotes in
   `check_prohibited_deps.py`'s table are data, and they stay.
5. **Line citations into rule files.** Six citations in
   `docs/retrospectives/2026-09-29-orchestrator-and-hooks-audit.md` point into
   `REQUIRED-READING.md` as it stood before the sweep (lines 194, 208, 261,
   271, 302 and 320 of that file). `python3 tools/check_citations.py` listed
   all six as at-risk while the sweep was unmerged. Now that the sweep is on
   `master` it reports nothing, because the tool can only flag a change while
   that change is on a branch. Read against today's file, some of the numbers
   still land on the right text (`.claude/REQUIRED-READING.md@3245958:100-126`, `.claude/REQUIRED-READING.md@3245958:118-126`), and some are off by
   a line: `.claude/REQUIRED-READING.md@3245958:49-54` now starts one line into the skills section at 48, and
   `.claude/REQUIRED-READING.md@3245958:110-112` misses the start of the "fresh yes" list at 109. Pinning (step 2)
   is what makes such citations durable. Three more in the same file point into
   hooks whose docstrings step 3 will shorten: `.claude/hooks/guard_governance.py@180b9ac:24-25`,
   `.claude/hooks/guard_governance.py@180b9ac:56-59`, and `.claude/hooks/guard_push.py@180b9ac:20-30`. Four point into rule files the
   migration edits: `docs/increments/05-noder.md:992` → `testing.md@3245958:301`;
   `docs/increments/05b-noder-driver.md:411` →
   `docs/increments/README.md@3245958:65-67`; and two in code comments,
   `tests/cpp/property/noding_generators.h:5` and
   `tests/cpp/property/prop_noding_no_crossings.cpp:98`, both → `testing.md@3245958:220`.
   That makes 13. The two in code comments are invisible to
   `check_citations.py`, which scans only `.md` (`SCAN_SUFFIXES`), and to the
   `*.md`-only grep this design first used. They were found with
   `git grep -nE '(PRINCIPLES|testing|CLAUDE|README)\.md:[0-9]'`, which covers
   every suffix.
6. **Subagents get the session's snapshot of `CLAUDE.md`** (5c). A subagent
   receives the `CLAUDE.md` hierarchy as the main session loaded it at its own
   start, which may be stale against every checkout. Step 5a records this
   mechanism in the generic "harness" section, next to the existing note that
   `SessionStart` reads the main checkout. The rule to write is: "after a rule
   file changes, start a new main session before relying on
   `CLAUDE.md`-delivered text in a spawned persona; a persona reads the rules
   through its pointer".

---

## 7. Migration plan for rasputin

Each step is one PR and leaves rasputin working. Its gates stay green, and
`python3 tools/check_citations.py` reports no *broken* citation. Every
*at-risk* entry is re-read as a quotation, as `REQUIRED-READING` requires.
Where a step edits a rule file that something cites by line, it first
converts that citation (the "Citations" column). A prose step touches rule
files, so `guard_governance` asks Ola on each file. A code step goes through
the loop: `@tester` red, `@developer` green, `@reviewer`. LOC figures are
**estimates** of production lines as `CLAUDE.md` §2 counts them.

| # | Step | Files | Estimate | Citations |
|---|---|---|---|---|
| 1 | **One statement per rule** (prose, early, small). Apply the "Fix" column of 1.5, except the PRINCIPLES rows, which step 5 covers. Rename E4 so its heading states no number. Resolve defect 1 once Ola rules | `.claude/agents/{architect,developer,orchestrator,reviewer,tester,perf}.md`; `.claude/skills/{computational-geometry,geospatial-data-formats,modern-cxx}/SKILL.md`; `CLAUDE.md` line 9 (edited in place, so the §N numbering is untouched); `testing.md` (the floor restated at 18, edited in place so the line count stays the same, and at 283-288); `docs/PRINCIPLES.md` (E4 heading); `docs/increments/README.md` (defect 1) | 0 LOC; 1 round + `@reviewer` | First convert `05-noder.md:992` → "`testing.md`, *Test layout conventions*", because the edits to `testing.md` at 283-288 may shift line 301. Step 1 does **not** move `testing.md@3245958:220`: line 18 is edited in place and every other edit is below 220. So the two `tests/cpp/property/` comment citations stay correct without touching `tests/`, which would need `@tester`; step 2 converts them. `docs/increments/README.md@3245958:65-67` is cited, but step 1 only edits below it. Nothing cites the agents or skills by line (checked with `git grep`, all suffixes) |
| 2 | **Citation tooling.** `check_citations.py` resolves `path@<rev>:<line>` (generalising the `legacy-archive` case) and ID citations (`file §ID`), and warns on unpinned line citations into governed files. Widen `SCAN_SUFFIXES` beyond `.md` to the code comments of `.py`, `.h`, `.hpp`, `.cpp`, `.cmake` and `CMakeLists.txt`. Then convert the other 12 citations of defect 5 (step 1 already converted the 13th): pin the historical ones to the commit each was written against; replace the two living ones in `tests/cpp/property/` with the heading citation "`testing.md`, *Frameworks*". `@tester` makes that edit, as part of the red step | `tools/check_citations.py`, `tests/python/test_check_citations.py`, `docs/retrospectives/2026-09-29-orchestrator-and-hooks-audit.md`, `docs/increments/05b-noder-driver.md`, `tests/cpp/property/noding_generators.h`, `tests/cpp/property/prop_noding_no_crossings.cpp` | ~80-120 LOC | After this step, every line citation into a rule file, in prose or in a code comment, is pinned or reported by the checker. The widened scan is what makes that claim checkable: run the checker after planting an unpinned `testing.md:<n>` in a `.cpp` comment, and it must report it |
| 3 | **History out of code docstrings** (defect 4). Move it verbatim into the 2026-09-30 history retrospective | `.claude/hooks/*.py`, `tools/{session_state,check_citations,check_prohibited_deps}.py`, `docs/retrospectives/2026-09-30-rule-file-history.md` | 0 LOC (docstrings do not count); 1 round | Safe after step 2, which pinned the audit's citations into the hooks |
| 4 | **The profile.** `.claude/profile.toml` and a stdlib reader, `tools/harness_profile.py`. `session_state`, `gates_after_commit`, `guard_governance` (adding `testing.md` and `.claude/skills/`, defect 3) and `check_citations` read it. Add the single-statement gate from 1.2 as a fast gate. `CLAUDE.md` §4 points at `[[gates]]` | the above, their tests, and `CLAUDE.md` §4 (edited in place) | ~200-250 LOC | `test_session_state.py` keeps passing unchanged, since `owner = "Ola"` renders `ASK OLA`. `.claude/settings.json` is not touched, so no fresh yes is needed |
| 5a | **Rule-text split, part 1.** Move `REQUIRED-READING`'s profile sections (skill names, stale artifacts including A5's `.pyc`, data paths) and the README's Acceptance and template-cross-product paragraphs into a new section 5 of `CLAUDE.md` (to be titled "Project rules"), appended at the end. Leave the README heading as a pointer. Add the note for defect 6. Re-test the `skills:` preload (defect 2) | `.claude/REQUIRED-READING.md`, `docs/increments/README.md`, `CLAUDE.md` (append only), and the files that cite "Acceptance: an increment that touches refine or mesh code" by name (`git grep -n 'Acceptance: an increment'`) | 0 LOC; 1 round | §1-§4 unchanged. Name-citations updated in the same PR |
| 5b | **Rule-text split, part 2.** Turn PRINCIPLES into an index (§4, after Ola rules). Add `harness:begin/end` regions to the personas and `CLAUDE.md`. Drop the "kept for citations" notes. Trim the skills per §5d, checking each line against the tree | `docs/PRINCIPLES.md`, `.claude/agents/*.md`, `CLAUDE.md`, `.claude/skills/*/SKILL.md` | 0 LOC; 1-2 rounds | IDs unchanged, so `guard_push.py`'s "E1" still resolves |
| 6 | **Extract the harness repo.** Build the copier template from rasputin's generic files and regions, with lock generation and `tools/check_harness.py`. Smoke test: `copier copy` into an empty repo with a toy profile, then check that the recap prints and that the guards ask on a planted push and a planted rule edit | a new repository (a publishing act: Ola) | ~150 LOC there | rasputin untouched |
| 7 | **rasputin adopts.** Run `copier copy --vcs-ref v0.1` over the tree. Acceptance: the diff is empty apart from `.copier-answers.yml`, `.claude/harness.lock` and `tools/check_harness.py`, which joins the fast gates | those three files and `.claude/profile.toml` (the gate entry) | ~100 LOC (vendored gate) | no rule text changes if the acceptance holds |
| 8 | **R-B guard** (6.1). First a probe: a hook that logs `agent_type` in one spawned `@tester` (a `.claude/settings.json` edit, so it needs a fresh yes). Then `guard_persona_paths.py`, with a suite that plants each forbidden write for each persona. Build it in the harness repo, and `copier update` it into rasputin | harness repo; then rasputin's `.claude/hooks/`, `.claude/settings.json` and `profile.paths` | ~150 LOC | — |
| 9 | **Ledger** (6.2), optional, after step 8 | guard table, `.claude/ledger/`, the orchestrator's generic region | ~50 LOC | — |

Steps 1-5 are useful on their own even if steps 6-9 never happen. They leave
rasputin's rules deduplicated, parameterised and free of history, which was
the original ask.

---

## 8. Open questions for Ola

1. **Carrier.** Should the generic layer live in a new repository as a copier
   template (for example a private `expertanalytics/agent-harness`), vendored
   into each project? *Recommend yes* (§2). Creating the repository and
   choosing its licence are yours.
2. **Owner's name in generic text.** Should the text be rendered per project
   ("Ola", `ASK OLA:`), or say "the owner"? *Recommend rendering.* It keeps
   your preference to be named, and rasputin's text and tests unchanged.
3. **PRINCIPLES.md.** Should it become an index (ID, title, home, check), with
   "last exercised" replaced by a `git grep` of the ID and retirement made
   event-driven? *Recommend yes* (§4). The stored field has not been updated
   since increment 5c.
4. **Parallel suites or one agent at a time** (defect 1). *Recommend D1 wins*
   and the README's "Independent suites run as parallel agents" goes, because
   that matches your standing instruction to dispatch serially.
5. **Gemini (c).** Should the persona pointer to `REQUIRED-READING.md` stay?
   *Recommend keeping it as one line.* §5c shows that a subagent receives
   the `CLAUDE.md` the main session loaded at its own start, which may be
   stale against every checkout. Only the pointer reads the rules as they are
   on the branch now.
6. **Ledger.** Should it be tracked, one file per branch, written only by the
   dispatcher (the main session and `@orchestrator`), and readable by the main
   session so it can relay it to you? Or would you rather run the main session
   as `claude --agent orchestrator`, so the guard sees an explicit identity?
   *Recommend the first*, which keeps "the main session is not the orchestrator
   persona" true.
7. **R-B timing.** Should the guard be built after step 4, since it needs
   `profile.paths`, starting with the `agent_type` probe (which needs your yes
   on `.claude/settings.json`)? *Recommend yes.*
8. **`@perf` in the generic layer.** Should the generic layer include it as an
   optional skeleton, off by default, with evidence and like-for-like
   comparison generic and the tool and triggers from the profile? *Recommend
   yes*, since every similar project will have a performance question.


## 9. Ola's rulings (2026-09-30)

All eight open questions are answered as recommended:

| # | Question | Ruling |
|---|---|---|
| 1 | Carrier | A copier template in a new repository (for example a private `expertanalytics/agent-harness`), vendored into each project. Creating the repository and choosing its licence remain Ola's acts, taken at step 6 |
| 2 | Owner's name | Rendered per project from the profile ("Ola", `ASK OLA:` in rasputin) |
| 3 | `docs/PRINCIPLES.md` | Becomes an index; "last exercised" is found by `git grep` of the ID, and retirement is event-driven (step 5b) |
| 4 | Parallel suites vs one agent at a time | One agent at a time wins; the README's "Independent suites run as parallel agents" goes (step 1) |
| 5 | Persona pointer to `REQUIRED-READING.md` | Kept, as one line |
| 6 | Ledger | Tracked, one file per branch, written only by the dispatcher, readable by the main session |
| 7 | R-B timing | After step 4, starting with the `agent_type` probe. The probe's `.claude/settings.json` edit still needs Ola's yes when it comes |
| 8 | `@perf` | In the generic layer as an optional skeleton, off by default |
