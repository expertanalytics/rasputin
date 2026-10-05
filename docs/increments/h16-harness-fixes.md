# Harness h16: guard fixes, a line counter, a scratch copy, brief fixes

Status: Ola ruled on §7 on 2026-10-05 (all three defaults) and on the afternoon questions (last section). PR A: green at `195b504` plus the test fix `6e52aeb` asked for by code review round 3, 287 net production lines by `tools/count_loc.py bc01cd8 a6b966e` against an estimate of 187 (the C++ scanner, the git plumbing and review's fixes were not priced); the rule lines of G5 and of R1 (but its *The harness* sentence, PR B's) written; code review round 4 asked for changes (two test docstrings), fixed in `a6b966e`; round 5 found nothing else; waits for Ola's yes to push, then CI. PR B: red `9cf533d` on `worktree-h16b`, green in progress.

Ola approved the items on 2026-10-05 (the main session's summary of his
rulings, not his words). He said this is the last harness increment before
a freeze of a few days, so the design is kept to what each item needs and
§5 names what should drop out. What the labels mean:

| Label | What it is | PR |
|---|---|---|
| G1 | `guard_push.py` asks before `git fetch` into a named ref and before `git replace` writes | B |
| G2 | `guard_push.py` asks before a git or gh command it does not know (an alias from user config) | B |
| G3 | scratchpad repositories are ordinary (night proposal N3), in two parts, G3a and G3b; G3c is the push half | B |
| G4 | guard import shadowing: `tools/` last on `sys.path`, stdlib-named `tools/` files governed | B |
| G5 | rule files written through the shell: one brief line | A |
| T1 | `tools/count_loc.py`, the one counter for `CLAUDE.md` §2 (night proposal N1) | A |
| T2 | `tools/scratch_copy.py <rev> <dir>` (night proposal N2), and evening P9 | A |
| T3 | `tools/brief.py`: note-file names that cannot collide; root Markdown files in `@architect`'s limit | A |
| R1 | rule lines: Monitor not `sleep`; a fallback kept for the last hour | A |

"The scratchpad" is the session's temporary directory,
`/private/tmp/claude-<uid>/<project>/<session>/scratchpad/`. "Governed"
means `guard_governance.py` asks before a write (`governed()` in that file).
Sources not on `master` are named by branch and commit:
`worktree-retro-1004c` holds the night retrospective
(`docs/retrospectives/2026-10-05-night.md`, a79f2d9) and the evening one
(`docs/retrospectives/2026-10-04-evening-merges-h11-h13-29.md`, c0e3451);
the brief said the evening one is on `master`, and it is not (`git branch -a
--contains c0e3451` lists only `worktree-retro-1004c`). `worktree-h12-design`
holds `docs/increments/h12-prose-fast-lane.md`.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature.* The sources each item rests on, read 2026-10-05:

- git 2.55, `git-config(1)`, `alias.*`: "To avoid confusion and troubles
  with script usage, aliases that hide existing Git commands are ignored
  except for deprecated commands." So an alias can carry any name that is
  not a current command, and can take a deprecated one (checked here:
  `git -c alias.whatchanged='rev-parse --short HEAD' whatchanged` printed a
  hash). G2 rests on it.
- `git-config(1)`, `url.<base>.pushInsteadOf`: "Any URL that starts with
  this value will not be pushed to; instead, it will be rewritten to start
  with <base>, and the resulting URL will be pushed to." G3c rests on it.
- `git-replace(1)`: "Typing "git replace" without arguments, also lists all
  replace refs"; `-l`/`--list` lists. G1 treats those forms as reads.
- `python3(1)`, `-I`: "In isolated mode sys.path contains neither the
  script's directory nor the user's site-packages directory"; `-P`: "Don't
  automatically prepend a potentially unsafe path to sys.path such as [...]
  the script's directory". Python's *The initialization of the sys.path
  module search path*: "The first entry in the module search path is the
  directory that contains the input script". G4 rests on these.
- Anthropic, *Writing effective tools for agents* (cited by the night
  retrospective for N1 and N2) and *Making Claude Code more secure and
  autonomous with sandboxing* (for N3); not re-read here.

What differs: the h12 design scrubs user config
(`GIT_CONFIG_GLOBAL=/dev/null GIT_CONFIG_NOSYSTEM=1`) so that it reads only
the governed `.git/config`. That is right for a tool that must not be
misled (T1, T2, G3a's lookup) and wrong for a guard that must predict what
a command will do, because the command reads the user config the scrub
hides. G2 and G3c therefore do not scrub; §2 says where each applies.

*Legacy.* Nothing. `git grep -l -i -e 'insteadOf' -e 'sys.path.insert' -e
'git archive' -e 'count_loc' legacy-archive -- legacy` returned no file
(exit 1); the legacy tree had no harness.

## 2. The items

Each item: what changes, the incident it prevents, the red test, its size.
Sizes are production lines by `CLAUDE.md` §2 (net), estimated; tests are
excluded.

### G1. `git fetch` into a named ref, and `git replace`

**Change.** In `guard_push.py`'s `segment_why`:

- `git fetch` or `git pull` with a positional other than the first (the
  repository) that contains `:` asks: "fetch writes a named ref". That is a
  refspec with a destination (`. <sha>:refs/remotes/origin/master`,
  `+a:b`). `git fetch origin`, `git fetch -q origin master` pass, as today.
- `git replace` asks ("replace refs change what git reads for an object")
  unless it is a listing: no arguments, or `-l`/`--list` (with or without a
  pattern), or `--format=…` with those. Every other form (create, `-d`,
  `-f`, `--edit`, `--graft`, `--convert-graft-file`) asks.

**Incident.** h12 design review round 2 (recorded in dfe3ad7 on
`worktree-h12-design`): in a scratch repository,
`git fetch . "<sha>:refs/remotes/origin/master"` repointed `origin/master`
at a planted commit, and `git replace` made `git archive` serve a planted
file; neither was asked. h12's own text (§3.3, "These checks do not depend
on any guard change") says the guard change makes the routes harder, not
h12's checker correct; that is all G1 claims.

**Red test** (`tests/python/test_guard_push.py`, parametrised over argv):
asks for `git fetch . abc:refs/remotes/origin/master`,
`git fetch origin +master:refs/heads/x`, `git pull . a:b`,
`git replace HEAD HEAD~1`, `git replace -d x`, `git replace --graft a b`,
`git replace --edit a`; passes `git fetch origin`, `git fetch -q origin master`,
`git fetch --all`, `git replace`, `git replace -l`, `git replace --list 'a*'`.
Pinned false positive: `git fetch --depth 1 git@github.com:a/b.git` asks
(the URL is the second positional after `1`); `--depth=1` does not.

**Size.** About 12 lines.

### G2. Git and gh aliases: user config changes what the guard sees

**Change.** `guard_push.py` judges a git command by its subcommand word, so
an alias hides what runs. Probe, run on this branch's base (bc01cd8):
`publishes(['git', '-c', 'alias.p=push', 'p', 'origin'])` and
`publishes(['git', 'p', 'origin'])` both return `[]`, and `segment_why`
returns `None`: an alias written into `~/.gitconfig` (not governed) or given
with `-c` pushes unasked. The fix: a git subcommand that is not a current
git command asks, with the reason "a git alias or a command the guard does
not know; it cannot see what it runs". "Current" is the output of
`git --list-cmds=main` less `git --list-cmds=deprecated` (181 and 2 names
with git 2.55), read once per hook run. Those two calls are not scrubbed:
they list commands, which no config changes. The same for gh: `words[1]`
outside a fixed set of gh's top-level commands asks. The set is what
`gh help` lists with gh 2.101 less its alias `co` (which a user can
redefine): `auth browse codespace discussion gist issue org pr project
release repo skill cache run workflow agent-task alias api attestation
completion config copilot extension gpg-key label licenses preview ruleset
search secret ssh-key status variable`. A gh extension is then asked
about too.

**Incident.** h12 design review round 3 (recorded in 2d77b6c on
`worktree-h12-design`): a user-level config file changed what git resolved
(`url.<x>.insteadOf` sent `ls-remote` and `fetch` elsewhere) with no
guarded command. That round's subject was the h12 checker, not the guards;
the alias route above is the same class turned on `guard_push.py`, found
while writing this design, with no incident of its own.

**Red test** (`test_guard_push.py`): asks for `git -c alias.p=push p origin`,
`git p origin`, `git whatchanged` (deprecated, so aliasable), `gh pm 12`,
`gh co 12`; passes `git status`, `git log -1`, `git worktree list`,
`gh pr view 12`, `gh api repos/x` (a GET). A test that `git --list-cmds=main`
contains `push`, so a git without `--list-cmds` fails loudly rather than
asking on everything.

**Size.** About 14 lines.

### G3. Scratchpad repositories are ordinary (N3)

Night retrospective §6, N3 (a79f2d9). Three parts; the third should drop
(§5).

**G3a, governance.** `governed(path)` returns False for an absolute path
whose real path (`os.path.realpath`, so `/tmp` resolves to `/private/tmp`
and a symlink to the repository's `.git` is followed) lies under a
scratchpad: the pattern
`^/private/tmp/claude-\d+/[^/]+/[^/]+/scratchpad(/|$)`. Any session's
scratchpad, not only the current one: all are temporary, and nothing the
harness reads lives there. A relative path is judged as today, since the
guard does not track `cd`.

**G3b, local git writes in a scratch repository.** `segment_why` returns
None for the local writes it now asks about (`config`, `remote`,
`symbolic-ref`, `update-ref`, G1's fetch and replace) when all hold:

- the git call has `-C <dir>` with `<dir>` absolute;
- `git -C <dir> rev-parse --absolute-git-dir --git-common-dir`, run with
  `GIT_CONFIG_GLOBAL=/dev/null GIT_CONFIG_NOSYSTEM=1`, succeeds, and both
  directories, resolved against `<dir>` and by `realpath`, lie under a
  scratchpad. This catches a scratch directory that is a linked worktree
  of the real repository (`git worktree add`), whose config, refs and
  replace refs are the real repository's;
- no `--git-dir`, `--work-tree`, `--global` or `--system` in the argv, and
  the command text contains no `GIT_` (a `GIT_DIR=` prefix overrides `-C`;
  `tools/shell_scan.py` strips leading assignments, so the text is checked).

`git config --file <path>` with `<path>` under a scratchpad (G3a's test)
passes too. Here the scrub is right: the lookup asks where the repository
is, and must not be steered by an `include` or other entry in a user file.

**G3c, a push to a scratch repository** (recommended to drop, §5). A push
passes only when the line parses to exactly one simple command, a
`git -C <scratch repository> push <dest> …` as in G3b; `<dest>` is written
literally as an absolute path or a `file://` URL whose real path lies under
a scratchpad (a remote name is never resolved); the argv has no `-c`,
`--repo`, `--receive-pack` or `--exec`; and
`git -C <dir> config --get-regexp '^url\.'` run **without** the scrub,
reading what the push will read, prints nothing. Any `url.*` entry at any
level voids the pass, because `pushInsteadOf` can send a scratch path to the
real remote. What it cannot see: a variable such as `GIT_CONFIG_GLOBAL`
exported by a shell startup file into the Bash tool's shell but not into the
hook's environment (h12 §5's residual).

**Incident.** The night of 2026-10-04/05, §2 of the night retrospective
(a79f2d9): four of seven `guard_unattended` refusals were scratch
repositories in the scratchpad (22:58:22 `git config user.*`,
`git remote add` and a push to a local bare repository, agent `a5b4f008`;
23:04:28 an append of `[url] insteadOf` to a scratch `.git/config` and
23:05:27 a push to a `file://` scratch remote, agent `a3cc5c9a`), and two of
them were then reached by another route. A refusal leaves no commit; the
retrospective's table is the record.

**Red test.** `test_guard_governance.py`: a Write of
`<scratchpad>/r/.git/config` passes; of `<scratchpad>/copy/CLAUDE.md`
passes; of `<scratchpad>/link/.git/config` where `link/.git` is a symlink to
a real repository's `.git` asks; of `.git/config` (relative) asks. A Bash
`cat >> <scratchpad>/r/.git/config` passes. `test_guard_push.py`, with a
real scratch repository made by the test under `tmp_path` and the pattern
pointed at it (a module constant the test patches): `git -C <it> config
user.name x` and `git -C <it> remote add o /x` pass; the same in a linked
worktree of a second repository asks; with `GIT_DIR=/x ` in front asks;
`git config user.name x` (no `-C`) asks. If G3c stays: `git -C <it> push
<scratch bare> master` passes; with a `GIT_CONFIG_GLOBAL` file (set in the
hook's environment) holding `url.<outside>.pushInsteadOf = <scratch bare>`
asks; `git -C <it> push origin master` (a remote name) asks;
`git -C <it> config url.x.insteadOf y && git -C <it> push <bare>` asks.

**Size.** G3a about 6 lines (a `scratchpad` helper, below, and the
exemption); G3b about 20; G3c about 25.

### G4. Guard import shadowing

**Change.** Two parts, both needed:

1. Every `sys.path.insert(0, <tools>)` becomes `sys.path.append(<tools>)`:
   `.claude/hooks/guard_push.py@bc01cd8:32`,
   `.claude/hooks/guard_governance.py@bc01cd8:44`,
   `.claude/hooks/guard_spawn.py@bc01cd8:80`,
   `.claude/hooks/guard_unattended.py@bc01cd8:76`,
   `tools/session_state.py@bc01cd8:48`, `tools/brief.py@bc01cd8:28`,
   `tools/away.py@bc01cd8:32` (seven lines changed, net 0). The standard library
   then wins over any file in `tools/`. The hooks' own directory,
   `.claude/hooks/`, stays first, and is governed by prefix.
2. A script run as `python3 tools/x.py` (the `SessionStart` hook, the main
   session's `brief.py`) has `tools/` first by Python's own rule, which no
   `append` changes. So `governed()` also returns True for a path whose
   component after a `tools` component, cut at the first `.`, is in
   `sys.stdlib_module_names` (`tools/ast.py`, `tools/json/__init__.py`,
   `tools/subprocess.cpython-314-darwin.so`). G3a's exemption comes first,
   so a probe in a scratch copy is not asked about.

**Incident.** h12 PR A review round 2 (recorded in 6403582 on
`worktree-h12-design`): a planted `tools/json.py` replaced the stdlib module
inside `tools/ci_changes.py`, fixed there with `python3 -I`. The reviewer's
lesson carried it to the gates; `@architect` then found it reaches the
hooks (main session transcript `806b4380`, line 7527, the `ASK OLA:` line
written at 07:00:27). Checked again on this branch's base in a scratch copy
of the two guards and `tools/`: with `tools/dataclasses.py` raising at
import, `guard_governance.py` given a Write of `CLAUDE.md` and
`guard_push.py` given `git push` both exited 1 with a traceback and no
decision, which Claude Code treats as a non-blocking error; with
`sys.path.append` in both, each printed its `ask`. No file in `tools/` or
`.claude/hooks/` is named after a stdlib module on any branch (the main
session's check, same transcript, line 7548).

**Red test.** `test_guard_governance.py`: `governed()` is True for
`tools/ast.py`, `tools/json/__init__.py`, `/abs/wt/tools/typing.py`; False
for `tools/count_loc.py` and `docs/x/ast.py`. A test per hook that copies
the hook and `tools/` into `tmp_path`, plants `tools/dataclasses.py` raising
at import, feeds a governed Write (or a push) and requires an `ask` on
stdout. A repository test: no file in `tools/` is named after a stdlib
module, so one approved by mistake still fails CI.

**Size.** About 6 lines net.

### G5. Rule files written through the shell

**Finding: the reported bypass did not happen.** The h12 `@architect`
(agent `a329f25e`, the run that committed 6403582) handed back "I wrote two
governed files with a Python heredoc, which the governance guard does not
see." Its transcript shows otherwise: the heredoc that wrote
`docs/increments/README.md` and `.claude/agents/reviewer.md` (line 344,
06:49:26 UTC) drew `guard_governance.py`'s `ask` naming both files (line
346), which was answered yes 71 seconds later (line 348); its reverts were
asked about too (lines 355, 364). Replayed on this branch's base,
`judge_bash` on that command returns `ask` for both files.

**What is true.** The Bash arm reads string literals, so a path built from
parts passes: on this branch's base, `judge_bash` returns None for
`Path('docs', 'increments', 'README.md').write_text(…)`,
`(Path('docs/increments') / 'README.md').write_text(…)`,
`open('CLAU' + 'DE.md', 'w')`, and `python3 /tmp/w.py`. The guard's own
docstring calls itself "a tripwire, not a sandbox"; no static reading of a
program closes this.

**Change (bound, not closed).** One line in `.claude/briefs/common.md`:
"Write a rule file with Edit or Write, never through the shell." The Edit and
Write arm judges the exact path, so the rule moves rule-file writes to the
arm that cannot miss. Prose, no red test. (G5b, joining literal path parts
in `tools/shell_scan.py`'s `candidates`, about 15 lines, is left out, §5.)

### T1. `tools/count_loc.py` (N1)

**Interface.**

    python3 tools/count_loc.py <base> [<head>]

`<head>` defaults to `HEAD`. The old side is `git merge-base <base> <head>`,
so `<base>` may be `origin/master` and the count is the PR's, as
`git diff <base>...<head>`. Output, tab-separated, one line per counted
file, then a total, then one `not counted:` line per changed file it skipped
with the reason (`tests/`, `docs/`, `not code`):

    src_python/tin_engine/gauge.py	117	0	117
    total	744	53	691
    not counted: tests/python/test_gauge.py (tests/)

Exit 0; exit 2 with one stderr line starting `count_loc:` on a git error.
Committed revisions only; the working tree is not read.

**Rule, from `CLAUDE.md` §2.** Added lines are the `+` ranges of
`git diff -U0` hunks, judged in the file at `<head>`; removed lines the `-`
ranges, judged at the old side. A line counts unless it is blank, a comment,
a docstring, or inside a raw literal's body. Per kind:

- Python (`.py`): `tokenize`, as the 29 PR 2 round-2 record states its
  method (cc52b8f): a line counts if it holds a token other than comments,
  newlines and indentation; docstrings (the first statement of a module,
  class or function, if a string; found with `ast`) do not count on any of
  their lines; any other string spanning lines counts on its first line
  only. A file that does not tokenize counts every non-blank line, with a
  stderr warning.
- C++ (`.h`, `.hpp`, `.cpp`, `.cc`, `.cxx`): a small scanner; a line counts
  if it has a character outside `//` and `/* */` comments and outside a raw
  string's body. A raw string `R"d(…)d"` counts on its opening line; its
  later lines count only for code after the closing `)d"`.
- CMake (`CMakeLists.txt`, `.cmake`) and shell (`.sh`): blank and
  `#`-comment lines do not count.
- Not counted: anything under `tests/` or `docs/`, and files of no kind
  above (Markdown, YAML, TOML, JSON, data). The `not counted:` lines show
  them, so a reviewer sees what was left out (question 2).

Every git call runs with `GIT_CONFIG_GLOBAL=/dev/null`,
`GIT_CONFIG_NOSYSTEM=1`, `GIT_NO_REPLACE_OBJECTS=1`, `--no-ext-diff
--no-color`: a user setting such as `diff.noprefix` changes the headers the
parser reads (h12 §3.3 names it). Renames are git's default detection.
The repository's own config is not neutralised by those variables, so the
diff also passes its settings explicitly: `--src-prefix=a/ --dst-prefix=b/`
and `--inter-hunk-context=0` (review round 1), and `--no-relative` (review
round 2: `diff.relative=true` in `.git/config`, run from a subdirectory,
would otherwise count only that subdirectory).

**Blueprint.** Pure functions, each tested alone, and one thin `main`:

    def kind_of(path: str) -> Kind | None
    def counted_lines(text: str, kind: Kind) -> frozenset[int]   # 1-based
    def hunks(diff: str) -> dict[str, tuple[list[int], list[int]]]  # path -> (new +lines, old -lines)
    def tally(...) -> list[Row]                                    # joins the three
    def main(argv) -> int                                          # git I/O only here

**Incident.** Night retrospective, lesson 3 (a79f2d9): eighteen runs wrote
their own counter; two personas got 558 and 577 from the same tree
(`a1fee960`'s lessons).

**Red test** (`tests/python/test_count_loc.py`): `counted_lines` cases per
kind (comment-only, blank, docstring of one and of three lines, a
multi-line non-docstring string, a C++ raw string over three lines with
`;` after it, `/* */` across lines with code after it); `hunks` on a
hand-written diff with a new file, a deletion and two hunks; a temporary
repository whose two commits give known counts; that a `GIT_CONFIG_GLOBAL`
file with `diff.noprefix=true` does not change the output; and two
recorded counts from real history, which need the full clone CI already
has (`fetch-depth: 0`): `count_loc.py 529613a 193079d` totals 589, 20,
569 (29 PR 2, review round 3), and `count_loc.py 9e666f4 9bb1723` totals
744, 53, 691 (29 PR 4, ef403bd's message and the round-6 review). If the
tool and a recorded count disagree, `@tester` reports which line differs
before anyone changes either.

**Size.** About 120 lines.

### T2. `tools/scratch_copy.py` (N2), and P9

**Interface.**

    python3 tools/scratch_copy.py <rev> <dir>

1. Refuses (exit 2, one stderr line) if `<dir>` exists and is not empty, or
   lies inside the worktree it is run from or inside the main checkout (the
   work tree of the common git dir, which holds `.claude/worktrees/`).
2. `git archive <rev> | tar -x -C <dir>`, with T1's three variables, so the
   copy is not a git work tree.
3. Copies the built `_core*.so` from the running worktree's
   `.venv/lib/python3.*/site-packages/tin_engine/` into
   `<dir>/src_python/tin_engine/`. None found: a stderr warning (Python-only
   suites still run). If `git diff --quiet <rev> HEAD -- include src
   bindings CMakeLists.txt` says the C++ differs, a stderr warning that the
   copied `_core` was built from `HEAD`.
4. Writes `<dir>/.scratch_copy/sitecustomize.py`, which drops the editable
   finder from `sys.meta_path`, and prints one line on stdout, the command
   that runs `pytest` against the copy: `cd <dir> && PYTHONPATH=<dir>/.scratch_copy:<dir>/src_python`
   and the worktree's `.venv/bin/python -c` with a program that only calls
   `pytest.main` on the remaining arguments (`tests/python/` by default;
   replace it to run some tests only). Because `PYTHONPATH` is inherited,
   the drop runs at every interpreter start, after the `.pth` file installed
   the finder, so `tin_engine` comes from the copy in pytest's process and
   in any Python child a test starts. Review round 1 caused this: the first
   version did the drop inside the `-c` program only, as
   `docs/benchmarks/2026-10-02/basin-memory-probe/blockprobe.py:23` does,
   and a child process reloaded the finder and imported the worktree's
   code. Limits: a child whose environment replaces `PYTHONPATH` (one built
   from scratch, as `tests/python/test_io_geotiff.py:1271` does) or that
   starts Python with `-I` or `-E` (both ignore `PYTHONPATH`) still imports
   the worktree's code, so a mutant run there can report a false survivor.

It also serves "run the new tests against the code before the change"
without `git stash`.

**P9** (evening retrospective, c0e3451): `test_hook_is_executable_in_the_checkout`
in `tests/python/test_settings_wiring.py` skips, with the reason stated,
when the tree is not a git work tree. `@tester`'s, in the red step.

**Incident.** Night retrospective lessons 1 and 2 and §3 (a79f2d9): six
`@tester` runs each built their own way round the editable finder; the
seven `test_settings_wiring` failures in a `git archive` copy (round 2's
"8 failed", c7d427a); five `git stash` uses in the shared stash.

**Red test** (`tests/python/test_scratch_copy.py`): the refusals; the copy
has no `.git`; the printed command, run with `tests/python/` replaced by
one generated test file asserting `tin_engine.__file__` lies under
`<dir>`, exits 0 (this is the check that the finder is dropped; it fails
if the drop is removed); the C++-differs warning on a revision before a
C++ change.

**Size.** About 60 lines.

### T3. `tools/brief.py`

**Note-file names.** Today `<persona>-<HHMMSS>.md`. Two briefs printed in
the same second for the same persona in different worktrees got the same
name: `architect-094230.md`, for `h15-ci` and for `h16`, on 2026-10-05
(main session transcript `912df417`, lines 147 and 155; no commit). The
name becomes `<persona>-<worktree>-<HHMMSS>.md`, `<worktree>` the
worktree directory's name; if that file exists, `-2`, `-3`, … are appended
before `.md`. `.claude/REQUIRED-READING.md` states the format and changes
with it. **Red test** (`test_brief.py`): two briefs with the clock frozen,
for `h15-ci` and `h16`, name different files; with the file already there,
the next name ends `-2.md`. About 6 lines.

**Root Markdown files.** The item is the evening retrospective's P1 (c0e3451): root Markdown files
(`README.md`, `INSTALL.md`, `NOTICE.md`, `testing.md`,
`project_structure.md`, `auto_catchments.md`, `parallel_refinement.md`)
are in no persona's write limit, and Ola ruled on 2026-10-05 that
`@architect` may edit them (recorded in ef403bd's message). `WRITES` still
lacks them: this run's own brief prints `@architect`'s limit without them.
The change: `@architect`'s entry gains "root *.md files". P1's second half,
a test that every tracked path falls in some limit, needs the limits as
patterns rather than prose, about 30 more lines; it is left out (§5).
Question 1 asks Ola whether this is the item. **Red test**: `brief.py
--persona architect` prints a limit naming root Markdown files. About 1 line.

### R1. Rule lines

- **Monitor, not `sleep`.** `.claude/briefs/common.md` gains: "Wait for a
  background run with the Monitor tool, never with `sleep`." Incident: night
  retrospective §2 (a79f2d9), three refusals by Claude Code's own `sleep`
  check (`ac6626cc` 00:48:08, `a9fbff4d` 01:54:49, `a79b3df8` 02:50:22),
  each an agent sleeping and then reading its log.
- **A fallback for the last hour.** `.claude/REQUIRED-READING.md`,
  *Unattended mode*, already requires a no-ruling fallback in `QUEUE:`
  before Ola leaves; it gains "and one such item stays queued for the
  window's last hour". Incident: night retrospective §1 and question 1
  (a79f2d9): from about 04:35 to 05:28 every queued item waited on Ola.
- **Pointers.** `tester.md`'s mutant paragraph replaces its `git archive`
  recipe and place list with "in a scratch copy made by
  `python3 tools/scratch_copy.py <rev> <dir>` in the session scratchpad,
  removed afterwards" (the night's cut C4, about 35 words fewer);
  `reviewer.md`'s LOC item and `CLAUDE.md` §2 name `tools/count_loc.py`;
  REQUIRED-READING's *The harness* names G1 to G4 in one sentence (PR B).

Prose; no red test. `@architect` writes them, in the PR that ships the
tool or guard they describe. `CLAUDE.md` changes, so the main session
restarts after PR A merges.

## 3. Shared pieces and boundaries

- **`tools/scratchpad.py`**, new, standard library only: the pattern and
  `def under(path: str) -> bool` (absolute paths only, by `realpath`).
  Imported by both guards (G3) and added to `GOVERNED`'s self-protecting
  set, as `shell_scan.py` is. About 10 lines.
- **`tools/count_loc.py`** joins `GOVERNED` too: it computes the arithmetic
  of a rule in `CLAUDE.md` §2, so a change to it is a change to the rule.
  `tools/scratch_copy.py` does not.
- The hooks stay pure apart from the git calls named above; the new git
  calls are G2's two command listings and G3b's (and G3c's) lookups, each
  with a fixed argv. The guards' `except Exception` → `deny` stays, so a
  failing lookup refuses rather than passes.

## 4. PR split and size

Under `CLAUDE.md` §2 (700 net per PR) everything fits one PR (about 250
lines, about 300 with G3c and G5b). It is split anyway, so the tools are
not held up by guard review rounds (h12's design took four):

| PR | Items | Production lines, about |
|---|---|---|
| A, tools | T1, T2 (+P9 test), T3, G5's line, R1 but its *The harness* sentence | 187 (count_loc 120, scratch_copy 60, brief.py 7) |
| B, guards | G1, G2, G3a, G3b, G4, `scratchpad.py`, the `GOVERNED` entries, R1's *The harness* sentence | 60 (G1 12, G2 14, G3a 6, G3b 20, G4 6, scratchpad 10, minus shared lines) |

Order: A first. Each PR runs red (`@tester`), green (`@developer`), review
(`@reviewer`). Neither touches refine or mesh code, so no `@perf` run. No
suite here is invariant-critical, so no mutation round. Writing the hooks,
`tools/` files in `GOVERNED` and the rule files is asked about at each
write, so both PRs are day work.

## 5. What should drop

- **G3c (the push to a scratch repository): drop, by default.** It is the
  one part that lets a push through the guard, its safety rests on the
  route list being complete, and h12's review found a new route in each of
  three rounds. The night's push refusals had a workaround that cost
  seconds (`git remote get-url --push`). Question 3.
- **G5b (joining literal path parts): drop.** The incident did not happen
  as reported, and what remains open after it stays open after it.
- **P1's coverage test: drop** from h16 (above).

Nothing else should drop: G4 is a live hole in the hooks, G1 and G2 are
small, and T1 and T2 remove the most repeated work of the night.

## 6. Residual, after h16

The guards remain tripwires. Not covered: programs that compute a path
(G5); a module planted in the user-writable Homebrew site-packages under a
`tools/` module's name, which with `append` would now win (h12 §5's
environment residual); shell startup files and `GIT_*` variables they
export; relative paths in a scratch repository (G3 needs `-C` and absolute
paths, and the refusal message says so).

## 7. Questions for Ola

1. **Which `brief.py` item did you approve as "a limit on root-file size it
   enforces"?** No source found states a size limit. Default: the
   root-Markdown item in T3 (your 2026-10-05 ruling that `@architect` edits
   root Markdown files, written into `brief.py`'s limits).
2. **Which files does the line counter count?** Default: code files (Python,
   C++, CMake, shell) outside `tests/` and `docs/`; workflow YAML,
   TOML and Markdown are listed as "not counted".
3. **Drop the push half of the scratchpad item (G3c) from h16?** Default:
   yes, drop it; scratch config and remote writes with `git -C` still pass.

## Ola's rulings

2026-10-05, Ola, verbatim: "yes, go with all three defaults". So: 1 the brief.py item is `@architect`'s write limit gaining root Markdown files (the earlier phrase "root-file limit" was the main session's shorthand; Ola: "I have no idea what a root-file limit would mean."); 2 the line counter counts Python, C++, CMake and shell outside `tests/` and `docs/`, and lists the rest as not counted; 3 G3c, the push to a scratchpad repository, drops out of h16.

## Review

### Round 1: `@reviewer`, code, early, PR A tools part, `bc01cd8..4db1eab`

`@reviewer`'s verdict, word for word:

**Code review, round 1 (early, tools part of PR A), 2026-10-05.** Range `bc01cd8..4db1eab` (design e9e3c46, rulings 093725b, red 11cee8e, partial green 4db1eab). Verdict: CHANGES REQUESTED.

LOC: 252 net production lines by `tools/count_loc.py bc01cd8 4db1eab`: `tools/count_loc.py` 195 and `tools/scratch_copy.py` 57. The workflow, the tests and the docs are listed as not counted. The 57 was recounted by hand, line by line. The estimate was 187 (count_loc 120, scratch_copy 60, brief.py 7). count_loc is 75 over, because the C++ scanner and the git plumbing were not priced; scratch_copy is 3 under. With brief.py's about 7 lines still to come, the total is about 259 against 700. The design names no split seam, and none is needed.

Packing: 11 `# fmt: skip` regions, 6 in count_loc and 5 in scratch_copy. Together they save 43 lines; ruff format would give 222 and 73. Each region keeps one call's arguments or one literal's entries on two lines instead of one per line, and each stays readable. They are packed for density only; none is needed for the ceiling.

The counter reproduces both recorded counts exactly:
- `529613a 193079d` gives 589, 20, 569.
- `9e666f4 9bb1723` gives 744, 53, 691.

On a range I chose, PR #163 (`45acf22^1..45acf22`), it gives 41, 20, 21. Both C++ files were checked by hand:
- `refine.hpp`: 10 added, 6 removed, 4 net.
- `lattice_mesh.hpp`: 31 added, 14 removed, 17 net.

Test strength: I planted 16 faults, one at a time, in a scratch copy. 14 were caught. Two survived:
- removing the explicit `--src-prefix`/`--dst-prefix`;
- dropping C++ character-literal handling.

`ruff` is clean and the prohibited-dependency gate passes. No red-step scaffolding remains. The workflow change gives the full clone (`fetch-depth: 0`) to the `python` job, which is the only job that runs the whole pytest suite; the pack is 40 MiB.

Citations: `check_citations.py` lists three at-risk citations, and all three are quoted review records of earlier revisions, so they stay.
- `h9-spawn-briefs.md:763` cites `test_brief.py:335` as the file stood at 4ee0328, where that line is a concurrency refusal. That is still true.
- `h11-ci-path-filter.md:451`'s `.github/workflows/main.yaml:306-307` [path written out in full by the main session so the citation gate resolves it] was already off on bc01cd8; the `CI result` name was at line 311 there and is at 314 now.
- `h10-merge-queue.md:183` is unaffected.

Blocking:
1. **[Ola]** There is no CI yet; the branch must be pushed and CI green. The 18 `test_brief.py` failures wait on the refused `tools/brief.py` change.
2. **[now]** `scratch_copy.py`: a child Python process imports the worktree's `tin_engine`, not the copy's. The finder drop lives only in the `-c` program, and a child reloads the editable finder at startup. Probe, with the main checkout's venv in a copy of 4db1eab:
   - parent: `…/scratchpad/copy1/src_python/tin_engine/__init__.py`
   - child: `/Users/skavhaug/projects/rasputin/src_python/tin_engine/__init__.py`

   Three suites spawn such children: `test_cli_mesh_geographic.py:886`, `test_features.py:521` and `test_io_geotiff.py:1271`. In a mutant run they would test the original code and report a false survivor. Fix: carry the drop into children. One route is a generated `sitecustomize.py` on `PYTHONPATH` in the printed command; a sitecustomize runs after the `.pth` file has installed the finder. `test_io_geotiff.py:1271` replaces the whole environment, so where the drop cannot reach, state the limit in the docstring. Add a test in which a child must see the copy.
3. **[now]** `scratch_copy.py:10-11` says "Append test paths to it, or keep the `tests/python/` it ends with". Appending keeps `tests/python/`, so the whole suite runs as well. It should say to replace `tests/python/`.
4. **[now]** `test_recorded_counts_from_real_history` fails in a `git archive` copy ("not a git repository"). That is the same class of failure P9 fixes for `test_settings_wiring.py`, here in a tool built for such copies. Skip it with a reason when the tree is not a git work tree, using P9's `is_work_tree_top`; keep the failure in a shallow clone.
5. **[now]** The status line of `docs/increments/h16-harness-fixes.md` still says "design … next the red step for PR A". Red is 11cee8e and partial green is 4db1eab. Record the measured 252 lines and the reason for the overrun.

Suggestions:
- **[now]** Strip `GIT_CONFIG_COUNT`/`GIT_CONFIG_PARAMETERS` from the environment the counter passes to git, and pass `--inter-hunk-context=0`. Without that, setting `diff.interHunkContext=50` through `GIT_CONFIG_*` in the environment turns PR #163's 41/20 into 63/42; the net stays 21.
- **[now]** Add a test with repo-level `diff.srcPrefix`/`diff.dstPrefix` set, and a C++ case of `'"'` followed by `/*` inside a string. Those are the two planted faults that survived.
- **[now]** Count `.pyi` as Python (`src_python/tin_engine/_core.pyi` is not counted).
- **[now]** `scratch_copy.py` refuses a target only inside the worktree it runs from. A target inside the main checkout, which contains `.claude/worktrees/`, passes; also refuse there.
- **[now]** A failure in `tar` prints more than one stderr line and leaves a partly filled directory.
- **[now]** `count_loc.py` and `scratch_copy.py` are tracked at 100755, while every other `tools/*.py` file is 100644. No rule decides; pick one.
- **[Ola]** A rename from an uncounted path into a counted one counts 0. Probe: a 40-line `docs/proto.py` renamed to `tools/proto.py` counts 0. `CLAUDE.md` §2's literal text (git's default rename detection) gives the same, so a change is a rule change.
- **[Ola]** `CLAUDE.md` §2 still says only "tests excluded", but the counter also leaves out `docs/` and non-code files (ruling 2). The pending §2 pointer line should state both.

Not pushed; no CI.

Taken in the recording commit: blocking item 5 (the status line).

### Round 2: `@reviewer`, code, early, PR A tools part, `4db1eab..731e2a0`

`@reviewer`'s verdict, word for word:

**Code review, round 2 (early, tools part of PR A), 2026-10-05.** Range `4db1eab..731e2a0` (round 1 recorded 1135e1d, red 3de6e1b, green 64e51a2, citation fix 731e2a0). Verdict: CHANGES REQUESTED.

LOC: 280 net production lines by `tools/count_loc.py bc01cd8 731e2a0`: `tools/count_loc.py` 199 and `tools/scratch_copy.py` 81. This round adds 28 net (43 added, 15 removed): count_loc 8/4/4 and scratch_copy 35/11/24. I recounted the scratch_copy part by hand from `git diff -U0` and got the same figures. The estimate was 187. count_loc is 79 over, from the C++ scanner and the git plumbing. scratch_copy is 21 over, from the sitecustomize, the main-checkout refusal and the cleanup after a failed tar; all three came from review and were not priced. With brief.py's roughly 7 lines still to come, the total is about 287 against 700. The design names no split seam, and none is needed.

Packing: there are now 15 `# fmt: skip` regions, 7 in count_loc and 8 in scratch_copy. Four are new this round:
- the `DIFF` tuple;
- the `git rev-parse` call in `_main_checkout`;
- the `tar` call in `_extract`;
- the `quoted` tuple.

Each one keeps one call's arguments, or one tuple, on two or three lines, and each stays readable. Together the regions save about 55 lines; ruff format would give about 229 and 106. They are packed for density only; none is needed for the ceiling.

Round 1's [now] items, each run:
- Blocking 2 (a child process imports the worktree's code). I repeated round 1's probe using the main checkout's venv and a copy of 731e2a0 made by the script. The pytest process and a child `sys.executable -c` both import `…/rv2/copy1/src_python/tin_engine/__init__.py`, and the child's `_core` is the copy's. In the control, the same program run without the printed `PYTHONPATH` imports `/Users/skavhaug/projects/rasputin/src_python/tin_engine/__init__.py` in both processes. Closed. One limit remains: a child started with `-I` still imports the worktree's code. The only such child in the suite is `test_hardening.py:247`, which loads `_core` by path and not as the package, so no test is affected.
- Blocking 3: the docstring now says to replace `tests/python/`. Closed.
- Blocking 4: in the copy, the 5 history tests skip with the reason ("… is not the top of a git work tree"). In a `--depth 1` clone of the branch they fail: 2 recorded-count cases and 3 PR #163 cases. Closed.
- Blocking 5: the status line was taken in 1135e1d. It is out of date again (see below).
- Config in the environment and hunk context: closed. Planting the removal of `--inter-hunk-context=0` is caught.
- Prefixes and character literal (round 1's two surviving faults): I planted both again in a scratch clone. The prefix fault fails 1 test and the character-literal fault fails 1 test. Closed.
- `.pyi`: closed; removing it is caught.
- Main checkout: closed; removing `_main_checkout()` from the refusal is caught.
- tar: closed. Removing the cleanup is caught, and so is a stderr longer than one line.
- Modes: closed. Both files are 100644, like every other file under `tools/`.

Test strength this round: I planted 10 faults, one at a time. 9 were caught. One survived: removing the stripping of config given in the environment (`ENV_CONFIG`). All 65 count_loc tests pass, because `--inter-hunk-context=0` alone overrides the key the tests use. The stripping does matter. Run from `src_python/` with `GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=diff.relative GIT_CONFIG_VALUE_0=true`, the counter gives 280 as written and `total 0 0 0` with the stripping removed.

Gates: 111 tests pass in `test_count_loc.py`, `test_scratch_copy.py` and `test_settings_wiring.py`. 347 tests pass and 15 skip in the copy, in the suites that run with the main venv's older `_core`. `ruff check` and `ruff format --check` are clean, and the prohibited-dependency gate passes. No red-step scaffolding remains.

Citations: `check_citations.py` exits 0 and lists 7 at-risk citations, re-read as quotations:
- The 3 that are new on this range all sit in this file's round-1 record, as the code stood at 4db1eab (`test_brief.py:335`, `.github/workflows/main.yaml:306` [path written out in full by `@architect` so the citation gate resolves it], `scratch_copy.py:10`). They stay.
- `h3-unattended-u1.md:759` names its revision (`as of 6a19357`). It stays.
- The other 3 are as in round 1.

The 731e2a0 edit to the round-1 record is accepted. It is needed (the bare path makes the gate exit 1 with "no such file"), it does not change the meaning, and it is bracketed and attributed.

Blocking:
1. **[Ola]** There is no CI yet; the branch must be pushed and CI green. The 18 `test_brief.py` failures wait on the refused `tools/brief.py` change.
2. **[now]** Design §2 T2 no longer describes the code, and `scratch_copy.py`'s docstring points to it as its spec. Step 1 says the script refuses a target only "inside the repository it is run from", but it now also refuses a target inside the main checkout. Step 4 says the command uses a `-c` program that "drops the editable finder … puts `<dir>/src_python` first on `sys.path`". The command now sets `PYTHONPATH` to a generated `sitecustomize.py` and to `<dir>/src_python`, and its `-c` program only calls `pytest.main`. Rewrite steps 1 and 4 to match, and say that review round 1 caused the change.
3. **[now]** The status line still says "partial green `4db1eab` (252 net …)" and "code review round 1 (early) asked for fixes, in progress". It should state the round-1 green (64e51a2), the 280 lines, and review round 2. The recording commit can take this.

Suggestions:
- **[now]** Add a test that fails when the environment-config stripping is removed. One case is `diff.relative=true` given through `GIT_CONFIG_COUNT`, with the counter run from a subdirectory.
- **[now]** The same key in the repository's own config is not neutralised. `diff.relative=true` in `.git/config`, with the counter run from a subdirectory, should count 0 by the same mechanism as above. I could not plant it: unattended mode refused the `git config` write in my scratch clone. Pass `--no-relative` in `DIFF`, or run git from the top level, and add the test.
- **[later]** The generated `sitecustomize.py` hides Homebrew's own `sitecustomize.py`, so in the copy's processes `sys.base_prefix` is not rewritten and the tk path is not added. No test I ran depends on either. Running the hidden file first, or stating the effect in the docstring, would make it explicit.
- **[later]** "In any Python process a test starts" leaves out children started with `-I` or `-E`, which ignore `PYTHONPATH`. The docstring names only children that replace the environment.
- **[later]** When tar fails on a missing target, `_extract` removes the target but keeps any parent folders that `mkdir(parents=True)` created.

Not pushed; no CI.

Taken in the recording commit: blocking item 3 (the status line).

### Round 3: `@reviewer`, code, early, PR A tools part, `731e2a0..6d35fbd`

`@reviewer`'s verdict, word for word:

**Code review, round 3 (early, tools part of PR A), 2026-10-05.** Range `731e2a0..6d35fbd` (round 2 recorded a61e848, design c3817b3, red 9a1f07a, green 6d35fbd). Verdict: CHANGES REQUESTED.

LOC: 280 net production lines by `tools/count_loc.py bc01cd8 6d35fbd`: `tools/count_loc.py` 199 and `tools/scratch_copy.py` 81. This round adds 0 net (1 added, 1 removed): the `DIFF` tuple's second line gains `"--no-relative"`. I checked this by hand against `git diff -U0`; the other 6 changed lines are comments. The estimate was 187. The overrun is as round 2 recorded it. The design names no split seam, and none is needed. No new `# fmt: skip` region.

Round 2's items, each run:
- Blocking 2 (design §2 T2 out of date): closed. Step 1 now names the main checkout. Step 4 names the generated `sitecustomize.py`, the `PYTHONPATH` and the `-c` program that only calls `pytest.main`, and it says review round 1 caused the change. I ran `scratch_copy.py 6d35fbd` into a scratch folder. The printed command and the `.scratch_copy/sitecustomize.py` it writes match step 4 exactly. A target in the main checkout and a target in the worktree are both refused, as step 1 says. With no built `_core` in this worktree, the step-3 warning is printed. The two citations the rewrite added (`test_io_geotiff.py:1271`, `blockprobe.py:23`) point at the subprocess call and at the line that drops the editable finder.
- Blocking 3 (status line): closed in a61e848.
- Suggestion: config in the repository (`diff.relative`): closed. I removed `--no-relative` by hand to check the test catches it, and the `repository config` case fails.
- Suggestion: a test that fails when the environment-config stripping is removed: **not closed**. I removed the stripping again (`_env` keeps every `GIT_CONFIG_*` key) and all 68 count_loc tests pass. The new `GIT_CONFIG_COUNT` and `GIT_CONFIG_PARAMETERS` cases fail only when `--no-relative` is removed as well: with both removed, 3 tests fail. `--no-relative` now overrides `diff.relative` however it is given, so the stripping is again caught by no test. The stripping still matters. With it removed, `GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=diff.algorithm GIT_CONFIG_VALUE_0=histogram` (or `patience`) turns `count_loc.py 9e666f4 9bb1723` from 744/53/691 into 746/55/691. As written, it stays 744/53/691.

Gates: 114 tests pass in `test_count_loc.py`, `test_scratch_copy.py` and `test_settings_wiring.py`. `ruff check` and `ruff format --check` are clean, and the prohibited-dependency gate passes. No red-step scaffolding remains.

Citations: `check_citations.py` exits 0 and lists 10 at-risk citations. 7 are as in round 2. The 3 new ones are on line 631: round 2's record quoting round 1's citations, as the code stood at 4db1eab. They stay.

Blocking:
1. **[Ola]** There is no CI yet; the branch must be pushed and CI green. The 18 `test_brief.py` failures wait on the refused `tools/brief.py` change.
2. **[now]** The comment above `RELATIVE` in `tests/python/test_count_loc.py` says "the environment forms are what the counter's stripping of `GIT_CONFIG_*` neutralises (no other test fails without it)". Since 6d35fbd that is false: `--no-relative` neutralises them too, and removing the stripping fails no test. Correct the comment. Then close round 2's suggestion with a test that fails when the stripping is removed. The direct way is a unit test that `_env()` drops `GIT_CONFIG_COUNT`, `GIT_CONFIG_KEY_n`, `GIT_CONFIG_VALUE_n` and `GIT_CONFIG_PARAMETERS`. Unlike an end-to-end key, it does not stop catching the removal when the next flag is pinned. A `diff.algorithm` case run through the counter is the end-to-end alternative.
3. **[now]** The status line still says "code review round 2 asked for fixes, in progress" and does not name red 9a1f07a and green 6d35fbd. The recording commit can take this.

Suggestions:
- **[now]** The same `diff.algorithm` key in the repository's own `.git/config` is probably not neutralised, because the stripping covers only the environment. In my probe above, added and removed changed and net did not. Passing `--diff-algorithm=myers` in `DIFF`, with a repository-config test, would pin it; myers is git's default, so the recorded counts stay. I could not test this: unattended mode refused the `git config` write in my scratch clone.
- **[later]** Round 2's three [later] items stand.

Not pushed; no CI.

Taken in the recording commit: blocking item 3 (the status line).

## Ola's rulings, 2026-10-05 afternoon

Ola, verbatim: "yes to all defaults, push both". So: `tools/brief.py` takes the two edits (note names `<persona>-<worktree>-<HHMMSS>.md`, root `*.md` files in `@architect`'s limit); `count_loc.py` pins `--diff-algorithm=myers`; a file renamed across the counted/uncounted boundary counts in full (added on the way in, removed on the way out); `CLAUDE.md` §2's pointer says that `docs/` and non-code files are not counted either. Red for the second and third: `f885887`.

### Round 4: `@reviewer`, code, PR A whole, `bc01cd8..095ab63`

CHANGES REQUESTED. 287 net production lines (`count_loc.py bc01cd8 095ab63`). Round 3's items closed; 219 harness tests pass; rule lines match the tools. Blocking: the docstring of `test_env_drops_config_given_in_the_environment` names `diff.algorithm` as unpinned, which `195b504` made false; no CI before the push. Suggestion: cut "At 98e31cd only the one added line counts" from the rename test's docstring.

### Round 5: `@reviewer`, code, PR A, `095ab63..a6b966e`

CHANGES REQUESTED, on the status line only (taken in the recording commit). 287 net production lines (`count_loc.py bc01cd8 a6b966e`). Round 4's docstring items closed. Green CI after the push makes it APPROVED with no further round.
