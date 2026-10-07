# Harness h18: a new worktree with its own venv, and a master-merge runner

Status: designed on `worktree-harness-tools` off `41bda81`; design review round 2 approved (below under "Review"); red at `e36a5dd6` (72 tests, all failing at load because neither tool exists yet), its pins ruled in §9; next `@developer`, green. One PR, both tools, about 275 net production lines (§5). §8's three questions were ruled by Ola on 2026-10-07: the defaults.

What this is. Two tools from the 2026-10-06 day retrospective
(`docs/retrospectives/2026-10-06-day-bottlenecks-and-merges.md`, table
"Proposals", rows T4 and T2, and the order note under it: "T4 and T2
together (the runner uses the worktree's venv)"). Ola said yes by default
on 2026-10-07 ("1-3 default"). Labels used here:

| Label | What it is |
|---|---|
| T4 | `tools/new_worktree.py`: makes a worktree with a real `.venv` of its own, so Python in that worktree runs that worktree's code |
| T2 | `tools/merge_master.py`: merges `origin/master` into a worktree's branch, stops for conflicts, runs every gate on the merged tree in that worktree's venv, and commits only from a message template |
| master merge | merging `origin/master` into a branch that is behind it (not the merge of a PR into master, which GitHub's merge queue does) |
| `_core` | the compiled C++ extension, `tin_engine._core` |
| the editable finder | `_editable_skbc_rasputin.py` in a venv's `site-packages`, which scikit-build-core's editable install uses to map `tin_engine.*` to one checkout's `src_python/` |
| guarded | `.claude/hooks/guard_push.py` or `guard_governance.py` asks Ola before the command runs |

Where the retrospective row and this file differ, §6 says so.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature.* Sources, read or run 2026-10-07:

- git 2.55, `git-merge(1)`: a merge stopped by a conflict or by
  `--no-commit` can be finished or abandoned with `git merge --abort` or
  `--continue`; "If there were uncommitted worktree changes present when
  the merge started, git merge --abort will in some cases be unable to
  reconstruct these changes. It is therefore recommended to always commit
  or stash your changes before running git merge." T2 refuses a dirty tree
  instead, because the stash is one stack for every worktree (the day
  retrospective, D5). `git-merge(1)` also says `--abort` "is equivalent to
  `git reset --merge` when MERGE_HEAD is present"; `guard_push.py` asks only
  before `reset --hard`, so that form is not guarded either way.
- `git-log(1)`, `--first-parent`: "follow only the first parent commit upon
  seeing a merge commit". On master that lists one merge commit per PR,
  which is what T2 prints and records.
- `git-worktree(1)`: `git worktree add [(-b | -B) <new-branch>] <path>
  [<commit-ish>]`.
- Claude Code's worktree guide (code.claude.com/docs/en/worktrees, as the
  night retrospective of 2026-10-06 quoted it, P6, not re-read here): "A
  worktree is a fresh checkout, so initialize your development environment
  there". The main session creates worktrees with `git worktree add`, so
  nothing initialises them today.
- Anthropic, "Building effective agents", as the day retrospective quoted
  it (section 6, not re-read here): "Poka-yoke your tools. Change the
  arguments so that it is harder to make mistakes." Both tools move steps
  that briefs repeat, and personas skip, into code that runs anyway.
- `uv pip install -e` with scikit-build-core, run for this design on a
  `git archive` copy of `41bda81` in the scratchpad, on AC power (§3.4
  gives the numbers). One finding fixed the design: `uv pip install -e
  "<path>[dev,codecs]"` fails ("Failed to parse: `` … Empty field is not
  allowed for PEP508"); run from the worktree as `-e ".[dev,codecs]"` it
  works.

What differs from the retrospective's rows: §6.

*Legacy.* Nothing; the legacy tree had no harness.

```
$ git grep -l -i -e 'worktree' -e 'merge --no-commit' -e 'uv venv' -e 'virtualenv' legacy-archive -- legacy
(no output; exit 1)
```

## 2. Blueprint

```
main session                       persona (one, end to end)
------------                       -------------------------
new_worktree.py <name>             merge_master.py <wt> [--no-build]
  git fetch origin master            preflight refusals
  git worktree add -b ...            git fetch origin master
  uv venv / uv pip install -e        print + record first-parent merges
  cmake configure build-pyext        git merge --no-ff --no-commit
  import check (paths printed)         conflict? -> exit 3, list files
        |                                 persona: Edit hunks, git add
        v                              merge_master.py <wt> --continue
  <wt>/.venv  <wt>/build-pyext  ---->    fast gates, [_core rebuild], pytest
                                         git commit -F <template>
```

Both tools are plain scripts in `tools/`, stdlib only, with `main(argv,
run)`: every subprocess goes through one function, `run(argv, cwd, *,
capture=True) -> subprocess.CompletedProcess[str]`, which the tests replace
(§3.6, §4.9). The default `run` returns exit 127 with "not found" on stderr
for a missing program, so "uv is not installed" is a refusal, not a
traceback. Neither imports `tin_engine`. A refusal is exit 2 with its
reason as the last line on stderr, starting `new_worktree:` or
`merge_master:` (a passed-through command's output may come before it), as
`tools/scratch_copy.py` does. Both go into `pyproject.toml`'s mypy `files`
list (strict), as `tools/shell_scan.py` is.

The worktree is always named by path, never taken from the current
directory: a subagent's shell returns to the session directory between
calls, which is how the after-commit gate checked the main checkout instead
of the merged worktree (D3 in the day retrospective).

## 3. T4: `tools/new_worktree.py`

### 3.1 Command line

```
python3 tools/new_worktree.py <name> [--from <rev>] [--python <version>]
python3 tools/new_worktree.py --existing <path> [--python <version>]
```

`<name>` gives the path `<main>/.claude/worktrees/<name>` and the branch
`worktree-<name>`, the convention every worktree in `git worktree list`
already follows. `<main>` is the parent of `git rev-parse
--path-format=absolute --git-common-dir`, so the tool run inside a worktree
does not nest the new one there. `--from` defaults to `origin/master`.
`--python` defaults to `3.14`: the main checkout's `.venv` is Python 3.14.7
(its `pyvenv.cfg`, checked 2026-10-07) and 3.14 is CI's newest leg (CI runs
3.12, 3.13 and 3.14).

`--existing <path>` sets up a worktree that already exists (64 of the 107
directories under `.claude/worktrees/` have no `.venv/bin/python`, this
design's own among them; counted 2026-10-07): it does steps 4 to 7
below, skipping each whose product exists, so it is safe to run again.

### 3.2 Steps

1. Refusals (§3.3).
2. If `--from` is `origin/master`: `git fetch origin master`.
3. `git worktree add -b worktree-<name> <main>/.claude/worktrees/<name> <from>`.
4. Unless `<wt>/.venv` exists: `uv venv --python <version> <wt>/.venv`.
5. Unless the step-7 check already passes on this venv (it imports
   `tin_engine` and its submodule from this worktree, not merely from
   somewhere): in `<wt>`,
   `uv pip install --python <wt>/.venv/bin/python -e ".[dev,codecs]" pybind11`.
   This is `INSTALL.md`'s developer install (`.[dev,codecs]`, so ruff,
   mypy, pytest and the codec tests are there), plus `pybind11` so step 6
   can find its CMake files. scikit-build-core builds `_core` with this
   venv's Python, by construction. On a venv whose finder points at another
   checkout, this install is meant to replace that finder with one for this
   worktree, so `--existing` repairs such a venv instead of skipping step 5
   and then refusing at step 7 (a dead end); step 7 checks the result, so
   if the install does not repair it, the refusal stands. `uv pip` writes no `uv.lock` (only `uv
   run` and `uv sync` do; `.gitignore` explains why that matters).
6. Unless `<wt>/build-pyext/CMakeCache.txt` exists: configure, do not build,

   ```
   cmake -S <wt> -B <wt>/build-pyext -DCMAKE_BUILD_TYPE=Release
         -DRASPUTIN_BUILD_PYTHON=ON -DRASPUTIN_BUILD_TESTS=OFF -DRASPUTIN_HARDENING=ON
         -DPYTHON_EXECUTABLE=<py> -DPython_EXECUTABLE=<py>
         -Dpybind11_DIR=<output of: <py> -m pybind11 --cmakedir>
   ```

   with `<py>` = `<wt>/.venv/bin/python`. This is the directory
   `.claude/REQUIRED-READING.md`'s rebuild recipe (`cmake --build
   build-pyext -j --target _core`, then `cp`) assumes. Both spellings of the
   Python variable are set: pybind11's CMake config runs in its classic
   FindPythonInterp mode (the configure run's CMake warning names it), which
   reads `PYTHON_EXECUTABLE` (checked: the main checkout's `build-pyext` has
   `PYTHON_EXECUTABLE=/opt/homebrew/bin/python3` while its
   `Python_EXECUTABLE` names the venv, the "wrong-Python module" lesson of
   the day retrospective). The hardening and build-type values match what
   the editable install builds (`pyproject.toml`'s
   `[tool.scikit-build.cmake.define]`).
7. The check, run with `cwd=<wt>`:
   `<py> -c "import tin_engine, tin_engine.cli, tin_engine._core as c; print(tin_engine.__file__); print(tin_engine.cli.__file__); print(c.__file__)"`,
   then `<py> -m ruff --version` and `<py> -m mypy --version`. It passes
   only if the first two paths lie under `<wt>/src_python/tin_engine/` and
   the third under `<wt>/.venv/`. Every path is compared after
   `Path.resolve()`, `<wt>` and `<main>` included: on macOS a temporary
   directory is `/var/...` as named and `/private/var/...` as Python reports
   a module's file. On a pass it prints one block:

   ```
   worktree  <wt>  (branch worktree-<name> at <short sha>)
   python    <wt>/.venv/bin/python  (Python 3.14.x)
   code      tin_engine and tin_engine.cli from <wt>/src_python
   _core     <path>
   rebuild   cmake --build <wt>/build-pyext -j --target _core   (configured, not built)
   ```

   A submodule is checked as well as the package because a venv whose
   finder points at another checkout routes every submodule there (night
   retrospective 2026-10-06, worktree-setup lessons).

The tool never deletes anything. A failure after step 3 leaves the
worktree, prints what exists, and names the command that finishes it
(`--existing <wt>`).

### 3.3 Refusals (exit 2, the reason the last line on stderr)

| When | Message |
|---|---|
| `<name>` is not `[a-z0-9][a-z0-9._-]*` | `new_worktree: name "<name>": use lower-case letters, digits, '.', '_' and '-' only` |
| the path exists | `new_worktree: <path> already exists; to give it a venv run: python3 tools/new_worktree.py --existing <path>` |
| the branch exists | `new_worktree: branch worktree-<name> already exists; pick another name` |
| `uv --version` or `cmake --version` exits non-zero | `new_worktree: <uv or cmake> is not installed or does not run: <its stderr line>` |
| the fetch fails | `new_worktree: could not fetch origin/master (<git's last line>); blocked on network: stop and hand back` |
| `--existing` path is not the top of a checkout of this repository | `new_worktree: <path> is not the top of a checkout of this repository` |
| step 3, `git worktree add`, fails | `new_worktree: git worktree add failed (exit <n>); its output is above. No worktree was made` |
| a later step's command fails | `new_worktree: <step> failed (exit <n>); its output is above. The worktree is kept; finish with: python3 tools/new_worktree.py --existing <wt>` |
| the check's paths are wrong | `new_worktree: <wt>/.venv imports <module> from <path>, not from <wt>; this venv runs another checkout's code` |

The step commands run with their output passed through (not captured), so
a failed build's own words are on screen.

### 3.4 Measured cost

On a `git archive` copy of `41bda81` in the scratchpad, AC power, uv
0.12.15, warm uv cache, 2026-10-07: `uv venv` 0.04 s; the editable install
with its `_core` build 17.3 s wall; the `build-pyext` configure 0.9 s. A
cold uv cache adds the downloads. The retrospective's "about a minute per
worktree" is an upper bound on this machine.

### 3.5 Guards

Checked by calling the guards' own functions (`guard_push.publishes` and
`segment_why`, `guard_governance.judge_bash`) on each command line, at
`41bda81` and at `worktree-h16b`'s head `c83dfe0e` (guard fix PR B, not
merged): `python3 tools/new_worktree.py x`, `--existing <path>`, and,
typed by hand, `git fetch origin master`, `git worktree add -b
worktree-x <path> origin/master`, `uv venv ...`, `uv pip install ...` and
`cmake -S ... -B ...`. All pass both guards. The tool writes no rule file
and no harness state. Nothing in it needs Ola at the keyboard, attended or
unattended.

### 3.6 Red tests (`tests/python/test_new_worktree.py`, `@tester`)

Marked `harness`; no network, no uv, no C++ build. A temporary repository
with a bare repository as its `origin`; `main([...], run=recorder)`, where
the recorder runs `git` for real and answers `uv`, `cmake` and the venv's
Python from a table (the check's three paths among them). The temporary
directory is passed to the tool unresolved, and the table answers with its
resolved form (`/private/var/...` on macOS), as Python would; the tests
compare paths after `Path.resolve()`.

1. A new worktree lands at `<tmp main>/.claude/worktrees/<name>` on branch
   `worktree-<name>` at `origin/master`'s commit; run from inside another
   worktree, it still lands under the main checkout.
2. The non-git calls, in order: `uv venv --python 3.14 <wt>/.venv`; `uv pip
   install --python <wt>/.venv/bin/python -e .[dev,codecs] pybind11` with
   `cwd=<wt>`; `cmake -S <wt> -B <wt>/build-pyext` carrying both
   `-DPYTHON_EXECUTABLE=<wt>/.venv/bin/python` and
   `-DPython_EXECUTABLE=<wt>/.venv/bin/python`; no `cmake --build`.
3. Each refusal of §3.3, with its message's first words, and nothing
   created for the refusals before step 3. A `git worktree add` that fails
   (`--from` a revision that does not exist) gives the "No worktree was
   made" message, and neither the path nor the branch exists.
4. The check: a table answer naming `<main>/src_python/tin_engine/cli.py`
   for the submodule is refused, naming both paths, and the worktree is
   kept.
5. `--existing` on a worktree that has `.venv` and `build-pyext`, whose
   check passes, makes no `uv venv`, no `uv pip install` and no `cmake`
   call. Whose check first names another checkout's `src_python` (the
   table answers so until the install has run): `uv pip install` is called,
   then the check again, and it passes.
6. A failing `uv pip install` exits 2 with the `--existing` hint and leaves
   the worktree.
7. Every argv the recorder saw passes `guard_push.publishes(argv) == []`
   and `guard_push.segment_why(argv) is None` (the hook loaded by path, as
   `tests/python/harness_fixtures.py` does): the tool runs nothing the push
   guard would ask about, so it is not a route around it.

### 3.7 Size

About 90 net production lines (the retrospective said about 40 with tests;
the difference is `--existing`, the refusals and the path check).

## 4. T2: `tools/merge_master.py`

### 4.1 Command line

```
python3 tools/merge_master.py <worktree> [--no-build]
python3 tools/merge_master.py <worktree> --continue --persona <name> --trailer "<line>" [--no-build]
python3 tools/merge_master.py <worktree> --abort
```

The first form starts the merge and, if git merges it cleanly, goes on as
`--continue` would; then `--persona` and `--trailer` are needed on the
first form too, and it refuses up front without them. `--trailer` is the
Co-Authored-By line from the persona's own system context; the tool cannot
know the model, so the persona passes it. `--persona` gives the subject
tag `(@<name>)`. `--no-build` is for a run whose brief says "You may not
build C++ in this run" (§4.4).

Exit codes: 0 committed, nothing to merge, or aborted; 2 refused; 3 stopped for
conflict resolution; 4 a gate is red, merge not committed.

The suite takes about 6 minutes (§4.8), so a persona runs the first form
and `--continue` with `run_in_background` and waits with Monitor.

### 4.2 State

One file, `<git-dir>/merge_master.json` (the worktree's own git dir,
`git rev-parse --git-dir`, never tracked): the `origin/master` commit being
merged, the merge base, the merge list as printed (§4.3 step 4), the PR
numbers, and the files that conflicted. The commit message is rendered
from it at commit time into `<git-dir>/MERGE_MASTER_MSG`. Both files are
removed after the commit and on `--abort`. The state file plus git's
`MERGE_HEAD` is "a merge this tool started".

### 4.3 Steps: start

1. Refusals (§4.6).
2. `git fetch origin master`.
3. If `git merge-base --is-ancestor origin/master HEAD`: print
   `merge_master: <branch> already contains origin/master (<short>); nothing to merge`
   and exit 0.
4. Print what comes in (recorded in the state file, below): `git merge-base HEAD origin/master`
   as the base, and `git log --merges --first-parent --oneline
   <base>..origin/master`. PR numbers come from subjects of the form
   `Merge pull request #N from ...`. This is the command's output, never a
   typed range (D9 in the day retrospective).

   Then, if `git diff --quiet <base> origin/master -- include src bindings
   lib CMakeLists.txt` says master changed C++ since the merge base: with
   `--no-build`, refuse (§4.6, the row "start with `--no-build`, and master
   changed C++ since the merge base"); with no `<wt>/build-pyext`, refuse
   (§4.6, the row "C++ came in and no `build-pyext`"). Both come before the
   merge, so a run that cannot rebuild `_core` never starts a merge it
   cannot finish.

   Only now, after every start-form check, is the state file (§4.2)
   written; a start that refuses or finds nothing to merge leaves none.
   The conflicted files are added to it after step 5.
5. `git merge --no-ff --no-commit --no-autostash --no-rerere-autoupdate
   origin/master`. `--no-autostash` because the stash is shared by every
   worktree; `--no-rerere-autoupdate` so a recorded resolution, if a user
   config turns rerere on, is never staged unseen.
6. Exit 0 from git: go to §4.4. Unmerged paths (`git diff --name-only
   --diff-filter=U`): record them, print the stop message, exit 3:

   ```
   merge_master: stopped for conflict resolution in <n> files:
     <file>
     <file>   (states rules: your Edit asks Ola; unattended, it is refused)
   Resolve hunk by hunk with Edit. Read each side with
     git show :2:<file>   (this branch)    git show :3:<file>   (origin/master)
   Never take a whole file with checkout --ours/--theirs: it drops the other
   side's clean hunks. Then git add each file (git rm one whose deletion
   you keep), and run
     python3 tools/merge_master.py <wt> --continue --persona <name> --trailer "<line>"
   or abandon the merge with --abort.
   ```

   "States rules" is `governed()` of `.claude/hooks/guard_governance.py`,
   loaded by path; if it cannot be loaded, the mark is left out.

The tool never resolves a conflict, never runs `checkout --ours/--theirs`,
`restore`, `stash`, `reset` or `rebase`, and never stages a file: the
persona's `git add` is the statement that a file is resolved.

### 4.4 Steps: `--continue` (and a clean start)

1. Refusals (§4.6), among them leftover conflict markers.
2. The fast gates, each with `cwd=<wt>`, each through `<wt>/.venv/bin/python`,
   all of them run, output passed through:
   `tools/check_prohibited_deps.py`, `tools/check_detria_boundary.py`,
   `tools/check_citations.py --base origin/master`, `-m ruff check .`,
   `-m ruff format --check .`, `-m mypy`. Any red: exit 4 after the last.
3. `_core`: if `git diff --quiet HEAD -- include src bindings lib
   CMakeLists.txt` says the merged tree's C++ differs from the branch's
   (master brought in C++), the suite needs a rebuilt `_core`. With
   `--no-build`, step 1 has already refused (§4.6), before any gate ran; on
   the first form that refusal came before the merge (§4.3 step 4), on
   `--continue` it leaves the merge in progress. Otherwise `cmake --build
   <wt>/build-pyext -j --target _core`, output passed through; a non-zero
   exit is exit 4 (a gate is red: pytest not run, nothing committed, the
   merge left in progress). Then `build-pyext` must hold exactly one
   `_core*.so` at its top, and its name must end with the venv's extension
   suffix (`<py> -c "import sysconfig;
   print(sysconfig.get_config_var('EXT_SUFFIX'))"`, e.g.
   `.cpython-314-darwin.so`); otherwise refuse (§4.6), since copying a
   `_core` built for another Python is the "wrong-Python module" of the day
   retrospective. That file is copied into
   `<wt>/.venv/lib/python3.*/site-packages/tin_engine/` and its time stamp
   touched (the recipe of `.claude/REQUIRED-READING.md`, "Stale
   artifacts"). No `build-pyext`: refused (§4.6, "C++ came in and no
   `build-pyext`"), on the first form before the merge (§4.3 step 4), on
   `--continue` in step 1, the merge left in progress.
4. The full suite: `<py> -m pytest -q`, with `pyproject.toml`'s options
   (coverage floor included, as CI). Red: exit 4.
5. Print `git diff --cached --stat HEAD` (what the merge changes on this
   branch) and `git diff --cached --stat MERGE_HEAD` (what this branch
   differs from master: the PR's own change set), the "diff against both
   parents" of the day retrospective's merge lessons.
6. Render the message and `git commit -F <git-dir>/MERGE_MASTER_MSG`;
   never `-m`, never `--amend`, never `--no-verify`. The message:

   ```
   Merge origin/master <short> into <branch>: brings in #199, #200 (@developer)

   Merge base <short>. First-parent merges brought in
   (git log --merges --first-parent --oneline <base>..<master>):
     ed125121 Merge pull request #199 from expertanalytics/worktree-landcover-speed
     ...
   Conflicts resolved by hand: a.py, b.md     (or: Conflicts: none)
   Changed by hand beyond git's merge: c.py   (line left out when empty)
   Gates on the merged tree, in <wt>/.venv: check_prohibited_deps,
   check_detria_boundary, check_citations --base origin/master, ruff check,
   ruff format --check, mypy, pytest.
   _core rebuilt: yes.     (or: _core rebuilt: no, no C++ came in.)

   Co-Authored-By: <as given>
   ```

   More than six PRs: the subject names six and "and N more". Without PR
   subjects: "brings in N commits". No body line starts with `#`.

   The hand-changed list is what the persona staged beyond git's own merge:
   `git merge-tree --write-tree HEAD MERGE_HEAD` gives the tree git's merge
   makes on its own (its first output line, printed with exit 1 when there
   are conflicts, conflicted files holding their markers), and `git diff
   --cached --name-only <that tree>` lists every path the index differs in.
   Paths in the conflict list go on the first line, the rest on the second.
   This reads the index and writes only objects, so the tool still stages
   nothing. Checked with git 2.55 on a scratch repository: a content
   conflict, a modify/delete conflict and a hand edit of a cleanly merged
   file all appear in that list, and nothing else.
7. Remove the two state files; print the new head and the subject.

The commit is made inside the tool, so the command line the hooks see is
`python3 tools/merge_master.py ...`, which contains neither `git commit`
nor `git merge`: `gates_after_commit.py` does not fire. The tool's own
gates are a superset of that hook's (it adds mypy and the suite, and runs
in the right tree).

### 4.5 `--abort`

`git merge --abort`, then remove the state files. Safe because the start
refused a dirty tree. "A merge in progress" here is `MERGE_HEAD`, whether
or not this tool started the merge (the start-form refusal for a merge in
progress names `--abort`). With no `MERGE_HEAD`: refuse (§4.6).

### 4.6 Refusals (exit 2, the reason the last line on stderr)

| When | Message |
|---|---|
| `<worktree>` is not the top of a checkout of this repository | `merge_master: <path> is not the top of a checkout of this repository` |
| HEAD is detached, or the branch is `master` | `merge_master: <path> is on <master or a detached HEAD>; run it on a worktree's own branch` |
| no `<wt>/.venv/bin/python` | `merge_master: <wt> has no .venv; run: python3 tools/new_worktree.py --existing <wt>` |
| start, and `git status --porcelain` is not empty | `merge_master: <wt> has uncommitted changes (<first file>); commit them first. Do not use git stash: every worktree shares one stash` |
| start, and a merge is already in progress | `merge_master: a merge is already in progress in <wt>; finish it with --continue or abandon it with --abort` |
| the first form or `--continue` without `--persona` or `--trailer` | `merge_master: <the missing option> is missing; starting or continuing a merge needs --persona and --trailer` |
| `--persona` not a persona of `tools/brief.py`, or a read-only one | `merge_master: --persona <name>: give one of architect, developer, orchestrator, perf, tester` (the writers, read from `brief.WRITES`, in its order) |
| `--trailer` does not match `^Co-Authored-By: [^<>]+ <[^<>@]+@[^<>]+>$` | `merge_master: --trailer must be the Co-Authored-By line from your system context` |
| the fetch fails | `merge_master: could not fetch origin/master (<git's last line>); blocked on network: stop and hand back` |
| `--continue`, and no state file | `merge_master: no merge started by this tool in <wt>; start one with: python3 tools/merge_master.py <wt>` |
| `--continue`, state file but no `MERGE_HEAD` | `merge_master: the merge in <wt> was left (MERGE_HEAD is gone: a checkout, reset or commit during it). Nothing was committed by this tool; stop and hand back` |
| `MERGE_HEAD` is not the recorded master commit | `merge_master: MERGE_HEAD is <a>, but this tool merged <b>; stop and hand back` |
| unmerged paths remain | `merge_master: still unresolved: <files>; resolve with Edit, then git add` |
| a recorded conflicted file has a line starting `<<<<<<< ` or `>>>>>>> ` | `merge_master: <file>:<line> still holds a conflict marker` |
| tracked changes not staged | `merge_master: <file> is changed but not added; git add it if it belongs to the merge` |
| `--abort`, and no `MERGE_HEAD` | `merge_master: no merge in progress in <wt>; nothing to abort` |
| start with `--no-build`, and master changed C++ since the merge base (§4.3 step 4) | `merge_master: origin/master changed C++ (<first file>); the suite needs a rebuilt _core and this run may not build C++. No merge was started; hand back` |
| `--continue` with `--no-build`, and C++ came in (§4.4 step 3) | `merge_master: origin/master changed C++ (<first file>); the suite needs a rebuilt _core and this run may not build C++. The merge is left in progress; hand back` |
| C++ came in and no `build-pyext` (start: §4.3 step 4; `--continue`: §4.4 step 1) | `merge_master: origin/master changed C++ (<first file>) and <wt> has no build-pyext to rebuild _core in; run: python3 tools/new_worktree.py --existing <wt>. ` then `No merge was started` (start) or `The merge is left in progress` (`--continue`) |
| after the `_core` build, `build-pyext` holds no `_core*.so`, more than one, or one without the venv's extension suffix | `merge_master: <wt>/build-pyext holds <the names, or no _core*.so>, not one _core<suffix>; remove build-pyext and run: python3 tools/new_worktree.py --existing <wt>. The merge is left in progress` |

A conflicted file that the resolution deleted (`git rm`, as a
modify/delete conflict may end) has no markers to search and is skipped.
Conflict markers are searched only in the files that conflicted, and only
the start and end markers: a line of `=======` is a Markdown heading
underline as often as it is a marker, and `git diff --check` would also
refuse whitespace errors that master's own text may carry.

### 4.7 Guards

The same check as §3.5, on `python3 tools/merge_master.py <wt>`, its
`--continue` and `--abort` forms, and, typed by hand, `git -C <wt> merge
--no-ff --no-commit origin/master`, `git -C <wt> commit -F <file>`, `git -C
<wt> merge --abort` and `git -C <wt> add CLAUDE.md`: all pass both guards,
at `41bda81` and at `c83dfe0e`. `git checkout --theirs CLAUDE.md` is asked
by `guard_governance` ("it writes CLAUDE.md"), which the tool never runs.

Three things follow.

- **The tool is not a route around the push guard, by test.** The hooks
  see only the line typed, not what a script runs, so a later edit adding
  `git push` to the tool would not be asked about. Red test 9 below fails
  on any argv the push guard would ask about. The two tools stay off
  `guard_governance.py`'s list (§8 question 3, ruled).
- **A merge writes rule files without asking, as it does today.** A clean
  merge that brings in master's changes to `CLAUDE.md` or a hook writes
  them through git, unseen by the guards, exactly as a hand-typed `git
  merge` does now (checked: `git merge` names no targets for
  `guard_governance`). That text is already on master, merged with Ola's
  yes. Not a new hole.
- **Where Ola is needed.** Only when a conflict falls in a rule file: the
  persona's Edit is asked about (attended) or refused and queued
  (unattended). The stop message marks those files, so the persona can
  `--abort` or leave the merge in progress and hand back an `ASK OLA:`
  line. The tool itself never needs Ola at the keyboard.

### 4.8 Measured cost

On the same scratch copy as §3.4: `mypy` 4.2 s; the full suite 6 min 3 s
(5683 passed, 31 skipped, of which 8 skip only because a `git archive` copy
is not a work tree); a `_core` rebuild in a configured `build-pyext` 16.3 s.

### 4.9 Red tests (`tests/python/test_merge_master.py`, `@tester`)

Marked `harness`; no network, no C++ build, no real suite. A temporary
repository with a bare `origin`, whose `master` gets first-parent merge
commits with subjects `Merge pull request #12 from x/y`; a worktree branch
off an older master commit; a stub `<wt>/.venv/bin/python` path for the
venv check. The recorder runs `git` for real and answers the gates and
`cmake` from a table (red or green per test), the venv's extension suffix
among them; for the `_core` tests the test puts a dummy
`build-pyext/_core<suffix>.so` in place, and the tool finds it by the glob
`_core*.so` and checks its suffix (the suffix differs between macOS and
Linux). The temporary repository is passed to the tool as the test made
it, unresolved; every path the tests compare is compared after
`Path.resolve()` (on macOS `/var` is `/private/var`).

1. Each refusal of §4.6 with its message's first words; no merge started
   by any start-form refusal, the `--no-build` one included (it comes
   before the merge, §4.3 step 4). The refusals that come after the merge
   (`--continue` with `--no-build`, the `_core*.so` check) leave it in
   progress: `MERGE_HEAD` and the state file still there.
2. Up to date: exit 0, no new commit, no state file.
3. Clean merge: the gates are called in §4.4's order with
   `<wt>/.venv/bin/python` and `cwd=<wt>`, pytest last; one commit with two
   parents, the old head and `origin/master`; its subject ends `(@developer)`
   and contains `#12`; its last line is the trailer as given; the recorded
   merge list equals `git log --merges --first-parent --oneline
   <base>..origin/master` run by the test; the state files are gone.
4. Conflict: exit 3; no commit; `MERGE_HEAD` present; the stop message
   lists the file; a conflicted `CLAUDE.md` carries the "states rules"
   mark; the file still holds its markers (the tool resolved nothing).
5. `--continue` with a marker left: refused, naming `file:line`. After the
   test resolves and `git add`s, and also edits and `git add`s a file that
   merged cleanly: committed, body "Conflicts resolved by hand: <file>"
   and "Changed by hand beyond git's merge: <the other file>"; with no such
   edit, no second line.

   Modify/delete: master deletes a file this branch changed. Exit 3, the
   stop message lists it; after the test `git rm`s it, `--continue`
   commits without tripping on the missing file, and the body lists it
   under "Conflicts resolved by hand".
6. A red fast gate: exit 4, every fast gate still called, pytest not
   called, no commit, `MERGE_HEAD` present, the gate's output reaches the
   tool's output unchanged. A red pytest: exit 4, no commit.
7. C++ on master (a file under `include/` changed): `cmake --build
   <wt>/build-pyext -j --target _core` is called after the fast gates and
   before pytest. With `--no-build` on the first form: refused, no build
   call, no merge started. With `--no-build` on `--continue` after a
   conflict: refused, no build call, merge left in progress. A build that
   exits non-zero: exit 4, pytest not called, no commit, merge left in
   progress. Two `_core*.so` files, or one with another suffix: refused,
   nothing copied.
8. `--continue` after the test deletes `MERGE_HEAD` (as a stray checkout
   would): refused with "the merge in <wt> was left"; nothing committed.
9. Every argv the recorder saw passes `guard_push.publishes(argv) == []`
   and `guard_push.segment_why(argv) is None`; none is `stash`, `reset`,
   `rebase`, `checkout`, `restore` or `push`; exactly one is `git commit`,
   and it carries `-F` and no `-m`.
10. `--abort`: `MERGE_HEAD` gone, the tree equal to the old head, the state
    files gone.

Not the invariant-critical suite of anything, so no mutation round
(`docs/increments/README.md`, "Cost constraints").

### 4.10 Size

About 185 net production lines (the retrospective, after the 2026-10-03
research, said about 120 with tests; the difference is the two-phase
command line, the `_core` rebuild with its suffix check, the hand-change
list and the refusal table).

## 5. PR and who does what

One PR, both tools, on a worktree branch made with `git worktree add` (T4
does not exist yet): about 275 net production lines by `python3
tools/count_loc.py <base> <head>`, under `CLAUDE.md` §2's limit. Not
refine or mesh code, so no `@perf` acceptance run.

- `@tester`: §3.6 and §4.9, red, committed before any tool exists.
- `@developer`: both tools, the two `pyproject.toml` mypy entries; touches
  no test file.
- `@reviewer`: code review, including the §3.5 and §4.7 guard check rerun
  on the branch, and one real run of `new_worktree.py` on a throwaway name
  (it needs the network and builds `_core`, so not in the suite), removed
  afterwards with `git worktree remove` and `git branch -D`.

Nothing in this PR writes a rule file: `tools/new_worktree.py` and
`tools/merge_master.py` are not on `guard_governance.py`'s list, and
neither are `tests/` and `pyproject.toml`. Ola is needed for the push and
the enqueue, as for every PR.

## 6. Where this differs from the retrospective's rows

- **T4's `--base <sha>`, a scratch venv for a base install, is left out**
  (§8 question 2, ruled: left out). Its two uses are covered: `tools/bench.py run --tree`
  builds an older tree into its own `build-bench/pkg`, and
  `tools/scratch_copy.py` runs the suites against an older revision in the
  worktree's interpreter. The lesson behind it, "base installs go in
  scratch venvs", is kept by T4 never installing anything but the
  worktree's own checkout.
- **T2's lineage.** The row calls T2 "stage 3, already ruled yes". Stage 3
  of `docs/research/2026-10-03-dispatcher-control.md` is a different tool,
  `tools/merge_queue.py`, for merging PRs into master (a lock, "no checks
  reported is not green"), whose job GitHub's merge queue now does (h10).
  T2's authority is Ola's yes of 2026-10-07, not stage 3.
- **"One persona runs it end to end"** meets the write limits: a conflict in
  `tests/` or `docs/` is outside `@developer`'s limit (`tools/brief.py`'s
  `WRITES`). Ruled by Ola: §8 question 1.
- **The gates** add `check_prohibited_deps.py` and
  `check_detria_boundary.py` to the row's list (both fast; both are
  `gates_after_commit.py` gates, and that hook does not fire on the tool's
  commit, §4.4), and
  rebuild `_core` when C++ came in, without which the suite tests new
  Python against an old extension.
- **Who builds `_core`.** The row says T4 builds `_core` "with
  `-DPYTHON_EXECUTABLE` set to that venv's Python". Here the editable
  install (§3.2 step 5) builds `_core`, with the venv's Python by
  construction, and `build-pyext` is only configured, with
  `-DPYTHON_EXECUTABLE` (step 6), for later rebuilds by the
  `REQUIRED-READING.md` recipe and by T2. The step-7 import check, which
  requires `_core` to load from `<wt>/.venv/`, is what shows the venv's
  `_core` is the one in use.

## 7. Follow-ups: rule and brief text these tools let us change (not in this PR)

The saving is in briefs (item 3), not in rule files: items 1 and 2 keep
their length and make the text name a route that sets up the venv; item 4
cuts nothing, since no rule line states that lesson today; items 5 and 6
add. Items 1, 2, 4, 5 and 6 are writes to governed files, so each asks Ola
(`guard_governance.py`).

1. `.claude/REQUIRED-READING.md`, "One session per working tree":
   "(`git worktree add`)" becomes "(`python3 tools/new_worktree.py
   <name>`)". Same length, but the one documented route then makes the
   venv.
2. `.claude/briefs/common.md`: "with your own build directory and venv"
   becomes "with its own `.venv` and `build-pyext` (made by
   `tools/new_worktree.py`)". About the same length; the brief then names
   a thing that exists instead of asking each persona to make one.
3. The main session's merge briefs (not rule text): the steps the 30c and
   audit D merge lessons asked briefs to spell out (`--no-commit`, `-F`
   with the trailer, no `--theirs`, ruff and the full suite, the
   `git log --merges --first-parent` range; the day retrospective,
   section 9) become one line, "run `python3
   tools/merge_master.py <wt>`". This is the largest saving, about 20 writer
   spawns a day of split merges on 2026-10-06 (the day retrospective, D8).
4. `.claude/REQUIRED-READING.md`, "Stale artifacts": the rebuild recipe
   stays (it is how a persona rebuilds outside a merge), and the "first
   build-pyext needs `-DPYTHON_EXECUTABLE`" lesson, which has no rule line
   today (`grep -n PYTHON_EXECUTABLE .claude/REQUIRED-READING.md` finds
   nothing), keeps needing none, because `build-pyext` now comes
   configured.
5. Question 1 is ruled (§8): one line in `CLAUDE.md` §3 or the h6 role table
   saying who runs a master merge. That adds a line; it replaces the
   per-merge split that D8 counts.
6. Optional, a check rather than a cut: `tools/brief.py` refuses a brief
   for a worktree with no `.venv/bin/python`, naming `new_worktree.py
   --existing`. Governed (`tools/brief.py` is on the list).

No persona file mentions venvs or master merges today (`grep -n -i -e
venv -e "merge master" -e "--theirs" .claude/agents/*.md` exits 1 with no
output), so none is shortened.

## 8. Questions for Ola: ruled

Ola, 2026-10-07: "defaults on the h18 questions". So:

1. **Who runs a master merge.** `@developer` runs `merge_master.py` end
   to end and may resolve a conflict in any file, inside its write limit or
   not: choosing between two already-reviewed versions hunk by hunk is not
   authoring. A test that fails after the merge but did not conflict (as in
   audit D's merge, where a cleanly merged test broke on master's changed
   helper) goes to `@tester` in the same worktree while the merge is still
   in progress; then `@developer` runs `--continue`.
2. **T4's `--base` scratch venv is left out.** `bench.py --tree` and
   `scratch_copy.py` cover its two uses (§6).
3. **The two tools are not added to `guard_governance.py`'s list.** Red
   test 9 of §4.9 and test 7 of §3.6 fail on any command the push guard
   would ask about, and review sees the failure.

## 9. The red step's pins, ruled

`@tester`'s red commit `e36a5dd6` (`tests/python/test_new_worktree.py`,
`tests/python/test_merge_master.py`, fixture
`tests/python/worktree_fixtures.py`) stated 19 assumptions beyond §3.6 and
§4.9. Each is ruled here; the sections above carry the ones that set text.
One changes a test (18); the rest keep the suite as committed.

1. **Kept.** Each tool is loaded from its copy in the temporary main
   checkout. "This repository" is the one holding the tool's own file
   (`Path(__file__).resolve()`, then its `git rev-parse --git-common-dir`),
   not the current directory (§2).
2. **Kept.** `main(argv: list[str], run: Run = run) -> int`, with `run` a
   keyword the tests pass; the exit code is returned (argparse's own usage
   errors may raise `SystemExit(2)`). A module-level `run(argv, cwd, *,
   capture=True)` is the default.
3. **Kept.** The tools run only `git`, `uv`, `cmake` and
   `<wt>/.venv/bin/python*`. Copying `_core`, touching it, globbing and
   reading files are done in Python (`shutil`, `Path`), never by `cp` or
   `touch` subprocesses.
4. **Kept.** The probe commands are exactly: `uv --version`, `cmake
   --version`, `<py> -m ruff --version`, `<py> -m mypy --version`, `<py>
   -m pybind11 --cmakedir`, the §3.2 step 7 check (one `-c` naming
   `tin_engine.cli` and `tin_engine._core`, three lines out), and, for
   `merge_master.py`, `<py> -c` printing `sysconfig`'s `EXT_SUFFIX`. The
   pass block's Python version comes from `<py> --version`.
5. **Kept.** The refusal is the last line on stderr (§2 and the §3.3 and
   §4.6 headings now say so); a refusal before any passed-through command
   (bad name, existing path, existing branch) is the only line.
6. **Kept.** `Alpha`, `a/b`, `_a`, `.hidden`, `a b` and the empty name are
   refused with §3.3's line exactly. A name the pattern allows but git
   refuses as a branch (`a..b`) fails at `git worktree add`: "No worktree
   was made".
7. **Kept.** The fetch is `git fetch origin master` (the argv ends `origin
   master`), on the first form only; `--existing` fetches nothing. The
   `uv` and `cmake` probes are step 1 on both forms, so `--existing` with a
   good venv and a configured `build-pyext` runs the probes and the check
   and nothing else. A failing install stops before the configure;
   `<step>` in its message is the command's name (`uv pip install`).
   `--existing` on a worktree with neither runs venv, install and
   configure.
8. **Kept; messages set.** A missing `--persona` or `--trailer` on the
   first form or on `--continue` is exit 2 with the new §4.6 row naming the
   missing option; `--abort` with no `MERGE_HEAD` is exit 2 with the new
   row "`--abort`, and no `MERGE_HEAD`" (§4.5).
9. **Kept.** The persona list is compared as a set; the message lists
   `brief.WRITES`'s writers in its order (§4.6 now says so).
10. **Kept.** The citations gate's argv ends `--base origin/master`; the
    build is `cmake --build <wt>/build-pyext -j --target _core`, nothing
    after `_core`. The gate scripts run as `<py> tools/<script>`.
11. **Kept.** The marker refusal names the first marker line of the first
    recorded conflicted file that still holds one, files in the order git
    listed them.
12. **Kept.** Both hashes in the `MERGE_HEAD` message are `git rev-parse
    --short` prefixes (7 or more characters).
13. **Kept; template pinned.** The body lines are literally `Conflicts:
    none`, `Conflicts resolved by hand: <a>, <b>` (comma and space),
    `Changed by hand beyond git's merge: <list>`, and the `_core` sentence
    on a line of its own, `_core rebuilt: yes.` or `_core rebuilt: no, no
    C++ came in.` (§4.4 step 6 now shows it so). Merge list lines are
    indented in the body and printed to stdout; the tests compare them
    trimmed. "brings in N commits" counts `git rev-list --first-parent
    --count <base>..<master>`.
14. **Kept.** The conflict stop message goes to stdout, all of it (it is
    instructions, not an error); one file per line, indented two spaces.
15. **Kept.** "Never stages" is tested as no `git add` or `git rm`, and no
    `stash`, `reset`, `rebase`, `checkout`, `restore` or `push`. The tool
    also runs no `update-index` or `read-tree`; `git merge --abort` and
    `git merge-tree --write-tree` are allowed.
16. **Kept.** A successful `--abort` exits 0 (§4.1's exit codes now say
    so).
17. **Kept.** The copy goes into the one directory matching
    `<wt>/.venv/lib/python3.*/site-packages/tin_engine/` (the stub venv's
    is `python3.12`); "touched" means its mtime is set to now after the
    copy (`Path.touch()`).
18. **Changed.** No `build-pyext` when C++ came in is a refusal with a set
    message (§4.6, "C++ came in and no `build-pyext`"), and on the first
    form it comes before the merge, beside the `--no-build` check (§4.3
    step 4), for the same reason. Test change, for `@tester`: in
    `test_cpp_from_master_without_build_pyext_is_refused_naming_new_worktree`,
    replace the head check with `assert_no_merge_started(setup, outcome)`
    and add `assert builds(outcome) == []`. The suite stays red either way.
    On `--continue` with `--no-build` and no `build-pyext`, the `--no-build`
    row wins (it comes first in §4.6).
19. **Kept.** "already contains origin/master … nothing to merge" goes to
    stdout. The gates, the build and pytest run with `capture=False`, so
    their output passes through unchanged.

## Review

Design review round 1, @reviewer: CHANGES REQUESTED, 2 blockers: (1) the h18 row at /Users/skavhaug/projects/rasputin/ROADMAP.md@360259c0:45 shifts ROADMAP.md:50 and :54, which 45 citations quote as rows 27 and 29 (/Users/skavhaug/projects/rasputin/ROADMAP.md@41bda81a:50, @41bda81a:54); move the row below line 54; (2) the §4.9 test 1 pin "no merge started by a start-form refusal" (/Users/skavhaug/projects/rasputin/docs/increments/h18-worktree-and-merge-tools.md@360259c0:448-449) contradicts the --no-build refusal, which on a clean start comes after the merge (@360259c0:330-331, :395, :468-469).

Design review round 2, @reviewer: APPROVED, range 360259c0..9e11868e (docs only, 0 production lines), both round 1 blockers answered: (1) the h18 row now sits at /Users/skavhaug/projects/rasputin/ROADMAP.md@9e11868e:70, and against 41bda81a the only change is that inserted line, so :50 and :54 are byte-identical to @41bda81a:50 and :54 (rows 27 and 29); (2) the --no-build check now runs before the merge on the first form (/Users/skavhaug/projects/rasputin/docs/increments/h18-worktree-and-merge-tools.md@9e11868e:313-316, :440), refuses after it on --continue (@9e11868e:441), and tests 1 and 7 pin both cases (@9e11868e:500-504, :528-535).
