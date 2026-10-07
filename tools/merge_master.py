#!/usr/bin/env python3
"""Merge origin/master into a worktree's branch, gate the result, commit from a template (h18 T2).

    python3 tools/merge_master.py <wt> --persona <name> --trailer "<line>" [--no-build]
    python3 tools/merge_master.py <wt> --continue --persona <name> --trailer "<line>" [--no-build]
    python3 tools/merge_master.py <wt> --abort

The first form fetches origin/master, prints the first-parent merges it brings
in and starts `git merge --no-ff --no-commit`. A conflict stops it (exit 3) for
the persona to resolve hunk by hunk and `git add`; `--continue` then runs every
gate and the suite in `<wt>/.venv`, rebuilds `_core` in `<wt>/build-pyext` if
C++ came in, and commits with `git commit -F`. The tool never resolves,
checks out, restores, stashes, resets, rebases or stages anything.

Exit 0 committed, nothing to merge, or aborted; 2 refused (the reason is the
last line on stderr); 3 stopped for conflict resolution; 4 a gate is red.
Spec: docs/increments/h18-worktree-and-merge-tools.md §4. Standard library only.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
from new_worktree import HERE, Refusal, Run, check_checkout, last_line, run

CXX = ("include", "src", "bindings", "lib", "CMakeLists.txt")
TRAILER = re.compile(r"^Co-Authored-By: [^<>]+ <[^<>@]+@[^<>]+>$")
PR = re.compile(r"Merge pull request #(\d+) from ")
GATES = (["tools/check_prohibited_deps.py"], ["tools/check_detria_boundary.py"],
         ["tools/check_citations.py", "--base", "origin/master"], ["-m", "ruff", "check", "."],
         ["-m", "ruff", "format", "--check", "."], ["-m", "mypy"])  # fmt: skip
GATE_NAMES = ("check_prohibited_deps, check_detria_boundary, check_citations --base origin/master,"
              " ruff check, ruff format --check, mypy, pytest")  # fmt: skip


def load(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(f"merge_master_{path.stem}", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # a dataclass in it looks its module up here
    spec.loader.exec_module(module)
    return module


class Tree:
    """One worktree, every git call made with `-C <wt>` through `run`."""

    def __init__(self, wt: Path, run: Run) -> None:
        self.wt, self.run, self.py = wt, run, str(wt / ".venv" / "bin" / "python")
        self.git_dir = Path(self.git("rev-parse", "--absolute-git-dir").stdout.strip())
        self.state, self.msg = self.git_dir / "merge_master.json", self.git_dir / "MERGE_MASTER_MSG"
        self.merge_head = self.git_dir / "MERGE_HEAD"

    def git(self, *args: str, capture: bool = True) -> subprocess.CompletedProcess[str]:
        return self.run(["git", "-C", str(self.wt), *args], capture=capture)

    def out(self, *args: str) -> str:
        return (self.git(*args).stdout or "").strip()

    def lines(self, *args: str) -> list[str]:
        return [line for line in (self.git(*args).stdout or "").splitlines() if line.strip()]

    def cxx(self, *against: str) -> str | None:
        """The first C++ file changed between `against`; None if none."""
        changed = self.lines("diff", "--name-only", *against, "--", *CXX)
        return changed[0] if changed else None

    def drop_state(self) -> None:
        self.state.unlink(missing_ok=True)
        self.msg.unlink(missing_ok=True)


def cxx_refusal(t: Tree, first: str, no_build: bool, ending: str) -> None:
    if no_build:
        raise Refusal(f"origin/master changed C++ ({first}); the suite needs a rebuilt _core and "
                      f"this run may not build C++. {ending}; hand back")  # fmt: skip
    if not (t.wt / "build-pyext").is_dir():
        raise Refusal(f"origin/master changed C++ ({first}) and {t.wt} has no build-pyext to "
                      f"rebuild _core in; run: python3 tools/new_worktree.py --existing {t.wt}. "
                      f"{ending}")  # fmt: skip


def start(t: Tree, args: argparse.Namespace) -> int:
    if t.merge_head.exists():
        raise Refusal(f"a merge is already in progress in {t.wt}; finish it with --continue or "
                      "abandon it with --abort")  # fmt: skip
    dirty = t.lines("status", "--porcelain")
    if dirty:
        raise Refusal(f"{t.wt} has uncommitted changes ({dirty[0][3:]}); commit them first. "
                      "Do not use git stash: every worktree shares one stash")  # fmt: skip
    fetch = t.git("fetch", "origin", "master")
    if fetch.returncode != 0:
        raise Refusal(f"could not fetch origin/master ({last_line(fetch)}); "
                      "blocked on network: stop and hand back")  # fmt: skip
    master = t.out("rev-parse", "origin/master")
    short = t.out("rev-parse", "--short", master)
    branch = t.out("rev-parse", "--abbrev-ref", "HEAD")
    if t.git("merge-base", "--is-ancestor", master, "HEAD").returncode == 0:
        print(f"merge_master: {branch} already contains origin/master ({short}); nothing to merge")
        return 0
    base = t.out("merge-base", "HEAD", master)
    merges = t.lines("log", "--merges", "--first-parent", "--oneline", f"{base}..{master}")
    print(f"merge_master: merge base {t.out('rev-parse', '--short', base)}; first-parent merges "
          "coming in:")  # fmt: skip
    print("".join(f"  {line}\n" for line in merges), end="")
    first = t.cxx(base, master)
    if first is not None:
        cxx_refusal(t, first, args.no_build, "No merge was started")
    count = t.out("rev-list", "--first-parent", "--count", f"{base}..{master}")
    state: dict[str, Any] = {"master": master, "base": base, "merges": merges, "conflicts": [],
                             "prs": [m[1] for m in map(PR.search, reversed(merges)) if m],
                             "count": count}  # fmt: skip
    t.state.write_text(json.dumps(state, indent=1))
    merge = t.git("merge", "--no-ff", "--no-commit", "--no-autostash", "--no-rerere-autoupdate",
                  "origin/master", capture=False)  # fmt: skip
    if merge.returncode == 0:
        return finish(t, args)
    state["conflicts"] = t.lines("diff", "--name-only", "--diff-filter=U")
    if not state["conflicts"] or not t.merge_head.exists():
        t.drop_state()
        raise Refusal(f"git merge failed (exit {merge.returncode}); its output is above")
    t.state.write_text(json.dumps(state, indent=1))
    stop(t, state["conflicts"])
    return 3


def stop(t: Tree, conflicts: list[str]) -> None:
    try:
        governed = load(HERE.parent / ".claude" / "hooks" / "guard_governance.py").governed
    except Exception:  # the mark is optional: a guard that does not load leaves it out
        governed = None
    mark = "   (states rules: your Edit asks Ola; unattended, it is refused)"
    listed = "".join(f"  {f}{mark if governed and governed(f) else ''}\n" for f in conflicts)
    print(f"merge_master: stopped for conflict resolution in {len(conflicts)} files:\n{listed}"
          "Resolve hunk by hunk with Edit. Read each side with\n"
          "  git show :2:<file>   (this branch)    git show :3:<file>   (origin/master)\n"
          "Never take a whole file with checkout --ours/--theirs: it drops the other\n"
          "side's clean hunks. Then git add each file (git rm one whose deletion\n"
          "you keep), and run\n"
          f"  python3 tools/merge_master.py {t.wt} --continue --persona <name> --trailer "
          "\"<line>\"\nor abandon the merge with --abort.")  # fmt: skip


def finish(t: Tree, args: argparse.Namespace) -> int:
    """§4.4: the refusals, the gates, `_core`, the suite, and the commit."""
    if not t.state.exists():
        raise Refusal(f"no merge started by this tool in {t.wt}; start one with: "
                      f"python3 tools/merge_master.py {t.wt}")  # fmt: skip
    if not t.merge_head.exists():
        raise Refusal(f"the merge in {t.wt} was left (MERGE_HEAD is gone: a checkout, reset or "
                      "commit during it). Nothing was committed by this tool; stop and "
                      "hand back")  # fmt: skip
    state = json.loads(t.state.read_text())
    found, merged = t.merge_head.read_text().strip(), state["master"]
    if found != merged:
        a, b = t.out("rev-parse", "--short", found), t.out("rev-parse", "--short", merged)
        raise Refusal(f"MERGE_HEAD is {a}, but this tool merged {b}; stop and hand back")
    unresolved = t.lines("diff", "--name-only", "--diff-filter=U")
    if unresolved:
        raise Refusal(f"still unresolved: {', '.join(unresolved)}; resolve with Edit, then git add")
    for name in state["conflicts"]:
        path = t.wt / name
        lines = path.read_text(errors="replace").splitlines() if path.is_file() else []
        for number, line in enumerate(lines, 1):
            if line.startswith(("<<<<<<< ", ">>>>>>> ")):
                raise Refusal(f"{name}:{number} still holds a conflict marker")
    unstaged = t.lines("diff", "--name-only")
    if unstaged:
        raise Refusal(f"{unstaged[0]} is changed but not added; git add it if it belongs to the "
                      "merge")  # fmt: skip
    first = t.cxx("--cached", "HEAD")
    if first is not None:
        cxx_refusal(t, first, args.no_build, "The merge is left in progress")
    red = [" ".join(g) for g in GATES if t.run([t.py, *g], t.wt, capture=False).returncode != 0]
    if red:
        print(f"merge_master: red: {'; '.join(red)}. Nothing committed; the merge is left in "
              "progress", file=sys.stderr)  # fmt: skip
        return 4
    if first is not None and not rebuild(t):
        return 4
    if t.run([t.py, "-m", "pytest", "-q"], t.wt, capture=False).returncode != 0:
        print("merge_master: pytest is red. Nothing committed; the merge is left in progress",
              file=sys.stderr)  # fmt: skip
        return 4
    t.git("diff", "--cached", "--stat", "HEAD", capture=False)
    t.git("diff", "--cached", "--stat", "MERGE_HEAD", capture=False)
    t.msg.write_text(message(t, args, state, rebuilt=first is not None))
    done = t.git("commit", "-F", str(t.msg), capture=False)
    if done.returncode != 0:
        raise Refusal(f"git commit failed (exit {done.returncode}); its output is above. The "
                      "merge is left in progress")  # fmt: skip
    t.drop_state()
    print(f"merge_master: committed {t.out('rev-parse', '--short', 'HEAD')} "
          f"{t.out('log', '-1', '--format=%s')}")  # fmt: skip
    return 0


def rebuild(t: Tree) -> bool:
    """§4.4 step 3: build `_core`, check its suffix, copy it into the venv and touch it."""
    build = t.wt / "build-pyext"
    if t.run(["cmake", "--build", str(build), "-j", "--target", "_core"],
             capture=False).returncode != 0:  # fmt: skip
        print("merge_master: the _core build is red. Nothing committed; the merge is left in "
              "progress", file=sys.stderr)  # fmt: skip
        return False
    suffix = t.run([t.py, "-c", "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"],
                   t.wt).stdout.strip()  # fmt: skip
    cores = sorted(p.name for p in build.glob("_core*.so"))
    packages = sorted(t.wt.glob(".venv/lib/python3.*/site-packages/tin_engine"))
    if cores != [f"_core{suffix}"] or len(packages) != 1:
        held = ", ".join(cores) or "no _core*.so"
        if len(packages) != 1:
            held += f" (and {t.wt}/.venv has {len(packages)} tin_engine package dirs)"
        raise Refusal(f"{build} holds {held}, not one _core{suffix}; remove build-pyext and run: "
                      f"python3 tools/new_worktree.py --existing {t.wt}. The merge is left in "
                      "progress")  # fmt: skip
    copied = packages[0] / cores[0]
    shutil.copyfile(build / cores[0], copied)
    copied.touch()
    return True


def message(t: Tree, args: argparse.Namespace, state: dict[str, Any], *, rebuilt: bool) -> str:
    prs, count = state["prs"], state["count"]
    more = f" and {len(prs) - 6} more" if len(prs) > 6 else ""
    brings = ", ".join(f"#{n}" for n in prs[:6]) + more if prs else f"{count} commits"
    short = t.out("rev-parse", "--short", state["master"])
    base = t.out("rev-parse", "--short", state["base"])
    branch = t.out("rev-parse", "--abbrev-ref", "HEAD")
    tree = t.out("merge-tree", "--write-tree", "HEAD", "MERGE_HEAD").splitlines()[0]
    conflicts = state["conflicts"]
    hand = [f for f in t.lines("diff", "--cached", "--name-only", tree) if f not in conflicts]
    lines = [f"Merge origin/master {short} into {branch}: brings in {brings} (@{args.persona})", "",
             f"Merge base {base}. First-parent merges brought in",
             f"(git log --merges --first-parent --oneline {base}..{short}):",
             *(f"  {m}" for m in state["merges"] or ["(none)"]),
             f"Conflicts resolved by hand: {', '.join(conflicts)}" if conflicts
             else "Conflicts: none",
             *([f"Changed by hand beyond git's merge: {', '.join(hand)}"] if hand else []),
             f"Gates on the merged tree, in {t.wt}/.venv:", f"{GATE_NAMES}.",
             "_core rebuilt: yes." if rebuilt else "_core rebuilt: no, no C++ came in.",
             "", args.trailer]  # fmt: skip
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None, run: Run = run) -> int:
    parser = argparse.ArgumentParser(prog="merge_master", description=__doc__.split("\n")[0])
    parser.add_argument("worktree", type=Path)
    form = parser.add_mutually_exclusive_group()
    form.add_argument("--continue", dest="resume", action="store_true")
    form.add_argument("--abort", action="store_true")
    parser.add_argument("--persona")
    parser.add_argument("--trailer")
    parser.add_argument("--no-build", action="store_true")
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)
    wt: Path = args.worktree
    try:
        if not args.abort:
            for option in ("persona", "trailer"):
                if getattr(args, option) is None:
                    raise Refusal(f"--{option} is missing; starting or continuing a merge needs "
                                  "--persona and --trailer")  # fmt: skip
            writes = load(HERE / "brief.py").WRITES
            writers = [name for name, places in writes.items() if places]
            if args.persona not in writers:
                raise Refusal(f"--persona {args.persona}: give one of {', '.join(writers)}")
            if not TRAILER.match(args.trailer):
                raise Refusal("--trailer must be the Co-Authored-By line from your system context")
        check_checkout(wt, run)
        t = Tree(wt, run)
        branch = t.out("rev-parse", "--abbrev-ref", "HEAD")
        if branch in ("HEAD", "master"):
            where = "master" if branch == "master" else "a detached HEAD"
            raise Refusal(f"{wt} is on {where}; run it on a worktree's own branch")
        if args.abort:
            if not t.merge_head.exists():
                raise Refusal(f"no merge in progress in {wt}; nothing to abort")
            done = t.git("merge", "--abort", capture=False)
            if done.returncode != 0:
                raise Refusal(f"git merge --abort failed (exit {done.returncode})")
            t.drop_state()
            return 0
        if not (wt / ".venv" / "bin" / "python").exists():
            raise Refusal(f"{wt} has no .venv; run: python3 tools/new_worktree.py --existing {wt}")
        return finish(t, args) if args.resume else start(t, args)
    except Refusal as exc:
        print(f"merge_master: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
