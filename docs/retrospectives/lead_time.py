"""Lead time per pull request, from git alone (Ola's ask of 2026-10-07).

For each PR merged into origin/master since a date (a local master can be stale), prints:
  start  the first commit on the PR's branch (author date)
  built  the last commit on the branch (author date): the work is done
  merged the merge commit on master (committer date)
  build  start -> built, in hours: design, tests, code, review, records
  wait   built -> merged, in hours: Ola's yes, the push, CI and the queue
Branches not merged yet are given by name; their "wait" runs to now.

    python3 docs/retrospectives/lead_time.py 2026-10-01 [branch ...]

Times are local (Oslo). Start counts from the branch's first commit, so a
design written on an earlier branch (20c's, on #207) is not in a later PR's
build time; per increment, use the first commit that touched its file:
    git log --all --reverse --format=%ad -- docs/increments/<file> | head -1
"""

import subprocess
import sys
import time


def git(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True, check=True).stdout


def span(base: str, tip: str) -> tuple[int, int, int]:
    stamps = [int(s) for s in git("log", "--no-merges", "--format=%at", f"{base}..{tip}").split()]
    return min(stamps), max(stamps), len(stamps)


def row(name: str, start: int, built: int, end: int, n: int, merged: str) -> None:
    clock = time.strftime("%m-%d %H:%M", time.localtime(start))
    build, wait = (built - start) / 3600, (end - built) / 3600
    print(f"{name:34} {clock}  {build:6.1f}  {wait:5.1f}  {build + wait:6.1f}  {n:3}  {merged}")


def main() -> None:
    since, pending = sys.argv[1], sys.argv[2:]
    print(f"{'PR':34} {'start':11}  {'build':>6}  {'wait':>5}  {'total':>6}  {'n':>3}  merged")
    # "--since=<date>" alone means that date at the current time of day
    log = git("log", "--first-parent", "origin/master", f"--since={since} 00:00", "--merges",
              "--grep=Merge pull request", "--format=%H %ct %s")  # fmt: skip
    for line in reversed(log.splitlines()):
        sha, ct, subject = line.split(" ", 2)
        name = subject.removeprefix("Merge pull request ").replace("from expertanalytics/", "")
        start, built, n = span(f"{sha}^1", f"{sha}^2")
        row(name, start, built, int(ct), n, time.strftime("%m-%d %H:%M", time.localtime(int(ct))))
    for branch in pending:
        base = git("merge-base", "origin/master", branch).strip()
        start, built, n = span(base, branch)
        row(branch, start, built, int(time.time()), n, "not yet")


if __name__ == "__main__":
    main()
