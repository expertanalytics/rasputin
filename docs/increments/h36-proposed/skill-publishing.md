---
name: publishing
description: Pushing, opening and merging pull requests, CI and the merge queue. Load before asking Ola to push or merge, before any push, PR or merge, and before calling a branch merge-ready.
---

# Publishing

* **One yes, one act.** Ola's yes covers one occurrence of the act it names:
  done once and the tree changed since, ask again. A yes recorded on disk is
  asked again after a context loss.
* **`@reviewer` first,** before a branch's first push and before a push that
  follows a verdict other than `APPROVED`, on any branch.
* **Conflicts before the question:** `git fetch`, then
  `git merge-tree --write-tree origin/master <branch>`, and against each
  branch queued ahead. Master conflicts: a writer runs
  `tools/merge_master.py`; a branch ahead conflicts: the question names the
  files.
* **Merge queue:** `gh pr merge <pr>` enqueues. Never `--admin`.
* **CI:** `CI result` passes a prose-only change (`tools/ci_changes.py`
  decides) on the governance gates alone; that green says nothing about code.
