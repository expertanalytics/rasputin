#!/bin/bash
# Size, trap-instruction count and sha256 of each built _core.
W=/Users/skavhaug/projects/rasputin/.claude/worktrees/agent-a0ff1ca8678bdc0c3
for d in "$@"; do
  f=$(ls "$W/$d"/_core*.so)
  printf '%s size=%s brk=%s sha256=%s\n' "$d" "$(stat -f %z "$f")" \
    "$(objdump -d --no-show-raw-insn "$f" | grep -cE '\bbrk\b')" \
    "$(shasum -a 256 "$f" | cut -d' ' -f1)"
done
