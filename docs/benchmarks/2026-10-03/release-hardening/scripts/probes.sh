#!/bin/bash
# The checks that the hardened builds are the objects measured, and that the
# define does what it says. Run from the worktree root after the four builds;
# writes nothing. Output committed as raw/probes.txt.
W=/Users/skavhaug/projects/rasputin/.claude/worktrees/agent-a0ff1ca8678bdc0c3
S=$(dirname "$0")
echo "## compile flags of _core per build (CMakeFiles/_core.dir/flags.make)"
for d in build-clang-plain build-clang-fast build-park-gplain build-bench; do
  echo "$d ($(cat "$W/$d/.variant")): $(grep '^CXX_FLAGS =' "$W/$d/CMakeFiles/_core.dir/flags.make")"
  echo "  compiler: $(grep '^CMAKE_CXX_COMPILER:' "$W/$d/CMakeCache.txt")"
done
echo
echo "## size, brk (trap) instructions, sha256 of each _core"
bash "$S/so_probe.sh" build-clang-plain build-clang-fast build-park-gplain build-bench
echo
echo "## libstdc++ assertion handler imported (nm -u | grep -c __glibcxx_assert_fail)"
for d in build-park-gplain build-bench; do
  echo "$d ($(cat "$W/$d/.variant")): $(nm -u "$W/$d"/_core*.so | grep -c glibcxx_assert_fail)"
done
echo
echo "## oob.cpp: one-past-the-end vector/span index, -O3 -DNDEBUG"
T=$(mktemp -d)
c++ -std=c++20 -O3 -DNDEBUG "$S/oob.cpp" -o "$T/a" && "$T/a"; echo "clang plain: exit $?"
c++ -std=c++20 -O3 -DNDEBUG -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_FAST "$S/oob.cpp" -o "$T/b" && "$T/b"; echo "clang FAST: exit $? (133 = SIGTRAP)"
g++-16 -std=c++20 -O3 -DNDEBUG "$S/oob.cpp" -o "$T/c" && "$T/c"; echo "gcc plain: exit $?"
g++-16 -std=c++20 -O3 -DNDEBUG -D_GLIBCXX_ASSERTIONS "$S/oob.cpp" -o "$T/d" 2>&1 && "$T/d" 2>&1; echo "gcc ASSERTIONS: exit $? (134 = SIGABRT)"
rm -rf "$T"
c++ --version | head -1; g++-16 --version | head -1
