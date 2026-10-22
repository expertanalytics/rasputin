# Increment 88: the parallel helper keeps its threads between calls

**Status: dropped, not built** (Ola's ruling, 2026-10-22). The request assumed
about 35,000 parallel calls a year from the private package's step loop; the
year makes 4,684 (35,000 was threads started, 8 a call; section 1), so a
persistent team would save only about 1 to 1.5 s a year. The evaluator's time
was per-receiver trigonometry, which the private package's increment 89 (the
sun evaluator speed-up) sped up instead. The design below is kept unchanged
because that increment cites it. Designed by `@architect`, 2026-10-21, on
master `3b21eb59`, branch `worktree-sun-threads`. Full path: it changes a C++
header that refine and the mesh scan run through, so `@perf`'s acceptance run
applies
(`docs/increments/README.md`, *Acceptance*). Number 88 because 87 is the
highest found in `docs/increments/` over every local branch, every
remote-tracking ref and every worktree under `.claude/worktrees/`, checked
2026-10-21; 86 and 87 are the private package's.

## 1. What is asked, and why

Ola said yes to question 45 on 2026-10-21: the core's parallel helper should
reuse its threads instead of starting new ones on every call. His words were
relayed by the main session and are not quoted here.

`parallel_util::for_each_block` (`include/terrain/parallel_util/chunks.hpp`)
starts `min(threads, blocks)` new `std::jthread`s on every call above its
inline threshold and joins them before returning. Refine calls it a few dozen
times a run, so that never mattered. The private package's step loop calls it
once per step and kernel: the sun-and-shadow evaluator, the snow kernel and
the new evaporation kernel each reach it through the same header, compiled
into their own extension module. Every caller gains from a change here
without changing its own code.

**What the evidence says, and what it does not.** From the year-run profile
(`git show ee9565d3:docs/benchmarks/2026-10-21/profile-year/README.md`, AC
power, M1 Max):

- macOS `sample` saw 14,241 distinct threads in 30 s of the 78 s step loop,
  about 37,000 over the loop if the rate held.
- The evaluator was called 4,684 times (cProfile), 3.95 ms a call. At 8
  threads a call that is about 37,000 threads: the "about 35,000" figure the
  request gave is threads started, not calls. Checked with
  `git show ee9565d3:docs/benchmarks/2026-10-21/profile-year/raw/pstats-listing.txt | grep evaluate_instant`.
- A start-and-join of 8 threads cost about 90 µs plus about 25 µs of join on
  this machine (`docs/benchmarks/2026-09-27/serial-profile/README.md`, lines
  229-231: 3.8 ms over 41 starts, join 1.0 ms; battery power).

So, an estimate, not measured: thread starts cost the evaluator about
4,684 × 115 µs ≈ 0.5 s of its 18.5 s, and the evaporation kernel (one call per
step since private increment 86, 8,784 steps) about 1 s more. The pool takes
back most of that, about 1-1.5 s of the year run. It does not make the
evaluator fast: its time is the per-receiver trigonometry the profile lists
(`atan2`, `tan`, `sincos`, `pow`), which is the private package's to fix. The
micro-benchmark in section 7 replaces this estimate with a measurement.

## 2. Prior art: legacy and literature

**Legacy.** Nothing to carry across:

```sh
$ git grep -n -i -e "thread_pool" -e "threadpool" -e "std::thread" -e "omp " -e "tbb" legacy-archive -- legacy
$ git grep -l -e "thread" legacy-archive -- legacy
```

Both return no files (the tag `a27db62d` holds 30 files under `legacy/`, so the
grep ran over a real tree). Increment 21's wider grep found the same
(`docs/increments/21-parallel-refine.md`, lines 64-70).

**Literature.** The helper is fork-join parallelism with dynamic loop
scheduling: OpenMP's `schedule(dynamic, b)` over a team of threads (OpenMP
Architecture Review Board, *OpenMP API Specification 5.2*, 2021, the `schedule`
clause of the worksharing-loop construct). A team that persists between
parallel regions is how every OpenMP runtime implements it (libgomp's and
LLVM libomp's "hot team"). This increment does the same: a persistent team,
woken per call, with the existing shared counter handing out blocks. No work
stealing (Blumofe and Leiserson, "Scheduling multithreaded computations by
work stealing", *JACM* 46(5), 1999; Lea, "A Java Fork/Join Framework", Java
Grande 2000): the loops here are flat, with no nested tasks to steal, and the
shared counter already balances them (increment 21a). No novelty is claimed.

The options, with what each costs us:

| option | what it is | cost to us |
|---|---|---|
| `std::jthread` per call (today) | start and join per call | about 115 µs per 8-thread call on the M1; nothing outlives a call |
| persistent fork-join team over `std::jthread` (chosen) | threads sleep on a condition variable between calls, the caller works as one of the team | one process-wide object per shared library, its shutdown and fork handling (section 4); no dependency |
| general task-queue pool (Williams, *C++ Concurrency in Action*, 2nd ed., 2019, ch. 9; `BS::thread_pool`) | a queue of type-erased tasks, a future per task | an allocation and a future per block, more than fork-join needs; nested waits deadlock unless tasks are run while waiting |
| oneTBB `task_arena` and `parallel_for` | work-stealing arena, a process-wide market of threads | a new compiled dependency in both wheels; its default partitioner changes the block partition (`simple_partitioner` with a grain keeps it); refused in increment 14 R7 for `std::execution::par` on the same ground |
| OpenMP `parallel for` | the hot team above, from the runtime | Apple Clang ships no OpenMP: `libomp` from Homebrew, and two OpenMP runtimes in one Python process (NumPy's BLAS may bring one) is a known crash; refused in increment 14 R7 |
| spinning workers (OpenMP's `OMP_WAIT_POLICY=active`) | workers spin between calls instead of sleeping | lower wake latency, but the private loop spends most of its time in single-threaded NumPy between calls; 7 spinning threads would take the cores and the battery that work needs. Not taken |

**A ruling this reverses.** Increment 14 R7 said "No pool: a pool is state that
outlives a call, which is the one thing that makes this hard to reason about"
(`docs/increments/14-adaptive-refinement.md`, lines 254-256), when calls
numbered in the tens. They now number in the thousands, which is why Ola's yes
was asked. What made a pool hard to reason about is kept out by section 4: no
result depends on which thread runs a block (already true, and the contract
says so), a busy or foreign pool is never waited on, and the state the pool
keeps is threads asleep, nothing a caller can see.

**Ideas file** (`rasputin_scratch/ideas/ideas.md`, grepped for thread, pool,
parallel, speed): no entry this design takes up or overlaps.

## 3. The blueprint

```
for_each_block(n, threads, schedule, fn)        chunks.hpp, template, unchanged signature
  ├─ resolve threads, block size b, blocks, workers = min(threads, blocks)   unchanged
  ├─ workers <= 1 or n < inline_below  ──► inline, ascending               unchanged
  ├─ drain: loop { k = next++ ; if k >= blocks stop ; run(k), catch into errors[k] }
  ├─ detail::pool().try_run({&drain}, workers - 1)
  │     true  ─► caller and workers-1 pool threads each ran drain; all returned
  │     false ─► spawn workers std::jthreads running drain, join      (today's path)
  └─ rethrow errors[lowest k]                                          unchanged
```

The partition, the counter, the per-block error slots and the rethrow order
are today's lines, untouched. Only who runs `drain` changes: today `workers`
fresh threads, after this the caller plus `workers - 1` pool threads.

### 3.1 The interface

A new header, `include/terrain/parallel_util/pool.hpp`, non-template, included
by `chunks.hpp`:

```cpp
namespace terrain::parallel_util {

// OS threads this helper has started in this shared library, pool and
// fallback together. For tests and benchmarks; monotonic.
[[nodiscard]] std::uint64_t thread_starts() noexcept;

namespace detail {

// One call's work: run(context) is called once on each participating thread
// and must not throw.
struct Job {
    void (*run)(void*) noexcept;
    void* context;
};

class Pool {
public:
    // Runs job on the calling thread and on `helpers` pool threads at once,
    // starting threads until the pool holds `helpers`, and returns once every
    // one has returned from job.run; their writes happen-before the return.
    // Returns false at once, having run nothing, when the pool cannot take the
    // call: another call is in it (a second calling thread, or a nested call
    // from inside a job), the process is a fork of the one that made the pool,
    // or the pool has shut down. A thread that cannot start throws
    // std::system_error, as std::jthread does today; the pool stays usable.
    bool try_run(Job job, std::size_t helpers);
};

Pool& pool() noexcept;  // this shared library's one pool, made on first use

}  // namespace detail
}  // namespace terrain::parallel_util
```

`thread_starts` is the one public addition, and the only way a test can see a
thread start; it costs one relaxed atomic increment per start.

### 3.2 The pool's mechanics

- **State:** a mutex, two condition variables (`wake` for the team, `done`
  for the caller), a grow-only `std::vector<std::jthread>`, the current
  `Job`, `helpers` and `pending` counts, a generation counter, a `stop` flag;
  outside the mutex, `std::atomic<bool> busy` and the owner's process ID.
- **Taking the pool:** `busy.compare_exchange_strong(false, true)`. A flag, not
  a `try_lock`, because a nested call comes from a thread that may already
  hold the pool, and `try_lock` by the owning thread is undefined. A failed
  exchange returns false; so does `getpid() != owner`.
- **A call:** grow the vector to `helpers` threads (each start counted); under
  the mutex store the job, set `helpers` and `pending = helpers`, bump the
  generation; `notify_all` on `wake`; run `job.run` on the caller; then wait on
  `done` until `pending == 0`; clear the job; release `busy`.
- **A pool thread** with index `i` waits on `wake` for a new generation or
  `stop`. On a new generation it runs the job only if `i < helpers`, then
  under the mutex decrements `pending` and notifies `done` at zero. Threads
  with `i >= helpers` go back to sleep, so `--threads 4` after an earlier
  `--threads 8` call still runs 4 at once.
- **Publication:** every write `fn` made on a pool thread is before that
  thread's decrement of `pending` under the mutex, and the caller reads
  `pending == 0` under the same mutex; that is the happens-before the join gives
  today, and one TSan can see.
- **Never waits on a busy pool.** The fallback is today's spawn-and-join, so a
  second calling thread (the catchment batch runs `delineate` in
  `asyncio.to_thread`, concurrently) or a nested call costs what it costs today
  and cannot deadlock on the pool.
- **No result depends on the pool.** Which thread runs a block was already free
  (the header says so, and `test_refinement_chunks_dynamic` does not pin it
  above the threshold); no caller uses `thread_local` (`grep -rn thread_local
  include src bindings` finds nothing). So refine, the slope and steepness
  passes and the private kernels give the same bits for every thread count, as
  T6 already pins.

### 3.3 Lifetime: one pool per shared library, shut down at exit

- `pool()` holds `static Pool* p = new Pool;`, never deleted, and a
  `static Reaper` beside it. The pointer outlives every static destructor, so a
  call from another static destructor after the reaper has run is safe: it
  finds the pool shut down and takes the fallback.
- **The reaper's destructor** runs at process exit (after `Py_Finalize`, when
  the C++ runtime runs static destructors; Python never unloads an extension
  module). In the owning process it takes `busy` (yielding until any in-flight
  call ends), sets `stop`, wakes the team and joins it, and leaves `busy` set
  for good. In a forked child it does nothing: the threads exist only in the
  parent, and joining their copies is undefined.
- **Per shared library.** `_core` and the private package's `_lsm` are each
  built with `pybind11_add_module`, which hides symbols, so each has its own
  pool; a process using both holds two sleeping teams. Where symbols are not
  hidden (a test binary, or Linux's unique inline statics), two libraries may
  share one pool; the busy flag keeps that correct.
- **Fork.** A child forked after the pool exists sees `getpid() != owner` and
  takes the fallback for good. Python 3.12 and later warns at `os.fork()` in a
  process with more than one thread; after this change that warning appears
  once any parallel call has run, where today it appears only during one.
  `tests/python/test_away.py` forks a pty and execs at once, so it is
  unaffected beyond the warning.
- POSIX only (`getpid` from `<unistd.h>`). CI builds on Ubuntu and macOS only.

## 4. Contract kept, line by line

| asked to keep | how |
|---|---|
| same block partition, bit-identical output for every thread count | partition, counter and error slots unchanged (section 3); T6 and T16 of `prop_refinement_refine` pin refine; section 6 adds a partition check through the pool |
| refine and mesh determinism | the same; no caller depends on the thread identity |
| the thread sanitizer passing | mutex-and-condition-variable publication; the new suite joins the TSan job's lists |
| `--threads` and `RASPUTIN_THREADS` | untouched: they set the `threads` argument, and a call uses `min(threads, blocks)` participants however large the pool has grown |
| no global state that breaks plugins or tests | one pool per shared library, never waited on when busy; nothing it keeps is visible to a caller |
| a clean shutdown at interpreter exit | the reaper (section 3.3); section 6, P7 |

The header comment of `chunks.hpp` ("No pool. Threads are created per call…")
and the `parallel_util` section of `project_structure.md` ("over
`std::jthread` created per call and joined before it returns, no pool") are
rewritten in the green commit to say what the code then does.

## 5. Left out

- **Retuning `inline_below` (256).** It came from a sweep with per-call starts
  (increment 21a); with cheap wakes a lower threshold may pay. Changing it
  changes which rounds go parallel, so refine's benchmark would move for a
  second reason. A follow-up, with its own sweep.
- **Shrinking the pool**, an idle timeout, or a cap below what callers ask
  for. Sleeping threads cost no CPU; the largest `threads` a process asked for
  bounds the team.
- **A Python view of `thread_starts`.** No binding change.
- **Making the evaluator itself fast.** The private package's work (section 1).

## 6. The tests that pin it

A new suite, `tests/cpp/unit/test_parallel_pool.cpp` (target
`test_parallel_pool`), added to the TSan job's build and run lists in
`.github/workflows/main.yaml`. Fixtures are synthetic: `n` items, each block
adding 1 to a per-index counter, sizes as given. Every case reads
`thread_starts()` as a difference, never as an absolute, so the cases pass in
any order and in one process.

- **P1, reuse.** One call with `n = 100,000`, `threads = 8`, block 0 (so 625
  blocks), then 200 more identical calls: the 200 start no thread
  (`thread_starts()` unchanged), and every index is hit exactly once in each.
- **P2, growth only to what is asked.** After a call at `threads = 3`, a call at
  `threads = 6` starts at most 3 more; a following call at `threads = 2` starts
  none.
- **P3, participants honour `threads`.** After a call at `threads = 8`, a call at
  `threads = 3`, `n = 100,000`, block 1,000, sees at most 3 distinct thread IDs
  calling `fn`, and with the blocks waiting for each other (as the existing
  "runs blocks concurrently" case does, 10 s timeout) exactly 3.
- **P4, partition through the pool.** The existing case "calls fn once per block
  of the partition", rerun twice in a row in one case for `threads` in
  {2, 7, 8}, `n` in {1,009, 100,000}: the second run goes through the warm pool
  and gives the same sorted call list.
- **P5, nested call.** `fn` itself calls `for_each_block(64, 4, {1, 0}, …)`; the
  outer `n = 1,000`, `threads = 4`, block 100: completes within 30 s (timed
  from a watchdog thread), every inner and outer index hit once.
- **P6, two calling threads.** Two `std::thread`s each make 50 calls with
  `n = 10,000`, `threads = 4`: every call's indices hit once, no hang (30 s
  watchdog). TSan runs it.
- **P7, exit.** Python: a subprocess imports `tin_engine._core`, runs one
  parallel refine of the smallest refine fixture at `threads = 4`, and exits;
  return code 0 within 60 s. A hang in the reaper fails it.
- **P8, fork.** A separate target `test_parallel_pool_fork`, not in the TSan
  job (TSan stops a child that starts threads after a multi-threaded fork):
  warm the pool at `threads = 4`, `fork()`, and in the child make one call at
  `threads = 4`, `n = 10,000` and `_exit(0)` if every index was hit once; the
  parent waits up to 30 s and requires status 0.
- **Existing, unchanged and required green:**
  `test_refinement_chunks_dynamic` (its exception-order cases now run through
  the pool), T6 and T16 of `prop_refinement_refine`, and the TSan job's whole
  list.

**Invariant-critical:** P1, P4 and the existing chunks suite. Mutation targets
for the red step's kill record: a pool thread that runs when `i >= helpers`
(P3 fails); `pending` decremented before the job ends, the caller returning
early (P1 and P4 fail under TSan or on counts); `busy` never released (P1
fails: every later call spawns); the fallback path removed, so a busy pool is
waited on (P5 or P6 hangs, caught by the watchdog).

## 7. Acceptance (`@perf`)

- **Refine must not regress.** `tools/bench.py run` on the 1 m benchmark and
  its thread sweep, against this branch's merge base measured back to back,
  same power state (`docs/increments/README.md`, *Acceptance*). A refine run
  makes a few dozen parallel calls, so the expected change is within noise;
  any slowdown outside the noise band fails.
- **The micro-benchmark.** A hidden Catch2 case, `[.][bench]`, in
  `test_parallel_pool.cpp`: 10,000 calls of `for_each_block(55,000, 8, {}, fn)`
  with a trivial `fn` (one add per index), wall clock per call, median of 5
  runs. It uses only `for_each_block`, so `@perf` compiles the same file
  against the merge base's header for the before figure. 55,000 is the year
  run's triangle count. Expected from section 1: before about 115 µs a call;
  after, the cost of a condition-variable wake, which this design does not
  guess. Recorded under `docs/benchmarks/<date>/`.
- **The private year run**, timed by `@perf` before and after: section 1's
  estimate is 1-1.5 s.

## 8. Size

Production code: `pool.hpp` about 75 counted lines, `chunks.hpp` about +10
(the `drain` lambda moves into a `Job`, the fallback keeps today's lines);
about 85 net. Tests (`tests/`) and the workflow YAML are not counted.
