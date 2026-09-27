"""chunksum.py FILES...: per file, scan totals from the R and C records."""
import sys
for path in sys.argv[1:]:
    rounds = []
    for line in open(path):
        if line[0] == "R":
            f = line.split(); rounds.append([float(f[3]), []])
        elif line[0] == "C":
            f = line.split(); rounds[-1][1].append((float(f[3]), float(f[4]), int(f[5])))
    scan = sum(r[0] for r in rounds)
    mx = sum(max(c[1] - c[0] for c in r[1]) for r in rounds)
    busy = sum(sum(c[1] - c[0] for c in r[1]) for r in rounds)
    mean = sum(sum(c[1] - c[0] for c in r[1]) / len(r[1]) for r in rounds)
    start = sum(max(c[0] for c in r[1]) for r in rounds)
    first = sum(min(c[0] for c in r[1]) for r in rounds)
    tail = sum(r[0] - max(c[1] for c in r[1]) for r in rounds)
    print(f"{path}\tscan_ms={1e3*scan:.1f}\tsum_max_chunk_ms={1e3*mx:.1f}\tsum_mean_chunk_ms={1e3*mean:.1f}\tthread_busy_ms={1e3*busy:.1f}\tfirst_start_ms={1e3*first:.2f}\tlast_start_ms={1e3*start:.2f}\tjoin_tail_ms={1e3*tail:.2f}")
