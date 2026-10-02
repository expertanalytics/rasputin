## Runs

| run | commit | dirty | power | pmset % | so sha256 | started |
|---|---|---|---|---|---|---|
| 15c-1-r1 | 6ab7ad5 | True | battery | 79 | 7dee5ad918d0 | 21:49:13 |
| 15c-1-r2 | 6ab7ad5 | True | battery | 78 | 7dee5ad918d0 | 21:54:02 |
| 15c-1-r3 | 6ab7ad5 | True | battery | 75 | 7dee5ad918d0 | 21:58:37 |
| 15c-1-r4 | 6ab7ad5 | True | battery | 74 | 7dee5ad918d0 | 22:01:45 |
| 15c-1-r5 | 6ab7ad5 | True | battery | 72 | 7dee5ad918d0 | 22:06:27 |
| base-a130f7c-r1 | a130f7c | False | battery | 81 | eb5f06e88a7c | 21:46:42 |
| base-a130f7c-r2 | a130f7c | False | battery | 79 | eb5f06e88a7c | 21:51:46 |
| base-a130f7c-r3 | a130f7c | False | battery | 76 | eb5f06e88a7c | 21:56:18 |
| base-a130f7c-r4 | a130f7c | False | battery | 73 | eb5f06e88a7c | 22:04:05 |
| base-a130f7c-r5 | a130f7c | False | battery | 71 | eb5f06e88a7c | 22:08:49 |

## Quality

| run | domain | worst angle | max degree | within tol | Delaunay checked / ambiguous / violations | mesh sha256 |
|---|---|---:|---:|---|---|---|
| 15c-1-r1 | quarter | 0.3955 | 18 | True | 641791 / 80136 / 0 | ccebf96a86c6c5e2 |
| 15c-1-r1 | tile | 0.6296 | 74 | True | 692056 / 94469 / 0 | 11741a81adfa17b3 |
| 15c-1-r2 | quarter | 0.3955 | 18 | True | 641791 / 80136 / 0 | ccebf96a86c6c5e2 |
| 15c-1-r2 | tile | 0.6296 | 74 | True | 692056 / 94469 / 0 | 11741a81adfa17b3 |
| 15c-1-r3 | quarter | 0.3955 | 18 | True | 641791 / 80136 / 0 | ccebf96a86c6c5e2 |
| 15c-1-r3 | tile | 0.6296 | 74 | True | 692056 / 94469 / 0 | 11741a81adfa17b3 |
| 15c-1-r4 | quarter | 0.3955 | 18 | True | 641791 / 80136 / 0 | ccebf96a86c6c5e2 |
| 15c-1-r4 | tile | 0.6296 | 74 | True | 692056 / 94469 / 0 | 11741a81adfa17b3 |
| 15c-1-r5 | quarter | 0.3955 | 18 | True | 641791 / 80136 / 0 | ccebf96a86c6c5e2 |
| 15c-1-r5 | tile | 0.6296 | 74 | True | 692056 / 94469 / 0 | 11741a81adfa17b3 |
| base-a130f7c-r1 | quarter | 0.3955 | 18 | True | 641791 / 80136 / 0 | ccebf96a86c6c5e2 |
| base-a130f7c-r1 | tile | 0.6296 | 74 | True | 692056 / 94469 / 0 | 11741a81adfa17b3 |
| base-a130f7c-r2 | quarter | 0.3955 | 18 | True | 641791 / 80136 / 0 | ccebf96a86c6c5e2 |
| base-a130f7c-r2 | tile | 0.6296 | 74 | True | 692056 / 94469 / 0 | 11741a81adfa17b3 |
| base-a130f7c-r3 | quarter | 0.3955 | 18 | True | 641791 / 80136 / 0 | ccebf96a86c6c5e2 |
| base-a130f7c-r3 | tile | 0.6296 | 74 | True | 692056 / 94469 / 0 | 11741a81adfa17b3 |
| base-a130f7c-r4 | quarter | 0.3955 | 18 | True | 641791 / 80136 / 0 | ccebf96a86c6c5e2 |
| base-a130f7c-r4 | tile | 0.6296 | 74 | True | 692056 / 94469 / 0 | 11741a81adfa17b3 |
| base-a130f7c-r5 | quarter | 0.3955 | 18 | True | 641791 / 80136 / 0 | ccebf96a86c6c5e2 |
| base-a130f7c-r5 | tile | 0.6296 | 74 | True | 692056 / 94469 / 0 | 11741a81adfa17b3 |

## Refine time, median of the per-run medians (5 base runs, 5 branch runs)

threads 0 is the CLI default (hardware_concurrency).

| domain | threads | base a130f7c (s) | base min-max | 15c-1 (s) | 15c-1 min-max | change |
|---|---:|---:|---|---:|---|---:|
| quarter | 0 | 0.1702 | 0.1652-0.1769 | 0.1696 | 0.1646-0.1758 | -0.3 % |
| quarter | 1 | 0.4723 | 0.4482-0.4837 | 0.4597 | 0.4232-0.4763 | -2.7 % |
| quarter | 2 | 0.2959 | 0.2801-0.2996 | 0.2946 | 0.2790-0.3058 | -0.4 % |
| quarter | 3 | 0.2382 | 0.2272-0.2427 | 0.2364 | 0.2271-0.2509 | -0.8 % |
| quarter | 4 | 0.2099 | 0.2046-0.2129 | 0.2077 | 0.1965-0.2140 | -1.0 % |
| quarter | 5 | 0.1914 | 0.1871-0.1970 | 0.1926 | 0.1818-0.1957 | +0.7 % |
| quarter | 6 | 0.1811 | 0.1778-0.1834 | 0.1803 | 0.1747-0.1819 | -0.5 % |
| quarter | 7 | 0.1740 | 0.1714-0.1758 | 0.1721 | 0.1691-0.1729 | -1.1 % |
| quarter | 8 | 0.1691 | 0.1634-0.1742 | 0.1712 | 0.1673-0.1776 | +1.2 % |
| quarter | 9 | 0.1714 | 0.1646-0.1740 | 0.1684 | 0.1627-0.1780 | -1.8 % |
| quarter | 10 | 0.1692 | 0.1646-0.1744 | 0.1725 | 0.1622-0.1867 | +2.0 % |
| quarter | 11 | 0.1706 | 0.1667-0.1735 | 0.1697 | 0.1630-0.1736 | -0.6 % |
| quarter | 12 | 0.1685 | 0.1620-0.1736 | 0.1708 | 0.1613-0.1767 | +1.4 % |
| quarter | 13 | 0.1706 | 0.1625-0.1744 | 0.1701 | 0.1618-0.1736 | -0.3 % |
| quarter | 14 | 0.1688 | 0.1650-0.1752 | 0.1695 | 0.1616-0.1722 | +0.5 % |
| quarter | 15 | 0.1688 | 0.1630-0.1735 | 0.1706 | 0.1625-0.1746 | +1.0 % |
| quarter | 16 | 0.1692 | 0.1643-0.1760 | 0.1678 | 0.1628-0.1739 | -0.9 % |
| quarter | 17 | 0.1699 | 0.1645-0.1737 | 0.1688 | 0.1641-0.1725 | -0.6 % |
| quarter | 18 | 0.1709 | 0.1652-0.1735 | 0.1698 | 0.1641-0.1741 | -0.7 % |
| quarter | 19 | 0.1722 | 0.1648-0.1750 | 0.1666 | 0.1627-0.1729 | -3.2 % |
| quarter | 20 | 0.1716 | 0.1653-0.1757 | 0.1692 | 0.1623-0.1742 | -1.4 % |
| tile | 0 | 0.1884 | 0.1869-0.1943 | 0.1952 | 0.1902-0.1967 | +3.6 % |
| tile | 1 | 0.5102 | 0.4748-0.5371 | 0.5174 | 0.5133-0.5444 | +1.4 % |
| tile | 2 | 0.3363 | 0.3145-0.3425 | 0.3343 | 0.3247-0.3421 | -0.6 % |
| tile | 3 | 0.2722 | 0.2637-0.2790 | 0.2720 | 0.2674-0.2921 | -0.1 % |
| tile | 4 | 0.2363 | 0.2310-0.2393 | 0.2385 | 0.2350-0.2443 | +0.9 % |
| tile | 5 | 0.2153 | 0.2101-0.2216 | 0.2187 | 0.2173-0.2287 | +1.6 % |
| tile | 6 | 0.2054 | 0.1994-0.2094 | 0.2077 | 0.2049-0.2117 | +1.1 % |
| tile | 7 | 0.1974 | 0.1925-0.2073 | 0.1969 | 0.1955-0.2005 | -0.3 % |
| tile | 8 | 0.1925 | 0.1857-0.1982 | 0.1921 | 0.1894-0.2011 | -0.2 % |
| tile | 9 | 0.1925 | 0.1849-0.1951 | 0.1934 | 0.1862-0.1980 | +0.5 % |
| tile | 10 | 0.1927 | 0.1846-0.1943 | 0.1917 | 0.1890-0.1955 | -0.5 % |
| tile | 11 | 0.1930 | 0.1846-0.1944 | 0.1928 | 0.1883-0.1974 | -0.1 % |
| tile | 12 | 0.1891 | 0.1855-0.1951 | 0.1953 | 0.1889-0.1965 | +3.3 % |
| tile | 13 | 0.1905 | 0.1900-0.1945 | 0.1909 | 0.1891-0.1952 | +0.2 % |
| tile | 14 | 0.1934 | 0.1880-0.1966 | 0.1923 | 0.1882-0.1999 | -0.6 % |
| tile | 15 | 0.1922 | 0.1864-0.2019 | 0.1899 | 0.1878-0.1990 | -1.2 % |
| tile | 16 | 0.1938 | 0.1909-0.2067 | 0.1904 | 0.1892-0.1991 | -1.8 % |
| tile | 17 | 0.1953 | 0.1849-0.1997 | 0.1938 | 0.1892-0.1958 | -0.8 % |
| tile | 18 | 0.1945 | 0.1853-0.1957 | 0.1937 | 0.1894-0.1985 | -0.5 % |
| tile | 19 | 0.1949 | 0.1880-0.1980 | 0.1905 | 0.1887-0.1990 | -2.3 % |
| tile | 20 | 0.1927 | 0.1854-0.1968 | 0.1949 | 0.1883-0.2046 | +1.1 % |

Over 42 cells: median change -0.39 %, range -3.2 % to +3.6 %, 0 cells above +5 %, 0 below -5 %.

## Per pair: median over the 21 thread cells of (branch / base - 1)

| pair | order | quarter | tile | verdict stored in the branch run.json |
|---|---|---:|---:|---|
| 1 | base first | +0.0 % | +1.0 % | REGRESSION: quarter refine_s[t=2] 0.2801 -> 0.2946 (+5.2 %); REGRESSION: quarter refine_s[t=10] 0.1646 -> 0.1867 (+13.4 %); REGRESSION: quarter refine_s[t=12] 0.1632 -> 0.1767 (+8.3 %) |
| 2 | base first | -2.3 % | +1.0 % | REGRESSION: tile refine_s[t=1] 0.4748 -> 0.5444 (+14.7 %); REGRESSION: tile refine_s[t=3] 0.2670 -> 0.2921 (+9.4 %); REGRESSION: tile refine_s[t=5] 0.2129 -> 0.2287 (+7.4 %); REGRESSION: tile refine_s[t=12] 0.1855 -> 0.1953 (+5.3 %) |
| 3 | base first | -0.8 % | +3.9 % | REGRESSION: tile refine_s[t=2] 0.3145 -> 0.3384 (+7.6 %); REGRESSION: tile refine_s[t=10] 0.1846 -> 0.1955 (+5.9 %); REGRESSION: tile refine_s[t=15] 0.1864 -> 0.1990 (+6.8 %); REGRESSION: tile refine_s[t=17] 0.1849 -> 0.1956 (+5.8 %); REGRESSION: tile refine_s[t=20] 0.1854 -> 0.2046 (+10.4 %) |
| 4 | branch first | -0.8 % | -2.7 % | NO BASELINE: no comparable stored run |
| 5 | branch first | +0.7 % | +1.4 % | NO BASELINE: no comparable stored run |

Pairs 4-5 ran branch first, so each branch run stored NO BASELINE; `bench.py compare` against its pair's base is in logs/pairs_reversed.log.

## Ceiling (pooled medians; 2026-09-26 reference 2.2x)

| domain | build | 1 -> 20 threads | best | at threads |
|---|---|---:|---:|---:|
| quarter | base a130f7c | 2.75x | 2.80x | 12 |
| quarter | 15c-1 | 2.72x | 2.76x | 19 |
| tile | base a130f7c | 2.65x | 2.70x | 12 |
| tile | 15c-1 | 2.65x | 2.72x | 15 |
Same build against itself (base-base and branch-branch run pairs, 840 cells): range -13.1 % to +13.1 %, 154 cells beyond +-5 %.
