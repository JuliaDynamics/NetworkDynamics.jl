# Benchmark Report for */var/folders/cx/04cwlxl9287_ndx30k6q3n1h0000gn/T/jl_lUxdLT/NetworkDynamics*

## Job Properties
* Time of benchmark: 26 Jul 2021 - 17:47
* Package commit: ebc9b4
* Julia commit: e76c9d
* Julia command flags: `-J/Users/hw/bin/julia/usr/lib/julia/sys.dylib`
* Environment variables: `JULIA_NUM_THREADS => 4`

## Results
Below is a table of this job's results, obtained by running the benchmarks.
The values listed in the `ID` column have the structure `[parent_group, child_group, ..., key]`, and can be used to
index into the BaseBenchmarks suite to retrieve the corresponding benchmarks.
The percentages accompanying time and memory values in the below table are noise tolerances. The "true"
time/memory value for a given benchmark is expected to fall within this percentage of the reported value.
An empty cell means that the value was zero.

| ID                                                   | time            | GC time | memory          | allocations |
|------------------------------------------------------|----------------:|--------:|----------------:|------------:|
| `["diffusion", "ode_edge", "assemble", "10"]`        |  59.000 μs (5%) |         |  50.97 KiB (1%) |         817 |
| `["diffusion", "ode_edge", "assemble", "100"]`       | 376.584 μs (5%) |         | 365.97 KiB (1%) |        5439 |
| `["diffusion", "ode_edge", "assemble", "1000"]`      |   3.647 ms (5%) |         |   3.50 MiB (1%) |       58760 |
| `["diffusion", "ode_edge", "assemble", "10000"]`     |  38.354 ms (5%) |         |  35.60 MiB (1%) |      637752 |
| `["diffusion", "ode_edge", "call", "10"]`            | 166.000 ns (5%) |         |                 |             |
| `["diffusion", "ode_edge", "call", "100"]`           |   1.833 μs (5%) |         |                 |             |
| `["diffusion", "ode_edge", "call", "1000"]`          |  24.042 μs (5%) |         |                 |             |
| `["diffusion", "ode_edge", "call", "10000"]`         | 318.792 μs (5%) |         |                 |             |
| `["diffusion", "ode_edge", "call_mt", "10"]`         |   8.917 μs (5%) |         |   4.62 KiB (1%) |          42 |
| `["diffusion", "ode_edge", "call_mt", "100"]`        |  10.709 μs (5%) |         |   4.62 KiB (1%) |          42 |
| `["diffusion", "ode_edge", "call_mt", "1000"]`       |  16.917 μs (5%) |         |   4.62 KiB (1%) |          42 |
| `["diffusion", "ode_edge", "call_mt", "10000"]`      |  91.583 μs (5%) |         |   4.62 KiB (1%) |          42 |
| `["diffusion", "static_edge", "assemble", "10"]`     |  58.083 μs (5%) |         |  49.14 KiB (1%) |         810 |
| `["diffusion", "static_edge", "assemble", "100"]`    | 370.458 μs (5%) |         | 350.56 KiB (1%) |        5432 |
| `["diffusion", "static_edge", "assemble", "1000"]`   |   3.649 ms (5%) |         |   3.36 MiB (1%) |       58744 |
| `["diffusion", "static_edge", "assemble", "10000"]`  |  38.442 ms (5%) |         |  34.13 MiB (1%) |      637738 |
| `["diffusion", "static_edge", "call", "10"]`         |  41.000 ns (5%) |         |                 |             |
| `["diffusion", "static_edge", "call", "100"]`        | 541.000 ns (5%) |         |                 |             |
| `["diffusion", "static_edge", "call", "1000"]`       |   6.958 μs (5%) |         |                 |             |
| `["diffusion", "static_edge", "call", "10000"]`      |  92.542 μs (5%) |         |                 |             |
| `["diffusion", "static_edge", "call_mt", "10"]`      |   8.375 μs (5%) |         |   4.73 KiB (1%) |          42 |
| `["diffusion", "static_edge", "call_mt", "100"]`     |   9.333 μs (5%) |         |   4.73 KiB (1%) |          42 |
| `["diffusion", "static_edge", "call_mt", "1000"]`    |  12.125 μs (5%) |         |   4.73 KiB (1%) |          42 |
| `["diffusion", "static_edge", "call_mt", "10000"]`   |  33.292 μs (5%) |         |   4.73 KiB (1%) |          42 |
| `["kuramoto", "heterogeneous", "assemble", "10"]`    |  62.625 μs (5%) |         |  50.11 KiB (1%) |         847 |
| `["kuramoto", "heterogeneous", "assemble", "100"]`   | 410.625 μs (5%) |         | 361.80 KiB (1%) |        5829 |
| `["kuramoto", "heterogeneous", "assemble", "1000"]`  |   3.987 ms (5%) |         |   3.48 MiB (1%) |       63489 |
| `["kuramoto", "heterogeneous", "assemble", "10000"]` |  42.039 ms (5%) |         |  35.32 MiB (1%) |      683008 |
| `["kuramoto", "heterogeneous", "call", "10"]`        | 583.000 ns (5%) |         |   1.09 KiB (1%) |          40 |
| `["kuramoto", "heterogeneous", "call", "100"]`       |   6.666 μs (5%) |         |  10.94 KiB (1%) |         400 |
| `["kuramoto", "heterogeneous", "call", "1000"]`      |  68.458 μs (5%) |         | 109.38 KiB (1%) |        4000 |
| `["kuramoto", "heterogeneous", "call", "10000"]`     | 702.917 μs (5%) |         |   1.07 MiB (1%) |       40000 |
| `["kuramoto", "heterogeneous", "call_mt", "10"]`     |  21.875 μs (5%) |         |   5.80 KiB (1%) |          81 |
| `["kuramoto", "heterogeneous", "call_mt", "100"]`    |  21.625 μs (5%) |         |  15.67 KiB (1%) |         442 |
| `["kuramoto", "heterogeneous", "call_mt", "1000"]`   |  67.125 μs (5%) |         | 114.11 KiB (1%) |        4042 |
| `["kuramoto", "heterogeneous", "call_mt", "10000"]`  | 251.125 μs (5%) |         |   1.07 MiB (1%) |       40043 |
| `["kuramoto", "homogeneous", "assemble", "10"]`      |  65.458 μs (5%) |         |  54.12 KiB (1%) |         892 |
| `["kuramoto", "homogeneous", "assemble", "100"]`     | 436.792 μs (5%) |         | 405.58 KiB (1%) |        6235 |
| `["kuramoto", "homogeneous", "assemble", "1000"]`    |   4.296 ms (5%) |         |   3.89 MiB (1%) |       68002 |
| `["kuramoto", "homogeneous", "assemble", "10000"]`   |  44.925 ms (5%) |         |  39.09 MiB (1%) |      728012 |
| `["kuramoto", "homogeneous", "call", "10"]`          | 208.000 ns (5%) |         |   64 bytes (1%) |           2 |
| `["kuramoto", "homogeneous", "call", "100"]`         |   2.042 μs (5%) |         |   64 bytes (1%) |           2 |
| `["kuramoto", "homogeneous", "call", "1000"]`        |  21.959 μs (5%) |         |   64 bytes (1%) |           2 |
| `["kuramoto", "homogeneous", "call", "10000"]`       | 267.334 μs (5%) |         |   64 bytes (1%) |           2 |
| `["kuramoto", "homogeneous", "call_mt", "10"]`       |  27.667 μs (5%) |         |   4.86 KiB (1%) |          44 |
| `["kuramoto", "homogeneous", "call_mt", "100"]`      |  22.708 μs (5%) |         |   4.86 KiB (1%) |          44 |
| `["kuramoto", "homogeneous", "call_mt", "1000"]`     |  54.750 μs (5%) |         |   4.86 KiB (1%) |          44 |
| `["kuramoto", "homogeneous", "call_mt", "10000"]`    | 132.042 μs (5%) |         |   4.89 KiB (1%) |          45 |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["diffusion", "ode_edge", "assemble"]`
- `["diffusion", "ode_edge", "call"]`
- `["diffusion", "ode_edge", "call_mt"]`
- `["diffusion", "static_edge", "assemble"]`
- `["diffusion", "static_edge", "call"]`
- `["diffusion", "static_edge", "call_mt"]`
- `["kuramoto", "heterogeneous", "assemble"]`
- `["kuramoto", "heterogeneous", "call"]`
- `["kuramoto", "heterogeneous", "call_mt"]`
- `["kuramoto", "homogeneous", "assemble"]`
- `["kuramoto", "homogeneous", "call"]`
- `["kuramoto", "homogeneous", "call_mt"]`

## Julia versioninfo
```
Julia Version 1.7.0-beta3
Commit e76c9dad42* (2021-07-07 08:12 UTC)
Platform Info:
  OS: macOS (arm64-apple-darwin20.5.0)
  uname: Darwin 20.5.0 Darwin Kernel Version 20.5.0: Sat May  8 05:10:31 PDT 2021; root:xnu-7195.121.3~9/RELEASE_ARM64_T8101 x86_64 i386
  CPU: Apple M1: 
              speed         user         nice          sys         idle          irq
       #1    24 MHz     465132 s          0 s     421564 s    4144316 s          0 s
       #2    24 MHz     436116 s          0 s     379290 s    4215496 s          0 s
       #3    24 MHz     397728 s          0 s     335802 s    4297371 s          0 s
       #4    24 MHz     367252 s          0 s     303627 s    4360022 s          0 s
       #5    24 MHz     246025 s          0 s      75884 s    4708991 s          0 s
       #6    24 MHz     170248 s          0 s      44494 s    4816159 s          0 s
       #7    24 MHz     122563 s          0 s      27951 s    4880386 s          0 s
       #8    24 MHz      92572 s          0 s      14593 s    4923735 s          0 s
       
  Memory: 16.0 GB (385.90625 MB free)
  Uptime: 1.067975e6 sec
  Load Avg:  1.8017578125  1.5673828125  1.34814453125
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-12.0.0 (ORCJIT, cyclone)
```