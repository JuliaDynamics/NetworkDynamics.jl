# Benchmark Report for */var/folders/cx/04cwlxl9287_ndx30k6q3n1h0000gn/T/jl_2hJZTs/NetworkDynamics*

## Job Properties
* Time of benchmark: 13 Aug 2021 - 18:59
* Package commit: a8c594
* Julia commit: 1b93d5
* Julia command flags: `-J/Applications/Julia-1.6.app/Contents/Resources/julia/lib/julia/sys.dylib`
* Environment variables: `JULIA_NUM_THREADS => 4`

## Results
Below is a table of this job's results, obtained by running the benchmarks.
The values listed in the `ID` column have the structure `[parent_group, child_group, ..., key]`, and can be used to
index into the BaseBenchmarks suite to retrieve the corresponding benchmarks.
The percentages accompanying time and memory values in the below table are noise tolerances. The "true"
time/memory value for a given benchmark is expected to fall within this percentage of the reported value.
An empty cell means that the value was zero.

| ID                                                   | time            | GC time   | memory          | allocations |
|------------------------------------------------------|----------------:|----------:|----------------:|------------:|
| `["diffusion", "ode_edge", "assemble", "10"]`        | 117.584 μs (5%) |           |  71.69 KiB (1%) |        1108 |
| `["diffusion", "ode_edge", "assemble", "100"]`       |   8.476 ms (5%) |           |   5.46 MiB (1%) |       88402 |
| `["diffusion", "ode_edge", "assemble", "1000"]`      | 946.734 ms (5%) | 29.063 ms | 524.06 MiB (1%) |     9036472 |
| `["diffusion", "ode_edge", "call", "10"]`            | 584.000 ns (5%) |           |                 |             |
| `["diffusion", "ode_edge", "call", "100"]`           |  77.583 μs (5%) |           |                 |             |
| `["diffusion", "ode_edge", "call", "1000"]`          |   8.371 ms (5%) |           |                 |             |
| `["diffusion", "ode_edge", "call_mt", "10"]`         |  31.042 μs (5%) |           |   3.91 KiB (1%) |          41 |
| `["diffusion", "ode_edge", "call_mt", "100"]`        |  52.750 μs (5%) |           |   3.97 KiB (1%) |          43 |
| `["diffusion", "ode_edge", "call_mt", "1000"]`       |   2.237 ms (5%) |           |   3.97 KiB (1%) |          43 |
| `["diffusion", "static_edge", "assemble", "10"]`     | 115.708 μs (5%) |           |  68.97 KiB (1%) |        1107 |
| `["diffusion", "static_edge", "assemble", "100"]`    |   8.548 ms (5%) |           |   5.13 MiB (1%) |       88399 |
| `["diffusion", "static_edge", "assemble", "1000"]`   | 959.874 ms (5%) | 32.461 ms | 491.16 MiB (1%) |     9036463 |
| `["diffusion", "static_edge", "call", "10"]`         | 208.000 ns (5%) |           |                 |             |
| `["diffusion", "static_edge", "call", "100"]`        |  27.042 μs (5%) |           |                 |             |
| `["diffusion", "static_edge", "call", "1000"]`       |   2.832 ms (5%) |           |                 |             |
| `["diffusion", "static_edge", "call_mt", "10"]`      |  24.208 μs (5%) |           |   3.53 KiB (1%) |          42 |
| `["diffusion", "static_edge", "call_mt", "100"]`     |  60.541 μs (5%) |           |   3.53 KiB (1%) |          42 |
| `["diffusion", "static_edge", "call_mt", "1000"]`    |   1.215 ms (5%) |           |   3.62 KiB (1%) |          45 |
| `["kuramoto", "heterogeneous", "assemble", "100"]`   | 599.833 μs (5%) |           | 364.41 KiB (1%) |        5969 |
| `["kuramoto", "heterogeneous", "assemble", "1000"]`  |   5.993 ms (5%) |           |   3.48 MiB (1%) |       64567 |
| `["kuramoto", "heterogeneous", "assemble", "10000"]` |  65.743 ms (5%) |           |  35.62 MiB (1%) |      693269 |
| `["kuramoto", "heterogeneous", "call", "100"]`       |   3.000 μs (5%) |           |  576 bytes (1%) |           8 |
| `["kuramoto", "heterogeneous", "call", "1000"]`      |  32.500 μs (5%) |           |  576 bytes (1%) |           8 |
| `["kuramoto", "heterogeneous", "call", "10000"]`     | 365.750 μs (5%) |           |  576 bytes (1%) |           8 |
| `["kuramoto", "heterogeneous", "call_mt", "100"]`    |  42.209 μs (5%) |           |   5.98 KiB (1%) |          71 |
| `["kuramoto", "heterogeneous", "call_mt", "1000"]`   |  46.916 μs (5%) |           |   6.02 KiB (1%) |          72 |
| `["kuramoto", "heterogeneous", "call_mt", "10000"]`  | 174.459 μs (5%) |           |   6.02 KiB (1%) |          72 |
| `["kuramoto", "homogeneous", "assemble", "100"]`     | 638.291 μs (5%) |           | 398.78 KiB (1%) |        6251 |
| `["kuramoto", "homogeneous", "assemble", "1000"]`    |   6.462 ms (5%) |           |   3.84 MiB (1%) |       68050 |
| `["kuramoto", "homogeneous", "assemble", "10000"]`   |  72.091 ms (5%) |           |  39.32 MiB (1%) |      728238 |
| `["kuramoto", "homogeneous", "call", "100"]`         |   2.625 μs (5%) |           |   64 bytes (1%) |           2 |
| `["kuramoto", "homogeneous", "call", "1000"]`        |  27.667 μs (5%) |           |   64 bytes (1%) |           2 |
| `["kuramoto", "homogeneous", "call", "10000"]`       | 331.458 μs (5%) |           |   64 bytes (1%) |           2 |
| `["kuramoto", "homogeneous", "call_mt", "100"]`      |  43.292 μs (5%) |           |   3.59 KiB (1%) |          43 |
| `["kuramoto", "homogeneous", "call_mt", "1000"]`     |  69.417 μs (5%) |           |   3.62 KiB (1%) |          44 |
| `["kuramoto", "homogeneous", "call_mt", "10000"]`    | 156.708 μs (5%) |           |   3.66 KiB (1%) |          45 |

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
Julia Version 1.6.2
Commit 1b93d53fc4 (2021-07-14 15:36 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin18.7.0)
  uname: Darwin 20.5.0 Darwin Kernel Version 20.5.0: Sat May  8 05:10:31 PDT 2021; root:xnu-7195.121.3~9/RELEASE_ARM64_T8101 x86_64 i386
  CPU: Apple M1: 
              speed         user         nice          sys         idle          irq
       #1    24 MHz     473711 s          0 s     398004 s    2913349 s          0 s
       #2    24 MHz     448116 s          0 s     362173 s    2974709 s          0 s
       #3    24 MHz     410586 s          0 s     326939 s    3047473 s          0 s
       #4    24 MHz     378113 s          0 s     297482 s    3109403 s          0 s
       #5    24 MHz     262559 s          0 s      83783 s    3438655 s          0 s
       #6    24 MHz     202835 s          0 s      52957 s    3529205 s          0 s
       #7    24 MHz     156097 s          0 s      39146 s    3589753 s          0 s
       #8    24 MHz     117402 s          0 s      20905 s    3646692 s          0 s
       
  Memory: 16.0 GB (231.88671875 MB free)
  Uptime: 720891.0 sec
  Load Avg:  1.98095703125  1.8681640625  1.91650390625
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, westmere)
```