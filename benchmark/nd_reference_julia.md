# Benchmark Report for */var/folders/cx/04cwlxl9287_ndx30k6q3n1h0000gn/T/jl_ztre65/NetworkDynamics*

## Job Properties
* Time of benchmark: 28 Jul 2021 - 12:52
* Package commit: 7b72fc
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

| ID                                                   | time            | GC time | memory          | allocations |
|------------------------------------------------------|----------------:|--------:|----------------:|------------:|
| `["diffusion", "ode_edge", "assemble", "10"]`        |  86.500 μs (5%) |         |  51.94 KiB (1%) |         822 |
| `["diffusion", "ode_edge", "assemble", "100"]`       | 545.334 μs (5%) |         | 360.78 KiB (1%) |        5449 |
| `["diffusion", "ode_edge", "assemble", "1000"]`      |   5.242 ms (5%) |         |   3.45 MiB (1%) |       58801 |
| `["diffusion", "ode_edge", "assemble", "10000"]`     |  54.266 ms (5%) |         |  35.42 MiB (1%) |      637968 |
| `["diffusion", "ode_edge", "call", "10"]`            | 291.000 ns (5%) |         |                 |             |
| `["diffusion", "ode_edge", "call", "100"]`           |   2.833 μs (5%) |         |                 |             |
| `["diffusion", "ode_edge", "call", "1000"]`          |  33.833 μs (5%) |         |                 |             |
| `["diffusion", "ode_edge", "call", "10000"]`         | 378.291 μs (5%) |         |                 |             |
| `["diffusion", "ode_edge", "call_mt", "10"]`         |  13.250 μs (5%) |         |   3.94 KiB (1%) |          42 |
| `["diffusion", "ode_edge", "call_mt", "100"]`        |  14.291 μs (5%) |         |   3.94 KiB (1%) |          42 |
| `["diffusion", "ode_edge", "call_mt", "1000"]`       |  22.208 μs (5%) |         |   3.94 KiB (1%) |          42 |
| `["diffusion", "ode_edge", "call_mt", "10000"]`      | 111.166 μs (5%) |         |   3.94 KiB (1%) |          42 |
| `["diffusion", "static_edge", "assemble", "10"]`     |  84.375 μs (5%) |         |  50.00 KiB (1%) |         815 |
| `["diffusion", "static_edge", "assemble", "100"]`    | 525.500 μs (5%) |         | 345.31 KiB (1%) |        5442 |
| `["diffusion", "static_edge", "assemble", "1000"]`   |   5.137 ms (5%) |         |   3.30 MiB (1%) |       58785 |
| `["diffusion", "static_edge", "assemble", "10000"]`  |  53.232 ms (5%) |         |  33.95 MiB (1%) |      637954 |
| `["diffusion", "static_edge", "call", "10"]`         |  83.000 ns (5%) |         |                 |             |
| `["diffusion", "static_edge", "call", "100"]`        | 792.000 ns (5%) |         |                 |             |
| `["diffusion", "static_edge", "call", "1000"]`       |  10.250 μs (5%) |         |                 |             |
| `["diffusion", "static_edge", "call", "10000"]`      | 132.583 μs (5%) |         |                 |             |
| `["diffusion", "static_edge", "call_mt", "10"]`      |  13.083 μs (5%) |         |   4.05 KiB (1%) |          42 |
| `["diffusion", "static_edge", "call_mt", "100"]`     |  12.500 μs (5%) |         |   4.05 KiB (1%) |          42 |
| `["diffusion", "static_edge", "call_mt", "1000"]`    |  17.042 μs (5%) |         |   4.05 KiB (1%) |          42 |
| `["diffusion", "static_edge", "call_mt", "10000"]`   |  41.834 μs (5%) |         |   4.05 KiB (1%) |          42 |
| `["kuramoto", "heterogeneous", "assemble", "10"]`    |  91.084 μs (5%) |         |  50.98 KiB (1%) |         852 |
| `["kuramoto", "heterogeneous", "assemble", "100"]`   | 585.875 μs (5%) |         | 358.61 KiB (1%) |        5840 |
| `["kuramoto", "heterogeneous", "assemble", "1000"]`  |   5.671 ms (5%) |         |   3.44 MiB (1%) |       63531 |
| `["kuramoto", "heterogeneous", "assemble", "10000"]` |  58.693 ms (5%) |         |  35.14 MiB (1%) |      683224 |
| `["kuramoto", "heterogeneous", "call", "10"]`        |   1.000 μs (5%) |         |   1.09 KiB (1%) |          40 |
| `["kuramoto", "heterogeneous", "call", "100"]`       |   9.291 μs (5%) |         |  10.94 KiB (1%) |         400 |
| `["kuramoto", "heterogeneous", "call", "1000"]`      |  94.375 μs (5%) |         | 109.38 KiB (1%) |        4000 |
| `["kuramoto", "heterogeneous", "call", "10000"]`     | 966.125 μs (5%) |         |   1.07 MiB (1%) |       40000 |
| `["kuramoto", "heterogeneous", "call_mt", "10"]`     |  38.667 μs (5%) |         |   5.14 KiB (1%) |          82 |
| `["kuramoto", "heterogeneous", "call_mt", "100"]`    |  47.458 μs (5%) |         |  14.95 KiB (1%) |         441 |
| `["kuramoto", "heterogeneous", "call_mt", "1000"]`   |  86.375 μs (5%) |         | 113.42 KiB (1%) |        4042 |
| `["kuramoto", "heterogeneous", "call_mt", "10000"]`  | 331.000 μs (5%) |         |   1.07 MiB (1%) |       40043 |
| `["kuramoto", "homogeneous", "assemble", "10"]`      |  95.500 μs (5%) |         |  55.25 KiB (1%) |         898 |
| `["kuramoto", "homogeneous", "assemble", "100"]`     | 631.042 μs (5%) |         | 396.77 KiB (1%) |        6245 |
| `["kuramoto", "homogeneous", "assemble", "1000"]`    |   6.138 ms (5%) |         |   3.82 MiB (1%) |       68043 |
| `["kuramoto", "homogeneous", "assemble", "10000"]`   |  64.223 ms (5%) |         |  39.16 MiB (1%) |      728229 |
| `["kuramoto", "homogeneous", "call", "10"]`          | 250.000 ns (5%) |         |   64 bytes (1%) |           2 |
| `["kuramoto", "homogeneous", "call", "100"]`         |   2.291 μs (5%) |         |   64 bytes (1%) |           2 |
| `["kuramoto", "homogeneous", "call", "1000"]`        |  23.500 μs (5%) |         |   64 bytes (1%) |           2 |
| `["kuramoto", "homogeneous", "call", "10000"]`       | 285.583 μs (5%) |         |   64 bytes (1%) |           2 |
| `["kuramoto", "homogeneous", "call_mt", "10"]`       |  45.583 μs (5%) |         |   4.17 KiB (1%) |          44 |
| `["kuramoto", "homogeneous", "call_mt", "100"]`      |  31.125 μs (5%) |         |   4.14 KiB (1%) |          43 |
| `["kuramoto", "homogeneous", "call_mt", "1000"]`     |  46.334 μs (5%) |         |   4.17 KiB (1%) |          44 |
| `["kuramoto", "homogeneous", "call_mt", "10000"]`    | 150.875 μs (5%) |         |   4.20 KiB (1%) |          45 |

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
       #1    24 MHz     504589 s          0 s     454498 s    4520678 s          0 s
       #2    24 MHz     473307 s          0 s     409305 s    4597024 s          0 s
       #3    24 MHz     431714 s          0 s     362431 s    4685489 s          0 s
       #4    24 MHz     398390 s          0 s     327595 s    4753649 s          0 s
       #5    24 MHz     271190 s          0 s      82968 s    5125474 s          0 s
       #6    24 MHz     184804 s          0 s      48091 s    5246740 s          0 s
       #7    24 MHz     133032 s          0 s      29822 s    5316779 s          0 s
       #8    24 MHz     101016 s          0 s      15760 s    5362857 s          0 s
       
  Memory: 16.0 GB (58.72265625 MB free)
  Uptime: 1.223098e6 sec
  Load Avg:  2.45458984375  2.08642578125  1.73779296875
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, westmere)
```