# Multi-Threading

Since version `0.3.0` multi-threading via the `Threads.@threads` macro is possible. This allows julia to integrate different nodes and edges in different threads, and can lead to significant performance gains on parallel architectures. To enable multi-threading call `network_dynamics` with the keyword argument `parallel=true`.

```julia
network_dynamics(vertices!, edges!, graph, parallel=true)
```

In order for this to take effect, multiple threads have to be available. This is achieved by setting the environment variable `JULIA_NUM_THREADS` **before** starting Julia.  To start Julia from a bash shell and with 4 threads use:
```
$ env JULIA_NUM_THREADS=4 julia
```

If you are using `Juno` for the `Atom` text editor `JULIA_NUM_THREADS` is set to the number of physical cores of your processor by default. This is also the number of threads we recommend to use.

!!! note
    The thread handling causes an overhead in the order of *20 μs* per call to the ODE function which might impair performance on small networks (<100 nodes) or on single core machines. In theses cases `network_dynamics` can be called without any additional arguments, since `parallel` defaults to `false`.



For more information on setting environment varibales see the [Julia documentation](https://docs.julialang.org/en/v1/manual/environment-variables/index.html#JULIA_NUM_THREADS-1).
