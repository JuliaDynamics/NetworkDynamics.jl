# Contributing to NetworkDynamics.jl

Thanks for your interest in NetworkDynamics.jl! Bug reports, questions, documentation fixes and
code are all welcome.

## Getting help

If you have a question about using the package, there are several places to ask:

- An [issue](https://github.com/JuliaDynamics/NetworkDynamics.jl/issues) on this repository.
  Questions are welcome there, not only bug reports.
- The [Julia Discourse](https://discourse.julialang.org/), category *Specific Domains > Modelling &
  Simulations*.
- The [Julia Slack](https://julialang.org/slack/), where you can message the maintainers directly.

Please have a look at the [documentation](https://juliadynamics.github.io/NetworkDynamics.jl/stable/)
first. The tutorials and examples cover many common use cases.

## Reporting bugs

Please open an [issue](https://github.com/JuliaDynamics/NetworkDynamics.jl/issues). A good bug
report contains

- a minimal example that reproduces the problem, ideally something we can copy into a fresh REPL,
- the full error message and stack trace, or a description of what you expected instead,
- the output of `versioninfo()` and `Pkg.status()`.

## Proposing changes

For anything larger than a small fix, please open an issue first and describe what you have in
mind. That way we can discuss the design before you spend time on the code.

Pull requests should

- do one thing, e.g. fix one bug or add one feature,
- come with tests for the new or fixed behavior,
- add a short entry to
  [`NEWS.md`](https://github.com/JuliaDynamics/NetworkDynamics.jl/blob/main/NEWS.md) if the
  change is visible to users,
- update the docstrings and documentation if the API changes.

Contributions written with the help of AI tools are welcome. You are responsible for the code
you submit, so please check and understand what the tool produced. Keep the changes local and
brief, and leave unrelated code alone.

## Running the tests

The package, the tests in `test/` and the documentation in `docs/` form a single Julia workspace,
so they share one resolved set of packages. The simplest way to run all tests is to start Julia in
the repository root with `julia --project=@.` and run
```julia-repl
pkg> test
```
The full test suite takes a while. The test files run in parallel on several worker processes, and
you can pass the names of individual test files to run only those:
```julia
using Pkg
Pkg.test("NetworkDynamics"; test_args=["callbacks_test"])
```
During development it is often easier to work in a REPL in the test environment and run test files
interactively. This needs [Revise.jl](https://github.com/timholy/Revise.jl) in your global
environment. Start Julia from the repository root with `julia --project=test` (on a fresh clone,
run `pkg> instantiate` once), load the shared test helpers and then include the file you are
working on:
```julia
using Revise
includet("test/testutils.jl") # utils needed by multiple tests
include("test/callbacks_test.jl")
```
The GPU tests only run if CUDA is available on your machine.

## Building the documentation

The documentation is built with Documenter.jl from the `docs/` environment. For working on the
docs locally, use the `localmake.jl` script:
```sh
julia docs/localmake.jl          # full build, serves the docs on localhost:8000
julia docs/localmake.jl --draft  # skips the examples and rerenders on every change
```
It keeps the Julia session alive between builds, so rebuilding after a change is fast. The script
needs Revise.jl and LiveServer.jl in your global environment.

Building all examples takes some time. Only the pages listed in `pages` in `docs/make.jl` are built,
so you can comment out the ones you are not working on to iterate faster.

## Maintenance and support

NetworkDynamics.jl is maintained by Hans Würfel at the Potsdam Institute for Climate Impact
Research (PIK). Frank Hellmann (PIK) leads the project and takes part in design decisions.
Issues and pull requests are answered on a best-effort basis, usually within a few days.
Development is partly funded through research projects, see the
[README](https://github.com/JuliaDynamics/NetworkDynamics.jl#funding) for details.
