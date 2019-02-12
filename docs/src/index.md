# NetworkDynamics.jl

## Overview

This package implements functions for defining and studying dynamics on networks.
The key construction is a callable struct compatible with the
DifferentialEquations.jl calling syntax.

```julia
nd = network_dynamics(nodes!, lines!, s_e, t_e)
nd(dx, x, p, t)
```

The key functions, or function arrays from which a network dynamics is
constructed are:

``$nodes![n](dx, x, [l]_s, [l]_t, p, t)$``
``$lines![e](dl, l, x_s, x_t, p, t)$``

Given edges $e$, and nodes $n$, as well as an orientation encoded by
the source function $s(e)$ and the target function $t(e)$
this implements the system of ODEs:

``$\frac{dx_n}{dt} = dx_n$``

``$\frac{dl_e}{dt} = dl_e$``

with $dx$ and $dl$ calculated by

``$[l]_s = [l_e \text{ if } s(e) = n]$``

``$[l]_t = [l_e \text{ if } t(e) = n]$``

``$nodes![n](dx_n, x_n, [l]_s, [l]_t, p_n, t)$``

``$lines![e](dl_e, l_e, x_{s(e)}, x_{t(e)}, p_e, t)$``

Something that relaxes to a diffusive network would for example be
implemented by


```julia
lines = (dl, l, x_1, x_2) -> dl .= 1000. * ((x_1 - x_2) - l)
nodes = (dx_n, x_n, l_s, l_t, p_n, t) -> dx_n .= f(x_n) - (sum(l_s) - sum(l_t))
```

## Static lines

For static line relations we similarly have:

```julia
sl_nd = static_lines_network_dynamics(nodes!, lines!, s_e, t_e)
sl_nd(dx, x, p, t)
```

With the convention for lines given by:

``$lines![e](l, x_s, x_t, p, t)$``

Given edges $e$, and nodes $n$, as well as an orientation encoded by
the source function $s(e)$ and the target function $t(e)$
this implements the system of ODEs:

``$\frac{dx_n}{dt} = dx_n$``

with $dx$ calculated by

``$lines![e](l_e, x_{s(e)}, x_{t(e)}, p_e, t)$``

``$[l]_s = [l_e \text{ if } s(e) = n]$``

``$[l]_t = [l_e \text{ if } t(e) = n]$``

``$nodes![n](dx_n, x_n, [l]_s, [l]_t, p_n, t)$``

A diffusive network would be implemented by

```julia
lines = (l, x_1, x_2) -> l .= x_1 - x_2
nodes = (dx_n, x_n, l_s, l_t, p_n, t) -> dx_n .= f(x_n) - (sum(l_s) - sum(l_t))
```

## Convenience functions for symbolic access to node variables

## Network DAEs
## Network SDEs
## Network DDEs

## API

```@autodocs
Modules = [NetworkDynamics]
```
