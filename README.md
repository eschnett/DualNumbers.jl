# DualNumbers.jl

Calculate partial derivatives via automatic differentiation.

This package is very similar to
[ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl) or
[Zygote.jl](https://github.com/FluxML/Zygote.jl). The difference is
that these two packages calculate the gradient of a function, i.e. its
derivative with respect to all the function's arguments. This package
here only evaluates a single partial derivative, which can be more
efficient in certain circumstances. If you are in doubt which package
to use, look at ForwardDiff or Zygote.

* [![Documenter](https://img.shields.io/badge/docs-dev-blue.svg)](https://eschnett.github.io/DualNumbers.jl/dev)
* [![GitHub
  CI](https://github.com/eschnett/DualNumbers.jl/workflows/CI/badge.svg)](https://github.com/eschnett/DualNumbers.jl/actions)
* [![Codecov](https://codecov.io/gh/eschnett/DualNumbers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/eschnett/DualNumbers.jl)

## Examples

Derivative of a function with a single (scalar) argument:

```Julia
julia> f′(x) = derivative(f, x)
```

Derivative of a function with a multiple arguments, passed as array:

```Julia
julia> f′(xs) = derivative(f, xs, i)
```
