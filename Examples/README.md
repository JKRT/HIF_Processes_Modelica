# Models produced by the translator

## SysWithoutProcess.jl
Does not use the new operator.
Solving this using RK-4 or Tsit will result in failure.

## SysWithProcess.jl
Uses the new operator, possible to simulate using explicit solvers.
It simulates the fastest using the Tsit5 solver.
