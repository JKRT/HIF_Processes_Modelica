# HIF_Processes_Modelica
Initial proof of concept of handling Zimmer theta operator in Modelica.

There are two instructions provided here, the first is the long install option.
These instructions, might require some familiarity with the Julia Language and its ecosystem.

To test models generated using the new code generator, without generating them on your own machine locally see the examples folder

## How to Install

### Using the new operator for the entire compiler
Note that this need to be done in order to generate files using this new operator.
However, since this is currently on an experimental branch it might be a bit cumbersome.

- Install OpenModelica.jl by following the instructions here:
  - https://github.com/JKRT/OM.jl
- Switch to the follwowing branch: https://github.com/JKRT/OMBackend.jl/tree/when-eq-improvements
  -It should be noted that this branch is currently under development and things might change, please email me if there are any issues.

## Additional dependencies
- DifferentialEquations.jl
- MetaGraphs.jl
- ModelingToolkit.jl
- Plots.jl
- CSV.jl
- DataFrames.jl
- JuliaFormatter.jl
- Graphs.jl
+ All assorted depencies described at https://github.com/JKRT/OM.jl

Once this is complete the file *hifProcessTest.jl* may be included.

To test the various models decribed below under *models*
See

```julia
function dynamicTest(;modelName = "TestThetaMethod.NonlinearCircuit.Test.Circuit1Dynamic")
  global file = "./Models/TestThetaMethod2.mo"
  simulateDynamic(;solver = :(FBDF()))
  sysInfo = experimentWithTeta(;modelName = modelName)
  global LATEST_SYS_INFO = sysInfo
  testRun(sysInfo; outputFileName = outputFileNameHelper(modelName))
#  simulate()
end

"""
  Runs codegen were theta is used in a component
"""
function dynamicTestTheta1(;
                          modelName = "TestThetaMethod.NonlinearCircuit.Test.ThetaCircuit1Dynamic")
  global file = "./Models/TestThetaMethod2.mo"
  simulateDynamicTheta(;solver = :(FBDF()))
  sysInfo = experimentWithTeta(;modelName = modelName)
  global LATEST_SYS_INFO = sysInfo
  testRun(sysInfo; outputFileName = outputFileNameHelper(modelName))
end

"""
  Runs codegen were theta is used at the top level
"""
function dynamicTestTheta2(;
                           modelName = "TestThetaMethod.NonlinearCircuit.Test.ThetaCircuit2Dynamic")
  global file = "./Models/TestThetaMethod2.mo"
  simulateDynamicTheta2(;solver = :(FBDF()))
  sysInfo = experimentWithTeta(;modelName = modelName)
  global LATEST_SYS_INFO = sysInfo
  testRun(sysInfo; outputFileName = outputFileNameHelper(modelName))
end
```

It should be noted that running this code evaluates and generate code just-in-time in the currently Julia environment.
Hence, if they are run twice they will overwrite eachother.

If you use the system this way the file name of the resulting models will be:
```
lowered_sys<qualified-modelica-path>.jl
```

The files should be included in your local directory.
One such file is included as a starter:

```
loweredSysTestThetaMethod__NonlinearCircuit__Test__ThetaCircuit2Dynamic.jl
```

It should be the same as the files under the Examples directory.


### Testing the model using the new operator
Two example models have been generated, along with a Modelica program to simulate the same system using OpenModelica.
They can be run and examined directely without running the steps above.

To run the OpenModelica reference file. Navigate to the folder and execute

```
omc buildDynamic.mos
````

The two generated models are called

- sysWithProcess.jl
- sysWithoutProcess.jl

These can be executed using the following set of commands, each in its separate terminal since they overwrite each others symbols:

```
> julia-1.9
include("sysWithProcess.jl")
using Plots
sol = simulateRK4()
plot(sol, idxs = (1))
```


```
> julia-1.9
include("sysWithoutProcess.jl")
using Plots
sol = simulate()
plot(sol, idxs = (1))
```

To compare timing between the two variants, one using theta and one without.
```julia
using BenchmarkTools
include("sysWithProcess.jl")
@benchmark simulateTsit5()
include("sysWithoutProcess.jl")
@benchmark simulate()
```

To experiment with each of these models using different solvers see the following functions:

```julia
function simulateRodas5(; tspan = (0.0, 1.0))
    prob = problemDef(; tspan = tspan)
    sol = solve(prob, Rodas5(autodiff = false); reltol = 1e-6, abstol = 1e-6)
    return sol
end

function simulateTsit5(; tspan = (0.0, 1.0))
    prob = problemDef(; tspan = tspan)
    sol = solve(prob, Tsit5(); reltol = 1e-6, abstol = 1e-6);
    return sol
end

function simulateRK4(; tspan = (0.0, 1.0))
    prob = problemDef(; tspan = tspan)
    sol = solve(prob, RK4(); reltol = 1e-6, abstol = 1e-6)
    return sol
end
```

Please note that the system that is not using theta, that is the system defined by sysWithoutProcess.jl explicit solvers such as RK-4 and Tsit5 will not work.
However, the process that uses Θ will run successfully for both RK4 and for TSIT5.


## The Modelica example
The example model in Modelica is available under `Models`.
The model used to produce the two examples is TestThetaMethod2.mo.

The two relevant models are:

```modelica
model ThetaCircuit1Dynamic
  extends Circuit1Static;
TestThetaMethod.ThetaCapacitator Cnl(C (displayUnit = "F")= 1e-12) annotation(
    Placement(visible = true, transformation(origin = {60, -14}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  equation
    connect(Cnl.n, ground.p) annotation(
    Line(points = {{60, -24}, {60, -30}, {0, -30}}, color = {0, 0, 255}));
    connect(diode.n, Cnl.p) annotation(
    Line(points = {{60, 10}, {60, -4}}, color = {0, 0, 255}));
  annotation(
    Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})),
__OpenModelica_simulationFlags(lv = "LOG_DASSL,LOG_SIMULATION,LOG_STATS", s = "dassl", cpu = "()"),
__OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=infoXmlOperations ");
end ThetaCircuit1Dynamic;
```

and

```modelica
model ThetaCircuit2Dynamic
  parameter Real THETA = 1.0;
  extends Circuit1Static;
 Modelica.Electrical.Analog.Basic.Capacitor Cp(C (displayUnit = "F")= 1e-12 * THETA) annotation(
    Placement(visible = true, transformation(origin = {60, -14}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  equation
    connect(Cp.n, ground.p) annotation(
    Line(points = {{60, -24}, {60, -30}, {0, -30}}, color = {0, 0, 255}));
    connect(diode.n, Cp.p) annotation(
    Line(points = {{60, 10}, {60, -4}}, color = {0, 0, 255}));
  annotation(
    Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})),
__OpenModelica_simulationFlags(lv = "LOG_DASSL,LOG_SIMULATION,LOG_STATS", s = "dassl", cpu = "()"),
__OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=infoXmlOperations ");
end ThetaCircuit2Dynamic;
```

The first uses theta as a part of a component while the second use it in the parameter modification.


## Translation notes
The
Modelica.Electrical.Analog.Semiconductors.exlin is not very well documented.
The implementation of that function is:

```
function exlin "Exponential function linearly continued for x > Maxexp"
  extends Modelica.Icons.Function;
  input Real x;
  input Real Maxexp;
  output Real z;
algorithm
  z := if x > Maxexp then exp(Maxexp)*(1 + x - Maxexp) else exp(x);
end exlin;
```
This function is inlined in the frontend.

### The file ThetaExperiment.jl
The file `ThetaExperiment.jl` is a temporary file showing MTK structure before
generation of the more low level representation.

## Algorithm high level overview
- Index reduction performed as before
- The Θ affected part of the system is extracted using dataflow analysis
- Construction of code for continuation solvers, those are to be available in G_z

### G_z
G_z consists of each separate infinitely fast process.
It is constructed by means of extracting the relevant structure using static analysis.

### Substeps
1) We treat Θ as a normal parameter and perform index reduction
2) The system should be examined post index reduction (This is currently not done)
3) Θ identify, the variables belonging to the infinitely fast processes and extract this code for continuation solvers

### Error handling and reporting
Examine the resulting system post index reduction and see if Θ has been applied in the correct way

- Θ is not allowed as an argument to functions
- Θ can be expressed as a factor of Θ to the power of {0, -1}, ensuring the use of Θ as a time constant

Either a normal time deriviative Θ raised to the power of {0} or infinitely fast, power of {-1}

### Extracting code for continuation solvers from the index-0 system
This step describes how the relevant structure is extracted from the system.
We extract it by looking at the dependencies of the original system, extracting those parts affected by Θ.


## Open questions based on Zimmers paper with artificial states
How should the step-size control of the main simulation loop be controlled w.r.t to the
convergence speed of the continuation solver?
Currently this is left as future work


## Ongoing work
The previous bullet is currently ongoing work
