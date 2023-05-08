#=
This program contains functions to simulate HIF processes
Part of this package implements a code generator
=#
using Revise
using Plots
using Test
using DataFrames
using CSV

import JuliaFormatter
import OM
import ModelingToolkit
import Graphs
import MetaGraphs
MTK = ModelingToolkit
#=
Change this to swap between the Francesco file and my own
=#
#global file::String = "./Models/NonLinearScaleable.mo"
#global file::String = "./Models/TestThetaMethod2.mo"
global file::String = "./Models/ThetaPaperExample.mo"
# Code for testinf HIF processes
function Static()
  res = OM.generateFlatModelica("TestThetaMethod.Test.Circuit1Dynamic", file; MSL = true)
  write("./FlatModels/PowerGridsSystem1.mo", res)
  solver = :(Rodas5())
end

function flattenDynamic()
  res = OM.generateFlatModelica("TestThetaMethod.NonlinearCircuit.Test.Circuit1Dynamic", file; MSL = true)
  write("./FlatModels/Dynamic.mo", res)
end

function simulate(modelName;startTime = 0.0, stopTime = 1.0, solver = :(FBDF()))
    sol = OM.simulate(modelName,
                      file;
                      startTime = startTime,
                      stopTime = stopTime,
                      MSL = true,
                      solver = solver);
  global LATEST_SOL = sol
  global LATEST_SYS = LATEST_SOL.prob.f.sys
  #global CURRENT_SUBSTITUTIONS = LATEST_SYS.substitutions.subs
  generateCSV(modelName * string("_", solver), sol)
end

function simulateNLP(;solver = Rodas5())
  simulate("NonLinearScalable"; startTime = 0.0, stopTime = 1.0, solver = solver)
end

function simulateWithRodas5()
  solver = :(Rodas5())
  simulateDynamic(solver = solver)
  simulateStatic(solver = solver)
end

function simulateWithRosenbrock23()
  solver = :(Rosenbrock23())
  simulateDynamic(solver = solver)
  simulateStatic(solver = solver)
end

function simulateDynamic(;solver = :(FBDF()))
  simulate("TestThetaMethod.NonlinearCircuit.Test.Circuit1Dynamic"; startTime = 0.0, stopTime = 1.0, solver = solver)
end

function simulateDynamicTheta(;solver = :(FBDF()))
  simulate("TestThetaMethod.NonlinearCircuit.Test.ThetaCircuit1Dynamic"; startTime = 0.0, stopTime = 1.0, solver = solver)
end

function simulateDynamicTheta2(;solver = :(FBDF()))
  simulate("TestThetaMethod.NonlinearCircuit.Test.ThetaCircuit2Dynamic"; startTime = 0.0, stopTime = 1.0, solver = solver)
end

function simulateStatic(;solver = :(FBDF()))
  simulate("TestThetaMethod.NonlinearCircuit.Test.Circuit1Static"; startTime = 0.0, stopTime = 1.0, solver = solver)
end

function simulateThetaStatic(;solver = :(FBDF()))
  simulate("TestThetaMethod.NonlinearCircuit.Test.ThetaCircuit1Static"; startTime = 0.0, stopTime = 1.0, solver = solver)
end

function simulateThetaPaperExample(;solver = :(FBDF()))
  simulate("ThetaPaperExample.M"; startTime = 0.0, stopTime = 1.0, solver = solver)
end

function writeDynamic()
  res = OM.generateFlatModelica("TestThetaMethod.NonlinearCircuit.Test.Circuit1Dynamic", file; MSL = true)
  write("./FlatModels/Dynamic.mo", res)
  OM.OMBackend.writeModelToFile("TestThetaMethod__NonlinearCircuit__Test__Circuit1Dynamic", "testTheta_Circuit1Dynamic.jl"; keepComments = false, keepBeginBlocks=false)
end


function writeDynamicTheta2()
  res = OM.generateFlatModelica("TestThetaMethod.NonlinearCircuit.Test.ThetaCircuit2Dynamic", file; MSL = true)
  write("./FlatModels/ThetaCircuit2Dynamic.mo", res)
  OM.OMBackend.writeModelToFile("TestThetaMethod__NonlinearCircuit__Test__ThetaCircuit2Dynamic", "testTheta_Circuit1Dynamic.jl"; keepComments = false, keepBeginBlocks=false)
end

function writeStatic()
  res = OM.generateFlatModelica("TestThetaMethod.NonlinearCircuit.Test.Circuit1Static", file; MSL = true)
  write("./FlatModels/Static.mo", res)
  OM.OMBackend.writeModelToFile("TestThetaMethod__NonlinearCircuit__Test__Circuit1Static", "testTheta_Circuit1Static.jl"; keepComments = false, keepBeginBlocks=false)
end

function writeStaticTheta()
  res = OM.generateFlatModelica("TestThetaMethod.NonlinearCircuit.Test.ThetaCircuit1Static", file; MSL = true)
  write("./FlatModels/ThetaStatic.mo", res)
  OM.OMBackend.writeModelToFile("TestThetaMethod__NonlinearCircuit__Test__ThetaCircuit1Static", "testTheta_ThetaCircuit1Static.jl"; keepComments = false, keepBeginBlocks=false)
end

function experimentWithTeta(;tspan  = (0.0, 1.0)
                            , modelName = "TestThetaMethod.NonlinearCircuit.Test.ThetaCircuit1Static")
  OM.translate(modelName, file; MSL = true)
  mModelName = replace(modelName, "." => "__")
  compiledModel = OM.OMBackend.getCompiledModel(mModelName)
  generateSimplifiedModel(compiledModel, modelName)
  println("Running testfile")
  @time runTestFile(tspan)
end

"""
  Generates a csv file from a plot
  s.t OMEdit can be used to check the results.
TODO: Make it so that the observed variables are also added to the resulting CSV file
"""
function generateCSV(modelName, sol)
  df1 = DataFrame(sol)
  rename!(df1, Dict(:timestamp => "time"))
  modelName = replace(modelName, "."=>"_")
  finalFileName = string(modelName,"_res.csv")
  CSV.write(finalFileName, df1)
end

"""
  Function to perform index reduction
"""
function systemSimplify(sys::ModelingToolkit.AbstractSystem, io = nothing; simplify = true, kwargs...)
  println("Calling local structural simplify")
  #sys = MTK.ode_order_lowering(sys)
  #sys = MTK.dae_index_lowering(sys)
  #sys = MTK.tearing(sys; simplify = simplify)
  #sys = MTK.tearing(sys; simplify = simplify)
  return MTK.structural_simplify(sys)
end

function generateSimplifiedModel(cModel, modelName)
  mModelName = replace(modelName, "." => "__")
  strippedModel = OM.OMBackend.CodeGeneration.stripBeginBlocks(cModel)
  @eval $strippedModel
  println(strippedModel)
  local modelConstructor = Meta.parse("OM.OMBackend.$(mModelName)Model()")
  (problem, callbacks, ivs, reducedSystem, tspan, pars, vars, irreductibles) = @eval $modelConstructor
  equationsBeforeOpt = ModelingToolkit.equations(reducedSystem)
  eqsStr = string(equationsBeforeOpt)
  eqsStr = replace(eqsStr, "Differential(t)" => "DER")
  eqsStr = replace(eqsStr, "(t)" => "")
  eqs = Meta.parse(eqsStr)
  initialValues = ivs
  tmpSystem = createSystem(reducedSystem, eqs, vars, pars, ivs, irreductibles)
  tmpSystem = OM.OMBackend.CodeGeneration.stripComments(tmpSystem)
  tmpSystem = OM.OMBackend.CodeGeneration.stripBeginBlocks(tmpSystem)
  modelStr::String = "$tmpSystem"
  writeFile("ThetaExperiment.jl", modelStr)
end

function writeFile(fileName, modelStr)
  local formattedResults = ""
  try
    formattedResults = JuliaFormatter.format_text(modelStr;
                                                  remove_extra_newlines = true,
                                                  always_use_return = false)
  catch e
    println("Exception during formatting!")
    println("String was:")
    println(modelStr)
    println("Exception was:")
    println(e)
    return
  end
  fdesc = open(fileName, "w")
  write(fdesc, formattedResults)
  close(fdesc)
end


function createSystem(reducedSystem, eqs, vars, pars, ivs, irreductibles = [])
    tmpSystem = quote
    begin
      import ModelingToolkit
      MTK = ModelingToolkit
    end

    struct SystemInfo
      system
      initialValues
      parameterMap
      equations
      variables
      parameters
      irreductibles
    end

    function system()
      MTK.@variables t
      DER = Differential(t)
      pars = MTK.@parameters($(Tuple(reducedSystem.ps)))
      vars = MTK.@variables($(Tuple(vars)))
      initialValues = $(Meta.parse(replace(string(ivs), "(t)" => "")))
      parameterMap = $(Meta.parse(replace(string(pars), "(t)" => "")))
      println("Marking irreductibles...")
      #= Make symbols available =#
      for var in vars
        varStr = replace(string(var), "(t)" => "")
        eval(:($(Symbol(varStr)) = $var))
      end
      for sym in $(irreductibles)
        eval(:($sym = SymbolicUtils.setmetadata($sym,
                                                ModelingToolkit.VariableIrreducible,
                                                true)))
      end

      println("Variables done")
      equations = $(eqs)
      system = MTK.ODESystem(
                equations,
                t,
                vars,
                pars;
                name = :($(Symbol("modelName"))),
      )
      return SystemInfo(system,
                        initialValues,
                        parameterMap,
                        equations,
                        vars,
                        pars,
                        $(if(isempty(irreductibles))
                            []
                          else
                            irreductibles
                          end))
    end
    end
end

function reverseLookupIdxToName(vHT::Dict{String, Int}, idx::Int)::String
  local vs = collect(values(vHT))
  local idxPos = findfirst(x->x==idx, vs)
  local ks = collect(keys(vHT))
  return ks[idxPos]
end

"""
  Returns the dependencies.
  That is the result of the equations before sorting.
  The index of the variables in ysHT are offeset by the #algebraic variables.
  The index of the variables in zsHT  are offeset by the #algebraic variables + #state variables
"""
function getEqDependencies(sys, xsHT, ysHT, zsHT)
  function getVariablesRecur(eq, xsHT, ysHT, zsHT)
    vStrs::Vector{String} = map(string, MTK.get_variables(eq))
    return vStrs
  end
  local algEqs = getTopSortEquations(sys)
  local stateNLEqs = MTK.equations(sys)
  local allEqs = vcat(algEqs, stateNLEqs)
  local eqsVarDeps = Vector{Int}[]
  for (idx,eq) in enumerate(allEqs)
    vStrs::Vector{String} = getVariablesRecur(eq, xsHT, ysHT, zsHT)
    deps = Int[]
    for vStr in vStrs
      if vStr in keys(xsHT)
        cIdx = xsHT[vStr]
        if idx != cIdx
          push!(deps, cIdx)
        end
      elseif vStr in keys(ysHT)
        cIdx = ysHT[vStr] + length(keys(xsHT))
        if idx != cIdx
          push!(deps, cIdx)
        end
      elseif vStr in keys(zsHT)
        push!(deps, zsHT[vStr] + length(keys(xsHT)) + length(keys(ysHT)))
      else
        throw("Unsupported system. Maybe parameters are not eliminated yet?")
      end
    end
    push!(eqsVarDeps, deps)
  end
  #= Postprocessing remove the self reference from afflicted variables =#

  return eqsVarDeps
end

"""
  Returns the sorted equations generated from MTK in symbolic form
"""
function getTopSortEquations(sys)
  local subs = ModelingToolkit.get_substitutions(sys)
  return deepcopy(subs.subs)
end

"""
  Maps the symbol to the index for the algebraic equations
"""
function mapSymbolToIdxAlgebraic(sys)::Dict{String, Int}
  local algEqs = getTopSortEquations(sys)
  local idx = Ref(0)
  local orderedSymIdxs = map(algEqs) do eq
    println(eq)
    idx.x = idx.x + 1
    string(eq.lhs), idx.x
  end
  ht = Dict(orderedSymIdxs)
  return ht
end


function createSymbolIntHT(sys)::Dict{Term{Real, Base.ImmutableDict{DataType, Any}}, Int}
  local algEqs = getTopSortEquations(sys)
  local idx = Ref(0)
  local orderedSymIdxs = map(algEqs) do eq
    idx.x = idx.x + 1
    eq.lhs, idx.x
  end
  ht = Dict(orderedSymIdxs)
  return ht
end

"""
  Same as mapSymbolToIdxAlgebraic
  but gives the parameter array
"""
function mapSymbolToIdxParameter(sys)
  pars = String[]
  for i in sys.ps
    push!(pars, string(i))
  end
  local idx = Ref(0)
  local symIdxPairs = map(pars) do p
    idx.x = idx.x + 1
    p,idx.x
  end
  ht = Dict(symIdxPairs)
  return ht
end

"""
  Map symbol to the index of state variables
"""
function mapSymbolToIdxState(sys)
  local stateVariables = states(sys)
  local idx = Ref(0)
  local symIdxPairs = map(stateVariables) do s
    idx.x = idx.x + 1
    string(s), idx.x
  end
  ht = Dict(symIdxPairs)
  return ht
end

"""
  Maps the index to a symbol
"""
function idxToSymbol(varSet)
  local idx = Ref(0)
  local orderedSymIdxs = map(varSet) do eq
    idx.x = idx.x + 1
    idx.x, string(eq.lhs)
  end
  ht = Dict(orderedSymIdxs)
  return ht
end

"""
  Maps symbol to variable name
"""
function symbolToName(sys)
  MTK.get_var_to_name(newProb.f.sys)
end

"""
  Replace the content of an equation s.t it maps to the correct index.
For NL Problems including the algebraic equations
"""
function replaceVariablesWithIdxNL(eq, algIdxMap, stateIdxMap, parIdxMap, nlAlgIdxMap, thetaHT)::String
  mtkVars = MTK.get_variables(eq)
  vars = map(mtkVars) do v
    string(v)
  end
  println(vars)
  local replacements = map(vars) do v
    println(v)
    if v in keys(nlAlgIdxMap)
      v, string("z[",nlAlgIdxMap[v] + maximum(values(algIdxMap)),"]")
    elseif v in keys(thetaHT)
      v, string("z[",thetaHT[v],"]")
    elseif v in keys(algIdxMap)
      v, string("z[",algIdxMap[v],"]")
    elseif v in keys(stateIdxMap)
      v, string("y[",stateIdxMap[v],"]")
    elseif v == "t" #= This case is left in but does not seem to occur=#
      "t", "t"
    else
      v, string("p[",parIdxMap[v],"]")
    end
  end
  local replacementMap = Dict(replacements)
  local eqAsStr = string(eq)
  for v in vars
    eqAsStr = replace(eqAsStr, v => replacementMap[v])
  end
  return eqAsStr
end

"""
Replace variable with the index and name of the corresponding array..
"""
function replaceVariablesWithIdx(eq,
                                 algIdxMap,
                                 stateIdxMap,
                                 parIdxMap,
                                 nlAlgIdxMap,
                                 thetaIdxMap = Dict{String, Int}())::String
  mtkVars = MTK.get_variables(eq)
  vars = map(mtkVars) do v
    string(v)
  end
  println(vars)
  local replacements = map(vars) do v
    println(v)
    if v in keys(thetaIdxMap)
      v, string("y[",thetaIdxMap[v],"]")
    elseif v in keys(nlAlgIdxMap)
      v, string("z[",nlAlgIdxMap[v],"]")
    elseif v in keys(algIdxMap)
      v, string("x[",algIdxMap[v],"]")
    elseif v in keys(stateIdxMap)
      v, string("y[",stateIdxMap[v],"]")
    elseif v == "t" #= This case is left in but does not seem to occur=#
      "t", "t"
    else
      v, string("p[",parIdxMap[v],"]")
    end
  end
  local replacementMap = Dict(replacements)
  local eqAsStr = string(eq)
  for v in vars
    eqAsStr = replace(eqAsStr, v => replacementMap[v])
  end
  return eqAsStr
end


function createInitialValuePairs(initialValues)
  map(initialValues) do iv
    (string(first(iv)), last(iv))
  end
end

function runTestFile(tspan; testFile = "ThetaExperiment.jl")
  include(testFile)
  sysInfo = system()
  reducedSystem = systemSimplify(sysInfo.system)
  @assign sysInfo.system = reducedSystem
  newProb = ODEProblem(reducedSystem, sysInfo.initialValues, tspan, sysInfo.parameterMap)
  global LATEST_PROB_FROM_CUSTOM_GEN = newProb
  #= Take FBDF =#
  sol = solve(newProb, Rodas5(); reltol = 1e-07, abstol = 1e-06)
  generateCSV("IFF_Custom", sol)
  return sysInfo
end

"""
This function should only be called from test run in its current position
"""
function createInitialValueAssignments(initialValuePairs, algIdxMap, stateIdxMap, parIdxMap, nlAlgIdxMap)
  assignments = map(initialValuePairs) do pair
    name,val = pair
    if name in keys(nlAlgIdxMap)
      string("z[",nlAlgIdxMap[name], "]", " = ", val)
    elseif name in keys(algIdxMap)
      string("x[",algIdxMap[name], "]", " = ", val)
    elseif name in keys(stateIdxMap)
      string("y[",stateIdxMap[name],"]", " = ", val)
    elseif name in keys(parIdxMap) #= This case is left in but does not seem to occur=#
      string("p[",parIdxMap[name],"]", " = ", val)
    else
    end
  end
end

"""
This function assumes that we have compiled the system in stage one
"""
function testRun(sysInfo; outputFileName = "loweredSys.jl")
  local sys = sysInfo.system
  local initialValues = sysInfo.initialValues
  local xsHT = mapSymbolToIdxAlgebraic(sys)
  local ysHT = mapSymbolToIdxState(sys)
  println(ysHT)
  psHt = mapSymbolToIdxParameter(sys)
  initialValuePairs = createInitialValuePairs(initialValues)
  algebraicEquations = getTopSortEquations(sys)
  stateEquations = MTK.equations(sys) #MTK.full_equations(sys)## #Get the complete equations for the nonlinear part
  nlEquations = filter(stateEquations) do eq
    if eq.lhs isa Int64 && eq.lhs == 0
      true
    else
      false
    end
  end
  differenceEquations = filter(stateEquations) do eq
    if eq.lhs isa Int64 && eq.lhs == 0
      false
    else
      true
    end
  end
  differenceVars = map(differenceEquations) do eq
    vars = Symbolics.get_variables(eq.lhs)
    @assert length(vars) == 1 "Variables should be 1. This needs to hold"
    string(first(vars))
  end
  #=
    Another HT for nl variables
  =#
  nlXsHt = Dict{String, Int}()
  @info "differenceVars" differenceVars
  nXsIndices = length(xsHT)
  offset = 1
  tearingXsIndiciesHT = Dict{Int, String}()
  #= Update the state variable maps =#
  for (k,v) in ysHT
    if k in differenceVars
      continue
    end
    #= Remove the special algebraic variable from the HT =#
    pop!(ysHT, k)
    #= We need to check if the index already exists in the HT =#
    newIdx = offset
#    xsHT[k] = newIdx
    nlXsHt[k] = newIdx
    offset += 1
  end
  global LATEST_ALG_INDICES = xsHT
  global LATEST_STATE_INDICES = ysHT
  global LATEST_NL_ALG_INDICES =  nlXsHt
  local thetaEquationsIdxs = Int[]
  local thetaIdxs = Int[]
  local thetaStates = String[]
  local thetaEquations = MTK.Equation[]
  local thetaHT = Dict{String, Int}()
  local oldDifferenceEquations = deepcopy(differenceEquations)
  #= Check if theta mode is relevant =#
  if "THETA(t)" in keys(xsHT)
    println("THETA MODE")
    (thetaSCCS, mg, thetaG) = extractThetaSCCS(sys, xsHT, ysHT, nlXsHt)
    println(thetaG)

    #= The theta induced state equations are to be removed from the differenceEquations =#
    for (i, thetaChild) in enumerate(last(thetaG))
      local thetaVarIdx = length(xsHT) + length(nlXsHt) + i
      println(thetaVarIdx)
      local stateIdx = thetaChild - nXsIndices
      println(stateIdx)
      local stateName = reverseLookupIdxToName(ysHT, stateIdx)
      push!(thetaEquationsIdxs, stateIdx)
      push!(thetaIdxs, thetaVarIdx)
      pop!(ysHT, stateName)
      #=
      Include the theta variable together with the algebraic variables.
      The index of theta is #Algebraic + #NL Algebraic
      =#
      thetaHT[stateName] = thetaVarIdx
      push!(thetaStates, stateName)
      #= Add constraint for the NL solver. That theta is 0=#
      local thetaEq = popat!(differenceEquations, stateIdx)
      thetaEq = MTK.Equation(0, thetaEq.rhs)
      push!(thetaEquations, thetaEq)
      println(thetaEquationsIdxs)
      println(xsHT)
      println(ysHT)
    end
    #=
    The indices of the original states need to be reordered.
    =#
    ysHT = adjustIndicesForTheta(thetaEquationsIdxs, ysHT)
    println(ysHT)
  else
    println("Codegen without THETA")
  end
  #=
    Code that is induced by theta is not to be solved anymore
  =#
  #= Create theta function =#
  local buffer = IOBuffer()
  println(buffer, "using NLsolve")
  println(buffer, "using DifferentialEquations")
  createThetaFunctions(buffer,
                       thetaEquationsIdxs,
                       thetaIdxs,
                       oldDifferenceEquations,
                       thetaStates,
                       xsHT,
                       ysHT,
                       psHt,
                       nlXsHt)
  #= End of preamble begin code generation =#
  #=
    The difference equations are the equations were the derivatives are solved.
  =#
  stateEquationStrs = map(differenceEquations) do eq
    println(eq)
    rwEq = replaceVariablesWithIdx(eq, xsHT, ysHT, psHt, nlXsHt, thetaHT)
    #= Some evil regexp code =#
    rwEq = replace(rwEq, r"Differential\(t\)\(y\[(\d+)\]\)" => s"dy[\1]")
    res = replace(rwEq, "~" => "=")
    res = "#STATE EQUATION $(eq)\n" * res
  end
  initialAssignments = filter(x->x!==nothing,
                              createInitialValueAssignments(initialValuePairs, xsHT, ysHT, psHt, nlXsHt))
  assigns = reduce(initialAssignments) do x,y
    string(x, "\n", y,"\n")
  end
  createHFunction(buffer, stateEquationStrs, thetaStates)
  createKFunction(buffer
                  ,sys
                  ,length(differenceEquations)
                  ,algebraicEquations
                  ,nlEquations
                  ,thetaEquations
                  ,xsHT
                  ,ysHT
                  ,psHt
                  ,nlXsHt
                  ,thetaHT)
  createLFunction(buffer
                  ,algebraicEquations::Vector{MTK.Equation}
                  ,xsHT
                  ,ysHT
                  ,psHt
                  ,nlXsHt
                  ,thetaHT)

  println(buffer, "\n\n#Problem Function")
  println(buffer, "function problemDef(;tspan = (0.0, 1.0))")
  #Parameters are assumed to have been removed
  println(buffer, "p::Vector{Float64} = zeros(1)")
  println(buffer, "x::Vector{Float64} = zeros($(length(keys(xsHT))))")
  println(buffer, "y::Vector{Float64} = zeros($(length(keys(ysHT))))")
  println(buffer, "z::Vector{Float64} = zeros($(length(keys(nlXsHt))))")
  println(buffer, "theta::Vector{Float64} = zeros($(length(keys(thetaHT))))")
  println(buffer, "aux::Vector{Vector{Float64}} = Vector{Vector{Float64}}[]")
  println(buffer, "#Assigning intial values")
  println(buffer, assigns)
  println(buffer, "tspan = (0.0,1.0)")
  println(buffer, "x = vcat(x, z, theta)")
  println(buffer, "push!(aux, x)")
  println(buffer, "push!(aux, y)")
  println(buffer, "push!(aux, p)")
  println(buffer, "l(aux)")
  println(buffer, "prob = ODEProblem(h, y, tspan, aux)")
  println(buffer, "end")
  println(buffer, "function simulate(;tspan=(0.0,1.0))")
  println(buffer, "prob = problemDef(;tspan=tspan)")
  println(buffer, "sol = solve(prob, Rodas5(autodiff=false); abstol = 1e-06)")
  println(buffer, "return sol")
  println(buffer, "end")
  #= Driver for TSIT5 =#
  println(buffer, "function simulateTsit5(;tspan=(0.0,1.0))")
  println(buffer, "prob = problemDef(;tspan=tspan)")
  println(buffer, "sol = solve(prob, Tsit5(); abstol = 1e-06)")
  println(buffer, "return sol")
  println(buffer, "end")
  #= Driver for RK4 (adaptive stepsize) =#
  println(buffer, "function simulateRK4(;tspan=(0.0,1.0))")
  println(buffer, "prob = problemDef(;tspan=tspan)")
  println(buffer, "sol = solve(prob, RK4(); abstol = 1e-06)")
  println(buffer, "return sol")
  println(buffer, "end")
  content = String(take!(buffer))
  writeFile(outputFileName, content)
end

"""
Create the L function.
Assigns additional initial values at time 0.
"""
function createLFunction(buffer
                         ,algebraicEquations::Vector{MTK.Equation}
                         ,xsHT
                         ,ysHT
                         ,psHt
                         ,nlXsHt
                         ,thetaHT)
  println(buffer, "function l(aux::Vector{Vector{Float64}})")
  println(buffer, "# Algebraics eq zero eqs #")
  println(buffer, "z = aux[1]")
  println(buffer, "y = aux[2]")
  println(buffer, "p = aux[3]")
  println(buffer, "t = last(p)")
  #= Mark the algebraic equations as residuals =#
  algebraicEquationsStrs = map(algebraicEquations) do eq
    rwEq = replaceVariablesWithIdxNL(eq, xsHT, ysHT, psHt, nlXsHt, thetaHT)
    replace(rwEq, "~" => "=")
    res = replace(rwEq, "~" => "=")
    res = "#ALG EQUATION $(eq)\n" * res
  end
  for aes in algebraicEquationsStrs
    println(buffer, aes)
  end
  println(buffer, "end")
end

"""
For NL Solve
"""
function createKFunction(buffer
                                ,sys
                                ,nDiffEqs
                                ,algebraicEquations::Vector{MTK.Equation}
                                ,nlEquations
                                ,thetaEquations::Vector{MTK.Equation}
                                ,xsHT
                                ,ysHT
                                ,psHt
                                ,nlXsHt
                                ,thetaHT)
  println("Calling K function NL Solve")
  local dxCounter::Int = nDiffEqs + length(thetaEquations) + 1
  @info "dxCounter" dxCounter
  local equationsDependencies = map(MTK.equation_dependencies(sys)) do eqDeps
    map(string, eqDeps)
  end
  @info "equationsDependencies" equationsDependencies
  algebraicEquationsStrs = map(algebraicEquations) do eq
    #= Mark the algebraic equations as residuals =#
    @assign eq.rhs = (eq.rhs - eq.lhs)
    @assign eq.lhs = 0
    rwEq = replaceVariablesWithIdxNL(eq, xsHT, ysHT, psHt, nlXsHt, thetaHT)
    replace(rwEq, "~" => "=")
    res = replace(rwEq, "~" => "=")
    res = "#ALG L EQUATION $(eq)\n" * res
  end
  algebraicStateEqStrs = map(nlEquations) do eq
    println(eq)
    rwEq = replaceVariablesWithIdxNL(eq, xsHT, ysHT, psHt, nlXsHt, thetaHT)
    #= Some evil regexp code =#
    rwEq = replace(rwEq, r"Differential\(t\)\(y\[(\d+)\]\)" => s"dy[\1]")
    res = replace(rwEq, "~" => "=")
    res = "#ALG NL EQUATION $(eq)\n" * res
  end
  thetaEqStrs = map(thetaEquations) do eq
    res = replaceVariablesWithIdxNL(eq, xsHT, ysHT, psHt, nlXsHt, thetaHT)
    res = "#THETA EQUATION $(eq)\n" * res
  end
  println(buffer, "#= The K function the z vector is the concatenation of the L and NL variables=#")
  println(buffer, "function k(oz::Vector{Float64}, z::Vector{Float64}, aux::Vector{Vector{Float64}})")
  println(buffer, "y = aux[2]")
  println(buffer, "p = aux[3]")
  println(buffer, "t = last(p)")
  local usedEqs::Set{Int} = Set{Int}()
  local unusedEqs::Set{String} = Set{String}()
  local nAlgebraicEquations = length(algebraicEquations)
  #= Write the algebraic equations =#
  local aLCounter = 1
  for aeqs in algebraicEquationsStrs
    println(buffer, replace(aeqs, "0 =" => "oz[$(aLCounter)] ="))
    push!(usedEqs, aLCounter)
    aLCounter += 1
  end
  @info "Before generation of the algebraic state equations" dxCounter
  for aes in algebraicStateEqStrs
    eqDeps = equationsDependencies[dxCounter]
    @info "Equation dependencies for $(aes) was: $(eqDeps)"
    local idx = 0
    local selected = "ERROR"
#    @assert length(eqDeps) == 0 "No dependencies for eq: $(aes)"
    for candidate in eqDeps
      @info "Equation dependencies" eqDeps
      @info "Looking for Candidate:" candidate
      @info "#algebraic:" nAlgebraicEquations
      if candidate in keys(nlXsHt)
        @info "Found candidate" nlXsHt[candidate]
        idx = nlXsHt[candidate] + nAlgebraicEquations
        selected = candidate
      end
    end
    @info "No dependencies found when matching. The tearing variables idx was:$(idx).\n Dependencies: $(eqDeps) \n $(aes)"
    tearingAlgIdx::Int = idx
    if tearingAlgIdx == 0
      println(buffer, "#= Solve $(selected) with idx $(tearingAlgIdx) =#")
      print(buffer, "#")
      push!(unusedEqs, replace(aes, "0 =" => "oz[$(tearingAlgIdx)] ="))
    else
      println(buffer, "#= Solve $(selected) with idx $(tearingAlgIdx) =#")
      push!(usedEqs, tearingAlgIdx)
      println(buffer, replace(aes, "0 =" => "oz[$(tearingAlgIdx)] ="))
    end
    dxCounter += 1
  end
  #=
  John 2023-03-18
  =#
#  fail()
  println(nlEquations)
  println(algebraicEquations)
  local sumAlgEqs = maximum(collect(values(nlXsHt))) + length(algebraicEquations)
  for idx in 1:sumAlgEqs
    if !(idx in usedEqs) && !(isempty(unusedEqs))
      push!(usedEqs, idx)
      res = pop!(unusedEqs)
      println(buffer, replace(res, "[0]"=>"[$(idx)]"))
    end
  end
  eqLen(x) = length(keys(x))
  l = eqLen(xsHT)  + eqLen(nlXsHt) + eqLen(thetaHT)
  #= Add the theta equations last =#
  for thetaEqStr in thetaEqStrs
    println(buffer, replace(thetaEqStr, "0 ~"=>"oz[$(l)] ="))
    l += 1
  end
  println(buffer,"end")
end

function createHFunction(buffer, stateEquationStrs, thetaStates = String[])
  println(buffer,
          "#= The H function. dy and y are Vector{Float64} and aux::Vector{Vector{Float64}} =#")
  println(buffer, "function h(dy, y, aux, t)")
  println(buffer, "local x = aux[1]")
  println(buffer, "aux[2] = y")
  println(buffer, "local p = aux[3]")
  println(buffer, "p[1] = t")
  for (i, ts) in enumerate(thetaStates)
    println(buffer, "#= IFF process for $(ts) =#")
    println(buffer, "thetaProcess$(i)(aux)")
  end
  println(buffer, "NLF! = (F, u) -> k(F, u, aux)")
  #=
  First the nonlinear problem need to be solved.
  The result is then fed into the set of state equations.
  =# #θ
  #= Solve NLP =#
  println(buffer, "x = (NLsolve.nlsolve(NLF!, x; autoscale = false, iterations = 1000, method=:newton)).zero")
  println(buffer,"#States#")
  for se in stateEquationStrs
    println(buffer, se)
  end
  println(buffer, "aux[1] = x")
  println(buffer, "aux[2] = y")
  println(buffer, "aux[3] = p")
  println(buffer, "p[length(p)] = t")
  println(buffer,"end")
end

"""
Function to create the theta function.
The thetaEquationIds vector contains the indices of the state equations that we generate infinitly fast suprocesses for.
The thetaIdxs is the index of the theta variable in the set of algebraic variables of the system.
"""
function createThetaFunctions(buffer,
                              thetaEquationsIdxs::Vector{Int},
                              thetaIdxs::Vector{Int},
                              differenceEquations::Vector{MTK.Equation},
                              thetaStates::Vector{String},
                              xsHT::Dict{String, Int},
                              ysHT::Dict{String, Int},
                              psHt::Dict{String, Int},
                              nlXsHt::Dict{String, Int})
  #= Generate nothing if we have no theta =#
  if isempty(thetaEquationsIdxs)
    return
  end
  for (i, eqIdx) in enumerate(thetaEquationsIdxs)
    println(buffer, "#= Theta function for: \"$(differenceEquations[thetaEquationsIdxs[i]])\" =#")
    println(buffer, "function thetaProcess$(i)(inAux)")
    println(buffer, "x = inAux[1]")
    println(buffer, "y = inAux[2]")
    println(buffer, "theta = [x[$(thetaIdxs[i])]]")
    println(buffer, "θy = [theta[1]]")
    println(buffer, "function F(dy, y, aux, t)")
    println(buffer, "x = aux[1]")
    println(buffer, "y = aux[2]")
    println(buffer, "p = aux[3]")
    println(buffer, "sT = last(p)")
    local thetaIdxMap = Dict{String, Int}()
    local thetaIdxMap[thetaStates[i]] = i
    local rwEq = replaceVariablesWithIdx(differenceEquations[eqIdx], xsHT, ysHT, psHt, nlXsHt, thetaIdxMap)
    rwEq = replace(rwEq, r"Differential\(t\)\(y\[(\d+)\]\) ~" => s"dy[\1] =")
    println(buffer, rwEq)
    println(buffer, "end")
    println(buffer, "prob = ODEProblem(F, θy, (0.0, 1.0), inAux)")
    println(buffer, "sol = DifferentialEquations.solve(prob, ImplicitEuler(autodiff = false))")
    println(buffer, "θy[1] = last(collect(Iterators.flatten(sol.u)))")
    println(buffer, "end")
  end
end

"""
  Extract the index of theta in the system.
  In the backend theta is assumed to have been renamed THETA everywhere.
  So it is really the same variable.
"""
function extractThetaIdx(xsHT::Dict{String, Int})
  local thetaStr = "THETA(t)"
  @assert thetaStr in keys(xsHT) "Error: THETA codegen without theta"
  local thetaIdx = xsHT[thetaStr]
  return thetaIdx
end

"""
Returns true if any of the variables in vars has an index that is present in one of the HTs.
"""
function containsIdx(vars
                     ,idx::Int
                     ,xsHT::Dict{String, Int}
                     ,ysHT::Dict{String, Int}
                     ,zsHT::Dict{String, Int})
  for v in vars
    if (v => idx) in xsHT
      return true
    elseif (v => idx) in ysHT
      return true
    elseif (v => idx) in zsHT
      return true
    else
      return false
    end
  end
end

"""
This function adjusts the indices w.r.t
the removed indices present in the thetaEquationsIdxs vector.
"""
function adjustIndicesForTheta(thetaEquationsIdxs::Vector{Int}, ysHT::Dict{String,Int})
  local sortedThetaEquationIdx = sort(thetaEquationsIdxs, rev = true)
  local yKeys = collect(keys(ysHT))
  local yVals = collect(values(ysHT))
  local newYsHt = Dict{String, Int}()
  for (i, removedStateIdx) in enumerate(sortedThetaEquationIdx)
    for j in 1:length(yVals)
      if yVals[j] > removedStateIdx
        yVals[j] = yVals[j] - 1
      end
    end
  end
  #= Create a new HT =#
  local yHT = Dict{String, Int}()
  for i in 1:length(yVals)
    push!(yHT, yKeys[i] => yVals[i])
  end
  return yHT
end

"""
  This functions extracts the equations affected by theta.
  The input is a MTK system, and hashtables from name to index for the algebraic variables,
  the state variables and the nonlinear variables.
The output is the graph along with the thetaSCC
"""
function extractThetaSCCS(sys::MTK.ODESystem,
                          xsHT::Dict{String, Int},
                          ysHT::Dict{String, Int},
                          zsHT::Dict{String, Int})
  #= Returns the variable idx of theta. =#
  local θidx = extractThetaIdx(xsHT)
  local algebraicEquations = getTopSortEquations(sys)
  local stateEquations = MTK.equations(sys)
  local θequations = Dict()
  local nAlgebraic = length(algebraicEquations)
  local nStates = length(keys(ysHT))
  local nNLs = length(zsHT)
  local nStateAndNLs = nStates + nNLs
  #=
  Get the variables for each equation and see which one is using theta.
  Providing a set of equations that is using theta.
  =#
  for (eqIdx, eq) in enumerate(vcat(algebraicEquations, stateEquations))
    vars = map(string, get_variables(eq))
    if containsIdx(vars, θidx, xsHT, ysHT, zsHT)
      push!(θequations, (eqIdx => eq))
    end
  end
  #= Get set of variables solved initially by the equation were Θ appears =#
  @info "The initial theta variables are:" θequations
  local dependentStates = Set{Int}()
  local dependentAlgs = Set{Int}()
  local dependentNLs = Set{Int}()
  for vIdx in keys(θequations)
    if vIdx > nAlgebraic
      local tmpVIdx = vIdx - (nStates + nNLs)
      if tmpVIdx - nNLs < nStates #State Variable
        push!(dependentStates, tmpVIdx)
      else
        push!(dependentNLs, tmpVIdx)
      end
    else
      push!(dependentAlgs, vIdx)
    end
  end
  @info "before returning the eq dependencies"
  local eqsDependencies = getEqDependencies(sys, xsHT, ysHT, zsHT)
  global EQ_DEPENDENCIES = eqsDependencies
  @info "dependentStates" dependentStates
  @info "dependentAlgs" dependentAlgs
  @info "dependentNLs" dependentNLs
  eqsDepsStr = dumpEquationDependencies(eqsDependencies, xsHT, ysHT, zsHT)
  println("Equation dependencies:\n" * eqsDepsStr)
  #= Create a meta graph, mg =#
  local mg = createGraph(eqsDependencies, xsHT, ysHT, zsHT)
  global MG = mg
  #=
  Extract the strongly connected components belonging to theta.
  Furthermore, we return the graph s.t we can generate code for the continuation solver
  =#
  local sccs::Vector{Vector{Int}} = Graphs.strongly_connected_components(mg)
  local thetaChildren::Vector{Int} = mg.graph.fadjlist[θidx]
  local thetaSCC::Vector{Vector{Int}} = Int[]
  #=
  This currently return all SCC.
  The intention is not to return only this but the theta induced subgraph

  =#
  for sc in sccs
    for tc in thetaChildren
      if tc in sc
        push!(thetaSCC, sc)
      end
    end
  end
  local flatThetaSCC = collect(Iterators.flatten(thetaSCC))
  #= Return the relevant scc components + the graph =#
  @assert length(thetaChildren) == 1 "Assuming one theta child for now: thetaChildren: $(thetaChildren)"
  return flatThetaSCC, mg, [θidx, thetaChildren]
end

"""
Maps a index to the corresponding ALG, S or NL idx.
"""
function idxIsState(idx, xsHT, ysHT, zsHT)
#  print("idxInState idx: $idx")
  nxs = length(keys(xsHT))
  nys = length(keys(ysHT))
  nzs = length(keys(zsHT))
#  println("nxs + nys = $(nxs + nys)")
  isstate = (idx > nxs) && (idx <=  nxs + nys)
  return isstate
end

function idxIsAlg(idx, xsHT, ysHT, zsHT)
  nxs = length(keys(xsHT))
  nys = length(keys(ysHT))
  nzs = length(keys(zsHT))
  return idx <= nxs
end

function idxIsNL(idx, xsHT, ysHT, zsHT)
  local alg = idxIsAlg(idx, xsHT, ysHT, zsHT)
  local state = idxIsState(idx, xsHT, ysHT, zsHT)
  return !(state || alg)
end

function getStateIdx(idx, xsHT, ysHT, zsHT)
  return idx - length(keys(xsHT))
end

function getAlgIdx(idx, xsHT, ysHT, zsHT)
  return idx - length(keys(zsHT))
end

function getNLIdx(idx, xsHT, ysHT, zsHT)
  return idx - length(keys(xsHT)) - length(keys(ysHT))
end


function dumpEquationDependencies(eqDependencies::Vector{Vector{Int}},
                                   xsHT,
                                   ysHT,
                                  zsHT)
  local nAlgs = length(keys(xsHT))
  local nStates = length(keys(ysHT))
  local nNLs = length(keys(zsHT))
  buffer = IOBuffer()
  println(buffer, "************************************************")
  println(buffer, "[")
  for (idx, eqDep) in enumerate(eqDependencies[1:length(keys(xsHT))])
    println(buffer, "\t[Algebraic Equation#$(idx)")
    for vDep in eqDep
      println(buffer, "\t\t " * string(vDep)
              * if idxIsState(vDep, xsHT, ysHT, zsHT)
                " State idx = $(getStateIdx(vDep, xsHT, ysHT, zsHT))"
              elseif idxIsNL(vDep, xsHT, ysHT, zsHT)
                " NL idx = $(getNLIdx(vDep, xsHT, ysHT, zsHT))"
              else
                " ALG idx = $(vDep)"
              end)
    end
    println(buffer, "\t]")
  end
  #= States =#
  for (idx, eqDep) in enumerate(eqDependencies[nAlgs + 1:nAlgs+nStates])
    println(buffer, "\t[State Equation#$(idx)")
    for vDep in eqDep
      println(buffer, "\t\t " * string(vDep)
              * if idxIsState(vDep, xsHT, ysHT, zsHT)
                " State idx = $(getStateIdx(vDep, xsHT, ysHT, zsHT))"
              elseif idxIsNL(vDep, xsHT, ysHT, zsHT)
                " NL idx = $(getNLIdx(vDep, xsHT, ysHT, zsHT))"
              else
                " ALG idx = $(vDep)"
              end)
    end
    println(buffer, "\t]")
  end
  #= NL Equations=#
  for (idx, eqDep) in enumerate(eqDependencies[nAlgs + nStates + 1:nAlgs + nStates + nNLs])
    println(buffer, "\t[NL Equation#$(idx)")
    for vDep in eqDep
      println(buffer, "\t\t " * string(vDep)
              * if idxIsState(vDep, xsHT, ysHT, zsHT)
                " State idx = $(getStateIdx(vDep, xsHT, ysHT, zsHT))"
              elseif idxIsNL(vDep, xsHT, ysHT, zsHT)
                " NL idx = $(getNLIdx(vDep, xsHT, ysHT, zsHT))"
              else
                " ALG idx = $(vDep)"
              end)
    end
    println(buffer, "\t]")
  end
  println(buffer, "]")
  println(buffer, "************************************************")
  return String(take!(buffer))
end

function testLoweredSys()
  include("loweredSys.jl")
  simulate()
end

function getVarIdxAttributeStr(idx::Int, xsHT, ysHT, zsHT)::String
  if idxIsState(idx, xsHT, ysHT, zsHT)
    return "STATE"
  elseif idxIsNL(idx, xsHT, ysHT, zsHT)
    return "NL"
  else
    return "ALG"
  end
end

function getVarNameAttribute(idx::Int, xsHT, ysHT, zsHT)::String
  if idxIsState(idx, xsHT, ysHT, zsHT)
    return collect(keys(ysHT))[idx - length(keys(xsHT))]
  elseif idxIsNL(idx, xsHT, ysHT, zsHT)
    return collect(keys(zsHT))[idx - length(keys(xsHT)) - length(keys(ysHT))]
  else
    return collect(keys(xsHT))[idx]
  end
end

"""
Create a graph using the Graphs library in Julia from the adjlst
TODO:
Write about the HTs
"""
function createGraph(adjLst::Vector{Vector{Int64}}, xsHT, ysHT, zsHT)
  local g = MetaGraphs.MetaDiGraph()
  for v in 1:length(adjLst)
    Graphs.add_vertex!(g)
  end
  for v in 1:length(adjLst)
    for e in adjLst[v]
      Graphs.add_edge!(g, v, e)
    end
  end
  #=Add meta information to the graph=#
  for v in 1:length(adjLst)
    MetaGraphs.set_prop!(g, v, :type, getVarIdxAttributeStr(v, xsHT, ysHT, zsHT))
    MetaGraphs.set_prop!(g, v, :name, getVarNameAttribute(v, xsHT, ysHT, zsHT))
  end
  return reverse(g)
end

function outputFileNameHelper(modelName)
  "loweredSys" * replace(modelName, "."=>"__") * ".jl"
end

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
  global file = "./Models/TestThetaMethod.mo"
  simulateDynamicTheta(;solver = :(FBDF()))
  sysInfo = experimentWithTeta(;modelName = modelName)
  global LATEST_SYS_INFO = sysInfo
  testRun(sysInfo; outputFileName = outputFileNameHelper(modelName))
#  simulate()
end

"""
  Runs codegen were theta is used at the top level
"""
function dynamicTestTheta2(;
                           modelName = "TestThetaMethod.NonlinearCircuit.Test.ThetaCircuit2Dynamic")
  global file = "./Models/TestThetaMethod.mo"
  simulateDynamicTheta2(;solver = :(FBDF()))
  sysInfo = experimentWithTeta(;modelName = modelName)
  global LATEST_SYS_INFO = sysInfo
  testRun(sysInfo; outputFileName = outputFileNameHelper(modelName))
#  simulate()
end

function nonLinearTest()
  modelName::String = "NonLinearScalable"
  global file = "./Models/NonLinearScaleable.mo"
  sysInfo = experimentWithTeta(;modelName = modelName)
  testRun(sysInfo, outputFileName = outputFileNameHelper(modelName))
#  simulate()
end

function dynamicTestThetaPaperExample(;
                                      modelName = "ThetaPaperExample.M")
  global file = "./Models/ThetaPaperExample.mo"
  simulateThetaPaperExample(;solver = :(FBDF()))
  sysInfo = experimentWithTeta(;modelName = modelName)
  global LATEST_SYS_INFO = sysInfo
  testRun(sysInfo, outputFileName = outputFileNameHelper(modelName))
#  simulate()
end
