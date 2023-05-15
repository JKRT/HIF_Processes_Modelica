using NLsolve
using DifferentialEquations
#= Theta function for: "Differential(t)(Cp_v(t)) ~ Cp_i(t) / (1.0e-12THETA(t))" =#
function thetaProcess1(inAux)
    x = inAux[1]
    y = inAux[2]
    theta = [x[30]]
    θy = [theta[1]]
    function F(dy, y, aux, t)
        x = aux[1]
        y = aux[2]
        p = aux[3]
        sT = last(p)
        dy[1] = x[24] / (1.0e-12x[4])
    end
    prob = ODEProblem(F, θy, (0.0, 1.0), inAux)
    sol = DifferentialEquations.solve(prob, ImplicitEuler(autodiff = false))
    θy[1] = last(collect(Iterators.flatten(sol.u)))
end
#= The H function. dy and y are Vector{Float64} and aux::Vector{Vector{Float64}} =#
function h(dy, y, aux, t)
    local x = aux[1]
    aux[2] = y
    local p = aux[3]
    p[1] = t
    #= IFF process for Cp_v(t) =#
    thetaProcess1(aux)
    nlp = NonlinearProblem{true}(k, x, aux)
    nlSol = solve(nlp)
    x = nlSol.u
    #States#
    #STATE EQUATION Differential(t)(C1_v(t)) ~ 9999.99999999999C1_i(t)
    dy[1] = 9999.99999999999x[25]
    #STATE EQUATION Differential(t)(ifCond1(t)) ~ -0.0
    dy[2] = -0.0
    #STATE EQUATION Differential(t)(ifCond2(t)) ~ -0.0
    dy[3] = -0.0
    #STATE EQUATION Differential(t)(ifCond3(t)) ~ -0.0
    dy[4] = -0.0
    #STATE EQUATION Differential(t)(ifCond4(t)) ~ -0.0
    dy[5] = -0.0
    #STATE EQUATION Differential(t)(ifCond5(t)) ~ -0.0
    dy[6] = -0.0
    aux[1] = x
    aux[2] = y
    aux[3] = p
    p[length(p)] = t
end
#= The K function the z vector is the concatenation of the L and NL variables=#
function k(oz::Vector{Float64}, z::Vector{Float64}, aux::Vector{Vector{Float64}})
    y = aux[2]
    p = aux[3]
    t = last(p)
    #ALG L EQUATION 0 ~ ifelse(ifCond5(t) == true, 2.3538526683702e17(25.0diode_v(t) - 39.0), exp(25.0diode_v(t))) - ifEq_tmp53(t)
    oz[1] =
        ifelse(y[6] == true, 2.3538526683702e17(25.0z[28] - 39.0), exp(25.0z[28])) - z[1]
    #ALG L EQUATION 0 ~ ifelse(ifCond4(t) == true, 2.3538526683702e17(25.0D3_v(t) - 39.0), exp(25.0D3_v(t))) - ifEq_tmp47(t)
    oz[2] =
        ifelse(y[5] == true, 2.3538526683702e17(25.0z[29] - 39.0), exp(25.0z[29])) - z[2]
    #ALG L EQUATION 0 ~ ifelse(ifCond1(t) == true, 0.0, 1.0) - ifEq_tmp18(t)
    oz[3] = ifelse(y[2] == true, 0.0, 1.0) - z[3]
    #ALG L EQUATION 0 ~ 1.0 - THETA(t)
    oz[4] = 1.0 - z[4]
    #ALG L EQUATION 0 ~ Cp_v(t) + diode_v(t) - D3_n_v(t)
    oz[5] = z[30] + z[28] - z[5]
    #ALG L EQUATION 0 ~ 293.15 - diode_T_heatPort(t)
    oz[6] = 293.15 - z[6]
    #ALG L EQUATION 0 ~ D3_n_v(t) + D3_v(t) - D1_n_v(t)
    oz[7] = z[5] + z[29] - z[7]
    #ALG L EQUATION 0 ~ 293.15 - D3_T_heatPort(t)
    oz[8] = 293.15 - z[8]
    #ALG L EQUATION 0 ~ 1.0e-9(ifEq_tmp47(t) - 1.0) + 1.0e-8D3_v(t) - diode_i(t)
    oz[9] = 1.0e-9(z[2] - 1.0) + 1.0e-8z[29] - z[9]
    #ALG L EQUATION 0 ~ diode_i(t)*diode_v(t) - diode_LossPower(t)
    oz[10] = z[9] * z[28] - z[10]
    #ALG L EQUATION 0 ~ D3_v(t)*diode_i(t) - D3_LossPower(t)
    oz[11] = z[29] * z[9] - z[11]
    #ALG L EQUATION 0 ~ D1_n_v(t) - C1_v(t) - D2_v(t)
    oz[12] = z[7] - y[1] - z[12]
    #ALG L EQUATION 0 ~ ifelse(ifCond3(t) == true, 2.3538526683702e17(25.0D2_v(t) - 39.0), exp(25.0D2_v(t))) - ifEq_tmp41(t)
    oz[13] =
        ifelse(y[4] == true, 2.3538526683702e17(25.0z[12] - 39.0), exp(25.0z[12])) - z[13]
    #ALG L EQUATION 0 ~ ifelse(ifCond2(t) == true, 2.3538526683702e17(-39.0 - 25.0D2_v(t)), exp(-25.0D2_v(t))) - ifEq_tmp35(t)
    oz[14] =
        ifelse(y[3] == true, 2.3538526683702e17(-39.0 - 25.0z[12]), exp(-25.0z[12])) - z[14]
    #ALG L EQUATION 0 ~ 293.15 - D2_T_heatPort(t)
    oz[15] = 293.15 - z[15]
    #ALG L EQUATION 0 ~ 1.0e-9(ifEq_tmp41(t) - 1.0) + 1.0e-8D2_v(t) - D2_i(t)
    oz[16] = 1.0e-9(z[13] - 1.0) + 1.0e-8z[12] - z[16]
    #ALG L EQUATION 0 ~ D2_i(t)*D2_v(t) - D2_LossPower(t)
    oz[17] = z[16] * z[12] - z[17]
    #ALG L EQUATION 0 ~ 293.15 - D1_T_heatPort(t)
    oz[18] = 293.15 - z[18]
    #ALG L EQUATION 0 ~ 1.0e-9(ifEq_tmp35(t) - 1.0) - D1_i(t) - 1.0e-8D2_v(t)
    oz[19] = 1.0e-9(z[14] - 1.0) - z[19] - 1.0e-8z[12]
    #ALG L EQUATION 0 ~ -D1_LossPower(t) - D1_i(t)*D2_v(t)
    oz[20] = -z[20] - z[19] * z[12]
    #ALG L EQUATION 0 ~ 300.15 - R1_T_heatPort(t)
    oz[21] = 300.15 - z[21]
    #ALG L EQUATION 0 ~ Cp_v(t)*R1_i(t) - R1_LossPower(t)
    oz[22] = z[30] * z[27] - z[22]
    #ALG L EQUATION 0 ~ 1000.0 - R1_R_actual(t)
    oz[23] = 1000.0 - z[23]
    #ALG L EQUATION 0 ~ diode_i(t) - Cp_i(t) - R1_i(t)
    oz[24] = z[9] - z[24] - z[27]
    #ALG L EQUATION 0 ~ D2_i(t) + ifEq_tmp18(t) - C1_i(t) - D1_i(t)
    oz[25] = z[16] + z[3] - z[25] - z[19]
    #ALG L EQUATION 0 ~ C1_i(t) + Cp_i(t) + R1_i(t) - ground_p_i(t) - ifEq_tmp18(t)
    oz[26] = z[25] + z[24] + z[27] - z[26] - z[3]
    #= Solve ERROR with idx 0 =#
    ##= Solve R1_i(t) with idx 27 =#
    #ALG NL EQUATION 0 ~ Cp_v(t) - R1_R_actual(t)*R1_i(t)
    oz[27] = z[30] - z[23] * z[27]
    #= Solve diode_v(t) with idx 28 =#
    #ALG NL EQUATION 0 ~ diode_i(t) - 1.0e-9(ifEq_tmp53(t) - 1.0) - 1.0e-8diode_v(t)
    oz[28] = z[9] - 1.0e-9(z[1] - 1.0) - 1.0e-8z[28]
    #ALG NL EQUATION 0 ~ D2_i(t) + diode_i(t) - D1_i(t)
    oz[29] = z[16] + z[9] - z[19]
    #THETA EQUATION oz[30] = Cp_i(t) / (1.0e-12THETA(t))
    oz[30] = z[24] / (1.0e-12z[4])
end
function l(aux::Vector{Vector{Float64}})
    # Algebraics eq zero eqs #
    z = aux[1]
    y = aux[2]
    p = aux[3]
    t = last(p)
    #ALG EQUATION ifEq_tmp53(t) ~ ifelse(ifCond5(t) == true, 2.3538526683702e17(25.0diode_v(t) - 39.0), exp(25.0diode_v(t)))
    z[1] = ifelse(y[6] == true, 2.3538526683702e17(25.0z[28] - 39.0), exp(25.0z[28]))
    #ALG EQUATION ifEq_tmp47(t) ~ ifelse(ifCond4(t) == true, 2.3538526683702e17(25.0D3_v(t) - 39.0), exp(25.0D3_v(t)))
    z[2] = ifelse(y[5] == true, 2.3538526683702e17(25.0z[29] - 39.0), exp(25.0z[29]))
    #ALG EQUATION ifEq_tmp18(t) ~ ifelse(ifCond1(t) == true, 0.0, 1.0)
    z[3] = ifelse(y[2] == true, 0.0, 1.0)
    #ALG EQUATION THETA(t) ~ 1.0
    z[4] = 1.0
    #ALG EQUATION D3_n_v(t) ~ Cp_v(t) + diode_v(t)
    z[5] = z[30] + z[28]
    #ALG EQUATION diode_T_heatPort(t) ~ 293.15
    z[6] = 293.15
    #ALG EQUATION D1_n_v(t) ~ D3_n_v(t) + D3_v(t)
    z[7] = z[5] + z[29]
    #ALG EQUATION D3_T_heatPort(t) ~ 293.15
    z[8] = 293.15
    #ALG EQUATION diode_i(t) ~ 1.0e-9(ifEq_tmp47(t) - 1.0) + 1.0e-8D3_v(t)
    z[9] = 1.0e-9(z[2] - 1.0) + 1.0e-8z[29]
    #ALG EQUATION diode_LossPower(t) ~ diode_i(t)*diode_v(t)
    z[10] = z[9] * z[28]
    #ALG EQUATION D3_LossPower(t) ~ D3_v(t)*diode_i(t)
    z[11] = z[29] * z[9]
    #ALG EQUATION D2_v(t) ~ D1_n_v(t) - C1_v(t)
    z[12] = z[7] - y[1]
    #ALG EQUATION ifEq_tmp41(t) ~ ifelse(ifCond3(t) == true, 2.3538526683702e17(25.0D2_v(t) - 39.0), exp(25.0D2_v(t)))
    z[13] = ifelse(y[4] == true, 2.3538526683702e17(25.0z[12] - 39.0), exp(25.0z[12]))
    #ALG EQUATION ifEq_tmp35(t) ~ ifelse(ifCond2(t) == true, 2.3538526683702e17(-39.0 - 25.0D2_v(t)), exp(-25.0D2_v(t)))
    z[14] = ifelse(y[3] == true, 2.3538526683702e17(-39.0 - 25.0z[12]), exp(-25.0z[12]))
    #ALG EQUATION D2_T_heatPort(t) ~ 293.15
    z[15] = 293.15
    #ALG EQUATION D2_i(t) ~ 1.0e-9(ifEq_tmp41(t) - 1.0) + 1.0e-8D2_v(t)
    z[16] = 1.0e-9(z[13] - 1.0) + 1.0e-8z[12]
    #ALG EQUATION D2_LossPower(t) ~ D2_i(t)*D2_v(t)
    z[17] = z[16] * z[12]
    #ALG EQUATION D1_T_heatPort(t) ~ 293.15
    z[18] = 293.15
    #ALG EQUATION D1_i(t) ~ 1.0e-9(ifEq_tmp35(t) - 1.0) - 1.0e-8D2_v(t)
    z[19] = 1.0e-9(z[14] - 1.0) - 1.0e-8z[12]
    #ALG EQUATION D1_LossPower(t) ~ -D1_i(t)*D2_v(t)
    z[20] = -z[19] * z[12]
    #ALG EQUATION R1_T_heatPort(t) ~ 300.15
    z[21] = 300.15
    #ALG EQUATION R1_LossPower(t) ~ Cp_v(t)*R1_i(t)
    z[22] = z[30] * z[27]
    #ALG EQUATION R1_R_actual(t) ~ 1000.0
    z[23] = 1000.0
    #ALG EQUATION Cp_i(t) ~ diode_i(t) - R1_i(t)
    z[24] = z[9] - z[27]
    #ALG EQUATION C1_i(t) ~ D2_i(t) + ifEq_tmp18(t) - D1_i(t)
    z[25] = z[16] + z[3] - z[19]
    #ALG EQUATION ground_p_i(t) ~ C1_i(t) + Cp_i(t) + R1_i(t) - ifEq_tmp18(t)
    z[26] = z[25] + z[24] + z[27] - z[3]
end

#Problem Function
function problemDef(; tspan = (0.0, 1.0))
    p::Vector{Float64} = zeros(1)
    x::Vector{Float64} = zeros(26)
    y::Vector{Float64} = zeros(6)
    z::Vector{Float64} = zeros(3)
    theta::Vector{Float64} = zeros(1)
    aux::Vector{Vector{Float64}} = Vector{Vector{Float64}}[]
    #Assigning intial values
    y[2] = 0.0
    y[3] = 0.0

    y[4] = 0.0

    y[5] = 0.0

    y[6] = 0.0

    x[4] = 0.0

    x[26] = 0.0

    x[25] = 0.0

    z[1] = 0.0

    x[22] = 0.0

    x[21] = 288.15

    x[23] = 0.0

    x[19] = 0.0

    x[7] = 0.0

    x[20] = 0.0

    x[18] = 288.15

    x[12] = 0.0

    x[16] = 0.0

    x[17] = 0.0

    x[15] = 288.15

    z[3] = 0.0

    x[5] = 0.0

    x[11] = 0.0

    x[8] = 288.15

    z[2] = 0.0

    x[9] = 0.0

    x[10] = 0.0

    x[6] = 288.15

    x[24] = 0.0

    x[3] = 0.0

    x[14] = 0.0

    x[13] = 0.0

    x[2] = 0.0

    x[1] = 0.0

    y[1] = 0.0

    tspan = (0.0, 1.0)
    x = vcat(x, z, theta)
    push!(aux, x)
    push!(aux, y)
    push!(aux, p)
    l(aux)
    prob = ODEProblem(h, y, tspan, aux)
end
function simulate(; tspan = (0.0, 1.0))
    prob = problemDef(; tspan = tspan)
    sol = solve(prob, Rodas5(autodiff = false); abstol = 1e-06)
    return sol
end
function simulateTsit5(; tspan = (0.0, 1.0))
    prob = problemDef(; tspan = tspan)
    sol = solve(prob, Tsit5(); abstol = 1e-06)
    return sol
end
function simulateRK4(; tspan = (0.0, 1.0))
    prob = problemDef(; tspan = tspan)
    sol = solve(prob, RK4(); abstol = 1e-06)
    return sol
end
