#cd("E:\\Julia test codes\\Cycle_new")
#include("Initial.jl")
#include("Props.jl")

function ΔT_RUP(Pᵢₕ, Tᵢₕ, Pₒₕ, Tₒₕ, Pᵢₗ, Tᵢₗ, Pₒₗ, Tₒₗ, n=10)
    # n is number of deltaT
    # 2_2 is workingfluid to workingfluid, 2_3 is workingfluid to Coolingfluid, and 1_2 is salt to working fluid
    hᵢₕ = PT2h(Pᵢₕ,Tᵢₕ) #Enthalpy inlet at high temperature
    hₒₕ = PT2h(Pₒₕ,Tₒₕ) #Enthalpy outlet at high temperature
    hᵢₗ = PT2h(Pᵢₗ,Tᵢₗ) #Enthalpy inlet at low temperature
    hₒₗ = PT2h(Pₒₗ,Tₒₗ) #Enthalpy outlet at low temperature

    Δhₕ = (hᵢₕ - hₒₕ)/(n) #Delta enthalpy in high temperature
    Δhₗ = (hₒₗ - hᵢₗ)/(n) #Delta enthalpy in low temperature
    ΔTᵢ = zeros(n+1) #Delta temperature field

    ### Fleid points: 1 = min; n = max
    hₕ = zeros(n+1) #High temperature field 
    hₕ[1] = hₒₕ
    hₕ[end] = hᵢₕ

    hₗ = zeros(n+1) #High temperature field 
    hₗ[1] = hᵢₗ
    hₗ[end] = hₒₗ

    Tₕ₋ₚₘₐₓ = zeros(n+1) #High temperature field with max pressure
    Tₕ₋ₚₘₐₓ[1] = Tₒₕ
    Tₕ₋ₚₘₐₓ[end] = Tᵢₕ

    Tₕ₋ₚₘᵢₙ = zeros(n+1) #High temperature field with min pressure
    Tₕ₋ₚₘᵢₙ[1] = Tₒₕ
    Tₕ₋ₚₘᵢₙ[end] = Tᵢₕ

    Tₗ₋ₚₘₐₓ = zeros(n+1) #Low temperature field with max pressure
    Tₗ₋ₚₘₐₓ[1] = Tᵢₗ
    Tₗ₋ₚₘₐₓ[end] = Tₒₗ

    Tₗ₋ₚₘᵢₙ = zeros(n+1) #Low temperature field with min pressure
    Tₗ₋ₚₘᵢₙ[1] = Tᵢₗ
    Tₗ₋ₚₘᵢₙ[end] = Tₒₗ

    ΔTᵢ[1] = Tₒₕ - Tᵢₗ #Initialize of deltaT
    ΔTᵢ[end] = Tᵢₕ - Tₒₗ
    
    for i = 2:n
        hₕ[i] = hₕ[i-1] + Δhₕ
        hₗ[i] = hₗ[i-1] + Δhₗ
        Tₕ₋ₚₘₐₓ[i] = Ph2T(Pᵢₕ,hₕ[i])
        Tₕ₋ₚₘᵢₙ[i] = Ph2T(Pₒₕ,hₕ[i])
        Tₗ₋ₚₘₐₓ[i] = Ph2T(Pᵢₗ,hₗ[i])
        Tₗ₋ₚₘᵢₙ[i] = Ph2T(Pₒₗ,hₗ[i])
        ΔTᵢ[i] = min(Tₕ₋ₚₘₐₓ[i],Tₕ₋ₚₘᵢₙ[i]) - max(Tₗ₋ₚₘₐₓ[i],Tₗ₋ₚₘᵢₙ[i])
    end
    ΔTₘᵢₙ = minimum(ΔTᵢ)
    return ΔTₘᵢₙ
end

#=
function ΔT_RUP(Pᵢₕ, Tᵢₕ, Pₒₕ, Tₒₕ, Pᵢₗ, Tᵢₗ, Pₒₗ, Tₒₗ, n=10)
    hᵢₕ = PT2h(Pᵢₕ, Tᵢₕ) # Enthalpy inlet at high temperature
    hₒₕ = PT2h(Pₒₕ, Tₒₕ) # Enthalpy outlet at high temperature
    hᵢₗ = PT2h(Pᵢₗ, Tᵢₗ) # Enthalpy inlet at low temperature
    hₒₗ = PT2h(Pₒₗ, Tₒₗ) # Enthalpy outlet at low temperature
    dh_h = hᵢₕ - hₒₕ
    dh_l = hₒₗ - hᵢₗ
    Ph = (Pᵢₕ + Pₒₕ) / 2
    Pl = (Pᵢₗ + Pₒₗ) / 2
    dt = zeros(n+1)
    for i = 1:n+1
        x = (i-1) / n
        dt[i] = Ph2T(Ph, x * dh_h + hₒₕ) - Ph2T(Pl, x * dh_l + hᵢₗ)
    end
    print(dt)
    DT = minimum(dt)
    return DT
end
=#
#=
    Pᵢₕ = 7600000.0
    Tᵢₕ = 837.8295388600052
    Pₒₕ = 7500000.0
    Tₒₕ = 380.59280967467015
    Pᵢₗ = 20000000
    Tᵢₗ = 372.59280967467015
    Pₒₗ = 19900000.0
    Tₒₗ = 757.762645492023
    ANS_RUP = ΔT_RUP(Pᵢₕ,Tᵢₕ,Pₒₕ,Tₒₕ,Pᵢₗ,Tᵢₗ,Pₒₗ,Tₒₗ) 
=#

#=
using JuMP
using NLopt,Ipopt

model = Model(NLopt.Optimizer)
set_optimizer_attribute(model, "algorithm", :LN_COBYLA)
register(model, :ΔT_RUP, 8, ΔT_RUP; autodiff = true)

@variable(model, 800 <= Tin_h <= 900)
@variable(model, 7e6 <= Pin_h <= 8e6)
#@variable(model, Pin == Pᵢ)
#@variable(model, Pout == Pₒ)
#@variable(model, η == η)
@NLobjective(model, Min, abs(ΔT_RUP(Pin_h, Tin_h, Pₒₕ, Tₒₕ, Pᵢₗ, Tᵢₗ, Pₒₗ, Tₒₗ)))
#set_start_value(T, 350)
#set_start_value(P, 7.5e6)
@time JuMP.optimize!(model)
print("\n")
@show value(Tin_h)
@show value(Pin_h)
@show objective_value(model)
Pin_h = value(Pin_h)
Tin_h = value(Tin_h)
ANS_1 = ΔT_RUP(Pin_h, Tin_h, Pₒₕ, Tₒₕ, Pᵢₗ, Tᵢₗ, Pₒₗ, Tₒₗ)
=#

function ΔT_CON(Pᵢₕ, Tᵢₕ, Pₒₕ, Tₒₕ, Pᵢₗ, Tᵢₗ, Pₒₗ, Tₒₗ, n=10)
     # n is number of deltaT
    # 2_2 is workingfluid to workingfluid, 2_3 is workingfluid to Coolingfluid, and 1_2 is salt to working fluid
    hᵢₕ = PT2h(Pᵢₕ,Tᵢₕ) #Enthalpy inlet at high temperature
    hₒₕ = PT2h(Pₒₕ,Tₒₕ) #Enthalpy outlet at high temperature
    hᵢₗ = PT2h_col(Pᵢₗ,Tᵢₗ) #Enthalpy inlet at low temperature
    hₒₗ = PT2h_col(Pₒₗ,Tₒₗ) #Enthalpy outlet at low temperature

    Δhₕ = (hᵢₕ - hₒₕ)/(n) #Delta enthalpy in high temperature
    Δhₗ = (hₒₗ - hᵢₗ)/(n) #Delta enthalpy in low temperature
    ΔTᵢ = zeros(n+1) #Delta temperature field

    ### Fleid points: 1 = min; n = max
    hₕ = zeros(n+1) #High temperature field 
    hₕ[1] = hₒₕ
    hₕ[end] = hᵢₕ

    hₗ = zeros(n+1) #High temperature field 
    hₗ[1] = hᵢₗ
    hₗ[end] = hₒₗ

    Tₕ₋ₚₘₐₓ = zeros(n+1) #High temperature field with max pressure
    Tₕ₋ₚₘₐₓ[1] = Tₒₕ
    Tₕ₋ₚₘₐₓ[end] = Tᵢₕ

    Tₕ₋ₚₘᵢₙ = zeros(n+1) #High temperature field with min pressure
    Tₕ₋ₚₘᵢₙ[1] = Tₒₕ
    Tₕ₋ₚₘᵢₙ[end] = Tᵢₕ

    Tₗ₋ₚₘₐₓ = zeros(n+1) #Low temperature field with max pressure
    Tₗ₋ₚₘₐₓ[1] = Tᵢₗ
    Tₗ₋ₚₘₐₓ[end] = Tₒₗ

    Tₗ₋ₚₘᵢₙ = zeros(n+1) #Low temperature field with min pressure
    Tₗ₋ₚₘᵢₙ[1] = Tᵢₗ
    Tₗ₋ₚₘᵢₙ[end] = Tₒₗ

    ΔTᵢ[1] = Tₒₕ - Tᵢₗ #Initialize of deltaT
    ΔTᵢ[end] = Tᵢₕ - Tₒₗ
    
    for i = 2:n
        hₕ[i] = hₕ[i-1] + Δhₕ
        hₗ[i] = hₗ[i-1] + Δhₗ
        Tₕ₋ₚₘₐₓ[i] = Ph2T(Pᵢₕ,hₕ[i])
        Tₕ₋ₚₘᵢₙ[i] = Ph2T(Pₒₕ,hₕ[i])
        Tₗ₋ₚₘₐₓ[i] = Ph2T_col(Pᵢₗ,hₗ[i])
        Tₗ₋ₚₘᵢₙ[i] = Ph2T_col(Pₒₗ,hₗ[i])
        ΔTᵢ[i] = min(Tₕ₋ₚₘₐₓ[i],Tₕ₋ₚₘᵢₙ[i]) - max(Tₗ₋ₚₘₐₓ[i],Tₗ₋ₚₘᵢₙ[i])
    end
    ΔTₘᵢₙ = minimum(ΔTᵢ)
    return ΔTₘᵢₙ
end

#=
    Pᵢₕ = 7500000.0
    Tᵢₕ = 380.59280982738767
    Pₒₕ = 7400000.0
    Tₒₕ = 305.99
    Pᵢₗ = 101225.0
    Tᵢₗ = 288.15
    Pₒₗ = 101325.0
    Tₒₗ = 323.33336313388935

    ANS_CON = ΔT_CON(Pᵢₕ,Tᵢₕ,Pₒₕ,Tₒₕ,Pᵢₗ,Tᵢₗ,Pₒₗ,Tₒₗ) 
=#