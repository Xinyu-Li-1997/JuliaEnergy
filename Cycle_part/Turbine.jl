# Sequential calculation of turbine
function TUR_ioh(hᵢ,Pᵢ,Pₒ,η)
    sᵢ = Ph2s(Pᵢ,hᵢ)
    sₒ₀ = sᵢ # Hypothesis of isentropic expansion 
    hₒ₀ = Ps2h(Pₒ,sₒ₀) # Enthalpy with 100% η
    hₒ = hᵢ - (hᵢ - hₒ₀)*η # Enthalpy real J/kg
    return hₒ
end 
function TUR_io(Tᵢ,Pᵢ,Pₒ,η)
    hᵢ = PT2h(Pᵢ,Tᵢ)
    sᵢ = PT2s(Pᵢ,Tᵢ)
    sₒ₀ = sᵢ # Hypothesis of isentropic expansion 
    hₒ₀ = Ps2h(Pₒ,sₒ₀) # Enthalpy with 100% η
    hₒ = hᵢ - (hᵢ - hₒ₀)*η # Enthalpy real J/kg
    Tₒ = Ph2T(Pₒ,hₒ)
    #return Tₒ
    return hₒ
end 