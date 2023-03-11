
function COM_ioh(hᵢ,Pᵢ,Pₒ,η)
    sᵢ = Ph2s(Pᵢ,hᵢ)
    sₒ₀ = sᵢ # Hypothesis of isentropic compress 
    hₒ₀ = Ps2h(Pₒ,sₒ₀) # Enthalpy with 100% η
    hₒ = hᵢ + (hₒ₀ - hᵢ)/η # Enthalpy real J/kg
    return hₒ
end 

function COM_io(Tᵢ,Pᵢ,Pₒ,η)
    hᵢ = PT2h(Pᵢ,Tᵢ)
    sᵢ = PT2s(Pᵢ,Tᵢ)
    sₒ₀ = sᵢ # Hypothesis of isentropic compress 
    hₒ₀ = Ps2h(Pₒ,sₒ₀) # Enthalpy with 100% η
    hₒ = hᵢ + (hₒ₀ - hᵢ)/η # Enthalpy real J/kg
    Tₒ = Ph2T(Pₒ,hₒ)
    #return Tₒ
    return hₒ
end 

function COM_oi(Tₒ,Pₒ,Pᵢ,η)
    hₒ = PT2h(Pₒ,Tₒ)
    sₒ = PT2s(Pₒ,Tₒ)
    sᵢ₀ = sₒ # Hypothesis of isentropic compress 
    hᵢ₀ = Ps2h(Pᵢ,sᵢ₀) # Enthalpy with 100% η
    #(hₒ - hᵢ₀)/(hₒ - hᵢ) = η
    #(hₒ - hᵢ₀)/η = (hₒ - hᵢ)
    hᵢ = hₒ - (hₒ - hᵢ₀)/η
    #return Tₒ
    return hᵢ
end 

