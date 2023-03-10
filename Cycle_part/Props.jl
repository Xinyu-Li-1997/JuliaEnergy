include("CO2.jl")
include("Air.jl")

## Workingfluid ##
function PT2h(P,T,ftype = Workingfluid)
    if ftype == "CO2"
        return CO2_PT2h.PT_h(P,T)
    elseif ftype == "Air"
        return Air_PT2h.PT_h(P,T)
    else
        print("Props ERROR!")
    end
end

function PT2s(P,T,ftype = Workingfluid)
    if ftype == "CO2"
        return CO2_PT2s.PT_s(P,T)
    elseif ftype == "Air"
        return Air_PT2s.PT_s(P,T)
    else
        print("Props ERROR!")
    end
end

function Ph2T(P,T,ftype = Workingfluid)
    if ftype == "CO2"
        return CO2_Ph2T.Ph_T(P,T)
    elseif ftype == "Air"
        return Air_Ph2T.Ph_T(P,T)
    else
        print("Props ERROR!")
    end
end

function Ph2s(P,T,ftype = Workingfluid)
    if ftype == "CO2"
        return CO2_Ph2s.Ph_s(P,T)
    elseif ftype == "Air"
        return Air_Ph2s.Ph_s(P,T)
    else
        print("Props ERROR!")
    end
end

function Ps2h(P,T,ftype = Workingfluid)
    if ftype == "CO2"
        return CO2_Ps2h.Ps_h(P,T)
    elseif ftype == "Air"
        return Air_Ps2h.Ps_h(P,T)
    else
        print("Props ERROR!")
    end
end

function Ps2T(P,T,ftype = Workingfluid)
    if ftype == "CO2"
        return CO2_Ps2T.Ps_T(P,T)
    elseif ftype == "Air"
        return Air_Ps2T.Ps_T(P,T)
    else
        print("Props ERROR!")
    end
end


## Coolingfluid ##

function PT2h_col(P,T,ftype = Coolingfluid)
    if ftype == "CO2"
        return CO2_PT2h.PT_h(P,T)
    elseif ftype == "Air"
        return Air_PT2h.PT_h(P,T)
    else
        print("Props ERROR!")
    end
end

function PT2s_col(P,T,ftype = Coolingfluid)
    if ftype == "CO2"
        return CO2_PT2s.PT_s(P,T)
    elseif ftype == "Air"
        return Air_PT2s.PT_s(P,T)
    else
        print("Props ERROR!")
    end
end

function Ph2T_col(P,T,ftype = Coolingfluid)
    if ftype == "CO2"
        return CO2_Ph2T.Ph_T(P,T)
    elseif ftype == "Air"
        return Air_Ph2T.Ph_T(P,T)
    else
        print("Props ERROR!")
    end
end

function Ph2s_col(P,T,ftype = Coolingfluid)
    if ftype == "CO2"
        return CO2_Ph2s.Ph_s(P,T)
    elseif ftype == "Air"
        return Air_Ph2s.Ph_s(P,T)
    else
        print("Props ERROR!")
    end
end

function Ps2h_col(P,T,ftype = Coolingfluid)
    if ftype == "CO2"
        return CO2_Ps2h.Ps_h(P,T)
    elseif ftype == "Air"
        return Air_Ps2h.Ps_h(P,T)
    else
        print("Props ERROR!")
    end
end

function Ps2T_col(P,T,ftype = Coolingfluid)
    if ftype == "CO2"
        return CO2_Ps2T.Ps_T(P,T)
    elseif ftype == "Air"
        return Air_Ps2T.Ps_T(P,T)
    else
        print("Props ERROR!")
    end
end
#=
function FLiBeh(T)
    return FLiBe.h(T)
end
=#
exergy(h,h0,s,s0,T0) = (h-h0)-T0*(s-s0)