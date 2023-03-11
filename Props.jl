include("CO2.jl")
include("Air.jl")
include("Salts.jl")

## CO2 ##
PT2h_CO2(P,T) = CO2_PT2h.PT_h(P,T)
PT2s_CO2(P,T) = CO2_PT2s.PT_s(P,T)
Ph2T_CO2(P,h) = CO2_Ph2T.Ph_T(P,h)
Ph2s_CO2(P,h) = CO2_Ph2s.Ph_s(P,h)
Ps2h_CO2(P,s) = CO2_Ps2h.Ps_h(P,s)
Ps2T_CO2(P,s) = CO2_Ps2T.Ps_T(P,s)

## Air ##
PT2h_Air(P,T) = Air_PT2h.PT_h(P,T)
PT2s_Air(P,T) = Air_PT2s.PT_s(P,T)
Ph2T_Air(P,h) = Air_Ph2T.Ph_T(P,h)
Ph2s_Air(P,h) = Air_Ph2s.Ph_s(P,h)
Ps2h_Air(P,s) = Air_Ps2h.Ps_h(P,s)
Ps2T_Air(P,s) = Air_Ps2T.Ps_T(P,s)

## FLiBe ##
T2h_FLiBe(T) = FLiBe.h(T)
T2s_FLiBe(T) = FLiBe.s(T)

exergy(h,h0,s,s0,T0) = (h-h0)-T0*(s-s0)