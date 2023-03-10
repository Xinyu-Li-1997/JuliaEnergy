## Working Directory ##
cd("G:\\JuMPEnergy\\Cycle_part") # 跳转目录
## Include Basic functions ##
include("Salts.jl")
include("Props.jl")
## Include module functions ##
include("Turbine.jl")
include("Compressor.jl")
include("APoint.jl")

Workingfluid = "CO2"
Coolingfluid = "Air"
## Boundary Initialize ##
UP_P = 22e6 # 压力上限
UP_T = 1200 # 温度上限 K
UP_W = 5000 # 循环流量上限 kg/s
UP_Wair = 5000 # 空气流量上限
UP_Yita = 1 # 热效率上限
DN_P = 7.5e6 # 压力下限
DN_T = 308 # 温度下限
DN_W = 10 # 循环流量下限  kg/s
DN_Wair = 10 # 空气流量下限
DN_Yita = 0 # 热效率下限
## Initialize ##
QA = 125e6 # 热功率
T_core_out = 973.15 # 堆芯出口温度 K
T_core_in = 923.15 # 堆芯入口温度 K
dT_RUP = 8 # 回热器端差
dT_Heater = 10 # 加热器端差
dT_Con = 10 # 冷却器端差
η_com = 0.85 # 压缩机绝热效率
η_tur = 0.9 # 透平等熵效率
dP_Heater = 0.2e6 # 加热器压降
dP_RUP_h = 0.1e6 # 回热器热侧压降
dP_RUP_l = 0.1e6 # 回热器冷侧压降
dP_Con = 0.1e6 # 冷却器热侧压降
dP_Air = 0.1e3 # 冷却器冷侧压降
P₀ = 101325 # 大气压力
T₀ = 288.15 # 大气温度
h₀_col = PT2h_col(P₀,T₀) # 大气焓
s₀_col = PT2s_col(P₀,T₀) # 大气熵
## Initial Calculate ##
# reactor #
h_core_in = FLiBe.h(T_core_in) # 堆芯入口焓
h_core_out = FLiBe.h(T_core_out) # 堆芯出口焓
h_core_0 = FLiBe.h(T₀) # 熔盐在大气压力下的死态焓
s_core_in = FLiBe.s(T_core_in)  # 堆芯入口熵
s_core_out = FLiBe.s(T_core_out) # 堆芯出口熵
s_core_0 = FLiBe.s(T₀) # 熔盐在大气压力下的死态熵
W_core = QA/(h_core_out - h_core_in) # 堆芯流量
# Condenser #
Tair_in = T₀ # 冷却器空气入口温度   
Pair_out = P₀ # 冷却器空气出口压力
Pair_in = Pair_out + dP_Air # 冷却器空气入口压力
## Loading JuMP and NLopt ##
using JuMP,Ipopt

# NLP solver
model = Model(Ipopt.Optimizer)
## Regists ##
register(model, :PT2h, 2, PT2h; autodiff = true) 
register(model, :PT2h_col, 2, PT2h_col; autodiff = true) 
register(model, :PT2s_col, 2, PT2h_col; autodiff = true) 
register(model, :TUR_io, 4, TUR_io; autodiff = true)
register(model, :COM_io, 4, TUR_io; autodiff = true)
register(model, :COM_oi, 4, COM_oi; autodiff = true)
register(model, :exergy, 5, exergy; autodiff = true)
## Variables Initialize ##
# 声明变量及其范围
@variables(model, begin 
    0 <= Q_com
    0 <= Q_tur
    Tair_in <= Tair_out <= T_core_in
    DN_Wair <= Wair <= UP_Wair
    DN_W <= W <= UP_W
    DN_P <= P[1:6] <= UP_P
    DN_T <= T[1:6] <= T_core_out
    DN_Yita <= Yita <= UP_Yita
end)
## Inequations ##
# 声明简单的不等式关系
@constraints(model,begin
    # Pressure Fields #
    P[4] == DN_P
    P[5] == UP_P
    #7.4e6 <= P[4] <= 8.5e6
    #20e6 <= P[5] <= 22e6
    # Temperature Fields #
    T[4] == DN_T
    #304 <= T[4] <= 320
    T[1] >= T[2]
    T[2] >= T[3]
    T[3] >= T[4]
    T[5] >= T[4]
    T[6] >= T[5]
    T[1] >= T[6]
    T[2] >= T[6]
    T[3] >= T[5]
    # Condenser #
    T[3] - Tair_in >= dT_Con
    T[4] - Tair_out >= dT_Con
    # Regenerator #
    T[2] - T[6] >= dT_RUP
    T[3] - T[5] >= dT_RUP
    # Reactor #
    T_core_out - T[1] >= dT_Heater
    T_core_in - T[6] >= dT_Heater
end)
## LEquations ##
# 声明压力线性等式关系
@constraints(model,begin
    # Pressure Fields #
    0 == P[3] - P[4] - dP_Con
    0 == P[5] - P[6] - dP_RUP_l
    0 == P[2] - P[3] - dP_RUP_h
    0 == P[6] - P[1] - dP_Heater
end)
## NLEquations ##
# 声明非线性方程组
@NLconstraints(model,begin
    # Condenser #
    0 == W*(PT2h(P[3],T[3]) - PT2h(P[4],T[4])) - Wair*(PT2h_col(Pair_out,Tair_out) - PT2h_col(Pair_in,Tair_in))
    # Compressor #
    0 == PT2h(P[4],T[4]) - COM_oi(T[5],P[5],P[4],η_com)
    0 == Q_com - W*(PT2h(P[5],T[5])- PT2h(P[4],T[4]))
    # Regenerator #
    0 == (PT2h(P[2],T[2]) - PT2h(P[3],T[3])) - ((PT2h(P[6],T[6]) - PT2h(P[5],T[5])))
    # Turbine #
    0 == (PT2h(P[2],T[2]) - TUR_io(T[1],P[1],P[2],η_tur))
    0 == Q_tur - W*(PT2h(P[1],T[1])- PT2h(P[2],T[2]))
    # Reactor #
    0 == QA - W*(PT2h(P[1],T[1])- PT2h(P[6],T[6]))
    # Yita #
    0 == Yita - (Q_tur - Q_com)/QA
end)

# exergy(h,h0,s,s0,T0) 㶲计算
DS = @NLexpression(model, 
    (W_core*s_core_in + Wair*PT2s_col(Pair_out,Tair_out)) - 
    (W_core*s_core_out + Wair*PT2s_col(Pair_in,Tair_in))
    )
DE = @NLexpression(model, 
    (W_core*exergy(h_core_out,h_core_0,s_core_out,s_core_0,T₀) + 
    Wair*exergy(PT2h_col(Pair_in,Tair_in),h₀_col,PT2s_col(Pair_in,Tair_in),s₀_col,T₀)) -
    (W_core*exergy(h_core_in,h_core_0,s_core_in,s_core_0,T₀) + 
    Wair*exergy(PT2h_col(Pair_out,Tair_out),h₀_col,PT2s_col(Pair_out,Tair_out),s₀_col,T₀))
)

## Optimize!! ##
@NLobjective(model,  Max, Yita) 
# 执行优化
JuMP.optimize!(model)

## Call back ## 
# 参数查询
Tair_out = value.(Tair_out)
Q_tur    = value.(Q_tur)
Q_com    = value.(Q_com)
Wair     = value.(Wair)
T        = value.(T)
P        = value.(P)
W        = value.(W)
Yita     = value.(Yita)
## Variables Soluting ##
Tmin = minimum(T)
Tmax = maximum(T)
Pmin = minimum(P)
Pmax = maximum(P)

# 校核夹点温差，可以不管这个
DT_RUP = ΔT_RUP(P[2], T[2], P[3], T[3], P[5], T[5], P[6], T[6])
DT_con = ΔT_CON(P[3], T[3], P[4], T[4], Pair_in, Tair_in, Pair_out, Tair_out)

# 输出压力平衡绝对误差
ERR_P = [
    P[3] - P[4] - dP_Con,
    P[5] - P[6] - dP_RUP_l,
    P[2] - P[3] - dP_RUP_h,
    P[6] - P[1] - dP_Heater
]

# 输出能量守恒绝对误差
ERR_EQ = [
    (W*(PT2h(P[3],T[3]) - PT2h(P[4],T[4])) - Wair*(PT2h_col(Pair_out,Tair_out) - PT2h_col(Pair_in,Tair_in)))/(W*(PT2h(P[3],T[3]) - PT2h(P[4],T[4]))),
    PT2h(P[5],T[5]) - COM_io(T[4],P[4],P[5],η_com),
    (Q_com - W*(PT2h(P[5],T[5])- PT2h(P[4],T[4])))/Q_com,
    (W*(PT2h(P[2],T[2]) - PT2h(P[3],T[3])) - W*((PT2h(P[6],T[6]) - PT2h(P[5],T[5]))))/(W*(PT2h(P[2],T[2]) - PT2h(P[3],T[3]))),
    PT2h(P[2],T[2]) - TUR_io(T[1],P[1],P[2],η_tur),
    (Q_tur - W*(PT2h(P[1],T[1])- PT2h(P[2],T[2])))/Q_tur,
    (QA - W*(PT2h(P[1],T[1])- PT2h(P[6],T[6])))/QA,
    Yita - (Q_tur - Q_com)/QA
]