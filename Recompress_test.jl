## Working Directory ##
cd("G:\\WorkingPath\\Julia test codes\\Cycle_组件化") # 跳转目录
## Include Basic functions ##
include("Salts.jl")
include("Props.jl")
## Include module functions ##
include("Turbine.jl")
include("Compressor.jl")
include("APoint.jl")
using JuMP, Ipopt

UP_P, DN_P, DN_T = [2.22606e7, 8.74496e6, 306.617]
cal_type = 1 # 选择优化方式，1和其他参数不同
## Initialize ##
Workingfluid = "CO2"
Coolingfluid = "Air"
QA = 125e6 # 热功率
T_core_out = 973.15 # 堆芯出口温度 K
T_core_in = 923.15  # 堆芯入口温度 K
dT_RUP = 10.0 # 回热器端差
dT_Heater = 10.0 # 加热器端差
dT_Con = 10.0 # 冷却器端差
η_com = 0.83 # 压缩机绝热效率
η_tur = 0.87 # 透平等熵效率
dP_Heater = 0.35e6 # 加热器压降
dP_RUP_h = 0.1e6 # 回热器热侧压降
dP_RUP_l = 0.1e6 # 回热器冷侧压降
dP_Con = 0.1e6 # 冷却器热侧压降
dP_Air = 0.1e3 # 冷却器冷侧压降
P₀ = 101325.0 # 大气压力
T₀ = 288.15 # 大气温度
## Initial Calculate ##
# Condenser #
Tair_in = T₀ # 冷却器空气入口温度   
Pair_out = P₀ # 冷却器空气出口压力
Pair_in = Pair_out + dP_Air # 冷却器空气入口压力
# NLP solver
model = Model(Ipopt.Optimizer)
## Regists ##
register(model, :PT2h, 2, PT2h; autodiff = true)
register(model, :PT2h_col, 2, PT2h_col; autodiff = true)
register(model, :PT2s_col, 2, PT2h_col; autodiff = true)
register(model, :TUR_ioh, 4, TUR_ioh; autodiff = true)
register(model, :COM_ioh, 4, COM_ioh; autodiff = true)
## Variables Initialize ##
@variables(model, begin
    Tair_in <= Tair_out <= T_core_in
    0 <= Wair
    0 <= W
    0 <= Wc
    DN_P <= P[1:10] <= UP_P
    DN_T <= T[1:10] <= T_core_out
    Q_com1 >= 0
    Q_com2 >= 0
    Q_tur >= 0
    0 <= Yita <= 1
end)
## Model Selcet ##
if cal_type == 1 # Cal model 采用定最高、最低压力和最低温度的方法来计算
    @variables(model, begin
        Pmin == DN_P
        Pmax == UP_P
        Tmin == DN_T

    end)
elseif cal_type == 2 # Optim model 采用变最高、最低压力和最低温度的方法来计算
    @variables(model, begin
        7.4e6 <= Pmin <= 10e6
        12e6 <= Pmax <= 22e6
        304 <= Tmin <= 320
    end)
    set_start_value(Pmin, DN_P)
    set_start_value(Pmax, UP_P)
    set_start_value(Tmin, DN_T)
else
    println("Type of cal ERROR")
    return 0.0
end
## Inequations ##
@constraints(model, begin
    # Flow constraints
    Wc <= W
    # Temperature Fields #
    T[5] == Tmin
    #0 == T[7] - T[8]

    T[4] >= T[5]
    T[3] >= T[4]
    T[2] >= T[3]
    T[1] >= T[2]

    T[6] >= T[5]
    T[7] >= T[6]
    T[8] >= T[4]
    T[10] >= T[9]
    # Condenser #
    T[5] - Tair_in >= dT_Con
    T[4] - Tair_out >= dT_Con
    # Regenerator #
    T[4] - T[6] == dT_RUP
    T[3] - T[7] >= dT_RUP
    T[3] - T[9] == dT_RUP
    T[2] - T[10] >= dT_RUP
    # Reactor #
    T_core_out - T[1] == dT_Heater
    T_core_in - T[10] >= dT_Heater
end)
## LEquations ##
@constraints(model, begin
    # Pressure Fields #
    P[5] == Pmin
    P[6] == Pmax
    # Pressure Fields #
    0 == P[4] - P[5] - dP_Con
    0 == P[6] - P[7] - dP_RUP_l
    0 == P[9] - P[10] - dP_RUP_l
    0 == P[3] - P[4] - dP_RUP_h
    0 == P[2] - P[3] - dP_RUP_h
    0 == P[10] - P[1] - dP_Heater
    0 == P[7] - P[8]
    0 == P[7] - P[9]
    Yita >= 0.4761
end)
## NLEquations-expression
@NLexpressions(model, begin
    # Condenser #
    ex_Condenser, Wc * (PT2h(P[4], T[4]) - PT2h(P[5], T[5])) - Wair * (PT2h_col(Pair_out, Tair_out) - PT2h_col(Pair_in, Tair_in))
    # Compressor1 #
    ex_Q_com1, Q_com1 - Wc * (PT2h(P[6], T[6]) - PT2h(P[5], T[5]))
    # Compressor2 #
    ex_Q_com2, Q_com2 - (W - Wc) * (PT2h(P[8], T[8]) - PT2h(P[4], T[4]))
    # Regenerator1 #
    ex_Regenerator1, (PT2h(P[2], T[2]) - PT2h(P[3], T[3])) - ((PT2h(P[10], T[10]) - PT2h(P[9], T[9])))
    # Regenerator2 #
    ex_Regenerator2, W * (PT2h(P[3], T[3]) - PT2h(P[4], T[4])) - Wc * ((PT2h(P[7], T[7]) - PT2h(P[6], T[6])))
    # Merge #
    ex_Merge, W * PT2h(P[9], T[9]) - (Wc * PT2h(P[7], T[7]) + (W - Wc) * PT2h(P[8], T[8]))
    # Turbine #
    ex_Q_tur, Q_tur - W * (PT2h(P[1], T[1]) - PT2h(P[2], T[2]))
    # Reactor #
    ex_Reactor, QA - W * (PT2h(P[1], T[1]) - PT2h(P[10], T[10]))
    # Dependent Variable
    ex_Yita, Yita - (Q_tur - Q_com1 - Q_com2) / QA
end)
## NLEquations ##
@NLconstraints(model, begin
    ex_Condenser == 0
    ex_Q_com1 == 0
    ex_Q_com2 == 0
    ex_Regenerator1 == 0
    ex_Regenerator2 == 0
    ex_Merge == 0
    ex_Q_tur == 0
    ex_Reactor == 0
    ex_Yita == 0
    PT2h(P[6], T[6]) - COM_ioh(PT2h(P[5], T[5]), P[5], P[6], η_com) == 0
    PT2h(P[8], T[8]) - COM_ioh(PT2h(P[4], T[4]), P[4], P[8], η_com) == 0
    PT2h(P[2], T[2]) - TUR_ioh(PT2h(P[1], T[1]), P[1], P[2], η_tur) == 0
end)
## Object Initialize ##
set_start_value(Yita, 1.0)
#set_start_value(Wc, 250)
## Optimize!! ##
@NLobjective(model, Max, Yita)
JuMP.optimize!(model)
## Call back ##
Tair_out = value.(Tair_out)
Wair = value.(Wair)
T = value.(T)
Tmin = value.(Tmin)
Pmin = value.(Pmin)
Pmax = value.(Pmax)
P = value.(P)
W = value.(W)
Wc = value.(Wc)
Q_tur = value.(Q_tur)
Q_com1 = value.(Q_com1)
Q_com2 = value.(Q_com2)
Yita = value.(Yita)
## Variables Soluting ##
ERR_P = [
    P[4] - P[5] - dP_Con,
    P[7] - P[6] - dP_RUP_l,
    P[10] - P[9] - dP_RUP_l,
    P[3] - P[4] - dP_RUP_h,
    P[2] - P[3] - dP_RUP_h,
    P[10] - P[1] - dP_Heater,
    P[7] - P[8],
    P[7] - P[9]
]
ex_Condenser = value(ex_Condenser)
ex_Q_com1 = value(ex_Q_com1)
ex_Q_com2 = value(ex_Q_com2)
ex_Regenerator1 = value(ex_Regenerator1)
ex_Regenerator2 = value(ex_Regenerator2)
ex_Q_tur = value(ex_Q_tur)
ex_Merge = value(ex_Merge)
ex_Reactor = value(ex_Reactor)
ex_Yita = value(ex_Yita)
ERR = [
    PT2h(P[6], T[6]) - COM_ioh(PT2h(P[5], T[5]), P[5], P[6], η_com),
    PT2h(P[8], T[8]) - COM_ioh(PT2h(P[4], T[4]), P[4], P[8], η_com),
    PT2h(P[2], T[2]) - TUR_ioh(PT2h(P[1], T[1]), P[1], P[2], η_tur),
    ex_Condenser,
    ex_Q_com1,
    ex_Q_com2,
    ex_Regenerator1,
    ex_Regenerator2,
    ex_Merge,
    ex_Q_tur,
    ex_Reactor,
    ex_Yita
]
