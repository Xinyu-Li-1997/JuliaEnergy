## Working Directory ##
cd("G:\\WorkingPath\\Julia test codes\\Cycle_组件化") # 跳转目录
## Include Basic functions ##
include("Salts.jl")
include("Props.jl")
## Include module functions ##
include("Turbine.jl")
include("Compressor.jl")
include("APoint.jl")
using JuMP, Ipopt  #,AmplNLWriter,Bonmin_jll,Couenne_jll,SHOT_jll
using BlackBoxOptim
UP_P, DN_P, DN_T = [2.0e7, 7.5e6, 305.0]
Workingfluid = "CO2"
Coolingfluid = "Air"
function Yita_thermal(UP_P, DN_P, DN_T, cal_type = 1)
    ## Initialize ##
    QA = 125e6
    T_core_out = 973.15 #max temp of _salt_ K
    T_core_in = 923.15 #min temp of _salt_ K
    dT_RUP = 10.0
    dT_Heater = 10.0
    dT_Con = 10.0
    η_com = 0.83
    η_tur = 0.87
    dP_Heater = 0.1e6
    dP_RUP_h = 0.1e6
    dP_RUP_l = 0.1e6
    dP_Con = 0.1e6
    dP_Air = 0.1e3
    P₀ = 101325.0
    T₀ = 288.15
    ## Initial Calculate ##
    # Condenser #
    Tair_in = T₀
    Pair_out = P₀
    Pair_in = Pair_out + dP_Air
    # NLP solver
    model = Model(Ipopt.Optimizer)
    set_silent(model) # 模型静默，不再打印结果，因为太长了
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
    if cal_type == 1 # Cal model
        @variables(model, begin
            Pmin == DN_P
            Pmax == UP_P
            Tmin == DN_T

        end)
    elseif cal_type == 2 # Optim model
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
    
    # 计算卡诺效率
    ηcannot = 1 - minimum(T) / maximum(T)
    println("ηcannot:", maximum(ηcannot))
    println("search:", [UP_P,DN_P,DN_T])

    # 热效率在0.1和卡诺效率之间，最大绝对误差小于10,存在最优解Yita(不存在则不会是浮点数)，没有求解失败
    if 0.1 < Yita < ηcannot && maximum(ERR) <= 10 && typeof(Yita) == Float64 && raw_status(model)!="Restoration_Failed"    
        println("maximum(ERR_EQ):", maximum(ERR))
        yita = Yita
        println("Yita = :", Yita, "\n")
        return -yita # 返回负的热效率，因为优化算法找最小值
    # 其他情况下，返回0。由于实际热效率的负数一定小于0，因此优化算法会在这里得到一个大的数，会继续计算。
    elseif typeof(Yita) == Nothing
        println("Yita is nothing!!")
        return 0.0
    else
        println("Other ERROR!")
        return 0.0
    end
end

# 封装优化函数，输入决策变量x=[UP_P, DN_P, DN_T]，输出热效率
function Optim_Yita(x)
    UP_P, DN_P, DN_T = x
    ANS = Yita_thermal(UP_P, DN_P, DN_T)
    if typeof(ANS) == Float64
        return ANS
    elseif typeof(ANS) == Nothing
        println("ERRORDATA:",x)
        return 0.1
    else
        return 0.1
    end    
end
# 设置优化上下限
UP_P_bound = (16.0e6, 24.0e6)
DN_P_bound = (7.4e6, 9.0e6)
DN_T_bound = (304.0, 320.0)
search = [UP_P_bound,DN_P_bound,DN_T_bound]
# 执行优化计算，MaxSteps = 50表示内迭代不超过50步
res2 = @time bboptimize(Optim_Yita; SearchRange = search, Method=:adaptive_de_rand_1_bin_radiuslimited, MaxSteps = 50)
# 获取最优决策变量
UP_P, DN_P, DN_T = best_candidate(res2)
# 获取最优解
OBj = best_fitness(res2)
