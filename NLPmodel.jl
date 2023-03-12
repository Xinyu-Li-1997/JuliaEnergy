## Include Basic functions ##
include("Props.jl")
## Include module functions ##
include("Turbomachine.jl")
include("APoint.jl")
include("Input.jl")
## Loading JuMP and NLopt ##
using JuMP,Ipopt,MadNLP
# NLP solver
model = Model(Ipopt.Optimizer)
#model = Model(()->MadNLP.Optimizer(print_level=MadNLP.INFO, max_iter=100))
## Regists ##
register(model, :sum, 1, sum; autodiff = true)
register(model, :T2h_FLiBe, 1, T2h_FLiBe; autodiff = true) 
register(model, :PT2h_CO2, 2, PT2h_CO2; autodiff = true)
register(model, :PT2s_CO2, 2, PT2s_CO2; autodiff = true)
register(model, :Ps2h_CO2, 2, Ps2h_CO2; autodiff = true)
register(model, :PT2h_Air, 2, PT2h_Air; autodiff = true) 
register(model, :PT2s_Air, 2, PT2s_Air; autodiff = true)
register(model, :Ps2h_Air, 2, Ps2h_Air; autodiff = true) 
register(model, :Ps2T_Air, 2, Ps2T_Air; autodiff = true)
register(model, :exergy, 5, exergy; autodiff = true)
## Variables Initialize ##
@variables(
        model,
        begin 
            T_min <= T[1:n_process] <= T_max # 温度
            p_min <= p[1:n_process] <= p_max # 压力
            300 <= W[1:n_process] # 流量
            0 <= Q_com[1:n_com] # 透平功
            0 <= Q_tur[1:n_tur] # 压缩功
            Q_HE[1:n_HE] # 换热量
            0 <= Q_heating[1:n_heatsource] # 热源输热量
            0 <= Yita <= 1 # 温度
            0 <= Q_heatall
            0 <= Q_comall
            0 <= Q_turboall
        end
    ) 

# 点约束函数
function Point_set(ID::Int64,value,type::String,GL::String)
    # ID 流股编号
    # value 数值
    # 数值类型 p,T,W
    # 大小关系 >,<,=
    if GL == "="
        if type == "p" @constraint(model,p[ID] == value) # 压力约束等式
        elseif type == "T" @constraint(model,T[ID] == value) # 温度约束等式
        elseif type == "W" @constraint(model,W[ID] == value) # 流量约束等式
        else print("ERROR type of point!")
        end
    elseif GL == ">"
        if type == "p" @constraint(model,value <= p[ID]) # 压力大于输入压力
        elseif type == "T" @constraint(model,value <= T[ID]) # 温度大于输入温度
        elseif type == "W" @constraint(model,value <= W[ID]) # 流量大于输入流量
        else print("ERROR type of point!")
        end
    elseif GL == "<"
        if type == "p" @constraint(model,p[ID] <= value) # 压力小于输入压力
        elseif type == "T" @constraint(model,T[ID] <= value) # 温度小于输入温度
        elseif type == "W" @constraint(model,W[ID] <= value) # 流量小于输入流量
        else print("ERROR type of point!")
        end
    else print(">,< or = must be input in Point_set!")   
    end
end
# 过程约束函数
function Process_set(ID::Int64,value,type::String,GL::String)
    # ID 流股编号
    # value 数值
    # 数值类型 Heat,Powerin,Powerout
    # 大小关系 >,<,=
    if GL == "="
        if type == "Heat" @constraint(model,Q_HE[ID] == value) # 传热量约束等式
        elseif type == "Powerin" @constraint(model,Q_com[ID] == value) # 输入功约束等式
        elseif type == "Powerout" @constraint(model,Q_tur[ID] == value) # 输出功约束等式
        else print("ERROR type of process!")
        end
    elseif GL == ">"
        if type == "Heat" @constraint(model,value <= Q_HE[ID]) # 传热量大于输入传热量
        elseif type == "Powerin" @constraint(model,value <= Q_com[ID]) # 输入功大于输入值
        elseif type == "Powerout" @constraint(model,value <= Q_tur[ID]) # 输出功大于输入值
        else print("ERROR type of process!")
        end
    elseif GL == "<"
        if type == "Heat" @constraint(model,Q_HE[ID] <= value) # 传热量小于输入传热量
        elseif type == "Powerin" @constraint(model,Q_com[ID] <= value) # 输入功小于输入值
        elseif type == "Powerout" @constraint(model,Q_tur[ID] <= value) # 输出功小于输入值
        else print("ERROR type of process!")
        end
    else print(">,< or = must be input in Process_set!")   
    end
end
# 热输运过程描述
function HeatProcess(inID::Int64,outID::Int64,dp,HEID,fld::String)
    # inID 入口流股编号
    # outID 出口流股编号
    # dp 压降
    # HEID 热过程编号
    # 流体类型
    @constraints(model,begin
        p[inID] - p[outID] == dp # 压降计算
        W[outID] - W[inID] == 0 # 流量守恒
    end)
    # 不同工质下的能量守恒计算
    if fld == "CO2"
        @NLconstraint(model, W[outID]*PT2h_CO2(p[outID], T[outID]) - W[inID]*PT2h_CO2(p[inID], T[inID]) == Q_HE[HEID])
    elseif fld == "Air"
        @NLconstraint(model, W[outID]*PT2h_Air(p[outID], T[outID]) - W[inID]*PT2h_Air(p[inID], T[inID]) == Q_HE[HEID])
    elseif fld == "FLiBe"
        @NLconstraint(model, W[outID]*T2h_FLiBe(T[outID]) - W[inID]*T2h_FLiBe(T[inID]) == Q_HE[HEID])
    else
        print("Props ERROR in HeatProcess!")
    end
end 
# 传热过程描述_端差
function Transfer_dT(HotinID::Int64,HotoutID::Int64,ColdinID::Int64,ColdoutID::Int64,HotHEID::Int64,ColdHEID::Int64,dT,GL::String)
    # HotinID 热入口流股编号
    # HotoutID 热出口流股编号
    # ColdinID 冷入口流股编号
    # ColdoutID 热出口流股编号
    # HotHEID 热侧过程编号
    # ColdHEID 冷侧过程编号
    # dT 端差
    # GL 大小关系 >,<,=
    @constraints(model,begin
        Q_HE[HotHEID] + Q_HE[ColdHEID] == 0 # 传热量守恒
        #T[HotinID] >= T[HotoutID] # 温差约束
        #T[ColdinID] <= T[ColdoutID] # 温差约束
    end)
    if GL == "="
        @constraint(model,min((T[HotinID] - T[ColdoutID]),(T[HotoutID] - T[ColdinID])) == dT) # 约束最小端差
    elseif GL == ">"
        @constraint(model,T[HotinID] - T[ColdoutID] >= dT) # 约束最小端差
        @constraint(model,T[HotoutID] - T[ColdinID] >= dT) # 约束最小端差
    elseif GL == "<"
        @constraint(model,T[HotinID] - T[ColdoutID] <= dT) # 约束最小端差
        @constraint(model,T[HotinID] - T[ColdoutID] >= 0) # 约束最小端差
        @constraint(model,T[HotoutID] - T[ColdinID] <= dT) # 约束最小端差
        @constraint(model,T[HotoutID] - T[ColdinID] >= 0) # 约束最小端差
    else
        print(">,< or = must be input in Heattransfer!") 
    end
end
# 压缩过程描述
function Compressor(inID::Int64,outID::Int64,comID::Int64,Yita,fld::String)
    # inID 入口流股编号
    # outID 出口流股编号
    # comID 压缩机编号
    # Yita 绝热效率
    # fld 流体类型
    @constraints(model,begin
        W[outID] - W[inID] == 0 # 流量守恒
        T[inID] <= T[outID] # 温差约束
    end)
    if fld == "CO2"
        @NLconstraints(model,begin
        (Ps2h_CO2(p[outID],PT2s_CO2(p[inID],T[inID])) - PT2h_CO2(p[inID],T[inID])) == (PT2h_CO2(p[outID],T[outID]) - PT2h_CO2(p[inID],T[inID]))*Yita # 不可逆过程
        Q_com[comID] == W[outID]*PT2h_CO2(p[outID],T[outID])- W[inID]*PT2h_CO2(p[inID],T[inID]) # 压缩功
        end)
    elseif fld == "Air"
        @NLconstraints(model,begin
        (Ps2h_Air(p[outID],PT2s_Air(p[inID],T[inID])) - PT2h_Air(p[inID],T[inID])) == (PT2h_Air(p[outID],T[outID]) - PT2h_Air(p[inID],T[inID]))*Yita # 不可逆过程
        Q_com[comID] == W[outID]*PT2h_Air(p[outID],T[outID])- W[inID]*PT2h_Air(p[inID],T[inID]) # 压缩功
        end)
    else
        print("Props ERROR in Compressor!")
    end
end
# 膨胀过程描述
function Turbine(inID::Int64,outID::Int64,turID::Int64,Yita,fld::String)
    # inID 入口流股编号
    # outID 出口流股编号
    # comID 压缩机编号
    # Yita 绝热效率
    # fld 流体类型
    @constraints(model,begin
        W[outID] - W[inID] == 0 # 流量守恒
        T[outID] <= T[inID] # 温差约束
    end)
    if fld == "CO2"
        @NLconstraints(model,begin
        (PT2h_CO2(p[inID],T[inID]) - PT2h_CO2(p[outID],T[outID])) == (PT2h_CO2(p[inID],T[inID]) - Ps2h_CO2(p[outID],PT2s_CO2(p[inID],T[inID])))*Yita # 不可逆过程
        Q_tur[turID] == -(W[outID]*PT2h_CO2(p[outID],T[outID])- W[inID]*PT2h_CO2(p[inID],T[inID])) # 膨胀功
        end)
    elseif fld == "Air"
        @NLconstraints(model,begin
        (PT2h_Air(p[inID],T[inID]) - PT2h_Air(p[outID],T[outID])) == (PT2h_Air(p[inID],T[inID]) - Ps2h_Air(p[outID],PT2s_Air(p[inID],T[inID])))*Yita # 不可逆过程
        Q_tur[turID] == -(W[outID]*PT2h_Air(p[outID],T[outID])- W[inID]*PT2h_Air(p[inID],T[inID])) # 膨胀功
        end)
    else
        print("Props ERROR in Turbine!")
    end
end
function Heater(heatID::Int64)
    # 计算加热量
    @constraint(model,Q_heating[heatID] + Q_HE[heatID] == 0) # 标记换热量
end
function Yitacal()
    # 计算最大热效率
    @constraint(model, Q_heatall == sum(Q_heating)) # 总加热量
    @constraint(model, Q_comall == sum(Q_com)) # 总压缩功
    @constraint(model, Q_turboall == sum(Q_tur)) # 总透平功
    @NLconstraint(model,Yita*Q_heatall == (Q_turboall - Q_comall)) # 标记换热量
    @NLobjective(model, Max, Yita)
end
