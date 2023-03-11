## Working Directory ##
#cd("G:\\JuMPEnergy\\Cycle_part") # 跳转目录
#cd(pwd()*"\\Cycle_part") # 跳转到当前目录
## Include Basic functions ##
include("Props.jl")
## Include module functions ##
include("Turbine.jl")
include("Compressor.jl")
include("APoint.jl")

## Boundary Initialize ##

p₀ = 101325 # 大气压力
T₀ = 288.15 # 大气温度

p_max_cycle = 22e6 # 压力上限
T_max_cycle = 1200 # 温度上限 K
W_max_cycle = 5000 # 循环流量上限 kg/s
p_min_cycle = 7.5e6 # 压力下限
T_min_cycle = 308 # 温度下限

p_max = max(p₀,p_max_cycle)
p_min = min(p₀,p_min_cycle)

T_max = max(T₀,T_max_cycle)
T_min = max(T₀,T_min_cycle)

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
h₀_col = PT2h_Air(p₀,T₀) # 大气焓
s₀_col = PT2s_Air(p₀,T₀) # 大气熵
## Initial Calculate ##
# reactor #
h_core_in = T2h_FLiBe(T_core_in) # 堆芯入口焓
h_core_out = T2h_FLiBe(T_core_out) # 堆芯出口焓
h_core_0 = T2h_FLiBe(T₀) # 熔盐在大气压力下的死态焓
s_core_in =T2s_FLiBe(T_core_in)  # 堆芯入口熵
s_core_out = T2s_FLiBe(T_core_out) # 堆芯出口熵
s_core_0 = T2s_FLiBe(T₀) # 熔盐在大气压力下的死态熵
W_core = QA/(h_core_out - h_core_in) # 堆芯流量
# Condenser #
Tair_in = T₀ # 冷却器空气入口温度   
Pair_out = p₀ # 冷却器空气出口压力
Pair_in = Pair_out + dP_Air # 冷却器空气入口压力

# 声明变量及其范围
n_process = 10 # 流股数量
n_com = 1 # 压缩机数量
n_tur = 1 # 透平数量
n_HE = 6 # 换热过程数量

## Loading JuMP and NLopt ##
using JuMP,Ipopt
# NLP solver
model = Model(Ipopt.Optimizer)
## Regists ##
register(model, :T2h_FLiBe, 1, T2h_FLiBe; autodiff = true) 
register(model, :PT2h_CO2, 2, PT2h_CO2; autodiff = true) 
register(model, :PT2h_Air, 2, PT2h_Air; autodiff = true) 
register(model, :PT2s_Air, 2, PT2s_Air; autodiff = true) 
register(model, :TUR_io, 4, TUR_io; autodiff = true)
register(model, :COM_io, 4, TUR_io; autodiff = true)
register(model, :COM_oi, 4, COM_oi; autodiff = true)
register(model, :exergy, 5, exergy; autodiff = true)
## Variables Initialize ##
@variables(
        model,
        begin 
            T_min <= T[1:n_process] <= T_max # 温度
            p_min <= p[1:n_process] <= p_max # 压力
            0 <= W[1:n_process] # 流量
            0 <= Q_com[1:n_com] # 透平功
            0 <= Q_tur[1:n_tur] # 压缩功
            Q_HE[1:n_HE] # 换热量
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
        p[outID] - p[inID] == dp # 压降计算
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
        print("Props ERROR!")
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
    @constraint(model,Q_HE[HotHEID] + Q_HE[ColdHEID] == 0)
    if GL == "="
        @constraint(model,min(T[HotinID] - T[ColdoutID],T[HotoutID] - T[ColdinID]) == dT) # 约束最小端差
    elseif GL == ">"
        @constraint(model,T[HotinID] - T[ColdoutID] >= dT) # 约束最小端差
        @constraint(model,T[HotoutID] - T[ColdinID] >= dT) # 约束最小端差
    elseif GL == "<"
        @constraint(model,T[HotinID] - T[ColdoutID] <= dT) # 约束最小端差
        @constraint(model,T[HotoutID] - T[ColdinID] <= dT) # 约束最小端差
    else
        print(">,< or = must be input in Heattransfer!") 
    end
end
### 建模操作 ###
## 压力点约束
Point_set(6,p_min_cycle,"p","=") # 循环最低压力
Point_set(7,p_max_cycle,"p","=") # 循环最高压力
Point_set(2,p₀,"p","=") # 热源热出口压力
Point_set(10,p₀,"p","=") # 环境压力
## 温度点约束
Point_set(6,T_min_cycle,"T","=") # 循环最低温度
Point_set(1,T_core_out,"T","=") # 热源热进口温度
Point_set(2,T_core_in,"T","=") # 热源热出口温度
## 流量点约束
Point_set(1,W_core,"W","=") # 热源热侧流量
## 换热过程约束
# 热源
HeatProcess(1,2,dP_Heater,1,"FLiBe") # 热源热侧压降
HeatProcess(8,3,dP_Heater,1,"CO2") # 热源冷侧压降
# 回热器
HeatProcess(4,5,dP_RUP_h,2,"CO2") # 回热器热侧压降
HeatProcess(7,8,dP_RUP_l,2,"CO2") # 回热器冷侧压降
# 冷却器
HeatProcess(5,6,dP_Con,3,"Air") # 冷却器热侧压降
HeatProcess(9,10,dP_Air,3,"Air") # 冷却器冷侧压降
