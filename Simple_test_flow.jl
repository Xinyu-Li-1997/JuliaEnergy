## Working Directory ##
#cd("G:\\JuMPEnergy\\Cycle_part") # 跳转目录
#cd(pwd()*"\\Cycle_part") # 跳转到当前目录
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

# 声明变量及其范围
n_process = 8 # 换热器数量
n_com = 1 # 压缩机数量
n_tur = 1 # 透平数量

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
@variables(
        model,
        begin 
            DN_T <= T_in[1:n_process] <= T_core_out # 进口温度
            DN_T <= T_out[1:n_process] <= T_core_out # 出口温度
            DN_P <= p_in[1:n_process] <= UP_P # 进口压力
            DN_P <= p_out[1:n_process] <= UP_P # 出口压力

            DN_W <= W[1:n_process] <= UP_W # 流量

            0 <= Q_com[1:n_com]
            0 <= Q_tur[1:n_tur]
        end
    ) 

function PressureLimit(ID,point,p) # 换热器压力约束函数
    if point == "in"
        @constraint(model,p_in[ID] == p)
    elseif point == "out"
        @constraint(model,p_out[ID] == p)
    else
        print("Only 'in','out' can be inputed")
    end
end

function TemperatureLimit(ID,point,T) # 换热器压力约束函数
    if point == "in"
        @constraint(model,T_in[ID] == T)
    elseif point == "out"
        @constraint(model,T_out[ID] == T)
    else
        print("Only 'in','out' can be inputed")
    end
end

function HeatTransfer()
    a = 1
    nothing
end

# 压力点约束
PressureLimit(5,"out",DN_P)
PressureLimit(7,"in",UP_P)
# 温度点约束
TemperatureLimit(5,"out",DN_T)
