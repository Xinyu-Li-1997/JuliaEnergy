include("NLPmodel.jl")
### 建模操作 ###
## 压力点约束
Point_set(6,p_min_cycle,"p","=") # 循环最低压力
Point_set(7,p_max_cycle,"p","=") # 循环最高压力
Point_set(2,dP_Heater,"p","=") # 热源热出口压力
Point_set(10,p₀,"p","=") # 环境压力
## 温度点约束
Point_set(6,T_min_cycle,"T","=") # 循环最低温度
Point_set(1,T_core_out,"T","=") # 热源热进口温度
Point_set(2,T_core_in,"T","=") # 热源热出口温度
Point_set(9,T₀,"T","=") # 冷端入口温度
Point_set(3,T_core_out - dT_Heater,"T","=") # 循环最高温度
## 流量点约束
Point_set(1,W_core,"W","=") # 热源热侧流量
Point_set(9,3000,"W","<") # 热源热侧流量
## 换热过程约束
# 热源
HeatProcess(1,2,dP_Heater,1,"FLiBe") # 热源热侧压降
HeatProcess(8,3,dP_Heater,2,"CO2") # 热源冷侧压降
Transfer_dT(1,2,8,3,1,2,dT_Heater,">") # 热源传热
# 回热器
HeatProcess(4,5,dP_RUP_h,3,"CO2") # 回热器热侧压降
HeatProcess(7,8,dP_RUP_l,4,"CO2") # 回热器冷侧压降
Transfer_dT(4,5,7,8,3,4,dT_RUP,">") # 回热器传热
# 冷却器
HeatProcess(5,6,dP_Con,6,"CO2") # 冷却器热侧压降
HeatProcess(9,10,dP_Air,5,"Air") # 冷却器冷侧压降
Transfer_dT(5,6,9,10,6,5,dT_Con,">") # 冷却器传热
## 变压过程约束
Compressor(6,7,1,η_com,"CO2")
Turbine(3,4,1,η_tur,"CO2")
## 热源标记
Heater(1)
## 计算热效率
set_start_value(Yita, 1.0)
#set_attributes(model, "tol" => 1e-4, "max_iter" => 1000)
Yitacal()
# 执行优化1
JuMP.optimize!(model)
print("热效率：",value.(Yita))