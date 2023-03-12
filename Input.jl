
## Boundary Initialize ##
p₀ = 101325 # 大气压力
T₀ = 288.15 # 大气温度
p_reactor = 2e6 # 反应堆操作压力
p_max_cycle = 22e6 # 压力上限
T_max_cycle = 1200 # 温度上限 K
W_max_cycle = 5000 # 循环流量上限 kg/s
p_min_cycle = 7.5e6 # 压力下限
T_min_cycle = 304.8 # 温度下限
p_max = max(p₀,p_max_cycle)
p_min = min(p₀,p_min_cycle)
T_max = max(T₀,T_max_cycle)
T_min = min(T₀,T_min_cycle)
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
## Initial Calculate ##
# reactor #
h_core_in = T2h_FLiBe(T_core_in) # 堆芯入口焓
h_core_out = T2h_FLiBe(T_core_out) # 堆芯出口焓
h_core_0 = T2h_FLiBe(T₀) # 熔盐在大气压力下的死态焓
s_core_in =T2s_FLiBe(T_core_in)  # 堆芯入口熵
s_core_out = T2s_FLiBe(T_core_out) # 堆芯出口熵
s_core_0 = T2s_FLiBe(T₀) # 熔盐在大气压力下的死态熵
W_core = QA/(h_core_out - h_core_in) # 堆芯流量
# 声明变量及其范围
n_process = 10 # 流股数量
n_com = 1 # 压缩机数量
n_tur = 1 # 透平数量
n_HE = 6 # 换热过程数量
n_heatsource = 1 # 热源数量
