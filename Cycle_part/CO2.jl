#cd("E:\\Julia test codes\\CO2")
module CO2_PT2h
    using DataFrames,XLSX
    #@pyimport CoolProp.CoolProp as CP
    ## Working Directory
    #cd("E:\\Julia test codes\\Spline")

    ## Reading PT of Spline
    data_PT = XLSX.readxlsx("CO2_PT.xlsx")
    P_data = (data_PT["Sheet1"][:])'
    T_data = (data_PT["Sheet2"][:])'
    n_P = length(P_data) #Number of P
    n_T = length(T_data) #Number of T

    ################## Reading coefficients of Spline and generate PT_? ################
    
    ######################## PT_h ########################
    data_coeffs = XLSX.readxlsx("CO2_PT2H.xlsx")
    coeffs = zeros(n_P,n_T-1,4)
    for i = 1:n_P
        coeffs[i,:,:] = data_coeffs[string("Sheet",string(i))][:]
    end

    # Spline in T and linear in P
    function PT_h(P,T)
        for j = 2:n_P # Select P
            if  P <= P_data[j] && P >= P_data[j-1]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            elseif P < P_data[1]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            elseif P > P_data[end]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            end
        end
    end

    #=
    test_P = [P for P in 7.37e6:0.1e5:7.8e6]
    test_T = [T for T in 300.0:0.5:320]
    ERR = ones(length(test_P),length(test_T))
    @sync @distributed for i=1:length(test_P)
        @distributed for j=1:length(test_T)
            ERR[i,j] = abs(PT_h(test_P[i],test_T[j]) - CP.PropsSI("H","P",test_P[i],"T",test_T[j],"CO2"))/CP.PropsSI("H","P",test_P[i],"T",test_T[j],"CO2")
        end
    end
    ERR_Max = findmax(ERR)
    ERR_Max_h = ERR_Max[1]
    ERR_Max_P = test_P[ERR_Max[2][1]]
    ERR_Max_T = test_T[ERR_Max[2][2]]
    ERR_Max_PT_h = [
        ERR_Max_h,
        ERR_Max_P,
        ERR_Max_T,
        CP.PropsSI("H","P",ERR_Max_P,"T",ERR_Max_T,"CO2")
        ]
    #print("Max_informations:\n",ERR_Max_Series,"\n")
    ######################################################
    export 
        PT_h,
        ERR_Max_PT_h
    =#
export PT_h
end

module CO2_PT2s
    using DataFrames,XLSX
    ## Working Directory
    #cd("E:\\Julia test codes\\Spline")

    ## Reading PT of Spline
    data_PT = XLSX.readxlsx("CO2_PT.xlsx")
    P_data = (data_PT["Sheet1"][:])'
    T_data = (data_PT["Sheet2"][:])'
    n_P = length(P_data) #Number of P
    n_T = length(T_data) #Number of T

    ################## Reading coefficients of Spline and generate PT_? ################
    
    ######################## PT_h ########################
    data_coeffs = XLSX.readxlsx("CO2_PT2S.xlsx")
    coeffs = zeros(n_P,n_T-1,4)
    for i = 1:n_P
        coeffs[i,:,:] = data_coeffs[string("Sheet",string(i))][:]
    end

    # Spline in T and linear in P
    function PT_s(P,T)
        for j = 2:n_P # Select P
            if  P <= P_data[j] && P >= P_data[j-1]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            elseif P < P_data[1]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            elseif P > P_data[end]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            end
        end
    end

    #=
    test_P = [P for P in 7.37e6:0.1e5:7.8e6]
    test_T = [T for T in 300.0:0.5:320]
    ERR = ones(length(test_P),length(test_T))
    @sync @distributed for i=1:length(test_P)
        @distributed for j=1:length(test_T)
            ERR[i,j] = abs(PT_s(test_P[i],test_T[j]) - CP.PropsSI("S","P",test_P[i],"T",test_T[j],"CO2"))/CP.PropsSI("S","P",test_P[i],"T",test_T[j],"CO2")
        end
    end
    ERR_Max = findmax(ERR)
    ERR_Max_h = ERR_Max[1]
    ERR_Max_P = test_P[ERR_Max[2][1]]
    ERR_Max_T = test_T[ERR_Max[2][2]]
    ERR_Max_PT_s = [
        ERR_Max_h,
        ERR_Max_P,
        ERR_Max_T,
        CP.PropsSI("S","P",ERR_Max_P,"T",ERR_Max_T,"CO2")
        ]
    #print("Max_informations:\n",ERR_Max_Series,"\n")
    ######################################################
    export 
        PT_s,
        ERR_Max_PT_s
        =#
export PT_s
end

module CO2_Ph2T
    using DataFrames,XLSX
    ## Working Directory
    #cd("E:\\Julia test codes\\Spline")

    ## Reading PT of Spline
    data_PT = XLSX.readxlsx("CO2_Ph.xlsx")
    P_data = (data_PT["Sheet1"][:])'
    T_data = (data_PT["Sheet2"][:])'
    n_P = length(P_data) #Number of P
    n_T = length(T_data) #Number of T

    ################## Reading coefficients of Spline and generate PT_? ################
    
    ######################## PT_h ########################
    data_coeffs = XLSX.readxlsx("CO2_Ph2T.xlsx")
    coeffs = zeros(n_P,n_T-1,4)
    for i = 1:n_P
        coeffs[i,:,:] = data_coeffs[string("Sheet",string(i))][:]
    end

    # Spline in T and linear in P
    function Ph_T(P,T)
        for j = 2:n_P # Select P
            if  P <= P_data[j] && P >= P_data[j-1]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            elseif P < P_data[1]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            elseif P > P_data[end]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            end
        end
    end

    #=
    test_P = [P for P in 7.37e6:0.1e5:7.8e6]
    test_T = [T for T in 0.14e6:0.01e6:1.3e6]
    ERR = ones(length(test_P),length(test_T))
    @sync @distributed for i=1:length(test_P)
        @distributed for j=1:length(test_T)
            ERR[i,j] = abs(Ph_T(test_P[i],test_T[j]) - CP.PropsSI("T","P",test_P[i],"H",test_T[j],"CO2"))/CP.PropsSI("T","P",test_P[i],"H",test_T[j],"CO2")
        end
    end
    ERR_Max = findmax(ERR)
    ERR_Max_h = ERR_Max[1]
    ERR_Max_P = test_P[ERR_Max[2][1]]
    ERR_Max_T = test_T[ERR_Max[2][2]]
    ERR_Max_Ph_T = [
        ERR_Max_h,
        ERR_Max_P,
        ERR_Max_T,
        CP.PropsSI("T","P",ERR_Max_P,"H",ERR_Max_T,"CO2")
        ]
    #print("Max_informations:\n",ERR_Max_Series,"\n")
    ######################################################
    export 
        Ph_T,
        ERR_Max_Ph_T
        =#
export Ph_T  
end

module CO2_Ph2s
    using DataFrames,XLSX
    ## Working Directory
    #cd("E:\\Julia test codes\\Spline")

    ## Reading PT of Spline
    data_PT = XLSX.readxlsx("CO2_Ph.xlsx")
    P_data = (data_PT["Sheet1"][:])'
    T_data = (data_PT["Sheet2"][:])'
    n_P = length(P_data) #Number of P
    n_T = length(T_data) #Number of T

    ################## Reading coefficients of Spline and generate PT_? ################
    
    ######################## PT_h ########################
    data_coeffs = XLSX.readxlsx("CO2_Ph2s.xlsx")
    coeffs = zeros(n_P,n_T-1,4)
    for i = 1:n_P
        coeffs[i,:,:] = data_coeffs[string("Sheet",string(i))][:]
    end

    # Spline in T and linear in P
    function Ph_s(P,T)
        for j = 2:n_P # Select P
            if  P <= P_data[j] && P >= P_data[j-1]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            elseif P < P_data[1]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            elseif P > P_data[end]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            end
        end
    end

    #=
    test_P = [P for P in 7.37e6:0.1e5:7.8e6]
    test_T = [T for T in 0.14e6:0.01e6:1.3e6]
    ERR = ones(length(test_P),length(test_T))
    @sync @distributed for i=1:length(test_P)
        @distributed for j=1:length(test_T)
            ERR[i,j] = abs(Ph_s(test_P[i],test_T[j]) - CP.PropsSI("S","P",test_P[i],"H",test_T[j],"CO2"))/CP.PropsSI("S","P",test_P[i],"H",test_T[j],"CO2")
        end
    end
    ERR_Max = findmax(ERR)
    ERR_Max_h = ERR_Max[1]
    ERR_Max_P = test_P[ERR_Max[2][1]]
    ERR_Max_T = test_T[ERR_Max[2][2]]
    ERR_Max_Ph_s = [
        ERR_Max_h,
        ERR_Max_P,
        ERR_Max_T,
        CP.PropsSI("S","P",ERR_Max_P,"H",ERR_Max_T,"CO2")
        ]
    #print("Max_informations:\n",ERR_Max_Series,"\n")
    ######################################################
    export 
        Ph_s,
        ERR_Max_Ph_s
    =#
export Ph_s
end

module CO2_Ps2h
    using DataFrames,XLSX
    ## Working Directory
    #cd("E:\\Julia test codes\\Spline")

    ## Reading PT of Spline
    data_PT = XLSX.readxlsx("CO2_Ps.xlsx")
    P_data = (data_PT["Sheet1"][:])'
    T_data = (data_PT["Sheet2"][:])'
    n_P = length(P_data) #Number of P
    n_T = length(T_data) #Number of T

    ################## Reading coefficients of Spline and generate PT_? ################
    
    ######################## PT_h ########################
    data_coeffs = XLSX.readxlsx("CO2_PS2H.xlsx")
    coeffs = zeros(n_P,n_T-1,4)
    for i = 1:n_P
        coeffs[i,:,:] = data_coeffs[string("Sheet",string(i))][:]
    end

    # Spline in T and linear in P
    function Ps_h(P,T)
        for j = 2:n_P # Select P
            if  P <= P_data[j] && P >= P_data[j-1]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            elseif P < P_data[1]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            elseif P > P_data[end]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            end
        end
    end

    #=
    test_P = [P for P in 7.37e6:0.1e5:7.8e6]
    test_T = [T for T in 700.0:100:3500.0]
    ERR = ones(length(test_P),length(test_T))
    @sync @distributed for i=1:length(test_P)
        @distributed for j=1:length(test_T)
            ERR[i,j] = abs(Ps_h(test_P[i],test_T[j]) - CP.PropsSI("H","P",test_P[i],"S",test_T[j],"CO2"))/CP.PropsSI("H","P",test_P[i],"S",test_T[j],"CO2")
        end
    end
    ERR_Max = findmax(ERR)
    ERR_Max_h = ERR_Max[1]
    ERR_Max_P = test_P[ERR_Max[2][1]]
    ERR_Max_T = test_T[ERR_Max[2][2]]
    ERR_Max_Ps_h = [
        ERR_Max_h,
        ERR_Max_P,
        ERR_Max_T,
        CP.PropsSI("H","P",ERR_Max_P,"S",ERR_Max_T,"CO2")
        ]
    #print("Max_informations:\n",ERR_Max_Series,"\n")
    ######################################################
    export 
        Ps_h,
        ERR_Max_Ps_h
        =#
export Ps_h
end

module CO2_Ps2T
    using DataFrames,XLSX
    ## Working Directory
    #cd("E:\\Julia test codes\\Spline")

    ## Reading PT of Spline
    data_PT = XLSX.readxlsx("CO2_Ps.xlsx")
    P_data = (data_PT["Sheet1"][:])'
    T_data = (data_PT["Sheet2"][:])'
    n_P = length(P_data) #Number of P
    n_T = length(T_data) #Number of T

    ################## Reading coefficients of Spline and generate PT_? ################
    
    ######################## PT_h ########################
    data_coeffs = XLSX.readxlsx("CO2_PS2T.xlsx")
    coeffs = zeros(n_P,n_T-1,4)
    for i = 1:n_P
        coeffs[i,:,:] = data_coeffs[string("Sheet",string(i))][:]
    end

    # Spline in T and linear in P
    function Ps_T(P,T)
        for j = 2:n_P # Select P
            if  P <= P_data[j] && P >= P_data[j-1]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            elseif P < P_data[1]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[1] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[1] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            elseif P > P_data[end]
                for i = 2:n_T
                    if T <= T_data[i] && T >= T_data[i-1]
                        S0 = (T - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T - T_data[i-1])^3*coeffs[j,i-1,1] + (T - T_data[i-1])^2*coeffs[j,i-1,2] + (T - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T < T_data[1]
                        S0 = (T_data[1] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[1] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[1] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[1] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    elseif T > T_data[end]
                        S0 = (T_data[end] - T_data[i-1])^3*coeffs[j-1,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j-1,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j-1,i-1,3] + coeffs[j-1,i-1,4]
                        S1 = (T_data[end] - T_data[i-1])^3*coeffs[j,i-1,1] + (T_data[end] - T_data[i-1])^2*coeffs[j,i-1,2] + (T_data[end] - T_data[i-1])*coeffs[j,i-1,3] + coeffs[j,i-1,4]
                        S = (P_data[end] - P_data[j])/(P_data[j-1] - P_data[j])*S0 + (P_data[end] - P_data[j-1])/(P_data[j] - P_data[j-1])*S1   
                        return S
                    end
                end
            end
        end
    end

    #=
    test_P = [P for P in 7.37e6:0.1e5:7.8e6]
    test_T = [T for T in 700.0:100:3500.0]
    ERR = ones(length(test_P),length(test_T))
    @sync @distributed for i=1:length(test_P)
        @distributed for j=1:length(test_T)
            ERR[i,j] = abs(Ps_T(test_P[i],test_T[j]) - CP.PropsSI("T","P",test_P[i],"S",test_T[j],"CO2"))/CP.PropsSI("T","P",test_P[i],"S",test_T[j],"CO2")
        end
    end
    ERR_Max = findmax(ERR)
    ERR_Max_h = ERR_Max[1]
    ERR_Max_P = test_P[ERR_Max[2][1]]
    ERR_Max_T = test_T[ERR_Max[2][2]]
    ERR_Max_Ps_T = [
        ERR_Max_h,
        ERR_Max_P,
        ERR_Max_T,
        CP.PropsSI("T","P",ERR_Max_P,"S",ERR_Max_T,"CO2")
        ]
    #print("Max_informations:\n",ERR_Max_Series,"\n")
    ######################################################
    export 
        Ps_T,
        ERR_Max_Ps_T
        =#
export Ps_T 
end
