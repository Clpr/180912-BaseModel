# Refactored Urban Model in JULIA 6.4, the main control chapter
# ==============================================================================
# I. Add storage folders to the loading path
tmp_LOADPATH = copy(LOAD_PATH)
if tmp_LOADPATH[end] != "./modules"
    push!(LOAD_PATH,"./modules")
end
# II. Import/Using Custom Modules
    # using PolicyFunctions # Functions about Policy Function Solving (DP) -> DEPRECATED
    using Policy_Analytical  # analytical household decision making for whole lifetime
    using Policy_Analytical_Retired  # analytical household decision making for only the retired
    using GeneProc # General Processes to deal with data, ploting, reporting and others (PARTLY DEPRECATED)
    using SteadyState # Steady State Searching
    using Transition # Transition Path Searching
    using PyPlot  # plotting
    using PlottingCustom  # customized plotting short-writing function
    using DataFrames  # to process results
    using CSV  # to read in DataFrames
# ==============================================================================
# Section 0: BASIC PARMS
# ==============================================================================
    const MAX_YEAR = 400
    const MAX_AGE = 80
    const RETIRE_AGE = 40
# ==============================================================================
# Section 1: DATA READING, ROOMING
# ==============================================================================
# 1. Declare empty data structures
    include("./modules/proc_VarsDeclare.jl")  # to include in main workspace, prepare memories for economic variables
# 2. Initilize parameters
    include("./modules/proc_ParsInit.jl") # to include in main workspace, initializes parameters (easier to dynamic testing the program)
# 3. Demography Data Reading
    # 3.1 Reading in Data
        # 3.1.1 Case: Base fertility
    mat_N,mat_Survival = func_DemoRead(MAX_YEAR::Int64,MAX_AGE::Int64,RETIRE_AGE::Int64,"./data/人口0913_base.csv")
    # 3.1.2 Case: High fertility
        # mat_N,mat_Survival = func_DemoRead(MAX_YEAR::Int64,MAX_AGE::Int64,RETIRE_AGE::Int64,"./data/人口0919_high_fertility.csv")
    # 3.1.3 Case: Low fertility
        # mat_N,mat_Survival = func_DemoRead(MAX_YEAR::Int64,MAX_AGE::Int64,RETIRE_AGE::Int64,"./data/人口0919_low_fertility.csv")
    # 3.2 Process survival probs
        mat_Survival[:,end] = 1.0
        mat_Survival[ mat_Survival .> 1.0 ] = 1.0
# 4. Age-dependent Data Reading
    vec_epsilon,mat_MA2MB = func_AgeDatRead(MAX_YEAR::Int64,MAX_AGE::Int64,RETIRE_AGE::Int64, "./data/DataByAge.csv"::String,colidx_epsilon=1,colidx_MA=2,colidx_MB=3)


# ==============================================================================
# Section I: DERIVATIVE CASES (MODIFIED ON INPUT DATA & PARS)
# ==============================================================================
# I.1 Case 1: Standard Case, using other papers' parameter set
    # doing nothing
# I.2 Case 2: No tech-growth Standard Case
    # vec_tech[:] = ones(MAX_YEAR)




# ==============================================================================
# Section 2: INITIAL STEADY STATE SEARCHING
# ==============================================================================
println("="^60);println("Section: Initial Steady State Search");println("="^60)
    # 1. Additional Pars
        YEAR = 1::Int64
    # 2. Guesses
        guess_r = 0.08  # net interest rate
        guess_L = 0.6  # labour
        guess_q = 0.15  # m2c ratio
    # 3. data packaging
        DATA_vec = [ mat_N[YEAR,:], mat_Survival[YEAR,:], mat_MA2MB[YEAR,:], vec_epsilon ]
        DATA_num = [ num_kappa, num_gamma, num_varrho, num_delta, num_alpha, num_mu, num_sigma,
            vec_beta[YEAR], vec_tech[YEAR], vec_doubleK[YEAR], vec_eta[YEAR], vec_theta[YEAR], vec_z[YEAR], vec_zeta[YEAR],
            vec_phiCoef[YEAR], vec_doubleA[YEAR], vec_doubleB[YEAR], vec_cpB[YEAR]  ]
    # 4. CALLING: initial steady state searching
        @time RESULT_vec, RESULT_num = SStateSearch(guess_r,guess_L,guess_q, MAX_AGE, RETIRE_AGE,
            DATA_vec::Array{Array{Float64,1},1}, DATA_num::Array{Float64,1} ,
            flag_diagnose = false, MagicNum = 2.0, TOL_SS = 1E-08, MAX_ITER = 1000, UPDATE_STEP = 0.7 )
    # 5. Result distributed
        mat_wprofiled[YEAR,:], mat_Lambda[YEAR,:], mat_a[YEAR,:], mat_Phi[YEAR,:],
            mat_c[YEAR,:], mat_l[YEAR,1:RETIRE_AGE], mat_M[YEAR,:], mat_MA[YEAR,:], mat_MB[YEAR,:] = RESULT_vec
        vec_K[YEAR], vec_L[YEAR], vec_Y[YEAR], vec_r[YEAR], vec_wmean[YEAR], vec_oCoef[YEAR], vec_TRc[YEAR], vec_TRw[YEAR],
            vec_gapUMP[YEAR], vec_D[YEAR], vec_D2Y[YEAR], vec_G[YEAR], vec_I[YEAR], vec_doubleP[YEAR] = RESULT_num
    # 6. Output the results to csv files
        Output_SteadyState(RESULT_vec,RESULT_num, path=".\\output\\", prefix="InitSS")


# ==============================================================================
# Section 3: FINAL STEADY STATE SEARCHING
# ==============================================================================
println("="^60);println("Section: Final Steady State Search");println("="^60)
    # 1. Additional Pars
        YEAR = MAX_YEAR::Int64
    # 2. Guesses
        guess_r = 0.03  # net interest rate
        guess_L = 0.6  # labour
        guess_q = 0.25  # m2c ratio
    # 3. data packages
        DATA_vec = [ mat_N[YEAR,:], mat_Survival[YEAR,:], mat_MA2MB[YEAR,:], vec_epsilon ]
        DATA_num = [ num_kappa, num_gamma, num_varrho, num_delta, num_alpha, num_mu, num_sigma,
            vec_beta[YEAR], vec_tech[YEAR], vec_doubleK[YEAR], vec_eta[YEAR], vec_theta[YEAR], vec_z[YEAR], vec_zeta[YEAR],
            vec_phiCoef[YEAR], vec_doubleA[YEAR], vec_doubleB[YEAR], vec_cpB[YEAR]  ]
    # 4. CALLING: final steady state searching
        @time RESULT_vec, RESULT_num = SStateSearch(guess_r,guess_L,guess_q, MAX_AGE, RETIRE_AGE,
            DATA_vec, DATA_num,
            flag_diagnose = false, MagicNum = 2.0, TOL_SS = 1E-08, MAX_ITER = 100, UPDATE_STEP = 0.75 )
    # 5. Result distributed
        mat_wprofiled[YEAR,:], mat_Lambda[YEAR,:], mat_a[YEAR,:], mat_Phi[YEAR,:],
            mat_c[YEAR,:], mat_l[YEAR,1:RETIRE_AGE], mat_M[YEAR,:], mat_MA[YEAR,:], mat_MB[YEAR,:] = RESULT_vec
        vec_K[YEAR], vec_L[YEAR], vec_Y[YEAR], vec_r[YEAR], vec_wmean[YEAR], vec_oCoef[YEAR], vec_TRc[YEAR], vec_TRw[YEAR],
            vec_gapUMP[YEAR], vec_D[YEAR], vec_D2Y[YEAR], vec_G[YEAR], vec_I[YEAR], vec_doubleP[YEAR] = RESULT_num
    # 6. Output the results to csv files
        Output_SteadyState(RESULT_vec,RESULT_num, path=".\\output\\", prefix="FinaSS")


# EX: define a q ratio path
vec_Q = Array(linspace(0.15,0.25,MAX_YEAR))
    # EX: 保存一套只有初始稳态和最终稳态的数据（transition输出同款），用于最后附加的环节
    snapshot_RESULT_mat = [copy(mat_a), copy(mat_Phi), copy(mat_c), copy(mat_l), copy(mat_M), copy(mat_MA), copy(mat_MB), copy(mat_wprofiled), copy(mat_Lambda)]
    snapshot_RESULT_vec = [copy(vec_Q), copy(vec_D), copy(vec_D2Y), copy(vec_G), copy(vec_I), copy(vec_K), copy(vec_L), copy(vec_r), copy(vec_wmean), copy(vec_TRc), copy(vec_TRw), copy(vec_gapUMP), copy(vec_U), copy(vec_Y), copy(vec_oCoef)]



# ==============================================================================
# Section 4: TRANSITION PATH SEARCHING
# ==============================================================================
println("="^60);println("Section: Transition Search");println("="^60)
    # 1. Data prepare
        DATA_mat = [ mat_a::Array{Float64,2}, mat_Phi, mat_c, mat_l, mat_M, mat_MA, mat_MB, mat_wprofiled, mat_Lambda, mat_N, mat_Survival, mat_MA2MB ]
        DATA_num = [ num_alpha::Float64, num_delta, num_gamma, num_kappa, num_mu, num_sigma, num_varrho ]
        DATA_vec = [ vec_Q, vec_D::Array{Float64,1}, vec_D2Y, vec_G, vec_I, vec_K, vec_L, vec_r, vec_wmean, vec_TRc, vec_TRw, vec_gapUMP, vec_U, vec_Y, vec_oCoef, vec_beta, vec_tech, vec_cpB, vec_doubleA, vec_doubleB, vec_doubleK, vec_doubleP, vec_eta, vec_theta, vec_z,vec_zeta, vec_phiCoef, vec_x, vec_epsilon ]
    # 2. Transition searching
        RESULT_mat,RESULT_vec = TransitionSearch(MAX_YEAR, MAX_AGE, RETIRE_AGE, DATA_mat, DATA_vec, DATA_num ,
            flag_diagnose = false, MagicNum = 2.0,
            TOL_SS = 1E-02, MAX_ITER = 1000,
            LogFilepath = ".\\output\\Log_Transition.txt", UPDATE_STEP = 0.5  )
    # 3. Result distributed
        mat_a, mat_Phi, mat_c, mat_l, mat_M, mat_MA, mat_MB, mat_wprofiled, mat_Lambda = copy(RESULT_mat)
        vec_Q, vec_D, vec_D2Y, vec_G, vec_I, vec_K, vec_L, vec_r, vec_wmean, vec_TRc, vec_TRw, vec_gapUMP, vec_U, vec_Y, vec_oCoef = copy(RESULT_vec)
    # 4. writes to csv files
        Output_Transition(RESULT_mat,RESULT_vec, path=".\\output\\", prefix="Trans")
    # 5. 读入profiling数据，绘制误差收敛图，第1列是loop，第
        mat_profile = readcsv("output/Log_Transition.txt")[2:end,:]
        figure(1)
            plot(mat_profile[10:end,1],mat_profile[10:end,3]); grid(true);
            plot(mat_profile[10:end,1],mat_profile[10:end,4]); title("Error by Loop"); xlabel("Loop"); ylabel("Rel Err");
            legend(["K","L"])
            tight_layout()
    # 6. (optional) profiling 性能profiling
        # Profile.clear()  # 清除所有profiling数据
        # @profile TransitionSearch(MAX_YEAR, MAX_AGE, RETIRE_AGE, DATA_mat, DATA_vec, DATA_num ,
        #     flag_diagnose = false, MagicNum = 2.0,
        #     TOL_SS = 1E-02, MAX_ITER = 1000,
        #     LogFilepath = ".\\output\\Log_Transition.txt", UPDATE_STEP = 0.5  )
        # # 打开一个文件流
        # fio = open("./ExtraScripts/性能Profile.txt","w")
        # # 写出数据
        # Profile.print(fio)
        # # 关闭文件流
        # close(fio)


# ==============================================================================
# Section 5: PLOTTING & DISPLAY
# ==============================================================================
# 5.0 get some useful statistics
vec_C = sum(mat_N .* mat_c, 2)  # 5.0.1 social total consumption
    # 5.0.2 total population
    vec_N = sum(mat_N,2)
    # 5.0.3 working-age population ratio
    vec_LabRatio = sum(mat_N[:,1:RETIRE_AGE],2) ./ sum(mat_N,2)
    # 5.0.4 GDP growth
    vec_Ygrow = vec_Y[2:end] ./ vec_Y[1:end-1] .- 1.0
    vec_Ygrow = vec_Ygrow .* 100
    # 5.0.5 social utility
    for t in 1:MAX_YEAR
            tmp_U = [ Policy_Analytical.U(mat_c[t,s],mat_l[t,s],vec_Q[t],num_gamma,num_alpha) for s in 1:MAX_AGE ]
            vec_U[t] = sum( tmp_U .* mat_N[t,:] )
    end

# 5.1 prepare year indexing & set time period we are interested in
IDX_YEAR = 1946:1946+MAX_YEAR-1
    IDX_YEAR_lag = 1946:1946+MAX_YEAR-2
    IDXRANGE = [ 65 , 265 ]  # NOTE: pls delete the 1st and the last year, 65 = 2010
    TRANGE = [  IDX_YEAR[IDXRANGE[1]] , IDX_YEAR[IDXRANGE[2]]  ]
    TRANGE_lag = [  IDX_YEAR_lag[IDXRANGE[1]] , IDX_YEAR_lag[IDXRANGE[2]]  ]

# 5.3 using quick functiosn in module PlottingCustom to display main variables
    # NOTE: SAMPLE OF USE: QuickSubplot_Custom1(3,3,5,IDX_YEAR,vec_Y,vec_N,"GDP","Year",["Population"],IDXRANGE,TRANGE)
SUBTITLE = "Case: Standard Case with base fertility"
    # SUBTITLE = "Case: Standard Case with high fertility"
    # SUBTITLE = "Case: Standard Case with low fertility"
    # SUBTITLE = "Case: No-tech-growth Standard Case with base fertility"
    # SUBTITLE = "Case: No-tech-growth Standard Case with high fertility"
    # SUBTITLE = "Case: No-tech-growth Standard Case with low fertility"

    # # 5.3.1 plotting: Panel 1 (Population)
    #     figure(2); suptitle(SUBTITLE)
    #         QuickSubplot_Custom1(3,4,1,IDX_YEAR,vec_Y,vec_N,"GDP (Y)","Year",["Population"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,2,IDX_YEAR,vec_K,vec_N,"Capital (K)","Year",["Population"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,3,IDX_YEAR,vec_L,vec_N,"Labour (L)","Year",["Population"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,4,IDX_YEAR,vec_C,vec_N,"Consumption (C)","Year",["Population"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,5,IDX_YEAR,vec_I,vec_N,"Investment (I)","Year",["Population"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,6,IDX_YEAR,vec_G,vec_N,"Gov Purchase (G)","Year",["Population"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,7,IDX_YEAR,vec_D,vec_N,"Gov Debt (D)","Year",["Population"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,8,IDX_YEAR,vec_gapUMP,vec_N,"Gap of Pool-Medical Account (LI)","Year",["Population"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,9,IDX_YEAR,vec_U ./ vec_N,vec_N,"Social Wealfare per capita (U)","Year",["Population"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,10,IDX_YEAR,vec_r,vec_N,"Net Interest Rates (r)","Year",["Population"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,11,IDX_YEAR,vec_wmean,vec_N,"Average Wage Level (\$\\bar{w}\$)","Year",["Population"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,12,IDX_YEAR_lag,vec_Ygrow,vec_N[1:end-1],"GDP Growth Rate (%)","Year",["Population"],IDXRANGE,TRANGE_lag,DELTA=0.1)
    #     subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.05, wspace=0.60, hspace=0.35 )
    # # 5.3.1 plotting: Panel 1 (Population)
    #     figure(3); suptitle(SUBTITLE)
    #         QuickSubplot_Custom1(3,4,1,IDX_YEAR,vec_Y,vec_LabRatio,"GDP (Y)","Year",["Labour Population Ratio"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,2,IDX_YEAR,vec_K,vec_LabRatio,"Capital (K)","Year",["Labour Population Ratio"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,3,IDX_YEAR,vec_L,vec_LabRatio,"Labour (L)","Year",["Labour Population Ratio"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,4,IDX_YEAR,vec_C,vec_LabRatio,"Consumption (C)","Year",["Labour Population Ratio"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,5,IDX_YEAR,vec_I,vec_LabRatio,"Investment (I)","Year",["Labour Population Ratio"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,6,IDX_YEAR,vec_G,vec_LabRatio,"Gov Purchase (G)","Year",["Labour Population Ratio"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,7,IDX_YEAR,vec_D,vec_LabRatio,"Gov Debt (D)","Year",["Labour Population Ratio"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,8,IDX_YEAR,vec_gapUMP,vec_LabRatio,"Gap of Pool-Medical Account (LI)","Year",["Labour Population Ratio"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,9,IDX_YEAR,vec_U ./ vec_N,vec_LabRatio,"Social Wealfare per capita (U)","Year",["Labour Population Ratio"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,10,IDX_YEAR,vec_r,vec_LabRatio,"Net Interest Rates (r)","Year",["Labour Population Ratio"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,11,IDX_YEAR,vec_wmean,vec_LabRatio,"Average Wage Level (\$\\bar{w}\$)","Year",["Labour Population Ratio"],IDXRANGE,TRANGE,DELTA=0.1)
    #         QuickSubplot_Custom1(3,4,12,IDX_YEAR_lag,vec_Ygrow,vec_LabRatio[1:end-1],"GDP Growth Rate (%)","Year",["Labour Population Ratio"],IDXRANGE,TRANGE_lag,DELTA=0.1)
    #     subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.05, wspace=0.60, hspace=0.35 )
    #
    # # 5.3.2 plotting: aggregated cross-sectional utility & per capita
    #     figure(4);
    #         tmp_idx = IDXRANGE[1]:IDXRANGE[2]
    #         title("Utility: Aggregated (left axis) v.s. per capita (right axis)")
    #         plot(IDX_YEAR[tmp_idx], vec_U[tmp_idx]);ylabel("Aggregated");
    #         legend(["social welfare aggregated by population"],loc="lower left")
    #         twinx();
    #         plot(IDX_YEAR[tmp_idx], (vec_U ./ vec_N)[tmp_idx], "--r" ); ylabel("per capita");
    #         legend(["social welfare per capita"])
    #
    # # 5.3.2 plotting: consumption per capita & leisure per capita dynamics
    #     figure(5); suptitle(SUBTITLE)
    #             tmp_idx = IDXRANGE[1]:IDXRANGE[2]
    #             tmp_c_percapita = vec_C ./ vec_N
    #             tmp_l_percapita = sum(mat_N .* mat_l,2)[:,1] ./ vec_N
    #             tmp_u_percapita = Policy_Analytical.U.(tmp_c_percapita,tmp_l_percapita,vec_Q,num_gamma,num_alpha)
    #         subplot(1,2,1)
    #             title("Consumption per capita (left axis) & Leisure per capita (right axis)")
    #             plot(IDX_YEAR[tmp_idx], tmp_c_percapita[tmp_idx]);ylabel("Consumption per capita");
    #             legend(["Consumption per capita"],loc="lower left")
    #             twinx();
    #             plot(IDX_YEAR[tmp_idx], tmp_l_percapita[tmp_idx], "--r" ); ylabel("Leisure per capita");
    #             legend(["Leisure per capita (retired considered)"])
    #         subplot(1,2,2)
    #             title("Utility per capita computed by \n Consumption per capita & Leisure per capita")
    #             plot(IDX_YEAR[tmp_idx],tmp_u_percapita[tmp_idx])
    #             xlabel("Year")
    #     subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.05, wspace=0.60, hspace=0.35 )
    #


# 校准============================================================
# 1. 核算数据准备
Calib_AcctGap = CSV.read("./data/城镇职工医保核算结果.csv",header=true)  # 读入核算的缺口数据
    chk = ismissing.(Calib_AcctGap[:Year])  # 处理缺失值
    Calib_AcctGap = Calib_AcctGap[.!chk,:]
    Calib_AcctGap[:AcctGap] = Calib_AcctGap[:Expenses] .- Calib_AcctGap[:Incomes]
    # 提取2010年到2018年的gap序列，计算增长率
    Cp_AcctGap = Calib_AcctGap[:AcctGap][2010 .<= Calib_AcctGap[:Year] .<= 2018]
    Cp_AcctGapGrowth = Cp_AcctGap[2:end] ./ Cp_AcctGap[1:end-1] .- 1
    # 然后从模拟结果中提取2010到2017年的gap序列并计算其增长率
    Cp_SimuGap = vec_gapUMP[65:73]
    Cp_SimuGapGrowth = Cp_SimuGap[2:end] ./ Cp_SimuGap[1:end-1] .- 1
# 2. GDP增长率准备
Calib_GDPGrow = CSV.read("./data/CalibGDP.csv",header=true)  # 读入GDP增长
    # 提取2010年到2017年的GDP增长率序列
    Cp_AcctGDPGrowth = collect(Float64, Calib_GDPGrow[:GDPGrowth][2010 .<= Calib_GDPGrow[:Year] .<= 2017] )
    # 提取2010年到2017年的模拟GDP增长序列
    Cp_SimuGDPGrowth = vec_Ygrow[65:72]
# 3. 绘图
figure(6,figsize=[10,4.5])
    suptitle("Calibration Indicators")
    subplot(1,2,1)  # gap -------------------
    title("Gap of Pool Account")
    plot(2010:2017,Cp_AcctGapGrowth,"--b"); grid(true); xlim([2009,2018]); xlabel("Year"); ylabel("Growth")
    plot(2010:2017,Cp_SimuGapGrowth,"r"); legend(["AcctGap","SimuGap"])
    subplot(1,2,2)  # GDP growth -------------------
    title("GDP Growth")
    plot(2010:2017,Cp_AcctGDPGrowth,"--b"); grid(true); xlim([2009,2018]); xlabel("Year"); ylabel("Growth")
    plot(2010:2017,Cp_SimuGDPGrowth ./ 100,"r"); legend(["AcctGDP","SimuGDP"])
    tight_layout()




# script ends
