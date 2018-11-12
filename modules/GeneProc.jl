__precompile__()
# GeneProc.jl - Doing jobs of general works
# FuncList:
#     1. func_VarsDeclare - declare empty variables
#     1. func_ParsInit - initialize parameters
#     1. func_DemoRead - reading in demography data
#     1. func_AgeDatRead - reading in age-indexed data
#     1. print_ChapterTitle - print chapter-style title onto the screen
# =============================================
module GeneProc
    export func_VarsDeclare,func_ParsInit,func_DemoRead,func_AgeDatRead, Output_SteadyState, Output_Transition
# ====================== MODULE BEGINS =========================================

    # =============== SUBFUNCS BEGINS =======================
    function func_VarsDeclare(MAX_YEAR::Int64,MAX_AGE::Int64,RETIRE_AGE::Int64)
        # func_VarsDeclare - declare & return a series of empty variables
        # ----------------------------------------
        # personal assets
        mat_a = zeros(MAX_YEAR,MAX_AGE)
        # labour factor
        vec_L = zeros(MAX_YEAR)
        # leisure
        mat_l = ones(MAX_YEAR,MAX_AGE)
        # capital
        vec_K = zeros(MAX_YEAR)
        # gdp
        vec_Y = zeros(MAX_YEAR)
        # interest rate
        vec_r = zeros(MAX_YEAR)
        # average wage
        vec_wmean = zeros(MAX_YEAR)
        # profiled wages
        mat_wprofiled = zeros(MAX_YEAR,RETIRE_AGE)
        # wage scaling coefficient
        vec_oCoef = zeros(MAX_YEAR)
        # cross-sectional social total utility/welfare
        vec_U = zeros(MAX_YEAR)
        # consumption
        mat_c = zeros(MAX_YEAR,MAX_AGE)
        # lifetime discount factor
        # NOTE: each row for those born in year t
        mat_V = ones(MAX_YEAR,MAX_AGE)
        # ACCOUNT: individual medical accounts
        mat_Phi = zeros(MAX_YEAR,MAX_AGE)
        # GAP: individual medical accounts
        mat_gapUMI = zeros(MAX_YEAR,MAX_AGE)
        # GAP: pooling medical accounts
        vec_gapUMP = zeros(MAX_YEAR)

        # government debt
        vec_D = zeros(MAX_YEAR)
        # goverment pruchase
        vec_G = zeros(MAX_YEAR)
        # investment
        vec_I = zeros(MAX_YEAR)
        # # bequests
        # mat_b = zeros(MAX_YEAR,MAX_AGE)
        # pension benefits
        mat_Lambda = zeros(MAX_YEAR,MAX_AGE-RETIRE_AGE)

        # wage tax revenue
        vec_TRw = zeros(MAX_YEAR)
        # consumption tax revenue
        vec_TRc = zeros(MAX_YEAR)

        # medical expenses (MA+MB after profiling)
        mat_M = zeros(MAX_YEAR,MAX_AGE)
        # outpatient expenses
        mat_MA = zeros(MAX_YEAR,MAX_AGE)
        # inpatient expenses
        mat_MB = zeros(MAX_YEAR,MAX_AGE)

        # D2Y RATIO
        vec_D2Y = zeros(MAX_YEAR)

        # transfer amount from firm contribution to individual medical account (retired pahse)
        vec_doubleP = zeros(MAX_YEAR)

        return mat_a,mat_Phi,mat_c,mat_l,vec_L,vec_K,vec_Y,vec_r,vec_wmean,mat_wprofiled,
            vec_D,vec_G,vec_I,vec_D2Y,mat_Lambda,mat_M,mat_MA,mat_MB,
            vec_oCoef,vec_gapUMP,vec_TRw,vec_TRc,vec_U,vec_doubleP

    end
    # =============== SUBFUNCS ENDS =========================

    # =============== SUBFUNCS BEGINS =======================
    function func_ParsInit(MAX_YEAR::Int64,MAX_AGE::Int64,RETIRE_AGE::Int64)
        # func_ParsInit - define parameters used in the programme
        # -------------------------------------------
        # technology
        vec_tech = ones(MAX_YEAR)
            # vec_tech = 1.0 * 1.01 .^ (0:(MAX_YEAR-1)  )
        TechGrowYear = 65 + 100
        vec_tech[1:TechGrowYear] = 1.0 * 1.01 .^ (0:TechGrowYear-1)
        vec_tech[TechGrowYear+1:end] = vec_tech[TechGrowYear]

        # capital income share (maybe 0.4~0.6 for China)
        vec_beta = 0.55 * ones(MAX_YEAR)
        # depreciation rate
        num_kappa = 0.05
        # inter-temporal consumption substitute elasticity
        num_gamma = 0.5
        # inter-temporal substitute elasticity (DEPRECATED)
        num_varrho = 0.80
        # discount rate (time preference)
        num_delta = 0.01
        # leisure substitution coefficient
        num_alpha = 1.5
        # elasticity of total medical expenses to gdp (with mark 'x')
        # \Delta M / M * Y / \Delta Y
        vec_x = 1.6 * ones(MAX_YEAR)

        # ratio: cap of D/Y
        vec_doubleK = 0.00 * ones(MAX_YEAR)

        # consumption tax
        num_mu = 0.1
        # wage tax
        num_sigma = 0.24
        # contribution: firm -> pension
        vec_eta = 0.2 * ones(MAX_YEAR)
        # contribution: individual -> pension
        vec_theta = 0.08 * ones(MAX_YEAR)
        # contribution: firm -> medical
        vec_zeta = 0.06 * ones(MAX_YEAR)
        # contribution: individual -> medical
        vec_phiCoef = 0.02 * ones(MAX_YEAR)
        # transfer: firm.medical -> individual accounts
        vec_doubleA = 0.3 * ones(MAX_YEAR)
        # transfer rate from firm contribution to medical to individual medical account (retired phase)
        vec_doubleB = 0.0 * ones(MAX_YEAR)

        # Collection rate of pension 养老金收缴率
        vec_z = 0.85 * ones(MAX_YEAR)

        # Compute total contributions on nominal wages
        # # pi
        # vec_pi = ( vec_eta + vec_theta ) ./ ( 1 + vec_eta + vec_zeta )
        # # pi^M
        # vec_piM = ( vec_zeta + vec_phi ) ./ ( 1 + vec_eta + vec_zeta )
        # # co-payment rate: outpatient
        # vec_cpA = 0.0 * ones(MAX_YEAR)
        # co-payment rate: inpatient
        vec_cpB = 0.3 * ones(MAX_YEAR)
        # m2c ratio ( q ) the ratio of total medical expenses to consumption
        # mat_m2c = ones(MAX_YEAR,MAX_AGE)

        return num_kappa,num_gamma,num_varrho,num_delta,num_alpha,num_mu,num_sigma,
            vec_beta,vec_tech,vec_x,vec_doubleK,vec_eta,vec_theta,vec_zeta,vec_phiCoef,vec_doubleA,vec_doubleB,vec_cpB, vec_z
    end
    # =============== SUBFUNCS ENDS =========================

    # =============== SUBFUNCS BEGINS =======================
    function func_DemoRead(MAX_YEAR::Int64,MAX_AGE::Int64,RETIRE_AGE::Int64,csvfilepath::String)
        # func_DemoRead - Reading in Demographic Data from given excel file
        # ----------------------------------------------
        mat_Nraw = readcsv(csvfilepath::String,header=false)

        # NOTE: 如果是“人口0913.csv”，则要做一下额外处理。因为这份数据是从0岁（第1列）到101岁的。所以要先截取第20岁（第21列）到最后
            # 1. 增加新的第1年，在原来第1年基础上修正成非增的稳定分布
            mat_Nraw2 = zeros(size(mat_Nraw)[1]+1,size(mat_Nraw)[2] )
            mat_Nraw2[2:end,:] = mat_Nraw
            mat_Nraw2[1,:] = copy(mat_Nraw2[2,:])
            # 1.1 adjust the new 1st year population to a non-increasing distribution
            for s in 2:size(mat_Nraw2)[2]
                if mat_Nraw2[1,s] > mat_Nraw2[1,s-1]
                    mat_Nraw2[1,s] = mat_Nraw2[1,s-1]
                end
            end
            # 1.2 raising the non-increasing distribution to make it also non-increasing when transferred to the 2nd year (old 1st year)
            mat_Nraw2[1,:] += findmax(abs.(mat_Nraw2[1,1:end-1] .- mat_Nraw2[2,2:end]))[1]
            # 2. 截取成年后（20岁开始）
            mat_Nraw = mat_Nraw2[:,21:end]


        # then go to our normal process
        mat_Nraw = mat_Nraw[1:MAX_YEAR+1,1:MAX_AGE+1]
        # make sure it is stably distributed (non-increasing) in year 0
        for idx_Age in 2:MAX_AGE+1
            if mat_Nraw[1,idx_Age]>mat_Nraw[1,idx_Age-1]
                mat_Nraw[1,idx_Age]=mat_Nraw[1,idx_Age-1]
            end
            if mat_Nraw[MAX_YEAR,idx_Age]>mat_Nraw[MAX_YEAR,idx_Age-1]
                mat_Nraw[MAX_YEAR,idx_Age]=mat_Nraw[MAX_YEAR,idx_Age-1]
            end
        end
        # Get Survival Prob
        mat_Survival = mat_Nraw[2:end,2:end]./mat_Nraw[1:end-1,1:end-1]
        mat_Survival[1,:] = mat_Nraw[1,2:end] ./ mat_Nraw[1,1:end-1]  # 单独处理初始稳态和最终稳态
        mat_Survival[1,:] = mat_Nraw[end-1,2:end] ./ mat_Nraw[end-1,1:end-1]
        # the last age survival probabilities set to 1.0
        mat_Survival[:,end] = 1.0
        mat_Survival[ mat_Survival .> 1.0 ] = 1.0
        mat_Survival[ mat_Survival .< 0.0 ] = 0.0
        # Cut the useless part
        mat_Nraw = mat_Nraw[1:end-1,1:end-1]
        # Population Normalization
        mat_N = mat_Nraw / sum(mat_Nraw[1,:])

        return mat_N::Array{Float64,2},mat_Survival::Array{Float64,2}

    end
    # =============== SUBFUNCS ENDS =========================

    # =============== SUBFUNCS BEGINS =======================
    function func_AgeDatRead(MAX_YEAR::Int64,MAX_AGE::Int64,RETIRE_AGE::Int64,csvfilepath::String ; colidx_epsilon=1,colidx_MA=2,colidx_MB=3)
        # func_AgeDatRead - Reading in Age-dependent data, like wage profiles & MA/MB
        # ---------------------------------------------
        dat_byAge,~ = readcsv(csvfilepath,header = true)
        # Wage profile
        vec_epsilon = dat_byAge[1:RETIRE_AGE,colidx_epsilon]
        vec_epsilon = vec_epsilon / vec_epsilon[1]
        # MA/MB Ratio (p ratio)
        mat_MA2MB = repmat(dat_byAge[1:MAX_AGE,colidx_MA]'./dat_byAge[1:MAX_AGE,colidx_MB]', MAX_YEAR, 1 )

        return vec_epsilon,mat_MA2MB
    end
    # =============== SUBFUNCS ENDS =========================

    # =============== SUBFUNCS BEGINS =======================
    """
        Output_SteadyState(RESULT_vec::Array{Array{Float64,1},1},RESULT_num::Array{Float64,1} ; path=".\\output\\", prefix="SteadyState")

    outputs the result of steady state searching to assinged path with assigned prefix in .csv format.
    The function output two files, with fixed postfixes "_num.csv" & "_vec.csv".
    The outputs can be easily used to plot or demostrate.

    **Inputs**
    1. RESULT_vec: a data package output by SStateSearch()
    1. RESULT_num: a data package output by SStateSearch()
    1. path: a string assigning the folder to output
    1. prefix: a string assigning the prefix of the output files

    **Outputs**
    1. nothing
    """
    function Output_SteadyState(RESULT_vec::Array{Array{Float64,1},1},RESULT_num::Array{Float64,1} ; path::String=".\\output\\", prefix="SteadyState")
        # unfold data
        wprofiled, Lambda, apath, Phipath, Cpath, Lpath, M, MA, MB = RESULT_vec
        Final_K, Final_L, Y, r, wmean, oCoef, TRc, TRw, gapUMP, D, D2Y, G, I, doubleP = RESULT_num
        # dictionaries
        dict_vec = Dict(
            "Profiled Wage" => wprofiled,
            "Pension Benefit" => Lambda,
            "Personal Wealth Distribution" => apath,
            "Individual Medical Account Distribution" => Phipath,
            "Consumption Path" => Cpath,
            "Leisure Path" => Lpath,
            "Total Medical Expenditure" => M,
            "Outpatient Expenditure" => MA,
            "Inpatient Expenditure" => MB
         )
        dict_num = Dict(
            "Capital" => Final_K,
            "Labour" => Final_L,
            "GDP/Output" => Y,
            "Net Interest Rate" => r,
            "Social Average Wage Level" => wmean,
            "Wage Scaling Factor" => oCoef,
            "Consumption Tax Revenue" => TRc,
            "Wage Tax Revenue" => TRw,
            "Gap of Pool Medical Account" => gapUMP,
            "Government Outstanding Debt" => D,
            "D/Y Ratio" => D2Y,
            "Government Purchase" => G,
            "Investment" => I,
            "Average Transfer Amount from Firm Medical Contribution to the Retired" => doubleP
        )
        # writes to csv files
        writecsv( string( path, prefix, "_vec.csv" ), dict_vec )
        writecsv( string( path, prefix, "_num.csv" ), dict_num )

        return nothing::Void
    end
    # =============== SUBFUNCS ENDS =========================

    # =============== SUBFUNCS BEGINS =======================
    """
        Output_Transition(RESULT_mat::Array{Array{Float64,2},1},RESULT_vec::Array{Float64,1} ; path::String=".\\output\\", prefix="Transition")

    outputs the result of transition searching to assinged path with assigned prefix in .csv format.
    The function output two files, with fixed postfixes "_(VarNames).csv" & "_VectorData.csv".
    The outputs can be easily used to plot or demostrate.

    **Inputs**
    1. RESULT_mat: a data package output by SStateSearch()
    1. RESULT_vec: a data package output by SStateSearch()
    1. path: a string assigning the folder to output
    1. prefix: a string assigning the prefix of the output files

    **Outputs**
    1. nothing
    """
    function Output_Transition(RESULT_mat::Array{Array{Float64,2},1},RESULT_vec::Array{Array{Float64,1},1} ; path::String=".\\output\\", prefix="Transition")
        # unfold data
        mat_a, mat_Phi, mat_c, mat_l, mat_M, mat_MA, mat_MB, mat_wprofiled, mat_Lambda = copy(RESULT_mat)
        vec_Q, vec_D, vec_D2Y, vec_G, vec_I, vec_K, vec_L, vec_r, vec_wmean, vec_TRc, vec_TRw, vec_gapUMP, vec_U, vec_Y, vec_oCoef = copy(RESULT_vec)
        # write matrices
        writecsv( string(path,prefix,"_", "a" ,".csv") , mat_a )
        writecsv( string(path,prefix,"_", "Phi" ,".csv") , mat_Phi )
        writecsv( string(path,prefix,"_", "Wealth" ,".csv") , mat_a .+ mat_Phi )
        writecsv( string(path,prefix,"_", "c" ,".csv") , mat_c )
        writecsv( string(path,prefix,"_", "l" ,".csv") , mat_l )
        writecsv( string(path,prefix,"_", "M" ,".csv") , mat_M )
        writecsv( string(path,prefix,"_", "MA" ,".csv") , mat_MA )
        writecsv( string(path,prefix,"_", "MB" ,".csv") , mat_MB )
        writecsv( string(path,prefix,"_", "wprofiled" ,".csv") , mat_wprofiled )
        writecsv( string(path,prefix,"_", "Lambda" ,".csv") , mat_Lambda )
        # dictionaries (vector length = MAX_YEAR)
        dict_vec = Dict(
            "Gov Outstanding Debt" => vec_D,
            "D/Y Ratio" => vec_D2Y,
            "Gov Purchase" => vec_G,
            "Investment" => vec_I,
            "Capital" => vec_K,
            "Labour" => vec_L,
            "Net Interest Rate" => vec_r,
            "Average Wage Level" => vec_wmean,
            "Consumption Tax Revenue" => vec_TRc,
            "Wage Tax Revenue" => vec_TRw,
            "Gap of Pool Medical Account" => vec_gapUMP,
            "Social Utility" => vec_U,
            "GDP/Output" => vec_Y,
            "Wage Scaling Coef" => vec_oCoef
        )
        # write to csv
        writecsv(string(path,prefix,"_VectorData.csv"),dict_vec)

        return nothing::Void
    end
    # =============== SUBFUNCS ENDS =========================





# ====================== MODULE ENDS ===========================================
end
# ====================== MODULE ENDS ===========================================
