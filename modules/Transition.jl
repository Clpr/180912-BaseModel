__precompile__()
"""
    Transition

Searches the transition path, based on given data
"""
module Transition
    using NumAlgo
    using GeneProc
    # using PolicyFunctions  # deprecated now
    using Policy_Analytical
    using Policy_Analytical_Retired

    export TransitionSearch
# ========================== MOD BEGIN =========================================
    """
        Profile_vec(vec)

    summaries vector data.
    """
    function Profile_vec(vec)
        return [ findmin(vec)[1], median(vec), mean(vec), findmax(vec)[1], vec[1], vec[end] ]
    end
    # ==========================================================================
    """
        TransitionSearch(MAX_YEAR::Int64, MAX_AGE::Int64, RETIRE_AGE::Int64,
            DATA_mat::Array{Array{Float64,2},1}, DATA_vec::Array{Array{Float64,1},1}, DATA_num::Array{Float64,1} ;
            flag_diagnose::Bool = true, MagicNum::Float64 = 2.0,
            TOL_SS::Float64 = 1E-05, MAX_ITER::Int64 = 200 , LogFilepath = "Log_perform.csv", UPDATE_STEP::Float64 = 0.3 )

    Searches the transition path when given data (init & final states have been filled).

    ##**Inputs**:
    1. MAX_YEAR: maximum year
    1. MAX_AGE: maximum age
    1. RETIRE_AGE: retire age
    1. DATA_mat: a DATA PACKAGE containing matrix data (in order)
    1. DATA_vec: a DATA PACKAGE containing vector data (in order)
    1. DATA_num: a DATA PACKAGE containing scalar data (in order)
    1. flag_diagnose: controls whether to print diagnose infos for each round
    1. MagicNum: a trick numebr used in Gauss-Seidel updating
    1. Arange: the range of policy function searching (DP)
    1. TOL_SS: relative tolerance of Gauss-Seidel algorithm
    1. TOL_POLICY: relative tolerance in policy function searching (DP)
    1. GRID_LEN: density of griding searching in policy function searching
    1. MAX_ITER: maximum iteration times of Gauu-Seidel algorithm
    1. LogFilepath: log file (any type) to write performace log (error, time-cost) to a csv file or csv-revertable file in UTF8 with titles
    1. UPDATE_STEP: relative step size to update K & L, in (0,1], similar to the category of "belief"

    ##**Outputs##:
    1. RESULT_mat: a DATA PACKAGE containing matrix data (in order)
    1. RESULT_vec: a DATA PACKAGE containing vector data (in order)

    ##**DATA PACKAGES**
    1. DATA_mat [len=12] (in order):
        1. mat_a [MAX_YEAR,MAX_AGE]: personal asset
        1. mat_Phi [MAX_YEAR,MAX_AGE]: individual medical account
        1. mat_c [MAX_YEAR,MAX_AGE]: consumption
        1. mat_l [MAX_YEAR,MAX_AGE]: leisure
        1. mat_M [MAX_YEAR,MAX_AGE]: total medical expenses
        1. mat_MA [MAX_YEAR,MAX_AGE]: outpatient expenses
        1. mat_MB [MAX_YEAR,MAX_AGE]: inpatient expenses
        1. mat_wprofiled [MAX_YEAR,RETIRE_AGE]: profiled wage level
        1. mat_Lambda [MAX_YEAR,MAX_AGE-RETIRE_AGE]: pension benefits
        1. (------- cutline ---------)
        1. mat_N [MAX_YEAR,MAX_AGE]: demography
        1. mat_Survival [MAX_YEAR,MAX_AGE]: survival prob
        1. mat_MA2MB [MAX_YEAR,MAX_AGE]: ratio of outpatient to inpatient
    1. DATA_vec [len=28] (in order):
        1. vec_D [MAX_YEAR]: government outstanding debt
        1. vec_D2Y [MAX_YEAR]: debt to GDP ratio
        1. vec_G [MAX_YEAR]: government purchase
        1. vec_I [MAX_YEAR]: investment
        1. vec_K [MAX_YEAR]: capital
        1. vec_L [MAX_YEAR]: labour
        1. vec_r [MAX_YEAR]: net interest rate
        1. vec_wmean [MAX_YEAR]: social average wage level
        1. vec_TRc [MAX_YEAR]: consumption tax revenue
        1. vec_TRw [MAX_YEAR]: wage tax revenue
        1. vec_gapUMP [MAX_YEAR]: gap of pooling medical account
        1. vec_U [MAX_YEAR]: social cross-sectional utility
        1. vec_Y [MAX_YEAR]: GDP
        1. vec_oCoef [MAX_YEAR]: wage scaling coefficient
        1. (------- cutline ---------)
        1. vec_beta [MAX_YEAR]: capital income share
        1. vec_tech [MAX_YEAR]: technology, the TFP
        1. vec_cpB [MAX_YEAR]: copayment rate of inpatient expenses
        1. vec_doubleA [MAX_YEAR]: transfer rate from firm contribution to individual medical account (working phase)
        1. vec_doubleB [MAX_YEAR]: transfer rate from firm contribution to individual medical account (retired phase)
        1. vec_doubleK [MAX_YEAR]: cap of D/Y ratio
        1. vec_doubleP [MAX_YEAR]: transfer amount from firm contribution to individual medical account (retired phase)
        1. vec_eta [MAX_YEAR]: firm contribution rate to pension
        1. vec_theta [MAX_YEAR]: personal contribution rate to pension
        1. vec_z [MAX_YEAR]: collection rate of pension
        1. vec_zeta [MAX_YEAR]: firm contribution rate to medical
        1. vec_phiCoef [MAX_YEAR]: personal contribution rate to medical
        1. vec_x [MAX_YEAR]: income elasticity of social medical expenses
        1. vec_epsilon [RETIRE_AGE]: wage profile coefficients/distribution
    1. DATA_num [len=7] (in order):
        1. Alpha: preference coef of leisure than Consumption
        1. Delta: utility discount rate
        1. Gamma: inter-temporal elasticity
        1. Kappa: depreciation rate of capital
        1. Mu: consumption tax rate
        1. Sigma: wage tax rate
        1. Varrho: intra-temporal consumption elasticity of leisure
    1. RESULT_mat [len=9] (in order):
        1. mat_a [MAX_YEAR,MAX_AGE]: personal asset
        1. mat_Phi [MAX_YEAR,MAX_AGE]: individual medical account
        1. mat_c [MAX_YEAR,MAX_AGE]: consumption
        1. mat_l [MAX_YEAR,MAX_AGE]: leisure
        1. mat_M [MAX_YEAR,MAX_AGE]: total medical expenses
        1. mat_MA [MAX_YEAR,MAX_AGE]: outpatient expenses
        1. mat_MB [MAX_YEAR,MAX_AGE]: inpatient expenses
        1. mat_wprofiled [MAX_YEAR,RETIRE_AGE]: profiled wage level
        1. mat_Lambda [MAX_YEAR,MAX_AGE-RETIRE_AGE]: pension benefits
    1. RESULT_vec [len=14] (in order):
        1. vec_Q [MAX_YEAR]: medical expenses to consumption ratio
        1. vec_D [MAX_YEAR]: government outstanding debt
        1. vec_D2Y [MAX_YEAR]: debt to GDP ratio
        1. vec_G [MAX_YEAR]: government purchase
        1. vec_I [MAX_YEAR]: investment
        1. vec_K [MAX_YEAR]: capital
        1. vec_L [MAX_YEAR]: labour
        1. vec_r [MAX_YEAR]: net interest rate
        1. vec_wmean [MAX_YEAR]: social average wage level
        1. vec_TRc [MAX_YEAR]: consumption tax revenue
        1. vec_TRw [MAX_YEAR]: wage tax revenue
        1. vec_gapUMP [MAX_YEAR]: gap of pooling medical account
        1. vec_U [MAX_YEAR]: social cross-sectional utility
        1. vec_Y [MAX_YEAR]: GDP
        1. vec_oCoef [MAX_YEAR]: wage scaling coefficient
    """
    function TransitionSearch(MAX_YEAR::Int64, MAX_AGE::Int64, RETIRE_AGE::Int64,
        DATA_mat::Array{Array{Float64,2},1}, DATA_vec::Array{Array{Float64,1},1}, DATA_num::Array{Float64,1} ;
        flag_diagnose::Bool = true, MagicNum::Float64 = 2.0,
        TOL_SS::Float64 = 1E-05, MAX_ITER::Int64 = 200 , LogFilepath = "Log_perform.csv", UPDATE_STEP::Float64 = 0.3 )
        # ----------------------------------------------------------------------
        # Section 1: DATA PKG VALIDATION & PROCESS 数据拆包与预处理
        # ----------------------------------------------------------------------
        # 1.1. Validation 合法性检查
            @assert((length(DATA_mat)==12)&&(length(DATA_vec)==29)&&(length(DATA_num)==7) , "Uncompatible DATA_mat, DATA_vec or DATA_num!"  )
        # 1.2. Unfolds data packages (call by reference) 拆包数据（派发引用）
            mat_a, mat_Phi, mat_c, mat_l, mat_M, mat_MA, mat_MB, mat_wprofiled, mat_Lambda, mat_N, mat_Survival, mat_MA2MB = DATA_mat
            Alpha, Delta, Gamma, Kappa, Mu, Sigma, Varrho = DATA_num
            vec_Q, vec_D, vec_D2Y, vec_G, vec_I, vec_K, vec_L, vec_r, vec_wmean, vec_TRc, vec_TRw, vec_gapUMP, vec_U, vec_Y, vec_oCoef, vec_beta, vec_tech, vec_cpB, vec_doubleA, vec_doubleB, vec_doubleK, vec_doubleP, vec_eta, vec_theta, vec_z, vec_zeta, vec_phiCoef, vec_x, vec_epsilon = DATA_vec
        # 1.3 FILL K,L,D,Q PATHS (LINEAR) as initial guesses 线性插值填充资本、劳动、闲暇、M2C比率等作为初始guesses
            # 1.3.1 using reference rather than copy
                vec_K[1:end] = Array(linspace(vec_K[1],vec_K[MAX_YEAR],MAX_YEAR) )
                vec_L[1:end] = Array(linspace(vec_L[1],vec_L[MAX_YEAR],MAX_YEAR) )
                mat_l[:,1:RETIRE_AGE] = repmat( (1-vec_L)/RETIRE_AGE , 1, RETIRE_AGE )
            # 1.3.2 specially define a mat for m/c (q) ratio
                # NOTE: convenient to select data, though more memories spent
                # vec_Q[1:end] = Array( linspace(mat_M[1,1]/mat_c[1,1],mat_M[MAX_YEAR,1]/mat_c[MAX_YEAR,1],MAX_YEAR) )
        # 1.4 Prepare empty vectors of K,L (using copy but not references)
            vec_K2 = copy(vec_K); vec_L2 = copy(vec_L)
        # 1.5 PREPARE INDEXES FOR INNER LOOPING (to make it clear) 准备索引缩写
            IDX_T = 2:MAX_YEAR-1
            AGE_R = RETIRE_AGE+1:MAX_AGE
            AGE_W = 1:RETIRE_AGE
        # 1.6 COPY INITIAL STEADY STATES 准备一份稳态的拷贝以防止搜索过程中数据读写对稳态的污染
            # NOTE: to avoid pollution raised by data I/O in looping
            copy_i_a = copy(mat_a[1,1:MAX_AGE]); copy_i_Phi = copy(mat_Phi[1,1:MAX_AGE])
            copy_i_c = copy(mat_c[1,1:MAX_AGE]); copy_i_l = copy(mat_l[1,1:MAX_AGE])
        # 1.7 COPY FINAL STEADY STATES
            copy_f_a = copy(mat_a[MAX_YEAR,1:MAX_AGE]); copy_f_Phi = copy(mat_Phi[MAX_YEAR,1:MAX_AGE])
            copy_f_c = copy(mat_c[MAX_YEAR,1:MAX_AGE]); copy_f_l = copy(mat_l[MAX_YEAR,1:MAX_AGE])
        # 1.8 PREPARE SHORT WRITINGS FOR SOCIAL SECURITY CONTRIBUTIONS
            vec_Pi = vec_z.*(vec_theta.+vec_eta)./(1.+vec_z.*vec_eta.+vec_zeta)

        # ----------------------------------------------------------------------
        # Section 2: LOOPING 迭代求解 Designed Gauss-Seidel Algorithm
        # ----------------------------------------------------------------------
        for idx_loop in 1:MAX_ITER
        # 0. print basic infomation & timing & validations
            println("Round: ",idx_loop)
            tic();
            # non-NaN validations
            @assert( !any(isnan.(vec_K)) , "NaN found in K!" )
            @assert( !any(isnan.(vec_L)) , "NaN found in L!" )

        # 1. Firm Department & Pension System 厂商部门&养老金 -------------------
            # NOTE: quite similar to SteadyState Searching, if no special notations
            for t in IDX_T
                # 1.0 GDP/Output
                    vec_Y[t]     = vec_tech[t] * vec_K[t]^vec_beta[t] * vec_L[t]^(1-vec_beta[t])
                # 1.1 interest rate (net investment returns)
                    vec_r[t]     = vec_beta[t]*vec_tech[t]*(vec_K[t]/vec_L[t])^(vec_beta[t]-1) - Kappa
                # 1.2 social average wage level
                    vec_wmean[t] = (1-vec_beta[t])*vec_tech[t]*(vec_K[t]/vec_L[t])^vec_beta[t]
                # 1.3 wage profiling
                    vec_oCoef[t] = vec_L[t] / sum(  mat_N[t,AGE_W].*vec_epsilon.*(1-mat_l[t,AGE_W])  )
                    mat_wprofiled[t,AGE_W] = vec_wmean[t]*vec_oCoef[t]*vec_epsilon
                # 1.4 pension benefits
                    avg_LambdaLevel = sum(  vec_Pi[t]*mat_wprofiled[t,AGE_W].*mat_N[t,AGE_W].*(1-mat_l[t,AGE_W])  ) / sum(  mat_N[t,AGE_R]  )
                    mat_Lambda[t,:] = avg_LambdaLevel
            end

        # 2. Household Department 家庭部门 --------------------------------------
            # NOTE: this section is the most challenging part in the algorithm;
            # we have several different types of decision to make:
            # 1. decision for those who alive in year 0 and will be still in the economy on the transition path; these people (usually) have non-zero wealth start point, and even start from years post-retired, we separately compute paths for the two types of people;
            # 2. decision for those who are born in transition years but will not meet the final steady state (die before MAX_YEAR); these people have complete dicision making problem, easy to play with;
            # 3. decision for those who have at least one year reached the final steady state (MAX_YEAR); we need to construct complete data seires for their whole life then make decision; complete problem though, most work on data operations & indexing computation
            # --------------------------------
            # 2.1 re-compute new paths for those alive in year 0 为初始稳态仍存活的人重新计算路径
                # 2.1.1 for those have not retired in year 0 开始时尚未退休的那部分人
                for s in AGE_W
                    # 2.1.1.1 different cohorts, different path lengths & start point positions
                        LEN = MAX_AGE - s + 1
                        LEN_W = RETIRE_AGE - s + 1
                        LEN_R = MAX_AGE - RETIRE_AGE
                    # 2.1.1.2 data series (we write each of them line by line for clearer debugging purpose)
                        Rpath = vec_r[1:LEN]
                        Fpath = 1-diag(mat_Survival,s-1)
                        Qpath = vec_Q[1:LEN]
                        Ppath = diag(mat_MA2MB,s-1)
                        cpBpath = vec_cpB[1:LEN]
                            Wpath = diag(mat_wprofiled,s-1)
                            phiCoefpath = vec_phiCoef[1:LEN_W]
                            Zetapath = vec_zeta[1:LEN_W]
                            Etapath = vec_eta[1:LEN_W]
                            Thetapath = vec_theta[1:LEN_W]
                            zpath = vec_z[1:LEN_W]
                            doubleApath = vec_doubleA[1:LEN_W]
                        Lambdapath = diag(mat_Lambda,0)
                        doublePpath = vec_doubleP[1:LEN_R]
                    # 2.1.1.3 data pacakges (Dict)
                        DATA = AmmoReload_DATA( LEN, Rpath, Fpath, Qpath, Ppath, cpBpath )
                        DATA_w = AmmoReload_DATA_w( LEN_W::Int64, Wpath, phiCoefpath, Zetapath, Etapath, Thetapath, zpath, doubleApath   )
                        DATA_r = AmmoReload_DATA_r( LEN, LEN_W, Lambdapath, doublePpath   )
                        PARS = AmmoReload_PARS( Mu, Sigma, Alpha, Gamma, Delta  )
                    # 2.1.1.4 compute new paths
                        apath, Phipath, Cpath, Lpath = PolicySolve_Analytical(mat_a[1,s],mat_Phi[1,s],
                            LEN, LEN_W, DATA, DATA_w, DATA_r, PARS   )
                    # 2.1.1.5 record the results (using diagwrite() in module GeneProc)
                        # NOTE: directly writing to the original variables but not copies (pls refer to the documentation of diagwrite())
                        diagwrite(mat_a,apath,offset=s-1)
                        diagwrite(mat_Phi,Phipath,offset=s-1)
                        diagwrite(mat_c,Cpath,offset=s-1)
                        diagwrite(mat_l,Lpath,offset=s-1)
                end
                # 2.1.2 for those have retired in year 0 开始时已经退休的那部分人
                for s in AGE_R
                    # 2.1.2.1 new path lengths
                        LEN = MAX_AGE - s + 1
                    # 2.1.2.2 data series
                        Rpath = vec_r[1:LEN]
                        Fpath = 1-diag(mat_Survival,s-1)
                        Qpath = vec_Q[1:LEN]
                        Ppath = diag(mat_MA2MB,s-1)
                        cpBpath = vec_cpB[1:LEN]
                            Lambdapath = diag(mat_Lambda,s-RETIRE_AGE-1)
                            doublePpath = vec_doubleP[1:LEN]
                    # 2.1.2.3 data packing
                        DATA = AmmoReload_Retired_DATA( MAX_AGE, s, Rpath, Fpath, Qpath, Ppath, cpBpath, Lambdapath, doublePpath  )
                        PARS = AmmoReload_Retired_PARS( Mu, Gamma, Delta  )
                    # 2.1.2.4 compute new paths
                        apath, Phipath, Cpath = PolicySolve_Analytical_Retired(mat_a[1,s],mat_Phi[1,s], MAX_AGE , s, DATA, PARS   )
                    # 2.1.2.5 record the results
                        diagwrite(mat_a,apath,offset=s-1)
                        diagwrite(mat_Phi,Phipath,offset=s-1)
                        diagwrite(mat_c,Cpath,offset=s-1)
                end
                # 2.1.3 info printing to report & monitor the algorithm process
                    println('\t',"+ Policy Done: for those alive in year 0")
                # 2.1.4 Refresh the initial steady state to avoid I/O pollution
                    mat_a[1,1:MAX_AGE] = copy_i_a; mat_Phi[1,1:MAX_AGE] = copy_i_Phi
                    mat_c[1,1:MAX_AGE] = copy_i_c; mat_l[1,1:MAX_AGE] = copy_i_l

            # 2.2 paths for those born on the transition path 为转轨路径上新出生的人计算决策路径
                # NOTE: but they won't reach the final steady state
                for t in 2:MAX_YEAR-MAX_AGE+1  # until the last generation who has just the very last year in the final steady state
                    # 2.2.1 new length marks
                        TMPIDX = t : t+MAX_AGE-1
                        TMPIDX_W = t : t+RETIRE_AGE-1;
                        TMPIDX_R = t+RETIRE_AGE : t+MAX_AGE-1
                    # 2.2.2 data series
                        Rpath = vec_r[TMPIDX]
                        Fpath = 1-diag(mat_Survival,1-t)
                        Qpath = vec_Q[TMPIDX]
                        Ppath = diag(mat_MA2MB,1-t)
                        cpBpath = vec_cpB[TMPIDX]
                            Wpath = diag(mat_wprofiled,1-t)
                            phiCoefpath = vec_phiCoef[TMPIDX_W]
                            Zetapath = vec_zeta[TMPIDX_W]
                            Etapath = vec_eta[TMPIDX_W]
                            Thetapath = vec_theta[TMPIDX_W]
                            zpath = vec_z[TMPIDX_W]
                            doubleApath = vec_doubleA[TMPIDX_W]
                        Lambdapath = diag(mat_Lambda,1-t)
                        doublePpath = vec_doubleP[TMPIDX_R]
                    # 2.2.3 data pacakges (Dict)
                        DATA = AmmoReload_DATA( MAX_AGE, Rpath, Fpath, Qpath, Ppath, cpBpath )
                        DATA_w = AmmoReload_DATA_w( RETIRE_AGE::Int64, Wpath, phiCoefpath, Zetapath, Etapath, Thetapath, zpath, doubleApath   )
                        DATA_r = AmmoReload_DATA_r( MAX_AGE, RETIRE_AGE, Lambdapath, doublePpath   )
                        PARS = AmmoReload_PARS( Mu, Sigma, Alpha, Gamma, Delta  )
                    # 2.2.4 compute new paths
                        apath, Phipath, Cpath, Lpath = PolicySolve_Analytical(0.0,0.0, MAX_AGE, RETIRE_AGE, DATA, DATA_w, DATA_r, PARS )
                    # 2.2.5 record the results
                        diagwrite(mat_a,apath,offset=1-t)
                        diagwrite(mat_Phi,Phipath,offset=1-t)
                        diagwrite(mat_c,Cpath,offset=1-t)
                        diagwrite(mat_l,Lpath,offset=1-t)
                end
                # 2.2.6 Refresh the final steady state
                    mat_a[MAX_YEAR,1:MAX_AGE] = copy_f_a; mat_Phi[MAX_YEAR,1:MAX_AGE] = copy_f_Phi
                    mat_c[MAX_YEAR,1:MAX_AGE] = copy_f_c; mat_l[MAX_YEAR,1:MAX_AGE] = copy_f_l
                # 2.2.7 info printing
                    println('\t',"+ Policy Done: for those born in transition")

            # 2.3 paths for those having more than 1 year reaching the final steady state 为那些至少有2年碰到了最终稳态的个体计算路径
                # NOTE: it means we need to expand our data; but pls remember, we still compute WHOLE path for them
                for t in MAX_YEAR-MAX_AGE+2 : MAX_YEAR-1
                    # 2.3.1 new indexes & locations
                        TMPLOC = t; TMPIDX = t:MAX_YEAR
                        TMPAGE = MAX_YEAR+1-t
                        TMPLOC_R = min( t+RETIRE_AGE , MAX_YEAR )
                    # 2.3.2 Data expansion & series (VecExpand() & VecExpand_ContentSpecific() in module NumAlgo)
                        Rpath = VecExpand(vec_r[TMPIDX],MAX_AGE,Direct="forward")
                        Fpath = VecExpand_ContentSpecific(1-diag(mat_Survival,1-t),(1-mat_Survival[MAX_YEAR,1:MAX_AGE]),MAX_AGE,Direct="forward"  )
                        Qpath = VecExpand(vec_Q[TMPIDX],MAX_AGE,Direct="forward")
                        Ppath = VecExpand_ContentSpecific(diag(mat_MA2MB,1-t),mat_MA2MB[MAX_YEAR,1:MAX_AGE],MAX_AGE,Direct="forward"  )
                        cpBpath = VecExpand(vec_cpB[TMPIDX],MAX_AGE,Direct="forward")
                            Wpath = VecExpand_ContentSpecific(diag(mat_wprofiled,1-t),mat_wprofiled[MAX_YEAR,1:RETIRE_AGE],RETIRE_AGE,Direct="forward"  )
                            phiCoefpath = VecExpand(vec_cpB[TMPIDX],RETIRE_AGE,Direct="forward")
                            Zetapath = VecExpand(vec_zeta[TMPIDX],RETIRE_AGE,Direct="forward")
                            Etapath = VecExpand(vec_eta[TMPIDX],RETIRE_AGE,Direct="forward")
                            Thetapath = VecExpand(vec_theta[TMPIDX],RETIRE_AGE,Direct="forward")
                            zpath = VecExpand(vec_z[TMPIDX],RETIRE_AGE,Direct="forward")
                            doubleApath = VecExpand(vec_doubleA[TMPIDX],RETIRE_AGE,Direct="forward")
                        if TMPAGE >= RETIRE_AGE
                            Lambdapath = VecExpand_ContentSpecific(diag(mat_Lambda,1-t),mat_Lambda[MAX_YEAR,:],MAX_AGE-RETIRE_AGE,Direct="forward"  )
                        else
                            Lambdapath = mat_Lambda[MAX_YEAR,:]
                        end
                        doublePpath = VecExpand(vec_doubleP[TMPIDX],RETIRE_AGE,Direct="forward")
                    # 2.3.3 Data packaging
                        DATA = AmmoReload_DATA( MAX_AGE, Rpath, Fpath, Qpath, Ppath, cpBpath )
                        DATA_w = AmmoReload_DATA_w( RETIRE_AGE, Wpath, phiCoefpath, Zetapath,
                                Etapath, Thetapath, zpath, doubleApath   )
                        DATA_r = AmmoReload_DATA_r( MAX_AGE, RETIRE_AGE, Lambdapath, doublePpath   )
                        PARS = AmmoReload_PARS( Mu, Sigma, Alpha, Gamma, Delta  )
                    # 2.3.4 compute new policy functions
                        apath, Phipath, Cpath, Lpath = PolicySolve_Analytical(0.0,0.0, MAX_AGE,RETIRE_AGE, DATA, DATA_w, DATA_r, PARS   )
                    # 2.3.5 record the results
                        diagwrite(mat_a,apath[1:TMPAGE]     ,offset=1-t)
                        diagwrite(mat_Phi,Phipath[1:TMPAGE] ,offset=1-t)
                        diagwrite(mat_c,Cpath[1:TMPAGE]     ,offset=1-t)
                        diagwrite(mat_l,Lpath[1:min(TMPAGE,RETIRE_AGE)]     ,offset=1-t)
                end
                # 2.3.6 Refresh the final steady state
                    mat_a[MAX_YEAR,1:MAX_AGE] = copy_f_a; mat_Phi[MAX_YEAR,1:MAX_AGE] = copy_f_Phi
                    mat_c[MAX_YEAR,1:MAX_AGE] = copy_f_c; mat_l[MAX_YEAR,1:MAX_AGE] = copy_f_l
                # 2.3.7 info printing
                    println('\t',"+ Policy Done: for those having reached the final steady state")


        # 3. Social Medical Security System 医保系统 ----------------------------
            # 3.1 update medical expenses
                for t in IDX_T
                    # 4.1 Medical expenses based on updated consumption & leisure
                        mat_M[t,:] = vec_Q[t]*mat_c[t,:]    # total expenses
                        mat_MA[t,:] = mat_M[t,:].*mat_MA2MB[t,:]./(1+mat_MA2MB[t,:])    # outpatient
                        mat_MB[t,:] = mat_M[t,:]./(1+mat_MA2MB[t,:])    # inpatient
                end
            # 3.2 Compute q (m/c) ratio through Y growth
                for t in IDX_T
                    # 3.2.1 Social total medical expenditure
                        Mtotal = (  vec_x[t-1]*(vec_Y[t]/vec_Y[t-1]-1) + 1  ) * sum( mat_N[t-1,:].*mat_M[t-1,:] )
                    # 3.2.2 Trans to expenses for everyone
                        # NOTE: there is a cap of 0.25 (European level)
                        # vec_Q[t] = max( 0.0, min(0.25  ,  Mtotal / sum( mat_N[t,:].*mat_c[t,:] )  )    )
                        vec_Q[t] = vec_Q[t]
                end

        # 4. Fiscal & Labour Accounting 财政&劳动力核算 --------------------------
            # 4.1 taxes, gaps & labour re-aggregation
                for t in IDX_T
                    # 4.1.1 Tax Revenue
                        vec_TRc[t] = Mu * sum(  mat_N[t,:].*mat_c[t,:]  )   # consumption tax
                        vec_TRw[t] = Sigma * sum(  mat_N[t,AGE_W].*mat_wprofiled[t,AGE_W].*(1-mat_l[t,AGE_W])  )    # wage tax
                    # 4.1.2 Pooling medical account gaps (income - expense)
                        vec_gapUMP[t] = sum( (1-vec_cpB[t])*mat_MB[t,:].*mat_N[t,:] ) -
                            ( (1-vec_doubleA[t]-vec_doubleB[t])*vec_zeta[t]/(1+vec_eta[t]+vec_zeta[t]) ) *
                            sum( mat_N[t,AGE_W].*mat_wprofiled[t,AGE_W].*(1-mat_l[t,AGE_W]) )
                    # 4.1.3 Aggregate Labour factor
                        vec_L2[t] = sum(  mat_N[t,AGE_W].*(1 - mat_l[t,AGE_W])  )
                    # 4.1.4 compute a new GDP with new labour and old capital (Gauss-Seidel style updating)
                        vec_Y[t] = vec_tech[t]*(vec_K[t]^vec_beta[t])*(vec_L2[t]^(1-vec_beta[t]))
                end
            # 4.2 government outstanding debt
                for t in IDX_T
                    # 4.2.1 through capital market equilibrium
                        vec_D[t] = vec_K[t] - sum(mat_N[t,:].*(mat_a[t,:].+mat_Phi[t,:]))
                    # 4.2.2 compute debt ratio & set a cap
                        vec_D2Y[t] = max( 0.0 , min( vec_doubleK[t] , vec_D[t] ./ vec_Y[t] ) )
                    # 4.2.3 update debt
                        vec_D[t] = vec_Y[t] * vec_D2Y[t]
                end

        # 5. Government Purchase, Investment & Capital 政府购买，投资&资本更新 -----------------
            # 5.1 government purchase from fiscal budget dynamics
                for t in IDX_T
                    vec_G[t] = vec_TRc[t] + vec_TRw[t] + vec_D[t+1] - vec_gapUMP[t] - vec_r[t]*vec_D[t]
                end
            # 5.2 get investment from good market clearing
                for t in IDX_T
                    vec_I[t] = vec_Y[t] - vec_G[t] - sum( mat_N[t,:].*mat_c[t,:] )
                end
            # 5.3 get capital from capital growth dynamics
                # NOTE: from the last year to the beginning
                for t in reverse(IDX_T)
                    vec_K2[t] = vec_K[t+1] - vec_I[t]
                    vec_K2[t] /= 1 - Kappa
                end


        # 6. Convergence Check 收敛检查 -------------------------------------------
            # 6.1 Compute maximum relative errors [K,L,Q]
                Err = [ findmax(abs.(vec_K2./vec_K-1))[1] , findmax(abs.(vec_L2./vec_L-1))[1]  ]
            # 6.2 Diagnose printing (using Profile_vec() to do basic descriptive statistics)
                if flag_diagnose
                    println('\t',"-"^50)
                    println('\t',"+ Diagnoses:")
                    println("\t\t","- format   : [min,  median,  mean,  max,  initSS,  finalSS]")
                    println("\t\t","- K current: ", Profile_vec(vec_K) )
                    println("\t\t","- K re-agg : ", Profile_vec(vec_K2) )
                    println("\t\t","- L current: ", Profile_vec(vec_L) )
                    println("\t\t","- L re-agg : ", Profile_vec(vec_L2) )
                    println("\t\t","- Q        : ", Profile_vec(vec_Q) )
                    println("\t\t","- Y        : ", Profile_vec(vec_Y) )
                    println("\t\t","- r        : ", Profile_vec(vec_r) )
                    println("\t\t","- w        : ", Profile_vec(vec_wmean) )
                    println("\t\t","- C        : ", Profile_vec( sum( mat_N.*mat_c, 2 ) ) )
                    println("\t\t","- I        : ", Profile_vec(vec_I) )
                    println("\t\t","- D        : ", Profile_vec(vec_D) )
                    println("\t\t","- D/Y      : ", Profile_vec(vec_D2Y) )
                    println("\t\t","- gapUMP   : ", Profile_vec(vec_gapUMP) )
                    println("\t\t","- TaxRev C : ", Profile_vec(vec_TRc) )
                    println("\t\t","- TaxRev W : ", Profile_vec(vec_TRw) )
                    println("\t\t","- Max a+Phi: ", findmax(mat_a.+mat_Phi)[1] )
                    println("\t\t","- Pension  : ", Profile_vec(mat_Lambda[:,1]) )
                end
            # 6.3 Check & Update 检查收敛并更新进入下一轮
                if all(Err .< TOL_SS)
                    println("\t+ Converged.\n\t+ Error [K,L,Q]: ",Err)
                    print("\t+ ");toc() # stop timing & print time cost
                    # packaging data
                    RESULT_mat = [mat_a, mat_Phi, mat_c, mat_l, mat_M, mat_MA, mat_MB, mat_wprofiled, mat_Lambda]
                    RESULT_vec = [vec_Q, vec_D, vec_D2Y, vec_G, vec_I, vec_K, vec_L, vec_r, vec_wmean, vec_TRc, vec_TRw, vec_gapUMP, vec_U, vec_Y, vec_oCoef]
                    # then return
                    return RESULT_mat, RESULT_vec
                elseif idx_loop==MAX_ITER
                    println("\t+ Max Iteration Reached.\n\t+ Error [K,L,Q]: ",Err)
                    print("\t+ ");toc()
                    RESULT_mat = [mat_a, mat_Phi, mat_c, mat_l, mat_M, mat_MA, mat_MB, mat_wprofiled, mat_Lambda]
                    RESULT_vec = [vec_Q, vec_D, vec_D2Y, vec_G, vec_I, vec_K, vec_L, vec_r, vec_wmean, vec_TRc, vec_TRw, vec_gapUMP, vec_U, vec_Y, vec_oCoef]
                    return RESULT_mat, RESULT_vec
                else
                    # 6.3.1 print errors
                        println("\t+ Error [K,L]: ",Err)
                    # 6.3.2 get elapsed time & print
                        print("\t+ ");flag_time = toc();
                    # 6.3.3 write logs (to a plain text file; the readers may change its postfix as .csv then use excel to insight it)
                        # NOTE: record: loop, time-cost of current loop, relative error of capital, relative error of labour
                        if idx_loop == 1
                            open(LogFilepath,"w") do file
                                err_K , err_L = Err
                                write(file, "Loop,Elapsed(s),Err_K,Err_L\n")
                                write(file, "$idx_loop,$flag_time,$err_K,$err_L\n")
                            end
                        else
                            open(LogFilepath,"a") do file
                                err_K , err_L = Err
                                write(file, "$idx_loop,$flag_time,$err_K,$err_L\n")
                            end
                        end
                    # 6.3.4 bounding L, Q
                        if any(vec_L2.<0.0) || !isreal(vec_L2)
                            vec_L2[ (vec_L2.<0.0) .| (vec_L2.!=real(vec_L2))  ] = 0.001
                        end
                    # 6.3.5 update L
                        vec_L[IDX_T] = vec_L[IDX_T] .+ UPDATE_STEP .* (vec_L2[IDX_T] .- vec_L[IDX_T])
                    # 6.3.6 set a floor for capital
                        tmp_floor = MagicNum * findmin(vec_L)[1]
                    # 6.3.7 bounding K
                        if any(vec_K2.<tmp_floor) || !isreal(vec_K2)
                            vec_K2[ (vec_K2.<0.0) .| (vec_K2.!=real(vec_K2))  ] = tmp_floor
                        end
                    # 6.3.8 update K
                        vec_K[IDX_T] = vec_K[IDX_T] .+ UPDATE_STEP .* (vec_K2[IDX_T] .- vec_K[IDX_T])
                    # 6.3.9 Diagnose (extra)
                        if flag_diagnose
                            println("\t\t","- K updated: ", Profile_vec(vec_K) )
                            println("\t\t","- L updated: ", Profile_vec(vec_L) )
                        end
                    # 6.3.10 Output snapshot of current loop (for monitoring & extra diagnoses)
                        RESULT_mat = [mat_a, mat_Phi, mat_c, mat_l, mat_M, mat_MA, mat_MB, mat_wprofiled, mat_Lambda]
                        RESULT_vec = [vec_Q, vec_D, vec_D2Y, vec_G, vec_I, vec_K, vec_L, vec_r, vec_wmean, vec_TRc, vec_TRw, vec_gapUMP, vec_U, vec_Y, vec_oCoef]
                        Output_Transition(RESULT_mat,RESULT_vec, path=".\\output\\", prefix="Trans")

                end # branch (convergence check) ends

        end # loop ends

        return nothing
    end
    # ==========================================================================



















# ========================== MOD END ===========================================
end
# ========================== MOD END ===========================================
#
