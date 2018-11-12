__precompile__()
"""
    SteadyState

Searches steady states based on given data.
"""
module SteadyState
    using NumAlgo
    using Policy_Analytical
    # using PolicyFunctions  # DP, deprecated

    export SStateSearch
# ========================== MOD BEGIN =========================================
    """
        SStateSearch(Rguess::Float64,Lguess::Float64,Qguess::Float64, MAX_AGE::Int64, RETIRE_AGE::Int64,
            DATA_vec::Array{Array{Float64,1},1}, DATA_num::Array{Float64,1} ;
            flag_diagnose::Bool = true, MagicNum::Float64 = 2.0, Arange::Array{Float64,1} = [0.0,25.0],
            TOL_SS::Float64 = 1E-05, TOL_POLICY::Float64 = 1E-08, GRID_LEN::Int64 = 50, MAX_ITER::Int64 = 200, , UPDATE_STEP::Float64 = 0.3 )

    Searches steady state, using Gauss-Seidel algorithm.

    ##Inputs:
    1. Rguess: a guess of interest rate to compute capital in the first round (for the purpose of faster convergence)
    1. Lguess: a guess of aggregated labour supply, in [0,1]
    1. Qguess: a CONSTANT level of the ratio of medical fee (outpatient + inpatient) to consumption
    1. MAX_AGE: maximum age
    1. RETIRE_AGE: retire age
    1. DATA_vec: a DATA PACKAGE of vector-like data, in a nested one-dimension array, in order
    1. DATA_num: a DATA PACKAGE of scalar data, in a one-dimension array, in order
    1. flag_diagnose: controls whether to print a short summary for each round
    1. MagicNum: a trick to control the algorithm, as the lowest level of K/L
    1. Arange: the range of policy function searching
    1. TOL_SS: relative tolerance of Gauss-Seidel algorithm, for convergence check
    1. TOL_POLICY: relative tolerance for the Golden Section searching in policy function solving (DP)
    1. GRID_LEN: the density of griding search in policy function solving (DP)
    1. MAX_ITER: maximum iteration times of Gauss-Seidel algorithm
    1. UPDATE_STEP: relative step size to update K & L, in (0,1], similar to the category of "belief"

    ##Outputs:
    1. RESULT_vec: a nested one-D array containing vector steady state variables (in order)
    1. RESULT_num: a nested one-D array containing scalar steady state variables (in order)

    ##**DATA PACKAGES**
    1. DATA_vec [len=4]
        1. N [len=MAX_AGE] : demography in steady state, distributes on age
        1. Survival [len=MAX_AGE] : survival prob by age
        1. MA2MB [len=MAX_AGE] : outpatient/inpatient ratio by age
        1. Epsilon [len=RETIRE_AGE] : wage profile coefficients
    1. DATA_num [len=18]
        1. Kappa : depreciation rate of capital
        1. Gamma : inter-temporal substitution elasticity
        1. Varrho : consumption elasticity of leisure
        1. Delta : discount rate of utility
        1. Alpha : preference coefficient of leisure than consumption
        1. Mu : consumption tax rate
        1. Sigma : wage tax rate
        1. Beta : capital income share
        1. Tech : TFP
        1. doubleK : cap of the ratio of D/Y
        1. Eta : firm contribution to pension
        1. Theta : personal contribution to pension
        1. Zeta : firm contribution to medical
        1. phiCoef : personal contribution to medical
        1. doubleA : transfer rate from firm contribution to individual medical account (working phase)
        1. doubleB : transfer rate from firm contribution to individual medical account (retired phase)
        1. cpB : copayment rate of inpatient fee
        1. z : collection rate of pension
    1. RESULT_vec [len=9]
        1. wprofiled [len=RETIRE_AGE] : profiled wage level
        1. Lambda [len=MAX_AGE-RETIRE_AGE] : pension benefits
        1. apath [len=MAX_AGE] : personal asset path
        1. Phipath [len=MAX_AGE] : individual medical account path
        1. Cpath [len=MAX_AGE] : consumption path
        1. Lpath [len=RETIRE_AGE] : anti-retire leisure path
        1. M [len=MAX_AGE] : total medical expense path
        1. MA [len=MAX_AGE] : outpatient expense path
        1. MB [len=MAX_AGE] : inpatient expense path
    1. RESULT_num [len=12]
        1. Final_K [num]: calibrated capital
        1. Final_L [num]: calibrated labour
        1. Y [num]: GDP
        1. r [num]: interest rate (net)
        1. wmean [num]: social average wage level
        1. oCoef [num]: wage scaling coefficient
        1. TRc [num]: total consumption tax revenue
        1. TRw [num]: total wage tax revenue
        1. gapUMP [num] : gap of the pooling medical account (positive for gap, negative for surplus)
        1. D [num] : government outstanding debt
        1. D2Y [num] : ratio of debt to GDP
        1. G [num] : government purchase
        1. I [num] : investment
        1. doubleP [num] : transfer AMOUNT from firm contribution to individual medical account (retired phase)
    """
    function SStateSearch(Rguess::Float64,Lguess::Float64,Qguess::Float64, MAX_AGE::Int64, RETIRE_AGE::Int64,
        DATA_vec::Array{Array{Float64,1},1}, DATA_num::Array{Float64,1} ;
        flag_diagnose::Bool = true, MagicNum::Float64 = 2.0,
        TOL_SS::Float64 = 1E-05, MAX_ITER::Int64 = 200, UPDATE_STEP::Float64 = 0.3 )
        # -------------------------
        # Timing
        time_begin = time()
        # ----------------------------------------------------------------------
        # Section 1: DATA PKG VALIDATION & PROCESS 数据拆包与预处理
        # ----------------------------------------------------------------------
        # 1.1 Length Check of input data packages
            @assert( (length(DATA_vec)==4)&&(length(DATA_num)==18) , "Uncompatible DATA packages")
        # 1.2 un-package data
            N, Survival, MA2MB, Epsilon = copy(DATA_vec)
            Kappa, Gamma, Varrho, Delta, Alpha, Mu, Sigma, Beta, Tech, doubleK, Eta, Theta, z, Zeta, phiCoef, doubleA, doubleB, cpB = copy(DATA_num)
            # essential validation for survival probs
            Survival[end] == 1.0 || begin Survival[end] = 1; warn("the survival prob in last age year should be 1.0 exactly"); end
        # 1.3 prepare memories for outputs 预留结果变量空间
            wprofiled, Lambda, apath, Phipath, Cpath, Lpath, M, MA, MB, Final_K, Final_L, Y, r, wmean, oCoef, TRc, TRw, gapUMP, D, D2Y, G, I, doubleP = Result_Spacing(MAX_AGE,RETIRE_AGE)
        # 1.4 data pre-process: fill initial leisure matrix with average 数据预处理（填充初始休闲分布）
            Lpath = ( (1-Lguess)/RETIRE_AGE)*ones(RETIRE_AGE)   # leisure path
        # 1.5 other variables in looping 迭代用变量
            Kguess=0.0::Float64  # capital (computed in the first round)
            Kguess2=0.0  # (re-aggregated capital)
            Lguess2=0.0  # (re-aggregated labour)
        # 1.6 short-writings of contribution rates
            # 1.6.1 total pension contribution rates
            Pi = z*(Theta+Eta)/(1+z*Eta+Zeta)
            # 1.6.2 total medical contribution rates
            PiM = (phiCoef+Zeta)/(1+z*Eta+Zeta)

        # ----------------------------------------------------------------------
        # Section 2: LOOPING 迭代搜索 designed Gauss-Seidel algorithm
        # ----------------------------------------------------------------------
        for idx_loop in 1:MAX_ITER
        # 1. printing info
            if flag_diagnose
                println("Round: ",idx_loop)
            end
        # 2. Validation of capital
            @assert( 0.0<=Lguess , "negative Labour!"  )

        # 3. Firm Department 厂商部门 -------------------------------------------
            # 3.1 interest rate, wage and GDP (output)
                if idx_loop == 1
                    # NOTE: for the 1st round, we use interest rate & labour to compute the guess of capital
                    r = Rguess
                    Kguess = Lguess*( (r+Kappa)/(Beta*Tech) ).^(1/(Beta-1))
                    @assert(Kguess>0,"negative K!")
                    wmean = (Rguess+Kappa)*(1/Beta-1)*Kguess/Lguess
                    Y = Tech.*(Kguess).^(Beta).*(Lguess).^(1-Beta)
                else
                    # NOTE: for the other loops, we directly start from K, L
                    Y = Tech.*(Kguess).^(Beta).*(Lguess).^(1-Beta)
                    r = Beta*Tech*(Kguess/Lguess)^(Beta-1) - Kappa
                    wmean = (1-Beta)*Tech*(Kguess/Lguess)^Beta
                end
            # 3.2 wage profiling
                # 3.2.1 a scaling factor for diff ages of agents
                oCoef = Lguess / sum(N[1:RETIRE_AGE].*Epsilon.*(1-Lpath) )
                # 3.2.2 profiling
                wprofiled = wmean*oCoef*Epsilon
            # 3.3 PAYG pension benefits
                # 3.3.1 total contribution from working people / total retired population
                Lambda = sum( Pi*wprofiled.*(1-Lpath[1:RETIRE_AGE]).*N[1:RETIRE_AGE] ) / sum( N[RETIRE_AGE+1:MAX_AGE] )
                # NOTE: here the benefit distributes even on ages
                Lambda = Lambda*ones(MAX_AGE-RETIRE_AGE)
            # 3.4 transfer level (amount) from firm medical contribution to those retired, pls refer to our paper
                # NOTE: 医保缴纳中对退休人群的转移支付, doubleB是从企业总缴纳部分中转移给当年退休人群的比例
                doubleP = doubleB * sum( Zeta/(1+Eta+Zeta)*wprofiled.*(1-Lpath).*N[1:RETIRE_AGE] )
                doubleP /= sum( N[RETIRE_AGE+1:MAX_AGE] )

        # 4. Household Department 家庭部门 --------------------------------------
            # 4.1 data packaging (using API in module Policy_Analytical, packaged as Dict structures)
                DATA = AmmoReload_DATA( MAX_AGE, r*ones(MAX_AGE), 1-Survival, Qguess*ones(MAX_AGE), MA2MB, cpB*ones(MAX_AGE) )
                DATA_w = AmmoReload_DATA_w( RETIRE_AGE, wprofiled, phiCoef*ones(RETIRE_AGE), Zeta*ones(RETIRE_AGE),
                        Eta*ones(RETIRE_AGE), Theta*ones(RETIRE_AGE), z*ones(RETIRE_AGE), doubleA*ones(RETIRE_AGE)   )
                DATA_r = AmmoReload_DATA_r( MAX_AGE, RETIRE_AGE, Lambda, doubleP*ones(MAX_AGE-RETIRE_AGE)   )
                PARS = AmmoReload_PARS( Mu, Sigma, Alpha, Gamma, Delta  )
            # 4.2 Solving, analytically, calling API in module Policy_Analytical 求解，得到解析解
                apath, Phipath, Cpath, Lpath = PolicySolve_Analytical(0.0, 0.0, MAX_AGE , RETIRE_AGE, DATA, DATA_w, DATA_r, PARS   )
            # 4.3 get medical expenses
                M = Cpath*Qguess
                MA = M.*MA2MB./(1+MA2MB); MB = M./(1+MA2MB)

        # 5. Fiscal & Social Security accounting 财政核算 ----------------------
            # 5.1 Tax Revenue (Consumption & Wage taxes)
                TRc = sum( Cpath.*N ) * Mu
                TRw = sum( wprofiled.*N[1:RETIRE_AGE].*(1-Lpath) ) * Sigma
            # 5.2 Pooling medical account gaps
                # NOTE: positive for gap, negative for surplus (benefits - incomes)
                gapUMP = sum( (1-cpB)*MB.*N ) - (1-doubleA-doubleB)*Zeta/(1+Eta+Zeta) * sum( N[1:RETIRE_AGE].*wprofiled.*(1-Lpath) )

        # 6. Update Labour & Gauss-Seidel style updating 更新劳动力 -------------
            # 6.1 Re-aggregate labour
                Lguess2 = sum( N[1:RETIRE_AGE].*(1-Lpath) )
            # 6.2 Update GDP thorugh production function with re-aggregated Labour factor
                # NOTE: the key characteristic of Gauss-Seidel iteration
                Y = Tech.*(Kguess).^(Beta).*(Lguess2).^(1-Beta)

        # 7. Update Capital 更新资本 --------------------------------------------
            # NOTE: there are two ways to the next loop: one using capital market clearing condition (Method A),
            # and the other using capital growth dynamics (Method B); the two are equivalent, I write them here both,
            # more convenient for the readers when modifying/generating the model design
            # 7.A Method A: updating through CAPITAL MARKET CLEARING CONDITION (资本市场均衡更新法)
                # 7.A.1 get investment I from capital growth dynamics
                    I = Kguess * Kappa
                # 7.A.2 get government purchase from good market clearing condition
                    # NOTE: GDP (Y) now contains new information from re-aggregated labour (6.2)
                    G = Y - I - sum( Cpath.*N )
                # 7.A.3 get government outstanding debt from fiscal budget constraint
                    D = (G + gapUMP - TRw - TRc) / (1-r)
                # 7.A.4 cross-sectional condition for gov-debt (0<=D/Y<=doubleK)
                    D = max( 0.0 , min(  D , Y*doubleK  ) )
                # 7.A.5 record the ratio D/Y
                    D2Y = D / Y
                # 7.A.6 update capital from capital market clearing (re-aggregation)
                    Kguess2 = sum( N.*(apath.+Phipath) ) + D
            # # 7.B Method B: updating through CAPITAL GROWTH DYNAMICS (资本增长动态更新法)
            #     # 7.B.1 get government outstanding debt from capital market clearing condition
            #         D = Kguess - sum(  N.*(apath.+Phipath)  )
            #     # 7.B.2 cross-sectional condition
            #         D = min( D , Y*doubleK )
            #     # 7.B.3 record D/Y ratio (important indicator in analysis)
            #         D2Y = D / Y
            #     # 7.B.4 get government purchase from fiscal budget constraint
            #         G = TRc + TRw + (1-r)*D - gapUMP
            #     # 7.B.5 get investment from good market clearing condition
            #         I = Y - sum(N.*Cpath) - G
            #     # 7.B.6 updating capital through capital growth dynamics
            #         Kguess2 = I / Kappa

            # 8. CONVERGENCE CHECK 收敛检查 -------------------------------------
                # 8.1 get relative errors of capital & labour (in digits)
                    Err = [abs(Kguess2/Kguess-1), abs(Lguess2/Lguess-1)]
                # 8.2 record initial values of current loop (using copy rather than references)
                    Kguess_init = copy(Kguess)
                    Lguess_init = copy(Lguess)
                # 8.3 checking & updating
                    if all(Err .< TOL_SS)   # if all endogeneous variables met the tolerance
                        # 8.3.1 get final K & L to output
                            Final_K = (Kguess+Kguess2)/2; Final_L = (Lguess+Lguess2)/2
                        # 8.3.2 Result data packaging
                            RESULT_vec = [ wprofiled, Lambda, apath, Phipath, Cpath, Lpath, M, MA, MB ]
                            RESULT_num = [ Final_K, Final_L, Y, r, wmean, oCoef, TRc, TRw, gapUMP, D, D2Y, G, I, doubleP ]
                        # 8.3.3 Diagnose printing
                            SummaryPrint("Converged.", MAX_AGE, RETIRE_AGE, N,
                                Kguess, Kguess2, Final_K, Lguess, Lguess2, Final_L, RESULT_vec, RESULT_num )
                        # 8.3.4 return, ends the function
                        return RESULT_vec::Array{Array{Float64,1},1}, RESULT_num::Array{Float64,1}
                    elseif idx_loop == MAX_ITER  # if not have found a solution in MAX_ITER iterations
                        # 8.3.1 get final K & L to output
                            Final_K = (Kguess+Kguess2)/2; Final_L = (Lguess+Lguess2)/2
                        # 8.3.2 Result data packaging
                            RESULT_vec = [ wprofiled, Lambda, apath, Phipath, Cpath, Lpath, M, MA, MB ]
                            RESULT_num = [ Final_K, Final_L, Y, r, wmean, oCoef, TRc, TRw, gapUMP, D, D2Y, G, I, doubleP ]
                        # 8.3.3 Diagnose printing
                            SummaryPrint("Max Iteration Reached.", MAX_AGE, RETIRE_AGE, N,
                                Kguess, Kguess2, Final_K, Lguess, Lguess2, Final_L, RESULT_vec, RESULT_num )
                        # 8.3.4 return, ends the function
                        return RESULT_vec::Array{Array{Float64,1},1}, RESULT_num::Array{Float64,1}
                    else  # if not converged, update guesses then goto next loop
                        # 8.3.1 Bounding check for L
                            (Lguess2<0 || ~isreal(Lguess2)) && begin warn("L floor touched at ",Lguess2); Lguess2 = 0.01; end
                        # 8.3.2 Update L for the next round
                            # NOTE: using UPDATE_STEP to control the relative step length to go
                            Lguess = Lguess + UPDATE_STEP*(Lguess2 - Lguess)
                        # 8.3.3 Set a floor for Kapital thorugh the Magic number
                            # NOTE: essential for convergence;
                            # the MagicNum should be large enough to avoid an IRREGULAR interest rate which may make the algorithm collapse in the error fluctuatings in beginning rounds
                            tmp_floor = MagicNum*Lguess
                        # 8.3.4 Update K for the next round, using floor to bound it
                            (Kguess2<tmp_floor || ~isreal(Kguess2)) && begin warn("K floor touched at ",Kguess2); Kguess2 = tmp_floor; end
                            Kguess = Kguess + UPDATE_STEP*(Kguess2 - Kguess)
                        # 8.3.5 Print reports (optional)
                            if flag_diagnose
                                # Result folding/packaging
                                RESULT_vec = [ wprofiled, Lambda, apath, Phipath, Cpath, Lpath, M, MA, MB ]
                                RESULT_num = [ Kguess, Lguess, Y, r, wmean, oCoef, TRc, TRw, gapUMP, D, D2Y, G, I, doubleP ]
                                SummaryPrint("Solving.", MAX_AGE, RETIRE_AGE, N,
                                    Kguess_init, Kguess2, Kguess, Lguess_init, Lguess2, Lguess, RESULT_vec, RESULT_num )
                            else
                                # if not, only prints a short error report line
                                println("\t+ [K, L] Rel Err: ",Err, " in round ",idx_loop)
                            end
                    end  # branch (convergence check) ends
        end  # loop ends

        # Timing ends
        if flag_diagnose
            println("(SS-Search) Time Elapsed: ",time()-time_begin," s")
        end

        return nothing::Void
    end
    # ==========================================================================
    """
        Result_Spacing(MAX_AGE,RETIRE_AGE)

    returns several empty vectors or 0.0 scalars
    """
    function Result_Spacing(MAX_AGE,RETIRE_AGE)
        wprofiled = zeros(RETIRE_AGE)
        Lambda = zeros(MAX_AGE-RETIRE_AGE)
        apath = zeros(MAX_AGE)
        Phipath = zeros(MAX_AGE)
        Cpath = zeros(MAX_AGE)
        Lpath = zeros(RETIRE_AGE)
        M = zeros(MAX_AGE)
        MA = zeros(MAX_AGE)
        MB = zeros(MAX_AGE)
        # ----------------------
        Final_K = 0.0; Final_L = 0.0; Y = 0.0; r = 0.0; wmean = 0.0; oCoef = 0.0; TRc = 0.0; TRw = 0.0
        gapUMP = 0.0; D = 0.0; D2Y = 0.0; G = 0.0; I = 0.0; doubleP = 0.0
        # -------------------------
        return wprofiled, Lambda, apath, Phipath, Cpath, Lpath, M, MA, MB, Final_K, Final_L, Y, r, wmean, oCoef, TRc, TRw, gapUMP, D, D2Y, G, I, doubleP
    end
    # ==========================================================================
    """
        SummaryPrint(Message::String, MAX_AGE::Int64, RETIRE_AGE::Int64, N::Array{Float64,1},
            Kguess_init::Float64, Kguess_reagg::Float64, Kguess_update::Float64,
            Lguess_init::Float64, Lguess_reagg::Float64, Lguess_update::Float64,
            RESULT_vec::Array{Array{Float64,1},1}, RESULT_num::Array{Float64,1} )

    Prints summary/Diagnoses.
    """
    function SummaryPrint(Message::String, MAX_AGE::Int64, RETIRE_AGE::Int64, N::Array{Float64,1},
        Kguess_init::Float64, Kguess_reagg::Float64, Kguess_update::Float64,
        Lguess_init::Float64, Lguess_reagg::Float64, Lguess_update::Float64,
        RESULT_vec::Array{Array{Float64,1},1}, RESULT_num::Array{Float64,1} )
        # ------------------------
        wprofiled, Lambda, apath, Phipath, Cpath, Lpath, M, MA, MB = RESULT_vec
        Final_K, Final_L, Y, r, wmean, oCoef, TRc, TRw, gapUMP, D, D2Y, G, I, doubleP = RESULT_num
        # ------------------------
        println("\t","-"^60)
        println("\t+ Status : ",Message)
        println("\t+ Capital: ",[Kguess_init,Kguess_reagg,Kguess_update] )
        println("\t+ Labour : ",[Lguess_init,Lguess_reagg,Lguess_update] )
        println("\t* Note   : Initial -> Reagg -> Updated")
        println("\t","."^40)
        println("\t+ GDP    : ",Y)
        println("\t+ Net r  : ",r)
        println("\t+ AvgWage: ",wmean)
        println("\t","."^40)
        println("\t+ CTaxRev: ",TRc)
        println("\t+ WTaxRev: ",TRw)
        println("\t+ D & D/Y: ",[D,D2Y])
        println("\t","."^40)
        println("\t+ Pen-Lev: ",mean(Lambda))
        println("\t+ gapUMP : ",gapUMP)
        # println("\t+ RepRat : ",sum(N[RETIRE_AGE+1:end].*Lambda)/sum(N[1:RETIRE_AGE].*wprofiled.*(1-Lpath))  )
        println("\t","."^40)
        println("\t+ C,I,G  : ",[sum(N.*Cpath),I,G])
        println("\t+ CIG->Y : ",[sum(N.*Cpath)/Y,I/Y,G/Y])
        println("\t","."^40)
        println("\t+ A      : ",sum(N.*apath))
        println("\t+ Phi    : ",sum(N.*Phipath))
        println("\t+ GovDebt: ",D)
        println("\t+ A&Phi&D: ",sum(N.*apath)+sum(N.*Phipath)+D)
        println("\t","."^40)
        println("\t+ Rel Err: ",[abs(Kguess_reagg/Kguess_init-1) , abs(Lguess_reagg/Lguess_init-1) ])
        println("\t* Note   : ","K  ,  L")
        println("\t","-"^40)
        println("\t","* Pen-Lev: Pension Benefits; gapUMP: Gap of Pool Medical Account")
        println("\t","* CIG->Y: ratio of C,I,G to GDP")
        println("\t","* max wealth found: ",findmax(apath+Phipath)[1])

        return nothing::Void
    end








































# ========================== MOD END ===========================================
end
# ========================== MOD END ===========================================
#
