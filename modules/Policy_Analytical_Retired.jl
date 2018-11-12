__precompile__()
"""
    Policy_Analytical_Retired

a modified module based on Policy_Analytical;
solves lifetime decisions for those start from their retired phase (usually those who are still live in year 1 when shock comes);
Please refer to "New_DepartableUtilityFunction.pdf" for mathematical contents;

The programmes were coded in Julia-0.6x version;
some deprecated features (like linspace()) should be modified in Julia 1.0 and later;
I marked them by "&&update_require&&";
the readers may use searching functions of their editors to easily locate where to update when they want to run the codes in later Julia versions;

Similar algorithm & treatment & adjustments in Policy_Analytical module;
pls refer to specific documentations of each function


by Tianhao Zhao
2018-9-10 \\ 1st draft
"""
module Policy_Analytical_Retired
    using Policy_Analytical
    using NumAlgo
    export PolicySolve_Analytical_Retired, AmmoReload_Retired_DATA, AmmoReload_Retired_PARS
# ==============================================================================
    """
        PolicySolve_Analytical_Retired(A0::Float64,Phi0::Float64,
            MAX_AGE::Int64 , START_AGE::Int64,
            DATA::Dict{String,Array{Float64,1}},
            PARS::Dict{String,Float64}    )

    Requires:
        1. A0: initial value of personal asset account when age s=START_AGE (beginning of the year)
        2. Phi0: initial value of individual medical account when age s=START_AGE (beginning of the year)
        3. MAX_AGE: maximum age limit (to death)
        4. START_AGE: age to start
        5. DATA, PARS: datasets of parameters and economic state variables (like interest rates and wage levels)
    Returns:
        1. apath::Array{Float64,1} [len=MAX_AGE]: path of personal asset
        1. Phipath::Array{Float64,1} [len=MAX_AGE]: path of individual medical account
        1. Cpath::Array{Float64,1} [len=MAX_AGE]: path of consumption

    There are what should be contained in DATA, DATA_w, DATA_r & PARS [format: VarName | DictKey : Definition]
    1. DATA [all elements have length MAX_AGE-START_AGE+1]
        1. Rpath | "Rpath": interest rates in each age year
        2. Fpath | "Fpath": mortality in each age year, asking for the last element (s=S) to be 0
        3. Qpath | "Qpath": m2c ratio (medical expenses / consumption) in each year
        4. Ppath | "Ppath": outpatient / inpatient expenses
        5. cpBpath | "cpBpath": inpatient copayment rates
        6. Lambdapath | "Lambdapath": pension benefits
        7. doublePpath | "doublePpath": transfer (amount) from firm contribution to those retired
    2. PARS [all elements are scalars]
        1. Mu | "Mu": consumption tax rate
        4. Gamma | "Gamma": inter-temporal substitution elasticity
        5. Delta | "Delta": utility discount rate

    We do not do length check inside the PolicySolve_Analytical();
    the check has been performed when using "AmmoReload_Retired_x" functions to construct the inputs;
    And, most validations of Level 0 data (original inputs) should also be pre-validated in "AmmoReload_Retired_x" functions;

    Depends extra:
        1. Bisection() in module: NumAlgo
    """
    function PolicySolve_Analytical_Retired(A0::Float64,Phi0::Float64,
        MAX_AGE::Int64 , START_AGE::Int64,
        DATA::Dict{String,Array{Float64,1}},
        PARS::Dict{String,Float64}    )
        # ---------------------------------------------------------------------
        # Section 0: Unload data & parameters from data collections (Level 0)
        # ---------------------------------------------------------------------
        # 0.1 re-define MAX_AGE by START_AGE
            MAX_AGE = MAX_AGE - START_AGE + 1
        # 0.1 Unload those len=MAX_AGE (work in the whole lifetime)
            Rpath = DATA["Rpath"] # interest rate
            Fpath = DATA["Fpath"] # mortality
            Qpath = DATA["Qpath"] # m2c ratio
            Ppath = DATA["Ppath"] # outpatient / inpatient
            cpBpath = DATA["cpBpath"] # inpatient copayment rate
            Lambdapath = DATA["Lambdapath"] # pension benefits (amount of money)
            doublePpath = DATA["doublePpath"] # transfer (amount) from firm contribution to those retired
        # 0.4 Unload those scalar parameters
            Mu = PARS["Mu"] # consumption tax rate
            Gamma = PARS["Gamma"] # inter-temporal substitution elasticity
            Delta = PARS["Delta"] # utility discount rate
        # 0.5 Prepare result arrays
            scriptA = zeros(MAX_AGE)  # total wealth path, = apath + Phipath
            apath = zeros(MAX_AGE)  # personal asset path
            Phipath = zeros(MAX_AGE)  # individual medical account path
            Cpath = zeros(MAX_AGE)  # consumption path
        # 0.5.1 initialization (essential, used in un-bounded solution)
            scriptA[1] = A0 + Phi0

        # ---------------------------------------------------------------------
        # Section 1: Special Validations of Level 0 data/short-writings
        # ---------------------------------------------------------------------
        # 1.1 Interest Rate cannot be -1
            all(Rpath .!= -1.0) || error(ErrorException("interest rates cannot not be -100%, or there will be infinity"))
        # 1.2 Utility discounting rate cannot be -1
            Delta != -1.0 || error(ErrorException("utility discounting rate (Delta) cannot be -100%, or there will be infinity"))
        # 1.3 The mortility of the very last age year should be 0, if not, correct it
            Fpath[end] == 0.0 || begin Fpath[end] == 0.0; warn("the mortality in last age year should always be 0.0 exactly") end

        # ---------------------------------------------------------------------
        # Section 2: Define Level 1 & Level 2 Short-Writings
        # ---------------------------------------------------------------------
        # 2.3 \scripta (Level 2)
            # scripta = 2.0 .- 1 ./ ( 1 .- Fpath )
            scripta = 1.0 .- Fpath
        # 2.5 \scriptd (Level 2)
            scriptd = Qpath .* ( Ppath .+ (1.0 .- cpBpath) ) ./ ( 1.0 .+ Ppath )
        # 2.7 \scriptg (Level 2)
            scriptg = -1.0 .* Qpath .* Ppath ./ ( 1.0 .+ Ppath )
        # 2.8 \scripth (Level 2) = \scriptd + \scriptg
            scripth = Qpath .* ( 1.0 .- cpBpath ) ./ ( 1.0 .+ Ppath )
        # 2.9 \scriptj (Level 2)
            scriptj = Lambdapath .+ doublePpath
        # 2.10 Discounting Factor (Money/Capital/Wealth/Asset)
            # NOTE: Julia supports any-level-precision floating computation, thus we can just write cumprod(), or (if you try to transplant the codes to other platforms) please seriously consider logarithm transformations to avoid accumulating floating errors (esp. in C++ and Python)
            # NOTE: to discount all cash flows to age 0 (year -1), not common but acceptable & convenient in mathematics also coding; no affact in path computation
            V = cumprod( 1.0 ./ ( 1.0 .+ Rpath ) )
        # 2.11 Mortality-adjusted Discounting Factor of cash flows
            # NOTE: we consider accidentaly deaths during the lifetime, which means a \scripta re-scaling multiplier for bequests in EACH year; all cash flows are divided by \scripta to there real values (after obtaining bequests); however, we found that the bequest effect, the multiplier, can be uniformly integrated into "V" we just defined, the discounting factor; it provides much convenience in mathematics & may has such a conception of "mortality-adjusted discounting factor"; more information, pls refer to our documentations and/or paper
            # NOTE: in practice, we overweite "V" rather than define a new V_adjusted; cauz "V" will never be solely used
            V = V ./ scripta
        # 2.11 Discounting Factor (Utility, considering mortalities)
            betatilde = ( 1.0 ./ (1.0 + Delta) ) .^ ( 0:(MAX_AGE-1) )
            betatilde = ( 1.0 .- Fpath ) .* betatilde
        # 2.12 Essential Validations
            # 2.12.1 all \scripth should be in the open range (0,1)
                # NOTE: or there will be numerical collapses in consumptions; however, in fact, according to the definition of \scripth, we have secured the condition in previous validations in "AmmoReload_DATA_w()" function; if much worse performance here, consider ignore/comment this validation process
                all( 0.0 .< scripth .< 1.0 ) || error(ErrorException("\scripth should be in the open range (0,1), or there will be numerical collapse of consumption path"))

        # ---------------------------------------------------------------------
        # Section 2.5: Check the case: MAX_AGE - START_AGE + 1 = 1
        # ---------------------------------------------------------------------
            # NOTE: if MAX_AGE-START_AGE+1==1, no need to do the following processes,
            # just use the inter-temporal budget constraint to get a single "c" then return
        if MAX_AGE == 1
            # a. use the inter-temporal budget in the last year to get consumption
                Cpath[1] = (1.0 + Rpath[1])*scriptA[1] + scriptj[1]
                Cpath[1] /= 1.0 - scripth[1]
            # b. return
                return [A0], [Phi0], Cpath
        end

        # ---------------------------------------------------------------------
        # Section 3: Define Level 3 Short-Writings (Mainly for Dynamics)
        # ---------------------------------------------------------------------
        # 3.0 Prepare empty vectors for dynamic-equation short-writings
            # 3.0.1 Dynamic component \scriptP{s,s+1}, s=1:S-1 in Euler Equation
            scriptP = zeros(MAX_AGE-1)
            # 3.0.2 Dynamic component \scriptQ{s,s+1}, s=1:S-1 in Euler Equation
            scriptQ = zeros(MAX_AGE-1)
        # 3.2 Filling the empty arrays
            for s in 1:MAX_AGE-1
                # NOTE: from s to s+1, locating at index = s
                # 3.2.1
                    scriptP[s] = (1+Rpath[s+1])/(1+Delta)
                    scriptP[s] *= (1-scripth[s])/(1-scripth[s+1])
                # 3.2.2
                    scriptQ[s] = (1-Qpath[s])/(1-Qpath[s+1])
            end
        # 3.3 Essential Valiations (IMPORTANT!) to make sure a bounded (c>0) consumption path
            all(scriptP .> 0) || error(ErrorException("invalid scriptP which leads to invalid consumption path"))
            all(scriptQ .> 0) || error(ErrorException("invalid scriptQ which leads to invalid consumption path"))

        # ---------------------------------------------------------------------
        # Section 4: Define Level 4 Short-Writings (Mainly for Accumulated Dynamics)
        # ---------------------------------------------------------------------
        # 4.1 \scriptT{1->s} [len=MAX_AGE][s=1:MAX_AGE]
            # NOTE: in the mathematics documentation, we define \scriptT{1->s} as the multiplier on initial wealth \scritpA{s=1}; it means: \scriptA{s} = \scriptT{1->s} * \scriptA{1}; naturally, the \scriptT{1->s} series should have MAX_AGE-1 elements, and works on \scriptA{2} to \scriptA{MAX_AGE}; however, for a clear & convenient & easy-understanding style, we define a MAX_AGE elements array here for \scriptT and set its first element as 1.0 (which means \scriptT{1->1}=1.0); thus, now we can UNIFORMLY write the accumulated dynamic relationship as \scriptA{s} = \scriptT{1->s} * \scriptA{1}; be careful & mind-clear when reading the documentation^& coding meantime;
            scriptT = ones(MAX_AGE)
            for s in 1:MAX_AGE-1
                scriptT[s+1] = scriptP[s]^Gamma * scriptQ[s]^(1-Gamma)
            end
            scriptT = cumprod(scriptT)
        # 4.2 Essential Validation to make sure bounded consumption (c>0)
            all( scriptT .> 0.0 ) || error(ErrorException("invalid scriptT which leads to invalid consumption path"))

        # ---------------------------------------------------------------------
        # Section 5: Define Level 5 Short-Writings (for unbounded c,l solving)
        # ---------------------------------------------------------------------
        # 5.1 \scriptX{1->MAX_AGE} (the left side of final equation)
            scriptX  = sum( V .* (1.0 .- scripth) .* scriptT )
        # 5.2 \scriptY{1->MAX_AGE} (the right side of final equation)
            scriptY  = V[1] * scriptA[1]
            scriptY += sum( V .* scriptj )
        # 5.3 Essential Validations to make sure consumptions finite & valid
            scriptX != 0.0  || error(ErrorException("zero \scriptX leads to infinite consumptions"))

        # ---------------------------------------------------------------------
        # Section 6: Get un-bounded solution and check whether the bounds are met
        # ---------------------------------------------------------------------
        # 6.1 Get unbounded solutions
            # 6.1.1 Get the first year consumption
                Cpath[1] = scriptY ./ scriptX
                Cpath[1] *= (1+Rpath[1])  # discount to s=1
            # 6.1.2 Use Euler Equation to generalize the path
                for s in 2:MAX_AGE; Cpath[s] = Cpath[1] .* scriptT[s]; end
        # 6.2 Check whether all consumptions meets the bounding
            all( 0.0 .<= Cpath .< Inf ) || error(ErrorException("negative or Inf consumption found in un-bounded solution"))
        # 6.4 Check whether the leisure path meets bounds
            flag_AdjustNeeded = false

        # ---------------------------------------------------------------------
        # Section 8: Get apath & Phipath
        # ---------------------------------------------------------------------
        # 8.1 Get personal asset path (apath)
            apath[1:end],apath_dead = GetFullPath_Retired_a( A0,
                    Cpath, Rpath, scriptd, scripta, Lambdapath   )

        # 8.2 Get individual medical account path (Phipath)
            Phipath[1:end],Phipath_dead,Gaps = GetFullPath_Retired_Phi( Phi0,
                    Cpath, Rpath, scriptg, scripta, doublePpath   )
            # gap covered by personal asset (apath)
            apath = apath .+ Gaps

        # 8.3 Get total wealth path (scriptA)
            # NOTE: scriptA_dead indicates the remaining wealth/capital when agent dies (at the end of the last age year)  # scriptA[1] = A0 + Phi0
            scriptA[1:end],scriptA_dead = GetFullPath_Retired_scriptA( A0 + Phi0, Cpath,
                        Rpath, scripth, scripta, scriptj  )

        # 8.4 Check whether the relationship "a+Phi=scriptA" met
            findmax( abs.( scriptA .- apath .- Phipath ) )[1] < (1E-06) || error(ErrorException("the relationship scriptA = a + Phi not met"))

        # ---------------------------------------------------------------------
        # Section 9: Search a c_1 which meets scriptA_dead == 0 (if scriptA_dead != 0)
        # ---------------------------------------------------------------------
        # 9.1 Tolerance to believe scriptA_dead = 0
            TOL = 1E-08
        # 9.2 Use bi-sectional method to search the zero point of scriptA_dead
            if abs( scriptA_dead - 0.0 ) > TOL
                # 9.2.0 Define a temp function to use c_1 to compute scriptA_dead (opjevt function)
                    objfunc(C0::Float64) = begin
                        # a. construct a Cpath by C0 input
                            tmp_Cpath = C0 * scriptT
                        # b. get the wealth when dead, then return it
                            tmp_scriptA,tmp_scriptA_dead = GetFullPath_Retired_scriptA( A0 + Phi0, tmp_Cpath, Rpath, scripth, scripta, scriptj  )
                            return tmp_scriptA_dead::Float64, tmp_scriptA::Array{Float64,1}
                    end
                # 9.2.1 An empty consumption path for searching
                    Cpath2 = copy(Cpath);
                # 9.2.2 CASE 1: scriptA_dead < 0 (over-consumption)
                    if scriptA_dead < 0.0
                        # 9.2.2.1 initialization
                            # NOTE: we've learn that when c_1 = 0, the scriptA_dead MUST be non-negative
                            C1 = 0.0; C2 = Cpath[1]
                        # 9.2.2.2 Bisection searching to find an initial c_1 which makes scriptA_dead == 0.0
                            # NOTE: Get solved initial consumption (s=1) (considering the namespace)
                            Cpath[1] = Bisection(objfunc,C1,C2,MAXLOOP=5000,TOL=1E-12)
                    elseif scriptA_dead > 0.0
                # 9.2.3 CASE 2: scriptA_dead > 0 (too-low-consumption)
                        # 9.2.3.1 initialization
                            # NOTE: we dont have a specific point that makes scriptA_dead < 0, so we search for one (cauz we know as long as consume more and more, we will finally get such a point)
                            # NOTE: (in chinese) 因为在给定一条收入时，死亡时刻资产scriptA_dead对初始消费c_1单调递减且凹（收入流不变而增加支出肯定会令死亡时刻资产减少），而我们又知道给定各种外生数据后几乎可以保证如果c_1=0那么scriptA_dead一定大于0
                            # NOTE: 显然可知目前的c_1仍然处在零点左边，那么我们要想办法找到任意一个（但最好是离零点很近）能够保证是零点右边的点C2就能用二分法搜索零点了
                            # NOTE: 那么为什么不成倍放大现在的c_1，以快速找到这样一个点呢？（凹函数下降速度加快）
                            C1 = Cpath[1]; C2 = 0.0 ; tmp_times = 2.0
                            while tmp_times > 0
                                C2 = tmp_times * C1;
                                if objfunc(C2)[1] < 0.0; break; else; tmp_times += 1.0; continue; end
                            end
                        # 9.2.3.2 Looping
                            Cpath[1] = Bisection(objfunc,C1,C2,MAXLOOP=5000,TOL=1E-12)
                    end  # branch ends (CASES)

                # 9.2.4 If sovled successfully
                    Cpath[1:end] = Cpath[1] .* scriptT
                    # get scriptA_dead and the whole path
                    scriptA[1:end],scriptA_dead = GetFullPath_Retired_scriptA( A0 + Phi0, Cpath, Rpath, scripth, scripta, scriptj  )
                    # get Phipath
                    Phipath[1:end],Phipath_dead,Gaps = GetFullPath_Retired_Phi( Phi0, Cpath, Rpath, scriptg, scripta, doublePpath   )
                    # get apath
                    apath[1:end] = scriptA .- Phipath
                # 9.2.5 Check scriptA >=0
                    # all(scriptA .>= 0.0) || error(ErrorException("Fail to solve a non-negative Wealth path, please try less generations!"))
                    abs(scriptA_dead) < 1E-3 || error(ErrorException(string("Fail to solve a non-negative Wealth path, please try less generations! ",scriptA_dead)  ))
                # 9.2.6 Finally, return
                    return apath,Phipath,Cpath
            end  # branch/search ends

        # Check scriptA >=0
            # all(scriptA .>= 0.0) || error(ErrorException("Fail to solve a non-negative Wealth path, please try less generations!"))
            abs(scriptA_dead) < 1E-3 || error(ErrorException(string("Fail to solve a non-negative Wealth path, please try less generations! ",scriptA_dead)  ))
        # Finally, return the results
            return apath,Phipath,Cpath

    end  # function ends


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """
    AmmoReload_Retired_DATA( MAX_AGE::Int64, START_AGE::Int64,
            Rpath::Array{Float64,1}, Fpath::Array{Float64,1},
            Qpath::Array{Float64,1}, Ppath::Array{Float64,1}, cpBpath::Array{Float64,1},
            Lambdapath::Array{Float64,1}, doublePpath::Array{Float64,1}  )

    Returns a dictionary DATA; also do validations before packaging
    """
    function AmmoReload_Retired_DATA( MAX_AGE::Int64, START_AGE::Int64,
        Rpath::Array{Float64,1}, Fpath::Array{Float64,1},
        Qpath::Array{Float64,1}, Ppath::Array{Float64,1}, cpBpath::Array{Float64,1},
        Lambdapath::Array{Float64,1}, doublePpath::Array{Float64,1}  )
        # --------------------------
        # 1. Validations
            # 1.0 re-define MAX_AGE
                MAX_AGE = MAX_AGE - START_AGE + 1
            # 1.1 Length Check
                @assert(MAX_AGE==length(Rpath)==length(Fpath)==length(Qpath)==length(Ppath)==length(cpBpath)==length(Lambdapath)==length(doublePpath),"Uncompatible MAX_AGE-length inputs")
            # 1.2 Range Check
                all( Rpath .!= -1.0 ) || error(ErrorException("interest rates cannot be equal to -100%, or there will be all zeros assets"))
                all( 0.0 .<= Fpath .< 1.0 ) || error(ErrorException("mortalities should be in the half open range [0,1)"))
                all( 0.0 .< Qpath .< 1.0 ) || error(ErrorException("m2c ratios should be in the open range (0,1)"))
                all( 0.0 .< Ppath ) || error(ErrorException("outpatient/inpatient expenses ratio should be in greater than 0"))
                all( 0.0 .<= cpBpath .<= 1.0 ) || error(ErrorException("inpatient copayment rates should be in the closed range [0,1]"))
                all( Lambdapath .>= 0.0 ) || error(ErrorException("pension benefits should be greater than or equal to 0"))
                all( doublePpath .>= 0.0 ) || error(ErrorException("medical transfer to the retired should be greater than or equal to 0"))
            # 1.3 Special Check
                Fpath[end] == 0.0 || begin Fpath[end] = 0.0; warn("the last element of Fpath should be exact 0.0") end
        # 2. Packaging
            DATA = Dict(
                "Rpath" => Rpath,
                "Fpath" => Fpath,
                "Qpath" => Qpath,
                "Ppath" => Ppath,
                "cpBpath" => cpBpath,
                "Lambdapath" => Lambdapath,
                "doublePpath" => doublePpath
            )
        # 3. Return
        return DATA::Dict{String,Array{Float64,1}}
    end  # function ends

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """
        AmmoReload_Retired_PARS( Mu::Float64, Gamma::Float64, Delta::Float64  )

    Returns a dictionary PARS; also do validations before packaging
    """
    function AmmoReload_Retired_PARS( Mu::Float64, Gamma::Float64, Delta::Float64  )
        # ----------------------
        # 1. Range Check
            Delta != -1.0 || error(ErrorException("Delta cannot be -1, or there will be infinite utility"))
            (0.0 < Gamma < 1.0) || error(ErrorException("Gamma should be in the open range (0,1)"))
        # 2. Packaging
            PARS = Dict(
                "Delta" => Delta,
                "Gamma" => Gamma,
                "Mu" => Mu,
            )
        # 3. Return
        return PARS::Dict{String,Float64}
    end  # function ends


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """
        GetFullPath_Retired_scriptA( scriptA0::Float64,
            Cpath::Array{Float64,1},
            Rpath::Array{Float64,1}, scripth::Array{Float64,1}, scripta::Array{Float64,1}, scriptj::Array{Float64,1}   )

    a subfunction, computes whole WEALTH path from s=START_AGE to s=S;
    returns a vector tmp_scriptA (copy, not reference), also a number scriptA_dead,
    which indicates the wealth when dead (used in bisection search);
    and, no validation cauz the function is only used inside the PolicySolve_Analytical_Retired();
    """
    function GetFullPath_Retired_scriptA( scriptA0::Float64,
        Cpath::Array{Float64,1},
        Rpath::Array{Float64,1}, scripth::Array{Float64,1}, scripta::Array{Float64,1}, scriptj::Array{Float64,1}   )
        # -----------------
        # 0. prepare empty vector for wealth (scriptA)
            MAX_AGE = length(Cpath)
            tmp_scriptA = zeros(MAX_AGE)
            tmp_scriptA[1] = scriptA0
        # 2. post-retire
            # NOTE: the loop is disigned for case: MAX_AGE-START_AGE>1, if MAX_AGE-START_AGE=1, only need to run section 3
            if MAX_AGE>1
                for s in 1:MAX_AGE-1
                    tmp_scriptA[s+1] = (1+Rpath[s])*tmp_scriptA[s] + scriptj[s] - (1-scripth[s])*Cpath[s]
                    tmp_scriptA[s+1] /= scripta[s]
                end
            end
        # 3. wealth at death
            tmp_scriptA_dead = (1+Rpath[end])*tmp_scriptA[end] + scriptj[end] - (1-scripth[end])*Cpath[end]
        # 4. returns
            return tmp_scriptA::Array{Float64,1}, tmp_scriptA_dead::Float64
    end  # function ends


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """
        GetFullPath_Retired_Phi( Phi0::Float64,
            Cpath::Array{Float64,1},
            Rpath::Array{Float64,1}, scriptg::Array{Float64,1}, scripta::Array{Float64,1},
            doublePpath::Array{Float64,1}   )

    a subfunction, computes whole INDIVIDUAL MEDICAL ACCOUNT path from s=START_AGE to s=S;
    returns a vector tmp_Phipath (copy, not reference), also a number tmp_Phipath_dead, and a path of g ap covered by personal asset (apath) where negative for gaps;
    which indicates the individual medical account when dead (used in bisection search);
    and, no validation cauz the function is only used inside the PolicySolve_Analytical_Retired();
    """
    function GetFullPath_Retired_Phi( Phi0::Float64,
        Cpath::Array{Float64,1},
        Rpath::Array{Float64,1}, scriptg::Array{Float64,1}, scripta::Array{Float64,1},
        doublePpath::Array{Float64,1}   )
        # -----------------
        # 0. prepare empty vector for wealth (scriptA)
            MAX_AGE = length(Cpath)
            tmp_Phipath = zeros(MAX_AGE)
            tmp_Phipath[1] = Phi0
        # 2. post-retire
            # NOTE: the loop is disigned for case: MAX_AGE-START_AGE>1, if MAX_AGE-START_AGE=1, only need to run section 3
            if MAX_AGE>1
                for s in 1:MAX_AGE-1
                    tmp_Phipath[s+1] = (1+Rpath[s])*tmp_Phipath[s] + doublePpath[s] + scriptg[s]*Cpath[s]
                    tmp_Phipath[s+1] /= scripta[s]
                end
            end
        # 3. wealth at death
            tmp_Phipath_dead = (1+Rpath[end])*tmp_Phipath[end] + doublePpath[end] + scriptg[end]*Cpath[end]
        # 4. Cover gaps by personal Aaset path (apath)
            Gaps = zeros(MAX_AGE)
            Gaps[ tmp_Phipath .< 0.0 ] = tmp_Phipath[ tmp_Phipath .< 0.0 ]
            tmp_Phipath[ tmp_Phipath .< 0.0 ] = 0.0
        # 5. returns
            return tmp_Phipath::Array{Float64,1}, tmp_Phipath_dead::Float64, Gaps::Array{Float64,1}
    end  # function ends

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """
        GetFullPath_a( A0::Float64,
            Cpath::Array{Float64,1}, Lpath::Array{Float64,1},
            Rpath::Array{Float64,1}, scriptd::Array{Float64,1}, scripta::Array{Float64,1},
            Wpath::Array{Float64,1}, scriptb::Array{Float64,1},
            Lambdapath::Array{Float64,1}   )

    a subfunction, computes whole PERSONAL ASSET path from s=1 to s=S;
    returns a vector tmp_apath (copy, not reference), also a number tmp_apath_dead;
    which indicates the personal asset when dead (used in bisection search);
    and, no validation cauz the function is only used inside the PolicySolve_Analytical();
    And, no gap from individual medical account (Phipath) considered;
    """
    function GetFullPath_Retired_a( A0::Float64,
        Cpath::Array{Float64,1},
        Rpath::Array{Float64,1}, scriptd::Array{Float64,1}, scripta::Array{Float64,1},
        Lambdapath::Array{Float64,1}   )
        # -----------------
        # 0. prepare empty vector for wealth (scriptA)
            MAX_AGE = length(Cpath)
            tmp_apath = zeros(MAX_AGE)
            tmp_apath[1] = A0
        # 2. post-retire
            # NOTE: the loop is disigned for case: MAX_AGE-START_AGE>1, if MAX_AGE-START_AGE=1, only need to run section 3
            if MAX_AGE>1
                for s in 1:MAX_AGE-1
                    tmp_apath[s+1] = (1+Rpath[s])*tmp_apath[s] + Lambdapath[s] - (1 - scriptd[s])*Cpath[s]
                    tmp_apath[s+1] /= scripta[s]
                end
            end
        # 3. wealth at death
            tmp_apath_dead = (1+Rpath[end])*tmp_apath[end] + Lambdapath[end] - (1 - scriptd[end])*Cpath[end]
        # 5. returns
            return tmp_apath::Array{Float64,1}, tmp_apath_dead::Float64
    end  # function ends




# ==============================================================================
end  # module ends
#
