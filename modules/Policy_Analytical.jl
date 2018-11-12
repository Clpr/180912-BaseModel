__precompile__()
"""
    Policy_Analytical

analytically make consumption (c) & leisure (l) decisions for household;
pls refer to "New_DepartableUtilityFunction.pdf" for mathematical contents;

The programmes were coded in Julia-0.6x version;
some deprecated features (like linspace()) should be modified in Julia 1.0 and later;
I marked them by "&&update_require&&";
the readers may use searching functions of their editors to easily locate where to update when they want to run the codes in later Julia versions;

## (The following section about the algorithm structure is written in Chinese (UTF-8)!)

这个算法的主要内容和各符号的约定、公式推导请参见New_DepartableUtilityFunction.pdf；
接下来主要说明的是PolicySolve_Analytical()函数是如何组织的，按照Section列出：
1. PolicySolve_Analytical()
    1. **Section 0**: 从输入参数/数据包中将外生参数、经济系统的各状态变量等取出来，并且定义+初始化存储结果用的结构
    2. **Section 1**: 进行基本的合法性检查（更初级的在AmmoReload_x系列函数构建输入数据包时已经进行了）
    3. **Section 2**: 定义Level1 & Level2的缩写变量（详见推导文档），这些缩写都是从原始输入直接计算的，多用于最后算出消费和闲暇后递推得到每期的资产、个人医保和财富，以及计算更高一级的缩写变量
    4. **Section 3**: 计算Level3的缩写变量，这些缩写都要用到前两级的缩写，这一层的缩写主要用于计算欧拉方程等一阶动态方程，以及生成更高级的缩写
    5. **Section 4**: 计算Level4的缩写变量，这部分其实只有一个scriptT乘子，用于从第一期消费直接得到各期消费，是欧拉方程的累乘
    6. **Section 5**: 计算Level5的缩写变量，这部分只有scriptX, scriptY，即最终Xc{1}=Y方程的左右两侧，用于求解第一期的消费（初始消费），有了初始消费，我们就可以通过scriptT乘子直接算出整条消费的路径，并且可以通过c-l Transfer Equation（消费-闲暇转换方程）得到每一期的闲暇
    7. **Section 6**: 从初始消费得到整条消费路径和闲暇路径
    8. **Section 7**: 因为我们上面得到的是没有加入边界条件（闲暇大于等于0小于等于1,而消费大于0已经在之前的合法性检查中保证了，详见数学文档）的最优路径，所以一旦存在有闲暇l超出了[0,1]的范围，就会执行此节。这一节将所有超出范围的点强行拉回0或1边界上得到修正后的闲暇路径，然后作为给定的闲暇来重新求解一条消费路径。由于欧拉方程告诉我们消费是单调递增的，所以这样处理是可以接受的。关于这节处理方法的合理性的讨论可以在数学文档中找到。
    9. **Section 8**: 利用跨期预算约束算出个人资产路径（apath）、个人医保账户路径（Phipath）以及个人总财富路径（scriptA=apath+Phipath）。我们将三个路径单独计算，然后验证是否满足相加的关系。这样做没有太大的额外计算开销，建议保留。如果满足条件：死亡时剩余财富为0（不留遗产），那么不执行Section 9而直接返回各路径结果。
    10. **Section 9**: 虽然我们理论上通过折现确实是已经满足了死亡时不留遗产的约束，但实际算起来发现还是满足不了，死亡时刻的财富（scriptA_dead）有的时候正，有的时候负。为此，我们假定前面一系列处理后得到的闲暇路径Lpath给定不再变动，从而死亡时刻财富scriptA_dead对初始消费c{1}一定单调递减（收入流给定，每期消费递增，so花的越多当然最后剩下的越少）。由于我们知道scriptA_dead(c{1}=0)一定大于等于0（绝大多数时候>0），那么一定存在唯一的一个零点real_c{1}>0。如果我们之前算得的scriptA_dead小于0，那么直接在[0,scriptA_dead]间用二分法搜索零点；如果之前算的scriptA_dead>0，那么首先一倍一倍c{1}，直到找到第一个使得scriptA_dead小于0的c{1}_2。（这个过程很快，一是因为我们之前的c{1}虽然不是零点但离零点通常不会太远，二是因为scriptA_dead(c{1})是一个凹函数，随着c{1}等距增加，scriptA_dead降到0以下的速度是逐渐加快的)。然后我们在[c{1},c{1}_2]之间用二分法搜索零点即可。搜索完成后，更新Cpath,apath,Phipath,scriptA并返回。这样一个算法虽然不一定是全局最优（因为conditional on Lpath，但Lpath was obtained not conditional on optimal Cpath），但仍然可以保证一个不错的局部最优


By Tianhao Zhao
2018-9-7 \\ 1st draft
"""
module Policy_Analytical
    using NumAlgo  # only uses Bisection()
    export PolicySolve_Analytical, AmmoReload_DATA, AmmoReload_DATA_w, AmmoReload_DATA_r, AmmoReload_PARS
# ==============================================================================
    """
        U(C::Float64,L::Float64,Q::Float64,Gamma::Float64,Alpha::Float64)

    Departable utility function, a dessert function to compute cross-sectional social wealfare;
    mathematics in latex:
        " \frac{1}{1-\gamma^{-1}} [ ((1-q)C)^{1-\gamma^{-1}} + \alpha l^{1-\gamma^{-1}} ] "

    Asking for:
        1. Gamma in (0,1) range
        2. q in (0,1) range
        3. Alpha is greater than 0
    But not validation for higher performance;
    Uses a minor amount to avoid infinity case;

    新采用的可消费、休闲可分离的CES效用函数；
    实际上是一个甜点函数，可以用于计算社会截面效用
    """
    function U(C::Float64,L::Float64,Q::Float64,Gamma::Float64,Alpha::Float64)
        U = 1/(1-1/Gamma) * ( ((1-Q)*C + 1E-6)^(1-1/Gamma) + Alpha * (L + 1E-6)^(1-1/Gamma) )
        return U
    end  # function ends


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """
        PolicySolve_Analytical(A0::Float64,Phi0::Float64,
            MAX_AGE::Int64 , RETIRE_AGE::Int64,
            DATA::Dict{String,Array{Float64,1}},
            DATA_w::Dict{String,Array{Float64,1}},
            DATA_::Dict{String,Array{Float64,1}},
            PARS::Dict{String,Float64}   )

    Analytically solves personal asset (a), individual medical account (Phi), consumption (c) and leisure before retirement (l);

    Requires:
        1. A0: initial value of personal asset account when age s=1
        2. Phi0: initial value of individual medical account when age s=1
        3. MAX_AGE: maximum age limit (to death)
        4. RETIRE_AGE: retirement age (the last year of working)
        5. DATA, DATA_w, DATA_r, PARS: datasets of parameters and economic state variables (like interest rates and wage levels)
    Returns:
        1. apath::Array{Float64,1} [len=MAX_AGE]: path of personal asset
        1. Phipath::Array{Float64,1} [len=MAX_AGE]: path of individual medical account
        1. Cpath::Array{Float64,1} [len=MAX_AGE]: path of consumption
        1. Lpath::Array{Float64,1} [len=RETIRE_AGE]: path of leisure

    In the old DP (dynamic programming, value function method) version,
    we use un-named lists (one-way array/vector) to manage data for a much greater performance;
    but now because we have got the analytical solutions,
    we use dictionaries to manage DATA, DATA_w, DATA_r and PARS;
    the structure has more flexibility in code modifying, generalization and updating when we do not need to care the code performance much;
    Specially, now do not need to care much in the order of variables;
    and, we provide several functions "AmmoReload_x" to module-lize the process (with essential validations),
    the functions can be easily used inside PolicySolve_Analytical(), also outside namespaces when calling the API;

    在之前动态规划（值函数）版本里面我们使用了大量一维数组来管理I/O以追求尽可能高的性能（因为DP太耗时了）；
    但是这意味着数据结构的脆弱性和算法的难以拓展，不利于阅读代码和后期修改模型；
    但是现在因为我们得到了解析的方法来计算路径，计算成本几乎可以忽略，
    所以可以使用拓展性更好的Dict结构来管理数据的I/O；
    虽然Dict的开销大约是数组的70倍左右，但仍然在单次求解中可以忽略（每次求解约0.0007秒）；
    由于我们要求Transition每一轮都有打印报告，所以Dict额外的开销相比于打印开销又可以完全忽略了；
    几秒钟的delay对于我们而言是完全可以接受的；
    此外，为了更方便地管理数据包，我们提供了配套的 "AmmoReload_x" 系列函数，对应DATA, DATA_w, DATA_r, PARS这四个不同种类的数据包；
    这个系列的函数用于将输入的数据转化为一个有标准键值的Dict，既用于PolicySolve_Analytical()内部，
    同时也是从外部调用这个求解API的必备套件。

    There are what should be contained in DATA, DATA_w, DATA_r & PARS [format: VarName | DictKey : Definition]
    1. DATA [all elements have length MAX_AGE]
        1. Rpath | "Rpath": interest rates in each age year
        2. Fpath | "Fpath": mortality in each age year, asking for the last element (s=S) to be 0
        3. Qpath | "Qpath": m2c ratio (medical expenses / consumption) in each year
        4. Ppath | "Ppath": outpatient / inpatient expenses
        5. cpBpath | "cpBpath": inpatient copayment rates
    2. DATA_w [all elements have length RETIRE_AGE]
        1. phiCoefpath | "phiCoefpath": personal contribution rate to medical
        2. Zetapath | "Zetapath": firm contribution rate to medical
        3. Etapath | "Etapath": firm contribution rate to pension
        4. Thetapath | "Thetapath": personal contribution rate to pension
        5. zpath |"zpath": collection rate of pension
        6. Wpath |"Wpath": wage level path
        7. doubleApath | "doubleApath": transfer (amount) from firm contribution to individual medical account
    3. DATA_r [all elements have length MAX_AGE - RETIRE_AGE]
        1. Lambdapath | "Lambdapath": pension benefits
        2. doublePpath | "doublePpath": transfer (amount) from firm contribution to those retired
    4. PARS [all elements are scalars]
        1. Mu | "Mu": consumption tax rate
        2. Sigma | "Sigma": wage tax rate
        3. Alpha | "Alpha": leisure preference than consumption
        4. Gamma | "Gamma": inter-temporal substitution elasticity
        5. Delta | "Delta": utility discount rate
        6. (NOT USED!) Varrho | "Varrho": consumption leisure substitution elasticity (DEPRECATED!!!)

    We do not do length check inside the PolicySolve_Analytical();
    the check has been performed when using "AmmoReload_x" functions to construct the inputs;
    And, most validations of Level 0 data (original inputs) should also be pre-validated in "AmmoReload_x" functions;

    Depends extra:
        1. Bisection() in module: NumAlgo
    """
    function PolicySolve_Analytical(A0::Float64,Phi0::Float64,
        MAX_AGE::Int64 , RETIRE_AGE::Int64,
        DATA::Dict{String,Array{Float64,1}},
        DATA_w::Dict{String,Array{Float64,1}},
        DATA_r::Dict{String,Array{Float64,1}},
        PARS::Dict{String,Float64}   )
        # ---------------------------------------------------------------------
        # Section 0: Unload data & parameters from data collections (Level 0)
        # ---------------------------------------------------------------------
        # 0.1 Unload those len=MAX_AGE (work in the whole lifetime)
            Rpath = DATA["Rpath"] # interest rate
            Fpath = DATA["Fpath"] # mortality
            Qpath = DATA["Qpath"] # m2c ratio
            Ppath = DATA["Ppath"] # outpatient / inpatient
            cpBpath = DATA["cpBpath"] # inpatient copayment rate
        # 0.2 Unload those len=RETIRE_AGE (only work in the working phase)
            Wpath = DATA_w["Wpath"] # wage level path
            phiCoefpath = DATA_w["phiCoefpath"] # personal contribution rate to medical
            Zetapath = DATA_w["Zetapath"] # firm contribution rate to medical
            Etapath = DATA_w["Etapath"] # firm contribution rate to pension
            Thetapath = DATA_w["Thetapath"] # persoanl contribution rate to pension
            zpath = DATA_w["zpath"] # collection rate of pension
            doubleApath = DATA_w["doubleApath"] # transfer (amount of money) from firm contribution to individual medical account
        # 0.3 Unload those len=MAX_AGE-RETIRE_AGE (only work in the retired phase)
            Lambdapath = DATA_r["Lambdapath"] # pension benefits (amount of money)
            doublePpath = DATA_r["doublePpath"] # transfer (amount) from firm contribution to those retired
        # 0.4 Unload those scalar parameters
            Mu = PARS["Mu"] # consumption tax rate
            Sigma = PARS["Sigma"] # wage tax rate
            Alpha = PARS["Alpha"] # leisure preference than consumption
            Gamma = PARS["Gamma"] # inter-temporal substitution elasticity
            Delta = PARS["Delta"] # utility discount rate
            # Varrho = PARS["Varrho"] # the elasticity has been deprecated in the new utility function
        # 0.5 Prepare result arrays
            scriptA = zeros(MAX_AGE)  # total wealth path, = apath + Phipath
            apath = zeros(MAX_AGE)  # personal asset path
            Phipath = zeros(MAX_AGE)  # individual medical account path
            Cpath = zeros(MAX_AGE)  # consumption path
            Lpath = zeros(RETIRE_AGE)  # leisure path
            # 0.5.1 initialization (essential, used in un-bounded solution)
            scriptA[1] = A0 + Phi0

        # ---------------------------------------------------------------------
        # Section 1: Special Validations of Level 0 data/short-writings
        # ---------------------------------------------------------------------
        # 1.0 MAX_AGE > 2, RETIRE_AGE < MAX_AGE (becasue the code cannot cook the case MAX_AGE=2 && RETIRE_AGE=1)
            (MAX_AGE>2 && 1<=RETIRE_AGE<MAX_AGE) || error(ErrorException("MAX_AGE > 2, RETIRE_AGE < MAX_AGE (becasue the code cannot cook the case MAX_AGE=2 && RETIRE_AGE=1)"))
        # 1.1 Interest Rate cannot be -1
            all(Rpath .!= -1.0) || error(ErrorException("interest rates cannot not be -100%, or there will be infinity"))
        # 1.2 Utility discounting rate cannot be -1
            Delta != -1.0 || error(ErrorException("utility discounting rate (Delta) cannot be -100%, or there will be infinity"))
        # 1.3 The mortility of the very last age year should be 0, if not, correct it
            Fpath[end] == 0.0 || begin Fpath[end] == 0.0; warn("the mortality in last age year should always be 0.0 exactly") end

        # ---------------------------------------------------------------------
        # Section 2: Define Level 1 & Level 2 Short-Writings
        # ---------------------------------------------------------------------
        # 2.1 Total pension contribution rates (Level 1)
            Pipath = zpath .* ( Thetapath .+ Etapath ) ./ ( 1.0 .+ zpath .* Etapath .+ Zetapath )
        # 2.2 Total medical contribution rates (Level 1)
            PiMpath = ( phiCoefpath .+ Zetapath ) ./ ( 1.0 .+ zpath .* Etapath .+ Zetapath )
        # 2.3 \scripta (Level 2)
            # scripta = 2.0 .- 1 ./ ( 1 .- Fpath )
            scripta = 1.0 .- Fpath
        # 2.4 \scriptb (Level 2)
            scriptb = 1.0 .- Sigma .- Pipath .- PiMpath
        # 2.5 \scriptd (Level 2)
            scriptd = Qpath .* ( Ppath .+ (1.0 .- cpBpath) ) ./ ( 1.0 .+ Ppath )
        # 2.6 \scriptf (Level 2)
            scriptf = ( phiCoefpath .+ doubleApath .* Zetapath ) ./ ( 1.0 .+ zpath .* Etapath .+ Zetapath )
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
            # 2.12.2 \scriptb + \scriptf should be greater than zero
                # NOTE: or all wage incomes are contributed/taxed, no left to consume
                all( (scriptb .+ scriptf) .> 0.0 ) || warn("\scriptb + \scriptf should be greater than zero, or all wage incomes are contributed/taxed, no left to consume")

        # ---------------------------------------------------------------------
        # Section 3: Define Level 3 Short-Writings (Mainly for Dynamics)
        # ---------------------------------------------------------------------
        # 3.0 Prepare empty vectors for dynamic-equation short-writings
            # NTOE: because the dynamics, all these vectors have length MAX_AGE-1 (S-1) or RETIRE_AGE-1 (Sr-1)
            # 3.0.1 Dynamic component \scriptP{s,s+1}, s=1:S-1 in Euler Equation
            scriptP = zeros(MAX_AGE-1)
            # 3.0.2 Dynamic component \scriptQ{s,s+1}, s=1:S-1 in Euler Equation
            scriptQ = zeros(MAX_AGE-1)
            # # 3.0.3 Dynamic component \scriptS{s,s+1}, s=1:Sr-1 in Leisure Dynamics
            #     # NOTE: deprecated parameter, it cannot deal with the case: RETIRE_AGE==1
            # scriptS = zeros(RETIRE_AGE-1)
        # 3.1 Con-temporal component \scriptR in consumption-leisure transfer equation (s=1:Sr-1)
            scriptR = (1.0 .- Qpath[1:RETIRE_AGE]) .* ( scriptb .+ scriptf ) .* Wpath ./ ( 1.0 .- scripth[1:RETIRE_AGE] )
            scriptR /= Alpha
        # 3.2 Filling the empty arrays
            for s in 1:MAX_AGE-1
                # NOTE: from s to s+1, locating at index = s
                # 3.2.1
                    scriptP[s] = (1+Rpath[s+1])/(1+Delta)
                    scriptP[s] *= (1-scripth[s])/(1-scripth[s+1])
                # 3.2.2
                    scriptQ[s] = (1-Qpath[s])/(1-Qpath[s+1])
            end
            # for s in 1:RETIRE_AGE-1
            #     # NOTE: from s to s+1, locating at index = s
            #     # 3.2.3
            #     scriptS[s] = (1+Rpath[s+1])/(1+Delta)
            #     scriptS[s] *= (scriptb[s]+scriptf[s])/(scriptb[s+1]+scriptf[s+1])
            #     scriptS[s] *= Wpath[s]/Wpath[s+1]
            # end
        # 3.3 Essential Valiations (IMPORTANT!) to make sure a bounded (c>0) consumption path
            all(scriptP .> 0) || error(ErrorException("invalid scriptP which leads to invalid consumption path"))
            all(scriptQ .> 0) || error(ErrorException("invalid scriptQ which leads to invalid consumption path"))
            all(scriptR .> 0) || error(ErrorException("invalid scriptR which leads to invalid consumption path"))
            # all(scriptS .> 0) || error(ErrorException("invalid scriptS which leads to invalid consumption path"))

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
            scriptX  = sum( V[1:RETIRE_AGE] .* (scriptb .+ scriptf) .* Wpath .* (1.0 .- Qpath[1:RETIRE_AGE]) .* scriptT[1:RETIRE_AGE] )
            scriptX += sum( V .* (1.0 .- scripth) .* scriptT )
        # 5.2 \scriptY{1->MAX_AGE} (the right side of final equation)
            scriptY  = sum( V[1:RETIRE_AGE] .* (scriptb .+ scriptf) .* Wpath )
            scriptY += V[1] * scriptA[1]
            scriptY += sum( V[RETIRE_AGE+1:end] .* scriptj )
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
        # 6.3 Get un-bounded leisure path through consumpton-leisure transfer equation
            for s in 1:RETIRE_AGE; Lpath[s] = Cpath[s] .* (1-Qpath[1]) / scriptR[1]^Gamma; end
        # 6.4 Check whether the leisure path meets bounds
            flag_AdjustNeeded = false
            all( 0.0 .<= Lpath .<= 1.0 ) || (flag_AdjustNeeded = true)

        # ---------------------------------------------------------------------
        # Section 7: Adjustment if the leisure path meets bounds
        # ---------------------------------------------------------------------
        # NOTE: the block runs only if flag_AdjustNeeded == true
        if flag_AdjustNeeded
            # 7.1 Find out which points are out of the bounds [0,1]
                loc_outL_0 = Lpath .<= 0.01
                loc_outL_1 = Lpath .>= 0.99
            # 7.2 Force the outside points back onto the bounds
                Lpath[loc_outL_0] = 0.01
                Lpath[loc_outL_1] = 0.99
            # 7.3 Define Level 6 short-writings
                # 7.3.1 \scriptX\tilde{1->MAX_AGE} (the new left side of final equation)
                    scriptX_tilde = sum( V .* (1.0 .- scripth) .* scriptT )
                # 7.3.2 \scriptY\tilde{1->MAX_AGE} (the new right side of final equation)
                    scriptY_tilde  = sum( V[1:RETIRE_AGE] .* (scriptb .+ scriptf) .* Wpath .* ( 1.0 .- Lpath ) )
                    scriptY_tilde += V[1] * scriptA[1]
                    scriptY_tilde += sum( V[RETIRE_AGE+1:end] .* scriptj )
            # 7.4 Get new consumption path still meets Euler Equation
                # 7.4.1 Get first-year consumption
                    Cpath[1] = scriptY_tilde ./ scriptX_tilde
                    Cpath[1] *= (1+Rpath[1])  # discount to s=1
                # 7.4.2 Generalize the path through Euler Equation
                    for s in 2:MAX_AGE; Cpath[s] = Cpath[1] .* scriptT[s]; end
            # 7.5 Validation
                all( 0.0 .<= Cpath .< Inf ) || error(ErrorException("negative or Inf consumption found in adjusted solution"))
        end

        # ---------------------------------------------------------------------
        # Section 8: Get apath & Phipath
        # ---------------------------------------------------------------------
        # 8.1 Get personal asset path (apath)
            apath[1:end],apath_dead = GetFullPath_a( A0,
                    Cpath, Lpath, Rpath, scriptd, scripta, Wpath, scriptb, Lambdapath   )

        # 8.2 Get individual medical account path (Phipath)
            Phipath[1:end],Phipath_dead,Gaps = GetFullPath_Phi( Phi0,
                    Cpath, Lpath, Rpath, scriptg, scripta, Wpath, scriptf, doublePpath   )
            # gap covered by personal asset (apath)
            apath = apath .+ Gaps

        # 8.3 Get total wealth path (scriptA)
            # NOTE: scriptA_dead indicates the remaining wealth/capital when agent dies (at the end of the last age year)  # scriptA[1] = A0 + Phi0
            scriptA[1:end],scriptA_dead = GetFullPath_scriptA( A0 + Phi0, Cpath, Lpath,
                        Rpath, scripth, scripta, Wpath, scriptb, scriptf, scriptj  )

        # 8.4 Check whether the relationship "a+Phi=scriptA" met
            findmax( abs.( scriptA .- apath .- Phipath ) )[1] < (1E-06) || error(ErrorException("the relationship scriptA = a + Phi not met"))

        # ---------------------------------------------------------------------
        # Section 9: Search a c_1 which meets scriptA_dead == 0 (if scriptA_dead != 0)
        # ---------------------------------------------------------------------
        # NOTE: keep Euler equation determined; keep leisure path determined (cauz we've adjusted it)
        # NOTE: use bisection method to search the zero point (scriptA_dead = 0)
        # 9.1 Tolerance to believe scriptA_dead = 0
            TOL = 1E-08
        # 9.2 Use bi-sectional method to search the zero point of scriptA_dead
            if abs( scriptA_dead - 0.0 ) > TOL
                # 9.2.0 Define a temp function to use c_1 to compute scriptA_dead (opjevt function)
                    objfunc(C0::Float64) = begin
                        # a. construct a Cpath by C0 input
                            tmp_Cpath = C0 * scriptT
                        # b. get the wealth when dead, then return it
                            tmp_scriptA,tmp_scriptA_dead = GetFullPath_scriptA( A0 + Phi0, tmp_Cpath, Lpath, Rpath, scripth, scripta, Wpath, scriptb, scriptf, scriptj  )
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
                            Cpath[1] = Bisection(objfunc,C1,C2,MAXLOOP=10000,TOL=1E-12)
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
                            Cpath[1] = Bisection(objfunc,C1,C2,MAXLOOP=10000,TOL=1E-12)
                    end  # branch ends (CASES)

                # 9.2.4 If sovled successfully
                    Cpath[1:end] = Cpath[1] .* scriptT
                    # get scriptA_dead and the whole path
                    scriptA[1:end],scriptA_dead = GetFullPath_scriptA( A0 + Phi0, Cpath, Lpath, Rpath, scripth, scripta, Wpath, scriptb, scriptf, scriptj  )
                    # get Phipath
                    Phipath[1:end],Phipath_dead,Gaps = GetFullPath_Phi( Phi0, Cpath, Lpath, Rpath, scriptg, scripta, Wpath, scriptf, doublePpath   )
                    # get apath
                    apath[1:end] = scriptA .- Phipath
                # 9.2.5 Check scriptA >=0
                    # all(scriptA .>= 0.0) || error(ErrorException("Fail to solve a non-negative Wealth path, please try less generations!"))
                    abs(scriptA_dead) < 1E-3 || error(ErrorException(string("Fail to solve a non-negative Wealth path, please try less generations! ",scriptA_dead)  ))
                # 9.2.6 Finally, return
                    return apath,Phipath,Cpath,Lpath
            end  # branch/search ends

        # Check scriptA >=0
            # all(scriptA .>= 0.0) || error(ErrorException("Fail to solve a non-negative Wealth path, please try less generations!"))
            abs(scriptA_dead) < 1E-3 || error(ErrorException(string("Fail to solve a non-negative Wealth path, please try less generations! ",scriptA_dead)  ))
        # Finally, return the results
            return apath,Phipath,Cpath,Lpath
    end  # function ends


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """
        AmmoReload_DATA( MAX_AGE::Int64,
            Rpath::Array{Float64,1}, Fpath::Array{Float64,1},
            Qpath::Array{Float64,1}, Ppath::Array{Float64,1}, cpBpath::Array{Float64,1} )

    Returns a dictionary DATA; also do validations before packaging
    """
    function AmmoReload_DATA( MAX_AGE::Int64,
        Rpath::Array{Float64,1}, Fpath::Array{Float64,1},
        Qpath::Array{Float64,1}, Ppath::Array{Float64,1}, cpBpath::Array{Float64,1} )
        # --------------------------
        # 1. Validations
            # 1.1 Length Check
                @assert(MAX_AGE==length(Rpath)==length(Fpath)==length(Qpath)==length(Ppath)==length(cpBpath),"Uncompatible MAX_AGE-length inputs")
            # 1.2 Range Check
                all( Rpath .!= -1.0 ) || error(ErrorException("interest rates cannot be equal to -100%, or there will be all zeros assets"))
                all( 0.0 .<= Fpath .< 1.0 ) || error(ErrorException("mortalities should be in the half open range [0,1)"))
                all( 0.0 .< Qpath .< 1.0 ) || error(ErrorException("m2c ratios should be in the open range (0,1)"))
                all( 0.0 .< Ppath ) || error(ErrorException("outpatient/inpatient expenses ratio should be in greater than 0"))
                all( 0.0 .<= cpBpath .<= 1.0 ) || error(ErrorException("inpatient copayment rates should be in the closed range [0,1]"))
            # 1.3 Special Check
                Fpath[end] == 0.0 || begin Fpath[end] = 0.0; warn("the last element of Fpath should be exact 0.0") end
        # 2. Packaging
            DATA = Dict(
                "Rpath" => Rpath,
                "Fpath" => Fpath,
                "Qpath" => Qpath,
                "Ppath" => Ppath,
                "cpBpath" => cpBpath
            )
        # 3. Return
        return DATA::Dict{String,Array{Float64,1}}
    end  # function ends


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """
        AmmoReload_DATA_w( RETIRE_AGE::Int64,
            Wpath::Array(Float64,1),
            phiCoefpath::Array(Float64,1), Zetapath::Array(Float64,1),
            Etapath::Array(Float64,1), Thetapath::Array(Float64,1),
            zpath::Array(Float64,1), doubleApath::Array(Float64,1)   )

    Returns a dictionary DATA_w; also do validations before packaging
    """
    function AmmoReload_DATA_w( RETIRE_AGE::Int64,
        Wpath::Array{Float64,1},
        phiCoefpath::Array{Float64,1}, Zetapath::Array{Float64,1},
        Etapath::Array{Float64,1}, Thetapath::Array{Float64,1},
        zpath::Array{Float64,1}, doubleApath::Array{Float64,1}   )
        # -----------------------
        # 1. Validations
            # 1.1 Length Check
            @assert(RETIRE_AGE==length(Wpath)==length(phiCoefpath)==length(Zetapath)==length(Etapath)==length(Thetapath)==length(doubleApath)==length(Etapath),"Uncompatible RETIRE_AGE-length inputs")
            # 1.2 Range Check
                all( Wpath .>= 0.0 ) || error(ErrorException("wage levels should be greater than 0"))
                all( 0.0 .< phiCoefpath .< 1.0 ) || error(ErrorException("phiCoefpath should be in the open range (0,1)"))
                all( 0.0 .< Zetapath .< 1.0 ) || error(ErrorException("Zetapath should be in the open range (0,1)"))
                all( 0.0 .< Etapath .< 1.0 ) || error(ErrorException("Etapath should be in the open range (0,1)"))
                all( 0.0 .< Thetapath .< 1.0 ) || error(ErrorException("Thetapath should be in the open range (0,1)"))
                all( 0.0 .< zpath .< 1.0 ) || error(ErrorException("zpath should be in the open range (0,1)"))
                all( 0.0 .< doubleApath .< 1.0 ) || error(ErrorException("doubleApath should be in the open range (0,1)"))
        # 2. Packaging
            DATA_w = Dict(
                "phiCoefpath" => phiCoefpath,
                "Zetapath" => Zetapath,
                "Etapath" => Etapath,
                "Thetapath" => Thetapath,
                "zpath" => zpath,
                "Wpath" => Wpath,
                "doubleApath" => doubleApath
            )
        # 3. Return
        return DATA_w::Dict{String,Array{Float64,1}}
    end  # function ends


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """
        AmmoReload_DATA_r( MAX_AGE::Int64, RETIRE_AGE::Int64,
            Lambdapath::Array{Float64,1}, doublePpath::Array{Float64,1}   )

    Returns a dictionary DATA_r; also do validations before packaging
    """
    function AmmoReload_DATA_r( MAX_AGE::Int64, RETIRE_AGE::Int64,
        Lambdapath::Array{Float64,1}, doublePpath::Array{Float64,1}   )
        # ------------------------
        # 1. Validations
            # 1.1 Length Check
            @assert((MAX_AGE-RETIRE_AGE)==length(Lambdapath)==length(doublePpath),"Uncompatible MAX_AGE - RETIRE_AGE length inputs")
            # 1.2 Range Check
                all( Lambdapath .>= 0.0 ) || error(ErrorException("pension benefits should be greater than or equal to 0"))
                all( doublePpath .>= 0.0 ) || error(ErrorException("medical transfer to the retired should be greater than or equal to 0"))
        # 2. Packaging
            DATA_r = Dict(
                "Lambdapath" => Lambdapath,
                "doublePpath" => doublePpath
            )
        # 3. Return
        return DATA_r::Dict{String,Array{Float64,1}}
    end  # function ends


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """
        AmmoReload_PARS( Mu::Float64, Sigma::Float64,
            Alpha::Float64, Gamma::Float64, Delta::Float64  )

    Returns a dictionary PARS; also do validations before packaging
    """
    function AmmoReload_PARS( Mu::Float64, Sigma::Float64,
        Alpha::Float64, Gamma::Float64, Delta::Float64  )
        # ----------------------
        # 1. Range Check
            Alpha > 0.0 || error(ErrorException("Alpha should be greater than 0"))
            Delta != -1.0 || error(ErrorException("Delta cannot be -1, or there will be infinite utility"))
            (0.0 < Gamma < 1.0) || error(ErrorException("Gamma should be in the open range (0,1)"))
            (0.0 < Sigma < 1.0) || error(ErrorException("Sigma should be in the open range (0,1)"))
        # 2. Packaging
            PARS = Dict(
                "Alpha" => Alpha,
                "Delta" => Delta,
                "Gamma" => Gamma,
                "Mu" => Mu,
                "Sigma" => Sigma
            )
        # 3. Return
        return PARS::Dict{String,Float64}
    end  # function ends


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """
        GetFullPath_scriptA( scriptA0::Float64,
            Cpath::Array{Float64,1}, Lpath::Array{Float64,1},
            Rpath::Array{Float64,1}, scripth::Array{Float64,1}, scripta::Array{Float64,1},
            Wpath::Array{Float64,1}, scriptb::Array{Float64,1}, scriptf::Array{Float64,1},
            scriptj::Array{Float64,1}   )

    a subfunction, computes whole WEALTH path from s=1 to s=S;
    returns a vector tmp_scriptA (copy, not reference), also a number scriptA_dead,
    which indicates the wealth when dead (used in bisection search);
    and, no validation cauz the function is only used inside the PolicySolve_Analytical();
    """
    function GetFullPath_scriptA( scriptA0::Float64,
        Cpath::Array{Float64,1}, Lpath::Array{Float64,1},
        Rpath::Array{Float64,1}, scripth::Array{Float64,1}, scripta::Array{Float64,1},
        Wpath::Array{Float64,1}, scriptb::Array{Float64,1}, scriptf::Array{Float64,1},
        scriptj::Array{Float64,1}   )
        # -----------------
        # 0. backup (to better understand the logic in PolicySolve_Analytical)
            # tmp_scriptA = zeros(MAX_AGE)
            # tmp_Cpath = C0 * scriptT
            # tmp_scriptA[1] = A0 + Phi0

        # 0. prepare empty vector for wealth (scriptA)
            MAX_AGE = length(Cpath); RETIRE_AGE = length(Lpath)
            tmp_scriptA = zeros(MAX_AGE)
            tmp_scriptA[1] = scriptA0
        # 1. anti-retire
            for s in 1:RETIRE_AGE
                tmp_scriptA[s+1] = (1+Rpath[s])*tmp_scriptA[s] + (scriptb[s]+scriptf[s])*Wpath[s]*(1-Lpath[s]) - (1-scripth[s])*Cpath[s]
                tmp_scriptA[s+1] /= scripta[s]
            end
        # 2. post-retire
            # NOTE: the loop is disigned for case: MAX_AGE-RETIRE_AGE>1, if MAX_AGE-RETIRE_AGE=1, only need to run section 3
            if MAX_AGE-RETIRE_AGE>1
                for s in RETIRE_AGE:MAX_AGE-1
                    tmp_scriptA[s+1] = (1+Rpath[s])*tmp_scriptA[s] + scriptj[s-RETIRE_AGE+1] - (1-scripth[s])*Cpath[s]
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
        GetFullPath_Phi( Phi0::Float64,
            Cpath::Array{Float64,1}, Lpath::Array{Float64,1},
            Rpath::Array{Float64,1}, scriptg::Array{Float64,1}, scripta::Array{Float64,1},
            Wpath::Array{Float64,1}, scriptf::Array{Float64,1},
            doublePpath::Array{Float64,1}   )

    a subfunction, computes whole INDIVIDUAL MEDICAL ACCOUNT path from s=1 to s=S;
    returns a vector tmp_Phipath (copy, not reference), also a number tmp_Phipath_dead, and a path of g ap covered by personal asset (apath) where negative for gaps;
    which indicates the individual medical account when dead (used in bisection search);
    and, no validation cauz the function is only used inside the PolicySolve_Analytical();
    """
    function GetFullPath_Phi( Phi0::Float64,
        Cpath::Array{Float64,1}, Lpath::Array{Float64,1},
        Rpath::Array{Float64,1}, scriptg::Array{Float64,1}, scripta::Array{Float64,1},
        Wpath::Array{Float64,1}, scriptf::Array{Float64,1},
        doublePpath::Array{Float64,1}   )
        # -----------------
        # 0. backup (to better understand the logic in PolicySolve_Analytical)
            # tmp_scriptA = zeros(MAX_AGE)
            # tmp_Cpath = C0 * scriptT
            # tmp_scriptA[1] = A0 + Phi0

        # 0. prepare empty vector for wealth (scriptA)
            MAX_AGE = length(Cpath); RETIRE_AGE = length(Lpath)
            tmp_Phipath = zeros(MAX_AGE)
            tmp_Phipath[1] = Phi0
        # 1. anti-retire
            for s in 1:RETIRE_AGE
                tmp_Phipath[s+1] = (1+Rpath[s])*tmp_Phipath[s] + scriptf[s]*Wpath[s]*(1-Lpath[s]) + scriptg[s]*Cpath[s]
                tmp_Phipath[s+1] /= scripta[s]
            end
        # 2. post-retire
            # NOTE: the loop is disigned for case: MAX_AGE-RETIRE_AGE>1, if MAX_AGE-RETIRE_AGE=1, only need to run section 3
            if MAX_AGE-RETIRE_AGE>1
                for s in RETIRE_AGE:MAX_AGE-1
                    tmp_Phipath[s+1] = (1+Rpath[s])*tmp_Phipath[s] + doublePpath[s-RETIRE_AGE+1] + scriptg[s]*Cpath[s]
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
    function GetFullPath_a( A0::Float64,
        Cpath::Array{Float64,1}, Lpath::Array{Float64,1},
        Rpath::Array{Float64,1}, scriptd::Array{Float64,1}, scripta::Array{Float64,1},
        Wpath::Array{Float64,1}, scriptb::Array{Float64,1},
        Lambdapath::Array{Float64,1}   )
        # -----------------
        # 0. backup (to better understand the logic in PolicySolve_Analytical)
            # tmp_scriptA = zeros(MAX_AGE)
            # tmp_Cpath = C0 * scriptT
            # tmp_scriptA[1] = A0 + Phi0

        # 0. prepare empty vector for wealth (scriptA)
            MAX_AGE = length(Cpath); RETIRE_AGE = length(Lpath)
            tmp_apath = zeros(MAX_AGE)
            tmp_apath[1] = A0
        # 1. anti-retire
            for s in 1:RETIRE_AGE
                tmp_apath[s+1] = (1+Rpath[s])*tmp_apath[s] + scriptb[s]*Wpath[s]*(1-Lpath[s]) - (1 - scriptd[s])*Cpath[s]
                tmp_apath[s+1] /= scripta[s]
            end
        # 2. post-retire
            # NOTE: the loop is disigned for case: MAX_AGE-RETIRE_AGE>1, if MAX_AGE-RETIRE_AGE=1, only need to run section 3
            if MAX_AGE-RETIRE_AGE>1
                for s in RETIRE_AGE:MAX_AGE-1
                    tmp_apath[s+1] = (1+Rpath[s])*tmp_apath[s] + Lambdapath[s-RETIRE_AGE+1] - (1 - scriptd[s])*Cpath[s]
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
