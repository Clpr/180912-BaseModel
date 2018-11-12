# 在确定了人口结构是影响稳态下居民生命期财富曲线和消费决策曲线形状的最主要因素后，
# 我们生成一系列不同形状的 非增 人口分布曲线（长度为MAX_AGE）和对应的死亡概率，
# 来讨论人口曲线形状&特殊特征（如断点、无死亡风险等）带来的对生命期财富曲线形状的影响
# ==============================================================
push!(LOAD_PATH,"./modules")
using PolicyFunctions # Functions about Policy Function Solving
using GeneProc # General Processes to deal with data, ploting, reporting and others
using SteadyState # Steady State Searching
using Transition # Transition Path Searching
using DataFrames # For data packing
# ================== BASIC PARMS ===============================================
const MAX_YEAR = 250
const MAX_AGE = 80
const RETIRE_AGE = 40
# const MAX_ITER = 200
const REL_TOL = 1E-05
# ================== DATA READING, ROOMING =====================================
# Declare empty data structures
mat_a,mat_Phi,mat_c,mat_l,vec_L,vec_K,vec_Y,vec_r,vec_wmean,mat_wprofiled,
vec_D,vec_G,vec_I,vec_D2Y,mat_Lambda,mat_M,mat_MA,mat_MB,
vec_oCoef,vec_gapUMP,vec_TRw,vec_TRc,vec_U,vec_doubleP = func_VarsDeclare(MAX_YEAR,MAX_AGE,RETIRE_AGE)
# Initilize parameters
num_kappa,num_gamma,num_varrho,num_delta,num_alpha,num_mu,num_sigma,
vec_beta,vec_tech,vec_x,vec_doubleK,vec_eta,vec_theta,vec_zeta,vec_phiCoef,vec_doubleA,vec_doubleB,vec_cpB,vec_z = func_ParsInit(MAX_YEAR,MAX_AGE,RETIRE_AGE)
# Demography Data Reading
mat_N,mat_Survival = func_DemoRead(MAX_YEAR::Int64,MAX_AGE::Int64,RETIRE_AGE::Int64,"./data/TotalPopu.csv")
# Age-dependent Data Reading
vec_epsilon,mat_MA2MB = func_AgeDatRead(MAX_YEAR::Int64,MAX_AGE::Int64,RETIRE_AGE::Int64, "./data/DataByAge.csv"::String,colidx_epsilon=1,colidx_MA=2,colidx_MB=3)

# ================== GRIDING SEARCH ON DIFFERENT TECHNOLOGY LEVELS =====================
# 接下来换用中文注释。
# 这个测试代码的作用是给定一组基本参数集（初始稳态参数集），然后给定不同水平的人口分布，观察带来的稳态下（财富、消费、休闲、利率等）的变化
# 在这种测试环境下，由于人口总量都被sclared到1，所以只考虑形状带来的问题
# -------------
# 最首先准备一个空矩阵
TEST_GRID = 7
DemogList = zeros(TEST_GRID,MAX_AGE)
SurvivalList = ones(TEST_GRID,MAX_AGE)
# 然后确定测试的曲线种类：
# 1. 完全平坦（even）=> 死亡率为0，无死亡风险，无因死亡带来的意外财产/遗产收入
DemogList[1,:] = 1/MAX_AGE
# 2. 直线下降（linear） => 凸死亡率/凹生存率，人口以线性下降至0
DemogList[2,:] = Array(linspace(1.0,0.01,MAX_AGE))
DemogList[2,:] /= sum(DemogList[2,:])
# 3. 凸函数下降（convex） => 人口以凸函数形式下降（先快速下降然后平坦），死亡率逐渐降低
DemogList[3,:] = 1./(1:MAX_AGE)
DemogList[3,:] /= sum(DemogList[3,:])
# 4. 凹函数下降（concave） => 人口以凹函数形式下降（先平坦然后快速下降），死亡率逐渐升高
DemogList[4,:] = log.(MAX_AGE:-1:1)
DemogList[4,:] /= sum(DemogList[4,:])
# 5. 断点（breaking point） => 两段平坦的水平直线间在某一年存在断崖式下跌，表明该岁存在突然的高死亡风险
DemogList[5,:] = append!( ones(30), 0.5*ones(MAX_AGE-30) )
DemogList[5,:] /= sum(DemogList[5,:])
# 6. 原本初始稳态的人口数据 => 该数据基本是直线下降，但在退休前存在两段平坦（even）和一段几乎是断点的形状
DemogList[6,:] = copy(mat_N[1,:])
# 7. 原本最终稳态的人口数据 => 该数据是平稳数据
DemogList[7,:] = copy(mat_N[end,:]./sum(mat_N[end,:]))
# ------------
# 当然，千万别忘了对应的生存概率
SurvivalList[:,1:end-1] = DemogList[:,2:end]./DemogList[:,1:end-1]
# ------------
# 然后决定观测的变量：
# 1. 个人财富
test_Wealth = zeros(TEST_GRID,MAX_AGE)
# 2. 消费
test_Consumption = zeros(TEST_GRID,MAX_AGE)
# 3. 休闲
test_Leisure = zeros(TEST_GRID,MAX_AGE)
# 4. 利率
test_r = zeros(TEST_GRID)
# 5. 替代率
test_rep = zeros(TEST_GRID)

# =======================================
# NOTE: 我们注意到不同的曲线形状通常要求不同的（Arange上限，GRID_LEN）组合才能保证收敛。为此我们给出推荐组合（精度均至少为1E-04，且保证100轮内收敛）
# 1. [30.0, 83]
# 2. [30.0, 83]
# 3. [40.0, 53]
# 4. [25.0, 41]
# 5. [30.0, 41]
# 6. [40.0, 59]
# 7. [40.0, 59]
# 读者可以轻易使用一个zip或list来自定义算法参数集以提高（复现）测试效率
# 下面的集合用于使用初始稳态技术系数时
ulimA = [60.0, 200.0, 80.0, 80.0, 30.0, 40.0, 40.0]   # Case 2, 3 无法收敛
gridA = [53,   57,    53,   41,   41,   59,   59  ]
# 下面的集合用于当使用最终稳态的技术系数时
# vec_tech[MAX_YEAR] = copy(vec_tech[end])  # 可选是否使用最后一年的技术系数
# ulimA = [400.0, 400.0, 2000.0, 350.0, 320.0, 900.0, 500.0] # 实际上此时CASE 3 并不收敛
# gridA = [157,   83,    117,    61,    111,   141,   71  ]

# 然后进行测试
for idx in 1:TEST_GRID
    idx = 4  # 用于手动调参，调整上面算法参数集合用

    YEAR = 1  # 统统以第1年为环境进行测试
    # Guesses
    guess_r = 0.10; guess_L = 0.6; guess_q = 0.15
    # Folding data packages
    DATA_vec = [ DemogList[idx,:], SurvivalList[idx,:], # 人口数据在这一行修改
        mat_MA2MB[YEAR,:], vec_epsilon ]

    DATA_num = [ num_kappa, num_gamma, num_varrho, num_delta, num_alpha, num_mu, num_sigma,
        vec_beta[YEAR], vec_tech[YEAR], vec_doubleK[YEAR], vec_eta[YEAR], vec_theta[YEAR], vec_z[YEAR], vec_zeta[YEAR],
        vec_phiCoef[YEAR], vec_doubleA[YEAR], vec_doubleB[YEAR], vec_cpB[YEAR]  ]
    # CALLING: initial steady state searching
    println("*"^40," ","Case: ",idx," ","*"^40)
    RESULT_vec, RESULT_num = SStateSearch(guess_r,guess_L,guess_q, MAX_AGE, RETIRE_AGE, DATA_vec, DATA_num,
        flag_diagnose = true, MagicNum = 2.0, Arange = [0.0, ulimA[idx]  ],   # 使用了特定的算法参数（Arange和GRID_LEN）来加快搜索
        TOL_SS = 1E-04, TOL_POLICY = 1E-08, GRID_LEN = gridA[idx] , MAX_ITER = 50 )
    # 拆包结果
    mat_wprofiled[YEAR,:], mat_Lambda[YEAR,:], mat_a[YEAR,:], mat_Phi[YEAR,:],
        mat_c[YEAR,:], mat_l[YEAR,1:RETIRE_AGE], mat_M[YEAR,:], mat_MA[YEAR,:], mat_MB[YEAR,:] = RESULT_vec
    vec_K[YEAR], vec_L[YEAR], vec_Y[YEAR], vec_r[YEAR], vec_wmean[YEAR], vec_oCoef[YEAR], vec_TRc[YEAR], vec_TRw[YEAR],
        vec_gapUMP[YEAR], vec_D[YEAR], vec_D2Y[YEAR], vec_G[YEAR], vec_I[YEAR], vec_doubleP[YEAR] = RESULT_num

    # 记录结果
    test_Wealth[idx,:] = mat_a[YEAR,:] .+ mat_Phi[YEAR,:]
    test_Consumption[idx,:] = mat_c[YEAR,:]
    test_Leisure[idx,:] = mat_l[YEAR,:]
    test_r[idx] = vec_r[YEAR]
    test_rep[idx] = mean(mat_Lambda[YEAR,:]) / sum( mat_wprofiled[YEAR,:].*( 1-mat_l[YEAR,1:RETIRE_AGE] ) ) * RETIRE_AGE

end
# 当然，别忘了先存出来，不然够蛋疼
writecsv("output\\测试_不同人口曲线形状\\DemographyUsed.csv",DemogList)
writecsv("output\\测试_不同人口曲线形状\\SurvivalProbUsed.csv",SurvivalList)
writecsv("output\\测试_不同人口曲线形状\\test_Wealth.csv",test_Wealth)
writecsv("output\\测试_不同人口曲线形状\\test_Consumption.csv",test_Consumption)
writecsv("output\\测试_不同人口曲线形状\\test_Leisure.csv",test_Leisure)
writecsv("output\\测试_不同人口曲线形状\\test_r.csv",test_r)
writecsv("output\\测试_不同人口曲线形状\\test_rep.csv",test_rep)


# ===================== 最后是绘图 ===========================
using PyPlot
# Case 1 ---------------
subplot(241)
plot(DemogList[1,:],"-")
legend(["Demography","Wealth"]); xlabel("Age")












#
