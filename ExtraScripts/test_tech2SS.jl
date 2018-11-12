# Different tech level, different steady states
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
# 这个测试代码的作用是给定一组基本参数集，然后给定不同水平的技术系数（TFP），观察技术系数带来的稳态下（财富、消费、休闲、利率等）
# 的变化。用以挑选最合适的TFP水平
# -------------
# 首先决定测试的griding density
TEST_GRID = 10
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

# 接下来确定测试的利率水平序列
input_tech = Array(linspace( 1.0, 3.0, TEST_GRID ))



# 然后开始测试初始稳态，并存储每一次的结果。由于Arange和GRID_LEN的选择会很大程度上影响最终能不能收敛（但一旦收敛，收敛到同一个结果），
# 所以比较推荐手动修改idx然后调试计算每一个tech水平对应的内容
idx = 10

    YEAR = MAX_YEAR
    # 替换新的技术系数水平
    vec_tech[YEAR] = input_tech[idx]
    # 使用最终稳态人口（可选）
    # mat_N[YEAR,:] = mat_N[MAX_YEAR,:]
    # Guesses
    guess_r = 0.10  # net interest rate
    guess_L = 0.6  # labour
    guess_q = 0.15  # m2c ratio
    # Folding data packages
    DATA_vec = [ mat_N[YEAR,:], mat_Survival[YEAR,:], mat_MA2MB[YEAR,:], vec_epsilon ]
    DATA_num = [ num_kappa, num_gamma, num_varrho, num_delta, num_alpha, num_mu, num_sigma,
        vec_beta[YEAR], vec_tech[YEAR], vec_doubleK[YEAR], vec_eta[YEAR], vec_theta[YEAR], vec_z[YEAR], vec_zeta[YEAR],
        vec_phiCoef[YEAR], vec_doubleA[YEAR], vec_doubleB[YEAR], vec_cpB[YEAR]  ]
    # CALLING: initial steady state searching
    RESULT_vec, RESULT_num = SStateSearch(guess_r,guess_L,guess_q, MAX_AGE, RETIRE_AGE, DATA_vec, DATA_num,
        flag_diagnose = true, MagicNum = 2.0, Arange = [0.0, 450.0],
        TOL_SS = 1E-04, TOL_POLICY = 1E-08, GRID_LEN = 277, MAX_ITER = 50 )
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



# =====================================================
# using Plots # Plotting, using gr() or plotly() backends
using PyPlot
# 计算完成后进行展示，使用GR() backend，Plots style
# 1. 首先展示wealth
plt_Wealth = plot( test_Wealth')
xlabel("Age"); title("Wealth by TFP"); legend([string("TFP=",round(x,2)) for x in input_tech]   )
grid("on"); xlim([1,MAX_AGE])
# 2. 然后展示Consumption
plt_Consumption = plot(test_Consumption')
xlabel("Age"); title("Consumption by TFP"); legend([string("TFP=",round(x,2)) for x in input_tech]   )
grid("on"); xlim([1,MAX_AGE])
# 3. 接下来是leisure
plt_Lesiure = plot(test_Leisure[:,1:RETIRE_AGE]')
xlabel("Age"); title("Leisure (before retire) by TFP"); legend([string("TFP=",round(x,2)) for x in input_tech]   )
grid("on"); xlim([1,RETIRE_AGE])
# 4. 还有利率
plt_r = plot( input_tech, test_r )
xlabel("TFP"); ylabel("Net interest rate");  title("TFP -> Interest rate")
grid("on")
# 5. 以及替代率
plt_rep = plot( input_tech, test_rep )
xlabel("TFP"); ylabel("Replacement rate");  title("TFP -> Replacement rate")
grid("on")





# 当然，别忘了先存出来，不然够蛋疼的
writecsv("output\\final_demog\\test_Wealth.csv",test_Wealth)
writecsv("output\\final_demog\\test_Consumption.csv",test_Consumption)
writecsv("output\\final_demog\\test_Leisure.csv",test_Leisure)
writecsv("output\\final_demog\\test_r.csv",test_r)
writecsv("output\\final_demog\\test_rep.csv",test_rep)

# # 后边用的时候就直接读进来，不重新算了
# test_Wealth = readcsv("output\\final_demog\\test_Wealth.csv")
# test_Consumption = readcsv("output\\final_demog\\test_Consumption.csv")
# test_Leisure = readcsv("output\\final_demog\\test_Leisure.csv")
# test_r = readcsv("output\\final_demog\\test_r.csv")
# test_rep = readcsv("output\\final_demog\\test_rep.csv")



#
