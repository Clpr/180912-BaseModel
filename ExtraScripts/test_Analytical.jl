# 测试解析解的新PolicySolve_Analytical()函数。引用Policy_Analytical.jl模块
# ========================================================
push!(LOAD_PATH,"./")
push!(LOAD_PATH,"./modules/")
using PyPlot  # plotting
using Policy_Analytical  # policy analytically
using Policy_Analytical_Retired  # policy analytically specially for the retired phase
# ========================================================
MAX_AGE = 80
RETIRE_AGE = 35
# ========================================================
    Rpath = 0.06 * ones(MAX_AGE)
    Fpath = 0.01 * ones(MAX_AGE); Fpath[end] = 0.0
    Qpath = 0.15 * ones(MAX_AGE)
    Ppath = 1.20 * ones(MAX_AGE)
    cpBpath = 0.30 * ones(MAX_AGE)
    # ----------------------------------------------------
    phiCoefpath = 0.02*ones(RETIRE_AGE)
    Zetapath = 0.06*ones(RETIRE_AGE)
    Etapath = 0.21*ones(RETIRE_AGE)
    Thetapath = 0.05*ones(RETIRE_AGE)
    Wpath = 1.50*ones(RETIRE_AGE)
    doubleApath = 0.30*ones(RETIRE_AGE)
    zpath = 0.90*ones(RETIRE_AGE)
    # ------------------------------
    Lambdapath = 0.02*ones(MAX_AGE-RETIRE_AGE)
    doublePpath = 0.01*ones(MAX_AGE-RETIRE_AGE)
    # ------------------------------
    Mu = 0.1; Sigma = 0.24; Alpha = 1.5; Gamma = 0.5; Varrho = 0.8; Delta = 0.015
    # -------------------------------
# ========================================================
DATA = AmmoReload_DATA( MAX_AGE, Rpath, Fpath, Qpath, Ppath, cpBpath )
DATA_w = AmmoReload_DATA_w( RETIRE_AGE, Wpath, phiCoefpath, Zetapath,
        Etapath, Thetapath, zpath, doubleApath   )
DATA_r = AmmoReload_DATA_r( MAX_AGE, RETIRE_AGE, Lambdapath, doublePpath   )
PARS = AmmoReload_PARS( Mu, Sigma, Alpha, Gamma, Delta  )
# ========================================================
# 1. 测试完整路径求解
A0 = 0.0; Phi0 = 0.0
@time apath, Phipath, Cpath, Lpath = PolicySolve_Analytical(A0,Phi0, MAX_AGE ,
    RETIRE_AGE, DATA, DATA_w, DATA_r, PARS   )

scriptA = apath .+ Phipath
plot(1:MAX_AGE,scriptA);title("Wealth")
xlabel("Age");grid(true)
# figure();plot(Cpath);title("C")
# figure();plot(Lpath);title("leisure")
# 2. 测试仅退休期路径求解



# 一口气400轮，80期测试，结果平均2.08~2.09秒，意味着2.5秒内完成1轮transition搜索
# @time for x in 1:400 PolicySolve_Analytical(A0,Phi0, MAX_AGE , RETIRE_AGE, DATA, DATA_w, DATA_r, PARS   ); end

# # ==========================================================
# # 接下来测试Analytical_Retired
# # 从上面直接截取数据而不是重新生成
# # 然后比较得到的路径和上边整条路径是否一致
# # ==========================================================
# START_AGE = RETIRE_AGE + 1   # START_AGE > RETIRE_AGE
#     # START_AGE = MAX_AGE - 1  # special case 1
#     # START_AGE = MAX_AGE  # special case 2
# # 1. 截取数据
# r_A0 = apath[START_AGE]
# r_Phi0 = Phipath[START_AGE]
#     r_Rpath = Rpath[START_AGE:end]
#     r_Fpath = Fpath[START_AGE:end]
#     r_Qpath = Qpath[START_AGE:end]
#     r_Ppath = Ppath[START_AGE:end]
#     r_cpBpath = cpBpath[START_AGE:end]
#     r_Lambdapath = Lambdapath[START_AGE-RETIRE_AGE:end]
#     r_doublePpath = doublePpath[START_AGE-RETIRE_AGE:end]
# # 2. 装载数据
# r_DATA = AmmoReload_Retired_DATA( MAX_AGE, START_AGE,
#             r_Rpath, r_Fpath,
#             r_Qpath, r_Ppath, r_cpBpath,
#             r_Lambdapath, r_doublePpath  )
# r_PARS = AmmoReload_Retired_PARS( Mu, Gamma, Delta  )
# # 3. 求解
# @time r_apath, r_Phipath, r_Cpath = PolicySolve_Analytical_Retired(r_A0,r_Phi0,
#         MAX_AGE , START_AGE, r_DATA, r_PARS   )
#
# # 4. 绘图比较（在之前的图上补充）
# plot(START_AGE:MAX_AGE,r_apath.+r_Phipath,"--")
# legend(["Whole Path","Retired"])







#
