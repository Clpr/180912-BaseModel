push!(LOAD_PATH,"./")
push!(LOAD_PATH,"./modules/")
using PyPlot  # plotting
using Policy_Analytical  # policy analytically
using Policy_Analytical_Retired  # policy analytically specially for the retired phase
using NumAlgo
using PolicyFunctions
using DataFrames  # 用于展示计算的折现福利
# ====================================================
MAX_AGE = 80
RETIRE_AGE = 35
# ========================================================
    Rpath = 0.12 * ones(MAX_AGE)
    Fpath = 0.01 * ones(MAX_AGE); Fpath[end] = 0.0
    Qpath = 0.15 * ones(MAX_AGE)
    Ppath = 1.20 * ones(MAX_AGE)
    cpBpath = 0.30 * ones(MAX_AGE)
    # ----------------------------------------------------
    phiCoefpath = 0.02*ones(RETIRE_AGE)
    Zetapath = 0.06*ones(RETIRE_AGE)
    Etapath = 0.21*ones(RETIRE_AGE)
    Thetapath = 0.05*ones(RETIRE_AGE)
    Wpath = 1.00*ones(RETIRE_AGE)
    doubleApath = 0.30*ones(RETIRE_AGE)
    zpath = 0.90*ones(RETIRE_AGE)
    # ------------------------------
    Lambdapath = 0.02*ones(MAX_AGE-RETIRE_AGE)
    doublePpath = 0.01*ones(MAX_AGE-RETIRE_AGE)
    # ------------------------------
    Mu = 0.1; Sigma = 0.24; Alpha = 1.5; Gamma = 0.5; Varrho = 0.8; Delta = 0.015
# ================================================================
A0 = 0.0; Phi0 = 0.0
# ============================================= 解析解测试
DATA = AmmoReload_DATA( MAX_AGE, Rpath, Fpath, Qpath, Ppath, cpBpath )
DATA_w = AmmoReload_DATA_w( RETIRE_AGE, Wpath, phiCoefpath, Zetapath,
        Etapath, Thetapath, zpath, doubleApath   )
DATA_r = AmmoReload_DATA_r( MAX_AGE, RETIRE_AGE, Lambdapath, doublePpath   )
PARS = AmmoReload_PARS( Mu, Sigma, Alpha, Gamma, Delta  )
# @time apath1, Phipath1, Cpath1, Lpath1 = PolicySolve_Analytical(A0,Phi0, MAX_AGE ,
    # RETIRE_AGE, DATA, DATA_w, DATA_r, PARS   )
# 写一个小函数在计算结果的同时返回求解时间
Timing_Analytical() = begin
    tic()
    apath, Phipath, Cpath, Lpath = PolicySolve_Analytical(A0,Phi0, MAX_AGE , RETIRE_AGE, DATA, DATA_w, DATA_r, PARS   )
    TIME = toc()
    return TIME, apath, Phipath, Cpath, Lpath
end
# 首先跑一次得到求解结果
~, apath0, Phipath0, Cpath0, Lpath0 = Timing_Analytical()
# 然后运行1000次得到求解时间的分布
Sample_Timing_Analytical = zeros(100)
for idx in 1:100
    tmp,~,~,~,~ = Timing_Analytical()
    Sample_Timing_Analytical[idx] = tmp
end
# 记录极端点界限，删除随机造成的离群点
tmpQunt = quantile(Sample_Timing_Analytical,0.98)
filter!( x -> x<tmpQunt, Sample_Timing_Analytical )


# ============================================== 动态规划测试
# 准备数据
DATA = [Rpath::Array{Float64,1},Fpath,Qpath,Ppath,cpBpath]
DATA_w = [phiCoefpath::Array{Float64,1},Zetapath,Etapath,Thetapath,zpath,Wpath,doubleApath]
DATA_r = [Lambdapath::Array{Float64,1},doublePpath]
PARS = [Mu::Float64,Sigma,Alpha,Gamma,Varrho,Delta]
# 定义一个快速函数
Timing_DP(Grid) = begin
    tic()
    apath, Phipath, Cpath, Lpath = PolicySolve(A0, Phi0, MAX_AGE, RETIRE_AGE, DATA, DATA_w, DATA_r, PARS, flag_print=false,
        GRID_LEN=Grid, TOL=1E-6, Arange=[0.0,45.0]);
    TIME = toc()
    return TIME, apath, Phipath, Cpath, Lpath
end
# 准备字典，盛放不同GRID_LEN下的DP的结果以及计时结果
tmpDict = Dict()
GridSet = [ 105, 115, 125, 135 ]
# 循环，先计算一次得到结果存储，然后计算多轮得到计时的统计数据
for idx in 1:length(GridSet)
    # 先计算一次得到当前GRID宽度下的结果
    ~,apath,Phipath,Cpath,Lpath = Timing_DP(GridSet[idx])
    # 再计算多次得到耗时的采样数据
    Sample_Timing_DP = [ Timing_DP(GridSet[idx])[1] for idx2 in 1:100   ]
    # 去掉因随机因素造成的极端值
    tmpQunt = quantile(Sample_Timing_DP,0.98)
    filter!( x -> x<tmpQunt, Sample_Timing_DP )
    # 将结果存储进入字典
    tmpDict[GridSet[idx]] = [ Sample_Timing_DP, apath,Phipath,Cpath,Lpath ]
end



# ================= 同屏绘图(差异性测试)
# 首先DP的结果（因为线会盖住底下的，而我们要突出analytical，所以那个最后画）
for k in GridSet
    subplot(1,3,1)  # 首先画资产
    plot( 1:MAX_AGE, tmpDict[k][2] .+ tmpDict[k][3] )
    subplot(1,3,2)  # 然后画消费
    plot( 1:MAX_AGE, tmpDict[k][4] )
    subplot(1,3,3)  # 最后画闲暇
    plot( 1:RETIRE_AGE, tmpDict[k][5][1:RETIRE_AGE] )
end
# 准备legend向量
tmpLegend = [ "DP Density="*string(x) for x in GridSet ]
push!(tmpLegend,"Analytical")
# 然后画analytical的结果，顺便加上图例等元素
    subplot(1,3,1)  # 首先画资产
    plot( 1:MAX_AGE, apath0 .+ Phipath0 )
        grid(true)
        xlabel("Age"); ylabel("Wealth (Capital)"); #legend(tmpLegend,loc="upper left")
    subplot(1,3,2)  # 然后画消费
    plot( 1:MAX_AGE, Cpath0 )
        grid(true)
        xlabel("Age"); ylabel("Consumption"); legend(tmpLegend,loc="upper left")
    subplot(1,3,3)  # 最后画闲暇
    plot( 1:RETIRE_AGE, Lpath0[1:RETIRE_AGE] )
        grid(true)
        xlabel("Age"); ylabel("Leisure"); #legend(tmpLegend,loc="upper left")


# ====================== 计算福利
# 定义一个小函数，用于快速计算折现效用，使用Policy_Analytical中的U函数
DiscU(Cpath,Lpath) = begin
    tmpDisc = cumprod( fill(1/(1+Delta),MAX_AGE) )  # 效用折现因子
    tmpUvec = 0.0
    for idx in 1:length(Cpath)
        tmpUvec += Policy_Analytical.U(Cpath[idx],Lpath[idx],Qpath[idx],Gamma,Alpha) * tmpDisc[idx]
    end
    return tmpUvec
end
# 为每种情形计算福利(DP)
tmpDict_U = Dict()
for idx in 1:length(GridSet)
    tmpKey = "DP Density="*string(GridSet[idx])
    # 提取C，L，组装完整的L
    tmpC = tmpDict[GridSet[idx]][4]
    tmpL = tmpDict[GridSet[idx]][5]
    append!(tmpL,ones(MAX_AGE-RETIRE_AGE))
    tmpVal = DiscU( tmpC, tmpL )
    # load to dictionary
    tmpDict_U[tmpKey] = tmpVal
end
# 手动计入解析解的情形
tmpC = Cpath0; tmpL = Lpath0
append!(tmpL,ones(MAX_AGE-RETIRE_AGE))
tmpDict_U["Analytical"] = DiscU( tmpC, tmpL )
# 制成数据框
tab_U = DataFrame(tmpDict_U)
print(tab_U)




# ================= 同屏绘图(性能比较测试)
# 直接画每种case的耗时的分布
# 首先DP的结果
subplot(2,3,1)  # GRID=105
plt[:hist]( tmpDict[105][1] ,density=true)
    grid(true); title("DP Density = 105")
    xlabel("Time Cost"); ylabel("Density"); #legend(tmpLegend,loc="upper left")
subplot(2,3,2)  # GRID=115
plt[:hist]( tmpDict[115][1] ,density=true)
    grid(true); title("DP Density = 115")
    xlabel("Time Cost"); ylabel("Density"); #legend(tmpLegend,loc="upper left")
subplot(2,3,3)  # GRID=125
plt[:hist]( tmpDict[125][1] ,density=true)
    grid(true); title("DP Density = 125")
    xlabel("Time Cost"); ylabel("Density"); #legend(tmpLegend,loc="upper left")
    subplot(2,3,4)  # GRID=135
plt[:hist]( tmpDict[135][1] ,density=true)
    grid(true); title("DP Density = 135")
    xlabel("Time Cost"); ylabel("Density"); #legend(tmpLegend,loc="upper left")
subplot(2,3,5)  # Analytical
plt[:hist]( Sample_Timing_Analytical ,density=true)
    grid(true); title("Adjusted Analytical")
    xlabel("Time Cost"); ylabel("Density"); #legend(tmpLegend,loc="upper left")






# ==================== 测试折现效用&时间成本随DP griding密度增加的规律，绘制函数
# 定义一个快速函数，用于计算给定grid的DP求解的时间和折现效用
func_DPtest_U_Time(Grid) = begin
    TIME, apath, Phipath, Cpath, Lpath = Timing_DP(Grid) # 求解
    append!(Lpath,ones(MAX_AGE-RETIRE_AGE))  # 拓展完整闲暇
    tmpU = DiscU(Cpath,Lpath)  # 计算折现效用
    return TIME,tmpU
end
# 准备grid序列
GridSet2 = Array(100:20:500)
# 准备字典
tmpDict_DPtest = Dict(
    "DP Density" => GridSet2,
    "Discounted Lifetime Utility" => zeros(length(GridSet2)),
    "Time Cost (s)" => zeros(length(GridSet2))
)
# 开始循环填充
for idx in 1:length(GridSet2)
    tmpDict_DPtest["Time Cost (s)"][idx],tmpDict_DPtest["Discounted Lifetime Utility"][idx] = func_DPtest_U_Time(GridSet2[idx])
end
# 制表，输出
tab_DPtest = DataFrame(tmpDict_DPtest)
writetable("./ExtraScripts/DPtest结果.csv",tab_DPtest)
# 绘图
plot( GridSet2, tmpDict_DPtest["Discounted Lifetime Utility"] )
    grid(true)
    xlabel("DP Density"); ylabel("Discounted Lifetime Utility")
twinx()
plot( GridSet2, tmpDict_DPtest["Time Cost (s)"], "r" )
    ylabel("Time Cost (s)")
legend(["Time Cost (s)"],loc="best")























#
