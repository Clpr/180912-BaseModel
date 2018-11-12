# PROC: initialize parameters for easy-modifying; more convenient for calibarations
# (only called from main.jl when MAX_YEAR, MAX_AGE & RETIRE_AGE defined)
# =================================================================================

# =================================================================================
# PART 1: Most-frequently used & modified in calibration
# =================================================================================
# 1.1 Technology (TFP) 技术系数（TFP）
    # NOTE: 年与坐标对应关系: idx = year - 1945
    vec_tech = ones(MAX_YEAR)
    # 0. 时间节点
        pt_1 = 35  # 1980
        pt_2 = 53  # 2008
        pt_3 = 73  # 2018
    # 1. 第1段，到1980年之前(35年)
        vec_tech[1:pt_1] = 1.0 .* 1.01 .^ (0:pt_1-1)
    # 2. 第2段，从1980到2008年(28年)
        Calib_TechGrow = linspace(0.04,0.045,pt_2-pt_1)
        vec_tech[pt_1+1:pt_2] = vec_tech[pt_1] .* cumprod(1 .+ Calib_TechGrow)
    # 3. 第3段，2008到2018年(共10年，长度为10)
        Calib_TechGrow = linspace(0.05,0.011,pt_3-pt_2)
        vec_tech[pt_2+1:pt_3] = vec_tech[pt_2] .* cumprod(1 .+ Calib_TechGrow)
    # 4. 第4段，2018年之后(增长50年)
        TechGrowYear = 50
        vec_tech[pt_3:pt_3+TechGrowYear] = vec_tech[pt_3] * 1.01 .^ (0:TechGrowYear)
        vec_tech[pt_3+TechGrowYear+1:end] = vec_tech[pt_3+TechGrowYear]



# 1.2 inter-temporal consumption substitute elasticity 跨期替代弹性
    num_gamma = 0.5

# 1.3 cap of government debt ratio (D/Y) 政府债务上限约束
    vec_doubleK = 0.00 * ones(MAX_YEAR)

# 1.4 capital income share 资本收入占比
    # NOTE: (maybe 0.4~0.6 for China, but finally calibrated at 0.55 according to literature)
    vec_beta = 0.55 * ones(MAX_YEAR)

# =================================================================================
# PART 2: Preferences, micro parameters & non-direct parameters
# =================================================================================
# 2.1 depreciation rate 折旧率
    num_kappa = 0.05

# 2.2 utility discounting rate 效用折现率
    num_delta = 0.01

# 2.3 leisure preference than consumption 闲暇相对消费的偏好
    num_alpha = 1.5

# 2.5 consumption tax 消费税
    num_mu = 0.1

# 2.6 wage tax 工资税
    num_sigma = 0.24

# 2.7 Collection rate of pension 养老金收缴率
    vec_z = 0.85 * ones(MAX_YEAR)

# 2.8 co-payment rate: inpatient 住院费用的自付比例
    vec_cpB = 0.3 * ones(MAX_YEAR)

# 2.9 contribution: firm -> pension 缴纳：企业 to 养老金
    vec_eta = 0.2 * ones(MAX_YEAR)

# 2.10 contribution: individual -> pension 缴纳：个人 to 养老金
    vec_theta = 0.08 * ones(MAX_YEAR)

# 2.11 contribution: firm -> medical 缴纳：企业 to 医保
    vec_zeta = 0.06 * ones(MAX_YEAR)

# 2.12 contribution: individual -> medical 缴纳：个人 to 养老金
    vec_phiCoef = 0.02 * ones(MAX_YEAR)

# 2.13 transfer: firm.medical -> individual accounts (just the contributors, working phase) 转移支付：医保企业缴纳中分给个人账户（工作期）的比例
    vec_doubleA = 0.3 * ones(MAX_YEAR)

# 2.14 transfer rate from firm contribution to medical to individual medical account (other alive people, retired phase) 转移支付：当年医保企业缴纳中分给当年存活退休人群个人账户的比例
    vec_doubleB = 0.0 * ones(MAX_YEAR)

# =================================================================================
# PART 3: Deprecated
# =================================================================================
# inter-temporal substitute elasticity (DEPRECATED)
    num_varrho = 0.80

# elasticity of total medical expenses to gdp (with mark 'x') (DEPRECATED)
    # \Delta M / M * Y / \Delta Y
    vec_x = 1.6 * ones(MAX_YEAR)

# Compute total contributions on nominal wages
# # pi (pension)
    # vec_pi = ( vec_eta + vec_theta ) ./ ( 1 + vec_eta + vec_zeta )
# # pi^M (medical)
    # vec_piM = ( vec_zeta + vec_phi ) ./ ( 1 + vec_eta + vec_zeta )
# # co-payment rate: outpatient
    # vec_cpA = 0.0 * ones(MAX_YEAR)
# # m2c ratio ( q ) the ratio of total medical expenses to consumption
    # mat_m2c = ones(MAX_YEAR,MAX_AGE)







#
