# PROC: prepare memories for economic variables
# (only called from main.jl when MAX_YEAR, MAX_AGE & RETIRE_AGE defined)
# ==============================================

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











#
