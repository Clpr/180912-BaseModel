__precompile__()
"""
    PolicyFunctions

Solves policy functions of personal asset (a) and individual-medical account (Phi).
Also contains Bellman functions.

"""
module PolicyFunctions
    using NumAlgo
    export PolicySolve, PolicySolve_r
# ========================== MOD BEGIN =========================================
    """
        U(C::Float64, L::Float64, Q::Float64, Alpha::Float64, Gamma::Float64, Varrho::Float64)

    Cross-sectional utility of the CES utility function

    Pars:
        1. C [num]: NET consumption, must be greater or equal to 0
        1. L [num]: leisure, must be in the range of [0,1]
        1. Q [num]: q ratio
        1. Alpha [num]: Leisure-consumption preference coef, usually be 1.5
        1. Gamma [num]: Inter-temporal consumption substitute elasticity, usually be various (0.25 in AK1987)
        1. Varrho [num]: Inter-temporal substitute elasticity, usually be 0.8

    Ret:
        1. ucrux: Cross-sectional utility value

    Depend:
        None

    NOTE:
        1. No validation for 0<=Q<=1, pls do it outside the function
        1. The function supports vector or matrix inputs of C & L, but Alpha, Gamma, Varrho MUST be scalar
        1. Two minor scalar are added to avoid Inf or -Inf when c,l are very close to zero
        1. The function corresponds to the CES utility function in 'Dynamic Fiscal Policy' by Auerbach and Kotlikoff, 1987
        1. No validation performed
    """
    function U(C::Float64, L::Float64, Q::Float64, Alpha::Float64, Gamma::Float64, Varrho::Float64)
        # ----------------------
        ucrux = (1-1/Gamma)*(  (  (1-Q)*C+(1E-03))^(1-1/Varrho) + Alpha*(L+(1E-03))^(1-1/Varrho)  )^( (1-1/Gamma)/(1-1/Varrho) )
        # Scaling the number out of beauty
        ucrux /= 1E+04
        return ucrux::Float64
    end
    # ==========================================================================
    """
        GetCL_r(Athis::Float64,Anext::Float64,PARS::Array{Float64,1})

    Compute c(s) for the retired when given A(s),A(s+1) and other data.
    where PARS contains (in order):
    1. R : interest rate
    1. F : mortality
    1. Q : m2c ratio
    1. P : outpatient / inpatient fee
    1. cpB : copayment ratio of inpatient
    1. Mu : consumption tax rate
    1. Lambda : pension benefit amount
    1. doubleP : post-retire contribution by firms (amount)
    """
    function GetCL_r(Athis::Float64,Anext::Float64,PARS::Array{Float64,1})
        # Unfold PARS
        R,F,Q,P,cpB,Mu,Lambda,doubleP = PARS
        # Compute short-write
        # here we add a minor const 1E-10 to avoid NaN in the very last year
        # (all left agents die at the death moment, but Anext = 0, so an nearly-infinite coefficient wont matter)
        Right = (1+R)*Athis - (2-1/(1-F+(1E-10)))*Anext + Lambda + doubleP
        Left = (1+Mu) - (1-cpB)/(1+P)*Q
        # Get C (non-negative check)
        Cthis = Right/Left < 0.0 ? -999.0 : Right/Left
        return Cthis::Float64
    end
    # ==========================================================================
    """
        GetCL_w(Athis::Float64,Anext::Float64,PARS::Array{Float64,1})

    Compute c(s) for the retired when given A(s),A(s+1) and other data.
    where PARS contains (in order):
    1. R : interest rate
    1. F : mortality
    1. Q : m2c ratio
    1. P : outpatient / inpatient fee
    1. cpB : copayment ratio of inpatient
    1. Mu : consumption tax rate
    1. Sigma : wage tax rate
    1. phiCoef : personal contribution to medical
    1. Zeta : firm contribution to medical
    1. Eta : firm contribution to pension
    1. Theta : personal contribution to pension
    1. z : collection rate of pension
    1. W : wage level
    1. doubleA : transfer rate from firm medical contribution amount to individual medical account
    1. Alpha : leisure preference than consumption
    1. Varrho : consumption-leisure substitution elasticity
    """
    function GetCL_w(Athis::Float64,Anext::Float64,PARS::Array{Float64,1})
        # unfold PARS
        R,F,Q,P,cpB,Mu,Sigma,phiCoef,Zeta,Eta,Theta,z,W,doubleA,Alpha,Varrho = PARS
        # Compute short writings
        Pi = z*(Theta+Eta)/(1+z*Eta+Zeta); PiM = (phiCoef+Zeta)/(1+z*Eta+Zeta)
        tmp_A = 1+Mu - (1-cpB)/(1+P)*Q
        tmp_B = (1-Sigma-Pi-PiM+ (phiCoef+doubleA*Zeta)/(1+z*Eta+Zeta) ) * W
        tmp_C = (1-Q)^(1-1/Varrho) / Alpha
        tmp_D = (1+R)*Athis - (2-1/(1-F+(1E-10)))*Anext
        # Get leisure
        Right_L = tmp_D + tmp_B
        Left_L = tmp_B + tmp_A * (tmp_A/tmp_B/tmp_C)^(1/(1-1/Varrho))
        Lthis = Right_L / Left_L
        # Get consumption
        Right_C = (tmp_A/tmp_B/tmp_C)^(1/(1-1/Varrho) )
        Cthis = Lthis * Right_C
        # Validation
        if Lthis < 0.0
            Lthis = 0.0
            Cthis = (tmp_D+tmp_B)/tmp_A
            if Cthis < 0.0 # CASE: no solution
                Cthis = -1.0
                return Cthis::Float64, Lthis::Float64
            end
            return Cthis::Float64, Lthis::Float64
        elseif Lthis>1.0
            Lthis = 1.0
            Cthis = tmp_D/tmp_A
            if Cthis < 0.0 # CASE: no solution
                Cthis = -1.0
                return Cthis::Float64, Lthis::Float64
            end
            return Cthis::Float64, Lthis::Float64
        end

        # CASE: everything is fine
        return Cthis::Float64, Lthis::Float64
    end
    # ==========================================================================
    """
        Bellman_dead(Athis::Float64, PARS::Array{Float64,1})

    Bellman equation for those at the beginning of the very last age year.

    PARS (in order):
    1. R : interest rate
    1. F : mortality
    1. Q : m2c ratio
    1. P : outpatient / inpatient fee
    1. cpB : copayment ratio of inpatient
    1. Mu : consumption tax rate
    1. Lambda : pension benefit amount
    1. doubleP : post-retire contribution by firms (amount)
    1. Alpha : leisure preference than consumption
    1. Gamma: inter-temporal substitution elasticity
    1. Varrho: consumption-leisure substitution elasticity
    """
    function Bellman_dead(Athis::Float64, PARS::Array{Float64,1})
        # Unfold PARS
        R,F,Q,P,cpB,Mu,Lambda,doubleP,Alpha,Gamma,Varrho = PARS
        # Fold a PARS for consumption computation
        PARS_C = [R,F,Q,P,cpB,Mu,Lambda,doubleP]
        # Compute consumption in the very last year
        Cthis = GetCL_r(Athis::Float64,0.0::Float64,PARS_C::Array{Float64,1})
        # Validation
        if Cthis < 0.0
            BVal = -6.66E20
        else
            BVal = U(Cthis, 1.0, Q, Alpha, Gamma, Varrho)
        end
        return BVal::Float64
    end
    # ==========================================================================
    """
        Bellman_r(Athis::Float64,Anext::Float64,AnextSeries::Array{Float64,1},VnextSeries::Array{Float64,1},PARS::Array{Float64,1})

    Bellman equation for the retired (the very last age year not included).

    PARS (in order):
    1. R : interest rate
    1. F : mortality
    1. Q : m2c ratio
    1. P : outpatient / inpatient fee
    1. cpB : copayment ratio of inpatient
    1. Mu : consumption tax rate
    1. Lambda : pension benefit amount
    1. doubleP : post-retire contribution by firms (amount)
    1. Alpha : leisure preference than consumption
    1. Gamma: inter-temporal substitution elasticity
    1. Varrho: consumption-leisure substitution elasticity
    """
    function Bellman_r(Athis::Float64,Anext::Float64,
        AnextSeries::Array{Float64,1},VnextSeries::Array{Float64,1},PARS::Array{Float64,1})
        # -------------------------
        # Unfold PARS
        R,F,Q,P,cpB,Mu,Lambda,doubleP,Alpha,Gamma,Varrho,Delta = PARS
        # Fold a PARS for consumption computation
        PARS_C = [R,F,Q,P,cpB,Mu,Lambda,doubleP]
        # Compute consumption in the very last year
        Cthis = GetCL_r(Athis::Float64,Anext::Float64,PARS_C::Array{Float64,1})
        # Validation
        if Cthis < 0.0
            BVal = -6.66E20
        else  # compute Bellman equation with linear interpolating
            # Discount Factor
            Beta = (1-F+(1E-10))/(1+Delta)
            BVal = U(Cthis, 1.0, Q, Alpha, Gamma, Varrho) + Beta*VFunc(Anext,AnextSeries,VnextSeries)
        end
        return BVal::Float64
    end
    # ==========================================================================
    """
        Bellman_w(Athis::Float64,Anext::Float64,AnextSeries::Array{Float64,1},VnextSeries::Array{Float64,1},PARS::Array{Float64,1})

    Bellman equation for the working agents.

    PARS (in order):
    1. R : interest rate
    1. F : mortality
    1. Q : m2c ratio
    1. P : outpatient / inpatient fee
    1. cpB : copayment ratio of inpatient
    1. Mu : consumption tax rate
    1. Sigma : wage tax rate
    1. phiCoef : personal contribution to medical
    1. Zeta : firm contribution to medical
    1. Eta : firm contribution to pension
    1. Theta : personal contribution to pension
    1. z : collection rate of pension
    1. W : wage level
    1. doubleA : transfer rate from firm medical contribution amount to individual medical account
    1. Alpha : leisure preference than consumption
    1. Varrho: consumption-leisure substitution elasticity
    1. Gamma: inter-temporal substitution elasticity
    """
    function Bellman_w(Athis::Float64,Anext::Float64,
        AnextSeries::Array{Float64,1},VnextSeries::Array{Float64,1},PARS::Array{Float64,1})
        # ----------------------
        # unfold PARS
        R,F,Q,P,cpB,Mu,Sigma,phiCoef,Zeta,Eta,Theta,z,W,doubleA,Alpha,Varrho,Delta,Gamma = PARS
        # fold a PARS for consumption, leisure computation
        PARS_CL = [R,F,Q,P,cpB,Mu,Sigma,phiCoef,Zeta,Eta,Theta,z,W,doubleA,Alpha,Varrho]
        # Get consumption & leisure
        Cthis,Lthis = GetCL_w(Athis::Float64,Anext::Float64,PARS_CL::Array{Float64,1})
        # Validation
        if Cthis < 0.0
            BVal = -6.66E20
        else  # compute Bellman equation with linear interpolating
            # Discount Factor
            Beta = (1-F+(1E-10))/(1+Delta)
            BVal = U(Cthis, 1.0, Q, Alpha, Gamma, Varrho) + Beta*VFunc(Anext,AnextSeries,VnextSeries)
        end
        return BVal::Float64
    end
    # ==========================================================================
    """
        PolicySolve(A0::Float64,Phi0::Float64,DATA::Array{Array{Float64,1},1},DATA_w::Array{Array{Float64,1},1},DATA_r::Array{Array{Float64,1},1}, PARS::Array{Float64,1} ;  flag_print::Bool=false, GRID_LEN::Int64=50, TOL::Float64=1E-8, Arange::Array{Float64,1}=[0.0,15.0])

    Solve optimal personal asset path, individual medical account path, consumption path and leisure path.

    In:
    1. A0 : initial personal asset (beginning of year)
    1. Phi0 : initial individual medical account (beginning of year)  (A0+Phi0 then we can get the initial personal wealth)
    1. MAX_AGE : maximum age (for indexing)
    1. RETIRE_AGE : retiring age (for indexing)
    1. DATA : a collect of vector which has a common length MAX_AGE
    1. DATA_w : a collect of vector which has a common length RETIRE_AGE
    1. DATA_r : a collect of vector which has a common length (MAX_AGE - RETIRE_AGE)
    1. PARS : a collect of scalars as parameters
    1. flag_print : whether to print timing information
    1. GRID_LEN : density of grid search
    1. TOL : tolerance of convergence

    Out:
    1. Apath::Array{Float64,1} : personal asset path, length = MAX_AGE
    1. Phipath::Array{Float64,1} : individual medical account path, length = MAX_AGE
    1. Cthis::Array{Float64,1} : consumption path, length = MAX_AGE
    1. Lthis::Array{Float64,1} : leisure path, length = RETIRE_AGE

    **EXTRA**
    **DATA** contains (in order):
    1. Rpath::Array{Float64,1} : interest rate path
    1. Fpath::Array{Float64,1} : mortality path
    1. Qpath::Array{Float64,1} : m2c ratio path
    1. Ppath::Array{Float64,1} : outpatient / inpatient fee path
    1. cpBpath::Array{Float64,1} : inpatient copayment rate path

    **DATA_w** contains (in order):
    1. phiCoefpath::Array{Float64,1} : personal contribution rate to medical
    1. Zetapath::Array{Float64,1} : firm contribution rate to medical
    1. Etapath::Array{Float64,1} : firm contribution rate to pension
    1. Thetapath::Array{Float64,1} : personal contribution rate to pension
    1. zpath::Array{Float64,1} : collection rate of pension
    1. Wpath::Array{Float64,1} : wage level path
    1. doubleApath::Array{Float64,1} : transfer (amount) from firm contribution to individual medical account
x
    **DATA_r** contains(in order):
    1. Lambdapath::Array{Float64,1} : pension benefits
    1. doublePpath::Array{Float64,1} : transfer (amount) from firm contribution to those retired

    **PARS** contains(in order):
    1. Mu : consumption tax rate
    1. Sigma : wage tax rate
    1. Alpha : leisure preference than consumption
    1. Gamma : inter-temporal substitution elasticity
    1. Varrho : consumption leisure substitution elasticity
    1. Delta : utility discount rate
    """
    function PolicySolve(A0::Float64,Phi0::Float64, MAX_AGE::Int64, RETIRE_AGE::Int64,
        DATA::Array{Array{Float64,1},1},DATA_w::Array{Array{Float64,1},1},DATA_r::Array{Array{Float64,1},1}, PARS::Array{Float64,1} ;
        flag_print::Bool=false, GRID_LEN::Int64=50, TOL::Float64=1E-8, Arange::Array{Float64,1}=[0.0,15.0])
        # ------------------------------
        # Timing
        time_begin = time()
        # Unfold data
        Rpath::Array{Float64,1},Fpath,Qpath,Ppath,cpBpath= DATA
        phiCoefpath::Array{Float64,1},Zetapath,Etapath,Thetapath,zpath,Wpath,doubleApath = DATA_w
        Lambdapath::Array{Float64,1},doublePpath = DATA_r
        Mu::Float64,Sigma,Alpha,Gamma,Varrho,Delta = PARS
        # Validation
        @assert(2==length(Arange),"Arange can only have two elements in PolicySolve()")
        @assert(MAX_AGE>RETIRE_AGE,"MAX_AGE must be greater than RETIRE_AGE in PolicySolve()")
        @assert(MAX_AGE==length(Rpath)==length(Fpath)==length(Qpath)==length(Ppath)==length(cpBpath),"Uncompatible MAX_AGE-length inputs in PolicySolve()")
        @assert(RETIRE_AGE==length(Wpath)==length(phiCoefpath)==length(Zetapath)==length(Etapath)==length(Thetapath)==length(doubleApath)==length(Etapath),"Uncompatible RETIRE_AGE-length inputs in PolicySolve()")
        @assert((MAX_AGE-RETIRE_AGE)==length(Lambdapath)==length(doublePpath),"Uncompatible Lambdapath in PolicySolve()")
        # Algorithm data spacing
        Aopt = zeros(GRID_LEN,MAX_AGE)  # CAUTION: Aopt here is the summation of personal asset and individual medical account!!!
        Copt = zeros(GRID_LEN,MAX_AGE); Lopt = zeros(GRID_LEN,RETIRE_AGE)
        Vr = zeros(GRID_LEN,MAX_AGE)   # Value function matrix
        Grid_A = Array(linspace(Arange[1],Arange[2],GRID_LEN))  # Griding vector of personal asset
        Apath = zeros(MAX_AGE)  # CAUTION: Aopt here is the summation of personal asset and individual medical account!!!
        apath = zeros(MAX_AGE)  # CAUTION: apath here is ONLY the personal asset!!!
        Phipath = zeros(MAX_AGE)  # CAUTION: Phipath here is ONLY the individual medical account!!!
        Cpath = zeros(MAX_AGE); Lpath = zeros(RETIRE_AGE)  # Final results
        # ----------------------
        # ----------------------
        # ------- BELLMAN: s=S ---------------
        idx_Age = MAX_AGE
        # Value function for s=S
        Aopt[:,idx_Age] = Grid_A
        PARS_DEAD = [  Rpath[idx_Age], Fpath[idx_Age], Qpath[idx_Age], Ppath[idx_Age], cpBpath[idx_Age],
            Mu, Lambdapath[idx_Age-RETIRE_AGE], doublePpath[idx_Age-RETIRE_AGE], Alpha, Gamma, Varrho, Delta  ]
        Vr[:,idx_Age] = [Bellman_dead( Grid_A[x], PARS_DEAD ) for x in 1:GRID_LEN]
        # Get optimal consumption
        Copt[:,idx_Age] = [GetCL_r( Grid_A[x], 0.0, PARS_DEAD[1:8] ) for x in 1:GRID_LEN]

        tmpV = zeros(GRID_LEN)  # a temporary value-function vector outside the loop to avoid repeatlly constructing
        # ------- BELLMAN: Sr<s<S ------------
        for idx_Age in MAX_AGE-1:-1:RETIRE_AGE+1 , idx_Grid in 1:GRID_LEN
                # ammo reloaded
                PARS_R = [  Rpath[idx_Age], Fpath[idx_Age], Qpath[idx_Age], Ppath[idx_Age], cpBpath[idx_Age],
                    Mu, Lambdapath[idx_Age-RETIRE_AGE], doublePpath[idx_Age-RETIRE_AGE], Alpha, Gamma, Varrho, Delta  ]
                afunc_Bellman_r(Anext::Float64) = Bellman_r(Grid_A[idx_Grid], Anext::Float64 ,Grid_A,Vr[:,idx_Age+1],PARS_R)
                # compute temporary value function for griding points
                for x in 1:GRID_LEN
                    tmpV[x] = afunc_Bellman_r(Grid_A[x])
                end
                # tmpV = [afunc_Bellman_r(Grid_A[x]) for x in 1:GRID_LEN  ]
                # find the peak
                ~,tmp_LocOptAs = findmax(tmpV)
                if tmp_LocOptAs==1
                    Aopt[idx_Grid,idx_Age] = Arange[1]
                elseif tmp_LocOptAs==GRID_LEN
                    Aopt[idx_Grid,idx_Age] = Arange[2]
                else
                    gold_a = Grid_A[tmp_LocOptAs-1] # pars for Golden Section Searching
                    gold_b = Grid_A[tmp_LocOptAs]
                    gold_c = Grid_A[tmp_LocOptAs+1]
                    # Golden Section Searching
                    Aopt[idx_Grid,idx_Age] = Golden(afunc_Bellman_r::Function, gold_a::Float64,gold_b::Float64,gold_c::Float64,TOL=TOL  )
                end
                # Get interpolated value function & consumption
                Vr[idx_Grid,idx_Age] = afunc_Bellman_r(Aopt[idx_Grid,idx_Age])
                Copt[idx_Grid,idx_Age] = GetCL_r( Grid_A[idx_Grid], Aopt[idx_Grid,idx_Age], PARS_R[1:8] )
        end
        # ----------- BELLMAN: 1<=s<=Sr (working phase) ------------
        for idx_Age in RETIRE_AGE:-1:1  ,  idx_Grid in 1:GRID_LEN
                PARS_W = [ Rpath[idx_Age], Fpath[idx_Age], Qpath[idx_Age], Ppath[idx_Age], cpBpath[idx_Age],
                    Mu, Sigma, phiCoefpath[idx_Age], Zetapath[idx_Age], Etapath[idx_Age],
                    Thetapath[idx_Age], zpath[idx_Age], Wpath[idx_Age], doubleApath[idx_Age], Alpha, Varrho, Delta, Gamma ]
                afunc_Bellman_w(Anext::Float64) = Bellman_w(Grid_A[idx_Grid],Anext::Float64,Grid_A,Vr[:,idx_Age+1],PARS_W)
                # compute temporary value function for griding points
                for x in 1:GRID_LEN
                    tmpV[x] = afunc_Bellman_w(Grid_A[x])
                end
                # tmpV = [afunc_Bellman_w(Grid_A[x]) for x in 1:GRID_LEN]
                # find the peak
                ~,tmp_LocOptAs = findmax(tmpV)
                if tmp_LocOptAs==1
                    Aopt[idx_Grid,idx_Age] = Arange[1]
                elseif tmp_LocOptAs==GRID_LEN
                    Aopt[idx_Grid,idx_Age] = Arange[2]
                else
                    gold_a = Grid_A[tmp_LocOptAs-1] # pars for Golden Section Searching
                    gold_b = Grid_A[tmp_LocOptAs]
                    gold_c = Grid_A[tmp_LocOptAs+1]
                    # Golden Section Searching
                    Aopt[idx_Grid,idx_Age] = Golden(afunc_Bellman_w::Function, gold_a::Float64,gold_b::Float64,gold_c::Float64,TOL=TOL  )
                end
                # Get interpolated value function & consumption
                Vr[idx_Grid,idx_Age] = afunc_Bellman_w(Aopt[idx_Grid,idx_Age])
                Copt[idx_Grid,idx_Age],Lopt[idx_Grid,idx_Age] = GetCL_w( Grid_A[idx_Grid], Aopt[idx_Grid,idx_Age], PARS_W[1:16] )
        end
        # ----------- EXTRACT OPTIMAL a+Phi path -------------------
        Apath[1] = A0 + Phi0
        # look back, interpolating & compute C,L
        for idx_Age in 1:RETIRE_AGE
            # uses computed value function & interpolates asset in the next year
            Apath[idx_Age+1] = APolate(Apath[idx_Age]::Float64, Aopt[:,idx_Age]::Array{Float64,1}, Arange::Array{Float64,1}  )
            PARS_W = [ Rpath[idx_Age], Fpath[idx_Age], Qpath[idx_Age], Ppath[idx_Age], cpBpath[idx_Age],
                Mu, Sigma, phiCoefpath[idx_Age], Zetapath[idx_Age], Etapath[idx_Age],
                Thetapath[idx_Age], zpath[idx_Age], Wpath[idx_Age], doubleApath[idx_Age], Alpha, Varrho ]
            Cpath[idx_Age],Lpath[idx_Age] = GetCL_w(Apath[idx_Age],Apath[idx_Age+1], PARS_W )
        end
        if MAX_AGE>RETIRE_AGE+1
            # cauz if MAX_AGE==RETIRE_AGE+1, no need to run the following loop
            for idx_Age in RETIRE_AGE+1:MAX_AGE-1
                PARS_R = [  Rpath[idx_Age], Fpath[idx_Age], Qpath[idx_Age], Ppath[idx_Age], cpBpath[idx_Age],
                    Mu, Lambdapath[idx_Age-RETIRE_AGE], doublePpath[idx_Age-RETIRE_AGE] ]
                Apath[idx_Age+1] = APolate(Apath[idx_Age]::Float64, Aopt[:,idx_Age]::Array{Float64,1}, Arange::Array{Float64,1}  )
                Cpath[idx_Age] = GetCL_r(Apath[idx_Age],Apath[idx_Age+1],PARS_R)
            end
        end
        # Get the consumption in the very last year
        idx_Age = MAX_AGE
        PARS_R = [  Rpath[idx_Age], Fpath[idx_Age], Qpath[idx_Age], Ppath[idx_Age], cpBpath[idx_Age],
            Mu, Lambdapath[idx_Age-RETIRE_AGE], doublePpath[idx_Age-RETIRE_AGE] ]
        Cpath[idx_Age] = GetCL_r(Apath[idx_Age],0.0::Float64,PARS_R)

        # ----- DEPART INDIVIDUAL-MEDICAL ACCOUNT FROM THE COMPREHENSIVE ASSET ACCOUNT ---------
        Phipath,~ = Depart_Phipath(Phi0,Cpath,Lpath, MAX_AGE, RETIRE_AGE, DATA,DATA_w,DATA_r, PARS )
        # get personal asset path
        apath = Apath .- Phipath

        # Timing ends
        if flag_print
            println("(Policy) Time Elasped: ",time()-time_begin," s")
        end

        return apath::Array{Float64,1},Phipath::Array{Float64,1},Cpath::Array{Float64,1},Lpath::Array{Float64,1}
    end
    # ==========================================================================
    """
        Depart_Phipath(Phi0::Float64,Cpath::Array{Float64,1},Lpath::Array{Float64,1}, MAX_AGE::Int64, RETIRE_AGE::Int64, DATA::Array{Array{Float64,1},1},DATA_w::Array{Array{Float64,1},1},DATA_r::Array{Array{Float64,1},1}, PARS::Array{Float64,1} )

    Given optimal consumption, leisure paths, computes individual medical account path.
    Most parameters are the same as PolicySolve()

    return a Phipath and a Gappath, cauz we have the constraint: individual medical account is non-negative
    """
    function Depart_Phipath(Phi0::Float64,Cpath::Array{Float64,1},Lpath::Array{Float64,1}, MAX_AGE::Int64, RETIRE_AGE::Int64,
        DATA::Array{Array{Float64,1},1},DATA_w::Array{Array{Float64,1},1},DATA_r::Array{Array{Float64,1},1}, PARS::Array{Float64,1} )
        # ----------------------
        Rpath::Array{Float64,1},Fpath,Qpath,Ppath,cpBpath= DATA
        phiCoefpath::Array{Float64,1},Zetapath,Etapath,Thetapath,zpath,Wpath,doubleApath = DATA_w
        Lambdapath::Array{Float64,1},doublePpath = DATA_r
        Mu::Float64,Sigma,Alpha,Gamma,Varrho,Delta = PARS
        # ------------------------
        Phipath = zeros(MAX_AGE); Gappath = zeros(MAX_AGE)
        # for working phase
        Phipath[1] = Phi0
        for idx_Age in 1:RETIRE_AGE
            Phipath[idx_Age+1] = (1+Rpath[idx_Age])*Phipath[idx_Age] + (phiCoefpath[idx_Age]+doubleApath[idx_Age]*Zetapath[idx_Age])/(1+zpath[idx_Age]*Etapath[idx_Age]+Zetapath[idx_Age])*Wpath[idx_Age]*(1.0-Lpath[idx_Age])
            Phipath[idx_Age+1] -= Qpath[idx_Age]*Ppath[idx_Age]/(1+Ppath[idx_Age])*Cpath[idx_Age]
            Phipath[idx_Age+1] /= 2-1/(1-Fpath[idx_Age])
        end
        if MAX_AGE>RETIRE_AGE+1
            # cauz if MAX_AGE==RETIRE_AGE+1, no need to run the following loop
            for idx_Age in RETIRE_AGE+1:MAX_AGE-1
                Phipath[idx_Age+1] = (1+Rpath[idx_Age])*Phipath[idx_Age] + doublePpath[idx_Age-RETIRE_AGE]
                Phipath[idx_Age+1] -= Qpath[idx_Age]*Ppath[idx_Age]/(1+Ppath[idx_Age])*Cpath[idx_Age]
                Phipath[idx_Age+1] /= 2-1/(1-Fpath[idx_Age])
            end
        else
            idx_Age = RETIRE_AGE+1
            Phipath[idx_Age+1] = (1+Rpath[idx_Age])*Phipath[idx_Age] + doublePpath[idx_Age-RETIRE_AGE]
            Phipath[idx_Age+1] -= Qpath[idx_Age]*Ppath[idx_Age]/(1+Ppath[idx_Age])*Cpath[idx_Age]
            Phipath[idx_Age+1] /= 2-1/(1-Fpath[idx_Age])
        end
        # ----------- Compute Gaps -----------------
        Gappath = -1.0 * (Phipath.<0.0) .* Phipath
        Phipath += Gappath  # gap covered

        return Phipath::Array{Float64,1}, Gappath::Array{Float64,1}
    end
    # ==========================================================================


    # ==========================================================================
    # SPECIAL-DESIGNED FUNCTIONS USED IN TRANSITION SEARCHING
    # ==========================================================================


    # ==========================================================================
    """
        PolicySolve_r(A0::Float64,Phi0::Float64, MAX_AGE::Int64
            DATA::Array{Array{Float64,1},1},DATA_r::Array{Array{Float64,1},1}, PARS::Array{Float64,1} ;
            flag_print::Bool=false, GRID_LEN::Int64=50, TOL::Float64=1E-8, Arange::Array{Float64,1}=[0.0,15.0])

    Solve Policy functions, only retired phase included.
    MAX_AGE = real MAX_AGE - real Start Age + 1, MAX_AGE is in the range 1,...,MAX_AGE-RETIRE_AGE.

    **DATA** contains (in order):
    1. Rpath::Array{Float64,1} : interest rate path
    1. Fpath::Array{Float64,1} : mortality path
    1. Qpath::Array{Float64,1} : m2c ratio path
    1. Ppath::Array{Float64,1} : outpatient / inpatient fee path
    1. cpBpath::Array{Float64,1} : inpatient copayment rate path

    **DATA_r** contains(in order):
    1. Lambdapath::Array{Float64,1} : pension benefits
    1. doublePpath::Array{Float64,1} : transfer (amount) from firm contribution to those retired

    **PARS** contains(in order):
    1. Mu : consumption tax rate
    1. Sigma : wage tax rate (though useless here, but remained for a universal & easy-coding style)
    1. Alpha : leisure preference than consumption
    1. Gamma : inter-temporal substitution elasticity
    1. Varrho : consumption leisure substitution elasticity
    1. Delta : utility discount rate
    """
    function PolicySolve_r(A0::Float64,Phi0::Float64, MAX_AGE::Int64 ,
        DATA::Array{Array{Float64,1},1},DATA_r::Array{Array{Float64,1},1}, PARS::Array{Float64,1} ;
        flag_print::Bool=false, GRID_LEN::Int64=50, TOL::Float64=1E-8, Arange::Array{Float64,1}=[0.0,15.0] )
        # ----------------------------------------------------
        # Timing
        time_begin = time()
        # Unfold data
        Rpath::Array{Float64,1},Fpath,Qpath,Ppath,cpBpath= DATA
        Lambdapath::Array{Float64,1},doublePpath = DATA_r
        Mu::Float64,Sigma,Alpha,Gamma,Varrho,Delta = PARS
        # Validation
        @assert(2==length(Arange),"Arange can only have two elements in PolicySolve_r()")
        @assert(MAX_AGE==length(Rpath)==length(Fpath)==length(Qpath)==length(Ppath)==length(cpBpath),"Uncompatible MAX_AGE-length inputs in PolicySolve_r()")
        @assert(MAX_AGE==length(Lambdapath)==length(doublePpath),"Uncompatible Lambdapath in PolicySolve_r()")
        # Algorithm data spacing
        Aopt = zeros(GRID_LEN,MAX_AGE)  # CAUTION: Aopt here is the summation of personal asset and individual medical account!!!
        Copt = zeros(GRID_LEN,MAX_AGE);
        Vr = zeros(GRID_LEN,MAX_AGE)   # Value function matrix
        Grid_A = Array(linspace(Arange[1],Arange[2],GRID_LEN))  # Griding vector of personal asset
        Apath = zeros(MAX_AGE)  # CAUTION: Aopt here is the summation of personal asset and individual medical account!!!
        apath = zeros(MAX_AGE)  # CAUTION: apath here is ONLY the personal asset!!!
        Phipath = zeros(MAX_AGE)  # CAUTION: Phipath here is ONLY the individual medical account!!!
        Cpath = zeros(MAX_AGE); # Final results
        # ----------------------
        # ----------------------
        # ------- BELLMAN: s=S ---------------
        idx_Age = MAX_AGE
        # Value function for s=S
        Aopt[:,idx_Age] = Grid_A
        PARS_DEAD = [  Rpath[idx_Age], Fpath[idx_Age], Qpath[idx_Age], Ppath[idx_Age], cpBpath[idx_Age],
            Mu, Lambdapath[idx_Age], doublePpath[idx_Age], Alpha, Gamma, Varrho, Delta  ]
        Vr[:,idx_Age] = [Bellman_dead( Grid_A[x], PARS_DEAD ) for x in 1:GRID_LEN]
        # Get optimal consumption
        Copt[:,idx_Age] = [GetCL_r( Grid_A[x], 0.0, PARS_DEAD[1:8] ) for x in 1:GRID_LEN]

        # CHECK & RETURN -------------
        if MAX_AGE == 1
            # rip off individual medical account
            apath[MAX_AGE] = A0; Phipath[MAX_AGE] = Phi0
            Cpath[MAX_AGE] = GetCL_r(A0,0.0,PARS_DEAD[1:8])
            return apath::Array{Float64,1}, Phipath::Array{Float64,1}, Cpath::Array{Float64,1}
        end

        # ------- BELLMAN: Sr<s<S ------------ (MAX_AGE>=2)
        for idx_Age in MAX_AGE-1:-1:1
            for idx_Grid in 1:GRID_LEN
                # ammo reloaded
                PARS_R = [  Rpath[idx_Age], Fpath[idx_Age], Qpath[idx_Age], Ppath[idx_Age], cpBpath[idx_Age],
                    Mu, Lambdapath[idx_Age], doublePpath[idx_Age], Alpha, Gamma, Varrho, Delta  ]
                afunc_Bellman_r(Anext::Float64) = Bellman_r(Grid_A[idx_Grid], Anext::Float64 ,Grid_A,Vr[:,idx_Age+1],PARS_R)
                # compute temporary value function for griding points
                tmpV = [afunc_Bellman_r(Grid_A[x]) for x in 1:GRID_LEN  ]
                # find the peak
                ~,tmp_LocOptAs = findmax(tmpV)
                if tmp_LocOptAs==1
                    Aopt[idx_Grid,idx_Age] = Arange[1]
                elseif tmp_LocOptAs==GRID_LEN
                    Aopt[idx_Grid,idx_Age] = Arange[2]
                else
                    gold_a = Grid_A[tmp_LocOptAs-1] # pars for Golden Section Searching
                    gold_b = Grid_A[tmp_LocOptAs]
                    gold_c = Grid_A[tmp_LocOptAs+1]
                    # Golden Section Searching
                    Aopt[idx_Grid,idx_Age] = Golden(afunc_Bellman_r::Function, gold_a::Float64,gold_b::Float64,gold_c::Float64,TOL=TOL  )
                end
                # Get interpolated value function & consumption
                Vr[idx_Grid,idx_Age] = afunc_Bellman_r(Aopt[idx_Grid,idx_Age])
                Copt[idx_Grid,idx_Age] = GetCL_r( Grid_A[idx_Grid], Aopt[idx_Grid,idx_Age], PARS_R[1:8] )
            end
        end

        # ----------- EXTRACT OPTIMAL a+Phi path -------------------
        Apath[1] = A0 + Phi0
        for idx_Age in 1:MAX_AGE-1
            PARS_R = [  Rpath[idx_Age], Fpath[idx_Age], Qpath[idx_Age], Ppath[idx_Age], cpBpath[idx_Age],
                Mu, Lambdapath[idx_Age], doublePpath[idx_Age] ]
            Apath[idx_Age+1] = APolate(Apath[idx_Age]::Float64, Aopt[:,idx_Age]::Array{Float64,1}, Arange::Array{Float64,1}  )
            Cpath[idx_Age] = GetCL_r(Apath[idx_Age],Apath[idx_Age+1],PARS_R)
        end
        # Get the consumption in the very last year
        idx_Age = MAX_AGE
        PARS_R = [  Rpath[idx_Age], Fpath[idx_Age], Qpath[idx_Age], Ppath[idx_Age], cpBpath[idx_Age],
            Mu, Lambdapath[idx_Age], doublePpath[idx_Age] ]
        Cpath[idx_Age] = GetCL_r(Apath[idx_Age],0.0::Float64,PARS_R)

        # ----- DEPART INDIVIDUAL-MEDICAL ACCOUNT FROM THE COMPREHENSIVE ASSET ACCOUNT ---------
        Phipath,~ = Depart_Phipath_r(Phi0,Cpath, MAX_AGE, DATA,DATA_r, PARS )
        # get personal asset path
        apath = Apath .- Phipath

        # Timing ends
        if flag_print
            println("(Policy) Time Elasped: ",time()-time_begin," s")
        end

        return apath::Array{Float64,1},Phipath::Array{Float64,1},Cpath::Array{Float64,1}
    end
    # ==========================================================================
    """
        Depart_Phipath_r(Phi0::Float64,Cpath::Array{Float64,1}, MAX_AGE::Int64,
            DATA::Array{Array{Float64,1},1},DATA_r::Array{Array{Float64,1},1}, PARS::Array{Float64,1} )

    Given optimal consumption, leisure paths, computes individual medical account path.
    Most parameters are the same as PolicySolve()

    return a Phipath and a Gappath, cauz we have the constraint: individual medical account is non-negative
    """
    function Depart_Phipath_r(Phi0::Float64,Cpath::Array{Float64,1}, MAX_AGE::Int64,
        DATA::Array{Array{Float64,1},1},DATA_r::Array{Array{Float64,1},1}, PARS::Array{Float64,1} )
        # ----------------------
        Rpath::Array{Float64,1},Fpath,Qpath,Ppath,cpBpath= DATA
        Lambdapath::Array{Float64,1},doublePpath = DATA_r
        Mu::Float64,Sigma,Alpha,Gamma,Varrho,Delta = PARS
        # ------------------------
        Phipath = zeros(MAX_AGE); Gappath = zeros(MAX_AGE)
        # for working phase
        Phipath[1] = Phi0
        if MAX_AGE>1
            # cauz if MAX_AGE==1, no need to run the following loop
            for idx_Age in 1:MAX_AGE-1
                Phipath[idx_Age+1] = (1+Rpath[idx_Age])*Phipath[idx_Age] + doublePpath[idx_Age]
                Phipath[idx_Age+1] -= Qpath[idx_Age]*Ppath[idx_Age]/(1+Ppath[idx_Age])*Cpath[idx_Age]
                Phipath[idx_Age+1] /= 2-1/(1-Fpath[idx_Age])
            end
        else
            idx_Age = 1
            Phipath[idx_Age+1] = (1+Rpath[idx_Age])*Phipath[idx_Age] + doublePpath[idx_Age]
            Phipath[idx_Age+1] -= Qpath[idx_Age]*Ppath[idx_Age]/(1+Ppath[idx_Age])*Cpath[idx_Age]
            Phipath[idx_Age+1] /= 2-1/(1-Fpath[idx_Age])
        end
        # ----------- Compute Gaps -----------------
        Gappath = -1.0 * (Phipath.<0.0) .* Phipath
        Phipath += Gappath  # gap covered

        return Phipath::Array{Float64,1}, Gappath::Array{Float64,1}
    end
    # ==========================================================================










# ========================== MOD END ===========================================
end
# ========================== MOD END ===========================================
#
