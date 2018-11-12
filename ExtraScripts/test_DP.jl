push!(LOAD_PATH,"./modules");
using NumAlgo
using PolicyFunctions

# reload("NumAlgo")
# reload("PolicyFunctions")

MAX_AGE = 80
RETIRE_AGE = 40

Rpath = 0.03*ones(MAX_AGE)
Fpath = 0.01*ones(MAX_AGE)
Qpath = 0.15*ones(MAX_AGE)
Ppath = 1.1*ones(MAX_AGE)
cpBpath = 0.3*ones(MAX_AGE)
# -----------------------------
phiCoefpath = 0.02*ones(RETIRE_AGE)
Zetapath = 0.06*ones(RETIRE_AGE)
Etapath = 0.21*ones(RETIRE_AGE)
Thetapath = 0.05*ones(RETIRE_AGE)
Wpath = 4.609*ones(RETIRE_AGE)
doubleApath = 0.3*ones(RETIRE_AGE)
zpath = 0.1*ones(RETIRE_AGE)
# ------------------------------
Lambdapath = 1.5*ones(MAX_AGE-RETIRE_AGE)
doublePpath = 0.01*ones(MAX_AGE-RETIRE_AGE)
# ------------------------------
Mu = 0.1; Sigma = 0.24; Alpha = 1.5; Gamma = 0.5; Varrho = 0.8; Delta = 0.015
# -------------------------------


DATA = [Rpath::Array{Float64,1},Fpath,Qpath,Ppath,cpBpath]
DATA_w = [phiCoefpath::Array{Float64,1},Zetapath,Etapath,Thetapath,zpath,Wpath,doubleApath]
DATA_r = [Lambdapath::Array{Float64,1},doublePpath]
PARS = [Mu::Float64,Sigma,Alpha,Gamma,Varrho,Delta]


# reload("NumAlgo")
# reload("PolicyFunctions")
apath = zeros(MAX_AGE); Phipath = zeros(MAX_AGE); Cpath = zeros(RETIRE_AGE); Lpath = zeros(RETIRE_AGE);


# for t in 1:5
#     @time for idx in 1:10
#             apath,Phipath,Cpath,Lpath = PolicySolve(0.0, 0.0, MAX_AGE, RETIRE_AGE,
#         DATA, DATA_w, DATA_r, PARS, flag_print=false, GRID_LEN=91, TOL=1E-8, Arange=[0.0,25.0]);
#     end
# end
# println("-"^20)
# for t in 1:5
#     @time for idx in 1:10
#             apath,Phipath,Cpath,Lpath = PolicySolve2(0.0, 0.0, MAX_AGE, RETIRE_AGE,
#         DATA, DATA_w, DATA_r, PARS, flag_print=false, GRID_LEN=91, TOL=1E-8, Arange=[0.0,25.0]);
#     end
# end

# -------------- profiling ------------------------
@time apath,Phipath,Cpath,Lpath = PolicySolve(0.0, 0.0, MAX_AGE, RETIRE_AGE, DATA, DATA_w, DATA_r, PARS, flag_print=false,
    GRID_LEN=5, TOL=1E-8, Arange=[0.0,25.0]);
# using ProfileView
# @profile PolicySolve2(0.0, 0.0, MAX_AGE, RETIRE_AGE, DATA, DATA_w, DATA_r, PARS, flag_print=false, GRID_LEN=91, TOL=1E-8, Arange=[0.0,25.0]);
# Profile.print()
# ProfileView.view();
















#
