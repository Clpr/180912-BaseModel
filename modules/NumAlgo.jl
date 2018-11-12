__precompile__()
"""
    NumAlgo

A module containing all numerical algorithm methods.
1. Golden : Golden Section Search
1. Apolate: interpolating A in a given A series (discrete)
1. VecExpand: Expand vector to a longer one (used in transition path compute)

"""
module NumAlgo
    export Golden, Bisection, APolate, VecExpand, VecExpand_ContentSpecific, VFunc, diagwrite
# ========================== MOD BEGIN =========================================
    """
        Golden(f::Function, LLB::Float64, LB::Float64, RB::Float64 ;  TOL::Float64=1E-08)

    Use Golden Method to get more precise solution (a) for a 1-D maximization problem

    Par:
        1. f [annonymous func]: target to search (a maximization problem), it MUST HAVE ONLY ONE PARAMETER f(x)
        1. LLB [num]: the (initial) very left bound of searching (initial)
        1. LB [num]: the (initial) left bound of searching
        1. RB [num]: the (initial) right bound of searching
        1. TOL [num]: tolerance, 1E-4 or 1E-5 is enough

    Ret:
        1. Xmin [num]: the solution

    Depend:
        1. func
    """
    function Golden(f::Function, LLB::Float64, LB::Float64, RB::Float64 ;  TOL::Float64=1E-08)
        # -------- DATA PROCESS ----------
        # NOTE: in Julia, the golden numebr is an integrated const
        r1 = golden - 1; r2 = 1 - r1
        x0 = LLB; x3 = RB
        if abs(RB-LB)<=abs(LB-LLB)
            x1 = LB; x2 = LB+r2*(RB-LB)
        else
            x2 = LB; x1 = LB-r2*(LB-LLB)
        end
        # initialization of function value (and turns a maximization to a minimization)
        f1 = - f(x1); f2 = - f(x2)
        # Searching
        for counter in 1:300
            if f2<f1
                x0=x1 # Update the very lower bound
                x1=x2; x2=r1*x1+r2*x3 # new left golden position
                f1=f2; f2= - f(x2)
            else
                x3=x2; x2=x1; x1=r1*x2+r2*x0; f2=f1; f1= - f(x1)
            end
            if abs(x3-x0)>TOL*(abs(x1)+abs(x2))
                break
            end
        end
        # Post-convergence
        if f1<=f2
            Xmin=x1
        else
            Xmin=x2
        end
        return Xmin::Float64
    end
    # ==========================================================================
    """
        Bisection(f::Function, X1::Float64, X2::Float64 ; MAXLOOP::Int64 = 5000, TOL::Float64 = 1E-12  )

    search ZERO-point for single input function f;
    if f returns multiple results, selects the 1st member;
    requires f(X1)*f(X2)<0 && X1 < X2;
    MAXLOOP decides the max loop to search, and TOL decides the tolerance of |X2 - X1|

    returns a Float64 "Xzero", the average of X2 & X1 in the last loop;
    if TOL not met, throws an error;
    """
    function Bisection(f::Function, X1::Float64, X2::Float64 ; MAXLOOP::Int64 = 5000, TOL::Float64 = 1E-12  )
        # 1. temp vars
            Xhalf = (X1+X2)/2
        # 2. initial evaluation
            F1 = f(X1)[1]; F2 = f(X2)[1]; Fhalf = f(Xhalf)[1]
        # 3. validation
            ( (X1<X2)&&(F1*F2<0) ) || throw(ErrorException("bisection search requires (X1<X2)&&(F1*F2<0)"))
        # 4. searching
            for idx in 1:MAXLOOP
                # 4.1 if tolarence met, returns the mid point
                    abs(X1-X2)<TOL && return (X1+X2)/2
                # 4.2 evaluation of this loop
                    Xhalf = (X1+X2)/2
                    F1 = f(X1)[1]; F2 = f(X2)[1]; Fhalf = f(Xhalf)[1]
                # 4.3 update
                    if F1*Fhalf < 0.0
                        X1 = X1; X2 = Xhalf
                        continue
                    elseif Fhalf*F2 < 0.0
                        X1 = Xhalf; X2 = X2
                        continue
                    else
                        if F1 == 0.0
                            return X1::Float64
                        elseif F2 == 0.0
                            return X2::Float64
                        else
                            return Xhalf::Float64
                        end
                    end
            end
        # 5. if MAXLOOP reached, throws an error
            throw(ErrorException(string("cannot find a zero point in ",MAXLOOP," rounds, please consider a larger MAXLOOP or check your object function"  )  ))
        return nothing
    end
    # ==========================================================================
    """
            APolate(Athis::Float64, ASeries::Array{Float64,1}, Arange::Array{Float64,1}  )

    APolate - interpolating policy functions for next period's asset

    Pars:
        1. Athis::Float64 - asset (non grid) value
        1. idx_Age::Int64 - age index
        1. ASeries::Array{Float64,1} - value function column
        1. Arange::Array{Float64,1} - range of asset, a 2-ele 1-D vector

    Ret:
        1. Apolated::Float64 - interpolated asset
    """
    function APolate(Athis::Float64, ASeries::Array{Float64,1}, Arange::Array{Float64,1}  )
        # ----------------------
        Amax = Arange[2]
        Na = length(ASeries)::Int64
        A0 = Athis/Amax * (Na-1) + 1
        n2 = floor(Int64,A0)
        if n2>(Na-1)
            n2 = Na-1
        end
        n1 = A0 - n2
        # Interpolating
        Apolated = (1-n1)*ASeries[n2] + n1*ASeries[n2+1]
        # Return
        return Apolated::Float64
    end
    # ==========================================================================
    """
        VecExpand(Vec::Array{Float64,1},Len::Int64; Direct="forward")

    Expand input vector to a longer one, where Len says the expanded length,
    and Direct= says which direction to expand ("forward" or "backward").
    Uses the nearest element to fill expanded positions.

    e.g. [1,2,3] --"forward"--> [1,2,3,3,3,3]

    [1,2,3] --"backward"--> [1,1,1,1,2,3]
    """
    function VecExpand(Vec::Array{Float64,1},Len::Int64; Direct="forward")
        Len0 = length(Vec)
        if Len0>=Len
            if Direct=="forward"
                return Vec[Len0-Len+1:end]
            elseif Direct=="backward"
                return Vec[1:Len]
            end
        elseif Direct=="forward"
            Vec1 = zeros(Len); Vec1[1:Len0]=Vec; Vec1[Len0+1:end]=Vec[end]
            return Vec1
        else # no diff backward & other cases
            Vec1 = zeros(Len); Vec1[Len-Len0+1:end]=Vec; Vec1[1:Len-Len0]=Vec[1]
            return Vec1
        end
    end
    # ==========================================================================
    """
        VecExpand_ContentSpecific(Vec::Array{Float64,1},ExtraVec::Array{Float64,1},Len::Int64; Direct="forward")

    Expand input vector to a longer one, where Len says the expanded length,
    and Direct= says which direction to expand ("forward" or "backward").
    Uses the corresponding elements in ExtraVec to fill blanks.

    ExtraVec should have the length Len.
    """
    function VecExpand_ContentSpecific(Vec::Array{Float64,1},ExtraVec::Array{Float64,1},Len::Int64; Direct="forward")
        Len0 = length(Vec)
        @assert(length(ExtraVec)==Len, "Uncompatible ExtraVec !"  )
        if Len0>=Len
            if Direct=="forward"
                return Vec[Len0-Len+1:end]
            elseif Direct=="backward"
                return Vec[1:Len]
            end
        elseif Direct=="forward"
            Vec1 = zeros(Len); Vec1[1:Len0]=Vec; Vec1[Len0+1:end]=ExtraVec[Len0+1:end]
            return Vec1
        else # no diff backward & other cases
            Vec1 = zeros(Len); Vec1[Len-Len0+1:end]=Vec; Vec1[1:Len-Len0]=ExtraVec[1:Len-Len0]
            return Vec1
        end
    end
    # ==========================================================================
    """
        VFunc(A::Float64, ASeries::Array{Float64,1}, VSeries::Array{Float64,1} )

    VFunc: Searching/Interpolating for (AFTER RETIRE) Value function v(s+1) given the concrete empirical series of f(a)=max v(a)

    Par:
        1. A [vec]: the asset a(s+1) to look up in the empirical function (s+1) [ASeries] -> [VSeries] (i.e. the concrete form of f(a)=max v(a)). The A neednt to belong to ASeries, but in the range [min(ASeries),max(ASeries)]
        2. ASeries [row-vec]: the concrete state space of A, increasingly sorted, and paired with elements in VSeries
        3. VSeries [row-vec,len=len(ASeries)]: the concrete function value of f(a)=max v(a) given an A in ASeries

    Ret:
        1. rv [vec]: the interpolated function value of f(A)=max v(A)

    Depend:
        None

    NOTE:
        0. The function now supports vectorized A
        1. The function uses simple linear interpolation
        2. A which is out of range [min(ASeries),max(ASeries)] will be bounded to min(ASeries) or max(ASeries)
        3. The function is usually used in the computation of v(s+1) in the Bellman Equation v(s) = max[ u(s) + beta*v(s+1) ]
        4. It is improved based on the gauss language programmes by Burkhard Heer in 2008. It now does not depend on global variables and the age index.
        5. Pls refer to: https://www.wiwi.uni-augsburg.de/vwl/maussner/dge_buch/dge_book_2ed/downloads_2nd/
    """
    function VFunc(A::Float64, ASeries::Array{Float64,1}, VSeries::Array{Float64,1} )
        # --------------------
        Len = length(ASeries)
        Amax,~ = findmax(ASeries)
        Amin,~ = findmin(ASeries)
        Loc = (A-Amin)/(Amax-Amin).*(Len-1)+1 # Locate A0 (continuously)
        Lb = floor(Int64,Loc) # The nearest lower bound to interpolate
        Dis = Loc - Lb # The distance between A's index and the Lower Bound (index)
        rv = 0
        # Those touching the lower bound
        if A<=Amin
            rv = VSeries[1] - (1-Loc).*(VSeries[2]-VSeries[1])
        elseif A>=Amax
        # The upper bound
            rv  = VSeries[end]
        else
        # The left, doing linear interpolation
            rv = (1-Dis).*VSeries[Lb] + Dis.*VSeries[Lb+1]
        end

        return rv::Float64
    end
    # ==========================================================================
    """
    diagwrite(MAT::Array{Float64,2}, VEC::Array{Float64,1} ; offset::Int64 = 0, RETURN::Bool = false)

    A MatLab also Python style function, both supported.
    writes a vector "VEC" into a matrix "MAT" diagonally, according to offset= (value starts from 0).
    two mode supported (default by reference, the python style):
        1. by reference (no return): just set RETURN=true. works in a Python style. the input "MAT" itself will be modified.
        1. by copy (return filled matrix): set RETURN=false. works in a Matlab style. the outside "MAT" will not be modified. but much slower (about 100x time-cost than by-reference)
    """
    function diagwrite(MAT::Array{Float64,2}, VEC::Array{Float64,1} ; offset::Int64 = 0, RETURN::Bool = false)
        # validation
        LEN = length(VEC)
        ROW,COL = size(MAT)
        ( (1-offset<=ROW)&&(offset<COL) ) || throw(ArgumentError("invalid offset= value, it must within the size of matrix to write"))
		# adjust LEN, to avoid too-long VEC (if VEC is too-long, only fill the MAT with beginning elements of VEC, if too-short, fill all elements of VEC)
        LEN = min( LEN, min(ROW,COL) )
        # if a return required, use copy but not reference
        RETURN && (MAT=copy(MAT))  
        # filling
        if offset == 0
            for idx in 1:LEN; MAT[idx,idx] = VEC[idx]; end
        elseif offset > 0
            for idx in 1:LEN; MAT[idx,offset+idx] = VEC[idx]; end
        elseif offset < 0
            for idx in 1:LEN; MAT[idx-offset,idx] = VEC[idx]; end
        end
        # return nothing if RETURN=false
        RETURN || return nothing
        # or, return the copy of MAT
        return MAT
    end
    # ==========================================================================





















# ========================== MOD END ===========================================
end
# ========================== MOD END ===========================================
#
