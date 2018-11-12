__precompile__()
"""
    OLGDatTypes



"""
module OLGDatTypes
    using PyPlot
    export YearSnapshot,  # a data package struct for Steady state result
           TransitionPath  # a data package struct for transition path solved
    # -------------------------------------------------------------------------------------
    """
        YearSeries(Series,StartYear::Int64)
    
    a structure to store year data or age data;
    """
    struct YearSeries
        # ----- attributes --------
        Years::Vector{Int}  # series of years
        Series::Vector  # data series
        # ---- constructor --------
        YearSeries(Series,StartYear::Int64) = begin
            Length = length(Series)
            return new( collect(Int,StartYear:(StartYear+Length-1)) ,Series  )
        end
    end
    # corresponding functions
    
    # -------------------------------------------------------------------------------------
    """
        YearAgeMatrix(Mat::Matrix,StartYear::Int,StartAge::Int)

    a matrix datatype for year * age matrix
    """
    struct YearAgeMatrix
        # ------ attributes --------
        Years::Vector{Int}  # year as rows
        Ages::Vector{Int}  # ages as columns
        Mat::Matrix  # data
        # ------ constructor -------
        YearAgeMatrix(Mat::Matrix,StartYear::Int,StartAge::Int) = begin
            T,S = size(Mat)
            return new(  collect(Int,StartYear:StartYear+T-1), collect(Int,StartAge:StartAge+S-1), Mat  )
        end
    end
    # corresponding functions

    # -------------------------------------------------------------------------------------
    """



    """
    mutable struct YearSnapshot
        # ------ scalars -------------
        Y::Real  # GDP
        K::Real  # capital
        L::Real  # labour
        r::Real  # net investment interest rate
        wmean::Real  # social average wage level
        Lambda::Real  # pension benefit level
        


        
    end
































end  # module ends
# ===============================================================================