#Packages
using Plots
using Interpolations
using BenchmarkTools
using ProgressMeter
using Parameters
using QuantEcon

#1. Parameter Section
@everywhere MyParam1 = @with_kw (
             T = 100, # Life-Horizon
             β = 0.98, # Discount Factor
             minc = 0.01, # Minimum Consumption Share
             r = 1.02, # Risk-free Rate
             σ = 2.0, # CRRA Coefficient
             σ_v = 1.06, # Std of Labor Income Shock
             μ_v = 0.05 # Mean of Labor Income Shock
             )

MyParam1c = MyParam1()
#typeof(MyParam1c)

#-> Gaussian Quadrature for normal distribution (discretize the shock and proceed more)











# Define the parameters as a dictionary
# TODO: how does speed compare to using a constant?
# TODO: create a Dict of various objects (params, objs, etc)
# TODO: https://docs.julialang.org/en/stable/manual/performance-tips/#tools-1
params                     = Dict{String, Float64}()
params["tol"]              = 1e-10               # max allowed error
params["minCons"]          = 1e-5                # min allowed consumption
params["r"]                = 1.0/0.95 - 1.0      # Interest rate
params["beta"]             = 0.95                # 1/(1+r) # Discount factor
params["gamma"]            = 1.5                 # Coefficient of relative risk aversion
params["gamma_mod"]        = 1.0-params["gamma"] # For speed, just do this once
params["startA"]           = 0.0                 # How much asset do people start life with
params["mu"]               = 0.0                 # mean of initial log income
params["sigma"]            = 0.25                # variance of innovations to log income
params["rho"]              = 0.75                # persistency of log income

# Constants
const interpMethod         = "linear"            # for now, I only allow linear option
const T                    = 60                  # Number of time period
const Tretire              = 45                  # Age at which retirement happens
const borrowingAllowed     = 0                   # allow borrowing
const isUncertainty        = 1                   # uncertain income (currently: only works if isUncertainty == 1)
const numPointsY           = 9                   # number of points in the income grid
const numPointsA           = 50                  # number of points in the discretised asset grid
const gridMethod           = "5logsteps"         # method to construct grid. One of equalsteps or 5logsteps
const normBnd              = 3                   # truncate the normal distrib: ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma
const numSims              = 10                  # How many individuals to simulate
const useEulerEquation     = false               # Solve the model using the euler equation?
const saveValue_inEE       = false               # When using euler equation to solve the model, do we want to compute EV? (Note: adds time due to interpolation)
const linearise            = false               # Whether to linearise the slope of EdU when using EE




function getGrid(minongrid, maxongrid, GridPoints, method)
    span = maxongrid - minongrid
    if method == "equalsteps"
        grid= linspace(minongrid, span, GridPoints)
    elseif method == "5logsteps"
        loggrid = linspace(log(1+log(1+log(1+log(1+log(1))))), log(1+log(1+log(1+log(1+log(1+span))))), GridPoints)
        grid = exp.(exp.(exp.(exp.(exp.(loggrid)-1)-1)-1)-1)-1
    end
end


function getIncomeGrid(params)

    # A function that returns:
    # 1. an income grid
    # 2. A Markovian transition matrix (Q) over income realisations
    # 3. A vector of minimum incomes in each year
    # 4. A vector of maximum incomes in each year

    #----------------------------------------#
    # Scenario where there is no uncertainty
    #----------------------------------------#
    if isUncertainty == 0
        y = exp( params["mu"] )                 # income is set equal to the exp of the log mean
        minInc = y
        maxInc = y
        Q = [1]                    # The transition matrix Q is simply a constant 1
                                   # with prob 1 each period income is 1

       #----------------------------------------#
       #Now get a matrix, T * numPointsY that holds the grid for each income in
       #each year. Do likewise with minimum and maximum income
       #----------------------------------------#
       Ygrid = repmat([y'], T, 1)
       minInc = repmat([minInc'], T, 1)
       maxInc = repmat([maxInc'], T, 1)

    #----------------------------------------#
    # Scenario where there is uncertainty - income draws are log normally distributed
    #----------------------------------------#

elseif isUncertainty == 1

           # First get the standard deviation of income (from sigma and rho)
           sig_inc = params["sigma"] / ((1- params["rho"]^2)^0.5)

           # Split the entire normal distribution into numPointsY sections that
           # are equiprobable. The output lNormDev gives the (numPointsY + 1)
           # points that bound the sections, the output ly gives the
           # (numPointsY) expected value in each section
           lNormDev, ly = getNormDev(params["mu"], sig_inc, normBnd, numPointsY );

           #---------------------#
           #Get transition matrix Q(i, j). The prob of income j in t+1
           #conditional on income i in t
           #---------------------#
           Q = zeros(numPointsY, numPointsY)              #initialise the transition matrix
           for i = 1:1:numPointsY
               for j = 1:1:numPointsY
                   hiDraw = lNormDev[j+1] - (1-params["rho"])*params["mu"] - params["rho"] * ly[i]; #highest innovation that will give us income j tomorrow
                   loDraw = lNormDev[j]   - (1-params["rho"])*params["mu"] - params["rho"] * ly[i]; #lowest  innovation that will give us income j tomorrow
                   Q[i,j] = stdnormcdf_manual(hiDraw/params["sigma"]) - stdnormcdf_manual(loDraw/params["sigma"]);
               end #j

               #Each of the rows of Q should add up to 1. But
               #due to the truncation of the normal distribution they will
               #not. So we divide through by the sum of elements in the row
               Q[i, :] = Q[i, :] ./ sum(Q[i, :])
           end #i

           y = exp.(ly);                      # Get y from log y
           minInc = exp(-normBnd * sig_inc); #Get the minimum income in each year
           maxInc = exp(normBnd * sig_inc);  #Get the maximum income in each year

           if (y[1] < 1e-4) || (y[ numPointsY ] > 1e5)
               warning("Combination of sigma and rho give a very high income variance. Numerical instability possible")
           end

           #----------------------------------------#
           #Now get a matrix, T * numPointsY that holds the grid for each income in
           #each year. Do likewise with minimum and maximum income
           #----------------------------------------#
           Ygrid = repmat(y', T, 1)
           minInc = repmat([minInc'], T, 1)
           maxInc = repmat([maxInc'], T, 1)

    end  # if isUncertainty == 0



    #----------------------------------------#
    #Replace these arrays with zeros for all years after retirement
    #----------------------------------------#

    if Tretire == 0         # no work (retired at birth)
       Ygrid[:, :] = 0
       minInc[:, :] = 0
       maxInc[:, :] = 0
    elseif (Tretire > 0) && (Tretire <=T)  #retire at some age
        Ygrid[Tretire:T, :] = 0
        minInc[Tretire:T, :] = 0
        maxInc[Tretire:T, :] = 0
    end

    return Ygrid, Q, minInc, maxInc

end
