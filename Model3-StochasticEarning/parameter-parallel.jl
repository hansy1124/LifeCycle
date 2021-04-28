#Parameter regarding economic problem
@everywhere MyParam1 = @with_kw (
             T = 15, # Life-Horizon
             reT = 10, # Retirement Age
             β = 0.98, # Discount Factor
             minc = 0.01, # Minimum Consumption Share
             r=1.03,
             σ = 2.0, # CRRA Coefficient
             σ_v = 0.001, # Std of Labor Income Shock
             μ_v = 0., # Mean of Labor Income Shock
             nv = 5# Number of grids for Income Shock
             )

#Parameter regarding numerical solution
@everywhere MyParam2 = @with_kw (
            nX = 500, # # of grid points for X
            minx = 0.1, # minimum x
            maxx = 100000., # maximum x
            nY = 100, # # of grid points for Y
            nA = 500 # # of grid points for A
            )

#Gaussian Hermite Quadrature for normal distribution
#http://fmwww.bc.edu/ec-p/software/Miranda/chapt6.pdf - 6.2
#https://www.wouterdenhaan.com/numerical/integrationslides.pdf - 24th slide

@everywhere function normal_gausshermite(μ, σ, n)
    x, w = gausshermite(n)
    x = √2σ*x .+ μ
    w = (1/√π)*w
    return x, w
end
# x - points / w - weights

#Life-Cycle Income Growth Rate(gbar)/Profile(lpf) Generator
@everywhere function generate_gbar(T, reT)
    lpf = [-(x-45.)^2 + 5000. for x in 1:T]
    lpf[reT:end] .= lpf[reT]*0.6
    append!(lpf,0.01)
    loglpf = log.(lpf)
    gbar = loglpf[2:end-1]-loglpf[1:end-2]
    append!(gbar,0.00)
    return gbar, lpf
end

#Life-Cycle Ygrid Generator
@everywhere function Ygrid_life(init_range, lpf, v, MyParam1c)
    #Initial Y grid
    minY=lpf[1]-init_range
    maxY=lpf[1]+init_range
    VmaxY=zeros(MyParam1c.T)
    VminY=zeros(MyParam1c.T)
    VmaxY[1] = maxY
    VminY[1] = minY
    #Before Retire
    for i in 2:1:MyParam1c.reT
        VmaxY[i] = VmaxY[i-1]*exp(gbar[i-1]+v[5])
        VminY[i] = VminY[i-1]*exp(gbar[i-1]+v[1])
    end
    #After Retire
    for i in MyParam1c.reT+1:1:MyParam1c.T
        VmaxY[i] = VmaxY[i-1]
        VminY[i] = VminY[i-1]
    end
    #After the death
    append!(VmaxY,0.00)
    append!(VminY,0.00)
    return VmaxY, VminY
end
