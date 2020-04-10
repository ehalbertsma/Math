######## 0. FOURIER SERIES #####################
# any periodic function f can be approximated by a fourier series: 
# f(x) = lim{n->infinity} a0/2 + Σ_n a_j*cos(jx) + b_j*sin(jx)
# this program generates the fourier series of a function truncated to n terms, which is approximately equal to f
# assumptions:
#   - we assume the input function f is periodic over (-pi,pi)
# instructions:
#   - write your input function in section 2. INPUT FUNCTION
#   - enter your number of iterations as N
# output:
#   - an array 'functionarray' containing the input function and each fourier series truncated from 2 to N
#   - a plot showing the input function and each fourier series truncated from 2 to N for comparison

# import packages
using Plots
gr()
Plots.GRBackend()


########### 1. FUNCTIONS ########################
function integr(f,xrange=Array{Float64}(-pi:0.01:pi))
    # inputs:
    #   function f
    #   array xrange containing the input domain of the function
    # outputs:
    #   left endpoint integral of f over xrange (step precision equal to that of xrange)

    I = 0
    dx = xrange[2]-xrange[1]
    for xi in xrange
        I+=f(xi)*dx
    end
    return I
end

function InnerProd(f1, f2, xrange=Array{Float64}(-pi:0.01:pi))
    # inputs:
    #   function f1
    #   function f2
    #   array xrange containing the integral domain of the function
    # outputs:
    #   integral inner product of f1 with f2 over the xrange, rounded to 4 places

    # integrand
    f3(x) = f1(x)*f2(x)
    # return integral
    return round(integr(f3,xrange), digits=4)
end

function SymmetryCheck(f, axis=0.0, dist=pi)
    # inputs:
    #   function 'f'
    #   float axis of symmetry 'axis'
    #   positive float 'dist', specifies how far from the axis you want to check for symmetry
    # outputs:
    #   boolean even, true if function is even
    #   boolean odd, true if function is odd

    even, odd = false, false
    
    xi = axis+0.01 # starting point shifted slightly in case the axis of symmetry is an asymptote
    eps = 0.1 # tolerance (this is both delta-x and the tolerance in asymmetry)
    
    # check for evenness (from the definition of evenness) from the axis of symmetry until the specified radius
    # function can be considered even as long as f(xi) and f(-xi) are within a tolerance
    while (abs(f(xi)-f(-xi)) < eps) && (xi<dist)
        xi+=eps
    end
    # if and only if there are no asymmetries, the loop will run to completion. so if it ran all the way, we flag it as even
    if abs(xi-dist) <= 0.1
        even = true
    end

    # only check for oddness if it's not even
    if !even
        xi = axis+0.01 # reset xi
        # search for symmetry-breaking points
        # check for oddness (from the definition of oddness) from the axis of symmetry until the specified radius
        while (abs(f(xi)+f(-xi)) < eps) && (xi<dist)
            xi+=eps
        end
        # if the loop runs to completion, there are no symmetry-breaking points
        if abs(xi-dist) <= 0.1
            odd = true
        end
    end

    return even, odd
end

function cosm(m::Int64)
    # inputs:
    #   m - a 64 bit integer
    # outputs:
    #   a function of x with output cos(mx)
    return function(x) return cos(m*x) end
end
function sinm(m::Int64)
    # inputs:
    #   m - a 64 bit integer
    # outputs:
    #   a function of x with output sin(mx)
    return function(x) return sin(m*x) end
end

function FourierCoefficients(f, N::Int,xrange=Array{Float64}(-pi:0.01:pi))
    # inputs:
    #   - a function "f"
    #   - number of iterations (truncation limit of the fourier series) "N"
    #   - array "xrange", the range of the periodic function f
    # outputs:
    #   - a0: the constant coefficient of the fourier series
    #   - aCof: array, the even coefficients of the fourier series (for the cos terms)
    #   - bCof: array, the odd coefficients of the fourier series (for the sin terms)

    aCof = Array{Float64}(undef, N) 
    bCof = Array{Float64}(undef, N) 
    
    unity(x) = 1

    # a_0 coefficient
    a0 = (1/pi)*InnerProd(f,unity,xrange)
    
    # a_i coefficients, 1<=i<=N
    for i in 1:1:N
        aCof[i] = (1/pi)*InnerProd(f,cosm(i),xrange)
    end

    # b_i coefficients, 1<=i<=N
    for i in 1:N
        bCof[i] = (1/pi)*InnerProd(f,sinm(i),xrange)
    end

    return a0, aCof, bCof
end

function FourierCoefficientsSym(f, N,xrange=Array{Float64}(-pi:0.01:pi))
    # This function is the same as FourierCoefficients(), except it takes symmetries into account.
    # As a result it runs 40% faster.
    # inputs:
    #   - a function "f"
    #   - number of iterations (truncation limit of the fourier series) "N"
    #   - array "xrange", the range of the periodic function f
    # outputs:
    #   - a0: the constant coefficient of the fourier series
    #   - aCof: array, the even coefficients of the fourier series (for the cos terms)
    #   - bCof: array, the odd coefficients of the fourier series (for the sin terms)

    a0 = 0.0
    aCof = zeros(Float64, N)
    bCof = zeros(Float64, N)

    # first, check for symmetry
    # 1. if f is even, all b coefficients vanish. only need to find a_n
    # 2. if f is odd, all a coefficients vanish. only need to find b_n
    # 3. if f is neither even nor odd, we must calculate. need to find a_n, b_n

    even,odd = SymmetryCheck(f)

    if !odd
        unity(x) = 1
        # a_0 coefficient is obtained by the inner product of f and 1
        a0 = (1/(2*pi))*InnerProd(f,unity,xrange)
        
        # a_i coefficients, 1<=i<=N
        for i in 1:1:N
            aCof[i] = (1/pi)*InnerProd(f,cosm(i),xrange)
        end
    end
    if !even
        # b_i coefficients, 1<=i<=N
        for i in 1:N
            bCof[i] = (1/pi)*InnerProd(f,sinm(i),xrange)
        end
    end

    return a0, aCof, bCof
end

function FourierSeries(f, N=25,sym=true)
    # This function calculates the fourier series to N terms for a given function.
    # inputs:
    #   - a function "f"
    #   - number of iterations (truncation limit of the fourier series) "N"
    #   - Boolean "sym" to specify whether we want to check for symmetries in f for faster computation
    # outputs:
    #   - function "fourier", which is the actual fourier series

    if sym
        @time a0, aCof, bCof = FourierCoefficientsSym(f, N)
    else
        @time a0, aCof, bCof = FourierCoefficients(f, N)
    end
    # println("Cos coeffs ",aCof)
    # println("Sin coeffs ",bCof)
    # println("a0 ",a0)

    fourier = function(x)

            cosses = [cosm(i)(x) for i in 1:1:N]
            sins = [sinm(i)(x) for i in 1:1:N]
            
            S(x) = a0 + sum(aCof.*cosses) + sum(bCof.*sins)
            return S(x)
        end
    
    return fourier

end

######## 2. INPUT FUNCTION #####################
# the function that we want to find a Fourier series of:
function f(x) return (x-1)^2 end
# truncate the fourier series at N (minimum = 2)
N = 25

############# 3. TEST ##########################
# - creates an array containing all 2 to N truncated fourier series

# initialize an array of functions to store your fourier series
functionarray = Array{Function}(undef,(N+1))

# place the input function as the first element of functionarray
functionarray[1] = f
# place the i-th truncated fourier series in the i-th element
for i in 2:(N+1)
    functionarray[i] = FourierSeries(f,i)
end

############# 4. PLOT ##########################
# - plots all 2 to N fourier series alongside the original function for comparison

# set the lower axis to contain fractions of pi instead of decimal values
xticcs = ([-π:π/2:π;], ["-\\pi","-\\pi/2","0","\\pi/2","\\pi"])
# plot
p1 = plot(functionarray, -pi, pi, leg=false, xticks=xticcs)
title!("Fourier Series approximation to $(N) terms")

# causes the plot to appear on-screen
gui()

# ensures the program doesn't stop early
print("Press enter to quit ")
readline()
