using Plots, LaTeXStrings, LinearAlgebra, Random
gr()
default(
    titlefont = (16, "times"), 
    legendfontsize = 12, 
    guidefont = (14), # Changes x and y axis label fonts
    linewidth = 2)

x_actual = collect(1:0.01:2.5);
y_actual = cos.(2π * x_actual) .* exp.(-x_actual);
plot(x_actual, y_actual, title=L"$f(x) = e^{-x}cos(2\pi x)$", legend=false)
xlabel!(L"x")
ylabel!(L"f(x)")

N = 50

# Allows us to get the same random numbers every time we run this cell
Random.seed!(12345678)

#= 
idx - random indices 
We use unique() because after flooring the numbers 
(rounding numbers to the nearest integer less than or equal to the number), 
we may have repeated indices
=#
idx = Int.(unique(floor.(rand(N) * length(x_actual))) .+ 1) 
N   = length(idx);

# Training input
x_measured = x_actual[idx]; 

# Noise from a scaled normal distribution
y_measured = y_actual[idx] + 0.02 * randn(N)

println("Each of the yellow dots is a measurement that we will 
    use for building an approximation of the function f(x) ")

plot(x_actual, y_actual, title="Noisy measurements from a known function",
    label=L"$f(x) = e^{-x}cos(2\pi x)$", xlabel=L"$x$", ylabel=L"$f(x)$")
scatter!(
    x_measured, 
    y_measured, 
    c=:orange, # set the color
    label="Noisy Measurements", 
    legend=:best) # automatically use best location to place legend in graph
N = length(x_measured)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Number of measurements
N = 44

# Example measured x values (replace this with actual data if available)
x_measured = rand(N)  # Replace with actual data for `x_measured`

# Initialize Phi as an N x 4 matrix with zeros
Phi = zeros(N, 4)

# Populate Phi row-by-row
for i in 1:N
    Phi[i, :] = [1, x_measured[i], x_measured[i]^2, x_measured[i]^3]
end
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#=
forwardsub(L, b)

It solves for x in an equation Lx = b, where L is lower triangular.
=#

function forwardsub(L, b)
    n = length(b)
    x = Vector{Float64}(undef, n); 
    x[1] = b[1]/L[1,1] 
    for i = 2:n 
        x[i] = (b[i]- (L[i,1:i-1])'*x[1:i-1] )/L[i,i] 
    end
    return x
end

#=
backwardsub(U, b)

It solves for x in an equation Ux = b, where U is upper triangular.
=#

function backwardsub(U, b)

    if minimum(abs.(diag(U))) < 1e-6
        println(" U is nearly singular. ")
        return x = NaN
    end
    
    n = length(b)
    x = Vector{Float64}(undef, n)
    x[n] = b[n] / U[n,n]
    for i = n-1:-1:1
        x[i] = (b[i] - (U[i,(i+1):n])' * x[(i+1):n]) / U[i,i]
    end
    
    return x    
end

# Phi::Array{Float64,2} simply checks that the input is a real matrix
# Y::Array{Float64,1} simply checks that the input is a real vector 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function least_squares_lu(Phi::Array{Float64,2}, Y::Array{Float64,1})
    # Step 1: Construct A and b
    A = Phi' * Phi
    b = Phi' * Y

    # Step 2: Perform LU decomposition with permutation
    LU_decomp = lu(A)
    L = LU_decomp.L
    U = LU_decomp.U
    P = LU_decomp.P  # Permutation matrix

    # Step 3: Solve the system using forward and backward substitution
    # Solve Lz = Pb
    z = forwardsub(L, P * b)
    # Solve U * alphaStar = z
    alphaStar = backwardsub(U, z)

    return alphaStar
end

^^^^^^^^^^^^^^^^^^^ 

# Calling the function and returning a_star
a_star = least_squares_lu(Phi,y_measured)

# Let's see the cofficients for our monomial fit model 
@show a1 = a_star[1];
@show a2 = a_star[2];
@show a3 = a_star[3];
@show a4 = a_star[4];
println("This is a_star!")
a_monomial = copy(a_star)

^^^^^^^^^^^^^^^^^^^^^^
#=
x is a scalar value (real number) or a vector of real numbers
a is a vector of parameters for the polynomial (a_star)
yHat is the polynomial approximation to the data at the point(s) x
=#
function yHatMonomials(x::Union{Real, AbstractVector{<:Real}}, a::Vector{Float64}=a_star)
    # Ensure 'a' has at least one coefficient
    if length(a) < 1
        error("Coefficient vector 'a' must have at least one value.")
    end
    
    # Evaluate the polynomial based on the type of x
    if isa(x, Real)
        # Scalar case: Compute yHat for scalar input x
        yHat = sum(a[i] * x^(i-1) for i in 1:length(a))
        return yHat
    elseif isa(x, AbstractVector{<:Real})
        # Vector case: Compute yHat for each element in x
        return [sum(a[i] * xi^(i-1) for i in 1:length(a)) for xi in x]
    else
        error("Input 'x' must be a scalar or a vector of real numbers.")
    end
end
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Test the implementation with the provided cases
is_it_correct_check1 = ((yHatMonomials(-1) - 47.694823671937684) < 1e-5) ? "Yes" : "No" 
is_it_correct_check2 = ((yHatMonomials(2) - 0.1058911899932653) < 1e-5) ? "Yes" : "No" 
is_it_correct_check3 = (norm(yHatMonomials.([1; 2]) - [0.579297437402257; 0.1058911899932653]) < 1e-5) ? "Yes" : "No" 

@show is_it_correct_check1;
@show is_it_correct_check2;
@show is_it_correct_check3;

check = 0

if is_it_correct_check1 == "Yes" && is_it_correct_check2 == "Yes" && is_it_correct_check3 == "Yes"
  check = 1
end

@assert check == 1

//////returns yes yes no for this **error**