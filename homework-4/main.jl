using LinearAlgebra
using Plots

# Definition of the SparseMatrixCSR{T} type
struct SparseMatrixCSR{T} <: AbstractMatrix{T}
    rows::Vector{Int64}
    cols::Vector{Int64}
    values::Vector{T}

    function SparseMatrixCSR{T}(rows::Vector{Int64}, cols::Vector{Int64}, values::Vector{T}) where T
        new{T}(rows, cols, values)
    end
end

# Implementing size and getindex for SparseMatrixCSR
import Base: size, getindex

function size(A::SparseMatrixCSR)
    return (maximum(A.rows), maximum(A.cols))
end

function getindex(A::SparseMatrixCSR{T}, i::Int, j::Int) where T
    for k in 1:length(A.rows)
        if A.rows[k] == i && A.cols[k] == j
            return A.values[k]
        end
    end
    return zero(T)
end

# Implementing matrix-vector multiplication
import LinearAlgebra: mul!

function mul!(out::AbstractVector, A::SparseMatrixCSR, x::AbstractVector)
    fill!(out, 0)
    for k in 1:length(A.values)
        out[A.rows[k]] += A.values[k] * x[A.cols[k]]
    end
    out
end

# Iterative solver function
function iterate(A, b, u, a; tol = 1e-3, max_iters = 1000)
    iterations = 0
    norm_r = Inf

    while true
        r = b - A * u
        norm_r = norm(r)
        if norm_r < tol || iterations >= max_iters
            break
        end
        u += a * r
        iterations += 1
    end

    return u, iterations, norm_r
end

# Function to create a tridiagonal matrix
function create_tridiagonal_matrix(n::Int)
    rows = Int64[]
    cols = Int64[]
    values = Float64[]

    h = 1.0 / n
    a = 0.5 * h^2

    for i = 1:n+1
        if i > 1
            push!(rows, i)
            push!(cols, i-1)
            push!(values, -a)
        end

        push!(rows, i)
        push!(cols, i)
        push!(values, 2*a)

        if i < n+1
            push!(rows, i)
            push!(cols, i+1)
            push!(values, -a)
        end
    end

    return SparseMatrixCSR{Float64}(rows, cols, values)
end

# Test function to demonstrate usage
function test_and_time(n::Int)
    println("Creating tridiagonal matrix...")
    A_sparse = create_tridiagonal_matrix(n)
    b = [cos(π * i * 1/n) for i in 0:n]
    initial_u = zeros(n+1)

    println("Running iterate with sparse matrix...")
    @time final_u_sparse, _, _ = iterate(A_sparse, b, initial_u, 0.5 * (1/n)^2)

    println("Plotting results...")
    p = plot(0:n, initial_u, label="Initial u", title="Initial vs Final Solutions u", xlabel="Index", ylabel="Value of u")
    plot!(p, 0:n, final_u_sparse, label="Final u (Sparse)")
    display(p)
    savefig(p, "plot_result.png")
end

# Main execution
println("Starting test...")
test_and_time(100) # Testing with matrix size 100
println("Test completed.")


# Create a test sparse matrix
n = 5
A_sparse = create_tridiagonal_matrix(n)

# Create vectors for testing
out = zeros(n + 1)
x = rand(n + 1)
b = [cos(π * i / n) for i in 0:n]
u = zeros(n + 1)
a = 0.5 * (1 / n)^2

# Check type stability for mul!(out, A, x)
@code_warntype mul!(out, A_sparse, x)

# Check type stability for iterate(A, b, u, a)
@code_warntype iterate(A_sparse, b, u, a)
