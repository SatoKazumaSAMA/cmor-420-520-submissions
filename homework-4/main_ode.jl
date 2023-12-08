using DifferentialEquations, LinearAlgebra, Plots, SparseArrays

struct SparseMatrixCSR{T} <: AbstractMatrix{T}
    rows::Vector{Int64}
    cols::Vector{Int64}
    values::Vector{T}
end

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

function mul!(y::Vector{T}, A::SparseMatrixCSR{T}, x::Vector{T}) where T
    fill!(y, zero(T))
    for k in 1:length(A.values)
        y[A.rows[k]] += A.values[k] * x[A.cols[k]]
    end
    return y
end

function create_tridiagonal_matrix(n::Int)
    rows = []
    cols = []
    values = []
    h = 1.0 / n
    a = -2.0 * h^2
    for i = 1:n
        push!(rows, i)
        push!(cols, i)
        push!(values, a)
        if i > 1
            push!(rows, i)
            push!(cols, i-1)
            push!(values, h^2)
        end
        if i < n
            push!(rows, i)
            push!(cols, i+1)
            push!(values, h^2)
        end
    end
    return SparseMatrixCSR{Float64}(rows, cols, values)
end

function rhs!(du, u, p, t)
    mul!(du, p.A, u)
    du .*= -1
    du .+= p.b
end

n = 100
h = 1.0 / n
A = create_tridiagonal_matrix(n)
b = [cos(π * i / n) for i in 1:n]  
u0 = zeros(n)

p = (A = A, b = b, h = h)

tspan = (0.0, 1.0)

eigen_estimation = integrator -> 2 / (integrator.p.h^2)

prob = ODEProblem(rhs!, u0, tspan, p)

solvers = [
    Tsit5(),
    SSPRK33(),
    ROCK4(eigen_est = eigen_estimation)
]

animations = Dict()

for solver in solvers
    sol = solve(prob, solver, dt = h^2, saveat = 0.1, abstol = 1e-6, reltol = 1e-3)
    
    anim = @animate for (i, u) in enumerate(sol.u)
        plot(u, ylims=(-1, 1), label="u(t)", title="Solver: $(typeof(solver)) at time t=$(round(sol.t[i], digits=2))")
    end
    
    animation_filename = "animation_$(typeof(solver)).gif"
    gif(anim, animation_filename, fps = 15)
    animations[typeof(solver)] = animation_filename

    println("Solver: ", typeof(solver))
    println("Number of accepted steps: ", sol.destats.naccept)
    println("Number of rejected steps: ", sol.destats.nreject)
    println("Number of function evaluations: ", sol.destats.nf)
    
    println("Solver: ", typeof(solver))
    println("Number of time steps: ", length(sol.t))
end

println("ODE solving and animations completed.")

