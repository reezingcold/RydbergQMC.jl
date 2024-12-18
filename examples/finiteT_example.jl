include("src/FiniteTemperatureSSE.jl")
using PyPlot
using PyFormattedStrings
using LinearAlgebra

β = 0.7
L = 11
x, y = 3., 1.2
C₆, d = 1000., 2.
V = C₆/d^6
Ω = V/y^6
δ = x * Ω
ϵ = 0.01
M = 300
Niter = 80000
lattice = rydberg_chain_obc(L, V)

function exact_diag(lattice::Lattice, β::Float64)::Float64
    H = Hamiltonian(lattice, Ω, δ)
    rho = exp(-β * H)
    Z = tr(rho)
    free_energy = tr(rho * H)/Z
    return free_energy
end

function qmc(lattice::Lattice, β::Float64)::Vector{Float64}
    @inline measure(cell::Cell)::Float64 = cell.n
    mc_result = finite_temperature_sse(
    lattice, 
    (β=β, Ω=Ω, δ=δ, ϵ=ϵ), 
    (M=M, Niter=Niter), 
    measure
    )
    return mc_result
end

function n_average(lst::Vector{Float64})::Vector{Float64}
    tmp = 0.
    tmp2 = []
    C = sum(construct_Cmat(lattice, δ, ϵ))
    for i in 1:length(lst)
        tmp = sum(lst[1:i])/i
        push!(tmp2, -(tmp)/β+C+Ω/2*L)
    end
    return tmp2
end

# running roughly 60s on an Apple M1-Max MacBook pro
@time begin
    plt.scatter(collect(1:Niter), n_average(qmc(lattice, β)), marker=".", s=3, label="QMC")
    plt.plot([1, Niter], repeat([exact_diag(lattice, β)], 2), color="orangered", label="ED")
    plt.xlabel("#sample")
    plt.ylabel("energy")
    plt.legend()
    plt.xlim(1, Niter)
    plt.savefig("mc_test.png", dpi=300)
    plt.close()
end
