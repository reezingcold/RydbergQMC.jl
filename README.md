# RydbergQMC.jl
A Julia package for doing quantum Monte Carlo in Rydberg atom systems.

## Introduction
This package is written based on [this paper](https://scipost.org/SciPostPhysCore.7.2.016/pdf). It is better to read the paper before using 
this package. 
The basic functionality of this package has been completed but many interfaces are not implemented yet so at this stage this package is not 
friendly to use. 

## Dependencies
This package requires [StatsBase.jl](https://juliastats.org/StatsBase.jl/stable/). Other packages, including 
[PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl), [PyFormattedStrings](https://github.com/JuliaAPlavin/PyFormattedStrings.jl) 
are used in examples but they are not necessary for using RydbergQMC.jl.

## Finite Temperature Stochastic Series Expansion (finite_temperature_sse)
Here is a basic example.
```Julia
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
```
The ground energy calculated by QMC is getting closer to ED result as the samples increase.
![finiteT_result](https://github.com/reezingcold/RydbergQMC.jl/blob/main/figures/mc_finiteT_test.png)

## Zero Temperature Stochastic Series Expansion (zero_temperature_sse)
Here is a basic example.
```Julia
include("src/ZeroTemperatureSSE.jl")
using PyPlot
using PyFormattedStrings
using Random
using LinearAlgebra
Random.seed!(114514)

L = 11

x, y = 2., 1.4
C₆, d = 1000., 2.
V = C₆/d^6
@show Ω = V/y^6
δ = x * Ω
lattice = rydberg_chain_obc(L, V)
ϵ = 0.001
M = 700
M *= 2
Niter = 80000

# about 3s
@time begin
H = Hamiltonian(lattice, Ω, δ)
eigsys = eigen(H)
ground_energy = eigsys.values[1]
end

C = sum(construct_Cmat(lattice, δ, ϵ))

function q(cell::Cell)::Float64
    return cell.n+1
end

# about 5min
@time mc_result = zero_temperature_sse(
    lattice, 
    (Ω=Ω, δ=δ, ϵ=ϵ), 
    (M=M, Niter=Niter), 
    q
)

function energy(x::Int64)::Float64
    aver_q = sum(mc_result[1:x])/x
    return (-L * (M+1) * Ω/2)/aver_q + C + Ω/2 * L
end

# about 4s
@time begin
plt.scatter(collect(1:Niter), energy.(collect(1:Niter)), marker=".", s=4)
# plt.scatter(collect(1:Niter), mc_result, marker=".", s=4)
plt.plot([1, Niter], [ground_energy, ground_energy], color="orangered")

plt.xlabel("#sample")
plt.ylabel("energy")
plt.savefig(f"mc0_nbar_x{x}_y{y}_M{M}_ϵ{ϵ}.png", dpi=300)
plt.close()
end
```
The ground energy calculated by Projective QMC is getting closer to ED result as the samples increase.
![zeroT_result](https://github.com/reezingcold/RydbergQMC.jl/blob/main/figures/mc_zeroT_test.png)






