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