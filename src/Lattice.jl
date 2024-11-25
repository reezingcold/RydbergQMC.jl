mutable struct Lattice
    # system size
    L::Int64
    # adjacency matrix
    M_ad::Matrix{Int64}
    # interaction strength matrix
    M_int::Matrix{Float64}
    # number of bonds
    Nbond::Int64
    # bonds and indices
    bonds::Vector{Tuple{Int64, Int64}}
end

Lattice(L::Int64) = Lattice(L, zeros(Int64, L, L), zeros(L, L), 0, [])

function add_bond!(lattice::Lattice, node1::Int64, node2::Int64, V12::Float64)
    lattice.M_ad[node1, node2] = 1
    lattice.M_ad[node2, node1] = 1
    lattice.M_int[node1, node2] = V12
    lattice.M_int[node2, node1] = V12
    lattice.Nbond += 1
    if node1 < node2
        push!(lattice.bonds, (node1, node2))
    else
        push!(lattice.bonds, (node2, node1))
    end
end

function check(lattice::Lattice)
    if lattice.Nbond == length(lattice.bonds) && lattice.Nbond == sum(lattice.M_ad)÷2
        return "all proper"
    else
        error("lattice constructed improperly")
    end
end

# function add_onsite_detuning(lattice::Lattice, node::Int64, V::Float64)
#     lattice.M_int[node, node] = V
# end

function adjacent_sites(lattice::Lattice, node1::Int64)::Vector{Int64}
    nodes_adjc = []
    for i in 1:lattice.L
        if lattice.M_ad[i, node1] == 1
            push!(nodes_adjc, i)
        end
    end
    return nodes_adjc
end

function uniform_sample_bond(lattice::Lattice)::Tuple{Int64, Int64}
    bond_ind = rand(1:lattice.Nbond)
    return lattice.bonds[bond_ind]
end

function Hamiltonian(lattice::Lattice, Ω::Float64, δ::Float64)::Matrix{Float64}
    L = lattice.L
    @assert L <= 12 "this function does not support system size greater than 12."
    qeye = [1. 0.; 0. 1.]
    n = [0. 0.; 0. 1.]
    sigmax = [0. 1.; 1. 0.]
    H = zeros(Float64, 2^L, 2^L)
    for i in 1:L
        cup = repeat([qeye], L)
        cup[i] = Ω/2 * sigmax - δ * n
        H += kron(cup...)
    end
    for i in 1:L
        for j in 1:L
            if lattice.M_ad[i, j] == 1
                cup = repeat([qeye], L)
                cup[i], cup[j] = n, n 
                H += lattice.M_int[i, j] * kron(cup...)
            end
        end
    end
    return H
end

function rydberg_chain_obc(L::Int64, V::Float64)::Lattice
    max_d::Float64=L+0.1
    chain_obc = Lattice(L)
    for i in 1:L-1
        for j in i+1:L
            if abs(i-j) < max_d
                add_bond!(chain_obc, i, j, V/(i-j)^6)
            end
        end
    end
    return chain_obc
end

function rydberg_square_obc(Lx::Int64, Ly::Int64, V::Float64)::Lattice
    max_d::Float64=Lx+Ly+0.1
    coo2ind(x::Int64, y::Int64)::Int64 = (y-1) * Lx + x
    ind2coo(s::Int64)::Tuple{Int64, Int64} = (s % Lx == 0) ? (Lx, s÷Lx) : (s%Lx, s÷Lx+1)
    L = Lx * Ly
    square_obc = Lattice(L)
    for i in 1:L-1
        for j in i+1:L
            ix, iy = ind2coo(i)
            jx, jy = ind2coo(j)
            dij = sqrt((ix-jx)^2+(iy-jy)^2)
            if dij < max_d
                add_bond!(square_obc, i, j, V/dij^6)
            end
        end
    end
    return square_obc
end