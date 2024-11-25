import StatsBase
include("Cell.jl")

function construct_Cmat(lattice::Lattice, δ::Float64, ϵ::Float64)::Matrix{Float64}
    L = lattice.L
    V = lattice.M_int
    δb = δ/(L-1)
    Cmat = zeros(L, L)
    for i in 1:L-1
        for j in i+1:L
            Cmat[i, j] += abs(min(0, δb, 2*δb-V[i, j])) + ϵ*abs(min(δb, 2*δb-V[i, j]))
        end
    end
    return Cmat
end

function construct_Wbond(lattice::Lattice, δ::Float64, Cmat::Matrix{Float64})::Vector{Matrix{Float64}}
    L = lattice.L
    V = lattice.M_int
    δb = δ/(L-1)
    W1 = Cmat
    W2 = Cmat .+ δb
    W3 = Cmat .+ δb
    W4 = -V + Cmat .+ 2*δb
    return [W1, W2, W3, W4]
end

function construct_Pmat(Wbond::Vector{Matrix{Float64}}, Ω::Float64)::Tuple{Matrix{Float64},Float64,Vector{Tuple{Int64,Int64}},Vector{Float64}}
    L = size(Wbond[1])[1]
    Pmat = zeros(L, L)
    mtcN = 0.
    vec_sites = []
    vec_weights = []
    for i in 1:L
        for j in i:L
            if i == j
                Pmat[i, j] = Ω/2
                mtcN += Ω/2
                push!(vec_sites, (i, j))
                push!(vec_weights, Ω/2)
            else
                Pmat[i, j] = max([Wbond[x][i, j] for x in 1:4]...)
                mtcN += Pmat[i, j]
                push!(vec_sites, (i, j))
                push!(vec_weights, Pmat[i, j])
            end
        end
    end
    return Pmat, mtcN, vec_sites, vec_weights
end

function zeroT_diagonal_update!(
    cell::Cell, 
    vec_sites::Vector{Tuple{Int64, Int64}}, 
    vec_weights::Vector{Float64}, 
    Pmat::Matrix{Float64}, 
    Ω::Float64, 
    Wbond::Vector{Matrix{Float64}}, 
    )::Cell
    for p in 1:cell.M
        if cell[p].type == "1,a" || cell[p].type == "1,b"
            while true
                loc1, loc2 = StatsBase.sample(vec_sites, StatsBase.pweights(vec_weights))
                if loc1 == loc2
                    if cell[p].type == "1,b"; cell.n += 1 end
                    cell[p].type = "1,a"
                    cell[p].sites = [loc1, loc2]
                    cell[p].value = Ω/2
                    break
                else # loc1 != loc2
                    q1, q2 = cell[p].in_stat[loc1], cell[p].in_stat[loc2]
                    Wactual = Wbond[t2s(q1, q2)][loc1, loc2]
                    Wsampled = Pmat[loc1, loc2]
                    if rand() < Wactual/Wsampled
                        if cell[p].type == "1,a"; cell.n -= 1 end
                        cell[p].type = "1,b"
                        cell[p].sites = [loc1, loc2]
                        cell[p].value = Wactual
                        break
                    end
                end
            end
        elseif cell[p].type == "-1,a"
            loc = cell[p].sites[1]
            cell[p].value = Ω/2
            if cell[p].in_stat[loc] == cell[p].out_stat[loc]
                # flip
                cell[p].out_stat[loc] = flip(cell[p].in_stat[loc])
            end
        else
            error("unexpected operator type")
        end
    end
    return cell
end

function zeroT_accept_prob_for_cluster_flip(p0::Int64, p1::Int64, phyind::Int64, 
    cell::Cell, 
    Ω::Float64, 
    Wbond::Vector{Matrix{Float64}})::Float64
    if abs(p0-p1) == 1
        return 1.
    else 
        W = prod(value.(cell[p0+1:p1-1]))
        Wp = 1.
        for p in p0+1:p1-1
            if cell[p].type == "1,b" && (phyind in cell[p].sites)
                Wp *= value(flip(cell[p], phyind, Ω, Wbond))
            end
        end
        return min(1, Wp/W)
    end
end

function zeroT_cluster_update_on_a_line!(cell::Cell, 
    p0::Int64, p1::Int64, phyind::Int64, 
    Ω::Float64, Wbond::Vector{Matrix{Float64}})::Cell
    if p0 == 0 && p1 == cell.M+1
        for p in 1:cell.M; flip!(cell[p], phyind, Ω, Wbond) end
    elseif p0 == 0 && p1 != cell.M+1
        cell[p1].type == "1,a" ? cell.n -= 1 : cell.n += 1
        for p in 1:p1-1; flip!(cell[p], phyind, Ω, Wbond) end
        off_diagonal_operator_flip!(cell[p1], phyind, "left")
    elseif p0 != 0 && p1 == cell.M+1
        cell[p0].type == "1,a" ? cell.n -= 1 : cell.n += 1
        for p in p0+1:cell.M; flip!(cell[p], phyind, Ω, Wbond) end
        off_diagonal_operator_flip!(cell[p0], phyind, "right")
    else
        cell[p0].type == "1,a" ? cell.n -= 1 : cell.n += 1
        cell[p1].type == "1,a" ? cell.n -= 1 : cell.n += 1
        off_diagonal_operator_flip!(cell[p0], phyind, "right")
        off_diagonal_operator_flip!(cell[p1], phyind, "left")
        for p in p0+1:p1-1; flip!(cell[p], phyind, Ω, Wbond) end
    end
    return cell
end


function zeroT_line_cluster_update!(cell::Cell, Ω::Float64, Wbond::Vector{Matrix{Float64}}
    )::Cell
    pvec_lst::Vector{Vector{Int64}} = repeat([[]], cell.L)
    for p in 1:cell.M
        if cell[p].type == "1,a" || cell[p].type == "-1,a"
            push!(pvec_lst[cell[p].sites[1]], p)
        end
    end
    for i in 1:cell.L
        Lp = length(pvec_lst[i])
        if Lp == 0
            if rand() < zeroT_accept_prob_for_cluster_flip(0, cell.M+1, i, cell, Ω, Wbond)
                zeroT_cluster_update_on_a_line!(cell, 0, cell.M+1, i, Ω, Wbond)
            end
        else
            if rand() < zeroT_accept_prob_for_cluster_flip(0, pvec_lst[i][1], i, cell, Ω, Wbond)
                zeroT_cluster_update_on_a_line!(cell, 0, pvec_lst[i][1], i, Ω, Wbond)
            end
            if rand() < zeroT_accept_prob_for_cluster_flip(pvec_lst[i][Lp], cell.M+1, i, cell, Ω, Wbond)
                zeroT_cluster_update_on_a_line!(cell, pvec_lst[i][Lp], cell.M+1, i, Ω, Wbond)
            end
            for j in 1:Lp-1
                if rand() < zeroT_accept_prob_for_cluster_flip(pvec_lst[i][j], pvec_lst[i][j+1], i, cell, Ω, Wbond)
                    zeroT_cluster_update_on_a_line!(cell, pvec_lst[i][j], pvec_lst[i][j+1], i, Ω, Wbond)
                end
            end
        end
    end
    return cell
end

function zero_temperature_sse(
    lattice::Lattice, 
    model_params::NamedTuple, 
    sampling_params::NamedTuple, 
    observe_func::Function
    )::Vector{Float64}
    M = sampling_params.M
    @assert M%2 == 0 "in zero-T simulation, M must be an even number"
    Niter = sampling_params.Niter
    Ω = model_params.Ω
    Cmat = construct_Cmat(lattice, model_params.δ, model_params.ϵ)
    Wbond = construct_Wbond(lattice, model_params.δ, Cmat)
    Pmat, mtcN, vec_sites, vec_weights = construct_Pmat(Wbond, Ω)
    cell = zeroT_randomCell(M, lattice, Ω, Wbond)
    @assert cell.type == "zero-T" "unexpected cell.type == $(cell.type)"
    observable_data = zeros(Niter)
    for j in 1:Niter
        zeroT_diagonal_update!(cell, vec_sites, vec_weights, Pmat, Ω, Wbond)
        zeroT_line_cluster_update!(cell, Ω, Wbond)
        observable_data[j] = observe_func(deepcopy(cell))
        # println(j, " | ", sum(observable_data)/j)
    end
    return observable_data
end

function zero_temperature_sse(
    lattice::Lattice, 
    model_params::NamedTuple, 
    sampling_params::NamedTuple, 
    observe_func_lst::Vector{Function}
    )::Matrix{Float64}
    M = sampling_params.M
    @assert M%2 == 0 "in zero-T simulation, M must be an even number"
    Niter = sampling_params.Niter
    Ω = model_params.Ω
    Cmat = construct_Cmat(lattice, model_params.δ, model_params.ϵ)
    Wbond = construct_Wbond(lattice, model_params.δ, Cmat)
    Pmat, mtcN, vec_sites, vec_weights = construct_Pmat(Wbond, Ω)
    cell = zeroT_randomCell(M, lattice, Ω, Wbond)
    @assert cell.type == "zero-T" "unexpected cell.type == $(cell.type)"
    flen = length(observe_func_lst)
    observable_data = zeros(Niter, flen)
    for j in 1:Niter
        zeroT_diagonal_update!(cell, vec_sites, vec_weights, Pmat, Ω, Wbond)
        zeroT_line_cluster_update!(cell, Ω, Wbond)
        cell_copy = deepcopy(cell)
        for i in 1:flen; observable_data[j, i] = observe_func_lst[i](cell_copy) end
        # println(j, " | ", Niter)
    end
    return observable_data
end
