include("Cell.jl")
import StatsBase

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

function finiteT_diagonal_update!(
    cell::Cell, 
    mtcN::Float64, 
    vec_sites::Vector{Tuple{Int64, Int64}}, 
    vec_weights::Vector{Float64}, 
    Pmat::Matrix{Float64}, 
    Wbond::Vector{Matrix{Float64}}, 
    model_params::NamedTuple
    )::Cell
    Ω, β = model_params.Ω, model_params.β
    for p in 1:cell.M
        if cell[p].type == "1,a" || cell[p].type == "1,b"
            if rand() < min((cell.M-cell.n+1)/(β*mtcN), 1)
                cell[p].type = "0,0"
                cell[p].sites = [1, 1]
                cell[p].value = 1.
                cell.n -= 1
            end
        elseif cell[p].type == "0,0"
            if rand() < min((β*mtcN)/(cell.M-cell.n), 1)
                loc1, loc2 = StatsBase.sample(vec_sites, StatsBase.pweights(vec_weights))
                if loc1 == loc2
                    cell[p].type = "1,a"
                    cell[p].sites = [loc1, loc2]
                    cell[p].value = Ω/2
                    cell.n += 1
                else # loc1 != loc2
                    q1, q2 = cell[p].in_stat[loc1], cell[p].in_stat[loc2]
                    Wactual = Wbond[t2s(q1, q2)][loc1, loc2]
                    Wsampled = Pmat[loc1, loc2]
                    if rand() < Wactual/Wsampled
                        cell[p].type = "1,b"
                        cell[p].sites = [loc1, loc2]
                        cell[p].value = Wactual
                        cell.n += 1
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

function accept_prob_for_cluster_flip(p0::Int64, 
    phyind::Int64, 
    cell::Cell, 
    Ω::Float64, 
    Wbond::Vector{Matrix{Float64}})::Float64
    W = prod(value.(cell[1:end]))/value(cell[p0])
    Wp = 1.
    for p in 1:cell.M
        if p == p0
            nothing
        else
            if cell[p].type == "1,b" && (phyind in cell[p].sites)
                Wp *= value(flip(cell[p], phyind, Ω, Wbond))
            end
        end
    end
    return min(1, Wp/W)
end

function accept_prob_for_cluster_flip(p0::Int64, p1::Int64, 
    phyind::Int64, 
    cell::Cell, 
    Ω::Float64, 
    Wbond::Vector{Matrix{Float64}}, 
    which::String="in")::Float64
    if which == "in"
        if abs(p0-p1) == 1
            return 1.
        else
            W = prod(value.(cell[p0+1:p1-1]))
            # values if flipped
            Wp = 1.
            for p in p0+1:p1-1
                if cell[p].type == "1,b" && (phyind in cell[p].sites)
                    Wp *= value(flip(cell[p], phyind, Ω, Wbond))
                end
            end
            return min(1, Wp/W)
        end
    elseif which == "out"
        if p0 == 1; Wl = 1.; else; Wl = prod(value.(cell[1:p0-1])) end
        if p1 == cell.M; Wr = 1.; else; Wr = prod(value.(cell[p1+1:cell.M])) end
        W = Wl * Wr
        # next calculate Wp_left and Wp_right
        Wp_left, Wp_right = 1., 1.
        for p in 1:p0-1
            if cell[p].type == "1,b" && (phyind in cell[p].sites)
                Wp_left *= value(flip(cell[p], phyind, Ω, Wbond))
            end
        end
        for p in p1+1:cell.M
            if cell[p].type == "1,b" && (phyind in cell[p].sites)
                Wp_right *= value(flip(cell[p], phyind, Ω, Wbond))
            end
        end
        Wp = Wp_left * Wp_right
        return min(1, Wp/W)
    else
        error("unsupported keyword.")
    end
end

function finiteT_cluster_update_on_a_line!(cell::Cell, 
    p0::Int64, p1::Int64, phyind::Int64, 
    Ω::Float64, Wbond::Vector{Matrix{Float64}}, which::String="in")::Cell
    if cell[p0].type != "-1,a" && cell[p0].type != "1,a"
        error("wrong cluster edge at $(p0), edge operator type is $(cell[p0].type)") 
    end
    if cell[p1].type != "-1,a" && cell[p1].type != "1,a"
        error("wrong cluster edge at $(p1), edge operator type is $(cell[p1].type)") 
    end
    if p0 == p1
        for p in 1:cell.M; flip!(cell[p], phyind, Ω, Wbond) end
    else
        if which == "in"
            off_diagonal_operator_flip!(cell[p0], phyind, "right")
            off_diagonal_operator_flip!(cell[p1], phyind, "left")
            for p in p0+1:p1-1; flip!(cell[p], phyind, Ω, Wbond) end
        elseif which == "out"
            off_diagonal_operator_flip!(cell[p0], phyind, "left")
            off_diagonal_operator_flip!(cell[p1], phyind, "right")
            for p in 1:p0-1; flip!(cell[p], phyind, Ω, Wbond) end
            for p in p1+1:cell.M; flip!(cell[p], phyind, Ω, Wbond) end
        else
            error("which parameter can only be in or out.")
        end
    end
    return cell
end

# function line_cluster_update!(cell::Cell, phyind::Int64, Ω::Float64, Wbond::Vector{Matrix{Float64}})::Cell
#     plst = []
#     for p in 1:cell.M
#         if phyind in cell[p].sites
#             if cell[p].type == "1,a" || cell[p].type == "-1,a"
#                 push!(plst, p)
#             end
#         end
#     end

#     Lp = length(plst)
#     if Lp == 0
#         if rand() < accept_prob_for_cluster_flip(0, cell.M+1, phyind, cell, Ω, Wbond)
#             for p in 1:cell.M; flip!(cell[p], phyind, Ω, Wbond) end
#         end
#     elseif Lp == 1
#         if rand() < accept_prob_for_cluster_flip(plst[1], phyind, cell, Ω, Wbond)
#             cluster_update_on_a_line!(cell, plst[1], plst[1], phyind, Ω, Wbond)
#         end
#     else
#         for i in 1:Lp-1
#             p0, p1 = plst[i], plst[i+1]
#             if rand() < accept_prob_for_cluster_flip(p0, p1, phyind, cell, Ω, Wbond)
#                 cluster_update_on_a_line!(cell, p0, p1, phyind, Ω, Wbond)
#             end
#         end
#         p0, p1 = plst[1], plst[Lp]
#         if rand() < accept_prob_for_cluster_flip(p0, p1, phyind, cell, Ω, Wbond, "out")
#             cluster_update_on_a_line!(cell, p0, p1, phyind, Ω, Wbond, "out")
#         end
#     end
#     return cell
# end

function finiteT_line_cluster_update!(cell::Cell, phyind::Int64, plst::Vector{Int64}, Ω::Float64, Wbond::Vector{Matrix{Float64}})::Cell
    Lp = length(plst)
    if Lp == 0
        if rand() < accept_prob_for_cluster_flip(0, cell.M+1, phyind, cell, Ω, Wbond)
            for p in 1:cell.M; flip!(cell[p], phyind, Ω, Wbond) end
        end
    elseif Lp == 1
        if rand() < accept_prob_for_cluster_flip(plst[1], phyind, cell, Ω, Wbond)
            finiteT_cluster_update_on_a_line!(cell, plst[1], plst[1], phyind, Ω, Wbond)
        end
    else
        for i in 1:Lp-1
            p0, p1 = plst[i], plst[i+1]
            if rand() < accept_prob_for_cluster_flip(p0, p1, phyind, cell, Ω, Wbond)
                finiteT_cluster_update_on_a_line!(cell, p0, p1, phyind, Ω, Wbond)
            end
        end
        p0, p1 = plst[1], plst[Lp]
        if rand() < accept_prob_for_cluster_flip(p0, p1, phyind, cell, Ω, Wbond, "out")
            finiteT_cluster_update_on_a_line!(cell, p0, p1, phyind, Ω, Wbond, "out")
        end
    end
    return cell
end

function finiteT_all_cluster_update!(cell::Cell, Ω::Float64, Wbond::Vector{Matrix{Float64}})::Cell
    pvec_lst::Vector{Vector{Int64}} = repeat([[]], cell.L)
    for p in 1:cell.M
        if cell[p].type == "1,a" || cell[p].type == "-1,a"
            push!(pvec_lst[cell[p].sites[1]], p)
        end
    end
    for i in 1:cell.L
        finiteT_line_cluster_update!(cell, i, pvec_lst[i], Ω, Wbond)
    end
    return cell
end

# model_params = (β=, Ω=, δ=, ϵ=)
# sampling_params = (M=, Niter=, )
function finite_temperature_sse(
    lattice::Lattice, 
    model_params::NamedTuple,
    sampling_params::NamedTuple, 
    observe_func::Function
    )::Vector{Float64}
    M = sampling_params.M
    Niter = sampling_params.Niter
    Ω = model_params.Ω
    Cmat = construct_Cmat(lattice, model_params.δ, model_params.ϵ)
    Wbond = construct_Wbond(lattice, model_params.δ, Cmat)
    Pmat, mtcN, vec_sites, vec_weights = construct_Pmat(Wbond, Ω)
    # cell = randomEmptyCell(M, lattice)
    cell = finiteT_randomCell(M, lattice, Ω, Wbond)
    @assert cell.type == "finite-T" "unexpected cell.type == $(cell.type)"
    observable_data = zeros(Niter)
    for j in 1:Niter
        finiteT_diagonal_update!(cell, mtcN, vec_sites, vec_weights, Pmat, Wbond, model_params)
        finiteT_all_cluster_update!(cell, Ω, Wbond)
        # line_cluster_update!(cell, j%lattice.L+lattice.L, Ω, Wbond)
        observable_data[j] = observe_func(deepcopy(cell))
    end
    return observable_data
end

function finite_temperature_sse(
    lattice::Lattice, 
    model_params::NamedTuple,
    sampling_params::NamedTuple, 
    observe_func_vec::Vector{Function}
    )::Matrix{Float64}
    M = sampling_params.M
    Niter = sampling_params.Niter
    Ω = model_params.Ω
    Cmat = construct_Cmat(lattice, model_params.δ, model_params.ϵ)
    Wbond = construct_Wbond(lattice, model_params.δ, Cmat)
    Pmat, mtcN, vec_sites, vec_weights = construct_Pmat(Wbond, Ω)
    # cell = randomEmptyCell(M, lattice)
    cell = randomCell(M, lattice, Ω, Wbond)
    flen = length(observe_func_vec)
    @assert cell.type == "finite-T" "unexpected cell.type == $(cell.type)"
    observation_matrix = zeros(Niter, flen)
    for j in 1:Niter
        finiteT_diagonal_update!(cell, mtcN, vec_sites, vec_weights, Pmat, Wbond, model_params)
        finiteT_all_cluster_update!(cell, Ω, Wbond)
        # line_cluster_update!(cell, j%lattice.L+lattice.L, Ω, Wbond)
        cell_copy = deepcopy(cell)
        for k in 1:flen; observation_matrix[j, k] = observe_func_vec[k](cell_copy) end
    end
    return observation_matrix
end