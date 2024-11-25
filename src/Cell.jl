include("Lattice.jl")

mutable struct CellSlice
    type::String
    sites::Vector{Int64}
    in_stat::Vector{Bool}
    out_stat::Vector{Bool}
    value::Float64
end

function calculate_value(slice::CellSlice, Ω::Float64, Wbond::Vector{Matrix{Float64}})::Float64
    value = 1.
    if slice.type == "0,0"
        i, j = slice.sites
        if i != j; error("sites error") end
        if slice.in_stat[i] != slice.out_stat[i]; error("evolution error") end
        value = 1.
    elseif slice.type == "-1,a"
        i, j = slice.sites
        if i != j; error("sites error") end
        if slice.in_stat[i] == slice.out_stat[i]; error("evolution error") end
        value = Ω/2
    elseif slice.type == "1,a"
        i, j = slice.sites
        if i != j; error("sites error") end
        if slice.in_stat[i] != slice.out_stat[i]; error("evolution error") end
        value = Ω/2
    elseif slice.type == "1,b"
        i, j = slice.sites
        if i == j; error("sites error") end
        if slice.in_stat[i] != slice.out_stat[i]; error("evolution error") end
        if slice.in_stat[j] != slice.out_stat[j]; error("evolution error") end
        if slice.in_stat[i] == 0 && slice.out_stat[j] == 0; value = Wbond[1][i, j]; 
        elseif slice.in_stat[i] == 0 && slice.out_stat[j] == 1; value = Wbond[2][i, j]; 
        elseif slice.in_stat[i] == 1 && slice.out_stat[j] == 0; value = Wbond[3][i, j]; 
        elseif slice.in_stat[i] == 1 && slice.out_stat[j] == 1; value = Wbond[4][i, j]; 
        else error("unknown state") end
    else
        error("unknown operator type.")
    end
    return value
end

function update_value!(slice::CellSlice, Ω::Float64, Wbond::Vector{Matrix{Float64}})
    slice.value = calculate_value(slice, Ω, Wbond)
end

function value(slice::CellSlice)::Float64
    return slice.value
end

function flip!(slice::CellSlice, i::Int64, Ω::Float64, Wbond::Vector{Matrix{Float64}})::CellSlice
    slice.in_stat[i] = flip(slice.in_stat[i])
    slice.out_stat[i] = flip(slice.out_stat[i])
    update_value!(slice, Ω, Wbond)
    return slice
end

function flip(slice::CellSlice, i::Int64, Ω::Float64, Wbond::Vector{Matrix{Float64}})::CellSlice
    slice_copy = deepcopy(slice)
    slice_copy.in_stat[i] = flip(slice_copy.in_stat[i])
    slice_copy.out_stat[i] = flip(slice_copy.out_stat[i])
    update_value!(slice_copy, Ω, Wbond)
    return slice_copy
end

function off_diagonal_operator_flip!(slice::CellSlice, i::Int64, side::String)::CellSlice
    if slice.type == "1,a"
        slice.type = "-1,a"
        if side == "left"
            slice.in_stat[i] = flip(slice.in_stat[i])
        elseif side == "right"
            slice.out_stat[i] = flip(slice.out_stat[i])
        else
            error("side can only be left or right.")
        end
    elseif slice.type == "-1,a"
        slice.type = "1,a"
        if side == "left"
            slice.in_stat[i] = flip(slice.in_stat[i])
        elseif side == "right"
            slice.out_stat[i] = flip(slice.out_stat[i])
        else
            error("side can only be left or right.")
        end
    else
        error("this function only acts on off-diagonal operator.")
    end
    return slice
end


@inline flip(x::Int64)::Int64 = 1-x
@inline flip(x::Int8)::Int8 = Int8(1)-x
@inline flip(x::Bool)::Bool = !x

mutable struct Cell
    type::String
    M::Int64
    L::Int64
    # number of non-identitys for finite-T and 
    # number of identitys for zero-T
    n::Int64
    slices::Vector{CellSlice}
end

function Base.getindex(cell::Cell, i::Int)::CellSlice
    return cell.slices[i]
end

function Base.getindex(cell::Cell, ii::UnitRange{Int64})::Vector{CellSlice}
    return cell.slices[ii]
end

function Base.setindex!(cell::Cell, slice::CellSlice, i::Int64)
    cell.slices[i] = slice
end

function Base.lastindex(cell::Cell)::Int64
    return cell.M
end

function check(cell::Cell)::String
    if cell.type == "finite-T"
        if cell.M != length(cell.slices); error("time steps dismatch") end
        if cell.L != length(cell.slices[1].in_stat); error("space dimensions dismatch") end
        m = 0
        for i in 1:cell.M
            if cell.slices[i].type == "0,0"
                m += 1
            end
            if i != 1 && i != cell.M
                if cell.slices[i].out_stat != cell.slices[i+1].in_stat; error("state propagates improperly at $(i)") end
            end
        end
        if cell.n != cell.M - m; error("number of non-identities may be wrong") end
        if cell.slices[1].in_stat != cell.slices[cell.M].out_stat; error("pbc in time direction broken") end
        return "all proper"
    elseif cell.type == "zero-T"
        return "Not implemented yet."
    else
        error("unknown cell type")
    end
end

function t2s(q1::Bool, q2::Bool)::Int8
    if q1 == 0 && q2 == 0; return 1;
    elseif q1 == 0 && q2 == 1; return 2;
    elseif q1 == 1 && q2 == 0; return 3;
    elseif q1 == 1 && q2 == 1; return 4;
    else error("unexpected state") end
end


function finiteT_randomCell(M::Int64, lattice::Lattice, Ω::Float64, Wbond::Vector{Matrix{Float64}})::Cell
    L = lattice.L
    state = Bool.(rand(0:1, L))
    randcell = Cell("finite-T", M, L, 0, [])
    for _ in 1:M
        oper_type = ["0,0","1,a","1,b"][rand(1:3)]
        if oper_type == "0,0"
            current_slice = CellSlice(oper_type, [1,1], deepcopy(state), deepcopy(state), 1.)
            push!(randcell.slices, current_slice)
        elseif oper_type == "1,a"
            loc = rand(1:L)
            current_slice = CellSlice(oper_type, [loc,loc], deepcopy(state), deepcopy(state), Ω/2)
            push!(randcell.slices, current_slice)
            randcell.n += 1
        else
            loc1, loc2 = uniform_sample_bond(lattice)
            val = Wbond[t2s(state[loc1], state[loc2])][loc1, loc2]
            current_slice = CellSlice(oper_type, [loc1,loc2], deepcopy(state), deepcopy(state), val)
            push!(randcell.slices, current_slice)
            randcell.n += 1
        end
    end
    check(randcell)
    return randcell
end

function randomEmptyCell(M::Int64, lattice::Lattice)::Cell
    L = lattice.L
    state = Bool.(rand(0:1, L))
    randemptycell = Cell("finite-T", M, L, 0, [])
    for _ in 1:M
        push!(randemptycell.slices, CellSlice("0,0", [1,1], deepcopy(state), deepcopy(state), 1.))
    end
    check(randemptycell)
    return randemptycell
end

function zeroT_randomCell(M::Int64, lattice::Lattice, Ω::Float64, Wbond::Vector{Matrix{Float64}})::Cell
    @assert M%2 == 0 "In zero T simulation, M must be an even number."
    L = lattice.L
    state = Bool.(rand(0:1, L))
    randcell = Cell("zero-T", M, L, 0, [])
    for _ in 1:M
        oper_type = ["1,a","1,b"][rand(1:2)]
        if oper_type == "1,a"
            loc = rand(1:L)
            current_slice = CellSlice(oper_type, [loc,loc], deepcopy(state), deepcopy(state), Ω/2)
            push!(randcell.slices, current_slice)
            randcell.n += 1
        else
            loc1, loc2 = uniform_sample_bond(lattice)
            val = Wbond[t2s(state[loc1], state[loc2])][loc1, loc2]
            current_slice = CellSlice(oper_type, [loc1,loc2], deepcopy(state), deepcopy(state), val)
            push!(randcell.slices, current_slice)
        end
    end
    return randcell
end

function Base.show(io::IO, ::MIME"text/plain", cell::Cell)
    println(io, "cell.type = $(cell.type), cell.M = $(cell.M), cell.L = $(cell.L)")
    for i in 1:cell.M
        print(io, "slices[$(i)]: ")
        print(io, "type = $(lpad(cell[i].type, 4)), sites = $(cell[i].sites), value = $(cell[i].value), ")
        print(io, "in_stat = $(cell[i].in_stat), ")
        println(io, "out_stat = $(cell[i].out_stat)")
    end
end
