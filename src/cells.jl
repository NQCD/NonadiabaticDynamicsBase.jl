using LinearAlgebra: norm, mul!
using Distances: evaluate, PeriodicEuclidean

export AbstractCell
export PeriodicCell
export InfiniteCell
export set_periodicity!
export set_vectors!
export apply_cell_boundaries!
export evaluate_periodic_distance
export check_atoms_in_cell

const periodic_distance = PeriodicEuclidean([1, 1, 1])

abstract type AbstractCell end

struct InfiniteCell <: AbstractCell end

"""
    PeriodicCell{T<:AbstractFloat} <: AbstractCell

Optionally periodic cell
"""
struct PeriodicCell{T<:AbstractFloat} <: AbstractCell
    vectors::Matrix{T}
    inverse::Matrix{T}
    periodicity::Vector{Bool}
    tmp_vector1::Vector{T}
    tmp_vector2::Vector{T}
    tmp_bools::Vector{Bool}
    function PeriodicCell{T}(vectors::AbstractMatrix, periodicity::Vector{Bool}) where {T}
        new{T}(vectors, inv(vectors), periodicity,
            zeros(size(vectors)[1]), zeros(size(vectors)[1]), zeros(Bool, size(vectors)[1]))
    end
end

Base.eltype(::PeriodicCell{T}) where {T} = T

function PeriodicCell(vectors::AbstractMatrix)
    vectors = austrip.(vectors)
    PeriodicCell{eltype(vectors)}(vectors, [true, true, true]) 
end

function PeriodicCell(vectors::AbstractMatrix{<:Integer})
    PeriodicCell{Float64}(vectors, [true, true, true]) 
end

function set_periodicity!(cell::PeriodicCell, periodicity::AbstractVector{Bool})
    cell.periodicity .= periodicity
end

function set_vectors!(cell::PeriodicCell, vectors::AbstractMatrix)
    cell.vectors .= vectors
    cell.inverse .= inv(cell.vectors)
end

function apply_cell_boundaries!(cell::PeriodicCell, R::AbstractMatrix)
    @views for i in axes(R, 2) # atoms
        apply_cell_boundaries!(cell, R[:,i])
    end
end
apply_cell_boundaries!(::InfiniteCell, ::AbstractArray) = nothing

function apply_cell_boundaries!(cell::PeriodicCell, R::AbstractVector)
    mul!(cell.tmp_vector1, cell.inverse, R)
    for j in axes(R, 1) # DoFs
        if cell.periodicity[j]
            cell.tmp_vector1[j] = mod(cell.tmp_vector1[j], 1)
        end
    end
    mul!(R, cell.vectors, cell.tmp_vector1)
end

"""
    check_atoms_in_cell(cell::PeriodicCell, R::AbstractMatrix)::Bool

True if all atoms are inside the cell, false otherwise.
"""
function check_atoms_in_cell(cell::PeriodicCell, R::AbstractMatrix)::Bool
    @views for i in axes(R, 2) # atoms
        mul!(cell.tmp_vector1, cell.inverse, R[:,i])
        @. cell.tmp_bools = (cell.tmp_vector1 > 1) | (cell.tmp_vector1 < 0)
        any(cell.tmp_bools) && return false
    end
    true
end

function evaluate_periodic_distance(cell::PeriodicCell, r1::AbstractVector, r2::AbstractVector)
    mul!(cell.tmp_vector1, cell.inverse, r1)
    mul!(cell.tmp_vector2, cell.inverse, r2)
    evaluate(periodic_distance, cell.tmp_vector1, cell.tmp_vector2)
end