using AtomsBase: AtomsBase
using StaticArraysCore: SVector

function NQCBase.Atoms(system::AtomsBase.AbstractSystem)
    return NQCBase.Atoms(AtomsBase.atomic_symbol(system))
end

function Cell(system::AtomsBase.AbstractSystem)
    if AtomsBase.isinfinite(system)
        return NQCBase.InfiniteCell()
    else
        box = AtomsBase.bounding_box(system)
        cell = PeriodicCell(reduce(hcat, box))
        NQCBase.set_periodicity!(cell, AtomsBase.periodicity(system))
        return cell
    end
end

function Position(system::AtomsBase.AbstractSystem)
    r = AtomsBase.position(system)
    output = zeros(AtomsBase.n_dimensions(system), Base.length(system))
    for i in axes(output, 2)
        for j in axes(output, 1)
            output[j, i] = austrip(r[i][j])
        end
    end
    return output
end

function Velocity(system::AtomsBase.AbstractSystem)
    v = AtomsBase.velocity(system)
    output = zeros(AtomsBase.n_dimensions(system), Base.length(system))
    for i in axes(output, 2)
        for j in axes(output, 1)
            output[j, i] = austrip(v[i][j])
        end
    end
    return output
end

function AtomsBase.bounding_box(cell::PeriodicCell)
    S = size(cell.vectors, 2)
    return SVector{S}(auconvert.(u"Å", vec) for vec in eachcol(cell.vectors))
end

function AtomsBase.boundary_conditions(cell::PeriodicCell)
    S = size(cell.vectors, 2)
    return SVector{S}(
        bc ? AtomsBase.Periodic() : AtomsBase.DirichletZero() for bc in cell.periodicity
    )
end
AtomsBase.isinfinite(::PeriodicCell) = false
AtomsBase.isinfinite(::InfiniteCell) = true

function System(atoms::NQCBase.Atoms, position::AbstractMatrix, cell::AbstractCell=InfiniteCell())
    output_atoms = AtomsBaseAtoms(atoms, position)
    return build_system(output_atoms, cell)
end

function System(atoms::NQCBase.Atoms, position::AbstractMatrix, velocity::AbstractMatrix, cell::AbstractCell=InfiniteCell())
    output_atoms = AtomsBaseAtoms(atoms, position, velocity)
    return build_system(output_atoms, cell)
end

function AtomsBaseAtoms(atoms::NQCBase.Atoms, position::AbstractMatrix)
    if length(atoms) != size(position, 2)
        @error atoms position
        error("The provided `Atoms` do not match the `position` array.")
    end

    output_atoms = AtomsBase.Atom[]
    sizehint!(output_atoms, length(atoms))
    for i in axes(position,2)
        r = auconvert.(u"Å", position[:,i])
        push!(output_atoms, AtomsBase.Atom(atoms.numbers[i], r))
    end
    return output_atoms
end

function AtomsBaseAtoms(atoms::NQCBase.Atoms, position::AbstractMatrix, velocity::AbstractMatrix)
    if length(atoms) != size(position, 2)
        @error atoms position
        error("The provided `Atoms` do not match the `position` array.")
    end
    if length(atoms) != size(velocity, 2)
        @error atoms velocity
        error("The provided `Atoms` do not match the `velocity` array.")
    end

    output_atoms = AtomsBase.Atom[]
    sizehint!(output_atoms, length(atoms))
    for i in axes(position,2)
        r = auconvert.(u"Å", position[:,i])
        v = auconvert.(u"Å/ps", velocity[:,i])
        push!(output_atoms, AtomsBase.Atom(atoms.numbers[i], r, v))
    end
    return output_atoms
end

function build_system(atoms, cell)
    if AtomsBase.isinfinite(cell)
        return AtomsBase.isolated_system(atoms)
    else
        box = AtomsBase.bounding_box(cell)
        bc = AtomsBase.boundary_conditions(cell)
        return AtomsBase.atomic_system(atoms, box, bc)
    end
end

function Trajectory(
    atoms::NQCBase.Atoms,
    position::Vector{<:AbstractMatrix},
    velocity::Vector{<:AbstractMatrix},
    cell::AbstractCell=InfiniteCell()
    )

    trajectory = AtomsBase.FlexibleSystem[]
    sizehint!(trajectory, length(position))

    for i in eachindex(position, velocity)
        push!(trajectory, System(atoms, position[i], velocity[i], cell))
    end

    return trajectory
end
