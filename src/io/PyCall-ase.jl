
using .PyCall
using Unitful, UnitfulAtomic

export convert_from_ase_atoms
export convert_to_ase_atoms

const ase = PyNULL()
copy!(ase, pyimport_conda("ase", "ase", "conda-forge"))

convert_to_ase_atoms(atoms::Atoms, R::Matrix) =
    ase.Atoms(positions=ustrip.(u"Å", R'u"bohr"), symbols=string.(atoms.types))

convert_to_ase_atoms(atoms::Atoms, R::Matrix, ::InfiniteCell) =
    convert_to_ase_atoms(atoms, R)

function convert_to_ase_atoms(atoms::Atoms, R::Matrix, cell::PeriodicCell)
    ase.Atoms(
        positions=ustrip.(u"Å", R'u"bohr"),
        cell=ustrip.(u"Å", cell.vectors'u"bohr"),
        symbols=string.(atoms.types),
        pbc=cell.periodicity)
end

function convert_to_ase_atoms(atoms::Atoms, R::Vector{<:Matrix}, cell::AbstractCell)
    convert_to_ase_atoms.(Ref(atoms), R, Ref(cell))
end

convert_from_ase_atoms(ase_atoms::PyObject) =
    Atoms(ase_atoms), positions(ase_atoms), Cell(ase_atoms)

Atoms(ase_atoms::PyObject) = Atoms{Float64}(Symbol.(ase_atoms.get_chemical_symbols()))

positions(ase_atoms::PyObject) = austrip.(ase_atoms.get_positions()'u"Å")

function Cell(ase_atoms::PyObject)
    if all(ase_atoms.cell.array .== 0)
        return InfiniteCell()
    else
        return PeriodicCell{Float64}(austrip.(ase_atoms.cell.array'u"Å"), ase_atoms.pbc)
    end
end
