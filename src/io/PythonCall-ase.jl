
using .PythonCall
using Unitful, UnitfulAtomic

export convert_from_ase_atoms
export convert_to_ase_atoms

const ase = pyimport("ase")

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

convert_from_ase_atoms(ase_atoms::Py) =
    Atoms(ase_atoms), positions(ase_atoms), Cell(ase_atoms)

Atoms(ase_atoms::Py) = Atoms{Float64}(Symbol.(PyList(ase_atoms.get_chemical_symbols())))

positions(ase_atoms::Py) = austrip.(PyArray(ase_atoms.get_positions())'u"Å")

function Cell(ase_atoms::Py)
    if all(PyArray(ase_atoms.cell.array) .== 0)
        return InfiniteCell()
    else
        return PeriodicCell{Float64}(austrip.(PyArray(ase_atoms.cell.array)'u"Å"), [Bool(x) for x in ase_atoms.pbc])
    end
end