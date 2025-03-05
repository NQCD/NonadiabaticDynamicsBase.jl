module NQCBase

using PeriodicTable
using Unitful, UnitfulAtomic
using Requires

include("unit_conversions.jl")
include("atoms.jl")
include("cells.jl")
include("io/extxyz.jl")

function __init__() # Conditional loading of IO extensions if said packages are available. Import these before NQCBase for the interfaces to load properly. 
    @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" @eval include("io/PyCall-ase.jl")
    @require PythonCall="6099a3de-0909-46bc-b1f4-468b9a2dfc0d" @eval include("io/PythonCall-ase.jl")
end

include("atoms_base.jl")
export Cell
export System
export Trajectory
export Position, Velocity

end
