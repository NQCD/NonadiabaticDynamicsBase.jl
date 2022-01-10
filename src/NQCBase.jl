module NQCBase

using PeriodicTable
using Unitful, UnitfulAtomic
using StaticArrays
using Requires

include("unit_conversions.jl")
include("atoms.jl")
include("cells.jl")
include("io/extxyz.jl")

function __init__()
    @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" @eval include("io/ase.jl")
end

end
