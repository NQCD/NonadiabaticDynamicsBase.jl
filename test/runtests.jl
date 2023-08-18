using NQCBase
using Test, SafeTestsets

@safetestset "Atoms tests" begin include("atoms.jl") end
@safetestset "Cells tests" begin include("cells.jl") end
@safetestset "ExtXYZ tests" begin include("io/extxyz.jl") end
@safetestset "ase tests" begin include("io/ase.jl") end
@safetestset "AtomsBase tests" begin include("atoms_base.jl") end
