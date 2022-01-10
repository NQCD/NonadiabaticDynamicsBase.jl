using NQCBase
using Test, SafeTestsets

@time @safetestset "Atoms tests" begin include("atoms.jl") end
@time @safetestset "Cells tests" begin include("cells.jl") end
@time @safetestset "ExtXYZ tests" begin include("io/extxyz.jl") end
@time @safetestset "ase tests" begin include("io/ase.jl") end
