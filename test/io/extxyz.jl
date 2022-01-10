
using Test
using NQCBase

atoms = Atoms([:H, :C, :O, :N])
cell = PeriodicCell(rand(3, 3) .* 10)
R = rand(3, 4) .* 10

@testset "to/from_extxyz_dict" begin
    dict = NQCBase.to_extxyz_dict(atoms, R, cell)
    @test dict["cell"] ≈ au_to_ang.(permutedims(cell.vectors, (2,1)))
    natoms, nR, ncell = NQCBase.from_extxyz_dict(dict)
    @test ncell.vectors ≈ cell.vectors
    @test ncell.inverse ≈ cell.inverse
end

@testset "Single frame" begin
    write_extxyz("output.xyz", atoms, R, cell)
    new_atoms, new_R, new_cell = read_extxyz("output.xyz")
    @test new_atoms == atoms
    @test new_cell.vectors ≈ cell.vectors
    @test new_cell.inverse ≈ cell.inverse
    @test new_cell.periodicity ≈ cell.periodicity
    @test new_R ≈ [R]
end

@testset "Multiple frames" begin
    atoms = Atoms([:H, :C, :O, :N])
    cell = PeriodicCell(rand(3, 3) .* 10)
    R = [rand(3, 4) .* 10 for _=1:100]
    write_extxyz("output.xyz", atoms, R, cell)
    new_atoms, new_R, new_cell = read_extxyz("output.xyz")
    @test new_atoms == atoms
    @test new_cell.vectors ≈ cell.vectors
    @test new_cell.inverse ≈ cell.inverse
    @test new_cell.periodicity ≈ cell.periodicity
    @test new_R ≈ R
end

rm("output.xyz")
