
function WriteTB(TB::AbstractTightBindModel, folder::AbstractString; mode = "w", format = "Cartesian", comment = "From LatticeModel.WriteTB.")
	WriteHR(HR(TB), joinpath(folder, "hr.dat"); mode, comment)
	WritePOSCAR(Cell(TB), joinpath(folder, "POSCAR"); mode, format, comment)
	WriteORBITAL(ORBITAL(TB), joinpath(folder, "centres.xyz"); mode, comment)
end
