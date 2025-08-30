
"""

"""
function Writewin(file::AbstractString; mode = "w", comment = "From LatticeModel.Writewin.",
	num_iter,
	num_print_cycles = 20,
	search_shells = 100,
	shell_list = nothing,
	dis_num_iter,
	dis_conv_tol = 1e-12,
	dis_conv_window = 5,
	exclude_bands = nothing,
	num_bands,
	dis_win_min = 0,
	dis_win_max = 0,
	dis_froz_min = 0,
	dis_froz_max = 0,
	num_wann,
	use_bloch_phases = false,
	spinors = false,
	guiding_centres = false,
	projections = nothing,
	wannier_plot = false,
	wannier_plot_format = "xcrysden",
	wannier_plot_supercell = nothing,
	bands_plot = false,
	bands_num_points = 30,
	bands_plot_format = "gnuplot xmgrace",
	kpoint_path = nothing,
	cell::Cell,
	grid,
	kpoints = nothing,
)

	if num_bands < num_wann
		error("num_bands should be equalto or larger than num_wann.")
	end

	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = "#" * comment * " " * now

	path = dirname(file)
	mkpath(path)
	file = open(file, mode)
	write(file, comment * "\n")

	println(
		file,
		"
#Iteration steps to get minimized WFs.
num_iter         = $num_iter
num_print_cycles = $num_print_cycles
	",
	)

	println(
		file,
		"
#Output file, can open them only after running vasp.
write_hr  = true
write_xyz = true
	",
	)

	println(
		file,
		"
#Reciprocal space shells used by finite difference.
search_shells = $search_shells
	",
	)
	if !isnothing(shell_list)
		#TODO
		error("To be continued.")
		println(file, "shell_list =	")
	end

	println(
		file,
		"
#disentanglement parameters
dis_num_iter    = $dis_num_iter
dis_conv_tol    = $dis_conv_tol
dis_conv_window = $dis_conv_window
	",
	)

	if !isnothing(exclude_bands)
		#TODO
		error("To be continued.")
	end
	println(
		file,
		"
num_bands     = $num_bands
dis_win_min   = $dis_win_min
dis_win_max   = $dis_win_max
dis_froz_min  = $dis_froz_min
dis_froz_max  = $dis_froz_max
	",
	)

	println(file, "num_wann  = $num_wann")
	if use_bloch_phases
		println(file, "use_bloch_phases  = true")
	end
	if spinors
		println(file, "spinors  = true")
	end
	if guiding_centres
		println(file, "guiding_centres  = true")
	end

	if !isnothing(projections)
		@printf(file, "begin projections\n")
		for projection in projections
			if projection[1] isa ReducedCoordinates
				@printf(file, "f=%.6f,%.6f,%.6f:%s\n", projection[1]..., projection[2])
			elseif projection[1] isa CartesianCoordinates
				@printf(file, "c=%.6f,%.6f,%.6f:%s\n", projection[1]..., projection[2])
			end
		end
		@printf(file, "end projections\n\n")
	end

	if wannier_plot
		#TODO
		error("To be continued.")
		println(
			file,
			"
#write out Wannier functions.
wannier_plot           = true
wannier_plot_format    = $wannier_plot_format
wannier_plot_supercell = $wannier_plot_supercell
	",
		)
	end

	if bands_plot
		println(
			file,
			"
#Bands plot, to be continued.
bands_plot         = true
bands_num_points   = $bands_num_points
bands_plot_format  = $bands_plot_format
begin kpoint_path",
		)
		nkpoints = length(kpoint_path)
		for i in 1:nkpoints-1
			@printf(file, "%-5s %9.6f %9.6f %9.6f ", kpoint_path[i][1], kpoint_path[i][2]...)
			@printf(file, "%-5s %9.6f %9.6f %9.6f \n", kpoint_path[i+1][1], kpoint_path[i+1][2]...)
		end
		println(file, "end kpoint_path\n")
	end

	@printf(file, "begin unit_cell_cart\n")
	for i in 1:3
		@printf(file, "%14.7f %14.7f %14.7f \n", cell.lattice[:, i]...)
	end
	@printf(file, "end unit_cell_cart\n\n")

	if typeof(cell).parameters[1] <: ReducedCoordinates
		@printf(file, "begin atoms_frac\n")
		for i in eachindex(cell.location)
			@printf(file, "%-5s %14.7f %14.7f %14.7f \n", cell.name[i], cell.location[i]...)
		end
		@printf(file, "end atoms_frac\n\n")
	elseif typeof(cell).parameters[1] <: CartesianCoordinates
		@printf(file, "begin atoms_cart\n")
		for i in eachindex(cell.location)
			@printf(file, "%-5s %14.7f %14.7f %14.7f \n", cell.name[i], cell.location[i]...)
		end
		@printf(file, "end atoms_cart\n\n")
	else
		error("Wrong cell.")
	end

	if grid isa MonkhorstPack
		kpoints = RedKgrid(grid).kdirect
	elseif grid isa RedKgrid
		kpoints = grid.kdirect
	else
		error("Wrong grid!")
	end

	@printf(file, "mp_grid = %6u %6u %6u \n", grid.kgrid_size...)
	@printf(file, "begin kpoints\n")
	for k in kpoints
		@printf(file, "%20.12f %20.12f %20.12f \n", k...)
	end
	@printf(file, "end kpoints\n")

	close(file)

	return nothing
end
