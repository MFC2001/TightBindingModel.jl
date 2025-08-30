struct Wannier90chk
	header::String
	num_bands::Int
	exclude_bands::Vector{Int}
	lattice::Lattice{Float64}
	rlattice::ReciprocalLattice{Float64}
	mp_grid::Vector{Int}
	kdirect::Vector{Vec3{Float64}}
	nntot::Int
	num_wann::Int
	checkpoint::String
	have_disentangled::Bool
	omega_invariant::Float64
	lwindow::Matrix{Bool}
	ndimwin::Vector{Int}
	u_matrix_opt::Array{ComplexF64, 3}
	u_matrix::Array{ComplexF64, 3}
	m_matrix::Array{ComplexF64, 4}
	wannier_centres::Vector{Vec3{Float64}}
	wannier_spreads::Vector{Float64}
end
function Readchk(file::AbstractString)
	file = open(file, "r")

	# 1. 读取头部信息 (33字符)
	header = read_record(file, String)

	# 2. 读取基本参数
	num_bands = Int(read_record(file, Int32))
	num_exclude_bands = Int(read_record(file, Int32))

	# 3. 读取排除能带列表
	exclude_bands = read_record(file, Int32, (num_exclude_bands))

	# 4. 读取晶格信息
	lattice = Lattice(read_record(file, Float64, (3, 3))')
	rlattice = ReciprocalLattice(read_record(file, Float64, (3, 3))')

	# 5. 读取k点信息
	num_kpts = Int(read_record(file, Int32))
	mp_grid = Int.(read_record(file, Int32, (3)))
	kdirect = read_record(file, Float64, (3, num_kpts))
	kdirect = map(Vec3{Float64}, eachcol(kdirect))

	# 6. 读取其他参数
	nntot = Int(read_record(file, Int32))
	num_wann = Int(read_record(file, Int32))
	checkpoint = read_record(file, String)
	have_disentangled = read_record(file, Bool)

	# 7. 处理disentangled相关数据
	if have_disentangled
		omega_invariant = read_record(file, Float64)
		lwindow = read_record(file, Bool, (num_bands, num_kpts))
		ndimwin = Int.(read_record(file, Int32, (num_kpts)))
		u_matrix_opt = read_record(file, ComplexF64, (num_bands, num_wann, num_kpts))
	else
		omega_invariant = 0
		lwindow = Matrix{Bool}(undef, 0, 0)
		ndimwin = Vector{Int}(undef, 0)
		u_matrix_opt = Array{ComplexF64, 3}(undef, 0, 0, 0)
	end

	# 8. 读取核心U矩阵
	u_matrix = read_record(file, ComplexF64, (num_wann, num_wann, num_kpts))

	# 9. 读取重叠矩阵m_matrix
	m_matrix = read_record(file, ComplexF64, (num_wann, num_wann, nntot, num_kpts))

	# 10. 读取Wannier中心
	wannier_centres = read_record(file, Float64, (3, num_wann))
	wannier_centres = map(Vec3{Float64}, eachcol(wannier_centres))

	wannier_spreads = read_record(file, Float64, (num_wann))

	close(file)

	# 返回结构化数据
	return Wannier90chk(header, num_bands, exclude_bands,
		lattice, rlattice,
		mp_grid, kdirect, nntot, num_wann,
		checkpoint,
		have_disentangled,
		omega_invariant, lwindow, ndimwin, u_matrix_opt,
		u_matrix, m_matrix,
		wannier_centres, wannier_spreads,
	)
end

function read_record(io::IOStream, ::Type{Bool}, dims)
	data = read_record(io)
	data = reinterpret(Int32, data)
	if length(data) ≠ prod(dims)
		error("record length mismatch: expected $(prod(dims)), got $(length(data))")
	end
	return reshape(data, dims) .≠ 0
end
function read_record(io::IOStream, T::Type, dims)
	data = read_record(io)
	data = reinterpret(T, data)
	if length(data) ≠ prod(dims)
		error("record length mismatch: expected $(prod(dims)), got $(length(data))")
	end
	return reshape(data, dims)
end
function read_record(io::IOStream, ::Type{String})::String
	data = read_record(io)
	return String(data)
end
function read_record(io::IOStream, ::Type{Bool})::Bool
	data = read_record(io)
	len_data = length(data)
	len_T = sizeof(Int32)
	if len_data ≠ len_T
		error("record type mismatch: expected $len_T bytes of Int32 for Bool, got $len_data bytes.")
	end
	data = reinterpret(Int32, data)[1]
	return data ≠ 0
end
function read_record(io::IOStream, T::Type)::T
	data = read_record(io)
	len_data = length(data)
	len_T = sizeof(T)
	if len_data ≠ len_T
		error("record type mismatch: expected $len_T bytes of $T, got $len_data bytes.")
	end
	data = reinterpret(T, data)[1]
	return data
end
function read_record(io::IOStream, nb::Integer; all = true)::Vector{UInt8}
	data = read_record(io)
	len_data = length(data)
	if len_data ≠ nb
		error("record length mismatch: expected $nb, got $len_data")
	end
	return data
end
function read_record(io::IOStream)::Vector{UInt8}
	nb_start = read(io, Int32)
	data = read(io, nb_start)
	nb_end = read(io, Int32)
	if nb_start ≠ nb_end
		error("record length mismatch: $nb_start vs $nb_end")
	end
	return data
end
