
function _KrylovKit_init()

	# 检查是否使用KrylovKit
	if isdefined(Main, :KrylovKit)
		include("./config.jl")
        include("./eigsolve.jl")
		@info "Loaded KrylovKit support."
	else
		@info "No KrylovKit support."
	end

end
_KrylovKit_init()
