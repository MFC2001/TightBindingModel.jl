
include("./integrate/integrate.jl")
include("ORBITAL2HR.jl")
include("./unstackPOSCAR.jl")

include("./TwoLayerPOSCAR2HR.jl")



function MoirePOSCAR2HR()

	#1：关于ORBITAL中的原子位置信息；
	#2：是用数据结构还是索引数组；
	#3：数据结构后期改动函数更方便，可以写关于数据结构的合并以及orbital2interlayerhr！至于会增加计算量，该函数总的计算量应该很难过大！

end
