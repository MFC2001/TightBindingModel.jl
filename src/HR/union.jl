# Base.union(hr,hrs...)


# function union end

# union(s, sets...) = union!(emptymutable(s, promote_eltype(s, sets...)), s, sets...)
# union(s::AbstractSet) = copy(s)

# const ∪ = union

# function union!(s::AbstractSet, sets...)
#     for x in sets
#         union!(s, x)
#     end
#     return s
# end

# max_values(::Type) = typemax(Int)
# max_values(T::Union{map(X -> Type{X}, BitIntegerSmall_types)...}) = 1 << (8*sizeof(T))
# # saturated addition to prevent overflow with typemax(Int)
# function max_values(T::Union)
#     a = max_values(T.a)::Int
#     b = max_values(T.b)::Int
#     return max(a, b, a + b)
# end
# max_values(::Type{Bool}) = 2
# max_values(::Type{Nothing}) = 1

# function union!(s::AbstractSet{T}, itr) where T
#     haslength(itr) && sizehint!(s, length(s) + Int(length(itr))::Int)
#     for x in itr
#         push!(s, x)
#         length(s) == max_values(T) && break
#     end
#     return s
# end