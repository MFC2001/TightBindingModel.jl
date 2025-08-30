
δ(i, j) = i == j ? 1 : 0

function kronecker_product(pairs::Vararg{Tuple{Int, Int}})
    for (i, j) in pairs
        i == j || return 0
    end
    return 1
end

macro kronecker_shortcircuit(pairs...)
    expr = 1
    for pair in reverse(pairs)
        (i, j) = pair.args
        expr = :(($i == $j) ? $expr : 0)
    end
    return expr
end