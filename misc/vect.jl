using StaticArrays
using LoopVectorization

function test(a)
    s = 0
    for i = 1:1000000
        s += fastMinimum(a)
        a = a .+ 0.1
    end
    return s
end
function vec(a)
    s = 0
    for i = 1:1000000
        s += vecMinimum(a)
        a = a .+ 0.1
    end
    return s
end

@inbounds function fastMinimum(a)
    ifelse(
        a[1] < a[2],
        ifelse(a[1] < a[3], a[1], a[3]),
        ifelse(a[2] < a[3], a[2], a[3]),
    )
end

function vecMinimum(a)
    min = a[1]
    for i = 2:3
        if a[i] < min
            min = a[i]
        end
    end
    return min
end