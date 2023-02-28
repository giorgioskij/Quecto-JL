
# const Vec3 = Vector{Float64}

function test(v1, v2)
    # v3 = Vec3(rand(), rand(), rand())
    v3 = rand(Float64, 3)
    for i = 1:10^7
        v3 = (v1 * v2) .+ v3
    end
    return v3
end

# v1 = Vec3(rand(), rand(), rand())
v1 = rand(Float64, 3)
# v2 = Vec3(rand(), rand(), rand())
v2 = rand(Float64, 3)

test(v1, v2)
@time v3 = test(v1, v2)

println(v3)