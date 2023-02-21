function randoms()
    c = 0
    for i = 1:10_000_000
        x, y = rand(Float32, 2)

        c = x + y
    end
    return c
end

using RandomNumbers

function randoms2()
    c = 0
    for i = 1:10_000_000
        x, y = rand(rng_xor, Float32, 2)

        c = x + y
    end
    return c
end

function test(rng)
    c = 0
    for i = 1:10_000_000
        # x, y = rand(rng, Float32, 2)
        x = rand(rng, Float32)
        y = rand(rng, Float32)

        c = x + y
    end
    return c
end

# mutable struct RngState
#     state::UInt64
#     inc::UInt64

#     RngState() = new(0x853c49e6748fea9b, 0xda3e39cb94b95bdb)
# end

# function advanceRng(rng::RngState)::UInt32
#     oldState::UInt64 = rng.state
#     rng.state = oldState * 0x5851f42d4c957f2d + rng.inc
#     xorShifted = UInt32(((oldState >> 0x12) ⊻ oldState) >> 0x1b)
#     rot = UInt32(oldstate >> 0x3b)

#     return (xorShifted >> rot) | (xorShifted << ((~rot + 0x1) & 0x1f))
# end