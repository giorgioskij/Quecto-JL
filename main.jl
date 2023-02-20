include("types.jl")
include("algebra.jl")
include("bvh.jl")
include("jtrace.jl")
include("loader.jl")

# set global paths
const baseDir = dirname(@__FILE__)
cd(baseDir)

# const scenePath = joinpath(baseDir, "03_texture/texture.json")
const scenePath = joinpath(baseDir, "02_matte/bunny.json")
# const scenePath = joinpath(baseDir, "04_envlight/envlight.json")

# set shader to use
const shader = shaderNormal
run(300, 200, 2)
