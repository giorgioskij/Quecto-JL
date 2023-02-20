include("types.jl")
include("algebra.jl")
include("bvh.jl")
include("jtrace.jl")
include("loader.jl")

# set global paths
const baseDir = dirname(@__FILE__)
cd(baseDir)
const scenePath = joinpath(baseDir, "03_texture/texture.json")

# set shader to use
const shader = shaderEyelight
# run(30, 20, 2)
