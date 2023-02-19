include("types.jl")
include("algebra.jl")
include("bvh.jl")
include("jtrace.jl")
include("loader.jl")

# set global paths
const baseDir = dirname(@__FILE__)
cd(baseDir)
const scenePath = joinpath(baseDir, "02_matte/matte.json")

# set shader to use
const shader = shaderEyelight
bvh, scene = run(30, 20, 2)
