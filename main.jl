include("types.jl")
include("algebra.jl")
include("bvh.jl")
include("jtrace.jl")
include("loader.jl")

# set global paths
const baseDir = dirname(@__FILE__)
cd(baseDir)

# const scenePath = joinpath(baseDir, "03_texture/texture.json")
const scenePath = joinpath(baseDir, "02_matte/matte.json")
# const scenePath = joinpath(baseDir, "02_matte/bunny.json")
# const scenePath = joinpath(baseDir, "04_envlight/envlight.json")
# const scenePath = joinpath(baseDir, "12_ecosys/ecosys.json")

# set shader to use
shader = shaderNormal
# const shader = shaderEyelight
run(shader, 1920, 1080, 2)
