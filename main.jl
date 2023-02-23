include("types.jl")
include("algebra.jl")
include("bvh.jl")
include("jtrace.jl")
include("loader.jl")

# set global paths
const baseDir = dirname(@__FILE__)
cd(baseDir)

#scenePath = "03_texture/texture.json"
# scenePath = "02_matte/matte.json"
#scenePath = "02_matte/bunny.json"
#scenePath = "04_envlight/envlight.json"
scenePath = "12_ecosys/ecosys.json"

# set shader to use
shader = shaderEyelight

run(scenePath, shader, 1920, 2)
