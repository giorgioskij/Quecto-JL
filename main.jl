include("types.jl")
include("jtrace.jl")
include("loader.jl")

# set current working directory to main project dir
const baseDir = dirname(@__FILE__)
cd(baseDir)
const scenePath = joinpath(baseDir, "02_matte/matte.json")

# set the shader to use
const shader = shaderNormal

run(300, 200, 2)
