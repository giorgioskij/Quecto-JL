include("types.jl")
include("jtrace.jl")
include("loader.jl")

# set current working directory to main project dir
const baseDir = dirname(@__FILE__)
cd(baseDir)
const scenePath = joinpath(baseDir, "02_matte/matte.json")

const sceneFile = "julia-pathtracer/02_matte/matte.json"
const shader = shaderEyelight
run(64, 64, 2)
