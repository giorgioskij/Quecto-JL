include("types.jl")
include("jtrace.jl")
include("loader.jl")

# set global paths
const baseDir = dirname(@__FILE__)
cd(baseDir)
const scenePath = joinpath(baseDir, "02_matte/matte.json")

# set shader to use
const shader = shaderEyelight

# render 
run(300, 200, 2)
