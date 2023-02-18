include("types.jl")
include("jtrace.jl")
include("loader.jl")


const sceneFile = "julia-pathtracer/02_matte/matte.json"
const shader = shaderEyelight
run(480, 320, 2)
