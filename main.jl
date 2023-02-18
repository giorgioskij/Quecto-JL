include("types.jl")
include("jtrace.jl")
include("loader.jl")


const sceneFile = "julia-pathtracer/02_matte/matte.json"
const shader = shaderColor
run(256, 256, 2)