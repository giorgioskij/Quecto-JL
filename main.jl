include("types.jl")
include("jtrace.jl")
include("loader.jl")


const sceneFile = "julia-pathtracer/02_matte/matte.json"
const shader = shaderNormal
run(300, 200, 2)