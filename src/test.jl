using Pkg
using Revise
using Jtrace
scenePath = "03_texture/texture.json"
# scenePath = "02_matte/matte.json"
#scenePath = "02_matte/bunny.json"
# scenePath = "04_envlight/envlight.json"
# scenePath = "12_ecosys/ecosys.json"
# scenePath = "materials1/materials1.json"

# @time trace(
#     scenePath = scenePath,
#     shader = "eyelight",
#     width = 1280,
#     samples = 10,
#     multithreaded = true,
# )
