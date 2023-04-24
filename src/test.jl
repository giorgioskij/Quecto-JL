using Pkg
using Revise
using Jtrace
# using Cthulhu
# using JET
# using ProfileView

# scenePath = "03_texture/texture.json"
# scenePath = "02_matte/matte.json"
#scenePath = "02_matte/bunny.json"
# scenePath = "04_envlight/envlight.json"
# scenePath = "12_ecosys/ecosys.json"
scenePath = "materials1/materials1.json"
# scenePath = "materials2/materials2.json"
# scenePath = "materials4/materials4.json"
# scenePath = "features2/features2.json"

trace(scenePath = scenePath, shader = "pathtrace", width = 1280, samples = 128)
