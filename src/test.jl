using Pkg
#Pkg.instantiate()
using Jtrace
using Revise
# using Cthulhu
# using JET
# using ProfileView

texture = "scenes/03_texture/texture.json"
matte = "scenes/02_matte/matte.json"
bunny = "scenes/02_matte/bunny.json"
envlight = "scenes/04_envlight/envlight.json"
ecosys = "scenes/12_ecosys/ecosys.json"
materials1 = "scenes/materials1/materials1.json"
materials2 = "scenes/materials2/materials2.json"
materials4 = "scenes/materials4/materials4.json"
scenePath = "scenes/features2/features2.json"
coffee = "scenes/coffee/coffee.json"
bathroom1 = "scenes/bathroom1/bathroom1.json"
livingroom1 = "scenes/livingroom1/livingroom1.json"
livingroom2 = "scenes/livingroom2/livingroom2.json"
livingroom3 = "scenes/livingroom3/livingroom3.json"
staircase1 = "scenes/staircase1/staircase1.json"
staircase2 = "scenes/staircase2/staircase2.json"
kitchen = "scenes/kitchen/kitchen.json"
classroom = "scenes/classroom/classroom.json"
bathroom2 = "scenes/bathroom2/bathroom2.json"

# trace(
#     scenePath = coffee,
#     shader = "volumetric",
#     resolution = 1280,
#     samples = 128,
#     filename = "coffee_1280_128_v.png",
# )

# trace(
#     scenePath = staircase1,
#     shader = "volumetric",
#     resolution = 1280,
#     samples = 4096,
#     filename = "staircase1_1280_4096_v.png",
# )

# trace(
#     scenePath = staircase2,
#     shader = "volumetric",
#     resolution = 1280,
#     samples = 4096,
#     filename = "staircase2_1280_4096_v.png",
# )
# trace(
#     scenePath = kitchen,
#     shader = "volumetric",
#     resolution = 1280,
#     samples = 4096,
#     filename = "kitchen_1280_4096_v.png",
# )
# trace(
#     scenePath = classroom,
#     shader = "volumetric",
#     resolution = 1280,
#     samples = 4096,
#     filename = "classroom_1280_4096_v.png",
# )
# trace(
#     scenePath = bathroom2,
#     shader = "volumetric",
#     resolution = 1280,
#     samples = 4096,
#     filename = "bathroom2_1280_4096_v.png",
# )

trace(
    scenePath = bathroom1,
    shader = "volumetric",
    resolution = 1280,
    samples = 128,
    filename = "bathroom1_1280_128_v.png",
)

# trace(
#     scenePath = livingroom2,
#     shader = "volumetric",
#     resolution = 1280,
#     samples = 8192,
#     filename = "livingroom2_1280_8192_v.png",
# )

# trace(
#     scenePath = livingroom3,
#     shader = "volumetric",
#     resolution = 1280,
#     samples = 8192,
#     filename = "livingroom3_1280_8192_v.png",
# )