module World

using StaticArrays: dot
using ..Algebra
using ..Types
using ..Types: SVec3f, SVec3i
using Images

export Scene,
    Instance,
    Shape,
    Texture,
    Material,
    Environment,
    Camera,
    evalCamera,
    sampleCamera,
    MaterialPoint,
    Subdiv

struct Camera
    frame::Frame
    lens::Float32
    film::Float32
    aspect::Float32
    focus::Float32
    aperture::Float32

    Camera() = new(
        Frame([1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]),
        0.05,
        0.036,
        1.5,
        10000,
        0,
    )
    Camera(frame) = new(frame, 0.05, 0.036, 1.5, 10000, 0)
    Camera(frame, lens, film, aspect, focus, aperture) =
        new(frame, lens, film, aspect, focus, aperture)
end

struct Instance
    frame::Frame
    shapeIndex::Int32
    materialIndex::Int32

    Instance() = new(Frame([1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]), -1, -1)
    Instance(frame, shapeIndex, materialIndex) =
        new(frame, shapeIndex, materialIndex)
end

struct Shape
    # element data
    quads::Vector{SVec4i}
    triangles::Vector{SVec3i}
    lines::Vector{SVec2i}
    points::Vector{Int32}

    # vertex data
    positions::Vector{SVec3f}
    normals::Vector{SVec3f}
    textureCoords::Vector{SVec2f}
    radius::Vector{Float32}
    #colors::Vector{SVec4f}
end

struct Environment
    frame::Frame
    emission::SVec3f
    emissionTex::Int32

    Environment(frame, emission, emissionTex) =
        new(frame, emission, emissionTex)
    Environment() = new(
        Frame(SVec3f(1, 0, 0), SVec3f(0, 1, 0), SVec3f(0, 0, 1), zeroSV3f),
        zeroSV3f,
        -1,
    )
end

struct Texture
    # image::Matrix{RGBA{N0f8}}
    # hdrImage::Matrix{RGB{N0f16}}

    image::Matrix{SVec4f}
    # hdrImage::Matrix{SVec4f}

    # some parameters from yocto
    # linear::Bool    # textures can be stored in linear or non-linear colorspace
    nearest::Bool
    clamp::Bool

    Texture(image, nearest, clamp) = new(image, nearest, clamp)
    Texture() = new(
        # Matrix{RGBA{N0f8}}(undef, 0, 0),
        # Matrix{RGB{N0f16}}(undef, 0, 0),
        Matrix{SVec4f}(undef, 0, 0),
        # Matrix{SVec4f}(undef, 0, 0),
        # false,
        false,
        false,
    )
end

struct Material
    type::String # TODO: create an enum of material types
    emission::SVec3f
    color::SVec3f
    metallic::Float32
    roughness::Float32
    ior::Float32
    trdepth::Float32
    scattering::SVec3f
    scanisotropy::Float32
    opacity::Float32
    emissionTex::Int32
    colorTex::Int32
    roughnessTex::Int32
    scatteringTex::Int32
    normalTex::Int32

    Material() = new(
        "matte", # TODO: create an enum of material types like 0->matte, 1->metal, 2->glass, etc.
        [0, 0, 0],
        [0, 0, 0],
        0,
        0,
        1.5f0,
        0.01f0,
        [0, 0, 0],
        0,
        1.0f0,
        -1,
        -1,
        -1,
        -1,
        -1,
    )
    Material(
        type,
        emission,
        color,
        metallic,
        roughness,
        ior,
        trdepth,
        scattering,
        scanisotropy,
        opacity,
        emissionTex,
        colorTex,
        roughnessTex,
        scatteringTex,
        normalTex,
    ) = new(
        type,
        emission,
        color,
        metallic,
        roughness,
        ior,
        trdepth,
        scattering,
        scanisotropy,
        opacity,
        emissionTex,
        colorTex,
        roughnessTex,
        scatteringTex,
        normalTex,
    )
end

struct Scene
    cameras::Vector{Camera}
    instances::Vector{Instance}
    shapes::Vector{Shape}
    textures::Vector{Texture}
    materials::Vector{Material}
    environments::Vector{Environment}
end

@fastmath function sampleCamera(
    camera::Camera,
    i::Int,
    j::Int,
    imwidth::Int,
    imheight::Int,
)::Ray
    f1 = rand(Float32)
    f2 = rand(Float32)

    u::Float32 = (i + f1) / imwidth
    v::Float32 = (j + f2) / imheight
    #u = clamp(u, 0, 1)
    #v = clamp(v, 0, 1)

    ray = evalCamera(camera, u, v)
    return ray
end

# takes a camera, image coordinates (uv) and generates a ray connecting them
@fastmath function evalCamera(camera::Camera, u::Float32, v::Float32)::Ray
    film::SVec2f =
        camera.aspect >= 1 ? SVec2f(camera.film, camera.film / camera.aspect) :
        SVec2f(camera.film * camera.aspect, camera.film)

    q = SVec3f(film[1] * (0.5f0 - u), film[2] * (v - 0.5f0), camera.lens)

    # ray direction through the lens center
    dc::SVec3f = -Algebra.norm(q)

    # generate random lens_uv
    uLens = rand(Float32)
    vLens = rand(Float32)

    uLens, vLens = sampleDisk(uLens, vLens)

    # point on the lens
    pointOnLens = SVec3f(
        uLens * camera.aperture / 2.0f0,
        vLens * camera.aperture / 2.0f0,
        0,
    )

    # point on the focus plane
    pointOnFocusPlane = dc * camera.focus / abs(dc[3])

    # correct ray direction to account for camera focusing
    direction = Algebra.norm(pointOnFocusPlane - pointOnLens)

    rayOrigin = transformPoint(camera.frame, pointOnLens)
    rayDirection = transformDirection(camera.frame, direction)

    return Ray(rayOrigin, rayDirection)
end

struct MaterialPoint
    type::String
    emission::SVec3f
    color::SVec3f
    opacity::Float32
    roughness::Float32
    metallic::Float32
    ior::Float32
    density::SVec3f
    scattering::SVec3f
    scanisotropy::Float32
    trdepth::Float32

    MaterialPoint(
        type = "matte",
        emission = zeroSV3f,
        color = zeroSV3f,
        opacity = 1,
        roughness = 0,
        metallic = 0,
        ior = 1,
        density = zeroSV3f,
        scattering = zeroSV3f,
        scanisotropy = 0,
        trdepth = 0.01f0,
    ) = new(
        type,
        emission,
        color,
        opacity,
        roughness,
        metallic,
        ior,
        density,
        scattering,
        scanisotropy,
        trdepth,
    )
end

struct Subdiv
    # face-varying primitives
    quadsPos::Vector{SVec4i}
    quadsNorm::Vector{SVec4i}
    quadsTexCoord::Vector{SVec4i}

    # vertex data
    positions::Vector{SVec3f}
    normals::Vector{SVec3f}
    textureCoords::Vector{SVec2f}

    # subdivision data
    subdivisions::Int32
    catmullclark::Bool
    smooth::Bool

    # displacement data
    displacement::Float32
    displacementTex::Int32

    # shape reference
    shapeIndex::Int32

    Subdiv(;
        quadsPos = Vector{SVec4i}[],
        quadsNorm = Vector{SVec4i}[],
        quadsTexCoord = Vector{SVec4i}[],
        positions = Vector{SVec3f}[],
        normals = Vector{SVec3f}[],
        textureCoords = Vector{SVec2f}[],
        subdivisions = 0,
        catmullclark = true,
        smooth = true,
        displacement = 0,
        displacementTex = -1,
        shapeIndex = -1,
    ) = new(
        quadsPos,
        quadsNorm,
        quadsTexCoord,
        positions,
        normals,
        textureCoords,
        subdivisions,
        catmullclark,
        smooth,
        displacement,
        displacementTex,
        shapeIndex,
    )
end

end