module World

using StaticArrays: dot
using ..Algebra
using ..Types
using Images

export Scene,
    Instance,
    Shape,
    Texture,
    Material,
    Environment,
    Camera,
    evalCamera,
    sampleCamera

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
    shapeIndex::Int64
    materialIndex::Int64

    Instance() = new(Frame([1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]), -1, -1)
    Instance(frame, shapeIndex, materialIndex) =
        new(frame, shapeIndex, materialIndex)
end

struct Shape
    # element data
    triangles::Matrix{Int64}
    quads::Matrix{Int64}

    # vertex data
    positions::Matrix{Float32}
    normals::Matrix{Float32}
    textureCoords::Matrix{Float32}
end

struct Environment
    frame::Frame
    emission::SVec3f
    emissionTex::Int

    Environment(frame, emission, emissionTex) =
        new(frame, emission, emissionTex)
    Environment() = new(
        Frame(
            SVec3f(1, 0, 0),
            SVec3f(0, 1, 0),
            SVec3f(0, 0, 1),
            SVec3f(0, 0, 0),
        ),
        SVec3f(0, 0, 0),
        -1,
    )
end

struct Texture
    image::Matrix{RGBA{N0f8}}
    hdrImage::Matrix{RGB{N0f16}}

    # some parameters from yocto
    linear::Bool    # textures can be stored in linear or non-linear colorspace
    nearest::Bool
    clamp::Bool

    Texture(image, hdrImage, linear, nearest, clamp) =
        new(image, hdrImage, linear, nearest, clamp)
    Texture() = new(
        Matrix{RGBA{N0f8}}(undef, 0, 0),
        Matrix{RGB{N0f16}}(undef, 0, 0),
        false,
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
    emissionTex::Int64
    colorTex::Int64
    roughnessTex::Int64
    scatteringTex::Int64
    normalTex::Int64

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

function sampleCamera(
    camera::Camera,
    i::Int64,
    j::Int64,
    imwidth::Int64,
    imheight::Int64,
)
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
function evalCamera(camera::Camera, u::Float32, v::Float32)
    film::SVec2f = (
        camera.aspect >= 1 ?
        SVec2f(camera.film, camera.film / camera.aspect) :
        SVec2f(camera.film * camera.aspect, camera.film)
    )

    q = SVec3f(film[1] * (0.5f0 - u), film[2] * (v - 0.5f0), camera.lens)

    # ray direction through the lens center
    dc::SVec3f = -Algebra.norm(q)

    # generate random lens_uv
    uLens = rand(Float32)
    vLens = rand(Float32)

    uLens, vLens = sampleDisk(uLens, vLens)

    # point on the lens
    pointOnLens =
        SVec3f(uLens * camera.aperture / 2.0, vLens * camera.aperture / 2.0, 0)

    # point on the focus plane
    pointOnFocusPlane = dc * camera.focus / abs(dc[3])

    # correct ray direction to account for camera focusing
    direction = Algebra.norm(pointOnFocusPlane - pointOnLens)

    ray_origin = transformPoint(camera.frame, pointOnLens)
    ray_direction = transformDirection(camera.frame, direction)

    return Ray(ray_origin, ray_direction)
end

end