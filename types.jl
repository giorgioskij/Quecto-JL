module Types
using StaticArrays
using Images

export SVec4f,
    SVec3f,
    SVec2f,
    SVec3i,
    Frame,
    Camera,
    Instance,
    Shape,
    Texture,
    Material,
    Environment,
    Scene,
    Ray,
    Intersection,
    Triangle,
    Quad

const SVec4f = SVector{4,Float32}
const SVec3f = SVector{3,Float32}
const SVec2f = SVector{2,Float32}
const SVec3i = SVector{3,Int}

struct Frame
    x::SVec3f
    y::SVec3f
    z::SVec3f
    o::SVec3f
end

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
    # for now a texture is a matrix of bytes
    image::Matrix{RGBA{N0f8}}

    # some parameters from yocto
    linear::Bool    # textures can be stored in linear or non-linear colorspace
    nearest::Bool
    clamp::Bool

    Texture(image, linear, nearest, clamp) = new(image, linear, nearest, clamp)
    Texture() = new(Matrix{RGBA{N0f8}}(undef, 0, 0), false, false, false)
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

struct Instance
    frame::Frame
    shapeIndex::Int64

    Instance() = new(Frame([1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]), -1)
    Instance(frame, shapeIndex) = new(frame, shapeIndex)
end

struct Scene
    cameras::Vector{Camera}
    instances::Vector{Instance}
    shapes::Vector{Shape}
    textures::Vector{Texture}
    materials::Vector{Material}
    environments::Vector{Environment}
end

struct Ray
    origin::SVec3f
    direction::SVec3f
    tmin::Float32
    tmax::Float32

    Ray(origin, direction, tmin, tmax) = new(origin, direction, tmin, tmax)
    Ray(origin, direction) = new(origin, direction, 0.0001f0, typemax(Float32))
end

struct Intersection
    hit::Bool
    instanceIndex::Int64
    elementIndex::Int64
    u::Float32
    v::Float32
    distance::Float32

    isTriangle::Bool

    Intersection(hit, instanceIndex, elementIndex, u, v, distance, isTriangle) =
        new(hit, instanceIndex, elementIndex, u, v, distance, isTriangle)
    Intersection(hit) = new(hit, -1, -1, 0, 0, 0, true)
end

struct Triangle
    x::SVec3f
    y::SVec3f
    z::SVec3f

    Triangle(x, y, z) = new(x, y, z)
end

struct Quad
    a::SVec3f
    b::SVec3f
    c::SVec3f
    d::SVec3f

    Quad(a, b, c, d) = new(a, b, c, d)
end

# end of module 
end