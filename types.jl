module Types
using StaticArrays

export SVec4f,
    SVec3f,
    SVec2f,
    Frame,
    Camera,
    Instance,
    Shape,
    Texture,
    Scene,
    Ray,
    Intersection,
    Triangle,
    Quad

const SVec4f = SVector{4,Float32}
const SVec3f = SVector{3,Float32}
const SVec2f = SVector{2,Float32}

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

struct Texture
    width::Int64
    hegiht::Int64

    # some parameters from yocto
    linear::Bool    # textures can be stored in linear or non-linear colorspace
    nearest::Bool
    clamp::Bool

    # a texture can be an array of float or byte pixels
    pixelsf::Vector{SVec4f}
    pixelsb::Vector{SVec4f}
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