module Types

export SVec3f, SVec2f, Frame, Camera, Instance, Shape, Scene, Ray, HitObject, Triangle


using StaticArrays

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
        0
    )
    Camera(frame) = new(frame, 0.05, 0.036, 1.5, 10000, 0)
    Camera(frame, lens, film, aspect, focus, aperture) = new(frame, lens, film, aspect, focus, aperture)
end

struct Shape
    # element data
    triangles::Matrix{Int32}
    # vertex data
    positions::Matrix{Float32}
    normals::Matrix{Float32}
    textureCoords::Matrix{Float32}
end

struct Instance
    frame::Frame
    shape::Int64

    Instance() = new(Frame([1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]), -1)
    Instance(frame, shape) = new(frame, shape)
end


struct Scene
    cameras::Vector{Camera}
    instances::Vector{Instance}
    shapes::Vector{Shape}

end

struct Ray
    origin::SVec3f
    direction::SVec3f
    tmin::Float32
    tmax::Float32

    Ray(origin, direction, tmin, tmax) = new(origin, direction, tmin, tmax)
    Ray(origin, direction) = new(origin, direction, 0.0001f0, typemax(Float32))
end

struct HitObject
    hit::Bool
    u::Float32
    v::Float32
    ray::Union{Ray,Nothing}
    HitObject(hit, u, v, ray) = new(hit, u, v, ray)
    HitObject(hit) = new(hit, -1, -1, nothing)
end

struct Triangle
    x::SVec3f
    y::SVec3f
    z::SVec3f
end

end