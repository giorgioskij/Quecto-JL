using StaticArrays, Images
using BenchmarkTools


const SVec3f = SVector{3,Float32}
const SVec2f = SVector{2,Float32}
# const Frame = SVector{4,Float32}


struct Frame
    x::SVec3f
    y::SVec3f
    z::SVec3f
    o::SVec3f
end

struct Sphere
    center::UInt32 # indexing Scene.points
    radius::Float32
    color::SVec3f # temp

    Sphere(center, radius) = new(center, radius, (0.925, 0.36, 0.38))
end

struct Scene
    points::Vector{SVec3f}
    spheres::Vector{Sphere}
end

struct Ray
    origin::SVec3f
    direction::SVec3f
    tmin::Float32
    tmax::Float32

    Ray(origin, direction, tmin, tmax) = new(origin, direction, tmin, tmax)
    Ray(origin, direction) = new(origin, direction, 0.0, 0.0)

end


struct HitObject
    hit::Bool
    sphereHit::Union{Sphere,Nothing}
    ray::Union{Ray,Nothing}

    HitObject(hit, sphereHit, ray) = new(hit, sphereHit, ray)
    HitObject(hit) = new(hit, nothing, nothing)

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

end



# main entry point to the program
function run(width, height)

    # reads params and initializes stuff

    # generate scene
    scene = generateScene()

    # generate camera
    camera = Camera()

    # generate empty starting image
    image = zeros(SVec3f, height, width)
    # image = zeros(Float32, height, width, 3)

    # call the function to trace samples
    traceSamples(image, scene, camera, width, height)

    # save the resulting image
    rgbImage = zeros(RGB, size(image))
    for i in 1:size(image)[1], j in 1:size(image)[2]
        # rgbImage[i, j] = RGB(image[i, j][1]^0.45, image[i, j][2]^0.45, image[i, j][3]^0.45)
        rgbImage[i, j] = RGB(image[i, j]...)
    end
    save("prova.png", rgbImage)

end

function generateScene()::Scene
    # generate point vector
    points = Vector{SVec3f}([[0, 0, -100]])

    # generate sphere
    spheres = Vector{Sphere}([Sphere(1, 5.0)])

    return Scene(points, spheres)
end

function traceSamples(image, scene, camera, imwidth, imheight)
    # loop over pixels
    # TODO: add threads
    for i in 1:size(image)[2]
        for j in 1:size(image)[1]
            color = traceSample(i, j, scene, camera, imwidth, imheight)
            image[j, i] += color
            # image[i, j, :] += color
        end
    end
end

function traceSample(i::Int,
    j::Int,
    scene::Scene,
    camera::Camera,
    imwidth::Int,
    imheight::Int
)::SVec3f

    # send a ray
    ray = sampleCamera(camera, i, j, imwidth, imheight)

    # call the shader
    radiance = shader(scene, ray, camera)

    return radiance

end

# function intersectEnvironment(ray)::SVec3f
#     # return {1, 1, 1}
#     return SVec3f([1, 1, 1])
# end

function shader(scene::Scene, ray::Ray, camera::Camera)::SVec3f
    background = SVec3f(0.105, 0.443, 0.90)

    hit::HitObject = hitSphere(ray, scene, scene.spheres[1])

    if !hit.hit
        color = background
        return color
    end

    ray::Ray = hit.ray
    sphere::Sphere = hit.sphereHit
    sphereCenter::SVec3f = scene.points[sphere.center]

    # compute coordinates of point hit
    pointHit::SVec3f = ray.origin + ray.tmin * ray.direction

    # compute normal: n = (p-c) / |p-c|
    normal = unitVector(pointHit - sphereCenter)

    radiance = 0.5 .* (normal .+ 1) .* sphere.color

    return radiance
end


function hitSphere(ray::Ray, scene::Scene, sphere::Sphere)::HitObject

    sphereCenter::SVec3f = scene.points[sphere.center]
    oc::SVec3f = ray.origin - sphereCenter
    a::Float32 = dot(ray.direction, ray.direction)
    b::Float32 = 2 * dot(oc, ray.direction)
    c::Float32 = dot(oc, oc) - sphere.radius^2

    delta = b^2 - 4 * a * c

    if delta < 0
        return HitObject(false)
    else
        hitDist::Float32 = (-b - sqrt(delta)) / (2 * a)

        if hitDist::Float32 < ray.tmin
            return HitObject(false)
        end

        new_ray = Ray(ray.origin, ray.direction, hitDist, ray.tmax)
        return HitObject(true, sphere, new_ray)
    end
end


function sampleCamera(camera::Camera, i::Int, j::Int, imwidth::Int, imheight::Int)
    u::Float32 = i / imwidth
    v::Float32 = j / imheight
    ray = evalCamera(camera, u, v)
    return ray
end

# takes a camera, image coordinates (uv) and generates a ray connecting them
function evalCamera(camera::Camera, u::Float32, v::Float32)
    film::SVec2f = (camera.aspect >= 1
                    ? SVec2f([camera.film, camera.film / camera.aspect])
                    : SVec2f([camera.film * camera.aspect, camera.film]))

    q = SVec3f([film[1] * (0.5f0 - u), film[2] * (v - 0.5f0), camera.lens])

    # ray direction through the lens center
    dc::SVec3f = -normalize(q)

    # point on the lens
    pointOnLens = SVec3f([
        u * camera.aperture / 2.0,
        v * camera.aperture / 2.0, 0
    ])

    # point on the focus plane
    pointOnFocusPlane = dc * camera.focus / abs(dc[3])

    # correct ray direction to account for camera focusing
    direction = normalize(pointOnFocusPlane - pointOnLens)


    ray_origin = transformPoint(camera.frame, pointOnLens)
    ray_direction = transformDirection(camera.frame, direction)

    return Ray(ray_origin, ray_direction)
end

function normalize(v::SVec3f)::SVec3f
    l = length(v)
    return (l != 0) ? v / l : v
end

function length(v::SVec3f)::Float32
    return sqrt(dot(v, v))
end

function transformPoint(frame::Frame, v::SVec3f)::SVec3f
    return frame.x * v[1] + frame.y * v[2] + frame.z * v[3] + frame.o
end

function transformVector(frame::Frame, v::SVec3f)::SVec3f
    return frame.x * v[1] + frame.y * v[2] + frame.z * v[3]
end

function transformDirection(frame::Frame, v::SVec3f)::SVec3f
    return normalize(transformVector(frame, v))
end

function unitVector(v::SVec3f)::SVec3f
    return v / length(v)
end


# println(Threads.nthreads())
run(1080, 720)




#TODO: benchmark iteration order
#TODO: try out if image is better srotolata or not
#TODO: check if explicit inlining is possible
#TODO: check argument passing