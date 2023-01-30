using JSON
using("jtrace.jl")
using PlyIO


struct Scene
    cameras::Vector{Camera}
    instances::Vector{Instance}
    shapes::Vector{Shape}
end

struct Instance
    frame::Frame
    shape::Int64

    Instance() = new( Frame([1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]), -1)
    Instance(frame, shape) = new(frame, shape)
end

struct Shape
    # element data

    # vertex data
    positions::Vector{SVec3f}
end

function loadJsonScene(filename)
    json = JSON.parsefile(filename)


    # if haskey(json, "asset")
    #     element = json["asset"]
    #     copyright = element["copyright"]
    # end

    # filenames
    shapeFilenames = Vector{String}(undef, 0)

    # CAMERAS
    if haskey(json, "cameras")
        group = json["cameras"]
        cameras = Vector{Camera}(undef, length(group))   
        # ignore camera_names
        defaultCamera = Camera()
        for (i, element) in enumerate(group)
            camera = Camera(
                get(element, "film", defaultCamera.film),
                get(element, "frame", defaultCamera.frame),
                get(element, "lens", defaultCamera.lens),
                get(element, "film", defaultCamera.film),
                get(element, "aspect", defaultCamera.aspect),
                get(element, "focus", defaultCamera.focus),
                get(element, "aperture", defaultCamera.aperture)
            )
            # ignore params "ortographics", "name", "lookat"
            cameras[i] = camera
        end

    # SHAPES
    if haskey(json, "shapes")
        group = json["shapes"]
        sizehint!(shapeFilenames, length(group))
        for (i, element) in enumerate(group)
            if !haskey(element, "uri")
                throw(MissingException("uri not present"))
            end
            shapeFilenames[i] = element["uri"]
        end
    end

    # INSTANCES
    if haskey(json, "instances")
        group = json["instances"]
        instances = Vector{Instance}(undef, legnth(group))

        defaultInstance =Instance()
        for (i, element) in enumerate(group)
            instance = Instance(
                get(element, "frame", defaultInstance.frame),
                get(element, "shape", defaultInstance.shape)
            )
            # ignore params "name", "material"
            instances[i] = instance
        end
    end

    # LOAD RESOURCES

    # load shapes
    shapes = Vector{Shape}(undef, length(shapeFilenames))

    for (i, filename) in enumerate(shapeFilenames)
        shape = loadShape(filename)
        shapes[i] = shape
    end

    # ignore load subdivs, textures


    
      

    # TODO: generate the scene from json        
    scene = Scene(cameras, instances, shapes)

    return scene
end


function loadShape(filename::String)
    ply = load_ply(filename)
    
end

loadJsonScene("julia-pathtracer/02_matte/matte.json")