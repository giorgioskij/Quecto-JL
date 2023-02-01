using JSON
using PlyIO
using .Types

function loadJsonScene(filename)
    json = JSON.parsefile(filename)


    # if haskey(json, "asset")
    #     element = json["asset"]
    #     copyright = element["copyright"]
    # end

    # CAMERAS
    if haskey(json, "cameras")
        group = json["cameras"]
        cameras = Vector{Camera}(undef, Base.length(group))
        # ignore camera_names
        defaultCamera = Camera()
        for (i, element) in enumerate(group)

            f = get(element, "frame", [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
            frame = Frame(f[1:3], f[4:6], f[7:9], f[10:12])

            camera = Camera(
                frame,
                # get(element, "frame", defaultCamera.frame),
                get(element, "lens", defaultCamera.lens),
                get(element, "film", defaultCamera.film),
                get(element, "aspect", defaultCamera.aspect),
                get(element, "focus", defaultCamera.focus),
                get(element, "aperture", defaultCamera.aperture)
            )
            # ignore params "ortographics", "name", "lookat"
            cameras[i] = camera
        end
    end

    # SHAPES
    # shapeFilenames::Vector{string}
    # filenames
    shapeFilenames = Vector{String}(undef, 0)
    if haskey(json, "shapes")
        group = json["shapes"]
        # sizehint!(shapeFilenames, length(group))
        shapeFilenames = Vector{String}(undef, Base.length(group))
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
        instances = Vector{Instance}(undef, Base.length(group))

        defaultInstance = Instance()
        for (i, element) in enumerate(group)

            f = get(element, "frame", [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
            frame = Frame(f[1:3], f[4:6], f[7:9], f[10:12])

            instance = Instance(
                frame,
                #get(element, "frame", defaultInstance.frame),
                get(element, "shape", defaultInstance.shape)
            )
            # ignore params "name", "material"
            instances[i] = instance
        end
    end

    # LOAD RESOURCES

    # load shapes
    shapes = Vector{Shape}(undef, Base.length(shapeFilenames))

    for (i, filename) in enumerate(shapeFilenames)
        fullFilename = "julia-pathtracer/02_matte/"
        shape = loadShape(fullFilename * filename)
        shapes[i] = shape
    end

    # ignore load subdivs, textures

    scene = Scene(cameras, instances, shapes)

    # println(scene)

    return scene
end


function loadShape(filename::String)
    ply = load_ply(filename)

    # check properties
    vertexProperties = plyname.(ply["vertex"].properties)

    # load positions
    if "x" in vertexProperties &&
       "y" in vertexProperties &&
       "z" in vertexProperties
        x = ply["vertex"]["x"]
        y = ply["vertex"]["y"]
        z = ply["vertex"]["z"]
        positions = hcat(x, y, z)
    else
        positions = Matrix{Float32}(undef, 0, 0)
    end

    # load normals
    if "nx" in vertexProperties &&
       "ny" in vertexProperties &&
       "nz" in vertexProperties
        nx = ply["vertex"]["nx"]
        ny = ply["vertex"]["ny"]
        nz = ply["vertex"]["nz"]
        normals = hcat(nx, ny, nz)
    else
        normals = Matrix{Float32}(undef, 0, 0)
    end

    # load texture coordinates
    if "u" in vertexProperties && "v" in vertexProperties
        u = ply["vertex"]["u"]
        v = ply["vertex"]["v"]
        textureCoords = hcat(u, v)
    else
        textureCoords = Matrix{Float32}(undef, 0, 0)
    end

    if size(ply["face"]["vertex_indices"][1]) != (3,)
        throw(MissingException("Only implemented triangles right now"))
    end
    # load triangles
    faces = reduce(hcat, (Vector(ply["face"]["vertex_indices"])))'

    # create shape
    return Shape(faces, positions, normals, textureCoords)
end

#scene = loadJsonScene("02_matte/bunny.json")