using JSON
using PlyIO
using .Types
using Images

function loadJsonScene(scenePath::String)::Scene
    json = JSON.parsefile(scenePath)

    # if haskey(json, "asset")
    #     element = json["asset"]
    #     copyright = element["copyright"]
    # end

    # CAMERAS
    if haskey(json, "cameras")
        group = json["cameras"]
        cameras = Vector{Camera}(undef, size(group, 1))
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
                get(element, "aperture", defaultCamera.aperture),
            )
            # ignore params "ortographics", "name", "lookat"
            cameras[i] = camera
        end
    end

    # TEXTURES
    textureFilenames = Vector{String}(undef, 0)
    textures = Vector{Texture}(undef, 0)
    if haskey(json, "textures")
        group = json["textures"]
        textureFilenames = Vector{String}(undef, size(group, 1))
        textures = Vector{Texture}(undef, size(group, 1))
        defaultTexture = Texture()
        for (i, element) in enumerate(group)
            if !haskey(element, "uri")
                throw(MissingException("uri not present"))
            end
            textureFilenames[i] = element["uri"]
            textures[i] = Texture(
                defaultTexture.image,
                # defaultTexture.hdrImage,
                # get(element, "linear", defaultTexture.linear),
                get(element, "nearest", defaultTexture.nearest),
                get(element, "clamp", defaultTexture.clamp),
            )
        end
    end

    # Materials
    materials = Vector{Material}(undef, 0)
    if haskey(json, "materials")
        group = json["materials"]
        materials = Vector{Material}(undef, size(group, 1))
        # ignore material_names
        defaultMaterial = Material()
        for (i, element) in enumerate(group)
            colorTex = get(element, "color_tex", defaultMaterial.colorTex)
            if colorTex != -1
                colorTex += 1
            end

            emissionTex =
                get(element, "emission_tex", defaultMaterial.emissionTex)
            if emissionTex != -1
                emissionTex += 1
            end

            roughnessTex =
                get(element, "roughness_tex", defaultMaterial.roughnessTex)
            if roughnessTex != -1
                roughnessTex += 1
            end

            scatteringTex =
                get(element, "scatterin_tex", defaultMaterial.scatteringTex)
            if scatteringTex != -1
                scatteringTex += 1
            end

            normalTex = get(element, "normal_tex", defaultMaterial.normalTex)
            if normalTex != -1
                normalTex += 1
            end

            material = Material(
                get(element, "type", defaultMaterial.type),
                get(element, "emission", defaultMaterial.emission),
                get(element, "color", defaultMaterial.color),
                get(element, "metallic", defaultMaterial.metallic),
                get(element, "roughness", defaultMaterial.roughness),
                get(element, "ior", defaultMaterial.ior),
                get(element, "trdepth", defaultMaterial.trdepth),
                get(element, "scattering", defaultMaterial.scattering),
                get(element, "scanisotropy", defaultMaterial.scanisotropy),
                get(element, "opacity", defaultMaterial.opacity),
                emissionTex,
                colorTex,
                roughnessTex,
                scatteringTex,
                normalTex,
            )
            materials[i] = material
        end
    end

    # SHAPES
    shapeFilenames = Vector{String}[]
    if haskey(json, "shapes")
        group = json["shapes"]
        shapeFilenames = Vector{String}(undef, size(group, 1))
        for (i, element) in enumerate(group)
            if !haskey(element, "uri")
                throw(MissingException("uri not present"))
            end
            shapeFilenames[i] = element["uri"]
        end
    end

    # SUBDIVS
    # subdivs = Vector{Subdiv}[]
    # subdivsFilenames = Vector{String}[]
    # if haskey(json, "subdivs")
    #     group = json["subdivs"]
    #     subdivs = Vector{Subdiv}(undef, size(group, 1))
    #     subdivsFilenames = Vector{String}(undef, size(group, 1))
    #     defaultSubdiv::Subdiv = Subdiv()
    #     for (i, element) in enumerate(group)
    #         if !haskey(element, "uri")
    #             throw(MissingException("uri not present"))
    #         end
    #         subdivsFilenames[i] = element["uri"]
    #         shapeIndex = get(element, "shape", defaultSubdiv.shapeIndex)
    #         subs = get(element, "subdivisions", defaultSubdiv.subdivisions)
    #         cmk = get(element, "catmullclark", defaultSubdiv.catmullclark)
    #         smooth = get(element, "smooth", defaultSubdiv.smooth)
    #         disp = get(element, "displacement", defaultSubdiv.displacement)
    #         dispTex =
    #             get(element, "displacement_tex", defaultSubdiv.displacementTex)

    #         if dispTex != -1
    #             dispTex += 1
    #         end
    #         if shapeIndex != -1
    #             shapeIndex += 1
    #         end

    #         subdivs[i] = Subdiv(
    #             shapeIndex = shapeIndex,
    #             subdivisions = subs,
    #             catmullclark = cmk,
    #             smooth = smooth,
    #             displacement = disp,
    #             displacementTex = dispTex,
    #         )
    #     end
    # end

    # INSTANCES
    instances = Vector{Instance}[]
    if haskey(json, "instances")
        group = json["instances"]
        instances = Vector{Instance}(undef, size(group, 1))

        defaultInstance = Instance()
        for (i, element) in enumerate(group)
            f = get(element, "frame", [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
            frame = Frame(f[1:3], f[4:6], f[7:9], f[10:12])

            # +1 because julia arrays starst at 1
            shapeIndex = get(element, "shape", defaultInstance.shapeIndex)
            if shapeIndex != -1
                shapeIndex += 1
            end
            materialIndex =
                get(element, "material", defaultInstance.materialIndex)
            if materialIndex != -1
                materialIndex += 1
            end

            instance = Instance(frame, shapeIndex, materialIndex)
            # ignore params "name", "material"
            instances[i] = instance
        end
    end

    # ENVIRONMENTS
    environments = Vector{Environment}(undef, 0)
    if haskey(json, "environments")
        group = json["environments"]
        environments = Vector{Environment}(undef, size(group, 1))

        defaultEnv = Environment()
        for (i, element) in enumerate(group)
            f = get(element, "frame", [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
            frame = Frame(f[1:3], f[4:6], f[7:9], f[10:12])

            emissionTex = get(element, "emission_tex", defaultEnv.emissionTex)
            if emissionTex != -1
                emissionTex += 1
            end
            environments[i] = Environment(
                frame,
                get(element, "emission", defaultEnv.emission),
                emissionTex,
            )
        end
    end

    # LOAD RESOURCES

    # load shapes
    shapes = Vector{Shape}(undef, size(shapeFilenames, 1))
    # no speedup with threading and collect for now, to try also with large scenes
    @inbounds Threads.@threads for (i, filename) in
                                   collect(enumerate(shapeFilenames))
        path = (dirname(scenePath), filename) |> joinpath
        shape = loadShape(path)
        shapes[i] = shape
    end

    # load subdivs
    # for (i, filename) in enumerate(subdivsFilenames)
    #     path = (dirname(scenePath), filename) |> joinpath
    # end

    # load textures
    # understand which textures are used as normal textures
    normalTexturesIndices = Set(map(m -> m.normalTex, materials))
    # no speedup with threading and collect for now, to try also with large scenes
    @inbounds Threads.@threads for (i, filename) in
                                   collect(enumerate(textureFilenames))
        if i in normalTexturesIndices
            convertToRgb = false
        else
            convertToRgb = true
        end

        path = (dirname(scenePath), filename) |> joinpath

        # check that extension is png
        extension = path[findlast(==('.'), path)+1:end]
        if lowercase(extension) == "png"
            image = loadTexturePng(path, convertToRgb)
            # hdrImage = Matrix{RGB{N0f16}}(undef, 0, 0)
            @inbounds textures[i] = Texture(
                image,
                # hdrImage,
                # false, # png is not linear
                textures[i].nearest,
                textures[i].clamp,
            )
        elseif lowercase(extension) == "hdr"
            image = loadTextureHdr(path)
            # image = Matrix{RGB{N0f8}}(undef, 0, 0)
            @inbounds textures[i] = Texture(
                image,
                # hdrImage,
                # true, # hdr is linear
                textures[i].nearest,
                textures[i].clamp,
            )
        else
            error("only png textures for now!")
        end
    end

    # ignore load subdivs

    scene = Scene(cameras, instances, shapes, textures, materials, environments)

    return scene
end

# loads a texture image in png format
# function loadTexturePng(filename::String)::Matrix{RGBA{N0f8}}
function loadTexturePng(filename::String, toRgb::Bool = true)::Matrix{SVec4f}
    image = load(filename)
    # undo srgb and transform into Vector{SVector}
    if toRgb
        @inbounds imageVector =
            map(x -> srgbToRgb(SVec4f(x.r, x.g, x.b, x.alpha)), image)
    else
        @inbounds imageVector = map(x -> SVec4f(x.r, x.g, x.b, x.alpha), image)
    end
    # reverse!(imageVector, dims = 1)
    return imageVector
end

# loads a texture image in hdr format
# function loadTextureHdr(filename::String)::Matrix{RGB{N0f16}}
function loadTextureHdr(filename::String)::Matrix{SVec4f}
    image = load(filename)
    # undo srgb and transform into Vector{SVector}
    # save("sky_test_pre.png", image)
    # FIX: julia applies srgb automatically when loading hdr
    @inbounds imageVector = map(x -> srgbToRgb(SVec4f(x.r, x.g, x.b, 1)), image)
    # save(
    #     "sky_test_post.png",
    #     map(x -> RGBA(x[1], x[2], x[3], x[4]), imageVector),
    # )
    return imageVector
end

function loadShape(filename::String)::Shape
    ply = load_ply(filename)

    # check properties
    vertexProperties = plyname.(ply["vertex"].properties)

    # load positions
    positions::Vector{SVec3f} = Vector{SVec3f}[]
    if "x" in vertexProperties &&
       "y" in vertexProperties &&
       "z" in vertexProperties
        x = ply["vertex"]["x"]
        y = ply["vertex"]["y"]
        z = ply["vertex"]["z"]

        positions = map(SVec3f, zip(x, y, z))
    end

    # load normals
    normals::Vector{SVec3f} = Vector{SVec3f}[]
    if "nx" in vertexProperties &&
       "ny" in vertexProperties &&
       "nz" in vertexProperties
        nx = ply["vertex"]["nx"]
        ny = ply["vertex"]["ny"]
        nz = ply["vertex"]["nz"]

        normals = map(SVec3f, zip(nx, ny, nz))
    end

    # load texture coordinates
    textureCoords::Vector{SVec2f} = Vector{SVec2f}[]
    if "u" in vertexProperties && "v" in vertexProperties
        u = ply["vertex"]["u"]
        v = ply["vertex"]["v"]
        # putting below lines make floor line disappear
        # elseif "s" in vertexProperties && "t" in vertexProperties
        #     u = ply["vertex"]["s"]
        #     v = ply["vertex"]["t"]

        # ignoring filpv value for now

        @inbounds textureCoords =
            map((x) -> SVec2f(x[1], 1 - x[2]), collect(zip(u, v)))
    end

    # ignoring load colors for now

    # load radius
    radius::Vector{Float32} = Vector{Float32}[]
    if "radius" in vertexProperties
        plyradius = ply["vertex"]["radius"]
        @inbounds radius = map(x -> x[1], plyradius)
    end

    triangles::Vector{SVec3i} = Vector{SVec3i}[]
    quads::Vector{SVec4i} = Vector{SVec4i}[]
    lines::Vector{SVec2i} = Vector{SVec2i}[]
    points::Vector{Int32} = Vector{Int32}[]

    if "face" in plyname.(ply.elements)
        faces = ply["face"]["vertex_indices"]

        # count triangles and quads
        @inbounds triCount =
            count(elem -> (size(elem[1])[end] == 3), eachrow(faces))
        @inbounds quadCount =
            count(elem -> (size(elem[1])[end] == 4), eachrow(faces))
        if triCount + quadCount != size(faces)[1]
            error("Only implemented triangles and quads as faces")
        end

        # load triangles and quads

        if quadCount > 0
            @inbounds quads = map(
                x ->
                    size(x[1])[end] == 3 ?
                    SVec4i(x[1][1] + 1, x[1][2] + 1, x[1][3] + 1, x[1][3] + 1) : SVec4i(x[1] .+ 1...),
                eachrow(faces),
            )
        elseif triCount > 0
            @inbounds triangles =
                map(x -> SVec3i((x[1] .+ 1)...), eachrow(faces))
        else
            error("Zero triangles and zero quads. No bueno?")
        end
    elseif "line" in plyname.(ply.elements)
        # plylines = ply["line"]["vertex_indices"]
        lines = map(
            x -> SVec2i((x[1] .+ 1)...),
            eachrow(ply["line"]["vertex_indices"]),
        )
    elseif "point" in plyname.(ply.elements)
        error("not sure how to load points")
        # points = map(x -> x[1] + 1, eachrow(ply["point"]["vertex_indices"]))
    end

    if (
        isempty(points) &&
        isempty(lines) &&
        isempty(triangles) &&
        isempty(quads)
    )
        error("No faces, lines, or points found")
    end

    # create shape
    return Shape(
        quads,
        triangles,
        lines,
        points,
        positions,
        normals,
        textureCoords,
        radius,
    )
end
