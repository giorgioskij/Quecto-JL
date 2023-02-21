# julia-pathtracer
A pathtracer in Julia inspired by Yocto-GL


## Yocto-GL structure
**ytrace.cpp:**  
main function `run`, loops over **samples**, calls `trace_samples`
stores empty image in `state.image`

**yocto_trace.cpp:**
- `trace_samples` loops over pixels, calls `trace_sample`. Writes over image
- `trace_sample` runs the selective shader (recursively). Returns nothing, 
    writes on the image

**ytrace.cpp:**   
now image is written, saves with `save_image`


## Scene Structure
- Vector{Int} points
- Vector{Sphere} spheres
- Vector{Triangles} tri
- Camera


Triangles = Vec3i
Spheres = {Int, Float}

## Schedule
- [x] Initialize flow and save first image 
- [x] Define structures for camera, frame, ray, etc...
- [x] Spheres as {index to point, radius}
- [x] Hit spheres with rays
- [x] Normal shading 
- [x] Antialiasing
- [x] Eyelight shading
- [x] Load triangles
- [ ] Evaluate normals
- [ ] Textures
- [ ] Environment maps
- [ ] Hit other stuff with rays


## Timings 


### Tests with shader normal on 02_matte/matte.json

(1920, 1080, 2): 6.5 s
(3000, 2000, 2): 2s
(300,200,2): 389.5 +- 42 ms

for this scene the construction of the bvh is roughly 230 ms
so only test much bigger than that are significative for rendering speed

Fixing randoms: (1920, 1080, 2) 6.5 -->  5.3 seconds

## Benchmarking
In a starting state, with no computation in the shader, representing the image
with SVectors as colors is slower.


## Long term todo:

Care about RNGs