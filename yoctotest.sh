#! /bin/zsh

~/dev/uni/computer-graphics/masters/yocto-gl/bin/ytrace \
    --scene ~/dev/uni/computer-graphics/masters/Jtrace/02_matte/matte.json \
    --samples 100 \
    --resolution 3000 \
    --output ~/dev/uni/computer-graphics/masters/Jtrace/out/yocto.png \
    --sampler eyelight
