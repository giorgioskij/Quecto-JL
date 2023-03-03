#! /bin/zsh

~/dev/uni/computer-graphics/masters/yocto-gl/bin/ytrace \
    --scene ~/dev/uni/computer-graphics/masters/Jtrace/04_envlight/envlight.json \
    --samples 32\
    --resolution 1920 \
    --output ~/dev/uni/computer-graphics/masters/Jtrace/out/yocto.png \
    --sampler naive
