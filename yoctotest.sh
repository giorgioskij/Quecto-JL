#! /bin/zsh

../yocto-gl/bin/ytrace \
    --scene 04_envlight/envlight.json \
    --samples 100 \
    --output out/yocto.png \
    --sampler eyelight 