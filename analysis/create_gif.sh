#!/bin/bash

path='plots/du10-z2/'

#ffmpeg -i $path/light_map_%04d.png -vf palettegen $path/palette.png
#ffmpeg -f image2 -framerate 15 -i $path/light_map_%04d.png -i $path/palette.png -lavfi paletteuse -loop 0 $path/du10-gas-short.gif

ffmpeg -y -f image2 -r 20 -i $path/light_map_%04d.png -b 2500k -vcodec libx264 -pix_fmt yuv420p -threads 4 -start_number 0 -vframes 3000 $path/z2-3Gyr.mp4

