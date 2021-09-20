#!/usr/bin/env vash

myuid=$(id -u $USER)
mygid=$(id -g $USER)
mygroup=$(id -g -n $USER)

docker build -t "gsplinespp" \
    --build-arg myuser="$USER" \
    --build-arg myuid="$myuid" \
    --build-arg mygroup="$mygroup" \
    --build-arg mygid="$mygid" \
    --no-cache  -f ./image.dockerfile .

exit 0
