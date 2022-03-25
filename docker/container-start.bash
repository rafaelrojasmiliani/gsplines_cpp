
main(){

    scriptdir=$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)
    if [ "$scriptdir" != "$(pwd)" ]; then
      echo "this script must be executed from $scriptdir".
      exit 1
    fi

    if lspci | grep -qi "vga .*nvidia" && \
        docker -D info 2>/dev/null | grep -qi "runtimes.* nvidia"; then
        DOCKER_NVIDIA_OPTIONS="\
            --runtime=nvidia \
            --env=NVIDIA_VISIBLE_DEVICES=${NVIDIA_VISIBLE_DEVICES:-all} \
            --env=NVIDIA_DRIVER_CAPABILITIES=${NVIDIA_DRIVER_CAPABILITIES:+$NVIDIA_DRIVER_CAPABILITIES,}all"
    else
        DOCKER_NVIDIA_OPTIONS=""
    fi

    XAUTH=/tmp/.docker.xauth
    touch $XAUTH
    xauth nlist $DISPLAY | sed -e 's/^..../ffff/' | xauth -f $XAUTH nmerge -

    DOCKER_VIDEO_OPTIONS="${DOCKER_NVIDIA_OPTIONS} --env=DISPLAY --env=QT_X11_NO_MITSHM=1 --env=XAUTHORITY=$XAUTH --volume=$XAUTH:$XAUTH --volume=/tmp/.X11-unix:/tmp/.X11-unix:rw"

    myuid=$(id -u $USER)
    mygid=$(id -g $USER)
    mygroup=$(id -g -n $USER)
    myuser="$USER"
    docker pull rafa606/cpp-vim
    docker run -it --rm \
        ${DOCKER_VIDEO_OPTIONS} \
        --volume $(pwd)/../:/workspace/src/tutorial: \
        --entrypoint="/bin/bash" \
        --privileged \
        "rafa606/cpp-vim" -c "addgroup --gid ${mygid} ${mygroup} --force-badname;  adduser --gecos \"\" --disabled-password  --uid ${myuid} --gid ${mygid} ${myuser} --force-badname ; usermod -a -G video ${myuser}; echo ${myuser} ALL=\(ALL\) NOPASSWD:ALL >> /etc/sudoers; sudo -EHu ${myuser}  bash"
}

main
