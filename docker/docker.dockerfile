ARG BASEIMAGE
FROM ${BASEIMAGE}
ARG ROS_DISTRO

RUN --mount=type=bind,source=./,target=/workspace,rw \
    cd /workspace \
    && mkdir build \
    && cd build \
    && source /opt/${ROS_DISTRO}/setup.bash \
    && cmake .. -DCMAKE_INSTALL_PREFIX=/usr \
    && make gsplines -j $(nproc) \
    && make install \
    && cd / \
    && rm -rf /gsplines
