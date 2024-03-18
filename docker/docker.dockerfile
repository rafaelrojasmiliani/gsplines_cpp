ARG BASEIMAGE
FROM ${BASEIMAGE}
ARG ROS_DISTRO

RUN --mount=type=bind,source=./,target=/workspace,rw \
    cd /workspace \
    && mkdir build \
    && cd build \
    && source /opt/ros/${ROS_DISTRO}/setup.bash \
    && cmake .. -DCMAKE_INSTALL_PREFIX=/usr  -DBUILD_TESTING=OFF \
    && make gsplines -j $(nproc) \
    && make install
