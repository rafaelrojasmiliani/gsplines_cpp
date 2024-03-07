ARG BASEIMAGE
FROM ${BASEIMAGE}

RUN --mount=type=bind,source=./,target=/workspace,rw \
    cd /workspace \
    && mkdir build \
    && cd build \
    && cmake .. -DCMAKE_INSTALL_PREFIX=/usr \
    && make gsplines -j $(nproc) \
    && make install \
    && cd / \
    && rm -rf /gsplines
