FROM ubuntu:22.04

RUN set -x && apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install \
                    -y --no-install-recommends -o Dpkg::Options::="--force-confnew" \
            build-essential \
            cmake \
            coinor-libipopt-dev \
            g++ \
            git \
            libeigen3-dev \
            python3-dev \
            python3-pip \
    && pip3 install \
            autopep8 \
            cmake-format \
            flake8 \
            flask \
            pip \
            sympy \
            numpy \
            scipy \
            pandas \
            tk \
            matplotlib \
    && git clone https://github.com/ethz-adrl/ifopt.git /ifopt \
    && cd /ifopt && mkdir build && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr .. \
    && make -j$(nproc) && make install \
    && cd / \
    && rm -rf /ifopt \
    && git clone --recursive https://github.com/rafaelrojasmiliani/gsplines_cpp.git /gsplines \
    && cd /gsplines && mkdir build &&  cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr .. \
    && make -j$(nproc) && make install \
    && cd / \
    && rm -rf /gsplines
