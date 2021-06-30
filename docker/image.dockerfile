# This file tells docker what image must be created
# in order to be ahble to test this library
FROM ubuntu:20.04


ENV TZ=Europe/Rome
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV UBUNTU_RELEASE=bionic


RUN apt-get update

# Install packages
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends -o Dpkg::Options::="--force-confnew" \
                    python3-pip git iputils-ping net-tools netcat screen build-essential lsb-release gnupg2 curl less \
                    python3-sympy coinor-libipopt-dev sudo valgrind \
                    build-essential pkg-config git \
                    liblapack-dev liblapack3 libopenblas-base libopenblas-dev \
                    libgfortran-7-dev cmake libgsl-dev gdb python3-tk libeigen3-dev \
                    libboost-math-dev

RUN pip3 install setuptools matplotlib Mosek scipy quadpy six cython tk

RUN sh -c 'echo "deb http://packages.ros.org/ros/ubuntu $(lsb_release -sc) main" > /etc/apt/sources.list.d/ros-latest.list'
RUN apt-key adv --keyserver 'hkp://keyserver.ubuntu.com:80' --recv-key C1CF6E31E6BADE8868B172B4F42ED6FBAB17C654
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends -o Dpkg::Options::="--force-confnew" \
                    ros-noetic-ifopt

# user handling
ARG myuser
ARG myuid
ARG mygroup
ARG mygid
ARG scriptdir
RUN addgroup --gid ${mygid} ${mygroup} --force-badname
RUN adduser --gecos "" --disabled-password  --uid ${myuid} --gid ${mygid} ${myuser} --force-badname
#add user to sudoers
RUN echo "${myuser} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
WORKDIR /gsplinespp

COPY vim_installation.bash /
RUN cd / && bash vim_installation.bash
COPY configfiles/vimrc /etc/vim/
COPY configfiles/ycm_extra_conf.py /etc/vim/
