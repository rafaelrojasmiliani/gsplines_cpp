# This file tells docker what image must be created
# in order to be ahble to test this library.
FROM rafa606/ros_noetic_vim

RUN git clone --recursive https://github.com/rafaelrojasmiliani/gsplines_cpp.git /gsplines && cd /gsplines && mkdir build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX=/usr && make && make install \
   && rm -rf /gsplines
