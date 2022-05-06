FROM rafa606/ros-galactic-vim
SHELL ["bash", "-c"]
RUN source /opt/ros/galactic/setup.bash && git clone --recursive https://github.com/rafaelrojasmiliani/gsplines_cpp.git /gsplines && cd /gsplines && mkdir build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX=/usr && make && make install \
   && rm -rf /gsplines
