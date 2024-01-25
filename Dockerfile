# syntax=docker/dockerfile:1

# ## Start with python3.9 on Debian base
FROM python:3.9-bookworm as base

# Install libopencv-dev and numpy in the base, as it is needed in the build and 
# final images
RUN apt-get update
RUN pip install numpy matplotlib

# ## Create a stage for building/compiling the application.
FROM base as build

RUN apt-get install -y libopencv-dev
RUN apt-get install -y cmake 

COPY . /opt
WORKDIR /opt
RUN mkdir build && cd build; \
    cmake .. ; \
    make

# ## Create a final stage for running your application.
FROM base AS final

# Copy the executable from the "build" stage.
COPY --from=build /opt/build/eSpecAnalysis /bin

# Install and copy the required OpenCV libraries
RUN apt-get install -y libopengl0 libqt5test5 libqt5opengl5

RUN mkdir /usr/local/lib/opencv
RUN export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

COPY --from=build /usr/lib/x86_64-linux-gnu/libopencv_core.so.* /usr/local/lib
COPY --from=build /usr/lib/x86_64-linux-gnu/libopencv_highgui.so.* /usr/local/lib
COPY --from=build /usr/lib/x86_64-linux-gnu/libopencv_imgcodecs.so.* /usr/local/lib
COPY --from=build /usr/lib/x86_64-linux-gnu/libopencv_calib3d.so.* /usr/local/lib
COPY --from=build /usr/lib/x86_64-linux-gnu/libopencv_imgproc.so.* /usr/local/lib
#COPY --from=build /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.* /usr/local/lib
#COPY --from=build /usr/lib/x86_64-linux-gnu/libQt5Gui.so.* /usr/local/lib
COPY --from=build /usr/lib/x86_64-linux-gnu/libgdcmMSFF.so.* /usr/local/lib
COPY --from=build /usr/lib/x86_64-linux-gnu/libgdcmDSED.so.* /usr/local/lib
COPY --from=build /usr/lib/x86_64-linux-gnu/libgdal.so.* /usr/local/lib
COPY --from=build /usr/lib/x86_64-linux-gnu/libopencv_features2d.so.* /usr/local/lib
COPY --from=build /usr/lib/x86_64-linux-gnu/libopencv_flann.so.* /usr/local/lib
COPY --from=build /usr/lib/x86_64-linux-gnu/libGLX.so.* /usr/local/lib
COPY --from=build /usr/lib/x86_64-linux-gnu/liblapack.so.* /usr/local/lib

# Create a non-privileged user that the app will run under.
# See https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#user
ARG UID=10001
RUN adduser \
    --disabled-password \
    --gecos "" \
    --home "/nonexistent" \
    --shell "/sbin/nologin" \
    --no-create-home \
    --uid "${UID}" \
    appuser
USER appuser

# Create folder for input and output data
WORKDIR /home/appuser
RUN mkdir data && \
    mkdir data/ePointing; \
    mkdir data/eScreenA; \
    mkdir data/eScreenB

# Copy configuration files
COPY config/settings.cfg /home/appuser
COPY config/Calibration /home/appuser/Calibration

# What the container should run when it is started.
ENTRYPOINT [ "/bin/eSpecAnalysis" ]
