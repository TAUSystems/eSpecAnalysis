# syntax=docker/dockerfile:1

# ## Start with python3.9 on Debian base
FROM python:3.9-slim-bookworm as base

# Install libopencv-dev and numpy in the base, as it is needed in the build and 
# final images
RUN apt-get update
RUN pip install numpy matplotlib

# ## Create a stage for building/compiling the application.
FROM base as build

RUN apt-get install -y libopencv-dev cmake g++

COPY . /opt
WORKDIR /opt
RUN mkdir build && cd build; \
    cmake .. ; \
    make

# ## Create a final stage for running your application.
FROM base AS final

# only install the libraries that are linked against
RUN apt-get install -y libopencv-core-dev libopencv-calib3d-dev \
                       libopencv-highgui-dev libopencv-imgproc-dev \
                       libopencv-features2d-dev

# Copy the executable from the "build" stage.
COPY --from=build /opt/build/eSpecAnalysis /bin

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

# Copy configuration files
COPY config/settings.cfg /home/appuser
COPY config/Calibration /home/appuser/Calibration

# What the container should run when it is started.
ENTRYPOINT [ "/bin/eSpecAnalysis" ]
