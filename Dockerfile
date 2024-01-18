# syntax=docker/dockerfile:1

FROM python:3.9-bookworm as base

################################################################################
# Create a stage for building/compiling the application.
#
# The following commands will leverage the "base" stage above to generate
# a "hello world" script and make it executable, but for a real application, you
# would issue a RUN command for your application's build process to generate the
# executable. For language-specific examples, take a look at the Dockerfiles in
# the Awesome Compose repository: https://github.com/docker/awesome-compose

# FROM base as build

RUN apt-get update
RUN apt-get install -y libopencv-dev cmake 
# mingw-w64 
RUN pip install numpy

# WORKDIR /tmp
# RUN git clone https://github.com/opencv/opencv.git && \
#     cd opencv; \
#     mkdir build && cd build; \
#     cmake \
#         -DCMAKE_SYSTEM_NAME=Windows \
#         -DCMAKE_INSTALL_PREFIX=${MINGW} \
#         -DCMAKE_FIND_ROOT_PATH=${MINGW} \
#         -DCMAKE_FIND_ROOT_PATH_MODE_PROGRAM=NEVER \
#         -DCMAKE_FIND_ROOT_PATH_MODE_LIBRARY=ONLY \
#         -DCMAKE_FIND_ROOT_PATH_MODE_INCLUDE=ONLY \
#         -DCMAKE_C_COMPILER=x86_64-w64-mingw32-gcc \
#         -DCMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++ \
#         -DCMAKE_RC_COMPILER=x86_64-w64-mingw32-windres \
#     .. ; \
#     make -j8; \
#     make install; 

COPY . /opt
WORKDIR /opt
RUN mkdir build && cd build; \
    cmake \
        # -DCMAKE_CXX_STANDARD=11 \
        # -DCMAKE_SYSTEM_NAME=Windows \
        # -DCMAKE_INSTALL_PREFIX=${MINGW} \
        # -DCMAKE_FIND_ROOT_PATH=${MINGW} \
        # -DCMAKE_FIND_ROOT_PATH_MODE_PROGRAM=NEVER \
        # -DCMAKE_FIND_ROOT_PATH_MODE_LIBRARY=ONLY \
        # -DCMAKE_FIND_ROOT_PATH_MODE_INCLUDE=ONLY \
        # -DCMAKE_C_COMPILER=x86_64-w64-mingw32-gcc \
        # -DCMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++ \
        # -DCMAKE_RC_COMPILER=x86_64-w64-mingw32-windres \
    .. ; \
    make

################################################################################
# Create a final stage for running your application.
#
# The following commands copy the output from the "build" stage above and tell
# the container runtime to execute it when the image is run. Ideally this stage
# contains the minimal runtime dependencies for the application as to produce
# the smallest image possible. This often means using a different and smaller
# image than the one used for building the application, but for illustrative
# purposes the "base" image is used here.

# FROM base AS final

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

# # Copy the executable from the "build" stage.
# COPY --from=build /home/* /home/

# # What the container should run when it is started.
# ENTRYPOINT [ "/bin/hello.sh" ]
