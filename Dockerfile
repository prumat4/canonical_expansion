FROM ubuntu:latest

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libgmp-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY . /app

RUN rm -rf build && mkdir build && cd build && cmake .. && make -j4

ENTRYPOINT ["./build/canonicalExpansionApp"]