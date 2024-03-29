# Canonical Expansion Project

## Short Description
This project implements an algorithm to find the canonical expansion of a given number. It incorporates various factorization methods, including trial divisions, Pollard's Rho method, and others, to decompose numbers into their prime factors.

## Required Software
- C++ compiler (GCC or Clang)
- CMake (version 3.10 or higher)
- Make
- Docker (for building and running within a container)
- GMP (GNU Multiple Precision Arithmetic Library)

## How to Use

### Building and Running the Executable
1. Clone the repository and navigate to the project directory.
2. Build the project using CMake and Make. Use the following commands:
    ```sh
    cmake . && make -j12
    ```
3. Run the executable with a number as an argument to find its canonical expansion. For example:
    ```sh
    ./canonicalExpansion 666
    ```

### Docker hub link:
https://hub.docker.com/repository/docker/omega111111/factorization/general
### Building and Running with Docker
1. Ensure Docker is installed and running on your system.
2. Build the Docker image from the Dockerfile provided in the project directory. Replace `your-image-name` with your preferred name for the Docker image:
    ```sh
    docker build -t your-image-name .
    ```
3. Once the image is built, run it using Docker. The following command runs the container and executes the canonical expansion algorithm with a specified number (e.g., 666):
    ```sh
    docker run your-image-name ./canonicalExpansion 666
    ```
Replace `your-image-name` with the name you used when building the Docker image.

By following these steps, you can compile and run the canonical expansion algorithm both directly on your machine and within a Docker container.
