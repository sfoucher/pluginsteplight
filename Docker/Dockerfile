# Use Ubuntu 22.04 as base image
FROM ubuntu:22.04

# Set environment variables to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Update package list and upgrade system packages to latest versions
RUN apt-get update && apt-get upgrade -y

# Install system dependencies required for building and running the application
RUN apt-get install -y \
    software-properties-common \
    build-essential \
    subversion \
    wget \
    curl \
    git \
    cmake \
    pkg-config \
    gcc \
    g++ \
    make \
    qtcreator \
    qtbase5-dev \
    qttools5-dev \
    qttools5-dev-tools \
    libeigen3-dev \
    libmuparser-dev \
    libgdal-dev \
    gdal-bin \
    libopencv-dev \
    python3-opencv \
    libpcl-dev \
    libboost-all-dev \
    libgsl-dev \
    gsl-bin \
    libqhull-dev \
    libflann-dev \
    liblaszip-dev \
    && rm -rf /var/lib/apt/lists/*

# Add Ubuntu toolchain repository for newer GCC versions
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
# Install GCC-11 and G++-11 compilers with updated libstdc++
RUN apt install -y gcc-11 g++-11 libstdc++6

# Enable Qt debug plugins for troubleshooting Qt-related issues
ENV QT_DEBUG_PLUGINS=1

# Copy the Computree installer to the container
COPY Computree_Installer_LINUX_6.0.392 .
# Make the installer executable
RUN chmod +x ./Computree_Installer_LINUX_6.0.392
# Run the Computree installer with automatic license acceptance and installation
RUN ./Computree_Installer_LINUX_6.0.392 --al --da -c install

# Create a non-root user for running Computree (security best practice)
RUN useradd -m -s /bin/bash computree_user && \
    chown -R computree_user:computree_user /opt/CompuTree

# Switch to non-root user for security
USER computree_user

# Set working directory to Computree installation directory
WORKDIR /opt/CompuTree

# Expose port 8080 for potential web interface or API
EXPOSE 8080

# Set the default command to start a bash shell
CMD ["bash"]
