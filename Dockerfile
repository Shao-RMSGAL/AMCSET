FROM ubuntu:20.04

LABEL maintainer="Nathaniel Thomas <nathaniel@swbell.net>"
LABEL description="Custom AMCSET development environment based on ubuntu-dev"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y \
        sudo \
        python3 \
        python3-pip \
        git \
        software-properties-common \
        curl \
        wget \
        coreutils \
        build-essential \
    && rm -rf /var/lib/apt/lists/*
     
RUN echo "developer ALL=(ALL) NOPASSWD: /usr/bin/chown, /usr/bin/dpkg" >> /etc/sudoers

RUN echo "root:rootpassword" | chpasswd \
    && useradd -m developer \
    && echo "developer:developer" | chpasswd \
    && adduser developer sudo

USER developer

WORKDIR /home/developer

COPY . ./AMCSET

RUN sudo chown -R developer:developer /home/developer/

ENV PATH="/home/developer/.local/bin:${PATH}"

COPY pip_requirements.txt /tmp/

# Download the Conan Debian package from the specified URL
RUN wget https://github.com/conan-io/conan/releases/download/2.4.0/conan-2.4.0-amd64.deb \
    && sudo dpkg -i conan-2.4.0-amd64.deb \
    && rm conan-2.4.0-amd64.deb 

RUN conan profile detect --force \
    && cd AMCSET \
    && conan install . \
    && conan install . -s build_type=Debug

RUN echo ". /home/developer/AMCSET/build/Debug/generators/conanbuild.sh" >> /home/developer/.bashrc

WORKDIR /home/developer/AMCSET/

CMD ["bash"]
