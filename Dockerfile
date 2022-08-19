FROM ubuntu:20.04

MAINTAINER jmonlong@ucsc.edu

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    wget \
    gcc \ 
    samtools \
    build-essential \
    bzip2 \
    git \
    sudo \
    less \
    pkg-config \
    cmake \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libbz2-dev \
    autoconf \
    apt-transport-https software-properties-common dirmngr gpg-agent \ 
    && rm -rf /var/lib/apt/lists/*

## install R
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    r-base \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('ggplot2', 'dplyr'))"

## install wambam
WORKDIR /build/wambam

ADD . /build/wambam/

RUN mkdir build && cd build && \
    cmake .. && make

ENV PATH=/build/wambam/build:$PATH

WORKDIR /home
