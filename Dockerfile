FROM ubuntu:20.04

MAINTAINER SergeDmi www.biophysics.fr

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y make cmake g++ libboost-all-dev git

RUN git clone https://github.com/SergeDmi/Plyssim.git Plyssim
WORKDIR Plyssim
RUN ./setup

