---
title: "Docker practice"
date: '2017-08-23'
slug: docker
categories: ["data at fingertips"]
tags: ["Bioinformatics","Docker"]
---

### Basic commands

Run a container from an image  
`docker run -it my_image`    
The `-it` instructs Docker to allocate a pseudo-TTY connected to the container’s stdin; creating an interactive bash shell in the container. 

Quit from a container with interactive bash shell  
`exit` or  
`ctrl+p` followed by `ctrl+q`

Stop a container  
`docker stop my_container`

Remove a container (only stopped containers can be removed)  
`docker rm my_container`

Stop and remove all containers  
`docker stop $(docker ps -a -q)`  
`docker rm $(docker ps -a -q)`

Remove images  
`docker rmi my_image`

List all exited containers  
`docker ps -aq -f status=exited`  

Remove stopped containers  
`docker ps -aq -f status=exited --no-trunc | xargs docker rm`

Remove dangling/untagged images  
`docker images -q --filter dangling=true | xargs docker rmi`

Remove containers created by a specific image  
`docker ps -a | grep 'my_image' | awk '{print $1}' | xargs docker rm`


### Using existing images from docker hub

Take one data processing tool, <a href=https://hub.docker.com/r/librecat/catmandu/ target="_blank">catmandu</a>, for example. 

Pull image from docker hub:  
`docker pull librecat/catmandu`

Start a container from the image, and mount a local directory onto /home/catmandu/Home directory inside the container:  
`docker run -v /home/yourdirectory:/home/catmandu/Home -it librecat/catmandu`  
The local directory should be writable for all users, which means, the permission triads look like "drwxrwxrwx". To change the permission of the directory, use `chmod 777 /home/yourdirectory` 

If you have a file named "test.json" in /home/yourdirectory/, and would like to convert it to CSV format, once inside the container, you should be able to find the same file under /home/catmandu/Home/. Then:  

`catmandu convert JSON to CSV < /home/catmandu/Home/test.json > /home/catmandu/Home/test.csv`

The converted file "test.csv" would appear in /home/yourdirectory/ as well.

### Customize dockerfile

Some example Dockerfile of commonly-used tools in bioinformatics:  

<a href=https://github.com/zhengh42/Dockerfiles target="_blank">zhengh42/Dockerfiles</a>

Take the Dockerfile of samtools/bcftools/htslib for example:
 
```
FROM ubuntu:16.04
MAINTAINER Hong Zheng <zhengh42@stanford.edu>

# Install dependencies
RUN apt-get update && apt-get install -y \
  curl \
  build-essential \
  libncurses-dev \
  zlib1g-dev \
  libz-dev \
  libbz2-dev \
  liblzma-dev

# Set up
ENV SAMTOOLS_RELEASE=1.5
ENV SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/
ENV BCFTOOLS_RELEASE=1.5
ENV BCFTOOLS_URL=https://github.com/samtools/bcftools/releases/download/
ENV HTSLIB_RELEASE=1.5
ENV HTSLIB_URL=https://github.com/samtools/htslib/releases/download/
ENV DEST_DIR=/opt/

# Download; untar & decompress; remove tr.bz2 file; compile & install samtools; remove unnecessary files
RUN curl -SLo ${DEST_DIR}/samtools-${SAMTOOLS_RELEASE}.tar.bz2 ${SAMTOOLS_URL}/${SAMTOOLS_RELEASE}/samtools-${SAMTOOLS_RELEASE}.tar.bz2 && \
  tar -xf ${DEST_DIR}/samtools-${SAMTOOLS_RELEASE}.tar.bz2 -C ${DEST_DIR} && \
  rm ${DEST_DIR}/samtools-${SAMTOOLS_RELEASE}.tar.bz2 && \
  cd ${DEST_DIR}/samtools-${SAMTOOLS_RELEASE} && \
  ./configure && \
  make && \
  make install && \
  rm -rf ${DEST_DIR}/samtools-${SAMTOOLS_RELEASE} && \
  curl -SLo ${DEST_DIR}/bcftools-${BCFTOOLS_RELEASE}.tar.bz2 ${BCFTOOLS_URL}/${BCFTOOLS_RELEASE}/bcftools-${BCFTOOLS_RELEASE}.tar.bz2 && \
  tar -xf ${DEST_DIR}/bcftools-${BCFTOOLS_RELEASE}.tar.bz2 -C ${DEST_DIR} && \
  rm ${DEST_DIR}/bcftools-${BCFTOOLS_RELEASE}.tar.bz2 && \
  cd ${DEST_DIR}/bcftools-${BCFTOOLS_RELEASE} && \
  ./configure && \
  make && \
  make install && \
  rm -rf ${DEST_DIR}/bcftools-${HTSLIB_RELEASE} && \
  curl -SLo ${DEST_DIR}/htslib-${HTSLIB_RELEASE}.tar.bz2 ${HTSLIB_URL}/${HTSLIB_RELEASE}/htslib-${HTSLIB_RELEASE}.tar.bz2 && \
  tar -xf ${DEST_DIR}/htslib-${HTSLIB_RELEASE}.tar.bz2 -C ${DEST_DIR} && \
  rm ${DEST_DIR}/htslib-${HTSLIB_RELEASE}.tar.bz2 && \
  cd ${DEST_DIR}/htslib-${HTSLIB_RELEASE} && \
  ./configure && \
  make && \
  make install && \
  rm -rf ${DEST_DIR}/htslib-${HTSLIB_RELEASE}

CMD ["/bin/bash"]
```

Build new image from the Dockerfile:  
`docker build -f samtools/1.5all/Dockerfile -t zhengh42/samtools:1.5all samtools/1.5all`

After that, an image named "zhengh42/samtools" would appear if we run `docker images`.

How to use this samtools from this image?

`docker run -v /home/yourdirectory:/mnt zhengh42/samtools:1.5all samtools view /mnt/yourfiles`

/home/yourdirectory stores your files, and it is mounted into /mnt inside docker containers.

