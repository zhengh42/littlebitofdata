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
`docker rm ***`

Stop and remove all containers  
`docker stop $(docker ps -a -q)`  
`docker rm $(docker ps -a -q)`

Remove images  
`docker rmi ***`

List all exited containers  
`docker ps -aq -f status=exited`  

Remove stopped containers  
`docker ps -aq -f status=exited --no-trunc | xargs docker rm`

Remove dangling/untagged images  
`docker images -q --filter dangling=true | xargs docker rmi`

Remove containers created after a specific container  
`docker ps --since *** -q | xargs docker rm`

Remove containers created before a specific container  
`docker ps --before *** -q | xargs docker rm`

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

Very useful example Dockerfile of commonly-used tools in bioinformatics:  

<a href=https://github.com/Duke-GCB/GCB-Dockerfiles target="_blank">GCB-Dockerfile</a>

Take the Dockerfile of bedtools for example:
 
```
FROM ubuntu:16.04
LABEL maintainer="john.bradley@duke.edu"

# picard requires java
RUN apt-get update && apt-get install -y \
  wget \
  openjdk-8-jre-headless

# Installs fastqc from compiled java distribution into /opt/FastQC
ENV PICARD_VERSION="2.10.7"
ENV PICARD_URL https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar

WORKDIR /opt/picard
RUN wget $PICARD_URL

CMD ["java", "-jar", "picard.jar"]
```

Build new image from the Dockerfile:  
`docker build -f picard/2.10.7/Dockerfile -t dukegcb/picard:2.10.7 picard/2.10.7`

After that, an image named "dukegcb/picard" would appear if we run `docker images`.

How to use picard tools such as CollectHsMetrics from this image? 

`docker run dukegcb/picard:2.10.7 java -jar picard.jar CollectHsMetrics [other arguments]`


