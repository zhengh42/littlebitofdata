---
title: "Docker practice"
date: '2017-08-23'
slug: docker
categories: ["data at fingertips"]
tags: ["Bioinformatics","Docker"]
---


Run a container from an image  
`docker run -it IAMGE`    
The `-it` instructs Docker to allocate a pseudo-TTY connected to the container’s stdin; creating an interactive bash shell in the container. 

Quit from a container with interactive bash shell  
`exit 13` or  
`ctrl+p` followed by `ctrl+q`


Remove images  
`docker rmi ***`

List all exited containers  
`docker ps -aq -f status=exited`  

Remove stopped containers  
`docker ps -aq --no-trunc | xargs docker rm`

Remove dangling/untagged images  
`docker images -q --filter dangling=true | xargs docker rmi`

Remove containers created after a specific container  
`docker ps --since *** -q | xargs docker rm`

Remove containers created before a specific container  
`docker ps --before *** -q | xargs docker rm`