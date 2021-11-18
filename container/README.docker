#See https://docs.docker.com/engine/reference/commandline/build/
#    https://stackify.com/docker-build-a-beginners-guide-to-building-docker-images/
#    https://docs.docker.com/engine/reference/commandline/exec/
#    https://developers.redhat.com/blog/2016/03/09/more-about-docker-images-size/


## mnxtools:0.0.3        is the built image target, usually (lowercase) name:version
## Dockerfile            is the Dockerfile, the commands used to build the image
## Dockerfile example: https://github.com/ElmerCSC/elmerfem/blob/devel/docker/elmerice.dockerfile
#                      https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/Dockerfile
#                      https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit/blob/master/Dockerfile

# Build Docker image
export MNXTOOLS_VERSION=0.0.3
docker build -t mnxtools:$MNXTOOLS_VERSION --no-cache=true --build-arg BUILD_DATE=$(date -u +'%Y-%m-%dT%H:%M:%SZ') --build-arg MNXTOOLS_VERSION=$MNXTOOLS_VERSION  -f Dockerfile .  2>&1 >MNXtools.Dockerfile.log

# List Docker local images (imported or built)
docker images

# Inspect images
docker inspect  mnxtools:$MNXTOOLS_VERSION


# Run bash in the Docker image
docker run --name mnxtools --rm -i -t mnxtools:$MNXTOOLS_VERSION bash
	# Mounting/binding a local repository (,readonly can be added to force readonly mounting)
	--mount type=bind,source=/software,target=/software

