# MNXtools

!!! This is work in progress - nothing below is yet completely true !!!

The code in this repository is primarily designed to build a dockerized application. The latest stable docker version can be obtained from dockerhub: https://hub.docker.com/r/sibswiss/mnxtools.   

see also: https://www.metanetx.org

## How to download MNXref cache
```
export MNXREF_VERSION=4.4;
cd cache/ && wget "https://www.metanetx.org/ftp/$MNXREF_VERSION/ChemSpace.bindump.gz" && gunzip ChemSpace.bindump.gz
```
