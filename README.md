[![DOI](https://zenodo.org/badge/DOI/10.1093/nar/gkaa992.svg)](https://dx.doi.org/10.1093/nar/gkaa992)
[![DOI](https://zenodo.org/badge/DOI/10.1093/nar/gkv1117.svg)](https://dx.doi.org/10.1093/nar/gkv1117)
[![DOI](https://zenodo.org/badge/DOI/10.1093/bioinformatics/btt036.svg)](https://dx.doi.org/10.1093/bioinformatics/btt036)
[![DOI](https://zenodo.org/badge/DOI/10.1093/bib/bbs058.svg)](https://dx.doi.org/10.1093/bib/bbs058)

[![Bluesky](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fpublic.api.bsky.app%2Fxrpc%2Fapp.bsky.actor.getProfile%2F%3Factor%3Dmetanetx.bsky.social&query=%24.followersCount&style=social&logo=bluesky&label=Follow%20%40MetaNetX)](https://bsky.app/profile/metanetx.bsky.social)
[![Mastodon](https://img.shields.io/mastodon/follow/110870430059853959?label=Follow%20%40MetaNetX)](https://mastodon.social/%40metanetx)

# MNXtools

!!! This is work in progress - nothing below is yet completely true !!!

The code in this repository is primarily designed to build a dockerized application. The latest stable docker version can be obtained from dockerhub: https://hub.docker.com/r/sibswiss/mnxtools.   

see also: https://www.metanetx.org

## How to download MNXref cache
```
export MNXREF_VERSION=4.4;
cd cache/ && wget "https://www.metanetx.org/ftp/$MNXREF_VERSION/ChemSpace.bindump.gz" && gunzip ChemSpace.bindump.gz
```
