The Docker container can be built from the _Dockerfile_ file provided here, or pulled from https://hub.docker.com/r/sibswiss/mnxtools

Once you have the mnxtools Docker container locally you can run it like this:
```
docker run --name mnxtools --rm -i -t mnxtools:tagname <MNXtools command>
```
- *tagname* is the container version
- *--name* assignes a name to the running container
- *--rm* automatically removes the container when it exits
- *-i* opens an interactive session with the container
- *-t* allocates a pseudo-TTY



Available **MNXtools commands** are:
```
SBML_to_MetaNetX.pl  [options] -sbml <SBML file to convert to MetaNetX TSV format> -outdir <dir>
convert_mnet.pl      [options] <in-dir> <out-dir>

```
