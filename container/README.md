The Docker container can be built from the _Dockerfile_ file provided here, or pulled from https://hub.docker.com/r/sibswiss/mnxtools

E.g.
```bash
docker pull sibswiss/mnxtools
```

Once you have the mnxtools Docker container locally you can run it like this:
```bash
export BINDING_PATH=/tmp
docker run --name mnxtools --mount type=bind,source=$BINDING_PATH,target=/mybinding --rm -i -t sibswiss/mnxtools <MNXtools command>
```
- *BINDING_PATH* is the outside world path mounted in the container
- */mybinding* is the outside world path within the container
- *--name* assignes a name to the running container
- *--rm* automatically removes the container when it exits
- *-i* opens an interactive session with the container
- *-t* allocates a pseudo-TTY



Available **MNXtools commands** are:
```
SBML_to_MetaNetX.pl  [options] -sbml <SBML file to convert to MetaNetX TSV format> -outdir <dir>
convert_mnet.pl      [options] <in-dir> <out-dir>

```
