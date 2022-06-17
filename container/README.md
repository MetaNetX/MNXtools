The Docker container can be built from the _Dockerfile_ file provided here, or pulled from https://hub.docker.com/r/sibswiss/mnxtools

E.g.
```bash
docker pull sibswiss/mnxtools
```

Once you have the mnxtools Docker container locally you can run it like this:
```bash
export BINDING_PATH=/tmp
docker run --mount type=bind,source=$BINDING_PATH,target=/mybinding --rm -i -t sibswiss/mnxtools
```
- *BINDING_PATH* is the outside world path mounted in the container
- */mybinding* is the outside world path within the container
- *--rm* automatically removes the container when it exits
- *-i* opens an interactive session with the container
- *-t* allocates a pseudo-TTY



```help``` lists the available **MNXtools commands**:
```bash
SBML_to_TSV.pl    [options] -sbml <SBML file to convert to MetaNetX TSV format> -outdir <dir>
map_mnet.pl       [options] <indir> <outdir>

```
