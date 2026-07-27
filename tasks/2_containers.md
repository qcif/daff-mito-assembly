# Building required containers

We need to build a mulled container for the RECRUIT process:

```sh
pipx install galaxy-tool-util
mulled-build build 'minimap2=2.31,seqtk=1.5,samtools=1.24'
docker tag \
    quay.io/biocontainers/mulled-v2-9278ceb357570ba6e25522f5a16e2a4d3ba61a68:74b9019c6ae38f81f95dd09e981151bb2d1028ee-0 \
    neoformit/daff-wf5-recruit:latest
docker push neoformit/daff-wf5-recruit:latest
```
