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

## COVERAGE_GATE (C2)

Originally `python=3.12,seqkit=2.13,seqtk=1.5`
([task 15](15_coverage_gate.md) §3). Task 25 added `minimap2` so the
gate can split the recruited pool by organelle panel before estimating
coverage ([spec §2.1.5](../../spec/02-stages.md#215-sibling-organelle-carry-over-in-the-estimate)):

```sh
mulled-build build 'python=3.12,seqkit=2.13,seqtk=1.5,minimap2=2.31'
docker tag \
    quay.io/biocontainers/mulled-v2-b98d675e1edf2fea65e6536a11db11ff302b7b11:552d068ed47a19fe375b369440dc98b7008e2fec-0 \
    neoformit/daff-wf5-coverage-gate:python3.12_seqkit2.13_seqtk1.5_minimap2-2.31
docker push neoformit/daff-wf5-coverage-gate:python3.12_seqkit2.13_seqtk1.5_minimap2-2.31
```

Pinned in [`conf/containers.config`](../../conf/containers.config) as
`@sha256:bbfec02eb9670d3bf8a0209bfd00a595f45fce9f3ff7189f7e38f5057cab07e8`.

Note: `mulled-build` drops an `involucro` binary into the working
directory — delete it, it is not a project file.
