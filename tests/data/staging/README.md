# Fixture staging

Landing zone for raw SRA / ENA downloads that feed the subsampled
fixtures in [`../`](../). Contents are git-ignored (see
[`../../../.gitignore`](../../../.gitignore)) except this README.

Fixtures should be regenerated from these sources using the recipes in
[`../README.md`](../README.md).

## Currently staged

### `SRR8306868_1.fastq.gz`

Feeds `animal.fastq.gz` and `animal_b.fastq.gz`.

| Field | Value |
|---|---|
| Organism | *Acyrthosiphon pisum* (Pea aphid) |
| Instrument platform | OXFORD_NANOPORE |
| Instrument model | MinION |
| Read count | 892,247 |
| Base count | 7,338,953,827 |
| Center name | SUB4910096 |
| Library layout | SINGLE |
| Library strategy | WGS |
| Library source | GENOMIC |
| Library name | Nanopore |
| Library selection | unspecified |
| ENA first public | 2019-05-02 |
| ENA last update | 2019-05-02 |

Download from ENA:

```
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR830/008/SRR8306868/SRR8306868_1.fastq.gz
```


### `SRR11315861_1.fastq.gz`

Feeds `plant.fastq.gz`.

| Field | Value |
|---|---|
| Organism | *Datura stramonium* (Jimsonweed) |
| Instrument platform | OXFORD_NANOPORE |
| Instrument model | MinION |
| Read count | 558,821 |
| Base count | 8,446,989,155 |
| Center name | SUB7154150 |
| Library layout | SINGLE |
| Library strategy | WGS |
| Library source | GENOMIC |
| Library name | Nanopore |
| Library selection | RANDOM |
| ENA first public | 2021-04-02 |
| ENA last update | 2021-04-02 |

Download from ENA:

```
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR113/061/SRR11315861/SRR11315861_1.fastq.gz
```
