# Can GeSeq Be Run Offline? Licensing Status and Open-Source Alternatives for a Public Chloroplast-Annotation Pipeline

> Not part of the spec, but useful for consideration when deciding how to
> build the plastid annotation workflow step.

## TL;DR
- **No public offline/standalone/dockerized/containerized GeSeq exists, and its source code is not publicly available.** GeSeq is a closed-source Java+PHP *web application* hosted only at chlorobox.mpimp-golm.mpg.de; there is no GitHub/GitLab repo, no Bioconda recipe, and no Docker/BioContainers image for GeSeq itself.
- **GeSeq is effectively proprietary ("all rights reserved" with an academic-use carve-out).** The CHLOROBOX Disclaimer states academic users may use it non-commercially, but that "the download of any component of the CHLOROBOX … requires an agreement" with Max-Planck-Innovation GmbH, and that any distribution needs the authors' prior written consent. It therefore cannot be legally redistributed, rewrapped, or reimplemented into a public pipeline without written permission, and it exposes no official API for automated submission.
- **For a public Nextflow pipeline, use the genuinely open-licensed alternatives instead: PGA (GPL-3.0) and Plann (BSD-3-Clause) are standalone, command-line, batch-capable plastome annotators**, and the open building blocks GeSeq itself wraps (BLAT, HMMER/nhmmer, Infernal, tRNAscan-SE, ARAGORN/ARWEN) are all separately available. Note that Chloë — the annotator GeSeq now embeds — is on public GitHub but carries **no license file (all rights reserved)**, so it is not redistributable either.

## Key Findings

1. **No offline GeSeq.** Exhaustive searching of GitHub, Bioconda, Docker Hub, Biostars/SEQanswers, and the literature found no downloadable, standalone, or containerized build of GeSeq. Every reference points back to the web server. In its own paper (Tillich et al. 2017, *Nucleic Acids Research* 45(W1):W6–W11, doi:10.1093/nar/gkx391), the tool is described as a **"web application … written in Java and PHP"** with "parallel (multi-threaded) job processing and user management" — i.e., a server application, not a distributable command-line tool.

2. **No public source code.** There is no GeSeq repository under the Max Planck Institute of Molecular Plant Physiology, under the authors (Michael Tillich, Stephan Greiner, Pascal Lehwark), or anywhere else. (A GitHub user handle "geseq" exists but belongs to an unrelated London software developer.)

3. **Proprietary licensing.** The CHLOROBOX Disclaimer §2 ("Limitation of use") states verbatim: *"Academic users may download the material offered on the site for their non-commercial use, but all copyright and other proprietary notices contained in the materials are to be retained. Non-academic/commercial, for-profit users may use the CHLOROBOX website and online services, but any other use, in particular the download of any component of the CHLOROBOX, requires an agreement. Please contact the Max-Planck-Innovation GmbH for the latter."* §5 adds that any duplication or distribution "beyond the scope of copyright law shall require the prior written consent of the author." The GeSeq *paper* is CC BY-NC 4.0, but that licenses the article text, not the software. Bottom line: GeSeq is not open source and cannot be legally bundled into a public pipeline.

4. **No API / anti-automation terms.** GeSeq supports "batch submission" only within its browser form; there is no documented REST API or programmatic endpoint. The only known automation approach is scraping the web form (as the CroP-BioDiv "irs_wrappers" project does — it notes GeSeq "sequence has to be processed with Web application and wrapper extracts IRs location from result GenBank file"). Given §2 and §7 (no commercial use of contact info) of the Disclaimer, automated/high-volume hammering of the server is legally and ethically problematic and against the spirit of "Citations keep this server running."

5. **What GeSeq wraps internally.** Per the paper and the "3rd Party Software" page: BLAT (standalone BLAT v.35), HMMER/nhmmer, Infernal, tRNAscan-SE (v1.3.1 and v2.0.7), ARAGORN (v1.2.38), ARWEN (v1.2.3), MUSCLE, TranslatorX, and OGDRAW for visualization. The current web form additionally offers two "3rd Party Stand-Alone Annotators": **Chloë v0.1.0** and **MFannot v1.34**. Most of these components are independently open source and could be assembled into an open pipeline.

## Details

### GeSeq itself
GeSeq (Tillich M, Lehwark P, Pellizzer T, Ulbricht-Jones ES, Fischer A, Bock R, Greiner S. 2017. *Nucleic Acids Research* 45:W6–W11; doi:10.1093/nar/gkx391) is a web application only. Its abstract states: *"We have developed the web application GeSeq (https://chlorobox.mpimp-golm.mpg.de/geseq.html) for the rapid and accurate annotation of organellar genome sequences, in particular chloroplast genomes."* Its annotation strategy is BLAT-based homology search against a curated/uploaded reference set (translated BLATx for CDS, BLATn for RNA/DNA features), plus optional profile-HMM (nhmmer) searches for CDS/rRNA, and de novo tRNA prediction. The output is GenBank/GFF3/JSON, visualized via OGDRAW. Recent papers cite "GeSeq v2.0.3," indicating the server is actively maintained, but there is no versioned downloadable release. GeSeq was created to fill a gap among the then-available tools — the paper notes Mitofy annotates only plant mitochondria, CpGAVAS and Verdant only chloroplasts, and DOGMA chloroplast and animal (but not plant) mitochondria.

The "anecdote" that GeSeq can be run offline is most likely a confusion with (a) the standalone third-party tools it wraps (BLAT, tRNAscan-SE, ARAGORN, HMMER, Chloë, MFannot — all of which do run locally), or (b) other command-line plastome annotators such as PGA and Plann. I found no evidence that a standalone GeSeq binary or container has ever been released to the public.

### Legal analysis for a public pipeline
- **Redistributing/reimplementing GeSeq:** Not permitted without a written agreement from Max-Planck-Innovation GmbH. A public Nextflow module that bundled GeSeq or its curated reference database would infringe.
- **Calling the public GeSeq web server from a pipeline:** No API exists; scraping the form for automated/high-volume use conflicts with the Disclaimer's use-limitation clauses. Not recommended for a published, reusable pipeline.
- **Practical conclusion:** Build the pipeline from openly licensed components.

### Recommended open-source, wrappable annotators

| Tool | Language | License | Automation fit |
|---|---|---|---|
| **PGA (Plastid Genome Annotator)** | Perl (needs BLAST+ 2.8.1+, Perl 5) | **GPL-3.0** | Standalone CLI, explicit batch annotation, reference-based; ideal Nextflow process |
| **Plann** | Perl (needs BLASTN + tbl2asn) | **BSD-3-Clause** | Locally executable CLI, explicitly designed for pipelines and GenBank submission |
| **Chloë** | Julia | **No license (all rights reserved)** | Runs locally and has a web API, but not redistributable; angiosperm-optimized |
| **GetOrganelle** | Python | **GPL-3.0** | Assembly (not annotation) — the upstream step; on Bioconda |
| **AGORA** | Python/PHP web app | Web app; "users can freely use" (no standalone OSS license found) | Web only |
| **CPGAVAS2** | Web server | Web only | Web only |

- **PGA** (Qu X-J, Moore MJ, Li D-Z, Yi T-S. 2019. *Plant Methods* 15:50, doi:10.1186/s13007-019-0435-7) is explicitly released under **GPL-3**; its Availability section states verbatim *"License: GPL-3"* and, for restrictions to use by non-academics, *"none."* Its abstract describes it as using *"reference plastomes as the query and unannotated target plastomes as the subject to locate genes, which we refer to as the reverse query-subject BLAST search approach. PGA accurately identifies gene and intron boundaries as well as intron loss."* GPL-3 permits building and publishing a public pipeline that invokes PGA (and even redistributing PGA), provided GPL terms are honored. Because Nextflow orchestrates PGA as a separate executable rather than linking to it, your pipeline code itself need not be GPL. A v2.0 (PGA2) is also available.
- **Plann** (Huang DI, Cronk QCB. 2015. *Applications in Plant Sciences* 3:1500026, doi:10.3732/apps.1500026) is **BSD-3-Clause** — the most permissive option. The LICENSE file (github.com/daisieh/plann) reads *"Copyright (c) 2015, Daisie Huang. All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met…"* (a standard 3-clause BSD text). Its paper's Conclusions state: *"Unlike Web-based annotation packages, Plann is a locally executable script that will accurately annotate a plastome sequence to a locally specified reference plastome. Because it executes from the command line, it is ready to use in other software pipelines and can be easily rerun as a draft plastome is improved."* It uses BLASTN and NCBI tbl2asn.
- **Chloë** (Ian Small lab, UWA; chloe.plastid.org; github.com/ian-small/Chloe.jl) is the modern alignment/XGBoost-based annotator GeSeq now embeds. It runs locally via Julia and offers a web API (e.g., `https://chloe.plastid.org/annotate-ncbi?ncid={accession}&force_circular=true`), but the repository contains **no LICENSE file and no license field in Project.toml** (confirmed: the GitHub sidebar shows no license entry, and authors are listed as Ian Small and Ian Castleden, version 0.1.18). By default this makes it "all rights reserved." You can run it, but you cannot legally redistribute/repackage it into a public pipeline without the authors' permission, and its README gates non-angiosperm use behind emailing the author (ian.small@uwa.edu.au).
- **GetOrganelle** (Jin et al. 2020, *Genome Biology* 21:241; GPL-3, on Bioconda via `conda install -c bioconda getorganelle`) and **NOVOPlasty** are the standard upstream *assembly* tools. A full open pipeline would typically chain GetOrganelle (assembly) → reorientation → PGA or Plann (annotation) → OGDRAW/visualization. A recent worked example is **CGAS** (Chloroplast Genome Analysis Suite; Abdullah et al. 2026, *iMetaOmics*), an **MIT-licensed** Python pipeline that already wraps fastp, GetOrganelle, and PGA end-to-end.

## Recommendations

1. **Do not distribute, rewrap, reimplement, or auto-submit to GeSeq in a public pipeline.** It is proprietary and has no API. If you specifically need GeSeq's output at scale, email Max-Planck-Innovation GmbH / the CHLOROBOX team to request a licensing agreement or a server-side batch arrangement — do not build a public tool around scraping the form.

2. **Build the public Nextflow pipeline on PGA (GPL-3.0) as the primary annotator, with Plann (BSD-3-Clause) as an alternative/complementary path.** Both are standalone, batch-capable, and unambiguously redistributable. Package each in its own container (Perl + BLAST+/tbl2asn) and pin versions. If you want a fully permissive license posture, prefer Plann and keep PGA as a separate opt-in module.

3. **Chain the full workflow from open components** to reproduce GeSeq-like functionality: GetOrganelle (assembly, GPL-3) → reorient to LSC-IR-SSC-IR → PGA/Plann (annotation) → optional tRNAscan-SE/ARAGORN for tRNAs → OGDRAW/visualization. This mirrors what GeSeq does internally using the same open building blocks. Consider studying the MIT-licensed CGAS pipeline as a reference implementation.

4. **If you want Chloë's accuracy specifically,** contact Ian Small (ian.small@uwa.edu.au) to obtain an explicit license before including it in a redistributable pipeline. Until then, treat it as a non-redistributable optional dependency the user installs themselves, or call its web API for low-volume use within its terms.

**Thresholds that would change this advice:** If MPI-MP later publishes GeSeq under an OSI-approved license or a Bioconda/BioContainers package appears, it could be wrapped directly. If Chloë adds a permissive LICENSE file, it becomes the best drop-in for GeSeq-style accuracy. If GeSeq exposes a documented API with terms permitting automation, remote invocation becomes viable.

## Caveats
- GeSeq is actively developed (the web form shows Chloë v0.1.0, MFannot v1.34, tRNAscan-SE v2.0.7, ARAGORN v1.2.38); the Disclaimer is dated "Chlorobox 2021," so licensing terms could change — verify current terms before relying on them.
- "No license" for Chloë is confirmed by the absence of a LICENSE file and license field on the ian-small/Chloe.jl GitHub repo; the authors may grant permission on request, and the `ian-small/chloe` clone URL resolves to the same repository.
- GPL-3 compliance for PGA/GetOrganelle requires honoring copyleft obligations for any *distributed modifications*; simply invoking them as separate executables from Nextflow is fine, but forking/modifying and redistributing them carries GPL obligations. If your own pipeline must be permissively licensed, keep GPL tools as separately-invoked containers, not as vendored/modified code.
- I could not find any official GeSeq API; the absence of documentation is itself the finding, but I cannot prove no private/undocumented endpoint exists.
- AGORA and CPGAVAS2 are web-only and were not confirmed to have standalone open-source distributions, so they are not suitable for embedding in a public offline pipeline.
