# Plant Mitogenome Annotation Tools: Legal & Technical Assessment for the Taxodactyl Nextflow Pipeline

## TL;DR
- **No single "drop-in" tool is both accurate for plant mitogenomes AND cleanly redistributable.** The most accurate plant-specific tool, **PMGA**, is AGPL-3.0 in its own code but its distributable container bundles **MAKER (academic/non-commercial only)** and **Vmatch (proprietary)**, so you cannot legally redistribute it; and **MITOS/MITOS2**, though MIT-licensed and standalone, is a *metazoan* annotator that mishandles the plant genetic code and cannot model introns.
- **The one legally-safe, standalone, plant-capable annotator is MFannot (GPL-3.0, GitHub `BFL-lab/Mfannot`)** — a Perl/BLAST/Exonerate/HMMER/RNAweasel pipeline that excels at the hardest part (Group I/II intron-containing genes) but is weaker on RNA-editing-created start/stop codons.
- **Recommendation:** wrap **MFannot (GPL-3.0)** as the structural core and chain it with GPL/permissively-licensed components — **tRNAscan-SE + ARAGORN** (tRNAs), **nhmmer/HMMER against a curated angiosperm reference/HMM set** (PCG/rRNA location), and optionally **Deepred-mt (MIT)** for C-to-U RNA-editing sites — mirroring how GeSeq internally chains BLAT+HMMER+tRNAscan-SE+ARAGORN. Treat PMGA and GeSeq as web-only "gold-standard" comparators you point users to, not code you bundle.

## Key Findings

**Tool-by-tool license & standalone status**

| Tool | Standalone? | License | Plant mito–capable? | Redistributable in public Nextflow pipeline? |
|---|---|---|---|---|
| **MFannot** | Yes (Perl, Unix; GitHub) | **GPL-3.0** | Yes (intron-aware; organelle) | **Yes** (GPL obligations) |
| **MITOS / MITOS2** | Yes (Bioconda, PyPI, Docker) | **MIT** | **No** (metazoan; no intron models) | Yes legally, but not fit for plants |
| **PMGA** | Yes (Singularity container) | AGPL-3.0 (own code) | **Yes — best accuracy** | **No** — container bundles MAKER + Vmatch |
| **Mitofy (original)** | Yes (Perl/macOS; web) | **No license stated** | Yes (seed plants) | Risky — no clear license on original |
| **Mitofy (`dsenalik` fork)** | Yes (Perl/Linux CGI) | **GPL-3.0** | Yes (seed plants) | **Yes** (GPL), but legacy/unmaintained |
| **AGORA** | **No (web-only)** | No OSS license ("free to use") | Partial (BLAST homology) | No |
| **GeSeq** | **No (web-only)** | Web service ToS | Yes (has Mitochondrial mode) | No (cannot bundle) |
| **Chloë (Chloe.jl)** | Yes (Julia) | **No license file (all-rights-reserved)** | **No — plastid only** (mito = "Emma") | No |
| **Emma** | Yes (Julia) | **No license file (all-rights-reserved)** | **No — metazoan mito** | N/A for plants |
| **fpma** | Yes (Rust binary) | **No license file (all-rights-reserved)** | Partial (angiosperm PCG location only) | **No** until a license is added |
| **DOGMA** | Legacy web | Legacy | **No** (plastid + animal mito only) | No |

**Accuracy benchmark (from the PMGA paper, errors as [genes, start/end, splice-site], lower is better):** On *Arabidopsis thaliana* — PMGA (1, 2, 0) ≪ GeSeq (18, 0, 10) < MFannot (21, 12, 6) < AGORA (42, 5, 16). On *Liriodendron tulipifera* — PMGA (1, 2, 0), GeSeq (16, 1, 10), MFannot (28, 20, 13), AGORA (45, 18, 20). PMGA is clearly the most accurate; GeSeq second; both are the tools you cannot bundle.

## Details

### 1. Mitofy
Mitofy (Alverson, Wei, Rice, Stern, Barry & Palmer 2010, *Mol Biol Evol* 27:1436–1448, DOI 10.1093/molbev/msq029) is the original plant-mitogenome-specific annotator. It is a Perl program (originally macOS-only, bundling macOS builds of tRNAscan-SE and BLAST+) hosted at `https://dogma.ccbb.utexas.edu/mitofy/`, which as of 2026 still serves a working landing page with "Mitofy Webserver / Mitofy Download / Mitofy Documentation" links. It uses BLAST against a curated set of angiosperm mitochondrial protein/rRNA references plus tRNAscan-SE, and requires substantial manual curation. It remains widely cited in 2024–2025 plant mitogenome papers (e.g. *Populus simonii*, *Blechnopsis orientalis*), usually via the public web server.

**Licensing:** the original distribution does not carry a stated open-source license, which is a redistribution risk. However, community forks clarify the situation: **`dsenalik/mitofy`** ("Fork of mitofy for linux web interface," running at `vcru.wisc.edu/cgi-bin/mitofy/mitofy.cgi`) is explicitly **GPL-3.0 licensed**, and `jianzuoyi/mitofyX` / `maclandrol/mitofyX` are further Perl forks that simplify installation. The GPL-3.0 Linux fork is standalone-capable and legally redistributable, but the codebase is effectively unmaintained legacy software and does not natively handle RNA editing, multipartite structures, or batch/multi-contig input.

### 2. MITOS / MITOS2
MITOS (Bernt et al. 2013) and MITOS2 (Donath et al. 2019) are hosted at `gitlab.com/Bernt/MITOS`, packaged on **Bioconda** (`mitos`, current v2.1.10, updated Dec 9 2025), **PyPI** (`pip install mitos`), and as **biocontainers Docker/Singularity images** on quay.io. The Bioconda recipe and Anaconda page both state the license is **MIT** — fully permissive and redistributable.

**However, MITOS is fundamentally a metazoan tool.** Its own Bioconda/GitLab summary reads verbatim: "MITOS is a tool for the annotation of metazoan mitochondrial genomes." The MFannot review (Lang et al. 2023) notes MITOS2 is "also tailored to animals" and, decisively, **"MITOS2 cannot model introns."** Plant mitogenomes are intron-rich (Group I and II introns in *cox1, nad1, nad2, nad4, nad5, nad7*, etc.), use a different genetic code from animals, and undergo C-to-U RNA editing. Running MITOS2 on a plant mitogenome will therefore miss intron-containing genes and misannotate boundaries. **Conclusion: MIT-licensed and easy to wrap, but not fit for purpose on plant mitogenomes.** (The 2024 tool **DeGeCI 1.1**, Fiedler, Bernt & Middendorf, *Bioinformatics Advances* 4(1):vbae072, is a related command-line + web de-Bruijn-graph annotator, but is again oriented to metazoan reference mitogenomes.)

### 3. PMGA — Plant Mitochondrial Genome Annotator
PMGA (Li et al., *Plant Communications* 6, published March 2025; DOI 10.1016/j.xplc.2024.101191; PMC11956084; from Chang Liu's IMPLAD group) is the newest and, on the published benchmarks, **the most accurate** plant-mitogenome annotator. It was purpose-built to overcome the three failings of prior tools: (1) it annotates **multiple chromosomes/contigs simultaneously**; (2) it correctly handles **splice sites / exon boundaries**; and (3) it uniquely annotates **advanced features — RNA editing sites (start-gain/stop-gain codons) and mitochondrial-plastid DNA insertions (MTPTs)**. Web server: `http://www.1kmpg.cn/pmga/`. Code: **`github.com/implad-bix/PMGA`**, licensed **AGPL-3.0**, with the runnable source distributed as a **Singularity container on figshare (DOI 10.6084/m9.figshare.27201798)**.

**The redistribution blocker.** PMGA is a meta-pipeline that internally calls MAKER, CPGAVAS2, BLAST, ARAGORN, tRNAscan-SE, MISA, TRF, Vmatch, Deepred-mt, and OGDRAW. Two of these dependencies cannot be freely redistributed:
- **MAKER** (Yandell Lab, `github.com/Yandell-Lab/maker`) is dual-licensed Artistic-2.0/GPL **but only for academic/non-profit use**; its notice states MAKER "is not available for commercial use without a license" (obtainable via the University of Utah TCO). This non-commercial field-of-use restriction is incompatible with redistributing an open pipeline that any user, including commercial ones, can pull.
- **Vmatch** (`vmatch.de`, S. Kurtz, University of Hamburg) is **proprietary**, distributed as binaries under a per-user signed license agreement with **no third-party redistribution grant**.

By contrast, PMGA's other flagged dependencies are fine: **Deepred-mt** (`github.com/aedera/deepredmt`, Edera et al. 2021, *Computers in Biology and Medicine* 136:104682) is **MIT**, and **TRF** (`github.com/Benson-Genomics-Lab/TRF`) has been **AGPL-3.0** since June 2020. But MAKER or Vmatch alone is sufficient to make the bundled PMGA container **non-redistributable in a public Nextflow pipeline.** You may point users to PMGA's own web server or ask them to install it themselves, but you cannot ship its container.

### 4. AGORA
AGORA (Jung, Kim, Jeong & Yi 2018, *Bioinformatics* 34:2661–2663, DOI 10.1093/bioinformatics/bty196) is a **web-only** application at `https://bigdata.dongguk.edu/gene_project/AGORA/`, implemented in Python + PHP. It handles **both mitochondrial and plastid** genomes across eukaryotes (including a "plant" organism category and exon–intron/IR handling) via BLAST homology to user-selectable references — so the recollection that it is plastid-focused is only half right; it is genuinely organelle-general. **But** it is web-only with no standalone source distribution and no open-source license (merely "users can freely use the software"), and it was the **least accurate** of the four tools in the PMGA benchmark. **Not usable for bundling.**

### 5. Other plant-capable / general organelle tools

**GeSeq** (Tillich et al. 2017, *NAR* 45:W6–W11) is web-only (`chlorobox.mpimp-golm.mpg.de/geseq.html`). It **does support mitochondrial annotation** — the interface exposes a "Mitochondrial" sequence type and multiple mitochondrial genetic-code tables — and it was designed for plants (plastids in particular). Internally it chains **BLAT + HMMER + tRNAscan-SE + ARAGORN/ARWEN**, and can call **Chloë** and **MFannot** as third-party annotators. It performed second-best in the PMGA benchmark. The GeSeq paper explicitly notes that, unlike DOGMA, it can "take RNA editing into account." **However, GeSeq is not distributed as source or a container**, so it cannot be wrapped/redistributed — it is a web service only. It is the reference implementation to imitate architecturally, not to bundle.

**Chloë (Chloe.jl)** (`github.com/ian-small/Chloe.jl`, Ian Small lab, Julia; web at `chloe.plastid.org`) is an XGBoost/alignment-based annotator **optimised for angiosperm chloroplast genomes only**. Its own documentation directs mitochondrial users elsewhere ("To annotate mitochondria go to Emma"), so **Chloë does not handle plant mitochondria**. Note that, as of the August 2026 repository state, **Chloe.jl carries no LICENSE file** (the GitHub Resources panel lists only a Readme; "No description, website, or topics provided") — so it is effectively **all-rights-reserved**, not GPL. It is standalone (Julia package + `chloe_references` repo), but that is irrelevant for the mitochondrial use case and it should not be assumed redistributable.

**Emma** (`github.com/ian-small/Emma`, Julia) is the Small lab's mitochondrial counterpart to Chloë, but its README states it is a **"Metazoan Mitochondrial Annotator … optimised for vertebrates, especially fish,"** using vertebrate/invertebrate NCBI translation tables. **It is not a plant tool** and carries **no license file** (treat as all-rights-reserved). Not applicable.

**fpma — Fast Plant Mito Annotation** (`github.com/tolkit/fpma`, Max Brown / Wellcome Sanger 2022; help header "\<Max Brown; Wellcome Sanger 2022\> Fast plant mito annotation (fmpa). Version: 0.1.2") is a small **standalone Rust binary** that runs `nhmmer` (HMMER3) with a bundled set of angiosperm HMM profiles — "a set of known genes (43 core + 31 tRNA)" — to locate mitochondrial genes and emit GFF3. It is genuinely plant-mitochondrion-specific and trivially installable (`cargo install`), and is used inside the tolkit `plant-organellome-assembly` (PMPAP) pipeline alongside `fppa` (its plastid sibling). Its limitations: it only performs **HMM-based gene location** — no intron modeling, no tRNA/rRNA structural prediction, no RNA-editing awareness. **Critically, the repository (v0.1.2, 27 commits) contains no LICENSE file and displays no license in its GitHub panels — it is effectively all-rights-reserved. Do not bundle it until an explicit open-source license is added**; its HMM profiles are nonetheless a useful design reference.

**DOGMA** (Wyman, Jansen & Boore 2004) is a legacy web annotator for **plastid and animal mitochondrial** genomes; it explicitly does **not** annotate plant mitogenomes and is superseded. Not applicable.

### 6. Why plant mitogenome annotation is harder than plastid annotation
The technical hazards that break naïve annotators, and how the tools cope:
- **Multipartite / recombining structure.** Plant mitogenomes span an enormous size range — from ~66 kb in the parasitic *Viscum scurruloideum* to 11.3 Mb in *Silene conica* (Sloan et al. 2012, *PLOS Biology*), with still-larger genomes reported since (11.7 Mb *Larix sibirica*; ~18.99 Mb *Cathaya argyrophylla*) — and frequently map to multiple circular/linear chromosomes rather than a single "master circle" (e.g. multi-chromosomal *Populus simonii*, *Pulsatilla patens*). Tools that assume a single circular molecule (MITOS2, Mitofy) require manual stitching; **PMGA is the only tool that natively concatenates and annotates multiple contigs**.
- **RNA editing (C-to-U).** Plant mitochondrial mRNAs are heavily edited — from 313 sites in *Populus* (Edera et al. 2018) up to 902 sites in *Pulsatilla patens*, where *nad4* alone carried 65 C-to-U changes (Szandar et al. 2022, *BMC Plant Biology*, DOI 10.1186/s12870-022-03492-1) — and editing frequently **creates start (ACG→AUG) and stop codons**. Annotators keying on genomic start/stop codons therefore misplace gene boundaries. This is exactly where MFannot loses points (12 and 20 start/end errors in the benchmark). **PMGA models start-gain/stop-gain editing (using Deepred-mt)**; GeSeq claims partial RNA-editing awareness.
- **Genetic code.** Plant mitochondria use the **standard genetic code (translation table 1)**, unlike animal mitochondria (tables 2/5). Tools hard-wired to animal mito codes (MITOS2 default, Emma) will mistranslate plant genes — you must ensure the standard table is selected.
- **Promiscuous DNA (MTPT/NUPT/NUMT-like insertions).** Plant mitogenomes accumulate large plastid-derived (MTPT) and nuclear-derived insertions, often near-complete plastid inverted repeats (as in *Pulsatilla*). These generate spurious plastid-gene hits. **PMGA explicitly annotates MTPTs (via CPGAVAS2)**; most other tools flag these only if the user curates.
- **Introns.** Plant mito genes carry numerous Group I and especially Group II introns, including *trans*-spliced ones. **MFannot is the strongest open tool here** (its whole design centers on Group I/II intron and ncRNA modelling via RNAweasel/ERPIN/Infernal + Exonerate); **MITOS2 cannot model introns at all.**

## Recommendations

**Stage 1 — Ship MFannot as the legally-safe structural core (now).**
Wrap **MFannot** (GPL-3.0, `github.com/BFL-lab/Mfannot` + `MFannot_data`, also GPL-3.0). It is standalone (Perl + BLAST + Exonerate + HMMER + RNAweasel), organelle-capable, and the strongest open tool for intron-containing plant mito genes. GPL-3.0 permits redistribution provided your pipeline complies with GPL obligations (make source available, preserve license notices; the pipeline itself must remain GPL-compatible). Before shipping, verify the license of MFannot's bundled RNAweasel/ERPIN component, since ERPIN's own terms are not GPL-obvious; the Infernal-based path (Infernal is BSD-licensed) is a safer substitute where available.

**Stage 2 — Chain complementary open components to match GeSeq's architecture.**
Because no single open tool matches PMGA's accuracy, replicate GeSeq's chaining strategy with only redistributable parts:
- **tRNAs:** tRNAscan-SE (GPLv3) + ARAGORN — both bundleable.
- **rRNAs / PCG location:** nhmmer/HMMER (BSD-licensed) against a curated angiosperm mitochondrial HMM/reference set. (`fpma`'s 43-core-gene + 31-tRNA HMM set is a good design template, but you must build/curate your own set or wait for `fpma` to add a license, since `fpma` itself is currently unlicensed and cannot be bundled.)
- **RNA-editing sites:** **Deepred-mt (MIT)** — freely bundleable — to predict C-to-U edits and rescue edit-created start/stop codons, the single biggest MFannot weakness.
- **MTPT/plastid-insertion masking:** a BLAST screen against a plastid reference before finalizing calls.

**Stage 3 — Offer PMGA and GeSeq as opt-in external services, not bundled code.**
Since PMGA is the most accurate tool but its container is non-redistributable (MAKER + Vmatch), and GeSeq is web-only, expose them as **optional external steps**: let users who have independently installed PMGA (or accept its web ToS) point Taxodactyl at it, and/or emit GeSeq/PMGA-web-submission-ready FASTA. This gives users best-in-class accuracy without you redistributing restricted code.

**Benchmarks / thresholds that would change this recommendation:**
- If the PMGA authors publish a **container with MAKER and Vmatch removed or swapped** for open equivalents (e.g. AUGUSTUS/miniprot instead of MAKER; a free repeat finder instead of Vmatch), PMGA becomes the recommended default — it is materially more accurate than everything else.
- If **`fpma` adds a permissive license** (MIT/BSD/Apache), promote it from design template to a bundled first-pass PCG locator.
- If a new 2025–2026 tool appears that is **both plant-mito-accurate and permissively licensed end-to-end**, prefer it; monitor the tolkit, Bernt/MITOS, and Small-lab repos.

## Caveats
- **PMGA licensing conclusion depends on the container contents.** The AGPL-3.0 on PMGA's own code is not the problem; the blockers are the bundled MAKER (academic/non-commercial only) and Vmatch (proprietary). If PMGA's figshare container ships only *wrappers* that call user-installed copies of MAKER/Vmatch (rather than the binaries themselves), the redistribution picture could differ — inspect the container before finalizing. MAKER's exact clause should be read from `github.com/Yandell-Lab/maker/blob/master/LICENSE`; the subagent verifying it could reach the reproduced license text but not the raw file body in-session.
- **Several repos lack explicit licenses** and must be treated as all-rights-reserved (not bundleable) until clarified: **Chloe.jl**, **Emma**, and **`fpma`** all show no LICENSE file as of August 2026. The original Mitofy also lacks a clear license; rely on the GPL-3.0 `dsenalik` fork instead.
- **GeSeq and AGORA are web services**; availability and terms can change, and neither offers redistributable source. Do not build hard pipeline dependencies on them.
- **Accuracy numbers** come from the PMGA authors' own benchmark on two genomes (*A. thaliana*, *L. tulipifera*), so there is a possible home-tool advantage; treat the relative ranking (PMGA > GeSeq > MFannot > AGORA) as indicative rather than definitive.
- MITOS2 remains excellent for **animal** mitogenomes and is MIT-licensed; the recommendation against it is strictly about plant unsuitability, not licensing.
- **TRF's** open-source AGPL-3.0 status applies to the current (2020+) GitHub release; pre-2020 TRF binaries carried a restrictive no-redistribution academic license, so confirm which TRF version any bundled component uses.
