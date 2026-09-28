# NTM *Mycolicibacterium fortuitum* complex genomics pipeline

Parameterized Bash scripts used for the genomic characterization of the
*Mycolicibacterium fortuitum* complex from post-surgical non-tuberculous
mycobacteria (NTM) isolates in Ecuador.

> **Associated publication.** [Ángel Sebastián Rodríguez-Pazmiño1, Henry Parra-Vera2, Greta Franco-Sotomayor3, Miguel Ángel García-Bereguiain1]. [2026]. [Genomic characterization of the Mycolicibacterium fortuitum complex from post-surgical infections in Ecuador reveals a clonal Mycolicibacterium houstonense/farcinogenes cluster]. *Microbiology Spectrum* [[volume/DOI — update on acceptance]].
> **Raw reads.** NCBI BioProject **PRJNA1524051** (per-isolate BioSample/SRA accessions in the paper's Supplementary Table).
> **Genome assemblies.** [[accessions — to be added upon deposition]].

These scripts implement the analysis exactly as described in the Methods. They
are analysis wrappers around established, third-party tools (cited in the paper);
they are **not** a stand-alone software package. Each script is self-documenting:
run it with `-h` for usage.

## Pipeline overview and run order

| Step | Script | Purpose |
|---|---|---|
| 0 | `build_type_panel.sh` | Build the anchored *M. fortuitum* complex **type-strain ANI panel** (downloaded by accession) shared by the downstream steps. |
| 1 | `ntm_diagnose.sh` | Read QC (fastp) → de novo assembly (Shovill/SPAdes) → size/GC diagnostics → fastANI vs the panel → verdict (mycobacterium / *Neisseria* / mixed). Triage of contaminated cultures. |
| 2 | `checkm2_ani.sh` | Genome quality (**CheckM2**) + all-vs-all fastANI with single-linkage clustering to flag **clonal** groups. |
| 3 | `ntm_type_resistome.sh` | Fine species assignment (fastANI vs type strains) + **resistome** (Bakta + AMRFinderPlus/Abricate + species-matched targeted BLAST). Use `-H` to supply a *houstonense* allele source and remove the fortuitum-centric bias. |
| 4 | `ntm_profiler_run.sh` | Orthogonal NTM species identification with **NTM-Profiler** (cross-check to ANI/dDDH). |
| 5 | `snippy_snp.sh` | **SNP** analysis of a clonal cluster (Snippy → snippy-core → snp-dists → IQ-TREE). |
| 6 | `verify_farci_houst.sh` | Taxonomy verification: *M. houstonense* vs *M. farcinogenes* conspecificity and reference-assembly integrity check. |
| 7 | `context_tree.sh` | Global **core-genome phylogeny** (NCBI download → QC/dereplication → Prokka → Panaroo → IQ-TREE), with parallel annotation and an outgroup. |

Digital DNA–DNA hybridization (dDDH) was obtained from the web service TYGS
(https://tygs.dsmz.de) and is not scripted here.

## Requirements

All dependencies are installable from [Bioconda](https://bioconda.github.io/).
See `environment.yml` (create with `conda env create -f environment.yml`).
**Software versions used in the paper must be pinned here before archiving** —
replace the `[[version]]` placeholders in `environment.yml` with the exact
versions reported in the Methods.

Some steps also require downloadable databases (Bakta DB, CheckM2 DB, Kraken2 DB,
NTM-Profiler DB); see each tool's documentation.

## Usage

Every script prints full help with `-h`, e.g.:

```bash
./scripts/build_type_panel.sh -o panel
./scripts/ntm_diagnose.sh -i reads/ -o diag -R panel -c -t 16
./scripts/context_tree.sh -g genomes/ -o tree -s RefSeq -D 0.001 \
    -O GCF_000015005.1 -t 16 -j 8
```

## Reproducibility notes

- The type-strain panel is built by **accession** (not by taxon) to avoid a
  mislabeled public assembly (GCF_000723385.1, deposited as *M. farcinogenes*
  but genomically *M. senegalense*); the correct *M. farcinogenes* type-strain
  assembly is GCA_025821245.
- Scripts are idempotent where practical (they skip completed steps), so
  interrupted runs can be resumed.

## AI-use disclosure

In accordance with the ASM Generative AI Policy, the authors disclose that a generative AI assistant (Anthropic Claude) was used to help draft, review, and debug the analysis scripts in this repository. All scripts were designed, executed, inspected, and validated by the authors, who take full responsibility for the code and for the results reported in the associated publication.

## License

Released under the MIT License (see `LICENSE`).

## Citation

If you use these scripts, please cite the associated publication (above) and this
archived repository (Zenodo DOI: [[to be added on deposit]]). See `CITATION.cff`.
