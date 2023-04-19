# covid-analysis
Sequencing and coronavirus strain identification of ONT samples.

This workflow begins with metadata and ONT fastqs. The reads are filtered using [artic](artic.readthedocs.io) guppyplex (for length), [kraken2](https://github.com/DerrickWood/kraken2/wiki/Manual) (to dehumanize), and [porechop](https://github.com/rrwick/Porechop) (to remove adapter/primer sequences). For more on filtration, see [Read Filtering](#read-filtering). Artic minion sequences the reads and looks for variants. Homopolish cleans up the consensus sequences which are analyzed by [pangolin](https://cov-lineages.org/resources/pangolin.html) and [nextclade](https://docs.nextstrain.org/projects/nextclade/en/latest/) to determine the strain associated with the sample. All of this is summed up in some CSV files at the end.

## Dependencies
User provided:
* [conda 3](https://docs.conda.io/en/latest/miniconda.html)

The following conda environments are required. If running via [./run_pipeline.sh](./run_pipeline.sh), these environments will be created automatically, but they can also be created by running [./prepare_envs.sh](prepare_envs.sh)
* env-nextflow (current version: `22.10.6.5843`, version must use DSL 2)
* env-artic
* env-kraken2
* env-porechop
* env-homopolish
* env-pangolin
* env-nextclade

## Installation
Clone this repository and create the conda environments:
```
git clone https://github.com/enviro-lab/covid-analysis
cd covid-analysis
# make sure the scripts are executable:
chmod +x *.sh *.py *.nf
# to install all conda environments, run
./prepare_envs.sh
 ## NOTE: ./run_pipeline.sh will attempt to do this ^^ for you, if needed
```

You'll need a kraken database which will be downloaded automatically if it doesn't already exist.
* If you already have a human db or want to set where it should be automatically installed, edit the `kraken_db` path in [./nextflow.config](nextflow.config).
* Setting up the human database requires [Kraken 2](https://github.com/DerrickWood/kraken2/wiki/Manual) to be installed.
    * (NOTE: This env should already exist from when you ran [./prepare_envs.sh](prepare_envs.sh))
* Use `kraken2` to download the human database.
  * This is done automatically the first time you run the workflow, but...
  * The download can also be done manually via:
    ```
    kraken_db=/path/to/db    # edit as needed
    threads=16               # edit as needed
    kraken2-build --db $kraken_db --download-taxonomy
    kraken2-build --db $kraken_db --download-library human
    kraken2-build --db $kraken_db --build --threads $threads
    kraken2-build --db $kraken_db --clean --threads $threads
    ```
* Note: if the --download library step fails with a message like: `rsync_from_ncbi.pl: unexpected FTP path (new server?) for https://ftp.ncbi.nlm.nih.gov/`, you can edit line 46 of `rsync_from_ncbi.pl` (which in my case appears at `./conda/env-kraken2/share/kraken2-2.1.2-3/libexec/rsync_from_ncbi.pl`).
  * change `^ftp` to `^https` (line 46)

Check out the `reportMap` variable in the params section of [./nextflow.config](nextflow.config). This determines how many reports and what kinds are produced at the end of the pipeline. Make sure the best option for you is the only one that's uncommented (comment out the other examples).

## Usage
This guide assumes you have a file structure like this, but the most important part is that you have a directory like `fastq_pass` containing individual files
```
Clinical-12-22-22-1-V2A-fastqs      # PLATE_DIR
├── Clinical-12-22-22-1-V2A.csv     # Contains sample metadata
└── fastq_pass                      # Contains ONT fastqs for each barcode
    ├── barcode01                   # Contains fastq(s) for barcode01
    |   ├── PAG76195_pass_barcode01_b3089c27_1.fastq.gz
    |   ├── PAG76195_pass_barcode01_b3089c27_2.fastq.gz
    |   └── PAG76195_pass_barcode01_b3089c27_n.fastq.gz
    ├── barcode02                   # Contains fastq(s) for barcode02
    └── barcodeNN                   # Contains fastq(s) for barcodeNN
```
To run, activate the nextflow environment and run analyzeReads.nf with a few necessary args like below.
```
# go to where the scripts are
cd covid-analysis

# using our wrapper script, run the pipeline
./run_pipeline.sh \
    --plate 12-22-22-1-V2A \
    --fastqs ${PLATE_DIR}/fastq_pass \
    --meta $PLATE_DIR/Clinical-12-22-22-1-V2A.csv \
    --out ${PLATE_DIR}/output
# with `--group` specified, this script will give your group modification access to the outputs once the workflow finishes

# or

# activate the env the run the pipeline
cd covid-analysis
conda activate conda/env-nextflow
nextflow ./analyzeReads.nf \
    --plate 12-22-22-1-V2A \
    --fastqs ${PLATE_DIR}/fastq_pass \
    --meta $PLATE_DIR/Clinical-12-22-22-1-V2A.csv \
    --out ${PLATE_DIR}/output
```

## Output
Once complete, your directory will look like this with a work directory used by nextflow and an output directory
```
Clinical-12-22-22-1-V2A-fastqs
├── Clinical-12-22-22-1-V2A.csv
├── fastq_pass
├── output                          # Contains all output
│   ├── kraken_trimmed              # Contains reads that have been trimmed and filtered of human content
│   ├── porechop_kraken_trimmed     # Contains reads that have been trimmed, dehumanized, and filtered by porechop
│   ├── raw_read_counts             # Lists sample_name and number of reads in original dataset
│   ├── read_counts                 # Lists sample_name and number of reads used in final analysis
│   ├── samples                     # Contains `artic minion` output
│   ├── nextclade                   # Contains extra `nextclade` output
│   └── overall                     # Contains csv output
|       ├── Sequencing-report-12-22-22-1-V2A-All.csv        # Combined report of metadata, strains, read counts
|       ├── Sequencing-report-12-22-22-1-V2A-Lab1.csv       # Same report but only where `Source Lab = Lab1`
|       ├── Sequencing-report-12-22-22-1-V2A-Lab2.csv       # Same report but only where `Source Lab = Lab2`
|       ├── consensus_12-22-22-1-V2A.fasta                  # Final consensus seqeuences for all samples
|       ├── lineage_report-12-22-22-1-V2A.csv               # Pangolin lineage report
|       └── nextclade_results-12-22-22-1-V2A.csv            # Nextclade lineage/variant report
└── work                            # Contains all intermediate pipeline outputs
```

## Inputs (command line arguments)
* `--plate`: A name used in the seqeuncing report filenames (preferrably unique)
* `--meta`: A csv containing metadata for each sample. All fields will be present in the sequencing report.
  * Required fields:
    * 'Seq ID':
      * Must be unique
      * This will appear as `Sample #` in the sequencing report
    * 'Primer Scheme'
      * used to decide `artic` inputs
    * 'Barcode'
     * maps sample to associated barcode directory
  * Recommended fields:
    * "Zip code", "Test date", "Sequence date", 
    * anything else
  * The following fields wil be renamed as shown (if present)
    * "Seq ID" --> "Sample #"
    * "Ct N1" --> "Ct N gene"
    * "Qubit reading" --> "post-ARTIC Qubit (ng/ul)"
    * "Zip Code" --> "Zip code"
    * "Sample date" --> "Test date"
* `--fastqs`: A directory containing individual barcode directories, 
  * each barcode directory contains 1 or more fastq file(s)
  * only the names mentioned in the 'Barcode' field of `--meta` file will be used
* `--kraken_out`: Path to your human kraken database
  * This can also be input as the environment variable HUMAN_KRAKEN_DATABASE
  * for installation steps, see [Read filtering](#Read-filtering) below
* `--out`: Where the output should go
  * it will be formatted as shown in the [output](#Output) section above

## Inputs (behind-the-scenes)

### Primer scheme info
In order to determine which primer scheme, file, and associated read lengths to use, that information is gathered from [./scheme_details.csv](./scheme_details.csv).

Summary of [./scheme_details.csv](./scheme_details.csv) fields:
* primer_scheme: how the scheme is labeled in the 'Primer Scheme' column of the `--meta` file
* scheme: the `scheme` value passed to `artic minion` (name of primer scheme)
* scheme_version: the `--scheme-version` value passed to `artic minion`
* scheme_dir: the `--scheme-directory` value passed to `artic minion`
* max: the `--max-length` value passed to `artic guppyplex`
* min: the `--min-length` value passed to `artic guppyplex`

To ensure artic can access any primer schemes you add, make sure the directory structure looks like this:
```
<scheme_dir>
├── <scheme>   # Contains 1 or more scheme version directories
│   └── <scheme_version>   # Contains scheme files
│       ├── <scheme>.primer.bed
│       └── <scheme>.reference.fasta
```
For more details about *.primer.bed files, see [artic's Primer Scheme page](https://artic.readthedocs.io/en/latest/primer-schemes/).

## Read filtering
This pipeline filters out reads in three steps.
1. `artic guppyplex` removes reads based on provided maximum and minimum values. These values are found in [./scheme_details.csv](./scheme_details.csv) and can be edited there, as desired.
2. `kraken2` filters out all human-like reads, leaving primarily SARS-CoV-2 reads.
   * Only the unclassified reads are retained, so if a different database is specified, make sure it doesn't have SARS-CoV-2 included.
3. `porechop` filters out adapter sequences, barcodes, and primer sequences.
   * After installing `Porechop` (via [./prepare_envs.sh](prepare_envs.sh)), additional adapters can be filtered out if added to `./Porechop/porechop/adapters.py`.