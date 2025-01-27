# ecSimulator

## Introduction
ecSimulator simulates ecDNA genome structures with user-specifiable parameters for SV frequencies, size, and mechanism
of origin. Users can invoke a read simulator on the simulated structures to generate simulated sequencing data.

It takes as input a reference genome fasta and a yaml configuration file.

The latest version of ecSimulator is **0.7.1**.

## Requirements and Installation
### Basic requirements
ecSimulator requires python >=3 and the `numpy` and `intervaltree` python libraries.

```shell
conda install numpy intervaltree  # or pip3 install numpy intervaltree
git clone https://github.com/jluebeck/ecSimulator.git
```

ecSimulator also requires the AmpliconArchitect data repo to be present, which contains reference genome builds and genome annotations. 

- Download AA data repositories and set environment variable AA_DATA_REPO:
  1. Go [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/) to locate data repo(s) of your choice and make note of the URL you want.
  2. Wget and set a bash environment variable AA_DATA_REPO to point to the data_repo directory:
      ```bash
      mkdir data_repo && cd data_repo
      wget [url of reference_build]
      tar zxf [reference_build].tar.gz
      # command below exports a bash variable which is the parent directory of the individual data repos
      echo export AA_DATA_REPO=$PWD/ >> ~/.bashrc 
      touch coverage.stats && chmod a+rw coverage.stats
      source ~/.bashrc
      ```

A reference genome fasta file from which to extract genome sequences is also required.
These are available in the [AmpliconArchitect data repo](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/) for users who may already have that tool installed.

### Optional dependencies for read simulation
ecSimulator enables users to simulate reads from the resulting structures, using Nanopore or Illumina read simulators.

#### Nanopore: NanoSim
[NanoSim](https://github.com/bcgsc/NanoSim) can be installed by running
```shell
conda install -c bioconda nanosim
```

We provide a pre-generated NanoSim read simulation model generated using input data from a Promethion with version 10.4.1 flow cell and Guppy as the basecaller (version dna_r10.4.1_e8.2_400bps_hac@v3.5.2).
No additional configuration of the simulator is required by the user.


#### Illumina: Mason2
[Mason2](https://github.com/seqan/seqan/tree/main/apps/mason2) can be installed by running
```shell
conda install -c bioconda mason
```



### Computing requirements
The memory an CPU requirements of ecSimulator are minimal, but please try to have more than 4Gb RAM available. The simulator should finish within 30s-1min for typical simulation runs.
Simulating thousands of structures may take slightly longer. Read simulation is more resource intensive, and timing varies by the tool used and depth of coverage given.

## Usage
### Focal amplification simulation
An example command to generate an ecDNA of episomal origin with at least two non-overlapping genomic segments
with default SV frequencies:
>`ecSimulator/src/ecSimulator.py --ref_name GRCh38 --ref_fasta /path/to/hg38.fa -o test`

Or with customized simulation parameters:
>`ecSimulator/src/ecSimulator.py --ref_name GRCh38 --ref_fasta /path/to/hg38.fa --config_file your_config.yaml -o test`

### Nanopore read simulation

To simulate Nanopore reads from the simulated focal amplification using NanoSim, you can do the following
>`ecSimulator/src/run_nanosim.py -t [threads] -o test_amplicon_1_nanosim --amplicon_fasta test_amplicon1.fasta --amplicon_coverage 1`

Check with `-h` to see other arguments for setting background regions, and adjusting other read simulation parameters.

### Illumina read simulation
To simulate Illumina reads from the simulated focal amplification using Mason, you can do the following
>`ecSimulator/src/run_mason.py -t [threads] -o test_amplicon_1_mason --amplicon_fasta test_amplicon1.fasta --amplicon_coverage 1`

Check with `-h` to see other arguments for setting background regions, and adjusting other read simulation parameters.


## Options
ecSimulator takes the following command-line arguments:

| Command              | Default               | Description                                                                                                                                 |
|----------------------|-----------------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| `--ref_name`         |                       | (Required) One of [`hg19, GRCh37, GRCh38, mm10`]                                                                                            | 
| `--ref_fasta`        |                       | (Required) Path to the reference genome fasta file. In most cases this may just be the fasta in the AA data repo.                           |
| `--config_file`      | `default_config.yaml` | (Optional) Path to ecSimulator run configuration yaml file. Defaults to the `default_config.yaml` file in the `ecSimulator/src/` directory. |
| `-o\--output_prefix` |                       | (Required) A prefix for the output files (e.g. `my_results_dir/sim_name`)                                                                   |
| `-n\--num_amplicons` | 1                     | (Optional) The number of amplicons to simulate with the settings given by the config file.                                                  |
| `-v\--version`       |                       | (Optional) Print the version number and exit.                                                                                               |

   
### ecSimulator run configuration YAML file
The ecDNA structure simulations performed by ecSimulator can be customized a number of ways. A default version is provided in 
`src/default_config.yaml` and will be used unless a different file is given for `--config_file`. The run configuration YAML file supports the following customizations:

```yaml
random_seed: 0  # sets a seed for the simulator, set to 'None' for random seed.
target_size: 2000000  # target size for the ecDNA size, in basepairs. Target size is approximate due to duplications and deletions.
origin: "episome"  # can be either episome, chromothripsis, or two-foldback.
mean_segment_size: 150000  # refers to the average distance between breakpoints.
min_segment_size: 1000  # minimum segment length allowed between breakpoints
num_breakpoints: "auto"  # number of breakpoints to assign inside the amplicon (approximate). Does not count initial breakpoints to form & circularize interval 
num_intervals: 2  # number of non-overlapping genomic regions to use for the amplicon. Breakpoints will be assigned within these larger intervals. Recommend setting to "auto" if origin is not "episome".
same_chromosome: False  # Allow intervals to be sampled from multiple different chromosomes (if num_intervals > 1).
allow_interval_reuse: True  # allow different amplicons from the run to re-use some of the same genomic coordinates.
overlap_bed: ""  # specify a path to a bed file of regions the amplicon must overlap.
viral_insertion: False  # create a hybrid human-viral ecDNA
viral_strain: "hpv16.fasta"  # Only used if viral_insertion is 'True'. Specify the name of the viral strain to be used from the oncoviruses directory.
sv_probs:  # probability this type of SV occurs during the iterative rearrangement process. Not mutually exclusive events.
  del: 0.6
  dup: 0.5
  inv: 0.4
  trans: 0.4  # translocation (moves segment to random location inside amplicon)
  fback: 0.05  # foldback, this stacks on top of the probability of independently getting both a dup and inv at the same time.

```

ecSimulator is designed to keep ecDNA from getting approximately 20% larger or smaller than the target size, so the actual duplication and deletion probabilities in the simulated ecDNA may be lower than specified by user. 

## Outputs

### Files created by ecSimulator
* Fasta file encoding the simulated ecDNA structure.
* AmpliconArchitect-formatted `_cycles.txt` file, encoding the order and orientation of the genomic segments.
* AmpliconArchitect-formatted `_graph.txt` file, encoding the genomic segments and SVs (breakpoint graph).

Coordinates reported in the output files are based.

### Intermediate structures
ecSimulator will also create a directory of the intermediate structures from the simulation, containing the same files as described above for the final amplicon, but also for each intermediate structure created during the simulation.
This can be useful if you later want to create mixtures of heterogeneous but highly similar amplicons. 

**Users can use the `cycles_file_to_fasta.py` script to convert any intermediate (or generally, any) AA-formatted `_cycles.txt` file to a fasta.**

### NanoSim outputs

The NanoSim process will create two fastq files (aligned.fastq and unaligned.fastq) and a file describing the locations of the errors in the reads. A complete description of these files is available [here](https://github.com/bcgsc/NanoSim#2-simulation-stage-1).

To combine the background and amplicon fastq files, users can simply do
>`cat [sample]_amplicon_unaligned_reads.fastq [sample]_amplicon_aligned_reads.fastq
[sample]_background_unaligned_reads.fastq [sample]_background_aligned_reads.fastq`
