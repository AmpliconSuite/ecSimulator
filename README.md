# ecSimulator

### Introduction
ecSimulator simulates ecDNA genome structures with user-specifiable parameters for SV frequencies, size, and mechanism
of origin. Users can invoke a read simulator on the simulated structures to generate simulated sequencing data.

It takes as input a reference genome fasta and a yaml configuration file.

The latest version of ecSimulator is **0.5.1**.

### Requirements and Installation
#### Basic requirements:
ecSimulator requires python3 and the `numpy` and `intervaltree` python libraries.

```shell
pip3 install numpy intervaltree  # or conda install, etc.
git clone https://github.com/jluebeck/ecSimulator.git
```

A reference genome fasta file for hg19, GRCh37 or GRCh38 (hg38) is also required.
These are available in the [AmpliconArchitect data repo](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/) for users who may already have that tool installed.

#### Optional dependencies for read simulation:

ecSimulator allows users to simulate nanopore reads from the resulting structures, using [NanoSim](https://github.com/bcgsc/NanoSim).
Nanosim can be installed by performing
```shell
conda install -c bioconda nanosim
```

#### Computing requirements:
The memory an CPU requirements of ecSimulator are minimal, but please try to have more than 4Gb RAM available. The simulator should finish within 30s-1min for typical simulation runs.
Simulating thousands of structures may take slightly longer.

### Usage
#### Focal amplification simulation:
An example command to generate an ecDNA of episomal origin with at least two non-overlapping genomic segments
with default SV frequencies:
```shell
ecSimulator/src/ecSimulator.py --ref_name GRCh38 --ref_fasta /path/to/hg38.fa -o test
```

Or with customized simulation parameters:
```shell
ecSimulator/src/ecSimulator.py --ref_name GRCh38 --ref_fasta /path/to/hg38.fa --config_file /path/to/your_config.yaml -o test
```

#### Nanopore read simulation:

To simulate Nanopore reads from the simulated focal amplification (using NanoSim), you can do the following
```shell
ecSimulator/src/run_nanosim.py -o test_amplicon_1_read_sim --amplicon_fasta test_amplicon1.fasta --amplicon_coverage 1
```

### Options
ecSimulator takes the following command-line arguments:

- `--ref_name [one of hg19, GRCh37, GRCh38]`: Name of reference build.

- `--ref_fasta [/path/to/ref.fa]`: Path to the reference genome fasta file (must match with `--ref_name`).

- `--config_file [/path/to/your_config.yaml]`: (Optional) Path to ecSimulator run configuration yaml file. Defaults to the `default_config.yaml` file in the `ecSimulator/src/` directory. 

- `-o | --output_prefix [your_prefix]`: Prefix for output files.

- `-n | --num_amplicons`: The number of distinct ecDNA amplicons to simulate using the given parameters. 

- `-v | --version`: Print version number and exit.
   
#### ecSimulator run configuration YAML file
The ecDNA structure simulations performed by ecSimulator can be customized a number of ways.
The run configuration YAML file supports the following customizations:

```yaml
random_seed: 0  # sets a seed for the simulator, set to 'None' for random seed.
target_size: 2000000  # target size for the ecDNA size, in basepairs. Target size is approximate due to duplications and deletions.
origin: "episome"  # type of ecDNA to simulate. Can be "episome", "chromothripsis", or "tst".
num_breakpoints: "auto"  # number of breakpoints to assign inside the amplicon (approximate). 
num_intervals: 2  # number of non-overlapping genomic regions to use for the amplicon. Breakpoints will be assigned within these larger intervals. Recommend setting to "auto" if origin is not "episome".
multichromosomal: True  # if num_intervals > 1, allow intervals to be sampled from multiple different chromosomes.
allow_interval_reuse: True  # allow different amplicons from the run to re-use some of the same genomic coordinates.
viral_insertion: False  # create a hybrid human-viral ecDNA
viral_strain: "hpv16.fasta"  # Only used if viral_insertion is 'True'. Specify the name of the viral strain to be used from the oncoviruses directory.
overlap_bed: ""  # specify a path to a bed file of regions the amplicon must overlap.
sv_probs:  # probability this type of SV occurs during the iterative rearrangement process.
  del: 0.6
  dup: 0.5
  inv: 0.4
  trans: 0.4  # translocation (moves segment to random location inside amplicon)
  fback: 0.05  # foldback, this stacks on top of the probability of independently getting both a dup and inv at the same time.

```

ecSimulator is designed to keep ecDNA from getting approximately 20% larger or smaller than the target size, so the actual duplication and deletion probabilities in the simulated ecDNA may be lower than specified by user. 

### Outputs
ecSimulator create the following:
* Fasta file encoding the simulated ecDNA structure.
* AmpliconArchitect `_cycles.txt` file, encoding the order and orientation of the genomic segments.
* AmpliconArchitect `_graph.txt` file, encoding the genomic segments and SVs (breakpoint graph).

The Nanosim process will create two fastq files (aligned.fastq and unaligned.fastq) and a file describing the locations of the errors in the reads. A complete description of these files is available [here](https://github.com/bcgsc/NanoSim#2-simulation-stage-1). 
