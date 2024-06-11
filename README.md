# Chimera_Buster


This package takes concensus fasta files from the MrHamer pipeline and eliminates chimeric reads by comparing the UMI sequences and finding any matches in the 5' or 3' UMIs and keeping the sequence that has the highest prevalence.

> Note: This code currently relies on a very specific setup of the concesus files to work properly. These files are found as follows:

> input_file_name = [MrHamer_folder]/clustering_consensus/[sample type folder]/clusters_consensus.fasta 

> input_prelim_name = [MrHamer_folder]/clustering/[sample type folder]/clusters_consensus.fasta 

## Installation



## Usage
    python3 CLI_chimera_buster_pd.py [options] input_final_name input_prelim_name output_file_prefix

positional arguments:
| Arguement | Function |
| ------ | ------ |
|input_file_name  |    Designates clusters_concesus.fasta file from the clustering_consensus folder to be filtered. This is required.|
|  input_prelim_name     |    Designates clusters_concesus.fasta file from the clustering folder to be filtered. This is required..|
|  output_file_prefix  |  Designates output file prefix. This is required.|

options:
| Arguement | Function |
| ------ | ------ |
|-h, --help |  show this help message and exit |
|-m int, --mismatch int|  Designates the maximum number of mismatched/indel bases allowed. Default is 1. |

## License

MIT - Copyright (c) 2024 Jessica Lauren Albert

