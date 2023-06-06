# Desmosome ODP
Code repository for the Desmosomal ODP project


Citation: Paper citation [doi]\(link to the paper\)


TODO IMAGE

## Layout
There are 3 main folders:  `input`, `scripts`, `results`. There is an explanation README for the folder-specific contents in each of the folders. The data that could not be uploaded here is uploaded in [Zenodo TODO](link) (the filtered final set of analyzed models for the main modeling run presented in the paper).

## Protocol Summary
TODO

## Versions and Requirements
All of the modeling and analysis was done in a multi-server setup with Linux Fedora using Bash scripts. Pre-processing (homology modeling, tomogram processing) was done in Windows 10. The python libraries (and their versions) used in the project are as follows:
1. matplotlib (3.6.2)
2. numpy (1.21.5)
3. scipy (1.9.3)
4. IMP (2.17)
5. mrcfile (1.4.3)
6. tqdm (4.62.0)
7. biopython (1.79)

Additional Linux software needed to run all the scripts in the repository:
1. [GNU Parallel](https://doi.org/10.5281/zenodo.3956817)
2. OpemMPI

## License
All the works are released under the CC-by-SA-4.0 License (see `LICENSE` file for more details)

