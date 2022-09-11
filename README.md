# mafSimplifier

A pipeline that simpifies a whole genome alignment and produces a file in MAF format. Beside simplification of WGA it allows for repetitive blocks resolution. Synteny blocks simplification is done with [maf2synteny](https://github.com/fenderglass/maf2synteny/blob/master/README.md) and the repetitive block resolution is heavily based on the method used in [Ragout](https://github.com/fenderglass/Ragout). 

## Dependencies
- [maf2synteny](https://github.com/fenderglass/maf2synteny/blob/master/README.md)

## Usage
```
usage: main.py [-h] [-o output_dir] [-r] [-b block_sizes] recipe_file

...

positional arguments:
  recipe_file           path to recipe file

optional arguments:
  -h, --help            show this help message and exit
  -o output_dir, --outdir output_dir
                        output directory
  -r, --resolve_repeats
                        enable repeat resolution algorithm
  -b block_sizes, --block_sizes block_sizes
                        list of minBlock parameters for maf2synteny (numbers
                        separated by commas)
```
A recipe file contains paths to the alignment file (in a line started with "a") and paths to fasta files with genomes' sequences (in lines starting with "g"). Lines starting with "#" are considered to be comments.
- a recipe file example

```
a /path/to/alignment_file.maf

g /path/to/genome/sequence/Dog.fa
g /path/to/genome/sequence/Mouse.fa
g /path/to/genome/sequence/Cow.fa
g /path/to/genome/sequence/Rat.fa
g /path/to/genome/sequence/Human.fa

```
