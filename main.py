import os
import argparse
import subprocess
from permutations import PermutationContainer
import toMaf

#configuration parameters
config_param =  {
            "maf2synteny" :
            [
                (30, 10),
                (100, 100),
                (500, 1000),
                (1000, 5000),
                (5000,15000)
            ]
        }


def parse_recipe(recipe_file):
    '''
    parses recipe file
    '''
    recipe_dict={}
    recipe_dict["genome_seqs"]=[]
    with open(recipe_file, "r") as f:

        for line in f:
            line=line.strip()

            if not line or line[0]=="#":
                continue

            elif line[0]=="g":
                recipe_dict["genome_seqs"].append(line[2:])

            elif line[0]=="a":
                recipe_dict["alignment_file"]=line[2:]

            elif line[0]=="b":
                recipe_dict["blocks"]=line[2:]

            else:
                raise Exception("Unknown parameter at the begining of the recipe file line")

    return recipe_dict


def make_params_file(params, out_file):
    '''
    creates file with parameters for m2s
    '''
    assert len(params)

    with open(out_file, "w") as f:
        for k, d in params:
            f.write("{0} {1}\n".format(k, d))


def run_m2s(maf_file,out_dir,min_blocks_list):
    '''
    Runs maf2synteny
    '''
    params_file = os.path.join(out_dir, "simpl_params.txt")
    make_params_file(config_param["maf2synteny"], params_file)

    M2S_EXEC="maf2synteny"

    cmdline=[M2S_EXEC,maf_file,"-o", out_dir, "-s",
            params_file,"-b", ",".join(map(str,min_blocks_list))]
    cmdline=[M2S_EXEC,maf_file,"-o", out_dir,"-b",",".join(map(str,min_blocks_list))]
    subprocess.run(cmdline)

    os.remove(params_file)

    return True


def main():
    parser = argparse.ArgumentParser(
                description="...")
    parser.add_argument("recipe", metavar="recipe_file",
                        help="path to recipe file")
    parser.add_argument("-o", "--outdir", dest="out_dir",default="pampuch-out",
                        metavar="output_dir",
                        help="output directory")
    parser.add_argument("-r","--resolve_repeats", action="store_true", default=False,
                        dest="resolve_repeats",
                        help="enable repeat resolution algorithm")
    parser.add_argument("-b", "--block_sizes", dest="block_sizes",default="5000",
                    metavar="block_sizes",
                    help="list of minBlock parameters for maf2synteny (numbers separated by commas)")
    args = parser.parse_args()


    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    recipe=parse_recipe(args.recipe)

    maf_file=recipe["alignment_file"]
    seq_files=recipe["genome_seqs"]
    block_sizes=args.block_sizes.split(",")

    run_m2s(maf_file,args.out_dir,block_sizes)


    for size in block_sizes: 
        coords_file = os.path.join(args.out_dir,
                                str(size),"blocks_coords.txt")
        perms=PermutationContainer(coords_file)
        if args.resolve_repeats:
            logfile_blocks = os.path.join(args.out_dir,"logfile_"+str(size)+".txt")
            perms.resolve_repeats(logfile_blocks)
        maf_file="WGA_out.maf"
        toMaf.convert_to_MAF(perms.perms, seq_files, args.out_dir,maf_file)


if __name__ == "__main__":
    main()
