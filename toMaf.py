
import os
import subprocess
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class Block_seq:
    '''
    Represents synteny block with additional data needed
    to recover sequence from a file
    '''
    def __init__(self, chromosome, chr_size, block, seq=None):
        self.chromosome = chromosome
        self.block = block
        self.size = chr_size
        self.seq = seq


    def find_in_genome(self,seq_dict):
        '''
        finds sequence of the chromosome fragment in synteny block
        '''

        chr_seq=seq_dict[self.chromosome].seq
        start=self.block.start
        end=self.block.end
        if self.block.sign>0:
            seq=chr_seq[start:end]
        else:
            seq=chr_seq.reverse_complement()[start:end]

        return SeqRecord(seq)

    def MAF_repr(self,seq,widths=(0,0,0,0)):
        '''
        returns string representing seguence in MAF formar
        '''
        if self.block.sign>0:
            sign="+"
        else:
            sign="-"

        s_line=f"s {self.chromosome : <{widths[0]+4}}{self.block.start : >{widths[1]}} {self.block.length() : >{widths[2]}} {sign} {self.size : >{widths[3]}} {seq}"

        return(s_line)



class Synteny_block:
    '''
    represents synteny block as a set of sequences with an id
    from coords file
    '''

    def __init__(self,id,block_seqs):
        self.id=id
        self.block_seqs=block_seqs

    def add_block_seq(self,block_seq):
        self.block_seqs.append(block_seq)

    def align_block(self,seq_dict,out_dir):
        '''
        returns an alignment of sequences in synteny block
        '''
        sequences=[]
        for block_seq in self.block_seqs:
            seq_record=block_seq.find_in_genome(seq_dict)
            sequences.append(seq_record)
            # block_seq.seq=seq_record.seq

        sequences_file=os.path.join(out_dir, "sequences")
        SeqIO.write(sequences, sequences_file, "fasta")

        alignment=get_alignment(sequences_file)
        os.remove(sequences_file)

        return alignment

    def add_seq_values(self,alignment):
        '''
        inserts an alignment sequence as a value of seq attribute of respective
        block_seq in synteny block
        '''
        sequences=parse_mafft_output(alignment)

        for index in range(len(sequences)):
            self.block_seqs[index].seq=sequences[index]

    def write_to_maf(self,file,out_dir,seq_dict):
        '''
        writes a synteny block into a MAF file
        '''
        alignment = self.align_block(seq_dict,out_dir)

        name_width=max([len(block_seq.chromosome) for block_seq in self.block_seqs])
        start_width=max([len(str(block_seq.block.start)) for block_seq in self.block_seqs])
        len_width=max([len(str(block_seq.block.length())) for block_seq in self.block_seqs])
        chr_size_width=max([len(str(block_seq.size)) for block_seq in self.block_seqs])

        widths=(name_width,start_width,len_width,chr_size_width)

        sequences=parse_mafft_output(alignment)

        file.write('a\n')
        for i in range(len(sequences)):
            block_seq=self.block_seqs[i]
            seq=sequences[i]
            line=block_seq.MAF_repr(seq,widths)
            file.write(line+"\n")
        file.write(line+"\n\n")

def get_blocks_info(permutations):
    '''
    gets data needed to recover synteny block alignments
    '''
    block_info={}
    for perm in permutations:
        chromosome=perm.genome_name+"."+perm.chr_name
        chr_size=perm.seq_len
        for block in perm.blocks:
            id=block.block_id
            if id in block_info.keys():
                block_info[id].append(Block_seq(chromosome,chr_size,block,None))
            else:
                block_info[id]=[Block_seq(chromosome,chr_size,block,None)]

    block_info=[Synteny_block(id,block_info) for id,block_info in block_info.items()]
    block_info.sort(key=lambda block: block.id, reverse=True)
    return block_info

def get_alignment(fasta_file):
    '''
    runs mafft to get an alignment of sequences from fasta file
    '''
    cmdline=["mafft",fasta_file]
    proc = subprocess.run(cmdline,capture_output=True,text=True)

    return proc.stdout

def parse_mafft_output(alignment):
    '''
    parses mafft output and returns list of sequences in the same order as
    respective block_seq in block_info[id]
    '''

    seq_list=re.findall(">[^>]+",alignment)
    sequences=[]

    for seq in seq_list:
        id_end=re.match(".+",seq).end()
        sequence=seq[id_end:].replace('\n','')
        if sequence:
            sequences.append(sequence)

    return(sequences)


def get_seq_dict(seq_files):
    '''
    creates dict {seq_id:SeqRecord} for all fasta files
    '''
    seq_dict={}
    for file in seq_files:
        record_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
        # record_dict={record.id:record.seq for record in records}
        seq_dict.update(record_dict)

    return seq_dict

def write_maf(block_info, outdir, file_name="pampuch-out.maf"):
    '''
    writes MAF file with data from block_info
    '''
    path=os.path.join(outdir,file_name)

    try:
        with open(path, "x") as file:

            file.write("##maf version=1\n")


            for block in block_info:
                block.write_to_maf(file)
    except:
        with open(path, "w") as file:

            file.write("##maf version=a\n")

            for block in block_info:
                block.write_to_maf(file)

def convert_to_MAF(permutations, seq_files,
                    outdir,file="WGA_out.maf"):
    '''
    saves permutations to file in MAF format
    '''
    file=os.path.join(outdir,file)

    #get Block_seq from permutations for easier access to data
    block_info=get_blocks_info(permutations)
    #parse genome sequences in fasta format
    seq_dict=get_seq_dict(seq_files)

    file_out=open(file,"w")
    file_out.write("##maf version=1")
    for block in block_info:
        # alignment=block.align_block(seq_dict,outdir)
        # block.add_seq_values(alignment)
        block.write_to_maf(file_out,outdir,seq_dict)
    
    file_out.close()


    # # now every block in block_info contains all necessery data to
    # #recreate it in maf format

    # write_maf(block_info, outdir,file)
