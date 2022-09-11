import networkx as nx
from itertools import chain, product, combinations
from datatypes import Block,Permutation
from collections import defaultdict
from copy import deepcopy

class Context:
    '''
    Represents context of a sequence in repetetive block
    '''
    def __init__(self, perm, pos, left, right):
        self.perm = perm
        self.pos = pos
        self.left = left
        self.right = right

    def __str__(self):
        return "({0}, {1}, {2}, {3})".format(self.perm.chr_name, self.pos,
                                             self.left, self.right)

    def equal(self, other):
        return self.right == other.right and self.left == other.left

class PermutationContainer:
    def __init__(self,block_coords_file):
        """
        Parses permutation files referenced from recipe and filters duplications
        """

        self.perms = _parse_blocks_coords(block_coords_file)

    def get_contexts(self, repeats):
        """
        Get repeats' contexts
        """
        WINDOW = 5

        contexts = defaultdict(list)
        for perm in self.perms:
            for pos in range(len(perm.blocks)):
                block = perm.blocks[pos]
                if block.block_id not in repeats:
                    continue

                left_start = max(0, pos - WINDOW)
                left_end = max(0, pos)
                #left_context = list(map(lambda b: b.signed_id() * block.sign,
                #                        perm.blocks[left_start:left_end]))
                left_context = [b.signed_id() * block.sign for b in
                                perm.blocks[left_start:left_end]]

                right_start = min(len(perm.blocks), pos + 1)
                right_end = min(len(perm.blocks), pos + WINDOW + 1)
                #right_context = list(map(lambda b: b.signed_id() * block.sign,
                #                         perm.blocks[right_start:right_end]))
                right_context = [b.signed_id() * block.sign for b in
                                 perm.blocks[right_start:right_end]]

                if block.sign < 0:
                    left_context, right_context = (right_context[::-1],
                                                   left_context[::-1])

                contexts[block.block_id].append(Context(perm, pos, left_context,
                                                        right_context))
        return contexts

    def find_repeats(self):
        """
        Returns a set of repetitive blocks
        """
        index = defaultdict(set)
        repeats = set()
        for perm in self.perms:
            for block in perm.blocks:
                if perm.genome_name in index[block.block_id]:
                    repeats.add(block.block_id)
                else:
                    index[block.block_id].add(perm.genome_name)
        return repeats

    def resolve_repeats(self,logfile):
        '''
        splits repetetive blocks by context similarity of sequences
        into smaller blocks (without repeats)
        '''
        repeats = self.find_repeats()
        max_block_id = 0
        for perm in self.perms:
            max_block_id = max(max_block_id,
                                max([b.block_id for b in perm.blocks]))

        contexts = self.get_contexts(repeats)
        f=open(logfile,"w")

        new_block_id=max_block_id+1

        for repeat_id in sorted(contexts):
            by_genome = defaultdict(list)
            for ctx in contexts[repeat_id]:
                by_genome[ctx.perm.genome_name].append(ctx)

            profiles = _split_into_profiles(by_genome, repeats)
            f.write(str(repeat_id)+" ->")
            removed=[]
            for prof in profiles:
                if len(prof)>1:
                    ids=[]
                    for ctx in prof:
                        old_block=ctx.perm.blocks[ctx.pos]
                        new_block=deepcopy(old_block)
                        new_block.block_id=new_block_id
                        ctx.perm.blocks.append(new_block)
                    f.write(" "+str(new_block_id)+",")
                    new_block_id+=1
                else:
                    # f.write(" removed")
                    ctx=prof[0]
                    perm=ctx.perm
                    block=perm.blocks[ctx.pos]
                    removed.append([perm.genome_name, perm.chr_name, block.sign, block.start, block.end])
            f.write("\n")
            f.write("removed: ")
            for seq in removed:
                f.write(" "+seq[0]+"."+seq[1]+" "+str(seq[2])+" "+str(seq[3])+"-"+str(seq[4])+"\n")
            f.write("\n")
        f.close()

        #removing repetetive blocks from permutations and permutations without blocks

        for perm in self.perms:
            without_repetetive=[b for b in perm.blocks if b.block_id not in repeats]
            if without_repetetive:
                perm.blocks=without_repetetive
            else:
                self.perms.remove(perm)

def _max_weight_matching(graph):
    edges = nx.max_weight_matching(graph, maxcardinality=True)
    unique_edges = set()
    for v1, v2 in edges:
        if not (v2, v1) in unique_edges:
            unique_edges.add((v1, v2))

    return list(unique_edges)

def _split_into_profiles(contexts_by_genome, repeats):
    """
    Given repeat contexts in each of genomes,
    joins them into "profiles" -- sets of matched contexts
    across different genomes
    """
    genomes = sorted(contexts_by_genome.keys(),
                    key=lambda x: len(contexts_by_genome.get(x)),reverse=True)

    profiles  = [[c] for c in contexts_by_genome[genomes[0]]]

    for genome in genomes[1:]:
        #finding a matching between existing profiles and a new genome
        genome_ctxs = contexts_by_genome[genome]
        graph = nx.Graph()
        for (pr_id, prof), (ctx_id, ctx) in product(enumerate(profiles),
                                                    enumerate(genome_ctxs)):
            node_prof = "profile" + str(pr_id)
            node_genome = "genome" + str(ctx_id)
            graph.add_node(node_prof, profile=True, prof=prof)
            graph.add_node(node_genome, profile=False, ctx=ctx)

            score = _profile_similarity(prof, ctx, repeats, same_len=True)
            if score > 0:
                graph.add_edge(node_prof, node_genome, weight=score)

        edges = _max_weight_matching(graph)
        for edge in edges:
            prof_node, genome_node = edge
            if graph.nodes[genome_node]["profile"]:
                prof_node, genome_node = genome_node, prof_node
            genome_ctxs.remove(graph.nodes[genome_node]["ctx"])

            graph.nodes[prof_node]["prof"].append(graph.nodes[genome_node]["ctx"])

        profiles+=[[i] for i in genome_ctxs]

    return profiles

def _parse_blocks_coords(filename):
    """
    Parses a file with blocks coords
    """
    perm_by_id = {}
    with open(filename, "r") as f:
        header = True
        for line in f:
            line = line.strip()
            if not line:
                continue

            if header:
                if line.startswith("Seq_id"):
                    continue

                if line.startswith("-"):
                    header = False
                    continue

                chr_id, chr_size, seq_name = line.split("\t")
                tokens = seq_name.split(".", 1)
                if len(tokens) != 2:
                    raise Exception("permutation ids in " + filename +
                                        " do not follow naming convention: " +
                                        "'genome.chromosome'")

                genome_name, chr_name = tokens
                perm_by_id[chr_id] = Permutation(genome_name, chr_name,
                                                 int(chr_size), [])

            else:
                if line.startswith("Seq_id") or line.startswith("-"):
                    continue

                if line.startswith("Block"):
                    block_id = int(line.split(" ")[1][1:])
                    continue

                seq_id, sign, start, end, _length = line.split("\t")
                if sign == "-":
                    start, end = end, start
                if int(end) < int(start):
                    raise Exception("Error in permutations file format")

                sign_num = 1 if sign == "+" else -1
                perm_by_id[seq_id].blocks.append(Block(block_id, sign_num,
                                                      int(start), int(end)))

    #blocks in each permutations are sorted so that they are in the same
    #order as sequences in chromosome
    for perm in perm_by_id.values():
        perm.blocks.sort(key=lambda b: b.start)

    #out_perms = list(filter(lambda b: len(b.blocks), perm_by_id.values()))
    out_perms = [b for b in  perm_by_id.values() if len(b.blocks)]
    if not len(out_perms):
        raise Exception("Permutations file is empty")
    return out_perms

def _context_similarity(ctx_ref, ctx_trg, repeats, same_len):
    """
    Compute similarity between two contexts
    """
    def alignment(ref, trg):
        """
        Computes global alignment
        """
        GAP = -2
        def match(a, b):
            mult = 1 if abs(a) in repeats or abs(b) in repeats else 2
            if a != b:
                return -mult
            else:
                return mult

        l1, l2 = len(ref) + 1, len(trg) + 1
        table = [[0 for _ in range(l2)] for _ in range(l1)]
        if same_len:
            for i in range(l1):
                table[i][0] = i * GAP
            for i in range(l2):
                table[0][i] = i * GAP

        for i, j in product(range(1, l1), range(1, l2)):
            table[i][j] = max(table[i-1][j] + GAP, table[i][j-1] + GAP,
                              table[i-1][j-1] + match(ref[i-1], trg[j-1]))
        return table[-1][-1]

    if len(ctx_trg.left) + len(ctx_trg.right) == 0:
        return 0

    left = alignment(ctx_ref.left, ctx_trg.left)
    right = alignment(ctx_ref.right[::-1], ctx_trg.right[::-1])
    #return float(left + right) / (len(ctx_trg.left) + len(ctx_trg.right))
    return left + right


def _profile_similarity(profile, genome_ctx, repeats, same_len):
    """
    Compute similarity of set of contexts vs one context
    """
    #scores = list(map(lambda c: _context_similarity(c, genome_ctx,
    #                                repeats, same_len), profile))
    scores = [_context_similarity(c, genome_ctx, repeats, same_len) for c in profile]
    return float(sum(scores)) / len(scores)
