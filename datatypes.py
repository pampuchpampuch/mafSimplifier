class Block:
    """
    Represents synteny block
    """
    def __init__(self, block_id, sign, start=None, end=None):
        self.block_id = block_id
        self.sign = sign
        self.start = start
        self.end = end

    def length(self):
        if self.start is None or self.end is None:
            return None

        assert self.end >= self.start
        return self.end - self.start

    def signed_id(self):
        return self.block_id * self.sign

class Permutation:
    """
    Represents signed permutation
    Stores info about a chromosome, including its fragments from
    synteny blocks
    """
    def __init__(self, genome_name, chr_name, seq_len, blocks):
        self.genome_name = genome_name
        self.chr_name = chr_name
        self.seq_start = 0
        self.seq_end = seq_len
        self.seq_len = seq_len
        self.blocks = blocks

    def length(self):
        assert self.seq_end > self.seq_start
        return self.seq_end - self.seq_start

    def name(self):
        if self.seq_start == 0 and self.seq_end == self.seq_len:
            return self.chr_name
        else:
            return "{0}[{1}:{2}]".format(self.chr_name, self.seq_start,
                                         self.seq_end)

    def __lt__(self, other):
        return repr(self) < repr(other)

    def __repr__(self):
        return ("[{0}, {1}, {2}, b:{3}, e:{4}]"
                    .format(self.genome_name, self.chr_name,
                            [b.signed_id() for b in self.blocks],
                            self.seq_start, self.seq_end))
