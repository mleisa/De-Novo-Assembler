class Sequence:
    def __init__(self, lines):
        self.name = lines[0].strip()[1:]
        self.bases = "".join([x.strip() for x in lines[1:]]).upper()

    def __str__(self):
        return self.name + ": " + self.bases[:20] + "..."

    def __repr__(self):
        return self.__str__()


class Read(Sequence):
    def get_seed(self, seedlength):
        return self.bases[:seedlength]

    def replace_kmers(self, replacements):
        for key, value in replacements.items():
            if self.bases.find(key) != -1:
                self.bases = self.bases.replace(key, value)


class Reference(Sequence):
    def __init__(self, lines):
        self.kmers = None
        super().__init__(lines)

    def calculate_kmers(self, kmersize):
        self.kmers = {}
        for pos in range(0, len(self.bases) - kmersize + 1):
            kmer = self.bases[pos:(pos + kmersize)]
            if kmer not in self.kmers:
                self.kmers[kmer] = []
            self.kmers[kmer] += [pos]

    def get_kmer_positions(self, kmer):
        if self.kmers is None or len(next(iter(self.kmers))) != len(kmer):
            self.calculate_kmers(len(kmer))
        if kmer not in self.kmers:
            return []
        return self.kmers[kmer]

    def count_mismatches(self, read, position):
        mismatches = 0
        for pos in range(position, position + len(read.bases)):
            if pos >= len(self.bases):
                break
            if read.bases[pos - position] != self.bases[pos]:
                mismatches += 1
        # Count every base of the read that goes out past the end of the reference as a mismatch
        mismatches += position + len(read.bases) - pos - 1
        return mismatches


class Mapping:
    def __init__(self, reference):
        self.reference = reference
        self.reads = {}

    def add_read(self, read, position):
        if position not in self.reads:
            self.reads[position] = []
        self.reads[position] += [read]

    def get_reads_at_position(self, position):
        if position not in self.reads:
            return []
        return self.reads[position]

    def __str__(self):
        res = ["Mapping to " + self.reference.name]
        for pos in self.reads:
            res += ["  " + str(len(self.reads[pos])) + " reads mapping at " + str(pos)]
        return "\n".join(res)


class SAMWriter:
    def __init__(self, mapping):
        self.mapping = mapping

    def write_mapping(self, filename):
        ref_name = self.mapping.reference.name[:10]
        ref_len = len(self.mapping.reference.bases)

        content = '@SQ\tSN:' + ref_name + '\tLN:' + str(ref_len) + '\n'

        for read in self.mapping.reads.items():
            read_name = read[1][0].name
            flag = '0'
            position = str(read[0] + 1)
            mapq = '255'
            cigar = str(len(read[1][0].bases)) + 'M'
            r_next = '*'
            p_next = '0'
            insert_len = '0'
            seq = read[1][0].bases
            quality = '*'

            content += read_name + '\t' + flag + '\t' + ref_name + '\t' + position + '\t' + mapq + '\t' + cigar + '\t' + r_next + '\t' + p_next + '\t' + insert_len + '\t' + seq + '\t' + quality
            content += '\n'

        with open(filename, 'a') as file:
            file.write(content)


class ReadPolisher:
    kmer_spectrum = {}

    def __init__(self, kmerlen):
        self.kmer_length = kmerlen

    def add_read(self, readseq):
        for pos in range(0, len(readseq) - self.kmer_length + 1):
            kmer = readseq[pos:(pos + self.kmer_length)]
            if kmer not in self.kmer_spectrum.keys():
                self.kmer_spectrum[kmer] = 1
            else:
                self.kmer_spectrum[kmer] += 1

    def get_replacements(self, minfreq):
        bases = ["A", "G", "C", "T"]
        kmers_to_check = [key for key, val in self.kmer_spectrum.items() if val < minfreq]
        possible_corrections = {}

        for kmer in kmers_to_check:
            for test_base in bases:
                for position in range(0, len(kmer)):
                    temp_kmer = kmer[0:position] + test_base + kmer[position + 1:]

                    if self.kmer_spectrum.get(temp_kmer, None) is not None and self.kmer_spectrum[temp_kmer] >= minfreq:
                        possible_corrections[kmer] = temp_kmer

        return possible_corrections


def read_fasta(fastafile, klassname):
    klass = globals()[klassname]
    f = open(fastafile, "r")
    readlines = []
    reads = []
    for line in f:
        if line[0] == '>' and len(readlines) != 0:
            reads += [klass(readlines)]
            readlines = []
        readlines += [line]
    reads += [klass(readlines)]
    f.close()
    return reads


def map_reads(reads, reference, kmersize, max_mismatches):
    mapping = Mapping(reference)
    reference.calculate_kmers(kmersize)
    for read in reads:
        seed = read.get_seed(kmersize)
        seed_positions = reference.get_kmer_positions(seed)
        for position in seed_positions:
            mismatches = reference.count_mismatches(read, position)
            if mismatches < max_mismatches:
                mapping.add_read(read, position)
    return mapping


def main():
    reads = read_fasta("data/fluA_reads.fasta", Read.__name__)
    reference = read_fasta("data/fluA.fasta", Reference.__name__)[0]
    mapping = map_reads(reads, reference, 8, 5)
    print("Mapping reads: " + len(mapping.reads))


if __name__ == "__main__":
    main()
