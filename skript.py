from mapper import read_fasta, ReadPolisher, Reference, Read, map_reads, SAMWriter


def get_replacement(reads, kmer_len, min_freq):
    polisher = ReadPolisher(kmer_len)

    for read in reads:
        polisher.add_read(read.bases)

    return polisher.get_replacements(min_freq)


def replace_reads(reads, replacements):
    for read in reads:
        read.replace_kmers(replacements)


def get_SAM(reads, reference, output, kmer_len, cutoff):
    reads = read_fasta(reads, Read.__name__)
    replacements = get_replacement(reads, kmer_len, cutoff)
    replace_reads(reads, replacements)
    reference = read_fasta(reference, Reference.__name__)[0]
    mapping = map_reads(reads, reference, kmer_len, cutoff)

    sam = SAMWriter(mapping)
    sam.write_mapping(output)


# get_SAM("data/fluA_reads.fasta", "data/fluA.fasta", "data/fluA_mapping.sam", 17, 3)
get_SAM("data/patient1.fasta", "data/rpoB.fasta", "data/patient1_12-2.sam", 12, 2)
get_SAM("data/patient2.fasta", "data/rpoB.fasta", "data/patient2_12-2.sam", 12, 2)
get_SAM("data/patient3.fasta", "data/rpoB.fasta", "data/patient3_12-2.sam", 12, 2)
get_SAM("data/patient4.fasta", "data/rpoB.fasta", "data/patient4_12-2.sam", 12, 2)
