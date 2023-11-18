import random

def create_reads_from_genome(genome, gap):
    reads = []
    for i in range(0, len(genome) - 100, gap):
        reads.append(genome[i:i+100])
    return reads

def mutate_nucleotide(nucleotide):
    """Return a different nucleotide than the one provided."""
    nucleotides = ['A', 'C', 'G', 'T']
    nucleotides.remove(nucleotide)
    return random.choice(nucleotides)

def create_reads_with_mutations_from_genome(genome, gap):
    reads = []
    for i in range(0, len(genome) - 100, gap):
        read = list(genome[i:i+100])
        
        # Introduce mutations at a 1% rate
        for j in range(len(read)):
            if random.random() < 0.1:  # 1% chance of mutation
                read[j] = mutate_nucleotide(read[j])
        
        reads.append(''.join(read))
    return reads

path = "genome.txt"
genome = ""
with open(path, 'r') as f:
    for line in f.readlines():
        genome += line.strip()

reads = create_reads_with_mutations_from_genome(genome, 5)
with open("tools/generated_reads.txt", 'w') as f:
    for read in reads:
        f.write(read + "\n")