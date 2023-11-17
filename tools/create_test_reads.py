def create_reads_from_genome(genome, gap):
    reads = []
    for i in range(0, len(genome) - 100, gap):
        reads.append(genome[i:i+100])
    return reads

path = "genome.txt"
genome = ""
with open(path, 'r') as f:
    for line in f.readlines():
        genome += line.strip()

reads = create_reads_from_genome(genome, 50)
with open("tools/generated_reads.txt", 'w') as f:
    for read in reads:
        f.write(read + "\n")