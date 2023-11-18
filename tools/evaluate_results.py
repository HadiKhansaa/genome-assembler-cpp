def calculate_similarity(genome1, genome2):
    if len(genome1) != len(genome2):
        print("genomes of different length")
    genome1 = genome1.upper()
    genome2 = genome2.upper()
    length = len(genome1)
    num_matches = sum(1 for base1, base2 in zip(genome1, genome2) if base1 == base2)
    similarity_score = (num_matches / length) * 100

    return similarity_score

genome = ""
assembeled_genome = ""
with open("genome.txt", 'r') as f:
    for line in f.readlines():
        genome+=line.strip()
f.close()

with open("assembled_genome.txt", 'r') as f:
    for line in f.readlines():
        assembeled_genome+=line.strip()
f.close()

print("Similarity: " + str(calculate_similarity(genome, assembeled_genome)) + "%")
