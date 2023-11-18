reference guided genome assembler based on bwt and suffix array
Usage:
1. place your reads (1 read per line) in reads.txt
2. place your reference genome in genome.txt
3. run command: g++ -O3 main.cpp libsais.c libsais16.c libsais64.c
4. file assembled_genome.txt will be generated with a contig per line
