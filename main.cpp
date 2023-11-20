#include <vector>
#include <string>
#include <iostream>
#include <chrono>
#include <unordered_map>
#include "search_BWT.cpp"
#include "construct_bwt_SA.cpp"
#include "construct_output_genome.cpp"
using namespace std;

int main(){

    //create bwt and sa
    pair<string, vector<int>> bwt_SA = construct_bwt_SA();
    string bwt = bwt_SA.first;
    vector<int> SA = bwt_SA.second;

    //get the genome and reads
    string genome;
    try {
        genome = read_genome("reference.txt");
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }
    //indexing
    // Start indexing timer
    auto start = chrono::high_resolution_clock::now();
    pair<vector<int>, vector<vector<int>>> c_occ = index_bwt(bwt);
    // Stop indexing timer
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "indexing time: " << duration.count()/1000.0 << " seconds" << endl;

    vector<int> first_occ = c_occ.first;
    vector<vector<int>> counts = c_occ.second;

    vector<unsigned long> unmatched_reads;
    
    //aligning
    // Start aligning timer
    start = chrono::high_resolution_clock::now();

    ReadFile readFile("reads.txt");
    vector<pair<int, unsigned long>> alignments;
    std::ifstream file("reads.txt");
    std::string read;
    unsigned long readId=1;

    // Creating index for the reads and searching for exact matches among them
    while(getline(file, read))
    {     
        
        std::transform(read.begin(), read.end(), read.begin(), ::toupper);

        vector<int> exact_hits = search_exact(bwt, read, SA, first_occ, counts);
        
        //if no exact match found for the current read, save it to align it later with mismatches
        if(exact_hits.size() == 0) { 
            unmatched_reads.push_back(readId);
        }
        
        // add the hits to alignments to be processed later into longer contigs
        for(auto hit : exact_hits){
            pair<int, unsigned long> alignment = {hit, readId};
            alignments.push_back(alignment);
        }
        readId++;
    }

    /* 
        creating longer contigs from the exact matches, 
        writing them into a file, and spotting the current gaps in the genome
    */
    std::ofstream contigTextFile("contigs.txt");
    pair<vector<pair<int,  unsigned long>>, vector<pair<int,int>>> contigs_gaps = alignMatches(genome.length(), alignments, contigTextFile, readFile);
    vector<pair<int, unsigned long>> contigs = contigs_gaps.first;
    vector<pair<int,int>> gaps = contigs_gaps.second;

    /*
        if we have gaps in the assembled genome and if we have reads that were not
        matched exactly, apply inexact matching to these reads
    */
    if(unmatched_reads.size()>0 && gaps.size()>0) {
        vector<pair<int, unsigned long>> alignments_with_mismatches;
        int i = 1;
        for(const auto& readId:unmatched_reads) {
            string read = readFile.getRead(readId);
                
            vector<int> inexact_hits = search_inexact(genome, bwt, read, SA, first_occ,counts,1);

            // add the hits to alignments to be used later to fill gaps
            for(int hit: inexact_hits) {
                pair<int, unsigned long> alignment = {hit, readId};
                alignments_with_mismatches.push_back(alignment);
            }

        }

        //if some reads where matched inexactly, try to fill the gaps with them
        vector<pair<int, unsigned long>> new_contigs = alignMismatchesWithGaps(alignments_with_mismatches,gaps,contigTextFile,readFile);
        
        // add the new contigs to the old ones
        for(const auto& pair:new_contigs) {
            contigs.push_back(pair);
        }
    }

    contigTextFile.close();

    // call alignMatches again on the combined contigs to form longer contigs if possible
    ReadFile contigFile("contigs.txt");
    std::ofstream contigTextFile2("contigs.txt", std::ios::app);
    if(unmatched_reads.size()>0 && gaps.size()>0) {
        contigs_gaps = alignMatches(genome.length(), contigs, contigTextFile2, contigFile);
        contigs = contigs_gaps.first;
        gaps = contigs_gaps.second;
    }

    contigTextFile2.close();

    // write the output contigs into a file 
    ReadFile contigFile2("contigs.txt");
    writeToFile("assembled_genome.txt", contigs, contigFile2);
    contigTextFile.close();

    // Stop aligning timer
    end = chrono::high_resolution_clock::now();

    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Aligning time: " << duration.count()/1000.0 << " seconds" << endl;

    // writeToFile("assembled_genome.txt", alignments);
    return 0;
}