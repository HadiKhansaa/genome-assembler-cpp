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
    vector<string> reads;
    try {
        genome = read_genome("genome.txt");
        reads = read_reads("reads.txt");
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

    vector<pair<int, string>> alignments;
    vector<string> unmatched_reads;
    
    //aligning
    // Start aligning timer
    start = chrono::high_resolution_clock::now();
    
    // processing reads for the first time, allowing no mismatches
    for(int i=0; i<reads.size(); i++)
    {
        string read = reads[i];

        vector<int> exact_hits = search_exact(bwt, read, SA, first_occ, counts);
        
        if(exact_hits.size() == 0) { //if no exact match found for the current read, save it to align it later with mismatches
            unmatched_reads.push_back(read);
        }
        
        for(auto hit : exact_hits){
            pair<int, string> alignment = {hit, read};
            alignments.push_back(alignment);
        }
    }

    pair<vector<pair<int, string>>, vector<pair<int,int>>> contigs_gaps = alignMatches(genome.length(), alignments);
    vector<pair<int, string>> contigs = contigs_gaps.first;
    vector<pair<int,int>> gaps = contigs_gaps.second;

    // // form longer contigs from the exact matches + array of gaps that exist for now
    // pair<vector<pair<int, string>>, vector<int>> contigs_gaps = alignMatches(genome.length(), alignments);
    // vector<pair<int, string>> contigs = contigs_gaps.first;
    // vector<int> gaps = contigs_gaps.second;


    // if we have gaps in the assembled genome and if we have reads that were not matched exactly, apply inexact matching to these reads
    if(unmatched_reads.size()>0 && gaps.size()>0) {
        vector<pair<int, string>> alignments_with_mismatches;
        int i = 1;
        for(const auto& read:unmatched_reads) {
            
            vector<int> inexact_hits = search_inexact(genome, bwt, read, SA, first_occ,counts,1);

            for(int hit: inexact_hits) {
                pair<int, string> alignment = {hit, read};
                alignments_with_mismatches.push_back(alignment);
            }

        }

        //if some reads where matched inexactly, try to fill the gaps with them
        vector<pair<int, string>> new_contigs = alignMismatchesWithGaps(alignments_with_mismatches,gaps);
        for(const auto& pair:new_contigs) {
            contigs.push_back(pair);
        }
        contigs_gaps = alignMatches(genome.length(), contigs);
        contigs = contigs_gaps.first;
        gaps = contigs_gaps.second;
    }

    // cout << "Contigs:" << endl;
    // for (const auto& pair : contigs) {
    //     cout << "Index: " << pair.first << ", String: " << pair.second << endl;
    // }

    // Printing the contents of the gaps vector
    // cout << "\nGaps:" << endl;
    // for(const auto& pair:gaps) {
    //     cout << "[" << pair.first << ", " << pair.second << "]\n";
    // }
    // cout << endl;

    // cout << "\nGaps:" << endl;
    // for (int gap : gaps) {
    //     cout << gap << " ";
    // }
    // cout << endl;

    writeToFile("assembled_genome.txt", contigs, genome.length());
    // writeToFile("assembled_genome.txt", contigs);
    

    // Stop aligning timer
    end = chrono::high_resolution_clock::now();

    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Aligning time: " << duration.count()/1000.0 << " seconds" << endl;

    // writeToFile("assembled_genome.txt", alignments);
    return 0;
}