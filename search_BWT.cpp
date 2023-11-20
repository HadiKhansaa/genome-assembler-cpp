#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <set>

using namespace std;

// map the character to its index
int charToIndex(char ch) {
    switch (ch) {
        case '$': return 0;
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case 'N': return 4;
        case 'T': return 5;
        default: return -1;
    }
}

/*
    Function to index the Burrows-Wheeler Transform (BWT) string.
    It creates two structures: a vector for the first occurrence of each character 
    and a matrix to count the occurrences of each character up to each position in the BWT string.
*/
pair<vector<int>, vector<vector<int>>> index_bwt(const string& bwt) {
    // Initialize the first occurrence array
    vector<int> first_occ(6, 0);

    // Count the number of occurrences of each character in the BWT string.
    // This loop populates the first_occ array with these counts.
    for (char i : bwt) {
        switch(i) {
            case '$': case '#': first_occ[0]++; break;
            case 'A': first_occ[1]++; break;
            case 'C': first_occ[2]++; break;
            case 'G': first_occ[3]++; break;
            case 'N': first_occ[4]++; break;
            case 'T': first_occ[5]++; break;
        }
    }

    // Adjust the first occurrence array to represent the actual positions where each character first appears.
    int num_A = first_occ[1]; //count of character 'A'
    first_occ[0] = 0; // character at position 0 in the first column of the BWM is always $
    first_occ[1] = 1; // character at position 1 in the first column of the BWM is always A
    for (int i = 2; i <= 5; ++i) {
        int temp = first_occ[i];
        // we can find the first occurence of a character from the position of the last occurence of the preceding character
        first_occ[i] = first_occ[i - 1] + num_A; 
        num_A = temp;
    }

    // Initialize the count matrix 
    // Each row represents a character, and each column represents a position in the BWT string.
    vector<vector<int>> Counts(6, vector<int>(bwt.size()+1, 0));
    // Populate the Counts matrix with the cumulative counts of each character.
    for (int i = 1; i < bwt.size()+1; ++i) {
        // Carry forward the counts from the previous position for all characters.
        if (i > 0) {
            Counts[0][i] = Counts[0][i - 1];
            Counts[1][i] = Counts[1][i - 1];
            Counts[2][i] = Counts[2][i - 1];
            Counts[3][i] = Counts[3][i - 1];
            Counts[4][i] = Counts[4][i - 1];
            Counts[5][i] = Counts[5][i - 1];
        }
        // Increment the count of the character found at the current position.
        switch(bwt[i-1]) {
            case '$': case '#': Counts[0][i]++; break;
            case 'A': Counts[1][i]++; break;
            case 'C': Counts[2][i]++; break;
            case 'G': Counts[3][i]++; break;
            case 'N': Counts[4][i]++; break;
            case 'T': Counts[5][i]++; break;
        }
    }

    return {first_occ, Counts};
}

// Function to find the start and end positions of a substring in the BWT string.
pair<int, int> search_in_bwt(const string& bwt, const string& read, const vector<int>& first_occ, const vector<vector<int>>& Counts) {
    int p = read.size();
    int i = p - 1;
    int character;

    // Convert the last character of the read to its index representation.
    character = charToIndex(read[p-1]);

    // Initialize start (sp) and end (ep) pointers for the substring search.
    int sp = first_occ[character];
    int ep = character == 5 ? bwt.size() - 1 : first_occ[character + 1] - 1;

    // Iterate backwards through the read, updating the sp and ep pointers.
    while(sp <= ep && i >= 1) {
        character = charToIndex(read[i-1]);
        sp = first_occ[character] + Counts[character][sp];
        ep = first_occ[character] + Counts[character][ep+1] - 1;
        i--;
    }
    
    return {sp, ep};
}

// Function to search the BWT for exact matches of a given read.
// Returns the positions in the genome using the suffix array.
vector<int> search_exact(const string& bwt, const string& read, const vector<int>& suffix_array, const vector<int>& first_occ, const vector<vector<int>>& Counts) {
    
    // Find the start and end positions of the read in the BWT string.
    pair<int, int> sp_ep = search_in_bwt(bwt, read, first_occ, Counts);
    int sp = sp_ep.first;
    int ep = sp_ep.second;

    // Collect all positions where the read matches exactly.
    vector<int> hits;
    for (int i = sp; i <= ep; ++i) {
        hits.push_back(suffix_array[i]);
    }

    return hits;
}

// Function to count mismatches between part of the genome and a read.
// Returns true if mismatches are within the threshold, false otherwise.
bool countMismatchesWithGenome(const string& sub_genome, const string& read, int d) {
    int e = 0;
    for(int i=0;i<sub_genome.length();i++) {
        if(sub_genome[i] != read[i]) {
            e++;
            if(e>d) return false;
        }
    }
    return true;
}

// Function to search the BWT string for inexact matches of a read allowing a certain number of mismatches
vector<int> search_inexact(const string& genome, const string& bwt, const string& read, const vector<int>& suffix_array, const vector<int>& first_occ, const vector<vector<int>>& Counts, int d) {
    int l = read.length(); 
    
    set<int> currOccs;
    vector<int> occs;
    vector<pair<string, int>> seeds;
    int k = l/(d+1);

    // Divide the read into (d+1) seeds for searching.
    for (int i = 0; i < d; ++i) {
        seeds.push_back({read.substr(i * k, k), i * k});
    }
    seeds.push_back({read.substr(d * k), d * k});

    // Search for exact matches of each seed in the BWT string.
    for (const auto& sead : seeds) {
        string p = sead.first;
        int offset = sead.second; 

        int top = 0;
        int bottom = bwt.length() - 1;
        int currIndex = p.length() - 1;

        // search to find positions where the seed matches exactly.
        while (top <= bottom) {
            if (currIndex >= 0) {
                char symbol = charToIndex(p[currIndex]);
                currIndex--;
                if (Counts.at(symbol)[bottom + 1] - Counts.at(symbol)[top] > 0) {
                    top = first_occ.at(symbol) + Counts.at(symbol)[top];
                    bottom = first_occ.at(symbol) + Counts.at(symbol)[bottom + 1] - 1;
                } else {
                    break; // No exact match for this seed.
                }
            } else {
                // Collect all positions where the seed matches exactly.
                for (int i = top; i <= bottom; ++i) {
                    int adjustedOcc = suffix_array[i] - offset; //position of exact match - offset = position of whole read
                    if (adjustedOcc >= 0) {
                        currOccs.insert(adjustedOcc);
                    }
                }
                break;
            }
        }
    }

    // Check if the mismatches are within the threshold for each occurrence.
    for(int occ: currOccs) {
        if(occ<genome.length() && countMismatchesWithGenome(genome.substr(occ,l),read,d)) {
            occs.push_back(occ);
        }
    }
    
    return occs;
}