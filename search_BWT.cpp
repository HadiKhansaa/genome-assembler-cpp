#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <set>

using namespace std;
// unordered_map<char, int> acgt_to_1234 = {{'$', 0}, {'A', 1}, {'C', 2}, {'G', 3}, {'N', 4}, {'T', 5}};

int charToIndex(char ch) {
    switch (ch) {
        case '$': return 0;
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case 'N': return 4;
        case 'T': return 5;
        default: return -1; // Default case for characters not in the mapping
    }
}


//function to index bwt
pair<vector<int>, vector<vector<int>>> index_bwt(const string& bwt) {
    vector<int> first_occ(6, 0);
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

    int num_A = first_occ[1];
    first_occ[0] = 0; first_occ[1] = 1;
    for (int i = 2; i <= 5; ++i) {
        int temp = first_occ[i];
        first_occ[i] = first_occ[i - 1] + num_A;
        num_A = temp;
    }

    vector<vector<int>> Counts(6, vector<int>(bwt.size()+1, 0));
    for (int i = 1; i < bwt.size()+1; ++i) {
        if (i > 0) {
            Counts[0][i] = Counts[0][i - 1];
            Counts[1][i] = Counts[1][i - 1];
            Counts[2][i] = Counts[2][i - 1];
            Counts[3][i] = Counts[3][i - 1];
            Counts[4][i] = Counts[4][i - 1];
            Counts[5][i] = Counts[5][i - 1];
        }

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

//function to return the start and end positions of the sub_string in the bwt
pair<int, int> search_in_bwt(const string& bwt, const string& read, const vector<int>& first_occ, const vector<vector<int>>& Counts) {
    int p = read.size();
    int i = p - 1;
    int character;

    character = charToIndex(read[p-1]);

    int sp = first_occ[character];
    int ep = character == 5 ? bwt.size() - 1 : first_occ[character + 1] - 1;

    while(sp <= ep && i >= 1) {
        character = charToIndex(read[i-1]);
        sp = first_occ[character] + Counts[character][sp];
        ep = first_occ[character] + Counts[character][ep+1] - 1;
        i--;
    }
    
    return {sp, ep};
}

//function to search bwt for posssible position hits of a read
vector<int> search_exact(const string& bwt, const string& read, const vector<int>& suffix_array, const vector<int>& first_occ, const vector<vector<int>>& Counts) {
    pair<int, int> sp_ep = search_in_bwt(bwt, read, first_occ, Counts);
    int sp = sp_ep.first;
    int ep = sp_ep.second;

    vector<int> hits;
    for (int i = sp; i <= ep; ++i) {
        hits.push_back(suffix_array[i]);
    }

    return hits;
}


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

vector<int> search_inexact(const string& genome, const string& bwt, const string& read, const vector<int>& suffix_array, const vector<int>& first_occ, const vector<vector<int>>& Counts, int d) {
    int l = read.length(); 
    
    set<int> currOccs;
    vector<int> occs;
    vector<pair<string, int>> seeds;
    int k = l/(d+1);

    // creating the d+1 seeds
    for (int i = 0; i < d; ++i) {
        seeds.push_back({read.substr(i * k, k), i * k});
    }
    seeds.push_back({read.substr(d * k), d * k});

    //matching the seeds until an exact match among them is found
    for (const auto& sead : seeds) {
        string p = sead.first;
        int offset = sead.second; 

        int top = 0;
        int bottom = bwt.length() - 1;
        int currIndex = p.length() - 1;

        while (top <= bottom) {
            if (currIndex >= 0) {
                char symbol = charToIndex(p[currIndex]);
                currIndex--;
                if (Counts.at(symbol)[bottom + 1] - Counts.at(symbol)[top] > 0) {
                    top = first_occ.at(symbol) + Counts.at(symbol)[top];
                    bottom = first_occ.at(symbol) + Counts.at(symbol)[bottom + 1] - 1;
                } else {
                    break; // no exact match for this seed
                }
            } else {
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

    //call countMismatchesWithGenome to check if the number of mismatches in the read <=d
    for(int occ: currOccs) {
        if(occ<genome.length() && countMismatchesWithGenome(genome.substr(occ,l),read,d)) {
            occs.push_back(occ);
        }
    }
    
    return occs;
}