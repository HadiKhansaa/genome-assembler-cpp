#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>

// std::unordered_map<char, int> acgt_to_1234 = {{'$', 0}, {'A', 1}, {'C', 2}, {'G', 3}, {'N', 4}, {'T', 5}};

//function to index bwt
std::pair<std::vector<int>, std::vector<std::vector<int>>> index_bwt(const std::string& bwt) {
    std::vector<int> c(6, 0);
    for (char i : bwt) {
        switch(i) {
            case '$': case '#': c[0]++; break;
            case 'A': c[1]++; break;
            case 'C': c[2]++; break;
            case 'G': c[3]++; break;
            case 'N': c[4]++; break;
            case 'T': c[5]++; break;
        }
    }

    int num_A = c[1];
    c[0] = 0; c[1] = 1;
    for (int i = 2; i <= 5; ++i) {
        int temp = c[i];
        c[i] = c[i - 1] + num_A;
        num_A = temp;
    }

    std::vector<std::vector<int>> Occ(6, std::vector<int>(bwt.size(), 0));
    for (int i = 0; i < bwt.size(); ++i) {
        if (i > 0) {
            Occ[0][i] = Occ[0][i - 1];
            Occ[1][i] = Occ[1][i - 1];
            Occ[2][i] = Occ[2][i - 1];
            Occ[3][i] = Occ[3][i - 1];
            Occ[4][i] = Occ[4][i - 1];
            Occ[5][i] = Occ[5][i - 1];
        }

        switch(bwt[i]) {
            case '$': case '#': Occ[0][i]++; break;
            case 'A': Occ[1][i]++; break;
            case 'C': Occ[2][i]++; break;
            case 'G': Occ[3][i]++; break;
            case 'N': Occ[4][i]++; break;
            case 'T': Occ[5][i]++; break;
        }
    }

    return {c, Occ};
}

//function to return the start and end positions of the sub_string in the bwt
std::pair<int, int> search_in_bwt(const std::string& bwt, const std::string& sub_string, const std::vector<int>& c, const std::vector<std::vector<int>>& Occ) {
    int p = sub_string.size();
    int i = p - 1;
    // int character = acgt_to_1234[sub_string[p - 1]];
    int character;

    switch (sub_string[p - 1]) {
        case '$': character = 0; break;
        case 'A': character = 1; break;
        case 'C':  character = 2; break;
        case 'G':  character = 3; break;
        case 'N':  character = 4; break;
        case 'T':  character = 5; break;
        default:  character = -1; // Default case for characters not in the mapping
    }

    

    int sp = c[character];
    int ep = character == 5 ? bwt.size() - 1 : c[character + 1] - 1;

    while(sp <= ep && i >= 1) {
        // character = acgt_to_1234[sub_string[i - 1]];
        switch (sub_string[i - 1]) {
            case '$': character = 0; break;
            case 'A': character = 1; break;
            case 'C':  character = 2; break;
            case 'G':  character = 3; break;
            case 'N':  character = 4; break;
            case 'T':  character = 5; break;
            default:  character = -1; // Default case for characters not in the mapping
        }
        sp = c[character] + Occ[character][sp - 1];
        ep = c[character] + Occ[character][ep] - 1;
        i--;
    }
    
    return {sp, ep};
}

//function to search bwt for posssible position hits of a read
std::vector<int> search(const std::string& bwt, const std::string& sub_string, const std::vector<int>& suffix_array, const std::vector<int>& c, const std::vector<std::vector<int>>& Occ) {
    std::pair<int, int> sp_ep = search_in_bwt(bwt, sub_string, c, Occ);
    int sp = sp_ep.first;
    int ep = sp_ep.second;

    std::vector<int> hits;
    for (int i = sp; i <= ep; ++i) {
        hits.push_back(suffix_array[i]);
    }

    return hits;
}
