#include <vector>
#include <string>
#include <unordered_map>

std::pair<std::vector<int>, std::vector<std::vector<int>>> index_bwt(const std::string& bwt) {
    std::unordered_map<char, int> charToIndex = {{'$', 0}, {'A', 1}, {'C', 2}, {'G', 3}, {'T', 4}};
    std::vector<int> C(5, 0);
    std::vector<std::vector<int>> Occ(5, std::vector<int>(bwt.size() + 1, 0));

    for (char ch : bwt) {
        C[charToIndex[ch]]++;
    }

    // Accumulate the counts
    for (size_t i = 1; i < C.size(); i++) {
        C[i] += C[i - 1];
    }

    // Fill Occ table
    for (size_t i = 0; i < bwt.size(); i++) {
        for (size_t j = 0; j < Occ.size(); j++) {
            Occ[j][i + 1] = Occ[j][i];
        }
        Occ[charToIndex[bwt[i]]][i + 1]++;
    }

    return {C, Occ};
}

std::vector<int> search_bwt(const std::string& bwt, const std::string& pattern, 
                            const std::vector<int>& suffixArray, 
                            const std::vector<int>& C, 
                            const std::vector<std::vector<int>>& Occ) {
    int n = bwt.size();
    std::unordered_map<char, int> charToIndex = {{'$', 0}, {'A', 1}, {'C', 2}, {'G', 3}, {'T', 4}};
    
    int l = 0, r = n - 1;
    for (int i = pattern.size() - 1; i >= 0; --i) {
        char c = pattern[i];
        l = C[charToIndex[c]] + Occ[charToIndex[c]][l];
        r = C[charToIndex[c]] + Occ[charToIndex[c]][r + 1] - 1;
        if (l > r) return {}; // no match
    }

    std::vector<int> positions;
    for (int i = l; i <= r; ++i) {
        positions.push_back(suffixArray[i]);
    }
    return positions;
}