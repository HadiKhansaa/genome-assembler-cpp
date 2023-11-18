#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <set>

using namespace std;

vector<pair<int, string>> alignMismatchesWithGaps(vector<pair<int, string>>& mismatches, vector<int>& gaps) {
    
    // Sort mismatches based on the starting position
    sort(mismatches.begin(), mismatches.end(), [](const auto& a, auto& b) {
        return a.first < b.first;
    });

    int i = -1;
    size_t j = 0;
    vector<pair<int, string>> new_contigs;

    while (i < static_cast<int>(mismatches.size()) - 1 && j < gaps.size()) {
        i++;
        int start = mismatches[i].first;
        string read_with_mismatches = mismatches[i].second;
        int end = start + static_cast<int>(read_with_mismatches.size()) - 1;

        if (gaps[j] > end) {
            continue;
        }

        while (j < gaps.size() && gaps[j] < start) {
            j++;
        }

        if (j >= gaps.size()) {
            break;
        }

        int overlap_start = gaps[j];
        
        while (j < gaps.size() && start <= gaps[j] && gaps[j] <= end) {
            j++;
        }
        
        int overlap_end = gaps[j - 1];

        // Construct new contig from the overlapping region
        int overlap_length = overlap_end - overlap_start + 1;
        string contig = read_with_mismatches.substr(overlap_start - start, overlap_length);
        new_contigs.push_back(make_pair(overlap_start, contig));
    }

    return new_contigs;
}

//returns longer contigs as pairs {pos, contig}, in addition to array of gaps.
pair<vector<pair<int, string>>, vector<int>> alignMatches(int genome_len, vector<pair<int, string>>& reads) {
    set<int> seenIndexes;
    int i = 0, first_index = 0;
    string contig = "";
    vector<pair<int, string>> contigs;
    vector<int> gaps;
    
    // Sort the reads based on their start position
    sort(reads.begin(), reads.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });
    
    for (const auto& read : reads) {
        // Only look at unique indexes
        if (seenIndexes.find(read.first) != seenIndexes.end()) {
            continue;
        }
        seenIndexes.insert(read.first);
        int start = 0;

        if (i >= read.first) { // Overlap
            start = i - read.first;
        } 
        else { // Gap
            if(contig != "") {
                contigs.push_back({first_index, contig});
            }
            contig = "";
            for (int j = i; j < read.first; ++j) {
                gaps.push_back(j);
            }
            i = read.first;
            first_index = i;
            start = 0;
        }

        // Interior
        if (start > read.second.length()) {
            continue;
        }
        contig += read.second.substr(start);
        
        i += read.second.length() - start;
    }
    if(contig!="") {
        contigs.push_back({first_index, contig});
    }
    
    for (int j = i; j < genome_len; ++j) {
        gaps.push_back(j);
    }
    return {contigs, gaps};
}

void writeToFile(const string& filePath, vector<pair<int, string>>& reads) {
    ofstream f(filePath);
    if (!f.is_open()) {
        throw runtime_error("Could not open file: " + filePath);
    }

    for (const auto& read : reads) {
        f << read.second << "\n";
    }

    f.close();
}
