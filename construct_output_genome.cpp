#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <set>

using namespace std;

// fill the gaps with inexact matches
vector<pair<int, string>> alignMismatchesWithGaps(vector<pair<int, string>>& mismatches, vector<pair<int,int>>& gaps) {
    // Sort mismatches based on the starting position
    sort(mismatches.begin(), mismatches.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });
    
    size_t i = 0;
    size_t j = 0;
    vector<pair<int, string>> new_contigs;

    while (i < mismatches.size() && j < gaps.size()) {
        int start = mismatches[i].first;
        string read_with_mismatches = mismatches[i].second;
        int end = start + static_cast<int>(read_with_mismatches.size()) - 1;

        // Skip reads that start and end before the current gap range
        if (gaps[j].first > end) {
            i++;
            continue;
        }

        // Skip gaps that end before the current read
        while (j < gaps.size() && gaps[j].second < start) {
            j++;
        }

        // If there are no more gaps, break out of the loop
        if (j >= gaps.size()) {
            break;
        }

        // Find overlaps with the current gap range
        while (j < gaps.size() && start <= gaps[j].second && end >= gaps[j].first) {
            int overlap_start = max(gaps[j].first, start);        
            int overlap_end = min(gaps[j].second, end);

            // If the overlap starts beyond the read, move to the next mismatch
            if (overlap_start - start >= read_with_mismatches.length()) {
                break;
            }

            // Construct new contig from the overlapping region
            int overlap_length = overlap_end - overlap_start + 1;
            string contig = read_with_mismatches.substr(overlap_start - start, overlap_length);
            new_contigs.push_back(make_pair(overlap_start, contig));

            // If the current read extends beyond the current gap, check the next gap
            if (end > gaps[j].second) {
                j++;
            } else {
                break; // The current read does not extend to the next gap
            }
        }

        // Move to the next read for the next iteration
        i++;
    }

    return new_contigs;
}

// align matches to form longer contigs and return the remaining gaps
pair<vector<pair<int, string>>, vector<pair<int, int>>> alignMatches(int genome_len, vector<pair<int, string>>& reads) {
    set<int> seenIndexes;
    int i = 0, first_index = 0;
    string contig = "";
    vector<pair<int, string>> contigs;
    vector<pair<int, int>> gaps;

    // Sort the reads based on their start position
    sort(reads.begin(), reads.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    // Check for a gap at the start
    if (!reads.empty() && reads.front().first > 0) {
        gaps.push_back({0, reads.front().first - 1});
    }
    
    for (const auto& read : reads) {
        // Only look at unique indexes
        if (seenIndexes.find(read.first) != seenIndexes.end()) {
            continue;
        }
        seenIndexes.insert(read.first);

        if (i >= read.first) { // Overlap
            int overlap = i - read.first;
            if (overlap < read.second.length()) {
                contig += read.second.substr(overlap);
                i += read.second.length() - overlap;
            }
        } 
        else { // Gap
            if(!contig.empty()) {
                contigs.push_back({first_index, contig});
                contig = "";
            }
            gaps.push_back({i, read.first - 1});
            first_index = read.first;
            contig += read.second;
            i = read.first + read.second.length();
        }
    }
    
    // Add the last contig if there is one
    if(!contig.empty()) {
        contigs.push_back({first_index, contig});
    }

    // Add the last gap if the end of the last read is before the end of the genome
    if(i < genome_len + 1) {
        gaps.push_back({i, genome_len});
    }

    return {contigs, gaps};
}


// vector<pair<int, string>> alignMismatchesWithGaps(vector<pair<int, string>>& mismatches, vector<int>& gaps) {
    
//     // Sort mismatches based on the starting position
//     sort(mismatches.begin(), mismatches.end(), [](const auto& a, auto& b) {
//         return a.first < b.first;
//     });

//     int i = -1;
//     size_t j = 0;
//     vector<pair<int, string>> new_contigs;

//     while (i < static_cast<int>(mismatches.size()) - 1 && j < gaps.size()) {
//         i++;
//         int start = mismatches[i].first;
//         string read_with_mismatches = mismatches[i].second;
//         int end = start + static_cast<int>(read_with_mismatches.size()) - 1;

//         if (gaps[j] > end) {
//             continue;
//         }

//         while (j < gaps.size() && gaps[j] < start) {
//             j++;
//         }

//         if (j >= gaps.size()) {
//             break;
//         }

//         int overlap_start = gaps[j];
        
//         while (j < gaps.size() && start <= gaps[j] && gaps[j] <= end) {
//             j++;
//         }
        
//         int overlap_end = gaps[j - 1];

//         // Construct new contig from the overlapping region
//         int overlap_length = overlap_end - overlap_start + 1;
//         if (overlap_start - start >= read_with_mismatches.length()) {
//             continue;
//         }
//         string contig = read_with_mismatches.substr(overlap_start - start, overlap_length);
//         new_contigs.push_back(make_pair(overlap_start, contig));
//     }

//     return new_contigs;
// }


// //returns longer contigs as pairs {pos, contig}, in addition to array of gaps.
// pair<vector<pair<int, string>>, vector<int>> alignMatches(int genome_len, vector<pair<int, string>>& reads) {
//     set<int> seenIndexes;
//     int i = 0, first_index = 0;
//     string contig = "";
//     vector<pair<int, string>> contigs;
//     vector<int> gaps;
//     // Sort the reads based on their start position
//     sort(reads.begin(), reads.end(), [](const auto& a, const auto& b) {
//         return a.first < b.first;
//     });
    
//     for (const auto& read : reads) {
//         // Only look at unique indexes
//         if (seenIndexes.find(read.first) != seenIndexes.end()) {
//             continue;
//         }
//         seenIndexes.insert(read.first);
//         int start = 0;

//         if (i >= read.first) { // Overlap
//             start = i - read.first;
//         } 
//         else { // Gap
//             if(contig != "") {
//                 contigs.push_back({first_index, contig});
//             }
//             contig = "";
//             for (int j = i; j < read.first; ++j) {
//                 gaps.push_back(j);
//             }
//             i = read.first;
//             first_index = i;
//             start = 0;
//         }

//         // Interior
//         if (start > read.second.length()) {
//             continue;
//         }
//         contig += read.second.substr(start);
        
//         i += read.second.length() - start;
//     }
    
//     if(contig!="") {
//         contigs.push_back({first_index, contig});
//     }
    
//     for (int j = i; j < genome_len; ++j) {
//         gaps.push_back(j);
//     }
//     return {contigs, gaps};
// }


// write each contig on a line
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

//write to file and fill inner gaps with -
void writeToFile(const string& filePath, vector<pair<int, string>>& reads ,int ln) {
    ofstream f(filePath);
    if (!f.is_open()) {
        throw runtime_error("Could not open file: " + filePath);
    }
    
    int prev = 0;
    for (const auto& read : reads) {
        while(prev!=read.first) {
            f << "-";
            prev++;
        }
        f << read.second;// << "\n";
        prev += read.second.length();
    }
    // if (prev < ln) {
    //     while(prev < ln) { // Change this line from prev <= ln to prev < ln
    //         f << "-";
    //         prev++;
    //     }
    // }
    f.close();
}