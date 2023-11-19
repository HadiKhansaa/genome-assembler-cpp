#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <set>
#include "index_reads.cpp"

using namespace std;

unsigned long currentContigId = 1;

// fill the gaps with inexact matches
vector<pair<int, unsigned long>> alignMismatchesWithGaps(vector<pair<int, unsigned long>>& mismatches, vector<pair<int,int>>& gaps, ofstream& contigFile, ReadFile& readFile) {
    // Sort mismatches based on the starting position
    sort(mismatches.begin(), mismatches.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });
    
    size_t i = 0;
    size_t j = 0;
    vector<pair<int, unsigned long>> new_contigs;

    while (i < mismatches.size() && j < gaps.size()) {
        int start = mismatches[i].first;
        string read_with_mismatches = readFile.getRead(mismatches[i].second);
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
            new_contigs.push_back(make_pair(overlap_start, currentContigId));
            currentContigId++;
            contigFile << contig << "\n";

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
pair<vector<pair<int, unsigned long>>, vector<pair<int, int>>> alignMatches(int genome_len, vector<pair<int, unsigned long>>& alignments, ofstream& contigFile, ReadFile& readFile) {
    int i = 0, first_index = 0;
    string contig = "";
    vector<pair<int, unsigned long>> contigs;
    vector<pair<int, int>> gaps;

    // Sort the alignments based on their start position
    sort(alignments.begin(), alignments.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    // Check for a gap at the start
    if (!alignments.empty() && alignments.front().first > 0) {
        gaps.push_back({0, alignments.front().first - 1});
    }
    
    for (const auto& p : alignments) {
        string read = readFile.getRead(p.second);
        if (i >= p.first) { // Overlap
            int overlap = i - p.first;
            if (overlap < read.length()) {
                contig += read.substr(overlap);
                i += read.length() - overlap;
            }
        } 
        else { // Gap
            if(!contig.empty()) {
                contigs.push_back({first_index, currentContigId});
                currentContigId++;
                contigFile << contig << "\n";
                contig = "";
            }
            gaps.push_back({i, p.first - 1});
            first_index = p.first;
            contig += read;
            i = p.first + read.length();
        }
    }
    
    // Add the last contig if there is one
    if(!contig.empty()) {
        contigs.push_back({first_index, currentContigId});
        currentContigId++;
        contigFile << contig << "\n";
    }

    // Add the last gap if the end of the last read is before the end of the genome
    if(i < genome_len + 1) {
        gaps.push_back({i, genome_len});
    }

    return {contigs, gaps};
}

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
void writeToFile(const string& filePath, vector<pair<int, unsigned long>>& contigs ,int ln, ReadFile& contigFile) {
    ofstream f(filePath);
    if (!f.is_open()) {
        throw runtime_error("Could not open file: " + filePath);
    }
    
    int prev = 0;
    for (const auto& contigPair : contigs) {
        string contig = contigFile.getRead(contigPair.second);
        while(prev!=contigPair.first) {
            f << "-";
            prev++;
        }
        f << contig;// << "\n";
        prev += contig.length();
    }
    // if (prev < ln) {
    //     while(prev < ln) { // Change this line from prev <= ln to prev < ln
    //         f << "-";
    //         prev++;
    //     }
    // }
    f.close();
}