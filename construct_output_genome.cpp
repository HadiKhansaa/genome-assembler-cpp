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
vector<pair<int, unsigned long>> alignMismatchesWithGaps(vector<pair<int, unsigned long>>& mismatchedReads, vector<pair<int,int>>& gaps, ofstream& contigFile, ReadFile& readFile) {
    
    // Sort mismatches based on the starting position
    sort(mismatchedReads.begin(), mismatchedReads.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });
    
    size_t i = 0;
    size_t j = 0;
    vector<pair<int, unsigned long>> new_contigs;

    // as long as there are more reads and more gaps to fill
    while (i < mismatchedReads.size() && j < gaps.size()) {
        int start = mismatchedReads[i].first;
        // access the read from the reads index
        string read_with_mismatches = readFile.getRead(mismatchedReads[i].second);
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

// align current contigs to form longer contigs and return the remaining gaps
pair<vector<pair<int, unsigned long>>, vector<pair<int, int>>> alignMatches(int genome_len, vector<pair<int, unsigned long>>& alignments, ofstream& contigFile, ReadFile& readFile, bool last_call) {
    int last_pos = 0, first_index = 0;
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
    
    for (const auto& alignPair : alignments) {
        // access the read/contig from the reads/contigs index
        string read = readFile.getRead(alignPair.second);

        // if the current contig overlaps with the previously aligned one
        if (last_pos >= alignPair.first) { 
            int overlap = last_pos - alignPair.first;
            if (overlap < read.length()) {
                contig += read.substr(overlap);
                last_pos += read.length() - overlap;
            }
        } 
        else { // if we have a gap
            if(!contig.empty()) {
                // add the so far created contig to the contigs file and add its position and index to contigs 
                contigs.push_back({first_index, currentContigId});
                currentContigId++;
                contigFile << contig << "\n";
                contig = "";
            }
            // add the gap to the gaps ranges 
            gaps.push_back({last_pos, alignPair.first - 1});
            first_index = alignPair.first;
            contig += read;
            last_pos = alignPair.first + read.length();
        }
    }
    
    // Add the last contig if there is one
    if(!contig.empty()) {
        contigs.push_back({first_index, currentContigId});
        currentContigId++;
        contigFile << contig << "\n";
    }

    // Add the last gap if the end of the last read is before the end of the genome
    if(last_pos < genome_len + 1) {
        gaps.push_back({last_pos, genome_len});
    }

    return {contigs, gaps};
}

// write each contig on a line
void writeToFile(const string& filePath, vector<pair<int, unsigned long>>& contigs, ReadFile& contigFile) {
    ofstream f(filePath);
    if (!f.is_open()) {
        throw runtime_error("Could not open file: " + filePath);
    }

    for (const auto& contigPair : contigs) {
        f << contigFile.getRead(contigPair.second) << "\n";
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

    f.close();
}