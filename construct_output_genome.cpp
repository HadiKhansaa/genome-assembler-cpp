#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <set>

using namespace std;

vector<pair<int, string>> extendContigs(int genome_len, vector<pair<int, string>>& contigs, vector<pair<int, string>>& readsWithMismatches) {
    vector<pair<int, string>> extendedContigs;

    // Sort the reads by their starting positions
    sort(readsWithMismatches.begin(), readsWithMismatches.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    int currentEnd = -1; // Tracks the end position of the last contig or read added (which is filled)
    size_t readIndex = 0; // Index for iterating through the reads

    for (const auto& contig : contigs) {

        int contigStart = contig.first;
        int contigEnd = contig.first + contig.second.length() - 1;

        // Add all reads that fit in the gap before this contig
        while (readIndex < readsWithMismatches.size() && readsWithMismatches[readIndex].first < contigStart) { //read exists before current contig
            auto& read = readsWithMismatches[readIndex];
            int readEnd = read.first + read.second.length() - 1;

            if (read.first > currentEnd || readEnd > currentEnd) { // if read starts after prev contig or if it overlaps with it but with remainder to the right
                // Trim the read from left and right if needed
                if (currentEnd >= read.first) { // if it overlaps with prev contig, trim read from left
                    read.second = read.second.substr(currentEnd - read.first + 1, contigStart - read.first);
                    read.first += currentEnd - read.first + 1;
                }
                else { // else no need to trim the read
                    read.second = read.second.substr(0, contigStart - read.first);
                }
                
                extendedContigs.push_back(read);
                currentEnd += read.second.length()-1;
            }
            ++readIndex;
        }

        // Add the contig if it extends beyond the last read or contig added
        if (contigStart > currentEnd) {
            extendedContigs.push_back(contig);
            currentEnd = contigEnd;
        }

        // Move to the next contig / gap
    }

    // Add any remaining reads that come after the last contig
    while (readIndex < readsWithMismatches.size()) {
        auto& read = readsWithMismatches[readIndex];
        int readEnd = read.first + read.second.length() - 1;

        if (read.first > currentEnd) {
            // Trim the read if it surpasses the genome length
            if (readEnd >= genome_len) {
                read.second = read.second.substr(0, genome_len - read.first);
                readEnd = genome_len - 1;
            }
            extendedContigs.push_back(read);
            currentEnd = readEnd;
        }
        ++readIndex;
    }

    return extendedContigs;
}


// align matches to form longer contigs and return the remaining gaps
vector<pair<int, string>> alignMatches(int genome_len, vector<pair<int, string>>& reads) {
    set<int> seenIndexes;
    int i = 0, first_index = 0;
    string contig = "";
    vector<pair<int, string>> contigs;
    // vector<pair<int, int>> gaps;

    // Sort the reads based on their start position
    sort(reads.begin(), reads.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });
    
    for (const auto& read : reads) {
        // Only look at unique indexes
        // what if the reads are of different lengths?
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
            first_index = read.first;
            contig += read.second;
            i = read.first + read.second.length();
        }
    }
    
    // Add the last contig if there is one
    if(!contig.empty()) {
        contigs.push_back({first_index, contig});
    }

    return contigs;
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