#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <set>

void writeToFile(const std::string& filePath, std::vector<std::pair<int, std::string>>& reads) {
    std::ofstream f(filePath);
    if (!f.is_open()) {
        throw std::runtime_error("Could not open file: " + filePath);
    }

    std::set<int> seenIndexes;
    int i = 0;

    // Sort the reads based on their start position
    std::sort(reads.begin(), reads.end(), [](const auto& a, const auto& b) {
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
        } else { // Gap
            f.write(std::string(read.first - i, '\n').c_str(), read.first - i);
            i = read.first;
        }

        // Interior
        if (start > static_cast<int>(read.second.length())) {
            continue;
        }

        f << read.second.substr(start);
        i += read.second.length() - start;
    }
}
