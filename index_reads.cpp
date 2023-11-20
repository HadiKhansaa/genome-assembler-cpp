#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <chrono>
#include <algorithm> 
#include <cctype>
using namespace std;

//create an index for a text file to avoid storing contigs and reads in memory
class ReadFile {
private:
    std::ifstream file;
    std::vector<unsigned long> index;
    std::string reads_path;

public:
    ReadFile(const std::string& reads_path) {
        this->reads_path = reads_path;
        // file.open(reads_path);
        file.open(reads_path, std::ios::in | std::ios::binary); // Open in binary mode
        if (!file) {
            std::cerr << "Unable to open file: " << reads_path << std::endl;
            exit(1);
        }
        createIndex();
    }

    ~ReadFile() {
        if (file.is_open()) {
            file.close();
        }
    }

    void createIndex() {
        std::string line;
        std::streampos position = 0;
        file.seekg(0);
        index.push_back(0);
        while (std::getline(file, line)) {
            position = file.tellg(); // Get the current position
            index.push_back(position);
        }

        // Reset the file pointer to the beginning of the file
        file.clear();
        file.seekg(0);

        //reopen file in normal mode
        file.close();
        file.open(reads_path);
    }


    std::string getRead(unsigned long readID) {
        readID -= 1; // Adjust for 0 indexing
        file.seekg(index[readID]);
        std::string line;
        std::getline(file, line);
        if (!line.empty() && std::islower(line[0]))
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
        return line;
    }

    int close() {
        if (file.is_open()) {
            file.close();
            return 1;
        }
        return 0;
    }
};

