#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <chrono>
using namespace std;

class ReadFile {
private:
    std::ifstream file;
    std::vector<unsigned long> index;

public:
    ReadFile(const std::string& reads_path) {
        file.open(reads_path);
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
        unsigned long position = 0;

        while (std::getline(file, line)) {
            index.push_back(position);
            position += line.length() + 2; // Increment position by the length of the line plus 2 for \r\n
        }

        // Reset the file pointer to the beginning of the file
        file.clear();
        file.seekg(0);
    }

    std::string getRead(unsigned long readID) {
        readID -= 1; // Adjust for 0 indexing
        file.seekg(index[readID]);
        std::string line;
        std::getline(file, line);

        return line;
    }

    unsigned long getNbOfLine() {
        return index.size();
    }
};

