#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>

std::vector<std::pair<double, double>> readXYPositions(const std::string& filename) {
    std::vector<std::pair<double, double>> positions;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    
    // Skip the header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        int col = 0;
        double xpos = 0.0, ypos = 0.0;

        while (std::getline(ss, token, ',')) {
            if (col == 1) xpos = std::stod(token);
            else if (col == 2) ypos = std::stod(token);
            col++;
        }

        positions.emplace_back(xpos, ypos);
    }

    return positions;
}

int main() {
    try {
        auto xyList = readXYPositions("bad_ecal_channels_04_05_25.csv");
        for (const auto& [x, y] : xyList) {
            std::cout << "x = " << x << ", y = " << y << "\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
}
