/*
  Class to read CSV files. 
  ----
  NOTE:
  1. It can read a CSV file with multiple rows and cols including one header row
  2. The data get stored in a map with the entries in the first columns as the keys
  3. There are methods to search the map based on the key to access all or any 
  particular column entry.
  4. There are methods to search based on the closest matched keys as well. For this
  to work properly, please make sure that the keys are sorted in the csv file.
  5. It has primarily been written to read EMFF values from look-up tables based on Q2.
*/
#ifndef LOOK_UP_TABLE_READER_H
#define LOOK_UP_TABLE_READER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <cmath>

class LookUpTableReader {
private:
  std::unordered_map<double, std::vector<double>> dataMap;
  std::vector<double> rowKeys; // Stores keys in order of appearance  

public:
  // Function to read CSV file and populate a map
  void readCSV(const std::string& filename)
  {
    std::ifstream file(filename);
    if (!file.is_open())
      {
	std::cerr << "Error opening file w/ lookup table " << filename << std::endl;
	exit(EXIT_FAILURE);
      }

    std::string line;
    double key;

    std::getline(file, line);  // Skip header line

    while (std::getline(file, line))
      {
	std::istringstream iss(line);
	std::string keyStr;
	std::getline(iss, keyStr, ',');  // Extract key (first column)
	key = std::stod(keyStr);  // Convert key from string to double

	double value;
	std::vector<double> values;

	while (iss >> value)
	  {
	    values.push_back(value);  // Extract values (remaining columns)
	    if (iss.peek() == ',') iss.ignore();  // Ignore comma
	  }

	dataMap[key] = values;  // Store key-value pair in map
	rowKeys.push_back(key); // Maintain row order
      }

    file.close();
  }

  // Function to search for values by key
  void SearchValueByKey(double key)
  {
    auto it = dataMap.find(key);
    if (it != dataMap.end())
      {
	std::cout << "Key: " << key << std::endl;
	std::cout << "Values: ";
	for (double value : it->second)
	  {
	    std::cout << value << " ";
	  }
	std::cout << std::endl;
      }
    else
      {
	std::cout << "Key not found: " << key << std::endl;
      }
  }

  // Function to find closest key in the map to the given key
  double findClosestKey(double targetKey)
  {
    double closestKey = dataMap.begin()->first;  // Initialize with the first key
    double minDiff = std::abs(targetKey - closestKey);

    for (const auto& pair : dataMap)
      {
	double currentKey = pair.first;
	double diff = std::abs(targetKey - currentKey);

	if (diff < minDiff)
	  {
	    minDiff = diff;
	    closestKey = currentKey;
	  }
      }

    return closestKey;
  }

  // Function to search for values by closest key
  void SearchClosestValueByKey(double targetKey)
  {
    double closestKey = findClosestKey(targetKey);

    auto it = dataMap.find(closestKey);
    if (it != dataMap.end())
      {
	std::cout << "Target Key: " << targetKey << std::endl;
	std::cout << "Closest Key: " << closestKey << std::endl;
	std::cout << "Values: ";
	for (double value : it->second)
	  {
	    std::cout << value << " ";
	  }
	std::cout << std::endl;
      }
    else
      {
	std::cout << "Closest key not found: " << closestKey << std::endl;
      }
  }

  // Function to search for values by closest key
  void SearchClosestValueByKey(double targetKey, int columnNo)
  {
    double closestKey = findClosestKey(targetKey);

    auto it = dataMap.find(closestKey);
    if (it != dataMap.end())
      {
	std::cout << "Target Key: " << targetKey << std::endl;
	std::cout << "Closest Key: " << closestKey << std::endl;
	std::cout << "Values: ";
	std::cout << dataMap[closestKey][columnNo];
	std::cout << std::endl;
      }
    else
      {
	std::cout << "Closest key not found: " << closestKey << std::endl;
      }
  }

  // Function to return values by key and column no
  double GetValueByKey(double Key, int columnNo)
  {
    auto it = dataMap.find(Key);
    if (it != dataMap.end())
      {
	return dataMap[Key][columnNo];
      }
    else
      {
	std::cerr << "Key not found: " << Key << std::endl;
	return -99999;
      }
  }
  
  // Function to return values by closest key and column no
  double GetClosestValueByKey(double targetKey, int columnNo)
  {
    double closestKey = findClosestKey(targetKey);

    auto it = dataMap.find(closestKey);
    if (it != dataMap.end())
      {
	return dataMap[closestKey][columnNo];
      }
    else
      {
	std::cerr << "Closest key not found: " << closestKey << std::endl;
	return -99999;
      }
  }

  // Function to return value by row index and column number
  // NOTE: Column count starts from 0, where 0th column is the key (i.e. the Q2).
  // For other methods, however, 0th column is the 2nd column (the one after the key)
  double GetValueByRowAndColumn(int row, int columnNo)
  {
    if (row < 0 || row >= rowKeys.size())
      {
	std::cerr << "Row index out of bounds: " << row << std::endl;
	return -99999;
      }

    double key = rowKeys[row];
    if (columnNo == 0)
      return key; // Return the key itself when columnNo is 0

    auto it = dataMap.find(key);
    if (it != dataMap.end() && columnNo > 0 && columnNo <= it->second.size())
      {
	return it->second[columnNo - 1]; // Adjust index to match data structure
      }
    else
      {
	std::cerr << "Invalid row or column: row " << row << " column " << columnNo << std::endl;
	return -99999;
      }
  }
};

#endif

// int main()
// {
//     LookUpTableReader reader;
//     std::string filename = "maps/ECAL_r_c_x_y_cpr.csv";
//     reader.readCSV(filename);

//     // Example: Search for values by key
//     std::cout << reader.GetValueByKey(4,2) << "\n"; 
//     //reader.SearchValueByKey(0.0000010115794543);
//     //reader.SearchValueByKey(0.0000010232929923);

//     // Example: Search for closest values by key
//     // reader.SearchClosestValueByKey(7.4,3);
//     //reader.SearchClosestValueByKey(0.0000010115794543);
//     //reader.SearchClosestValueByKey(0.0000010232929923);

//     return 0;
// }
