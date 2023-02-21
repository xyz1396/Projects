#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
namespace fs = std::filesystem;
using namespace std;

struct alignas(64) sipPSM
{
    string fileName;
    vector<int> scanNumbers;
    vector<int> parentCharges;
    vector<double> measuredParentMasses;
    vector<double> calculatedParentMasses;
    string scanType;
    string searchName;
    string scoringFunction;
    vector<int> ranks;
    vector<float> scores;
    vector<string> identifiedPeptides;
    vector<string> originalPeptides;
    vector<string> proteinNames;
};

class sipFileReader
{
private:
public:
    string workingPath;
    vector<string> sipFileNames;
    vector<sipPSM> sipPSMs;
    vector<string> tokens;
    sipPSM currentSipPSM;
    string currentFilePath;
    fstream sipFileStream;
    // for single file
    sipFileReader();
    // for multi files;
    sipFileReader(string mWorkingPath);
    ~sipFileReader();
    void splitString(const string &mString);
    // fill vectors in currentSipPSM
    void fillVectors();
    void readOneFile(string sipFileName);
    void readAllFiles();
};
