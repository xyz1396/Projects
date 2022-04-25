#include "sipFileReader.h"

sipFileReader::sipFileReader()
{
}

// get all .sip files' full path in workingPath
sipFileReader::sipFileReader(string mWorkingPath) : workingPath(mWorkingPath)
{
    string ext(".sip");
    // recursive_directory_iterator can get .sip file in all child path
    for (auto &p : fs::directory_iterator(workingPath))
    {
        if (p.path().extension() == ext)
            sipFileNames.push_back(p.path().string());
    }
    if (sipFileNames.size() == 0)
        cout << "no .sip file was found in " << workingPath << " !" << endl;
}

sipFileReader::~sipFileReader()
{
}

void sipFileReader::splitString(const string &mString)
{
    string sep = "\t";
    size_t start = 0;
    size_t end = mString.find(sep);
    tokens.clear();
    while (end != std::string::npos)
    {
        tokens.push_back(mString.substr(start, end - start));
        start = end + sep.length();
        end = mString.find(sep, start);
    }
    tokens.push_back(mString.substr(start));
}

void sipFileReader::fillVectors()
{
    currentSipPSM.scanNumbers.push_back(stoi(tokens[1]));
    currentSipPSM.parentCharges.push_back(stoi(tokens[2]));
    currentSipPSM.measuredParentMasses.push_back(stod(tokens[3]));
    currentSipPSM.calculatedParentMasses.push_back(stod(tokens[4]));
    currentSipPSM.ranks.push_back(stoi(tokens[8]));
    currentSipPSM.scores.push_back(stof(tokens[9]));
    currentSipPSM.identifiedPeptides.push_back(tokens[10]);
    currentSipPSM.originalPeptides.push_back(tokens[11]);
    currentSipPSM.proteinNames.push_back(tokens[12]);
}

void sipFileReader::readOneFile(string sipFileName)
{
    setlocale(LC_ALL, "C");
    ios_base::sync_with_stdio(false);
    if (fs::exists(sipFileName))
    {
        sipFileStream.open(sipFileName.c_str(), ios::in);
        if (!sipFileStream.is_open())
        {
            cout << "Cannot open " << sipFileName << endl;
        }
        string currentLine;
        // ignore annotation in sip result file
        while (!sipFileStream.eof())
        {
            getline(sipFileStream, currentLine);
            if (currentLine[0] != '#')
                break;
        }
        // read first line
        if (!sipFileStream.eof())
        {
            getline(sipFileStream, currentLine);
            splitString(currentLine);
            currentSipPSM.fileName = tokens[0];
            currentSipPSM.scanType = tokens[5];
            currentSipPSM.searchName = tokens[6];
            currentSipPSM.scoringFunction = tokens[7];
            fillVectors();
        }
        // read more lines
        // ignore the last /n
        while (sipFileStream.peek() != EOF)
        {
            getline(sipFileStream, currentLine);
            splitString(currentLine);
            fillVectors();
        }
        if (sipFileStream.is_open())
            sipFileStream.close();
    }
    else
    {
        cout << sipFileName << " does not exists" << endl;
    }
    
}

void sipFileReader::readAllFiles()
{
    for (string sipFileName : sipFileNames)
    {
        currentSipPSM = sipPSM();
        readOneFile(sipFileName);
        sipPSMs.push_back(currentSipPSM);
    }
}