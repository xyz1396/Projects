#pragma once
#include "sipFileReader.h"
#include <unordered_map>
#include <algorithm>
#include <unordered_set>
#include <tuple>

struct alignas(64) PSMinfo
{
    string ftFileName;
    int scanNumber;
    int parentCharge;
    float retentionTime;
    // save top 5 PSM with highest scores of one scan for filtering
    float bestScores[5]={0};
    int parentCharge[5];
    string searchName[5];
    bool isDecoy[5];
    double measuredParentMasses[5];
    double calculatedParentMasses[5];
    string identifiedPepSeqs[5];
    string originalPepSeqs[5];
    string pepLengths[5];
    string proNames[5];
    string proCounts[5];
    // for accurate isotopic abundance calculation
    vector<double> precursorIsotopicMasses;
    vector<double> precursorIsotopicIntensities;
    PSMinfo(string &mSearchName, string &ftFileName, sipPSM &psm);
};

class PSMsFiltrator
{
private:
public:
    float FDRthreshold;
    vector<string> tokens;
    // <ftFileName+scanNumber, PSMinfo>
    unordered_map<string, PSMinfo> PSMsMap;
    // identified peptides with different charge states
    // we always don't select singly charged precursor in expriment
    // charge 2: 3 : >3 is nearly 100:10:1
    // score of precursor with high charge is low due to resolution limit
    // peptidesCharge2 contains peptides with charge <=2
    // <score,isDecoy,pepSeq>
    vector<tuple<float, bool, string>> peptidesCharge2;
    vector<tuple<float, bool, string>> peptidesCharge3;
    vector<tuple<float, bool, string>> peptidesChargeLargerThan3;
    int decoyCountCharge2, decoyCountCharge3, decoyCountChargeLargerThan3;
    int pepCountCharge2, pepCountCharge3, pepCountChargeLargerThan3;
    float scoreThresholdCharge2, scoreThresholdCharge3, scoreThresholdChargeLargerThan3;
    PSMsFiltrator(const vector<sipPSM> &sipPSMs, float mFDRthreshold);
    ~PSMsFiltrator();
    void splitString(const string mString);
    bool detectDecoy(const string &proteinName);
    void getRentionTime(string &workPath);
    void writePercolatorTSV();
    void readPercolatorTSV();
    void fillPeptidesCharge();
    void sortPeptideAllScore();
    // <decoyCount,pepCount,scoreThreshold>
    tuple<size_t, size_t, float> getDecoyCountScoreThreshold();
    void filterPSMsMap();
    void precursorIsotopicMassesIntensities(string &workPath);
    void writeFilteredPSMs();
};
