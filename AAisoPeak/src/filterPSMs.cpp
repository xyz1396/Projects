#include "lib/PSMsFiltrator.h"
#include <Rcpp.h>
using namespace Rcpp;

//' getUnfilteredPSMs
//' @param workingPath a full path with .sip files in it
//' @return a dataframe of unique PSMs and whether it is decoy sequence
//' @export
// [[Rcpp::export]]
DataFrame getUnfilteredPSMs(CharacterVector workingPath)
{
    sipFileReader reader(as<string>(workingPath));
    reader.readAllFiles();
    PSMsFiltrator filtrator(reader.sipPSMs, 0.01);
    vector<string> psmIDs;
    vector<string> ftFileNames;
    vector<int> scanNumbers;
    vector<float> retentionTimes;
    vector<float> bestScores;
    vector<int> parentCharges;
    vector<string> searchNames;
    vector<bool> isDecoys;
    vector<double> measuredParentMasses;
    vector<double> calculatedParentMasses;
    vector<string> identifiedPepSeqs;
    vector<string> originalPepSeqs;
    vector<int> pepLengths;
    vector<string> proNames;
    vector<int> proCounts;
    for (auto psmIX : filtrator.PSMsMap)
    {
        for (size_t i = 0; i < 5; i++)
        {
            if (psmIX.second.bestScores[i] > 0)
            {
                psmIDs.push_back(psmIX.first);
                ftFileNames.push_back(psmIX.second.ftFileName);
                scanNumbers.push_back(psmIX.second.scanNumber);
                retentionTimes.push_back(psmIX.second.retentionTime);
                bestScores.push_back(psmIX.second.bestScores[i]);
                parentCharges.push_back(psmIX.second.parentCharges[i]);
                searchNames.push_back(psmIX.second.searchNames[i]);
                isDecoys.push_back(psmIX.second.isDecoys[i]);
                measuredParentMasses.push_back(psmIX.second.measuredParentMasses[i]);
                calculatedParentMasses.push_back(psmIX.second.calculatedParentMasses[i]);
                identifiedPepSeqs.push_back(psmIX.second.identifiedPepSeqs[i]);
                originalPepSeqs.push_back(psmIX.second.originalPepSeqs[i]);
                pepLengths.push_back(psmIX.second.pepLengths[i]);
                proNames.push_back(psmIX.second.proNames[i]);
                proCounts.push_back(psmIX.second.proCounts[i]);
            }
            else
                break;
        }
    }
    return DataFrame::create(Named("psmIDs") = move(psmIDs),
                             _("ftFileNames") = move(ftFileNames),
                             _["scanNumbers"] = move(scanNumbers),
                             _["retentionTimes"] = move(retentionTimes),
                             _["bestScores"] = move(bestScores),
                             _["parentCharges"] = move(parentCharges),
                             _["searchNames"] = move(searchNames),
                             _["isDecoys"] = move(isDecoys),
                             _["measuredParentMasses"] = move(measuredParentMasses),
                             _["calculatedParentMasses"] = move(calculatedParentMasses),
                             _["identifiedPepSeqs"] = move(identifiedPepSeqs),
                             _["originalPepSeqs"] = move(originalPepSeqs),
                             _["pepLengths"] = move(pepLengths),
                             _["proNames"] = move(proNames),
                             _["proCounts"] = move(proCounts));
}