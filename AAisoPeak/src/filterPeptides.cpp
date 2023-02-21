#include "lib/peptidesFiltrator.h"
#include <Rcpp.h>
using namespace Rcpp;

//' getUnfilteredPeptides
//' @param workingPath a full path with .sip files in it
//' @return a dataframe of unique peptides and whether it is decoy sequence
//' @export
// [[Rcpp::export]]
DataFrame getUnfilteredPeptides(CharacterVector workingPath)
{
    sipFileReader reader(as<string>(workingPath));
    reader.readAllFiles();
    peptidesFiltrator filtrator(reader.sipPSMs, 0.01);
    filtrator.peptideMap.merge(filtrator.peptideMapCharge2);
    filtrator.peptideMap.merge(filtrator.peptideMapCharge3);
    filtrator.peptideMap.merge(filtrator.peptideMapChargeLargerThan3);
    size_t pepCount = filtrator.peptideMap.size();
    vector<string> identifiedPeptides(pepCount);
    vector<float> bestScores(pepCount);
    vector<bool> isDecoy(pepCount);
    size_t i = 0;
    for (auto pepIX : filtrator.peptideMap)
    {
        identifiedPeptides[i] = pepIX.first;
        bestScores[i] = pepIX.second.bestScore;
        isDecoy[i] = pepIX.second.isDecoy;
        i++;
    }
    return DataFrame::create(Named("identifiedPeptides") = move(identifiedPeptides),
                             _["bestScores"] = move(bestScores),
                             _["isDecoy"] = move(isDecoy));
}

//' getFilterThreshold
//' @param workingPath a full path with .sip files in it
//' @return a dataframe about filter threshold and FDR results
//' @export
// [[Rcpp::export]]
DataFrame getFilterThreshold(CharacterVector workingPath, NumericVector OverallThreshold)
{
    sipFileReader reader(as<string>(workingPath));
    reader.readAllFiles();
    peptidesFiltrator filtrator(reader.sipPSMs, as<float>(OverallThreshold));
    filtrator.filterPeptideMap();
    vector<int> decoyCount{filtrator.decoyCountCharge2, filtrator.decoyCountCharge3,
                           filtrator.decoyCountChargeLargerThan3};
    vector<int> pepCount{filtrator.pepCountCharge2, filtrator.pepCountCharge3,
                         filtrator.pepCountChargeLargerThan3};
    vector<float> scoreThreshold{filtrator.scoreThresholdCharge2, filtrator.scoreThresholdCharge3,
                                 filtrator.scoreThresholdChargeLargerThan3};
    return DataFrame::create(
        Named("decoyCount") = decoyCount,
        _["pepCount"] = pepCount,
        _["scoreThreshold"] = scoreThreshold);
}