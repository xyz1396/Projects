#include "lib/peptide.h"
#include "lib/ms2scan.h"
#include "lib/initSIP.h"
#include <Rcpp.h>
using namespace Rcpp;

//' scorePSM
//' @param realMZ mz vector in MS2 scan
//' @param realIntensity intensity vector in MS2 scan
//' @param pepSeq a string of peptide
//' @param Atom "C13" or "N15"
//' @param Prob its SIP abundance (0.0~1.0)
//' @return a score of this PSM
//' @export
// [[Rcpp::export]]
double scorePSM(const NumericVector &realMZ, const NumericVector &realIntensity,
                const NumericVector &realCharge, const String &pepSeq, const String &Atom, double Prob)
{
    // read default config
    string config = get_extdata();
    ProNovoConfig::setFilename(config);
    // compute residue mass and prob again
    computeResidueMassIntensityAgain(Atom, Prob);
    Peptide myPep;
    string sOriginalPeptide, sProteinName;
    int ibeginPos;
    double dPeptideMass;
    char cIdentifyPrefix, cIdentifySuffix, cOriginalPrefix, cOriginalSuffix;
    myPep.setPeptide(pepSeq, sOriginalPeptide, sProteinName,
                     ibeginPos, dPeptideMass, cIdentifyPrefix,
                     cIdentifySuffix, cOriginalPrefix, cOriginalSuffix);
    map<char, double> mapResidueMass;
    // mapResidueMass is not used in this function
    myPep.preprocessing(true, mapResidueMass);
    // myPep.calculateIsotope(pepSeq, mapResidueMass);
    MS2Scan myScan;
    myScan.isMS1HighRes = true;
    myScan.vdIntensity = as<vector<double>>(realIntensity);
    myScan.vdMZ = as<vector<double>>(realMZ);
    myScan.viCharge = as<vector<int>>(realCharge);
    myScan.preprocess();
    myScan.scoreWeightSumHighMS2(&myPep);
    return myScan.vpWeightSumTopPeptides[0]->dScore;
}
