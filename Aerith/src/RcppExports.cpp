// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getUnfilteredPSMs
DataFrame getUnfilteredPSMs(String sipPath, String ftPath, size_t topN);
RcppExport SEXP _Aerith_getUnfilteredPSMs(SEXP sipPathSEXP, SEXP ftPathSEXP, SEXP topNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type sipPath(sipPathSEXP);
    Rcpp::traits::input_parameter< String >::type ftPath(ftPathSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    rcpp_result_gen = Rcpp::wrap(getUnfilteredPSMs(sipPath, ftPath, topN));
    return rcpp_result_gen;
END_RCPP
}
// getUnfilteredPeptides
DataFrame getUnfilteredPeptides(CharacterVector workingPath);
RcppExport SEXP _Aerith_getUnfilteredPeptides(SEXP workingPathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type workingPath(workingPathSEXP);
    rcpp_result_gen = Rcpp::wrap(getUnfilteredPeptides(workingPath));
    return rcpp_result_gen;
END_RCPP
}
// getFilterThreshold
DataFrame getFilterThreshold(CharacterVector workingPath, NumericVector OverallThreshold);
RcppExport SEXP _Aerith_getFilterThreshold(SEXP workingPathSEXP, SEXP OverallThresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type workingPath(workingPathSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type OverallThreshold(OverallThresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(getFilterThreshold(workingPath, OverallThreshold));
    return rcpp_result_gen;
END_RCPP
}
// getFilterThresholdTopPSMs
List getFilterThresholdTopPSMs(CharacterVector workingPath, NumericVector OverallThreshold, size_t topN);
RcppExport SEXP _Aerith_getFilterThresholdTopPSMs(SEXP workingPathSEXP, SEXP OverallThresholdSEXP, SEXP topNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type workingPath(workingPathSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type OverallThreshold(OverallThresholdSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    rcpp_result_gen = Rcpp::wrap(getFilterThresholdTopPSMs(workingPath, OverallThreshold, topN));
    return rcpp_result_gen;
END_RCPP
}
// generateOneCFG
bool generateOneCFG(String cfgPath, String outPath, String element, int pct, int center, int width);
RcppExport SEXP _Aerith_generateOneCFG(SEXP cfgPathSEXP, SEXP outPathSEXP, SEXP elementSEXP, SEXP pctSEXP, SEXP centerSEXP, SEXP widthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type cfgPath(cfgPathSEXP);
    Rcpp::traits::input_parameter< String >::type outPath(outPathSEXP);
    Rcpp::traits::input_parameter< String >::type element(elementSEXP);
    Rcpp::traits::input_parameter< int >::type pct(pctSEXP);
    Rcpp::traits::input_parameter< int >::type center(centerSEXP);
    Rcpp::traits::input_parameter< int >::type width(widthSEXP);
    rcpp_result_gen = Rcpp::wrap(generateOneCFG(cfgPath, outPath, element, pct, center, width));
    return rcpp_result_gen;
END_RCPP
}
// generateCFGs
bool generateCFGs(String cfgPath, String outPath, String element);
RcppExport SEXP _Aerith_generateCFGs(SEXP cfgPathSEXP, SEXP outPathSEXP, SEXP elementSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type cfgPath(cfgPathSEXP);
    Rcpp::traits::input_parameter< String >::type outPath(outPathSEXP);
    Rcpp::traits::input_parameter< String >::type element(elementSEXP);
    rcpp_result_gen = Rcpp::wrap(generateCFGs(cfgPath, outPath, element));
    return rcpp_result_gen;
END_RCPP
}
// precursor_peak_calculator
DataFrame precursor_peak_calculator(CharacterVector AAstr);
RcppExport SEXP _Aerith_precursor_peak_calculator(SEXP AAstrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type AAstr(AAstrSEXP);
    rcpp_result_gen = Rcpp::wrap(precursor_peak_calculator(AAstr));
    return rcpp_result_gen;
END_RCPP
}
// residue_peak_calculator_DIY
DataFrame residue_peak_calculator_DIY(String residue, String Atom, double Prob);
RcppExport SEXP _Aerith_residue_peak_calculator_DIY(SEXP residueSEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type residue(residueSEXP);
    Rcpp::traits::input_parameter< String >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< double >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(residue_peak_calculator_DIY(residue, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// precursor_peak_calculator_DIY
DataFrame precursor_peak_calculator_DIY(CharacterVector AAstr, CharacterVector Atom, NumericVector Prob);
RcppExport SEXP _Aerith_precursor_peak_calculator_DIY(SEXP AAstrSEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type AAstr(AAstrSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(precursor_peak_calculator_DIY(AAstr, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// precursor_peak_calculator_DIY_averagine
List precursor_peak_calculator_DIY_averagine(StringVector AAstrs, String Atom, double Prob);
RcppExport SEXP _Aerith_precursor_peak_calculator_DIY_averagine(SEXP AAstrsSEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type AAstrs(AAstrsSEXP);
    Rcpp::traits::input_parameter< String >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< double >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(precursor_peak_calculator_DIY_averagine(AAstrs, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// BYion_peak_calculator_DIY
DataFrame BYion_peak_calculator_DIY(CharacterVector AAstr, CharacterVector Atom, NumericVector Prob);
RcppExport SEXP _Aerith_BYion_peak_calculator_DIY(SEXP AAstrSEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type AAstr(AAstrSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(BYion_peak_calculator_DIY(AAstr, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// readOneScanMS2
List readOneScanMS2(CharacterVector ftFile, NumericVector scanCount);
RcppExport SEXP _Aerith_readOneScanMS2(SEXP ftFileSEXP, SEXP scanCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type ftFile(ftFileSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scanCount(scanCountSEXP);
    rcpp_result_gen = Rcpp::wrap(readOneScanMS2(ftFile, scanCount));
    return rcpp_result_gen;
END_RCPP
}
// readOneScanMS1
List readOneScanMS1(CharacterVector ftFile, NumericVector scanCount);
RcppExport SEXP _Aerith_readOneScanMS1(SEXP ftFileSEXP, SEXP scanCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type ftFile(ftFileSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scanCount(scanCountSEXP);
    rcpp_result_gen = Rcpp::wrap(readOneScanMS1(ftFile, scanCount));
    return rcpp_result_gen;
END_RCPP
}
// readFTheader
List readFTheader(CharacterVector ftFile);
RcppExport SEXP _Aerith_readFTheader(SEXP ftFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type ftFile(ftFileSEXP);
    rcpp_result_gen = Rcpp::wrap(readFTheader(ftFile));
    return rcpp_result_gen;
END_RCPP
}
// readScansMS1
List readScansMS1(CharacterVector ftFile, NumericVector scanCount);
RcppExport SEXP _Aerith_readScansMS1(SEXP ftFileSEXP, SEXP scanCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type ftFile(ftFileSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scanCount(scanCountSEXP);
    rcpp_result_gen = Rcpp::wrap(readScansMS1(ftFile, scanCount));
    return rcpp_result_gen;
END_RCPP
}
// readAllScanMS1
List readAllScanMS1(CharacterVector ftFile);
RcppExport SEXP _Aerith_readAllScanMS1(SEXP ftFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type ftFile(ftFileSEXP);
    rcpp_result_gen = Rcpp::wrap(readAllScanMS1(ftFile));
    return rcpp_result_gen;
END_RCPP
}
// readScansMS2
List readScansMS2(CharacterVector ftFile, NumericVector scanCount);
RcppExport SEXP _Aerith_readScansMS2(SEXP ftFileSEXP, SEXP scanCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type ftFile(ftFileSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scanCount(scanCountSEXP);
    rcpp_result_gen = Rcpp::wrap(readScansMS2(ftFile, scanCount));
    return rcpp_result_gen;
END_RCPP
}
// readAllScanMS2
List readAllScanMS2(CharacterVector ftFile);
RcppExport SEXP _Aerith_readAllScanMS2(SEXP ftFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type ftFile(ftFileSEXP);
    rcpp_result_gen = Rcpp::wrap(readAllScanMS2(ftFile));
    return rcpp_result_gen;
END_RCPP
}
// readSip
List readSip(CharacterVector sipFile);
RcppExport SEXP _Aerith_readSip(SEXP sipFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type sipFile(sipFileSEXP);
    rcpp_result_gen = Rcpp::wrap(readSip(sipFile));
    return rcpp_result_gen;
END_RCPP
}
// readSips
List readSips(CharacterVector workingPath);
RcppExport SEXP _Aerith_readSips(SEXP workingPathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type workingPath(workingPathSEXP);
    rcpp_result_gen = Rcpp::wrap(readSips(workingPath));
    return rcpp_result_gen;
END_RCPP
}
// readFilesScansTopPSMs
DataFrame readFilesScansTopPSMs(CharacterVector workingPath, size_t topN);
RcppExport SEXP _Aerith_readFilesScansTopPSMs(SEXP workingPathSEXP, SEXP topNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type workingPath(workingPathSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    rcpp_result_gen = Rcpp::wrap(readFilesScansTopPSMs(workingPath, topN));
    return rcpp_result_gen;
END_RCPP
}
// readFilesScansTopPSMsFromOneFT2
DataFrame readFilesScansTopPSMsFromOneFT2(String workingPath, String pattern, size_t topN);
RcppExport SEXP _Aerith_readFilesScansTopPSMsFromOneFT2(SEXP workingPathSEXP, SEXP patternSEXP, SEXP topNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type workingPath(workingPathSEXP);
    Rcpp::traits::input_parameter< String >::type pattern(patternSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    rcpp_result_gen = Rcpp::wrap(readFilesScansTopPSMsFromOneFT2(workingPath, pattern, topN));
    return rcpp_result_gen;
END_RCPP
}
// scoreIntensity
double scoreIntensity(const bool observed, const double realIntensity, const double expectedIntensity, const String& Atom, double Prob);
RcppExport SEXP _Aerith_scoreIntensity(SEXP observedSEXP, SEXP realIntensitySEXP, SEXP expectedIntensitySEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const bool >::type observed(observedSEXP);
    Rcpp::traits::input_parameter< const double >::type realIntensity(realIntensitySEXP);
    Rcpp::traits::input_parameter< const double >::type expectedIntensity(expectedIntensitySEXP);
    Rcpp::traits::input_parameter< const String& >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< double >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(scoreIntensity(observed, realIntensity, expectedIntensity, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// scoreIntensityByCE
double scoreIntensityByCE(const NumericVector& expectedIntensity, const NumericVector& observedIntensity);
RcppExport SEXP _Aerith_scoreIntensityByCE(SEXP expectedIntensitySEXP, SEXP observedIntensitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type expectedIntensity(expectedIntensitySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type observedIntensity(observedIntensitySEXP);
    rcpp_result_gen = Rcpp::wrap(scoreIntensityByCE(expectedIntensity, observedIntensity));
    return rcpp_result_gen;
END_RCPP
}
// scorePSM
double scorePSM(const NumericVector& realMZ, const NumericVector& realIntensity, const NumericVector& realCharge, const String& pepSeq, const String& Atom, double Prob);
RcppExport SEXP _Aerith_scorePSM(SEXP realMZSEXP, SEXP realIntensitySEXP, SEXP realChargeSEXP, SEXP pepSeqSEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type realMZ(realMZSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type realIntensity(realIntensitySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type realCharge(realChargeSEXP);
    Rcpp::traits::input_parameter< const String& >::type pepSeq(pepSeqSEXP);
    Rcpp::traits::input_parameter< const String& >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< double >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(scorePSM(realMZ, realIntensity, realCharge, pepSeq, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// scorePSMold
double scorePSMold(const NumericVector& realMZ, const NumericVector& realIntensity, const NumericVector& realCharge, const String& pepSeq, const String& Atom, double Prob);
RcppExport SEXP _Aerith_scorePSMold(SEXP realMZSEXP, SEXP realIntensitySEXP, SEXP realChargeSEXP, SEXP pepSeqSEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type realMZ(realMZSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type realIntensity(realIntensitySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type realCharge(realChargeSEXP);
    Rcpp::traits::input_parameter< const String& >::type pepSeq(pepSeqSEXP);
    Rcpp::traits::input_parameter< const String& >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< double >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(scorePSMold(realMZ, realIntensity, realCharge, pepSeq, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// calc_sum
double calc_sum(NumericVector x);
RcppExport SEXP _Aerith_calc_sum(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_sum(x));
    return rcpp_result_gen;
END_RCPP
}
// test_ftFileReader
void test_ftFileReader(CharacterVector ftFile);
RcppExport SEXP _Aerith_test_ftFileReader(SEXP ftFileSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type ftFile(ftFileSEXP);
    test_ftFileReader(ftFile);
    return R_NilValue;
END_RCPP
}
// rankyfify
NumericVector rankyfify(NumericVector a);
RcppExport SEXP _Aerith_rankyfify(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(rankyfify(a));
    return rcpp_result_gen;
END_RCPP
}
// denoiseOneMS2ScanHasCharge
List denoiseOneMS2ScanHasCharge(List scanList, float window, float step, float threshold);
RcppExport SEXP _Aerith_denoiseOneMS2ScanHasCharge(SEXP scanListSEXP, SEXP windowSEXP, SEXP stepSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type scanList(scanListSEXP);
    Rcpp::traits::input_parameter< float >::type window(windowSEXP);
    Rcpp::traits::input_parameter< float >::type step(stepSEXP);
    Rcpp::traits::input_parameter< float >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(denoiseOneMS2ScanHasCharge(scanList, window, step, threshold));
    return rcpp_result_gen;
END_RCPP
}
// writeAllScanMS2
bool writeAllScanMS2(List header, List scansList, CharacterVector ftFile);
RcppExport SEXP _Aerith_writeAllScanMS2(SEXP headerSEXP, SEXP scansListSEXP, SEXP ftFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type header(headerSEXP);
    Rcpp::traits::input_parameter< List >::type scansList(scansListSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type ftFile(ftFileSEXP);
    rcpp_result_gen = Rcpp::wrap(writeAllScanMS2(header, scansList, ftFile));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Aerith_getUnfilteredPSMs", (DL_FUNC) &_Aerith_getUnfilteredPSMs, 3},
    {"_Aerith_getUnfilteredPeptides", (DL_FUNC) &_Aerith_getUnfilteredPeptides, 1},
    {"_Aerith_getFilterThreshold", (DL_FUNC) &_Aerith_getFilterThreshold, 2},
    {"_Aerith_getFilterThresholdTopPSMs", (DL_FUNC) &_Aerith_getFilterThresholdTopPSMs, 3},
    {"_Aerith_generateOneCFG", (DL_FUNC) &_Aerith_generateOneCFG, 6},
    {"_Aerith_generateCFGs", (DL_FUNC) &_Aerith_generateCFGs, 3},
    {"_Aerith_precursor_peak_calculator", (DL_FUNC) &_Aerith_precursor_peak_calculator, 1},
    {"_Aerith_residue_peak_calculator_DIY", (DL_FUNC) &_Aerith_residue_peak_calculator_DIY, 3},
    {"_Aerith_precursor_peak_calculator_DIY", (DL_FUNC) &_Aerith_precursor_peak_calculator_DIY, 3},
    {"_Aerith_precursor_peak_calculator_DIY_averagine", (DL_FUNC) &_Aerith_precursor_peak_calculator_DIY_averagine, 3},
    {"_Aerith_BYion_peak_calculator_DIY", (DL_FUNC) &_Aerith_BYion_peak_calculator_DIY, 3},
    {"_Aerith_readOneScanMS2", (DL_FUNC) &_Aerith_readOneScanMS2, 2},
    {"_Aerith_readOneScanMS1", (DL_FUNC) &_Aerith_readOneScanMS1, 2},
    {"_Aerith_readFTheader", (DL_FUNC) &_Aerith_readFTheader, 1},
    {"_Aerith_readScansMS1", (DL_FUNC) &_Aerith_readScansMS1, 2},
    {"_Aerith_readAllScanMS1", (DL_FUNC) &_Aerith_readAllScanMS1, 1},
    {"_Aerith_readScansMS2", (DL_FUNC) &_Aerith_readScansMS2, 2},
    {"_Aerith_readAllScanMS2", (DL_FUNC) &_Aerith_readAllScanMS2, 1},
    {"_Aerith_readSip", (DL_FUNC) &_Aerith_readSip, 1},
    {"_Aerith_readSips", (DL_FUNC) &_Aerith_readSips, 1},
    {"_Aerith_readFilesScansTopPSMs", (DL_FUNC) &_Aerith_readFilesScansTopPSMs, 2},
    {"_Aerith_readFilesScansTopPSMsFromOneFT2", (DL_FUNC) &_Aerith_readFilesScansTopPSMsFromOneFT2, 3},
    {"_Aerith_scoreIntensity", (DL_FUNC) &_Aerith_scoreIntensity, 5},
    {"_Aerith_scoreIntensityByCE", (DL_FUNC) &_Aerith_scoreIntensityByCE, 2},
    {"_Aerith_scorePSM", (DL_FUNC) &_Aerith_scorePSM, 6},
    {"_Aerith_scorePSMold", (DL_FUNC) &_Aerith_scorePSMold, 6},
    {"_Aerith_calc_sum", (DL_FUNC) &_Aerith_calc_sum, 1},
    {"_Aerith_test_ftFileReader", (DL_FUNC) &_Aerith_test_ftFileReader, 1},
    {"_Aerith_rankyfify", (DL_FUNC) &_Aerith_rankyfify, 1},
    {"_Aerith_denoiseOneMS2ScanHasCharge", (DL_FUNC) &_Aerith_denoiseOneMS2ScanHasCharge, 4},
    {"_Aerith_writeAllScanMS2", (DL_FUNC) &_Aerith_writeAllScanMS2, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_Aerith(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
