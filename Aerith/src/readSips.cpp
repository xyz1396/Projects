#include "lib/sipFileReader.h"
#include <Rcpp.h>
using namespace Rcpp;

//' readSip
//' @param sipFile a .sip file's full path
//' @export
// [[Rcpp::export]]
List readSip(CharacterVector sipFile)
{
    sipFileReader reader;
    reader.readOneFile(as<string>(sipFile));
    DataFrame psmDf = DataFrame::create(Named("scanNumbers") = reader.currentSipPSM.scanNumbers,
                                        _["parentCharges"] = reader.currentSipPSM.parentCharges,
                                        _["measuredParentMasses"] = reader.currentSipPSM.measuredParentMasses,
                                        _["calculatedParentMasses"] = reader.currentSipPSM.calculatedParentMasses,
                                        _["scores"] = reader.currentSipPSM.scores,
                                        _["ranks"] = reader.currentSipPSM.ranks,
                                        _["identifiedPeptides"] = reader.currentSipPSM.identifiedPeptides,
                                        _["originalPeptides"] = reader.currentSipPSM.originalPeptides,
                                        _["proteinNames"] = reader.currentSipPSM.proteinNames);
    List mSipList = List::create(Named("fileName") = reader.currentSipPSM.fileName,
                                 _["scanType"] = reader.currentSipPSM.scanType,
                                 _["searchName"] = reader.currentSipPSM.searchName,
                                 _["scoringFunction"] = reader.currentSipPSM.scoringFunction,
                                 _["PSM"] = psmDf);
    return mSipList;
}

//' readSips
//' @param workingPath a full path with .sip files in it
//' @export
// [[Rcpp::export]]
List readSips(CharacterVector workingPath)
{
    sipFileReader reader(as<string>(workingPath));
    reader.readAllFiles();
    List psmList(reader.sipPSMs.size());
    sipPSM *mPSM;
    for (size_t i = 0; i < reader.sipPSMs.size(); i++)
    {
        mPSM = &reader.sipPSMs[i];
        // use move to speed up copy vector
        DataFrame psmDf = DataFrame::create(Named("scanNumbers") = move(mPSM->scanNumbers),
                                            _["parentCharges"] = move(mPSM->parentCharges),
                                            _["measuredParentMasses"] = move(mPSM->measuredParentMasses),
                                            _["calculatedParentMasses"] = move(mPSM->calculatedParentMasses),
                                            _["scores"] = move(mPSM->scores),
                                            _["ranks"] = move(mPSM->ranks),
                                            _["identifiedPeptides"] = move(mPSM->identifiedPeptides),
                                            _["originalPeptides"] = move(mPSM->originalPeptides),
                                            _["proteinNames"] = move(mPSM->proteinNames));
        List mSipList = List::create(Named("fileName") = mPSM->fileName,
                                     _["scanType"] = mPSM->scanType,
                                     _["searchName"] = mPSM->searchName,
                                     _["scoringFunction"] = mPSM->scoringFunction,
                                     _["PSM"] = psmDf);
        psmList[i] = move(mSipList);
    }
    return psmList;
}

//' readFilesScansTopPSMs read each scan's top PSMs from multiple .sip files
//' @param workingPath a full path with .sip files in it
//' @param topN store top N PSMs of each scan of one .FT file
//' @export
// [[Rcpp::export]]
DataFrame readFilesScansTopPSMs(CharacterVector workingPath, size_t topN)
{
    sipFileReader reader(as<string>(workingPath));
    reader.topN = topN;
    reader.readAllFilesTopPSMs();
    sipPSM topPSMs = reader.convertFilesScansTopPSMs();
    DataFrame psmDf = DataFrame::create(Named("fileNames") = move(topPSMs.fileNames),
                                        _["scanNumbers"] = move(topPSMs.scanNumbers),
                                        _["parentCharges"] = move(topPSMs.parentCharges),
                                        _["measuredParentMasses"] = move(topPSMs.measuredParentMasses),
                                        _["calculatedParentMasses"] = move(topPSMs.calculatedParentMasses),
                                        _["searchNames"] = move(topPSMs.searchNames),
                                        _["scores"] = move(topPSMs.scores),
                                        _["identifiedPeptides"] = move(topPSMs.identifiedPeptides),
                                        _["originalPeptides"] = move(topPSMs.originalPeptides),
                                        _["proteinNames"] = move(topPSMs.proteinNames));
    return psmDf;
}