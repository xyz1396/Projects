#include "lib/ftFileReader.h"
#include <Rcpp.h>

using namespace Rcpp;

//' readScansMS1
//' @param ftFile a ft1 file's full path
//' @param scanNumber the scanNumber th scan
//' @export
// [[Rcpp::export]]
List readScansMS1(CharacterVector ftFile, NumericVector scanCount)
{
    ftFileReader reader(as<string>(ftFile));
    size_t scanCountInt = as<int>(scanCount);
    reader.readScans(scanCountInt);
    // in case there is not enough scan in ft file
    scanCountInt = scanCountInt < reader.Scans.size() ? scanCountInt : reader.Scans.size();
    List scanList(scanCountInt);
    Scan *mScan;
    for (size_t i = 0; i < scanCountInt; i++)
    {
        mScan = &reader.Scans[i];
        DataFrame peakDf = DataFrame::create(Named("mz") = mScan->mz,
                                             _["intensity"] = mScan->intensity,
                                             _["resolution"] = mScan->resolution,
                                             _["baseLine"] = mScan->baseLine,
                                             _["noise"] = mScan->noise,
                                             _["charge"] = mScan->charge);
        List mScanList = List::create(Named("scanNumber") = mScan->scanNumber,
                                      _["retentionTime"] = mScan->retentionTime,
                                      _["peaks"] = peakDf);
        scanList[i] = mScanList;
    }

    return scanList;
}

//' readAllScanMS1
//' @param ftFile a ft1 file's full path
//' @export
// [[Rcpp::export]]
List readAllScanMS1(CharacterVector ftFile)
{
    ftFileReader reader(as<string>(ftFile));
    reader.readAllScan();
    List scanList(reader.Scans.size());
    Scan *mScan;
    for (size_t i = 0; i < reader.Scans.size(); i++)
    {
        mScan = &reader.Scans[i];
        DataFrame peakDf = DataFrame::create(Named("mz") = mScan->mz,
                                             _["intensity"] = mScan->intensity,
                                             _["resolution"] = mScan->resolution,
                                             _["baseLine"] = mScan->baseLine,
                                             _["noise"] = mScan->noise,
                                             _["charge"] = mScan->charge);
        List mScanList = List::create(Named("scanNumber") = mScan->scanNumber,
                                      _["retentionTime"] = mScan->retentionTime,
                                      _["peaks"] = peakDf);
        scanList[i] = mScanList;
    }
    return scanList;
}

//' readScansMS2
//' @param ftFile a ft2 file's full path
//' @param scanNumber the scanNumber th scan
//' @export
// [[Rcpp::export]]
List readScansMS2(CharacterVector ftFile, NumericVector scanCount)
{
    ftFileReader reader(as<string>(ftFile));
    size_t scanCountInt = as<int>(scanCount);
    reader.readScans(scanCountInt);
    // in case there is not enough scan in ft file
    scanCountInt = scanCountInt < reader.Scans.size() ? scanCountInt : reader.Scans.size();
    List scanList(scanCountInt);
    Scan *mScan;
    for (size_t i = 0; i < scanCountInt; i++)
    {
        mScan = &reader.Scans[i];
        DataFrame peakDf = DataFrame::create(Named("mz") = mScan->mz,
                                             _["intensity"] = mScan->intensity,
                                             _["resolution"] = mScan->resolution,
                                             _["baseLine"] = mScan->baseLine,
                                             _["noise"] = mScan->noise,
                                             _["charge"] = mScan->charge);
        List mScanList = List::create(Named("scanNumber") = mScan->scanNumber,
                                      _["retentionTime"] = mScan->retentionTime,
                                      _["precursorScanNumber"] = mScan->precursorScanNumber,
                                      _["precursorMz"] = mScan->precursorMz,
                                      _["precursorCharge"] = mScan->precursorCharge,
                                      _["peaks"] = peakDf);
        scanList[i] = mScanList;
    }
    return scanList;
}

//' readAllScanMS2
//' @param ftFile a ft1 file's full path
//' @export
// [[Rcpp::export]]
List readAllScanMS2(CharacterVector ftFile)
{
    ftFileReader reader(as<string>(ftFile));
    reader.readAllScan();
    List scanList(reader.Scans.size());
    Scan *mScan;
    for (size_t i = 0; i < reader.Scans.size(); i++)
    {
        mScan = &reader.Scans[i];
        DataFrame peakDf = DataFrame::create(Named("mz") = mScan->mz,
                                             _["intensity"] = mScan->intensity,
                                             _["resolution"] = mScan->resolution,
                                             _["baseLine"] = mScan->baseLine,
                                             _["noise"] = mScan->noise,
                                             _["charge"] = mScan->charge);
        List mScanList = List::create(Named("scanNumber") = mScan->scanNumber,
                                      _["retentionTime"] = mScan->retentionTime,
                                      _["precursorScanNumber"] = mScan->precursorScanNumber,
                                      _["precursorMz"] = mScan->precursorMz,
                                      _["precursorCharge"] = mScan->precursorCharge,
                                      _["peaks"] = peakDf);
        scanList[i] = mScanList;
    }
    return scanList;
}
