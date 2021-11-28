# include "ms2scanvector.h"
using namespace std;

string sFT2filename="/mnt/d/projects/SiprosTest/test/Ecoli_OU.ms2";
string sOutputDirectory="./";
string sConfigFilename="";
bool bScreenOutput=true;

int main()
{
    MS2ScanVector * pMainMS2ScanVector = new MS2ScanVector(sFT2filename, sOutputDirectory, sConfigFilename, bScreenOutput);
    if (pMainMS2ScanVector->loadMassData()) 
    cout << "Reading MS2 scan file " << sFT2filename << " succeed" << endl;
}

