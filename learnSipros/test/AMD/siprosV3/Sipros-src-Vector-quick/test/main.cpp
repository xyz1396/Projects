#include "isotopologue.h"
#include "proNovoConfig.h"
#include "ms2scanvector.h"

void test_loadFT()
{
    cout << "\n"
         << "test_loadFT" << endl;
    MS2ScanVector *pMainMS2ScanVector = new MS2ScanVector(
        "timeCompare/AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2",
        "timeCompare",
        "timeCompare/SiproConfig.N15_0Pct.cfg",
        true);
    cout << "load " << pMainMS2ScanVector->loadFT2file() << endl;
    cout << "loaded " << pMainMS2ScanVector->vpAllMS2Scans.size() << " scans" << endl;
    // pMainMS2ScanVector->vpAllMS2Scans[66].print();
    // pMainMS2ScanVector->vpAllMS2ScanPtrs[66]->print();
    bool testPassed = true;
    for (size_t i = 6666; i < 6686; i++)
        cout << pMainMS2ScanVector->vpAllMS2ScanPtrs[i]->dParentNeutralMass << endl;
    for (size_t i = 0; i < pMainMS2ScanVector->vpAllMS2Scans[66].vdIntensity.size(); i++)
    {
        if (pMainMS2ScanVector->vpAllMS2Scans[66].vdIntensity[i] != pMainMS2ScanVector->vpAllMS2ScanPtrs[66]->vdIntensity[i])
        {
            testPassed = false;
            break;
        }
    }
    if (testPassed)
        cout << "test_loadFT passed" << endl;
    delete pMainMS2ScanVector;
}

void test_constructor()
{
    cout << "test_constructor" << endl;
    vector<double> vM = {200, 250};
    vector<double> vP = {0.5, 0.5};
    IsotopeDistribution a(vM, vP);
    a.print();
    IsotopeDistribution b;
    vector<double> vM2 = {300, 400, 500};
    vector<double> vP2 = {0.3, 0.4, 0.3};
    IsotopeDistribution c(vM2, vP2);
    c = a;
    c.print();
    a.print();
    b = a;
    b.print();
}

void test_filter()
{
    Isotopologue mIsotopologue;
    vector<double> vM = {200, 201, 202, 203};
    vector<double> vP = {1e-10, 0.5, 0.4, 0.003};
    IsotopeDistribution a(vM, vP);
    cout << "\n"
         << "test filter" << endl;
    cout << a.vProb.size() << endl;
    a.filterProbCutoff(0.01);
    cout << a.vProb.size() << endl;
    a.print();
}

void test_sum()
{
    cout << "\n"
         << "test_sum" << endl;
    Isotopologue mIsotopologue;
    vector<double> vM = {200, 250};
    vector<double> vP = {0.5, 0.5};
    IsotopeDistribution a(vM, vP);
    IsotopeDistribution b;
    vector<double> vM2 = {300, 400, 500};
    vector<double> vP2 = {0.3, 0.4, 0.3};
    IsotopeDistribution c(vM2, vP2);
    b = mIsotopologue.sum(a, c);
    b.print();
    b = mIsotopologue.sum(b, c);
    b.print();
}

void test_multiply()
{
    cout << "\n"
         << "test_multiply" << endl;
    Isotopologue mIsotopologue;
    vector<double> vM = {200, 201, 202};
    vector<double> vP = {0.97, 0.02, 0.01};
    IsotopeDistribution a(vM, vP);
    IsotopeDistribution b = mIsotopologue.multiply(a, 3);
    b.print();
    b = mIsotopologue.multiply(a, -3);
    b.print();
}

void test_computeProductIon()
{
    cout << "\n"
         << "test_computeProductIon" << endl;
    cout << "load " << ProNovoConfig::setFilename("timeCompare/SiproConfig.N15_0Pct.cfg") << endl;
    // cout << "\n"
    //      << ProNovoConfig::getProtonMass() << endl;
    vector<vector<double>> mvvdYionMass;
    vector<vector<double>> mvvdYionProb;
    vector<vector<double>> mvvdBionMass;
    vector<vector<double>> mvvdBionProb;
    ProNovoConfig::configIsotopologue.computeProductIon("[SKJHGHGHGHGHGHGHGC]", mvvdYionMass, mvvdYionProb,
                                                        mvvdBionMass, mvvdBionProb);
    cout << "Mass " << '\t' << "Inten" << endl;
    for (size_t i = 0; i < mvvdBionMass.size(); i++)
    {
        for (size_t j = 0; j < mvvdBionMass[i].size(); j++)
        {
            cout << setprecision(8) << mvvdBionMass[i][j] << '\t' << mvvdBionProb[i][j] << endl;
        }
    }
}

int main(int argc, char const *argv[])
{
    test_loadFT();
    test_constructor();
    test_filter();
    test_sum();
    test_multiply();
    test_computeProductIon();
    return 0;
}
