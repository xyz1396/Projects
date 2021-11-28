#include "proNovoConfig.h"
#include <mkl.h>

int mklConv1D(
    const double h[], int inch,
    const double x[], int incx,
    double y[], int incy,
    int nh, int nx, int iy0, int ny)
{
    int status;
    VSLConvTaskPtr task;
    vslsConvNewTask1D(&task, VSL_CONV_MODE_DIRECT, nh, nx, ny);
    vslConvSetStart(task, &iy0);
    status = vsldConvExec1D(task, h, inch, x, incx, y, incy);
    vslConvDeleteTask(&task);
    return status;
}

IsotopeDistribution mklSum(const IsotopeDistribution &distribution0, const IsotopeDistribution &distribution1)
{
    double ProbabilityCutoff_local = 0.000001;
    IsotopeDistribution sumDistribution;
    int nd0 = distribution0.vProb.size();
    int nd1 = distribution0.vProb.size();
    const double *d0vProb = distribution0.vProb.data();
    const double *d1vProb = distribution1.vProb.data();
    int nd2 = nd0 + nd1 - 1;
    double *d2vProb = (double *)mkl_malloc((nd2) * sizeof(double), 64);
    mklConv1D(d1vProb, 1, d0vProb, 1, d2vProb, 1, nd1, nd0, 0, nd2);
    double tmpMass = distribution0.vMass[0] + distribution0.vMass[0];
    double neutronMass = 1;
    for (int i = 0; i < nd2; i++)
    {
        sumDistribution.vProb.push_back(d2vProb[i]);
        sumDistribution.vMass.push_back(tmpMass);
        tmpMass += neutronMass;
    }
    mkl_free(d2vProb);
    return sumDistribution;
}

int main(int argc, char const *argv[])
{
    string config = "/mnt/d/projects/SIPmkl/test/SiprosConfigN15SIP.cfg";
    if (ProNovoConfig::setFilename(config))
        cout << ProNovoConfig::getSearchName() << endl;
    cout << ProNovoConfig::getFASTAfilename() << endl;
    // the original Isotopologue object from configfile
    // ProNovoConfig::configIsotopologue;
    vector<IsotopeDistribution> nature = ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution;
    for (IsotopeDistribution element : nature)
    {
        element.print();
        int n = element.vMass.size();
        cout << "neutronMass:" << endl;
        // for (int i = 1; i < n; i++)
        // {
        //     cout << element.vMass[i] - element.vMass[i - 1] << endl;
        // }
        double *vMass = element.vMass.data();
        for (int i = 1; i < n; i++)
        {
            cout << vMass[i] - vMass[i - 1] << endl;
        }
    }
    IsotopeDistribution myIso;
    ProNovoConfig::configIsotopologue.computeIsotopicDistribution("KRSKCH", myIso);
    cout << "Pep Iso Mass:" << endl;
    for (int i = 0; i < myIso.vMass.size(); i++)
    {
        cout << myIso.vMass[i] << " " << myIso.vProb[i] << endl;
    }
    return 0;
}