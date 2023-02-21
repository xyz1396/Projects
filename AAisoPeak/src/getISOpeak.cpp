#include "lib/proNovoConfig.h"
#include <Rcpp.h>

using namespace Rcpp;

string get_extdata() {
  Environment base("package:base");
  Function sys_file = base["system.file"];
  StringVector resVector =
      sys_file("src/lib", "SiprosConfigN15SIP.cfg", _["package"] = "AAisoPeak");
  string res = as<std::string>(resVector);
  return res;
}

//' Simple peak calculator of natural isotopic distribution
//' @param AAstr a CharacterVector
//' @export
// [[Rcpp::export]]
DataFrame peak_calculator(CharacterVector AAstr) {
  if (AAstr.length() > 1)
    Rcerr << "only one string one time" << endl;
  string config = get_extdata();
  ProNovoConfig::setFilename(config);
  IsotopeDistribution myIso;
  string AAstr_str = as<std::string>(AAstr);
  ProNovoConfig::configIsotopologue.computeIsotopicDistribution(AAstr_str,
                                                                myIso);
  DataFrame df =
      DataFrame::create(Named("Mass") = myIso.vMass, _["Prob"] = myIso.vProb);
  return df;
}

//' Simple peak calculator of user defined isotopic distribution
//' @param AAstr a CharacterVector
//' @param Atom a CharacterVector C13 or N15
//' @param Prob a NumericVector for its abundance
//' @export
// [[Rcpp::export]]
DataFrame peak_calculator_DIY(CharacterVector AAstr, CharacterVector Atom,
                              NumericVector Prob) {
  // check and convert input
  if (AAstr.length() > 1 || Atom.length() > 1 || Prob.length() > 1)
    Rcerr << "only one string one time" << endl;
  string Atom_str = as<string>(Atom);
  double Prob_d = as<double>(Prob);
  if (Prob_d < 0 || Prob_d > 1)
    Rcout << "Wrong isotopic percentage" << endl;

  // read default config
  string config = get_extdata();
  ProNovoConfig::setFilename(config);

  // change Prob
  if (Atom_str == "C13") {
    ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].vProb[0] =
        1.0 - Prob_d;
    ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].vProb[1] =
        Prob_d;
  } else if (Atom_str == "N15") {
    ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].vProb[0] =
        1.0 - Prob_d;
    ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].vProb[1] =
        Prob_d;
  } else
    Rcerr << "this element is not support" << endl;
  ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].print();
  ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].print();

  // compute residue mass and prob again
  map<string, vector<int>>::iterator ResidueIter;
  IsotopeDistribution tempIsotopeDistribution;
  for (ResidueIter =
           ProNovoConfig::configIsotopologue.mResidueAtomicComposition.begin();
       ResidueIter !=
       ProNovoConfig::configIsotopologue.mResidueAtomicComposition.end();
       ResidueIter++) {
    if (!ProNovoConfig::configIsotopologue.computeIsotopicDistribution(
            ResidueIter->second, tempIsotopeDistribution)) {
      cerr << "ERROR: cannot calculate the isotopic distribution for residue "
           << ResidueIter->first << endl;
      return false;
    }
    ProNovoConfig::configIsotopologue
        .vResidueIsotopicDistribution[ResidueIter->first] =
        tempIsotopeDistribution;
  }

  IsotopeDistribution myIso;
  string AAstr_str = as<std::string>(AAstr);
  ProNovoConfig::configIsotopologue.computeIsotopicDistribution(AAstr_str,
                                                                myIso);
  DataFrame df =
      DataFrame::create(Named("Mass") = myIso.vMass, _["Prob"] = myIso.vProb);
  return df;
}
