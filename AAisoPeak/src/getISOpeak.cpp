#include "lib/proNovoConfig.h"
#include <Rcpp.h>

using namespace Rcpp;

string get_extdata()
{
	Environment base("package:base");
	Function sys_file = base["system.file"];
	StringVector resVector =
		sys_file("src/lib", "SiprosConfigN15SIP.cfg", _["package"] = "AAisoPeak");
	string res = as<std::string>(resVector);
	return res;
}

void computeResidueMassIntensityAgain(const string Atom_str, double Prob_d)
{
	// change Prob
	if (Atom_str == "C13")
	{
		ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].vProb[0] =
			1.0 - Prob_d;
		ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].vProb[1] =
			Prob_d;
	}
	else if (Atom_str == "N15")
	{
		ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].vProb[0] =
			1.0 - Prob_d;
		ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].vProb[1] =
			Prob_d;
	}
	else
		Rcerr << "this element is not support" << endl;
	// compute residue mass and prob again
	map<string, vector<int>>::iterator ResidueIter;
	IsotopeDistribution tempIsotopeDistribution;
	for (ResidueIter =
			 ProNovoConfig::configIsotopologue.mResidueAtomicComposition.begin();
		 ResidueIter !=
		 ProNovoConfig::configIsotopologue.mResidueAtomicComposition.end();
		 ResidueIter++)
	{
		if (!ProNovoConfig::configIsotopologue.computeIsotopicDistribution(
				ResidueIter->second, tempIsotopeDistribution))
		{
			Rcerr << "ERROR: cannot calculate the isotopic distribution for residue "
				  << ResidueIter->first << endl;
		}
		ProNovoConfig::configIsotopologue
			.vResidueIsotopicDistribution[ResidueIter->first] =
			tempIsotopeDistribution;
	}
}

//' Simple peak calculator of natural isotopic distribution
//' @param AAstr a CharacterVector
//' @export
// [[Rcpp::export]]
DataFrame precursor_peak_calculator(CharacterVector AAstr)
{
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

//' Simple peak calculator of user defined isotopic distribution of one peptide
//' @param AAstr a CharacterVector
//' @param Atom a CharacterVector C13 or N15
//' @param Prob a NumericVector for its abundance
//' @export
// [[Rcpp::export]]
DataFrame precursor_peak_calculator_DIY(CharacterVector AAstr, CharacterVector Atom,
										NumericVector Prob)
{
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
	// compute residue mass and prob again
	computeResidueMassIntensityAgain(Atom_str, Prob_d);
	// test if prob change take effct
	// ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].print();
	// ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].print();
	IsotopeDistribution myIso;
	string AAstr_str = as<std::string>(AAstr);
	ProNovoConfig::configIsotopologue.computeIsotopicDistribution(AAstr_str,
																  myIso);
	DataFrame df =
		DataFrame::create(Named("Mass") = myIso.vMass, _["Prob"] = myIso.vProb);
	return df;
}

//' peak calculator of B Y ione from of one peptide using user defined isotopic distribution
//' @param AAstr a CharacterVector
//' @param Atom a CharacterVector C13 or N15
//' @param Prob a NumericVector for its abundance
//' @export
// [[Rcpp::export]]
DataFrame BYion_peak_calculator_DIY(CharacterVector AAstr, CharacterVector Atom,
									NumericVector Prob)
{
	// check and convert input
	if (AAstr.length() > 1 || Atom.length() > 1 || Prob.length() > 1)
		Rcerr << "only one string one time" << endl;
	string AAstr_str = as<std::string>(AAstr);
	string Atom_str = as<string>(Atom);
	double Prob_d = as<double>(Prob);
	if (Prob_d < 0 || Prob_d > 1)
		Rcout << "Wrong isotopic percentage" << endl;
	// read default config
	string config = get_extdata();
	ProNovoConfig::setFilename(config);
	// compute residue mass and prob again
	computeResidueMassIntensityAgain(Atom_str, Prob_d);
	// AA string format is [AAKRCI]
	AAstr_str = "[" + AAstr_str + "]";
	vector<vector<double>> vvdYionMass, vvdYionProb, vvdBionMass, vvdBionProb;
	ProNovoConfig::configIsotopologue.computeProductIon(AAstr_str, vvdYionMass,
														vvdYionProb, vvdBionMass, vvdBionProb);
	vector<double> masses, probs;
	vector<string> kinds;
	for (size_t i = 0; i < vvdYionMass.size(); i++)
	{
		for (size_t j = 0; j < vvdYionMass[i].size(); j++)
		{
			masses.push_back(vvdYionMass[i][j]);
			probs.push_back(vvdYionProb[i][j]);
			kinds.push_back("Y" + to_string(i + 1));
		}
	}
	for (size_t i = 0; i < vvdBionMass.size(); i++)
	{
		for (size_t j = 0; j < vvdBionMass[i].size(); j++)
		{
			masses.push_back(vvdBionMass[i][j]);
			probs.push_back(vvdBionProb[i][j]);
			kinds.push_back("B" + to_string(i + 1));
		}
	}
	DataFrame df =
		DataFrame::create(Named("Mass") = move(masses),
						  _["Prob"] = move(probs), _["Kind"] = move(kinds));
	return df;
}
