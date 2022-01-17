#include "isotopologue.h"

IsotopeDistribution::IsotopeDistribution() : aSize(0), massArray(nullptr), probArray(nullptr)
{
}

IsotopeDistribution::IsotopeDistribution(const IsotopeDistribution &mIsotopeDistribution)
{
	aSize = mIsotopeDistribution.aSize;
	massArray = new double[aSize];
	probArray = new double[aSize];
	copy(mIsotopeDistribution.massArray, mIsotopeDistribution.massArray + aSize, massArray);
	copy(mIsotopeDistribution.probArray, mIsotopeDistribution.probArray + aSize, probArray);
}

IsotopeDistribution::IsotopeDistribution(const vector<double> &vItsMass, const vector<double> &vItsProb)
{
	aSize = vItsMass.size();
	massArray = new double[aSize];
	probArray = new double[aSize];
	copy(vItsMass.begin(), vItsMass.end(), massArray);
	copy(vItsProb.begin(), vItsProb.end(), probArray);
}

IsotopeDistribution::IsotopeDistribution(int vItsSize) : aSize(vItsSize),
														 massArray(new double[vItsSize]),
														 probArray(new double[vItsSize]) {}

IsotopeDistribution &IsotopeDistribution::operator=(const IsotopeDistribution &mIsotopeDistribution)
{
	if (this == &mIsotopeDistribution)
		return *this;
	aSize = mIsotopeDistribution.aSize;
	delete[] massArray;
	delete[] probArray;
	massArray = new double[aSize];
	probArray = new double[aSize];
	copy(mIsotopeDistribution.massArray, mIsotopeDistribution.massArray + aSize, massArray);
	copy(mIsotopeDistribution.probArray, mIsotopeDistribution.probArray + aSize, probArray);
	return *this;
}

IsotopeDistribution::~IsotopeDistribution()
{
	delete[] massArray;
	delete[] probArray;
}

void IsotopeDistribution::print()
{
	cout << "Mass " << '\t' << "Inten" << endl;
	for (unsigned int i = 0; i < aSize; i++)
	{
		cout << setprecision(8) << massArray[i] << '\t' << probArray[i] << endl;
	}
}

double IsotopeDistribution::getMostAbundantMass()
{
	double dMaxProb = 0;
	double dMass = 0;
	for (unsigned int i = 0; i < aSize; ++i)
	{
		if (dMaxProb < probArray[i])
		{
			dMaxProb = probArray[i];
			dMass = massArray[i];
		}
	}
	return dMass;
}

double IsotopeDistribution::getAverageMass()
{
	double dSumProb = 0;
	double dSumMass = 0;
	for (unsigned int i = 0; i < aSize; ++i)
	{
		dSumProb = dSumProb + probArray[i];
		dSumMass = dSumMass + probArray[i] * massArray[i];
	}

	if (dSumProb <= 0)
		return 1.0;

	return (dSumMass / dSumProb);
}

void IsotopeDistribution::filterProbCutoff(double dProbCutoff)
{
	int iteBegin = 0;
	while (iteBegin != aSize)
	{
		if (probArray[iteBegin] > dProbCutoff)
		{
			break;
		}
		iteBegin++;
	}

	int iteEnd = aSize;
	while (iteEnd != 0)
	{
		if (probArray[iteEnd - 1] > dProbCutoff)
		{
			break;
		}
		iteEnd--;
	}
	if (iteBegin != 0 || iteEnd != aSize)
	{
		aSize = iteEnd - iteBegin;
		double *massArrayNew = new double[aSize];
		double *probArrayNew = new double[aSize];
		copy(massArray + iteBegin, massArray + iteEnd, massArrayNew);
		copy(probArray + iteBegin, probArray + iteEnd, probArrayNew);
		delete[] massArray;
		delete[] probArray;
		massArray = massArrayNew;
		probArray = probArrayNew;
	}
}

void IsotopeDistribution::normalize()
{
	double sumProb = 0;
	for (int i = 0; i < aSize; ++i)
	{
		sumProb += probArray[i];
	}

	if (sumProb <= 0)
	{
		cerr << "Error: Sum of distribution is zero" << endl;
		exit(1);
	}

	for (int i = 0; i < aSize; ++i)
	{
		probArray[i] = probArray[i] / sumProb;
	}
}

double IsotopeDistribution::getLowestMass()
{
	return *min_element(massArray, massArray + aSize);
}

Isotopologue::Isotopologue() : MassPrecision(0.01), ProbabilityCutoff(0.000000001)
{
}

Isotopologue::~Isotopologue()
{
	// destructor
}

bool Isotopologue::setupIsotopologue(const string &sTable, const string &AtomNameInput)
{
	// CHONPS
	AtomName = AtomNameInput;
	AtomNumber = AtomName.size();

	istringstream issStream(sTable);
	string sResidue;
	vector<int> viAtomVector;
	int iNumber;
	unsigned int i;

	// parse out the RESIDUE_ATOMIC_COMPOSITION table
	while (!(issStream.eof()))
	{
		// each row is expected to start with the residue name, following by 6 numbers for the natuual CHONPS
		issStream >> sResidue;
		if (sResidue == "")
			continue;
		viAtomVector.clear();
		viAtomVector.reserve(AtomNumber);
		for (i = 0; i < AtomNumber; ++i)
		{
			if (issStream.eof())
			{
				// this row doesn't have 6 fields
				cerr << "ERROR:  the RESIDUE_ATOMIC_COMPOSITION table in ProNovoConfig is not correct!" << endl;
				return false;
			}
			issStream >> iNumber;
			viAtomVector.push_back(iNumber);
		}
		// add this row into the mResidueAtomicComposition table
		mResidueAtomicComposition[sResidue] = viAtomVector;
		sResidue = "";
	}

	// push 6 empty IsotopeDistributions into vAtomIsotopicDistribution
	vAtomIsotopicDistribution.reserve(AtomNumber);
	for (i = 0; i < (AtomNumber); ++i)
	{
		IsotopeDistribution TempDistribution;
		vAtomIsotopicDistribution.push_back(TempDistribution);
	}

	// variables to be passed as reference to ProNovoConfig::getAtomIsotopicComposition
	// to receive its return value
	vector<double> vdMassTemp;
	vector<double> vdNaturalCompositionTemp;

	// the isotopic distribution is pushed into vAtomIsotopicDistribution in the order of
	// natural CHONPS
	for (i = 0; i < AtomName.size(); ++i)
	{
		if (!ProNovoConfig::getAtomIsotopicComposition(
				AtomName[i],
				vdMassTemp,
				vdNaturalCompositionTemp))
		{
			cerr << "ERROR: cannot retrieve isotopic composition for atom " << AtomName[i] << " from ProNovoConfig" << endl;
			return false;
		}
		vAtomIsotopicDistribution[i] = IsotopeDistribution(vdMassTemp,vdNaturalCompositionTemp);
	}

	// calculate Isotopic distribution for all residues
	map<string, vector<int>>::iterator ResidueIter;
	IsotopeDistribution tempIsotopeDistribution;
	for (ResidueIter = mResidueAtomicComposition.begin(); ResidueIter != mResidueAtomicComposition.end(); ResidueIter++)
	{
		if (!computeIsotopicDistribution(ResidueIter->second, tempIsotopeDistribution))
		{
			cerr << "ERROR: cannot calculate the isotopic distribution for residue " << ResidueIter->first << endl;
			return false;
		}

		vResidueIsotopicDistribution[ResidueIter->first] = tempIsotopeDistribution;

		//		cout << "Residue " << ResidueIter->first << endl;
		//		tempIsotopeDistribution.print();
	}

	return true;
}

double Isotopologue::computeMostAbundantMass(string sSequence)
{
	IsotopeDistribution tempIsotopeDistribution;
	if (!computeIsotopicDistribution(sSequence, tempIsotopeDistribution))
		return 0;
	else
		return tempIsotopeDistribution.getMostAbundantMass();
}

double Isotopologue::computeAverageMass(string sSequence)
{
	IsotopeDistribution tempIsotopeDistribution;
	if (!computeIsotopicDistribution(sSequence, tempIsotopeDistribution))
		return 0;
	else
		return tempIsotopeDistribution.getAverageMass();
}

double Isotopologue::computeMonoisotopicMass(string sSequence)
{
	IsotopeDistribution tempIsotopeDistribution;
	if (!computeIsotopicDistribution(sSequence, tempIsotopeDistribution))
		return 0;
	else
		return tempIsotopeDistribution.getLowestMass();
}

bool Isotopologue::getSingleResidueMostAbundantMasses(vector<string> &vsResidues, vector<double> &vdMostAbundantMasses, double &dTerminusMassN, double &dTerminusMassC)
{
	vsResidues.clear();
	vdMostAbundantMasses.clear();

	map<string, IsotopeDistribution>::iterator ResidueIter;
	string sCurrentResidue;
	IsotopeDistribution currentDistribution;
	double dCurrentMostAbundantMasses;

	// for single amino acid

	for (ResidueIter = vResidueIsotopicDistribution.begin(); ResidueIter != vResidueIsotopicDistribution.end(); ResidueIter++)
	{
		sCurrentResidue = ResidueIter->first;
		currentDistribution = ResidueIter->second;
		dCurrentMostAbundantMasses = currentDistribution.getMostAbundantMass();
		if (sCurrentResidue.size() == 1)
		{
			vsResidues.push_back(sCurrentResidue);
			vdMostAbundantMasses.push_back(dCurrentMostAbundantMasses);
			//			cout << sCurrentResidue << "  " << dCurrentMostAbundantMasses << endl;
		}
		else if (sCurrentResidue == "NTerm" || sCurrentResidue == "Nterm")
		{
			dTerminusMassN = dCurrentMostAbundantMasses;
		}
		else if (sCurrentResidue == "CTerm" || sCurrentResidue == "Cterm")
		{
			dTerminusMassC = dCurrentMostAbundantMasses;
		}
		else
		{
			cerr << "ERROR: Cannot recognize the configuration for residue " << sCurrentResidue << endl;
		}
	}

	unsigned int i;

	// bubble sort the list by mass
	unsigned int n = vsResidues.size();
	unsigned int pass;
	double dCurrentMass;
	for (pass = 1; pass < n; pass++)
	{ // count how many times
		// This next loop becomes shorter and shorter
		for (i = 0; i < n - pass; i++)
		{
			if (vdMostAbundantMasses[i] > vdMostAbundantMasses[i + 1])
			{
				// exchange
				dCurrentMass = vdMostAbundantMasses[i];
				sCurrentResidue = vsResidues[i];

				vdMostAbundantMasses[i] = vdMostAbundantMasses[i + 1];
				vsResidues[i] = vsResidues[i + 1];

				vdMostAbundantMasses[i + 1] = dCurrentMass;
				vsResidues[i + 1] = sCurrentResidue;
			}
		}
	}

	return true;
}

bool Isotopologue::computeIsotopicDistribution(string sSequence, IsotopeDistribution &myIsotopeDistribution)
{
	IsotopeDistribution sumDistribution;
	IsotopeDistribution currentDistribution;
	map<string, IsotopeDistribution>::iterator ResidueIter;

	ResidueIter = vResidueIsotopicDistribution.find("Nterm");
	if (ResidueIter != vResidueIsotopicDistribution.end())
	{
		currentDistribution = ResidueIter->second;
		sumDistribution = currentDistribution;
	}
	else
	{
		cerr << "ERROR: can't find the N-terminus" << endl;
		return false;
	}

	ResidueIter = vResidueIsotopicDistribution.find("Cterm");
	if (ResidueIter != vResidueIsotopicDistribution.end())
	{
		currentDistribution = ResidueIter->second;
		sumDistribution = sum(currentDistribution, sumDistribution);
	}
	else
	{
		cerr << "ERROR: can't find the C-terminus" << endl;
		return false;
	}

	// add up all residues's isotopic distribution
	for (unsigned int j = 0; j < sSequence.length(); j++)
	{
		string currentResidue = sSequence.substr(j, 1);
		ResidueIter = vResidueIsotopicDistribution.find(currentResidue);
		if (ResidueIter != vResidueIsotopicDistribution.end())
		{
			currentDistribution = ResidueIter->second;
			sumDistribution = sum(currentDistribution, sumDistribution);
		}
		else
		{
			cerr << "ERROR: can't find the residue " << currentResidue << endl;
			return false;
		}
	}

	myIsotopeDistribution = sumDistribution;

	return true;
}

bool Isotopologue::computeProductIon(string sSequence,
									 vector<vector<double>> &vvdYionMass,
									 vector<vector<double>> &vvdYionProb,
									 vector<vector<double>> &vvdBionMass,
									 vector<vector<double>> &vvdBionProb)
{
	// get the mass for a proton
	double dProtonMass = ProNovoConfig::getProtonMass();
	unsigned int i = 0;

	//	if(!isalpha(sSequence[0]))
	//	{
	//		cout << "ERROR: First character in a peptide sequence can't be a PTM." << endl;
	//		return false;
	//	}

	int iPeptideLength = 0;
	for (i = 0; i < sSequence.length(); ++i)
	{
		if (isalpha(sSequence[i]))
		{
			iPeptideLength = iPeptideLength + 1;
		}
	}

	if (iPeptideLength < ProNovoConfig::getMinPeptideLength())
	{
		cerr << "ERROR: Peptide sequence is too short " << sSequence << endl;
		return false;
	}

	//	cout << "sSequence = " << sSequence << endl;

	vvdYionMass.clear();
	vvdYionProb.clear();
	vvdBionMass.clear();
	vvdBionProb.clear();

	vvdYionMass.reserve(iPeptideLength);
	vvdYionProb.reserve(iPeptideLength);
	vvdBionMass.reserve(iPeptideLength);
	vvdBionProb.reserve(iPeptideLength);

	map<string, IsotopeDistribution>::iterator ResidueIter;

	vector<IsotopeDistribution> vResidueDistribution;
	IsotopeDistribution sumDistribution;
	IsotopeDistribution currentDistribution;

	string currentResidue;
	string currentPTM;

	if (sSequence[0] != '[')
	{
		cerr << "ERROR: First character in a peptide sequence must be [." << endl;
		return false;
	}

	unsigned int iStartResidueIndex = 1;
	ResidueIter = vResidueIsotopicDistribution.find("Nterm");
	if (ResidueIter != vResidueIsotopicDistribution.end())
	{
		currentDistribution = ResidueIter->second;

		if (!isalpha(sSequence[1]))
		{
			iStartResidueIndex = 2;
			currentPTM = sSequence[1];
			ResidueIter = vResidueIsotopicDistribution.find(currentPTM);
			if (ResidueIter == vResidueIsotopicDistribution.end())
			{
				cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
				return false;
			}

			currentDistribution = sum(currentDistribution, ResidueIter->second);
		}

		vResidueDistribution.push_back(currentDistribution);
	}
	else
	{
		cerr << "ERROR: can't find the N-terminus" << endl;
		return false;
	}

	for (i = iStartResidueIndex; i < sSequence.length(); i++)
	{
		if (sSequence[i] == ']')
		{
			break;
		}

		if (!isalpha(sSequence[i]))
		{
			cerr << "ERROR: One residue can only have one PTM (Up to only one symbol after an amino acid)" << endl;
			return false;
		}
		currentResidue = sSequence.substr(i, 1);
		ResidueIter = vResidueIsotopicDistribution.find(currentResidue);
		if (ResidueIter == vResidueIsotopicDistribution.end())
		{
			cerr << "ERROR: cannot find this residue in the config file. " << currentResidue << endl;
			return false;
		}
		currentDistribution = ResidueIter->second;
		if (i + 1 < sSequence.length())
		{
			if (!isalpha(sSequence[i + 1]) && sSequence[i + 1] != ']') // this residue is modified
			{
				currentPTM = sSequence[i + 1];
				ResidueIter = vResidueIsotopicDistribution.find(currentPTM);
				if (ResidueIter == vResidueIsotopicDistribution.end())
				{
					cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
					return false;
				}

				// this PTM can substract mass from the residue
				// the currentDistribution could even be negative.
				currentDistribution = sum(currentDistribution, ResidueIter->second);

				i = i + 1; // increment index to the next residue
			}
		}

		vResidueDistribution.push_back(currentDistribution);
	}

	ResidueIter = vResidueIsotopicDistribution.find("Cterm");
	if (ResidueIter != vResidueIsotopicDistribution.end())
	{
		currentDistribution = ResidueIter->second;

		if (i + 1 < sSequence.length())
		{
			if (!isalpha(sSequence[i + 1]))
			{
				currentPTM = sSequence[i + 1];
				ResidueIter = vResidueIsotopicDistribution.find(currentPTM);
				if (ResidueIter == vResidueIsotopicDistribution.end())
				{
					cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
					return false;
				}

				currentDistribution = sum(currentDistribution, ResidueIter->second);
			}
		}

		vResidueDistribution.push_back(currentDistribution);
	}
	else
	{
		cerr << "ERROR: can't find the C-terminus" << endl;
		return false;
	}

	// compute B-ion series

	vector<IsotopeDistribution> vBionDistribution;
	vBionDistribution.reserve(iPeptideLength);

	// start with the N-terminus distibution, which should the first element
	currentDistribution = vResidueDistribution.front();
	sumDistribution = currentDistribution;

	int j;
	for (j = 1; j < iPeptideLength; j++)
	{
		currentDistribution = vResidueDistribution[j];
		sumDistribution = sum(currentDistribution, sumDistribution);
		vBionDistribution.push_back(sumDistribution);
	}
	for (i = 0; i < vBionDistribution.size(); i++)
	{
		vvdBionMass.push_back(vector<double>(vBionDistribution[i].massArray,
											 vBionDistribution[i].massArray + vBionDistribution[i].aSize));
		vvdBionProb.push_back(vector<double>(vBionDistribution[i].probArray,
											 vBionDistribution[i].probArray + vBionDistribution[i].aSize));
		//	cout << "b " << i << endl;
		//	vBionDistribution[i].print();
	}

	// compute Y-ion series
	vector<IsotopeDistribution> vYionDistribution;
	vYionDistribution.reserve(iPeptideLength);

	// start with the C-terminus distibution, which should the last element
	currentDistribution = vResidueDistribution.back();
	sumDistribution = currentDistribution;

	for (j = iPeptideLength; j > 1; j--)
	{
		currentDistribution = vResidueDistribution[j];
		sumDistribution = sum(currentDistribution, sumDistribution);
		vYionDistribution.push_back(sumDistribution);
	}

	for (i = 0; i < vYionDistribution.size(); i++)
	{
		vvdYionMass.push_back(vector<double>(vYionDistribution[i].massArray,
											 vYionDistribution[i].massArray + vYionDistribution[i].aSize));
		vvdYionProb.push_back(vector<double>(vYionDistribution[i].probArray,
											 vYionDistribution[i].probArray + vYionDistribution[i].aSize));
		//	cout << "y " << i << endl;
		//	vYionDistribution[i].print();
	}

	// change the masses to correct for the proton transfer during peptide bond cleavage
	unsigned int n;
	unsigned int m;
	for (n = 0; n < vvdYionMass.size(); ++n)
	{
		for (m = 0; m < vvdYionMass[n].size(); ++m)
		{
			vvdYionMass[n][m] += dProtonMass;
		}
	}

	for (n = 0; n < vvdBionMass.size(); ++n)
	{
		for (m = 0; m < vvdBionMass[n].size(); ++m)
		{
			vvdBionMass[n][m] -= dProtonMass;
		}
	}

	// cout << "Y		B" << endl;
	// for (n = 0; n < vvdBionMass.size(); ++n)
	// {
	// 	cout << vvdYionMass[n][0] << "\t\t" << vvdBionMass[n][0] << endl;
	// }

	return true;
}

bool Isotopologue::computeIsotopicDistribution(vector<int> AtomicComposition, IsotopeDistribution &myIsotopeDistribution)
{
	IsotopeDistribution sumDistribution;
	IsotopeDistribution currentAtomDistribution;
	currentAtomDistribution = multiply(vAtomIsotopicDistribution[0], AtomicComposition[0]);
	sumDistribution = currentAtomDistribution;

	for (unsigned int i = 1; i < AtomNumber; i++)
	{
		currentAtomDistribution = multiply(vAtomIsotopicDistribution[i], AtomicComposition[i]);
		sumDistribution = sum(currentAtomDistribution, sumDistribution);
	}
	myIsotopeDistribution = sumDistribution;

	return true;
}

bool Isotopologue::computeAtomicComposition(string sSequence, vector<int> &myAtomicComposition)
{
	vector<int> AtomicComposition;
	vector<int> CurrentComposition;
	unsigned int i;
	map<string, vector<int>>::iterator ResidueIter;

	for (i = 0; i < AtomNumber; i++)
		AtomicComposition.push_back(0);

	ResidueIter = mResidueAtomicComposition.find("Nterm");
	if (ResidueIter != mResidueAtomicComposition.end())
	{
		CurrentComposition = ResidueIter->second;
		for (i = 0; i < AtomNumber; i++)
			AtomicComposition[i] = CurrentComposition[i];
	}
	else
	{
		cerr << "ERROR: can't find the atomic composition for the N-terminus" << endl;
		return false;
	}

	ResidueIter = mResidueAtomicComposition.find("Cterm");
	if (ResidueIter != mResidueAtomicComposition.end())
	{
		CurrentComposition = ResidueIter->second;
		for (i = 0; i < AtomNumber; i++)
			AtomicComposition[i] += CurrentComposition[i];
	}
	else
	{
		cerr << "ERROR: can't find the atomic composition for the C-terminus" << endl;
		return false;
	}

	for (unsigned int j = 0; j < sSequence.length(); j++)
	{
		string currentResidue = sSequence.substr(j, 1);
		ResidueIter = mResidueAtomicComposition.find(currentResidue);
		if (ResidueIter != mResidueAtomicComposition.end())
		{
			CurrentComposition = ResidueIter->second;
			for (i = 0; i < AtomNumber; i++)
				AtomicComposition[i] += CurrentComposition[i];
		}
		else
		{
			cerr << "ERROR: can't find the atomic composition for residue/PTM: " << currentResidue << endl;
			return false;
		}
	}

	myAtomicComposition = AtomicComposition;
	return true;
}

IsotopeDistribution Isotopologue::sum(const IsotopeDistribution &distribution0, const IsotopeDistribution &distribution1)
{
	double ProbabilityCutoff_local = 0.000001;

	double currentMass;
	double currentProb;
	int iSizeDistribution0 = distribution0.aSize;
	int iSizeDistribution1 = distribution1.aSize;
	int iCount = 0;
	double dSum = 0;
	IsotopeDistribution sumDistribution(iSizeDistribution0 + iSizeDistribution1 - 1);
	for (int k = 0; k < iSizeDistribution0 + iSizeDistribution1 - 1; k++)
	{
		double sumweight = 0, summass = 0;
		int start = k < (iSizeDistribution1 - 1) ? 0 : k - iSizeDistribution1 + 1; // max(0, k-f_n+1)
		int end = k < (iSizeDistribution0 - 1) ? k : iSizeDistribution0 - 1;	   // min(g_n - 1, k)
		iCount = 0;
		dSum = 0;
		for (int i = start; i <= end; i++)
		{
			double weight = distribution0.probArray[i] * distribution1.probArray[k - i];
			double mass = distribution0.massArray[i] + distribution1.massArray[k - i];
			sumweight += weight;
			summass += weight * mass;
			iCount += 1;
			dSum += mass;
		}
		if (sumweight == 0)
		{
			currentMass = dSum / ((double)iCount);
		}
		else
		{
			currentMass = summass / sumweight;
		}
		currentProb = sumweight;
		// sumDistribution.vMass.push_back(currentMass);
		// sumDistribution.vProb.push_back(currentProb);
		sumDistribution.massArray[k] = currentMass;
		sumDistribution.probArray[k] = currentProb;
	}
	int iSizeSumDistribution;
	int i;

	//prune small probabilities
	sumDistribution.filterProbCutoff(ProbabilityCutoff_local);

	// normalize the probability space to 1
	sumDistribution.normalize();

	return sumDistribution;
}

IsotopeDistribution Isotopologue::multiply(const IsotopeDistribution &distribution0, int count)
{
	if (count == 1)
		return distribution0;
	IsotopeDistribution productDistribution(1);
	productDistribution.massArray[0] = 0.0;
	productDistribution.probArray[0] = 1.0;
	if (count < 0)
	{
		IsotopeDistribution negativeDistribution = distribution0;
		for (unsigned int n = 0; n < negativeDistribution.aSize; n++)
		{
			negativeDistribution.massArray[n] = -negativeDistribution.massArray[n];
		}

		for (int i = 0; i < abs(count); ++i)
		{
			productDistribution = sum(productDistribution, negativeDistribution);
		}
	}
	else
	{
		for (int i = 0; i < count; ++i)
		{
			productDistribution = sum(productDistribution, distribution0);
		}
	}
	return productDistribution;
}

void Isotopologue::shiftMass(IsotopeDistribution &distribution0, double dMass)
{
	for (unsigned int i = 0; i < distribution0.aSize; ++i)
	{
		distribution0.massArray[i] = distribution0.massArray[i] + dMass;
	}
}
