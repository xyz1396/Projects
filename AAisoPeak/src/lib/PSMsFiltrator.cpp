#include "PSMsFiltrator.h"

PSMsFiltrator::PSMsFiltrator(const vector<sipPSM> &sipPSMs, float mFDRthreshold)
	: FDRthreshold(mFDRthreshold)
{
	for (sipPSM psm : sipPSMs)
	{
		for (size_t i = 0; i < psm.scanNumbers.size(); i++)
		{
			string psmID = psm.fileName + to_string(psm.scanNumbers[i]);
			auto psmIX = PSMsMap.find(psmID);
			if (psmIX != PSMsMap.end())
				fillPSMinfo(psmIX->second, psm, i);
			else
			{
				PSMinfo temp;
				temp.ftFileName = psm.fileName;
				temp.scanNumber = psm.scanNumbers[i];
				fillPSMinfo(temp, psm, i);
				PSMsMap.insert({psmID, temp});
			}
		}
	}
}

PSMsFiltrator::~PSMsFiltrator()
{
}

void PSMsFiltrator::splitString(const string mString)
{
	string sep = ",";
	size_t start = 0;
	size_t end = mString.find(sep);
	tokens.clear();
	while (end != std::string::npos)
	{
		tokens.push_back(mString.substr(start, end - start));
		start = end + sep.length();
		end = mString.find(sep, start);
	}
	tokens.push_back(mString.substr(start));
}

pair<int, bool> PSMsFiltrator::detectProDecoy(const string &proteinName)
{
	// remove {} out of protein names then split it
	splitString(proteinName.substr(1, proteinName.size() - 2));
	// if one protein is not decoy the peptide is not decoy
	for (string token : tokens)
	{
		if (token.substr(0, 4) != "Rev_")
			return {tokens.size(), false};
	}
	return {tokens.size(), true};
}

void PSMsFiltrator::fillPSMinfo(PSMinfo &mPSMinfo, const sipPSM &psm, const int psmIX)
{
	for (size_t i = 0; i < 5; i++)
	{
		if (psm.scores[psmIX] > mPSMinfo.bestScores[i])
		{
			mPSMinfo.bestScores[i] = psm.scores[psmIX];
			mPSMinfo.calculatedParentMasses[i] = psm.calculatedParentMasses[psmIX];
			mPSMinfo.identifiedPepSeqs[i] = psm.identifiedPeptides[psmIX];
			tie(mPSMinfo.proCounts[i], mPSMinfo.isDecoys[i]) = detectProDecoy(psm.proteinNames[psmIX]);
			mPSMinfo.measuredParentMasses[i] = psm.measuredParentMasses[psmIX];
			mPSMinfo.originalPepSeqs[i] = psm.originalPeptides[psmIX];
			mPSMinfo.parentCharges[i] = psm.parentCharges[psmIX];
			mPSMinfo.searchNames[i] = psm.searchName;
			mPSMinfo.proNames[i] = psm.proteinNames[psmIX];
			mPSMinfo.pepLengths[i] = psm.originalPeptides[psmIX].length() - 4;
			break;
		}
	}
}

void PSMsFiltrator::getRentionTime(string &workPath) {}
void PSMsFiltrator::writePercolatorTSV() {}
void PSMsFiltrator::readPercolatorTSV() {}
void PSMsFiltrator::fillPeptidesCharge() {}
void PSMsFiltrator::sortPeptideAllScore() {}

tuple<size_t, size_t, float> PSMsFiltrator::getDecoyCountScoreThreshold()
{
	return {0, 0, 0};
}

void PSMsFiltrator::filterPSMsMap() {}
void PSMsFiltrator::precursorIsotopicMassesIntensities(string &workPath) {}
void PSMsFiltrator::writeFilteredPSMs() {}