/******************************************************************************
 * $Id: StNbdFitMaker.h,v 1.2 2012/04/25 05:15:30 hmasui Exp $
 * $Log: StNbdFitMaker.h,v $
 * Revision 1.2  2012/04/25 05:15:30  hmasui
 * Added centrality calculation. Real data histogram name can be now specified in ReadData() function
 *
 ******************************************************************************/
#ifndef StNbdFitMaker_HH
#define StNbdFitMaker_HH
#include "TH2.h"
#include "StNegativeBinomial.h"
class StNbdFitMaker
{
public:
	StNbdFitMaker();
	virtual ~StNbdFitMaker(); /// Default destructor
	void ReadData(TFile *, TFile *);
	void Fit(const Int_t nevents = 1000, const TString outputFileName = "");
	/// Find minimum chi2/NDF in (npp, k, efficiency) space
	Int_t Scan(const Int_t nevents,
			   const Int_t nppbin, const Double_t nppmin, const Double_t nppmax,
			   const Int_t kbin, const Double_t kmin, const Double_t kmax,
			   const Int_t xbin, const Double_t xmin, const Double_t xmax,
			   const Int_t effbin, const Double_t effmin, const Double_t effmax,
			   const Double_t triggerbias = 1.0, const Bool_t isConstEfficiency = kTRUE);
	/// Set parameters
	void SetParameters(const Double_t npp, const Double_t k, const Double_t x,
					   const Double_t efficiency, const Double_t triggerbias, const Bool_t isConstEfficiency);
	/// Set minimum multiplicity cuts to avoid inefficiency (default is M>50)
	void SetMinimumMultiplicityCut(const Double_t cut);
	void SetOutputName(const std::string &name) { OutputName = name; };

private:
	TFile *OutputFile;
	std::string OutputName;
	Double_t MinChi2Ndof;
	Int_t MinChi2Index;
	std::string MinChi2Name;
	Double_t GetNormalization(const TH1 &h1, const TH1 &h2,
							  const Double_t min, const Double_t max) const; // Get normalization factor in (min, max)
	Double_t CalculateChi2(const TH1 &hdata, const TH1 &hfunc,
						   const Double_t minimumMultiplicityCut); // Get chi2 from data and func
	// Data members
	StNegativeBinomial *mNBinomial; /// Negative binomial distribution
	TH1 *hRefMult;
	TH2 *hNcollNpart;
	std::vector<TH1 *> hRefMultSim;
	Int_t mNData;					  /// Number of data points used in fit
	Double_t mMinimumMultiplicityCut; /// Minimum multiplicity cut for fitting
};
#endif
