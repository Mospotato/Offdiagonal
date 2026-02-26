#define STNegativeBinomial_cxx
#include <assert.h>
#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TH3.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"

#include "StNegativeBinomial.h"
#include "StNbdFitMaker.h"
//____________________________________________________________________________________________________
// Default constructor
StNbdFitMaker::StNbdFitMaker()
{
	mNBinomial = 0;
	hRefMult = 0;
	MinChi2Ndof = 1e9;
	mNData = 0;
	OutputName = "Debug.root";
	TH1::SetDefaultSumw2();
	mMinimumMultiplicityCut = 100.0; // >50 by default. Can be changed by setMinimumMultiplicityCut(const Double_t)
}

//____________________________________________________________________________________________________
// Default destructor
StNbdFitMaker::~StNbdFitMaker()
{
	if (mNBinomial)
		delete mNBinomial;
}

//____________________________________________________________________________________________________
Double_t StNbdFitMaker::GetNormalization(const TH1 &h1, const TH1 &h2, const Double_t min, const Double_t max) const
{
	// Get normalization factor in (min, max)
	Double_t numerator = 0.0;
	Double_t denominator = 0.0;
	const Int_t mulminbin = h1.GetXaxis()->FindBin(min);
	const Int_t mulmaxbin = h1.GetXaxis()->FindBin(max);

	for (Int_t mul = mulminbin; mul < mulmaxbin; mul++)
	{
		const Double_t n1 = h1.GetBinContent(mul + 1);
		const Double_t n1Error = h1.GetBinError(mul + 1);

		if (n1 == 0.0 || n1Error == 0.0)
			continue;

		const Double_t n2 = h2.GetBinContent(mul + 1);

		numerator += n1 * n2 / (n1Error * n1Error);
		denominator += n2 * n2 / (n1Error * n1Error);
	}

	const Double_t norm = (denominator != 0.0) ? numerator / denominator : 1.0;

	// printf("StNbdFitMaker::GetNormalization  (min, max, norm) = (%5d, %5d, %1.7f)\n", (Int_t)min, (Int_t)max, norm);
	return norm;
}

//____________________________________________________________________________________________________
Double_t StNbdFitMaker::CalculateChi2(const TH1 &hdata, const TH1 &hfunc, const Double_t minimumMultiplicityCut)
{
	/// Calculate chi2 from data and func
	mNData = 0;

	Double_t chi2 = 0.0;
	for (Int_t ix = 1; ix <= hdata.GetNbinsX(); ix++)
	{
		// Lower multiplicity cut off
		const Double_t mult = hdata.GetBinCenter(ix);
		if (mult < minimumMultiplicityCut)
			continue;

		// Check data points
		const Double_t dataError = hdata.GetBinError(ix);
		if (dataError == 0.0)
			continue;
		// Check sim points
		const Double_t simError = hfunc.GetBinError(ix);
		if (simError == 0.0)
			continue;
		const Double_t data = hdata.GetBinContent(ix);
		if (data == 0.0)
		{
			continue;
		}
		
		// Combine errors in quadrature
		const Double_t totError = TMath::Power((dataError * dataError + simError * simError), 0.5);

		// Calculate chi2
		const Double_t func = hfunc.GetBinContent(ix);
		const Double_t delta = (data - func) / totError;
		chi2 += delta * delta;
		mNData++;
	}

	printf("StNbdFitMaker::CalculateChi2: (eff, npp, x, k) = (%.3f, %1.5f, %1.5f, %1.5f)  chi2/ndof = %1.3f/%d = %1.3f\n",mNBinomial->GetEfficiency(), mNBinomial->GetNpp(), mNBinomial->GetX(), mNBinomial->GetK(), chi2, mNData - 2, chi2 / static_cast<Double_t>(mNData - 2.));
	if (std::isnan(chi2) || std::isinf(chi2))
	{
		return 1e9;
	}
	return chi2;
}

//____________________________________________________________________________________________________
void StNbdFitMaker::SetParameters(const Double_t npp, const Double_t k, const Double_t x,
								  const Double_t efficiency, const Double_t triggerbias, const Bool_t isConstEfficiency)
{
	if (mNBinomial)
		delete mNBinomial;

	mNBinomial = new StNegativeBinomial(npp, k, x, efficiency, triggerbias, isConstEfficiency);
}

//____________________________________________________________________________________________________
void StNbdFitMaker::SetMinimumMultiplicityCut(const Double_t cut)
{
	mMinimumMultiplicityCut = cut;
	printf("StNbdFitMaker::setMinimumMultiplicityCut Set low multiplicity cut off : M > %f\n", cut);
}

//____________________________________________________________________________________________________
void StNbdFitMaker::ReadData(TFile *fData, TFile *fGlauber)
{
	hRefMult = (TH1F *)fData->Get("RefMult");
	if (!hRefMult)
	{
		Error("StNbdFitMaker::readData", "RefMult doesn't exist");
		assert(hRefMult);
	}
	hRefMult->SetLineColor(1);
	hNcollNpart = (TH2D *)fGlauber->Get("hNcoll_Npart");
	if (!hNcollNpart)
	{
		Error("StNbdFitMaker::readData", "hNcoll_Npart doesn't exist");
		assert(hNcollNpart);
	}
}

void StNbdFitMaker::Fit(const Int_t nevents, TString name) // zaochen
{
	hRefMultSim.emplace_back(new TH1D(name.Data(), "", hRefMult->GetNbinsX(), hRefMult->GetXaxis()->GetXmin(), hRefMult->GetXaxis()->GetXmax()));
	TH1 *mhRefMultSim = hRefMultSim.back();

	Int_t ievent = 0;
	while (ievent < nevents)
	{
		Double_t npart, ncoll;
		hNcollNpart->GetRandom2(npart, ncoll);
		const Bool_t isNpartNcollOk = (npart >= 2 && ncoll >= 1);
		if (!isNpartNcollOk)
			continue;

		const Int_t multiplicity = mNBinomial->GetMultiplicity(npart, static_cast<Int_t>(ncoll));
		mhRefMultSim->Fill(multiplicity);
		ievent++;
	}

	// Normalization
	const Double_t norm = GetNormalization(*hRefMult, *mhRefMultSim, mMinimumMultiplicityCut, hRefMult->GetXaxis()->GetXmax());
	mhRefMultSim->Scale(norm);

	// Get chi2
	const Double_t Chi2overNdof = CalculateChi2(*hRefMult, *mhRefMultSim, mMinimumMultiplicityCut) / (mNData - 2.0);
	if (fabs(Chi2overNdof - 1.) < fabs(MinChi2Ndof - 1.))
	{
		MinChi2Ndof = Chi2overNdof;
		MinChi2Index = hRefMultSim.size() - 1;
		MinChi2Name = name;
	}

	mhRefMultSim->SetXTitle("Refmult (MC)");
	mhRefMultSim->SetYTitle("Count");
	mhRefMultSim->SetTitle(name.Data());
}

//____________________________________________________________________________________________________
Int_t StNbdFitMaker::Scan(const Int_t nevents,
						  const Int_t nppbin, const Double_t nppmin, const Double_t nppmax,
						  const Int_t kbin, const Double_t kmin, const Double_t kmax,
						  const Int_t xbin, const Double_t xmin, const Double_t xmax,
						  const Int_t effbin, const Double_t effmin, const Double_t effmax,
						  const Double_t triggerbias, const Bool_t isConstEfficiency)
{
	OutputFile = TFile::Open(OutputName.c_str(), "recreate");
	OutputFile->mkdir("Simulations");
	OutputFile->cd("Simulations");
	if (!OutputFile || !OutputFile->IsOpen())
	{
		Error("StNbdFitMaker::Scan", "Cannot open output file");
		assert(OutputFile);
	}
	if (!hRefMult)
	{
		Error("StNbdFitMaker::Fit", "hRefMult doesn't exist");
		assert(hRefMult);
	}
	if (!hNcollNpart)
	{
		Error("StNbdFitMaker::Fit", "hNcoll_Npart doesn't exist");
		assert(hNcollNpart);
	}
	const Double_t effstep = (effbin == 1) ? 0 : (effmax - effmin) / static_cast<Double_t>(effbin - 1);
	const Double_t nppstep = (nppbin == 1) ? 0 : (nppmax - nppmin) / static_cast<Double_t>(nppbin - 1);
	const Double_t kstep = (kbin == 1) ? 0 : (kmax - kmin) / static_cast<Double_t>(kbin - 1);
	const Double_t xstep = (xbin == 1) ? 0 : (xmax - xmin) / static_cast<Double_t>(xbin - 1);
	printf("StNbdFitMaker::Scan  eff: %d, effStep = %1.3f, effMin = %1.3f, effMax = %1.3f\n", effbin, effstep, effmin, effmax);
	printf("StNbdFitMaker::Scan  nBin: %d, nppStep = %1.3f, nppMin = %1.3f, nppMax = %1.3f\n", nppbin, nppstep, nppmin, nppmax);
	printf("StNbdFitMaker::Scan  kBin: %d, kStep = %1.3f, kMin = %1.3f, kMax = %1.3f\n", kbin, kstep, kmin, kmax);
	printf("StNbdFitMaker::Scan  xBin: %d, xStep = %1.3f, xMin = %1.3f, xMax = %1.3f\n", xbin, xstep, xmin, xmax);
	for (Int_t ie = 0; ie < effbin; ie++)
	{
		const Double_t efficiency = effmin + effstep * ie;
		for (Int_t ix = 0; ix < nppbin; ix++)
		{
			const Double_t npp = nppmin + nppstep * ix;
			for (Int_t iy = 0; iy < kbin; iy++)
			{
				const Double_t k = kmin + kstep * iy;
				for (Int_t iz = 0; iz < xbin; iz++)
				{
					const Double_t x = xmin + xstep * iz;
					SetParameters(npp, k, x, efficiency, triggerbias, isConstEfficiency);
					Fit(nevents, Form("npp%1.5f_k%1.5f_x%1.5f_eff%1.5f", npp, k, x, efficiency));
				}
			}
		}
	}

	printf("\nFinish!\nStNbdFitMaker::Scan  Minimum chi2/ndof = %1.3f for %s\n", MinChi2Ndof, MinChi2Name.c_str());
	mNBinomial->GetNbd()->SetDirectory(OutputFile);
	OutputFile->Write();
	OutputFile->cd();
	hRefMult->Write();
	hRefMultSim[MinChi2Index]->SetName("RefMultSim");
	hRefMultSim[MinChi2Index]->Write();

	TH1 *hRatio = (TH1D *)hRefMultSim[MinChi2Index]->Clone();
	hRatio->SetName("hRatio");
	hRatio->Divide(hRefMult);
	hRatio->SetYTitle("MC/data");

	hRatio->SetMinimum(0);
	hRatio->SetMaximum(2);
	hRatio->Write();
	// OutputFile->Close();
	printf("Output file %s is saved\n", OutputName.c_str());
	return 0;
}