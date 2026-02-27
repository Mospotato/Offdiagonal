#define ReadUrQMD_cxx
#include "ReadUrQMD.h"
#include "../../Dependency/ThreadPool.h"
#include "../../Dependency/Dictionary.h"
#include "../../Dependency/Configuration.h"
#include "EPDSimulator.h"
ReadUrQMD::ReadUrQMD(Int_t Energy, Bool_t flag, std::string jobid)
{
   JobId = jobid;
   FXTFlag = flag;
   auto &Config = Configuration::getInstance();
   Config.SetCentrality(Energy);
   if (FXTFlag)
   {
      Double_t gamma = Energy / (2 * 0.938);
      Double_t beta_abs = sqrt(1 - (1 / pow(gamma, 2)));
      BoostFactor = {0, 0, -beta_abs};
      BeamRapidity = Config.BeamRapidityMap[Energy];
      Info("ReadUrQMD", "FXT Mode: Beam Rapidity: %f", BeamRapidity);
   }
}

void ReadUrQMD::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
   auto &Config = Configuration::getInstance();
   auto *Dict = Dictionary::getInstance();
   //=====================================================================
   OutProfile = new TFile(Form("%s.root", JobId.c_str()), "RECREATE");
   hNpart = new TH1D("Npart", "Npart", 400, 0.5, 400.5);
   hRefMult = new TH1D("RefMult", "Reference Multiplicity", Config.CentVec[0], 0.5, Config.CentVec[0] + 0.5);
   hMultNpart = new TH2D("Mult2Npart", ";RefMult;Npart", Config.CentVec[0], 0.5, Config.CentVec[0] + 0.5, 400, 0.5, 400.5);
   hMultImpact = new TH2D("Mult2Impact", ";RefMult;Impact", Config.CentVec[0], 0.5, Config.CentVec[0] + 0.5, 160, 0, 16);
   ConstituentMap.resize(Config.CentVec[0]);
   for (Int_t RefMult = 0; RefMult < Config.CentVec.front(); RefMult++)
   {
      ConstituentMap[RefMult].resize(Config.nAcceptance);
      for (Int_t iAcc = 0; iAcc < Config.nAcceptance; iAcc++)
      {
         ConstituentMap[RefMult][iAcc].resize(Dict->nConstituent, 0);
      }
   }
   Info("Begin", "Initialization Completed!");
}

void ReadUrQMD::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();
}

Double_t ReadUrQMD::GetRapidity(TLorentzVector &vector)
{
   Double_t rapidity = FXTFlag ? BeamRapidity + vector.Rapidity() : vector.Rapidity();
   return rapidity;
}

Bool_t ReadUrQMD::IsEtaValid(TLorentzVector &vector)
{
   Double_t Eta = vector.Eta();
   if (FXTFlag)
   {
      return Eta > -2. && Eta < 0.;
   }
   else
   {
      return fabs(Eta) < 1.0;
   }
   return false;
}

Bool_t ReadUrQMD::Process(Long64_t entry)
{
   fChain->GetTree()->GetEntry(entry);
   if (entry != 0 && entry % 1000 == 0)
      Info("Process", "Processing Event %lld", entry);
   if (!Npart)
      return kTRUE;
   auto &PDG = PDGData::getInstance();
   auto &Config = Configuration::getInstance();
   // Int_t RefMult = 0;
   TLorentzVector FourMomentum;
   std::vector<std::map<Int_t, Int_t>> ParticleMap(Config.nAcceptance);
   auto &EPDSim = EPDSimulator::getInstance();
   for (Int_t iTrk = 0; iTrk < mul; iTrk++)
   {
      Int_t apid = abs(pid[iTrk]);
      ParticleEntry *ParticleData = PDG.PDGMap[apid];
      if (!ParticleData)
      {
         continue;
      }
      FourMomentum.SetXYZM(px[iTrk], py[iTrk], pz[iTrk], ParticleData->Mass);
      if (FXTFlag)
      {
         FourMomentum.Boost(BoostFactor);
      }
      // if (apid == 211 && IsEtaValid(FourMomentum))
      // {
      //    RefMult++;
      //    continue;
      // }
      EPDSim.FillEPD(FourMomentum);
      Double_t rapidity = GetRapidity(FourMomentum);
      if (!Config.IsAccepted(FourMomentum.Pt(), rapidity))
      {
         continue;
      }
      std::vector<int> AcceptanceVec = Config.GetAcceptanceVec(FourMomentum.Pt(), rapidity);
      Int_t PID = apid == 3212 ? 3122 : apid;
      Int_t sign = pid[iTrk] > 0 ? 1 : -1;
      PID *= sign;
      for (auto iAcc : AcceptanceVec)
      {
         ParticleMap[iAcc][PID]++;
         if (ParticleData->Baryon)
            ParticleMap[iAcc][1] += sign * ParticleData->Baryon;
         if (ParticleData->Strangeness)
            ParticleMap[iAcc][2] += sign * ParticleData->Strangeness;
         if (ParticleData->Charge)
            ParticleMap[iAcc][3] += sign * ParticleData->Charge;
      }
   }
   auto nMIP = EPDSim.GetnMIP();
   hNpart->Fill(Npart);
   hRefMult->Fill(nMIP);
   hMultNpart->Fill(nMIP, Npart);
   hMultImpact->Fill(nMIP, b);
   GetMoment(nMIP, ParticleMap);
   EPDSim.Reset();
   return kTRUE;
}

void ReadUrQMD::GetMoment(Int_t EPDMult, const std::vector<std::map<Int_t, Int_t>> &ParticleMap)
{
   auto &Config = Configuration::getInstance();
   if (EPDMult < Config.CentVec.back() || EPDMult >= Config.CentVec.front())
   {
      return;
   }
   Int_t &nAccpetance = Config.nAcceptance;
#if 0
   auto &pool = ThreadPool::getInstance();
   for (Int_t iAcc = 0; iAcc < nAccpetance; iAcc++)
   {
      pool.enqueue([iAcc, EPDMult, &ParticleMap, this]
                             { 
                              getMoment(true, ParticleMap[iAcc], ConstituentMap[EPDMult][iAcc]);});
   }
   pool.waitForCompletion();
#else
   for (Int_t iAcc = 0; iAcc < nAccpetance; iAcc++)
   {
      getMoment(ParticleMap[iAcc], ConstituentMap[EPDMult][iAcc]);
   }
#endif
   return;
}

void ReadUrQMD::getMoment(const std::map<Int_t, Int_t> &ParticleMap, std::vector<Double_t> &Moment)
{
   auto *Dict = Dictionary::getInstance();
   auto &constituent = Dict->Constituent;
   Int_t ibin = 0;
   for (const auto &component : constituent)
   {
      Double_t moment = 1.;
      for (const auto &pair : component.Internal)
      {
         if (ParticleMap.find(pair.first) == ParticleMap.end())
         {
            moment = 0;
            break;
         }
         moment *= std::pow(ParticleMap.at(pair.first), pair.second);
      }
      Moment[ibin] += moment;
      ibin++;
   }
}

void ReadUrQMD::SlaveTerminate()
{
}

void ReadUrQMD::Terminate()
{
   Info("Terminate", "Writing Files!");
   auto &Config = Configuration::getInstance();
   auto *Dict = Dictionary::getInstance();
   Int_t &nAcceptance = Config.nAcceptance;
   auto pComponent = new TProfile3D("Moment", "", Config.CentVec[0], 0.5, Config.CentVec[0] + 0.5, Dict->nConstituent, 0.5, Dict->nConstituent + 0.5, Config.nAcceptance, -0.5, Config.nAcceptance - 0.5);
   for (Int_t EPDMult = 0; EPDMult < Config.CentVec.front(); EPDMult++)
   {
      Double_t nEntries = hRefMult->GetBinContent(EPDMult);
      for (Int_t iAcc = 0; iAcc < nAcceptance; iAcc++)
      {
         for (Int_t iComponent = 0; iComponent < Dict->nConstituent; iComponent++)
         {
            Int_t binNumber = pComponent->GetBin(EPDMult, iComponent + 1, iAcc + 1);
            pComponent->SetBinContent(binNumber, ConstituentMap[EPDMult][iAcc][iComponent]);
            pComponent->SetBinEntries(binNumber, nEntries);
         }
      }
   }
   OutProfile->Write();
   delete OutProfile;
   return;
}