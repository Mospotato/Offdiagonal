///////////////////////////////
// Main Program
///////////////////////////////
#include <iostream>
#include <fstream>
#include <stdlib.h> //atoi
#include <unistd.h> //getopt
#include <cstdlib>

#include <TROOT.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>

#include "ReadUrQMD.h"
int main(int argc, char **argv)
{

  // command line parameter defaults
  Int_t Energy = 0;
  std::string filelist = "test.list";
  std::string Jobid = "DEBUG";
  Bool_t FXTFlag = false;
  int opt;
  std::string USAGE = "Usage: ./run -e Energy -f filelist -n jobid\n";
  while ((opt = getopt(argc, argv, "e:f:n:F")) != -1)
  {
    switch (opt)
    {
    case 'e':
      Energy = std::stoi(optarg);
      break;
    case 'f':
      filelist = optarg;
      break;
    case 'n':
      Jobid = optarg;
      break;
    case 'F':
      FXTFlag = true;
      break;
    default:
      std::cerr << USAGE;
      return 1;
    }
  }
  if (!Energy)
  {
    std::cerr << "Please specify the energy with -s option\n";
    return 1;
  }
  std::cout << "===========================================\n           Reading Files                   \n===========================================\n";
  const Char_t *TreeName = "urqmd";
  TChain chain(TreeName);
  Int_t sum = 0;
  std::ifstream list(filelist.c_str());
  std::string InputFile;
  while (std::getline(list, InputFile))
  {
    TFile *tempFile = TFile::Open(InputFile.c_str());
    if (!tempFile || tempFile->IsZombie())
    {
      std::cout << "Invalid: " << InputFile << '\n';
      continue;
    }
    TTree *tempTree;
    tempFile->GetObject(TreeName, tempTree);
    if (tempTree == NULL)
    {
      std::cout << "No Tree is Found in " << InputFile << '\n';
      continue;
    }
    std::cout << "Read in file: " << InputFile << " Entries: " << tempTree->GetEntries() << '\n';
    chain.Add(InputFile.c_str());
    sum++;
    delete tempFile;
  }
  list.close();

  // print out the number of events
  std::cout << "===========================================\n"
            << "       " << sum << " file(s) is(are) merged      \n"
            << "===========================================\n";

  std::cout << "===========================================\n"
            << "           Analyzing  Files                   "
            << "\n===========================================\n";

  Long64_t runevt = chain.GetEntries();
  std::cout << "Entries: " << runevt << std::endl;

  ReadUrQMD reader(Energy, FXTFlag, Jobid);
  chain.Process(&reader, "", runevt);
  return 0;
}
