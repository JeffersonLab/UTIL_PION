#include <TProof.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

void Q0p375W2p2Left1()
{
  TChain ch("T");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8544_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8545_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8546_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8547_-1.root");

  TProof *proof = TProof::Open("workers=8");
  //proof->SetProgressDialog(0);  
  ch.SetProof();
  ch.Process("/home/cdaq/hallc-online/hallc_replay_lt/UTIL_PION/scripts_Yield/PionYield.C+","1");
  proof->Close();
  
}
