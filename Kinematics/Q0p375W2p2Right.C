#include <TProof.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

void Q0p375W2p2Right()
{
  TChain ch("T");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8704_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8705_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8706_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8707_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8708_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8709_-1.root");

   
     
  TProof *proof = TProof::Open("workers=8");
  //proof->SetProgressDialog(0);  
  ch.SetProof();
  ch.Process("/home/cdaq/hallc-online/hallc_replay_lt/UTIL_PION/scripts_Yield/PionYield.C+","1");
  proof->Close();
  
}
