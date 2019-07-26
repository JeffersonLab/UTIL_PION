#include <TProof.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

void Q1p45W2p02Left()
{
  TChain ch("T");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8671_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8672_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8673_-1.root"); 
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8674_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8675_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8676_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8677_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8678_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8679_-1.root");
 

  TProof *proof = TProof::Open("workers=8");
  //proof->SetProgressDialog(0);  
  ch.SetProof();
  ch.Process("/home/cdaq/hallc-online/hallc_replay_lt/UTIL_PION/scripts_Yield/PionYield.C+","1");
  proof->Close();
  
}
