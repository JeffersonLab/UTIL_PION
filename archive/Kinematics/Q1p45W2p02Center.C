#include <TProof.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

void Q1p45W2p02Center()
{
  TChain ch("T");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8583_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8584_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8661_-1.root"); 
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8662_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8663_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8664_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8665_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8666_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8667_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8668_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8669_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8670_-1.root");
 
  TProof *proof = TProof::Open("workers=8");
  //proof->SetProgressDialog(0);  
  ch.SetProof();
  ch.Process("/home/cdaq/hallc-online/hallc_replay_lt/UTIL_PION/scripts_Yield/PionYield.C+","1");
  proof->Close();
  
}
