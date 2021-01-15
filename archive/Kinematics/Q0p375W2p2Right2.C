#include <TProof.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

void Q0p375W2p2Right2()
//void combine_files()
{

  //Combine Runs of Kin. Setting: Q2 = 0.375,  W = 2.2,  Middle Epsilon
  TChain ch("T");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8585_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8586_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8587_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8603_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8604_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8605_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8606_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8607_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8608_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8609_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8610_-1.root");
  ch.Add("/net/cdaq/cdaql1data/cdaq/hallc-online/ROOTfiles/PionLT_coin_replay_production_8611_-1.root");
   
  TProof *proof = TProof::Open("workers=8");
  //proof->SetProgressDialog(0);  
  ch.SetProof();
  ch.Process("/home/cdaq/hallc-online/hallc_replay_lt/UTIL_PION/scripts_Yield/PionYield.C+","1");
  proof->Close();
  
}
