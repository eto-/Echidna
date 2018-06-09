{
  std::cout << ">>> Welcome in the RootEchidna world" << std::endl;
  gSystem->Load("libHist");
  gSystem->Load("libPhysics");
  gROOT->LoadMacro("rootechidna.so");
  gROOT->LoadMacro("SetAliases.C");
//  gSystem->SetAclicMode(TSystem::kDebug);
  #include <stdexcept>
}

