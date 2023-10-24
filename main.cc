// -*- C++ -*-

#include <G4MTRunManager.hh>
#include <G4RunManager.hh>
#include <G4RunManagerFactory.hh>
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"
#include <G4UIExecutive.hh>
#include "G4VisExecutive.hh"
#include <QGSP_BERT.hh>

#include "ConfMan.hh"
#include "S2SAnaManager.hh"
#include "S2SActionInitialization.hh"
#include "S2SDetectorConstruction.hh"
#include "GetNumberFromKernelEntropyPool.hh"
#include "S2SRunAction.hh"
// #include "S2SPhysicsList.hh"
#include "S2SSteppingAction.hh"
#include "S2SPrimaryGeneratorAction.hh"
#include "S2SEventAction.hh"

enum EArgv { kProcess, kConfFile, kOutFile, kG4Macro, kArgc };

//_____________________________________________________________________________
int
main(int argc, char** argv)
{
  if(argc != kArgc-1 && argc != kArgc){
    G4cout << "Usage: " << argv[kProcess]
           << " [ConfFile] [OutputName] (G4Macro)" << G4endl;
    return EXIT_SUCCESS;
  }

  auto& confMan = ConfMan::GetInstance();
  if(!confMan.Initialize(argv[kConfFile])){
    return EXIT_FAILURE;
  }
  auto& anaMan = S2SAnaManager::GetInstance();
  anaMan.SetFileName(argv[kOutFile]);

  auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetUserInitialization(new S2SDetectorConstruction);
  runManager->SetUserInitialization(new QGSP_BERT);
  // runManager-> SetUserInitialization(new S2SPhysicsList);
  runManager->SetUserInitialization(new S2SActionInitialization);
  runManager->Initialize();

  auto visManager = new G4VisExecutive;
  visManager->SetVerboseLevel(0);
  visManager->Initialize();

  auto uiManager = G4UImanager::GetUIpointer();
  if(argc == kArgc-1){
    auto ui = new G4UIExecutive(argc, argv);
    uiManager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }else{
    G4String command = "/control/execute ";
    G4String macroFile(argv[kG4Macro]);
    uiManager->ApplyCommand(command+macroFile);
  }

  delete visManager;
  delete runManager;
  return EXIT_SUCCESS;
}
