/*
  "S2S.cc"
  
  A Monte Carlo simulation for S-2S spectrometer.
  Original code was made by S.Kanatsuki.
  Modified by Toshiyuki Gogami , 10Nov2014 
*/

#ifdef G4UI_USE_TCSH
#include "G4UItcsh.hh"
#endif

#include "G4RunManager.hh"
#include "G4UIterminal.hh"
#include "S2SDetectorConstruction.hh"
//#include "DetectorConstructionStsGas.hh"
#include "GetNumberFromKernelEntropyPool.hh"
#include "S2SRunAction.hh"
#include "S2SPhysicsList.hh"
#include "S2SSteppingAction.hh"
#include "S2SPrimaryGeneratorAction.hh"
#include "S2SEventAction.hh"
#include "S2SAnalysis.hh"
#include "ConfMan.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

using namespace CLHEP;

// ##################################
// ##### Main #######################
int main(int argc, char** argv) {//##
// ##################################

//  ====== Choose the Random Engine ========= //
//HepRandom::setTheEngine(new RanecuEngine);
//for(;;) G4cout<<mrad<<G4endl;
  //CLHEP::HepRandom::setTheSeed(time(NULL));
  int initSeed=GetIntFromKernelEntropyPool()&0x7FFFFFFF;
  CLHEP::HepRandom::setTheSeed(initSeed);
  int startSeed=CLHEP::HepRandom::getTheSeed();
  G4cout << "*** Initial Seed = " << startSeed << G4endl;
  CLHEP::HepRandom::showEngineStatus();
  //HepRandom::setTheEngine(new RanecuEngine);
  

  // run manager
  G4RunManager* runManager= new G4RunManager;  
  G4cout << G4endl;

  // ==== Initialize parameters ====================
  //G4double momentum = 1300.0;
  G4double distance = 600.0;
  //G4String confFile = "Conf/sts.conf.default";
  G4String confFile = "./conf/param.conf";
  G4String histFile = "../analysis/root/test.root";
  G4String dataFile = "../analysis/dat/test.dat";
  G4String macFile = "./macro/run.mac";
  if (argc==2){
    macFile = argv[1];
  }
  else if(argc == 3){
    macFile  = argv[1];
    confFile = argv[2];
  }
  else if(argc == 5){
    macFile  = argv[1];
    confFile = argv[2];
    histFile = argv[3];
    dataFile = argv[4];
  }
  else{
    G4cout << G4endl;
    G4cout << "Please type: S2S (macro file) (configuration file) (outputroot) (outputdat)" << G4endl;
    G4cout << G4endl;
    return 99;
  }
//  else if(argc == 8){
//    distance = atof(argv[3]);
//    momentum = atof(argv[4]);
//    confFile = argv[5];
//    histFile = argv[6];
//    dataFile = argv[7];
//  }
//  else if(argc!=1){
//    G4cerr<<"bin gun run T2Q1 mom conf hist data"<<G4endl;
//    return -1;
//  }

  // ===== Configuration Manager ===================
  ConfMan * confManager;
  if(argc==5) confManager= new ConfMan( confFile, histFile);
  else confManager= new ConfMan( confFile);
  confManager->Initialize();

  // ===== Analysis Manager =======================
  //S2SAnalysis * anaMan = new S2SAnalysis( histFile );
  S2SAnalysis * anaMan;
  if(argc==5) anaMan = new S2SAnalysis( histFile );
  else anaMan = new S2SAnalysis( confManager );
  anaMan->SetDataFile( dataFile );

  // ======== Detector Construction ==========================
  //runManager-> SetUserInitialization( new S2SDetectorConstruction() );
  runManager-> SetUserInitialization( new S2SDetectorConstruction(confManager) );
  
  // ============ Particles and Physics processes ===============
  //runManager-> SetUserInitialization(new S2SPhysicsList);
  runManager-> SetUserInitialization( new S2SPhysicsList(confManager) );
  
  // ======= Primary Generator =============
  //runManager-> SetUserAction(new S2SPrimaryGeneratorAction(anaMan,momentum,distance));
  //runManager-> SetUserAction(new S2SPrimaryGeneratorAction(anaMan,confManager,distance));
  runManager-> SetUserAction( new S2SPrimaryGeneratorAction(anaMan,confManager) );

  // ====== Set User Acttions =====================
  runManager-> SetUserAction(new S2SEventAction(anaMan));
  // runManager-> SetUserAction(new EventAction);
  runManager-> SetUserAction(new S2SSteppingAction);
  G4UserRunAction* run_action = new S2SRunAction(anaMan);
  runManager->SetUserAction(run_action);

  // ======== initialize visualization package ==========
#ifdef G4VIS_USE
  G4VisManager* visManager= new G4VisExecutive;
  visManager-> Initialize();
  G4cout << G4endl;
#endif

  // =========== User Session ========================
  runManager->Initialize();

  if(argc==1) { // interactive session, if no arguments given
#ifdef G4UI_USE_TCSH
    // ~~~~~~ tcsh-like ~~~~~~~~
    G4UItcsh* tcsh= new 
      G4UItcsh("s-2s(%s)[%/]:");
    G4UIterminal* session= new G4UIterminal(tcsh);
#else
    // ~~~~~~ csh-like ~~~~~~~~
    G4UIterminal* session= new G4UIterminal();
    session-> SetPrompt("s-2s(%s)[%/]:");
#endif
    session-> SessionStart();
    delete session;
  } 
  // ~~~~~~ batch mode ~~~~~~~
  else { 
    G4UImanager* UImanager= G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    //G4String fileName1 = argv[1];
    //G4String fileName2 = argv[2];
    //UImanager-> ApplyCommand(command+fileName1);
    //UImanager-> ApplyCommand(command+fileName2);
    UImanager-> ApplyCommand(command+macFile);
  }

  
  anaMan->Terminate();
  //  G4cout<<"!!! anaMan->Terminate() !!!"<<G4endl;
  
  // ========== terminating... =======================
#ifdef G4VIS_USE
  delete visManager;
#endif
  //  G4cout<<"!!! visManager !!!"<<G4endl;
  delete anaMan;
  //  G4cout<<"!!! anaMan delete !!!"<<G4endl;
  delete runManager;  
  confManager->ShowParam();
  confManager->OutputLog();
  
  return 0;
}
