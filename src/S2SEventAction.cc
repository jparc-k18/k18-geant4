// -*- C++ -*-

#include "S2SEventAction.hh"

#include <iomanip>
#include <time.h>

#include <G4Event.hh>
#include <G4SDManager.hh>
#include <G4Trajectory.hh>
#include <G4TrajectoryContainer.hh>
#include <G4UnitsTable.hh>
#include <G4VTrajectory.hh>
#include <G4VVisManager.hh>
#include <Randomize.hh>

#include "FuncName.hh"
#include "S2SAnaManager.hh"

namespace
{
auto& anaMan = S2SAnaManager::GetInstance();
}

//_____________________________________________________________________________
S2SEventAction::S2SEventAction()
{
}

//_____________________________________________________________________________
S2SEventAction::~S2SEventAction()
{
}

//_____________________________________________________________________________
void
S2SEventAction::BeginOfEventAction(const G4Event* evt)
{
  anaMan.BeginOfEvent(evt);
}

//_____________________________________________________________________________
void
S2SEventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id = evt->GetEventID();
  if(event_id % 1000 == 0){
    G4cout << FUNC_NAME << " " << event_id << G4endl;
  }
  anaMan.EndOfEvent( evt );

  // auto trajectoryContainer = evt->GetTrajectoryContainer();
  // if(trajectoryContainer && G4VVisManager::GetConcreteInstance()){
  //   G4int n_trajectories = trajectoryContainer->entries();
  //   for(G4int i=0; i<n_trajectories; ++i){
  //     auto trj = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
  //     trj->DrawTrajectory();
  //   }
  // }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
