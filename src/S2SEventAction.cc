/*
  "S2SEventAction.cc"

  Modified by Toshiyuki Gogami on 20Nov2014.
*/


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "S2SEventAction.hh"
//#include "s2sCounterHit.hh"
#include "S2SAnaManager.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SDManager.hh"
#include "DCHit.hh"
#include "Randomize.hh"
#include <iomanip>
#include <time.h>

namespace
{
using namespace CLHEP;
auto& anaMan = S2SAnaManager::GetInstance();
}

S2SEventAction::S2SEventAction()
  : drawFlag("all"), printModulo(100), DCcolID(-1)
{
  start = time(NULL);
  time(&start);
}

S2SEventAction::~S2SEventAction()
{
}

void S2SEventAction::BeginOfEventAction(const G4Event* evt)
{
  anaMan.BeginOfEvent( evt );
}

void S2SEventAction::EndOfEventAction(const G4Event* evt)
{
  // ~~~~~~~~ Periodic printing ~~~~~~~~~
  G4int event_id = evt->GetEventID();
  char anatime[100];
  if((int)event_id % 10000 == 0 ){
    end = time(NULL);
    time(&end);
    sprintf( anatime,"%.0f Sec",difftime(end,start) );
    G4cout << "Event ID : " << evt->GetEventID()
	   << " ( " << anatime << " ) " << G4endl;
  }

  anaMan.EndOfEvent( evt );

  // extract the trajectories and draw them

  // You can get a default drawing without this code by using, e.g.,
  // /vis/scene/add/trajectories 1000
  // The code here adds sophistication under control of drawFlag.

  // See comments in G4VTrajectory::DrawTrajectory for the
  // interpretation of the argument, 1000.

  //
  //  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  //  if (pVisManager)
  //    {
  //      G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  //      G4int n_trajectories = 0;
  //      if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  //
  //      for (G4int i=0; i<n_trajectories; i++)
  //        { G4VTrajectory* trj = ((*(evt->GetTrajectoryContainer()))[i]);
  //	  //trj->DrawTrajectory(1);
  //	  /*
  //	    if (drawFlag == "all") pVisManager->Draw(*trj,1000);
  //	    else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
  //	    pVisManager->Draw(*trj,1000);
  //	    else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.))
  //	    pVisManager->Draw(*trj,1000);
  //	  */
  //        }
  //    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
