
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: ExN03DetectorConstruction.cc,v 1.23 2006/06/29 17:48:54 gunter Exp $
// GEANT4 tag $Name: geant4-08-03-patch-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstructionStsGas.hh"
#include "MagnetConstant.hh"
#include "StsField.hh"
#include "MaterialList.hh"
#include "RadDeg.hh"
#include "DCSD.hh"
#include "SlSD.hh"
#include "TOFSD.hh"
#include "ACSD.hh"
#include "WCSD.hh"
//#include "Area.hh"

#include "ConfMan.hh"
#include "DCGeomMan.hh"

#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4SimpleRunge.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "G4SDManager.hh"
//#include "CalorimeterSD.hh"
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstructionStsGas::DetectorConstructionStsGas()
  :solidWorld(0),logicWorld(0),physiWorld(0),mList_(0)
{
}

DetectorConstructionStsGas::~DetectorConstructionStsGas()
{// delete detectorMessenger;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MaterialList* DetectorConstructionStsGas::DefineMaterial()
{
  if(mList_) delete mList_;
  return new MaterialList();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstructionStsGas::Construct()
{
  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //Define Materials
  mList_ = DefineMaterial();


  /////////     
  // World
  /////////
  solidWorld = new G4Box("World", 5.0*m, 5.0*m, 5.0*m);

  logicWorld = new G4LogicalVolume(solidWorld, mList_->Vacuum, "World");

  physiWorld = new G4PVPlacement(0, G4ThreeVector(),
				 logicWorld, "World", 0, false, 0);

  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

  //Magnetic Field
  MakeField();

  //Construct s-2s magnets
  ConstructQ1(physiWorld);
  ConstructQ2(physiWorld);
  ConstructD1(physiWorld);
  
//   //Make s-2s position detectors
  MakePositionDetector(physiWorld);

//   //Make s-2s counters
  MakeTOFCounter(physiWorld);
//   MakeAerogelCounter(physiWorld);
//   MakeWaterCounter(physiWorld);

  //Make slits
  //  MakeSlits(physiWorld);

  //always return the physical World
  return physiWorld;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void DetectorConstructionStsGas::MakeField()
{

  ConfMan *confMan = ConfMan::GetConfManager();
  G4cout<<"called Conf in MakeField"<<G4endl;

  StsField *field = new StsField(confMan->GetFieldMap());
  G4cout<<"call field in Makefield"<<G4endl;

  G4FieldManager *fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->SetDetectorField(field);
  fieldManager->CreateChordFinder(field);

//   G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs(field);
//   G4MagIntegratorStepper *pStepper = new G4SimpleRunge(fEquation);
//   G4ChordFinder *pChordFinder = new G4ChordFinder(field, 1.0E-3*mm, pStepper);
//   fieldManager->SetChordFinder(pChordFinder);

  fieldManager->GetChordFinder()->SetDeltaChord( 1.0E-3*mm);
  //fieldManager->SetDeltaIntersection( 1.0E-6*mm );
  //fieldManager->SetDeltaOneStep( 1.0E-4*mm );
  //fieldManager->SetMaximumEpsilonStep( 1.0E-3 );
}


//////////
/// Q1 ///
//////////

void DetectorConstructionStsGas::ConstructQ1(G4VPhysicalVolume *pMother)
{

  G4double x, y, z;
  G4Material *PoleMater = mList_->Fe;
  G4Material *Q1GapMater = mList_->HeGas;//Vacuum;
  
  G4double tmpB1 = 193.75*mm;
  G4double tmpPoleR = 180.84*mm;
  G4double Q1a0 = 155*mm;

  G4RotationMatrix rotQ1Box;
  rotQ1Box.rotateZ(45*deg);
  G4double Q1Box3V = 900+1000/sqrt(2);

  G4Box *solQ1Box = new G4Box("solidQ1Pole", 2400*mm/2, 2400*mm/2, 880*mm/2);
  G4Box *tmpQ1Box = new G4Box("tmp",2000*mm/2, 2000*mm/2, 2000*mm/2);
  G4SubtractionSolid *solQ1Pole1 = new G4SubtractionSolid("tmp", solQ1Box, tmpQ1Box,
							  G4Transform3D(rotQ1Box, G4ThreeVector(Q1Box3V,Q1Box3V,0)));
  G4SubtractionSolid *solQ1Pole2 = new G4SubtractionSolid("tmp", solQ1Pole1, tmpQ1Box,
							  G4Transform3D(rotQ1Box, G4ThreeVector(-Q1Box3V,Q1Box3V,0)));
  G4SubtractionSolid *solQ1Pole3 = new G4SubtractionSolid("tmp", solQ1Pole2, tmpQ1Box,
							  G4Transform3D(rotQ1Box, G4ThreeVector(Q1Box3V,-Q1Box3V,0)));
  G4SubtractionSolid *solQ1Pole4 = new G4SubtractionSolid("tmp", solQ1Pole3, tmpQ1Box,
							  G4Transform3D(rotQ1Box, G4ThreeVector(-Q1Box3V,-Q1Box3V,0)));

  G4Box *tmpBox1 = new G4Box("tmp", tmpB1, tmpB1, 882*mm/2);
  G4Tubs *tmpSolPole = new G4Tubs("tmp",
				 0*m, tmpPoleR, 884*mm/2,
				 0*degree, 360*degree );

  G4double tmpPoleCent = (tmpPoleR+Q1a0)/sqrt(2);
  G4SubtractionSolid *tmpSolGap1 = new G4SubtractionSolid("tmp", tmpBox1, tmpSolPole,
							0, G4ThreeVector(+tmpPoleCent, +tmpPoleCent, 0) );
  G4SubtractionSolid *tmpSolGap2 = new G4SubtractionSolid("tmp", tmpSolGap1, tmpSolPole,
							0, G4ThreeVector(-tmpPoleCent, +tmpPoleCent, 0) );
  G4SubtractionSolid *tmpSolGap3 = new G4SubtractionSolid("tmp", tmpSolGap2, tmpSolPole,
							0, G4ThreeVector(-tmpPoleCent, -tmpPoleCent, 0) );
  G4SubtractionSolid *tmpSolGap4 = new G4SubtractionSolid("tmp", tmpSolGap3, tmpSolPole,
							0, G4ThreeVector(+tmpPoleCent, -tmpPoleCent, 0) );

  G4Trd *tmpSolTrd = new G4Trd("tmp", 37, 62, 882/2., 882/2., (257-193.75)/2);
  G4RotationMatrix tmpGapRot1;
  tmpGapRot1.rotateX(90.*deg);

  G4UnionSolid *tmpSolGap5 = new G4UnionSolid("tmp", tmpSolGap4, tmpSolTrd, 
					      G4Transform3D(tmpGapRot1,G4ThreeVector(0,(257+193.75)/2,0)));
  tmpGapRot1.rotateZ(90*deg);
  G4UnionSolid *tmpSolGap6 = new G4UnionSolid("tmp", tmpSolGap5, tmpSolTrd, 
					      G4Transform3D(tmpGapRot1,G4ThreeVector(-(257+193.75)/2,0,0)));
  tmpGapRot1.rotateZ(90*deg);
  G4UnionSolid *tmpSolGap7 = new G4UnionSolid("tmp", tmpSolGap6, tmpSolTrd, 
					      G4Transform3D(tmpGapRot1,G4ThreeVector(0,-(257+193.75)/2,0)));
  tmpGapRot1.rotateZ(90*deg);
  G4UnionSolid *tmpSolGap8 = new G4UnionSolid("tmp", tmpSolGap7, tmpSolTrd, 
					      G4Transform3D(tmpGapRot1,G4ThreeVector((257+193.75)/2,0,0)));

  G4Box *tmpBox2 = new G4Box("tmp", 37., 11./2, 882./2);
  G4RotationMatrix tmpGapRot2;

  G4UnionSolid *tmpSolGap9 = new G4UnionSolid("tmp", tmpSolGap8, tmpBox2, 
					      G4Transform3D(tmpGapRot2,G4ThreeVector(0,5.5+257,0)));
  tmpGapRot2.rotateZ(90.*deg);
  G4UnionSolid *tmpSolGap10 = new G4UnionSolid("tmp", tmpSolGap9, tmpBox2, 
					       G4Transform3D(tmpGapRot2,G4ThreeVector(-5.5-257,0,0)));
  tmpGapRot2.rotateZ(90.*deg);
  G4UnionSolid *tmpSolGap11 = new G4UnionSolid("tmp", tmpSolGap10, tmpBox2, 
					      G4Transform3D(tmpGapRot2,G4ThreeVector(0,-5.5-257,0)));
  tmpGapRot2.rotateZ(90.*deg);
  G4UnionSolid *tmpSolGap12 = new G4UnionSolid("tmp", tmpSolGap11, tmpBox2, 
					       G4Transform3D(tmpGapRot2,G4ThreeVector(5.5+257,0,0)));

  //  G4SubtractionSolid *solQ1Pole5 = new G4SubtractionSolid("solQ1Pole", solQ1Pole4, tmpSolPolePart12);

  /////////////////////////////////////////////////////////
  ////////////

  //  G4Trd *tmpSolTrd = new G4Trd("tmp", 37, 62, 882/2., 882/2., (257-193.75)/2);
  //  G4RotationMatrix tmpGapRot1;
  tmpGapRot1.rotateX(90.*deg);

  G4SubtractionSolid *tmpSolPolePart5 = new G4SubtractionSolid("tmp", solQ1Pole4, tmpSolTrd, 
					      G4Transform3D(tmpGapRot1,G4ThreeVector(0,(257+193.75)/2,0)));
  tmpGapRot1.rotateZ(90*deg);
  G4SubtractionSolid *tmpSolPolePart6 = new G4SubtractionSolid("tmp", tmpSolPolePart5, tmpSolTrd, 
					      G4Transform3D(tmpGapRot1,G4ThreeVector(-(257+193.75)/2,0,0)));
  tmpGapRot1.rotateZ(90*deg);
  G4SubtractionSolid *tmpSolPolePart7 = new G4SubtractionSolid("tmp", tmpSolPolePart6, tmpSolTrd, 
					      G4Transform3D(tmpGapRot1,G4ThreeVector(0,-(257+193.75)/2,0)));
  tmpGapRot1.rotateZ(90*deg);
  G4SubtractionSolid *tmpSolPolePart8 = new G4SubtractionSolid("tmp", tmpSolPolePart7, tmpSolTrd, 
					      G4Transform3D(tmpGapRot1,G4ThreeVector((257+193.75)/2,0,0)));


  //  G4Box *tmpBox2 = new G4Box("tmp", 37., 11./2, 882./2);
  //  G4RotationMatrix tmpGapRot2;
  tmpGapRot2.rotateZ(90.*deg);

  G4SubtractionSolid *tmpSolPolePart9 = new G4SubtractionSolid("tmp", tmpSolPolePart8, tmpBox2, 
					      G4Transform3D(tmpGapRot2,G4ThreeVector(0,5.5+257,0)));
  tmpGapRot2.rotateZ(90.*deg);
  G4SubtractionSolid *tmpSolPolePart10 = new G4SubtractionSolid("tmp", tmpSolPolePart9, tmpBox2, 
					       G4Transform3D(tmpGapRot2,G4ThreeVector(-5.5-257,0,0)));
  tmpGapRot2.rotateZ(90.*deg);
  G4SubtractionSolid *tmpSolPolePart11 = new G4SubtractionSolid("tmp", tmpSolPolePart10, tmpBox2, 
					      G4Transform3D(tmpGapRot2,G4ThreeVector(0,-5.5-257,0)));
  tmpGapRot2.rotateZ(90.*deg);
  G4SubtractionSolid *tmpSolPolePart12 = new G4SubtractionSolid("tmp", tmpSolPolePart11, tmpBox2, 
					       G4Transform3D(tmpGapRot2,G4ThreeVector(5.5+257,0,0)));

  //  G4SubtractionSolid *solQ1Pole5 = new G4SubtractionSolid("solQ1Pole", solQ1Pole4, tmpBox1);
  G4SubtractionSolid *solQ1Pole5 = new G4SubtractionSolid("solQ1Pole", tmpSolPolePart12, tmpSolGap4);
  //////////////
  /////////////////////////////////////////////////////////

  G4LogicalVolume *logQ1Pole = new G4LogicalVolume(solQ1Pole5, PoleMater, "logQ1Pole");
  G4LogicalVolume *logQ1Gap = new G4LogicalVolume(tmpSolGap12, Q1GapMater, "logQ1Gap");
  //  G4LogicalVolume *logQ1Gap = new G4LogicalVolume(tmpBox1, Q1GapMater, "logQ1Gap");

  G4RotationMatrix rotQ1Pole;
  rotQ1Pole.rotateZ(90*deg);
  rotQ1Pole.rotateY(90*deg);
  x = -(rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2 + Q2z + driftL1 + Q1z/2.);
  y = 0.*m;
  z = 0.*m;
  G4ThreeVector gloPosQ1(x, y, z);

  //  G4VPhysicalVolume* physQ1Pole = 
  new G4PVPlacement( G4Transform3D(rotQ1Pole,gloPosQ1),
		     "physQ1Pole", logQ1Pole, pMother, false, 0 );
  //  G4VPhysicalVolume* physQ1Gap = 
  new G4PVPlacement( G4Transform3D(rotQ1Pole,gloPosQ1),
		     "physQ1Gap", logQ1Gap, pMother, false, 0 );
  
  logQ1Pole->SetVisAttributes(G4VisAttributes::Invisible);//G4VisAttributes(true,G4Colour(0.0, 1.0, 0.0)));
  logQ1Gap->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 1.0, 0.0)));

}

//////////
/// Q2 ///
//////////

void DetectorConstructionStsGas::ConstructQ2(G4VPhysicalVolume *pMother)
{

  G4double x, y, z;
  G4Material *PoleMater = mList_->Fe;
  G4Material *Q2GapMater = mList_->HeGas;

  G4double tmpB2=243.61*mm;
  G4double tmpPoleR2 = 219.28*mm;
  G4double Q2a0 = 180*mm;

  G4Box *solQ2Box = new G4Box("solidQ2Pole", 1050*mm, 770*mm, 540*mm/2);

  G4Box *tmpBox1 = new G4Box("tmp", tmpB2, tmpB2, 542*mm/2);
  G4Tubs *tmpSolPole = new G4Tubs("tmp",
				 0*m, tmpPoleR2, 544*mm/2,
				 0*degree, 360*degree );

  G4double tmpPoleCent = (tmpPoleR2+Q2a0)/sqrt(2);
  G4SubtractionSolid *tmpSolGap1 = new G4SubtractionSolid("tmp", tmpBox1, tmpSolPole,
							0, G4ThreeVector(+tmpPoleCent, +tmpPoleCent, 0) );
  G4SubtractionSolid *tmpSolGap2 = new G4SubtractionSolid("tmp", tmpSolGap1, tmpSolPole,
							0, G4ThreeVector(-tmpPoleCent, +tmpPoleCent, 0) );
  G4SubtractionSolid *tmpSolGap3 = new G4SubtractionSolid("tmp", tmpSolGap2, tmpSolPole,
							0, G4ThreeVector(-tmpPoleCent, -tmpPoleCent, 0) );
  G4SubtractionSolid *tmpSolGap4 = new G4SubtractionSolid("tmp", tmpSolGap3, tmpSolPole,
							0, G4ThreeVector(+tmpPoleCent, -tmpPoleCent, 0) );

  G4Trd *tmpSolTrdY = new G4Trd("tmp", 53, 66.5, 542/2., 542/2., (249.85-243.61)/2);
  G4RotationMatrix tmpGapRot1Y;
  tmpGapRot1Y.rotateX(90.*deg);

  G4UnionSolid *tmpSolGap5 = new G4UnionSolid("tmp", tmpSolGap4, tmpSolTrdY, 
					      G4Transform3D(tmpGapRot1Y,G4ThreeVector(0,(249.85+243.61)/2,0)));
  tmpGapRot1Y.rotateZ(180*deg);
  G4UnionSolid *tmpSolGap6 = new G4UnionSolid("tmp", tmpSolGap5, tmpSolTrdY, 
					      G4Transform3D(tmpGapRot1Y,G4ThreeVector(0,-(249.85+243.61)/2,0)));

  G4Trd *tmpSolTrdX = new G4Trd("tmp", 42.6316, 66.5, 542/2., 542/2., (380-243.61)/2);
  G4RotationMatrix tmpGapRot1X;
  tmpGapRot1X.rotateX(90.*deg);

  tmpGapRot1X.rotateZ(-90.*deg);
  G4UnionSolid *tmpSolGap7 = new G4UnionSolid("tmp", tmpSolGap6, tmpSolTrdX, 
					      G4Transform3D(tmpGapRot1X,G4ThreeVector((380+243.61)/2,0,0)));
  tmpGapRot1X.rotateZ(180*deg);
  G4UnionSolid *tmpSolGap8 = new G4UnionSolid("tmp", tmpSolGap7, tmpSolTrdX, 
					      G4Transform3D(tmpGapRot1X,G4ThreeVector(-(380+243.61)/2,0,0)));

  G4Trd *tmpSolTrdX2 = new G4Trd("tmp", 33, 42.6316, 542/2., 542/2., (465.85-380)/2);
  G4RotationMatrix tmpGapRot1X2;
  tmpGapRot1X2.rotateX(90.*deg);

  tmpGapRot1X2.rotateZ(-90.*deg);
  G4UnionSolid *tmpSolGap9 = new G4UnionSolid("tmp", tmpSolGap8, tmpSolTrdX2, 
					      G4Transform3D(tmpGapRot1X2,G4ThreeVector((465.85+380)/2,0,0)));
  tmpGapRot1X2.rotateZ(180*deg);
  G4UnionSolid *tmpSolGap10 = new G4UnionSolid("tmp", tmpSolGap9, tmpSolTrdX2, 
					      G4Transform3D(tmpGapRot1X2,G4ThreeVector(-(465.85+380)/2,0,0)));

  ///Gap Box///
  G4Box *tmpBox2Y = new G4Box("tmp", 53., (303-249.85)/2, 542./2);
  G4RotationMatrix tmpGapRot2Y;

  G4UnionSolid *tmpSolGap11 = new G4UnionSolid("tmp", tmpSolGap10, tmpBox2Y, 
					      G4Transform3D(tmpGapRot2Y,G4ThreeVector(0,(303+249.85)/2,0)));
  tmpGapRot2Y.rotateZ(180.*deg);
  G4UnionSolid *tmpSolGap12 = new G4UnionSolid("tmp", tmpSolGap11, tmpBox2Y, 
					       G4Transform3D(tmpGapRot2Y,G4ThreeVector(0,-(303+249.85)/2,0)));

  G4Box *tmpBox2X = new G4Box("tmp", (593-465.85)/2, 33., 542./2);
  G4RotationMatrix tmpGapRot2X;

  G4UnionSolid *tmpSolGap13 = new G4UnionSolid("tmp", tmpSolGap12, tmpBox2X, 
					       G4Transform3D(tmpGapRot2X,G4ThreeVector((593+465.85)/2,0,0)));
  G4UnionSolid *tmpSolGap14 = new G4UnionSolid("tmp", tmpSolGap13, tmpBox2X, 
					       G4Transform3D(tmpGapRot2X,G4ThreeVector(-(593+465.85)/2,0,0)));

  G4SubtractionSolid *solQ2Pole5 = new G4SubtractionSolid("solQ2Pole", solQ2Box, tmpSolGap14);

  G4LogicalVolume *logQ2Pole = new G4LogicalVolume(solQ2Pole5, PoleMater, "logQ2Pole");
  G4LogicalVolume *logQ2Gap = new G4LogicalVolume(tmpSolGap14, Q2GapMater, "logQ2Gap");

  G4RotationMatrix rotQ2Pole;
  rotQ2Pole.rotateZ(90*deg);
  rotQ2Pole.rotateY(90*deg);
  x = -(rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2 + Q2z/2);
  y = 0.*m;
  z = 0.*m;
  G4ThreeVector gloPosQ2(x, y, z);

  //  G4VPhysicalVolume* physQ2Pole = 
    new G4PVPlacement( G4Transform3D(rotQ2Pole,gloPosQ2),
		       "physQ2Pole", logQ2Pole, pMother, false, 0 );
  //  G4VPhysicalVolume* physQ2Gap = 
    new G4PVPlacement( G4Transform3D(rotQ2Pole,gloPosQ2),
 		       "physQ2Gap", logQ2Gap, pMother, false, 0 );
    
  logQ2Pole->SetVisAttributes(G4VisAttributes(true,G4Colour(0.0, 1.0, 0.0)));
  logQ2Gap->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 1.0, 0.0)));

}


void DetectorConstructionStsGas::ConstructD1(G4VPhysicalVolume *pMother)
{

  G4double x, y, z;
  G4Material *PoleMater = mList_->Fe;
  G4Material *D1GapMater = mList_->HeGas;

  // D magnet surface
  G4Tubs *solD1Tub = new G4Tubs( "solD1Tub", 
				 0,  Dfr2, DfHalfGap,
				 0.*degree, bendAngleDf*degree );
  // D magnet
  G4Tubs *solD1Gap = new G4Tubs( "solD1Gap", 
				 Dr1,  Dr2,  DHalfGap,
				 0.*degree, bendAngleD*degree );

  G4Box *tmpD1 = new G4Box("tmp", 1000, 1000, 1000);
  G4RotationMatrix tmpD1rot;
  tmpD1rot.rotateZ(35*deg);
  G4double tmpD1r = 1380-1000;
  G4ThreeVector tmpD13Y(tmpD1r*cos(35*deg),tmpD1r*sin(35*deg),0);

  // D magnet yoke
  G4SubtractionSolid *tmpsolD1Pole = new G4SubtractionSolid( "tmp", solD1Tub, solD1Gap );
  G4SubtractionSolid *solD1Pole = new G4SubtractionSolid( "solidDy", tmpsolD1Pole, tmpD1, G4Transform3D(tmpD1rot, tmpD13Y) );

  G4LogicalVolume *logD1Gap = new G4LogicalVolume(solD1Gap, D1GapMater, "logD1Gap");
  G4LogicalVolume *logD1Pole = new G4LogicalVolume(solD1Pole, PoleMater, "logD1Pole");

  G4RotationMatrix rotD1;
  rotD1.rotateZ(-90.*deg);

  x = -rhoD/cos(bendAngleD/2.*degree) * sin((bendAngleD/2.)*degree);
  y = rhoD/cos(bendAngleD/2.*degree) * cos((bendAngleD/2.)*degree);
  z = 0.*m;
  G4ThreeVector gloPosD1(x,y,z);

  //  G4VPhysicalVolume* physD1Pole = 
  new G4PVPlacement( G4Transform3D(rotD1, gloPosD1),
		     "physD1Pole", logD1Pole, pMother, false, 0 );
//  G4VPhysicalVolume* physD1Gap = 
  new G4PVPlacement( G4Transform3D(rotD1, gloPosD1),
		     "physD1Gap", logD1Gap, pMother/*physD1Pole*/, false, 0 );

  // D magnet endguard
  G4Box *solEnd1 = new G4Box("solEnd1",1880*mm/2.,1600*mm/2.,76*mm/2.);
  G4Box *solEnd2 = new G4Box("solEnd2",810*mm/2.,330*mm/2.,80*mm/2.);
  G4SubtractionSolid *solEnd = new G4SubtractionSolid("solEnd", solEnd1, solEnd2);

  G4LogicalVolume *logEnd = new G4LogicalVolume(solEnd, mList_->Fe, "logEnd");

  x = rhoD*tan(bendAngleD/2.*Deg2Rad) + 240*mm + 76*mm/2.;
  y = 0.*m;
  z = 0.*m;
  G4ThreeVector gloPosEnd(x, y, z);
  gloPosEnd.rotateZ(70.*deg);

  G4RotationMatrix rotEnd;
  rotEnd.rotateX(90.*deg);
  rotEnd.rotateZ(-20.*deg);

  //  G4VPhysicalVolume *physEnd = 
    new G4PVPlacement(G4Transform3D(rotEnd, gloPosEnd),
						 "physEnd", logEnd, pMother, false, 0);

  logD1Pole->SetVisAttributes(G4VisAttributes(true,G4Colour(0.0, 0.0, 1.0)));
  logD1Gap->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 1.0, 0.0)));
  logEnd->SetVisAttributes(G4VisAttributes(true,G4Colour(0.0, 0.0, 1.0)));
}



void DetectorConstructionStsGas::MakePositionDetector(G4VPhysicalVolume *pMother)
{
  //  bool flagInvisible = false;
  char name[100];

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  // Material
  G4Material *DCuAreaMater = mList_->Ar80IsoButane20Gas;
  G4Material *DCAreaMater = mList_->Ar50Ethane50Gas;
  G4Material *DCLayerMater = mList_->W;//Ar80IsoButane20Gas;
  G4Material *DCuboxMater  = mList_->Mylar;
  G4Material *DCboxMater  = mList_->Mylar;

  G4double DCuAreaX = 1000.0;//mm
  G4double DCuAreaY = 1000.0;//mm
  G4double DCuAreaZ = 80.0;//mm
  G4double DCuBoxX = 1100.0;//mm
  G4double DCuBoxY = 1100.0;//mm
  G4double DCuBoxZ = 80.2;//mm

  ///////////////////// BC1(x,x,v,v,u,u)
  G4Box *solidBc1Box =
    new G4Box( "BC1Box", DCuBoxZ/2.0*mm, DCuBoxX/2.0*mm, DCuBoxY/2.0*mm );
  G4Box *solidBc1Area =
    new G4Box( "BC1Area", DCuAreaZ/2.0*mm, DCuAreaX/2.0*mm, DCuAreaY/2.0*mm );
  G4Box *solidBc1Layer =
    new G4Box( "BC1Layer", 0.0001/2.0*mm, DCuAreaX/2.0*mm, DCuAreaY/2.0*mm );

  G4LogicalVolume *logBc1Box =
    new G4LogicalVolume( solidBc1Box,   DCuboxMater,  "BC1Box",   0, 0, 0 );
  G4LogicalVolume *logBc1Area =
    new G4LogicalVolume( solidBc1Area,  DCuAreaMater,  "BC1Area",  0, 0, 0 );
  G4LogicalVolume *logBc1Layer =
    new G4LogicalVolume( solidBc1Layer, DCLayerMater,  "BC1Layer", 0, 0, 0 );

  G4int lnum_bc1;
  lnum_bc1 = geomMan.GetDetectorId("BC1-u-1");
  G4ThreeVector posBc1G1 = geomMan.GetGlobalPosition( lnum_bc1 );
  lnum_bc1 = geomMan.GetDetectorId("BC1-u-2");
  G4ThreeVector posBc1G2 = geomMan.GetGlobalPosition( lnum_bc1 );
  G4ThreeVector posBc1G = (posBc1G1+posBc1G2)/2.0;

  G4double RotateAngleBc11;
  G4double RotateAngleBc12;
  RotateAngleBc12 = geomMan.GetRotAngle2( lnum_bc1 );
  RotateAngleBc11 = geomMan.GetRotAngle1( lnum_bc1 );
  G4RotationMatrix RMBc1;
  //  RMBc1.rotateX( RotateAngleBc11*degree );
  RMBc1.rotateZ( (RotateAngleBc12-90)*degree );

  //  G4RotationMatrix rotBc1;

  G4VPhysicalVolume *physBc1Box =
    new G4PVPlacement( G4Transform3D(RMBc1, posBc1G), 
		       "BC1Box", logBc1Box, pMother, false, 0 );
  G4VPhysicalVolume *physBc1Area =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       "BC1Area", logBc1Area, physBc1Box, false, 0 );

  char Bc1Name[6][100] = {"BC1-v-1", "BC1-v-2", "BC1-u-1", "BC1-u-2", "BC1-x-1", "BC1-x-2"};

  G4VPhysicalVolume *physBc1Layer[6];
  for (int i=0; i<6; i++) {

    lnum_bc1 = geomMan.GetDetectorId(Bc1Name[i]);
    G4ThreeVector posBc1L1 = geomMan.GetGlobalPosition( lnum_bc1 );
    G4double posBc1L1D = (posBc1L1-posBc1G).mag();

    //   G4double posBc1L1D = 10. + (i-3).*20;//mm
    if (i<3)
      posBc1L1D *= -1.; 

    physBc1Layer[i] =
      new G4PVPlacement( 0, G4ThreeVector( posBc1L1D*mm, 0.0*mm, 0.0*mm ),
			 name, logBc1Layer, physBc1Area, false,  101+i );
  }


  ///////////////////// BC2(x,x,v,v,u,u)
  G4Box *solidBc2Box =
    new G4Box( "BC2Box", DCuBoxZ/2.0*mm, DCuBoxX/2.0*mm, DCuBoxY/2.0*mm );
  G4Box *solidBc2Area =
    new G4Box( "BC2Area", DCuAreaZ/2.0*mm, DCuAreaX/2.0*mm, DCuAreaY/2.0*mm );
  G4Box *solidBc2Layer =
    new G4Box( "BC2Layer", 0.0001/2.0*mm, DCuAreaX/2.0*mm, DCuAreaY/2.0*mm );

  G4LogicalVolume *logBc2Box =
    new G4LogicalVolume( solidBc2Box,   DCuboxMater,  "BC2Box",   0, 0, 0 );
  G4LogicalVolume *logBc2Area =
    new G4LogicalVolume( solidBc2Area,  DCuAreaMater,  "BC2Area",  0, 0, 0 );
  G4LogicalVolume *logBc2Layer =
    new G4LogicalVolume( solidBc2Layer, DCLayerMater,  "BC2Layer", 0, 0, 0 );

  G4int lnum_bc2;
  lnum_bc2 = geomMan.GetDetectorId("BC2-u-1");
  G4ThreeVector posBc2G1 = geomMan.GetGlobalPosition( lnum_bc2 );
  lnum_bc2 = geomMan.GetDetectorId("BC2-u-2");
  G4ThreeVector posBc2G2 = geomMan.GetGlobalPosition( lnum_bc2 );
  G4ThreeVector posBc2G = (posBc2G1+posBc2G2)/2.0;

  G4double RotateAngleBc21;
  G4double RotateAngleBc22;
  RotateAngleBc21 = geomMan.GetRotAngle1( lnum_bc2 );
  RotateAngleBc22 = geomMan.GetRotAngle2( lnum_bc2 );
  G4RotationMatrix RMBc2;
  //  RMBc2.rotateX( RotateAngleBc21*degree );
  RMBc2.rotateZ( (RotateAngleBc22-90)*degree );

  G4VPhysicalVolume *physBc2Box =
    new G4PVPlacement( G4Transform3D(RMBc2, posBc2G), 
		       "BC2Box", logBc2Box, pMother, false, 0 );
  G4VPhysicalVolume *physBc2Area =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       "BC2Area", logBc2Area, physBc2Box, false, 0 );

  char Bc2Name[6][100] = {"BC2-v-1", "BC2-v-2", "BC2-u-1", "BC2-u-2", "BC2-x-1", "BC2-x-2"};

  G4VPhysicalVolume *physBc2Layer[6];
  for (int i=0; i<6; i++) {

    lnum_bc2 = geomMan.GetDetectorId(Bc2Name[i]);
    G4ThreeVector posBc2L1 = geomMan.GetGlobalPosition( lnum_bc2 );
    G4double posBc2L1D = (posBc2L1-posBc2G).mag();

    if (i<3)
      posBc2L1D *= -1.; 

    physBc2Layer[i] =
      new G4PVPlacement( 0, G4ThreeVector( posBc2L1D*mm, 0.0*mm, 0.0*mm ),
			 name, logBc2Layer, physBc2Area, false,  107+i );
  }

  logBc1Box->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc1Area->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc1Layer->SetVisAttributes(G4VisAttributes::Invisible);
  logBc2Box->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc2Area->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc2Layer->SetVisAttributes(G4VisAttributes::Invisible);

  G4double DCAreaX = 2000;//mm
  G4double DCAreaY = 1000;//mm
  G4double DCAreaZ = 100;//mm
  G4double DCBoxX = 2100;//mm
  G4double DCBoxY = 1100;//mm
  G4double DCBoxZ = 100.2;//mm

  ///////////////////// BC3(x,x,y,y)
  G4Box *solidBc3Box =
    new G4Box( "BC3Box", DCBoxZ/2.0*mm, DCBoxX/2.0*mm, DCBoxY/2.0*mm );
  G4Box *solidBc3Area =
    new G4Box( "BC3Area", DCAreaZ/2.0*mm, DCAreaX/2.0*mm, DCAreaY/2.0*mm );
  G4Box *solidBc3Layer =
    new G4Box( "BC3Layer", 0.0001/2.0*mm, DCAreaX/2.0*mm,  DCAreaY/2.0*mm );

  G4LogicalVolume *logBc3Box =
    new G4LogicalVolume( solidBc3Box,   DCboxMater,  "BC3Box",   0, 0, 0 );
  G4LogicalVolume *logBc3Area =
    new G4LogicalVolume( solidBc3Area,  DCAreaMater,  "BC3Area",  0, 0, 0 );
  G4LogicalVolume *logBc3Layer =
    new G4LogicalVolume( solidBc3Layer, DCLayerMater,  "BC3Layer", 0, 0, 0 );

  G4int lnum_bc3;
  lnum_bc3 = geomMan.GetDetectorId("BC3-x-2");
  G4ThreeVector posBc3G1 = geomMan.GetGlobalPosition( lnum_bc3 );
//   posBc3G1 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//   posBc3G1.rotateZ(bendAngleD*degree);
  lnum_bc3 = geomMan.GetDetectorId("BC3-y-1");
  G4ThreeVector posBc3G2 = geomMan.GetGlobalPosition( lnum_bc3 );
//   posBc3G2 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//   posBc3G2.rotateZ(bendAngleD*degree);
  G4ThreeVector posBc3G = (posBc3G1+posBc3G2)/2.0;

  G4double RotateAngleBc31;
  G4double RotateAngleBc32;
  RotateAngleBc31 = geomMan.GetRotAngle1( lnum_bc3 );
  RotateAngleBc32 = geomMan.GetRotAngle2( lnum_bc3 );
  G4RotationMatrix RMBc3;
  //  RMBc3.rotateX( RotateAngleBc31*degree );
  RMBc3.rotateZ( (RotateAngleBc32-90)*degree );

  G4VPhysicalVolume *physBc3Box =
    new G4PVPlacement( G4Transform3D(RMBc3, posBc3G), 
		       "BC3Box", logBc3Box, pMother, false, 0 );
  G4VPhysicalVolume *physBc3Area =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       "BC3Area", logBc3Area, physBc3Box, false, 0 );

  char Bc3Name[4][100] = {"BC3-x-1", "BC3-x-2", "BC3-y-1", "BC3-y-2"};

  G4VPhysicalVolume *physBc3Layer[4];
  for (int i=0; i<4; i++) {

    lnum_bc3 = geomMan.GetDetectorId(Bc3Name[i]);
    G4ThreeVector posBc3L1 = geomMan.GetGlobalPosition( lnum_bc3 );
//     posBc3L1 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//     posBc3L1.rotateZ(bendAngleD*degree);
    G4double posBc3L1D = (posBc3L1-posBc3G).mag();

    //G4double posBc3L1D = 10. + (i-3)*20.;

    if (i<2)
      posBc3L1D *= -1.; 
 
    //    sprintf(name, "BC3Layer%d", i+1);
    physBc3Layer[i] =
      new G4PVPlacement( 0, G4ThreeVector( posBc3L1D*mm, 0.0*mm, 0.0*mm ),
			 Bc3Name[i], logBc3Layer, physBc3Area, false,  113+i );
  }

  logBc3Box->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc3Area->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc3Layer->SetVisAttributes(G4VisAttributes::Invisible);


  ///////////////////// BC4(x,x,x,x,x,x)
  G4Box *solidBc4Box =
    new G4Box( "BC4Box", DCBoxZ/2.0*mm, DCBoxX/2.0*mm, DCBoxY/2.0*mm );
  G4Box *solidBc4Area =
    new G4Box( "BC4Area", DCAreaZ/2.0*mm, DCAreaX/2.0*mm, DCAreaY/2.0*mm );
  G4Box *solidBc4Layer =
    new G4Box( "BC4Layer", 0.0001/2.0*mm, DCAreaX/2.0*mm,  DCAreaY/2.0*mm );

  G4LogicalVolume *logBc4Box =
    new G4LogicalVolume( solidBc4Box,   DCboxMater,  "BC4Box",   0, 0, 0 );
  G4LogicalVolume *logBc4Area =
    new G4LogicalVolume( solidBc4Area,  DCAreaMater,  "BC4Area",  0, 0, 0 );
  G4LogicalVolume *logBc4Layer =
    new G4LogicalVolume( solidBc4Layer, DCLayerMater,  "BC4Layer", 0, 0, 0 );

  G4int lnum_bc4;
  lnum_bc4 = geomMan.GetDetectorId("BC4-x-3");
  G4ThreeVector posBc4G1 = geomMan.GetGlobalPosition( lnum_bc4 );
//   posBc4G1 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//   posBc4G1.rotateZ(bendAngleD*degree);
  lnum_bc4 = geomMan.GetDetectorId("BC4-x-4");
  G4ThreeVector posBc4G2 = geomMan.GetGlobalPosition( lnum_bc4 );
//   posBc4G2 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//   posBc4G2.rotateZ(bendAngleD*degree);
  G4ThreeVector posBc4G = (posBc4G1+posBc4G2)/2.0;

  G4double RotateAngleBc41;
  G4double RotateAngleBc42;
  RotateAngleBc41 = geomMan.GetRotAngle1( lnum_bc4 );
  RotateAngleBc42 = geomMan.GetRotAngle2( lnum_bc4 );
  G4RotationMatrix RMBc4;
  //  RMBc4.rotateX( RotateAngleBc41*degree );
  RMBc4.rotateZ( (RotateAngleBc42-90)*degree );

  //  G4ThreeVector posBc4(3.6*m,0*m,0*m);
  //  posBc4G.rotateZ(bendAngleD*degree);

  //  G4RotationMatrix RMBc4;
  //  RMBc4.rotateZ(bendAngleD*degree);

  G4VPhysicalVolume *physBc4Box =
    new G4PVPlacement( G4Transform3D(RMBc4, posBc4G), 
		       "BC4Box", logBc4Box, pMother, false, 0 );
  G4VPhysicalVolume *physBc4Area =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       "BC4Area", logBc4Area, physBc4Box, false, 0 );
  G4VPhysicalVolume *physBc4Layer[6];

  char Bc4Name[6][100] = {"BC4-x-1", "BC4-x-2", "BC4-x-3", "BC4-x-4", "BC4-x-5", "BC4-x-6"};

  for (int i=0; i<6; i++) {

    lnum_bc4 = geomMan.GetDetectorId(Bc4Name[i]);
    G4ThreeVector posBc4L1 = geomMan.GetGlobalPosition( lnum_bc4 );
//     posBc4L1 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//     posBc4L1.rotateZ(bendAngleD*degree);
    G4double posBc4L1D = (posBc4L1-posBc4G).mag();

    if (i<3)
      posBc4L1D *= -1.; 

    physBc4Layer[i] =
      new G4PVPlacement( 0, G4ThreeVector( posBc4L1D*mm, 0.0*mm, 0.0*mm ),
			 Bc4Name[i], logBc4Layer, physBc4Area, false,  117+i );
  }

  logBc4Box->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc4Area->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  //  logBc4Layer->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc4Layer->SetVisAttributes(G4VisAttributes::Invisible);

  ///////////////////// Bc5(y,y,y,y,y,y)
  G4Box *solidBc5Box =
    new G4Box( "BC5Box", DCBoxZ/2.0*mm, DCBoxX/2.0*mm, DCBoxY/2.0*mm );
  G4Box *solidBc5Area =
    new G4Box( "BC5Area", DCAreaZ/2.0*mm, DCAreaX/2.0*mm, DCAreaY/2.0*mm );
  G4Box *solidBc5Layer =
    new G4Box( "BC5Layer", 0.0001/2.0*mm, DCAreaX/2.0*mm,  DCAreaY/2.0*mm );

  G4LogicalVolume *logBc5Box =
    new G4LogicalVolume( solidBc5Box,   DCboxMater,  "BC5Box",   0, 0, 0 );
  G4LogicalVolume *logBc5Area =
    new G4LogicalVolume( solidBc5Area,  DCAreaMater,  "BC5Area",  0, 0, 0 );
  G4LogicalVolume *logBc5Layer =
    new G4LogicalVolume( solidBc5Layer, DCLayerMater,  "BC5Layer", 0, 0, 0 );

  G4int lnum_bc5;
  lnum_bc5 = geomMan.GetDetectorId("BC5-y-3");
  G4ThreeVector posBc5G1 = geomMan.GetGlobalPosition( lnum_bc5 );
//   posBc5G1 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//   posBc5G1.rotateZ(bendAngleD*degree);
  lnum_bc5 = geomMan.GetDetectorId("BC5-y-4");
  G4ThreeVector posBc5G2 = geomMan.GetGlobalPosition( lnum_bc5 );
//   posBc5G2 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//   posBc5G2.rotateZ(bendAngleD*degree);
  G4ThreeVector posBc5G = (posBc5G1+posBc5G2)/2.0;

  G4double RotateAngleBc51;
  G4double RotateAngleBc52;
  RotateAngleBc51 = geomMan.GetRotAngle1( lnum_bc5 );
  RotateAngleBc52 = geomMan.GetRotAngle2( lnum_bc5 );
  G4RotationMatrix RMBc5;
  //  RMBc5.rotateX( RotateAngleBc51*degree );
  RMBc5.rotateZ( (RotateAngleBc52-90)*degree );

  //  G4ThreeVector posBc5(3.6*m,0*m,0*m);
  //  posBc5G.rotateZ(bendAngleD*degree);

  //  G4RotationMatrix RMBc5;
  //  RMBc5.rotateZ(bendAngleD*degree);

  G4VPhysicalVolume *physBc5Box =
    new G4PVPlacement( G4Transform3D(RMBc5, posBc5G), 
		       "BC5Box", logBc5Box, pMother, false, 0 );
  G4VPhysicalVolume *physBc5Area =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       "BC5Area", logBc5Area, physBc5Box, false, 0 );
  G4VPhysicalVolume *physBc5Layer[6];

  char Bc5Name[6][100] = {"BC5-y-1", "BC5-y-2", "BC5-y-3", "BC5-y-4", "BC5-y-5", "BC5-y-6"};

  for (int i=0; i<6; i++) {

    lnum_bc5 = geomMan.GetDetectorId(Bc5Name[i]);
    G4ThreeVector posBc5L1 = geomMan.GetGlobalPosition( lnum_bc5 );
//     posBc5L1 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//     posBc5L1.rotateZ(bendAngleD*degree);
    G4double posBc5L1D = (posBc5L1-posBc5G).mag();

    if (i<3)
      posBc5L1D *= -1.; 

    physBc5Layer[i] =
      new G4PVPlacement( 0, G4ThreeVector( posBc5L1D*mm, 0.0*mm, 0.0*mm ),
			 Bc5Name[i], logBc5Layer, physBc5Area, false,  123+i );
  }

  logBc5Box->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc5Area->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  //  logBc5Layer->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc5Layer->SetVisAttributes(G4VisAttributes::Invisible);


  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  DCSD *dcSD = new DCSD("BcSD");
  SDMan->AddNewDetector(dcSD);


  ///////////////////////////////////////////////////
  //////////////Other Tracking Detector

  //HeBag between
  //target and BC1 0
  //BC1 and BC2    1
  //BC2 and Q1     2
  //Q1 and Q2      3
  //Q2 and D1      4
  //D1 and BC3     5
  //BC3 and BC4    6
  G4Box *solidHeBag[7];
  G4LogicalVolume *logHeBag[7];
  G4VPhysicalVolume *physHeBag[7];
  G4Material *HeBagMater[7]={mList_->HeGas,//0
			     mList_->HeGas,//1
			     mList_->HeGas,//2
			     mList_->HeGas,//3
			     mList_->HeGas,//4
			     mList_->HeGas,//5
			     mList_->HeGas};//6

  G4double tgt=(rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2 + Q2z + driftL1 + Q1z + driftL0);
  G4double bc1=posBc1G.mag();
  G4double bc2=posBc2G.mag();
  G4double Q1=rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2 + Q2z + driftL1 + Q1z/2.;
  G4double Q2=rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2 + Q2z/2.;
  G4double bc3=posBc3G.mag();
  G4double bc4=posBc4G.mag();
  G4double D1=rhoD*tan(bendAngleD/2.*Deg2Rad);

  G4ThreeVector HeBagG[7];
  HeBagG[0]=G4ThreeVector(-(tgt+(bc1+DCuBoxZ/2.))/2.,0,0);
  HeBagG[1]=G4ThreeVector(-((bc1-DCuBoxZ/2.)+(bc2+DCuBoxZ/2.))/2.,0,0);
  HeBagG[2]=G4ThreeVector(-((bc2-DCuBoxZ/2.)+(Q1+Q1z/2.))/2.,0,0);
  HeBagG[3]=G4ThreeVector(-((Q1-Q1z/2.)+(Q2+Q2z/2.))/2.,0,0);
  HeBagG[4]=G4ThreeVector(-((Q2-Q2z/2.)+D1)/2.,0,0);
  HeBagG[5]=G4ThreeVector((D1+(bc3-DCBoxZ/2.))/2.,0,0);
  HeBagG[6]=G4ThreeVector(((bc3+DCBoxZ/2.)+(bc4-DCBoxZ/2.))/2.,0,0);
  G4double thick[7];
  thick[0]=tgt-(bc1+DCuBoxZ/2.);
  thick[1]=(bc1-DCuBoxZ/2.)-(bc2+DCuBoxZ/2.);
  thick[2]=(bc2-DCuBoxZ/2.)-(Q1+Q1z/2.);
  thick[3]=(Q1-Q1z/2.)-(Q2+Q2z/2.);
  thick[4]=(Q2-Q2z/2.)-rhoD*tan(bendAngleD/2.*Deg2Rad);
  thick[5]=-rhoD*tan(bendAngleD/2.*Deg2Rad)+(bc3-DCBoxZ/2.);
  thick[6]=-(bc3+DCBoxZ/2.)+(bc4-DCBoxZ/2.);

  HeBagG[5].rotateZ(bendAngleD*degree);
  HeBagG[6].rotateZ(bendAngleD*degree);
  G4RotationMatrix RMHeBag[7];
  RMHeBag[5].rotateZ(bendAngleD*degree);
  RMHeBag[6].rotateZ(bendAngleD*degree);
  
  char hename[7][100]={"HeBag0","HeBag1","HeBag2","HeBag3","HeBag4","HeBag5","HeBag6"};

  for(int i=0;i<7;i++){
    if(i!=0||1){
      if(i<2){
	solidHeBag[i] =  new G4Box( hename[i], thick[i]/2.,   DCuAreaX/2.*mm, DCuAreaY/2.*mm );
	G4cout<<i<<" "<<thick[i]<<G4endl;
      }
      else{
	solidHeBag[i] =  new G4Box( hename[i], thick[i]/2.,   DCAreaX/2.*mm, DCAreaY/2.*mm );
	G4cout<<i<<" "<<thick[i]<<G4endl;
      }
      logHeBag[i] = new G4LogicalVolume( solidHeBag[i], HeBagMater[i], hename[i], 0, 0, 0);
      physHeBag[i] = new G4PVPlacement( G4Transform3D(RMHeBag[i], HeBagG[i]), 
      					hename[i], logHeBag[i], pMother, false, 0 );
    }
  }


  logBc1Layer->SetSensitiveDetector(dcSD);
  logBc2Layer->SetSensitiveDetector(dcSD);
  logBc3Layer->SetSensitiveDetector(dcSD);
  logBc4Layer->SetSensitiveDetector(dcSD);
  logBc5Layer->SetSensitiveDetector(dcSD);

}//Make PositionDetectors

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionStsGas::MakeSlits(G4VPhysicalVolume *pMother)
{

  G4Material *SlitMater = mList_->Scin;

  G4Box *solidSlit = new G4Box("solidSlit",2.0/2.*m, 2.0/2.*m, 0.000001/2.*mm);
  G4LogicalVolume *logicSlit = new G4LogicalVolume(solidSlit, SlitMater, "logicSlit");
  G4VPhysicalVolume *physSlit[8];

  G4RotationMatrix rotQ;
  rotQ.rotateZ(90.*deg);
  rotQ.rotateY(90.*deg);
  G4ThreeVector gloPosSlit[8];
  double x,y,z;
  y = 0;
  z = 0;

  x = -(rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2 + Q2z + driftL1 + Q1z/2.);
  gloPosSlit[0]=G4ThreeVector(x-Q1z/2.-5, y, z);
  gloPosSlit[1]=G4ThreeVector(x+Q1z/2.+5, y, z);

  x = -(rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2 + Q2z/2.);
  gloPosSlit[2]=G4ThreeVector(x-Q2z/2.-5, y, z);
  gloPosSlit[3]=G4ThreeVector(x+Q2z/2.+5, y, z);

  gloPosSlit[4]=G4ThreeVector(-rhoD*tan(bendAngleD/2.*Deg2Rad)-5, 0, 0);
  gloPosSlit[5]=G4ThreeVector( rhoD*tan(bendAngleD/2.*Deg2Rad)+5, 0, 0);
  gloPosSlit[6]=G4ThreeVector( rhoD*tan(bendAngleD/2.*Deg2Rad)+240+76, 0, 0);
  gloPosSlit[7]=G4ThreeVector( rhoD*tan(bendAngleD/2.*Deg2Rad)+500+1000+1000, 0, 0);

  gloPosSlit[5].rotateZ(bendAngleD*Deg2Rad);
  gloPosSlit[6].rotateZ(bendAngleD*Deg2Rad);
  gloPosSlit[7].rotateZ(bendAngleD*Deg2Rad);
  G4RotationMatrix rotSlit5;
  rotSlit5.rotateZ(90.*degree);
  rotSlit5.rotateY(90.*degree);
  rotSlit5.rotateZ((bendAngleD)*Deg2Rad);
  G4RotationMatrix rotSlit6;
  rotSlit6.rotateZ(90.*degree);
  rotSlit6.rotateY(90.*degree);
  rotSlit6.rotateZ((bendAngleD)*Deg2Rad);
  G4RotationMatrix rotSlit7;
  rotSlit7.rotateZ(90.*degree);
  rotSlit7.rotateY(90.*degree);
  rotSlit7.rotateZ((bendAngleD)*Deg2Rad);

  char name[8][100]={"slitNo.1","slitNo.2","slitNo.3","slitNo.4","slitNo.5","slitNo.6","slitNo.7","slitNo.8"};

//   physSlit[0] = 
//     new G4PVPlacement( G4Transform3D(rotQ, gloPosSlit[0]),
// 		       "physSlit[0]", logicSlit, pMother, false, 0 );
  for(int i=0;i<5;i++){
    physSlit[i] = 
      new G4PVPlacement( G4Transform3D(rotQ, gloPosSlit[i]),
			 name[i], logicSlit, pMother, false, i );
  }
  for(int i=5;i<8;i++){
    physSlit[i] = 
      new G4PVPlacement( G4Transform3D(rotSlit5, gloPosSlit[i]),
			 name[i], logicSlit, pMother, false, i );
  }
//   physSlit[6] = 
//     new G4PVPlacement( G4Transform3D(rotSlit6, gloPosSlit[6]),
// 		       "physSlit[6]", logicSlit, pMother, false, 6 );
//   physSlit[7] = 
//     new G4PVPlacement( G4Transform3D(rotSlit7, gloPosSlit[7]),
// 		       "physSlit[7]", logicSlit, pMother, false, 7 );
  
  logicSlit->SetVisAttributes(G4VisAttributes(true,G4Colour(0.0,1.0,0.0)));//::Invisible);

  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  SlSD *slSD = new SlSD("SlSD");
  SDMan->AddNewDetector(slSD);

  logicSlit->SetSensitiveDetector(slSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void DetectorConstructionStsGas::MakeTOFCounter(G4VPhysicalVolume *pMother)
{

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  //TOF
  const int SegNumTOF = 25;
  double TOFX = 70.*mm;
  double TOFY = 20.*mm;
  double TOFZ = 1000.*mm;

  // Material
  G4Material *TOFMater    = mList_->Scin;

  G4ThreeVector localPosTOF[SegNumTOF];
  G4ThreeVector globalPosTOF[SegNumTOF];

  localPosTOF[12] = G4ThreeVector(0, 0, 0);
  globalPosTOF[12] = geomMan.Local2GlobalPos(61, localPosTOF[12]);

  //  G4double TOFpos = 3800.;
//   localPosTOF[8] = G4ThreeVector(0, 0, 0);
//   globalPosTOF[8] = G4ThreeVector(TOFpos, 0, 0);

  for (int i=SegNumTOF/2-1; i>=0; i--) {
    double localx = localPosTOF[i+1].x() - TOFX;
    double localz = 0;

    localPosTOF[i] = G4ThreeVector(localx, 0, localz);
    globalPosTOF[i] = geomMan.Local2GlobalPos(61, localPosTOF[i]);
    //globalPosTOF[i] = G4ThreeVector(TOFpos+localz, localx, 0);
  }

  for (int i=SegNumTOF/2+1; i<=SegNumTOF; i++) {
    double localx = localPosTOF[i-1].x() + TOFX;
    double localz = 0;

    localPosTOF[i] = G4ThreeVector(localx, 0, localz);
    globalPosTOF[i] = geomMan.Local2GlobalPos(61, localPosTOF[i]);
    //globalPosTOF[i] = G4ThreeVector(TOFpos+localz, localx, 0);
  }

  G4double RotateAngleTOF = 70;//geomMan.GetRotAngle2( lnum_bh1 );
  G4RotationMatrix RMTOF; 
  RMTOF.rotateZ( (RotateAngleTOF)*degree );

  G4Box *solidTOF = new G4Box( "solidTOF", TOFY/2.*mm,  TOFX/2.*mm,  TOFZ/2.*mm );
  G4LogicalVolume *logTOF = new G4LogicalVolume( solidTOF, TOFMater, "logicTOF", 0, 0, 0 );
  G4VPhysicalVolume *physTOF[SegNumTOF];
  char buf[100];

  for(int i=0; i<SegNumTOF; i++) {
  
    sprintf(buf, "TOF_%d", i); 
    //    globalPosTOF[i].rotateZ(70*degree);
    physTOF[i] =
      new G4PVPlacement( G4Transform3D(RMTOF, globalPosTOF[i]),
                       buf, logTOF, pMother, false, i );
  }

  logTOF->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 1.0)));

  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  TOFSD *tofSD = new TOFSD("TOFSD");
  SDMan->AddNewDetector(tofSD);
  logTOF->SetSensitiveDetector(tofSD);

}//MakeTOFCounter

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstructionStsGas::MakeAerogelCounter(G4VPhysicalVolume* pMother){

  G4Material *ACMater = mList_->Vacuum;//Aerogel;

  G4Box *solidACBox =
    new G4Box( "ACBox",   510.0/2.0*mm, 1200.0/2.0*mm, 600.0/2.0*mm );
  G4Box *solidACArea =
    new G4Box( "ACArea",  500.0/2.0*mm, 1100.0/2.0*mm, 550.0/2.0*mm );
  G4Box *solidACRad =
    new G4Box( "ACRad",   100.0/2.0*mm, 1100.0/2.0*mm, 550.0/2.0*mm );
  G4Box *solidACLayer =
    new G4Box( "ACLayer", 0.01/2.0*mm,  1100.0/2.0*mm, 550.0/2.0*mm );
//   G4SubtractionSolid *solidACFrame =
//     new G4SubtractionSolid( "ACFrame", solidACBox, solidACArea );

  G4LogicalVolume *logACBox =
    new G4LogicalVolume( solidACBox,   mList_->Vacuum/*Air*/,  "ACBox",   0, 0, 0 );
//   G4LogicalVolume *logACFrame =
//     new G4LogicalVolume( solidACFrame, mList_->Air,  "ACFrame", 0, 0, 0 );
  G4LogicalVolume *logACArea =
    new G4LogicalVolume( solidACArea,  mList_->Vacuum/*Air*/,  "ACArea",  0, 0, 0 );
  G4LogicalVolume *logACRad =
    new G4LogicalVolume( solidACRad,   ACMater,  "ACRad", 0, 0, 0 );
  G4LogicalVolume *logACLayer =
    new G4LogicalVolume( solidACLayer, ACMater,  "ACLayer", 0, 0, 0 );

  /*
  G4int lnum_ac;
  lnum_ac = geomMan.GetDetectorId("AC-v-1");
  G4ThreeVector posACG1 = geomMan.GetGlobalPosition( lnum_ac );
  lnum_ac = geomMan.GetDetectorId("AC-v-2");
  G4ThreeVector posACG2 = geomMan.GetGlobalPosition( lnum_ac );
  G4ThreeVector posACG = (posACG1+posACG2)/2.0;

  G4double RotateAngleAC;
  RotateAngleAC = geomMan.GetRotAngle2( lnum_ac );

  G4RotationMatrix RMAC; RMAC.rotateZ( RotateAngleAC*degree );
  */

  G4ThreeVector posAC(4.2*m,0*m,0*m);
  posAC.rotateZ(bendAngleD*degree);
  G4RotationMatrix RMAC;
  RMAC.rotateZ(bendAngleD*degree);

  G4VPhysicalVolume *physACBox =
    new G4PVPlacement( G4Transform3D(RMAC,posAC), 
		       "ACBox", logACBox, pMother, false, 0 );
//   G4VPhysicalVolume *physACFrame =
//     new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
//                        "ACFrame", logACFrame, physACBox, false, 0 );
  G4VPhysicalVolume *physACArea =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       "ACArea", logACArea, physACBox, false, 0 );
  G4VPhysicalVolume *physACRad =
    new G4PVPlacement( 0, G4ThreeVector( -200.0*mm, 0.0*mm, 0.0*mm ),
                       "ACRad", logACRad, physACArea, false, 0 );
  //  G4VPhysicalVolume *physACLayer =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       "ACLayer", logACLayer, physACRad, false, 0 );


  //char ACName[6]/*[100] = {*/ = "AC";//"AC-x-1", "AC-x-2", "AC-v-1", "AC-v-2", "AC-u-1", "AC-u-2"};

  //  G4VPhysicalVolume *physACLayer[6];
  /*
  for (int i=0; i<6; i++) {
    
    lnum_ac = geomMan.GetDetectorId(ACName[i]);
    G4ThreeVector posACL1 = geomMan.GetGlobalPosition( lnum_ac );
    G4double posACL1D = (posACL1-posACG).mag();
    
    //    G4double posACL1D = 0.5 + (i-3)*1.0;
    
    if (i>=3)
      posACL1D *= -1.; 
    
    sprintf(name, "ACLayer%d", i+1);
    physACLayer[i] =
      new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, posACL1D*mm, 0.0*mm ),
			 name, logACLayer, physACArea, false,  113+i );
  }
  */

  logACBox->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 1.0)));
//  logACFrame->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logACArea->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 1.0)));
  logACRad->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 1.0)));
  logACLayer->SetVisAttributes(G4VisAttributes::Invisible);

  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  ACSD *acSD = new ACSD("AC");
  SDMan->AddNewDetector(acSD);
  logACLayer->SetSensitiveDetector(acSD);

}//MakeAerogelCounter

//MakeWaterCherenkovCounter
void DetectorConstructionStsGas::MakeWaterCounter(G4VPhysicalVolume* pMother){

  const int SegNumWC = 16;
  double WCframeX = 80.*mm;
  double WCframeY = 156.*mm;
  double WCframeZ = 506.*mm;
  double WCradX = 74.*mm;
  double WCradY = 150.*mm;
  double WCradZ = 500.*mm;

  double zdiff   = WCframeX/2;

  // Material
  G4Material *WCMater    = mList_->Vacuum;//Water;//mList_->Scin;

  G4ThreeVector localPosWC[SegNumWC];
  G4ThreeVector globalPosWC[SegNumWC];

  //  localPosWC[5] = G4ThreeVector(0, 0, -zdiff);
  //  globalPosWC[5] = geomMan.Local2GlobalPos(lnum_bh1, localPosWC[5]);

  G4double WCpos = 4600.;
  G4double WCtpos = 0.;
  G4int WCstart = 7;
  localPosWC[WCstart] = G4ThreeVector(WCtpos, 0, 0);
  globalPosWC[WCstart] = G4ThreeVector(WCpos+zdiff, WCtpos, 0);

  for (int i=WCstart-1; i>=0; i--) {
    double localx = localPosWC[i+1].x() - WCframeY/2;
    double localz;
    if(i%2 == 0)  localz = -zdiff;
    else  localz = zdiff;

    localPosWC[i] = G4ThreeVector(localx, 0, localz);
    //    globalPosWC[i] = geomMan.Local2GlobalPos(lnum_bh1, localPosWC[i]);
    globalPosWC[i] = G4ThreeVector(WCpos+localz, localx, 0);
  }

  for (int i=WCstart+1; i<=SegNumWC; i++) {
    double localx = localPosWC[i-1].x() + WCframeY/2;
    double localz;
    if (i%2 == 0)
      localz = -zdiff;
    else
      localz = zdiff;

    localPosWC[i] = G4ThreeVector(localx, 0, localz);
    //    globalPosWC[i] = geomMan.Local2GlobalPos(lnum_bh1, localPosWC[i]);
    globalPosWC[i] = G4ThreeVector(WCpos+localz, localx, 0);
  }

  char buf[100];

  G4double RotateAngleWC = 70;//geomMan.GetRotAngle2( lnum_bh1 );
  G4RotationMatrix RMWC; 
  RMWC.rotateZ( RotateAngleWC*degree );

  G4Box *solidWCBox = new G4Box( "solidWCBox", WCframeX/2.*mm,  WCframeY/2.*mm,  WCframeZ/2.*mm );
  G4Box *solidWCArea = new G4Box( "solidWCArea", WCradX/2.*mm,  WCradY/2.*mm,  WCradZ/2.*mm );
  //  G4SubtractionSolid *solidWCFrame = new G4SubtractionSolid( "solidWCFrame", solidWCBox,  solidWCArea );
  G4Box *solidWCLayer = new G4Box( "solidWCLayer", 0.01/2.*mm,  WCradY/2.*mm,  WCradZ/2.*mm );

  G4LogicalVolume *logWCBox = new G4LogicalVolume( solidWCBox, WCMater, "logicWCBox", 0, 0, 0 );
  G4LogicalVolume *logWCArea = new G4LogicalVolume( solidWCArea, WCMater, "logicWCArea", 0, 0, 0 );
  //  G4LogicalVolume *logWCFrame = new G4LogicalVolume( solidWCFrame, WCMater, "logicWCFrame", 0, 0, 0 );
  G4LogicalVolume *logWCLayer = new G4LogicalVolume( solidWCLayer, WCMater, "logicWCLayer", 0, 0, 0 );

  G4VPhysicalVolume *physWCBox[SegNumWC];
  G4VPhysicalVolume *physWCArea[SegNumWC];
  //  G4VPhysicalVolume *physWCFrame[SegNumWC];
  G4VPhysicalVolume *physWCLayer[SegNumWC];

  for (int i=0; i<SegNumWC; i++) {
  
    globalPosWC[i].rotateZ(70*degree);

    sprintf(buf, "WCBox_%d", i); 
    physWCBox[i] =
      new G4PVPlacement( G4Transform3D(RMWC,globalPosWC[i]), 
			 buf, logWCBox, pMother, false, 0 );
//     sprintf(buf, "WCFrame_%d", i); 
//     physWCFrame[i] =
//       new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
// 			 buf, logWCFrame, physWCBox[i], false, 0 );
    sprintf(buf, "WCArea_%d", i); 
    physWCArea[i] =
      new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
			 buf, logWCArea, physWCBox[i], false, 0 );
    sprintf(buf, "WCLayer_%d", i); 
    physWCLayer[i] =
      new G4PVPlacement( 0, G4ThreeVector( -10.0*mm, 0.0*mm, 0.0*mm ),
			 buf, logWCLayer, physWCArea[i], false, i );

    logWCBox->SetVisAttributes(G4VisAttributes::Invisible);
    //    logWCFrame->SetVisAttributes(G4VisAttributes(true,G4Colour(0.0, 1.0, 1.0)));
    logWCArea->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 1.0)));
    logWCLayer->SetVisAttributes(G4VisAttributes::Invisible);
  }

  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  WCSD *wcSD = new WCSD("WC");
  SDMan->AddNewDetector(wcSD);
  logWCLayer->SetSensitiveDetector(wcSD);

}//MakeWaterCounter

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DetectorConstructionStsGas::
IsVolumeStopper( G4VPhysicalVolume *physVol ) const
{
  G4String name = physVol->GetName();
  if( name=="physQ1" || name=="physQ2" || name=="physDy" ) return true;
  else
    return false;
}
