/*
  "S2SDetectorConstruction.hh"
  
  This was originally coded by S.Kanatsuki.
  Modified by Toshi Gogami , 11Nov2014.
*/


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
// $Id: ExN03DetectorConstruction.hh,v 1.7 2006/06/29 17:48:32 gunter Exp $
// GEANT4 tag $Name: geant4-08-03-patch-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef S2SDetectorConstruction_h
#define S2SDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "ConfMan.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;

struct MaterialList;
//class s2sDetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class S2SDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    S2SDetectorConstruction();
  S2SDetectorConstruction(ConfMan*);
   ~S2SDetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

  private:
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World

     G4UniformMagField* magField;      //pointer to the magnetic field
     
  //     s2sDetectorMessenger* detectorMessenger;  //pointer to the Messenger
      
  private:
  //  G4Material* defaultMaterial;
  G4Material* DMaterial;

  MaterialList* mList_;
  MaterialList* DefineMaterial();
  //  void DefineMaterials();

  void MakeField();
  //     void ConstructCalorimeter(G4VPhysicalVolume* pMother);     
  void ConstructTarget (G4VPhysicalVolume* pMother);
  void ConstructQ1(G4VPhysicalVolume* pMother);
  void ConstructQ2(G4VPhysicalVolume* pMother);
  void ConstructD1(G4VPhysicalVolume* pMother);

  void SetTarget(G4VPhysicalVolume* pMother);
  void MakePositionDetector(G4VPhysicalVolume* pMother);
  void MakeTOFCounter(G4VPhysicalVolume* pMother);
  void MakeAerogelCounter(G4VPhysicalVolume* pMother);
  void MakeWaterCounter(G4VPhysicalVolume* pMother);
  void MakeSlits(G4VPhysicalVolume* pMother);
  G4bool IsVolumeStopper( G4VPhysicalVolume *physVol ) const ;
  ConfMan* confman;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

