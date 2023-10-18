// -*- C++ -*-

#ifndef MATERIAL_LIST_HH
#define MATERIAL_LIST_HH

#include <map>

#include <G4String.hh>

class G4Element;
class G4Material;

//_____________________________________________________________________________
struct MaterialList
{
public:
  static G4String ClassName();
  MaterialList();
  ~MaterialList();

  std::map<G4String, G4Element*> element_map;
  std::map<G4String, G4Material*> material_map;
  G4Material* at(const G4String& key) const;

  G4Element *elH;
  G4Element *elHe;
  G4Element *elLi;
  G4Element *elB10;
  G4Element *elB11;
  G4Element *elC;
  G4Element *elN;
  G4Element *elO;
  G4Element *elNa;
  G4Element *elAl;
  G4Element *elSi;
  G4Element *elP;
  G4Element *elS;
  G4Element *elAr;
  G4Element *elTi;
  G4Element *elCr;
  G4Element *elMn;
  G4Element *elFe;
  G4Element *elNi;
  G4Element *elCu;
  G4Element *elZn;
  G4Element *elGe;
  G4Element *elMo;
  G4Element *elI;
  G4Element *elCs;
  G4Element *elW;
  G4Element *elPt;
  G4Element *elPb;
  G4Element *elBi;

  G4Material *HeLiq;
  G4Material *Li;
  G4Material *Be9;
  G4Material *B10;
  G4Material *B11;
  G4Material *C;
  G4Material *Al;
  G4Material *Si;
  G4Material *ArGas;
  G4Material *Ti;
  G4Material *Ni;
  G4Material *Cu;
  G4Material *Ge;
  G4Material *W;
  G4Material *Pt;
  G4Material *Pb;

  G4Material *Water;
  G4Material *BGO;
  G4Material *PWO;
  G4Material *CsI;
  G4Material *NaI;
  G4Material *SUS316L;
  G4Material *Aerogel;
  G4Material *Scin;
  G4Material *Polyethylene;
  G4Material *Acrylic;
  G4Material *LiO;
  G4Material *LiN;
  G4Material *Brass;
  G4Material *PET;
  G4Material *LH;
  G4Material *LD;
  G4Material *MethaneGas;
  G4Material *EthaneGas;
  G4Material *IsoButaneGas;
  G4Material *P10Gas;
  G4Material *Ar50Ethane50Gas;
  G4Material *Ar80IsoButane20Gas;

private:
  MaterialList(const MaterialList&);
  MaterialList& operator=(const MaterialList&);

};

//_____________________________________________________________________________
inline G4String
MaterialList::ClassName()
{
  static const G4String s_name("MaterialList");
  return s_name;
}

#endif
