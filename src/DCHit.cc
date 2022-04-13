/*
  DCHit.cc
  2007/4  K.Shirotori
*/

#include "DCHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

#include <iomanip>

DCHit::DCHit()
  : layerID_(0), wireID_(0), time_(-100.0), trackNo_(0)
{}

G4Allocator<DCHit> DCHitAllocator;

// void DCHit::Draw() const
// {
//   G4VVisManager *pVisManager = G4VVisManager::GetConcreteInstance();
//   if( pVisManager ){
//     G4Circle circle( pos_ );
//     circle.SetScreenSize( 3.0 );              // in pixels
//     circle.SetFillStyle( G4Circle::filled );
//     G4Colour colour( 0., 1., 0. );            // green
//     G4VisAttributes attribs( colour );
//     circle.SetVisAttributes( attribs );
//     pVisManager->Draw( circle );
//   }
// }

// void DCHit::Print() const
// {
//   G4cout << "DCHit " << std::setw(2) << layerID_ << "-"
// 	 << std::setw(3) << wireID_ << " > Time: "
// 	 << time_/ns << " [ns] Pos: ( " << xl_/mm << " , "
// 	 << yl_/mm << " )" << G4endl;
// }
