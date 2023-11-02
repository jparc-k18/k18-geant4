/*
  FieldMan.cc
*/

#include "FieldMan.hh"
#include "SksFieldMap.hh"
#include "FieldElements.hh"
#include "RadDeg.hh"

FieldMan *FieldMan::fMan_ = 0;

const double Delta = 0.1;

FieldMan::FieldMan()
  : Sksmap_(0)
{}

FieldMan::~FieldMan()
{
  delete Sksmap_;
}

FieldMan & FieldMan::GetInstance( void )
{
  if( !fMan_ ){
    fMan_ = new FieldMan();
  }
  return *fMan_;
}

bool FieldMan::Initialize( void )
{
  SksFieldMap *pSks = new SksFieldMap( filename_.c_str() );
  if( pSks ){
    delete Sksmap_;
    Sksmap_ = pSks;
  }
  if( Sksmap_ )
    return Sksmap_->Initialize();
  else
    return false;
}

ThreeVector FieldMan::GetField( const ThreeVector & position ) const
{
  ThreeVector field( 0.0, 0.0, 0.0 );
  if( Sksmap_ ){
    double p[3], b[3], pt[3];
    p[0]=position.x(); p[1]=position.y(); p[2]=position.z(); 
    //std::cout<<"[GetField] position"<<position<<std::endl;

    //Map2 region -> conversion
    ThreeVector x(p[0],p[1],p[2]);
    ThreeVector MagOrg(600/tan(55*Deg2Rad),-600,0);
    if( p[1]>-tan(70*Deg2Rad)*p[0] ){
      x = x - MagOrg;
      x.rotateZ(-70*Deg2Rad);
      x = MagOrg + x;
      pt[0]=x.x(); pt[1]=x.y(); pt[2]=x.z();
    }
    else{
      pt[0]=x.x(); pt[1]=x.y(); pt[2]=x.z();
    }      

    //conversion finish
    //std::cout<<"[GetField] positionx"<<x<<std::endl;

    if( Sksmap_->GetFieldValue( pt, b ) ){
      field.setX(b[0]); field.setY(b[1]); field.setZ(b[2]); 
    }
    if( p[1]>-tan(70*Deg2Rad)*p[0] ){
      field.rotateZ(70*Deg2Rad);
    }
    //std::cout<<"[GetField] b("<<b[0]<<","<<b[1]<<","<<b[2]<<")"<<std::endl;
    //std::cout<<"[GetField] field"<<field<<std::endl;

  }

#if 1
  FEIterator end = elemList_.end();
  for( FEIterator itr=elemList_.begin(); itr!=end; ++itr ){
    if( (*itr)->ExistField( position ) )
      field += (*itr)->GetField( position );
  }   
#endif

  return field;
}

ThreeVector FieldMan::GetdBdX( const ThreeVector & position ) const
{
  ThreeVector p1 = position + ThreeVector( Delta, 0., 0. );
  ThreeVector p2 = position - ThreeVector( Delta, 0., 0. );

  ThreeVector B1 = GetField( p1 );
  ThreeVector B2 = GetField( p2 );

 
  return 0.5*(B1-B2)/Delta;
} 

ThreeVector FieldMan::GetdBdY( const ThreeVector & position ) const
{
  ThreeVector p1 = position + ThreeVector( 0., Delta,  0. );
  ThreeVector p2 = position - ThreeVector( 0., Delta, 0. );

  ThreeVector B1 = GetField( p1 );
  ThreeVector B2 = GetField( p2 );

  return 0.5*(B1-B2)/Delta;
} 

ThreeVector FieldMan::GetdBdZ( const ThreeVector & position ) const
{
  ThreeVector p1 = position + ThreeVector( 0., 0., Delta );
  ThreeVector p2 = position - ThreeVector( 0., 0., Delta );

  ThreeVector B1 = GetField( p1 );
  ThreeVector B2 = GetField( p2 );

  return 0.5*(B1-B2)/Delta;
} 

void FieldMan::cleanupElementsList( void )
{
  elemList_.clear();
}

void FieldMan::AddElement( FieldElements *elem )
{
  elemList_.push_back( elem );
}


double FieldMan::StepSize( const ThreeVector &position,
			   double defStepSize,
			   double MinStepSize ) const
{
  double d=fabs(defStepSize);
  double s=defStepSize/d;
  MinStepSize = fabs(MinStepSize);

  bool flag = true;
  FEIterator end = elemList_.end();
  while ( flag && d>MinStepSize ){
    for( FEIterator itr=elemList_.begin(); itr!=end; ++itr ){
      if( (*itr)->checkRegion( position, d) != FERSurface ){
	flag=false;
	break;
      }
    }
    d *= 0.5;
  }
  if( flag ) return s*MinStepSize;
  else       return s*d;
}
