// -*- C++ -*-

#include "DebugCounter.hh"

#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>

#include <globals.hh>

#include "FuncName.hh"

namespace debug
{

//_____________________________________________________________________________
ObjectCounter::ObjectCounter()
  : m_map()
{
}

//_____________________________________________________________________________
ObjectCounter::~ObjectCounter()
{
}

//_____________________________________________________________________________
void
ObjectCounter::check([[maybe_unused]]const G4String& arg) const
{
#ifdef MemoryLeak
  G4bool has_leak = false;
  ObjectIter itr, end=m_map.end();
  for(itr=m_map.begin(); itr!=end; ++itr){
    if(itr->second!=0) has_leak = true;
  }
  if(has_leak)
    print(arg+" "+FUNC_NAME);
#endif
}

//_____________________________________________________________________________
void
ObjectCounter::print(const G4String& arg) const
{
  G4cout << "#DCounter " << FUNC_NAME << " " << arg << G4endl;
  ObjectIter itr, end=m_map.end();
  for(itr=m_map.begin(); itr!=end; ++itr){
    G4cout << std::setw(20) << std::left
		<< itr->first << " : " << itr->second << G4endl;
  }
}

}
