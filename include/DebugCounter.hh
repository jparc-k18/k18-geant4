// -*- C++ -*-

#ifndef DEBUG_COUNTER_HH
#define DEBUG_COUNTER_HH

#include <map>

#include <G4String.hh>

//_____________________________________________________________________________
namespace debug
{
class ObjectCounter
{
public:
  static const G4String& ClassName();
  static ObjectCounter& GetInstance();
  ~ObjectCounter();

private:
  ObjectCounter();
  ObjectCounter(const ObjectCounter&);
  ObjectCounter& operator =(const ObjectCounter&);

private:
  using ObjectMap = std::map<G4String, G4int>;
  using ObjectIter = ObjectMap::const_iterator;
  ObjectMap m_map;

public:
  void check(const G4String& arg="") const;
  void print(const G4String& arg="") const;

public:
  static void decrease(const G4String& key);
  static void increase(const G4String& key);
};

//_____________________________________________________________________________
inline const G4String&
ObjectCounter::ClassName()
{
  static G4String s_name("ObjectCounter");
  return s_name;
}

//_____________________________________________________________________________
inline ObjectCounter&
ObjectCounter::GetInstance()
{
  static ObjectCounter s_instance;
  return s_instance;
}

//_____________________________________________________________________________
inline void
ObjectCounter::decrease([[maybe_unused]]const G4String& key)
{
#ifdef MemoryLeak
  --(GetInstance().m_map[key]);
#endif
}

//_____________________________________________________________________________
inline void
ObjectCounter::increase([[maybe_unused]]const G4String& key)
{
#ifdef MemoryLeak
  ++(GetInstance().m_map[key]);
#endif
}

}

#endif
