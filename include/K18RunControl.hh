#ifndef K18_RUN_CONTROL_HH
#define K18_RUN_CONTROL_HH

#include <cerrno>
#include <cstdlib>
#include <limits>

#include <G4Event.hh>
#include <G4ios.hh>
#include <globals.hh>

namespace K18RunControl
{
// Parse optional integer controls strictly so malformed batch settings fail
// before they can alter event numbering or random streams.
inline G4bool
ReadEnvironmentLong(const char* name, G4long& value)
{
  const char* raw = std::getenv(name);
  if(!raw || !*raw)
    return false;

  errno = 0;
  char* end = nullptr;
  const long parsed = std::strtol(raw, &end, 10);
  if(errno == ERANGE || end == raw || !end || *end != '\0'){
    G4cerr << name << " must be an integer; got '" << raw << "'" << G4endl;
    std::exit(EXIT_FAILURE);
  }
  value = static_cast<G4long>(parsed);
  return true;
}

inline G4long
CheckedAdd(G4long lhs, G4long rhs, const char* label)
{
  const G4long min = std::numeric_limits<G4long>::min();
  const G4long max = std::numeric_limits<G4long>::max();
  if((rhs > 0 && lhs > max-rhs) || (rhs < 0 && lhs < min-rhs)){
    G4cerr << label << " exceeds the G4long range" << G4endl;
    std::exit(EXIT_FAILURE);
  }
  return lhs+rhs;
}

inline G4long
GlobalEventIndex(const G4Event* event)
{
  if(!event){
    G4cerr << "cannot compute a global index for a null event" << G4endl;
    std::exit(EXIT_FAILURE);
  }
  G4long offset = 0;
  ReadEnvironmentLong("K18G4_EVENT_OFFSET", offset);
  return CheckedAdd(offset, event->GetEventID(), "global event index");
}
}

#endif
