// -*- C++ -*-

#include "ConfMan.hh"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <libgen.h>
#include <unistd.h>

#include <TString.h>

#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"
#include "HistMan.hh"

//_____________________________________________________________________________
ConfMan::ConfMan()
  : m_conf_key("CONF"),
    m_conf_dir(),
    m_conf_buf()
{
}

//_____________________________________________________________________________
ConfMan::~ConfMan()
{
}

//_____________________________________________________________________________
G4bool
ConfMan::Initialize(const G4String& file_name)
{
  static const std::string funcname = "[ConfMan::Initialize]";

  m_file[m_conf_key] = file_name;
  m_conf_dir = ::dirname(const_cast<char*>(m_file[m_conf_key].data()));
  m_conf_buf.clear();
  m_conf_buf += "\n";

  std::ifstream ifs(file_name);
  if(!ifs.is_open()){
    G4cerr << FUNC_NAME << " file open fail : " << file_name << G4endl;
    return false;
  }

  G4String line;
  while(ifs.good() && std::getline(ifs, line)){
    if(line.empty() || line[0]=='#') continue;
    m_conf_buf += line + "\n";
    std::istringstream iss(line);
    G4String key, val;
    iss >> key >> val;
    if(key.empty() || val.empty())
      continue;
    if(key.back() == ':')
      key.pop_back();
    G4cout << " key = "   << std::setw(20) << std::left << key
	   << " value = " << std::setw(30) << std::left << val
	   << G4endl;

    m_file[key]   = FilePath(val);
    m_string[key] = val;
    m_double[key] = std::strtod(val, nullptr);
    m_int[key]    = std::strtol(val, nullptr, 10);
    m_bool[key]   = static_cast<G4bool>(std::strtol(val, nullptr, 10));
  }

  return InitializeParameterFiles();
}

//_____________________________________________________________________________
bool
ConfMan::InitializeParameterFiles()
{
  return (true
          && InitializeParameter<DCGeomMan>("DCGEO")
          // && InitializeParameter<BeamMan>("BEAM")
          && InitializeParameter<DetSizeMan>("DSIZE")
          && InitializeParameter<HistMan>("HIST")
          // && InitializeParameter<JamMan>("JAM")
          // && InitializeParameter<IncMan>("INC")
          );
}

//_____________________________________________________________________________
G4String
ConfMan::FilePath(const G4String& src) const
{
  std::ifstream tmp(src);
  if (tmp.good())
    return src;
  else
    return m_conf_dir + "/" + src;
}
