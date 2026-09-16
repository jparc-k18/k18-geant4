// -*- C++ -*-

#ifndef ConfMan_h
#define ConfMan_h 1

#include <exception>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <G4String.hh>
#include <G4ios.hh>

#include <TMath.h>

//_____________________________________________________________________________
class ConfMan
{
public:
  static G4String ClassName();
  static ConfMan& GetInstance();
  ~ConfMan();
private:
  ConfMan();
  ConfMan(const ConfMan&);
  ConfMan& operator =(const ConfMan&);

private:
  using StrList = std::map<G4String, G4String>;
  using DoubleList = std::map<G4String, G4double>;
  using IntList = std::map<G4String, G4int>;
  using BoolList = std::map<G4String, G4bool>;
  G4String   m_conf_key;
  G4String   m_conf_dir;
  G4String   m_conf_buf;
  StrList    m_file;
  StrList    m_string;
  DoubleList m_double;
  IntList    m_int;
  BoolList   m_bool;
  std::vector<G4String> m_config_warnings;
  std::vector<G4String> m_config_errors;

public:
  G4bool Initialize(const G4String& file_name);
  const G4String& ConfBuf() const { return m_conf_buf; }
  const G4String& ConfPath() const { return m_file.at(m_conf_key); }
  template <typename T>
  static const T& Get(const G4String& key);

private:
  G4String FilePath(const G4String& src) const;
  G4bool   ValidateConfiguration();
  void     PrintConfigurationDiagnostics() const;
  void     StoreValue(const G4String& key, const G4String& value);
  void     EraseValue(const G4String& key);
  G4bool   InitializeParameterFiles();
  template <typename T>
  G4bool    InitializeParameter();
  template <typename T>
  G4bool    InitializeParameter(const G4String& key);
  template <typename T>
  G4bool    InitializeParameter(const G4String& key1,
                                const G4String& key2);
  G4bool   ShowResult(G4bool status, const G4String& name) const;
};

//_____________________________________________________________________________
inline G4String
ConfMan::ClassName()
{
  static const G4String s_name("ConfMan");
  return s_name;
}

//_____________________________________________________________________________
inline ConfMan&
ConfMan::GetInstance()
{
  static ConfMan s_instance;
  return s_instance;
}

//_____________________________________________________________________________
template <>
inline const G4String&
ConfMan::Get<G4String>(const G4String& key)
{
  return GetInstance().m_string[key];
}

//_____________________________________________________________________________
template <>
inline const G4double&
ConfMan::Get<G4double>(const G4String& key)
{
  return GetInstance().m_double[key];
}

//_____________________________________________________________________________
template <>
inline const G4int&
ConfMan::Get<G4int>(const G4String& key)
{
  return GetInstance().m_int[key];
}

//_____________________________________________________________________________
template <>
inline const G4bool&
ConfMan::Get<G4bool>(const G4String& key)
{
  return GetInstance().m_bool[key];
}

//_____________________________________________________________________________
inline G4bool
ConfMan::ShowResult(G4bool status, const G4String& name) const
{
  if(status)
    std::cout << std::setw(20) << std::left
	      << " ["+name+"]"
	      << "-> Initialized" << std::endl;
  else
    std::cout << std::setw(20) << std::left
	      << " ["+name+"]"
	      << "-> Failed" << std::endl;
  return status;
}

//_____________________________________________________________________________
template <typename T>
inline G4bool
ConfMan::InitializeParameter()
{
  try {
    return
      ShowResult(T::GetInstance().Initialize(),
                 T::GetInstance().ClassName());
  } catch(const std::exception& error) {
    G4cerr << "#E [ConfMan::InitializeParameter] "
           << T::GetInstance().ClassName() << ": " << error.what()
           << G4endl;
  } catch(...) {
    G4cerr << "#E [ConfMan::InitializeParameter] "
           << T::GetInstance().ClassName() << ": unknown exception"
           << G4endl;
  }
  return ShowResult(false, T::GetInstance().ClassName());
}

//_____________________________________________________________________________
template <typename T>
inline G4bool
ConfMan::InitializeParameter(const G4String& key)
{
  try {
    return
      ShowResult(T::GetInstance().Initialize(m_file[key]),
                 T::GetInstance().ClassName());
  } catch(const std::exception& error) {
    G4cerr << "#E [ConfMan::InitializeParameter] "
           << T::GetInstance().ClassName() << ": " << error.what()
           << G4endl;
  } catch(...) {
    G4cerr << "#E [ConfMan::InitializeParameter] "
           << T::GetInstance().ClassName() << ": unknown exception"
           << G4endl;
  }
  return ShowResult(false, T::GetInstance().ClassName());
}

//_____________________________________________________________________________
template <typename T>
inline G4bool
ConfMan::InitializeParameter(const G4String& key1,
                             const G4String& key2)
{
  try {
    return
      ShowResult(T::GetInstance().Initialize(m_file[key1],
                                             m_file[key2]),
                 T::GetInstance().ClassName());
  } catch(const std::exception& error) {
    G4cerr << "#E [ConfMan::InitializeParameter] "
           << T::GetInstance().ClassName() << ": " << error.what()
           << G4endl;
  } catch(...) {
    G4cerr << "#E [ConfMan::InitializeParameter] "
           << T::GetInstance().ClassName() << ": unknown exception"
           << G4endl;
  }
  return ShowResult(false, T::GetInstance().ClassName());
}

#endif
