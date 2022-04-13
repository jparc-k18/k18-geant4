/*
  K18PhysicsList.cc
  2007/4  K.Shirotori
*/

#include "globals.hh"
#include "s2sPhysicsList.hh"
//#include "Transportation.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Ions.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4KL3DecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4Material.hh"
#include "G4OpticalPhoton.hh"
#include "G4Ions.hh"
#include "G4ios.hh"
#include <iomanip>

#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

//#include "ConfMan.hh"

s2sPhysicsList::s2sPhysicsList()
  : G4VUserPhysicsList()
{
  defaultCutValue = 0.1*cm;
}

s2sPhysicsList::~s2sPhysicsList()
{
}

void s2sPhysicsList::ConstructParticle()
{
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructHeavyIon();
}

void s2sPhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  // gamma
  G4Gamma::GammaDefinition();
  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

void s2sPhysicsList::ConstructLeptons()
{
  // leptons
  //  e+/-
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  // mu+/-
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  // nu_e
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  // nu_mu
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

void s2sPhysicsList::ConstructMesons()
{
  //  mesons
  //    light mesons
  // Pion +/- 0
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  // Kain +/- 0
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
  // Eta
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();

  ConstructStableMesons();

  //For Kaon Decay
  G4DecayTable *decayTable;
  G4VDecayChannel *mode;
  G4ParticleDefinition *particle;

#if 1
  //K->Pi+ Pi- Pi-
  particle 
    = new G4ParticleDefinition(
			       "kaon1-",    0.493677*GeV, 5.315e-14*MeV,    -1.*eplus, 
			       0,              -1,             0,          
			       1,              -1,             0,             
			       "meson",         0,             0,            -321,
			       false,          0.0,           NULL,
			       false,         "kaon" );

  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel( "kaon1-", 1.0, 3, "pi+", "pi-", "pi-" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  //K->e- Pi0 nue
  particle 
    = new G4ParticleDefinition (            
				"kaon2-",    0.493677*GeV, 5.315e-14*MeV,    -1.*eplus, 
				0,              -1,             0,          
				1,              -1,             0,             
				"meson",         0,             0,           -321,
				false,         0.0,            NULL,
				false,        "kaon" );
  
  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel("kaon2-",1.0, 3,
				       "e-","pi0", "anti_nu_e");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);  

  //K->Mu- Pi0 nuMu
  particle 
    = new G4ParticleDefinition (            
				"kaon3-",    0.493677*GeV, 5.315e-14*MeV,    -1.*eplus, 
				0,              -1,             0,          
				1,              -1,             0,             
				"meson",         0,             0,            -321,
				false,         0.0,            NULL,
				false,        "kaon" );
  
  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel("kaon3-",1.0, 3,
				       "mu-","pi0", "anti_nu_mu");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  //K->Pi- Pi0 Pi0
  particle 
    = new G4ParticleDefinition (            
				"kaon4-",    0.493677*GeV, 5.315e-14*MeV,    -1.*eplus, 
				0,              -1,             0,          
				1,              -1,             0,             
				"meson",         0,             0,            -321,
				false,         0.0,            NULL,
				false,        "kaon" );
  
  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel("kaon4-",1.0, 3,
				       "pi-","pi0", "pi0");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);
#endif

}

void s2sPhysicsList::ConstructBaryons()
{
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();

  G4Deuteron::DeuteronDefinition();

  ConstructStableHyperons();
  ConstructDecayingHyperons();


  // Construct tagged proton
  G4Ions* particle;
  particle = new G4Ions(
		"proton_d",   0.9382723*GeV,       0.0*MeV,       eplus, 
		    1,              +1,             0,          
		    1,              +1,             0,             
	     "baryon",               0,            +1,        2212,
		 true,            -1.0,          NULL,
		false,       "nucleon",             0,
                  0.0 
		);

}

void s2sPhysicsList::ConstructHeavyIon()
{
  //Ions
  // Deuteron
  G4Deuteron::DeuteronDefinition();
  // Triton
  G4Triton::TritonDefinition();
  // Alpha
  G4Alpha::AlphaDefinition();
}

void s2sPhysicsList::ConstructProcess()
{

  //  AddTransportationSks();

  //  ConfMan *confMan = ConfMan::GetConfManager();
  
  //  if (confMan->EMFlag() == 1)
    ConstructEM();
  
    //  if (confMan->DecayFlag() == 1)
    //    ConstructDecay();
  //ConstructHadronic();
}

///////Transportation
/*
void s2sPhysicsList::AddTransportationSks()
{
  Transportation* theTransportationProcess= new Transportation();

  // loop over all particles in G4ParticleTable
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if ( !particle->IsShortLived() ) {
      // Add transportation process for all particles other than  "shortlived"
      if ( pmanager == 0) {
        // Error !! no process manager
        G4Exception("s2sPhysicsList::AddTransportation : no process manager!");
      } 
      else {
        // add transportation with ordering = ( -1, "first", "first" )
        pmanager->AddProcess(theTransportationProcess);
        pmanager->SetProcessOrderingToFirst(theTransportationProcess,
                                            idxAlongStep);
        pmanager ->SetProcessOrderingToFirst(theTransportationProcess,
                                             idxPostStep);
      }
    }
    else {
      // shortlived particle case
    }
  }
}
*/
///////Cut
void s2sPhysicsList::SetCuts()
{
  // Suppress error message int case e/gamma/proton do not exist
  G4int temp = GetVerboseLevel();
  // Retrive verbose level
  SetVerboseLevel(temp);

  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");

  //SetCutsWithDefault();  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition *particle=theParticleIterator->value();
    particle->SetApplyCutsFlag( true );
    //////////////////////////////////////////////////////////////////////
    //    G4cout << particle->GetParticleName() << " ==> ApplyCutFlag = "
    //     << particle->GetApplyCutsFlag() << G4endl;
    //////////////////////////////////////////////////////////////////////
  }
}

///////EM
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

#include "G4UserSpecialCuts.hh"

void s2sPhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pManager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if( particleName == "gamma" ){
      pManager->AddProcess( new G4UserSpecialCuts(), -1, -1, 1 );    
      pManager->AddDiscreteProcess( new G4PhotoElectricEffect() );
      pManager->AddDiscreteProcess( new G4ComptonScattering() );
      pManager->AddDiscreteProcess( new G4GammaConversion() );
    }
    else if( particleName == "e-" ){
      pManager->AddProcess( new G4UserSpecialCuts(),    -1, -1, 1 );
      pManager->AddProcess( new G4MultipleScattering(), -1,  1, 2 );
      pManager->AddProcess( new G4eIonisation(),        -1,  2, 3 );
      pManager->AddProcess( new G4eBremsstrahlung(),    -1, -1, 4 );
    }
    else if( particleName == "e+" ){
      pManager->AddProcess( new G4UserSpecialCuts(),    -1, -1, 1 );
      pManager->AddProcess( new G4MultipleScattering(), -1,  1, 2 );
      pManager->AddProcess( new G4eIonisation(),        -1,  2, 3 );
      pManager->AddProcess( new G4eBremsstrahlung(),    -1, -1, 4 );
      pManager->AddProcess( new G4eplusAnnihilation(),   1, -1, 5 );
    }
    else if( particleName == "mu+" || particleName == "mu-" ){
      pManager->AddProcess( new G4MultipleScattering(), -1,  1, 1 );
      pManager->AddProcess( new G4MuIonisation(),       -1,  2, 2 ); 
      pManager->AddProcess( new G4MuBremsstrahlung(),   -1, -1, 3 );
      pManager->AddProcess( new G4MuPairProduction(),   -1, -1, 4 );
    }
    else if( !particle->IsShortLived() &&
             particle->GetPDGCharge()!=0.0 &&
	     !( particleName=="chargedgeantino"
		|| particleName=="antichargedgeantino")){
      pManager->AddProcess( new G4UserSpecialCuts(),      -1, -1, 1 );
      pManager->AddProcess( new G4MultipleScattering(),   -1,  1, 2 );
      pManager->AddProcess( new G4hIonisation(),          -1,  2, 3 );
    }
  }
}

///////Hadron
#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4LElastic.hh"

//Pion
#include "G4PionPlusInelasticProcess.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4LEPionMinusInelastic.hh"

//Kaon
#include "G4KaonPlusInelasticProcess.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4LEKaonZeroLInelastic.hh"

//Nucleon
#include "G4ProtonInelasticProcess.hh"
#include "G4LEProtonInelastic.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4LENeutronInelastic.hh"

void s2sPhysicsList:: ConstructHadronic()
{
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  G4LElastic* theElasticModel = new G4LElastic;
  theElasticProcess->RegisterMe(theElasticModel);
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    
    //Pion      
    if( particleName == "pi+" ) {
      pmanager->AddProcess( theElasticProcess );
      G4PionPlusInelasticProcess* theInelasticProcess 
	= new G4PionPlusInelasticProcess( "inelastic" );
      G4LEPionPlusInelastic* theLEInelasticModel 
	= new G4LEPionPlusInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "pi-" ) {
      pmanager->AddProcess( theElasticProcess );
      G4PionMinusInelasticProcess* theInelasticProcess 
	= new G4PionMinusInelasticProcess( "inelastic" );
      G4LEPionMinusInelastic* theLEInelasticModel 
	= new G4LEPionMinusInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    //Kaon
    else if( particleName == "kaon+" ) {
      pmanager->AddProcess( theElasticProcess );
      G4KaonPlusInelasticProcess* theInelasticProcess 
	= new G4KaonPlusInelasticProcess( "inelastic" );
      G4LEKaonPlusInelastic* theLEInelasticModel 
	= new G4LEKaonPlusInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "kaon-" ) {
      pmanager->AddProcess( theElasticProcess );
      G4KaonMinusInelasticProcess* theInelasticProcess 
	= new G4KaonMinusInelasticProcess( "inelastic" );
      G4LEKaonMinusInelastic* theLEInelasticModel 
	= new G4LEKaonMinusInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "kaon0S" ) {
      pmanager->AddProcess( theElasticProcess );
      G4KaonZeroSInelasticProcess* theInelasticProcess 
	= new G4KaonZeroSInelasticProcess( "inelastic" );
      G4LEKaonZeroSInelastic* theLEInelasticModel 
	= new G4LEKaonZeroSInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "kaon0L" ) {
      pmanager->AddProcess( theElasticProcess );
      G4KaonZeroLInelasticProcess* theInelasticProcess 
	= new G4KaonZeroLInelasticProcess( "inelastic" );
      G4LEKaonZeroLInelastic* theLEInelasticModel 
	= new G4LEKaonZeroLInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    //Nucleon
    else if( particleName == "proton" ) {
      pmanager->AddProcess( theElasticProcess );
      G4ProtonInelasticProcess* theInelasticProcess 
	= new G4ProtonInelasticProcess( "inelastic" );
      G4LEProtonInelastic* theLEInelasticModel 
	= new G4LEProtonInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "neutron" ) {
      pmanager->AddProcess( theElasticProcess );
      G4NeutronInelasticProcess* theInelasticProcess 
	= new G4NeutronInelasticProcess( "inelastic" );
      G4LENeutronInelastic* theLEInelasticModel 
	= new G4LENeutronInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
}

///////Decay
#include "G4Decay.hh"

void s2sPhysicsList::ConstructDecay()
{
  G4Decay *theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pManager = particle->GetProcessManager();
    if( theDecayProcess->IsApplicable(*particle) ){
      pManager->AddProcess( theDecayProcess );
      pManager->SetProcessOrdering( theDecayProcess, idxPostStep );
      pManager->SetProcessOrdering( theDecayProcess, idxAtRest );
    }
  }
}


void s2sPhysicsList::ConstructStableMesons()
{
  //  G4DecayTable* decayTable;
  //  G4VDecayChannel* mode;
  G4ParticleDefinition* particle;


  // skaon- non-decay kaon-
  particle = new G4ParticleDefinition(
               "skaon-",    0.493677*GeV,   5.315e-14*MeV,    -1.*eplus,
                    0,              -1,             0,
                    1,              -1,             0,
              "meson",               0,             0,        -321,
	       true,       0,          NULL);

  // skaon+ non-decay kaon+
  particle = new G4ParticleDefinition(
                 "skaon+",    0.493677*GeV,  5.315e-14*MeV,    +1.*eplus,
                    0,              -1,             0,
                    1,              +1,             0,
              "meson",               0,             0,         321,
		 true,       0,          NULL);


  // spi+ non-decay kaon-
  particle = new G4ParticleDefinition(
               "spi+",    0.1395700*GeV,  2.5284e-14*MeV,    +1.*eplus,
                    0,              -1,             0,
                    2,              +2,            -1,
              "meson",               0,             0,         211,
	       true,       0,          NULL);

  // spi- non-decay kaon-
  particle = new G4ParticleDefinition(
               "spi-",    0.1395700*GeV, 2.5284e-14*MeV,    -1.*eplus,
                    0,              -1,             0,
                    2,              -2,            -1,
              "meson",               0,             0,        -211,
	       true,               0,          NULL);

}

void s2sPhysicsList::ConstructStableHyperons()
{
  G4DecayTable* decayTable;
  G4VDecayChannel* mode;
  G4ParticleDefinition* particle;


  // ssigma+ non-decay sigma+
  particle = new G4ParticleDefinition(
           "ssigma+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   true,                0,          NULL);


  // sigma1+  decay only to sigma+ -> pi+ neutron channel
  particle = new G4ParticleDefinition(
           "sigma1+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   false,                0.0799*ns,          NULL);
  
  decayTable =  new G4DecayTable();
  // sigma+ -> neutron + pi+
  mode = new G4PhaseSpaceDecayChannel("sigma1+",1.0,2,"neutron","pi+");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // sigma2+  decay only to sigma+ -> pi0 proton channel
  particle = new G4ParticleDefinition(
           "sigma2+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   false,                0.0799*ns,          NULL);
  
  decayTable =  new G4DecayTable();
  // sigma+ -> proton + pi0
  mode = new G4PhaseSpaceDecayChannel("sigma2+", 1.0,2,"proton","pi0");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

}      

void s2sPhysicsList::ConstructDecayingHyperons()
{
  G4DecayTable* decayTable;
  G4VDecayChannel* mode;
  G4ParticleDefinition* particle;


  // usigma- life time is 0
  particle = new G4ParticleDefinition(
	    "usigma-",    1.19744*GeV,  4.45e-12*MeV,    -1*eplus,
                    1,              +1,             0,
                    2,              -2,             0,
             "baryon",               0,            +1,        3112,
	    false,       0.0*ns,          NULL);

 //create Decay Table 
  decayTable = new G4DecayTable();
  // create decay channels
  // sigma- -> neutron + pi-
  mode = new G4PhaseSpaceDecayChannel("usigma-",1.00,2,"neutron","pi-");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);


  // usigma+ life time is 0
  particle = new G4ParticleDefinition(
           "usigma+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   false,                0,          NULL);

 //create Decay Table
  decayTable = new G4DecayTable();

  // create decay channels
  G4VDecayChannel** modetable = new G4VDecayChannel*[2];
  // sigma+ -> proton + pi0
  modetable[0] = new G4PhaseSpaceDecayChannel("usigma+",0.516,2,"proton","pi0");
  // sigma+ -> neutron + pi+
  modetable[1] = new G4PhaseSpaceDecayChannel("usigma+",0.483,2,"neutron","pi+");

  for (G4int index=0; index <2; index++ ) decayTable->Insert(modetable[index]);
  delete [] modetable;
  particle->SetDecayTable(decayTable);
  
  // usigma1+  decay only to sigma+ -> pi+ neutron channel
  particle = new G4ParticleDefinition(
           "usigma1+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   false,                0.0*ns,          NULL);
  
  decayTable =  new G4DecayTable();
  // usigma+ -> neutron + pi+
  mode = new G4PhaseSpaceDecayChannel("usigma1+",1.0,2,"neutron","pi+");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // usigma2+  decay only to sigma+ -> pi0 proton channel
  particle = new G4ParticleDefinition(
           "usigma2+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   false,                0.0*ns,          NULL);
  
  decayTable =  new G4DecayTable();
  // sigma+ -> proton + pi0
  mode = new G4PhaseSpaceDecayChannel("usigma2+", 1.0,2,"proton","pi0");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

}      
