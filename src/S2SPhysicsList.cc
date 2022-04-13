/*
  "S2SPhysicsList.cc"

  Toshiyuki Gogami , 10Nov2014
  ( Modified from Ref. N03 )
*/

//#include "globals.hh"
#include "S2SPhysicsList.hh"
#include "S2STransportation.hh"
#include "ConfMan.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
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
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonConstructor.hh"


S2SPhysicsList::S2SPhysicsList()
  : G4VUserPhysicsList(), 
    defaultCutValue(0.1*cm),
    emflag(1), decayflag(0), hadronflag(0)
{
  defaultCutValue = 0.1*cm;
  emflag     = 1;
  decayflag  = 0;
  hadronflag = 0;
}

S2SPhysicsList::S2SPhysicsList(ConfMan* conf)
  : G4VUserPhysicsList(),
    defaultCutValue(0.1*cm),
    emflag(1), decayflag(0), hadronflag(0)
{
  defaultCutValue = 0.1*cm;
  emflag     = conf->GetEMFlag();
  decayflag  = conf->GetDecayFlag();
  hadronflag = conf->GetHadronFlag();
}

S2SPhysicsList::~S2SPhysicsList()
{}

void S2SPhysicsList::ConstructParticle()
{
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  //ConstructHeavyIon();
}

void S2SPhysicsList::ConstructBosons()
{
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();
  
  ///  // pseudo-particles
  ///  G4Geantino::GeantinoDefinition();
  ///  G4ChargedGeantino::ChargedGeantinoDefinition();
  ///  // gamma
  ///  G4Gamma::GammaDefinition();
  ///  // optical photon
  ///  G4OpticalPhoton::OpticalPhotonDefinition();
}

void S2SPhysicsList::ConstructLeptons()
{
  
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();
  
  ///  // leptons
  ///  //  e+/-
  ///  G4Electron::ElectronDefinition();
  ///  G4Positron::PositronDefinition();
  ///  // mu+/-
  ///  G4MuonPlus::MuonPlusDefinition();
  ///  G4MuonMinus::MuonMinusDefinition();
  ///  // nu_e
  ///  G4NeutrinoE::NeutrinoEDefinition();
  ///  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  ///  // nu_mu
  ///  G4NeutrinoMu::NeutrinoMuDefinition();
  ///  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

void S2SPhysicsList::ConstructMesons()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();
  
  
  // mesons
  // light mesons
  // Pion +/- 0
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  // Kaon +/- 0
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
  // Eta
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  
  //For Kaon Decay
  
  G4DecayTable *decayTable;
  G4VDecayChannel *mode;
  G4ParticleDefinition *particle;

  /*  
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
      = new G4ParticleDefinition ( "kaon4-",    0.493677*GeV, 5.315e-14*MeV,    -1.*eplus, 
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
  */

}

void S2SPhysicsList::ConstructBaryons()
{
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
  
  //  ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();
}

/*void S2SPhysicsList::ConstructHeavyIon()
{
  //Ions
  // Deuteron
  G4Deuteron::DeuteronDefinition();
  // Triton
  G4Triton::TritonDefinition();
  // Alpha
  G4Alpha::AlphaDefinition();
}*/

#include "S2SMinEkineCuts.hh"
#include "G4UserSpecialCuts.hh"
#include "G4StepLimiter.hh"

void S2SPhysicsList::ConstructProcess()
{
  G4StepLimiter* StepLimit = new G4StepLimiter();
  G4UserSpecialCuts* UserCuts = new G4UserSpecialCuts();
  AddTransportation();
  auto theParticleIterator=GetParticleIterator(); //added for ver.10.03
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    pmanager->AddDiscreteProcess(new S2SMinEkineCuts());
    pmanager->AddDiscreteProcess( StepLimit );
    pmanager->AddDiscreteProcess( UserCuts );
  }
  G4int EMFlag    = emflag;
  G4int DecayFlag = decayflag;
  G4int HadronFlag= hadronflag;

  if (EMFlag!=0){ConstructEM();}
  if (DecayFlag!=0){ConstructDecay();}
  if (HadronFlag!=0){ConstructHadronic();}
}

/////Transportation
void S2SPhysicsList::AddTransportation()
{
  S2STransportation* theTransportationProcess= new S2STransportation();

  // loop over all particles in G4ParticleTable
  auto theParticleIterator=GetParticleIterator(); //added for ver.10.03
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if ( !particle->IsShortLived() ) {
      // Add transportation process for all particles other than  "shortlived"
      if ( pmanager == 0) {
        // Error !! no process manager
        //G4Exception("S2SPhysicsList::AddTransportation : no process manager!");
	G4cout << "S2SPhysicsList::AddTransportation : no process manager!" << G4endl;
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

///////Cut
void S2SPhysicsList::SetCuts()
{
  // Suppress error message int case e/gamma/proton do not exist
  G4int temp = GetVerboseLevel();
  // Retrive verbose level
  SetVerboseLevel(temp);

  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
  
  //special for low energy physics
  G4double lowlimit=250*eV;  
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit,100.*GeV);

  auto theParticleIterator=GetParticleIterator(); //added for ver.10.03
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


#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4MollerBhabhaModel.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"

#include "G4UserSpecialCuts.hh"

void S2SPhysicsList::ConstructEM()
{
  auto theParticleIterator=GetParticleIterator(); //added for ver.10.03
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pManager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    
    if( particleName == "gamma" ){
      pManager->AddProcess( new G4UserSpecialCuts(), -1, -1, 1 );    
      pManager->AddDiscreteProcess( new G4PhotoElectricEffect );
      pManager->AddDiscreteProcess( new G4ComptonScattering );
      pManager->AddDiscreteProcess( new G4GammaConversion );
    }
    else if( particleName == "e-" ){
      pManager->AddProcess( new G4UserSpecialCuts,    -1, -1, 1 );
      pManager->AddProcess( new G4eMultipleScattering, -1, 1, 2 );
      pManager->AddProcess( new G4eIonisation,        -1,  2, 3 );
      pManager->AddProcess( new G4eBremsstrahlung,    -1, -1, 4 );
      //pManager->AddProcess( new G4MollerBhabhaModel,   -1,  3, 5 );//added by gogami
    }
    else if( particleName == "e+" ){
      pManager->AddProcess( new G4UserSpecialCuts,    -1, -1, 1 );
      pManager->AddProcess( new G4eMultipleScattering, -1, 1, 2 );
      pManager->AddProcess( new G4eIonisation,        -1,  2, 3 );
      pManager->AddProcess( new G4eBremsstrahlung,    -1,  3, 4 );
      pManager->AddProcess( new G4eplusAnnihilation,   0, -1, 5 );
      //pManager->AddProcess( new G4MollerBhabhaModel,   -1,  4, 6 );//added by gogami
    }
    else if( particleName == "mu+" || particleName == "mu-" ){
      pManager->AddProcess( new G4MuMultipleScattering, -1,  1, 1 );
      pManager->AddProcess( new G4MuIonisation,       -1,  2, 2 ); 
      pManager->AddProcess( new G4MuBremsstrahlung,   -1,  3, 3 );
      pManager->AddProcess( new G4MuPairProduction,   -1,  4, 4 );
    }
    else if( particleName == "proton" ||
	     particleName == "pi-" ||
	     particleName == "pi+" || 
	     particleName == "xi-" ||
	     particleName == "xi0" ) {
      pManager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pManager->AddProcess(new G4hIonisation,         -1, 2, 2);
      pManager->AddProcess(new G4hBremsstrahlung,     -1, 3, 3);
      pManager->AddProcess(new G4hPairProduction,     -1, 4, 4);
    }
    else if ( particleName == "alpha" ||
	      particleName == "He3" ){ 
      pManager->AddProcess( new G4hMultipleScattering,  -1,  1, 1 );
      pManager->AddProcess( new G4ionIonisation,        -1,  2, 2 );
    }
    else if ( particleName == "GenericIon" ){
      pManager->AddProcess( new G4hMultipleScattering,  -1,  1, 1 );
      pManager->AddProcess( new G4ionIonisation,        -1,  2, 2 );
    }
    else if( (!particle->IsShortLived()) &&
	     (particle->GetPDGCharge() != 0.0) &&
	     (particle->GetParticleName() != "chargedgeantino")){
      pManager->AddProcess( new G4hMultipleScattering,  -1,  1, 1 );
      pManager->AddProcess( new G4hIonisation,          -1,  2, 2 );
    }
    /*else if( !particle->IsShortLived() &&
      particle->GetPDGCharge()!=0.0 &&
      !( particleName=="chargedgeantino"
      || particleName=="antichargedgeantino")){
      pManager->AddProcess( new G4UserSpecialCuts(),      -1, -1, 1 );
      pManager->AddProcess( new G4MultipleScattering(),   -1,  1, 2 );
      pManager->AddProcess( new G4hIonisation(),          -1,  2, 3 );
      }*/
  }
}

///////Hadron
#include "G4HadronElasticProcess.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"
//#include "G4HadronInelasticProcess.hh"
//#include "G4LElastic.hh"

////Pion
//#include "G4PionPlusInelasticProcess.hh"
//#include "G4LEPionPlusInelastic.hh"
//#include "G4PionMinusInelasticProcess.hh"
//#include "G4LEPionMinusInelastic.hh"
//
////Kaon
//#include "G4KaonPlusInelasticProcess.hh"
//#include "G4LEKaonPlusInelastic.hh"
//#include "G4KaonMinusInelasticProcess.hh"
//#include "G4LEKaonMinusInelastic.hh"
//#include "G4KaonZeroSInelasticProcess.hh"
//#include "G4LEKaonZeroSInelastic.hh"
//#include "G4KaonZeroLInelasticProcess.hh"
//#include "G4LEKaonZeroLInelastic.hh"
//
////Nucleon
//#include "G4ProtonInelasticProcess.hh"
//#include "G4LEProtonInelastic.hh"
//#include "G4NeutronInelasticProcess.hh"
//#include "G4LENeutronInelastic.hh"
// Inelastic processes:
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

// High energy FTFP model and Bertini cascade
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4CascadeInterface.hh"

// Cross sections
#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4CrossSectionElastic.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4AntiNuclElastic.hh"

#include "G4CrossSectionInelastic.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
//#include "G4GGNuclNuclCrossSection.hh"

#include "G4HadronElastic.hh"
#include "G4HadronCaptureProcess.hh"

// Neutron high-precision models: <20 MeV
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"

// Stopping processes
#include "G4PiMinusAbsorptionBertini.hh"
#include "G4KaonMinusAbsorptionBertini.hh"
#include "G4AntiProtonAbsorptionFritiof.hh"


void S2SPhysicsList:: ConstructHadronic()
{
  //Elastic models
  const G4double elastic_elimitPi = 1.0*GeV;
  
  G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  G4HadronElastic* elastic_lhep1 = new G4HadronElastic();
  elastic_lhep1->SetMaxEnergy( elastic_elimitPi );
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE(); 
  elastic_he->SetMinEnergy( elastic_elimitPi );
  
  // Inelastic scattering
  const G4double theFTFMin0 =    0.0*GeV;
  const G4double theFTFMin1 =    4.0*GeV;
  const G4double theFTFMax =   100.0*TeV;
  const G4double theBERTMin0 =   0.0*GeV;
  const G4double theBERTMin1 =  19.0*MeV;
  const G4double theBERTMax =    5.0*GeV;
  //const G4double theHPMin =      0.0*GeV;
  //const G4double theHPMax =     20.0*MeV;
  
  G4FTFModel * theStringModel = new G4FTFModel;
  G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay( new G4LundStringFragmentation );
  theStringModel->SetFragmentationModel( theStringDecay );
  G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
  G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

  G4TheoFSGenerator * theFTFModel0 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel0->SetHighEnergyGenerator( theStringModel );
  theFTFModel0->SetTransport( theCascade );
  theFTFModel0->SetMinEnergy( theFTFMin0 );
  theFTFModel0->SetMaxEnergy( theFTFMax );

  G4TheoFSGenerator * theFTFModel1 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel1->SetHighEnergyGenerator( theStringModel );
  theFTFModel1->SetTransport( theCascade );
  theFTFModel1->SetMinEnergy( theFTFMin1 );
  theFTFModel1->SetMaxEnergy( theFTFMax );

  G4CascadeInterface * theBERTModel0 = new G4CascadeInterface;
  theBERTModel0->SetMinEnergy( theBERTMin0 );
  theBERTModel0->SetMaxEnergy( theBERTMax );

  G4CascadeInterface * theBERTModel1 = new G4CascadeInterface;
  theBERTModel1->SetMinEnergy( theBERTMin1 );
  theBERTModel1->SetMaxEnergy( theBERTMax );

  G4VCrossSectionDataSet * thePiData = new G4CrossSectionPairGG( new G4PiNuclearCrossSection, 91*GeV );
  G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  //G4VCrossSectionDataSet * theGGNuclNuclData = G4CrossSectionDataSetRegistry::Instance()->
  //GetCrossSectionDataSet(G4GGNuclNuclCrossSection::Default_Name());
  
  auto theParticleIterator=GetParticleIterator(); //added for ver.10.03
   theParticleIterator->reset();
  while ((*theParticleIterator)()) 
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "pi+") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  //Inelastic scattering
	  G4PionPlusInelasticProcess* theInelasticProcess = 
	    new G4PionPlusInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( thePiData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	} 

      else if (particleName == "pi-") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  //Inelastic scattering
	  G4PionMinusInelasticProcess* theInelasticProcess = 
	    new G4PionMinusInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( thePiData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
	  //Absorption
	  pmanager->AddRestProcess(new G4PiMinusAbsorptionBertini, ordDefault);
	}
      
      else if (particleName == "kaon+") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering	
	  G4KaonPlusInelasticProcess* theInelasticProcess = 
	    new G4KaonPlusInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      
      else if (particleName == "kaon0S") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering	 
	  G4KaonZeroSInelasticProcess* theInelasticProcess = 
	    new G4KaonZeroSInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
	}

      else if (particleName == "kaon0L") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  // Inelastic scattering
	  G4KaonZeroLInelasticProcess* theInelasticProcess = 
	    new G4KaonZeroLInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 ); 
	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
	}

      else if (particleName == "kaon-") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4KaonMinusInelasticProcess* theInelasticProcess = 
	    new G4KaonMinusInelasticProcess("inelastic");	
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	  pmanager->AddRestProcess(new G4KaonMinusAbsorptionBertini, ordDefault);
	}

      else if (particleName == "proton") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->
					GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name()));
          theElasticProcess->RegisterMe( elastic_chip );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  // Inelastic scattering
	  G4ProtonInelasticProcess* theInelasticProcess = 
	    new G4ProtonInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      else if (particleName == "anti_proton") 
	{
	  // Elastic scattering
          const G4double elastic_elimitAntiNuc = 100.0*CLHEP::MeV;
          G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
          elastic_anuc->SetMinEnergy( elastic_elimitAntiNuc );
          G4CrossSectionElastic* elastic_anucxs = new G4CrossSectionElastic( elastic_anuc->GetComponentCrossSection() );
          G4HadronElastic* elastic_lhep2 = new G4HadronElastic();
          elastic_lhep2->SetMaxEnergy( elastic_elimitAntiNuc );
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( elastic_anucxs );
          theElasticProcess->RegisterMe( elastic_lhep2 );
          theElasticProcess->RegisterMe( elastic_anuc );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  // Inelastic scattering
	  G4AntiProtonInelasticProcess* theInelasticProcess = 
	    new G4AntiProtonInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( theAntiNucleonData );
	  theInelasticProcess->RegisterMe( theFTFModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	  // Absorption
	  pmanager->AddRestProcess(new G4AntiProtonAbsorptionFritiof, ordDefault);
	}

//      else if (particleName == "neutron") {
//	// elastic scattering
//	G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
//        theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name()));
//        G4HadronElastic* elastic_neutronChipsModel = new G4ChipsElasticModel();
//	elastic_neutronChipsModel->SetMinEnergy( 19.0*CLHEP::MeV );
//        theElasticProcess->RegisterMe( elastic_neutronChipsModel );
//	G4NeutronHPElastic * theElasticNeutronHP = new G4NeutronHPElastic;
//        theElasticNeutronHP->SetMinEnergy( theHPMin );
//        theElasticNeutronHP->SetMaxEnergy( theHPMax );
//	theElasticProcess->RegisterMe( theElasticNeutronHP );
//	theElasticProcess->AddDataSet( new G4NeutronHPElasticData );
//	pmanager->AddDiscreteProcess( theElasticProcess );
//	// inelastic scattering		
//	G4NeutronInelasticProcess* theInelasticProcess =
//	  new G4NeutronInelasticProcess("inelastic");
//	theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Neutron::Neutron() ) );
//	theInelasticProcess->RegisterMe( theFTFModel1 );
//        theInelasticProcess->RegisterMe( theBERTModel1 );
//	G4NeutronHPInelastic * theNeutronInelasticHPModel = new G4NeutronHPInelastic;
//        theNeutronInelasticHPModel->SetMinEnergy( theHPMin );
//        theNeutronInelasticHPModel->SetMaxEnergy( theHPMax );
//	theInelasticProcess->RegisterMe( theNeutronInelasticHPModel );
//	theInelasticProcess->AddDataSet( new G4NeutronHPInelasticData );
//	pmanager->AddDiscreteProcess(theInelasticProcess);
//	// capture
//	G4HadronCaptureProcess* theCaptureProcess =
//	  new G4HadronCaptureProcess;
//	G4NeutronHPCapture * theLENeutronCaptureModel = new G4NeutronHPCapture;
//	theLENeutronCaptureModel->SetMinEnergy(theHPMin);
//	theLENeutronCaptureModel->SetMaxEnergy(theHPMax);
//	theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
//	theCaptureProcess->AddDataSet( new G4NeutronHPCaptureData);
//	pmanager->AddDiscreteProcess(theCaptureProcess);
//
//      }
/////      else if (particleName == "anti_neutron") 
/////	{
/////	  // Elastic scattering
/////          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
/////          theElasticProcess->RegisterMe( elastic_lhep0 );
/////	  pmanager->AddDiscreteProcess( theElasticProcess );
/////          // Inelastic scattering (include annihilation on-fly)
/////	  G4AntiNeutronInelasticProcess* theInelasticProcess = 
/////	    new G4AntiNeutronInelasticProcess("inelastic");
/////	  theInelasticProcess->AddDataSet( theAntiNucleonData );
/////	  theInelasticProcess->RegisterMe( theFTFModel0 );
/////	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
/////	}
/////
/////      else if (particleName == "deuteron") 
/////	{
/////	  // Elastic scattering
/////          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
/////          theElasticProcess->RegisterMe( elastic_lhep0 );
/////	  pmanager->AddDiscreteProcess( theElasticProcess );
/////          // Inelastic scattering
/////	  G4DeuteronInelasticProcess* theInelasticProcess = 
/////	    new G4DeuteronInelasticProcess("inelastic");
/////	  theInelasticProcess->AddDataSet( theGGNuclNuclData );
/////	  theInelasticProcess->RegisterMe( theFTFModel1 );
/////            theInelasticProcess->RegisterMe( theBERTModel0 );
/////	  pmanager->AddDiscreteProcess( theInelasticProcess );
/////	}
/////      
/////      else if (particleName == "triton") 
/////	{
/////	  // Elastic scattering
/////          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
/////          theElasticProcess->RegisterMe( elastic_lhep0 );
/////	  pmanager->AddDiscreteProcess( theElasticProcess );
/////          // Inelastic scattering
/////	  G4TritonInelasticProcess* theInelasticProcess = 
/////	    new G4TritonInelasticProcess("inelastic");
/////	  theInelasticProcess->AddDataSet( theGGNuclNuclData );
/////	  theInelasticProcess->RegisterMe( theFTFModel1 );
/////          theInelasticProcess->RegisterMe( theBERTModel0 );
/////	  pmanager->AddDiscreteProcess( theInelasticProcess );
/////	}
/////      else if (particleName == "alpha") 
/////	{
/////	  // Elastic scattering
/////          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
/////          theElasticProcess->RegisterMe( elastic_lhep0 );
/////	  pmanager->AddDiscreteProcess( theElasticProcess );
/////          // Inelastic scattering
/////	  G4AlphaInelasticProcess* theInelasticProcess = 
/////	    new G4AlphaInelasticProcess("inelastic");	 
/////          theInelasticProcess->AddDataSet( theGGNuclNuclData );
/////	  theInelasticProcess->RegisterMe( theFTFModel1 );
/////          theInelasticProcess->RegisterMe( theBERTModel0 );
/////	  pmanager->AddDiscreteProcess( theInelasticProcess );
/////	}
    }


}

///////Decay
#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"

void S2SPhysicsList::ConstructDecay()
{
  G4Decay *theDecayProcess = new G4Decay();
  auto theParticleIterator=GetParticleIterator(); //added for ver.10.03
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pManager = particle->GetProcessManager();
    if( theDecayProcess->IsApplicable(*particle) && !particle->IsShortLived()){
      pManager->AddProcess( theDecayProcess );
      pManager->SetProcessOrdering( theDecayProcess, idxPostStep );
      pManager->SetProcessOrdering( theDecayProcess, idxAtRest );
    }
  }
  G4int VerboseLevel = 1;
  // Declare radioactive decay to the GenericIon in the IonTable.
  const G4IonTable *theIonTable = 
    G4ParticleTable::GetParticleTable()->GetIonTable();
  G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();
  
  for (G4int i=0; i<theIonTable->Entries(); i++) 
    {
      G4String particleName = theIonTable->GetParticle(i)->GetParticleName();
      G4String particleType = theIonTable->GetParticle(i)->GetParticleType();
      
      if (particleName == "GenericIon") 
	{
	  G4ProcessManager* pmanager = 
	    theIonTable->GetParticle(i)->GetProcessManager();
	  pmanager->SetVerboseLevel(VerboseLevel);
	  pmanager ->AddProcess(theRadioactiveDecay);
	  pmanager ->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
	  pmanager ->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
	} 
    }
}
