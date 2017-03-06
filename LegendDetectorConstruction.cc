///////////////////////////////////////////////////////////////////////////
// This code implementation is the intellectual property of the
// ton-scale 0vbb in Germanium collaboration. It is based on Geant4, an
// intellectual property of the RD44 GEANT4 collaboration.
//
// *********************
//
// Neither the authors of this software system, nor their employing
// institutes, nor the agencies providing financial support for this
// work make any representation or warranty, express or implied,
// regrding tis sotware system or assue any lability for its use.
// By copying, distrbuting or modifying the Proam (or any work based
// on he Program) yu indicate your accptance of ts statement,
// and all its terms.
/// \file LegendDetectorConstruction.cc
/// \brief Implementation of the LegendDetectorConstruction class

#include "LegendDetectorConstruction.hh"
#include "LegendDetectorMessenger.hh"
//#include "LegendScintSD.hh"
//#include "LegendPMTSD.hh"

//#include "LegendWLSSlab.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"
//#include "G4LogicalBorderSurface.hh"

#include "math.h"

#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4Colour.hh"
#include "G4MaterialPropertyVector.hh"
//#include "G4Material.hh"
//
// ROOT includes
#include "TFile.h"
#include "TH1F.h"

#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
using namespace CLHEP;
const G4double LegendDetectorConstruction::LambdaE = twopi *1.973269602e-16 * m * GeV;

LegendDetectorConstruction::LegendDetectorConstruction()
: G4VUserDetectorConstruction()
{
  detectorMessenger = new LegendDetectorMessenger(this);
	detector_type = "GERDA";
	innerVessel_FillMaterial = "ArgonLiquid";//"NitrogenGas";
	checkOverlaps = true;//false;

    
  G4String pathFile = "External_data/tpbGhemann.root";
  TFile *tpbFile = new TFile(pathFile.data());
  if (!tpbFile ) 
    G4cout<<" LegendDetectorConstruction ERROR:: file " << pathFile << " not found " << G4endl;
  else
    G4cout<<" LegendDetectorConstruction INFO:: file " << pathFile << " opened " << G4endl;
  fTPBspec=NULL;
  tpbFile->GetObject("tpbGhemann",fTPBspec);
  if (!fTPBspec ) 
    G4cout<<" LegendDetectorConstruction ERROR:: not graph tpbBhemann in file " << pathFile <<G4endl;
  else 
    G4cout<<" LegendDetectorConstruction info tpbBhemann graph found " <<G4endl;

   G4cout<<" LegendDetectorConstruction constant LambdaE =" << LambdaE << G4endl;
  
  // create directory 
  fDir = LegendAnalysis::Instance()->topDir()->mkdir("detec");
  fDir->cd();
  G4cout<<" LegendDetectorAction working root directory  is  " << G4endl;  
  gDirectory->pwd();
  G4cout << " ... " << G4endl;
  G4double LowE =1.4*eV;//885.6013 2.4796*eV;//500 nm
  G4double HighE = 12.3984*eV;//100 nm
  hDetecWLSPhotonE = new TH1F("DetecWLSPhotonE"," photon energy in WLS",1000,LowE,HighE);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LegendDetectorConstruction::~LegendDetectorConstruction()
{
// delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* LegendDetectorConstruction::Construct()
{

	#include "LegendDetectorMaterials.icc"
  ArgonOpticalProperties();
//	#include "LegendScintillationDefinitions.icc"
	////////////////////////////////////////////////////////////////////////////////////////
	//
  // World
  //
  solid_World = new G4Box("sol_World",50*m,50*m,30*m);
  logical_World = new G4LogicalVolume(solid_World,mat_air,"log_World");
	logical_World->SetVisAttributes (G4VisAttributes::Invisible);
  physical_World = new G4PVPlacement(0,G4ThreeVector(),logical_World,"phy_World",0,false,0,checkOverlaps);



	solid_Rock = new G4Box("sol_Rock",50*m,50*m,30*m);
	solid_Lab = new G4Box("sol_Lab",35*m,10*m,4*m);
	G4SubtractionSolid *solid_Rock2 = new G4SubtractionSolid("sol_Rock2", solid_Rock, solid_Lab ,0 , G4ThreeVector(-25*m,0,10.5*m));
	G4Tubs* solid_CutOut = new G4Tubs("sol_CutOut",0, 6.50001*m ,6.50001*m, 0, 2*M_PI);
	G4SubtractionSolid *solid_Rock3 = new G4SubtractionSolid("sol_Rock2", solid_Rock2, solid_CutOut ,0 , G4ThreeVector(0,0,0));


  logical_Rock = new G4LogicalVolume(solid_Rock3,mat_Rock,"log_Rock");
	logical_Rock->SetVisAttributes ( new G4VisAttributes(G4Colour(0.1,0.1,0.7,0.5) ) );//(0.7, 0.7, 0.7, 0.5) )); //grey 50% transparent
  physical_Rock = new G4PVPlacement(0,G4ThreeVector(),logical_Rock,"phy_Rock",logical_World,false,0,checkOverlaps);


  //TODO do not forget about adding back in the Ge Array

  G4double innerR_cryo = 0.95*m;
  G4double delta = 0.00001*m;

	//inner Vessel
	G4double innerVessel_Z[6] = {1*m,innerR_cryo+delta,innerR_cryo,-innerR_cryo,-innerR_cryo-delta,-1*m};
	G4double innerVessel_RMin[6] = {0*m,0*m,innerR_cryo,innerR_cryo,0*m,0*m};
	G4double innerVessel_RMax[6] = {1*m,1*m,1*m,1*m,1*m,1*m};

	solid_innerVessel = new G4Polycone("sol_innerVessel", 0, 2*M_PI,6,innerVessel_Z,innerVessel_RMin,innerVessel_RMax);
	logical_innerVessel = new G4LogicalVolume(solid_innerVessel, mat_Cu, "log_innerVessel" );
	logical_innerVessel->SetVisAttributes ( new G4VisAttributes(G4Colour(0.62, 0.3, 0.2,0.7) ));
	physical_innerVessel = new G4PVPlacement(0,G4ThreeVector(0,0,0),logical_innerVessel,"phy_innerVessel",logical_World,false,0,checkOverlaps);
  
  //Inner Vessel Optical properties
  //
  static const G4int NUMENTRIES_LAr = 69;
  G4double refl = .96;
  G4double effncy = 0.;
  G4double ee;
  G4double LAr_Energy[NUMENTRIES_LAr];
  G4double LArHighE = LambdaE / (115*nanometer);
  G4double LArLowE = LambdaE / (650*nanometer); 
  G4double de = ((LArHighE - LArLowE) / ((G4double)(NUMENTRIES_LAr-1)));
  G4double reflectivity_skin[NUMENTRIES_LAr];// = {refl,refl,refl,refl,refl,refl,refl};
  G4double efficiency_skin[NUMENTRIES_LAr];// = {effncy, effncy, effncy, effncy, effncy, effncy, effncy};
  
  for (int ji = 0; ji < NUMENTRIES_LAr; ji++)
  {
      ee=LArLowE+ ((G4double)ji) * de;
      LAr_Energy[ji] = ee;
      reflectivity_skin[ji] = refl;
      efficiency_skin[ji] = effncy;
      //G4cout<<LAr_Energy[ji]<<G4endl;
  }
  G4MaterialPropertiesTable* cryoHousing = new G4MaterialPropertiesTable();
  
  cryoHousing->AddProperty("REFLECTIVITY", LAr_Energy, reflectivity_skin,NUMENTRIES_LAr);
  cryoHousing->AddProperty("EFFICIENCY", LAr_Energy, efficiency_skin, NUMENTRIES_LAr);
  
  G4OpticalSurface* cryoHousingSurface = new G4OpticalSurface("HousingSurface",unified,polished,dielectric_metal);
  cryoHousingSurface->SetMaterialPropertiesTable(cryoHousing);

  skin_copper = new G4LogicalSkinSurface("CU_Cyro_Surf",logical_copperShield, cryoHousingSurface);
  //fill gas
  solid_fillGas = new G4Tubs("sol_fillGas", 0, innerR_cryo ,innerR_cryo, 0, 2*M_PI);
  logical_fillGas = new G4LogicalVolume(solid_fillGas, mat_fillGas, "log_fillGas" );
  logical_fillGas->SetVisAttributes ( new G4VisAttributes(G4Colour(0.7,0.1,0.1,0.5) ) );//0.5, 0.5, 0.5, 0.5) )); //grey 50% transparent
  physical_fillGas = new G4PVPlacement(0,G4ThreeVector(0,0,0),logical_fillGas,"phy_fillGas",logical_World,false,0,checkOverlaps);

  #include "Detector_MJDStyle.icc"

  //Pmts  
  G4double innerR_pmt = 0.*cm;
  G4double outerR_pmt = 2.3*cm;
  G4double height_pmt =0.03175*m;
  G4double startAngle_pmt = 0.*deg;
  G4double spanningAngle_pmt = 360.*deg;

  solid_Pmt = new G4Tubs("solid_pmt",innerR_pmt,outerR_pmt,height_pmt,startAngle_pmt,spanningAngle_pmt);
  solid_Photocath = new G4Tubs("solid_photocath",innerR_pmt,outerR_pmt,height_pmt/2,startAngle_pmt,spanningAngle_pmt);
 
  logical_Pmt = new G4LogicalVolume(solid_Pmt,G4Material::GetMaterial("Glass"),"logical_Pmt");
  logical_Photocath = new G4LogicalVolume(solid_Photocath,G4Material::GetMaterial("Al"),"logical_Photocath_log");
 
  physical_Photocath = new G4PVPlacement(0,G4ThreeVector(0,0,-height_pmt/2),logical_Photocath,"photocath",logical_Pmt,false,0,checkOverlaps);
  physical_PMT = new G4PVPlacement(0,G4ThreeVector(0,0,-innerR_cryo+height_pmt),logical_Pmt,"phys_Pmt",logical_fillGas,false,0,checkOverlaps);

  logical_Pmt->SetVisAttributes ( new G4VisAttributes(G4Colour(0.6,0.1,0.7) ) );
 
  //fPMTGlassOptSurface defined in LegendDetectorMaterials.icc
  new G4LogicalSkinSurface("PMTGlass_surf",logical_Pmt,fPMTGlassOptSurface);
 
  //Photocathode surface
  G4double photocath_EFF[NUMENTRIES_LAr];//={1.,1.,1.,1.,1.,1.,1.}; //Enables 'detection' of photons
  G4double photocath_ReR[NUMENTRIES_LAr];//={1.92,1.92,1.92,1.92,1.92,1.92,1.92};
  G4double photocath_ImR[NUMENTRIES_LAr];//={1.69,1.69,1.69,1.69,1.69,1.69,1.69};
   for (int ij = 0; ij < NUMENTRIES_LAr; ij++){
     photocath_EFF[ij] = 1.;//Enables 'detection' of photons
     photocath_ReR[ij] = 1.92;
     photocath_ImR[ij] = 1.69;
   }
  G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
  
  photocath_mt->AddProperty("EFFICIENCY",LAr_Energy,photocath_EFF,NUMENTRIES_LAr);
  photocath_mt->AddProperty("REALRINDEX",LAr_Energy,photocath_ReR,NUMENTRIES_LAr);
  photocath_mt->AddProperty("IMAGINARYRINDEX",LAr_Energy,photocath_ImR,NUMENTRIES_LAr);
  
  G4OpticalSurface* photocath_opsurf = new G4OpticalSurface("photocath_opsurf",glisur,polished,dielectric_metal);
  
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);
  skin_photocath = new G4LogicalSkinSurface("photocath_surf",logical_Photocath,photocath_opsurf);
  // add WLS 
  // -- WLS: TPB (Tetraphenyl butadiene)
  // --M.Gold from Gehmann et al plot

   fTPB = G4Material::GetMaterial("TPB", false);
   if (fTPB == 0) {
     G4NistManager* nist = G4NistManager::Instance();
     G4Element* elementH = nist->FindOrBuildElement("H");
     G4Element* elementC = nist->FindOrBuildElement("C");
     fTPB= new G4Material("TPB", 1*g/cm3, 2, kStateSolid);
     fTPB->AddElement(elementH, 22);
     fTPB->AddElement(elementC, 28);
   }

   // Now attach the optical properties to it.
   // Build table with photon energies
   
   const G4int numTPB =500;// 63;;
   G4double PPSCHighETPB = LambdaE /(350*nanometer);
   G4double PPSCLowETPB = LambdaE /(650*nanometer);//(650*nanometer); //598
   G4double deeTPB = ((PPSCHighETPB - PPSCLowETPB) / ((G4double)(numTPB-1)));
   G4double LAr_SCPPTPB[numTPB];
   for (G4int ji = 0; ji < numTPB; ji++)  LAr_SCPPTPB[ji]=PPSCLowETPB+ ((G4double)ji) * deeTPB;
   //outFile->cd();
 //G4cout<< PPSCHighETPB <<" = PPSCHighETPB"<<PPSCLowETPB<<" PPSCLowETPB"<<G4endl;
   G4double WLS_absorption[numTPB];
   G4double WLS_emission[numTPB];
   G4double Refraction[numTPB];
   tpbTable = new G4MaterialPropertiesTable();
  //cheesy way to include tables from https://arxiv.org/abs/1104.3259
  //#include "fLAr_SCPPTPB.icc"
  //#include "fTPBspec.icc"
  for (G4int ji=0;ji < numTPB; ji++) {
    //convert to nm to eV
    //LAr_SCPPTPB[ji] = LambdaE/LAr_SCPPTPB[ji];
    Refraction[ji] = 1.6; //this is just a guess
     // Should the TPB shift the Cherenkov light?
     // This makes a tail starting at 128 until the visible.
     if (LAr_SCPPTPB[ji] > 3.31*eV)
       // For the moment set it to always absorb photons
       WLS_absorption[ji] = 0.001*nm; //absorbs UV (always)
     else
       WLS_absorption[ji] = 10.5*m; //otherwise transparent
       WLS_emission[ji] = TPBEmissionSpectrum(LAr_SCPPTPB[ji]);//fTPBspec[ji];
       hDetecWLSPhotonE->SetBinContent(ji,WLS_emission[ji]);
       //G4cout<<" WLS emission "<<LAr_SCPPTPB[ji]<<", "<<WLS_emission[ji]<<G4endl;
   }
   tpbTable->AddProperty("RINDEX",LAr_SCPPTPB,Refraction,numTPB);
   tpbTable->AddProperty("WLSABSLENGTH",LAr_SCPPTPB,WLS_absorption,numTPB);
   tpbTable->AddProperty("WLSCOMPONENT",LAr_SCPPTPB,WLS_emission,numTPB);
   // From WArP
   tpbTable->AddConstProperty("WLSTIMECONSTANT", 0.01*ns);
   G4double WLSyield = 1.2;
   tpbTable->AddConstProperty("WLSMEANNUMBERPHOTONS",WLSyield);
   fTPB->SetMaterialPropertiesTable(tpbTable);


  //********************* Wave Length Shifters (WLS) ******/
  //In Lxe example, they create a class WLS. 
  //This seems overly eloquent, thus contruct 
  //in the detector construction seems fine
  G4double height_WLS =innerR_cryo;//0.05*m;
  G4double startAngle_WLS = 0.*deg;
  G4double spanningAngle_WLS = 360.*deg;
  G4double innerR_WLS = innerR_cryo/2;//0.*m;
  G4double outerR_WLS = innerR_cryo/2+5.*cm;//2 +delta*100;
  G4Tubs* solid_ScintSlab = new G4Tubs("Slab",innerR_WLS,outerR_WLS,height_WLS,startAngle_WLS,spanningAngle_WLS);
  logical_ScintSlab = new G4LogicalVolume(solid_ScintSlab,G4Material::GetMaterial("Polystyrene"),"ScintSlab",0,0,0);
  physical_ScintSlab = new G4PVPlacement(0,G4ThreeVector(0,0,0/*innerR_cryo-height_WLS/2*/),
                                         logical_ScintSlab,"phy_ScintSlab",logical_fillGas,false,0,checkOverlaps);
   /**** WLS skin ****/

   // Define a rough optical surface to be used in the interface between WLS and LAr
   // 50% roughness in the surface
   // This surface will be attached between the WLS and the LAr in all instances
   G4double roughness = 0.5;
   fWLSoptSurf = new G4OpticalSurface("WLS_rough_surf",glisur,ground,dielectric_dielectric,roughness);
   fWLSoptSurf->SetMaterialPropertiesTable(tpbTable);
   fSkin_WLS = new G4LogicalSkinSurface("WLS_Surf",logical_ScintSlab,fWLSoptSurf);

  return physical_World;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//  ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//   ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LegendDetectorConstruction::ConstructSDandField()
{
  //This top method is taken from the LXe example
  if (!Scint_SD.Get()) {
		G4cout<<"Construction /Legend/scintSD"<<G4endl;
		LegendScintSD* scint_SD = new LegendScintSD("/LegendDet/scintSD");//intialize with some G4string, to call later
		Scint_SD.Put(scint_SD);
	}
	G4SDManager::GetSDMpointer()->AddNewDetector(Scint_SD.Get());
  SetSensitiveDetector(logical_fillGas,Scint_SD.Get());
  
  // PMT SD
  if (!Pmt_SD.Get()) {
    //Created here so it exists as pmts are being placed
    G4cout << "Construction /Legend/pmtSD" << G4endl;
    LegendPMTSD* pmt_SD = new LegendPMTSD("/LegendDet/pmtSD");
    Pmt_SD.Put(pmt_SD);

    pmt_SD->InitPMTs(1); //let pmtSD know # of pmts
    pmt_SD->SetPmtPosition(G4ThreeVector(0,0,-0.95*m+0.03175*m));
  }
  G4SDManager::GetSDMpointer()->AddNewDetector(Pmt_SD.Get());

}
/// methods imported from the mpiklarge class
/// optical properties of lar in several places
void LegendDetectorConstruction::ArgonOpticalProperties()
{
	//Taken from home/gold/MaGe-master/munichteststand/src/GEMPIKLArGe.cc
  static const G4int NUMENTRIES = 69;
  const G4int num = 69;
  static const G4double temp = 88.5*kelvin;
 // static const G4double LambdaE = twopi *1.973269602e-16 * m * GeV;
  G4double scint_yield = 28120./MeV;//40000./MeV; //sono 40000

  G4int ji;
  G4double e;
  G4double ee;

  G4double PPCKOVHighE = LambdaE / (115*nanometer);
  G4double PPCKOVLowE = LambdaE / (650*nanometer); 
  G4double de = ((PPCKOVHighE - PPCKOVLowE) / ((G4double)(NUMENTRIES-1)));
  G4double LArAbsLength = 3*m; //just a number. ICARUS says it is negligible over
  // liquid argon (LAr)  
  G4double LAr_PPCK[(NUMENTRIES)];
  G4double LAr_RIND[(NUMENTRIES)];
  G4double LAr_RAYL[(NUMENTRIES)];
  G4double LAr_ABSL[(NUMENTRIES)];
  
  G4double lar_absl_xuv = 60*cm;
	G4double lar_absl_vis = 1000*m;
 
  for (ji = 0; ji < NUMENTRIES; ji++)
    {
      e = PPCKOVLowE + ((G4double)ji) * de;
      LAr_PPCK[ji] = e;
      LAr_RIND[ji] = LArRefIndex((LambdaE / e));
      LAr_RAYL[ji] = LArRayLength((LambdaE / e), temp);
      //LAr_ABSL[ji] = LArAbsLength;
      
      if (((LambdaE / e)/nm) < 200.0) {
	    	  LAr_ABSL[ji] = lar_absl_xuv;
	    } else {
	    	 LAr_ABSL[ji] = lar_absl_vis;
	    }
    }

  G4double PPSCHighE = LambdaE /(115*nanometer);
  G4double PPSCLowE = LambdaE /(136*nanometer);
  G4double dee = ((PPSCHighE - PPSCLowE) / ((G4double)(num-1)));
  G4double LAr_SCIN[num];
  G4double LAr_SCPP[num];
  for (ji = 0; ji < num; ji++)
    {
      ee=PPSCLowE+ ((G4double)ji) * dee;
      LAr_SCPP[ji]=ee;
      LAr_SCIN[ji]=ArScintillationSpectrum((LambdaE/ee)/nanometer);
      //G4cout<<LAr_SCPP[ji]<<G4endl;
    }

  G4MaterialPropertiesTable* LAr_mt = new G4MaterialPropertiesTable();

  LAr_mt->AddProperty("RINDEX",        LAr_PPCK, LAr_RIND, NUMENTRIES);
  LAr_mt->AddProperty("RAYLEIGH",      LAr_PPCK, LAr_RAYL, NUMENTRIES);
  LAr_mt->AddProperty("ABSLENGTH",     LAr_PPCK, LAr_ABSL, NUMENTRIES);
  if ( (LAr_SCPP[0] >= PPCKOVLowE) &&
       (LAr_SCPP[(sizeof(LAr_SCPP)/sizeof(G4double) - 1)] <= PPCKOVHighE) )
    {
      LAr_mt->AddProperty("FASTCOMPONENT",LAr_SCPP,LAr_SCIN,num);
      LAr_mt->AddProperty("SLOWCOMPONENT",LAr_SCPP,LAr_SCIN,num);
    }
  LAr_mt->AddConstProperty("SCINTILLATIONYIELD",scint_yield);
  LAr_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
  LAr_mt->AddConstProperty("FASTTIMECONSTANT", 5.95*ns);//6.*ns);
  LAr_mt->AddConstProperty("SLOWTIMECONSTANT",922*ns);//1590.*ns);
  LAr_mt->AddConstProperty("YIELDRATIO",0.23);
// G4cout<<LAr_SCIN<<G4endl;
  G4double fano = 0.11;
  LAr_mt->AddConstProperty("RESOLUTIONSCALE",fano); 
  mat_ArLiq->SetMaterialPropertiesTable(LAr_mt); // G4Material defined in Detector_Materials.icc
  mat_ArLiq->GetIonisation()->SetBirksConstant(5.1748e-4*cm/MeV);//0.0144*mm/MeV);
 
}

G4double LegendDetectorConstruction::LArEpsilon(const G4double lambda)
{
  G4double epsilon;
  if (lambda < 110*nanometer) return 1.0e4; // lambda MUST be > 110.0 nm
  epsilon = lambda / micrometer; // switch to micrometers
  epsilon = 1.0 / (epsilon * epsilon); // 1 / (lambda)^2
  epsilon = 1.2055e-2 * ( 0.2075 / (91.012 - epsilon) +
                          0.0415 / (87.892 - epsilon) +
                          4.3330 / (214.02 - epsilon) );
  epsilon *= (8./12.); // Bideau-Sellmeier -> Clausius-Mossotti
  G4double LArRho = 1.396*g/cm3;
  G4double GArRho = 1.66e-03*g/cm3;
  epsilon *= (LArRho / GArRho); // density correction (Ar gas -> LAr liquid)
  if (epsilon < 0.0 || epsilon > 0.999999) return 4.0e6;
  epsilon = (1.0 + 2.0 * epsilon) / (1.0 - epsilon); // solve Clausius-Mossotti
  return epsilon;
}

G4double LegendDetectorConstruction::LArRefIndex(const G4double lambda)
{
  return ( sqrt(LArEpsilon(lambda)) ); // square root of dielectric constant
}
G4double LegendDetectorConstruction::LArRayLength(const G4double lambda,const
				   G4double temp)
{
  G4double dyne = 1.0e-5*newton;
  static const G4double LArKT = 2.18e-10 * cm2/dyne; // LAr isothermal compressibility
  static const G4double k = 1.380658e-23 * joule/kelvin; // the Boltzmann constant
  G4double h;
  h = LArEpsilon(lambda);
  if (h < 1.00000001) h = 1.00000001; // just a precaution
  h = (h - 1.0) * (h + 2.0); // the "dielectric constant" dependance
  h *= h; // take the square
  h *= LArKT * temp * k; // compressibility * temp * Boltzmann constant
  h /= lambda * lambda * lambda * lambda; // (lambda)^4
  h *= 9.18704494231105429; // (2 * Pi / 3)^3
  if ( h < (1.0 / (10.0 * km)) ) h = 1.0 / (10.0 * km); // just a precaution
  if ( h > (1.0 / (0.1 * nanometer)) ) h = 1.0 / (0.1 * nanometer); // just a precaution
  return ( 1.0 / h );
}
G4double LegendDetectorConstruction::ArScintillationSpectrum(const G4double kk)
{
  G4double waveL;
  waveL =exp(-0.5*((kk-128.0)/(2.929))*((kk-128.0)/(2.929)));
  return waveL;
}
#include "G4RunManager.hh"

void LegendDetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LegendDetectorConstruction::SetOverlapsCheck(G4bool f_check)
{
	checkOverlaps = f_check;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LegendDetectorConstruction::SetFillGas(G4String f_gas)
{
	innerVessel_FillMaterial = f_gas;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LegendDetectorConstruction::SetShieldStyle(G4String f_type)
{
	detector_type = f_type;
}
