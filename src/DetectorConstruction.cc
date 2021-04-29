//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include <fstream>

using namespace std;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
	: G4VUserDetectorConstruction()
{	// World
	world_half_side = 100*mm;

	// FATGEM
	fatgem_radius   = 20.*mm;
	fatgem_thickness= 2.3*mm; // half thickness has to be given!!!
	hole_radius		= 1.*mm;
	hole_pitch		= 5.*mm;
	hole_voff		= sqrt((hole_pitch*hole_pitch) - (hole_pitch/2*hole_pitch/2));

	// Photomultipliers
	PMT_radius     = 12.5*mm;
	PMT_thickness  = 2.5*mm;  // half thickness has to be given!!!
	phc_thickness  = 0.05*mm; // half thickness has to be given!!!

	// source
	source_side_x = 12.*mm;
	source_side_y = 6.*mm;
	source_side_z = 0.5*mm;
	kapton_thickness = 0.035*mm; // half thickness has to be given!!!
	al_thickness = 0.025*mm; // half thickness has to be given!!!

	// Drifts
	drift_1    		= 15.*mm;
	drift_2			= 20.5*mm;

	// Teflon base
	base_radius    = 50.*mm;

	// booleans
	disable_FATGEM_base = false;
	disable_PMT_base   = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() { ; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	// Option to switch on/off checking of volumes overlaps
	G4bool checkOverlaps = true;
	G4NistManager* nist = G4NistManager::Instance();
	
	// ------------ Generate & Add Material Properties Table ------------

	// Xe ------------------------------------------
	// gas density
	G4double pressure = 9.*bar;
	G4double temperature = 298.15 * kelvin;
	G4Material* matXe = nist->ConstructNewGasMaterial("GXe", "G4_Xe", temperature, pressure);
	G4double density_xe = (pressure / atmosphere) * 131.29 / (temperature / kelvin * 82.058); // g/cm^3

	// refractive index calculations
	const G4int nXe_entries = 100;
	G4double nXe_energy[nXe_entries];
	G4double refractiveIndex_Xe[nXe_entries];
	G4double energy_start = 1.8 * eV;
	G4double energy_finish = 8.9 * eV;

	for (G4int i = 0; i < nXe_entries; i++) {
		nXe_energy[i] = (energy_start + i * (energy_finish - energy_start) / nXe_entries);


		// Formula for the refractive index taken from
		// A. Baldini et al., "Liquid Xe scintillation calorimetry 
		// and Xe optical properties", arXiv:physics/0401072v1 [physics.ins-det]

		// The Lorentz-Lorenz equation (also known as Clausius-Mossotti equation)
		// relates the refractive index of a fluid with its density:
		// (n^2 - 1) / (n^2 + 2) = - A · d_M,     (1)
		// where n is the refractive index, d_M is the molar density and
		// A is the first refractivity viral coefficient:
		// A(E) = \sum_i^3 P_i / (E^2 - E_i^2),   (2)
		// with:
		G4double P[3] = { 71.23, 77.75, 1384.89 }; // [eV^3 cm3 / mole]
		G4double E[3] = { 8.4, 8.81, 13.2 };       // [eV]

		// Note.- Equation (1) has, actually, a sign difference with respect 
		// to the one appearing in the reference. Otherwise, it yields values
		// for the refractive index below 1.

		// Let's calculate the virial coefficient.
		// We won't use the implicit system of units of Geant4 because
		// it results in loss of numerical precision.

		G4double energy_ite = nXe_energy[i] / eV;

		G4double virial = 0.;
		for (G4int j = 0; j < 3; j++)
			virial = virial + P[j] / (energy_ite * energy_ite - E[j] * E[j]);

		G4double mol_density = density_xe / 131.29;
		G4double alpha = virial * mol_density;

		// Isolating now the n2 from equation (1) and taking the square root
		refractiveIndex_Xe[i] = (1. - 2 * alpha) / (1. + alpha);

		if (refractiveIndex_Xe[i] < 1.) {
			// "Non-physical refractive index for energy "
			refractiveIndex_Xe[i] = 1.;
		}

	}
	assert(sizeof(refractiveIndex_Xe) == sizeof(nXe_energy));

	G4MaterialPropertiesTable* Xe_MPT1 = new G4MaterialPropertiesTable();

	Xe_MPT1->AddProperty("RINDEX", nXe_energy, refractiveIndex_Xe, nXe_entries)
		->SetSpline(true);

	G4cout << "Xe G4MaterialPropertiesTable" << G4endl;

	matXe->SetMaterialPropertiesTable(Xe_MPT1);


	// PMT window ----------------------------------------
	G4Material* PMT_window_mat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

	// define window refractive index (source: janis.com/Libraries/Window_Transmissions/FusedSilicaUVGrade_SiO2_TransmissionCurveDataSheet.sflb.ashx)
	G4double photonEnergy[] = { 7.3074*eV, 6.7149*eV, 6.2113*eV, 5.7941*eV, 4.4319*eV, 4.1121*eV, 3.4035*eV,
								3.0704*eV, 2.8505*eV, 2.2748*eV, 2.1141*eV, 2.108*eV, 1.9296*eV, 1.8928*eV, 1.441*eV };
	
	G4double refractiveIndex1[] = { 1.6150, 1.5750, 1.5500, 1.5337, 1.4940, 1.4872, 1.4745,
									1.4696, 1.4666, 1.4601, 1.4585, 1.4584, 1.4567, 1.4564, 1.4525 };

	const G4int nEntries = sizeof(photonEnergy) / sizeof(G4double);

	G4MaterialPropertiesTable* pmtMPT1 = new G4MaterialPropertiesTable();

	pmtMPT1->AddProperty("RINDEX", photonEnergy, refractiveIndex1, nEntries);
	PMT_window_mat->SetMaterialPropertiesTable(pmtMPT1);


	// PMT photocathode ----------------------------------
	G4Material* PMT_phc_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");


	// Plexiglass (for FATGEM and source) ----------------
	G4Material* PMMA_mat = nist->FindOrBuildMaterial("G4_PLEXIGLASS");


	// Teflon bases   ------------------------------------
	G4Material* teflon_mat = nist->FindOrBuildMaterial("G4_TEFLON");

	const G4int num3 = 4;
	// source: C. Silva et al. “Reflectance of Polytetrafluoroethylene (PTFE) for Xenon Scintillation Light”,
	// J. Appl. Phys. 107 (2010) 064902
	G4double ephoton_teflon[num3] = { 2.21*eV , 3.95*eV, 4.87*eV, 7.3*eV };
	G4double reflectivity_teflon[num3] = { 0.98, 0.93, 0.85, 0.61};
	//G4double reflectivity_teflon[num3] = { 1., 1., 1., 1.};

	// Kapton layer at source ---------------------------
	G4Material* kapton_mat = nist->FindOrBuildMaterial("G4_KAPTON");

	// Aluminum layer at source ---------------------------
	G4Material* al_mat = nist->FindOrBuildMaterial("G4_Al");

	//PEN poly(ethylene 2,6-naphthalate)
	G4int ncomponents, z;
	G4double a, density=1.36*g/cm3;
	G4String name, symbol;
	G4Element* H  = new G4Element(name="Hydrogen",  symbol="H" , z=1,  a=1.00794*g/mole );  
	G4Element* C  = new G4Element(name="Carbon",    symbol="C",  z=6,  a=12.0107*g/mole );
	G4Element* O  = new G4Element(name="Oxygen",    symbol="O",  z=8,  a=15.9994*g/mole );
	G4Material* fPEN = new G4Material (name="PEN",density,ncomponents=3,kStateSolid);
	fPEN->AddElement(C,14);
	fPEN->AddElement(H,10);
	fPEN->AddElement(O,4);
	G4Material* fPENfoil = new  G4Material(name="PEN",density, fPEN);

	G4MaterialPropertiesTable *myPEN = new G4MaterialPropertiesTable();
	G4MaterialPropertiesTable *myPENfoil = new G4MaterialPropertiesTable();
	const G4double hc = h_Planck*c_light/eV/nm;

	G4double PEN_ENE[1000], PEN_EMISSION_VAL[1000], PEN_ABSORPTION_VAL[1000];
	G4int dim = 0 ;
	G4double myene, myene2, myvalue, PENNorma = 0 ;
	ifstream fpen_emission("pen_emission_spectrum.dat");
	if ( !fpen_emission.is_open())
        	cerr << "ERROR: Could not open PEN emission file" <<endl;
	while(!fpen_emission.eof()) {
		fpen_emission >> myene >> myvalue ;
		if(fpen_emission.eof()) break;
		PEN_ENE[dim] = myene*eV;
		if(hc/ (PEN_ENE[dim]/nm) > 600) {
			PEN_EMISSION_VAL[dim] = 0 ;
          	} else  {
			PEN_EMISSION_VAL[dim] = myvalue;
			PENNorma += myvalue ;
		}
		dim++;
	}
	fpen_emission.close();
	myPEN->AddProperty("WLSCOMPONENT",PEN_ENE,  PEN_EMISSION_VAL,       dim);
	myPENfoil->AddProperty("WLSCOMPONENT",PEN_ENE,  PEN_EMISSION_VAL,       dim);

	G4double PEN_ENE2[1000], PEN_WLS_ABSORPTION_VAL[1000], PEN_RINDEX[1000], PEN_RAYL[1000], PEN_RAYL_FOIL[1000];
	dim = 0 ;
	ifstream fpen_absorption("pen_absorption_length.dat");
	if ( !fpen_absorption.is_open())
		cerr << "ERROR: Could not open PEN absorption file" <<endl;

	while(!fpen_absorption.eof()) {
          fpen_absorption >> myene2 >> myvalue ;
          if(fpen_absorption.eof()) break;
          PEN_ENE2[dim] = myene2*eV;
	  // PEN_WLS_ABSORPTION_VAL[dim] = myvalue;
          if (myene2 > 1240./250.)  { // transition point at 250 nm
                // in VUV normal absorption and scattering disabled (WLS absorption will dominate)
		  PEN_WLS_ABSORPTION_VAL[dim] = 100*nm;
                  PEN_ABSORPTION_VAL[dim]  = 1000.*m;
                  PEN_RINDEX[dim]          = 1.75;
                  PEN_RAYL[dim]            = 1000.*m;
		  PEN_RAYL_FOIL[dim]       = 1000.*m;
          }
          else {
                // in visible absorption and scattering meaningful
		  PEN_WLS_ABSORPTION_VAL[dim] = 1000*m;
                  PEN_ABSORPTION_VAL[dim]  = 1.*cm; // Efremenko et al.
                  PEN_RINDEX[dim]          = 1.75;
                  PEN_RAYL[dim]            = 10*cm; // educated guess
		  PEN_RAYL_FOIL[dim]       = 150.*micrometer; // educated guess
          }
          dim++;
	}
	fpen_absorption.close();

	myPEN->AddProperty("WLSABSLENGTH",PEN_ENE2, PEN_WLS_ABSORPTION_VAL, dim);
	myPEN->AddProperty("ABSLENGTH",   PEN_ENE2, PEN_ABSORPTION_VAL,     dim);
	myPEN->AddProperty("RINDEX",      PEN_ENE2, PEN_RINDEX,             dim);
	myPEN->AddProperty("RAYLEIGH",    PEN_ENE2, PEN_RAYL,               dim);
	myPEN->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.091);
	myPEN->AddConstProperty("WLSTIMECONSTANT",20.*ns);
	myPEN->AddConstProperty("WLSEFFICIENCY",0.091); // https://arxiv.org/abs/2103.03232
	fPEN->SetMaterialPropertiesTable(myPEN);

	myPENfoil->AddProperty("WLSABSLENGTH",PEN_ENE2, PEN_WLS_ABSORPTION_VAL, dim);
	myPENfoil->AddProperty("ABSLENGTH",   PEN_ENE2, PEN_ABSORPTION_VAL,     dim);
	myPENfoil->AddProperty("RINDEX",      PEN_ENE2, PEN_RINDEX,             dim);
	myPENfoil->AddProperty("RAYLEIGH",    PEN_ENE2, PEN_RAYL_FOIL,          dim);
	myPENfoil->AddConstProperty("WLSMEANNUMBERPHOTONS", 1.);
	myPENfoil->AddConstProperty("WLSTIMECONSTANT",20.*ns);
	myPENfoil->AddConstProperty("WLSEFFICIENCY",0.34); // https://arxiv.org/abs/2103.03232
	fPENfoil->SetMaterialPropertiesTable(myPENfoil);

  // ----------- colors ------------

	G4VisAttributes * blue = new G4VisAttributes(G4Colour(0. ,0.8 ,0.9, 0.6));
	blue -> SetVisibility(true);
	blue -> SetForceSolid(true);
	G4VisAttributes * black = new G4VisAttributes(G4Colour(0.2 ,0.2 ,0.2));
	black -> SetVisibility(true);
	black -> SetForceSolid(true);
	G4VisAttributes * red = new G4VisAttributes(G4Colour(0.3 ,0.1 ,0., 0.5));
	red -> SetVisibility(true);
	red -> SetForceSolid(true);
	G4VisAttributes * white = new G4VisAttributes(G4Colour(0.95, 0.95, 0.95, 0.5));
	white -> SetVisibility(true);
	white -> SetForceSolid(true);
	G4VisAttributes * yellow = new G4VisAttributes(G4Colour(1. ,1. ,0.));
	yellow -> SetVisibility(true);
	yellow -> SetForceSolid(true);
	G4VisAttributes * gray = new G4VisAttributes(G4Colour(0.7 ,0.7 ,0.7));
	gray -> SetVisibility(true);
	gray -> SetForceSolid(true);

  // ------------- Volumes --------------

    // The experimental Hall
	G4Box* expHall_box = new G4Box("World", world_half_side, world_half_side, world_half_side);
	G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box, matXe, "World", 0, 0, 0);
	G4VPhysicalVolume* expHall_phys = new G4PVPlacement(0, G4ThreeVector(), expHall_log, "World", 0, false, 0, checkOverlaps);


	// PMT
	// window
	G4Tubs* PMT_win = new G4Tubs("PMT_win", 0.*cm, PMT_radius, PMT_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* pmtWindow_log = new G4LogicalVolume(PMT_win, PMT_window_mat, "PMT_win", 0, 0, 0);
	G4VPhysicalVolume* pmtWindow_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -drift_2 -fatgem_thickness - PMT_thickness), pmtWindow_log, "PMT_win", expHall_log, false, 0, checkOverlaps);
	pmtWindow_log -> SetVisAttributes(blue);
	
	// photocathode
	G4Tubs* PMT_phc = new G4Tubs("PMT_phCa", 0.*cm, PMT_radius, phc_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* pmtPhc_log = new G4LogicalVolume(PMT_phc, PMT_phc_mat, "PMT_phCa", 0, 0, 0);
	G4VPhysicalVolume* pmtPhc_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -drift_2 -fatgem_thickness - 2*PMT_thickness - phc_thickness), pmtPhc_log, "PMT_phCa", expHall_log, false, 0, checkOverlaps);
	pmtPhc_log -> SetVisAttributes(black);


	// Source
	G4Box* source = new G4Box("source", source_side_x, source_side_y, source_side_z);
	G4LogicalVolume* source_log = new G4LogicalVolume(source, PMMA_mat, "source", 0, 0, 0);
	G4VPhysicalVolume* source_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, drift_1 +fatgem_thickness +2*al_thickness +2*kapton_thickness +source_side_z), source_log, "source", expHall_log, false, 0, checkOverlaps);
	source_log -> SetVisAttributes(blue);

	// kapton layer
	G4Tubs* kap_layer = new G4Tubs("kap_layer", 0.*cm, base_radius, kapton_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* kap_layer_log = new G4LogicalVolume(kap_layer, kapton_mat, "kap_layer", 0, 0, 0);
	G4VPhysicalVolume* kap_layer_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, drift_1 +fatgem_thickness +2*al_thickness +kapton_thickness), kap_layer_log, "kap_layer", expHall_log, false, 0, checkOverlaps);
	kap_layer_log -> SetVisAttributes(red);	

	// aluminum layer
	G4Tubs* al_layer = new G4Tubs("al_layer", 0.*cm, base_radius, al_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* al_layer_log = new G4LogicalVolume(al_layer, al_mat, "al_layer", 0, 0, 0);
	G4VPhysicalVolume* al_layer_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, drift_1 +fatgem_thickness +al_thickness), al_layer_log, "al_layer", expHall_log, false, 0, checkOverlaps);
	al_layer_log -> SetVisAttributes(gray);	


	// FATGEM
	// volume
	G4Tubs* fatgem = new G4Tubs("fatgem", 0.*cm, fatgem_radius, fatgem_thickness, 0. * deg, 360. * deg);
	//G4LogicalVolume* fatgem_log = new G4LogicalVolume(fatgem, PMMA_mat, "fatgem", 0, 0, 0);
	G4LogicalVolume* fatgem_log = new G4LogicalVolume(fatgem, fPEN, "fatgem", 0, 0, 0);
	G4VPhysicalVolume* fatgem_phys = new G4PVPlacement(0, G4ThreeVector(), fatgem_log, "fatgem", expHall_log, false, 0, checkOverlaps);
	fatgem_log -> SetVisAttributes(red);

	G4Tubs* hole = new G4Tubs("hole", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* hole_log = new G4LogicalVolume(hole, matXe, "hole", 0, 0, 0);

	G4double hole_x[31]={0,-hole_pitch,hole_pitch,-hole_pitch*2,hole_pitch*2,
			     hole_pitch/2.,-hole_pitch/2.,hole_pitch*1.5,-hole_pitch*1.5,hole_pitch*2.5,-hole_pitch*2.5,
			     hole_pitch/2.,-hole_pitch/2.,hole_pitch*1.5,-hole_pitch*1.5,hole_pitch*2.5,-hole_pitch*2.5,
			     0,-hole_pitch,hole_pitch,-hole_pitch*2,hole_pitch*2,
			     0,-hole_pitch,hole_pitch,-hole_pitch*2,hole_pitch*2,
			     hole_pitch/2.,-hole_pitch/2.,
			     hole_pitch/2.,-hole_pitch/2.};
	G4double hole_y[31]= {0,0,0,0,0,
			      -hole_voff,-hole_voff,-hole_voff,-hole_voff,-hole_voff,-hole_voff,
			      hole_voff,hole_voff,hole_voff,hole_voff,hole_voff,hole_voff,
			      -2*hole_voff,-2*hole_voff,-2*hole_voff,-2*hole_voff,-2*hole_voff,
			      2*hole_voff,2*hole_voff,2*hole_voff,2*hole_voff,2*hole_voff,
			       -3*hole_voff,-3*hole_voff,
			      3*hole_voff,3*hole_voff};

	G4VPhysicalVolume* hole_phys[31];
	for(int i=0; i<31; i++)
	  hole_phys[i] = new G4PVPlacement(0, G4ThreeVector(hole_x[i], hole_y[i], 0), hole_log, "hole", fatgem_log, false, i, checkOverlaps);
	
	// optical surface 
	G4OpticalSurface* OS_fatgem_hole = new G4OpticalSurface("fatgem_hole_surface");

	for(int i=0; i<31; i++) {
		new G4LogicalBorderSurface("fatgem_hole_surface", fatgem_phys, hole_phys[i], OS_fatgem_hole);
		new G4LogicalBorderSurface("fatgem_hole_surface", hole_phys[i], fatgem_phys, OS_fatgem_hole);
	}

	OS_fatgem_hole->SetType(dielectric_dielectric);
	OS_fatgem_hole->SetFinish(ground);
	OS_fatgem_hole->SetModel(glisur);
	OS_fatgem_hole->SetPolish(0.1);

	const G4int num = 11;
	G4double ephoton_fatgem[num] = { 0.12*eV, 0.31*eV, 0.62*eV, 1.24*eV, 1.77*eV, 2.067*eV, 2.48*eV, 3.1*eV, 4.13*eV, 6.2*eV,  8.9*eV };
	G4double reflectivity_fatgem[num] = { 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. };
	G4MaterialPropertiesTable* R_fatgem = new G4MaterialPropertiesTable();
	R_fatgem->AddProperty("REFLECTIVITY", ephoton_fatgem, reflectivity_fatgem, num);
	OS_fatgem_hole->SetMaterialPropertiesTable(R_fatgem);

	G4OpticalSurface* OS_clevios = new G4OpticalSurface("clevios_surface");
	new G4LogicalBorderSurface("clevios_surface", fatgem_phys, expHall_phys, OS_clevios);
	new G4LogicalBorderSurface("clevios_surface", expHall_phys, fatgem_phys, OS_clevios);

	OS_clevios->SetType(dielectric_dielectric);
        OS_clevios->SetFinish(polished);
        OS_clevios->SetModel(glisur);
        OS_clevios->SetPolish(1.);

	//const G4int num = 11;
        G4double ephoton_clevios[num] = { 0.12*eV, 0.31*eV, 0.62*eV, 1.24*eV, 1.77*eV, 2.067*eV, 2.48*eV, 3.1*eV, 4.13*eV, 6.2*eV,  8.9*eV };
	// 1% loss in VIS, 20% loss in VUV
        G4double reflectivity_clevios[num] = { 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.8, 0.8 };
	G4MaterialPropertiesTable* R_clevios = new G4MaterialPropertiesTable();
	R_clevios->AddProperty("REFLECTIVITY", ephoton_clevios, reflectivity_clevios, num);
	OS_clevios->SetMaterialPropertiesTable(R_clevios);

	// FATGEM base
	if (!disable_FATGEM_base) {

		// volume
		G4Tubs* FATGEM_base = new G4Tubs("FATGEM_base", fatgem_radius, base_radius, fatgem_thickness, 0. * deg, 360. * deg);
		G4LogicalVolume* FATGEM_base_log	= new G4LogicalVolume(FATGEM_base, teflon_mat, "FATGEM_base", 0, 0, 0);
		G4VPhysicalVolume* FATGEM_base_phys = new G4PVPlacement(0, G4ThreeVector(), FATGEM_base_log, "cathode", expHall_log, false, 0, checkOverlaps);
		FATGEM_base_log  -> SetVisAttributes(yellow);

		// optical surface
		G4OpticalSurface* op_FATGEM_base = new G4OpticalSurface("FATGEM_base_Surface");
		G4LogicalSkinSurface* FATGEM_base_Surface = new G4LogicalSkinSurface("FATGEM_base_Surface", FATGEM_base_log, op_FATGEM_base);

		G4MaterialPropertiesTable* FATGEM_base_ST2 = new G4MaterialPropertiesTable();
		FATGEM_base_ST2->AddProperty("REFLECTIVITY", ephoton_teflon, reflectivity_teflon, num3);

		op_FATGEM_base->SetMaterialPropertiesTable(FATGEM_base_ST2);
		op_FATGEM_base->SetType(dielectric_metal);
		op_FATGEM_base->SetFinish(ground);
		op_FATGEM_base->SetPolish(0.2);
		op_FATGEM_base->SetModel(glisur);
	}


	// PMT base
	if (!disable_PMT_base) {

		// volume
		G4Tubs* PMT_base = new G4Tubs("PMT_base", PMT_radius, base_radius, PMT_thickness, 0. * deg, 360. * deg);
		G4LogicalVolume* PMT_base_log	= new G4LogicalVolume(PMT_base, teflon_mat, "PMT_base", 0, 0, 0);
		G4VPhysicalVolume* PMT_base_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -drift_2 -fatgem_thickness - PMT_thickness), PMT_base_log, "PMT_base", expHall_log, false, 0, checkOverlaps);
		PMT_base_log  -> SetVisAttributes(yellow);

		// optical surface
		G4OpticalSurface* op_PMT_base = new G4OpticalSurface("PMT_base_Surface");
		G4LogicalSkinSurface* PMT_base_Surface = new G4LogicalSkinSurface("PMT_base_Surface", PMT_base_log, op_PMT_base);

		G4MaterialPropertiesTable* PMT_base_ST2 = new G4MaterialPropertiesTable();
		PMT_base_ST2->AddProperty("REFLECTIVITY", ephoton_teflon, reflectivity_teflon, num3);

		op_PMT_base->SetMaterialPropertiesTable(PMT_base_ST2);
		op_PMT_base->SetType(dielectric_metal);
		op_PMT_base->SetFinish(ground);
		op_PMT_base->SetPolish(0.2);
		op_PMT_base->SetModel(glisur);
	}


	// ------------- Surfaces --------------

	// PMT Quartz Window
	G4OpticalSurface* opPMT_window = new G4OpticalSurface("PMT_Surface");
	opPMT_window->SetType(dielectric_LUT);
	opPMT_window->SetFinish(polishedair);
	opPMT_window->SetModel(LUT);

	G4LogicalBorderSurface* PMT_Surface = new G4LogicalBorderSurface("PMT_BorderSurface1", pmtWindow_phys, expHall_phys, opPMT_window);
	G4OpticalSurface* PMTopticalSurface = dynamic_cast <G4OpticalSurface*>
		(PMT_Surface->GetSurface(pmtWindow_phys, expHall_phys)->
			GetSurfaceProperty());
	if (PMTopticalSurface) PMTopticalSurface->DumpInfo();

	// photocathode
	const G4int num2 = 2;
	G4double ephoton_abs[num2] = { 1.0 * eV , 10.0 * eV };
	G4double efficiency[num2] = { 1.0 , 1.0 };
	G4double reflectivity_ph[num2] = { 0.0, 0.0 };

	G4OpticalSurface* opPhCa = new G4OpticalSurface("PhCa_Surface");

	G4MaterialPropertiesTable* PhCaST2 = new G4MaterialPropertiesTable();
	PhCaST2->AddProperty("EFFICIENCY", ephoton_abs, efficiency, num2);
	PhCaST2->AddProperty("REFLECTIVITY", ephoton_abs, reflectivity_ph, num2);
	
	opPhCa->SetMaterialPropertiesTable(PhCaST2);
	opPhCa->SetType(dielectric_metal);
	opPhCa->SetFinish(polished);
	opPhCa->SetModel(glisur);

	G4LogicalSkinSurface* PhCa_Surface = new G4LogicalSkinSurface("PhCa_SkinSurface", pmtPhc_log, opPhCa);
	
	//always return the physical World
	return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
