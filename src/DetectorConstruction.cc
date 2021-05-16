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
	//drift_2			= 20.5*mm;
	drift_2			= 3.*mm;
//	drift_2			= 15.*mm;

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
	G4int configuration = 6;
	
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

  G4MaterialPropertiesTable *myAcrylic = new G4MaterialPropertiesTable();

  const int myAcrSize = 504;
  G4double myAcrAbsEnergy2[myAcrSize];
  G4double myAcrAbsEnergy[myAcrSize] = {
    60,200,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,
    303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,
    331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,
    359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,
    388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,
    416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,
    444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,
    472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,
    500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,
    529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,
    557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,
    585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,
    613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,
    641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,
    670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,686,687,688,689,690,691,692,693,694,695,696,697,
    698,699,700,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,720,721,722,723,724,725,
    726,727,728,729,730,731,732,733,734,735,736,737,738,739,740,741,742,743,744,745,746,747,748,749,750,751,752,753,
    754,755,756,757,758,759,760,761,762,763,764,765,766,767,768,769,770,771,772,773,774,775,776,777,778,779,780,800
  };

  G4double myAcrAbs2[myAcrSize];
  G4double myAcrAbs[myAcrSize] = {
    100*nm,100*nm,1.32284,1.29312,1.30753,1.30943,1.26653,1.32372,1.29873,1.3092,1.30371,1.30577,1.32406,1.31243,1.3171,1.28569,
    1.30494,1.33508,1.30388,1.29962,1.29142,1.30584,1.28672,1.31488,1.30954,1.31865,1.29685,1.28557,1.3144,1.32283,1.31425,
    1.27107,1.302,1.29515,1.29604,1.31832,1.27713,1.33705,1.31362,1.34722,1.34076,1.33069,1.35506,1.35238,1.33769,1.33566,1.35435,
    1.32866,1.35414,1.34693,1.40201,1.3828,1.39343,1.31766,1.30547,1.34563,1.31455,1.34873,1.32728,1.281,1.35277,1.31982,1.30124,
    1.28945,1.29276,1.2827,1.28362,1.29432,1.25258,1.21353,1.27814,1.26223,1.22418,1.23389,1.19746,1.20099,1.23181,1.2066,1.21319,
    1.21229,1.22824,1.2102,1.20421,1.17812,1.19581,1.23024,1.15775,1.21839,1.20776,1.21069,1.1799,1.20016,1.1395,1.19189,1.16832,
    1.13224,1.24984,1.50017,1.77389,1.78046,2.09596,1.93195,2.85164,3.44997,2.62686,4.97381,5.97033,14.0861,18.7694,22.3729,23.5664,
    24.2125,24.918,27.6302,32.0202,38.1495,45.9201,55.8979,68.4572,83.8411,102.648,125.635,153.884,187.777,229.237,278.717,339.222,
    410.576,495.219,597.237,707.002,839.051,980.058,1143.41,1321,1510.09,1705.82,1932.72,2165.05,2359.05,2580.08,2771.07,3033.05,
    3167.62,3409.31,3479.23,3707.19,3851.96,3897.3,4225.34,4205.58,4315.72,4419.21,4657.92,4658.72,4759.31,4806.68,4802.77,
    5000.05,5121.6,5309.11,5333.83,5338.58,5643.52,5467.12,5694.31,5779.82,5594.7,5879.86,5981.84,5787.35,6262.51,6254.24,6416.64,
    6962.87,6881.92,7197.05,7142.11,7293.65,7876.71,8372.81,7681.55,7740.36,7990.87,8104.59,8363.88,8592.66,8136.1,8629.33,8565.14,
    8687.34,9020.28,9114.16,8924.67,9189.4,9277.85,9158.99,9574.41,10854.5,9460.01,9664.68,10045.7,10371.2,9946.08,10623.3,11275.5,
    10485.8,11033.6,11213.6,11160.2,11517.3,11466.1,11184,11771.9,11520.6,12506,11235.2,12208.5,12694.5,11598.4,13119.1,11726.9,
    12782.3,12282.6,13379.2,14005.3,12249.8,12855.6,12727.9,12390.6,13738.5,14171.4,13698.8,14075.3,14030.4,15346.9,13350.7,14807,
    14633.6,17460.4,15319,15551.3,15958,14883.8,14163.4,15687.7,14723.8,15681.3,16460.4,15074.8,15287.5,15808.1,16863.4,15432.6,
    16162.6,16028.7,17023.9,16573.2,14573.7,15343.5,16263,14356.3,14277.8,15545.9,16271.4,14849.5,14983.5,15728.8,15986.8,14954.6,
    16380.3,14752.6,14499.7,16345.6,18455.6,13926.8,17307.8,15290.8,15605.9,16007.6,17090.7,18384.3,16610.8,16896,17407.5,19442.5,
    19325.7,20769.7,18066.3,17648.5,18863.4,22099.9,20030.3,23828.9,19256,20513.2,26932.9,20529.7,18444.9,20507.2,20803.8,23749.3,
    20872.5,21161.2,20857.3,21194.2,22768.2,20472.9,18433.5,21932.4,19893.6,21459.4,19480.1,19474.4,22644.7,25062.3,20185.8,
    21164,20515.7,21679.3,23751.5,21387,24228.8,22827.1,22586.3,21646.6,19780.9,18207.4,17148.8,17317.4,14297.8,15752.4,15363.6,
    14105.6,12582.5,12050,11068.4,10582.3,10052.6,9658.73,9156.67,9611.99,8740.63,9141.28,8545.84,8675.79,8096.21,8927.55,8966.51,
    9515.41,8602.26,9440.08,10283.8,10265.4,10888.3,10528.8,12258.1,12131.1,14388.4,15769.1,16420.3,17624.4,15974.7,19390.5,21176.3,
    19823.3,24277.3,20107.8,21017.3,26572.9,24786.3,23206.9,25619.9,29196.8,22952.4,27438.1,22454,30586.4,20477.8,22649.7,22041.4,
    20035.7,17730.3,17692.6,20030.4,19537,16425.3,16144.3,19343.1,15152.6,14117.4,14447.6,15711.6,13035.8,14160.7,15208.4,14528.4,
    13168.9,13102.1,20340.7,14727,14351,13941.1,13676.7,14920.8,14700,13243,12718.4,15343.7,13408.6,13693.9,12562.9,13435.8,
    13291.9,12706.8,13482.3,12490.9,12314.1,11564,10328.3,10482,9638.82,9286.56,8971.19,8252.47,7627.67,7387.04,6738.8,6368.22,6173.61,
    5617.03,5641.01,5206.5,4730.11,4560.32,4190.11,3948.63,3640.44,3263.39,3022.28,2781.11,2542.86,2304.26,2125.74,1978.46,
    1820.81,1705.02,1617.41,1552.42,1469.53,1420.46,1389.11,1373.85,1360.76,1356.43,1357.77,1376.89,1393.14,1407.4,1442.59,1487.05,
    1538.95,1622.01,1698.69,1792.43,1935.75,2089.58,2274.06,2455.86,2639.89,2934.61,3027.92,3422.37,3714.02,4134.33,4542.13,5160.33,
    5363.87,5977.23,6442.42,7013.28,7784.87,8675.49,8942.09,9518.33,9760.12,10260.8,9922.78,10979.6,10000.5,10391.3,9134.34,
    9521.78,9595.54,9179.65,8279.66,8223.8,8020.19,6951.69,6875.43,6592.61,6159.45,6023.03,6023.03
  };
  for (int i=myAcrSize-1; i>=0; i--) {
    myAcrAbsEnergy2[i] = h_Planck*c_light / (nm*myAcrAbsEnergy[myAcrSize-1-i]);
    myAcrAbs2[i] = myAcrAbs[myAcrSize-1-i]*mm;
  }

  G4double RINDEX_value1[53] =  { 60 , 200 , 300 , 310 , 320 , 330 , 340 ,       350 , 360 , 370 , 380 , 390 , 400 , 410 , 420 , 430 , 440 , 450 , 460 , 470 , 480 , 490 , 500 , 510 , 520 , 530 , 540 ,       550 , 560 , 570 , 580 , 590 , 600 , 610 , 620 , 630 , 640 , 650 , 660 , 670 , 680 , 690 , 700 , 710 , 720 , 730 , 740 ,       750 , 760 , 770 , 780 , 790 , 800  } ;
  G4double RINDEX_value2[53] = { 1.65, 1.65, 1.527 , 1.524 , 1.521 , 1.519 , 1.516 , 1.514 , 1.512 , 1.510 , 1.509 , 1.507 , 1.506 , 1.505 , 1.503 , 1.502 , 1.501 , 1.500 , 1.499 , 1.499 , 1.498 , 1.497 , 1.496 , 1.496 , 1.495 , 1.494 , 1.494 , 1.493 , 1.493 , 1.492 , 1.492 , 1.491 , 1.491 , 1.490 , 1.490 , 1.490 , 1.489 , 1.489 , 1.488 , 1.488 , 1.488 , 1.488 , 1.487 , 1.487 , 1.487 , 1.486 , 1.486 , 1.486 , 1.486 , 1.485 , 1.485 , 1.485 , 1.485 } ;

  G4double myacr_RINDEX_ene[53], RINDEX_value[53] ;
  for (int i=52;i>=0; i--) {
    myacr_RINDEX_ene[i] =  h_Planck*c_light/ (RINDEX_value1[52-i]*nm);
    RINDEX_value[i] = RINDEX_value2[52-i];
  }

  myAcrylic->AddProperty("RINDEX",  myacr_RINDEX_ene  , RINDEX_value,   53);
  myAcrylic->AddProperty("ABSLENGTH", myAcrAbsEnergy2, myAcrAbs2 , myAcrSize);
  PMMA_mat->SetMaterialPropertiesTable(myAcrylic);


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
	G4Material* fPEN = new G4Material (name="PENQ51",density,ncomponents=3,kStateSolid);
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
		  PEN_WLS_ABSORPTION_VAL[dim] = 1.*nm;
                  PEN_ABSORPTION_VAL[dim]  = 1000.*m;
                  PEN_RINDEX[dim]          = 1.75;
                  PEN_RAYL[dim]            = 1000.*m;
		  PEN_RAYL_FOIL[dim]       = 1000.*m;
          }
          else {
                // in visible absorption and scattering meaningful
		  PEN_WLS_ABSORPTION_VAL[dim] = 1000*m;
                  if(configuration==1 || configuration==9 || configuration==10 || configuration==14) PEN_ABSORPTION_VAL[dim]  = 5.*cm; // Efremenko et al.
                  else PEN_ABSORPTION_VAL[dim]  = 1.*cm; // Efremenko et al.
                  PEN_RINDEX[dim]          = 1.75;
                  PEN_RAYL[dim]            = 10*cm; // educated guess
		  PEN_RAYL_FOIL[dim]       = 150.*micrometer; // educated guess
          }
          dim++;
	}
	fpen_absorption.close();

	myPEN->AddProperty("WLSABSLENGTH",PEN_ENE2, PEN_WLS_ABSORPTION_VAL, dim);

	if(configuration==6 || configuration==7 || configuration==11) { // using PMMA for the tile material, except with WLS property to emulate the TPB coating
        	myPEN->AddProperty("RINDEX",  myacr_RINDEX_ene  , RINDEX_value,   53);
        	myPEN->AddProperty("ABSLENGTH", myAcrAbsEnergy2, myAcrAbs2 , myAcrSize);
	} else {
		myPEN->AddProperty("RINDEX",      PEN_ENE2, PEN_RINDEX,             dim);
		myPEN->AddProperty("ABSLENGTH",   PEN_ENE2, PEN_ABSORPTION_VAL,     dim);
		myPEN->AddProperty("RAYLEIGH",    PEN_ENE2, PEN_RAYL,               dim);
	}

	G4double holeWLSE=1.;
	if(configuration==5 || configuration==8 || configuration==9 || configuration==10 || configuration==12 || configuration==14) holeWLSE=0.091;
	else if (configuration==4) holeWLSE=0.5;
	else if (configuration==2 || configuration==13) holeWLSE=0;

	myPEN->AddConstProperty("WLSMEANNUMBERPHOTONS", holeWLSE);
	myPEN->AddConstProperty("WLSTIMECONSTANT",20.*ns);
	fPEN->SetMaterialPropertiesTable(myPEN);

	myPENfoil->AddProperty("WLSABSLENGTH",PEN_ENE2, PEN_WLS_ABSORPTION_VAL, dim);
	myPENfoil->AddProperty("ABSLENGTH",   PEN_ENE2, PEN_ABSORPTION_VAL,     dim);
	myPENfoil->AddProperty("RINDEX",      PEN_ENE2, PEN_RINDEX,             dim);
	myPENfoil->AddProperty("RAYLEIGH",    PEN_ENE2, PEN_RAYL_FOIL,          dim);

	G4double foilWLSE=1.;
	if(configuration==5 || configuration==10 || configuration==13) foilWLSE=0.34;
	else if(configuration==4) foilWLSE=0.5;
	else if(configuration==2) foilWLSE=0.;
	
	myPENfoil->AddConstProperty("WLSMEANNUMBERPHOTONS", foilWLSE);
	myPENfoil->AddConstProperty("WLSTIMECONSTANT",20.*ns);
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
//	G4VPhysicalVolume* pmtWindow_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -drift_2 -fatgem_thickness - PMT_thickness), pmtWindow_log, "PMT_win", expHall_log, false, 0, checkOverlaps);
	pmtWindow_log -> SetVisAttributes(blue);
	
	// photocathode
	G4Tubs* PMT_phc = new G4Tubs("PMT_phCa", 0.*cm, PMT_radius, phc_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* pmtPhc_log = new G4LogicalVolume(PMT_phc, PMT_phc_mat, "PMT_phCa", 0, 0, 0);
//	G4VPhysicalVolume* pmtPhc_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -drift_2 -fatgem_thickness - 2*PMT_thickness - phc_thickness), pmtPhc_log, "PMT_phCa", expHall_log, false, 0, checkOverlaps);
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
	G4Tubs* fatgem_wls = new G4Tubs("fatgem_wls", 0.*cm, fatgem_radius, 25.*micrometer/2., 0. * deg, 360. * deg);
	//G4LogicalVolume* fatgem_log = new G4LogicalVolume(fatgem, PMMA_mat, "fatgem", 0, 0, 0);
	G4LogicalVolume* fatgem_log = new G4LogicalVolume(fatgem, fPEN, "fatgem", 0, 0, 0);
	G4LogicalVolume* fatgemwls_log = new G4LogicalVolume(fatgem_wls, fPENfoil, "fatgem_wls", 0, 0, 0);
	G4VPhysicalVolume* fatgem_phys = new G4PVPlacement(0, G4ThreeVector(), fatgem_log, "fatgem", expHall_log, false, 0, checkOverlaps);
	G4VPhysicalVolume* fatgemwls_phys=NULL;
	if(configuration==5 || configuration==6 || configuration==7 || configuration==10 || configuration==11 || configuration==13) 
		fatgemwls_phys = new G4PVPlacement(0, G4ThreeVector(0,0,fatgem_thickness-25.*micrometer/2.), fatgemwls_log, "fatgem_wls", fatgem_log, false, 0, checkOverlaps);
	fatgem_log -> SetVisAttributes(red);

	G4Tubs* hole = new G4Tubs("hole", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4Tubs* hole_wls = new G4Tubs("hole_wls", 0, hole_radius, 25.*micrometer/2., 0. * deg, 360. * deg);
	G4LogicalVolume* hole_log = new G4LogicalVolume(hole, matXe, "hole", 0, 0, 0);
	G4LogicalVolume* holewls_log = new G4LogicalVolume(hole_wls, matXe, "hole_wls", 0, 0, 0);
	hole_log->SetVisAttributes(white);
	
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

	G4VPhysicalVolume *hole_phys[31], *holewls_phys[31];
	for(int i=0; i<31; i++) {
	  hole_phys[i] = new G4PVPlacement(0, G4ThreeVector(hole_x[i], hole_y[i], 0), hole_log, "hole", fatgem_log, true, i, checkOverlaps);
	  if (configuration==5 || configuration==6 || configuration==7 || configuration==10 || configuration==11 || configuration==13) holewls_phys[i] = new G4PVPlacement(0, G4ThreeVector(hole_x[i], hole_y[i], 0), holewls_log, "hole_wls", fatgemwls_log, true, i, checkOverlaps);
	}
	
	// optical surface 
	G4OpticalSurface* OS_fatgem_hole = new G4OpticalSurface("fatgem_hole_surface");

	for(int i=0; i<31; i++) {
		new G4LogicalBorderSurface("fatgem_hole_surface", fatgem_phys, hole_phys[i], OS_fatgem_hole);
		new G4LogicalBorderSurface("fatgem_hole_surface", hole_phys[i], fatgem_phys, OS_fatgem_hole);
		if (configuration==5 || configuration==6 || configuration==7 || configuration==10 || configuration==11 || configuration==13) {
			new G4LogicalBorderSurface("fatgem_hole_wls_surface", fatgemwls_phys, holewls_phys[i], OS_fatgem_hole);
			new G4LogicalBorderSurface("fatgem_hole_wls_surface", holewls_phys[i], fatgemwls_phys, OS_fatgem_hole);
		}
	}

	OS_fatgem_hole->SetType(dielectric_dielectric);
	OS_fatgem_hole->SetFinish(ground);
	OS_fatgem_hole->SetModel(glisur);
	OS_fatgem_hole->SetPolish(0.1);

	const G4int num = 11;
	G4double ephoton_fatgem[num] = { 0.12*eV, 0.31*eV, 0.62*eV, 1.24*eV, 1.77*eV, 2.067*eV, 2.48*eV, 3.1*eV, 4.13*eV, 6.2*eV,  8.9*eV };
	G4double reflectivity_fatgem[num] = { 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999 };
	//G4double transmittance_fatgem[num] = { .9, .9, .9, .9, .9, .9, .9, .9, .9, .9, .9 };
	G4MaterialPropertiesTable* R_fatgem = new G4MaterialPropertiesTable();
	R_fatgem->AddProperty("REFLECTIVITY", ephoton_fatgem, reflectivity_fatgem, num);
	// line below commented out to rely on rindex
	//	R_fatgem->AddProperty("TRANSMITTANCE", ephoton_fatgem, transmittance_fatgem, num);
	OS_fatgem_hole->SetMaterialPropertiesTable(R_fatgem);

	G4OpticalSurface* OS_clevios = new G4OpticalSurface("clevios_surface");
	G4OpticalSurface* OS_clevios_outgoing = new G4OpticalSurface("clevios_surface_outgoing");
	G4OpticalSurface* OS_esr = new G4OpticalSurface("esr_surface");

	if (configuration==5 || configuration==6 || configuration==7 || configuration==10 || configuration==11 || configuration==13)	{
		new G4LogicalBorderSurface("clevios_surface_outgoing", fatgemwls_phys, expHall_phys, OS_clevios_outgoing);
		new G4LogicalBorderSurface("clevios_surface", expHall_phys, fatgemwls_phys, OS_clevios);
	}
	if (configuration==11) {
		new G4LogicalBorderSurface("ESR_surface", fatgem_phys, fatgemwls_phys, OS_esr);
		new G4LogicalBorderSurface("ESR_surface", fatgemwls_phys, fatgem_phys, OS_esr);
	}

	new G4LogicalBorderSurface("clevios_surface_outgoing", fatgem_phys, expHall_phys, OS_clevios_outgoing);
	new G4LogicalBorderSurface("clevios_surface", expHall_phys, fatgem_phys, OS_clevios);

	OS_clevios->SetType(dielectric_dielectric);
        OS_clevios->SetFinish(ground);
        OS_clevios->SetModel(glisur);
        OS_clevios->SetPolish(0.95);

	OS_clevios_outgoing->SetType(dielectric_dielectric);
        OS_clevios_outgoing->SetFinish(ground);
        OS_clevios_outgoing->SetModel(glisur);
        OS_clevios_outgoing->SetPolish(0.95);

	OS_esr->SetType(dielectric_metal);
        OS_esr->SetFinish(polished);
        OS_esr->SetModel(glisur);
        OS_esr->SetPolish(1.);

	//const G4int num = 11;
        G4double ephoton_clevios[num] = { 0.12*eV, 0.31*eV, 0.62*eV, 1.24*eV, 1.77*eV, 2.067*eV, 2.48*eV, 3.1*eV, 4.13*eV, 6.2*eV,  8.9*eV };
	// 1% loss in VIS, 20% loss in VUV

        G4double reflectivity_clevios[num] = { 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.75, 0.75 };
        //G4double reflectivity_clevios[num] = { 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99999, 0.99999 };
	//G4double transmittance_clevios[num] = { .9, .9, .9, .9, .9, .9, .9, .9, .9, .9, .9 };

        G4double reflectivity_esr[num] = { 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98 };

	// assuming actual SS_R of 55% for visible and 5% for VUV and 29% coverage with the mesh
	// [inspired by Figure 1 of Karlsson and Ribbing, J. Appl. Phys. Vol 53. No. 9 Sept. 1982]
	// for VIS 0.29*(1.0-0.55) is then lost -> in case of dielectric_dielectric implementation must be then used as (1-R)
        G4double reflectivity_cu[num] =  { 0.87, 0.87, 0.87, 0.87, 0.87, 0.87, 0.87, 0.87, 0.87, 0.28, 0.28 };

	// this is for photons entering the plate. 0.71% hits the plastic, and 4% of them is reflected: 0.96*0.71=0.68
	G4double transmittance_cu[num] = { 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68 };

	G4MaterialPropertiesTable* R_clevios = new G4MaterialPropertiesTable();
	G4MaterialPropertiesTable* R_clevios_outgoing = new G4MaterialPropertiesTable();

	if (configuration==7 || configuration==12 || configuration==14) {
		R_clevios_outgoing->AddProperty("REFLECTIVITY", ephoton_clevios, reflectivity_cu, num);
		R_clevios->AddProperty("REFLECTIVITY", ephoton_clevios, transmittance_cu, num);
	} else {
		R_clevios->AddProperty("REFLECTIVITY", ephoton_clevios, reflectivity_clevios, num);
		R_clevios_outgoing->AddProperty("REFLECTIVITY", ephoton_clevios, reflectivity_clevios, num);
	}

	// line below commented out to rely on rindex
	// R_clevios->AddProperty("TRANSMITTANCE", ephoton_clevios, transmittance_clevios, num);
	OS_clevios->SetMaterialPropertiesTable(R_clevios);
	OS_clevios_outgoing->SetMaterialPropertiesTable(R_clevios_outgoing);

	G4MaterialPropertiesTable* R_esr = new G4MaterialPropertiesTable();
	R_esr->AddProperty("REFLECTIVITY", ephoton_clevios, reflectivity_esr, num);
	OS_esr->SetMaterialPropertiesTable(R_esr);

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
		//G4Tubs* PMT_base = new G4Tubs("PMT_base", PMT_radius, base_radius, PMT_thickness, 0. * deg, 360. * deg);
		G4Tubs* PMT_base = new G4Tubs("PMT_base", 0., base_radius, PMT_thickness, 0. * deg, 360. * deg);
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
/*
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
*/	
	//always return the physical World
	return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
