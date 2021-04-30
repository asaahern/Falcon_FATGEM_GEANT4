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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "Randomize.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include "CLHEP/Units/PhysicalConstants.h"


using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(), 
   fParticleGun(0),fCounter(1)
{
    // define particle type and energy
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("opticalphoton");

    // create particle gun
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);
    fParticleGun->SetParticleDefinition(particle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    fParticleGun->GeneratePrimaryVertex(anEvent);

    // define initial position:
    G4double pitch  = 5.*mm;
    G4double fatgem_thickness= 4.6*mm;
    G4double x0     = 0;
    G4double y0     = 0;
    G4double sigma  = 0.33*mm; //sigma = (0.15, 0.21, 0.28, 0.32, 0.33) for R = (0.5, 0.75, 1, 1.25, 1.5)

    // distribute photons in different holes, overwriting initial values of x0, y0
    G4double n = 0;//G4UniformRand();
    if (n < 0.33){
        x0  = 0;
        y0  = 0;
    }
    else if(n<0.66){
        x0  = -pitch;
        y0  = 0;  
    }
    else {
        x0  = pitch/2;
        y0  = sqrt((pitch*pitch) - (pitch/2*pitch/2));  
    } 
 
    G4double posX   = G4RandGauss::shoot(x0, sigma);
    G4double posY   = G4RandGauss::shoot(y0, sigma);  
    G4double posZ   = (G4UniformRand()-0.5)*fatgem_thickness;

    if(fCounter>1000000) {posX=posY=0; posZ=17.3;} 

    fParticleGun->SetParticlePosition(G4ThreeVector(posX, posY, posZ));

    G4double wavelength = G4RandGauss::shoot(170.*nm, 6.*nm);//172 * nm;   // Xe scintillation
    G4double energy = h_Planck * c_light / wavelength;
    fParticleGun->SetParticleEnergy(energy);

    // define momentum direction: isotropic and random
    G4double theta, phi;
    G4ThreeVector momentumDirection;
    if(fCounter>1000000){
    	do {
		theta  = std::acos(-G4UniformRand())*rad;
	        phi    = G4UniformRand() * 360.0*deg;
        	momentumDirection.setRThetaPhi(1.,theta,phi);
	    } while (momentumDirection.angle(G4ThreeVector(0,0,-1))>atan(20./17.3));

    } else {
	theta  = std::acos(1.-2.*G4UniformRand())*rad;
        phi    = G4UniformRand() * 360.0*deg;
       	momentumDirection.setRThetaPhi(1.,theta,phi);
    }

    fParticleGun->SetParticleMomentumDirection(momentumDirection);
    fCounter++;
    //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1)); // launch all photons directed towards base

    // define polarization direction: random but contained in the plane perpendicular to the momentum direction of each photon
    G4ThreeVector normal(1., 0., 0.);
    G4ThreeVector product = normal.cross(momentumDirection);
    G4double modul2 = product * product;

    G4ThreeVector e_perpend(0., 0., 1.);
    if (modul2 > 0.) e_perpend = (1. / std::sqrt(modul2)) * product;
    G4ThreeVector e_paralle = e_perpend.cross(momentumDirection);

    G4double angle = G4UniformRand() * 360.0 * deg;
    G4ThreeVector polar = std::cos(angle) * e_paralle + std::sin(angle) * e_perpend;
    fParticleGun->SetParticlePolarization(polar);

}

