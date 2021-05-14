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

#include "SteppingAction.hh"
#include "G4SystemOfUnits.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "RunAction.hh"
#include "CLHEP/Units/PhysicalConstants.h"


#include <fstream>
#include <iostream>
#include <sstream>

using namespace CLHEP;
using namespace std;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* runAction)
: G4UserSteppingAction(), fRunAction(runAction)
{
  myfile.open("out.txt", ios::out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ ; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{

    G4StepPoint* endPoint = step->GetPostStepPoint();
    if (endPoint->GetStepStatus() == fGeomBoundary) {
	if (abs(endPoint->GetPosition().z()) > 5 ) {
		if (G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()>1000000) 
			myfile << "S1 ";
		else
			myfile << "S2 ";

		if(endPoint->GetPosition().z()>5)
			myfile << "TPC ";
		else
			myfile << "PMT ";

		myfile << endPoint->GetPosition().x() << " " << endPoint->GetPosition().y() << " ";
	
		if (h_Planck*c_light/step->GetTrack()->GetKineticEnergy()/nm>250)
			myfile << "VIS" << G4endl;
		else
			myfile << "VUV" << G4endl;

		step->GetTrack()->SetTrackStatus(fStopAndKill );
	}
        const G4String touchableName =
            endPoint->GetTouchable()->GetVolume()->GetName();
        //G4cout << touchableName << G4endl;
        if (touchableName == "PMT_phCa") {
            fRunAction->AddEdep(1);
        }
    }
}
