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
#include "Par04DetectorMessenger.hh"
#include "Par04DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04DetectorMessenger::Par04DetectorMessenger(
  Par04DetectorConstruction* aDetector)
  : G4UImessenger()
  , fDetector(aDetector)
{
  fExampleDir = new G4UIdirectory("/Par04/");
  fExampleDir->SetGuidance("UI commands specific to this example");

  fDetectorDir = new G4UIdirectory("/Par04/detector/");
  fDetectorDir->SetGuidance("Detector construction UI commands");

  fPrintCmd = new G4UIcmdWithoutParameter("/Par04/detector/print", this);
  fPrintCmd->SetGuidance("Print current settings.");
  
  fNbLayersCmd =
    new G4UIcmdWithAnInteger("/Par04/detector/setNbOfLayers", this);
  fNbLayersCmd->SetGuidance("Set number of layers.");
  fNbLayersCmd->SetParameterName("NbLayers", false);
  fNbLayersCmd->SetRange("NbLayers>0");
  fNbLayersCmd->AvailableForStates(G4State_PreInit);
  fNbLayersCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04DetectorMessenger::~Par04DetectorMessenger()
{
  delete fPrintCmd;
  delete fNbLayersCmd;
  delete fDetectorDir;
  delete fExampleDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04DetectorMessenger::SetNewValue(G4UIcommand* aCommand,
                                         G4String aNewValue)
{
  if(aCommand == fPrintCmd)
  {
    fDetector->Print();
  }
  else if(aCommand == fNbLayersCmd)
  {
    fDetector->SetNbOfLayers(fNbLayersCmd->GetNewIntValue(aNewValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String Par04DetectorMessenger::GetCurrentValue(G4UIcommand* aCommand)
{
  G4String cv;

  if(aCommand == fNbLayersCmd)
  {
    cv = fNbLayersCmd->ConvertToString(fDetector->GetNbOfLayers());
  }
  return cv;
}
