#include "ForcedPairConversion.hh"

// SLIC
#include "MCParticleManager.hh"

// Geant4
#include "G4DynamicParticle.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Gamma.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "Randomize.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "G4VEmModel.hh"
#include "G4VEmProcess.hh"
#include "G4VPhysicalVolume.hh"

// LCIO
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"

namespace slic {

// Start with layer 1 only; retune/retry other layers later by changing these
// strings. Each entry names the upstream (first, unshielded) sensor for one
// half: top -> axial, bottom -> stereo (see header comment for why these
// differ by half). Matching only the sensor volumes (not the halfmodule
// wrapper) means the forced conversion always lands in actual silicon.
const std::vector<std::string> ForcedPairConversion::kTargetVolumes = {
    "module_L1t_halfmodule_axial_sensor",
    "module_L1b_halfmodule_stereo_sensor",
};

G4VEmProcess *ForcedPairConversion::findConversionProcess(const G4Step *step) {
  if (convProcess_ != nullptr) {
    return convProcess_;
  }

  G4ProcessManager *pm = step->GetTrack()->GetDefinition()->GetProcessManager();
  G4ProcessVector *processes = pm->GetProcessList();
  for (G4int i = 0; i < processes->size(); i++) {
    if ((*processes)[i]->GetProcessName() == "conv") {
      convProcess_ = static_cast<G4VEmProcess *>((*processes)[i]);
      break;
    }
  }
  return convProcess_;
}

void ForcedPairConversion::stepping(const G4Step *step) {

  if (hasConverted_) {
    // Primary photon has already converted (naturally or forced) this
    // event; everything else (including any secondary photons produced by
    // the resulting e+/e-) is left to propagate under normal physics.
    return;
  }

  G4Track *track = step->GetTrack();

  // Only the primary (generator-level) photon is ever a candidate for
  // forcing. Any photon with a parent (e.g. from eBrem on a forced e+/e-)
  // is skipped here, so a forced conversion can never cascade downstream.
  if (track->GetDefinition() != G4Gamma::Gamma() || track->GetParentID() != 0) {
    return;
  }

  // If the photon really converted on its own in this step, we're done:
  // nothing to force, and the event is unweighted (weight stays at 1).
  const std::vector<const G4Track *> *secondaries = step->GetSecondaryInCurrentStep();
  if (secondaries->size() > 0 &&
      (*secondaries)[0]->GetCreatorProcess() != nullptr &&
      (*secondaries)[0]->GetCreatorProcess()->GetProcessName() == "conv") {
    hasConverted_ = true;
    return;
  }

  G4VPhysicalVolume *volume = track->GetVolume();
  if (volume == nullptr) {
    return;
  }
  const G4String &volumeName = volume->GetName();
  bool inTargetSensor = false;
  for (const std::string &target : kTargetVolumes) {
    if (volumeName.find(target) != std::string::npos) {
      inTargetSensor = true;
      break;
    }
  }
  if (!inTargetSensor) {
    // Not yet in (or already past) the forcing sensor.
    return;
  }

  // Only force at the far boundary of the layer, i.e. once we know the
  // photon transited the whole slab without converting naturally.
  if (step->GetPostStepPoint()->GetStepStatus() != fGeomBoundary) {
    return;
  }

  G4VEmProcess *conv = findConversionProcess(step);
  if (conv == nullptr) {
    return;
  }

  const G4double kinEnergy = track->GetKineticEnergy();
  const G4MaterialCutsCouple *couple = track->GetMaterialCutsCouple();

  // Local conversion probability for this slab: 1 - exp(-t/lambda), with
  // 1/lambda given by the real cross section per unit volume. This is the
  // weight the event is worth once we force the conversion to happen here.
  const G4double sigma = conv->CrossSectionPerVolume(kinEnergy, couple);
  const G4double stepLength = step->GetStepLength();
  const G4double weight = 1.0 - std::exp(-sigma * stepLength);

  // Sample the actual e+/e- kinematics from the same model Geant4 would
  // have used, then hand the secondaries to the track stack and kill the
  // photon -- the net effect is a conversion at this point, but driven by
  // the real model rather than hand-rolled kinematics.
  size_t idx = (size_t)couple->GetIndex();
  G4VEmModel *model = conv->SelectModelForMaterial(kinEnergy, idx);

  std::vector<G4DynamicParticle *> secondaryParticles;
  model->SampleSecondaries(&secondaryParticles, couple, track->GetDynamicParticle());

  // Place the conversion at a random depth within the slab rather than
  // always at the exit face, so the resulting e+/e- can deposit charge in
  // this sensor too. The conversion probability per slab is small, so a
  // uniform draw along the chord is an adequate stand-in for the (more
  // correct but barely different) exponential-attenuation distribution.
  const G4ThreeVector &entryPosition = step->GetPreStepPoint()->GetPosition();
  const G4ThreeVector &exitPosition = step->GetPostStepPoint()->GetPosition();
  const G4double entryTime = step->GetPreStepPoint()->GetGlobalTime();
  const G4double exitTime = step->GetPostStepPoint()->GetGlobalTime();

  const G4double u = G4UniformRand();
  const G4ThreeVector position = entryPosition + u * (exitPosition - entryPosition);
  const G4double globalTime = entryTime + u * (exitTime - entryTime);

  // Push onto the current step's own secondary bucket (the same G4TrackVector
  // that G4VEmModel::SampleSecondaries() and friends use internally, exposed
  // via G4Step::GetfSecondary()/NewSecondaryVector()). This is the native
  // Geant4 hook for adding secondaries during a stepping action: after this
  // UserSteppingAction call returns, G4SteppingManager/G4TrackingManager
  // harvest this same vector (via GimmeSecondaries()) and hand it to
  // G4EventManager::StackTracks(), which is what assigns each track a real
  // (non-zero) track ID via its internal counter. Calling StackTracks()
  // directly isn't an option -- it's a private method -- and pushing onto
  // the G4StackManager directly (the original approach) bypasses ID
  // assignment entirely, leaving every forced secondary with the default
  // trackID of 0. That in turn breaks LCIO's SimTrackerHit->MCParticle
  // lookup once more than one forced conversion occurs in a run.
  G4TrackVector *secondaryBucket = const_cast<G4Step *>(step)->GetfSecondary();
  for (G4DynamicParticle *dp : secondaryParticles) {
    G4Track *secondaryTrack = new G4Track(dp, globalTime, position);
    secondaryTrack->SetParentID(track->GetTrackID());
    secondaryBucket->push_back(secondaryTrack);
  }

  track->SetTrackStatus(fStopAndKill);

  hasConverted_ = true;
  ++forcedCount_;

  // Stash the weight on the MCParticle collection under the key that
  // LcioManager::createLCEvent() reads at end-of-event
  // (LcioManager.cc:331), since by the time any plugin's endEvent() hook
  // runs the LCEvent has already been created and its weight already read.
  LCCollectionVec *mcpColl = MCParticleManager::instance()->getMCParticleCollection();
  if (mcpColl != nullptr) {
    mcpColl->parameters().setValue("_weight", (float)weight);
  }

  if (verbose_ > 1) {
    G4cout << "[ ForcedPairConversion ]: forced conversion in " << volumeName
           << " with weight " << weight << G4endl;
  }
}

void ForcedPairConversion::endEvent(const G4Event *) {
  hasConverted_ = false;
}

void ForcedPairConversion::endRun(const G4Run *) {
  G4cout << "[ ForcedPairConversion ]: forced conversions in " << forcedCount_ << " events" << G4endl;
}

} // namespace slic
DECLARE_PLUGIN(slic, ForcedPairConversion)
