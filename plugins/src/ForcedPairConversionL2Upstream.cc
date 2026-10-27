#include "ForcedPairConversionL2Upstream.hh"

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

// Upstream (first) sensor per half of Layer 2 -- same axial/stereo-vs-
// layer-parity convention as Layer 1 (top: axial, bottom: stereo). Matching
// only the sensor volumes (not the halfmodule wrapper) means the forced
// conversion always lands in actual silicon.
const std::vector<std::string> ForcedPairConversionL2Upstream::kTargetVolumes = {
    "module_L2t_halfmodule_axial_sensor",
    "module_L2b_halfmodule_stereo_sensor",
};

// All four Layer 1 sensor volumes (both the upstream and downstream sensor
// of each half) -- the union of ForcedPairConversion::kTargetVolumes and
// ForcedPairConversionL1Downstream::kTargetVolumes. Natural conversion is
// suppressed here so every event's primary photon survives the whole of
// Layer 1 to reach the Layer 2 target above instead of some fraction
// converting naturally in L1 first.
const std::vector<std::string> ForcedPairConversionL2Upstream::kUpstreamVolumes = {
    "module_L1t_halfmodule_axial_sensor",
    "module_L1b_halfmodule_stereo_sensor",
    "module_L1t_halfmodule_stereo_sensor",
    "module_L1b_halfmodule_axial_sensor",
};

namespace {
bool matchesAnyVolume(const G4String &name, const std::vector<std::string> &targets) {
  for (const std::string &target : targets) {
    if (name.find(target) != std::string::npos) {
      return true;
    }
  }
  return false;
}
} // namespace

G4VEmProcess *ForcedPairConversionL2Upstream::findConversionProcess(const G4Step *step) {
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

void ForcedPairConversionL2Upstream::stepping(const G4Step *step) {

  G4Track *track = step->GetTrack();

  // Only the primary (generator-level) photon is ever a candidate for
  // forcing, or relevant to the Layer-1 deactivation toggle below. Any
  // photon with a parent (e.g. from eBrem on a forced e+/e-) is skipped
  // here, so a forced conversion can never cascade downstream.
  if (track->GetDefinition() != G4Gamma::Gamma() || track->GetParentID() != 0) {
    return;
  }

  // Suppress natural conversion while the primary photon transits any of
  // Layer 1's four sensor volumes, and restore it the moment it leaves the
  // last one, so every event's photon survives all of Layer 1 intact to
  // reach the Layer 2 target instead of some fraction converting naturally
  // along the way. Checked ahead of the hasConverted_ guard below since
  // this must run on every step of the primary photon's track, not just
  // while it's still a forcing candidate.
  //
  // Deliberately keyed off GetNextVolume() (the volume the upcoming step
  // will be in), not GetVolume() (the volume the just-finished step was
  // in): stepping() runs after a step's physics already executed, so by the
  // time GetVolume() reports "inside a Layer 1 sensor", any natural
  // conversion during that very step has already happened -- one step too
  // late to suppress it. GetNextVolume() lets the toggle take effect before
  // the photon's first step inside the sensor.
  G4VPhysicalVolume *volume = track->GetVolume();
  G4VPhysicalVolume *nextVolume = track->GetNextVolume();
  const bool enteringOrInUpstream = nextVolume != nullptr && matchesAnyVolume(nextVolume->GetName(), kUpstreamVolumes);
  if (enteringOrInUpstream != convDeactivated_) {
    G4VEmProcess *conv = findConversionProcess(step);
    if (conv != nullptr) {
      track->GetDefinition()->GetProcessManager()->SetProcessActivation(conv, !enteringOrInUpstream);
      convDeactivated_ = enteringOrInUpstream;
    }
  }

  if (hasConverted_) {
    // Primary photon has already converted (naturally or forced) this
    // event; everything else (including any secondary photons produced by
    // the resulting e+/e-) is left to propagate under normal physics.
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

  if (volume == nullptr) {
    return;
  }
  const G4String &volumeName = volume->GetName();
  if (!matchesAnyVolume(volumeName, kTargetVolumes)) {
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
    G4cout << "[ ForcedPairConversionL2Upstream ]: forced conversion in " << volumeName
           << " with weight " << weight << G4endl;
  }
}

void ForcedPairConversionL2Upstream::endEvent(const G4Event *) {
  // Defensive restore: normally the stepping action re-activates "conv" the
  // moment the primary photon leaves the last Layer-1 volume, but if that
  // track is ever killed by some other process while still inside one
  // (e.g. photoelectric absorption), no further step would run the restore
  // -- leaving conv deactivated for every other gamma, potentially leaking
  // into the next event.
  if (convDeactivated_ && convProcess_ != nullptr) {
    G4Gamma::Gamma()->GetProcessManager()->SetProcessActivation(convProcess_, true);
  }
  convDeactivated_ = false;
  hasConverted_ = false;
}

void ForcedPairConversionL2Upstream::endRun(const G4Run *) {
  G4cout << "[ ForcedPairConversionL2Upstream ]: forced conversions in " << forcedCount_ << " events" << G4endl;
}

} // namespace slic
DECLARE_PLUGIN(slic, ForcedPairConversionL2Upstream)
