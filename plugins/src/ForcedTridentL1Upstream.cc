#include "ForcedTridentL1Upstream.hh"

// SLIC
#include "MCParticleManager.hh"

// Geant4
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4ePairProduction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4LorentzVector.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4Positron.hh"
#include "Randomize.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "G4VEmModel.hh"
#include "G4VPhysicalVolume.hh"

// LCIO
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"

namespace slic {

// Start with L1 upstream only, same sensors used by
// ForcedPairConversionL1Upstream: for the top half the axial sensor is
// upstream, for the bottom half the stereo sensor is upstream (see that
// plugin's header comment for why these differ by half).
const std::vector<std::string> ForcedTridentL1Upstream::kTargetVolumes = {
    "module_L1t_halfmodule_axial_sensor",
    "module_L1b_halfmodule_stereo_sensor",
};

void ForcedTridentL1Upstream::beginRun(const G4Run *) {
  if (initialized_) {
    return;
  }

  // G4ePairProduction is not part of QGSP_BERT's default electron physics.
  // This plugin owns a private instance purely as a sampling engine --
  // CrossSectionPerVolume()/SelectModelForMaterial()/SampleSecondaries()
  // are all called manually from stepping() below, exactly as the existing
  // ForcedPairConversionL1Upstream plugin calls the real "conv" process.
  // Deliberately NOT attached via AddDiscreteProcess(): doing so makes
  // Geant4's own per-step dispatch (G4SteppingManager->
  // PostStepGetPhysicalInteractionLength->ComputeLambdaForScaledEnergy)
  // invoke this process on every electron step, which segfaults because
  // the lambda-table infrastructure that path depends on is never built
  // (only PreparePhysicsTable(), not BuildPhysicsTable(), is called here,
  // deliberately, to avoid touching shared dE/dx tables -- see below). No
  // shared physics-list code is touched either way.
  pairProd_ = new G4ePairProduction();

  // PreparePhysicsTable() (not BuildPhysicsTable(), which would rebuild
  // dE/dx and range tables for every registered electron energy-loss
  // process) is sufficient to construct and initialise the underlying
  // G4MuPairProductionModel so that CrossSectionPerVolume(),
  // SelectModelForMaterial() and SampleSecondaries() all work correctly
  // below -- CrossSectionPerVolume() is intended to fall back to a direct
  // model evaluation when no lambda table has been built. However,
  // PreparePhysicsTable() itself leaves theLambdaTable pointing at a
  // non-null but still-empty G4PhysicsTable (allocated by
  // G4PhysicsTableHelper::PreparePhysicsTable(), never populated since
  // BuildPhysicsTable() is deliberately never called), which otherwise
  // makes CrossSectionPerVolume() wrongly take the lambda-table lookup
  // branch and segfault dereferencing an empty table. Null it out so the
  // intended direct-model-evaluation fallback path is actually taken.
  pairProd_->PreparePhysicsTable(*G4Electron::Electron());
  pairProd_->SetLambdaTable(nullptr);

  initialized_ = true;
}

void ForcedTridentL1Upstream::stepping(const G4Step *step) {

  if (hasInteracted_) {
    // Primary beam electron has already been forced to interact this
    // event; everything else is left to propagate under normal physics.
    return;
  }

  G4Track *track = step->GetTrack();

  // Only the primary (generator-level) beam electron is ever a candidate
  // for forcing.
  if (track->GetDefinition() != G4Electron::Electron() || track->GetParentID() != 0) {
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

  const G4double kinEnergy = track->GetKineticEnergy();
  const G4double logKinEnergy = track->GetDynamicParticle()->GetLogKineticEnergy();
  const G4MaterialCutsCouple *couple = track->GetMaterialCutsCouple();

  // Local interaction probability for this slab: 1 - exp(-t/lambda), with
  // 1/lambda given by the real cross section per unit volume. This is the
  // weight the event is worth once we force the trident to happen here.
  // The 3-arg overload is called explicitly (mirroring
  // G4VEnergyLossProcess::MeanFreePath()'s own usage) because the 2-arg
  // overload, while declared in the installed Geant4 10.6.p01 header, has
  // no corresponding definition/export in the installed shared library.
  const G4double sigma = pairProd_->CrossSectionPerVolume(kinEnergy, couple, logKinEnergy);
  const G4double stepLength = step->GetStepLength();
  const G4double weight = 1.0 - std::exp(-sigma * stepLength);

  // Sample the e+/e- pair kinematics from the same model Geant4 would have
  // used. SampleSecondaries() only appends the two pair secondaries to the
  // output vector -- the degraded primary's kinematics are written by the
  // model to the process's own internal (protected) particle-change object,
  // which isn't reachable from here, so the scattered primary electron's
  // final 4-momentum is instead reconstructed below from 4-momentum
  // conservation using the two sampled secondaries.
  size_t idx = (size_t)couple->GetIndex();
  G4VEmModel *model = pairProd_->SelectModelForMaterial(kinEnergy, idx);

  std::vector<G4DynamicParticle *> secondaryParticles;
  model->SampleSecondaries(&secondaryParticles, couple, track->GetDynamicParticle());

  if (secondaryParticles.size() != 2) {
    // Nothing physical was sampled; leave the primary untouched.
    return;
  }

  const G4LorentzVector primaryIn4Momentum = track->GetDynamicParticle()->Get4Momentum();
  G4LorentzVector secondaries4Momentum(0., 0., 0., 0.);
  for (G4DynamicParticle *dp : secondaryParticles) {
    secondaries4Momentum += dp->Get4Momentum();
  }
  // Nuclear recoil carries away the (small) remaining momentum imbalance;
  // neglecting it here and assigning everything else to the scattered
  // primary electron is the same order of approximation already accepted
  // for the forced pair-conversion plugins' vertex placement.
  const G4LorentzVector scatteredElectron4Momentum = primaryIn4Momentum - secondaries4Momentum;

  G4DynamicParticle *scatteredElectron =
      new G4DynamicParticle(G4Electron::Electron(), scatteredElectron4Momentum);

  // Place the vertex at a random depth within the slab rather than always
  // at the exit face, so the resulting tracks can deposit charge in this
  // sensor too -- same approximation already accepted for the forced
  // pair-conversion plugins' "randpos" variant.
  const G4ThreeVector &entryPosition = step->GetPreStepPoint()->GetPosition();
  const G4ThreeVector &exitPosition = step->GetPostStepPoint()->GetPosition();
  const G4double entryTime = step->GetPreStepPoint()->GetGlobalTime();
  const G4double exitTime = step->GetPostStepPoint()->GetGlobalTime();

  const G4double u = G4UniformRand();
  const G4ThreeVector position = entryPosition + u * (exitPosition - entryPosition);
  const G4double globalTime = entryTime + u * (exitTime - entryTime);

  // Push onto the current step's own secondary bucket -- the native Geant4
  // hook for adding secondaries during a stepping action (see
  // ForcedPairConversionL1Upstream.cc for the detailed rationale: this is
  // what gets each forced track a real, non-zero track ID via
  // G4EventManager::StackTracks()).
  G4TrackVector *secondaryBucket = const_cast<G4Step *>(step)->GetfSecondary();

  G4Track *scatteredElectronTrack = new G4Track(scatteredElectron, globalTime, position);
  scatteredElectronTrack->SetParentID(track->GetTrackID());
  secondaryBucket->push_back(scatteredElectronTrack);

  for (G4DynamicParticle *dp : secondaryParticles) {
    G4Track *secondaryTrack = new G4Track(dp, globalTime, position);
    secondaryTrack->SetParentID(track->GetTrackID());
    secondaryBucket->push_back(secondaryTrack);
  }

  track->SetTrackStatus(fStopAndKill);

  hasInteracted_ = true;
  ++forcedCount_;

  // Stash the weight on the MCParticle collection under the key that
  // LcioManager::createLCEvent() reads at end-of-event. Unlike the
  // stdhep-driven ForcedPairConversion* plugins, this plugin's GPS-driven
  // gun never calls MCParticleManager::generateEvent() (GPS bypasses
  // MCParticleManager entirely), so the collection doesn't exist yet at
  // this point in the event and must be created here. This is safe because
  // TrackManager::saveTrackSummaries() now determines its geant4Generator
  // flag from the actual event source type (EventSourceManager), not from
  // whether the collection already exists -- see that file.
  MCParticleManager::instance()->createMCParticleCollection();
  LCCollectionVec *mcpColl = MCParticleManager::instance()->getMCParticleCollection();
  mcpColl->parameters().setValue("_weight", (float)weight);

  if (verbose_ > 1) {
    G4cout << "[ ForcedTridentL1Upstream ]: forced trident in " << volumeName << " with weight "
           << weight << G4endl;
  }
}

void ForcedTridentL1Upstream::endEvent(const G4Event *) {
  hasInteracted_ = false;
}

void ForcedTridentL1Upstream::endRun(const G4Run *) {
  G4cout << "[ ForcedTridentL1Upstream ]: forced tridents in " << forcedCount_ << " events"
         << G4endl;
}

} // namespace slic
DECLARE_PLUGIN(slic, ForcedTridentL1Upstream)
