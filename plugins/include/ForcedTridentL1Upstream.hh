#ifndef PLUGINS_FORCED_TRIDENT_L1_UPSTREAM_H
#define PLUGINS_FORCED_TRIDENT_L1_UPSTREAM_H

//~~~ SLIC ~~~//
#include "UserActionPlugin.hh"

// STL
#include <string>
#include <vector>

// Forward declarations
class G4Step;
class G4Run;
class G4ePairProduction;

namespace slic {

/**
 * Class defining a UserAction that forces a primary (generator-level) beam
 * electron to undergo a genuine one-step trident interaction,
 * e- + N -> e- + e+ + e- + N, in a chosen tracker sensor, via Geant4's
 * G4ePairProduction process (direct nuclear-field e+e- pair production off
 * a charged lepton). This process scales linearly with target thickness and
 * dominates the two-step "pseudo-trident" (eBrem followed by pair
 * conversion, which scales quadratically) in thin silicon sensor layers.
 *
 * G4ePairProduction is not part of the default QGSP_BERT electron physics.
 * This plugin constructs and initializes (via PreparePhysicsTable()) its
 * own private instance of the process at beginRun(), purely as a sampling
 * engine for manual CrossSectionPerVolume()/SampleSecondaries() calls in
 * stepping() -- it is never attached to G4Electron::Electron()'s process
 * manager, so Geant4's own per-step dispatch never invokes it (that path
 * requires lambda tables this plugin never builds). No shared physics-list
 * code is touched.
 *
 * A compensating event weight equal to the real local interaction
 * probability for the traversed slab is written to the MCParticle
 * collection so the forced sample remains correctly normalized. For this
 * plugin's GPS-driven gun, correct weight propagation also required a
 * small fix in TrackManager::saveTrackSummaries() (see that file) so that
 * pre-creating the MCParticle collection here doesn't defeat its
 * geant4Generator detection.
 *
 * Only the primary beam electron (track parent ID 0) is ever forced. The
 * forced interaction kills the primary and injects three new tracks (the
 * degraded scattered electron and the sampled e+/e- pair), all with a
 * nonzero parent ID, so the forcing can never cascade.
 */
class ForcedTridentL1Upstream : public UserActionPlugin {

public:
  /// Default constructor
  ForcedTridentL1Upstream() = default;

  /// Destructor
  ~ForcedTridentL1Upstream() = default;

  /// @return The name of this plugin
  std::string getName() { return "ForcedTridentL1Upstream"; }

  /// @return True if this plugin has a stepping action
  bool hasSteppingAction() final override { return true; }

  /// @return True if this plugin has an event action
  bool hasEventAction() final override { return true; }

  /// @return True if this plugin has a run action
  bool hasRunAction() final override { return true; }

  /**
   * Run action: construct and initialize this plugin's private
   * G4ePairProduction instance on G4Electron::Electron(). Guarded so it
   * only runs once even across multiple runs in one job.
   *
   * @param[in] run Geant4 run object.
   */
  void beginRun(const G4Run *run) final override;

  /**
   * Stepping action: watches the primary beam electron, and once it is
   * inside the target tracker sensor, forces the trident interaction.
   *
   * @param[in] step Geant4 step
   */
  void stepping(const G4Step *step) final override;

  /**
   * Reset the per-event forced flag.
   *
   * @param event Geant4 event object.
   */
  void endEvent(const G4Event *) final override;

  /**
   * Print a summary of how many events were forced over the run.
   *
   * @param run Geant4 run object.
   */
  void endRun(const G4Run *run) final override;

private:
  /// Substrings identifying the individual Layer 1 sensor volumes where the
  /// trident is forced -- same upstream (first, unshielded) sensors used by
  /// ForcedPairConversionL1Upstream, kept for consistency.
  static const std::vector<std::string> kTargetVolumes;

  /// This plugin's own, privately-owned G4ePairProduction process instance.
  /// Constructed and attached to G4Electron::Electron() in beginRun().
  G4ePairProduction *pairProd_{nullptr};

  /// Whether beginRun() has already run its one-time setup.
  bool initialized_{false};

  /// Whether the primary beam electron has already been forced to
  /// interact in this event.
  bool hasInteracted_{false};

  /// Count of events in which a forced trident occurred, for the
  /// end-of-run summary.
  long forcedCount_{0};

}; // ForcedTridentL1Upstream
} // namespace slic

#endif // PLUGINS_FORCED_TRIDENT_L1_UPSTREAM_H
