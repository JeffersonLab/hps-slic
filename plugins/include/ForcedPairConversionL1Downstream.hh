#ifndef PLUGINS_FORCED_PAIR_CONVERSION_L1_DOWNSTREAM_H
#define PLUGINS_FORCED_PAIR_CONVERSION_L1_DOWNSTREAM_H

//~~~ SLIC ~~~//
#include "UserActionPlugin.hh"

// STL
#include <string>
#include <vector>

// Forward declarations
class G4Step;
class G4VEmProcess;

namespace slic {

/**
 * Class defining a UserAction that forces the primary (generator-level) WAB
 * photon to pair-convert in the chosen tracker layer's downstream (second,
 * shielded) sensor, rather than relying on its natural (very small)
 * conversion probability. A compensating event weight equal to the real
 * local conversion probability is written to the MCParticle collection so
 * the forced sample remains correctly normalized.
 *
 * This is the downstream-sensor counterpart to ForcedPairConversion, which
 * targets the upstream (first, unshielded) sensor of the same Layer 1 pair.
 * The two plugins are otherwise identical.
 *
 * Only the primary photon (track parent ID 0, no creator process) is ever
 * forced. Any photon produced downstream of the forced conversion (e.g. by
 * eBrem on the forced e+/e-) has a nonzero parent ID and is left alone, so
 * the conversion cannot cascade into later layers.
 *
 * Natural conversion is also temporarily deactivated for the gamma "conv"
 * process while the primary photon transits the upstream (first) sensor of
 * the pair, so that every event reaches the downstream target instead of
 * some fraction being "wasted" on a natural upstream conversion. The
 * deactivation is global to all gammas (G4's process activation flag is per
 * particle-type, not per-track), but is only held for the brief window the
 * primary photon is inside the upstream sensor, then restored.
 */
class ForcedPairConversionL1Downstream : public UserActionPlugin {

public:
  /// Default constructor
  ForcedPairConversionL1Downstream() = default;

  /// Destructor
  ~ForcedPairConversionL1Downstream() = default;

  /// @return The name of this plugin
  std::string getName() { return "ForcedPairConversionL1Downstream"; }

  /// @return True if this plugin has a stepping action
  bool hasSteppingAction() final override { return true; }

  /// @return True if this plugin has an event action
  bool hasEventAction() final override { return true; }

  /// @return True if this plugin has a run action
  bool hasRunAction() final override { return true; }

  /**
   * Stepping action: watches the primary photon, and once it reaches the
   * exit boundary of the target tracker layer without having converted
   * naturally, forces the conversion right there.
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
  /**
   * Find and cache the gamma conversion ("conv") process registered on the
   * photon's process manager.
   */
  G4VEmProcess *findConversionProcess(const G4Step *step);

  /// Substrings identifying the individual Layer 1 sensor volumes where the
  /// conversion is forced. Layer 1's axial and stereo sensors are separate,
  /// non-stacked sibling volumes (not two slabs of one "module"), so the
  /// downstream (second, shielded) sensor differs by half: for the top half
  /// the stereo sensor is downstream, for the bottom half the axial sensor
  /// is downstream (see HpsSiSensor's axial/stereo-vs-layer-parity
  /// convention) -- the mirror image of ForcedPairConversion's upstream
  /// targeting. Each substring matches both "..._sensor_volume" (full
  /// silicon slab) and "..._sensor_active_volume" (its sensitive
  /// sub-volume), but not the enclosing non-sensitive
  /// "..._halfmodule_{axial,stereo}_volume" wrapper.
  static const std::vector<std::string> kTargetVolumes;

  /// Substrings identifying the upstream (first, unshielded) sensor of the
  /// same Layer 1 pair -- i.e. ForcedPairConversion::kTargetVolumes. Natural
  /// conversion is deactivated while the primary photon transits these
  /// volumes, so it always survives to reach the downstream target below.
  static const std::vector<std::string> kUpstreamVolumes;

  /// Cached pointer to the gamma conversion ("conv") process, resolved lazily
  /// on first use since the process manager isn't available at construction.
  G4VEmProcess *convProcess_{nullptr};

  /// Whether the primary photon has already been forced to convert, or has
  /// already converted naturally, in this event.
  bool hasConverted_{false};

  /// Whether the "conv" process is currently deactivated because the primary
  /// photon is inside the upstream sensor right now. Tracked explicitly
  /// rather than queried back from the process (G4VProcess has no public
  /// getter for its own activation flag) so it can be restored exactly once,
  /// on the step where the primary photon leaves the upstream volume.
  bool convDeactivated_{false};

  /// Count of events in which a forced (non-natural) conversion occurred,
  /// for the end-of-run summary.
  long forcedCount_{0};

}; // ForcedPairConversionL1Downstream
} // namespace slic

#endif // PLUGINS_FORCED_PAIR_CONVERSION_L1_DOWNSTREAM_H
