#ifndef PLUGINS_FORCED_PAIR_CONVERSION_L2_UPSTREAM_H
#define PLUGINS_FORCED_PAIR_CONVERSION_L2_UPSTREAM_H

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
 * photon to pair-convert in the upstream (first) sensor of Layer~2, rather
 * than relying on its natural (very small) conversion probability. A
 * compensating event weight equal to the real local conversion probability
 * is written to the MCParticle collection so the forced sample remains
 * correctly normalized.
 *
 * This is the Layer-2 counterpart to ForcedPairConversion (Layer-1 upstream)
 * and ForcedPairConversionL1Downstream (Layer-1 downstream). The three
 * plugins are otherwise identical in structure.
 *
 * Only the primary photon (track parent ID 0, no creator process) is ever
 * forced. Any photon produced downstream of the forced conversion (e.g. by
 * eBrem on the forced e+/e-) has a nonzero parent ID and is left alone, so
 * the conversion cannot cascade into later layers.
 *
 * Natural conversion is also temporarily deactivated for the gamma "conv"
 * process while the primary photon transits any of Layer~1's four sensor
 * volumes (both the upstream and downstream sensor of each half), so that
 * every event's photon survives the whole of Layer~1 intact and reaches the
 * Layer~2 target instead of some fraction being "wasted" on a natural
 * Layer-1 conversion. The deactivation is global to all gammas (G4's process
 * activation flag is per particle-type, not per-track), but is only held for
 * the brief window the primary photon is inside one of those four volumes,
 * then restored.
 */
class ForcedPairConversionL2Upstream : public UserActionPlugin {

public:
  /// Default constructor
  ForcedPairConversionL2Upstream() = default;

  /// Destructor
  ~ForcedPairConversionL2Upstream() = default;

  /// @return The name of this plugin
  std::string getName() { return "ForcedPairConversionL2Upstream"; }

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

  /// Substrings identifying the individual Layer 2 sensor volumes where the
  /// conversion is forced. Layer 2's axial and stereo sensors are separate,
  /// non-stacked sibling volumes (not two slabs of one "module"), so the
  /// upstream (first) sensor differs by half: for the top half the axial
  /// sensor is upstream (layer index 3 in the LCDD physvolid table), for the
  /// bottom half the stereo sensor is upstream -- the same axial/stereo-vs-
  /// layer-parity convention as Layer 1 (see HpsSiSensor). Each substring
  /// matches both "..._sensor_volume" (full silicon slab) and
  /// "..._sensor_active_volume" (its sensitive sub-volume), but not the
  /// enclosing non-sensitive "..._halfmodule_{axial,stereo}_volume" wrapper.
  static const std::vector<std::string> kTargetVolumes;

  /// Substrings identifying all four Layer 1 sensor volumes (both the
  /// upstream and downstream sensor of each half) -- the union of
  /// ForcedPairConversion::kTargetVolumes and
  /// ForcedPairConversionL1Downstream::kTargetVolumes. Natural conversion is
  /// deactivated while the primary photon transits any of these volumes, so
  /// it always survives the whole of Layer 1 to reach the Layer 2 target
  /// below.
  static const std::vector<std::string> kUpstreamVolumes;

  /// Cached pointer to the gamma conversion ("conv") process, resolved lazily
  /// on first use since the process manager isn't available at construction.
  G4VEmProcess *convProcess_{nullptr};

  /// Whether the primary photon has already been forced to convert, or has
  /// already converted naturally, in this event.
  bool hasConverted_{false};

  /// Whether the "conv" process is currently deactivated because the primary
  /// photon is inside one of the Layer 1 sensors right now. Tracked
  /// explicitly rather than queried back from the process (G4VProcess has no
  /// public getter for its own activation flag) so it can be restored
  /// exactly once, on the step where the primary photon leaves Layer 1.
  bool convDeactivated_{false};

  /// Count of events in which a forced (non-natural) conversion occurred,
  /// for the end-of-run summary.
  long forcedCount_{0};

}; // ForcedPairConversionL2Upstream
} // namespace slic

#endif // PLUGINS_FORCED_PAIR_CONVERSION_L2_UPSTREAM_H
