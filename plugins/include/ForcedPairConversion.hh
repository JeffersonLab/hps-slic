#ifndef PLUGINS_FORCED_PAIR_CONVERSION_H
#define PLUGINS_FORCED_PAIR_CONVERSION_H

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
 * photon to pair-convert in a chosen tracker layer, rather than relying on
 * its natural (very small) conversion probability. A compensating event
 * weight equal to the real local conversion probability is written to the
 * MCParticle collection so the forced sample remains correctly normalized.
 *
 * Only the primary photon (track parent ID 0, no creator process) is ever
 * forced. Any photon produced downstream of the forced conversion (e.g. by
 * eBrem on the forced e+/e-) has a nonzero parent ID and is left alone, so
 * the conversion cannot cascade into later layers.
 */
class ForcedPairConversion : public UserActionPlugin {

public:
  /// Default constructor
  ForcedPairConversion() = default;

  /// Destructor
  ~ForcedPairConversion() = default;

  /// @return The name of this plugin
  std::string getName() { return "ForcedPairConversion"; }

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
  /// upstream (first, unshielded) sensor differs by half: for the top half
  /// the axial sensor is upstream, for the bottom half the stereo sensor is
  /// upstream (see HpsSiSensor's axial/stereo-vs-layer-parity convention).
  /// Each substring matches both "..._sensor_volume" (full silicon slab) and
  /// "..._sensor_active_volume" (its sensitive sub-volume), but not the
  /// enclosing non-sensitive "..._halfmodule_{axial,stereo}_volume" wrapper.
  static const std::vector<std::string> kTargetVolumes;

  /// Cached pointer to the gamma conversion ("conv") process, resolved lazily
  /// on first use since the process manager isn't available at construction.
  G4VEmProcess *convProcess_{nullptr};

  /// Whether the primary photon has already been forced to convert, or has
  /// already converted naturally, in this event.
  bool hasConverted_{false};

  /// Count of events in which a forced (non-natural) conversion occurred,
  /// for the end-of-run summary.
  long forcedCount_{0};

}; // ForcedPairConversion
} // namespace slic

#endif // PLUGINS_FORCED_PAIR_CONVERSION_H
