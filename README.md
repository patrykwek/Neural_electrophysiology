# Neural_electrophysiology

MATLAB analysis code for the B.A. thesis *Striking the Right Note: The Role of the Dorsomedial Striatum in Learning Sensorimotor Associations Underlying Sequential Behavior* (Wekwejt, Harvard College, 2026).

Rats performed a cue-guided three-lever "piano" task. Analyses combine targeted DMS lesion behavior with chronic in vivo tetrode recordings in the DMS.

## Layout

- `lesions/` — DMS-lesion behavioral analysis (Fig 3)
- `ephys/` — chronic recordings, organized by figure:
  - `00_preprocessing/` — trial construction, behavior, Weibull fit
  - `fig02_unit_classification/` — MSN/FSI waveform classification
  - `fig04_population_magnitude/` — whole-trial FR and variability
  - `fig05_temporal_reliability/` — split-half reliability
  - `fig06_session_reorganization/` — cross-session similarity
  - `fig07_glm_task_epochs/` — event-kernel GLM
  - `fig08_epoch_modulation/` — learning-related epoch changes
  - `fig09_decoding_learning/` — cue-type decoding across training
  - `fig10_cue_response_dissociation/` — press-decoding on incorrect trials

Requires MATLAB R2024b. MIT licensed.
