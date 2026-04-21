# Neural_electrophysiology

MATLAB analysis code for the B.A. thesis *Striking the Right Note: The Role of the Dorsomedial Striatum in Learning Sensorimotor Associations Underlying Sequential Behavior* (Wekwejt, Harvard College, 2026).

Many natural behaviors, from pressing piano keys to producing speech, rely on flexibly combining discrete motor actions into cue-guided sequences. This study asks whether the dorsomedial striatum (DMS) is necessary for learning such sensorimotor associations, and how DMS population activity evolves with training. Rats were trained on a three-lever "piano" task in which LED cues instructed lever presses, and DMS contributions were probed with excitotoxic lesions and chronic in vivo tetrode recordings.

DMS-lesioned animals failed to acquire cue-action mappings despite preserved response latency, establishing DMS necessity for associative learning rather than basic movement execution. In recordings, aggregate population magnitude and variability remained stable across training, while within-session temporal reliability and across-session population similarity increased with learning. Event-kernel modulation reweighted toward action and outcome epochs, and cue-type decoding became preferentially aligned to the impending action during late learning.

## Layout

- `lesions/`: DMS-lesion behavioral analysis (Fig 3)
- `ephys/`: chronic recordings, organized by figure:
  - `00_preprocessing/`: trial construction, behavior, Weibull fit
  - `fig02_unit_classification/`: MSN/FSI waveform classification
  - `fig04_population_magnitude/`: whole-trial FR and variability
  - `fig05_temporal_reliability/`: split-half reliability
  - `fig06_session_reorganization/`: cross-session similarity
  - `fig07_glm_task_epochs/`: event-kernel GLM
  - `fig08_epoch_modulation/`: learning-related epoch changes
  - `fig09_decoding_learning/`: cue-type decoding across training
  - `fig10_cue_response_dissociation/`: press-decoding on incorrect trials

Requires MATLAB R2024b. MIT licensed.
