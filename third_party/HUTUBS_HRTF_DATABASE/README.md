# HUTUBS HRTF database (optional, not shipped)

This folder holds the HUTUBS HRTF database, which is only needed for binaural rendering with
`hrtf_database = HUTUBS` in the `.ini` file. The database is **not** distributed with the framework.
The default HRTF set (`hrtf_database = FABIAN`) is shipped in `third_party/FABIAN_HRTF_DATABASE_v4/`.

## Installation

1. Download `HRIRs.zip` from the HUTUBS database: https://doi.org/10.14279/depositonce-8487
2. Extract it into this folder, so that the files are found as, e.g.,
   `third_party/HUTUBS_HRTF_DATABASE/HRIRs/pp26_SHcoefficients_measured.mat`
3. Run `startup_auralization.m` again (it adds this folder and its subfolders to the MATLAB path).

Only the files `ppN_SHcoefficients_measured.mat` (about 1.2 MB each, N = 1 ... 96) are used. The
other files (`*_HRIRs_*.sofa`, `*_SHcoefficients_simulated.mat`, `Figures/`) can be deleted, or
only the subjects needed can be kept. The contents of this folder are ignored by git (except this file).

## Notes

- Subjects 1 and 96 are repeated measurements of the FABIAN head-and-torso simulator; subjects 22
  and 88 are repeated measurements of the same human subject. HUTUBS pp1/pp96 are not identical to
  `hrtf_database = FABIAN` (different measurement system), see the `.ini` files and the references below.
- The HUTUBS HRIRs are obtained with `AKhrirDatabase` (AKtools), which needs `AKhrirDatabase.m`,
  `AKfractionalDelayCyclic.m` and `AKkaiser.m` in `third_party/AKtools/`.

## License and references

The HUTUBS database is licensed under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/). If you
use it, cite:

- F. Brinkmann, M. Dinakaran, R. Pelzer, P. Grosche, D. Voss, and S. Weinzierl, A cross-evaluated
  database of measured and simulated HRTFs including 3D head meshes, anthropometric features, and
  headphone impulse responses, J. Audio Eng. Soc., vol. 67, no. 9, pp. 705–718, 2019.
  DOI: [10.17743/jaes.2019.0024](https://doi.org/10.17743/jaes.2019.0024)
- F. Brinkmann, M. Dinakaran, R. Pelzer, J. J. Wohlgemuth, F. Seipel, and S. Weinzierl, The HUTUBS
  head-related transfer function (HRTF) database, Data repository, 2019.
  DOI: [10.14279/depositonce-8487](https://doi.org/10.14279/depositonce-8487)
