# Introduction - UTIL_PION

- The structure of this REPO has been changed to match that of UTIL_KAONLT and UTIL_PROTON

  - scripts/        : contains all analysis scripts
  - batch/          : for running batch jobs to farm
  - archive/        : archive files that may be useful in the future, old TProof scripts and online running scripts
  - config/         : replay configuration files (e.g. DEF-files, PARAM)

- Please contact Stephen Kay (stephen.kay@uregina.ca) or Richard Trotta (trotta@cua.edu) or Muhammad Junaid (mjo147@uregina.ca) for more information.

# Initial Setup

- Before running ensure you are using ROOT 6.18.04 -

  - source /apps/root/6.18.04/setroot_CUE.csh 

  - This assumes you're running this on the JLab iFarm

- You should also make sure you have the relevant packages, in particular, if you have not done so previously.
, execute -

  - pip install --user root_numpy --user

  - pip install --user root_pandas --user

- Before running any scripts, execute the following -

  - cp -r bin/python/ltsep ~/.local/lib/python3.4/site-packages/

  - If you're not on the farm, copy the ltsep package to wherever your local packages for python 3.4 are.

  - After you copy this package into place, you can run the sym link setup script.

- Sym links are required for many of the analysis scripts to function. Check the pathing called in the analysis scripts carefully and examine the sym link script in this directory for more information -

  - UTIL_SymLinkSetup.sh

# UTIL_PION - Scientific Motivation

The pion occupies a special place in nature as one of the lightest hadrons, with one valence quark, and one valence antiquark. Small as it might be, the pion is also responsible for the long range character of the strong interaction that binds the atomic nucleus together. If chiral symmetry, which implies that Dirac Fermions are massless, were an exact global symmetry of strong interactions then pions would be massless. Through gluon-quark interaction and by explicit inclusion of light quark masses, chiral symmetry of massless QCD undergoes explicit symmetry breaking, thus giving the pion its mass. This puts the pion at the core of the mechanism that dynamically generates all the mass of the hadrons and makes it a crucial element in understanding hadron structure. Globally, this experiment will aim to confirm the potential of pion measurements both for studies of the pion structure itself and of the 3D structure of the proton, in terms of spatial imaging (tomography). In particular, E12-07-105 will probe if the measurements to map the spatial extension of the charged pion can be utilized to enable 3D spatial tomography of light quarks.

The E12-07-105 experiment is an exclusive measurement of the L/T separated pion electroproduction cross section, it aims to probe conditions for factorization of deep exclusive measurements for charged pions in GPD studies and the pion form factor. It will make precision measurements of the L/T separated pion electroproduction cross sections to the highest achievable value of Q2 at the 12 GeV Jefferson Lab, ~10 GeV2. Before considering the extraction of GPDs from pion electroproduction data, we must confirm that the charged pion longitudinal cross section is more dominant than the transverse cross section as expected from GPD models. We will test the longitudinal cross section dominance in charged pion production by conducting a measurement of the Q2 dependence of the
π+ longitudinal and transverse cross sections. A significant longitudinal response may be indicative of the realisation of the scaling expectation of the GPD formalism for charged pion electroproduction. The results may be useful in constraining the longitudinal backgrounds in the extraction of the pion form factor from pion electroproduction data. We will measure forward π+ electroproduction by detecting the the pion produced in coincidence with the scattered electron using the Hall C High Momentum Spectrometer (HMS) and Super High Momentum Spectrometer (SHMS). We will extract the cross sections using the Rosenbluth technique. In additon, measurements in non-parallel kinematics will allow for simultaneous extraction of the interference terms and serve as a test of t-dependence of the π+ cross section.
