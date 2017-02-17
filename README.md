# halla-xem-analysis

This is the single-arm version of SIMC for the HRS Spectrometers.

Author: Barak Schmookler

Contents
--------

* Generates events uniformly in delta, phi (y'), theta (x'), and along the target. Input data can be taken from text file instead.
* There are four versions. The "phase_space_hrs" subdirectory contains the 6GeV spectrometer models. The "phase_space_sos" subdirectory also contains the 6GeV spectrometer models, but with the SOS quad replacing Q1. The "phase_space_spring16" and "phase_space_fall16" subdirectory contain models of the spectrometers as they were in Spring and Fall 2016, respectively.

Additional Information
----------------------

* Some information can be found [here](https://hallaweb.jlab.org/wiki/index.php/Simulation_using_SIMC).
* Currently, the 'carbon_foil.inp' file is the only up-to-date input file in each subdirectory.
* A paw ntuple output file is created in the "worksim" subdirectory. It can be converted to a ROOT file using the h2root utility.
* I plan on addition some more info on what is written out to the output file at a later date.
