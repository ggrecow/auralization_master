
# Framework for aircraft noise auralization 

This repository contains MATLAB scripts that compose a framework for aircraft noise auralization. 

## Overview

The **auralization framework** can be understood as a tool that transforms sound source descriptions into sound pressure signals at desired receiver positions. To achieve that, the following is done: 

- `Sound source`: predictions of an operational aircraft are provided. Therefore, the framework does not have any control over the models used to describe the sound source. It is tailored to work with input data from the well-established **Parametric Aircraft Noise Analysis Module** (PANAM) [1], an in-house software from the **German Aerospace Center** (DLR). PANAM provides inputs related to the aircraft trajectory and its associated performance and sound emissions. The sound emissions are provided in discrete time-steps, for which SPLs associated with the particular receiver positions are provided for different componential sound sources, namely:`Airframe noise`, `Engine noise`, `Fan tonal noise` and `Buzzsaw noise`. **Important:** sound emissions are provided already including Doppler effect. 

- `Sound synthesis`: sound source descriptions are transformed into sound pressure signals in time domain using different techniques from [2]. Signals containing tonal noise from `Fan tonal noise` and `Buzzsaw noise` partial sound sources are obtained using additive synthesis. Signals containing broadband noise from `Airframe noise`, `Engine noise` partial sound sources are obtained using granular synthesis. At the end of the synthesis procedure, sound signals from all partial sound sources are summed to obtain the overall aircraft noise.       

- `Propagation effects`: sound propagation effects through an inhomogeneous and moving atmosphere are obtained based on ray-tracing simulations in combination with different models. The
open-source **Atmospheric Ray Tracing (ART)** software [3] is used to simulate sound paths connecting the source and receiver positions, also known as eigenrays. Direct and 1st order reflected eigenrays are considered, for which inputs related to the aircraft trajectory and receiver positions, as well as atmospheric conditions (i.e. temperature, wind, relative humidity, and pressure) are necessary. The ART accounts for sound refraction and advection caused by sound speed gradients and wind, respectively. Based on the simulated eigenrays, and their associated propagation distance, the following propagation effects are applied: spreading loss, atmospheric absorption, and reflection factor. Finally, atmospheric transfer functions are
established by combining the aforementioned propagation effects and the phase-shift caused by the propagation times of the direct and reflected sound paths. Thereby, the overall sound propagation effect is described in frequency domain for each combination of source and receiver positions considered within a discrete flight trajectory.

- `Receiver model`: the simpliest case is that of an omnidirectional microphone (with 0 dB gain for all frequencies). This is the default case, for which monoaural diotic .WAV files are obtained at the receiver positions. A human receiver can also be considered using **Head-Related Transfer Functions (HRTF)** to account for sound propagation effects caused by the interaction of sound with human body parts. In this case, binaural sound signals are rendered using HRTFs from the FABIAN database. In the default case, binaural signals are rendered considering that the receiver is looking towards the sound source.  

**References**

[1] L. Bertsch, Noise Prediction within Conceptual Aircraft Design. Doctoral thesis, Technische Universität Braunschweig, 2013. DOI: [10.34912/n0is3-d3sign](https://doi.org/10.34912/n0is3-d3sign).

[2] M. P. Allen, Analysis and synthesis of aircraft engine fan noise for use in psychoacoustic studies. Master thesis, Virginia Polytechnic Institute and State University, 2012

[3] P. Schäfer and M. Vorländer, Atmospheric ray tracing: an efficient, open-source framework for finding eigenrays in a stratified, moving medium, Acta Acust., vol. 5, p. 26, 2021. DOI: [10.1051/aacus/2021018](https://doi.org/10.1051/aacus/2021018)

[4]  F. Brinkmann, A. Lindau, S. Weinzierl, G. Geissler, S. van de Par, M. Müller-Trapet, R. Opdam, and M. Vorländer, The FABIAN head-related transfer function data base, Data repository,” 2020. DOI: [10.14279/depositonce-5718.5](https://doi.org/10.14279/depositonce-5718.5).

# References

The auralization framework is comprehensively described and used in:

> G. Felix Greco, Sound quality analysis of virtual aircraft prototypes: framework development and application. Doctoral thesis, Technische Universität Braunschweig (to be published)

