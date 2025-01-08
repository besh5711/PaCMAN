The scripts contained in this folder enable the user to compare the extended 
Ptychographic Iterative Engine (ePIE) and PaCMAN (Partial Coherence and 
Monochromatization Algorithm with noise). The code was designed in MATLAB 2021b 
and also verified in 2023a, so for older (or newer) releases minor changes may be 
needed. Additionally, access to the following toolboxes are needed to run various
functions:

1. Image Processing Toolbox for the roipoly function (for easily defining parasitic
   scattering regions.
2. Parallel Computing Toolbox for the the parfor, parpool, and gpuArray 
   functions.
3. Statistics and Machine Learning Toolbox for easily adding Poisson noise with the
   poissrnd function.


If you use (Ms.) PaCMAN, please cite our paper "Robust Broadband Ptychography
Algorithms for High-Harmonic Soft X-Ray Supercontinua" (2025).

####################################################################################
SINGLE-WAVELENGTH

Below is a code tree for reference, in sequential order from top-to-bottom:

compareSingle   ----    generateData    ----    makeProbe
                                        ----    makeScanGrid
                                        ----    makeObject
                                        ----    makeSpectrum    
                                        ----    plotESWA
                                        ----    addNoise
                                        ----    generateC

                ----    ePIE            ----    iterationPlot

                ----    monochromatize  ----    BiCGSTAB
                                        ----    CGLS

                ----    PaCMAN          ----    iterationPlot

The main script is compareSingle. Here, you define all the variables. After defining
these variables, the script then runs the function generateData, which is responsible
for creating the measured diffraction patterns used in ePIE or PaCMAN. We begin by
making a Gaussian probe with a hard aperture. With the probe, we make a square scan
grid with random offsets. Then, we make the object from "cameraman" (modulus) and
"peppers" (phase). Next, we make a Gaussian spectrum, which is then used to generate
the transformation matrix C for monochromatization. We then compute the diffraction
patterns at each position and wavelength, and add the appropriate Gaussian/Poisson
noise and parasitic scattering background, which gives the final measured broadband
diffraction patterns.

Next, we run ePIE reconstructions. This function is entirely self-contained, with
the exception of an optional function iterationPlot which plots the object modulus
at each iteration to monitor reconstruction quality.

After ePIE, we monochromatize the measured diffraction patterns. We can solve
b = Cm using Conjugate Gradient Least Squares (CGLS) or Biconjugate Gradient 
Stabilized (BiCGSTAB). Soft thresholding or other sorts of denoising can help 
improve the accuracy of the monochromatization.

Finally, we run PaCMAN on the monochromatized diffraction patterns. Again,
iterationPlot can be used to monitor reconstruction quality.

####################################################################################
MULTI-WAVELENGTH (work in progress)

Next, we compare Ptychographic Information Multiplexing (PIM) and the Multi-spectral
Partial Coherence Mitigation Algorithm with Noise (Ms. PaCMAN).

Below is a code tree for reference, in sequential order from top-to-bottom:

compareAlgos    ----    generateData_multi    ----    makeProbe
                                              ----    makeScanGrid
                                              ----    makeObject
                                              ----    makeSpectrum
                                              ----    addNoise

                ----    PIM                   ----    iterationPlot

                ----    Ms_PaCMAN             ----    iterationPlot