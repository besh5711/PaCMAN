MONOCHROMATIZATION

The scripts contained in this folder enable the user to compare the extended 
Ptychographic Iterative Engine (ePIE) and PaCMAN (Partial Coherence and 
Monochromatization Algorithm with noise). The code was designed in MATLAB 2021b 
and also verified in 2023a, so for older (or newer) releases minor changes may be 
needed. Additionally, access to the Image Processing Toolbox (for the roipoly 
function) and Parallel Computing Toolbox (for the parfor, parpool, and gpuArray 
functions) are required.

Below is a code tree for reference, in sequential order from top-to-bottom:

compareSingle   ----    generateData    ----    makeProbe
                                        ----    makeScanGrid
                                        ----    makeObject
                                        ----    makeSpectrum        ----    generateC
                                        ----    plotESWA
                                        ----    addNoise

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
MULTI-WAVELENGTH

Next, we compare Ptychographic Information Multiplexing (PIM) and the Multi-spectral
Partial Coherence Mitigation Algorithm with Noise (Ms. PaCMAN).

Below is a code tree for reference, in sequential order from top-to-bottom:

compareAlgos    ----    generateData    ----    makeProbe
                                        ----    makeScanGrid
                                        ----    makeObject
                                        ----    makeSpectrum        ----    generateC
                                        ----    addNoise
                                        ----    addBackground

                ----    ePIE            ----    iterationPlot

                ----    PaCMAN          ----    iterationPlot