This code package enables the user to compare the extended Ptychographic
Iterative Engine (ePIE) to the Partial Coherence and Monochromatization
Algorithm with Noise (PaCMAN), which is done in the script compareSingle.
We can also compare multi-wavelength algorithms that have wavelength-
dependent object/probe, namely Ptychographic Information Multiplexing (PIM)
and the Multi-spectral Partial Coherence Mitigation Algorithm with Noise 
(Ms. PaCMAN). This is done in the script compareMulti.

The code was designed in MATLAB 2021b and also verified in 2023a, so for 
older (or newer) releases minor changes may be needed. Additionally, access
to the following toolboxes are needed to run various functions:

1. Image Processing Toolbox for the roipoly function (for easily defining
   parasitic scattering regions.
2. Parallel Computing Toolbox for the the parfor, parpool, and gpuArray 
   functions.
3. Statistics and Machine Learning Toolbox for easily adding Poisson noise 
   with the poissrnd function.

Finally, if you use (Ms.) PaCMAN, please reference the following work:

Benjamin Shearer, Henry Kapteyn, Iona Binnie, Nicholas W. Jenkins, and 
Margaret Murnane, "Robust broadband ptychography algorithms for high-
harmonic soft X-ray supercontinua," Opt. Express 33, 717-735 (2025).

If you have questions about the code or spot an error, feel free to email 
me at besh5711@colorado.edu.

###########################################################################
SINGLE-WAVELENGTH

Below is a code tree for reference, in sequential order from top-to-bottom:

compareSingle   ----    generateData    ----    makeSpectrum
                                        ----    makeProbe
                                        ----    makeScanGrid
                                        ----    makeObject
                                        ----        
                                        ----    plotESWA
                                        ----    addNoise
                                        ----    generateC

                ----    ePIE            ----    iterationPlot

                ----    monochromatize  ----    BiCGSTAB
                                        ----    CGLS

                ----    PaCMAN          ----    iterationPlot

                ----    NRMSE

The main script is compareSingle. Here, you define all the variables. After 
defining these variables, the script then runs the function generateData, 
which is responsible for creating the measured diffraction patterns used in 
ePIE or PaCMAN. We begin by making a Gaussian probe with a hard aperture. 
With the probe, we make a square scan grid with random offsets. Then, we 
make the object from "cameraman" (modulus) and "peppers" (phase). Next, we 
make a Gaussian spectrum, which is then used to generate the transformation 
matrix C for monochromatization. We then compute the diffraction patterns 
at each position and wavelength, and add the appropriate Gaussian/Poisson
noise and parasitic scattering background, which gives the final measured 
broadband diffraction patterns.

Next, we run ePIE reconstructions. This function is entirely self-contained, 
with the exception of an optional function iterationPlot which plots the 
object modulus at each iteration to monitor reconstruction quality.

After ePIE, we monochromatize the measured diffraction patterns. We can 
solve b = Cm using Conjugate Gradient Least Squares (CGLS) or Biconjugate 
Gradient Stabilized (BiCGSTAB). Soft thresholding or other sorts of 
denoising can help improve the accuracy of the monochromatization.

Finally, we run PaCMAN on the monochromatized diffraction patterns. Again,
iterationPlot can be used to monitor reconstruction quality.

At the end of the script, we calculate the normalized root-mean-square 
error (NRMSE) for each of the three reconstructions and plot them.

###########################################################################
MULTI-WAVELENGTH

Next, we compare Ptychographic Information Multiplexing (PIM) and the 
Multi-spectral Partial Coherence Mitigation Algorithm with Noise (Ms. 
PaCMAN).

Below is a code tree for reference, in sequential order from top-to-bottom:

compareMulti    ----    generateData_multi    ----    makeSpectrum
                                              ----    makeProbe
                                              ----    makeScanGrid
                                              ----    makeObject
                                              ----    addNoise

                ----    PIM                   ----    iterationPlot_multi

                ----    Ms_PaCMAN             ----    iterationPlot_multi

                ----    NRMSE

The script compareMulti is very similar to compareSingle. The main 
difference is that the probes and objects are now functions of wavelength.
When we run generateData_multi, I assume that the object does not change as
a function of wavelength, which means that the exact objects and probes are 
simply scaled (e.g. interpolated to a different grid size) versions of the 
center wavelength. Of course this doesn't have to be true, which is the
primary reason why we would use PIM/Ms. PaCMAN over ePIE/PaCMAN. When we 
set the number of wavelengths used to generate the data equal to the number 
of wavelengths used to reconstruct the objects/probes, we get the best 
convergence, as we'd expect (e.g. multi-wavelength algorithms are excellent
for comb-like spectra).


