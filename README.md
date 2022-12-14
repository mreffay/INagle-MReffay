Application tutorial
IrŠne Nagle ? Myriam Reffay

 Installation 
       TensioX can be installed and downloaded from http://www.msc.univ-paris-diderot.fr/spip.php?rubrique274&lang=en. The application is coded in MATLAB and the source code can be obtained at the same location. 
       No MATLAB license is required to run the application but MATLAB Runtime is required and is installed simultaneously with the application. 

 Introduction and principle 
       This application is designed to measure the surface tension and the Young modulus of magnetic multicellular spheroids by flattening them with an external magnet as described in [1]. 
       Briefly, iron oxide superparamagnetic nanoparticles (?-Fe2O3) obtained via Massart?s procedure are incorporated into the cells then spheroids are formed by magnetic molding (overnight formation). Finally, spheroids are flattened with an external magnet and the side profile of each spheroid is monitored with a camera. The equilibrium shape of the spheroid is determined by the competition between surface tension and magnetic forces. Surface tension is measured by fitting the spheroid side profile at equilibrium while Young modulus is determined thanks to the radius of the contact zone using Hertz theory. 
    To obtain the measurements, four elements are needed:
 A picture of the spheroid at t_0 (the initial spheroid needs to be spherical to ensure the accuracy of the measurements).
 A picture of the spheroid at t_f when the equilibrium shape of the spheroid under flattening is reached.
 The scale factor of the imaging system in m/pixels.
 The magnetic force per unit of volume applied on the spheroid f_v=M_vžgrad(B), with M_v the magnetic moment per unit of volume (that can be measured using VSM, Vibrating Sample Magnetometry, measurements for example) and grad(B) the magnetic field gradient of the external magnet at the position of the spheroid. 

 Step-by-step tutorial

 Open TensioX.

 Select the spheroid image at t=t_0.

 Then select image of spheroid at t=t_f when the equilibrium shape is reached.

 Select folder where results will be saved.

 Spheroid image at t=t_0 opens and several windows open asking to select successively the left side, the right side, the apex and the bottom of the spheroid. 

For each point, press OK then click on the left, right, apex or bottom of the initial spheroid (as shown on the below figure) using the crosshair cursor, then press ENTER.

From this, the initial volume of the spheroid is estimated


 Spheroid image at t=t_f opens and several windows open asking to select successively the left side, the right side, the apex and the bottom of the spheroid.

For each point, press OK then click on the left, right, apex or bottom of the flattened spheroid (as shown on the below figure) using the crosshair cursor, then press ENTER. 

From this, the theoretical profile is optimized by minimizing the width the height and the volume of the spheroid, to obtain the surface tension of the spheroid. 

 Finally, two windows open asking to select the right and left limit of the contact area of the flattened spheroid.

As previously, for each point, press OK then click on the left or right of the contact area of the flattened spheroid (as shown on the below figure) using the crosshair cursor, then press ENTER. 









 Enter the image scale in m/pixel, the measured value of the magnetic force per unit of volume in N/m3 and an estimated value of gamma in mN/m. If the estimated value of gamma is unknown, leave the default value. Then press OK. 








 Two windows open containing the computed surface tension and Young modulus values and an overlay of the flattened spheroid image and the fitted profile for the corresponding surface tension. A good agreement between the experimental and the computed profile (in red) is required to guarantee the reliability of the surface tension results. 






Both files are saved in the previously selected folder. The computed surface tension and the Young modulus are saved in a results.txt file.














 Detailed description of the MATLAB code 
    Surface tension and Young modulus are determined as described in [2]. The following paragraphs describe what the MATLAB code contains.  
 Surface tension measurement
    Theoretical profiles are obtained by resolving numerically (ode45 function) the classical Laplace equation of capillarity describing the mechanical equilibrium conditions for two homogeneous fluids separated by an interface and in non-wetting conditions. To adjust the numerical profile to the experimental profile, the quadratic error e on the height h, the width w and the volume V of the flattened spheroid is minimized with respect to the curvature at apex of the flattened spheroid (b) and  c=  f_v/?  the capillary constant of the system.
e=(1-h_th/h_exp )^2+(1-w_th/w_exp )^2+(1-V_th/V_exp )^2
h_exp, w_exp are computed from the extracted values on the flattened spheroid image while V_exp is computed from the extracted values on the initial spheroid image (volume of a sphere). 
    The initial parameters for the minimization are taken such as b_0=1/R_0  with 2R_0 the average between the width and the height of the initial spheroid and c_0=  f_v/?_0   with f_v and ?_0 the magnetic force per unit of volume and the estimated surface tension respectively, entered manually. The quadratic error e is minimized with the function fminsearch, to obtain the capillary constant of the system from which the surface tension ? is deduced. The final numerical profile and the image of the flattened spheroid are superimposed to check for the reliability of the result.

 Young modulus measurement 
    The Young modulus E is computed using Hertz theory for an elastic sphere with an initial radius R_0 which gives E=((1-?^2 )  ? M_(v ) grad(B)  ?R_0?^4)/L     where ? stands for the Poisson ratio (? = 1/2) and L for the radius of the contact zone computed from the extracted values of the contact area of the flattened spheroid image. 

 Cautions and remarks 

 The initial shape of the spheroid at t_0 has to be spherical to give reliable measurements. 
 The substrate on which the spheroid is has to be flat and the camera has to be correctly aligned to provide accurate measurements. 
 This application can be adapted to other systems than multicellular aggregates such as any type of viscoelastic fluid or material in non-wetting conditions with respect to the substrate.   

 Authors, copyright, distribution policy
This application was written by Irène Nagle. 
Copyright ¸ 2021 Irène Nagle and Myriam Reffay, University of Paris. 
This application is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 or any later version.
This application is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this application.  If not, see http://www.gnu.org/licenses/.
Contact : 
IrŠne Nagle, irene.nagle@u-paris.fr, Myriam Reffay, myriam.reffay@u-paris.fr 

References: 
[1]      Irène Nagle, Florence Delort, Sylvie H‚non, Claire Wilhelm, Sabrina Batonnet-Pichon, Myriam Reffay, The importance of intermediate filaments in the shape maintenance of myoblast model tissues, eLife, 11:e76409, 2022
[2] Francois Mazuel, Myriam Reffay, Vicard Du, Jean-Claude Bacri, Jean-Paul Rieu, Claire Wilhelm, Magnetic flattening of stem-cell spheroids indicates a size-dependent elastocapillary transition, PRL, 114(9), 2015

