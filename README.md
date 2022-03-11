# SoilMirMod
Software for modeling multi-angular near-specular reflectance spectra of soiled solar mirrors

Installation on Linux (tested on Manjaro distro)
1) copy the source file in a folder and create a dedicated working folder where store experimental input-data
2) customize the pathfile QString pathRoot and QString fileExp at row 52 and 52 of soilmirmoc.cpp
3) install Qt5 (included QtCreator) and the following libraries: qwt , cminpack , blas , cblas
4) open SoilMirMod.pro by QtCreator and configure the project
5) compile

Quick user-guide
Arrange the spectra got far a given specimen as the following example

artificial soiled HairSpray#1 coupon <-----COMMENT
thetas 8 8 15 15 30 30 60 60 <------ incidence angle for clean and soiled pair of spectra
wl(nm)  Rh_clean    Rh_soil Rns1_clean  Rns1_soil Rns2_clean  Rns2_soil Rns3_clean  Rns3_soil <--------columns
320	0.0780618	0.0871	0.0780618	0.06531105969259	0.0783686	0.06546166532275	0.126388	0.117609475407
325	0.188149	0.1749	0.188149	0.1453782157093	0.190443	0.1541748468521	0.241975	0.2239781813101
(...)

Run the software or by QtCreator o from terminal, in the build directory, with the command ./SoilMirMod
Load the desired file
Manually adjust the parameters to keep confidence
Enable for fit those most relevant and push "Best Fit in"
limit the wavelength range excluding the most noisy parts
adjust the vertical scale
