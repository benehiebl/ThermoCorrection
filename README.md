# ThermoCorrection
version: 1.3
Author: Bene Hiebl
Date: 20.12.2021

An atmospheric and emissivity correction for georeferenced thermal infrared images from close range sensing
based on Caselles et al. (1996), Kodimalar et al. (2020), Wiecek (2011) and Minkina et al. (2016).

Input requirements:
•	Georeferenced thermo raster images (.tif) + timeinfo
•	DSM
•	FVC
•	Red Band (optional)
•	Camera location
•	Vertical gradients air temperature and humidity
•	Emissivity for soil/vegetation
•	Radiance air, downwelling radiance (optional)

Define Input Parameters in Correction_Parameters.txt
For further infos see correction_documentation.pdf
