# ThermoCorrection
version: 1.3
Author: Bene Hiebl
Date: 20.12.2021

An atmospheric and emissivity correction for georeferenced thermal infrared images from close range sensing
based on Caselles et al. (1996), Kodimalar et al. (2020), Wiecek (2011) and Minkina et al. (2016).

Input requirements:
•	Georeferenced thermo raster images (.tif) + timeinfo (.csv)
•	DSM (.tif)
•	FVC (.tif)
•	Red Band RGB image (optional, snow detection) (.tif)
•	Camera location
•	Emissivity for soil/vegetation
•	Vertical gradients of air temperature and humidity (.csv)
•	Radiance air, downwelling radiance (optional)

Define Input Parameters in Correction_Parameters.txt
For further infos see correction_documentation.pdf
