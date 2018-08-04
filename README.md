# Barebone-fault-PSHA

The purpose of this project is to develop an open source barebone probabilistic seismic hazard analysis (PSHA) code for fault sources. Some facts of code are listed as follows. 

  +	It follows the classical formulation of PSHA as described in McGuire (1995) and Pagani et al. (2014)
  +	Stirling approach is adopted for a dipping fault which bends along strike.
  +	For a bending rupture, Rx is taken as the one for the rupture segment with minimum Rrup (closest distance to the rupture plane)
  +	Fault or rupture geometry convention follows Stein and Wysession(2009)
  +	Distance conversion is based on geometries of hypothetic ruptures.
  + Strike direction is defined as first point toward ending point in the surface trace
  + It is working with some of the test cases in the PEER PSHA code verification project

  + It is developed under Windows 7 using free IDE of Code::Blocks and free compiler of GNU Fortran distributed with MSYS2. 
  + It has modular coding structure for easy extension.
  + It is built on many existing works. These works are cited in the comments in source code.
  + It uses keyword driven free format input.
  + A relatively generic ground motion prediction equation (GMPE) interface is signed for easy plugin of new GMPEs. Only two GMPEs (Sadigh et al. (1997) and Chiou & Youngs (2014)) are currently available.

# References

  + McGuire, R. K. (1995) “Probabilistic seismic hazard analysis and design earthquakes: closing the loop”, Bulletin of the Seismological Society of America, 85(5), pp. 1275-1284.
  + Stein, S., & Wysession, M. (2009). An introduction to seismology, earthquakes, and earth structure. John Wiley & Sons.
  + Pagani, M., Monelli, D., Weatherill, G. A. and Garcia, J. (2014). The OpenQuake-engine Book: Hazard. Global Earthquake Model (GEM) Technical Report 2014-08, doi: 10.13117/- GEM.OPENQUAKE.TR2014.08, 67 pages.
  + Sadigh at al. (1997) “Attenuation Relationships for Shallow Crustal Earthquakes Based on California Strong Motion Data”, Seismological Research Letters, 68(1), pp. 180-189
  + Chiou, B. S. J., & Youngs, R. R. (2014) “Update of the Chiou and Youngs NGA model for the average horizontal component of peak ground motion and response spectra”, Earthquake Spectra, 30(3), pp. 1117-1153.
