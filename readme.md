MuMaP (MUltiphae MAterial Properties) is a F2003-based code to calculate the
melt volume fraction from seismic signature.

There are three main modules:
global.f90 (utilities)
microgeodynamics.f90 (rock physics calculatios)
regional.f90 (data input-output utilities, object definitions)

The input data should be placed in the library named 'input'

Input format
##########################################################################
The format of the input data should be  the following:
Lat Lon LVL 410 Rnorm T1 ... Tn
where Lat and Lon give locations of each data point, LVL is the depth to LVL
410 is the depth to the top of the TZ and T1...Tn are temperature anomalies
calculated separately from the transition zone thickness
##########################################################################       	   	 

Output format
##########################################################################
Output is generated in the folder titled 'data' and a typical output file
has the format:
# Latitude, Longitude,  LVL(km), 410(km), Potential T (K),
# Vs_obs(m/s),  Vs_sol(m/s),Melt fraction,
# K/K0, G/G0,Residual, Amplitude
# where xx is the dihedral angle
# yy is the basalt percent in solid
# zzzz is reference potential temperature
##########################################################################

Module microgeodynamics.f90
##########################################################################
!> This module contains  a number of mineral and rock physics
  !!utility routines. There are a number of derived types defined
  !! also. fit contains the data for polynomial fitting parameters.
  !! It is used by functions to convert temperature into Vs, Vp,
  !! and density. Type melt contains parameters for use in EOS
  !! for melt. Composition contains information on the bulk
  !! composition of the reference mantle and the potential
  !! temperature. Type mantle contains the information from
  !! from mineral physics on the physical properties of the
  !! reference mantle for a given value of the variable composition.
  !! Finally, unit cell contains the information on the seismic
  !! signature of the reference mantle (subscript sol), the observed
  !! seismic signature (subscript obs), and the effective
  !! values (subscript eff). If the difference between reference and
  !! observed values can be explained by melting, then the
  !! melt fraction and melt fraction dependent contiguity
  !! are calculated. Temperature is the potential (not actual)
  !! temperature for the unit cell. For example, it can be 
  !! calculated from the depression/uplift of the transition zone
  !! boundaries in later modules. See the description of each
  !! function within the function body.
  !! Created by Saswata Hier-Majumder, August, 2012. !<
##########################################################################
Module regional.f90

##########################################################################
!> This module contains data types and functions applicable
  !!to the region of interest. Each location in the region
  !!is assigned a unit cell, which contains the physical properties.
  !!The parameters for the equation of state of the melt phase are
  !!stored in the derived type MAGMA. It is assumed that only one
  !!kind of melt exists in the region. Three dimensional location
  !!of each point in the region are stored in the array LOC. Typically
  !!this data is stored as lat,lon, depth (km). The variable composition
  !!is also regional, it contains the basalt fraction
  !!of the bulk ocmposition and the average potential temperature
  !!of the region.
  !! Created by Saswata Hier-Majumder August, 2012.!<
#############################################################################  