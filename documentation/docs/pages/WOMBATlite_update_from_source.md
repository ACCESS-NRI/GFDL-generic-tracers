# Description of the WOMBATlite ocean biogeochemical model
## Subroutine - "update_from_source"

`!         (\___/)  .-.   .-. .--. .-..-..---.  .--. .-----.            !`\
`!         / o o \  : :.-.: :: ,. :: '' :: .; :: .; :'-. .-'            !`\
`!        (   "   ) : :: :: :: :: :: .. ::   .':    :  : :              !`\
`!         \__ __/  : '' '' ;: :; :: :; :: .; :: :: :  : :              !`\
`!                   '.,'.,' '.__.':_;:_;:___.':_;:_;  :_;              !`\
`!                                                                      !`\
`!  World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT)  !`


The subroutine `generic_WOMBATlite_update_from_source` is the heart of the World Ocean Model of Biogeochemistry And Trophic‑dynamics. 
Its purpose is to apply biological source–sink terms to ocean tracers (nutrients, phytoplankton, zooplankton, detritus, iron and carbon pools) 
on each NPZD time‑step. The subroutine is documented internally by a list of numbered steps (see code comments). These steps are:

1. Light attenuation through the water column.
2. Nutrient limitation of phytoplankton.
3. Temperature‑dependent heterotrophy and remineralisation.
4. Light limitation of phytoplankton.
5. Growth of chlorophyll.
6. Phytoplankton uptake of iron.
7. Iron chemistry (precipitation, scavenging and coagulation).
8. Mortality scalings and grazing.
9. CaCO3 calculations.
10. Sources and sinks.
11. Tracer tendencies (update tracer concentrations).
12. Check for conservation of mass.

Below is a step‑by‑step explanation of each section together with the key equations. Variable names in grey follow the Fortran code, while 
variable names in math font are pointers to the equations; i,j,k refer to horizontal and vertical indices; square brackets denote units.


### 1. Light attenuation through the water column

Photosynthetically available radiation (PAR) is split into blue, green and red wavelengths. Surface incoming shortwave radiation is multiplied by 
0.43 to return the photosynthetically active radiation (PAR, [W m<sup>-2</sup>) of total shortwave radiation, and is then split evenly into each 
of blue, green and red light bands. At the top (`par_bgr_top(k,b)`, $PAR^{top}$) and mid‑point (`par_bgr_mid(k,b)`, $PAR^{mid}$) of each layer `k` 
the subroutine calculates the downward irradiance by exponential decay of each band `b` through the layer thickness (`dzt(i,j,k)`, $\Delta z$, [m]) using 
band‑specific attenuation coefficients (`ek_bgr(k,b)`, $ex_{bgr}$, [m<sup>-1</sup>]) obtained from a look‑up table (`zbgr`). The look-up table relates 
chlorophyll concentration [mg m<sup>-3</sup>] to the extinction rate of blue, green and red light through the water column. For example, 
the PAR in the blue band (`b=1`) at the top of level k is computed as

$PAR^{top}(k,1) = PAR^{top}(k-1,1) * e^{(-ex_{bgr}(k-1,1) * \Delta z(k-1))}$

The irradiance in the red band (`b=3`) at the mid point of layer `k`, for example, is equal to 

$PAR^{mid}(k,3) = PAR^{mid}(k-1,3) * e^{(-0.5*(ex_{bgr}(k-1,3) * \Delta z(k-1) + ex_{bgr}(k,3) * \Delta z(k)))}$

The total PAR available to phytoplantkon is assumed to be the sum of the blue, green and red bands. Because we assume that phytoplankton are 
homogenously distributed within a layer `k`, but we do not assume that light is homogenously distributed within that layer, we solve for the 
PAR that is seen by the average phytoplankton within that cell ('radbio', $PAR$, [W m<sup>-2</sup}], where

$PAR(k) = \sum_{b=1}^3 \frac{(PAR^{top}(k,b) - PAR^{top}(k+1,b))}{ex_{bgr}(k,b) * \Delta z(k)}$

This ensures phytoplankton growth in the model responds to the mean light they experience in the cell, not just light at one point.

The euphotic depth is defined as the depth where `radbio` falls below the 1% threshold of incidient shortwave radiation. 


### 2. Nutrient limitation of phytoplankton

 
 The maximum potential growth rate for phytoplankton, phy_mumax(i,j,k), is prescribed by an Eppley curve based on water temperature Temp(i,j,k):

\mu_{\max} = \mu_{ ext{ref}} imes \expigl(r_{ ext{Eppley}} (T- T_{ ext{ref}})igr),

where mu_ref and r_Eppley are constants.
