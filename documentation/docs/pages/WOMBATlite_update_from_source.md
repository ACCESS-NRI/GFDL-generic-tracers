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
3. Temperature‑dependent autotrophy and heterotrophy.
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
PAR that is seen by the average phytoplankton within that cell ('radbio', $PAR$, [W m<sup>-2</sup>]), where

$PAR(k) = \sum_{b=1}^3 \frac{(PAR^{top}(k,b) - PAR^{top}(k+1,b))}{ex_{bgr}(k,b) * \Delta z(k)}$

This ensures phytoplankton growth in the model responds to the mean light they experience in the cell, not just light at one point. See Eq. 19 from [Baird et al. (2020)](https://gmd.copernicus.org/articles/13/4503/2020/).

The euphotic depth is defined as the depth where `radbio` falls below the 1% threshold of incidient shortwave radiation and is therefore soley dependent on the concentration of chlorophyll within the water column.


### 2. Nutrient limitation of phytoplankton

At the start of each vertical loop the code computes biomass of phytoplankton (`biophy`, $B_{phy}$, [mmol C m<sup>-3</sup>]). Phytoplankton biomass 
is used to scale how both nitrate (`biono3`, $NO_{3}$, [mmol NO<sub>3</sub> m<sup>-3</sup>]) and dissolved iron (`biofer`, $dFe$, [µmol dFe m<sup>-3</sup>]) affect the growth of phytoplankton. Using compilations of marine phytoplankton and zooplankton communities, [Wickman et al. (2024)](https://www.science.org/doi/10.1126/science.adk6901) show that the nutrient affinity, $aff$, of a phytoplankton cell is related to its volume, $V$, via

$aff = V^{-0.57}$

Additionally, the authors demonstrate that the volume of the average phytoplankton cell is related to the density (i.e., concentration) of phytoplankton via

$V = (B_{phy})^{0.65}$
 
when combining panels c and f of their Figure 1. This then relates the affinity of an average cell to the concentration of phytoplankton biomass as

$aff = (B_{phy})^{-0.37}$

With this information, we allow the half-saturation terms for nitrate (`phy_kni(i,j,k)`, $K_{phy}^{N}$, [mmol N m<sup>-3</sup>]) and dissolved iron  (`phy_kfe(i,j,k)`, $K_{phy}^{Fe}$, [µmol dFe m<sup>-3</sup>]) uptake to vary as a function of phytoplankton biomass concentration. We set reference values for the half-saturation coefficient of nitrate (`phykn`, $K_{phy}^{N,0}$, [mmol N m<sup>-3</sup>]) and dissolved iron (`phykf`, $K_{phy}^{Fe,0}$, [µmol dFe m<sup>-3</sup>]) as input parameters to the model, and also set a threshold phytoplankton concentration (`phybiot`, $B_{phy}^{thresh}$, [mmol C m<sup>-3</sup>]) beneath which cell size cannot decrease and affinity can no longer increase. At this minimum, where affinity is maximised, the half-saturation coefficients are bounded to be 10% of their reference values.

$K_{phy}^{N} = K_{phy}^{N,0} * \max(0.1, \max(0.0, (B_{phy}-B_{phy}^{thresh}))^{0.37} )$
$K_{phy}^{Fe} = K_{phy}^{Fe,0} * \max(0.1, \max(0.0, (B_{phy}-B_{phy}^{thresh}))^{0.37} )$

Limitation of phytoplankton growth by nitrate (`phy_lnit(i,j,k)`, $L_{phy}^{N}$), [dimensionless]) then follows the Monod equation:

$L_{phy}^{N} = \frac{NO_3}{NO_3 + K_{phy}^{N}}$
 
where $NO_3$ is the ambient nitrate concentration (`biono3`, $NO_3$, [mmol N m<sup>-3</sup>]).

For iron the model employs an internal quota approach ([Droop, 1983](https://www.degruyterbrill.com/document/doi/10.1515/botm.1983.26.3.99/html)). Phytoplankton have a minimum iron quota (`phy_minqfe`, $Q_{phy}^{-Fe:C}$, [mol/mol]) and an optimal quota for growth (`phy_optqfe`, $Q_{phy}^{*Fe:C}$, [mol/mol]). The minimum iron quota, $Q_{phy}^{-Fe:C}$, is dependent on the chlorophyll content of the cell and the degree of nitrogen limitation according to

$$
\begin{align}
Q_{phy}^{-Fe:C} = & 0.00167 / 55.85 * Q_{phy}^{Chl:C} * 12 \\
&  + 1.21 \times 10^{-5} * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * L_{phy}^{N} \\
&  + 1.15 \times 10^{-4} * 14.0 / 55.85 / 7.625 * 0.5 * L_{phy}^{N}
\end{align}
$$

The minimum iron requirements of the cell for growth increases as the ratio of chlorophyll to carbon increases (`phy_chlc`, $Q_{phy}^{Chl:C}$, [mol/mol]) and as the cell becomes less nitrate limited (i.e., performs more nitrate reduction). This formulation and the coefficients applied to chlorophyll content and nitrate use derive from [Flynn & Hipkin (1999)](https://onlinelibrary.wiley.com/doi/10.1046/j.1529-8817.1999.3561171.x).

The Fe limitation factor (`phy_lfer(i,j,k)`, $L_{phy}^{Fe}$), [dimensionless]) is computed from the present Fe:C quota of the phytoplankton cells (`phy_Fe2C`, $Q_{phy}^{Fe:C}$, [mol/mol]) relative to these quotas.

$L_{phy}^{Fe} = \max(0.0, \min(1.0, \frac{ Q_{phy}^{Fe:C} - Q_{phy}^{-Fe:C} }{Q_{phy}^{*Fe:C}} )) $

If the cell is Fe‑replete with a quota that exceeds the minimum quota by as much as the optimal quota, then Fe does not limit growth ($L_{phy}^{Fe}$ = 1). If the cell is Fe‑deplete with a quota equal to the minimum quota, then the growth rate is reduced to zero. The optimal quota ($Q_{phy}^{*Fe:C}$) is therefore a measure of how much excess Fe is required to allow unrestricted growth.


### 3. Temperature‑dependent autotrophy and heterotrophy

**Autotrophy.**
The maximum potential growth rate for phytoplankton (`phy_mumax(i,j,k)`, $\mu_{phy}^{max}$, [day<sup>-1</sup>]) is prescribed by the temperature-dependent Eppley curve ([Eppley, 1972](https://spo.nmfs.noaa.gov/content/temperature-and-phytoplankton-growth-sea)). This formulation scales a reference growth rate at 0ºC via a power-law scaling with temperature (`Temp(i,j,k)`, $T$, [ºC]).

$\mu_{phy}^{max} = \mu_{phy}^{0^{\circ}C} \cdot β_{auto}^{(T)}$

In the above, both $\mu_{phy}^{0ºC}$ and $β_{auto}$ are reference values input to the model.

**Heterotrophy.**
Heterotrophic processes include mortality of phytoplankton and zooplankton, grazing rates of zooplankton and the remineralisation rate of detritus in the water column and sediments. These processes are scaled similarly to autotrophy, where some reference rate at 0ºC ($\mu_{het}^{0^{\circ}C}$, [day<sup>-1</sup>]) is multiplied by a power-law with temperature ($β_{hete}$).

$\mu_{het} = \mu_{het}^{0^{\circ}C} \cdot β_{hete}^{(T)}$


### 4. Light limitation of phytoplankton

Phytoplankton growth is limited by light through a photosynthesis–irradiance (P–I) relationship that links cellular chlorophyll content and available photosynthetically active radiation ('radbio', $PAR$, [W m<sup>-2</sup>]).

First, The initial slope of the P–I curve, (`phy_pisl`, $\alpha_{phy}$, [day<sup>-1</sup> (W m<sup>-2</sup>)<sup>-1</sup>]), determines how efficiently phytoplankton convert light into carbon fixation. It is scaled by the cellular chlorophyll-to-carbon ratio (`phy_chlc`, $Q_{phy}^{Chl:C}$, [mg/mg]).

$\alpha_{phy} = \max(\alpha_{phy}^{Chl} \cdot Q_{phy}^{Chl:C} \ , \ \alpha_{phy}^{Chl} \cdot Q_{phy}^{-Chl:C}$ )

where $\alpha_{phy}^{Chl}$ is the photosynthetic efficiency per unit chlorophyll in units of [day<sup>-1</sup> (W m<sup>-2</sup>)<sup>-1</sup> (mg/mg)<sup>-1</sup>] and $Q_{phy}^{-Chl:C}$ is the minimum chlorophyll to carbon ratio of the cell. This constraint prevents photosynthesis from collapsing unrealistically at low chlorophyll concentrations.

Second, light limitation (`phy_lpar(i,j,k)`, $L_{phy}^{PAR}$), [dimensionless]) is calculated using an exponential P–I formulation.

$L_{phy}^{PAR} = 1 - e^{- \alpha_{phy} PAR }$

At low irradiance (PAR), growth increases approximately linearly with light, while at high irradiance photosynthesis asymptotically saturates. Photoinhibition is not included in this formulation.


### 5. Growth of chlorophyll
