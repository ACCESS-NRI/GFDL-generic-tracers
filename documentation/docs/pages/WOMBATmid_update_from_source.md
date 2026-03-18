# Description of the WOMBATmid ocean biogeochemical model
## Subroutine - "update_from_source"

`!         (\___/)  .-.   .-. .--. .-..-..---.  .--. .-----.            !`\
`!         / o o \  : :.-.: :: ,. :: '' :: .; :: .; :'-. .-'            !`\
`!        (   "   ) : :: :: :: :: :: .. ::   .':    :  : :              !`\
`!         \__ __/  : '' '' ;: :; :: :; :: .; :: :: :  : :              !`\
`!                   '.,'.,' '.__.':_;:_;:___.':_;:_;  :_;              !`\
`!                                                                      !`\
`!  World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT)  !` 

---

The subroutine `generic_WOMBATmid_update_from_source` is the heart of the World Ocean Model of Biogeochemistry And Trophic‑dynamics. 
Its purpose is to apply biological source–sink terms to ocean tracers (nutrients, phytoplankton, zooplankton, bacteria, particulate
detritus, dissolved and particulate iron, dissolved organics, alkalinity, nitrous oxide, oxygen and carbon pools) at each time‑step. 
The subroutine is documented internally by a list of numbered steps (see code comments). These steps are:

1. Light attenuation through the water column.
2. Nutrient limitation of phytoplankton.
3. Temperature-dependent metabolism and POC-->DOC.
4. Light limitation of phytoplankton.
5. Realized growth rate of phytoplankton.
6. Dissolved organic carbon release by phytoplankton.
7. Growth of chlorophyll.
8. Phytoplankton uptake of iron.
9. Phytoplankton uptake of silicic acid.
10. Iron chemistry (scavenging, coagulation, dissolution).
11. Biogenic silica dissolution. 
12. Mortality terms.
13. Zooplankton grazing, egestion, excretion and assimilation.
14. Calcium carbonate production and dissolution.
15. Implicit nitrogen fixation.
16. Facultative bacterial heterotrophy.
17. Chemoautotrophy.
18. Nominal oxidation state of dissolved organic carbon.
19. Tracer tendencies.
20. Check for conservation of mass.
21. Additional operations on tracers.
22. Sinking rate of particulates.
23. Sedimentary processes.

Below is a step‑by‑step explanation of each section together with the key equations. Variable names in `grey` follow the Fortran code, while 
variable names in $math font$ are pointers to the equations; `i,j,k` refer to horizontal and vertical indices; [square brackets] denote units. 
If a variable is without i,j,k dimensions, this variable is held as a scalar and not an array.

The model carries tracers in [mol kg-1]. That is, moles of solute/tracer per kilogram of seawater (i.e., molarity). Some calculations herein are performed by converting tracers to units of [mmol m-3] or in the case of dissolved iron [µmol m-3]. However, we stress that all tracer tendency terms are converted back to [mol kg-1 s-1] when sources and sinks are applied.
---


### 1. Light attenuation through the water column.

Photosynthetically available radiation (PAR) is split into blue, green and red wavelengths. The incoming visible (photosynthetically available) short wave radiation flux (PAR, [W m<sup>-2</sup>]) is received from the physical model, and is then split evenly into each of blue, green and red light bands.

At the top (`par_bgr_top(k,b)`, $PAR^{top}$) and mid‑point (`par_bgr_mid(k,b)`, $PAR^{mid}$) of each layer `k` we calculate the downward irradiance by exponential decay of each band `b` through the layer thickness (`dzt(i,j,k)`, $\Delta z$, [m]) using band‑specific attenuation coefficients (`ek_bgr(k,b)`, $ex_{bgr}$, [m<sup>-1</sup>]). These attenuation coefficients are related to the concentration of chlorophyll (`chl`, [mg m<sup>-3</sup>]), organic detritus (`ndet`, [mg N m<sup>-3</sup>]) and calcium carbonate (`carb`, [kg m<sup>-3</sup>]) in the water column. For chlorophyll, attentuation coefficients for each of blue, green and red light are retrieved from the look-up table of [Morel & Maritorena (2001)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2000JC000319) (their Table 2) that explicitly relates chlorophyll concentration to attenutation rates and accounts for the packaging effect of chlorophyll in larger cells. For organic detritus, attenutation coefficients for blue, green and red light are taken from [Dutkiewicz et al. (2015)](https://bg.copernicus.org/articles/12/4447/2015/bg-12-4447-2015.html) (their Fig. 1b), while for calcium carbonate we take the coefficients defined in [Soja-Wozniak et al. (2019)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JC014998). For both detritus and calcium carbonate, these studies provide concentration normalized attenutation coefficients, which must be multipled against concentrations to retrieve the correct units of [m<sup>-1</sup>].

Because WOMBAT-mid has two forms of phytoplankton (nanophytoplankton and microphytoplankton) with their own chlorophyll quotas and two forms of particulate detritus (small and large), we sum both chlorophyll pools and particulate detritus pools to return the total chlorophyll and the total particulate detritus. 

As an example, the PAR in the blue band (`b=1`) at the top of level k is computed as

$PAR^{top}(k,1) = PAR^{top}(k-1,1) * e^{(-ex_{bgr}(k-1,1) * \Delta z(k-1))}$

where the total attenutation rate of blue light in the grid cell above `k` is the sum of attenuation due to all particulates in that grid cell, which includes chlorophyll, detritus and calcium carbonate:

$ex_{bgr}(k-1,1) = ex_{chl}(k-1,1) + ex_{det}(k-1,1) + ex_{CaCO_3}(k-1,1)$

where
- $ex_{chl}(k-1,1)$ is the attenuation rate of blue light (`b=1`) in the overlying grid cell (`k=k-1`) due to chlorophyll (`zbgr(2,ichl)`, [m<sup>-1</sup>])
- $ex_{det}(k-1,1)$ is the attenuation rate of blue light (`b=1`) in the overlying grid cell (`k=k-1`) due to detritus (`ndet * dbgr(1)`, [m<sup>-1</sup>])
- $ex_{CaCO_3}(k-1,1)$ is the attenuation rate of blue light (`b=1`) in the overlying grid cell (`k=k-1`) due to calcium carbonate (`carb * cbgr(1)`, [m<sup>-1</sup>])

The irradiance in the red band (`b=3`) at the mid point of layer `k`, in contrast, is equal to 

$PAR^{mid}(k,3) = PAR^{mid}(k-1,3) * e^{(-0.5*(ex_{bgr}(k-1,3) * \Delta z(k-1) + ex_{bgr}(k,3) * \Delta z(k)))}$

where
- $PAR^{mid}(k-1,3)$ is the red light (`b=3`) at the mid-point of the overlying grid cell (`par_bgr_mid(k-1,3)`, [W m<sup>-2</sup>])
- $ex_{bgr}(k-1,3)$ is the total attenuation of red light (`b=3`) in the overlying grid cell (`ek_bgr(k-1,3)`, [m<sup>-1</sup>])
- $ex_{bgr}(k,3)$ is the total attenuation of red light (`b=3`) in the current grid cell (`ek_bgr(k,3)`, [m<sup>-1</sup>])
- $\Delta z(k-1)$ and $\Delta z(k)$ are the grid cell thicknesses of the overlying and current grid cells (`dzt(i,j,k)`, [m])

The total PAR available to phytoplantkon is assumed to be the sum of the blue, green and red bands. Because we assume that phytoplankton are 
homogenously distributed within a layer `k`, but we do not assume that light is homogenously distributed within that layer, we solve for the 
PAR that is seen by the average phytoplankton within that cell (`radbio`, $PAR$, [W m<sup>-2</sup>])

$PAR(k) = \sum_{b=1}^3 \dfrac{(PAR^{top}(k,b) - PAR^{top}(k+1,b))}{ex_{bgr}(k,b) * \Delta z(k)}$

where
- $PAR^{top}(k,b)$ is the incoming photosynthetically active radiation at the top of grid cell `k` and light band `b` (`par_bgr_top(k,b)`, [W m<sup>-2</sup>])
- $ex_{bgr}(k,b)$ is the attenuation rate of light band `b` in grid cell `k` (`ek_bgr(k,b)`, [m<sup>-1</sup>])
- $\Delta z(k)$ is the grid cell thickness of grid cell `k` (`dzt(i,j,k)`, [m])

This ensures phytoplankton growth in the model responds to the mean light they experience in the cell, not just light at one point. See Eq. 19 from [Baird et al. (2020)](https://gmd.copernicus.org/articles/13/4503/2020/).

The euphotic depth (`zeuphot(i,j)`, [m]) is defined as the depth where `radbio` falls below the 1% threshold of incidient shortwave radiation or below 0.01 W m<sup>-2</sup>, whichever is shallower.

---

### 2. Nutrient limitation of phytoplankton.

At the start of each vertical loop `k=1` through `k=kmax` the code computes the biomass of nano-phytoplankton (`biophy`, $B_{np}$, [mmol C m<sup>-3</sup>]) and micro-phytoplankton (`biodia`, $B_{mp}$, [mmol C m<sup>-3</sup>]). Phytoplankton biomass is used to scale how nitrogen in the form of nitrate (`biono3`, $NO_{3}$, [mmol N m<sup>-3</sup>]) and ammonium (`bionh4`, $NH_{4}$, [mmol N m<sup>-3</sup>]), dissolved iron (`biofer`, $dFe$, [µmol dFe m<sup>-3</sup>]) and silicic acid in the case of micro-phytoplankton (`biosil`, $Si(OH)_{4}$, [mmol S m<sup>-3</sup>]) affect the growth of phytoplankton. Using compilations of marine phytoplankton and zooplankton communities, [Wickman et al. (2024)](https://www.science.org/doi/10.1126/science.adk6901) show that the nutrient affinity, $aff$, of a phytoplankton cell is related to its volume, $V$, via

$aff = V^{-0.57}$

Additionally, the authors demonstrate that the volume of the average phytoplankton cell is related to the density (i.e., concentration) of phytoplankton via

$V = (B_{phy})^{0.65}$
 
when combining panels c and f of their Figure 1. This then relates the affinity of an average cell to the concentration of phytoplankton biomass as

$aff = (B_{phy})^{-0.37}$

With this information, we allow the half-saturation terms for nitrogen (`phy_kni(i,j,k)`, $K_{np}^{N}$, [mmol N m<sup>-3</sup>]; `dia_kni(i,j,k)`, $K_{mp}^{N}$, [mmol N m<sup>-3</sup>]), dissolved iron  (`phy_kfe(i,j,k)`, $K_{np}^{Fe}$, [µmol dFe m<sup>-3</sup>]; `dia_kfe(i,j,k)`, $K_{mp}^{Fe}$, [µmol dFe m<sup>-3</sup>]) and silic acid (`dia_ksi(i,j,k)`, $K_{mp}^{Si}$, [µmol Si m<sup>-3</sup>]) uptake to vary as a function of phytoplankton biomass concentration. We set reference values for the half-saturation coefficient of nitrogen (`phykn`, $K_{np}^{N,0}$, [mmol N m<sup>-3</sup>]; `diakn`, $K_{mp}^{N,0}$, [mmol N m<sup>-3</sup>]), dissolved iron (`phykf`, $K_{np}^{Fe,0}$, [µmol dFe m<sup>-3</sup>]; `diakf`, $K_{mp}^{Fe,0}$, [µmol dFe m<sup>-3</sup>]) and silicic acid (`diaks`, $K_{mp}^{Si,0}$, [µmol Si m<sup>-3</sup>]) as input parameters to the model, and also set thresholds of nano-phytoplankton concentration (`phybiot`, $B_{np}^{thresh}$, [mmol C m<sup>-3</sup>]) and micro-phytoplankton concentration (`diabiot`, $B_{mp}^{thresh}$, [mmol C m<sup>-3</sup>]) beneath which cell size cannot decrease and affinity can no longer increase. At this minimum, where affinity is maximised, the half-saturation coefficients are bounded to be 10% of their reference values.

$K_{np}^{N} = K_{np}^{N,0} * \max(0.1, \max(0.0, (B_{np}-B_{np}^{thresh}))^{0.37} )$\
$K_{np}^{Fe} = K_{np}^{Fe,0} * \max(0.1, \max(0.0, (B_{np}-B_{np}^{thresh}))^{0.37} )$\
$K_{mp}^{N} = K_{mp}^{N,0} * \max(0.1, \max(0.0, (B_{mp}-B_{mp}^{thresh}))^{0.37} )$\
$K_{mp}^{Fe} = K_{mp}^{Fe,0} * \max(0.1, \max(0.0, (B_{mp}-B_{mp}^{thresh}))^{0.37} )$\
$K_{mp}^{Si} = K_{mp}^{Si,0} * \max(0.1, \max(0.0, (B_{mp}-B_{mp}^{thresh}))^{0.37} )$

where
- $K_{np}^{N}$ and $K_{mp}^{N}$ are the half-saturation coefficients for nitrogen uptake by nano- and micro-phytoplankton (`phy_kni(i,j,k)` and `dia_kni(i,j,k)`, [mmol N m<sup>-3</sup>])
- $K_{np}^{Fe}$ and $K_{mp}^{Fe}$ are the half-saturation coefficients for iron uptake by nano- and micro-phytoplankton (`phy_kfe(i,j,k)` and `dia_kfe(i,j,k)`, [µmol Fe m<sup>-3</sup>])
- $K_{mp}^{Si}$ are the half-saturation coefficients for nitrogen uptake by nano- and micro-phytoplankton (`dia_ksi(i,j,k)`, [mmol Si m<sup>-3</sup>])


**Limitation of phytoplankton growth by nitrogen** (`phy_lnit(i,j,k)`, $L_{np}^{N}$), [dimensionless]; `dia_lnit(i,j,k)`, $L_{mp}^{N}$), [dimensionless]) is split between ammonium (`phy_lnh4(i,j,k)`, $L_{np}^{NH_4}$), [dimensionless]; `dia_lnh4(i,j,k)`, $L_{mp}^{NH_4}$), [dimensionless]) and nitrate (`phy_lno3(i,j,k)`, $L_{np}^{NO_3}$), [dimensionless]; `dia_lno3(i,j,k)`, $L_{mp}^{NO_3}$), [dimensionless]). Phytoplankton preferentially consume and grow on ammonium because it is most efficiently converted to glutamate for biomass synthesis, while nitrate must be first reduced within the cell ([Dortch, 1990](https://www.jstor.org/stable/24842258)). To represent this preference, we follow [Buchanan et al., 2025](https://bg.copernicus.org/articles/22/4865/2025/) who assert a 5-fold preference of phytoplankton for ammonium over nitrate and show that this reproduces preferences of ammonium-fueled growth in ocean field data.

$l_{np}^{NH_4} = \dfrac{NH_4}{NH_4 + K_{np}^{N}}$\
$l_{np}^{NO_3} = \dfrac{NO_3}{NO_3 + K_{np}^{N}}$\
$l_{np}^{N} = \dfrac{NH_4 + NO_3}{NH_4 + NO_3 + K_{np}^{N}}$\
$L_{np}^{NH_4} = \dfrac{5 \cdot l_{np}^{N} l_{np}^{NH_4}}{l_{np}^{NO_3} + 5 \cdot l_{np}^{NH_4}}$\
$L_{np}^{NO_3} = \dfrac{l_{np}^{N} l_{np}^{NO_3}}{l_{np}^{NO_3} + 5 \cdot l_{np}^{NH_4}}$\
$L_{np}^{N} = L_{np}^{NH_4} + L_{np}^{NO_3}$

where
- $NH_4$ is the in situ concentration of ammonium (`bionh4`, [mmol N m<sup>-3</sup>]) 
- $NO_3$ is the in situ concentration of nitrate (`biono3`, [mmol N m<sup>-3</sup>]) 
- $l_{np}^{NH_4}$ is the limitation term of nano-phytoplankton growth on ammonium before preferencing (`phy_limnh4`, [dimensionless]) 
- $l_{np}^{NO_3}$ is the limitation term of nano-phytoplankton growth on nitrate before preferencing (`phy_limno3`, [dimensionless]) 
- $l_{np}^{N}$ is the limitation term of nano-phytoplankton growth on nitrogen before preferencing (`phy_limdin`, [dimensionless]) 
- $L_{np}^{NH_4}$ is the limitation term of nano-phytoplankton growth on ammonium (`phy_lnh4(i,j,k)`, [dimensionless]) 
- $L_{np}^{NO_3}$ is the limitation term of nano-phytoplankton growth on nitrate (`phy_lno3(i,j,k)`, [dimensionless]) 
- $L_{np}^{N}$ is the limitation term of nano-phytoplankton growth on nitrogen (`phy_lnit(i,j,k)`, [dimensionless]) 


The same set of equations are applied to micro-phytoplankton:

$l_{mp}^{NH_4} = \dfrac{NH_4}{NH_4 + K_{mp}^{N}}$\
$l_{mp}^{NO_3} = \dfrac{NO_3}{NO_3 + K_{mp}^{N}}$\
$l_{mp}^{N} = \dfrac{NH_4 + NO_3}{NH_4 + NO_3 + K_{mp}^{N}}$\
$L_{mp}^{NH_4} = \dfrac{5 \cdot l_{mp}^{N} l_{mp}^{NH_4}}{l_{mp}^{NO_3} + 5 \cdot l_{mp}^{NH_4}}$\
$L_{mp}^{NO_3} = \dfrac{l_{mp}^{N} l_{mp}^{NO_3}}{l_{mp}^{NO_3} + 5 \cdot l_{mp}^{NH_4}}$\
$L_{mp}^{N} = L_{mp}^{NH_4} + L_{mp}^{NO_3}$

where
- $NH_4$ is the in situ concentration of ammonium (`bionh4`, [mmol N m<sup>-3</sup>]) 
- $NO_3$ is the in situ concentration of nitrate (`biono3`, [mmol N m<sup>-3</sup>]) 
- $l_{mp}^{NH_4}$ is the limitation term of micro-phytoplankton growth on ammonium before preferencing (`dia_limnh4`, [dimensionless]) 
- $l_{mp}^{NO_3}$ is the limitation term of micro-phytoplankton growth on nitrate before preferencing (`dia_limno3`, [dimensionless]) 
- $l_{mp}^{N}$ is the limitation term of micro-phytoplankton growth on nitrogen before preferencing (`dia_limdin`, [dimensionless]) 
- $L_{mp}^{NH_4}$ is the limitation term of micro-phytoplankton growth on ammonium (`dia_lnh4(i,j,k)`, [dimensionless]) 
- $L_{mp}^{NO_3}$ is the limitation term of micro-phytoplankton growth on nitrate (`dia_lno3(i,j,k)`, [dimensionless]) 
- $L_{mp}^{N}$ is the limitation term of micro-phytoplankton growth on nitrogen (`dia_lnit(i,j,k)`, [dimensionless]) 


Note that although phytoplankton prefer $NH_4$ over $NO_3$, as $NO_3$ becomes more abundant than $NH_4$ the $L_{mp}^{NO_3}$ term begins to exceed the $L_{mp}^{NH_4}$ term such that phytoplankton switch from regenerated production ($NH_4$-based) to new production ($NO_3$-based). This reproduces the known switch of phytoplankton from regenerated to new production that is observed in the real ocean ([Dugdale & Goering, 1967](https://doi.org/10.4319/lo.1967.12.2.0196), [Buchanan et al., 2025](https://bg.copernicus.org/articles/22/4865/2025/)). Furthermore, if $K_{mp}^{N}$ > $K_{np}^{N}$, this ensures that (i) micro-phytoplankton are less competitive for $NH_4$ than nano-phytoplankton at any concentration and (ii) micro-phytoplankton growth is greater than nano-phytoplankton under abundant $NO_3$, which is consistent with theory and observations ([Fawcett et al., 2011](https://doi.org/10.1038/ngeo1265), [Glibert et al., 2016](https://doi.org/10.1002/lno.10203))

**Limitation of phytoplankton growth by iron** follows an internal quota approach ([Droop, 1983](https://www.degruyterbrill.com/document/doi/10.1515/botm.1983.26.3.99/html)). Phytoplankton have a minimum iron quota (`phy_minqfe`, $Q_{np}^{-Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]; `dia_minqfe`, $Q_{mp}^{-Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]) and an optimal quota for growth (`phy_optqfe`, $Q_{np}^{*Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]; `dia_optqfe`, $Q_{mp}^{*Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]). The minimum iron quota, $Q_{np}^{-Fe:C}$ and $Q_{mp}^{-Fe:C}$, is dependent on three terms that each correspond to the iron required by photosystems, respiration and nitrate reduction ([Flynn & Hipkin, 1999](https://onlinelibrary.wiley.com/doi/10.1046/j.1529-8817.1999.3561171.x)):

$$
\begin{align}
Q_{np}^{-Fe:C} = & 0.00167 / 55.85 * Q_{np}^{Chl:C} * 12 \\
&  + 1.21 \times 10^{-5} * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * L_{np}^{N} \\
&  + 1.15 \times 10^{-4} * 14.0 / 55.85 / 7.625 * 0.5 * L_{np}^{NO_3}
\end{align}
$$

$$
\begin{align}
Q_{mp}^{-Fe:C} = & 0.00167 / 55.85 * Q_{mp}^{Chl:C} * 12 \\
&  + 1.21 \times 10^{-5} * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * L_{mp}^{N} \\
&  + 1.15 \times 10^{-4} * 14.0 / 55.85 / 7.625 * 0.5 * L_{mp}^{NO_3}
\end{align}
$$

The first term reflects the amount of iron required for photosystems I and II. `0.00167/55.85` is equivalent to the grams of Fe per gram of chlorophyll divided by the grams of Fe per mol Fe, giving mol Fe per gram chlorophyll. This term is multipled by the chlorophyll to carbon ratio of the phytoplantkon cell (`phy_chlc`, $Q_{np}^{Chl:C}$, [mol C (mol C)<sup>-1</sup>]; `dia_chlc`, $Q_{mp}^{Chl:C}$, [mol C (mol C)<sup>-1</sup>]) and grams of C per mol C, returning mol Fe per mol C. At a healthy chlorophyll:C ratio of 0.03, this term returns an Fe:C ratio of roughly 10 µmol:mol, which reproduces well known requirements of phytoplankton cells ([Morel, Rueter & Price, 1991](https://www.jstor.org/stable/43924569)). The second term, representing the respiratory iron requirement, is derived from [Flynn & Hipkin (1999)](https://onlinelibrary.wiley.com/doi/10.1046/j.1529-8817.1999.3561171.x) who estimated 1.21 $\times 10^{-5}$ grams Fe per gram N assimilated into the cell, which is converted to mol Fe per mol C with 14 g N per mol N divided by 55.85 g Fe per mol Fe $\times$ 7.625 mol C per mol N. This second term assumes that respiration is reduced as growth becomes more limited by available nitrogen (`phy_lnit(i,j,k)`, $L_{np}^{N}$, [dimensionless]; `dia_lnit(i,j,k)`, $L_{mp}^{N}$, [dimensionless]). Finally, the third term represents the iron required by nitrate/nitrite reduction. Nitrate assimilation requires roughly 1.8-fold more iron than ammonia assimilation ([Raven, 1988](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/j.1469-8137.1988.tb04196.x)). [Flynn & Hipkin (1999)](https://onlinelibrary.wiley.com/doi/10.1046/j.1529-8817.1999.3561171.x) estimated a demand of 1.15 $\times 10^{-4}$ g Fe per mol NO$_3$ reduced, which is accounted for by the nitrate limitation term (`phy_lno3(i,j,k)`, $L_{np}^{NO_3}$, [dimensionless]; `dia_lno3(i,j,k)`, $L_{mp}^{NO_3}$, [dimensionless])). Note that the `1.5` is designed to account for dark respiration (i.e., respiration when the cells are not growing) and the `0.5` refers to the fact that during cell division the cell must reinstate half of its Fe reserves. 

The Fe limitation factor (`phy_lfer(i,j,k)`, $L_{np}^{Fe}$, [dimensionless]; `dia_lfer(i,j,k)`, $L_{mp}^{Fe}$, [dimensionless]) is then computed from the present Fe:C quota of the phytoplankton cells (`phy_Fe2C`, $Q_{np}^{Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]; `dia_Fe2C`, $Q_{mp}^{Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]) relative to the minimum and optimal quotas.

$L_{np}^{Fe} = \max\left(0.0, \min\left(1.0, \dfrac{ Q_{np}^{Fe:C} - Q_{np}^{-Fe:C} }{Q_{np}^{*Fe:C}} \right)\right)$

where
- $Q_{np}^{-Fe:C}$ is the minimum Fe:C quota of the nano-phytoplankton cell (`phy_minqfe`, [mol Fe (mol C)<sup>-1</sup>])
- $Q_{np}^{*Fe:C}$ is the optimal Fe:C quota of the nano-phytoplankton cell (`phyoptqf`, [mol Fe (mol C)<sup>-1</sup>])
- $Q_{np}^{Fe:C}$ is the in situ Fe:C quota of the nano-phytoplankton cell (`phy_Fe2C`, [mol Fe (mol C)<sup>-1</sup>])

$L_{mp}^{Fe} = \max\left(0.0, \min\left(1.0, \dfrac{ Q_{mp}^{Fe:C} - Q_{mp}^{-Fe:C} }{Q_{mp}^{*Fe:C}} \right)\right)$

where
- $Q_{mp}^{-Fe:C}$ is the minimum Fe:C quota of the micro-phytoplankton cell (`dia_minqfe`, [mol Fe (mol C)<sup>-1</sup>])
- $Q_{mp}^{*Fe:C}$ is the optimal Fe:C quota of the micro-phytoplankton cell (`diaoptqf`, [mol Fe (mol C)<sup>-1</sup>])
- $Q_{mp}^{Fe:C}$ is the in situ Fe:C quota of the micro-phytoplankton cell (`dia_Fe2C`, [mol Fe (mol C)<sup>-1</sup>])

If the cell is Fe‑replete with a quota that exceeds the minimum quota by as much as the optimal quota, then Fe does not limit growth ($L_{np}^{Fe}$ = 1; $L_{mp}^{Fe}$ = 1). If the cell is Fe‑deplete with a quota equal to or less than the minimum quota, then the growth rate is reduced to zero. The optimal quota ($Q_{np}^{*Fe:C}$; $Q_{mp}^{*Fe:C}$) is therefore a measure of how much excess Fe is required to allow unrestricted growth.

**Limitation of micro-phytoplankton growth by silicic acid** is computed as a gating constraint on division via:

$L_{mp}^{Si} = \min\left( 1.0, \max\left( 0.0, \dfrac{ Q_{mp}^{Si:C} - Q_{mp}^{-Si:C} }{ Q_{mp}^{*Si:C} - Q_{mp}^{-Si:C} }  \right) \right)$

where
- $Q_{mp}^{-Si:C}$ is the minimum Si:C quota of the micro-phytoplankton cell (`diaminqs`, [mol Si (mol C)<sup>-1</sup>])
- $Q_{mp}^{*Si:C}$ is the optimal Si:C quota of the micro-phytoplankton cell (`diaoptqs`, [mol Si (mol C)<sup>-1</sup>])

This formulation treats silicification as linearly limiting to growth between the minimum and optimal quotas. Above the optimal quota silica limitation does not exist. This reflects evidence that diatoms division is structurally constrained by silica until a threshold reserve is reached, at which point division can proceed ([Martin-Jézéquel, Hildebrand & Brzezinski, 2003](https://onlinelibrary.wiley.com/doi/full/10.1046/j.1529-8817.2000.00019.x)). This treatment is also supported by weak or even negative relationships between Si:C quotas and growth rates of marine diatoms ([María Mejía et al., 2013](https://www.sciencedirect.com/science/article/pii/S001670371300344X?via%3Dihub)) and is consistent with the apparent increase in Si:C quotas under Fe-limited growth ([Hutchins & Bruland, 1998](https://www.nature.com/articles/31203), [Takeda, 1998](https://www.nature.com/articles/31674)), which suggests that Si:C quotas can be decoupled from growth. 

---


### 3. Temperature-dependent metabolism and POM-->DOM.

**Autotrophy**\
The maximum potential growth rate for nano-phytoplankton (`phy_mumax(i,j,k)`, $\mu_{np}^{max}$, [day<sup>-1</sup>]) and micro-phytoplankton (`dia_mumax(i,j,k)`, $\mu_{mp}^{max}$, [day<sup>-1</sup>]) is prescribed by the temperature-dependent Eppley curve ([Eppley, 1972](https://spo.nmfs.noaa.gov/content/temperature-and-phytoplankton-growth-sea)). This formulation scales a reference growth rate at 0ºC via a power-law scaling with temperature (`Temp(i,j,k)`, $T$, [ºC]).

$\mu_{np}^{max} = \mu_{np}^{0^{\circ}C} \cdot (β_{np})^{T}$

$\mu_{mp}^{max} = \mu_{mp}^{0^{\circ}C} \cdot (β_{mp})^{T}$

where
- $\mu_{np}^{0^{\circ}}C$ is the rate of nano-phytoplankton growth at 0ºC (`abioa_phy`, [s<sup>-1</sup>])
- $β_{np}$ is the base temperature-sensitivity coefficient for autotrophy by nano-phytoplankton (`bbioa_phy`, [dimenionless])
- $\mu_{mp}^{0^{\circ}}C$ is the rate of micro-phytoplankton growth at 0ºC (`abioa_dia`, [s<sup>-1</sup>])
- $β_{np}$ is the base temperature-sensitivity coefficient for autotrophy by micro-phytoplankton (`bbioa_dia`, [dimenionless])
- $T$ is in situ water temperature (`Temp(i,j,k)`, [ºC])

In the above, $\mu_{np}^{0ºC}$, $\mu_{mp}^{0ºC}$, $β_{np}$ and $β_{mp}$ are reference values input to the model at run time. This allows the user to configure nano-phytoplankton and micro-phytoplankton with different maximum potential growth rates and different sensitivities to temperature ([Anderson et al., 2021](https://www.nature.com/articles/s41467-021-26651-8)).

**Heterotrophy**\
Heterotrophic processes include mortality of ecosystem functional types, grazing rates of zooplankton, growth rates of heterotrophic bacteria consuming dissovled organic matter (DOC and DON) and the hydrolysation rate of particulate detritus in the water column and sediments. These processes are scaled similarly to autotrophy, where some reference rate at 0ºC ($\mu_{het}^{0ºC}$, [<sup>s-1</sup>]) is multiplied by a power-law with temperature ($β_{hete}$). Each heterotrophic process has a different $\mu_{het}^{0ºC}$ value and we expand on this later under the mortality, grazing and bacterial heterotrophy sections. However, the basic formulation for scaling heterotrophic metabolisms with temperature takes the form:

$\mu_{het} = \mu_{het}^{0ºC} \left(β_{hete}\right)^{T}$

where
- $\mu_{het}^{0ºC}$ is the rate of some heterotrophic metabolism at 0ºC ([s<sup>-1</sup>])
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless])
- $T$ is the in situ temperature of seawater (`Temp(i,j,k)`, [ºC]) 

In the code, the combined term $\left(β_{hete}\right)^{T}$ is saved as `fbc`. See sections below for further details on heterotrophic metabolisms.

**POM --> DOM**
WOMBAT-mid considers the hydrolysation of sinking particulate organic matter (POM) into suspended dissolved organic matter (DOM), which occurs before the remineralisation of the DOM by heterotrophic bacteria. The hydrolysation rate of small sinking organic detritus (`detremi(i,j,k)`, $\Gamma_{sd}^{&rarr; C}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) and large sinking organic detritus (`bdetremi(i,j,k)`, $\Gamma_{ld}^{&rarr; C}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) is computed as:

$\Gamma_{sd}^{&rarr; C} = \Gamma_{sd}^{0ºC} \left(β_{hete}\right)^{T} \left(B_{sd}^{C}\right)^{2}$

$\Gamma_{ld}^{&rarr; C} = \Gamma_{ld}^{0ºC} \left(β_{hete}\right)^{T} \left(B_{ld}^{C}\right)^{2}$

where
- $\Gamma_{sd}^{0ºC} = \Gamma_{ld}^{0ºC}$ is the base hydrolysation rate of sinking detritus at 0ºC (`detlrem`, [(mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>])
- $\left(β_{hete}\right)^{T}$ is the temperature-dependent scaler of heterotrophic metabolism (`fbc`, [dimenionless])
- $B_{sd}^{C}$ and $B_{ld}^{C}$ are the in situ concentrations of small and large sinking organic detritus (`biodet`; `biobdet`, [mmol C m<sup>-3</sup>])

WOMBAT-mid also carries a distinct dissolved organic nitrogen tracer (`f_don(i,j,k)`, $B_{DOM}^{N}$, [mol N kg<sup>-1</sup>]) that receives material during the hydrolysation of particulate organics:

$\Gamma_{sd}^{&rarr; N} = \Gamma_{sd}^{&rarr; B_{DOM}^{C}} \dfrac{16}{122}$

$\Gamma_{ld}^{&rarr; N} = \Gamma_{ld}^{&rarr; B_{DOM}^{C}} \dfrac{16}{122}$

where
- $\dfrac{16}{122}$ is the ratio of nitrogen to carbon in organic material ([mol N (mol C)<sup>-1</sup>])
- $\Gamma_{sd}^{&rarr; C}$ and $\Gamma_{ld}^{&rarr; C}$ are the rates of hydrolysation of small and large particulate organic carbon (`detremi(i,j,k)`; `bdetremi(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])

---


### 4. Light limitation of phytoplankton

Phytoplankton growth is limited by light through a photosynthesis–irradiance (P–I) relationship that links cellular chlorophyll content and photosynthetically available radiation (`radbio`, $PAR$, [W m<sup>-2</sup>]).

First, The initial slope of the P–I curve, (`phy_pisl`, $\alpha_{np}$, [(W m<sup>-2</sup>)<sup>-1</sup>]; `dia_pisl`, $\alpha_{mp}$, [(W m<sup>-2</sup>)<sup>-1</sup>]), determines how efficiently phytoplankton convert light into carbon fixation. It is scaled by the cellular chlorophyll-to-carbon ratio (`phy_chlc`, $Q_{np}^{Chl:C}$, [mol C (mol C)<sup>-1</sup>]; `dia_chlc`, $Q_{mp}^{Chl:C}$, [mol C (mol C)<sup>-1</sup>]).

$\alpha_{np} = \max(\alpha_{np}^{Chl} \cdot Q_{np}^{Chl:C} \ , \ \alpha_{np}^{Chl} \cdot Q_{np}^{-Chl:C})$\
$\alpha_{mp} = \max(\alpha_{mp}^{Chl} \cdot Q_{mp}^{Chl:C} \ , \ \alpha_{mp}^{Chl} \cdot Q_{mp}^{-Chl:C})$ 

where 
- $\alpha_{np}^{Chl}$ is the photosynthetic efficiency per unit chlorophyll in nano-phytoplankton (`alphabio_phy`, [(W m<sup>-2</sup>)<sup>-1</sup> (mol C (mol C)<sup>-1</sup>)<sup>-1</sup>])
- $\alpha_{mp}^{Chl}$ is the photosynthetic efficiency per unit chlorophyll in micro-phytoplankton (`alphabio_dia`, [(W m<sup>-2</sup>)<sup>-1</sup> (mol C (mol C)<sup>-1</sup>)<sup>-1</sup>])
- $Q_{np}^{-Chl:C}$ is the minimum chlorophyll to carbon ratio of nano-phytoplankton cells (`phyminqc`, [mol C (mol C)<sup>-1</sup>])
- $Q_{mp}^{-Chl:C}$ is the minimum chlorophyll to carbon ratio of micro-phytoplankton cells (`diaminqc`, [mol C (mol C)<sup>-1</sup>])
- $Q_{np}^{Chl:C}$ is the in situ chlorophyll to carbon ratio of nano-phytoplankton cells (`phy_chlc`, [mol C (mol C)<sup>-1</sup>])
- $Q_{mp}^{Chl:C}$ is the in situ chlorophyll to carbon ratio of micro-phytoplankton cells (`dia_chlc`, [mol C (mol C)<sup>-1</sup>])

This constraint prevents photosynthesis from collapsing unrealistically at low chlorophyll concentrations. These values are parameter inputs at run time and can differ between nano-phytoplankton and micro-phytoplankton ([Edwards et al., 2015](https://doi.org/10.1002/lno.10033), [Litchman 2022](https://link.springer.com/chapter/10.1007/978-3-030-92499-7_1)).

Second, light limitation (`phy_lpar(i,j,k)`, $L_{np}^{PAR}$), [dimensionless]; `dia_lpar(i,j,k)`, $L_{mp}^{PAR}$), [dimensionless]) is calculated using an exponential P–I formulation.

$L_{np}^{PAR} = 1 - e^{- \alpha_{np} PAR }$ \
$L_{mp}^{PAR} = 1 - e^{- \alpha_{mp} PAR }$

where
- $PAR$ is the downwelling photosynthetically available radiation (`radbio`, [W m<sup>-2</sup>])

At low irradiance ($PAR$), growth increases approximately linearly with light, while at high irradiance photosynthesis asymptotically saturates. We do not account for photoinhibition at very high irradiances.

---


### 5. Realized growth rate of phytoplankton.

Realized growth of nano-phytoplankton (`phy_mu(i,j,k)`, $\mu_{np}$, [s<sup>-1</sup>]) and micro-phytoplankton (`dia_mu(i,j,k)`, $\mu_{mp}$, [s<sup>-1</sup>]) is calculated as:

$\mu_{np} = \mu_{np}^{max} L_{np}^{PAR} \min(L_{np}^{N}, L_{np}^{Fe})$\
$\mu_{mp} = \mu_{mp}^{max} L_{mp}^{PAR} \min(L_{mp}^{N}, L_{mp}^{Fe}) L_{mp}^{Si}$

where
- $\mu_{np}^{max}$ is the maximum potential rate of carbon fixation by nano-phytoplankton (`phy_mumax`, [s<sup>-1</sup>])
- $L_{np}^{PAR}$ is the growth limiter by light of nano-phytoplankton (`phy_lpar(i,j,k)`, [dimensionless])
- $L_{np}^{N}$ is the growth limiter by nitrogen of nano-phytoplankton (`phy_lnit(i,j,k)`, [dimensionless])
- $L_{np}^{Fe}$ is the growth limiter by iron of nano-phytoplankton (`phy_lfer(i,j,k)`, [dimensionless])
- $\mu_{np}^{max}$ is the maximum potential rate of carbon fixation by micro-phytoplankton (`dia_mumax`, [s<sup>-1</sup>])
- $L_{mp}^{PAR}$ is the growth limiter by light of micro-phytoplankton (`dia_lpar(i,j,k)`, [dimensionless])
- $L_{mp}^{N}$ is the growth limiter by nitrogen of micro-phytoplankton (`dia_lnit(i,j,k)`, [dimensionless])
- $L_{mp}^{Fe}$ is the growth limiter by iron of micro-phytoplankton (`dia_lfer(i,j,k)`, [dimensionless])
- $L_{mp}^{Si}$ is the growth limiter by silicic acid of micro-phytoplankton (`dia_lsil(i,j,k)`, [dimensionless])

Liebig's law of the minimum ([Liebig, 1840](https://archive.org/details/organicchemistry00liebrich/mode/2up), [Blackman, 1905](https://doi.org/10.1093/oxfordjournals.aob.a089000)) is applied to resources that are required for biomass synthesis (N and Fe). For micro-phytoplankton, their growth is additionally restricted by silica limitation applied outside of Liebig's law because we treat silica limitation (`dia_lsil(i,j,k)`, $L_{mp}^{Si}$, [dimensionless]) as a structural threshold, rather than as a metabolic throttle (see below).

Carbon fixation by phytoplankton is then calculated as:

$\mu_{np}^{&larr; C} = \mu_{np} B_{np}^{C}$\
$\mu_{mp}^{&larr; C} = \mu_{mp} B_{mp}^{C}$

where 
- $\mu_{np}^{&larr; C}$ is the realized rate of carbon biomass growth by nano-phytoplankton (`phygrow(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $\mu_{mp}^{&larr; C}$ is the realized rate of carbon biomass growth by micro-phytoplankton (`diagrow(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $B_{np}^{C}$ is the in situ concentration of nano-phytoplankton biomass (`f_phy(i,j,k)`, [mol C kg<sup>-1</sup>])
- $B_{mp}^{C}$ is the in situ concentration of micro-phytoplankton biomass (`f_dia(i,j,k)`, [mol C kg<sup>-1</sup>])

---


### 6. Dissolved organic carbon release by phytoplankton.

We implement the overflow hypothesis ([Fogg, 1983](https://doi.org/10.1515/botm.1983.26.1.3); [Hansell & Carlson, 2014](https://books.google.com.au/books?id=7iKOAwAAQBAJ&lpg=PP1&ots=kzkdHuHMF_&dq=Carlson%20Hansell%202014%20doi&lr&pg=PP1#v=onepage&q&f=false)), which posits that phytoplankton can exude their assimilated carbon as dissolved organic carbon (DOC) in high light, low nutrient conditions. We thus account for a phytoplankton-mediated creation of DOC from dissolved inorganic carbon (DIC) via:

$\mu_{np}^{&rarr; DOC} = \min\left( f_{overflow} \mu_{np}^{totalC}, \max\left(0.02 \cdot \mu_{np}^{totalC}, \mu_{np}^{totalC} - \mu_{np}^{&larr; C}\right) \right)$\
$\mu_{mp}^{&rarr; DOC} = \min\left( f_{overflow} \mu_{mp}^{totalC}, \max\left(0.02 \cdot \mu_{mp}^{totalC}, \mu_{mp}^{totalC} - \mu_{mp}^{&larr; C}\right) \right)$

where 
- $\mu_{np}^{&rarr; DOC}$ is the overflow production of DOC by nano-phytoplankton (`phydoc(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $\mu_{mp}^{&rarr; DOC}$ is the overflow production of DOC by micro-phytoplankton (`diadoc(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $\mu_{np}^{totalC}$ is the total carbon fixation rate of nano-phytoplankton (`zval`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $\mu_{mp}^{totalC}$ is the total carbon fixation rate rate of micro-phytoplankton (`zval`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $\mu_{np}^{&larr; C}$ is the realized biomass growth rate of nano-phytoplankton (`phygrow(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $\mu_{mp}^{&larr; C}$ is the realized biomass growth rate of micro-phytoplankton (`diagrow(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $f_{overflow}$ is the maximum fraction total carbon fixation that goes to DOC exudation (`overflow`, [dimenionless])

The total carbon fixation rate of phytoplankton type $p$ is 

$\mu_{p}^{totalC} = \mu_{p}^{&rarr; DOC} + \mu_{p}^{&larr; C} = \mu_{p}^{max} L_{p}^{PAR}$

This formulation is derived from the idea that DOC exudation occurs as a result of the difference between carbon fixation capacity, which is bounded by light, and biosynthesis, which is bounded by light and nutrient resources. Since [Thornton (2014)](https://doi.org/10.1080/09670262.2013.875596) identified that as much as 50% of total phytoplankton carbon fixation can be routed to DOC exudation, we cap DOC exudation at $f_{overflow}$ of total carbon fixation, which is set as to a default of 0.5. We also set a hard bound that 2% of total carbon fixation must at minimum go to DOC production based on the findings of [Bjørnsen (1988)](https://doi.org/10.4319/lo.1988.33.1.0151) who identified that even the healthiest cells lose a small fraction of their assimilated carbon as DOC via passive diffusion across the cell membrane.

---


### 7. Growth of chlorophyll

This step diagnoses the **rate of chlorophyll production** as a function of mixed-layer light, the phytoplankton growth rate and nutrient availability. The structure is consistent with the [Geider, MacIntyre & Kana (1997)](https://doi.org/10.3354/meps148187) formulation that relaxes the chlorophyll-to-carbon ratio towards an optimal value that supports photosynthetic growth under prevailing light and nutrient conditions.

We first solve for the optimal chlorophyll-to-carbon ratio (`phy_chlc`, $Q_{np}^{*Chl:C}$, [mol C (mol C)<sup>-1</sup>]; `dia_chlc`, $Q_{mp}^{*Chl:C}$, [mol C (mol C)<sup>-1</sup>]), which is diagnosed as the ratio required to support maximal photosynthetic carbon fixation under the ambient mean light level in the mixed layer, while accounting for nutrient limitation of biosynthesis:

$Q_{np}^{*Chl:C} = \dfrac{Q_{np}^{+Chl:C}}{1 + \dfrac{\alpha_{np} PAR_{MLD} Q_{np}^{+Chl:C}}{2 \mu_{np}^{max} \min \left(L_{np}^{N}, L_{np}^{Fe} \right) }}$\
$Q_{mp}^{*Chl:C} = \dfrac{Q_{mp}^{+Chl:C}}{1 + \dfrac{\alpha_{mp} PAR_{MLD} Q_{mp}^{+Chl:C}}{2 \mu_{mp}^{max} \min \left(L_{mp}^{N}, L_{mp}^{Fe} \right) }}$

where
- $Q_{np}^{+Chl:C}$ and $Q_{mp}^{+Chl:C}$ are the maximum allowable chlorophyll-to-carbon ratios (`phymaxqc`; `diamaxqc`, [mol C (mol C)<sup>-1</sup>])
- $\alpha_{np}$ and $\alpha_{mp}$ are the chlorophyll-specific initial slopes of the P–I curve (`alphabio_phy`; `alphabio_dia`, [(W m<sup>-2</sup>)<sup>-1</sup> (mol C (mol C)<sup>-1</sup>)<sup>-1</sup>])
- $PAR_{MLD}$ is mean photosynthetically available radiation over the mixed layer (`radmld(i,j,k)`, [W m<sup>-2</sup>])
- $\mu_{np}^{max}$ and $\mu_{mp}^{max}$ are the temperature-dependent maximum phytoplankton growth rates (`phy_mumax(i,j,k)`; `dia_mumax(i,j,k)`, [s<sup>-1</sup>] )
- $L_{np}^{N}$ and $L_{np}^{Fe}$ are the nano-phytoplankton limitation factors for growth on N and Fe (`phy_lnit(i,j,k)`; `phy_lfer(i,j,k)`, [dimensionless])
- $L_{mp}^{N}$ and $L_{mp}^{Fe}$ are the micro-phytoplankton limitation factors for growth on N and Fe (`dia_lnit(i,j,k)`; `dia_lfer(i,j,k)`, [dimensionless])

We set a floor for the minimum chlorophyll-to-carbon ratio of phytoplankton via:

$Q_{np}^{*Chl:C} = \min \left( Q_{np}^{*Chl:C}, Q_{np}^{-Chl:C} \right)$\
$Q_{mp}^{*Chl:C} = \min \left( Q_{mp}^{*Chl:C}, Q_{mp}^{-Chl:C} \right)$

where
- $Q_{np}^{-Chl:C}$ and $Q_{mp}^{-Chl:C}$ are the minimum allowable chlorophyll-to-carbon ratios (`phyminqc`; `diaminqc`, [mol C (mol C)<sup>-1</sup>])

Growth of chlorophyll by nano-phytoplankton and micro-phytoplankton (`pchl_mu(i,j,k)`; `dchl_mu(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) is then calculated as:

$\mu_{np}^{&larr; Chl} = \mu_{np} B_{np}^{Chl} + \dfrac{ Q_{np}^{*Chl:C} - Q_{np}^{Chl:C} }{\tau^{Chl}} \cdot B_{np}^{C}$\
$\mu_{mp}^{&larr; Chl} = \mu_{mp} B_{mp}^{Chl} + \dfrac{ Q_{mp}^{*Chl:C} - Q_{mp}^{Chl:C} }{\tau^{Chl}} \cdot B_{mp}^{C}$

where
- $Q_{np}^{Chl:C}$ and $Q_{mp}^{Chl:C}$ are the in-situ chlorophyll-to-carbon ratios (`phy_chlc`; `dia_chlc`, [mol C (mol C)<sup>-1</sup>])
- $B_{np}^{Chl}$ and $B_{mp}^{Chl}$ are the in-situ concentations of phytoplankton chlorophyll (`f_pchl(i,j,k)`; `f_dchl(i,j,k)`, [mol kg<sup>-1</sup>])
- $\mu_{np}$ and $\mu_{mp}$ are the realized growth rates of phytoplankton (`phy_mu(i,j,k)`; `dia_mu(i,j,k)`, [s<sup>-1</sup>] )
- $\tau^{Chl}$ is the timescale over which chlorophyll synthesis occurs within the cell (`chltau`, [s])
- $B_{np}^{C}$ and $B_{mp}^{C}$ are the in-situ concentations of phytoplankton carbon (`f_phy(i,j,k)`; `f_dia(i,j,k)`, [mol kg<sup>-1</sup>])

This formulation elevates chlorophyll-to-carbon ratios in low light and supresses synthesis when nutrients are low. $\tau^{Chl}$ is an input parameter at run time and should ideally be less than the doubling time of phytplankton given that phytoplankton can internally regulate their chlorophyll stores at rates greater than their overall growth.

---


### 8. Phytoplankton uptake of iron

Like chlorophyll, the iron content of phytoplankton is explicitly tracked as a tracer in WOMBAT-mid. First, a maximum quota is found based on the maximum Fe:C ratio of the phytoplankton type:

$B_{np}^{+Fe} = B_{np}^{C} Q_{np}^{+Fe:C}$\
$B_{mp}^{+Fe} = B_{mp}^{C} Q_{mp}^{+Fe:C}$

where
- $B_{np}^{+Fe}$ and $B_{mp}^{+Fe}$ are the maximum Fe quotas of the nano-phytoplankton and micro-phytoplankton cells (`phy_maxqfe`; `dia_maxqfe`, [mmol Fe m-3])
- $B_{np}^{C}$ and $B_{mp}^{C}$ are the in situ concentrations of nano-phytoplankton and micro-phytoplankton (`biophy`; `biodia`, [mmol C m<sup>-3</sup>])
- $Q_{np}^{+Fe:C}$ and $Q_{mp}^{+Fe:C}$ are the maximum Fe:C ratios of nano-phytoplankton and micro-phytoplankton cells (`phymaxqf`; `diamaxqf`, [mol Fe (mol C)-1])

Following [Aumont et al. (2015)](https://gmd.copernicus.org/articles/8/2465/2015/), this rate is scaled by three terms relating to (i) michaelis-menten type affinity for dFe, (ii) up-regulation of dFe uptake representing investment in transporters when cell quotas are limiting to growth, and (iii) down regulation of dFe uptake associated with enriched cellular quotas.

($i_{np}$) $\dfrac{dFe}{dFe + K_{np}^{Fe}}$ \
($i_{mp}$) $\dfrac{dFe}{dFe + K_{mp}^{Fe}}$

($ii_{np}$) $4 - \dfrac{4.5 L_{np}^{Fe}}{0.5 + L_{np}^{Fe}}$ \
($ii_{mp}$) $4 - \dfrac{4.5 L_{mp}^{Fe}}{0.5 + L_{mp}^{Fe}}$

($iii_{np}$) $\max\left(0, 1 - \dfrac{B_{np}^{Fe} / B_{np}^{+Fe}}{\left|1.05 - B_{np}^{Fe} / B_{np}^{+Fe}\right|} \right)$

($iii_{mp}$) $\max\left(0, 1 - \dfrac{B_{mp}^{Fe} / B_{mp}^{+Fe}}{\left|1.05 - B_{mp}^{Fe} / B_{mp}^{+Fe}\right|} \right)$

where
- $dFe$ is the in situ dissolved iron concentration (`biofer`, [µmol Fe m<sup>-3</sup>])
- $K_{np}^{Fe}$ and $K_{mp}^{Fe}$ are the half-saturation coefficients for dFe uptake by nano-phytoplankton and micro-phytoplankton (`phy_kfe(i,j,k)`; `dia_kfe(i,j,k)`, [µmol Fe <sup>m-3</sup>])
- $L_{np}^{Fe}$ and $L_{mp}^{Fe}$ are the growth limiters of nano-phytoplankton and micro-phytoplankton by iron (`phy_lfer(i,j,k)`; `dia_lfer(i,j,k)`, [dimensionless])
- $B_{np}^{Fe}$ and $B_{mp}^{Fe}$ are the in situ Fe quotas of nano-phytoplankton and micro-phytoplankton cells (`biophyfe`; `biodiafe`, [mmol Fe <sup>m-3</sup>])
- $B_{np}^{+Fe}$ and $B_{mp}^{+Fe}$ are the maximum Fe quotas of nano-phytoplankton and micro-phytoplankton cells (`phy_maxqfe`; `dia_maxqfe`, [mmol Fe <sup>m-3</sup>])

Note that we additionally include a fourth term that decreases the maximum dFe uptake of a cell under light limitation. This is informed by slower uptake of Fe by cells grown in darkness compared to those grown in light by roughly 10-fold ([Strzepek et al., 2025](https://doi.org/10.1093/ismejo/wraf015)), which may be due to physiological stimulation of Fe uptake machinery or photoreduction of ligand-bound iron complexes ([Kong et al., 2023](https://doi.org/10.1002/lno.12331); [Maldonado et al., 2005](https://doi.org/10.1029/2005GB002481)), or possibly a combination of both. To obtain a 10-fold relative increase in Fe uptake rates under light, we applied the following term:

($iv_{np}$) $\max\left(0.01, L_{np}^{PAR}\right)^{0.5}$ \
($iv_{mp}$) $\max\left(0.01, L_{mp}^{PAR}\right)^{0.5}$

where
- $L_{np}^{PAR}$ and $L_{mp}^{PAR}$ are the growth limiters of nano-phytoplankton and micro-phytoplankton by light (`phy_lpar(i,j,k)`; `dia_lpar(i,j,k)`, [dimensionless])

Under very low light, this fourth term reduces maximum potential Fe uptake by 10-fold than what it otherwise would be. All four terms are dimensionless and are designed to scale dissolved iron uptake either up or down. Dissolved iron uptake by nano-phytoplankton and micro-phytoplankton (`phy_dfeupt(i,j,k)`; `dia_dfeupt(i,j,k)`, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]) is then calculated as:

$\mu_{np}^{&larr; dFe} = \mu_{np}^{max} B_{np}^{+Fe} \cdot (i_{np}) \cdot (ii_{np}) \cdot (iii_{np}) \cdot (iv_{np})$ \
$\mu_{mp}^{&larr; dFe} = \mu_{mp}^{max} B_{mp}^{+Fe} \cdot (i_{mp}) \cdot (ii_{mp}) \cdot (iii_{mp}) \cdot (iv_{mp})$

where 
- $\mu_{np}^{&larr; dFe}$ and $\mu_{mp}^{&larr; dFe}$ are the realized uptake rate of dissolved iron by nano-phytoplankton and micro-phytoplankton (`phy_dfeupt(i,j,k)`; `dia_dfeupt(i,j,k)`, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>])
- $\mu_{np}^{max}$ and $\mu_{mp}^{max}$ are the maximum potential growth rates of nano-phytoplankton and micro-phytoplankton (`phy_mumax(i,j,k)`; `dia_mumax(i,j,k)`, [s<sup>-1</sup>])
- $B_{np}^{+Fe}$ and $B_{mp}^{+Fe}$ are the maximum Fe quotas of nano-phytoplankton and micro-phytoplankton cells (`phy_maxqfe`; `dia_maxqfe`, [mol Fe kg<sup>-1</sup>])

---


### 9. Phytoplankton uptake of silicic acid.

Like chlorophyll and iron, the silicon content of micro-phytoplankton is explicitly tracked as a tracer in WOMBAT-mid. Uptake of silicic acid by micro-phytoplankton (`dia_silupt(i,j,k)`, $\mu_{mp}^{&larr; Si}$, [mol Si kg<sup>-1</sup> s<sup>-1</sup>]) is scaled by two terms relating to (i) michaelis-menten type affinity for Si(OH)<sub>4</sub> and (ii) down regulation of Si(OH)<sub>4</sub> uptake associated with enriched cellular quotas. 

(i) 

$\dfrac{Si(OH)<sub>4</sub>}{Si(OH)<sub>4</sub> + K_{mp}^{Si}}$

(ii) 
$\left(\max\left(0.0, \dfrac{Q_{mp}^{Si:C} - Q_{mp}^{-Si:C}}{Q_{mp}^{+Si:C} - Q_{mp}^{-Si:C}} \right)\right)^{0.5}$

where 
- $Si(OH)_{4}$ is the in situ silicic acid concentration (`biosil`, [mmol Si m<sup>-3</sup>])
- $K_{mp}^{Si}$ is the half-saturation coefficient for siliic acid uptake by micro-phytoplankton (`dia_ksi(i,j,k)`, [mmol Si m<sup>-3</sup>])
- $Q_{mp}^{Si:C}$ is the in situ Si:C ratios of micro-phytoplankton cells (`dia_Si2C`, [mol Si (mol C)-1])
- $Q_{mp}^{+Si:C}$ is the maximum Si:C ratios of micro-phytoplankton cells (`diamaxqs`, [mol Si (mol C)-1])
- $Q_{mp}^{-Si:C}$ is the minimum Si:C ratios of micro-phytoplankton cells (`diaminqs`, [mol Si (mol C)-1])

Uptake is then calculated as

$\mu_{mp}^{&larr; Si} = \max \left(0, V_{mp}^{Si} B_{mp}^{C} \cdot (i) \cdot (ii) \right)$

where
- $V_{mp}^{Si}$ is the maximum uptake rate of silicon to carbon by a micro-phytoplankton cell (`diaVmaxs`, [mol Si (mol C)<sup>-1</sup> s<sup>-1</sup>])
- $B_{mp}^{C}$ is the in situ concentration of micro-phytoplankton carbon biomass (`f_dia(i,j,k)`, [mol C kg<sup>-1</sup>])

Unlike iron uptake, we do not include upregulation terms for silicic acid uptake. This is on the basis that highly silicified diatoms are caused by slow growth rather than increased/luxury uptake. Both light-limited and iron-limited diatoms show increases in their Si:C content by roughly 3-fold and it is suggested that this due to decoupling of biogenic silica precipitation from slowing carbon fixation ([Liu et al., 2016](https://doi.org/10.3389/fmars.2016.00089); [Hutchins & Bruland 1998](https://doi.org/10.1038/31203); [Takeda 1998](https://doi.org/10.1038/31674)). To properly decouple silicification from biomass growth we therefore make $V_{mp}^{Si}$ temperature-independent, which ensures that polar diatoms have a tendency towards heavier silicification then tropical diatoms ([Baines et al., 2010](https://doi.org/10.1029/2010GB003856)). 

---


### 10. Iron chemistry (scavenging, coagulation, dissolution).

Treatment of dissolved iron (`f_fe(i,j,k)`, $dFe$, mol kg<sup>-1</sup>) follows a combination of [Aumont et al. (2015)](https://gmd.copernicus.org/articles/8/2465/2015/) and [Tagliabue et al. (2023)](https://www.nature.com/articles/s41586-023-06210-5). Our calculations involve:
1. Solving for the distinct pools of dissolved iron: free iron, ligand-bound iron and colloidal iron.
2. Computing scavenging of free iron to authigenic sinking phases.
3. Computing coagulation of colloidal iron to authigenic sinking phases.
4. Computing dissolution of authigenic sinking phases back to dissolved iron.

We first estimate the **solubility of free Fe from Fe<sup>3+</sup>** in solution using temperature, pH and salinity using the thermodynamic equilibrium equations of [Liu & Millero (2002)](https://www.sciencedirect.com/science/article/abs/pii/S0304420301000743). 

$T_K = \max(5.0, T) + 273.15$ \
$T_K^{-1} = \dfrac{1}{T_K}$ \
$I_{S} = \dfrac{19.924,S}{1000 - 1.005,S}$ \

Solubility constants:

$Fe_{sol1} = 10^{\left(-13.486 - 0.1856\sqrt{I_S} + 0.3073 I_S + 5254,T_K^{-1}\right)}$ \
$Fe_{sol2} = 10^{\left(2.517 - 0.8885\sqrt{I_S} + 0.2139 I_S - 1320,T_K^{-1}\right)}$ \
$Fe_{sol3} = 10^{\left(0.4511 - 0.3305\sqrt{I_S} - 1996,T_K^{-1}\right)}$ \
$Fe_{sol4} = 10^{\left(-0.2965 - 0.7881\sqrt{I_S} - 4086,T_K^{-1}\right)}$ \
$Fe_{sol5} = 10^{\left(4.4466 - 0.8505\sqrt{I_S} - 7980,T_K^{-1}\right)}$

Final Fe(III) solubility:

$Fe_{sol} = Fe_{sol1}\left([H^+]^3 + Fe_{sol2}[H^+]^2 + Fe_{sol3}[H^+] + Fe_{sol4} + \dfrac{Fe_{sol5}}{[H^+]}\right)\times10^{9}$

where 
- $T_{K}$ is in situ water temperature (`ztemk`, [ºK])
- $I_{S}$ is a salinity coefficient (`zval`, [dimenionless])
- $[H^+]$ is in situ hydrogen ion concentration (`hp`, [mol L<sup>-1</sup>])
- $Fe_{sol}$ is the final estimated solubility of dissolved iron in seawater (`fe3sol`, [nmol Fe kg<sup>-1</sup>]).

Next we **estimate the concentration of colloidal iron** in solution following [Tagliabue et al. 2023](https://www.nature.com/articles/s41586-023-06210-5). Colloidal dissolved Fe (`fecol(i,j,k)`, $dFe_{col}$, [mmol Fe m<sup>-3</sup>]) is whatever exceeds the inorganic solubility ceiling (`fe3sol`, $dFe_{sol}$, [mmol Fe m<sup>-3</sup>]), but we enforce a hard minimum that colloids are at least 10% of total dissolved Fe (`biofer`, $dFe$, [mmol Fe m<sup>-3</sup>]).

$dFe_{col} = \max\left(0.1 dFe\, \ dFe - dFe_{sol} \right)$

Following solving for colloidal Fe, we **partition the remaining dissolved Fe into ligand-bound and free iron**. To do so, we find the remaining dissolved iron not in colloidal form (`fe_sfe`, $dFe_{sFe}$, [mmol Fe m<sup>-3</sup>]), 

$dFe_{sFe} = \max\left(0.0,\ dFe - dFe_{col} \right)$

Partitioning is done using a standard quadratic form that drops out of mass balance + equilibrium for 1:1 complexation with a single ligand class. To do so, we need the temperature- and light-dependent conditional stability constant for Fe–ligand complexation (`fe_keq`, $Fe_{Keq}$, [nmol Fe kg<sup>-1</sup>]). The temperature dependency comes from [Volker & Tagliabue (2015)](https://doi.org/10.1016/j.marchem.2014.11.008). The light-dependency accounts for the photoreduction of photoreactive ligands, which was identified to reduce the conditional stability constant of aquachelin by 0.7 log10 units ([Barbeau et al., 2001](https://doi.org/10.1038/35096545); [Vraspir & Butler, 2009](https://doi.org/10.1146/annurev.marine.010908.163712)):

$Fe_{Keq} = 10^{ \left(17.27 - 1565.7 T_K^{-1} \right) }\times 10^{-9}$, 

After finding $Fe_{Keq}$ we solve for the free dissolved Fe concentration (`feIII`, $dFe_{free}$, [nmol Fe kg<sup>-1</sup>]):

$z = 1.0 + [Ligand] \cdot Fe_{Keq} - dFe_{sFe}\cdot Fe_{Keq}$ \
$Fe_{free} = \dfrac{-z + \sqrt{z^2 + 4.0 Fe_{Keq} dFe_{sFe}}}{2 Fe_{Keq} + \varepsilon}$
$Fe_{free} = \max\left(0,\ \min(dFe_{free}, dFe_{sFe})\right)$

where
- $[Ligand]$ is the in situ concentration of iron-binding ligands (`ligand`, [nmol kg<sup>-1</sup>]).

Whatever soluble dissolved iron is not present as inorganic free iron is assigned to ligand-bound dissolved iron:

$dFe_{lig} = dFe_{sFe} - dFe_{free}$

Now that we have separated the dissolved Fe pool into its subcomponents of free, ligand-bound and colloidal Fe, we solve for scavenging of free iron and coagulation of colloidal, both of which remove dissolved iron and transfer these to two sinking authigenic particles. These authigenic sinking particles included a small, slowly sinking type (`f_afe(i,j,k)`, $Fe_{sA}$, [mol Fe kg<sup>-1</sup>]) and a large, fast sinking type (`f_bafe(i,j,k)`, $Fe_{lA}$, [mol Fe kg<sup>-1</sup>]). Both scavenging and colloidal coagulation are the major sinks of dissolved iron outside of phytoplankton uptake and this dissolved iron is transferred to the authigenics.

**Scavenging:** \
Scavenging of dissolved iron specifically affects free iron, is accelerated by the presence of particles in the water column and routes this iron to two sinking authigenic phases. Total scavenging of dissolved iron (`fescaven(i,j,k)`, $Sc_{dFe}^{&rarr;}$, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]) is calculated as

$Sc_{dFe}^{&rarr;} = dFe_{free} \left(10^{-7} + \gamma_{dFe}^{scav} \cdot B_{particles} \right)$ \

where
- $dFe_{free}$ is the in situ concentration of dissolved free iron (`feIII(i,j,k)`, [nmol Fe kg<sup>-1</sup>])
- $\gamma_{dFe}^{scav}$ is the rate constant of scavenging (`kscav_dfe`, [(mmol m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>])
- $B_{particles}^{M}$ is the in situ concentration of detrital particles in the water column (`partic`, [mmol m<sup>-3</sup>])

$B_{particles}^{M} = 2 \cdot B_{sd}^{C} + 2 \cdot B_{ld}^{C} + 2 \cdot B_{ld}^{Si} + 8.3 \cdot B_{CaCO_3}^{C}$

where
- $B_{sd}^{C}$ is the in situ concentration of small organic carbon detritus (`biodet`, [mmol C m<sup>-3</sup>])
- $B_{ld}^{C}$ is the in situ concentration of large organic carbon detritus (`biobdet`, [mmol C m<sup>-3</sup>])
- $B_{ld}^{Si}$ is the in situ concentration of biogenic silica detritus (`biodet`, [mmol Si m<sup>-3</sup>])
- $B_{CaCO_3}^{C}$ is the in situ concentration of calcium carbonate detritus (`biocaco3`, [mmol C m<sup>-3</sup>])

Organic carbon-based particle types $B_{sd}^{C}$ and $B_{ld}^{C}$ are multipled by 2 assuming that carbon represents half the mass of the particle, $B_{ld}^{Si}$ is multipled by 2 assuming that it represents biogenic silica with a molecular mass of 60 g mol<sup>-1</sup>, and inorganic carbon-based particles $B_{CaCO_3}^{C}$ is multipled by 8.3 since the moleculate weight of calcium carbonate is 100 g mol<sup>-1</sup>. 

Total scavenging ($Sc_{dFe}^{&rarr;}$) of free iron is then broken into two parts: scavenging to small authigenic particles (`fescaafe(i,j,k)`, $Sc_{dFe}^{&rarr; Fe_{sA}}$, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]) and scavenging to large authigenic particles (`fescabafe(i,j,k)`, $Sc_{dFe}^{&rarr; Fe_{lA}}$, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]). 

$Sc_{dFe}^{&rarr; Fe_{sA}} = Sc_{dFe}^{&rarr;} \cdot \dfrac{ 2 \cdot B_{sd}^{C} + 8.3 \cdot B_{CaCO_3}^{C} }{ B_{particles}^{M} }$\
$Sc_{dFe}^{&rarr; Fe_{lA}} = Sc_{dFe}^{&rarr;} \cdot \dfrac{ 2 \cdot B_{ld}^{C} + 2 \cdot B_{ld}^{Si} }{ B_{particles}^{M} }$


**Coagulation:** \
Similarly to scavenging of free iron, coagulation routes dissolved iron to two sinking authigenic phases. However, coagulation acts on the colloidal fraction of dissolved iron ([Tagliabue et al., 2023](https://www.nature.com/articles/s41586-023-06210-5)). Rates of coagulation of colloidal iron to small, slowly sinking authigenic iron (`fecoag2afe(i,j,k)`, $Co_{dFe}^{&rarr; Fe_{sA}}$, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]) and large, fast sinking authigenic iron (`fecoag2bafe(i,j,k)`, $Co_{dFe}^{&rarr; Fe_{lA}}$, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]) follow the form:

$Co_{dFe}^{&rarr; Fe_{sA}} = dFe_{col} \gamma_{dFe}^{coag} \cdot S_{coag}^{sA} $\
$Co_{dFe}^{&rarr; Fe_{lA}} = dFe_{col} \gamma_{dFe}^{coag} \cdot S_{coag}^{lA} $

where
- $dFe_{col}$ is the in situ concentration of dissolved colloidal iron (`fecol(i,j,k)`, [mol Fe kg<sup>-1</sup>])
- $\gamma_{dFe}^{coag}$ is the iron coagulation rate constant (`kcoag_dfe`, [(mmol m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>])
- $S_{coag}^{sA}$ and $S_{coag}^{lA}$ are scaling coefficients to decelerate or accelerate coagulation of small and large particles (`zval`, [mmol C m<sup>-3</sup>])

The coagulation scaling coefficients are themselves dependent on the concentrations of dissolved organic carbon, particulate organic carbon, phytoplankton biomass and the rate of mixing. For small particle coagulation:

$S_{coag}^{sA} = H_{mix} \left(12 \cdot F_{coag} B_{DOM}^{C} + 9 \cdot B_{sd}^{C}\right) + 2.5 \cdot B_{sd}^{C} + 128 \cdot F_{coag} B_{DOC}^{C} + 725 \cdot B_{sd}^{C}$ \
$F_{coag} = \dfrac{B_{np}^{C} + B_{mp}^{C}}{B_{np}^{C} + B_{mp}^{C} + 0.03}$ 

where
- $H_{mix}$ is a heavyside step function that is equalt to 1 in the mixed layer and 0.01 beneath the mixed layer (`shear`, [dimensionless])
- $F_{coag}$ is a phytoplankton concentration dependent coagulation factor (`biof`, [dimensionless]) 
- $B_{np}^{C}$ and $B_{mp}^{C}$ are the concentrations of nano- and micro-phytoplankton biomass (`biophy`; `biodia`, [mmol C m<sup>-3</sup>]) 
- $B_{DOM}^{C}$ is the concentration of dissolved organic matter in carbon (`biodoc`, [mmol C m<sup>-3</sup>])
- $B_{sd}^{C}$ is the concentration of small organic detrital particles (`biodet`, [mmol C m<sup>-3</sup>])

For large particle coagulation:

$S_{coag}^{lA} = \left(2 \cdot H_{mix} + 1.37 \right) B_{ld}^{C} + 1.94 \cdot B_{ld}^{C}$

where
- $H_{mix}$ is a heavyside step function that is equalt to 1 in the mixed layer and 0.01 beneath the mixed layer (`shear`, [dimensionless])
- $B_{ld}^{C}$ is the concentration of large organic detrital particles (`biobdet`, [mmol C m<sup>-3</sup>])

Together, these terms implement a biologically mediated coagulation pathway in which iron removal from the dissolved pool is tightly coupled to ecosystem state. The formulation reflects the central conclusion of [Tagliabue et al. (2023)](https://www.nature.com/articles/s41586-023-06210-5): that iron cycling is not governed solely by inorganic chemistry, but is strongly regulated by biological activity, organic matter dynamics, and particle ecology across the upper ocean. 


**Dissolution:** \
Small, slow sinking authigenic (`f_afe(i,j,k)`, $Fe_{sA}$, [mol Fe kg<sup>-1</sup>]) and a large, fast sinking authigenic iron (`f_bafe(i,j,k)`, $Fe_{lA}$, [mol Fe kg<sup>-1</sup>]) are returned back to the dissolved iron phase through reductive processes or complexation with ligands ([Tagliabue et al., 2023](https://www.nature.com/articles/s41586-023-06210-5)). We represent this process simply via dissolution rate cofficients:

$D_{sA}^{&rarr; dFe} = Fe_{sA} \gamma_{sA}^{diss}$\
$D_{lA}^{&rarr; dFe} = Fe_{lA} \gamma_{lA}^{diss}$

---


### 11. Biogenic silica dissolution. 

**Silicic acid equilibrium concentration**\
To determine the rate of biogenic silica dissolution we must first determine the equilibrium concentration of silicic acid ($Si(OH)_{4}$) in seawater. To do so, we solve for this equilibrium concentration via thermodynamic first-principles:

$K_{Si(OH)_{4}}(T,P) = \dfrac{\gamma_{Si(OH)_{4}^{0}} \cdot [Si(OH)_{4}]^{eq}}{(a_{H_{2}O})^{2}}$

where
- $K_{Si(OH)_{4}}(T,P)$ is the thermodynamic equilibrium constant in seawater at a given temperature and pressure (`K_am_silica`, [mol Si kg<sup>-1</sup>])
- $\gamma_{Si(OH)_{4}^{0}}$ is the activity ratio of $Si(OH)_{4}$ in seawater (`gamma0`, [dimensionless])
- $[Si(OH)_{4}]^{eq}$ is the equilibrium concentration of $Si(OH)_{4}$ (`sileqc(i,j,k)`, [mol Si kg<sup>-1</sup>])
- $a_{H_{2}O}$ is the activity of seawater (`alphaH2O`, [dimensionless])

The equation is rearranged such that:

$[Si(OH)_{4}]^{eq} = \dfrac{K_{Si(OH)_{4}}(T,P) \cdot (a_{H_{2}O})^{2}}{\gamma_{Si(OH)_{4}^{0}}}$

The activity of seawater is slightly less than 1 due to dissolved salts lowering its chemical potential and so we set $a_{H_{2}O}$ equal to 0.999 ([IOC, SCOR & IAPSO, 2010](https://www.teos-10.org/pubs/TEOS-10_Manual.pdf)). For $\gamma_{Si(OH)_{4}^{0}}$ we follow [Savenko 2014](https://doi.org/10.1134/S0001437014020222) who demonstrated that the solubility of $Si(OH)_{4}$ decreases predictably with salinity according to

$\gamma_{Si(OH)_{4}^{0}} = 1 + 0.0053 \cdot S - 0.000034 \cdot S^{2}$

where
- $S$ is the in situ salinity of seawater (`Salt(i,j,k)`, [psu])

For $K_{Si(OH)_{4}}(T,P)$ we follow the derivation of [Gunnarsson & Arnórsson (2000)](https://doi.org/10.1016/S0016-7037(99)00426-3) who relate the thermodynamic equilibrium constant of $Si(OH)_{4}$ to variations in temperature at a constant pressure of 1 bar ($P^{1}$):

$K(T,P^{1}) = 10^{ -8.476 - \dfrac{485.24}{T_{K}} - 2.268 \times 10^{-6} \cdot (T_{K})^{2} + 3.068 log_{10}(T_{K})}$

where
- $T_{K}$ is the in situ temperature of seawater ([`zval`, [ºK]])

We add a classic pressure correction to $K(T,P^{1})$ to retrieve $K(T,P)$ of the form:

$K(T,P) = K(T,P^{1}) \cdot e^{- \dfrac{\Delta V^{0}}{RT_{K}}P}$

where
- $\Delta V^{0}$ is the partial molal volume change (`deltaV0`, [m<sup>3</sup> mol<sup>-1</sup>])
- $R$ is the universal gas constant (`Rgas`, [J ºK<sup>-1</sup> mol<sup>-1</sup>])
- $T_{K}$ is the in situ temperature of seawater ([`zval`, [ºK]])
- $P$ is the in situ pressure [`zm(i,j,k) * 1.0e4`, [Pa]]

We obtain $\Delta V^{0}$ from [Willey, 1982](https://doi.org/10.1016/0016-7037(82)90015-1) and [Loucaides et al., 2012](https://doi.org/10.1016/j.marchem.2012.04.003) of roughly -9.0 [cm<sup>3</sup> mol<sup>-1</sup>], which we convert to [m<sup>3</sup> mol<sup>-1</sup>] by multiplying by $10^{6}$. The negative value of $\Delta V^{0}$ implies an increase in dissolution of silica at higher pressures. These value return equilibrium concentrations of silicic acid on the order of 1000 to 1800 mmol m<sup>-3</sup>. Temperature increases are the largest control, while pressure increases from the surface to the ocean bottom increase solubility by 15-20%. 


**Biogenic silica dissolution**\
Biogenic silica is only considered associated with the large type of sinking particulate organic matter and dissolution (`bsidiss(i,j,k)`, $D_{B_{ld}^{Si}}^{&rarr; Si}$, [mol Si kg<sup>-1</sup> s<sup>-1</sup>]) occurs via

$D_{B_{ld}^{Si}}^{&rarr; Si} = d_{B_{ld}^{Si}} \cdot B_{ld}^{Si}$

where
- $d_{B_{ld}^{Si}}$ is the rate of biogenic silica dissolution (`disssi(i,j,k)`, [s<sup>-1</sup>])
- $B_{ld}^{Si}$ is the in situ concentration of biogenic silica (`f_bdetsi(i,j,k)`, [mol Si kg<sup>-1</sup>])

We treat the dissolution rate of biogenic silica ($d_{B_{ld}^{Si}}$) as dependent on three conditions: the degree of undersaturation ([Rickert, 2000](https://epic.awi.de/id/eprint/26530/1/BerPolarforsch2000351.pdf); [Van Cappellen & Qiu, 1997](https://doi.org/10.1016/S0967-0645(96)00112-9); [Van Cappellen et al., 2002](https://doi.org/10.1029/2001GB001431)), in situ temperature ([Kamatani, 1982](https://doi.org/10.1007/BF00393146); [Greenwood et al., 2005](https://doi.org/10.1007/s10498-004-9515-y)) and the activity of heterotrophic microbes ([Bidle & Azam, 1999](https://doi.org/10.1038/17351); [Bidle et al., 2003](https://doi.org/10.4319/lo.2003.48.5.1855)). To account for these conditions, we formulate the rate of dissolution of biogenic silica (`disssi(i,j,k)`, $d_{Si}$, [s<sup>-1</sup>]) as the product of three terms:

$d_{B_{ld}^{Si}} = d_{B_{ld}^{Si}}^{T} \cdot S_{B_{ld}^{Si}}^{Sat} \cdot S_{B_{ld}^{Si}}^{bio}$

where
- $d_{B_{ld}^{Si}}^{T}$ is the temperature-dependent rate of dissolution (`disssi_temp`, [s<sup>-1</sup>])
- $S_{B_{ld}^{Si}}^{Sat}$ is a scaling factor that decelerates dissolution as the in situ concentration approachs the equilibrium concentration (`disssi_usat`, [dimenionless])
- $S_{B_{ld}^{Si}}^{bio}$ is a scaling factor that accelerates dissolution in the presence of heterotrophic bacterial biomass (`disssi_bact`, [dimenionless])

First, we solve for $d_{B_{ld}^{Si}}^{T}$. [Kamatani, 1982](https://doi.org/10.1007/BF00393146) measured dissolution rates of biogenic silica collected in Tokyo Bay between 8ºC and 30ºC and identified that these roughly obeyed the equation:

$d_{B_{ld}^{Si}}^{T} = \dfrac{e^{\alpha + β T}}{3600}$

where 
- $\alpha$ is a species-dependent dissolution intercept that ranges between -7.35 and -10.38. We set $\alpha$ = -8.0
- $β$ is the slope common to all species and is equal to 0.0833
- $T$ is the in situ temperature of seawater ([`Temp(i,j,k)`, [ºC]])
- $3600$ converts the rate from [hour<sup>-1</sup>] to [s<sup>-1</sup>]

Next, we apply scaling terms that either decelerate or accelerate dissolution. Given that equilibrium concentrations of Si(OH)<sub>4</sub> vary between 1000 to 1800 mmol m<sup>-3</sup> in the ocean, while actual in situ concentrations rarely exceed 200 mmol m<sup>-3</sup>, $Si(OH)_{4}$ is always undersaturated. We therefore assume that $Si(OH)_{4}$ is highly undersaturated everywhere in the ocean. According to [Van Cappellen et al., 2002](https://doi.org/10.1029/2001GB001431) "Detailed kinetic studies of biogenic silica dissolution conducted in flow-through reactors demonstrate that at very high degrees of undersaturation the dissolution kinetics switch from a linear dependence on the degree of undersaturation to an exponential one". Hence, we apply equation 2.13 from [Rickert, 2000](https://epic.awi.de/id/eprint/26530/1/BerPolarforsch2000351.pdf):

$S_{B_{ld}^{Si}}^{Sat} = 1 - \left(\dfrac{[Si(OH)_{4}]}{[Si(OH)_{4}]^{eq}}\right)^{2}$
 
where
- $[Si(OH)_{4}]$ is the in situ concentration of $Si(OH)_{4}$ (`f_sil(i,j,k)`, [mol Si kg<sup>-1</sup>])
- $[Si(OH)_{4}]^{eq}$ is the equilibrium concentration of $Si(OH)_{4}$ (`sileqc(i,j,k)`, [mol Si kg<sup>-1</sup>])
- we use an exponent of $2$ informed by the organic carbon-rich biogenic silica dissolution kinetics reported in Table 3.4 of [Rickert, 2000](https://epic.awi.de/id/eprint/26530/1/BerPolarforsch2000351.pdf)

The scaling term associated with activity of heterotrophic bacteria is informed by substantial evidence. According to [Rickert et al., 2002](https://doi.org/10.1016/S0016-7037(01)00757-8) "The removal of organic or inorganic coatings enhance the reactivity by at least an order of magnitude". Order of magnitude increases in silica dissolution have been reported for diatom frustules in contact with bacteria ([Bidle & Azam, 1999](https://doi.org/10.1038/17351)), while anti-biotic treatments to mesocosms off Monterey Bay caused silica dissolution to reduced by 50% ([Bidle et al., 2003](https://doi.org/10.4319/lo.2003.48.5.1855)). We represent this bacterially-induced stimulation of dissolution with

$S_{B_{ld}^{Si}}^{bio} = 1 + F_{B_{ld}^{Si}}^{bac} \cdot \dfrac{B_{bac}^{C}}{B_{bac}^{C} + K_{B_{ld}^{Si}}^{bac}} $

where
- $F_{B_{ld}^{Si}}^{bac}$ is the factor increase in dissolution caused by peak bacterial biomass (`bsi_fbac`, [dimenionless])
- $K_{B_{ld}^{Si}}^{bac}$ is the half-saturation coefficient for stimulation of silica dissolution in the presence of bacterial biomass (`bsi_kbac`, [mmol C m<sup>-1</sup>])
- $B_{bac}^{C}$ is the in situ concentration of bacterial biomass (`biobac1` + `biobac2`, [mmol C m<sup>-1</sup>])


---

### 12. Mortality terms

Mortality of ecological functional types are affected by both linear ($\gamma$) and quadratic ($\Gamma$) terms. Linear terms are per-capita losses associated with the costs of basal metabolism. Quadratic, and thus density-dependent losses, are associated with disease, aggregation and coagulation, viruses, infection and canabalism. None of these processes are represented explicitly within the model, so we represent them implicitly.

**Linear losses** of nano-phytoplankton ($_{np}$), micro-phytoplankton ($_{mp}$), micro-zooplankton ($_{mz}$), meso-zooplankton ($_{Mz}$), facultative nitrate-reducing bacteria ($_{b1}$), facultative nitrous oxide-reducing bacteria ($_{b2}$) and ammonia oxidizing archaea ($_{aoa}$) in [mol kg<sup>-1</sup> s<sup>-1</sup>] are modelled as

$\gamma_{np}^{&rarr; C} = \gamma_{np}^{0ºC} (β_{hete})^{T} B_{np}^{C}$\
$\gamma_{mp}^{&rarr; C} = \gamma_{mp}^{0ºC} (β_{hete})^{T} B_{mp}^{C}$\
$\gamma_{mz}^{&rarr; C} = \gamma_{mz}^{0ºC} (β_{hete})^{T} F_{mz}^{\gamma} B_{mz}^{C}$\
$\gamma_{Mz}^{&rarr; C} = \gamma_{Mz}^{0ºC} (β_{hete})^{T} F_{Mz}^{\gamma} B_{Mz}^{C}$\
$\gamma_{b1}^{&rarr; C} = \gamma_{b1}^{0ºC} (β_{hete})^{T} B_{b1}^{C}$\
$\gamma_{b2}^{&rarr; C} = \gamma_{b2}^{0ºC} (β_{hete})^{T} B_{b2}^{C}$\
$\gamma_{aoa}^{&rarr; C} = \gamma_{aoa}^{0ºC} (β_{hete})^{T} B_{aoa}^{C}$

where
- $\gamma_{np}^{0ºC}$ is the rate of linear mortality of nano-phytoplankton at 0ºC (`phylmor`, [s<sup>-1</sup>])
- $\gamma_{mp}^{0ºC}$ is the rate of linear mortality of micro-phytoplankton at 0ºC (`dialmor`, [s<sup>-1</sup>])
- $\gamma_{mz}^{0ºC}$ is the rate of linear mortality of micro-zooplankton at 0ºC (`zoolmor`, [s<sup>-1</sup>])
- $\gamma_{Mz}^{0ºC}$ is the rate of linear mortality of meso-zooplankton at 0ºC (`meslmor`, [s<sup>-1</sup>])
- $\gamma_{b1}^{0ºC}$ is the rate of linear mortality of facultative $NO_{3}$-reducing bacteria at 0ºC (`bac1lmor`, [s<sup>-1</sup>])
- $\gamma_{b2}^{0ºC}$ is the rate of linear mortality of facultative $N_{2}O$-reducing bacteria at 0ºC (`bac2lmor`, [s<sup>-1</sup>])
- $\gamma_{aoa}^{0ºC}$ is the rate of linear mortality of ammonia oxidizing archaea at 0ºC (`aoalmor`, [s<sup>-1</sup>])
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless])
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC])
- $B_{np}^{C}$ is the concentration of nano-phytoplankton carbon biomass (`f_phy(i,j,k)`, [mol kg<sup>-1</sup>])
- $B_{mp}^{C}$ is the concentration of micro-phytoplankton carbon biomass (`f_dia(i,j,k)`, [mol kg<sup>-1</sup>])
- $B_{mz}^{C}$ is the concentration of micro-zooplankton carbon biomass (`f_zoo(i,j,k)`, [mol kg<sup>-1</sup>])
- $B_{Mz}^{C}$ is the concentration of meso-zooplankton carbon biomass (`f_mes(i,j,k)`, [mol kg<sup>-1</sup>])
- $B_{b1}^{C}$ is the concentration of facultative $NO_{3}$-reducing bacteria carbon biomass (`f_bac1(i,j,k)`, [mol kg<sup>-1</sup>])
- $B_{b2}^{C}$ is the concentration of facultative $N_{2}$O-reducing bacteria carbon biomass (`f_bac2(i,j,k)`, [mol kg<sup>-1</sup>])
- $B_{aoa}^{C}$ is the concentration of ammonia oxidizing archaea carbon biomass (`f_aoa(i,j,k)`, [mol kg<sup>-1</sup>])
- $S_{mz}^{\gamma}$ is a scaling factor that reduces micro-zooplankton linear mortality at low biomass (`zoo_slmor`, [dimenionless])
- $S_{Mz}^{\gamma}$ is a scaling factor that reduces meso-zooplankton linear mortality at low biomass (`mes_slmor`, [dimenionless])

In the above, we scale down **linear mortality** of micro- and meso-zooplannkton when their populations are very small small with

$S_{mz}^{\gamma} = \dfrac{B_{mz}^{C}}{B_{mz}^{C} + K_{mz}^{\gamma}}$\
$S_{Mz}^{\gamma} = \dfrac{B_{Mz}^{C}}{B_{Mz}^{C} + K_{Mz}^{\gamma}}$

where
- $B_{mz}^{C}$ is the concentration of micro-zooplankton carbon biomass (`biozoo`, [mmol C m<sup>-1</sup>])
- $B_{Mz}^{C}$ is the concentration of meso-zooplankton carbon biomass (`biomes`, [mmol C m<sup>-1</sup>])
- $K_{mz}^{\gamma}$ is the half-saturation coefficient for scaling down linear mortality losses for micro-zooplankton (`zookz`, [mmol C m<sup>-1</sup>])
- $K_{Mz}^{\gamma}$ is the half-saturation coefficient for scaling down linear mortality losses for meso-zooplankton (`meskz`, [mmol C m<sup>-1</sup>])


**Quadratic losses** of nano-phytoplankton ($_{np}$), micro-phytoplankton ($_{mp}$), micro-zooplankton ($_{mz}$), meso-zooplankton ($_{Mz}$), facultative nitrate-reducing bacteria ($_{b1}$), facultative nitrous oxide-reducing bacteria ($_{b2}$) and ammonia oxidizing archaea ($_{aoa}$) in [mol kg<sup>-1</sup> s<sup>-1</sup>] are modelled as

$\Gamma_{np}^{&rarr; C} = \Gamma_{np}^{0ºC} (β_{hete})^{T} B_{np}^{C}$\
$\Gamma_{mp}^{&rarr; C} = \Gamma_{mp}^{0ºC} (β_{hete})^{T} B_{mp}^{C}$\
$\Gamma_{mz}^{&rarr; C} = \Gamma_{mz}^{0ºC} (β_{hete})^{T} F_{mz}^{\gamma} B_{mz}^{C}$\
$\Gamma_{Mz}^{&rarr; C} = \Gamma_{Mz}^{0ºC} (β_{hete})^{T} F_{Mz}^{\gamma} B_{Mz}^{C}$\
$\Gamma_{b1}^{&rarr; C} = \Gamma_{b1}^{0ºC} (β_{hete})^{T} B_{b1}^{C}$\
$\Gamma_{b2}^{&rarr; C} = \Gamma_{b2}^{0ºC} (β_{hete})^{T} B_{b2}^{C}$\
$\Gamma_{aoa}^{&rarr; C} = \Gamma_{aoa}^{0ºC} (β_{hete})^{T} B_{aoa}^{C}$

where
- $\Gamma_{np}^{0ºC}$ is the rate of quadratic mortality of nano-phytoplankton at 0ºC (`phyqmor`, [(mol C kg<sup>-1</sup>)<sup>-1</sup> s<sup>-1</sup>])
- $\Gamma_{mp}^{0ºC}$ is the rate of quadratic mortality of micro-phytoplankton at 0ºC (`diaqmor`, [(mol C kg<sup>-1</sup>)<sup>-1</sup> s<sup>-1</sup>])
- $\Gamma_{mz}^{0ºC}$ is the rate of quadratic mortality of micro-zooplankton at 0ºC (`zooqmor`, [(mol C kg<sup>-1</sup>)<sup>-1</sup> s<sup>-1</sup>])
- $\Gamma_{Mz}^{0ºC}$ is the rate of quadratic mortality of meso-zooplankton at 0ºC (`mesqmor`, [(mol C kg<sup>-1</sup>)<sup>-1</sup> s<sup>-1</sup>])
- $\Gamma_{b1}^{0ºC}$ is the rate of quadratic mortality of facultative $NO_{3}$-reducing bacteria at 0ºC (`bac1qmor`, [(mol C kg<sup>-1</sup>)<sup>-1</sup> s<sup>-1</sup>])
- $\Gamma_{b2}^{0ºC}$ is the rate of quadratic mortality of facultative $N_{2}O$-reducing bacteria at 0ºC (`bac2qmor`, [(mol C kg<sup>-1</sup>)<sup>-1</sup> s<sup>-1</sup>])
- $\Gamma_{aoa}^{0ºC}$ is the rate of quadratic mortality of ammonia oxidizing archaea at 0ºC (`aoaqmor`, [(mol C kg<sup>-1</sup>)<sup>-1</sup> s<sup>-1</sup>])
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless])
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC])
- $B_{np}^{C}$ is the concentration of nano-phytoplankton carbon biomass (`f_phy(i,j,k)`, [mol kg<sup>-1</sup>])
- $B_{mp}^{C}$ is the concentration of micro-phytoplankton carbon biomass (`f_dia(i,j,k)`, [mol kg<sup>-1</sup>])
- $B_{mz}^{C}$ is the concentration of micro-zooplankton carbon biomass (`f_zoo(i,j,k)`, [mol kg<sup>-1</sup>])
- $B_{Mz}^{C}$ is the concentration of meso-zooplankton carbon biomass (`f_mes(i,j,k)`, [mol kg<sup>-1</sup>])
- $B_{b1}^{C}$ is the concentration of facultative $NO_{3}$-reducing bacteria carbon biomass (`f_bac1(i,j,k)`, [mol kg<sup>-1</sup>])
- $B_{b2}^{C}$ is the concentration of facultative $N_{2}$O-reducing bacteria carbon biomass (`f_bac2(i,j,k)`, [mol kg<sup>-1</sup>])
- $B_{aoa}^{C}$ is the concentration of ammonia oxidizing archaea carbon biomass (`f_aoa(i,j,k)`, [mol kg<sup>-1</sup>])

---


### 13. Zooplankton grazing, egestion, excretion and assimilation.

**Grazing** by micro-zooplankton (`g_zoo`, $g_{mz}$, [s<sup>-1</sup>]) and meso-zooplankton (`g_mes`, $g_{Mz}$, [s<sup>-1</sup>]) is computed using a Holling Type III functional response [Holling, 1959](https://doi.org/10.4039/Ent91385-7), where:

$g_{mz} = \dfrac{\mu_{mz}^{max} (β_{hete})^{T} \sum_{i} \left(\varepsilon_{mz}^{i} \left(\phi_{mz}^{i} B_{i}^{C}\right)^{2}\right)}{\mu_{mz}^{max} (β_{hete})^{T} + \sum_{i} \left( \varepsilon_{mz}^{i} \phi_{mz}^{i} (B_{i}^{C})^{2}\right)}$\
$g_{Mz} = \dfrac{\mu_{Mz}^{max} (β_{hete})^{T} \sum_{i} \left(\varepsilon_{Mz}^{i} \left(\phi_{Mz}^{i} B_{i}^{C}\right)^{2}\right)}{\mu_{Mz}^{max} (β_{hete})^{T} + \sum_{i} \left( \varepsilon_{Mz}^{i} \left(\phi_{Mz}^{i} B_{i}^{C}\right)^{2}\right)}$

where 
- $\mu_{mz}^{max}$ is the maximum rate of micro-zooplankton grazing at 0ºC (`zoogmax`, [s<sup>-1</sup>])
- $\mu_{Mz}^{max}$ is the maximum rate of meso-zooplankton grazing at 0ºC (`mesgmax`, [s<sup>-1</sup>])
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless])
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC])
- $B_{i}^{C}$ is the concentration of prey type $i$ carbon biomass ([mmol C m<sup>-3</sup>])
- $\phi_{mz}^{i}$ is the relative prey preference of micro-zooplankton for prey type $i$ ([dimenionless])
- $\phi_{Mz}^{i}$ is the relative prey preference of meso-zooplankton for prey type $i$ ([dimenionless])
- $\varepsilon_{mz}^{i}$ is the prey capture rate coefficient of micro-zooplankton for prey type $i$ ([(mmol C m<sup>-3</sup>)<sup>-2</sup>])
- $\varepsilon_{Mz}^{i}$ is the prey capture rate coefficient of meso-zooplankton for prey type $i$ ([(mmol C m<sup>-3</sup>)<sup>-2</sup>])

and where

$\sum_{i} \phi_{mz}^{i} = \sum_{i} \phi_{Mz}^{i} = 1.0$

This formulation suppresses grazing at low prey biomass ($B_{i}^{C}$) due to reduced encounter and clearance rates, accelerates grazing at intermediate prey biomass as zooplankton effectively learn and switch to available prey, and saturates at high prey biomass due to handling-time limitation ([Gentleman and Neuheimer, 2008](https://doi.org/10.1093/plankt/fbn078); Rohr et al., [2022](https://doi.org/10.1016/j.pocean.2022.102878), [2024](https://doi.org/10.1029/2023GL107732)). This choice increases ecosystem stability and prolongs phytoplankton blooms relative to a Type II formulation.

The application of the temperature-dependent maximum growth rate in both the numerator and denominator makes this grazing formula unique [(Rohr et al., 2023)](https://www.nature.com/articles/s43247-023-00871-w) and equivalent to a disk formulation, rather than a Michaelis–Menten formulation [(Rohr et al., 2022)](https://doi.org/10.1016/j.pocean.2022.102878). Practically, this amplifies grazing in warmer climes, but to a lesser extent than other formulations that apply the temperature amplification ($(β_{hete})^{T}$) only in the numerator [(Rohr et al., 2023)](https://www.nature.com/articles/s43247-023-00871-w). This dampens the effect that variations in temperature have on grazing activity, amplifying the effect of $\varepsilon^{i}$ and aligning with observations that the ratio of grazing to phytoplankton growth varies little between tropical and polar climes [(Calbet and Landry, 2004)](https://doi.org/10.4319/lo.2004.49.1.0051). Theoretically, this assumes some evolutionary adaptation to account for the physiological effects of temperature across environmental niches, such that the efficiency of prey capture and handling becomes more important to grazers than metabolic constraints due to temperature.

The normalized prey preferences (i.e., dietary fractions) are further modified by prey switching prior to computation of total prey biomass ([Gentleman et al., 2003](https://doi.org/10.1016/j.dsr2.2003.07.001)) such that

$\phi_{z}^{i} = \left( \phi_{z}^{i} B_{i}^{C} \right)^{s_{z}}$ \

where
- $\phi_{z}^{i}$ is the relative prey preference of zooplankton type $z$ for prey type $i$
- $B_{i}^{C}$ is the concentration of prey type $i$ in carbon biomass
- $s_{z}$ is the prey-switching exponent of zooplankton type $z$ (`zoopreyswitch`; `mespreyswitch`)

When $s_{z} < 1$, zooplankton feed equally across all prey items irrespective of availability \ 
When $s_{z} = 1$, zooplankton feed according to pre-defined dietary fractions \
When $s_{z} > 1$, zooplankton exhibit prey-switching and feed disproportionately on most abundant prey \

Again, prey preferences are normalized to ensure 

$\sum_{i} \phi_{mz}^{i} = \sum_{i} \phi_{Mz}^{i} = 1.0$

The community average prey capture rate coefficients of micro-zooplankton (`zooeps(i,j,k)`, $\varepsilon_{mz}$, [(mmol C m<sup>-3</sup>)<sup>-2</sup>]) and meso-zooplankton (`meseps(i,j,k)`, $\varepsilon_{Mz}$, [(mmol C m<sup>-3</sup>)<sup>-2</sup>]) vary as a function of the prey biomasses and the consequential variations in prey preferences associated with prey-switching, which is consistent with the prey-dependent behaviour described by [Rohr et al. (2024)](doi.org/10.1029/2023GL107732). 

Total grazing of biomass by micro-zooplankton ([mol C kg<sup>-1</sup> day<sup>-1</sup>]) is therefore

$g_{mz}^{&larr; C} = g_{mz} B_{mz}^{C}$\
$g_{Mz}^{&larr; C} = g_{Mz} B_{Mz}^{C}$

where
- $g_{mz}$ is the total specific rate of grazing of micro-zooplankton (`g_zoo`, [s<sup>-1</sup])
- $g_{Mz}$ is the total specific rate of grazing of meso-zooplankton (`g_mes`, [s<sup>-1</sup])
- $B_{mz}^{C}$ is the in situ concentration of micro-zooplankton carbon biomass (`f_zoo(i,j,k)`, [mol C kg<sup>-1</sup>])
- $B_{Mz}^{C}$ is the in situ concentration of meso-zooplankton carbon biomass (`f_mes(i,j,k)`, [mol C kg<sup>-1</sup>])

Total grazing of prey can also be expressed as the sum of individual prey type consumption:

$g_{mz}^{&larr; C} = g_{mz}^{&larr; B_{np}^{C}} + g_{mz}^{&larr; B_{mp}^{C}} + g_{mz}^{&larr; B_{sd}^{C}} + g_{mz}^{&larr; B_{b1}^{C}} + g_{mz}^{&larr; B_{b2}^{C}} + g_{mz}^{&larr; B_{aoa}^{C}}$\
$g_{Mz}^{&larr; C} = g_{Mz}^{&larr; B_{np}^{C}} + g_{Mz}^{&larr; B_{mp}^{C}} + g_{Mz}^{&larr; B_{sd}^{C}} + g_{Mz}^{&larr; B_{ld}^{C}} + g_{Mz}^{&larr; B_{b1}^{C}} + g_{Mz}^{&larr; B_{b2}^{C}} + g_{Mz}^{&larr; B_{aoa}^{C}} + g_{Mz}^{&larr; B_{mz}^{C}}$

In this formulation, consumption of each prey item $i$ in [mol C kg<sup>-1</sup>] is equal to:

$g_{z}^{&larr; B_{i}^{C}} = g_{z} B_{z}^{C} \cdot \dfrac{\varepsilon_{mz}^{i} \left(\phi_{mz}^{i} B_{i}^{C} \right)^{2}}{\sum_{i} \varepsilon_{mz}^{i} \left(\phi_{mz}^{i} B_{i}^{C} \right)^{2}}$

Thus:
- $g_{mz}^{&larr; B_{np}^{C}}$ is the grazing rate of nano-phytoplankton by micro-zooplankton (`zoograzphy(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{mz}^{&larr; B_{mp}^{C}}$ is the grazing rate of micro-phytoplankton by micro-zooplankton (`zoograzdia(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{mz}^{&larr; B_{sd}^{C}}$ is the grazing rate of small particulate detritus by micro-zooplankton (`zoograzdet(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{mz}^{&larr; B_{b1}^{C}}$ is the grazing rate of facultative $NO_{3}$-reducing bacteria by micro-zooplankton (`zoograzbac1(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{mz}^{&larr; B_{b2}^{C}}$ is the grazing rate of facultative $N_{2}O$-reducing bacteria by micro-zooplankton (`zoograzbac2(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{mz}^{&larr; B_{aoa}^{C}}$ is the grazing rate of ammonia oxidizing archaea by micro-zooplankton (`zoograzaoa(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{Mz}^{&larr; B_{np}^{C}}$ is the grazing rate of nano-phytoplankton by meso-zooplankton (`mesgrazphy(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{Mz}^{&larr; B_{mp}^{C}}$ is the grazing rate of micro-phytoplankton by meso-zooplankton (`mesgrazdia(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{Mz}^{&larr; B_{sd}^{C}}$ is the grazing rate of small particulate detritus by meso-zooplankton (`mesgrazdet(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{Mz}^{&larr; B_{ld}^{C}}$ is the grazing rate of large particulate detritus by meso-zooplankton (`mesgrazbdet(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{Mz}^{&larr; B_{b1}^{C}}$ is the grazing rate of facultative $NO_{3}$-reducing bacteria by meso-zooplankton (`mesgrazbac1(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{Mz}^{&larr; B_{b2}^{C}}$ is the grazing rate of facultative $N_{2}O$-reducing bacteria by meso-zooplankton (`mesgrazbac2(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{Mz}^{&larr; B_{aoa}^{C}}$ is the grazing rate of ammonia oxidizing archaea by meso-zooplankton (`mesgrazaoa(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{Mz}^{&larr; B_{mz}^{C}}$ is the grazing rate of micro-zooplankton by meso-zooplankton (`mesgrazzoo(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])


**Zooplankton egestion, excretion and assimilation** are then calculated assuming static assimilation coefficients. Grazed biomass is routed to either egestion or ingestion via an ingestion coefficient ($\lambda^{C}$, [mol C (mol C)<sup>-1</sup>]), with the egested fraction being equal to $1.0 - \lambda^{C}$. The biomass that is ingested is then split between assimilation and excretion based on an assimilation coefficient ($\eta^{C}$, [mol C (mol C)<sup>-1</sup>]) with the excreted fraction being equal to $1.0 - \eta^{C}$. Egestion ($E$), excretion ($X$) and assimilation ($A$) of organic carbon due to grazing of prey type $i$ by zooplankton type $z$ are:

$E_{z}^{&larr; B_{i}^{C}} = g_{z}^{&larr; B_{i}^{C}} \left(1 - \lambda_{z}^{C} \right)$\
$X_{z}^{&larr; B_{i}^{C}} = g_{z}^{&larr; B_{i}^{C}} \lambda_{z}^{C} \left(1 - \eta_{z}^{C} \right)$\
$A_{z}^{&larr; B_{i}^{C}} = g_{z}^{&larr; B_{i}^{C}} \lambda_{z}^{C} \eta_{z}^{C}$

where 
- $E_{z}^{&larr; B_{i}^{C}}$ is the rate of egestion of carbon biomass by zooplankton type $z$ feeding on prey type $i$ ([mol C kg<sup>-1</sup>])
- $X_{z}^{&larr; B_{i}^{C}}$ is the rate of excretion of carbon biomass by zooplankton type $z$ feeding on prey type $i$ ([mol C kg<sup>-1</sup>])
- $A_{z}^{&larr; B_{i}^{C}}$ is the rate of assimilation of carbon biomass by zooplankton type $z$ feeding on prey type $i$ ([mol C kg<sup>-1</sup>])
- $g_{z}^{&larr; B_{i}^{C}}$ is the grazing rate of zooplankton type $z$ on prey type $i$ ([mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $\lambda_{z}^{C}$ is the fraction of prey carbon biomass that is ingested by zooplankton type $z$ (`zooCingest`; `mesCingest`, [mol C (mol C)<sup>-1</sup>])
- $\eta_{z}^{C}$ is the fraction of ingested prey carbon biomass that is assimilated by zooplankton type $z$ (`zooCassim`; `mesCassim`, [mol C (mol C)<sup>-1</sup>])

Total egestion, excretion and assimilation or carbon are therefore:

$E_{z}^{&larr; C} = g_{z}^{&larr; C} \left(1 - \lambda_{z}^{C} \right)$\
$X_{z}^{&larr; C} = g_{z}^{&larr; C} \lambda_{z}^{C} \left(1 - \eta_{z}^{C} \right)$\
$A_{z}^{&larr; C} = g_{z}^{&larr; C} \lambda_{z}^{C} \eta_{z}^{C}$


Because we track both carbon and iron through the ecosystem components, we assign unique ingestion and assimilation coefficients to carbon and iron. This separation of ingestion and assimilation coefficients for iron and carbon follows [Le Mézo & Galbraith (2021)](https://doi.org/10.1002/lno.11597). For iron, we apply unique ingestion ($\lambda^{Fe}$, [mol Fe (mol Fe)<sup>-1</sup>]) and assimilation coefficients ($\eta^{Fe}$, [mol Fe (mol Fe)<sup>-1</sup>]). [Le Mézo & Galbraith (2021)](https://doi.org/10.1002/lno.11597) show that if $\lambda^{Fe} << \lambda^{C}$ then egestion is enriched in Fe:C, and it follows that $\eta^{Fe} >> \eta^{C}$ so that zooplankton can absorb sufficient iron from their prey. Consequently:

$E_{z}^{&larr; B_{i}^{Fe}} = g_{z}^{&larr; B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \left(1 - \lambda_{z}^{Fe} \right)$\
$X_{z}^{&larr; B_{i}^{Fe}} = g_{z}^{&larr; B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \lambda_{z}^{Fe} \left(1 - \eta_{z}^{Fe} \right)$\
$A_{z}^{&larr; B_{i}^{Fe}} = g_{z}^{&larr; B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \lambda_{z}^{Fe} \eta_{z}^{Fe}$

where 
- $E_{z}^{&larr; B_{i}^{Fe}}$ is the rate of egestion of iron biomass by zooplankton type $z$ feeding on prey type $i$ ([mol Fe kg<sup>-1</sup>])
- $X_{z}^{&larr; B_{i}^{Fe}}$ is the rate of excretion of iron biomass by zooplankton type $z$ feeding on prey type $i$ ([mol Fe kg<sup>-1</sup>])
- $A_{z}^{&larr; B_{i}^{Fe}}$ is the rate of assimilation of iron biomass by zooplankton type $z$ feeding on prey type $i$ ([mol Fe kg<sup>-1</sup>])
- $g_{z}^{&larr; B_{i}^{C}}$ is the grazing rate of zooplankton type $z$ on prey type $i$ ([mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $\dfrac{B_{i}^{Fe}}{B_{i}^{C}}$ is the Fe:C ratio of prey type $i$ ([mol Fe (mol C)<sup>-1</sup>])
- $\lambda_{z}^{Fe}$ is the fraction of prey iron biomass that is ingested by zooplankton type $z$ (`zooFeingest`; `mesFeingest`, [mol Fe (mol Fe)<sup>-1</sup>])
- $\eta_{z}^{Fe}$ is the fraction of ingested prey iron biomass that is assimilated by zooplankton type $z$ (`zooFeassim`; `mesFeassim`, [mol Fe (mol Fe)<sup>-1</sup>])

Total egestion, excretion and assimilation or iron are therefore:

$E_{z}^{&larr; Fe} = \sum_{i} \left( g_{z}^{&larr; B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \right) \cdot \left(1 - \lambda_{z}^{Fe} \right)$\
$X_{z}^{&larr; Fe} = \sum_{i} \left( g_{z}^{&larr; B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \right) \cdot \lambda_{z}^{Fe} \left(1 - \eta_{z}^{Fe} \right)$\
$A_{z}^{&larr; Fe} = \sum_{i} \left( g_{z}^{&larr; B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \right) \cdot \lambda_{z}^{Fe} \eta_{z}^{Fe}$


**Excretion of nitrogen**\
For zooplankton preying on phytoplankton, other zooplankton and detritus, their excretion of nitrogen will be equal to:

$X_{z}^{&larr; B_{i}^{N}} = X_{z}^{&larr; B_{i}^{C}} \dfrac{16}{122}$

However, since both micro-zooplankton and meso-zooplankton consume heterotrophic bacteria and ammonia oxidizing arcaheal types, which have different C:N ratios to other ecosystem biomass components, we must also compute the specific excretion of $NH_4$ and $B_{DOM}^{N}$ by zooplankton when feeding on these types. Since these types are richer in N than the other prey types, zooplankton excrete more $NH_4$ and $B_{DOM}^{N}$ when bacteria and archaea represent a greater proportion of their diet ([Sterner & Elser, 2002](https://press.princeton.edu/books/ebook/9781400885695/ecological-stoichiometry-pdf)). Total excretion of nitrogen from bacterial/archaeal type $i$ by zooplankton type $z$ is as follows:

$X_{z}^{&larr; B_{i}^{N}} = g_{z}^{&larr; B_{i}^{C}} \dfrac{1}{R_{i}^{C:N}} - \dfrac{A_{z}^{&larr; B_{i}^{C}} + E_{z}^{&larr; B_{i}^{C}}}{R_{z}^{C:N}}$


---


### 14. Calcium carbonate production and dissolution.

**Dynamic $CaCO_3$ production and dissolution**

When $CaCO_3$ dynamics are enabled (`do_caco3_dynamics = .true.`), the model computes both particulate inorganic carbon production (via the PIC:POC ratio) and $CaCO_3$ dissolution rates as functions of carbonate chemistry, temperature, and organic matter availability.

**Production** of $CaCO_3$ in WOMBAT-mid comes from five sources: (1) density-dependent mortality of nano-phytoplankton (i.e., coccolithophorids), (2) density-dependent mortality of micro-zooplankton (i.e., foraminifera), (3) micro-zooplankton egestion of grazed nano-phytoplankton, (4) meso-zooplankton egestion of grazed nano-phytoplankton, and (5) meso-zooplankton egestion of grazed micro-zooplankton. Each term is multiplied by the particulate inorganic to organic carbon production ratio (`pic2poc`, $PIC:POC$, [mol/mol]) to return a rate of $CaCO_3$ production in mol C kg<sup>-1</sup> s<sup>-1</sup>.

(1) $P_{CaCO_3}^{\Gamma_{np}^{C}} = \Gamma_{np}^{&rarr; C} \cdot PIC:POC$\
(2) $P_{CaCO_3}^{\Gamma_{mz}^{C}} = \Gamma_{mp}^{&rarr; C} \cdot PIC:POC$\
(3) $P_{CaCO_3}^{g_{mz}^{&larr; B_{np}^{C}}} = g_{mz}^{&larr; B_{np}^{C}} \cdot PIC:POC \left(1 - F_{gut}\right)$\
(4) $P_{CaCO_3}^{g_{Mz}^{&larr; B_{np}^{C}}} = g_{Mz}^{&larr; B_{np}^{C}} \cdot PIC:POC \left(1 - F_{gut}\right)$\
(5) $P_{CaCO_3}^{g_{Mz}^{&larr; B_{mz}^{C}}} = g_{Mz}^{&larr; B_{mz}^{C}} \cdot PIC:POC \left(1 - F_{gut}\right)$

where
- $\Gamma_{np}^{&rarr; C}$ is the quadratic (density-dependent) loss rate of nano-phytoplankton biomass (`phymorq`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $\Gamma_{mz}^{&rarr; C}$ is the quadratic (density-dependent) loss rate of micro-zooplankton biomass (`zoomorq`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{mz}^{&larr; B_{np}^{C}}$ is the grazing rate of nano-phytoplankton by micro-zooplankton (`zoograzphy(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{Mz}^{&larr; B_{np}^{C}}$ is the grazing rate of nano-phytoplankton by meso-zooplankton (`mesgrazphy(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{Mz}^{&larr; B_{mz}^{C}}$ is the grazing rate of micro-zooplankton by meso-zooplankton (`mesgrazzoo(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $F_{gut}$ is the fraction of $CaCO_3$ that is dissolved within zooplankton guts (`fgutdiss`, [mol C (mol C)<sup>-1</sup>])

In the above, the $PIC:POC$ ratio is formulated as 
 
$PIC:POC = \min \left( 0.3,  \left( f_{\text{inorg}} + 10^{-3 + 4.31 \times 10^{-6} \left( \dfrac{[HCO_3^-]}{[H^+]} \right)} \right) F_T \right)$

where
- $f_{inorg}$ is the background PIC:POC ratio (f_inorg, [mol C (mol C)-1])
- $[HCO_{3}^{-}]$ is the concentration of bicarbonate ions (hco3, [mol kg-1])
- $[H^{+}]$ is concentration of free hydrogen ions (htotal(i,j,k), [µmol kg-1])
- $F_{T}$ is a temperature-dependent suppression term and if defined by $F_{T} = 0.55 + 0.45 \cdot \tanh\left(T - 4 \right)$ 

This formulation of $PIC:POC$ is therefore a function of the substrate–inhibitor ratio between bicarbonate and free hydrogen ions (`hco3 / htotal(i,j,k)`, $\dfrac{[HCO_{3}^{-}]}{[H^{+}]}$, [mol µmol<sup>-1</sup>]), following [Lehmann & Bach (2025)](https://www.nature.com/articles/s41561-025-01644-0). This reflects the sensitivity of calcification to carbonate system speciation, which is nonlinearly enhanced with an increasing substrate-inhibitor ratio. Moreover, the $F_{T}$ term strongly reduces production in cold waters, enforcing near-zero calcification below approximately 3°C consistent with observations of _Emiliania huxleyi_ growth limits in polar environments ([Fielding, 2013](https://doi.org/10.4319/lo.2013.58.2.0663)). Finally, we also cap the $PIC:POC$ ratio at an upper bound of 0.3 to prevent unrealistically high inorganic carbon production and accord with the highest measured ratios in the ocean.

**Dissolution** of $CaCO_3$ is computed as the sum of five contributions: 

(1) undersaturation-driven dissolution of calcite (`caldiss(i,j,k)`, $D_{CaCO_3}^{\Omega_{cal}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
(2) undersaturation-driven dissolution of aragonite (`aradiss(i,j,k)`, $D_{CaCO_3}^{\Omega_{ara}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
(3) biologically-mediated dissolution associated with degredation of small detrital organic matter (`pocdiss(i,j,k)`, $D_{CaCO_3}^{\Gamma_{sd}^{&rarr; C}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
(4) dissolution within micro-zooplankton during their digestion of detrital aggregates (`zoodiss(i,j,k)`, $D_{CaCO_3}^{g_{mz}^{&larr; B_{sd}^{C}}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
(5) dissolution within meso-zooplankton during their digestion of detrital aggregates (`mesdiss(i,j,k)`, $D_{CaCO_3}^{g_{Mz}^{&larr; B_{sd}^{C}}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>])

Total $CaCO_3$ dissolution is:

$D_{CaCO_3} = D_{CaCO_3}^{\Omega_{cal}} + D_{CaCO_3}^{\Omega_{ara}} + D_{CaCO_3}^{\Gamma_{sd}^{&rarr; C}} + D_{CaCO_3}^{g_{mz}^{&larr; B_{sd}^{C}}} + D_{CaCO_3}^{g_{Mz}^{&larr; B_{sd}^{C}}}$

This formulation, at leas the first three terms, follows [Kwon et al. (2024)](https://www.science.org/doi/full/10.1126/sciadv.adl0779).

$D_{CaCO_3}^{\Omega_{cal}} = d_{CaCO_3}^{\Omega_{cal}} \max\left(0,  1 - \Omega_{cal}\right)^{2.2} B_{CaCO_3}^{C}$\ 
$D_{CaCO_3}^{\Omega_{ara}} = d_{CaCO_3}^{\Omega_{ara}} \max\left(0,  1 - \Omega_{ara}\right)^{1.5} B_{CaCO_3}^{C}$\
$D_{CaCO_3}^{\Gamma_{sd}^{&rarr; C}} = d_{CaCO_3}^{\Gamma_{sd}} \Gamma_{sd}^{&rarr; C} B_{CaCO_3}^{C}$

where
- $\Omega_{cal}$ is the saturation state of calcite (`omega_cal(i,j,k)`, [dimenionless])
- $\Omega_{ara}$ is the saturation state of calcite (`omega_ara(i,j,k)`, [dimenionless])
- $d_{CaCO_3}^{\Omega_{cal}}$ is the reference dissolution rate constant for calcite (`disscal`, [s<sup>-1</sup>]) 
- $d_{CaCO_3}^{\Omega_{ara}}$ is the reference dissolution rate constant for aragonite (`dissara`, [s<sup>-1</sup>]) 
- $d_{CaCO_3}^{\Gamma_{sd}}$ is the reference dissolution rate constant per unit of small detrital organic carbon remineralised (`dissdet`, [(mmol C m<sup>-3</sup>)<sup>-1</sup>]) 
- $\Gamma_{sd}^{&rarr; C}$ is the in situ remineralisation rate of small detrital organic carbon (`detremi(i,j,k)`, [mmol C m<sup>-3</sup> s<sup>-1</sup>])
- $B_{CaCO_3}^{C}$ is the in situ concentration of $CaCO_3$ in carbon units (`f_caco3(i,j,k)`, [mol C kg<sup>-1</sup>])

For $D_{CaCO_3}^{\Omega_{cal}}$ and $D_{CaCO_3}^{\Omega_{ara}}$, dissolution is activated only under undersaturated conditions ($\Omega_{cal} < 1$; $\Omega_{ara} < 1$) and increases nonlinearly with increasing undersaturation. In contrast, $D_{CaCO_3}^{\Gamma_{sd}^{&rarr; C}}$ represents shallow water dissolution due to reducing microenvironments. In this scenario, $\Omega_{cal}$ and $\Omega_{ara}$ tend to be > 1 ([Sulpis et al., 2021](https://doi.org/10.1038/s41561-021-00743-y)) but dissolution nonetheless occurs in microenvironments enriched in $CO_{2}^{*}$ due to heterotrophic activity ([Borer et al., 2026](https://doi.org/10.1073/pnas.2510025123)).

The fourth and fifth terms (`zoodiss(i,j,k)`; `mesdiss(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) represent dissolution of $CaCO_3$ during zooplankton digestion of detrital particulates. 

$D_{CaCO_3}^{g_{mz}^{&larr; B_{sd}^{&larr; C}}} = g_{mz}^{&larr; B_{sd}^{C}} F_{gut} \dfrac{B_{CaCO_3}^{C}}{B_{sd}^{C}}$\
$D_{CaCO_3}^{g_{Mz}^{&larr; B_{sd}^{&larr; C}}} = g_{Mz}^{&larr; B_{sd}^{C}} F_{gut} \dfrac{B_{CaCO_3}^{C}}{B_{sd}^{C}}$

where
- $g_{mz}^{&larr; B_{sd}^{C}}$ is the grazing rate of small particulate detritus by micro-zooplankton (`zoograzdet(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $g_{Mz}^{&larr; B_{sd}^{C}}$ is the grazing rate of small particulate detritus by meso-zooplankton (`mesgrazdet(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $F_{gut}$ is the fraction of $CaCO_3$ that is dissolved within zooplankton guts (`fgutdiss`, [mol C (mol C)<sup>-1</sup>])
- $\dfrac{B_{CaCO_3}^{C}}{B_{sd}^{C}}$ is the in situ ratio of $CaCO_3$ to small organic carbon detritus (`biocaco3/biodet`, [mol C (mol C)<sup>-1</sup>])

Here we note that the processing of $CaCO_3$ by zooplankton grazing is treated differently to processing of organic carbon. For organic carbon, we route the biomass between zooplankton biomass (assimilation), inorganic nutrients (excretion) and particulate detritus (egestion). For $CaCO_3$ consumption by both micro-zooplankton and meso-zooplankton the $CaCO_3$ is not assimilated since it does not contain nitrogen or other key elements for biosynthesis, and so is only routed between excretion to DIC and alkalinity or goes undissolved and remains $CaCO_3$ that sinks through the water column. This is supported by the fact that micro- and meso-zooplankton may dissolve 92±7% and 38-73% of coccolithophore calcite during feeding, respectively ([Smith et al., 2024](https://doi.org/10.1126/sciadv.adr5453); [White et al., 2018](https://doi.org/10.1038/s41598-018-28073-x); [Harris, 1994](https://doi.org/10.1007/BF00347540)), and that the remainder is excreted and not assimilated ([Mayers et al., 2020](https://doi.org/10.3389/fmars.2020.569896)).

**Static $CaCO_{3}$ production and dissolution**

When $CaCO_3$ dynamics are disabled (`do_caco3_dynamics = .false.`), the model uses a static PIC:POC ratio (`f_inorg + 0.025`, [mol C (mol C)<sup>-1</sup>]) and $CaCO_3$ dissolution rate (`caco3lrem`, [s<sup>-1</sup>]). These are set as input parameters to the model.

---


### 15. Implicit nitrogen fixation.

Because we do not consider diazotrophs as an explicit phytoplankton functional type, we represent the fixation of nitrogen implicitly using a simple parameterization dependent on temperature, nutrient and light availability. The equation for new nitrogen (specifically $NH_4$) added via diazotrophy is:

$\mu_{diazo}^{&rarr; NH_4} = \mu_{diazo}^{max} \left(1 - L_{np}^{N} \right) \cdot \min\left(L_{diazo}^{Fe}, L_{diazo}^{PAR}\right) R_{diazo}^{N:C} \cdot 1 \times 10^{-6}$ 

where
- $\mu_{diazo}^{max}$ is the temperature-dependent maximum growth rate of diazotrophs (`trimumax(i,j,k)`, [s<sup>-1</sup>])
- $L_{np}^{N}$ is the limitation term of nano-phytoplankton growth on nitrogen (`phy_lnit(i,j,k)`, [dimensionless]) 
- $L_{diazo}^{Fe}$ is the limitation term of diazotroph growth on iron (`tri_lfer(i,j,k)`, [dimensionless]) 
- $L_{diazo}^{PAR}$ is the limitation term of diazotroph growth on light (`tri_lpar(i,j,k)`, [dimensionless]) 
- $R_{diazo}^{N:C}$ is the ratio of N:C within diazotrophic biomass (`trin2c`, [mol N (mol C)<sup>-1</sup>])
- $1 \times 10^{-6}$ is a conversion factor to mol kg<sup>-1</sup>. 

The temperature-dependent maximum growth rate ($\mu_{diazo}^{max}$) is taken directly from [Wrightson et al. (2022)]( https://doi.org/10.1111/gcb.16399) who based their formulation on the work of [Jiang et al. (2018)](https://doi.org/10.1038/s41558-018-021):

$\mu_{diazo}^{max} = \left( -0.000399(T)^{3} + 0.02685(T)^{2} - 0.555T + 3.633 \right) \dfrac{1}{86400}$ 

where
- $T$ is in situ water temperature (`Temp(i,j,k)`, [ºC]) and we only consider $T > 15.8$ºC
- $\dfrac{1}{86400}$ converts their formula from units of [day<sup>-1</sup>] to [s<sup>-1</sup>]

The iron and light limitation terms are as follows:

$L_{diazo}^{Fe} = \dfrac{dFe}{dFe + K_{diazo}^{Fe}}$\
$L_{diazo}^{PAR} = 1 - e^{- \alpha_{diazo} PAR }$

where
- $dFe$ is the in situ concentration of dissolved iron (`biofer`, [nmol Fe kg<sup>-1</sup>])
- $K_{diazo}^{Fe}$ is the half-saturation coefficient for uptake of dissolved iron by diazotrophs (`trikf`, [nmol Fe kg<sup>-1</sup>])
- $\alpha_{diazo}$ is the chlorophyll-adjusted slope of the photosynthesis-irradience curve of diazotrophs (`alphabio_tri * trichlc`, [(W m<sup>-2</sup>)<sup>-1</sup>])
- $PAR$ is the downwelling photosynthetically available radiation (`radbio`, [W m<sup>-2</sup>])

---


### 16. Facultative bacterial heterotrophy.

We remineralise dissolved organic matter into inorganic consituents via the activity of two facultative bacterial heterotrophs. These bacterial heterotrophs oxidise dissolved organic carbon ($B_{DOM}^{C}$) and reduce dissolved oxygen ($O_2$). However, we consider these bacterial types, which relfect the traits of the ubiquitous SAR11, to be faculatively anaerobic ([Zumft, 1997](https://doi.org/10.1128/mmbr.61.4.533-616.1997); [Tsementzi et al., 2016](https://doi.org/10.1038/nature19068)). This means that they can shift their metabolism to using either nitrate ($NO_3$) or nitrous oxide ($N_{2}O$) as alternative electron acceptors when $O_2$ is limiting. These populations of heterotrophic bacteria also assimilate dissolved organic nitrogen ($DON$) and ammonium ($NH_4$) and dissolved iron ($dFe$) to support biosynthesis. By taking up $NH_4$ and $dFe$ bacteria compete directly with phytoplankton, consistent with prior observations ([Kirchman, 1994](https://www.jstor.org/stable/4251383); [Tortell et al., 1996](https://doi.org/10.1038/383330a0); [Kirchman & Wheeler, 1998](https://doi.org/10.1016/S0967-0637(97)00075-7); [Fourquez et al., 2015](https://doi.org/10.5194/bg-12-1893-2015); [Deng et al., 2021](https://doi.org/10.1002/lno.11883); [Strzepek et al., 2025](https://doi.org/10.1093/ismejo/wraf015)).

Our formulation of heterotrophic bacterial growth follows that developed by [Zakem et al. (2020)](https://doi.org/10.1038/s41396-019-0523-8) and subsequently expanded in [Sun et al. (2024)](https://doi.org/10.1073/pnas.2417421121) and [Buchanan et al. (2025)](https://www.science.org/doi/full/10.1126/science.ado0742). In these studies, the realized biomass growth rate (integration of carbon into biomass) of bacterial functional type $b$ (`bac1grow(i,j,k)`; `bac2grow(i,j,k)`, $\mu_{b}^{&larr; C}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) is defined by:

$\mu_{b}^{&larr; C} = \max\left(\mu_{b}^{aer}, \mu_{b}^{ana} \right) (β_{hete})^{T} B_{b}^{C}$

where
- $\mu_{b}^{aer}$ is the realized growth rate due to aerobic metabolism (`bac_muaer`, [s<sup>-1</sup>])
- $\mu_{b}^{ana}$ is the realized growth rate due to anaerobic metabolism (`bac_muana`, [s<sup>-1</sup>])
- $(β_{hete})^{T}$ is the temperature-dependent scaling on heterotrophic metabolism (`fbc`, [dimensionless])
- $B_{b}^{C}$ is the in situ concentration of bacterial functional type $b$ (`f_bac1(i,j,k)`; `f_bac2(i,j,k)`, [mol C kg<sup>-1</sup>])

Thus, whichever of aerobic and anaerobic metabolism offers the greatest growth rate will be chosen as the means by which bacteria grow. We must define what controls aerobic ($\mu_{b}^{aer}$) and anaerobic ($\mu_{b}^{ana}$) growth. In the case of both aerobic and anaerobic metabolisms we compute these growth rates as the minimum of three rates associated with four essential resources: dissolved organic carbon ($B_{DOM}^{C}$), nitrogen ($N$), dissolved iron ($dFe$) and the electron acceptor ($EA$).

$\mu_{b}^{aer} = \min \left(\mu_{b}^{aer(DOC)}, \mu_{b}^{aer(N)}, \mu_{b}^{aer(dFe)}, \mu_{b}^{aer(EA)} \right)$\
$\mu_{b}^{ana} = \min \left(\mu_{b}^{ana(DOC)}, \mu_{b}^{ana(N)}, \mu_{b}^{ana(dFe)}, \mu_{b}^{ana(EA)} \right)$

For aerobic growth, these resource-specific growth rates are calculated as:

$\mu_{b}^{aer(DOC)} = V_{b}^{DOC} y_{b}^{aer(DOC)}$\
$\mu_{b}^{aer(N)} = \left(V_{b}^{DON} + V_{b}^{NH_{4}}\right) y_{b}^{aer(N)}$\
$\mu_{b}^{aer(dFe)} = V_{b}^{dFe} y_{b}^{aer(Fe)}$\
$\mu_{b}^{aer(EA)} = V_{b}^{O_{2}} y_{b}^{O_{2}}$

where
- $V_{b}^{DOC}$ is the maximum uptake rate of dissolved organic carbon by bacterial type $b$ (`bac_Vdoc`, [mmol C m<sup>-3</sup> s<sup>-1</sup>])
- $V_{b}^{DON}$ is the maximum uptake rate of dissolved organic nitrogen by bacterial type $b$ (`bac1_Vdon`; `bac2_Vdon`, [mmol N m<sup>-3</sup> s<sup>-1</sup>])
- $V_{b}^{NH_4}$ is the maximum uptake rate of ammonium by bacterial type $b$ (`bac1_Vnh4`; `bac2_Vnh4`, [mmol N m<sup>-3</sup> s<sup>-1</sup>])
- $V_{b}^{dFe}$ is the maximum uptake rate of dissolved iron by bacterial type $b$ (`bac_VdFe`, [mmol Fe m<sup>-3</sup> s<sup>-1</sup>])
- $V_{b}^{O_2}$ is the maximum uptake rate of dissolved oxygen by bacterial type $b$ (`bac_Voxy`, [mmol C m<sup>-3</sup> s<sup>-1</sup>])
- $y_{b}^{aer(DOC)}$ is the aerobic biomass growth yield on dissolved organic carbon by bacterial type $b$ (`bac1_ydoc(i,j,k)`; `bac2_ydoc(i,j,k)`, [mol C (mol C)<sup>-1</sup>])
- $y_{b}^{aer(N)}$ is the aerobic biomass growth yield on nitrogen by bacterial type $b$ (`bac1_ydonC`; `bac2_ydonC`, [mol C (mol N)<sup>-1</sup>])
- $y_{b}^{aer(Fe)}$ is the aerobic biomass growth yield on iron by bacterial type $b$ (`bac1_C2Fe`; `bac2_C2Fe`, [mol C (mol Fe)<sup>-1</sup>])
- $y_{b}^{O_2}$ is the biomass growth yield on dissolved oxygen by bacterial type $b$ (`bac1_yoxyC`; `bac2_yoxyC`, [mol C (mol $O_{2}$)<sup>-1</sup>])

For anaerobic growth, these resource-specific growth rates are calculated in the same manner, except that the yields associated with growth are altered to reflect anaerobic metabolism. Also, the electron acceptor is no longer oxygen and is now either $NO_3$ for the $NO_3$-reducing bacteria (`f_bac1(i,j,k)`, $b1$) or $N_{2}O$ for the $N_{2}O$-reducing bacteria (`f_bac2(i,j,k)`, $b2$):

$\mu_{b}^{ana(DOC)} = V_{b}^{DOC} y_{b}^{ana(DOC)}$\
$\mu_{b}^{ana(N)} = \left(V_{b}^{DON} + V_{b}^{NH_{4}}\right) y_{b}^{ana(N)}$\
$\mu_{b}^{ana(dFe)} = V_{b}^{dFe} y_{b}^{ana(Fe)}$\
$\mu_{b1}^{NO_{3}} = V_{b}^{NO_{3}} y_{b}^{NO_{3}}$
$\mu_{b2}^{N_{2}O} = V_{b}^{N_{2}O} y_{b}^{N_{2}O}$

where
- $V_{b}^{NO_3}$ is the maximum uptake rate of nitrate by bacterial type $b1$ (`bac_Vno3`, [mmol N m<sup>-3</sup> s<sup>-1</sup>])
- $V_{b}^{N_{2}O}$ is the maximum uptake rate of nitrous oxide by bacterial type $b2$ (`bac_Vn2o`, [mmol N m<sup>-3</sup> s<sup>-1</sup>])
- $y_{b}^{ana(DOC)}$ is the anaerobic biomass growth yield on dissolved organic carbon by bacterial type $b$ (`bac1_yanaC`; `bac2_yanaC`, [mol C (mol C)<sup>-1</sup>])
- $y_{b}^{ana(N)}$ is the anaerobic biomass growth yield on nitrogen by bacterial type $b$ (`bac1_ydonC * bacanapen`; `bac2_ydonC * bacanapen`, [mol C (mol N)<sup>-1</sup>])
- $y_{b}^{ana(Fe)}$ is the anaerobic biomass growth yield on iron by bacterial type $b$ (`bac1_C2Fe * bacanapen`; `bac2_C2Fe * bacanapen`, [mol C (mol Fe)<sup>-1</sup>])
- $y_{b}^{NO_3}$ is the biomass growth yield on nitrate by bacterial type $b1$ (`bac1_yno3C`, [mol C (mol $NO_{3}$)<sup>-1</sup>])
- $y_{b}^{N_{2}O}$ is the biomass growth yield on nitrous oxide by bacterial type $b2$ (`bac2_yn2o`, [mol C (mol $N_{2}O$)<sup>-1</sup>])

Whether aerobic or anaerobic metabolism results in higher growth is therefore dependent on the differences in substrate **uptake rates ($V_{b}$)** and the substrate-specific **biomass yields ($y_{b}$)**. It is crucial to estimate both quantities. 

**Uptake rates**\
Uptake rates of reductant ($DOC$), addition resources ($DON$, $NH_4$ and $dFe$) and oxidant (electron acceptors) are calculated as:

$V_{b}^{DOC} = V_{b}^{max,DOC} \cdot \dfrac{B_{DOM}^{C}}{B_{DOM}^{C} + K_{b}^{DOC}}$\
$V_{b}^{DON} = V_{b}^{max,DON} \cdot \dfrac{B_{DOM}^{N}}{B_{DOM}^{N} + K_{b}^{DON}}$\
$V_{b}^{NH_4} = V_{b}^{max,NH_4} \cdot \dfrac{NH_4}{NH_4 + K_{b}^{NH_4}}$\
$V_{b}^{dFe} = V_{b}^{max,dFe} \cdot \dfrac{dFe}{dFe + K_{b}^{dFe}}$\
$V_{b}^{O_{2}} = \rho_{b}^{O_2} \cdot O_{2}$\
$V_{b}^{NO_{3}} = V_{b}^{max,NO_{3}} \cdot \dfrac{NO_{3}}{NO_{3} + K_{b}^{NO_{3}}}$\
$V_{b}^{N_{2}O} = \rho_{b}^{N_{2}O} \cdot N_{2}O$

where
- $V_{b}^{max,DOC}$ is the maximum rate of DOC uptake by bacterial functional type $b$ (`bac1_Vmax_doc`; `bac2_Vmax_doc`, [mmol C m<sup>-3</sup> s<sup>-1</sup>])
- $V_{b}^{max,DON}$ is the maximum rate of DON uptake by bacterial functional type $b$ (`bac1_Vmax_don`; `bac2_Vmax_don`, [mmol N m<sup>-3</sup> s<sup>-1</sup>])
- $V_{b}^{max,NH_{4}}$ is the maximum rate of $NH_{4}$ uptake by bacterial functional type $b$ (`bac1_Vmax_nh4`; `bac2_Vmax_nh4`, [mmol N m<sup>-3</sup> s<sup>-1</sup>])
- $V_{b}^{max,dFe}$ is the maximum rate of $dFe$ uptake by bacterial functional type $b$ (`bac1_Vmax_dFe`; `bac2_Vmax_dFe`, [mmol Fe m<sup>-3</sup> s<sup>-1</sup>])
- $V_{b}^{max,NO_{3}}$ is the maximum rate of $NO_{3}$ uptake by bacterial functional type $b$ (`bac1_Vmax_no3`, [mmol N m<sup>-3</sup> s<sup>-1</sup>])
- $\rho_{b}^{O_{2}}$ is the diffusive uptake limit of $O_2$ by bacterial functional type $b$ (`bac1_poxy`; `bac2_poxy`, [(mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>])
- $\rho_{b}^{N_{2}O}$ is the diffusive uptake limit of $N_{2}O$ by bacterial functional type $b$ (`bac2_pn2o`, [(mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>])
- $K_{b}^{DOC}$ is the half-saturation coefficient for uptake of $DOC$ by bacterial functional type $b$ (`bac1_kdoc`; `bac2_kdoc`, [mmol C m<sup>-3</sup>])
- $K_{b}^{DON}$ is the half-saturation coefficient for uptake of $DON$ by bacterial functional type $b$ (`bac1_kdon`; `bac2_kdon`, [mmol N m<sup>-3</sup>])
- $K_{b}^{NH_{4}}$ is the half-saturation coefficient for uptake of $NH_4$ by bacterial functional type $b$ (`bac1_knh4`; `bac2_knh4`, [mmol N m<sup>-3</sup>])
- $K_{b}^{dFe}$ is the half-saturation coefficient for uptake of $dFe$ by bacterial functional type $b$ (`bac1_kfer`; `bac2_kfer`, [µmol Fe m<sup>-3</sup>])
- $K_{b}^{NO_{3}}$ is the half-saturation coefficient for uptake of $NO_3$ by bacterial functional type $b$ (`bac1_kno3`, [mmol N m<sup>-3</sup>])
- $B_{DOM}^{C}$ is the in situ concentration of dissolved organic carbon (`biodoc`, [mmol C m<sup>-3</sup>])
- $B_{DOM}^{N}$ is the in situ concentration of dissolved organic nitrogen (`biodon`, [mmol N m<sup>-3</sup>])
- $NH_{4}$ is the in situ concentration of $NH_{4}$ (`bionh4`, [mmol N m<sup>-3</sup>])
- $dFe$ is the in situ concentration of $dFe$ (`biofer`, [µmol Fe m<sup>-3</sup>])
- $NO_{3}$ is the in situ concentration of $NO_{3}$$ (`biono3`, [mmol N m<sup>-3</sup>])
- $O_{2}$ is the in situ concentration of $O_{2}$$ (`biooxy`, [mmol $O_{2}$ m<sup>-3</sup>])
- $N_{2}O$ is the in situ concentration of $N_{2}O$ (`bion2o`, [mmol $N_{2}O$ m<sup>-3</sup>])


**Biomass yields**\
We further expand on bacterial heterotrophic dynamics by also integrating the findings of [Wang & Kuzyakov (2023)](https://doi.org/10.1111/gcb.16925) to solve for biomass yields as a function of the nominal oxidation state of carbon (NOSC) of DOM (`f_nosdoc(i,j,k)`, $DOM^{NOSC}$, [dimensionless]). [Wang & Kuzyakov (2023)](https://doi.org/10.1111/gcb.16925) conceptually link the carbon use efficiency (i.e., yield) of bacteria to the NOSC. They identify a theoretical positive relationship between yield and NOSC and suggest that more oxidized compounds offer a greater energy of content per carbon atom because more reduced compounds require greater processing costs. To account for this dynamic, we scale the biomass yield of heterotrophic bacteria (`bac_ydon(i,j,k)`, $y_{b}^{DON}$, [mol N biomass (mol DON)<sup>-1</sup>]) as:

$y_{b}^{DON} = \max\left( y_{b}^{min(DON)}, \min\left( y_{b}^{max(DON)}, y_{b}^{min(DON)} + DOM^{NOSC} \cdot \left(y_{b}^{max(DON)} - y_{b}^{min(DON)} \right) \right) \right)$

where
- $y_{b}^{min(DON)}$ is the minimum biomass yield of bacterial functional type $b$ growing on DON (`bac_ydonmin`, [mol N biomass (mol DON)<sup>-1</sup>]) 
- $y_{b}^{max(DON)}$ is the maximum biomass yield of bacterial functional type $b$ growing on DON (`bac_ydonmax`, [mol N biomass (mol DON)<sup>-1</sup>]) 
- $DOM^{NOSC}$ is the in situ nominal oxidation state of dissolved organic carbon that is normalized to vary between 0 (most reduced) and 1 (most oxidised) (`f_nosdoc(i,j,k)`, [dimensionless])

From this base biomass yield on N, we can compute growth yields on DOC (`bac1_ydoc(i,j,k)`; `bac2_ydoc(i,j,k)`, $y_{b}^{DOC}$, [mol C biomass (mol DOC)<sup>-1</sup>]), for growth on $O_2$ [mol C biomass (mol DOC)<sup>-1</sup>] and for anaerobic growth on alternative electron acceptors. We do so by first finding the electron potential (`e_dom`; `e_bac`, $\kappa$) per mole of N of the bacterial biomass and the in situ DOM from basic stoichiometry ([Zakem et al., 2020](https://doi.org/10.1038/s41396-019-0523-8)):

$\kappa_{b} = 4 \cdot R_{b}^{C:N} + 1 \cdot R_{b}^{H:N} - 2 \cdot R_{b}^{O:N} - 3$\
$\kappa_{DOM} = 4 \cdot \dfrac{B_{DOM}^{C}}{B_{DOM}^{N}} + 1 \cdot R_{DOM}^{H:N} - 2 \cdot R_{DOM}^{O:N} - 3$

where
- $\dfrac{B_{DOM}^{C}}{B_{DOM}^{N}}$ is the in situ carbon:nitrogen ratio of DOM (`1/dom_N2C`, [mol C (mol N)<sup>-1</sup>])
- $R_{b}^{C:N}$, $R_{b}^{H:N}$ and $R_{b}^{O:N}$ are the carbon:nitroge, hydrogen:nitrogen and oxygen:nitrogen ratios of bacterial functional type $b$ ([mol (mol N)<sup>-1</sup>])
- $R_{DOM}^{H:N}$ and $R_{DOM}^{O:N}$ are the hydrogen:nitrogen and oxygen:nitrogen ratios of DOM ([mol (mol N)<sup>-1</sup>])

These ratios associated with the stoichimometry of bacterial biomass and dissolved organic matter are constants and are informed by reported values from the literature. For dissolved organic matter we set $R_{DOM}^{H:N} = 10.9$ and $R_{DOM}^{O:N} = 2.6$ ([Anderson et al., 1995](https://doi.org/10.1016/0967-0637(95)00072-E). For bacterial biomass we set $R_{b}^{C:N} = 5$, $R_{b}^{H:N} = 7$ and $R_{b}^{O:N} = 2$ ([Zimmerman et al., 2014](https://doi.org/10.1111/1462-2920.12329)). Therefore

$\kappa_{b} = = 20.0$\
$\kappa_{DOM} = 4 \cdot \dfrac{B_{DOM}^{C}}{B_{DOM}^{N}} + 2.7$

Using these electron potentials ($\kappa$) we can solve for the fraction of electrons used for biomass synthesis in both aerobic (`f_bac`, $f_{b}^{aer(\kappa)}$, [dimenionless]) and anaerobic (`f_bac`, $f_{b}^{ana(\kappa)}$, [dimenionless]) growth with

$f_{b}^{aer(\kappa)} = \min\left(0.8, y_{b}^{DON} \dfrac{\kappa_{b}}{\kappa_{DOM}} \right)$
$f_{b}^{ana(\kappa)} = f_{b}^{aer(\kappa)} P_{ana}$

where 
- $P_{ana}$ is a small penalty on growth yields due to anaerobic metabolism (`bacanapen`, [dimensionless]) ([Zakem et al., 2020](https://doi.org/10.1038/s41396-019-0523-8); [Sun et al., 2024](https://doi.org/10.1073/pnas.2417421121))

Now that we have the fraction of electrons that are partitioned towards biomass synthesis for both aerobic and anaerobic metabolisms, we can compute all growth yields:

$y_{b}^{aer(DON)} = y_{b}^{DON} \cdot R_{b}^{C:N}$
$y_{b}^{aer(DOC)} = y_{b}^{DON} \cdot \dfrac{R_{b}^{C:N}}{\dfrac{B_{DOM}^{C}}{B_{DOM}^{N}}}$
$y_{b}^{O_2} = \dfrac{\dfrac{f_{b}^{aer(\kappa)}}{\kappa_{b}}}{\dfrac{1 - f_{b}^{aer(\kappa)}}{4}} \cdot R_{b}^{C:N}$\
$y_{b}^{ana(DOC)} = y_{b}^{DON} P_{ana} \cdot \dfrac{R_{b}^{C:N}}{\dfrac{B_{DOM}^{C}}{B_{DOM}^{N}}}$
$y_{b}^{NO_{3} &rarr; N_{2}O} = \dfrac{\dfrac{f_{b}^{ana(\kappa)}}{\kappa_{b}}}{\dfrac{1 - f_{b}^{ana(\kappa)}}{4}} \cdot R_{b}^{C:N}$\
$y_{b}^{N_{2}O &rarr; N_{2}} = \dfrac{\dfrac{f_{b}^{ana(\kappa)}}{\kappa_{b}}}{\dfrac{1 - f_{b}^{ana(\kappa)}}{1}} \cdot R_{b}^{C:N}$

where
- $y_{b}^{aer(DON)}$ is the aerobic growth yield of bacterial functional type $b$ on DON (`bac1_ydonC`; `bac2_ydon`, [mol C biomass (mol DON)<sup>-1</sup>])
- $y_{b}^{aer(DOC)}$ is the aerobic growth yield of bacterial functional type $b$ on DOC (`bac1_ydoc(i,j,k)`; `bac2_ydoc(i,j,k)`, [mol C biomass (mol DOC)<sup>-1</sup>])
- $y_{b}^{ana(DOC)}$ is the anaerobic growth yield of bacterial functional type $b$ on DOC (`bac1_yana`; `bac2_yana`, [mol C biomass (mol DOC)<sup>-1</sup>])
- $y_{b}^{O_2}$ is the aerobic growth yield of bacterial functional type $b$ on $O_2$ (`bac1_yoxyC`; `bac2_yoxyC`, [mol C biomass (mol $O_2$)<sup>-1</sup>])
- $y_{b}^{NO_{3} &rarr; N_{2}O}$ is the anaerobic growth yield of bacterial functional type $b$ on $NO_3$ (`bac1_yno3C`, [mol C biomass (mol $NO_3$)<sup>-1</sup>])
- $y_{b}^{N_{2}O &rarr; N_{2}}$ is the anaerobic growth yield of bacterial functional type $b$ on $N_{2}O$ (`bac2_yn2oC`, [mol C biomass (mol $N_{2}O$)<sup>-1</sup>])


**Consumption of substrates**\
Heterotrophic bacterial growth consumes $B_{DOM}^{C}$, $B_{DOM}^{N}$, $NH_4$, $dFe$, $O_2$, $NO_3$ and $N_{2}O$. The consumption of a resource $R$ by bacterial functional type $b$ is calculated as:

$\mu_{b}^{&larr; R} = \dfrac{\mu_{b}^{&larr; C}}{y_{b}^{R}}$

where
- $\mu_{b}^{&larr; C}$ is the realized biomass growth rate (i.e., integration of carbon into biomass) of bacterial functional type $b$ ([mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $y_{b}^{R}$ is the growth yield of bacterial functional type $b$ on resource $R$ ([mol C biomass (mol R)<sup>-1</sup>]) 

We must be careful to accommodate the different growth yields asssociated with aerobic and anaerobic growth since these bacterial types are facultatively anaerobic, which means that:

$\mu_{b}^{&larr; R} = \dfrac{\mu_{b}^{&larr; C}}{y_{b}^{aer(R)}} \left(1 - f_{ana} \right) + \dfrac{\mu_{b}^{&larr; C}}{y_{b}^{ana(R)}} f_{ana}$

where
- $y_{b}^{aer(R)}$ is the aerobic growth yield of bacterial functional type $b$ on resource $R$ ([mol C biomass (mol R)<sup>-1</sup>]) 
- $y_{b}^{ana(R)}$ is the anaerobic growth yield of bacterial functional type $b$ on resource $R$ ([mol C biomass (mol R)<sup>-1</sup>]) 
- $f_{ana}$ is the fraction of growth that is supported by anaerobic metabolism (`bac1_fanaer(i,j,k)`; `bac2_fanaer(i,j,k)`, [dimenionless])

$B_{DOM}^{C}$ uptake (`doc1remi(i,j,k)`; `doc2remi(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]), for example, is calculated as:

$\mu_{b}^{&larr; B_{DOM}^{C}} = \dfrac{\mu_{b}^{&larr; C}}{y_{b}^{aer(DOC)}} \left(1 - f_{ana} \right) + \dfrac{\mu_{b}^{&larr; C}}{y_{b}^{ana(DOC)}} f_{ana}$

However, both $B_{DOM}^{N}$ (`don1remi(i,j,k)`; `don2remi(i,j,k)`, [mol N kg<sup>-1</sup> s<sup>-1</sup>]) and $NH_4$ (`bac1unh4(i,j,k)`; `bac2unh4(i,j,k)`, [mol N kg<sup>-1</sup> s<sup>-1</sup>]) serve the bacterial demand for nitrogen and so we must split uptake between these two resources. We do so using the relative $V_{b}^{DON}$ and $V_{b}^{NH_4}$ which depend on the prescribed affinities of bacteria for these resources and the in situ concentrations of $DON$ and $NH_4$:

$\mu_{b}^{&larr; B_{DOM}^{N}} = \left(\dfrac{\mu_{b}^{&larr; C}}{y_{b}^{aer(DON)}} \left(1 - f_{ana} \right) + \dfrac{\mu_{b}^{&larr; C}}{y_{b}^{ana(DON)}} f_{ana} \right) \dfrac{V_{b}^{DON}}{V_{b}^{DON} + V_{b}^{NH_4}}$\
$\mu_{b}^{&larr; NH_4} = \left(\dfrac{\mu_{b}^{&larr; C}}{y_{b}^{aer(DON)}} \left(1 - f_{ana} \right) + \dfrac{\mu_{b}^{&larr; C}}{y_{b}^{ana(DON)}} f_{ana} \right) \dfrac{V_{b}^{NH_4}}{V_{b}^{DON} + V_{b}^{NH_4}}$

Since bacteria release $NH_4$ back to the environment, computing the uptake of $NH_4$ is only for housekeeping purposes and we only consider uptake of $DON$ in tracer tendency equations (see below).

Consumption of the electron acceptors $O_2$ (`bac1resp(i,j,k)`; `bac2resp(i,j,k)`, [mol $O_2$ kg<sup>-1</sup> s<sup>-1</sup>]), $NO_3$ (`bac1deni(i,j,k)`, [mol N kg<sup>-1</sup> s<sup>-1</sup>]) and $N_{2}O$ (`bac2deni(i,j,k)`, [mol $N_{2}O$ kg<sup>-1</sup> s<sup>-1</sup>]) are also calculated only considering aerobic or anaerobic metabolism:

$\mu_{b}^{&larr; O_2} = \dfrac{\mu_{b}^{&larr; C}}{y_{b}^{O_2}} \left(1 - f_{ana} \right) 
$\mu_{b1}^{&larr; NO_3} = \dfrac{\mu_{b1}^{&larr; C}}{y_{b1}^{NO_3}} f_{ana}$
$\mu_{b2}^{&larr; N_{2}O} = \dfrac{\mu_{b2}^{&larr; C}}{y_{b2}^{N_{2}O}} f_{ana}$

Consumption of dissolved iron is simply:

$\mu_{b}^{&larr; dFe} = \dfrac{\mu_{b}^{&larr; C}}{R_{b}^{C:Fe}}$

where
- $R_{b}^{C:Fe}$ is the ratio of carbon to iron in the biomass of the bacterial functional type (`bac1_C2Fe`; `bac2_C2Fe`, [mol C (mol Fe)<sup>-1</sup>])


**Production of substrates**\
Heterotrophic bacterial growth produces $DIC$, $NH_4$ and $N_{2}O$. Because we carry ecosystem components in units of carbon, the production of $DIC$ by bacterial functional type $b$ is calculated as:

$\mu_{b}^{&rarr; DIC} = \mu_{b}^{&larr; C} \left(\dfrac{1}{y_{b}^{DOC}} - 1\right)$

As before, we must accommodate differences in yields between aerobic and anaerobic growth:

$\mu_{b}^{&rarr; DIC} = \mu_{b}^{&larr; C} \left(\dfrac{1}{y_{b}^{(aer)DOC}} - 1\right)\left(1 - f_{ana}\right) + \mu_{b}^{&larr; C} \left(\dfrac{1}{y_{b}^{(ana)DOC}} - 1\right) f_{ana}$

Meanwhile, production of $NH_4$ requires knowledge of how much $DON$ was assimilated by the cell. Assimilated $DON$ is converted into $NH_4$ for biosynthesis and any excess is exuded into the environment. Release of $NH_4$ is therefore: 

$\mu_{b}^{&rarr; NH_4} = \mu_{b}^{&larr; B_{DOM}^{N}} - \dfrac{\mu_{b}^{&larr; C}}{R_{b}^{C:N}}$

In the scenario where $V_{b}^{DON} << V_{b}^{NH_4}$, either due to very low concentrations of $DON$ or poor affinity of the bacteria for $DON$ uptake, $\mu_{b}^{&rarr; NH_4}$ will be negative. In this case, this negative release of $NH_4$ functions as a consumption of $NH_4$.

Finally, the $NO_3$-reducing bacterial functional type produces $N_{2}O$. We carry $N_{2}O$ as an explicit tracer and in units of mol $N_{2}O$ kg<sup>-1</sup>, such that at every conversion of $NO_3$ to $N_{2}O$ we must account for the necessity of 2 mol $NO_3$ for every 1 mol $N_{2}O$. Production of $N_{2}O$ via $NO_3$ reduction is calcualted as:

$\mu_{b1}^{&rarr; N_{2}O} = \dfrac{\mu_{b1}^{&larr; NO_3}}{2}$

---


### 17. Chemoautotrophy.

We consider two forms of chemoautotrophy carried out by two distinct forms of microbes: ammonia oxidizing archaea and anaerobic ammonia oxidizing (anammox) bacteria. Ammonia oxidizing archaea are considered explicitly within WOMBAT-mid (`f_aoa(i,j,k)`, [mol C kg<sup>-1</sup>]), while anammox bacteria are considered implicitly and therefore do not have varying biomasses (i.e., we only compute rates of anammox).

**Ammonia oxidizing archaea**\
Growth of ammonia oxidizing archaea (`aoagrow(i,j,k)`, $\mu_{aoa}^{C}$, [mol C kg<sup>-1</sup>]) is defined similarly to other microbes:

$\mu_{aoa}^{C} = \mu_{aoa} B_{aoa}^{C}$

where
- $\mu_{aoa}$ is the realized growth rate of ammonia oxidizing archaea (`aoa_mu(i,j,k)`, [s<sup>-1</sup])
- $B_{aoa}^{C}$ is the in situ concentration of carbon biomass of ammonia oxidizing archaea (`f_aoa(i,j,k)`, [mol C kg<sup>-1</sup>])

The realized growth rate, $\mu_{aoa}$, is the minimum growth achievable on oxygen and ammonium:

$\mu_{aoa} = \min\left(\mu_{aoa}^{NH_4}, \mu_{aoa}^{O_2}\right)$\
$\mu_{aoa}^{NH_4} = \mu_{aoa}^{max} \dfrac{NH_4}{NH_4 + K_{aoa}^{NH_4}}$\
$\mu_{aoa}^{O_2} = \dfrac{\rho_{aoa}^{O_2} O_2}{y_{aoa}^{O_2}}$

where
- $\mu_{aoa}^{max}$ is a temperature-dependent maximum growh rate of ammonia oxidizing archaea (`aoa_mumax(i,j,k)`, [s<sup>-1</sup])
- $K_{aoa}^{NH_4}$ is the half-saturation coefficient for uptake of $NH_4$ by ammonia oxidizing archaea (`aoa_knh4`, [mmol N m<sup>-3</sup>])
- $\rho_{aoa}^{O_2}$ is the diffusive uptake limit of $O_2$ by ammonia oxidizing archaea (`aoa_poxy`, [(mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>])
- $y_{aoa}^{O_2}$ is the aerobic growth yield of ammonia oxidizing archaea on $O_2$ (`aoa_yoxy`, [mol C biomass (mol $O_2$)<sup>-1</sup>])
- $NH_{4}$ is the in situ concentration of $NH_{4}$ (`bionh4`, [mmol N m<sup>-3</sup>])
- $O_{2}$ is the in situ concentration of $O_{2}$$ (`biooxy`, [mmol $O_{2}$ m<sup>-3</sup>])

The temperature-dependent maximum growth rate of ammonia oxidizing archaea is informed by the cultures of [Qin et al. (2015)](https://doi.org/10.1073/pnas.1501568112):

$\mu_{aoa}^{max} = \dfrac{\max\left(0.2, 0.029 \cdot T - 0.147 \right)}{86400}$

where
- $T$ is the in situ temperature of seawater (`Temp(i,j,k)`, [ºC])

In reality, ammonia oxidizing archaea perform the first step of the nitrification process by oxidizing ammonia through to nitrite. However, in WOMBAT-mid we do not consider nitrite oxidizing bacteria that then complete the second step of the nitrification process to produce nitrate. Hence, in this version of WOMBAT-mid we consider ammonia oxidizing archaea to perform full nitrification and oxidize $NH_4$ direclty to $NO_3$. Consumption of $NH_4$ (`ammox(i,j,k)`, [mol N kg<sup>-1</sup> s<sup>-1</sup>]) and $O_2$ (`aoaresp(i,j,k)`, [mol $O_2$ kg<sup>-1</sup> s<sup>-1</sup>]) are calculated as:

$\mu_{aoa}^{&larr; NH_4} = \dfrac{\mu_{aoa}^{C}}{y_aoa^{NH_4}}$\
$\mu_{aoa}^{&larr; O_2} = \dfrac{\mu_{aoa}^{C}}{y_aoa^{O_2}}$

where
- $y_aoa^{NH_4}$ is the aerobic growth yield of ammonia oxidizing archaea on $NH_4$ (`aoa_ynh4`, [mol C biomass (mol $NH_4$)<sup>-1</sup>])
- $y_{aoa}^{O_2}$ is the aerobic growth yield of ammonia oxidizing archaea on $O_2$ (`aoa_yoxy`, [mol C biomass (mol $O_2$)<sup>-1</sup>])

Ammonia ozidizing archaea produce both $NO_3$ and a small amount of $N_{2}O$ as they oxidize $NH_4$. The production of $N_{2}O$ is informed by numerous studies that show a relationship between the rate of ammonia oxidation, the ambient oxygen concentration and the production of $N_{2}O$ ([Goreau et al., 1980](https://doi.org/10.1128/aem.40.3.526-532.1980); [Santoro et al., 2011](https://www.science.org/doi/full/10.1126/science.1208239); [Qin et al., 2017](https://doi.org/10.1111/1758-2229.12525); [Ji et al., 2018](https://doi.org/10.1029/2018GB005887); [Frey et al., 2023](https://doi.org/10.1002/lno.12283)). We compute these production terms as:

$\mu_{aoa}^{&rarr; NO_3} = \mu_{aoa}^{&larr; NH_4} - \dfrac{\mu_{aoa}^{C}}{R_{aoa}^{C:N}} - \mu_{aoa}^{C} \cdot 2 p_{aoa}^{N_{2}O}$\
$\mu_{aoa}^{&rarr; N_{2}O} = \mu_{aoa}^{C} p_{aoa}^{N_{2}O}$

where
- $R_{aoa}^{C:N}$ is the ratio of carbon to nitrogen within the biomass of ammonia oxidizing archaea (`aoa_C2N`, [mol C (mol N)<sup>-1</sup>])
- $p_{aoa}^{N_{2}O}$ is the production yield of $N_{2}O$ by ammonia oxidizing archaea (`aoa_yn2o`, [mol $N_{2}O$ (mol C biomass)<sup>-1</sup>])

Determining $p_{aoa}^{N_{2}O}$ is key to determining the fraction of $NH_4$ that is routed to $N_{2}O$ during ammonia oxidation. For this we use the oxygen-dependent relationship identified by [Frey et al. (2023)](https://doi.org/10.1002/lno.12283) who found a maximum of 3% per mol of $NO_2$ produced:

$p_{aoa}^{\%N_{2}O} = \min\left(0.03, \dfrac{0.002}{O_2} + p_{aoa}^{min(N_{2}O)} \right)$

Because [Frey et al. (2023)](https://doi.org/10.1002/lno.12283) give the production yield of $N_{2}O$ in terms of % per mol $NO_2$ produced, we convert this through to units of [mol $N_{2}O$ (mol C biomass)<sup>-1</sup>] via:

$p_{aoa}^{N_{2}O} = \dfrac{p_{aoa}^{\%N_{2}O} \cdot \left(y_{aoa}^{NH_4} - dfrac{1}{R_{aoa}^{C:N}} \right)}{2 \cdot p_{aoa}^{\%N_{2}O} + 1}$ 

since

$a\cdot NH_4 + b \cdot O_2 &rarr; c \cdot B_{aoa}^{C} + d \cdot N_{2}O + e \cdot NO_{3}$\
$Y = \% N_{2}O produced per NO_{3} produced$\
$d = \dfrac{\left(a - c\right) \cdot Y}{2 \cdot Y + 1}$


**Anaerobic ammonia oxidizing (Anammox) bacteria**\
Anammox bacteria are considered to be an implicit population within WOMBAT-mid and we do not track variations in their biomass. Rather then computing growth of anammox bacteria we therefore compute rates of anammox, which convert $NH_4$ to $N_2$. As with $N_{2}O$-reducing heterotrophic bacteria, this nitrogen is then permanently lost from the ocean. When `do_anammox = .true.`, we perform this metabolism as:

$\mu_{aox}^{NH_4 &rarr; N_2} = \mu_{aox}^{max} \left(β_{hete}\right)^{T} f_{ana} L_{aox}^{NH_4}$

- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless])
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC])
- $f_{ana}$ is the fraction of growth that is supported by anaerobic metabolism (`bac1_fanaer(i,j,k)`, [dimenionless])
- $L_{aox}^{NH_4}$ is the growth limiter of anammox associated with $NH_4$ availability (`aox_lnh4(i,j,k)`, [dimensionless])

Note that anammox is considered to be present only when anaerobic metabolisms are ocurring. While anammox bacteria can perform anammox in oxygenated and deoxygenated environments, this metabolism is only appreciably measured in deoxygenated environments due to reduced competition with ammonia oxidizing archaea for a limited supply of $NH_4$. Because we do not resolve this competition explicitly, we apply $f_{ana}$ here. The growth limiter due to ammonium availability is a simple michealis-menten limitation function:

$L_{aox}^{NH_4} = \dfrac{NH_4}{NH_4 + K_{aox}^{NH_4}}$

where
- $K_{aoa}^{NH_4}$ is the half-saturation coefficient for uptake of $NH_4$ by anammox bacteria (`aox_knh4`, [mmol N m<sup>-3</sup>])
- $NH_{4}$ is the in situ concentration of $NH_{4}$ (`bionh4`, [mmol N m<sup>-3</sup>])

---


### 18. Nominal oxidation state of dissolved organic carbon.

In addition to DOC and DON, dissolved organic matter is also affected by changes to the nominal oxidation state of carbon, which we carry as a tracer (`f_nosdoc(i,j,k)`, $DOM^{NOSC}$, [dimenionless]). The nominal oxidation state of carbon (NOSC) is estimated by [La Rowe & Van Cappellen (2011)](https://doi.org/10.1016/j.gca.2011.01.020) from stoichiometry:

$NOSC = 4 - \dfrac{\left( 4 \cdot C + H - 3 \cdot N - 2 \cdot O + 5\cdot P - 2 \cdot S\right)}{C}$

This characteristic describes how oxidized or reduced the average carbon atoms are within organic molecules. The NOSC acts as a measure of the energetic potential of organic matter oxidation, where the Gibbs energy released is correlated to the NOSC ([La Rowe & Van Cappellen, 2011](https://doi.org/10.1016/j.gca.2011.01.020)). The more reduced the molecule, the more negative the NOSC value and the more energy is held within the material. The scale varies from -4 to +4, with $CH_4$ representing fully reduced carbon and $CO_2$ representing fully oxidized carbon. 

In WOMBAT-mid, we do not carry information of changes in the stoichiometry of marine organic matter, either dissolved or particulate. We therefore estimate how certian processes affect the NOSC of DOM ($DOM^{DOSC}$). To do so, we apply changes in $DOM^{NOSC}$ occur as a result of changes to the $B_{DOM}^{C}$ pool. Sources of $B_{DOM}^{C}$ alter $DOM^{NOSC}$ via:

$\dfrac{\partial DOM^{NOSC}}{\partial t} = \dfrac{\partial B_{DOM}^{C}}{\partial t} \dfrac{NOSC_{source} - NOSC}{B_{DOM}^{C}}$

where
- $B_{DOM}^{C}$ is the in situ concentration of dissolved organic carbon (`f_doc(i,j,k)`, [mol C kg<sup>-1</sup>])
- $DOM^{NOSC}$ is the in situ value of the normalized nominal oxidation state of carbon (`f_nosdoc(i,j,k)`, [dimensionless])
- $DOM_{source}^{NOSC}$ is the source value of $DOM^{NOSC}$ added to the $B_{DOM}^{C}$ pool during a given process ([dimenionless])

Sources of $B_{DOM}^{C}$ include (1) phytoplankton overflow production (`nosdoc_overflow(i,j,k)`, $\Delta DOM_{overflow}^{NOSC}$, [NOSC s<sup>-1</sup>]), (2) excretion by zooplankton (`nosdoc_excretion(i,j,k)`, $\Delta DOM_{excretion}^{NOSC}$, [NOSC s<sup>-1</sup>]), (3) phytoplankton cell death and lysis (`nosdoc_phylysis(i,j,k)`, $\Delta DOM_{photolyse}^{NOSC}$, [NOSC s<sup>-1</sup>]), (4) bacterial cell death and lysis (`nosdoc_baclysis(i,j,k)`, $\Delta DOM_{bacterlyse}^{NOSC}$, [NOSC s<sup>-1</sup>]), and (5) hydrolysation of particulate organic matter (`nosdoc_dethydro(i,j,k)`, $\Delta DOM_{dethydro}^{NOSC}$, [NOSC s<sup>-1</sup>]). The equations describing the effect on $DOM^{NOSC}$ of each process are:

(1) $\Delta DOM_{overflow}^{NOSC} = \left(\sum_{p} \mu_{p}^{&rarr; DOC} \right) \cdot \dfrac{NOSC_{overflow} - DOM^{NOSC}}{B_{DOM}^{C}}$\
(2) $\Delta DOM_{excretion}^{NOSC} = \left(\sum_{z,i} \left( X_{z}^{&rarr; i^{C}} \cdot f_{z}^{X &rarr; DOM} \right) \right) \cdot \dfrac{NOSC_{excretion} - DOM^{NOSC}}{B_{DOM}^{C}}$\
(3) $\Delta DOM_{photolyse}^{NOSC} = \left(\sum_{p} \gamma_{p}^{&rarr; C} \right) \cdot \dfrac{NOSC_{phytolyse} - DOM^{NOSC}}{B_{DOM}^{C}}$\
(4) $\Delta DOM_{bacterlyse}^{NOSC} = \left(\sum_{b} \gamma_{b}^{&rarr; C} \right) \cdot \dfrac{NOSC_{bacterlyse} - DOM^{NOSC}}{B_{DOM}^{C}}$\
(5) $\Delta DOM_{dethydro}^{NOSC} = \left(\sum_{d} \Gamma_{d}^{&rarr; C} \right) \cdot \dfrac{NOSC_{dethydro} - DOM^{NOSC}}{B_{DOM}^{C}}$

where
- $\mu_{p}^{&rarr; DOC}$ is the overflow production of DOC by phytoplankton type $p$ ([mol C kg<sup>-1</sup>])
- $X_{z}^{&rarr; i^{C}}$ is the excretion of carbon by zooplankton type $z$ feeding on prey type $i$ ([mol C kg<sup>-1</sup>])
- $f_{z}^{X &rarr; DOM}$ is the fraction of excreted carbon that is directed to DOC by zooplankton type $z$ ([mol C (mol C)<sup>-1</sup>])
- $\gamma_{p}^{&rarr; C}$ is the linear mortality rate of phytoplankton type $p$ ([mol C kg<sup>-1</sup>])
- $\gamma_{b}^{&rarr; C}$ is the linear mortality rate of bacteria type $b$ (includes ammonia oxidizing archaea) ([mol C kg<sup>-1</sup>])
- $\Gamma_{d}^{&rarr; C}$ is the hydrolysation rate of particulate carbon type $d$ ([mol C kg<sup>-1</sup>])
- $NOSC_{overflow}$ is the NOSC associated with overflow production of DOC by phytoplankton (`noscphyover`, [dimenionless])
- $NOSC_{excretion}$ is the NOSC associated with excretion of DOC by zooplankton (`nosczooexcr`, [dimenionless])
- $NOSC_{phytolyse}$ is the NOSC associated with lysis of cell contents of phytoplankton (`noscphylyse`, [dimenionless])
- $NOSC_{bacterlyse}$ is the NOSC associated with lysis of cell contents of bacteria and archaea (`noscbaclyse`, [dimenionless])
- $NOSC_{dethydro}$ is the NOSC associated with hydrolysation of sinking particulate detritus (`noscdethydr`, [dimenionless])
- $DOM^{NOSC}$ is the in situ normalized value of NOSC (`f_nosdoc(i,j,k)`, [dimenionless])
- $B_{DOM}^{C}$ is the in situ concentration of DOC (`f_doc(i,j,k)`, [mol C kg<sup>-1</sup>])

In WOMBAT-mid the only sink of $B_{DOM}^{C}$ is consumption by bacteria (`doc1remi(i,j,k)`; `doc2remi(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]). We acknowledge that abiotic photooxidation of DOC does occur in sunlit waters and can be a substantial contributor to total DOC oxidation to DIC ([Cory et al., 2014](https://www.science.org/doi/full/10.1126/science.1253119); [Aarnos et al., 2018](https://doi.org/10.1002/2017GB005698)). Moreover, the photo-oxidation of DOC is shown to make recalcitrant, high-molecular weight molecules more bioavailable to microbes ([Gonsior et al., 2014](https://doi.org/10.1016/j.marchem.2014.04.002)), which seems an important process to capture. However, in this version of WOMBAT-mid we include only bacteria-mediated breakdown of DOC to DIC and its effects on NOSC and leave photo-oxidation for a later version.

Consumption of $B_{DOM}^{C}$ by bacteria alters $DOM^{NOSC}$ (`nosdoc_docconsu(i,j,k)`, $\Delta DOM_{DOCconsume}^{NOSC}$, [NOSC s<sup>-1</sup>]) via:

$\Delta DOM_{DOCconsume}^{NOSC} = \left(\sum_{b} \mu_{b}^{&larr; B_{DOM}^{C}} \right) \cdot \dfrac{- NOSC_{DOCconsume}}{B_{DOM}^{C}}$\

where
- $\mu_{b}^{&larr; B_{DOM}^{C}}$ is the consumption of DOC by bacterial type $b$ (`doc1remi(i,j,k)`; `doc2remi(i,j,k)`, [mol C kg<sup>-1</sup>])
- $NOSC_{DOCconsume}$ is the effect (i.e., offset) to existing NOSC associated with bacterial reworking of DOC (`noscdocproc`, [dimenionless])
- $B_{DOM}^{C}$ is the in situ concentration of DOC (`f_doc(i,j,k)`, [mol C kg<sup>-1</sup>])

---


### 19. Tracer tendencies

**Nitrate** (`f_no3(i,j,k)`, $NO_3$, [mol N kg<sup>-1</sup>])

$\dfrac{\Delta NO_3}{\Delta t} = \mu_{aoa}^{&rarr; NO_3} 
                               - \mu_{b1}^{&larr; NO_3} 
                               - \left( \mu_{np}^{&larr; C} \dfrac{L_{np}^{NO_3}}{L_{np}^{N}} 
                                      + \mu_{mp}^{&larr; C} \dfrac{L_{mp}^{NO_3}}{L_{mp}^{N}} \right) \cdot \dfrac{16}{122}$

**Ammonium** (`f_nh4(i,j,k)`, $NH_4$, [mol N kg<sup>-1</sup>])

$\dfrac{\Delta NH_4}{\Delta t} = \left( X_{mz}^{&larr; B_{np}^{N}}
                                      + X_{mz}^{&larr; B_{mp}^{N}}
                                      + X_{mz}^{&larr; B_{sd}^{N}}
                                      + X_{mz}^{&larr; B_{b1}^{N}}
                                      + X_{mz}^{&larr; B_{b2}^{N}}
                                      + X_{mz}^{&larr; B_{aoa}^{N}} \right) \left(1 - f_{mz}^{X &rarr; DOM} \right)
                               + \left( X_{Mz}^{&larr; B_{np}^{N}}
                                      + X_{Mz}^{&larr; B_{mp}^{N}}
                                      + X_{Mz}^{&larr; B_{sd}^{N}}
                                      + X_{Mz}^{&larr; B_{ld}^{N}}
                                      + X_{Mz}^{&larr; B_{mz}^{N}}
                                      + X_{Mz}^{&larr; B_{b1}^{N}}
                                      + X_{Mz}^{&larr; B_{b2}^{N}}
                                      + X_{Mz}^{&larr; B_{aoa}^{N}} \right) \cdot \left(1 - f_{Mz}^{X &rarr; DOM} \right)
                               + \left( \gamma_{mz}^{&rarr; C} 
                                      + \gamma_{Mz}^{&rarr; C} \right) \cdot \dfrac{16}{122}
                               + \mu_{diazo}^{&rarr; NH_4}
                               + \mu_{b1}^{&rarr; NH_4} 
                               + \mu_{b2}^{&rarr; NH_4}
                               - \mu_{aox}^{NH_4 &rarr; N_2}
                               - \mu_{aoa}^{&larr; NH_4}
                               - \left( \mu_{np}^{&larr; C} \dfrac{L_{np}^{NH_4}}{L_{np}^{N}} 
                                      + \mu_{mp}^{&larr; C} \dfrac{L_{mp}^{NH_4}}{L_{mp}^{N}} \right) \cdot \dfrac{16}{122}$

**Silicic acid** (`f_sil(i,j,k)`, $Si(OH)_{4}$, [mol Si kg<sup>-1</sup>])

$\dfrac{\Delta Si(OH)_{4}}{\Delta t} = \left( \gamma_{mp}^{&rarr; C} 
                                            + g_{mz}^{&larr; B_{mp}^{C}} \right) \cdot Q_{mp}^{Si:C}
                                     + D_{B_{ld}^{Si}}^{&rarr; Si}
                                     - \mu_{mp}^{&larr; Si}$


**Nitrous oxide** (`f_n2o(i,j,k)`, $N_{2}O$, [mol $N_{2}O$ kg<sup>-1</sup>])

$\dfrac{\Delta N_{2}O}{\Delta t} = \mu_{aoa}^{&rarr; N_{2}O} 
                                 + \mu_{b1}^{&rarr; N_{2}O}
                                 - \mu_{b2}^{&larr; N_{2}O}$


**Oxygen** (`f_o2(i,j,k)`, $O_2$, [mol O<sub>2</sub> kg<sup>-1</sup>])

$\dfrac{\Delta O_2}{\Delta t} = \left( X_{mz}^{&larr; C} \left(1 - f_{mz}^{X &rarr; DOM} \right)
                                     + X_{Mz}^{&larr; C} \left(1 - f_{Mz}^{X &rarr; DOM} \right)
                                     + \gamma_{mz}^{&rarr; C} 
                                     + \gamma_{Mz}^{&rarr; C}
                                     - \mu_{np}^{&larr; C} 
                                     - \mu_{mp}^{&larr; C} \right) \dfrac{-132}{122} 
                              - \mu_{b1}^{&larr; O_2} 
                              - \mu_{b2}^{&larr; O_2} 
                              - \mu_{aoa}^{&larr; O_2} \right)$
 

**Dissolved iron** (`f_fe(i,j,k)`, $dFe$, [mol Fe kg<sup>-1</sup>])

$\dfrac{\Delta dFe}{\Delta t} = \Gamma_{sd}^{&rarr; C} Q_{sd}^{Fe:C} 
                              + \Gamma_{ld}^{&rarr; C} Q_{ld}^{Fe:C} 
                              + \gamma_{np}^{&rarr; C} Q_{np}^{Fe:C} 
                              + \gamma_{mp}^{&rarr; C} Q_{mp}^{Fe:C}
                              + \gamma_{mz}^{&rarr; C} Q_{mz}^{Fe:C} 
                              + \gamma_{Mz}^{&rarr; C} Q_{Mz}^{Fe:C}
                              + \dfrac{\gamma_{b1}^{&rarr; C} + \Gamma_{b1}^{&rarr; C}}{R_{b1}^{C:Fe}}
                              + \dfrac{\gamma_{b2}^{&rarr; C} + \Gamma_{b2}^{&rarr; C}}{R_{b2}^{C:Fe}}
                              + \dfrac{\gamma_{aoa}^{&rarr; C} + \Gamma_{aoa}^{&rarr; C}}{R_{aoa}^{C:Fe}}
                              + X_{mz}^{&larr; Fe} 
                              + X_{Mz}^{&larr; Fe}
                              + D_{Fe_{sA}}^{&rarr; dFe} 
                              + D_{Fe_{lA}}^{&rarr; dFe}
                              - \mu_{np}^{&larr; dFe} 
                              - \mu_{mp}^{&larr; dFe}
                              - \mu_{b1}^{&larr; dFe}
                              - \mu_{b2}^{&larr; dFe}
                              - \mu_{aoa}^{&larr; dFe}
                              - Sc_{dFe}^{&rarr; Fe_{sA}} 
                              - Sc_{dFe}^{&rarr; Fe_{lA}}
                              - Co_{dFe}^{&rarr; Fe_{sA}} 
                              - Co_{dFe}^{&rarr; Fe_{lA}}$


**Small authigenic iron** (`f_afe(i,j,k)`, $Fe_{sA}$, [mol Fe kg<sup>-1</sup>])

$\dfrac{\Delta Fe_{sA}}{\Delta t} = Sc_{dFe}^{&rarr; Fe_{sA}}
                                  + Co_{dFe}^{&rarr; Fe_{sA}}
                                  - D_{Fe_{sA}}^{&rarr; dFe}$


**Large authigenic iron** (`f_bafe(i,j,k)`, $Fe_{lA}$, [mol Fe kg<sup>-1</sup>])

$\dfrac{\Delta Fe_{lA}}{\Delta t} = Sc_{dFe}^{&rarr; Fe_{lA}}
                                  + Co_{dFe}^{&rarr; Fe_{lA}}
                                  - D_{Fe_{lA}}^{&rarr; dFe}$


**Nano-phytoplankton** (`f_phy(i,j,k)`, $B_{np}^{C}$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta B_{np}^{C}}{\Delta t} = \mu_{np}^{&larr; C} 
                                     - \Gamma_{np}^{&rarr; C} 
                                     - \gamma_{np}^{&rarr; C} 
                                     - g_{mz}^{&larr; B_{np}^{C}} 
                                     - g_{Mz}^{&larr; B_{np}^{C}}$ 


**Nano-phytoplankton chlorophyll** (`f_pchl(i,j,k)`, $B_{np}^{chl}$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta B_{np}^{Chl}}{\Delta t} = \mu_{np}^{&larr; Chl} 
                                       - \left( \Gamma_{np}^{&rarr; C} 
                                              + \gamma_{np}^{&rarr; C} 
                                              + g_{mz}^{&larr; B_{np}^{C}} 
                                              + g_{Mz}^{&larr; B_{np}^{C}} \right) \cdot Q_{np}^{Chl:C}$ 

**Nano-phytoplankton iron** (`f_phyfe(i,j,k)`, $B_{np}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$\dfrac{\Delta B_{np}^{Fe}}{\Delta t} = \mu_{np}^{&larr; dFe} 
                                      - \left( \Gamma_{np}^{&rarr; C} 
                                             + \gamma_{np}^{&rarr; C} 
                                             + g_{mz}^{&larr; B_{np}^{C}} 
                                             + g_{Mz}^{&larr; B_{np}^{C}} \right) \cdot Q_{np}^{Fe:C}$ 


**Micro-phytoplankton** (`f_dia(i,j,k)`, $B_{mp}^{C}$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta B_{mp}^{C}}{\Delta t} = \mu_{mp}^{&larr; C} 
                                     - \Gamma_{mp}^{&rarr; C} 
                                     - \gamma_{mp}^{&rarr; C} 
                                     - g_{mz}^{&larr; B_{mp}^{C}} 
                                     - g_{Mz}^{&larr; B_{mp}^{C}} $ 


**Micro-phytoplankton chlorophyll** (`f_dchl(i,j,k)`, $B_{mp}^{chl}$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta B_{mp}^{Chl}}{\Delta t} = \mu_{mp}^{&larr; Chl} 
                                       - \left( \Gamma_{mp}^{&rarr; C} 
                                              + \gamma_{mp}^{&rarr; C} 
                                              + g_{mz}^{&larr; B_{mp}^{C}} 
                                              + g_{Mz}^{&larr; B_{mp}^{C}} \right) \cdot Q_{mp}^{Chl:C}$ 

**Micro-phytoplankton iron** (`f_diafe(i,j,k)`, $B_{mp}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$\dfrac{\Delta B_{mp}^{Fe}}{\Delta t} = \mu_{mp}^{&larr; dFe} 
                                      - \left( \Gamma_{mp}^{&rarr; C} 
                                             + \gamma_{mp}^{&rarr; C} 
                                             + g_{mz}^{&larr; B_{mp}^{C}} 
                                             + g_{Mz}^{&larr; B_{mp}^{C}} \right) \cdot Q_{mp}^{Fe:C}$ 

**Micro-phytoplankton silica** (`f_diasi(i,j,k)`, $B_{mp}^{Si}$, [mol Si kg<sup>-1</sup>])

$\dfrac{\Delta B_{mp}^{Si}}{\Delta t} = \mu_{mp}^{&larr; Si}
                                      - \left( \Gamma_{mp}^{&rarr; C} 
                                             + \gamma_{mp}^{&rarr; C} 
                                             + g_{mz}^{&larr; B_{mp}^{C}} 
                                             + g_{Mz}^{&larr; B_{mp}^{C}} \right) \cdot Q_{mp}^{Si:C}$ 


**Micro-zooplankton** (`f_zoo(i,j,k)`, $B_{mz}^{C}$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta B_{mz}^{C}}{\Delta t} = \A_{mz}^{&larr; C} 
                                     - \Gamma_{mz}^{&rarr; C} 
                                     - \gamma_{mz}^{&rarr; C}
                                     - g_{Mz}^{&larr; B_{mz}^{C}}$ 


**Micro-zooplankton iron** (`f_zoofe(i,j,k)`, $B_{mz}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$\dfrac{\Delta B_{mz}^{Fe}}{\Delta t} = A_{mz}^{&larr; Fe}
                                      - \left( \Gamma_{mz}^{&rarr; C} 
                                             - \gamma_{mz}^{&rarr; C} 
                                             - g_{Mz}^{&larr; B_{mz}^{C}} \right) \cdot Q_{mz}^{Fe:C}$ 


**Meso-zooplankton** (`f_mes(i,j,k)`, $B_{Mz}^{C}$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta B_{Mz}^{C}}{\Delta t} = \A_{Mz}^{&larr; C}
                                     - \Gamma_{Mz}^{&rarr; C} 
                                     - \gamma_{Mz}^{&rarr; C}$ 


**Meso-zooplankton iron** (`f_mesfe(i,j,k)`, $B_{Mz}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$\dfrac{\Delta B_{Mz}^{Fe}}{\Delta t} = A_{Mz}^{&larr; Fe}
                                      - \left( \Gamma_{Mz}^{&rarr; C} 
                                             - \gamma_{Mz}^{&rarr; C} \right) \cdot Q_{Mz}^{Fe:C}$ 


**Small detritus** (`f_det(i,j,k)`, $B_{sd}^{C}$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta B_{sd}^{C}}{\Delta t} = E_{mz}^{&larr; C}
                                     + \Gamma_{np}^{&rarr; C} 
                                     + \Gamma_{mz}^{&rarr; C} 
                                     - g_{mz}^{&larr; B_{sd}^{C}}
                                     - g_{Mz}^{&larr; B_{sd}^{C}}
                                     - \Gamma_{sd}^{&rarr; C}$ 


**Small detritus iron** (`f_detfe(i,j,k)`, $B_{sd}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$\dfrac{\Delta B_{sd}^{Fe}}{\Delta t} = E_{mz}^{&larr; Fe}
                                      +\Gamma_{np}^{&rarr; C} Q_{np}^{Fe:C} 
                                      + \Gamma_{mz}^{&rarr; C} Q_{mz}^{Fe:C}
                                      - \left( g_{mz}^{&larr; B_{sd}^{C}}
                                             + g_{Mz}^{&larr; B_{sd}^{C}}
                                             + \Gamma_{sd}^{&rarr; C} \right) Q_{sd}^{Fe:C}$ 


**Large detritus** (`f_bdet(i,j,k)`, $B_{ld}^{C}$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta B_{ld}^{C}}{\Delta t} = E_{Mz}^{&larr; C}
                                     + \Gamma_{mp}^{&rarr; C} 
                                     + \Gamma_{Mz}^{&rarr; C} 
                                     - g_{Mz}^{&larr; B_{ld}^{C}}
                                     - \Gamma_{ld}^{&rarr; C}$ 


**Large detritus iron** (`f_bdetfe(i,j,k)`, $B_{ld}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$\dfrac{\Delta B_{ld}^{Fe}}{\Delta t} = E_{Mz}^{&larr; Fe}
                                      + \Gamma_{mp}^{&rarr; C} Q_{mp}^{Fe:C} 
                                      + \Gamma_{Mz}^{&rarr; C} Q_{Mz}^{Fe:C}
                                      - \left( g_{Mz}^{&rarr; B_{ld}^{C}}
                                             + \Gamma_{ld}^{&rarr; C} \right) Q_{ld}^{Fe:C}$ 


**Large detritus silicon** (`f_bdetsi(i,j,k)`, $B_{ld}^{Si}$, [mol Si kg<sup>-1</sup>])

$\dfrac{\Delta B_{ld}^{Si}}{\Delta t} = \left( \Gamma_{mp}^{&rarr; C} 
                                             + g_{Mz}^{&larr; B_{mp}^{C}} \right) Q_{mp}^{Si:C}  
                                      - D_{B_{ld}^{Si}}^{&rarr; Si}$ 


**Faculative $NO_3$-reducing heterotrophic bacteria** (`f_bac1(i,j,k)`, $B_{b1}^{C}$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta B_{b1}^{C}}{\Delta t} = \mu_{b1}^{&larr; C} 
                                     - g_{mz}^{&larr; B_{b1}^{C}}
                                     - g_{Mz}^{&larr; B_{b1}^{C}}
                                     - \gamma_{b1}^{&rarr; C} 
                                     - \Gamma_{b1}^{&rarr; C}$ 


**Faculative $N_{2}O$-reducing heterotrophic bacteria** (`f_bac2(i,j,k)`, $B_{b2}^{C}$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta B_{b1}^{C}}{\Delta t} = \mu_{b2}^{&larr; C} 
                                     - g_{mz}^{&larr; B_{b2}^{C}}
                                     - g_{Mz}^{&larr; B_{b2}^{C}}
                                     - \gamma_{b2}^{&rarr; C} 
                                     - \Gamma_{b2}^{&rarr; C}$ 


**Ammonia oxidizing archaea** (`f_aoa(i,j,k)`, $B_{aoa}^{C}$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta B_{aoa}^{C}}{\Delta t} = \mu_{aoa}^{&larr; C} 
                                      - g_{mz}^{&larr; B_{aoa}^{C}}
                                      - g_{Mz}^{&larr; B_{aoa}^{C}}
                                      - \gamma_{aoa}^{&rarr; C} 
                                      - \Gamma_{aoa}^{&rarr; C}$ 


**Dissolved organic carbon** (`f_doc(i,j,k)`, $B_{DOM}^{C}$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta B_{DOM}^{C}}{\Delta t} = \mu_{np}^{&rarr; DOC} 
                                      + \mu_{mp}^{&rarr; DOC}
                                      + \Gamma_{sd}^{&rarr; C}
                                      + \Gamma_{ld}^{&rarr; C}
                                      + \gamma_{np}^{&rarr; C}
                                      + \gamma_{mp}^{&rarr; C}
                                      + \Gamma_{b1}^{&rarr; C}
                                      + \gamma_{b1}^{&rarr; C}
                                      + \Gamma_{b2}^{&rarr; C}
                                      + \gamma_{b2}^{&rarr; C}
                                      + \Gamma_{aoa}^{&rarr; C}
                                      + \gamma_{aoa}^{&rarr; C}
                                      + X_{mz}^{&larr; C} f_{mz}^{X &rarr; DOM} 
                                      + X_{Mz}^{&larr; C} f_{Mz}^{X &rarr; DOM} 
                                      - \mu_{b1}^{&larr; B_{DOM}^{C}} 
                                      - \mu_{b2}^{&larr; B_{DOM}^{C}}$ 


**Dissolved organic nitrogen** (`f_don(i,j,k)`, $B_{DOM}^{N}$, [mol N kg<sup>-1</sup>])

$\dfrac{\Delta B_{DOM}^{N}}{\Delta t} = \left( \Gamma_{b1}^{&rarr; C} + \gamma_{b1}^{&rarr; C}
                                             + \Gamma_{b2}^{&rarr; C} + \gamma_{b2}^{&rarr; C}
                                             + \Gamma_{aoa}^{&rarr; C} + \gamma_{aoa}^{&rarr; C}
                                             + X_{mz}^{&larr; B_{b1}^{C}}
                                             + X_{mz}^{&larr; B_{b2}^{C}}
                                             + X_{mz}^{&larr; B_{aoa}^{C}}
                                             + X_{Mz}^{&larr; B_{b1}^{C}}
                                             + X_{Mz}^{&larr; B_{b2}^{C}}
                                             + X_{Mz}^{&larr; B_{aoa}^{C}} \right) \cdot \dfrac{1}{5}
                                      + \left( \Gamma_{sd}^{&rarr; C} 
                                             + \Gamma_{ld}^{&rarr; C}
                                             + \gamma_{np}^{&rarr; C}
                                             + \gamma_{mp}^{&rarr; C}
                                             + \left( X_{mz}^{&larr; B_{np}^{C}}
                                                    + X_{mz}^{&larr; B_{mp}^{C}}
                                                    + X_{mz}^{&larr; B_{sd}^{C}} \right) f_{mz}^{X &rarr; DOM}
                                             + \left( X_{Mz}^{&larr; B_{np}^{C}}
                                                    + X_{Mz}^{&larr; B_{mp}^{C}}
                                                    + X_{Mz}^{&larr; B_{sd}^{C}}
                                                    + X_{Mz}^{&larr; B_{ld}^{C}}
                                                    + X_{Mz}^{&larr; B_{mz}^{C}} \right) f_{Mz}^{X &rarr; DOM} \right) \cdot \dfrac{16}{122}
                                      - \mu_{b1}^{&larr; B_{DOM}^{N}} 
                                      - \mu_{b2}^{&larr; B_{DOM}^{N}}$ 


**Nominal oxidation state of dissolved organiccarbon** (`f_nosdoc(i,j,k)`, $DOM^{NOSC}$, [dimenionless])

$\dfrac{\Delta DOM^{NOSC}}{\Delta t} = \Delta DOM_{overflow}^{NOSC}
                                     + \Delta DOM_{excretion}^{NOSC}
                                     + \Delta DOM_{photolyse}^{NOSC}
                                     + \Delta DOM_{bacterlyse}^{NOSC}
                                     + \Delta DOM_{dethydro}^{NOSC}
                                     + \Delta DOM_{DOCconsume}^{NOSC}$


**Calcium Carbonate** (`f_caco3(i,j,k)`, $CaCO_3$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta CaCO_3}{\Delta t} = P_{CaCO_3}^{\Gamma_{np}^{&rarr; C}}
                                 + P_{CaCO_3}^{\Gamma_{mz}^{&rarr; C}}
                                 + P_{CaCO_3}^{g_{mz}^{&larr; B_{np}^{C}}}
                                 + P_{CaCO_3}^{g_{Mz}^{&larr; B_{np}^{C}}}
                                 + P_{CaCO_3}^{g_{Mz}^{&larr; B_{mz}^{C}}}
                                 - D_{CaCO_3}^{\Omega_{cal}}
                                 - D_{CaCO_3}^{\Omega_{ara}}
                                 - D_{CaCO_3}^{\Gamma_{sd}^{&rarr; C}}
                                 - D_{CaCO_3}^{g_{mz}^{&larr; B_{sd}^{C}}}
                                 - D_{CaCO_3}^{g_{Mz}^{&larr; B_{sd}^{C}}}$


**Dissolved Inorganic Carbon** (`f_dic(i,j,k)`, $DIC$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta DIC}{\Delta t} = X_{mz}^{&larr; C} \left(1 - f_{mz}^{X &rarr; DOM} \right)
                              + X_{Mz}^{&larr; C} \left(1 - f_{Mz}^{X &rarr; DOM} \right)
                              + \gamma_{mz}^{&rarr; C}
                              + \gamma_{Mz}^{&rarr; C}
                              + \mu_{b1}^{&rarr; DIC}
                              + \mu_{b2}^{&rarr; DIC}
                              + D_{CaCO_3}^{\Omega_{cal}}
                              + D_{CaCO_3}^{\Omega_{ara}}
                              + D_{CaCO_3}^{\Gamma_{sd}^{&rarr; C}}
                              + D_{CaCO_3}^{g_{mz}^{&larr; B_{sd}^{C}}}
                              + D_{CaCO_3}^{g_{Mz}^{&larr; B_{sd}^{C}}}
                              - \mu_{np}^{&larr; C}
                              - \mu_{mp}^{&larr; C}
                              - \mu_{aoa}^{&larr; C}
                              - \mu_{np}^{&rarr; DOC}
                              - \mu_{mp}^{&rarr; DOC}
                              - P_{CaCO_3}^{\Gamma_{np}^{&rarr; C}}
                              - P_{CaCO_3}^{\Gamma_{mz}^{&rarr; C}}
                              - P_{CaCO_3}^{g_{mz}^{&larr; B_{np}^{C}}}
                              - P_{CaCO_3}^{g_{Mz}^{&larr; B_{np}^{C}}}
                              - P_{CaCO_3}^{g_{Mz}^{&larr; B_{mz}^{C}}}$


**Alkalinity** (`f_alk(i,j,k)`, $Alk$, [mol Eq kg<sup>-1</sup>])

$\dfrac{\Delta Alk}{\Delta t} = \left( \mu_{np}^{&larr; C} \cdot \dfrac{L_{np}^{NO_3}}{L_{np}^{N}} 
                                     + \mu_{mp}^{&larr; C} \cdot \dfrac{L_{mp}^{NO_3}}{L_{mp}^{N}} 
                                     - \mu_{np}^{&larr; C} \cdot \dfrac{L_{np}^{NH_4}}{L_{np}^{N}} 
                                     - \mu_{mp}^{&larr; C} \cdot \dfrac{L_{mp}^{NH_4}}{L_{mp}^{N}} 
                                     + \gamma_{mz}^{&rarr; C} 
                                     + \gamma_{Mz}^{&rarr; C} \right) \cdot \dfrac{16}{122}
                              + \left( X_{mz}^{&larr; B_{np}^{N}}
                                     + X_{mz}^{&larr; B_{mp}^{N}}
                                     + X_{mz}^{&larr; B_{sd}^{N}}
                                     + X_{mz}^{&larr; B_{b1}^{N}}
                                     + X_{mz}^{&larr; B_{b2}^{N}}
                                     + X_{mz}^{&larr; B_{aoa}^{N}} \right) \left(1 - f_{mz}^{X &rarr; DOM} \right)
                              + \left( X_{Mz}^{&larr; B_{np}^{N}}
                                     + X_{Mz}^{&larr; B_{mp}^{N}}
                                     + X_{Mz}^{&larr; B_{sd}^{N}}
                                     + X_{Mz}^{&larr; B_{ld}^{N}}
                                     + X_{Mz}^{&larr; B_{mz}^{N}}
                                     + X_{Mz}^{&larr; B_{b1}^{N}}
                                     + X_{Mz}^{&larr; B_{b2}^{N}}
                                     + X_{Mz}^{&larr; B_{aoa}^{N}} \right) \cdot \left(1 - f_{Mz}^{X &rarr; DOM} \right)
                              + \mu_{b1}^{&rarr; NH_4} 
                              + \mu_{b2}^{&rarr; NH_4}
                              + \mu_{b1}^{&larr; NO_{3}}
                              - \mu_{aoa}^{&larr; NH_{4}}
                              - \mu_{aoa}^{&rarr; NO_{3}}
                              - \mu_{aox}^{NH_4 &rarr; N_2}
                              + \left( D_{CaCO_3}^{\Omega_{cal}}
                                     + D_{CaCO_3}^{\Omega_{ara}}
                                     + D_{CaCO_3}^{\Gamma_{sd}^{&rarr; C}}
                                     + D_{CaCO_3}^{g_{mz}^{&larr; B_{sd}^{C}}}
                                     + D_{CaCO_3}^{g_{Mz}^{&larr; B_{sd}^{C}}}
                                     - P_{CaCO_3}^{\Gamma_{np}^{&rarr; C}}
                                     - P_{CaCO_3}^{\Gamma_{mz}^{&rarr; C}}
                                     - P_{CaCO_3}^{g_{mz}^{&larr; B_{np}^{C}}}
                                     - P_{CaCO_3}^{g_{Mz}^{&larr; B_{np}^{C}}}
                                     - P_{CaCO_3}^{g_{Mz}^{&larr; B_{mz}^{C}}} \right) \cdot 2$


---


### 20. Check for conservation of mass

When checks for the conservation of mass is enabled (`do_check_n_conserve = .true.` or `do_check_c_conserve = .true.` or `do_check_si_conserve = .true.`), the model will calculate the budget of nitrogen or carbon or silicon before and after the ecosystem equations have completed. This checks that the ecosystem equations detailed above have indeed conserved the mass of both nitrogen and carbon within the ocean. In WOMBAT-mid, both nitrogen and carbon and silicon should be perfectly conserved during ecosystem cycling. The only exception to this is for nitrogen, where if any of `do_nitrogen_fixation = .true.`, `do_anammox = .true.`, `do_wc_denitrification = .true.` or `do_benthic_denitrification = .true.` then the model does not and should not be expected to conserve nitrogen.

---

### 21. Additional operations on tracers

**First**, dissolved iron concentrations are set to equal 1 nM everywhere where the depth of the water column is less than 200 metres deep. WOMBAT-mid is not considered to be a model of the coastal ocean, but rather a model of the global pelagic ocean. Given that coastal waters are not limited in dissolved iron due to substantial interactions with sediments and exchange with the land, we universally set the dissolved iron concentration in these waters to 1 nM.

**Second**, if dissolved iron concentrations dip below that measureable by operational detection limits, we reset these concentrations to this minimum  (`zfermin`, $[dFe]^{min}$, [µmol m<sup>-3</sup>]):

$[dFe]^{min} = \min\left( \max\left( 0.003 \cdot [NO_3]^{2}, 0.005 \right), 0.007 \right)$

where:
- $[NO_3]$ is the ambient nitrate concentration in units of [mmol m<sup>-3</sup>] 

This resetting of minimum dFe concentration comes directly from the PISCES ocean model and functions essentially as a constant source of dFe to the ocean when surface concentrations are drawn down to near zero values.

---


### 22. Sinking rate of particulates.

WOMBAT-mid functions with a spatially variable sinking rate of organic detritus (`f_det(i,j,k)`; `f_bdet(i,j,k)`), calcium carbonate (`f_caco3(i,j,k)`) and biogenic silica (`f_bdetsi(i,j,k)`). Sinking of organic iron (`f_detfe(i,j,k)`; `f_bdetfe(i,j,k)`)occurs at the same rate as their respective organic particulate carbon types, while small and large authigenic iron particles (`f_afe(i,j,k)`; `f_bafe(i,j,k)`) sink at their own unique rates. The algorithm to compute sinking rates functions by computing:

1. the average radii of particles in the community;
2. the seawater dynamic viscosity (if `do_viscous_sinking =.true.`);
3. the effect of mineral ballasting, specifically $CaCO_3$ and biogenic silica, on particulate excess density;
4. Rubey's equation for sinking rates of particles.

This approach is inspired by the study of [Dinauer et al. (2022)](https://doi.org/10.1029/2021GB007131). We deal with each of these steps below.

**Average radii of particulates**\
We first estimate the average radius of sinking particles belonging to small and large detritus. Nano-phytoplankton and micro-zooplankton contribute to the small detritus pool, and as such variations in the mean size of these plankton types determine the mean radius of small particles. Similarly, micro-phytoplankton and meso-zooplankton contribute to the large detritus pool, and their sizes determine the mean radius of large particles. According to [Wickman et al. (2024)]((https://doi.org/10.1126/science.adk6901)), the average volume, $V$, of phytoplankton, $p$, in the marine community scales with the biomass density according to:

$V = \left(B_{p}^{C}\right)^{0.65}$

We can relate the radius $r$ in units of µm to the volume of phytoplankton cells via:

$r = 0.5 \cdot \left(\dfrac{6}{\pi} \left(B_{p}^{C}\right)^{0.65} \right)^{\dfrac{1}{3}}$

Therefore, the average radius of nano-phytoplankton (`rad_phy`, [m]) and micro-phytoplankton (`rad_dia`, [m]) is equal to:

$r_{np} = r_{np}^{-} \cdot 0.5 \cdot \left(\dfrac{6}{\pi} \left(B_{np}^{C}\right)^{0.65} \right)^{\dfrac{1}{3}} \cdot 1 \times 10^{-6}$\
$r_{mp} = r_{mp}^{-} \cdot 0.5 \cdot \left(\dfrac{6}{\pi} \left(B_{mp}^{C}\right)^{0.65} \right)^{\dfrac{1}{3}} \cdot 1 \times 10^{-6}$

which simplifies to:

$r_{np} = r_{np}^{-} \left(B_{np}^{C}\right)^{0.217} \cdot 0.62 \times 10^{-6}$\
$r_{mp} = r_{mp}^{-} \left(B_{mp}^{C}\right)^{0.217} \cdot 0.62 \times 10^{-6}$

For micro-zooplankton, we use a relationship presented by [Menden-Deuer & Lessard (2000)](https://doi.org/10.4319/lo.2000.45.3.0569) who identified that the carbon concentration of diverse protists, including heterotrophic dinoflagellates and other micro-zooplankton, was related to their cell volume to the power of 0.939. Hence, we estimate the radius of micro-zooplankton (`rad_zoo`, [m]) from their carbon biomass concentration by inverting this exponent:

$r_{mz} = r_{mz}^{-} \cdot 0.5 \cdot \left(\dfrac{6}{\pi} \left(B_{mz}^{C}\right)^{1.065} \right)^{\dfrac{1}{3}} \cdot 1 \times 10^{-6}$

which simplifies to:

$r_{mz} = r_{mz}^{-} \left(B_{mz}^{C}\right)^{0.355} \cdot 0.62 \times 10^{-6}$

For meso-zooplankton, we assume that the dry carbon biomass scales with the body length to the power of 3, such that $B_{Mz}^{C} \prop L^{3}$ ([Uye, 1982](https://doi.org/10.1007/BF02110286); [Lehette & Hernandez-Leon, 2009](https://doi.org/10.4319/lom.2009.7.304)). Thus, $L \prop \left(B_{Mz}^{C}\right)^{\dfrac{1}{3}}$, and:

$r_{Mz} = r_{Mz}^{-} \cdot 0.5 \cdot \left(B_{Mz}^{C}\right)^{\dfrac{1}{3}} \cdot 1 \times 10^{-6}$

We determine the mean radius of small (`rad_det`, [m]) and large particulate detritus (`rad_bdet`, [m]) by taking the biomass-weighted means of each plankton functional type:

$r_{s} = \dfrac{B_{np}^{C} r_{np} + B_{mz}^{C} r_{mz}}{B_{np}^{C} + B_{mz}^{C}}\
$r_{l} = \dfrac{B_{mp}^{C} r_{mp} + B_{Mz}^{C} r_{Mz}}{B_{mp}^{C} + B_{Mz}^{C}}$

where
- $B_{np}^{C}$ is the concentration of nano-phytoplankton biomass at the surface of the water column (`biophy`, [mmol C m<sup>-3</sup>])
- $B_{mp}^{C}$ is the concentration of micro-phytoplankton biomass at the surface of the water column (`biodia`, [mmol C m<sup>-3</sup>])
- $B_{mz}^{C}$ is the concentration of micro-zooplankton biomass at the surface of the water column (`biozoo`, [mmol C m<sup>-3</sup>])
- $B_{Mz}^{C}$ is the concentration of meso-zooplankton biomass at the surface of the water column (`biomes`, [mmol C m<sup>-3</sup>])


**Seawater dynamic viscosity**\
If `do_viscous_sinking = .true.`, we calculate the dynamic viscosity of the in situ seawater (`dynvis_sw(i,j,k)`, $\eta_{sw}$, [kg m<sup>-1</sup> s<sup>-1</sup>]) that particulates are sinking through. This involves three steps and is dependent on temperature, salinity and pressure.

The dynamic viscosity of seawater at atmospheric pressure ($\eta_{sw}^{1atm}$) is described by equations 22 and 23 from [Sharqawy et al. (2010)](https://doi.org/10.5004/dwt.2010.1079), which are based on [Isdale et al. (1972)](https://doi.org/10.1016/S0011-9164(00)80002-8):

$\eta_{sw}^{1atm} = \eta_{w}^{1atm} \left(1 + A \dfrac{S}{1000} + B \left(\dfrac{S}{1000}\right)^{2} \right)$

where

$A = 1.541 + 1.998 \times 10^{-2} \cdot T - 9.52 \times 10^{-5} \cdot T^{2}$\
$B = 7.974 - 7.561 \times 10^{-2} \cdot T + 4.724 \times 10^{-4} \cdot T^{2}$\
$\eta_{w}^{1atm} = 4.2844 \times 10^{-5} + \left(0.157 \left( T + 64.993 \right)^{2} - 91.296 \right)^{-1}$

and where
- T is in situ seawater temperature (`Temp(i,j,k)`, [ºC])
- S is in situ seawater salinity (`Salt(i,j,k)`, [g kg<sup>-1</sup>])

After calculating $\eta_{sw}^{1atm}$, we must correct for pressure changes in the water column. This is done by calculating the effect of pressure and temperature changes on the dynamic viscosity of pure water ($\eta_{w}$), and then applying this correction to our estimate of seawater dynamic viscosity at atmospheric pressure, such that:

$\eta_{sw} = \eta_{sw}^{1atm} \dfrac{\eta_{w}}{\eta_{w}^{1atm}}$

We solve for $\eta_{w}$, the dynamic viscosity of pure water corrected for pressure effects, by following the [IAPWS (2008)](https://iapws.org/documents/release/viscosity). This formulation requires multiple steps. First, we estimate the density of pure water, $\rho_{w}$, at a given temperature $T$ and pressure $P$ using equation 14 from the UNESCO EOS-80:

$P_{MPa} = \left(101325 + 9.81 * 1025.0 * z\right) 1 \times 10^{-6}$\
$\rho_{0} = 999.842594 + 6.793952 \times 10^{-2} \cdot T - 9.095290 \times 10^{-3} \cdot T^{2} + 1.001685 \times 10^{-4} \cdot T^{3} - 1.120083 \times 10^{-6} \cdot T^{4} + 6.536332 \times 10^{-9} \cdot T^{5}$\
$\rho_{w} = \dfrac{\rho_{0}}{1 - \dfrac{P_{MPa} - 0.101325}{2.2 \times 10^{3}}}$\

where
- $P_{MPa}$ is the in situ pressure at a given depth [`P_MPa`, [MPa]]
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC])

Next, we solve for the dynamic viscosity at the dilute gas-limit, $\hat{\eta_{0}}$, detailed in equation 11 in the [IAPWS (2008)](https://iapws.org/documents/release/viscosity) and with $H_{i}$ coefficients detailed in their Table 1.

$\hat{\eta_{0}} = \dfrac{100 \sqrt{\hat{T}}}{\sum_{i=0}^{3} \dfrac{H_{i}}{\left(\hat{T}\right)^{i}}}$\

where
- $\hat{T} = \dfrac{T + 273.15}{647.096}$ (`T_hat`, [dimenionless])

Next, we solve for the contribution of finite density to dynamic viscosity, $\hat{\eta_{1}}$, detailed in equation 12 in the [IAPWS (2008)](https://iapws.org/documents/release/viscosity) and with $H_{ij}$ coefficients detailed in their Table 2.

$\hat{\eta_{1}} = e^{\left[ \hat{\rho} \sum_{i=0}^{5}\left( \dfrac{1}{\hat{T} - 1} \right)^{i} \sum_{j=0}^{6} H_{ij}\left(\hat{\rho} - 1\right)^{j} \right]}$
 
where
- $\hat{\rho} = \dfrac{\rho_{w}}{322}$ (`rho_hat`, [dimensionless])

Finally, we compute the density-corrected pure water dynamic viscosity:

$\eta_{w} = \left(\hat{\eta_{0}} \cdot \hat{\eta_{1}} \right) \eta^{*}$

where
- $\eta^{*} = 1 \times 10^{-6}$ (`mu_star`, [kg m<sup>-1</sup> s<sup>-1</sup>])

which we apply above to calculate the dynamic viscosity of seawater ($\eta_{sw}$) for a given temperature, salinity and pressure. We note that it is expected that the dynamic viscosity of water actually decreases with increasing pressure at low temperatures ([Percy W. Bridgman (1925)](https://api.semanticscholar.org/CorpusID:27500909)). 


**Mineral ballasting and excess density**\
WOMBAT-mid explicitly considers small organic carbon, large aggregates of organic carbon, $CaCO_3$ and biogenic silica. Each of these particulate types have unique densities. We compute the mass of each particulate type in [kg m</sup>-3</sup>]:

$M_{sd} = B_{sd}^{C} \cdot 1 \times 10{-6} \cdot \dfrac{12}{0.4}$\
$M_{ld} = B_{ld}^{C} \cdot 1 \times 10{-6} \cdot \dfrac{12}{0.4}$\
$M_{CaCO_3} = B_{CaCO_3}^{C} \cdot 1 \times 10{-6} \cdot 100$\
$M_{BSi} = B_{ld}^{Si} \cdot 1 \times 10{-6} \cdot 60$

where
- $B_{sd}^{C}$ is the in situ concentration of small particulate organic carbon (`biodet`, [mmol C m<sup>-3</sup>])
- $B_{ld}^{C}$ is the in situ concentration of large particulate organic carbon (`biobdet`, [mmol C m<sup>-3</sup>])
- $B_{CaCO_3}^{C}$ is the in situ concentration of calcium carbonate carbon (`biocaco3`, [mmol C m<sup>-3</sup>])
- $B_{ld}^{Si}$ is the in situ concentration of biogenic silica (`biobdetsi`, [mmol Si m<sup>-3</sup>])
- $\dfrac{12}{0.4}$ is the g (mol C)<sup>-1</sup> and assuming that 40% of the total biomass of particulate organics is carbon.
- $100$ is the g (mol C)<sup>-1</sup> of calcium carbonate.
- $60$ is the g (mol Si)<sup>-1</sup> of biogenic silica.

We consider $CaCO_3$ to be part of the small sinking particulates because, although more dense than organic matter, $CaCO_3$ particles tend to be smaller than organic aggregates and sink at a slower rate ([De La Rocha & Passow, 2007](https://doi.org/10.1016/j.dsr2.2007.01.004); [Zhang et al., 2018](https://doi.org/10.5194/bg-15-4759-2018)). Furthermore, the shedding of coccoliths by coccolithophores, which are near-neutrally bouyant, also contributes to a slower mean sinking speed of $CaCO_3$ ([Balch et al., 2009](https://doi.org/10.1029/2008JC004902)). In contrast, we consider biogenic silica to be part of the large sinking particulates:

$M_{s} = M_{sd} + M_{CaCO_3}$\
$M_{l} = M_{ld} + M_{BSi}$

And we compute the harmonic mean density of the small particulates (`rho_small`, $\rho_{s}$, [kg m<sup>-3</sup>]) and large particulates (`rho_large`, $\rho_{l}$, [kg m<sup>-3</sup>]) weighted by mass fractions, which accounts for the fact that less dense mass fractions account for greater volume within aggregates:

$\rho_{s} = \dfrac{1}{\dfrac{M_{sd} / M_{s}}{\rho_{orgC}} + \dfrac{M_{CaCO_3} / M_{s}}{\rho_{CaCO_3}} }$\
$\rho_{l} = \dfrac{1}{\dfrac{M_{ld} / M_{l}}{\rho_{orgC}} + \dfrac{M_{Bsi} / M_{l}}{\rho_{Bsi}} }$

where
- $\rho_{orgC}$ is the density of organic carbon particles (`detrho`, [kg m<sup>-3</sup>])
- $\rho_{CaCO_3}$ is the density of calcium carbonate particles (`caco3rho`, [kg m<sup>-3</sup>])
- $\rho_{BSi}$ is the density of biogenic silica particles (`bsirho`, [kg m<sup>-3</sup>])

Finally, we incorporate estimates of particle porosity to their density:

$\rho_{s} = \left(1 - p_{s}\right) \rho_{s} + p_{s} \rho_{sw} $\
$\rho_{l} = \left(1 - p_{l}\right) \rho_{l} + p_{l} \rho_{sw} $

where
- $p_{s}$ is the porosity of small particles (`detphi`, [dimensionless])
- $p_{l}$ is the porosity of large particles (`bdetphi`, [dimensionless])
- $\rho_{sw}$ is the density of seawater, which we set here to a constant 1025 (kg m<sup>-1</sup>)


**Rubey's equation**\
Rubey's equation ([Rubey, 1933](https://doi.org/10.2475/ajs.s5-25.148.325)) combines the radius of particles, their excess density relative to fluid and the dynamic viscosity of that fluid to compute the sinking rate of particles. We find the sinking rate of small (`wsink1(k)`, [m s<sup>-1</sup>]) and large particles (`wsink2(k)`, [m s<sup>-1</sup>]) using Rubey's equation:

$\omega_{s} = \dfrac{ \sqrt{\dfrac{4}{3} 9.8 \cdot \rho_{sw} \left(\rho_{s} - \rho_{sw} \right) \left(r_{s}\right)^{3} + 9 \left(\eta_{sw}\right)^{2}} - 3 \eta_{sw} }{\rho_{sw} \cdot r_{s}}$
$\omega_{l} = \dfrac{ \sqrt{\dfrac{4}{3} 9.8 \cdot \rho_{sw} \left(\rho_{l} - \rho_{sw} \right) \left(r_{l}\right)^{3} + 9 \left(\eta_{sw}\right)^{2}} - 3 \eta_{sw} }{\rho_{sw} \cdot r_{l}}$

where
- $r_{s}$ and $r_{l}$ are the mean radii of small and large particles (`rad_det`; `rad_bdet`; [m])
- $\eta_{sw}$ is the dynamic viscosity of seawater at in situ temperature, salinity and pressure (`dynvis_sw(i,j,k)`, [kg m<sup>-1</sup> s<sup>-1</sup>])
- $\rho_{s}$ and $\rho_{l}$ are the harmonic mean densities of small and large particles (`rho_small`; `rho_large`, [kg m<sup>-3</sup>]) 
- $\rho_{sw}$ is the density of seawater, which we set here to a constant 1025 (kg m<sup>-1</sup>)

Our approach therefore considers mineral ballasting on particle excess density, particle size and the viscosity of fluid in determining sinking rates. This allows for "an environmentally dependent, space-varying $\omega_{s}$ and $\omega_{l}$" ([Dinauer et al., 2022](https://doi.org/10.1029/2021GB007131)).


---


### 23. Sedimentary processes.

WOMBAT-mid tracks the accumulation of organic detrital carbon (`p_det_sediment(i,j)`, $B_{det,sed}^{C}$, [mol C m<sup>-2</sup>]), organic detrital iron (`p_detfe_sediment(i,j)`, $B_{det,sed}^{Fe}$, [mol Fe m<sup>-2</sup>]), organic detrital silica (`p_detsi_sediment(i,j)`, $B_{det,sed}^{Si}$, [mol Si m<sup>-2</sup>]) and $CaCO_3$ (`p_caco3_sediment(i,j)`, $B_{CaCO_3,sed}^{C}$, [mol C m<sup>-2</sup>]) within sedimentary pools. The organic pools contribute to bottom fluxes of dissolved organic carbon (DOC), dissolved organic nitrogen (DON), dissolved inorganic carbon (DIC), dissolved iron (dFe), silicic acid ($Si(OH)_{4}$), oxygen ($O_2$) and alkalinity (Alk). 


**Organics**\
Remineralisation of organic carbon ($\gamma_{sed}^{&rarr; C}$) produces DOC and DON and removes $O_2$. Remineralisation of organic iron produces dFe and remineralisation of organic silica produces silicic acid. Ratios of nitrogen to carbon and oxygen to carbon are static at 16:122 and 132:122.

$\gamma_{sed}^{&rarr; DOC} = \gamma_{sed}^{0^{\circ}C} (β_{hete})^{T} B_{sed}^{C}$\
$\gamma_{sed}^{&rarr; DON} = \gamma_{sed}^{&rarr; DOC} R^{N:C}$\
$\gamma_{sed}^{&larr; O_2} = \gamma_{sed}^{&rarr; DOC} R^{O_2:C}$\
$\gamma_{sed}^{&rarr; dFe} = \gamma_{sed}^{0^{\circ}C} (β_{hete})^{T} B_{sed}^{Fe}$

where
- $\gamma_{sed}^{0^{\circ}C}$ is a base rate of organic matter hydrolysation at 0ºC in the sediments (`detlrem_sed`, [s<sup>-1</sup>])
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless])
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC])
- $B_{sed}^{C}$ is the concentration of organic carbon in the sediment pool (`p_det_sediment(i,j,1)`, [mol C m<sup>-2</sup>])
- $B_{sed}^{Fe}$ is the concentration of organic iron in the sediment pool (`p_detfe_sediment(i,j,1)`, [mol Fe m<sup>-2</sup>])
- $R^{N:C}$ is the static Redfield ratio of nitrogen to carbon in the organic matter (`16/122)` [mol N (mol C)<sup>-1</sup>])
- $R^{O_2:C}$ is the static Redfield ratio of dissolved oxygen to carbon required to hydrolyse organic matter (`132/122)` [mol $O_2$ (mol C)<sup>-1</sup>])


**Dissolution of biogenic silica**\
With regard to the dissolution of sedimentary biogenic silica, we compute it in the same way as how it is computed in the water column:

$\gamma_{sed}^{&rarr; Si} = \left(d_{sed^{Si}}^{T} \cdot S_{sed^{Si}}^{Sat} \cdot S_{sed^{Si}}^{bio}\right) B_{sed}^{Si}$

where
- $B_{sed}^{Si}$ is the concentration of organic silicon in the sediment pool (`p_detsi_sediment(i,j,1)`, [mol Si m<sup>-2</sup>])
- $d_{sed^{Si}}^{T}$ is the temperature-dependent rate of dissolution (`disssi_temp`, [s<sup>-1</sup>])
- $S_{sed^{Si}}^{Sat}$ is a scaling factor that decelerates dissolution as the in situ concentration approachs the equilibrium concentration (`disssi_usat`, [dimenionless])
- $S_{sed^{Si}}^{bio}$ is a scaling factor that accelerates dissolution in the presence of heterotrophic bacterial biomass (`disssi_bact`, [dimenionless])

Please refer to the description above for the equations that describe these terms.


**Dissolution of $CaCO_3$**\
Dissolution of $CaCO_3$ produces DIC and alkalinity. The sedimenary $CaCO_3$ pool is considered as entirely calcite. If `do_caco3_dynamics = .true.`, then sedimentary dissolution is controlled by bottom water temperature and an estimate of the pore-water calcite saturation state ($\Omega_{cal,sed}$):

$\D_{CaCO_{3},sed} = d_{CaCO_{3},sed} (β_{hete})^{T} \left(1 - \Omega_{cal,sed}\right)^{4.5}$

where
- $d_{CaCO_{3},sed}$ is a base rate of dissolution in units of [day<sup>-1</sup>], and
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless])
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC])
- $\Omega_{cal,sed}$ is the calcite saturation state within sedimentary pore waters (`sedomega_cal(i,j)`, [dimensionless])

The $\Omega_{cal,sed}$ is calculated using the `mocsy` package for solving carbonate chemistry of seawater ([Orr & Epitalon, 2015](https://doi.org/10.5194/gmd-8-485-2015)). These routines require Alk and DIC as inputs, along with nutrient concentrations and temperature and salinity of bottom waters. For DIC, we chose to sum the water column concentration of DIC and the organic carbon content of the sediment to approximate the interstitial (i.e., porewater) DIC concentration. We assume that the organic carbon content of the sediment (`p_det_sediment`), which is held in units of in [mol m<sup>-2</sup>] is relevant over 10 centimeters, and therefore can be automatically converted to [mol m<sup>-3</sup>] via division by 0.1. With this assumption these arrays can be added together and approximates the reducing conditions of organic-rich sediments, which have lower $\Omega_{cal,sed}$. This ensures a greater rate of $CaCO_3$ dissolution within the sediment as organic matter accumulates.

However, if `do_caco3_dynamics = .false.`, then dissolution of $CaCO_3$ in the sediments proceeds according to a constant assumed $\Omega_{cal,sed}$ of 0.2. We are aware that such a low $\Omega_{cal,sed}$ would not occur in real sediments given the buffering effect of dissolving $CaCO_3$ and the subsequent release of alkalinity. However, in the absence of feedbacks between organic carbon remineralisation and $CaCO_3$ dissolution, we assert a low $\Omega_{cal,sed}$ to ensure that sufficient $CaCO_3$ is dissolved back into the water column.


**Benthic denitrification**\
We also consider the consumption of $NO_3$ via benthic denitrification. When `do_benthic_denitrification = .true.`, a portion of the particulate organic matter within the sediments that is hydrolysed to DOC and DON is performed anaerobically (i.e., using $NO_3$ as the electron acceptor). Unlike this process in the water column, which is performed by bacterial metabolism, we estimate this process using an empirical parameterization from [Bohlen et al. (2012)](https://doi.org/10.1029/2011GB004198):

$\gamma_{sed}^{&larr; NO_3} = \gamma_{sed}^{&rarr; DOC} \min\left(0.9 \dfrac{94}{122}, \left(0.083 + 0.21 \cdot 0.98^{O_2 - NO_3} \right) \right)$

where
- $0.9$ is a hard upper limit stating that 90% of organic matter hydrolysation can potentially be performed anaerobically via denitrification
- $\dfrac{94}{122}$ is the stoichiometry of nitrate demand per mol of organic carbon hydrolysed ([Paulmier et al., 2009](https://doi.org/10.5194/bg-6-923-2009))
- $O_2$ is the bottom water concentration of dissolved oxygen (mmol m<sup>-3</sup>)
- $NO_3$ is the bottom water concentration of nitrate (mmol m<sup>-3</sup>)

and where the fraction of organic matter that is hydrolysed via denitrification is equal to:

$f_{sed}^{denit} = \gamma_{sed}^{&larr; NO_3} \dfrac{\dfrac{122}{94}}{\gamma_{sed}^{&rarr; DOC}}$


**Tendencies from sediment processes**\
Overall bottom fluxes of tracers are:

$\dfrac{\Delta DOC}{\Delta t} = \gamma_{det,sed}^{&rarr; DOC}$\
$\dfrac{\Delta DON}{\Delta t} = \gamma_{det,sed}^{&rarr; DON}$\
$\dfrac{\Delta NO_3}{\Delta t} = \gamma_{det,sed}^{&larr; NO_3}$\
$\dfrac{\Delta O_2}{\Delta t} = \gamma_{det,sed}^{&larr; O_2} \left(1 - f_{sed}^{denit}\right)$\
$\dfrac{\Delta Si}{\Delta t} = \gamma_{det,sed}^{&rarr; Si}$\
$\dfrac{\Delta dFe}{\Delta t} = \gamma_{det,sed}^{&rarr; dFe}$\
$\dfrac{\Delta DIC}{\Delta t} = \D_{CaCO_{3},sed}$\
$\dfrac{\Delta Alk}{\Delta t} = 2 \cdot \D_{CaCO_{3},sed}$

---

