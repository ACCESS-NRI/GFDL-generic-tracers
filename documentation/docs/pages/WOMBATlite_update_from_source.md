# Description of the WOMBATlite ocean biogeochemical model
## Subroutine - "update_from_source"

`!         (\___/)  .-.   .-. .--. .-..-..---.  .--. .-----.            !`\
`!         / o o \  : :.-.: :: ,. :: '' :: .; :: .; :'-. .-'            !`\
`!        (   "   ) : :: :: :: :: :: .. ::   .':    :  : :              !`\
`!         \__ __/  : '' '' ;: :; :: :; :: .; :: :: :  : :              !`\
`!                   '.,'.,' '.__.':_;:_;:___.':_;:_;  :_;              !`\
`!                                                                      !`\
`!  World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT)  !`

---

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

$PAR(k) = \sum_{b=1}^3 \dfrac{(PAR^{top}(k,b) - PAR^{top}(k+1,b))}{ex_{bgr}(k,b) * \Delta z(k)}$

This ensures phytoplankton growth in the model responds to the mean light they experience in the cell, not just light at one point. See Eq. 19 from [Baird et al. (2020)](https://gmd.copernicus.org/articles/13/4503/2020/).

The euphotic depth is defined as the depth where `radbio` falls below the 1% threshold of incidient shortwave radiation and is therefore soley dependent on the concentration of chlorophyll within the water column.

---
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

$L_{phy}^{N} = \dfrac{NO_3}{NO_3 + K_{phy}^{N}}$
 
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

$L_{phy}^{Fe} = \max(0.0, \min(1.0, \dfrac{ Q_{phy}^{Fe:C} - Q_{phy}^{-Fe:C} }{Q_{phy}^{*Fe:C}} )) $

If the cell is Fe‑replete with a quota that exceeds the minimum quota by as much as the optimal quota, then Fe does not limit growth ($L_{phy}^{Fe}$ = 1). If the cell is Fe‑deplete with a quota equal to the minimum quota, then the growth rate is reduced to zero. The optimal quota ($Q_{phy}^{*Fe:C}$) is therefore a measure of how much excess Fe is required to allow unrestricted growth.

---
### 3. Temperature‑dependent autotrophy and heterotrophy

**Autotrophy.**
The maximum potential growth rate for phytoplankton (`phy_mumax(i,j,k)`, $\mu_{phy}^{max}$, [day<sup>-1</sup>]) is prescribed by the temperature-dependent Eppley curve ([Eppley, 1972](https://spo.nmfs.noaa.gov/content/temperature-and-phytoplankton-growth-sea)). This formulation scales a reference growth rate at 0ºC via a power-law scaling with temperature (`Temp(i,j,k)`, $T$, [ºC]).

$\mu_{phy}^{max} = \mu_{phy}^{0^{\circ}C} \cdot β_{auto}^{(T)}$

In the above, both $\mu_{phy}^{0ºC}$ and $β_{auto}$ are reference values input to the model.

**Heterotrophy.**
Heterotrophic processes include mortality of phytoplankton and zooplankton, grazing rates of zooplankton and the remineralisation rate of detritus in the water column and sediments. These processes are scaled similarly to autotrophy, where some reference rate at 0ºC ($\mu_{het}^{0^{\circ}C}$, [day<sup>-1</sup>]) is multiplied by a power-law with temperature ($β_{hete}$).

$\mu_{het} = \mu_{het}^{0^{\circ}C} \cdot β_{hete}^{(T)}$

---
### 4. Light limitation of phytoplankton

Phytoplankton growth is limited by light through a photosynthesis–irradiance (P–I) relationship that links cellular chlorophyll content and available photosynthetically active radiation ('radbio', $PAR$, [W m<sup>-2</sup>]).

First, The initial slope of the P–I curve, (`phy_pisl`, $\alpha_{phy}$, [day<sup>-1</sup> (W m<sup>-2</sup>)<sup>-1</sup>]), determines how efficiently phytoplankton convert light into carbon fixation. It is scaled by the cellular chlorophyll-to-carbon ratio (`phy_chlc`, $Q_{phy}^{Chl:C}$, [mg/mg]).

$\alpha_{phy} = \max(\alpha_{phy}^{Chl} \cdot Q_{phy}^{Chl:C} \ , \ \alpha_{phy}^{Chl} \cdot Q_{phy}^{-Chl:C}$ )

where $\alpha_{phy}^{Chl}$ is the photosynthetic efficiency per unit chlorophyll in units of [day<sup>-1</sup> (W m<sup>-2</sup>)<sup>-1</sup> (mg/mg)<sup>-1</sup>] and $Q_{phy}^{-Chl:C}$ is the minimum chlorophyll to carbon ratio of the cell. This constraint prevents photosynthesis from collapsing unrealistically at low chlorophyll concentrations.

Second, light limitation (`phy_lpar(i,j,k)`, $L_{phy}^{PAR}$), [dimensionless]) is calculated using an exponential P–I formulation.

$L_{phy}^{PAR} = 1 - e^{- \alpha_{phy} PAR }$

At low irradiance (PAR), growth increases approximately linearly with light, while at high irradiance photosynthesis asymptotically saturates. Photoinhibition is not included in this formulation.

Realized growth of phytoplankton is then calculated as:

$\mu_{phy} = \mu_{phy}^{max} L_{phy}^{PAR} \min(L_{phy}^{N}, L_{phy}^{Fe})$

---
### 5. Growth of chlorophyll

This step diagnoses the **rate of chlorophyll production** as a function of mixed-layer light, the phytoplankton growth rate, and nutrient availability. The structure is consistent with common **photoacclimation / variable chlorophyll-to-carbon** approaches that accelerate or decelerate chlorophyll growth relative to carbon growth.

A chlorophyll-specific initial P-I slope (`pchl_pisl`, $\alpha_{chl}$, [day<sup>-1</sup> (W m<sup>-2</sup>)<sup>-1</sup>]) is constructed by scaling the photosynthetic carbon-specific P-I slope ($\alpha_{phy}$):

$\alpha_{chl} = \dfrac{ \alpha_{phy} }{ \mu_{phy}^{max} (1 - \min(L_{phy}^{N}, L_{phy}^{Fe})) + \varepsilon }$

where $\varepsilon$ is a small constant preventing division by zero. Above, warm temperatures that elevate $\mu_{phy}^{max}$ and nutrient-stress decrease the chlorophyll-specific P-I slope, while nutrient-rich conditions steepen it. There is strong evidence that both nitrogen and iron limitation limit chlorophyll-to-carbon ratios of phytoplankton cells, while warm waters make enzymatic carbon fixation more efficient and lower light harvesting demand.

Then, light limitation (`pchl_lpar(i,j,k)`, $L_{phy}^{PAR}$), [dimensionless]) is calculated using an exponential P–I formulation of the same form as carbon-specific light limitation.

$L_{chl}^{PAR} = 1 - e^{- \alpha_{chl} PAR_{MLD} }$

where:

- $PAR_{MLD} =$ `radmld` is the average downwelling irradience in the mixed layer and ensures that chlorophyll growth is amplified relative to carbon growth in the lower part of the mixed layer, but suppressed in the upper part.

Before solving the growth rate of chlorophyll, we determine the minimum and optimal rates of chlorophyll production in mg m<sup>-3</sup> day<sup>-1</sup>, which are realised from the already known phytoplankton carbon growth rate:

$\mu_{chl}^{min} = Q_{phy}^{-Chl:C} \mu_{phy} B_{phy}^{C} 12$

$\mu_{chl}^{opt} = Q_{phy}^{*Chl:C} \mu_{phy} B_{phy}^{C} 12$

where:

- $Q_{phy}^{-Chl:C} =$ `phyminqc` is the minimum chlorophyll quota,
- $Q_{phy}^{*Chl:C} =$ `phyoptqc` is the optimal chlorophyll quota,
- $\mu_{phy} =$ `phy_mu` is the realised phytoplankton growth rate,
- $B_{phy}^{C} 12$ is phytoplankton carbon concentration in mg m<sup>-3</sup> (implemented as `biophy * 12.0`).


Chlorophyll production is then computed as

$\mu_{chl}^{\delta} = (\mu_{chl}^{opt} - \mu_{chl}^{min}) \cdot L_{chl}^{PAR} \min(L_{phy}^{N}, L_{phy}^{Fe})$

$\mu_{chl} = \phi ( \mu_{chl}^{min} + \dfrac{ \mu_{chl}^{\delta} }{ \alpha_{phy} PAR_{MLD} } )$

where:

- $\phi = \dfrac{PAR_{MLD}}{PAR_{MLD} + K_{chl}^{PAR}}$ and reduces growth during extended periods of low light, such as the polar winter;
- $K_{chl}^{PAR}$ is a half-saturation light scale for chlorophyll growth (`chlkWm2`, [W m<sup>-2<\sup>]).
- $\mu_{chl}^{\delta}$ is an additional growth of chlorophyll above what is required to maintain the minimum possible quota.

---

### 6. Phytoplankton uptake of iron

Like chlorophyll, the iron content of phytoplankton is explicitly tracked as a tracer in WOMBAT-lite. First, a maximum quota is found dependent on the maximum quota of Fe within the phytoplankton type (`phymaxqf`, $Q_{phy}^{^Chl:C}$, mol/mol) and the phytoplankton concentration in the water column:

$B_{phy}^{"Fe} = B_{phy}^{C} Q_{phy}^{"Chl:C}$

Following [Aumont et al. (2015)](https://gmd.copernicus.org/articles/8/2465/2015/), this rate is scaled by three terms relating to (i) michaelis-menten type affinity for dFe, (ii) up-regulation of dFe uptake representing investment in transporters when cell quotas are limiting to growth, and (iii) down regulation of dFe uptake associated with enriched cellular quotas:

(i) $\dfrac{dFe}{dFe + K_{phy}^{Fe}}$

(ii) $4 - \dfrac{4.5 L_{phy}^{Fe}}{0.5 + L_{phy}^{Fe}}$

(iii) $\max\left(0, 1 - \dfrac{B_{phy}^{Fe} / B_{phy}^{"Fe}}{|1.05 - B_{phy}^{Fe} / B_{phy}^{"Fe}|} \right)$

dFe uptake by phytoplankton is then calculated as

$\mu_{phy}^{Fe} = \mu_{phy}^{max} B_{phy}^{"Fe} \max\left(0.2, L_{phy}^{PAR} \cdot L_{phy}^{N}\right) \cdot (i) \cdot (ii) \cdot (iii)$

where the maximum dFe uptake is decreased to 20% of its maximum potential rate at night and when nitrogen is limited. The iron to carbon ratios of phytoplankton are passed to zooplankton and detritus and are also tracked in these pools.

---
### 7. Iron chemistry (precipitation, scavenging and coagulation)

Treatment of dissolved iron (dFe, µmol m-3) follows [Aumont et al. (2015)](https://gmd.copernicus.org/articles/8/2465/2015/) and [Tagliabue et al. (2023)](https://www.nature.com/articles/s41586-023-06210-5). 

We first estimate the **solubility of free Fe from Fe<sup>3+</sup>** in solution using temperature, pH and salinity using the thermodynamic equilibrium equations of [Liu & Millero (2002)](https://www.sciencedirect.com/science/article/abs/pii/S0304420301000743). 

$T_K = \max(5.0, T) + 273.15$ \
$T_K^{-1} = \dfrac{1}{T_K}$ \
$I_{S} = \dfrac{19.924,S}{1000 - 1.005,S}$

Solubility constants:

$Fe_{sol1} = 10^{\left(-13.486 - 0.1856\sqrt{I_S} + 0.3073 I_S + 5254,T_K^{-1}\right)}$ \
$Fe_{sol2} = 10^{\left(2.517 - 0.8885\sqrt{I_S} + 0.2139 I_S - 1320,T_K^{-1}\right)}$ \
$Fe_{sol3} = 10^{\left(0.4511 - 0.3305\sqrt{I_S} - 1996,T_K^{-1}\right)}$ \
$Fe_{sol4} = 10^{\left(-0.2965 - 0.7881\sqrt{I_S} - 4086,T_K^{-1}\right)}$ \
$Fe_{sol5} = 10^{\left(4.4466 - 0.8505\sqrt{I_S} - 7980,T_K^{-1}\right)}$

Final Fe(III) solubility:

$Fe_{sol} = Fe_{sol1}\left([H^+]^3 + Fe_{sol2}[H^+]^2 + Fe_{sol3}[H^+] + Fe_{sol4} + \dfrac{Fe_{sol5}}{[H^+]}\right)\times10^{9}$

where $[H^+]$ is in units of mol/L and $Fe_{sol}$ is in units of mmol/m<sup>3</sup>.

Next we **estimate the concentration of colloidal iron** in solution following [Tagliabue et al. 2023](https://www.nature.com/articles/s41586-023-06210-5). Colloidal Fe (`fecol(i,j,k)`, $Fe_{col}$, [mmol Fe m<sup>-3</sup>]) is whatever exceeds the inorganic solubility ceiling (`fe3sol`, $Fe_{sol}$, [mmol Fe m<sup>-3</sup>]), but we enforce a hard minimum that colloids are at least 10% of total dissolved Fe (`biofer`, $dFe$, [mmol Fe m<sup>-3</sup>]).

$Fe_{col} = \max\left(0.1 dFe\, \ dFe - Fe_{sol} \right)$

Following solving for colloidal Fe, we **partition the remaining dissolved Fe into ligand-bound and free iron**. To do so, we find the remain dissolved iron not in colloidal form (`fe_sfe`, $Fe_{sFe}$, [mmol Fe m<sup>-3</sup>]), 

$Fe_{sFe} = \max\left(0.0,\ dFe - Fe_{col} \right)$

Partitioning is done using a standard quadratic form that drops out of mass balance + equilibrium for 1:1 complexation with a single ligand class. To do so, we need the temperature-dependent conditional stability constant for Fe–ligand complexation ($Fe_{Keq}$), 
 
$Fe_{Keq} = 10^{ \left(17.27 - 1565.7 T_K^{-1} \right) }\times 10^{-9}$, 

and solving for free Fe (`feIII`, $Fe_{free}$, [mmol Fe m<sup>-3</sup>]):

$z = 1.0 + [Ligand] \cdot Fe_{Keq} - Fe_{sFe}\cdot Fe_{Keq}$ \
$Fe_{free} = \dfrac{-z + \sqrt{z^2 + 4.0 Fe_{Keq} Fe_{sFe}}}{2 Fe_{Keq} + \varepsilon}$
$Fe_{free} = \max\left(0,\ \min(Fe_{free},\ Fe_{sFe})\right)$

Whatever soluble Fe is not present as inorganic Fe' is assigned to ligand-bound Fe:

$Fe_{lig} = Fe_{sFe} - Fe_{free}$

Now that we have separated the dissolved Fe pool into its subcomponents of free, ligand-bound and colloidal Fe, we **solve for precipiation of nanoparticles, scavenging and coagulation** of dissolved Fe, all of which remove dFe from the water column.

Precipiation:
$Fe_{nanop}^{&rarr;} = \max\left(0.0,\ Fe_{free} - Fe_{sFe}\right)\, \gamma_{Fe}^{nano}$

Scavenging:
$Fe_{scav}^{&rarr;} = Fe_{free},\ (10^{-7} + \gamma_{Fe}^{scav} \cdot (B^{C}_{det} + B^{C}_{CaCO_3}))$ \
$Fe_{scav}^{&rarr;Det} = Fe_{scav}^{&rarr;} \cdot \dfrac{ B_{det}^{C} }{ B_{det}^{C} + B_{CaCO_3}^{C} }$

Coagulation:


---
### 8. Mortality scalings and grazing

---
### 9. CaCO3 calculations

---
### 10. Sources and Sinks

---
### 11. Tracer tendencies

---
### 12. Check for conservation of mass
