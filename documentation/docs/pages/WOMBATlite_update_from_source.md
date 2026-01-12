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
at each tracer time‑step. The subroutine is documented internally by a list of numbered steps (see code comments). These steps are:

1. Light attenuation through the water column.
2. Nutrient limitation of phytoplankton.
3. Temperature‑dependent autotrophy and heterotrophy.
4. Light limitation of phytoplankton.
5. Growth of chlorophyll.
6. Phytoplankton uptake of iron.
7. Iron chemistry (precipitation, scavenging and coagulation).
8. Mortality and remineralisation
9. Zooplankton grazing.
10. CaCO3 calculations.
11. Tracer tendencies (update tracer concentrations).
12. Check for conservation of mass.

Below is a step‑by‑step explanation of each section together with the key equations. Variable names in grey follow the Fortran code, while 
variable names in math font are pointers to the equations; i,j,k refer to horizontal and vertical indices; square brackets denote units. If a variable is without i,j,k dimensions, this variable is held as a scalar and not an array.

---


### 1. Light attenuation through the water column

Photosynthetically available radiation (PAR) is split into blue, green and red wavelengths. Surface incoming shortwave radiation is multiplied by 0.43 to return the photosynthetically active radiation (PAR, [W m<sup>-2</sup>) of total shortwave radiation, and is then split evenly into each of blue, green and red light bands. 

At the top (`par_bgr_top(k,b)`, $PAR^{top}$) and mid‑point (`par_bgr_mid(k,b)`, $PAR^{mid}$) of each layer `k` we calculate the downward irradiance by exponential decay of each band `b` through the layer thickness (`dzt(i,j,k)`, $\Delta z$, [m]) using band‑specific attenuation coefficients (`ek_bgr(k,b)`, $ex_{bgr}$, [m<sup>-1</sup>]). These attenuation coefficients are related to the concentration of chlorophyll ([mg m<sup>-3</sup>]), organic detritus ([mg N m<sup>-3</sup>]) and calcium carbonate ([kg m<sup>-3</sup>]) in the water column. For chlorophyll, attentuation coefficients for each of blue, green and red light are retrieved from the look-up table of [Morel & Maritorena (2001)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2000JC000319) (their Table 2) that explicitly relates chlorophyll concentration to attenutation rates and accounts for the packaging effect of chlorophyll in larger cells. For organic detritus, attenutation coefficients for blue, green and red light are taken from [Dutkiewicz et al. (2015)](https://bg.copernicus.org/articles/12/4447/2015/bg-12-4447-2015.html) (their Fig. 1b), while for calcium carbonate we take the coefficients defined in [Soja-Wozniak et al. (2019)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JC014998). For both detritus and calcium carbonate, these studies provide concentration normalized attenutation coefficients, which must be multipled against concentrations to retrieve the correct units of [m<sup>-1</sup>].

As an example, the PAR in the blue band (`b=1`) at the top of level k is computed as

$PAR^{top}(k,1) = PAR^{top}(k-1,1) * e^{(-ex_{bgr}(k-1,1) * \Delta z(k-1))}$

where the total attenutation rate of blue light in the grid cell above `k` is the sum of attenuation due to all particulates in that grid cell:

$ex_{bgr}(k-1,1) = ex_{chl}(k-1,1) + ex_{det}(k-1,1) + ex_{CaCO_3}(k-1,1)$

The irradiance in the red band (`b=3`) at the mid point of layer `k`, in contrast, is equal to 

$PAR^{mid}(k,3) = PAR^{mid}(k-1,3) * e^{(-0.5*(ex_{bgr}(k-1,3) * \Delta z(k-1) + ex_{bgr}(k,3) * \Delta z(k)))}$

The total PAR available to phytoplantkon is assumed to be the sum of the blue, green and red bands. Because we assume that phytoplankton are 
homogenously distributed within a layer `k`, but we do not assume that light is homogenously distributed within that layer, we solve for the 
PAR that is seen by the average phytoplankton within that cell ('radbio', $PAR$, [W m<sup>-2</sup>]), where

$PAR(k) = \sum_{b=1}^3 \dfrac{(PAR^{top}(k,b) - PAR^{top}(k+1,b))}{ex_{bgr}(k,b) * \Delta z(k)}$

This ensures phytoplankton growth in the model responds to the mean light they experience in the cell, not just light at one point. See Eq. 19 from [Baird et al. (2020)](https://gmd.copernicus.org/articles/13/4503/2020/).

The euphotic depth is defined as the depth where `radbio` falls below the 1% threshold of incidient shortwave radiation and is therefore soley dependent on the concentration of chlorophyll within the water column.

---


### 2. Nutrient limitation of phytoplankton

At the start of each vertical loop the code computes biomass of phytoplankton (`biophy`, $B_{phy}$, [mmol C m<sup>-3</sup>]). Phytoplankton biomass is used to scale how both nitrate (`biono3`, $NO_{3}$, [mmol NO<sub>3</sub> m<sup>-3</sup>]) and dissolved iron (`biofer`, $dFe$, [µmol dFe m<sup>-3</sup>]) affect the growth of phytoplankton. Using compilations of marine phytoplankton and zooplankton communities, [Wickman et al. (2024)](https://www.science.org/doi/10.1126/science.adk6901) show that the nutrient affinity, $aff$, of a phytoplankton cell is related to its volume, $V$, via

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

$\mu_{phy}^{max} = \mu_{phy}^{0^{\circ}C} \cdot (β_{auto})^{T}$

In the above, both $\mu_{phy}^{0ºC}$ and $β_{auto}$ are reference values input to the model and control how productive the ocean is.

**Heterotrophy.**
Heterotrophic processes include mortality of phytoplankton and zooplankton, grazing rates of zooplankton and the remineralisation rate of detritus in the water column and sediments. These processes are scaled similarly to autotrophy, where some reference rate at 0ºC ($\mu_{het}^{0^{\circ}C}$, [day<sup>-1</sup>]) is multiplied by a power-law with temperature ($β_{hete}$).

$\mu_{het} = \mu_{het}^{0^{\circ}C} \cdot (β_{hete})^{T}$

See sections below for further details on mortality and grazing terms.

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

and 

$\mu_{phy}^{C} = \mu_{phy} B_{phy}^{C}$

where $\mu_{phy}^{C}$ is in units of [mmol C m<sup>-3</sup> day<sup>-1</sup>].


---


### 5. Growth of chlorophyll

This step diagnoses the **rate of chlorophyll production** as a function of mixed-layer light, the phytoplankton growth rate and nutrient availability. The structure is consistent with common **photoacclimation / variable chlorophyll-to-carbon** approaches that accelerate or decelerate chlorophyll growth relative to carbon growth.

A chlorophyll-specific initial P-I slope (`pchl_pisl`, $\alpha_{chl}$, [day<sup>-1</sup> (W m<sup>-2</sup>)<sup>-1</sup>]) is constructed by scaling the photosynthetic carbon-specific P-I slope ($\alpha_{phy}$) by particular environmental conditions:

$\alpha_{chl} = \dfrac{ \alpha_{phy} }{ \mu_{phy}^{max} e^{\left(-\min\left(L_{phy}^{N}, L_{phy}^{Fe} \right)\right) } }$

Above, warm temperatures that elevate $\mu_{phy}^{max}$ decrease the chlorophyll-specific P-I slope, and nutrient-stress also decreases this slope. Meanwhile, cold and nutrient-rich conditions steepen it. Such conditions are found in at the bottom of euphotic zones. There is strong evidence that both nitrogen and iron limitation limit chlorophyll-to-carbon ratios of phytoplankton cells [](), while warm waters make enzymatic carbon fixation more efficient and lower light harvesting demand, thereby decreasing Chl:C ratios []().

After solving for the chlorophyll-specific P-I slope, we calculate light limitation (`pchl_lpar(i,j,k)`, $L_{phy}^{PAR}$), [dimensionless]) using an exponential P–I formulation of the same form as carbon-specific light limitation,

$L_{chl}^{PAR} = 1 - e^{- \alpha_{chl} PAR_{MLD} }$

but here we use $PAR_{MLD}$ rather than $PAR$. In this case, $PAR_{MLD}$ (`radmld`, [W m<sup>-2</sup>]) is the average downwelling irradience in the mixed layer. This ensures that chlorophyll growth is amplified relative to carbon growth in the lower part of the mixed layer, but suppressed in the upper part. As such, chlorophyll growth is able to be accelerated relative to carbon growth by both (1) scaling $\alpha_{chl}$ and (2) using $PAR_{MLD}$ rather than $PAR$. For waters beneath the mixed layer, $PAR_{MLD}$ = $PAR$, so for these waters chlorophyll to carbon ratios are only affected by scaling $\alpha_{chl}$.

Before solving the growth rate of chlorophyll, we determine the minimum and maximum rates of chlorophyll production in mg m<sup>-3</sup> day<sup>-1</sup>, which are realised from the already known phytoplankton carbon growth rate:

$\mu_{phy}^{-chl} = Q_{phy}^{-Chl:C} \mu_{phy}^{C} \cdot 12$

$\mu_{phy}^{+chl} = Q_{phy}^{+Chl:C} \mu_{phy}^{C} \cdot 12$

where:

- $Q_{phy}^{-Chl:C} =$ `phyminqc` is the minimum chlorophyll quota,
- $Q_{phy}^{+Chl:C} =$ `phymaxqc` is the maximum chlorophyll quota,
- $\mu_{phy}^{C} =$ is the realised phytoplankton growth rate in units of [mmol C m<sup>-3</sup> day<sup>-1</sup>].

These are the rates of chlorophyll growth required to support phytoplankton growth at both a minimum and maximum cellular quota. Chlorophyll-specific light limitation is then used to determine the extent to which optimal chlorophyll growth can be achieved:

$\delta_{phy}^{chl} = \phi \left(\mu_{phy}^{+chl} - \mu_{phy}^{-chl}\right) \cdot L_{chl}^{PAR}$

$\mu_{phy}^{chl} = \mu_{phy}^{-chl} + \delta_{phy}^{chl} $

where:

- $\phi = \dfrac{PAR_{MLD}(k=1)}{PAR_{MLD}(k=1) + K_{chl}^{PAR}}$ and reduces growth during extended periods of low light, such as the polar winter;
- $K_{chl}^{PAR}$ is a half-saturation light scale for chlorophyll growth (`chlkWm2`, [W m<sup>-2</sup>]).
- $\delta_{phy}^{chl}$ is the additional growth of chlorophyll above what is required to maintain the minimum quota.

---


### 6. Phytoplankton uptake of iron

Like chlorophyll, the iron content of phytoplankton is explicitly tracked as a tracer in WOMBAT-lite. First, a maximum quota is found dependent on the maximum quota of Fe within the phytoplankton type (`phymaxqf`, $Q_{phy}^{Chl:C}$, [mol/mol]) and the phytoplankton concentration in the water column:

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

Now that we have separated the dissolved Fe pool into its subcomponents of free, ligand-bound and colloidal Fe, we solve for precipiation of nanoparticles, scavenging and coagulation of dissolved Fe, all of which remove dFe from the water column. These are the major sinks outside of phytoplankton uptake.

**Precipiation:** \
$Fe_{nanop}^{&rarr;} = \max\left(0.0,\ Fe_{free} - Fe_{sFe}\right)\, \gamma_{Fe}^{nano}$

**Scavenging:** \
$Fe_{scav}^{&rarr;} = Fe_{free} \left(10^{-7} + \gamma_{Fe}^{scav} \cdot (B_{det}^{C} + B_{CaCO_3}^{C}) \right)$ \
$Fe_{scav}^{&rarr;det} = Fe_{scav}^{&rarr;} \cdot \dfrac{ B_{det}^{C} }{ B_{det}^{C} + B_{CaCO_3}^{C} }$

**Coagulation:** \
$Fe_{coag}^{&rarr;det} = Fe_{col} \Gamma_{Fe}^{coag}$

When the grid cell is within the mixed layer:
- $\Gamma_{Fe}^{coag} = \left(\ \ \ \ \ \ \ \ \  (12 \cdot F_{coag} [DOC] + 9 \cdot B_{det}^{C}) + 2.5 \cdot B_{det}^{C} + 128 \cdot F_{coag} [DOC] + 725 \cdot B_{det}^{C} \right) \gamma_{Fe}^{coag}$

When the grid cell is beneath the mixed layer:
- $\Gamma_{Fe}^{coag} = \left( 0.01 \cdot (12 \cdot F_{coag} [DOC] + 9 \cdot B_{det}^{C}) + 2.5 \cdot B_{det}^{C} + 128 \cdot F_{coag} [DOC] + 725 \cdot B_{det}^{C} \right) \gamma_{Fe}^{coag}$

And where:
- $F_{coag}$ is the phytoplankton-dependent coagulation factor <br> $F_{coag}= \max\left(\dfrac{1}{3},\ \dfrac{B_{phy}^{C}}{B_{phy}^{C} + 0.03}\right)$
- $[DOC]$ is a proxy estimate of the concentration of dissolved organic carbon <br> $[DOC] = 10 + 40 \cdot\left(1 - \min(L_{phy}^{N},\ L_{phy}^{Fe})\right) $
- $\gamma_{Fe}^{coag}$ is the iron coagulation rate constant
- $B_{phy}^{C}$ is the concentration of phytoplankton carbon biomass
- $B_{det}^{C}$ is the concentration of detrital carbon biomass

Together, these terms implement a biologically mediated coagulation pathway in which iron removal from the dissolved pool is tightly coupled to ecosystem state. The formulation reflects the central conclusion of [Tagliabue et al. (2023)](https://www.nature.com/articles/s41586-023-06210-5): that iron cycling is not governed solely by inorganic chemistry, but is strongly regulated by biological activity, organic matter dynamics, and particle ecology across the upper ocean.

---


### 8. Mortality and remineralisation 

Mortality of phytoplankton and zooplankton are affected by both linear ($\gamma$) and quadratic ($\Gamma$) terms. Linear terms are per-capita losses associated with the costs of basal metabolism. Quadratic, and thus density-dependent losses, are associated with disease, aggregation and coagulation, viruses, infection and canabalism. None of these processes are represented explicitly within the model, so we represent them implicitly.

**Linear losses** for phytoplankton and zooplankton in [mmol m<sup>-3</sup> s<sup>-1</sup>] are modelled as

$\gamma_{phy}^{C} = \gamma_{phy}^{0^{\circ}C} (β_{hete})^{T} B_{phy}^{C}$

$\gamma_{zoo}^{C} = \gamma_{zoo}^{0^{\circ}C} (β_{hete})^{T} F_{zoo}^{\gamma} B_{zoo}^{C}$

In the above, we scale down **linear zooplankton mortality** when zooplankton biomass is small, such that

$F_{zoo}^{\gamma} = \dfrac{B_{zoo}^{C}}{B_{zoo}^{C} + K_{zoo}^{\gamma}}$,

where
- $B_{zoo}^{C}$ is the concentration of zooplankton carbon biomass
- $K_{zoo}^{\gamma}$ is the half-saturation coefficient for scaling down linear mortality losses


**Quadratic losses** of phytoplankton and zooplankton in [mmol m<sup>-3</sup> s<sup>-1</sup>] are modelled as

$\Gamma_{phy}^{C} = \Gamma_{phy}^{0^{\circ}C} (β_{hete})^{T} \left(B_{phy}^{C}\right)^{2}$

$\Gamma_{zoo}^{C} = \Gamma_{zoo}^{0^{\circ}C} (β_{hete})^{T} \left(B_{zoo}^{C}\right)^{2}$

where
- $\Gamma^{0^{\circ}C}$ is the reference rate of biomass loss in units of [(mmol C m<sup>-3</sup>)<sup>-1</sup> day<sup>-1</sup>].


**Remineralisation** of detritus is only affected by a quadratic, density-dependent loss term,

$\Gamma_{det}^{C} = \Gamma_{det}^{0^{\circ}C} (β_{hete})^{T} \left(B_{det}^{C}\right)^{2}$

since hydrolyzation of organic detritus is performed by an heterotrophic bacterial population that is not explicitly resolved in the model.


---


### 9. Zooplankton grazing, growth, excretion and egestion

**Grazing by zooplankton** (`g_npz`, $g_{zoo}$, [day<sup>-1</sup>]) is computed using a Holling Type III functional response [Holling, 1959](https://doi.org/10.4039/Ent91385-7), where:

$g_{zoo} = \dfrac{\mu_{zoo}^{max} (β_{hete})^{T} \varepsilon (B_{prey}^{C})^{2}}{\mu_{zoo}^{max} (β_{hete})^{T} + \varepsilon (B_{prey}^{C})^{2}}$

This formulation suppresses grazing at very low prey biomass ($B_{prey}^{C}$) due to reduced encounter and clearance rates, accelerates grazing at intermediate prey biomass as zooplankton effectively learn and switch to available prey, and saturates at high prey biomass due to handling-time limitation ([Gentleman and Neuheimer, 2008](https://doi.org/10.1093/plankt/fbn078); Rohr et al., [2022](https://doi.org/10.1016/j.pocean.2022.102878), [2024](https://doi.org/10.1029/2023GL107732)). This choice increases ecosystem stability and prolongs phytoplankton blooms relative to a Type II formulation.

The application of $g_{zoo}^{max} (β_{hete})^{T}$ in both the numerator and denominator makes this grazing formula unique [(Rohr et al., 2023)](https://www.nature.com/articles/s43247-023-00871-w) and equivalent to a disk formulation, rather than a Michaelis–Menten formulation [(Rohr et al., 2022)](https://doi.org/10.1016/j.pocean.2022.102878). Practically, this amplifies grazing in warmer climes, but to a lesser extent than other formulations that apply the temperature amplification ($(β_{hete})^{T}$) only in the numerator [(Rohr et al., 2023)](https://www.nature.com/articles/s43247-023-00871-w). This dampens the effect that variations in temperature have on grazing activity, amplifying the effect of $\varepsilon and aligning with observations that the ratio of grazing to phytoplankton growth varies little between tropical and polar climes [(Calbet and Landry, 2004)](https://doi.org/10.4319/lo.2004.49.1.0051). Theoretically, this assumes some evolutionary adaptation to account for the physiological effects of temperature across environmental niches, such that the efficiency of prey capture and handling becomes more important to grazers than metabolic constraints due to temperature.

The total prey biomass available to zooplankton is defined as a preference-weighted sum of phytoplankton and detritus:

$B_{prey}^{C} = \phi_{zoo}^{phy} B_{phy}^{C} + \phi_{zoo}^{det} B_{det}^{C}$

where $B_{phy}^{C}$ and $B_{det}^{C}$ are phytoplankton and detrital carbon biomass, respectively, and $\phi_{zoo}$ terms define relative grazing preferences for these prey items.

The prey capture rate coefficient, $\varepsilon$ (`zooeps(i,j,k)`, $\varepsilon$, [(mmol C m<sup>-3</sup>)<sup>-1</sup>]),  is allowed to vary as a function of prey biomass, following the prey-dependent behaviour described by [Rohr et al. (2024)](doi.org/10.1029/2023GL107732). This reflects a transition from microzooplankton-like feeding with higher prey capture rate coefficients at low prey biomass to mesozooplankton-like feeding with lower prey capture rate coefficients at high prey biomass.

A prey-dependent scaling factor is defined as

$F_{prey} = e^{\left(-B_{prey}^{C} \varepsilon_{shift} \right)}$

and the effective capture rate coefficient, $\varepsilon$ is then computed as

$\varepsilon = \varepsilon_{\min} + (\varepsilon_{\max} - \varepsilon_{\min}) F_{prey}$

At low prey biomass, $\varepsilon \rightarrow \varepsilon_{\max}$, enhancing grazing efficiency. At high prey biomass, $\varepsilon \rightarrow \varepsilon_{\min}$, reducing capture efficiency as handling time and feeding mode are more ineffective on average in a community with relatively more mesozooplankton.

Zooplankton total grazing of biomass ([mmol C m<sup>-3</sup> day<sup>-1</sup>]) is therefore

$g_{zoo}^{C} = g_{zoo} B_{zoo}^{C} = \left(g_{zoo}^{&rarr; phy} + g_{zoo}^{&rarr; det}\right) B_{zoo}^{C} $

where:
- $g_{zoo}^{&rarr; phy} = g_{zoo} \dfrac{\phi_{zoo}^{phy} B_{phy}^{C}}{B_{prey}^{C}}$ and is the proportion of zooplankton grazing of phytoplankton
- $g_{zoo}^{&rarr; det} = g_{zoo} \dfrac{\phi_{zoo}^{det} B_{det}^{C}}{B_{prey}^{C}}$ and is the proportion of zooplankton grazing of detritus


**Zooplankton growth**, or biomass accumulation ($\mu_{zoo}^{C}$), is then calculated assuming a static assimilation coefficient ($\lambda$). $\lambda$ is the fraction of ingested material that is assimilated into zooplankton biomass and thus contributes to their growth.

$\mu_{zoo}^{C} = g_{zoo}^{C} \lambda$


**Zooplankton excretion and egestion** are then calculated from the unassimilated material based on an excretion coefficient ($\eta$). $\eta$ determines the fraction of unassimilated material that is routed towards excretion ($\eta_{zoo}^{C}$), rather than egestion ($E_{zoo}^{C}$).

$\eta_{zoo}^{C} = g_{zoo}^{C} \left(1 - \lambda\right)\eta$

$E_{zoo}^{C} = g_{zoo}^{C} \left(1 - \lambda\right) \left(1-\eta\right)$


---


### 10. CaCO3 calculations

**Dynamic CaCO$_3$ production and dissolution**

When $CaCO_3$ dynamics are enabled (`do_caco3_dynamics = .true.`), the model computes both particulate inorganic carbon production (via the PIC:POC ratio) and $CaCO_3$ dissolution rates as functions of carbonate chemistry, temperature, and organic matter availability.

The particulate inorganic to organic carbon production ratio (`pic2poc`, $PIC:POC$, [mol/mol]) is formulated as a function of the substrate–inhibitor ratio between bicarbonate and free hydrogen ions (`hco3 / htotal(i,j,k)`, $\dfrac{[HCO_3^-]}{[H^+]}$, [mol/µmol]), following [Lehmann & Bach (2025)](https://www.nature.com/articles/s41561-025-01644-0). This reflects the sensitivity of calcification to carbonate system speciation.

Bicarbonate concentration is diagnosed as $[HCO_3^-] = DIC - [CO_3^{2-}] - [CO_2^*]$ where $DIC$ is total dissolved inorganic carbon and $CO_2^*$ is dissolved $CO_2$. The $PIC:POC$ ratio is then computed as

$PIC:POC = \min \left( 0.3,  \left( f_{\text{inorg}} + 10^{-3 + 4.31 \times 10^{-6} \left( \dfrac{[HCO_3^-]}{[H^+]} \right)} \right) F_T \right)$

where $f_{\text{inorg}}$ is a low background inorganic fraction and $[H^+]$ is total free hydrogen ion concentration. The exponential term captures the nonlinear enhancement of calcification with increasing $[HCO_3^-]/[H^+]$. Calcification is further modulated by a temperature-dependent suppression term,

$F_T = 0.55 + 0.45 \tanh (T - 4)$

which strongly reduces $CaCO_3$ production in cold waters. This formulation enforces near-zero calcification below approximately 3 °C, consistent with observations of _Emiliania huxleyi_ growth limits in polar environments [(Fielding, 2013)](https://doi.org/10.4319/lo.2013.58.2.0663). The $PIC:POC$ ratio is also capped at an upper bound of 0.3 to prevent unrealistically high inorganic carbon production and accord with the highest measured ratios in the ocean.


Dissolution of $CaCO_3$ (`dissrat`) is computed as the sum of three contributions: undersaturation-driven dissolution of (i) calcite and (ii) aragonite, and (iii) biologically mediated dissolution associated with degredation of detrital organic matter. This formulation follows [Kwon et al. (2024)](https://www.science.org/doi/full/10.1126/sciadv.adl0779).

Calcite dissolution is given by

$D_{\text{cal}} = d_{\text{cal}} \max(0,  1 - \Omega_{\text{cal}})^{2.2}$

and aragonite dissolution by

$D_{\text{ara}} = d_{\text{ara}} \max(0,; 1 - \Omega_{\text{ara}})^{1.5}$

where $\Omega_{\text{cal}}$ and $\Omega_{\text{ara}}$ are the saturation states of calcite and aragonite, respectively, and $d_{\text{cal}}$ and $d_{\text{ara}}$ are reference dissolution rate constants in units of [day<sup>-1</sup>]. Dissolution is activated only under undersaturated conditions ($\Omega < 1$) and increases nonlinearly with increasing undersaturation.

An additional detritus-associated dissolution term is included (`diss_det`, $D_{\text{det}}$, [day<sup>-1</sup>]):

$D_{\text{det}} = d_{\text{det}} \Gamma_{det}$

where $\Gamma_{det}$ is the local remineralisation rate of organic matter ([mmol m<sup>-3</sup> day<sup>-1</sup>]. This term represents shallow water $CaCO_3$ dissolution due to reducing microenvironments. In this scenario, $\Omega$ of calcite and aragonite in the water column tend to be > 1 [(Sulpis et al. (2021)](https://www.nature.com/articles/s41561-021-00743-y) but dissolution nonetheless occurs in microenvironments enriched in $CO_2^*$.

The rate of $CaCO_3$ dissolution is then calculated by summing these three dissolution terms and applying them to the biomass of $CaCO_3$ in the water column ($B_{CaCO_3}$): 

$D_{\text{CaCO}3} = \left( D{\text{cal}} + D_{\text{ara}} + D_{\text{det}} \right) B_{CaCO_3}$ 


**Static CaCO$_3$ production and dissolution**

When $CaCO_3$ dynamics are disabled (`do_caco3_dynamics = .false.`), the model uses a static PIC:POC ratio and $CaCO_3$ dissolution rate. These are set as input parameters to the model.

---


### 11. Tracer tendencies

**Nitrate** (`f_no3(i,j,k)`, $NO_3$, [mol N kg<sup>-1</sup>])

$\dfrac{\Delta NO_3}{\Delta t} = \left(\Gamma_{det}^{C} + \gamma_{zoo}^{C} + \gamma_{phy}^{C} + \eta_{zoo}^{C} - \mu_{phy}^{C} \right) \cdot R^{N:C}$

**Oxygen** (`f_o2(i,j,k)`, $O_2$, [mol O<sub>2</sub> kg<sup>-1</sup>])

$\dfrac{\Delta O_2}{\Delta t} = \left( \mu_{phy}^{C} - \Gamma_{det}^{C} - \gamma_{zoo}^{C} - \gamma_{phy}^{C} - \eta_{zoo}^{C} \right) \cdot R^{O_2:C}$

**Dissolved iron** (`f_fe(i,j,k)`, $dFe$, [mol Fe kg<sup>-1</sup>])

$\dfrac{\Delta dFe}{\Delta t} = \Gamma_{det}^{C} Q_{det}^{Fe:C} + \gamma_{phy}^{C} Q_{phy}^{Fe:C} + \gamma_{zoo}^{C} Q_{zoo}^{Fe:C} +
\left( g_{zoo}^{&rarr; phy} \cdot Q_{phy}^{Fe:C} + g_{zoo}^{&rarr; det} \cdot Q_{det}^{Fe:C} \right) \left(1 - \lambda \right) \cdot \eta - 
\mu_{phy}^{Fe} - Fe_{nanop}^{&rarr;} - Fe_{scav}^{&rarr;} - Fe_{coag}^{&rarr;det} $ 

**Phytoplankton** (`f_phy(i,j,k)`, $phy$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta phy}{\Delta t} = \mu_{phy}^{C} - \Gamma_{phy}^{C} - \gamma_{phy}^{C} - g_{zoo}^{&rarr; phy}$ 

**Phytoplankton chlorophyll** (`f_pchl(i,j,k)`, $chl$, [mol Chl kg<sup>-1</sup>])

$\dfrac{\Delta chl}{\Delta t} = \mu_{phy}^{chl} - \left( \Gamma_{phy}^{C} + \gamma_{phy}^{C} + g_{zoo}^{&rarr; phy} \right) \cdot Q_{phy}^{Chl:C}$ 

**Phytoplankton iron** (`f_phyfe(i,j,k)`, $phy^{Fe}$, [mol Fe kg<sup>-1</sup>])

$\dfrac{\Delta phy^{Fe}}{\Delta t} = \mu_{phy}^{Fe} - \left( \Gamma_{phy}^{C} + \gamma_{phy}^{C} + g_{zoo}^{&rarr; phy} \right) \cdot Q_{phy}^{Fe:C}$ 

**Zoooplankton** (`f_zoo(i,j,k)`, $zoo$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta zoo}{\Delta t} = \mu_{zoo}^{C} - \Gamma_{zoo}^{C} - \gamma_{zoo}^{C}$ 

**Zoooplankton iron** (`f_zoofe(i,j,k)`, $zoo^{Fe}$, [mol Fe kg<sup>-1</sup>])

$\dfrac{\Delta zoo^{Fe}}{\Delta t} = \left( g_{zoo}^{&rarr; phy} \cdot Q_{phy}^{Fe:C} + g_{zoo}^{&rarr; det} \cdot Q_{det}^{Fe:C} \right) \lambda - \left( \Gamma_{zoo}^{C} - \gamma_{zoo}^{C} \right) \cdot Q_{zoo}^{Fe:C}$ 

**Detritus** (`f_det(i,j,k)`, $det$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta det}{\Delta t} = \Gamma_{phy}^{C} + \Gamma_{zoo}^{C} + E_{zoo}^{C} - \left( \Gamma_{det}^{C} + g_{zoo}^{&rarr; det} \right)$ 

**Detritus iron** (`f_detfe(i,j,k)`, $det^{Fe}$, [mol Fe kg<sup>-1</sup>])

$\dfrac{\Delta det^{Fe}}{\Delta t} = \Gamma_{phy}^{C} Q_{phy}^{Fe:C} + \Gamma_{zoo}^{C} Q_{zoo}^{Fe:C} + 
\left( g_{zoo}^{&rarr; phy} \cdot Q_{phy}^{Fe:C} + g_{zoo}^{&rarr; det} \cdot Q_{det}^{Fe:C} \right) \left(1 - \lambda \right) \left(1 - \eta \right) +
Fe_{scav}^{&rarr;det} + Fe_{coag}^{&rarr;det} -
\left( \Gamma_{det}^{C} + g_{zoo}^{&rarr; det} \right) \cdot Q_{det}^{Fe:C}$ 

**Calcium Carbonate** (`f_caco3(i,j,k)`, $CaCO_3$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta CaCO_3}{\Delta t} = \left( \Gamma_{zoo}^{C} + \Gamma_{phy}^{C} \right) PIC:POC + g_{zoo}^{&rarr; phy} \left(1 - F_{zoo}^{diss} \right) PIC:POC - D_{\text{CaCO}3}$ 

**Dissolved Inorganic Carbon** (`f_dic(i,j,k)`, $DIC$, [mol C kg<sup>-1</sup>])

$\dfrac{\Delta DIC}{\Delta t} = \Gamma_{det}^{C} + \gamma_{zoo}^{C} + \gamma_{phy}^{C} + \eta_{zoo}^{C} - \mu_{phy}^{C} - \dfrac{\Delta CaCO_3}{\Delta t}$ 

**Alkalinity** (`f_alk(i,j,k)`, $Alk$, [mol Eq kg<sup>-1</sup>])

$\dfrac{\Delta Alk}{\Delta t} = - \dfrac{\Delta NO_3}{\Delta t} - 2 \cdot \dfrac{\Delta CaCO_3}{\Delta t}$ 


---


### 12. Check for conservation of mass

When checks for the conservation of mass is enabled (`do_check_n_conserve = .true.` or `do_check_c_conserve = .true.`), the model will calculate the budget of nitrogen or carbon both before and after the ecosystem equations have completed. This checks that the ecosystem equations detailed above have indeed conserved the mass of both nitrogen and carbon within the ocean. 

---

### 13. Additional operations on tracers

**First**, dissolved iron concentrations are set to equal 1 nM everywhere where the depth of the water column is less than 200 metres deep. WOMBAT-lite is not considered to be a model of the coastal ocean, but rather a model of the global pelagic ocean. Given that coastal waters are not limited in dissolved iron due to substantial interactions with sediments and exchange with the land, we universally set the dissolved iron concentration in these waters to 1 nM.

**Second**, if dissolved iron concentrations dip below that measureable by operational detection limits, we reset these concentrations to this minimum  (`zfermin`, $[dFe]^{min}$, [µmol m<sup>-3</sup>]):

$[dFe]^{min} = \min\left( \max\left( 0.003 \cdot [NO_3]^{2}, 0.005 \right), 0.007 \right)$

where:
- [NO_3] is the ambient nitrate concentration in units of [mmol m<sup>-3</sup>] 

**Third**, in the absence of riverine additions of tracers, we apply a correction to add tracer where substantial dilution occurs. Both alkalinity ($Alk$) and dissolved inorganic carbon ($DIC$) are reset to minimum seawater concentrations. These concentrations are set in the input parameter namelist.

---


### 14. Calculate sinking rate of particulates

WOMBAT-lite functions with a spatially variable sinking rate of both organic and inorganic (i.e., $CaCO_3$) particulate matter. 

**Organic detritus**

We consider that base sinking rates of detritus (`wdetbio`, $\omega_{det}^{0}$, [m day<sup>-1</sup>]) are varied as a function of phytoplankton concentration (in a similar fashion to half-saturation coefficients described earlier), as well as the fraction of particulate matter that is $CaCO_3$. This approach is taken to emulate observations of varying sinking speeds [(Riley et al., 2012)](https://doi.org/10.1029/2011GB004085) and because such variations may be strongly dependent on phytoplankton community composition [(Bach et al., 2016)](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016GB005372).

In accordance with a more general Navier–Stokes drag equation and using a compilation of particle sinking speeds, [Cael et al. (2021)](https://doi.org/10.1029/2020GL091771) identified that the sinking velocity of detrital particles ($\omega_{det}$) in [m day<sup>-1</sup>]) is proportional to their diameter raised to the power of roughly 0.63, such that

$\omega_{det} \propto d^{0.63}$.

Knowing that 

$d = 6 V \pi \dfrac{1}{3}$

and given that the average volume of phytoplankton cells can be approximated by $V = (B_{phy}^{C})^{0.65}$ [(Wickman et al., 2024)](https://doi.org/10.1126/science.adk6901), we can relate $\omega_{det}$ to the biomass concentration of phytoplankton multiplied by the scaler $\omega_{det}^{0}$

$\omega_{det} = \omega_{det}^{0} \cdot \left( B_{phy}^{C} \right)^{0.21}$

This formula is identical to that presented by [Cael et al. (2021)](https://doi.org/10.1029/2020GL091771) in their Eq. (3), with the exception that we have related sinking rates to the biomass concentration of phytoplankton ($B_{phy}^{C}$) by assuming that $V = (B_{phy}^{C})^{0.65}$ based on marine phytoplankton data [(Wickman et al., 2024)](https://doi.org/10.1126/science.adk6901).

As phytoplankton concentrations are negligible beneath the euphotic zone we use $B_{phy}^{C}$ only in the uppermost grid cell (k = 1). This assumes that the sinking velocities of marine aggregates can be related to phytoplankton community composition at the surface ([Bach et al., 2016](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016GB005372); [Iversen and Lampitt, 2020](https://doi.org/10.1016/j.pocean.2020.102445)), which varies more horizontally across the ocean than vertically. Moreover, because we do not include dissolved/suspended organic matter as a tracer in WOMBAT-lite, we must also account for the large fraction of organics that are suspended and thus neutrally buoyant in the gyres. As such, we include a phytoplankton biomass threshold (`phybiot`, $B_{phy}^{thresh}$, [mmol C m<sup>-3</sup>]) above which sinking accelerates and beneath which any produced detritus emulates dissolved (neutrally buoyant) organic matter:

$\omega_{det} = \omega_{det}^{0} \cdot \max \left( 0.0, B_{phy}^{C,k=1} - B_{phy}^{thresh} \right)^{0.21}$

Sinking speeds are then accelerated depending on the fraction of $CaCO_3$ within the particulate mass by

$\omega_{det} = \omega_{det} + 10 \cdot \min \left( 1, \dfrac{B_{CaCO_3}^{C}}{B_{CaCO_3}^{C} + B_{phy}^{C}} \right)$

where we chose a maximum increase of 10 metres per day based on the more modest effect observed in mesocosm experiments ([Bach et al., 2016](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016GB005372).

Finally, we apply a linear increase to sinking speeds with depth to ensure that the trend in the concentration of detritus with depth exhibits a power law behavior, which is widely observed ([Berelson, 2001](https://doi.org/10.1016/S0967-0645(01)00102-3); [Martin et al., 1987](https://doi.org/10.1016/0198-0149(87)90086-0)), thought to be associated with a greater attenuation of more slowly sinking particles, and shows better performance than a constant sinking rate in models ([Tjiputra et al., 2020](https://gmd.copernicus.org/articles/13/2393/2020/)). This is applied after the previous equation as:

$\omega_{det} = \omega_{det} + \max\left(0.0, \dfrac{z}{5000} \cdot \left(\omega_{det}^{max} - \omega_{det} \right)\right)$


**Calcium carbonate**

The sinking rate of $CaCO_3$ is considered to always be a fraction of the sinking rate of organic detritus, such that

$\omega_{CaCO_3} = \omega_{det} \cdot \dfrac{\omega_{CaCO_3}^{0}}{\omega_{det}^{0}}$,

and where

$\omega_{CaCO_3}^{0} < \omega_{det}^{0}$

Although more dense than organic matter, $CaCO_3$ particles tend to be smaller than organic aggregates and sink at a slower rate ([De La Rocha & Passow, 2007](https://doi.org/10.1016/j.dsr2.2007.01.004); [Zhang et al., 2018](https://doi.org/10.5194/bg-15-4759-2018)). Furthermore, the shedding of coccoliths by coccolithophores, which are near-neutrally bouyant, also contributes to a slower mean sinking speed of $CaCO_3$ ([Balch et al., 2009](https://doi.org/10.1029/2008JC004902)).

The degree to which $CaCO_3$ particles sink more slowly than organic detritus is controlled by the user through altering the base sinking speed ($\omega_{CaCO_3}^{0}$) relative to the base sinking speed of organics ($\omega_{det}^{0}$).


---


### 15. Sedimentary remineralisation of tracers

WOMBAT-lite tracks the accumulation of organic carbon, organic iron and $CaCO_3$ within a sedimentary pool. The organic carbon pool, and organic nitrogen by extension using the model's static N:C ratio, contributes to bottom fluxes of dissolved inorganic carbon (DIC), nitrate ($NO_3$), oxygen ($O_2$) and alkalinity (Alk). Remineralisation of organic carbon produces DIC and $NO_3$, but removes $O_2$ and Alk. Ratios of nitrogen to carbon and oxygen to carbon are static at 16:122 and 172:122.

Alk and DIC are in turn affected by dissolution of the sedimenary $CaCO_3$ pool, which is considered as entirely calcite. Sedimentary dissolution is controlled by bottom water temperature and an estimate of the pore-water calcite saturation state ($\Omega_{cal}^{sed}$):

$\D_{cal}^{sed} = d_{cal}^{sed} (β_{hete})^{T} \max(1 - \Omega_{cal}^{+sed},  1 - \Omega_{cal}^{sed})^{4.5}$

where
- $d_{cal}^{sed}$ is a base rate of dissolution in units of [day<sup>-1</sup>], and
- $\Omega_{cal}^{+sed}$ is the maximum possible saturation state within sediment pore water.

The $Omega_{cal}^{sed}$ is calculated using the `mocsy` package for solving carbonate chemistry of seawater. For the sediments, we chose to sum the water column concentration of DIC and the organic carbon content of the sediment to approximate the interstitial (i.e., porewater) DIC concentration. We assume that the organic carbon content of the sediment (p_det_sediment) in mol/m2 is relevant over one meter, and therefore can be automatically converted to mol/m3 and then subsequently converted through the mol/kg using Rho_0. With this assumption these arrays can be added together. We add these arrays together to simulate the reducing conditions of organic-rich sediments, and to calculate a lower omega for calcite, which ensures greater rates of dissolution of CaCO3 within the sediment as organic matter accumulates.

---

### 16. Conserve tracers following burial?

