# Description of the WOMBATlite ocean biogeochemical model

```
        (\___/)  .-.   .-. .--. .-..-..---.  .--. .-----.
        / o o \  : :.-.: :: ,. :: '' :: .; :: .; :'-. .-'
       (   "   ) : :: :: :: :: :: .. ::   .':    :  : :
        \__ __/  : '' '' ;: :; :: :; :: .; :: :: :  : :
                  '.,'.,' '.__.':_;:_;:___.':_;:_;  :_;

 World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT)
```

_Contact Pearse J. Buchanan for any questions_

_Pearse.Buchanan@csiro.au_

---

## Tracers

The following are the active tracers in WOMBAT-lite

| Tracer           | Code name  | Description                                 | Units                     | Default on? |
|------------------|------------|---------------------------------------------|---------------------------|-------------|
| O<sub>2</sub>    | `p_o2`     | Dissolved oxygen                            | mol O2 kg<sup>-1</sup>    | Yes         |
| NO<sub>3</sub>   | `p_no3`    | Nitrate                                     | mol N kg<sup>-1</sup>     | Yes         |
| dFe              | `p_fe`     | Dissolved iron                              | mol Fe kg<sup>-1</sup>    | Yes         |
| Phy              | `p_phy`    | Phytoplankton                               | mol C kg<sup>-1</sup>     | Yes         |
| Zoo              | `p_zoo`    | Zooplankton                                 | mol C kg<sup>-1</sup>     | Yes         |
| Det              | `p_det`    | Detritus                                    | mol C kg<sup>-1</sup>     | Yes         |
| Chl              | `p_pchl`   | Phytoplankton chlorophyll content           | mol C kg<sup>-1</sup>     | Yes         |
| PhyFe            | `p_phyfe`  | Phytoplankton iron content                  | mol Fe kg<sup>-1</sup>    | Yes         |
| ZooFe            | `p_zoofe`  | Zoooplankton iron content                   | mol Fe kg<sup>-1</sup>    | Yes         |
| DetFe            | `p_detfe`  | Detritus iron content                       | mol Fe kg<sup>-1</sup>    | Yes         |
| DIC              | `p_dic`    | Dissolved inorganic carbon                  | mol C kg<sup>-1</sup>     | Yes         |
| Alk              | `p_alk`    | Dissolved alkalinity                        | mol Eq kg<sup>-1</sup>    | Yes         |
| CaCO<sub>3</sub> | `p_caco3`  | Calcium carbonate                           | mol C kg<sup>-1</sup>     | Yes         |
| DICp             | -          | Preformed dissolved inorganic carbon        | mol C kg<sup>-1</sup>     | No          |
| DICr             | `p_dicr`   | Remineralised dissolved inorganic carbon    | mol C kg<sup>-1</sup>     | No          |

---

## Logical controls

The following are logical statements within the `input.nml` namelist file that can be switched to TRUE or FALSE at runtime. 

| Logical                        | Description                                                                     | Default setting  |
|--------------------------------|---------------------------------------------------------------------------------|------------------|
| `do_caco3_dynamics`            | Production and dissolution of CaCO3 depends on carbon system state              | .true.           |
| `do_colloidal_shunt`           | Fraction of dissolved iron is colloids that coagulate onto sinking material     | .true.           |
| `do_two_ligands`               | Complex soluble iron using two ligands (weak + strong) rather than one          | .false.          |
| `do_burial`                    | Permanently bury a fraction of sinking detrital material into the sediments     | .false.          |
| `do_nitrogen_fixation`         | Allow implicit nitrogen fixing popoulation to add NO<sub>3</sub> to ocean       | .false.          |
| `do_benthic_denitrification`   | Allow implicit benthic bacterial population to remove NO<sub>3</sub> from ocean | .false.          |
| `do_tracer_dicp`               | Carry preformed dissolved inorganic carbon (dicp) as a tracer                   | .false.          |
| `do_tracer_dicr`               | Carry remineralised dissolved inorganic carbon (dicr) as a tracer               | .false.          |
| `do_check_n_conserve`          | Checks that the ecosystem calculations are conserving the mass of nitrogen      | .false.          |
| `do_check_c_conserve`          | Checks that the ecosystem calculations are conserving the mass of carbon        | .false.          |

We note that when `do_two_ligands` is set to `.true.`, the `ligK` diagnostic variable reflects the binding strength of the strong ligand. However, when `do_two_ligands` is set to `.false.`, this diagnostic (`ligK`) reflects the binding strength of the bulk ligand pool.

---

## Diagnostic outputs

The following are all **2D** diagnostic output variables from WOMBAT-lite.

| Diagnostic        | Description                                                                                          | Units                                |
| ----------------- | ---------------------------------------------------------------------------------------------------- | ------------------------------------ |
| `pco2`            | Surface aqueous partial pressure of CO₂                                                              | µatm                                 |
| `det_sed_remin`   | Rate of remineralisation of detritus in accumulated sediment                                         | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `det_sed_depst`   | Rate of deposition of detritus to sediment at base of water column                                   | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `det_sed_denit`   | Rate of denitrification (NO<sub>3</sub> consumption) in accumulated sediment                         | mol N m<sup>-2</sup> s<sup>-1</sup>  |
| `fbury`           | Fraction of deposited detritus permanently buried beneath sediment                                   | dimensionless                        |
| `fdenit`          | Fraction of detritus remineralised via denitrification                                               | dimensionless                        |
| `detfe_sed_remin` | Rate of remineralisation of detrital iron in accumulated sediment                                    | mol Fe m<sup>-2</sup> s<sup>-1</sup> |
| `detfe_sed_depst` | Rate of deposition of detrital iron to sediment at base of water column                              | mol Fe m<sup>-2</sup> s<sup>-1</sup> |
| `caco3_sed_remin` | Rate of remineralisation of CaCO₃ in accumulated sediment                                            | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `caco3_sed_depst` | Rate of deposition of CaCO₃ to sediment at base of water column                                      | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `zeuphot`         | Depth of the euphotic zone (1% incident light)                                                       | m                                    |
| `seddep`          | Depth of the bottom layer                                                                            | m                                    |
| `sedmask`         | Mask of active sediment points                                                                       | dimensionless                        |
| `sedtemp`         | Temperature in the bottom layer                                                                      | °C                                   |
| `sedsalt`         | Salinity in the bottom layer                                                                         | psu                                  |
| `sedno3`          | Nitrate concentration in the bottom layer                                                            | mol N kg<sup>-1</sup>                |
| `sedo2`           | Dissolved oxygen concentration in the bottom layer                                                   | mol O2 kg<sup>-1</sup>               |
| `seddic`          | Dissolved inorganic carbon concentration in the bottom layer                                         | mol C kg<sup>-1</sup>                |
| `sedalk`          | Alkalinity concentration in the bottom layer                                                         | mol Eq kg<sup>-1</sup>               |
| `sedhtotal`       | H<sup>+</sup> ion concentration in the bottom layer                                                  | mol H<sup>+</sup>kg<sup>-1</sup>     |
| `sedco3`          | CO₃<sup>2−</sup> ion concentration in the bottom layer                                               | mol C kg<sup>-1</sup>                |
| `sedomega_cal`    | Calcite saturation state in the bottom layer                                                         | dimensionless                        |
| `o2_stf`          | Surface flux of dissolved oxygen into ocean                                                          | mol O2 m<sup>-2</sup> s<sup>-1</sup> |
| `no3_stf`         | Surface flux of nitrate into ocean                                                                   | mol N m<sup>-2</sup> s<sup>-1</sup>  |
| `fe_stf`          | Surface flux of dissolved iron into ocean                                                            | mol Fe m<sup>-2</sup> s<sup>-1</sup> |
| `det_stf`         | Surface flux of detritus into ocean                                                                  | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `dic_stf`         | Surface flux of dissolved inorganic carbon into ocean                                                | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `dicp_stf`        | Surface flux of preformed dissolved inorganic carbon into ocean                                      | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `alk_stf`         | Surface flux of alkalinity into ocean                                                                | mol Eq m<sup>-2</sup> s<sup>-1</sup> |
| `no3_vstf`        | Virtual flux of nitrate into ocean due to salinity restoring/correction                              | mol N m<sup>-2</sup> s<sup>-1</sup>  |
| `dic_vstf`        | Virtual flux of dissolved inorganic carbon into ocean due to salinity restoring/correction           | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `dicp_vstf`       | Virtual flux of preformed dissolved inorganic carbon into ocean due to salinity restoring/correction | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `alk_vstf`        | Virtual flux of alkalinity into ocean due to salinity restoring/correction                           | mol Eq m<sup>-2</sup> s<sup>-1</sup> |
| `o2_btf`          | Bottom flux of dissolved oxygen into ocean                                                           | mol O2 m<sup>-2</sup> s<sup>-1</sup> |
| `no3_btf`         | Bottom flux of nitrate into ocean                                                                    | mol N m<sup>-2</sup> s<sup>-1</sup>  |
| `fe_btf`          | Bottom flux of dissolved iron into ocean                                                             | mol Fe m<sup>-2</sup> s<sup>-1</sup> |
| `dic_btf`         | Bottom flux of dissolved inorganic carbon into ocean                                                 | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `dicr_btf`        | Bottom flux of preformed dissolved inorganic carbon into ocean                                       | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `alk_btf`         | Bottom flux of alkalinity into ocean                                                                 | mol Eq m<sup>-2</sup> s<sup>-1</sup> |


The following are all **3D** diagnostic output variables from WOMBAT-lite.


| Diagnostic    | Description                                                            | Units                                            |
| ------------- | ---------------------------------------------------------------------- | ------------------------------------------------ |
| `htotal`      | Concentration of H<sup>+</sup> ion                                     | mol kg<sup>-1</sup>                              |
| `omega_ara`   | Saturation state of aragonite                                          | dimensionless                                    |
| `omega_cal`   | Saturation state of calcite                                            | dimensionless                                    |
| `co3`         | Carbonate ion concentration                                            | mol kg<sup>-1</sup>                              |
| `co2_star`    | CO2* (CO2(g) + H2CO3)) concentration                                   | mol kg<sup>-1</sup>                              |
| `radbio`      | Photosynthetically active radiation available for phytoplankton growth | W m<sup>-2</sup>                                 |
| `radmid`      | Photosynthetically active radiation at centre point of grid cell       | W m<sup>-2</sup>                                 |
| `radmld`      | Photosynthetically active radiation averaged in mixed layer            | W m<sup>-2</sup>                                 |
| `npp3d`       | Net primary productivity                                               | mol kg<sup>-1</sup> s<sup>-1</sup>               |
| `zsp3d`       | Zooplankton secondary productivity                                     | mol kg<sup>-1</sup> s<sup>-1</sup>               |
| `phy_mumax`   | Maximum growth rate of phytoplankton                                   | s<sup>-1</sup>                                   |
| `phy_mu`      | Realised growth rate of phytoplankton                                  | s<sup>-1</sup>                                   |
| `pchl_mu`     | Realised growth rate of phytoplankton chlorophyll                      | mol kg<sup>-1</sup> s<sup>-1</sup>               |
| `phy_lpar`    | Limitation of phytoplankton by light                                   | dimensionless                                    |
| `phy_kni`     | Half-saturation coefficient of nitrogen uptake by phytoplankton        | mmol m<sup>-3</sup>                              |
| `phy_kfe`     | Half-saturation coefficient of iron uptake by phytoplankton            | µmol m<sup>-3</sup>                              |
| `phy_lnit`    | Limitation of phytoplankton by nitrogen                                | dimensionless                                    |
| `phy_lfer`    | Limitation of phytoplankton by iron                                    | dimensionless                                    |
| `phy_dfeupt`  | Uptake of dFe by phytoplankton                                         | mol kg<sup>-1</sup> s<sup>-1</sup>               |
| `tri_mumax`   | Maximum growth rate of trichodesmium                                   | s<sup>-1</sup>                                   |
| `tri_lpar`    | Limitation of trichodesmium by light                                   | dimensionless                                    |
| `tri_lfer`    | Limitation of trichodesmium by iron                                    | dimensionless                                    |
| `nitrfix`     | Implicit nitrogen fixation rate (NO<sub>3</sub> production)            | mol N kg<sup>-1</sup> s<sup>-1</sup>             |
| `feIII`       | Free iron (Fe<sup>3+</sup>)                                            | mol kg<sup>-1</sup>                              |
| `ligK`        | Ligand stability constant                                              | L mol<sup>-1</sup>                               |
| `felig`       | Ligand-bound dissolved iron                                            | mol kg<sup>-1</sup>                              |
| `fecol`       | Colloidal dissolved iron                                               | mol kg<sup>-1</sup>                              |
| `feprecip`    | Precipitation of free Fe onto nanoparticles                            | mol kg<sup>-1</sup> s<sup>-1</sup>               |
| `fescaven`    | Scavenging of free Fe onto detritus (organic + inorganic)              | mol kg<sup>-1</sup> s<sup>-1</sup>               |
| `fescadet`    | Scavenging of free Fe onto organic detritus                            | mol kg<sup>-1</sup> s<sup>-1</sup>               |
| `fecoag2det`  | Coagulation of colloidal dFe onto detritus                             | mol kg<sup>-1</sup> s<sup>-1</sup>               |
| `fesources`   | Total source of dFe in water column                                    | mol kg<sup>-1</sup> s<sup>-1</sup>               |
| `fesinks`     | Total sink of dFe in water column                                      | mol kg<sup>-1</sup> s<sup>-1</sup>               |
| `phy_feupreg` | Factor up regulation of dFe uptake by phytoplankton                    | dimensionless                                    |
| `phy_fedoreg` | Factor down regulation of dFe uptake by phytoplankton                  | dimensionless                                    |
| `phygrow`     | Growth of phytoplankton                                                | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `phymorl`     | Linear mortality of phytoplankton                                      | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `phymorq`     | Quadratic mortality of phytoplankton                                   | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `zooeps`      | Zooplankton prey capture rate coefficient                              | (mmol C m<sup>-3</sup>)<sup>-2</sup> s<sup>-1</sup> |
| `zooprefphy`  | Dietary fraction of phytoplankton in zooplankton grazing               | dimensionless                                    |
| `zooprefdet`  | Dietary fraction of detritus in zooplankton grazing                    | dimensionless                                    |
| `zoograzphy`  | Grazing rate of zooplankton on phytoplankton                           | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `zoograzdet`  | Grazing rate of zooplankton on detritus                                | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `zoomorl`     | Linear mortality of zooplankton                                        | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `zoomorq`     | Quadratic mortality of zooplankton                                     | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `zooexcrphy`  | Excretion rate of zooplankton eating phytoplankton                     | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `zooexcrdet`  | Excretion rate of zooplankton eating detritus                          | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `zooassiphy`  | Assimilation into biomass of zooplankton feeding on phytoplankton      | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `zooassidet`  | Assimilation into biomass of zooplankton feeding on detritus           | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `zooegesphy`  | Egestion of zooplankton feeding on phytoplankton                       | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `zooegesdet`  | Egestion of zooplankton feeding on detritus                            | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `reminr`      | Rate of remineralisation                                               | s<sup>-1</sup>                                   |
| `detremi`     | Remineralisation of detritus                                           | mol C kg<sup>-1</sup> s<sup>-1</sup>             |
| `pic2poc`     | Inorganic (CaCO3) to organic carbon ratio                              | dimensionless                                    |
| `dissratcal`  | Dissolution rate of calcite CaCO3                                      | s<sup>-1</sup>                                   |
| `dissratara`  | Dissolution rate of aragonite CaCO3                                    | s<sup>-1</sup>                                   |
| `dissratpoc`  | Dissolution rate of CaCO3 due to POC remin                             | s<sup>-1</sup>                                   |
| `zoodiss`     | Dissolution of CaCO3 due to zooplankton grazing                        | mol CaCO3 kg<sup>-1</sup> s<sup>-1</sup>         |
| `caldiss`     | Dissolution of calcite CaCO3                                           | mol CaCO3 kg<sup>-1</sup> s<sup>-1</sup>         |
| `aradiss`     | Dissolution of aragonite CaCO3                                         | mol CaCO3 kg<sup>-1</sup> s<sup>-1</sup>         |
| `pocdiss`     | Dissolution of CaCO3 due to POC remin                                  | mol CaCO3 kg<sup>-1</sup> s<sup>-1</sup>         |
| `det_vmove`   | Sinking rate of detritus                                               | m s<sup>-1</sup>                                 |
| `detfe_vmove` | Sinking rate of detrital iron                                          | m s<sup>-1</sup>                                 |
| `caco3_vmove` | Sinking rate of CaCO3                                                  | m s<sup>-1</sup>                                 |

---

## Subroutine - "update_from_source"

The subroutine `generic_WOMBATlite_update_from_source` is the heart of the World Ocean Model of Biogeochemistry And Trophic‑dynamics. Its purpose is to apply biological source–sink terms to ocean tracers (nutrients, phytoplankton, zooplankton, detritus, iron and carbon pools) at each tracer time‑step. The subroutine is documented internally by a list of numbered steps (see code comments). These steps are:

1. Light attenuation through the water column.
2. Nutrient limitation of phytoplankton.
3. Temperature‑dependent autotrophy and heterotrophy.
4. Light limitation of phytoplankton.
5. Realized growth rate of phytoplankton.
6. Synthesis of chlorophyll.
7. Phytoplankton uptake of iron.
8. Iron chemistry.
9. Mortality and remineralisation.
10. Zooplankton grazing, growth, excretion and egestion.
11. CaCO3 calculations.
12. Implicit nitrogen fixation
13. Tracer tendencies (update tracer concentrations).
14. Check for conservation of mass.
15. Additional operations on tracers.
16. Sinking rate of particulates.
17. Sedimentary processes.

Below is a step‑by‑step explanation of each section together with the key equations. Variable names in grey follow the Fortran code, while variable names in math font are pointers to the equations; i,j,k refer to horizontal and vertical indices; square brackets denote units. If a variable is without i,j,k dimensions, this variable is held as a scalar and not an array.

The model carries tracers in [mol kg<sup>-1</sup>]. That is, moles of solute/tracer per kilogram of seawater (i.e., molality). Some calculations herein are performed by converting tracers to units of [mmol m<sup>-3</sup>] or in the case of dissolved iron [µmol m<sup>-3</sup>]. However, we stress that all tracer tendency terms are converted back to [mol kg<sup>-1</sup> s<sup>-1</sup>] when sources and sinks are applied.

---

### Parameter set and default values

| Parameter          | Description                                                                 | Value              |
|--------------------|-----------------------------------------------------------------------------|--------------------|
| `alphabio`         | Phytoplankton initial slope of P–I curve [(W/m²)⁻¹ (mg/mg)⁻¹]               | 3.0                |
| `abioa`            | Autotrophy maximum growth rate parameter a [s⁻¹]                            | 1.0/86400.0        |
| `bbioa`            | Autotrophy maximum growth rate parameter b [1]                              | 1.050              |
| `bbioh`            | Heterotrophy maximum growth rate parameter b [1]                            | 1.060              |
| `phykn`            | Phytoplankton half saturation constant for nitrogen uptake [mmolN/m³]       | 2.0                |
| `phykf`            | Phytoplankton half saturation constant for iron uptake [µmolFe/m³]          | 1.0                |
| `phyminqc`         | Phytoplankton minimum chlorophyll:C quota [mg/mg]                           | 0.004              |
| `phymaxqc`         | Phytoplankton maximum chlorophyll:C quota [mg/mg]                           | 0.060              |
| `phytauqc`         | Phytoplankton timescale for chlorophyll:C adjustment [s]                    | 86400.0            |
| `phyoptqf`         | Phytoplankton optimal Fe:C quota [mol/mol]                                  | 10e-6              |
| `phymaxqf`         | Phytoplankton maximum Fe:C quota [mol/mol]                                  | 50e-6              |
| `phylmor`          | Phytoplankton linear mortality rate constant [s⁻¹]                          | 0.005/86400.0      |
| `phyqmor`          | Phytoplankton quadratic mortality rate constant [(mmolC/m³)⁻¹ s⁻¹]          | 0.05/86400.0       |
| `phybiot`          | Phytoplankton biomass threshold [mmolC/m³]                                  | 0.6                |
| `alphabio_tri`     | Trichodesmium initial slope of P–I curve [(W/m²)⁻¹ (mg/mg)⁻¹]               | 1.8                |
| `trikf`            | Trichodesmium half saturation constant for iron uptake [µmolFe/m³]          | 0.125              |
| `trichlc`          | Trichodesmium chlorophyll to carbon ratio [mol Chl (mol C)<sup>-1</sup>]    | 0.01               |
| `trin2c`           | Trichodesmium nitrogen to carbon ratio [molN (mol C)<sup>-1</sup>]          | 50/300             |
| `zooCingest`       | Zooplankton ingestion efficiency of carbon [molC/molC]                      | 0.8                |
| `zooCassim`        | Zooplankton assimilation efficiency of carbon [molC/molC]                   | 0.3                |
| `zooFeingest`      | Zooplankton ingestion efficiency of iron [molFe/molFe]                      | 0.2                |
| `zooFeassim`       | Zooplankton assimilation efficiency iron [molFe/molFe]                      | 0.9                |
| `fgutdiss`         | CaCO₃ dissolution efficiency in zooplankton guts [molC/molC]                | 0.75               |
| `zookz`            | Half-saturation coefficient for zooplankton mortality [mmolC/m³]            | 0.25               |
| `zoogmax`          | Zooplankton maximum grazing rate [s⁻¹]                                      | 3.0/86400.0        |
| `zooepsmin`        | Zooplankton minimum prey capture rate [m⁶/mmol²/s]                          | 0.005/86400.0      |
| `zooepsmax`        | Zooplankton maximum prey capture rate [m⁶/mmol²/s]                          | 0.25/86400.0       |
| `zooepsrat`        | Zooplankton transition rate of epsilon [(mmolC/m³)⁻¹]                       | 0.1                |
| `zoopreyswitch`    | Zooplankton prey switching exponent [dimenionless]                          | 1.8                |
| `zprefphy`         | Zooplankton preference for phytoplankton [dimensionless]                    | 1.0                |
| `zprefdet`         | Zooplankton preference for detritus [dimensionless]                         | 0.50               |
| `zoolmor`          | Zooplankton respiration rate [s⁻¹]                                          | 0.0025/86400.0     |
| `zooqmor`          | Zooplankton quadratic mortality [(mmolC/m³)⁻¹ s⁻¹]                          | 0.8/86400.0        |
| `detlrem`          | Detritus remineralisation rate [(mmolC/m³)⁻¹ s⁻¹]                           | 0.5/86400.0        |
| `wdetbio`          | Base detritus sinking rate [m/s]                                            | 25.0/86400.0       |
| `wdetmax`          | Maximum detritus sinking rate [m/s]                                         | 42.0/86400.0       |
| `wcaco3`           | CaCO₃ sinking rate [m/s]                                                    | 12.5/86400.0       |
| `detlrem_sed`      | Sediment detritus remineralisation [s⁻¹]                                    | 0.01/86400.0       |
| `caco3lrem`        | Base CaCO₃ remineralisation [s⁻¹]                                           | 0.01/86400.0       |
| `caco3lrem_sed`    | Sediment CaCO₃ remineralisation [s⁻¹]                                       | 0.01/86400.0       |
| `omegamax_sed`     | Omega ceiling in sediments [dimenionless]                                   | 0.8                |
| `f_inorg`          | Base inorganic fraction of CaCO₃ within detritus [molC/molC]                | 0.045              |
| `disscal`          | Calcite dissolution factor [s⁻¹]                                            | 0.10/86400.0       |
| `dissara`          | Aragonite dissolution factor [s⁻¹]                                          | 0.10/86400.0       |
| `dissdet`          | Dissolution factor from detritus remineralisation [molC/molC]               | 0.200              |
| `ligW`             | Weak ligand background concentration [µmol/m³]                              | 1.7                |
| `ligS`             | Strong ligand background concentration [µmol/m³]                            | 0.4                |
| `dfefloor`         | Minimum dissolved Fe concentration [µmol/m³]                                | 0.05               |
| `knano_dfe`        | Fe nanoparticle precipitation rate [s⁻¹]                                    | 0.1/86400.0        |
| `kscav_dfe`        | Fe scavenging rate [(mmol/m³)⁻¹ s⁻¹]                                        | 0.01/86400         |
| `kcoag_dfe`        | Fe coagulation rate [(mmolC/m³)⁻¹ s⁻¹]                                      | 1e-5/86400         |
| `kagg_col`         | Colloidal Fe aggregation rate [s⁻¹]                                         | 0.1/86400.0        |
| `kagg_kcol`        | Half-saturation for colloidal Fe aggregation [µmolFe/m³]                    | 2.0                |
| `bottom_thickness` | Bottom layer thickness [m]                                                  | 1.0                |

---


### 1. Light attenuation through the water column.

Photosynthetically available radiation (PAR) is split into blue, green and red wavelengths. The incoming visible (photosynthetically available) short wave radiation flux (PAR, [W m<sup>-2</sup>]) is received from the physical model, and is then split evenly into each of blue, green and red light bands. 

At the top (`par_bgr_top(k,b)`, $PAR^{top}$) and mid‑point (`par_bgr_mid(k,b)`, $PAR^{mid}$) of each layer `k` we calculate the downward irradiance by exponential decay of each band `b` through the layer thickness (`dzt(i,j,k)`, $\Delta z$, [m]) using band‑specific attenuation coefficients. These attenuation coefficients are related to the concentration of chlorophyll (`chl`, [mg m<sup>-3</sup>]), organic detritus (`ndet`, [mg N m<sup>-3</sup>]) and calcium carbonate (`carb`, [kg m<sup>-3</sup>]) in the water column. 

For chlorophyll, attenuation coefficients for each of blue, green and red light (`zbgr(ichl,b)`, [m<sup>-1</sup>]) are retrieved from the look-up table of [Morel & Maritorena (2001)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2000JC000319) (their Table 2) that explicitly relates chlorophyll concentration to attenuation rates and accounts for the packaging effect of chlorophyll in larger cells. Within `zbgr(ichl,b)`, `ichl` is an integer that corresponds to a particular band of chlorophyll concentration, with increasing chlorophyll concentrations associated with increasing attenuation.  

For organic detritus, attenuation coefficients for blue, green and red light (`dbgr(b)`, [(mg N m<sup>-3</sup>)<sup>-1</sup> m<sup>-1</sup>]) are taken from [Dutkiewicz et al. (2015)](https://bg.copernicus.org/articles/12/4447/2015/bg-12-4447-2015.html) (their Fig. 1b), while for calcium carbonate (`cbgr(b)`, [(kg CaCO<sub>3</sub> m<sup>-3</sup>)<sup>-1</sup>m<sup>-1</sup>]) we take the coefficients defined in [Soja-Wozniak et al. (2019)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JC014998). For both detritus and calcium carbonate, these studies provide concentration-normalized attenuation coefficients, which must be multiplied against concentrations to retrieve the correct units of [m<sup>-1</sup>].

As an example, the PAR in the blue band (`b=1`) at the top of level k is computed as

$$
\begin{align}
PAR^{top}(k,1) =& \quad PAR^{top}(k-1,1) \cdot e^{\left(-ex_{bgr}(k-1,1) \cdot \Delta z(k-1)\right)}
\end{align}
$$

where the total attenutation rate of blue light in the grid cell above `k` is the sum of attenuation due to all particulates in that grid cell, which includes chlorophyll, detritus and calcium carbonate:

$$
\begin{align}
ex_{bgr}(k-1,1) =& \quad ex_{chl}(k-1,1) + ex_{det}(k-1,1) + ex_{CaCO_3}(k-1,1)
\end{align}
$$

_where_ <br>
- $ex_{chl}(k-1,1)$ is the attenuation rate of blue light (`b=1`) in the overlying grid cell (`k-1`) due to chlorophyll (`zbgr(2,ichl)`, [m<sup>-1</sup>]) <br>
- $ex_{det}(k-1,1)$ is the attenuation rate of blue light (`b=1`) in the overlying grid cell (`k-1`) due to detritus (`ndet * dbgr(1)`, [m<sup>-1</sup>]) <br>
- $ex_{CaCO_3}(k-1,1)$ is the attenuation rate of blue light (`b=1`) in the overlying grid cell (`k-1`) due to calcium carbonate (`carb * cbgr(1)`, [m<sup>-1</sup>]) <br>

The irradiance in the red band (`b=3`) at the mid point of layer `k`, in contrast, is equal to 

$$
\begin{align}
PAR^{mid}(k,3) =& \quad PAR^{mid}(k-1,3) \cdot e^{\left(-0.5\left(ex_{bgr}(k-1,3) \cdot \Delta z(k-1) + ex_{bgr}(k,3) \cdot \Delta z(k)\right)\right)}
\end{align}
$$

_where_ <br>
- $PAR^{mid}(k-1,3)$ is the red light (`b=3`) at the mid-point of the overlying grid cell (`par_bgr_mid(k-1,3)`, [W m<sup>-2</sup>]) <br>
- $ex_{bgr}(k-1,3)$ is the total attenuation of red light (`b=3`) in the overlying grid cell (`ek_bgr(k-1,3)`, [m<sup>-1</sup>]) <br>
- $ex_{bgr}(k,3)$ is the total attenuation of red light (`b=3`) in the current grid cell (`ek_bgr(k,3)`, [m<sup>-1</sup>]) <br>
- $\Delta z(k-1)$ and $\Delta z(k)$ are the grid cell thicknesses of the overlying and current grid cells (`dzt(i,j,k)`, [m]) <br>

The total PAR available to phytoplantkon is assumed to be the sum of the blue, green and red bands. Because we assume that phytoplankton are homogenously distributed within a layer `k`, but we do not assume that light is homogenously distributed within that layer, we solve for the PAR that is seen by the average phytoplankton within that grid cell (`radbio`, $PAR$, [W m<sup>-2</sup>])

$$
\begin{align}
PAR(k) =& \quad \sum_{b=1}^3 \dfrac{PAR^{top}(k,b) - PAR^{top}(k+1,b)}{ex_{bgr}(k,b) \cdot \Delta z(k)}
\end{align}
$$

_where_ <br>
- $PAR^{top}(k,b)$ is the incoming photosynthetically active radiation at the top of grid cell `k` and light band `b` (`par_bgr_top(k,b)`, [W m<sup>-2</sup>]) <br>
- $ex_{bgr}(k,b)$ is the attenuation rate of light band `b` in grid cell `k` (`ek_bgr(k,b)`, [m<sup>-1</sup>]) <br>
- $\Delta z(k)$ is the grid cell thickness of grid cell `k` (`dzt(i,j,k)`, [m]) <br>

This ensures phytoplankton growth in the model responds to the mean light they experience in the grid cell. See Eq. 19 from [Baird et al. (2020)](https://gmd.copernicus.org/articles/13/4503/2020/).

The euphotic depth (`zeuphot(i,j)`, [m]) is defined as the depth where `radbio` falls below the 1% threshold of incidient shortwave radiation or below 0.01 W m<sup>-2</sup>, whichever is shallower.

---


### 2. Nutrient limitation of phytoplankton.

At the start of each vertical loop the code computes biomass of phytoplankton (`phy_mmolm3`, $B_{phy}^{C}$, [mmol C m<sup>-3</sup>]). Phytoplankton biomass is used to scale how both nitrate (`no3_mmolm3`, $NO_{3}$, [mmol NO<sub>3</sub> m<sup>-3</sup>]) and dissolved iron (`fe_umolm3`, $dFe$, [µmol Fe m<sup>-3</sup>]) affect the growth of phytoplankton. Using compilations of marine phytoplankton and zooplankton communities, [Wickman et al. (2024)](https://www.science.org/doi/10.1126/science.adk6901) show that the nutrient affinity, $\mathit{aff}$, of a phytoplankton cell is related to its volume, $V$, via

$$
\begin{align}
\mathit{aff} =& \quad V^{-0.57}
\end{align}
$$

Additionally, the authors demonstrate that the volume of the average phytoplankton cell is related to the density (i.e., concentration) of phytoplankton via

$$
\begin{align}
V =& \quad (B_{phy}^{C})^{0.65}
\end{align}
$$
 
when combining panels c and f of their Figure 1. This then relates the affinity of an average cell to the concentration of phytoplankton biomass as

$$
\begin{align}
\mathit{aff} =& \quad (B_{phy}^{C})^{-0.37}
\end{align}
$$

With this information, we allow the half-saturation terms for nitrate (`phy_kni(i,j,k)`, $K_{phy}^{N}$, [mmol N m<sup>-3</sup>]) and dissolved iron  (`phy_kfe(i,j,k)`, $K_{phy}^{Fe}$, [µmol Fe m<sup>-3</sup>]) uptake to vary as a function of phytoplankton biomass concentration. We set reference values for the half-saturation coefficient of nitrate (`phykn`, $K_{phy}^{N,0}$, [mmol N m<sup>-3</sup>]) and dissolved iron (`phykf`, $K_{phy}^{Fe,0}$, [µmol dFe m<sup>-3</sup>]) as input parameters to the model, and also set a threshold phytoplankton concentration (`phybiot`, $B_{phy}^{thresh}$, [mmol C m<sup>-3</sup>]) beneath which cell size cannot decrease and affinity can no longer increase. At this minimum, where affinity is maximised, the half-saturation coefficients are bounded to be 10% of their reference values.

$$
\begin{align}
K_{phy}^{N} =& \quad K_{phy}^{N,0} \cdot \max\left(0.1, \max\left(0.0, B_{phy}-B_{phy}^{thresh}\right)^{0.37} \right) \\
K_{phy}^{Fe} =& \quad K_{phy}^{Fe,0} \cdot \max\left(0.1, \max\left(0.0, B_{phy}-B_{phy}^{thresh}\right)^{0.37} \right)
\end{align}
$$

**Limitation of phytoplankton growth by nitrate** (`phy_lnit(i,j,k)`, $L_{phy}^{N}$), [dimensionless]) then follows the Monod equation:

$$
\begin{align}
L_{phy}^{N} =& \quad \dfrac{NO_3}{NO_3 + K_{phy}^{N}}
\end{align}
$$
 
_where _ <br>
- $NO_3$ is the ambient nitrate concentration (`no3_mmolm3`, $NO_3$, [mmol N m<sup>-3</sup>]) <br>
- $K_{phy}^{N}$ is the michaelis-menten half-saturation coefficient (`phy_kni(i,j,k)`, [mmol N m<sup>-3</sup>]) <br>

**Limitation of phytoplankton growth by iron** follows an internal quota approach ([Droop, 1983](https://www.degruyterbrill.com/document/doi/10.1515/botm.1983.26.3.99/html)). Phytoplankton have a minimum iron quota (`phy_minqfe`, $Q_{phy}^{-Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]) and an optimal quota for growth (`phyoptqf`, $Q_{phy}^{*Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]). The minimum iron quota, $Q_{phy}^{-Fe:C}$, is dependent on the chlorophyll content of the cell and the degree of nitrogen limitation according to

$$
\begin{align}
Q_{phy}^{-Fe:C} =& \quad \dfrac{0.00167}{55.85} \cdot \max\left( Q_{phy}^{Chl:C}, Q_{phy}^{-Chl:C} \right) \cdot 12 \\
& \quad + 1.21 \times 10^{-5} \cdot \dfrac{14}{55.85 \cdot 7.625} \cdot 0.5 \cdot 1.5 \cdot L_{phy}^{N} \\
& \quad + 1.15 \times 10^{-4} \cdot \dfrac{14}{55.85 \cdot 7.625} \cdot 0.5 \cdot L_{phy}^{N}
\end{align}
$$

The first term reflects the amount of iron required for photosystems I and II. $\dfrac{0.00167}{55.85}$ is equivalent to the grams of Fe per gram of chlorophyll divided by the grams of Fe per mol Fe, giving mol Fe per gram chlorophyll. This term is multipled by the chlorophyll to carbon ratio of the phytoplantkon cell (`phy_chlc`, $Q_{phy}^{Chl:C}$, [mol C (mol C)<sup>-1</sup>]), taking into account the minimum possible chlorophyll to carbon ratio ($Q_{phy}^{-Chl:C}$), and grams of C per mol C, returning mol Fe per mol C. At a healthy chlorophyll:C ratio of 0.03, this term returns an Fe:C ratio of roughly 10 µmol:mol, which reproduces well known requirements of phytoplankton cells ([Morel, Rueter & Price, 1991](https://www.jstor.org/stable/43924569)). The second term, representing the respiratory iron requirement, is derived from [Flynn & Hipkin (1999)](https://onlinelibrary.wiley.com/doi/10.1046/j.1529-8817.1999.3561171.x) who estimated 1.21 $\times 10^{-5}$ grams Fe per gram N assimilated into the cell, which is converted to mol Fe per mol C with 14 g N per mol N divided by 55.85 g Fe per mol Fe $\times$ 7.625 mol C per mol N. This second term assumes that respiration is reduced as growth becomes more limited by available nitrogen (`phy_lnit(i,j,k)`, $L_{phy}^{N}$, [dimensionless]). Finally, the third term represents the iron required by nitrate/nitrite reduction. Nitrate assimilation requires roughly 1.8$\times$ more iron than ammonia assimilation ([Raven, 1988](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/j.1469-8137.1988.tb04196.x)). [Flynn & Hipkin (1999)](https://onlinelibrary.wiley.com/doi/10.1046/j.1529-8817.1999.3561171.x) estimated a demand of 1.15 $\times 10^{-4}$ g Fe per mol $NO_3$ reduced, which is accounted for by the nitrogen limitation term. Note that the 1.5 in the equation for $Q_{phy}^{-Fe:C}$ is designed to account for dark respiration (i.e., respiration when the cells are not growing) and the 0.5 refers to the fact that during cell division the cell must reinstate half of its Fe reserves.

The Fe limitation factor (`phy_lfer(i,j,k)`, $L_{phy}^{Fe}$) is then computed from the present Fe:C quota of the phytoplankton cells (`phy_Fe2C`, $Q_{phy}^{Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]) relative to the minimum and optimal quotas.

$$
\begin{align}
L_{phy}^{Fe} =& \quad \max\left(0.0, \min\left(1.0, \dfrac{ Q_{phy}^{Fe:C} - Q_{phy}^{-Fe:C} }{Q_{phy}^{*Fe:C}} \right)\right)
\end{align}
$$

If the cell is Fe‑replete with a quota that exceeds the minimum quota by as much as the optimal quota, then Fe does not limit growth ($L_{phy}^{Fe}$ = 1). If the cell is Fe‑deplete with a quota equal to or less than the minimum quota, then the growth rate is reduced to zero. The optimal quota ($Q_{phy}^{*Fe:C}$) is therefore a measure of how much excess Fe is required to allow unrestricted growth. The minimum iron requirements of the cell for growth increases as the ratio of chlorophyll to carbon increases (`phy_chlc`, $Q_{phy}^{Chl:C}$, [mol C (mol C)<sup>-1</sup>]), as respiration increases and as the cell uses more nitrate. This formulation and the coefficients applied to chlorophyll content and nitrate use derive from [Flynn & Hipkin (1999)](https://onlinelibrary.wiley.com/doi/10.1046/j.1529-8817.1999.3561171.x).


---


### 3. Temperature‑dependent autotrophy and heterotrophy.

**Autotrophy**

The maximum potential growth rate for phytoplankton (`phy_mumax(i,j,k)`, $\mu_{phy}^{max}$, [s<sup>-1</sup>]) is prescribed by the temperature-dependent Eppley curve ([Eppley, 1972](https://spo.nmfs.noaa.gov/content/temperature-and-phytoplankton-growth-sea)). This formulation scales a reference growth rate at 0ºC via a power-law scaling with temperature (`Temp(i,j,k)`, $T$, [ºC]).

$$
\begin{align}
\mu_{phy}^{max} =& \quad \mu_{phy}^{0^{\circ}C} \cdot \left(β_{auto}\right)^{T}
\end{align}
$$

_where_ <br>
- $\mu_{phy}^{0ºC}$ is the rate of phytoplankton growth at 0ºC (`abioa`, [s<sup>-1</sup>]) <br>
- $β_{auto}$ is the base temperature-sensitivity coefficient for autotrophy (`bbioa`, [dimenionless]) <br>

In the above, both $\mu_{phy}^{0ºC}$ and $β_{auto}$ are reference values input to the model and control how productive the ocean is and can therefore be modified to reflect known variations in the response of different phytoplankton types to temperature ([Anderson et al., 2021](https://www.nature.com/articles/s41467-021-26651-8)).

**Heterotrophy**

Heterotrophic processes include mortality of phytoplankton and zooplankton, grazing rates of zooplankton and the remineralisation rate of detritus in the water column and sediments. These processes are scaled similarly to autotrophy, where some reference rate at 0ºC ($\mu_{het}^{0^{\circ}C}$, [s<sup>-1</sup>]) is multiplied by a power-law with temperature ($β_{hete}$). Each heterotrophic process has a different $\mu_{het}^{0^{\circ}C}$ value and we expand on this later under mortality and grazing sections. However, the basic formulation for scaling heterotrophic metabolisms with temperature takes the form:

$$
\begin{align}
\mu_{het} =& \quad \mu_{het}^{0^{\circ}C} \cdot \left(β_{hete}\right)^{T}
\end{align}
$$

_where_  <br>
- $\mu_{het}^{0ºC}$ is the rate of some heterotrophic metabolism at 0ºC ([s<sup>-1</sup>]) <br>
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless]) <br>

See sections below for further details on heterotrophic metabolisms of mortality and grazing.

---


### 4. Light limitation of phytoplankton.

Phytoplankton growth is limited by light through a photosynthesis–irradiance (P–I) relationship that links cellular chlorophyll content and available photosynthetically active radiation (`radbio`, $PAR$, [W m<sup>-2</sup>]).

First, The initial slope of the P–I curve, (`phy_pisl`, $\alpha_{phy}$, [s<sup>-1</sup> (W m<sup>-2</sup>)<sup>-1</sup>]), determines how efficiently phytoplankton convert light into carbon fixation. It is scaled by the cellular chlorophyll-to-carbon ratio (`phy_chlc`, $Q_{phy}^{Chl:C}$, [mol C (mol C)<sup>-1</sup>]).

$$
\begin{align}
\alpha_{phy} =& \quad \max\left( Q_{phy}^{Chl:C}\ , \ Q_{phy}^{-Chl:C} \right) \alpha_{phy}^{Chl}
\end{align}
$$

_where_ <br>
- $\alpha_{phy}^{Chl}$ is the photosynthetic efficiency per unit chlorophyll (`alphabio`, [s<sup>-1</sup> (W m<sup>-2</sup>)<sup>-1</sup> (mg/mg)<sup>-1</sup>]) <br>
- $Q_{phy}^{Chl:C}$ is the in situ chlorophyll to carbon ratio (`phy_chlc`, [mol C (mol C)<sup>-1</sup>]) <br>
- $Q_{phy}^{-Chl:C}$ is the minimum chlorophyll to carbon ratio of the cell (`phyminqc`, [mol C (mol C)<sup>-1</sup>]) <br>

This constraint prevents photosynthesis from collapsing unrealistically at low chlorophyll concentrations. These values are parameter inputs at run time and can differ across phytoplankton taxa ([Edwards et al., 2015](https://doi.org/10.1002/lno.10033), [Litchman 2022](https://link.springer.com/chapter/10.1007/978-3-030-92499-7_1)).

Second, light limitation (`phy_lpar(i,j,k)`, $L_{phy}^{PAR}$), [dimensionless]) is calculated using an exponential P–I formulation.

$$
\begin{align}
L_{phy}^{PAR} =& \quad 1 - e^{\left(- \alpha_{phy} PAR \right)}
\end{align}
$$

_where_ <br>
- $PAR$ is the total photosynthetically available radiation (`radbio`, [W m<sup>-2</sup>]) <br>

At low irradiance (PAR), growth increases approximately linearly with light, while at high irradiance photosynthesis asymptotically saturates. Photoinhibition is not included in this formulation.

---


### 5. Realized growth rate of phytoplankton. 

Realized growth of phytoplankton (`phy_mu(i,j,k)`, $\mu_{phy}$, [s<sup>-1</sup>]) is calculated as:

$$
\begin{align}
\mu_{phy} =& \quad \mu_{phy}^{max} L_{phy}^{PAR} \min\left(L_{phy}^{N}, L_{phy}^{Fe}\right)
\end{align}
$$

_where_ <br>
- $\mu_{phy}^{max}$ is the maximum potential rate of carbon fixation (`phy_mumax`, [s<sup>-1</sup>]) <br>
- $L_{phy}^{PAR}$ is the growth limiter by light (`phy_lpar(i,j,k)`, [dimensionless]) <br>
- $L_{phy}^{N}$ is the growth limiter by nitrogen (`phy_lnit(i,j,k)`, [dimensionless]) <br>
- $L_{phy}^{Fe}$ is the growth limiter by iron (`phy_lfer(i,j,k)`, [dimensionless]) <br>

We apply Liebig's law of the minimum ([Liebig, 1840](https://archive.org/details/organicchemistry00liebrich/mode/2up), [Blackman, 1905](https://doi.org/10.1093/oxfordjournals.aob.a089000)) to resources that are required for biomass synthesis (N and Fe), while light is considered multiplicative because it is an energy supply constraint that powers nutrient aquisition and biomass synthesis.

Carbon fixation by phytoplankton (`phygrow(i,j,k)`, $\mu_{phy}^{\leftarrow C}$, [mmol C m<sup>-3</sup> s<sup>-1</sup>]) is then calculated as:

$$
\begin{align}
\mu_{phy}^{\leftarrow C} =& \quad \mu_{phy} B_{phy}^{C}
\end{align}
$$

_where_ <br>
- $\mu_{phy}^{\leftarrow C}$ is the realized rate of carbon biomass growth by phytoplankton (`phygrow(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $B_{phy}^{C}$ is the carbon concentration of phytoplankton (`phy_mmolm3`, [mmol C m<sup>-3</sup>]) <br>

---


### 6. Synthesis of chlorophyll.

This step diagnoses the **rate of chlorophyll synthesis** as a function of mixed-layer light, the phytoplankton growth rate and nutrient availability. The structure is consistent with the [Geider, MacIntyre & Kana (1997)](https://doi.org/10.3354/meps148187) formulation that relaxes the chlorophyll-to-carbon ratio towards an optimal value that supports photosynthetic growth under prevailing light and nutrient conditions.

We first solve for the optimal chlorophyll-to-carbon ratio (`theta_opt`, $Q_{phy}^{*Chl:C}$, [mol C (mol C)<sup>-1</sup>]), which is diagnosed as the ratio required to support maximal photosynthetic carbon fixation under the ambient mean light level in the mixed layer, while accounting for nutrient limitation:

$$
\begin{align}
Q_{phy}^{*Chl:C} =& \quad \dfrac{Q_{phy}^{+Chl:C}}{1 + \dfrac{\alpha_{phy} PAR_{MLD} Q_{phy}^{+Chl:C}}{2 \mu_{phy}^{max} \min \left(L_{phy}^{N}\ , \ L_{phy}^{Fe} \right) }}
\end{align}
$$

_where_ <br>
- $Q_{phy}^{+Chl:C}$ is the maximum allowable chlorophyll-to-carbon ratio (`phymaxqc`, [mol C (mol C)<sup>-1</sup>]) <br>
- $\alpha_{phy}^{Chl}$ is the chlorophyll-specific initial slope of the P–I curve (`alphabio`, [s<sup>-1</sup> (W m<sup>-2</sup>)<sup>-1</sup> (mg/mg)<sup>-1</sup>]) <br>
- $PAR_{MLD}$ is mean photosynthetically available radiation over the mixed layer (`radmld`, [W m<sup>-2></sup>]) <br>
- $\mu_{phy}^{max}$ is the temperature-dependent maximum phytoplankton growth rate (`phy_mumax`, [s<sup>-1</sup>]) <br>
- $L_{phy}^{N}$ and $L_{phy}^{Fe}$ are the nutrient limitation factors for growth on N and Fe (`phy_lnit(i,j,k)`, [dimensionless]) (`phy_lfer(i,j,k)`, [dimensionless]) <br>

We set a floor for the minimum chlorophyll-to-carbon ratio of phytoplankton via:

$$
\begin{align}
Q_{phy}^{*Chl:C} =& \quad \min \left( Q_{phy}^{*Chl:C}, Q_{phy}^{-Chl:C} \right)
\end{align}
$$

_where_ <br>
- $Q_{phy}^{-Chl:C}$ is the minimum allowable chlorophyll-to-carbon ratio (`phyminqc`, [mol C (mol C)<sup>-1</sup>]) <br>

Growth of chlorophyll (`pchl_mu(i,j,k)`, $\mu_{phy}^{\leftarrow Chl}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) is then calculated as:

$$
\begin{align}
\mu_{phy}^{\leftarrow Chl} =& \quad \mu_{phy} B_{phy}^{Chl} + \dfrac{ Q_{phy}^{*Chl:C} - Q_{phy}^{Chl:C} }{\tau^{Chl}} \cdot B_{phy}^{C} \\
\end{align}
$$

_where_ <br>
- $Q_{phy}^{Chl:C}$ is the in-situ chlorophyll-to-carbon ratio (`phy_chlc`, [mol C (mol C)<sup>-1</sup>]) <br>
- $B_{phy}^{Chl}$ is the in-situ concentation of phytoplankton chlorophyll (`p_pchl(i,j,k,tau)`, [mol C kg<sup>-1</sup>]) <br>
- $\mu_{phy}$ is the realized growth rate of phytoplankton (`phy_mu(i,j,k)`, [s<sup>-1</sup>]) <br>
- $\tau^{Chl}$ is the timescale over which chlorophyll synthesis occurs within the cell (`chltau`, [s]) <br>
- $B_{phy}^{C}$ is the in-situ concentation of phytoplankton carbon (`p_phy(i,j,k,tau)`, [mol C kg<sup>-1</sup>]) <br>

This formulation elevates chlorophyll-to-carbon ratios in low light and supresses synthesis when nutrients are low. $\tau_{phy}^{Chl}$ is an input parameter at run time and should ideally be less than the doubling time of phytplankton given that phytoplankton can internally regulate their chlorophyll stores at rates greater than their overall growth.

---


### 7. Phytoplankton uptake of iron.

Like chlorophyll, the iron content of phytoplankton is explicitly tracked as a tracer in WOMBAT-lite. First, a maximum Fe quota is found dependent on the maximum quota of Fe within the phytoplankton type and the phytoplankton concentration in the water column:

$$
\begin{align}
B_{phy}^{+Fe} =& \quad B_{phy}^{C} Q_{phy}^{+Fe:C}
\end{align}
$$

_where_ <br>
- $B_{phy}^{+Fe}$ is the maximum Fe quota of the cell (`phy_maxqfe`, [mmol Fe m<sup>-3</sup>]) <br>
- $Q_{phy}^{+Fe:C}$ is the maximum Fe:C ratio of the cell (`phymaxqf`, [mol Fe (mol C)<sup>-1</sup>]) <br>
- $B_{phy}^{C}$ is the in-situ concentation of phytoplankton carbon (`p_phy(i,j,k,tau)`, [mmol C m<sup>-3</sup>]) <br>

Following [Aumont et al. (2015)](https://gmd.copernicus.org/articles/8/2465/2015/), this rate is scaled by three terms relating to (i) the michaelis-menten type affinity for dFe, (ii) up-regulation of dFe uptake representing investment in transporters when cell quotas are limiting to growth, and (iii) down regulation of dFe uptake associated with enriched cellular quotas:

$$
\begin{align}
(i) & \qquad \dfrac{dFe}{dFe + K_{phy}^{Fe}} \\
(ii) & \qquad 4 - \dfrac{4.5 L_{phy}^{Fe}}{0.5 + L_{phy}^{Fe}} \\
(iii) & \qquad \max\left(0, 1 - \dfrac{B_{phy}^{Fe} / B_{phy}^{+Fe}}{\left|1.05 - B_{phy}^{Fe} / B_{phy}^{+Fe}\right|} \right) \\
\end{align}
$$

_where_ <br>
- $dFe$ is the in situ dissolved iron concentation (`fe_umolm3`, [µmol Fe m<sup>-3</sup>]) <br>
- $K_{phy}^{Fe}$ is the half-saturation coefficient for dFe uptake by phytoplankton (`phy_kfe(i,j,k)`, [µmol Fe m<sup>-3</sup>]) <br>
- $L_{phy}^{Fe}$ is the growth limiter by iron (`phy_lfer(i,j,k)`, [dimensionless]) <br>
- $B_{phy}^{Fe}$ is the in situ Fe quota of the cell (`phyfe_mmolm3`, [mmol Fe m<sup>-3</sup>]) <br>
- $B_{phy}^{+Fe}$ is the maximum Fe quota of the cell (`phy_maxqfe`, [mmol Fe m<sup>-3</sup>]) <br>

Note that we additionally include a fourth term that decreases the maximum dFe uptake of a cell under light limitation. This is informed by slower uptake of Fe by cells grown in darkness compared to those grown in light by roughly 10-fold ([Strzepek et al., 2025](https://doi.org/10.1093/ismejo/wraf015)), which may be due to physiological stimulation of Fe uptake machinery or photoreduction of ligand-bound iron complexes ([Kong et al., 2023](https://doi.org/10.1002/lno.12331); [Maldonado et al., 2005](https://doi.org/10.1029/2005GB002481)), or possibly a combination of both. To obtain a 10-fold relative increase in Fe uptake rates under light, we applied the following term:

$$
\begin{align}
(iv) & \qquad \max\left(0.01, L_{phy}^{PAR}\right)^{0.5}
\end{align}
$$

_where_ <br>
- $L_{phy}^{PAR}$ is the growth limiter by light (`phy_lpar(i,j,k)`, [dimensionless]) <br>

Under very low light, this fourth term reduces maximum potential Fe uptake by 10-fold than what it otherwise would be. All four terms are dimensionless and are designed to scale dissolved iron uptake either up or down. Dissolved iron uptake by phytoplankton (`phy_dfeupt(i,j,k)`, $\mu_{phy}^{\leftarrow dFe}$, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]) is then calculated as:

$$
\begin{align}
\mu_{phy}^{\leftarrow dFe} =& \quad \mu_{phy}^{max} B_{phy}^{+Fe} \cdot (i) \cdot (ii) \cdot (iii) \cdot (iv)
\end{align}
$$

_where_ <br>
- $\mu_{phy}^{\leftarrow dFe}$ is the realized uptake rate of dissolved iron by phytoplankton (`phy_dfeupt(i,j,k)`, [mmol Fe m<sup>-3</sup>]) <br>
- $\mu_{phy}^{max}$ is the temperature-dependent maximum phytoplankton growth rate (`phy_mumax`, [s<sup>-1</sup>]) <br>


---


### 8. Iron chemistry.

Treatment of dissolved iron (`fe_umolm3`, $dFe$, [nmol Fe kg<sup>-1</sup>]) follows a combination of [Aumont et al. (2015)](https://gmd.copernicus.org/articles/8/2465/2015/) and [Tagliabue et al. (2023)](https://www.nature.com/articles/s41586-023-06210-5). Our calculations involve: <br>
1. Solving for the distinct pools of dissolved iron: free iron, ligand-bound iron and colloidal iron. <br>
2. Computing precipitation of free iron into nanoparticles that are permanently lost. <br>
3. Computing scavenging of free iron onto sinking organic particles. <br>
4. Computing coagulation of colloidal iron onto sinking organic particles. <br>

_NOTE: WOMBAT-lite differs from WOMBAT-mid in that sinking authigenic pools of iron are not resolved._

We first estimate the **solubility of free Fe from Fe<sup>3+</sup>** in solution using temperature, pH and salinity using the thermodynamic equilibrium equations of [Liu & Millero (2002)](https://www.sciencedirect.com/science/article/abs/pii/S0304420301000743). 

$$
\begin{align}
T_K =& \quad \max(5.0, T) + 273.15 \\
(T_K)^{-1} =& \quad \dfrac{1}{T_K} \\
I_{S} =& \quad \dfrac{19.924 S}{1000 - 1.005 S}
\end{align}
$$

Solubility constants:

$$
\begin{align}
Fe_{sol1} =& \quad 10^{\left(-13.486 - 0.1856\sqrt{I_S} + 0.3073 I_S + 5254,\ (T_K)^{-1}\right)} \\
Fe_{sol2} =& \quad 10^{\left(2.517 - 0.8885\sqrt{I_S} + 0.2139 I_S - 1320,\ (T_K)^{-1}\right)} \\
Fe_{sol3} =& \quad 10^{\left(0.4511 - 0.3305\sqrt{I_S} - 1996,\ (T_K)^{-1}\right)} \\
Fe_{sol4} =& \quad 10^{\left(-0.2965 - 0.7881\sqrt{I_S} - 4086,\ (T_K)^{-1}\right)} \\
Fe_{sol5} =& \quad 10^{\left(4.4466 - 0.8505\sqrt{I_S} - 7980,\ (T_K)^{-1}\right)}
\end{align}
$$

Final Fe(III) solubility:

$$
\begin{align}
dFe_{sol} =& \quad Fe_{sol1}\left([H^+]^3 + Fe_{sol2}[H^+]^2 + Fe_{sol3}[H^+] + Fe_{sol4} + \dfrac{Fe_{sol5}}{[H^+]}\right)\times10^{9}
\end{align}
$$

_where_ <br>
- $T_K$ is in situ water temperature (`ztemk`, [ºK]) <br>
- $I_{S}$ is a salinity coefficient (`zval`, [dimenionless]) <br>
- $[H^+]$ is in situ hydrogen ion concentration (`hp`, [mol L<sup>-1</sup>]) <br>
- $\times10^{9}$ converts [mol Fe kg<sup>-1</sup>] to [nmol Fe kg<sup>-1</sup>] <br>
- $dFe_{sol}$ is the final estimated solubility of dissolve iron in seawater (`fe3sol`, [nmol Fe kg<sup>-1</sup>]). <br>

Next we **estimate the concentration of colloidal iron** in solution following [Tagliabue et al. 2023](https://www.nature.com/articles/s41586-023-06210-5) in the case that `do_colloidal_shunt == .true.`. If `do_colloidal_shunt == .false.` we consider no dissolved Fe to be in colloidal form. Colloidal dissolved Fe (`fecol(i,j,k)`, $dFe_{col}$, [mmol Fe m<sup>-3</sup>]) is whatever exceeds the inorganic solubility ceiling (`fe3sol`, $dFe_{sol}$, [mmol Fe m<sup>-3</sup>]), but we enforce a hard minimum that colloids are at least 10% of total dissolved Fe (`fe_umolm3`, $dFe$, [mmol Fe m<sup>-3</sup>]).

$$
\begin{align}
dFe_{col} =& \quad \max\left(0.1 dFe, \ dFe - dFe_{sol} \right)
\end{align}
$$

Following solving for colloidal dFe, we **partition the remaining dissolved Fe into ligand-bound and free iron**. To do so, we find the remaining dissolved iron not in colloidal form (`fe_sfe`, $dFe_{sFe}$, [mmol Fe m<sup>-3</sup>]), 

$$
\begin{align}
dFe_{sFe} =& \quad \max\left(0.0,\ dFe - dFe_{col} \right)
\end{align}
$$

Partitioning of iron between free and ligand-bound forms is done using one of two approaches.

When `do_two_ligands == .false.`, we use a single ligand class and solve for the equilibrium fractionation between ligand-bound and free iron using a standard quadratic form. When `do_two_ligands == .true.`, we assume complexation of iron by a weak and a strong ligand and therefore solve for the equilibrium fractionation between free iron, weakly ligand-bound iron and strongly ligand-bound iron via an iterative root solver.

In either case, we first determine the conditional stability constant(s) of the ligand(s). In the case of `do_two_ligands == .true.`, we solve for the stability constant of a strong ligand (`ligK(i,j,k)`, $Lig_{s}^{K}$, [kg mol<sup>-1</sup>]) and then consider the stability constant of a weak ligand (`ligW_K`) to be a constant offset equal to -1.5 log<sub>10</sub> units based on [Gledhill & Buck (2012)](https://doi.org/10.3389/fmicb.2012.00069). In the case of `do_two_ligands == .false.`, we solve for the stability constant of the strong (`ligK(i,j,k)`) and weak ligands (`ligW_K`), but take the concentration-weighted average binding strength to get the bulk ligand binding strength. 

The stability constant (`ligK(i,j,k)`, $Lig_{s}^{K}$, [kg mol<sup>-1</sup>]) is known to vary with the environmental conditions. In WOMBAT-lite, we consider the effect of temperature, light, pH and the concentration of labile DOC on the binding strength. The temperature dependency comes from [Volker & Tagliabue (2015)](https://doi.org/10.1016/j.marchem.2014.11.008) and warmer waters increase binding strength. The light-dependency accounts for the photoreduction of photoreactive ligands, which was identified to reduce the conditional stability constant of aquachelin by 0.7 log<sub>10</sub> units ([Barbeau et al., 2001](https://doi.org/10.1038/35096545); [Vraspir & Butler, 2009](https://doi.org/10.1146/annurev.marine.010908.163712)). The pH and DOC concentration dependency comes from [Ye et al. (2020)](https://doi.org/10.1029/2019GB006425) and increases binding strength at lower pH and higher concentrations of DOC.

$$
\begin{align}
Lig_{s}^{K} =& \quad 10^{-9} \cdot \bigg( 10^{ \left(17.27 - 1565.7 \left(T_K\right)^{-1} \right)} \\
             & \qquad  10^{\left(-0.7 \dfrac{PAR}{PAR + 10}\right)} \\
             & \qquad 10^{\left(-0.0002 [DOC]^{2} + 0.034 [DOC] - 1.67 \cdot pH + 24.36\right)} \bigg)
\end{align}
$$

_where_ <br>
- $T_K$ is in situ water temperature (`ztemk`, [ºK]) <br>
- $PAR$ is the total photosynthetically available radiation (`radbio`, [W m<sup>-2</sup>]) <br>
- pH is the in situ pH <br>
- $[DOC]$ is an empirical concentration of DOC and is equal to $40 + 40 \left( 1 - \min\left(L_{phy}^{N} \ , \ L_{phy}^{Fe}\right) \right)$ (`biodoc`, [mmol m<sup>-3</sup>]) <br>

After finding $Lig_{s}^{K}$ we solve for the free dissolved Fe concentration (`feIII`, $dFe_{free}$, [nmol Fe kg<sup>-1</sup>]) via the analytic method when `do_two_ligands == .false.`:

$$
\begin{align}
z =& \quad 1.0 + [Ligand] \cdot Lig_{bulk}^{K} - dFe_{sFe}\cdot Lig_{bulk}^{K} \\
dFe_{free} =& \quad \dfrac{-z + \sqrt{z^2 + 4.0 Lig_{bulk}^{K} dFe_{sFe}}}{2 Lig_{bulk}^{K} + \varepsilon} \\
dFe_{free} =& \quad \max\left(0, \min(dFe_{free}, dFe_{sFe})\right)
\end{align}
$$

_where_ <br>
- $[Ligand]$ is the in situ concentration of bulk ligands and in this case, where `do_two_ligands == .false.`, is equal to the sum of weak and strong ligand concentrations (`ligW` + `ligS`, [nmol kg<sup>-1</sup>]) <br>
- $Lig_{bulk}^{K}$ is the conditional stability constant of bulk ligands and in this case, where `do_two_ligands == .false.`, is equal to $Lig_{s}^{K} \cdot 10^{-0.5}$ <br>

In the case of `do_two_ligands == .true.`, we solve for (`feIII`, $dFe_{free}$, [nmol Fe kg<sup>-1</sup>]) via the iterative method. For this approach, we know that:

$$
\begin{align}
dFe_{free} =& \quad dFe_{sFe} - \sum_{i=1}^{2} \left( \dfrac{Lig_{i}^{K} dFe_{free} [Lig_{i}]}{1 + Lig_{i}^{K} dFe_{free}} \right) \\    
\end{align}
$$

and we seek the root of the residual of free iron ($R(dFe_{free})$) defined as:

$$
\begin{align}
R(dFe_{free}) =& \quad dFe_{sFe} - \sum_{i=1}^{2} \left( \dfrac{Lig_{i}^{K} dFe_{free} [Lig_{i}]}{1 + Lig_{i}^{K} dFe_{free}} \right)  - dFe_{free} \\    
\end{align}
$$

To do so, we apply Newton-Raphson iteration using an initial guess that is the maximum of two limiting-case approximations — one accurate when ligands are unsaturated (low $dFe_{sFe}$) and one accurate when ligands are saturated (high $dFe_{sFe}$):

$$
\begin{align}
dFe_{free}^{0} =& \quad \max\left( \dfrac{dFe_{sFe}}{1 + Lig_{s}^{K}[Lig_{s}] + Lig_{w}^{K}[Lig_{w}]}, \quad dFe_{sFe} - [Lig_{s}] - [Lig_{w}] \right)
\end{align}
$$

In the low-Fe regime the first term is close to the true $dFe_{free}$ while the second term is near zero or negative; in the high-Fe regime the second term is close to the true $dFe_{free}$ while the first term is very small. Taking the maximum selects the informative approximation in each regime. Both terms underestimate $dFe_{free}$ in their respective asymptotic regimes, and after clamping to $[0, dFe_{sFe}]$ the maximum provides a valid lower bound that ensures Newton–Raphson starts from a physically meaningful value.

Whatever soluble Fe is not present as inorganic free iron is assigned to ligand-bound Fe:

$$
\begin{align}
dFe_{lig} =& \quad dFe_{sFe} - dFe_{free}
\end{align}
$$

Now that we have separated the dissolved Fe pool into its subcomponents of free, ligand-bound and colloidal Fe all in units of [nmol Fe kg<sup>-1</sup>], we solve for precipiation of nanoparticles, scavenging and coagulation of dissolved Fe, all of which remove dFe from the water column. These are the major sinks outside of phytoplankton uptake.


**Precipitation:**

Precipitation of dissolved iron specifically affects free iron when the concentration is greater than that deemed soluble. This iron is permanently lost from the water column. This only occurs when `do_colloidal_shunt == .false.`.

$$
\begin{align}
Pr_{dFe}^{\rightarrow} =& \quad \max\left(0.0, dFe_{free} - dFe_{sol}\right) \gamma_{Fe}^{nano}
\end{align}
$$

_where_ <br>
- $dFe_{free}$ is the concentration of free dissolved iron (`feIII(i,j,k)`, [nmol Fe kg<sup>-1</sup>]) <br>
- $dFe_{sol}$ is the final estimated solubility of dissolve iron in seawater (`fe3sol`, [nmol Fe kg<sup>-1</sup>]) <br>
- $\gamma_{Fe}^{nano}$ is the rate constant of nanoparticle precipitation (`knano_dfe`, [s<sup>-1</sup>]) <br>
- $Pr_{dFe}^{\rightarrow}$ is the rate of loss of dFe via nanoparticle precipitation (`feprecip(i,j,k)`, [nmol Fe kg<sup>-1</sup> s<sup>-1</sup>]) <br>

**Scavenging:**

Scavenging of dissolved iron specifically affects free iron, is accelerated by the presence of particles in the water column ([Ye et al., 2011](https://bg.copernicus.org/articles/8/2107/2011/); [Tagliabue et al., 2019](https://doi.org/10.1038/s41467-019-12775-5)) and we route this iron to sinking organic particles in WOMBAT-lite since we do not represent sinking authigenic iron particles.

$$
\begin{align}
Sc_{dFe}^{\rightarrow} =& \quad dFe_{free} \left(10^{-7} + \gamma_{dFe}^{scav} \cdot B_{particles}^{M} \right)
\end{align}
$$

_where_ <br>
- $Sc_{dFe}^{\rightarrow}$ is the rate of loss of dFe via scavenging onto all particulates (`fescaven(i,j,k)`, [nmol Fe kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $dFe_{free}$ is the concentration of free dissolved iron (`feIII(i,j,k)`, [nmol Fe kg<sup>-1</sup>]) <br>
- $\gamma_{Fe}^{scav}$ is the rate constant of scavenging (`kscav_dfe`, [(mmol m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $B_{particles}^{M}$ is the in situ mass concentration of detrital particles in the water column (`partic`, [mmol m<sup>-3</sup>]) <br>

Organic carbon-based detritus $B_{det}^{C}$ is multipled by 2 assuming that carbon represents half the mass of the particle and inorganic carbon-based particles $B_{CaCO_3}^{C}$ is multipled by 8.3 since the molecular weight of calcium carbonate is 100 g mol<sup>-1</sup>. 

$$
\begin{align}
B_{particles}^{M} =& \quad 2 \cdot B_{det}^{C} + 8.3 \cdot B_{CaCO_3}^{C}
\end{align}
$$

_where_ <br>
- $B_{det}^{C}$ is the concentration of particulate organic carbon (`det_mmolm3`, [mmol C m<sup>-3</sup>]) <br>
- $B_{CaCO_3}^{C}$ is the concentration of particulate calcium carbonate in units of carbon (`caco3_mmolm3`, [mmol C m<sup>-3</sup>]) <br>

Total scavenging ($Sc_{dFe}^{\rightarrow}$) of free iron is broken into the part that attaches to organic detritus (`fescadet(i,j,k)`, $Sc_{dFe}^{\rightarrow B_{det}^{Fe}}$, [nmol Fe kg<sup>-1</sup> s<sup>-1</sup>]):

$$
\begin{align}
Sc_{dFe}^{\rightarrow B_{det}^{Fe}} =& \quad Sc_{dFe}^{\rightarrow} \cdot \dfrac{ 2 \cdot B_{det}^{C} }{ B_{particles}^{M} } \\
\end{align}
$$


**Coagulation:**

Similar to scavenging of free iron, coagulation routes dissolved iron to sinking organic particles (since we do not represent sinking authigenic particles). However, coagulation acts on the colloidal fraction of dissolved iron ([Tagliabue et al., 2023)](https://www.nature.com/articles/s41586-023-06210-5)). The rate of coagulation of colloidal iron to sinking organic particles (`fecoag2det(i,j,k)`, $Co_{dFe}^{\rightarrow B_{det}^{Fe}}$, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]) is of the form:

$$
\begin{align}
Co_{dFe}^{\rightarrow B_{det}^{Fe}} =& \quad dFe_{col} \gamma_{dFe}^{coag} \cdot S_{coag}
\end{align}
$$

_where_ <br>
- $dFe_{col}$ is the in situ concentration of colloidal iron (`fecol(i,j,k)`, [nmol Fe kg<sup>-1</sup>]) <br>
- $\gamma_{dFe}^{coag}$ is the iron coagulation rate constant (`kcoag_dfe`, [(mmol m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $S_{coag}$ is the scaling coefficient to decelerate or accelerate coagulation (`zval`, [mmol C m<sup>-3</sup>]) <br>

The coagulation scaling coefficient is dependent on the concentrations of dissolved organic carbon, particulate organic carbon, phytoplankton biomass and the rate of mixing via

$$
\begin{align}
S_{coag} =& \quad H_{mix} \left(12 \cdot F_{coag} [DOC] + 9.05 \cdot B_{det}^{C}\right) \\
          & \quad + \left( 2.49 \cdot B_{det}^{C} + 128 \cdot F_{coag} [DOC] + 725 \cdot B_{det}^{C} \right) \\
          & \quad + \left(\gamma_{dFe}^{agg} \cdot \dfrac{\left(dFe_{col}\right)^{4}}{\left(dFe_{col}\right)^{4} + \left(K_{dFe}^{agg}\right)^{4}} \right) \\ 
F_{coag} =& \quad \dfrac{B_{phy}^{C}}{B_{phy}^{C} + 0.03}
\end{align}
$$

_where_ <br>
- $H_{mix}$ is a Heaviside step function that is equalt to 1 in the mixed layer and 0.01 beneath the mixed layer (`shear`, [dimensionless]) <br>
- $F_{coag}$ is a phytoplankton concentration dependent coagulation factor (`biof`, [dimensionless])  <br>
- $B_{phy}^{C}$ is the concentrations of phytoplankton biomass (`phy_mmolm3`, [mmol C m<sup>-3</sup>])  <br>
- $[DOC]$ is an empirical concentration of DOC and is equal to $40 + 40 \left( 1 - \min\left(L_{phy}^{N} \ , \ L_{phy}^{Fe}\right) \right)$ (`biodoc`, [mmol m<sup>-3</sup>]) <br>
- $B_{det}^{C}$ is the concentration of organic detrital particles (`det_mmolm3`, [mmol C m<sup>-3</sup>]) <br>
- $\gamma_{dFe}^{agg}$ is the colloidal iron aggregation rate constant (`kagg_col`, [s<sup>-1</sup>]) <br>
- $K_{dFe}^{agg}$ is the half-saturation coefficient for colloidal iron aggregation (`kagg_kcol`, [µmol m<sup>-3</sup>]) <br>


Together, these terms implement a biologically mediated coagulation pathway in which iron removal from the dissolved pool is tightly coupled to ecosystem state. The formulation reflects the central conclusion of [Tagliabue et al. (2023)](https://www.nature.com/articles/s41586-023-06210-5): that iron cycling is not governed solely by inorganic chemistry, but is strongly regulated by biological activity, organic matter dynamics, and particle ecology across the upper ocean. However, we note the absence of slowly sinking authigenic phases of iron in WOMBAT-lite that are represented by [Tagliabue et al. (2023)](https://www.nature.com/articles/s41586-023-06210-5). Instead, we note that colloidal coagulation passes dissolved iron directly to sinking organic detritus.

---


### 9. Mortality and remineralisation.

Mortality of phytoplankton and zooplankton are affected by both linear ($\gamma$) and quadratic ($\Gamma$) terms. Linear terms are per-capita losses associated with the costs of basal metabolism. Quadratic, and thus density-dependent losses, are associated with disease, aggregation and coagulation, viruses, infection and canabalism. None of these processes are represented explicitly within the model, so we represent them implicitly.

**Linear losses** for phytoplankton and zooplankton in [mmol C m<sup>-3</sup> s<sup>-1</sup>] are modelled as

$$
\begin{align}
\gamma_{phy}^{\rightarrow C} =& \quad \gamma_{phy}^{0^{\circ}C} \left(β_{hete}\right)^{T} B_{phy}^{C} \\
\gamma_{zoo}^{\rightarrow C} =& \quad \gamma_{zoo}^{0^{\circ}C} \left(β_{hete}\right)^{T} F_{zoo}^{\gamma} B_{zoo}^{C}
\end{align}
$$

In the above, we scale down **linear zooplankton mortality** when zooplankton biomass is small, such that

$$
\begin{align}
F_{zoo}^{\gamma} =& \quad \dfrac{B_{zoo}^{C}}{B_{zoo}^{C} + K_{zoo}^{\gamma}}
\end{align}
$$

_where_ <br>
- $\gamma_{phy}^{\rightarrow C}$ is the linear loss rate of phytoplankton biomass (`phymorl`, [mmol C m<sup>-3</sup> s<sup>-1</sup>]) <br>
- $\gamma_{zoo}^{\rightarrow C}$ is the linear loss rate of zooplankton biomass (`zoomorl`, [mmol C m<sup>-3</sup> s<sup>-1</sup>]) <br>
- $\gamma_{phy}^{0^{\circ}C}$ is the linear loss rate of phytoplankton biomass at 0ºC (`phylmor`, [s<sup>-1</sup>]) <br>
- $\gamma_{zoo}^{0^{\circ}C}$ is the linear loss rate of zooplankton biomass at 0ºC (`zoolmor`, [s<sup>-1</sup>]) <br>
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless]) <br>
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC]) <br>
- $B_{phy}^{C}$ is the concentration of phytoplankton carbon biomass (`phy_mmolm3`, [mmol C m<sup>-3</sup>]) <br>
- $B_{zoo}^{C}$ is the concentration of zooplankton carbon biomass (`zoo_mmolm3`, [mmol C m<sup>-3</sup>]) <br>
- $K_{zoo}^{\gamma}$ is the half-saturation coefficient for scaling down linear mortality losses (`zookz`, [mmolC m <sup>-3</sup>]) <br>


**Quadratic losses** of phytoplankton and zooplankton in [mmol m<sup>-3</sup> s<sup>-1</sup>] are modelled as

$$
\begin{align}
\Gamma_{phy}^{\rightarrow C} =& \quad \Gamma_{phy}^{0^{\circ}C} \left(β_{hete}\right)^{T} \left(B_{phy}^{C}\right)^{2} \\
\Gamma_{zoo}^{\rightarrow C} =& \quad \Gamma_{zoo}^{0^{\circ}C} \left(β_{hete}\right)^{T} \left(B_{zoo}^{C}\right)^{2}
\end{align}
$$

_where_ <br>
- $\Gamma_{phy}^{\rightarrow C}$ is the quadratic (density-dependent) loss rate of phytoplankton biomass (`phymorq`, [mmol C m<sup>-3</sup> s<sup>-1</sup>]) <br>
- $\Gamma_{zoo}^{\rightarrow C}$ is the quadratic (density-dependent) loss rate of zooplankton biomass (`zoomorq`, [mmol C m<sup>-3</sup> s<sup>-1</sup>]) <br>
- $\Gamma_{phy}^{0^{\circ}C}$ is the quadratic (density-dependent) loss rate of phytoplankton biomass at 0ºC (`phyqmor`, [(mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\Gamma_{zoo}^{0^{\circ}C}$ is the quadratic (density-dependent) loss rate of zooplankton biomass at 0ºC (`zooqmor`, [(mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless]) <br>
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC]) <br>
- $B_{phy}^{C}$ is the concentration of phytoplankton carbon biomass (`phy_mmolm3`, [mmol C m<sup>-3</sup>]) <br>
- $B_{zoo}^{C}$ is the concentration of zooplankton carbon biomass (`zoo_mmolm3`, [mmol C m<sup>-3</sup>]) <br>


**Remineralisation** of detritus is only affected by a quadratic, density-dependent loss term,

$$
\begin{align}
\Gamma_{det}^{\rightarrow C} =& \quad \Gamma_{det}^{0^{\circ}C} \left(β_{hete}\right)^{T} \left(B_{det}^{C}\right)^{2}
\end{align}
$$

_where_ <br>
- $\Gamma_{det}^{\rightarrow C}$ is the quadratic (density-dependent) loss rate of particulate organic detritus (`detremi`, [mmol C m<sup>-3</sup> s<sup>-1</sup>]) <br>
- $\Gamma_{det}^{0^{\circ}C}$ is the quadratic (density-dependent) loss rate of particulate organic detritus at 0ºC (`detlrem`, [(mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless]) <br>
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC]) <br>
- $B_{det}^{C}$ is the concentration of particulate organic detritus (`det_mmolm3`, [mmol C m<sup>-3</sup>]) <br>

since hydrolyzation of organic detritus is performed by an heterotrophic bacterial population that is not explicitly resolved in the model and their acitivity is density-dependent.

---


### 10. Zooplankton grazing, egestion, excretion and assimilation.

**Grazing by zooplankton** (`g_zoo`, $g_{zoo}$, [s<sup>-1</sup>]) is computed using a Holling Type III functional response [Holling, 1959](https://doi.org/10.4039/Ent91385-7):

$$
\begin{align}
g_{zoo} =& \quad \dfrac{\mu_{zoo}^{max} \left(β_{hete}\right)^{T} \cdot L_{zoo}^{O_2} \varepsilon \left(B_{prey}^{C}\right)^{2}}{\mu_{zoo}^{max} \left(β_{hete}\right)^{T} + \varepsilon \left(B_{prey}^{C}\right)^{2}}
\end{align}
$$

_where_ <br>
- $\mu_{zoo}^{max}$ is the maximum rate of zooplankton grazing at 0ºC (`zoogmax`, [s<sup>-1</sup>]) <br>
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless]) <br>
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC]) <br>
- $L_{zoo}^{O_2}$ is a limiter of grazing in low oxygen conditions (`zoo_o2lim`, [dimensionless]) <br>
- $B_{prey}^{C}$ is the concentration of prey biomass (`zooprey`, [mmol C m<sup>-3</sup>]) <br>
- $\varepsilon$ is the prey capture rate coefficient (`zooeps(i,j,k)`, [(mmol C m<sup>-3</sup>)<sup>-2</sup>]) <br>

We apply an oxygen limitation to grazing at low oxygen concentrations (`zoo_o2lim`, $L_{zoo}^{O_2}$, [dimensionless]) based on the review of [Medina et al. (2017)](https://doi.org/10.3389/fmars.2017.00105) who found that from 5 µM to undetectable concentrations of oxygen, protist consumption of a bacterial population decreased from 28% to 13% of the population. Although they still see some consumption at undetectable oxygen concentrations, we scale grazing down to zero at an oxygen concentration of zero to avoid unrealistic negative oxygen concentrationsdue to grazing at zero oxygen:

$$
\begin{align}
L_{zoo}^{O_2} =& \quad \left(1 - e^{\left(-\dfrac{O_2}{10}\right)}\right)
\end{align}
$$

Total grazing of biomass by zooplankton ([mol C kg<sup>-1</sup> day<sup>-1</sup>]) is therefore

$$
\begin{align}
g_{zoo}^{\leftarrow C} =& \quad g_{zoo} B_{zoo}^{C} 
\end{align}
$$

_where_ <br>
- $g_{zoo}$ is the total specific rate of grazing of zooplankton (`g_zoo`, [s<sup>-1</sup>]) <br>
- $B_{zoo}^{C}$ is the in situ concentration of zooplankton carbon biomass (`p_zoo(i,j,k,tau)`, [mol C kg<sup>-1</sup>]) <br>

This formulation suppresses grazing at very low prey biomass ($B_{prey}^{C}$) due to reduced encounter and clearance rates, accelerates grazing at intermediate prey biomass as zooplankton effectively learn and switch to available prey, and saturates at high prey biomass due to handling-time limitation ([Gentleman and Neuheimer, 2008](https://doi.org/10.1093/plankt/fbn078); Rohr et al., [2022](https://doi.org/10.1016/j.pocean.2022.102878), [2024](https://doi.org/10.1029/2023GL107732)). This choice increases ecosystem stability and prolongs phytoplankton blooms relative to a Type II formulation.

The application of $g_{zoo}^{max} (β_{hete})^{T}$ in both the numerator and denominator makes this grazing formula unique [(Rohr et al., 2023)](https://www.nature.com/articles/s43247-023-00871-w) and equivalent to a disk formulation, rather than a Michaelis–Menten formulation [(Rohr et al., 2022)](https://doi.org/10.1016/j.pocean.2022.102878). Practically, this amplifies grazing in warmer climes, but to a lesser extent than other formulations that apply the temperature amplification ($(β_{hete})^{T}$) only in the numerator [(Rohr et al., 2023)](https://www.nature.com/articles/s43247-023-00871-w). This dampens the effect that variations in temperature have on grazing activity, amplifying the effect of $\varepsilon$ and aligning with observations that the ratio of grazing to phytoplankton growth varies little between tropical and polar climes [(Calbet and Landry, 2004)](https://doi.org/10.4319/lo.2004.49.1.0051). Theoretically, this assumes some evolutionary adaptation to account for the physiological effects of temperature across environmental niches, such that the efficiency of prey capture and handling becomes more important to grazers than metabolic constraints due to temperature.

The total prey biomass available to zooplankton is defined as a preference-weighted sum of phytoplankton and detritus that is normalized to reflect explicit dietary fractions ([Gentleman et al., 2003](https://doi.org/10.1016/j.dsr2.2003.07.001)):

$$
\begin{align}
B_{prey}^{C} =& \quad \phi_{zoo}^{phy} B_{phy}^{C} + \phi_{zoo}^{det} B_{det}^{C}
\end{align}
$$

_where_ <br>
- $B_{phy}^{C}$ is the concentration of phytoplankton biomass (`phy_mmolm3`, [mmol C m<sup>-3</sup>]) <br>
- $B_{det}^{C}$ is the concentration of particulate organic detritus (`det_mmolm3`, [mmol C m<sup>-3</sup>]) <br>
- $\phi_{zoo}^{phy}$ is the relative preference of zooplankton grazing on phytoplankton (`zooprefphy(i,j,k)`, [dimensionless]) <br>
- $\phi_{zoo}^{det}$ is the relative preference of zooplankton grazing on particulate detritus (`zooprefdet(i,j,k)`, [dimensionless]) <br>

and where:

$$
\begin{align}
\phi_{zoo}^{phy} + \phi_{zoo}^{det} =& \quad 1
\end{align}
$$

The normalized prey preferences (i.e., dietary fractions) are further modified by prey switching prior to computation of total prey biomass ([Gentleman et al., 2003](https://doi.org/10.1016/j.dsr2.2003.07.001)) such that

$$
\begin{align}
\phi_{zoo}^{phy} =& \quad \left( \phi_{zoo}^{phy} B_{phy}^{C} \right)^{s_{zoo}} \\
\phi_{zoo}^{det} =& \quad \left( \phi_{zoo}^{det} B_{det}^{C} \right)^{s_{zoo}}
\end{align}
$$

_where_ <br>
- $\phi_{zoo}^{phy}$ and $\phi_{zoo}^{det}$ are the relative prey preference of zooplankton for phytoplankton and detritus (`zooprefphy(i,j,k)`; `zooprefdet(i,j,k)`, [dimensionless]) <br>
- $B_{phy}^{C}$ and $B_{det}^{C}$ are the concentrations of phytoplankton and detritus in carbon biomass (`phy_mmolm3`; `det_mmolm3`, [mmol C m<sup>-3</sup>])<br>
- $s_{zoo}$ is the prey-switching exponent of zooplankton (`zoopreyswitch`) <br>

When $s_{zoo} < 1$, zooplankton feed equally across all prey items irrespective of availability  <br>
When $s_{zoo} = 1$, zooplankton feed according to pre-defined dietary fractions  <br>
When $s_{zoo} > 1$, zooplankton exhibit prey-switching and feed disproportionately on most abundant prey  <br>

Again, prey preferences are normalized to ensure that $\phi_{zoo}^{phy} + \phi_{zoo}^{det} = 1$.

The prey capture rate coefficient, $\varepsilon$ (`zooeps(i,j,k)`, $\varepsilon$, [(mmol C m<sup>-3</sup>)<sup>-2</sup>]), is allowed to vary as a function of prey biomass, following the prey-dependent behaviour described by [Rohr et al. (2024)](doi.org/10.1029/2023GL107732). This reflects a transition from microzooplankton-like feeding with higher prey capture rate coefficients at low prey biomass to mesozooplankton-like feeding with lower prey capture rate coefficients at high prey biomass.

A prey-dependent scaling factor (`g_peffect`, $F_{prey}$, [dimensionless]) is defined as

$$
\begin{align}
F_{prey} =& \quad e^{\left(-B_{prey}^{C} \varepsilon_{shift} \right)}
\end{align}
$$

and the effective capture rate coefficient, (`zooeps`, $\varepsilon$, [(mmol C m<sup>-3</sup>)<sup>-2</sup>]) is then computed as

$$
\begin{align}
\varepsilon =& \quad \varepsilon_{\min} + \left(\varepsilon_{\max} - \varepsilon_{\min}\right) F_{prey}
\end{align}
$$

_where_ <br>
- $\varepsilon_{shift}$ is the rate at which the prey capture efficiency transitions from micro- to meso-zooplankton (`zooepsrat`, [(mmol C m<sup>-3</sup>)<sup>-1</sup>]) <br>
- $\varepsilon_{\min}$ is the prey capture rate coefficient of a zooplankton community dominated by meso-zooplankton (`zooepsmin`, [(mmol C m<sup>-3</sup>)<sup>-2</sup>]) <br>
- $\varepsilon_{\max}$ is the prey capture rate coefficient of a zooplankton community dominated by micro-zooplankton (`zooepsmax`, [(mmol C m<sup>-3</sup>)<sup>-2</sup>]) <br>

At low prey biomass, $\varepsilon \rightarrow \varepsilon_{\max}$, enhancing grazing efficiency. At high prey biomass, $\varepsilon \rightarrow \varepsilon_{\min}$, reducing capture efficiency as handling time and feeding mode are more ineffective on average in a community with relatively more mesozooplankton.

Total grazing of prey can also be expressed as the sum of individual prey type consumption:

$$
\begin{align}
g_{zoo}^{\leftarrow C} =& \quad g_{zoo}^{\leftarrow B_{phy}^{C}} + g_{zoo}^{\leftarrow B_{det}^{C}}
\end{align}
$$

_where_ <br>
- $g_{zoo}^{\leftarrow B_{phy}^{C}} = g_{zoo} \dfrac{\phi_{zoo}^{phy} B_{phy}^{C}}{B_{prey}^{C}}$ and is the proportion of zooplankton grazing of phytoplankton (`zoograzphy`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{zoo}^{\leftarrow B_{det}^{C}} = g_{zoo} \dfrac{\phi_{zoo}^{det} B_{det}^{C}}{B_{prey}^{C}}$ and is the proportion of zooplankton grazing of detritus (`zoograzdet`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>

In this formulation, consumption of each prey item $i$ in [mol C kg<sup>-1</sup>] can also be expressed as:

$$
\begin{align}
g_{zoo}^{\leftarrow B_{i}^{C}} =& \quad g_{zoo} B_{zoo}^{C} \cdot \varepsilon \dfrac{\left(\phi_{zoo}^{i} B_{i}^{C} \right)^{2}}{\sum_{i} \left(\phi_{zoo}^{i} B_{i}^{C} \right)^{2}}
\end{align}
$$


**Zooplankton egestion, excretion and assimilation** are then calculated assuming static assimilation coefficients. Grazed biomass is routed to either egestion or ingestion via an ingestion coefficient ($\lambda^{C}$, [mol C (mol C)<sup>-1</sup>]), with the egested fraction being equal to $1.0 - \lambda^{C}$. The biomass that is ingested is then split between assimilation and excretion based on an assimilation coefficient ($\eta^{C}$, [mol C (mol C)<sup>-1</sup>]) with the excreted fraction being equal to $1.0 - \eta^{C}$. Egestion ($E$), excretion ($X$) and assimilation ($A$) of organic carbon due to grazing of prey type $i$ by zooplankton are calculated as:

$$
\begin{align}
E_{zoo}^{\leftarrow B_{i}^{C}} =& \quad g_{zoo}^{\leftarrow B_{i}^{C}} \left(1 - \lambda^{C} \right) \\
X_{zoo}^{\leftarrow B_{i}^{C}} =& \quad g_{zoo}^{\leftarrow B_{i}^{C}} \lambda^{C} \left(1 - \eta^{C} \right) \\
A_{zoo}^{\leftarrow B_{i}^{C}} =& \quad g_{zoo}^{\leftarrow B_{i}^{C}} \lambda^{C} \eta^{C}
\end{align}
$$

_where_ <br>
- $E_{zoo}^{\leftarrow B_{i}^{C}}$ is the rate of egestion of carbon biomass by zooplankton feeding on prey type $i$ ([mol C kg<sup>-1</sup>]) <br>
- $X_{zoo}^{\leftarrow B_{i}^{C}}$ is the rate of excretion of carbon biomass by zooplankton feeding on prey type $i$ ([mol C kg<sup>-1</sup>]) <br>
- $A_{zoo}^{\leftarrow B_{i}^{C}}$ is the rate of assimilation of carbon biomass by zooplankton feeding on prey type $i$ ([mol C kg<sup>-1</sup>]) <br>
- $g_{zoo}^{\leftarrow B_{i}^{C}}$ is the grazing rate of zooplankton on prey type $i$ ([mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\lambda^{C}$ is the fraction of prey carbon biomass that is ingested by zooplankton (`zooCingest`, [mol C (mol C)<sup>-1</sup>]) <br>
- $\eta^{C}$ is the fraction of ingested prey carbon biomass that is assimilated by zooplankton (`zooCassim`, [mol C (mol C)<sup>-1</sup>]) <br>

Total egestion, excretion and assimilation or carbon are therefore:

$$
\begin{align}
E_{zoo}^{\leftarrow C} =& \quad g_{zoo}^{\leftarrow C} \left(1 - \lambda^{C} \right) \\
X_{zoo}^{\leftarrow C} =& \quad g_{zoo}^{\leftarrow C} \lambda^{C} \left(1 - \eta^{C} \right) \\
A_{zoo}^{\leftarrow C} =& \quad g_{zoo}^{\leftarrow C} \lambda^{C} \eta^{C}
\end{align}
$$

Because we track both carbon and iron through the ecosystem components, we assign unique ingestion and assimilation coefficients to carbon and iron. This separation of ingestion and assimilation coefficients for iron and carbon follows [Le Mézo & Galbraith (2021)](https://doi.org/10.1002/lno.11597). For iron, we apply unique ingestion ($\lambda^{Fe}$, [mol Fe (mol Fe)<sup>-1</sup>]) and assimilation coefficients ($\eta^{Fe}$, [mol Fe (mol Fe)<sup>-1</sup>]). [Le Mézo & Galbraith (2021)](https://doi.org/10.1002/lno.11597) show that if $\lambda^{Fe} << \lambda^{C}$ then egestion is enriched in Fe:C, and it follows that $\eta^{Fe} >> \eta^{C}$ so that zooplankton can absorb sufficient iron from their prey. Consequently:

$$
\begin{align}
E_{zoo}^{\leftarrow B_{i}^{Fe}} =& \quad g_{zoo}^{\leftarrow B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \left(1 - \lambda^{Fe} \right) \\
X_{zoo}^{\leftarrow B_{i}^{Fe}} =& \quad g_{zoo}^{\leftarrow B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \lambda^{Fe} \left(1 - \eta^{Fe} \right) \\
A_{zoo}^{\leftarrow B_{i}^{Fe}} =& \quad g_{zoo}^{\leftarrow B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \lambda^{Fe} \eta^{Fe}
\end{align}
$$

_where_ <br>
- $E_{zoo}^{\leftarrow B_{i}^{Fe}}$ is the rate of egestion of iron biomass by zooplankton feeding on prey type $i$ ([mol Fe kg<sup>-1</sup>]) <br>
- $X_{zoo}^{\leftarrow B_{i}^{Fe}}$ is the rate of excretion of iron biomass by zooplankton feeding on prey type $i$ ([mol Fe kg<sup>-1</sup>]) <br>
- $A_{zoo}^{\leftarrow B_{i}^{Fe}}$ is the rate of assimilation of iron biomass by zooplankton feeding on prey type $i$ ([mol Fe kg<sup>-1</sup>]) <br>
- $g_{zoo}^{\leftarrow B_{i}^{C}}$ is the grazing rate of zooplankton on prey type $i$ ([mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\dfrac{B_{i}^{Fe}}{B_{i}^{C}}$ is the Fe:C ratio of prey type $i$ ([mol Fe (mol C)<sup>-1</sup>]) <br>
- $\lambda^{Fe}$ is the fraction of prey iron biomass that is ingested by zooplankton (`zooFeingest`, [mol Fe (mol Fe)<sup>-1</sup>]) <br>
- $\eta^{Fe}$ is the fraction of ingested prey iron biomass that is assimilated by zooplankton (`zooFeassim`, [mol Fe (mol Fe)<sup>-1</sup>]) <br>

Total egestion, excretion and assimilation or iron are therefore:

$$
\begin{align}
E_{zoo}^{\leftarrow Fe} =& \quad \sum_{i} \left( g_{zoo}^{\leftarrow B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \right) \cdot \left(1 - \lambda^{Fe} \right) \\
X_{zoo}^{\leftarrow Fe} =& \quad \sum_{i} \left( g_{zoo}^{\leftarrow B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \right) \cdot \lambda^{Fe} \left(1 - \eta^{Fe} \right) \\
A_{zoo}^{\leftarrow Fe} =& \quad \sum_{i} \left( g_{zoo}^{\leftarrow B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \right) \cdot \lambda^{Fe} \eta^{Fe}
\end{align}
$$


---


### 11. CaCO3 calculations.

**Dynamic $CaCO_3$ production and dissolution**

When $CaCO_3$ dynamics are enabled (`do_caco3_dynamics = .true.`), the model computes both particulate inorganic carbon production (via the PIC:POC ratio) and $CaCO_3$ dissolution rates as functions of carbonate chemistry, temperature, and organic matter availability.


**Production** of $CaCO_3$ in WOMBAT-lite comes from three sources: (1) density-dependent phytoplankton mortality ($P_{CaCO_3}^{\Gamma phy}$), (2) density-dependent zooplankton mortality ($P_{CaCO_3}^{\Gamma zoo}$) and (3) zooplankton egestion of grazed phytoplankton ($P_{CaCO_3}^{g_{zoo}^{\leftarrow phy}}$). Each term is multipled by the particulate inorganic to organic carbon production ratio (`pic2poc(i,j,k)`, $PIC:POC$, [mol C (mol C)<sup>-1</sup>]) to achieve a rate of $CaCO_3$ production [mol C kg<sup>-1</sup> s<sup>-1</sup>]:

$$
\begin{align}
(1) & \qquad P_{CaCO_3}^{\Gamma_{phy}^{\rightarrow C}} = \quad \Gamma_{phy}^{\rightarrow C} \cdot PIC:POC \\
(2) & \qquad P_{CaCO_3}^{\Gamma_{phy}^{\rightarrow C}} = \quad \Gamma_{zoo}^{\rightarrow C} \cdot PIC:POC \\
(3) & \qquad P_{CaCO_3}^{g_{zoo}^{\leftarrow B_{phy}^{C}}} = \quad g_{zoo}^{\leftarrow B_{phy}^{C}} \cdot PIC:POC \cdot \left(1 - F_{gut} \right)
\end{align}
$$

_where_ <br>
- $\Gamma_{phy}^{\rightarrow C}$ is the quadratic (density-dependent) loss rate of phytoplankton biomass (`phymorq`, [mmol C m<sup>-3</sup> s<sup>-1</sup>]) <br>
- $\Gamma_{zoo}^{\rightarrow C}$ is the quadratic (density-dependent) loss rate of zooplankton biomass (`zoomorq`, [mmol C m<sup>-3</sup> s<sup>-1</sup>]) <br>
- $g_{zoo}^{\leftarrow B_{phy}^{C}}$ is the phytoplankton grazed by zooplankton (`zoograzphy`, [mol C kg<sup>-1</sup> s<sup>-1</sup>])
- $PIC:POC$ is the particulate inorganic (i.e., $CaCO_3$) to particulate organic ratio (`pic2poc(i,j,k)`, [mol C (mol C)<sup>-1</sup>]) <br>
- $F_{gut}$ is the fraction of $CaCO_3$ that is dissolved within zooplankton guts (`fgutdiss`, [mol C (mol C)<sup>-1</sup>]) <br>

In the above, the $PIC:POC$ ratio is formulated as:

$$
\begin{align}
PIC:POC =& \quad \min \left( 0.3,  \left( f_{inorg} + 10^{-3 + 4.31 \times 10^{-6} \left( \dfrac{[HCO_3^-]}{[H^+]} \right)} \right) F_T \right)
\end{align}
$$

_where_ <br>
- $f_{inorg}$ is the background PIC:POC ratio (`f_inorg`, [mol C (mol C)<sup>-1</sup>]) <br>
- $[HCO_{3}^{-}]$ is the concentration of bicarbonate ions (`hco3`, [mol kg<sup>-1</sup>]) <br>
- $[H^{+}]$ is concentration of free hydrogen ions (`htotal(i,j,k)`, [µmol kg<sup>-1</sup>]) <br>
- $F_{T}$ is a temperature-dependent suppression term and if defined by $F_T = 0.55 + 0.45 \cdot \tanh \left(T - 4\right)$ <br>

This formulation of $PIC:POC$ is therefore a function of the substrate–inhibitor ratio between bicarbonate and free hydrogen ions (`hco3 / htotal(i,j,k)`, $\dfrac{[HCO_3^-]}{[H^+]}$, [mol µmol<sup>-1</sup>]), following [Lehmann & Bach (2025)](https://www.nature.com/articles/s41561-025-01644-0). This reflects the sensitivity of calcification to carbonate system speciation, which is nonlinearly enhanced with increasing $[HCO_3^-]/[H^+]$. Moreover, the $F_{T}$ term strongly reduces $CaCO_3$ production in cold waters, enforcing near-zero calcification below approximately 3 °C consistent with observations of _Emiliania huxleyi_ growth limits in polar environments [(Fielding, 2013)](https://doi.org/10.4319/lo.2013.58.2.0663). Finally, we also cap the $PIC:POC$ ratio at an upper bound of 0.3 to prevent unrealistically high inorganic carbon production and accord with the highest measured ratios in the ocean.


**Dissolution** of $CaCO_3$ is computed as the sum of four contributions: 

(1) undersaturation-driven dissolution of calcite (`caldiss(i,j,k)`, $D_{CaCO_3}^{\Omega_{cal}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
(2) undersaturation-driven dissolution of aragonite (`aradiss(i,j,k)`, $D_{CaCO_3}^{\Omega_{ara}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
(3) biologically-mediated dissolution associated with degredation of detrital organic matter (`pocdiss(i,j,k)`, $D_{CaCO_3}^{\Gamma_{det}^{\rightarrow C}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
(4) dissolution within zooplankton during their digestion of detrital aggregates (`zoodiss(i,j,k)`, $D_{CaCO_3}^{g_{zoo}^{\leftarrow B_{det}^{C}}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>

Total $CaCO_3$ dissolution is:

$$
\begin{align}
D_{CaCO_3} =& \quad D_{CaCO_3}^{\Omega_{cal}} + D_{CaCO_3}^{\Omega_{ara}} + D_{CaCO_3}^{\Gamma_{det}^{\rightarrow C}} + D_{CaCO_3}^{g_{zoo}^{\leftarrow B_{det}^{C}}} 
\end{align}
$$

The first three terms follow [Kwon et al. (2024)](https://www.science.org/doi/full/10.1126/sciadv.adl0779):

$$
\begin{align}
(1) & \qquad D_{CaCO_3}^{\Omega_{cal}} = \quad d_{CaCO_3}^{\Omega_{cal}} \max\left(0,  1 - \Omega_{cal}\right)^{2.2} B_{CaCO_3}^{C} \\ 
(2) & \qquad D_{CaCO_3}^{\Omega_{ara}} = \quad d_{CaCO_3}^{\Omega_{ara}} \max\left(0,  1 - \Omega_{ara}\right)^{1.5} B_{CaCO_3}^{C} \\
(3) & \qquad D_{CaCO_3}^{\Gamma_{det}^{\rightarrow C}} = \quad d_{CaCO_3}^{\Gamma_{det}} \Gamma_{det}^{\rightarrow C} B_{CaCO_3}^{C}
\end{align}
$$

_where_ <br>
- $\Omega_{cal}$ is the saturation state of calcite (`omega_cal(i,j,k)`, [dimenionless]) <br>
- $\Omega_{ara}$ is the saturation state of aragonite (`omega_ara(i,j,k)`, [dimenionless]) <br>
- $d_{CaCO_3}^{\Omega_{cal}}$ is the reference dissolution rate constant for calcite (`disscal`, [s<sup>-1</sup>])  <br>
- $d_{CaCO_3}^{\Omega_{ara}}$ is the reference dissolution rate constant for aragonite (`dissara`, [s<sup>-1</sup>])  <br>
- $d_{CaCO_3}^{\Gamma_{det}}$ is the reference dissolution rate constant per unit of small detrital organic carbon remineralised (`dissdet`, [(mmol C m<sup>-3</sup>)<sup>-1</sup>])  <br>
- $\Gamma_{det}^{\rightarrow C}$ is the in situ remineralisation rate of small detrital organic carbon (`detremi(i,j,k)`, [mmol C m<sup>-3</sup> s<sup>-1</sup>]) <br>
- $B_{CaCO_3}^{C}$ is the in situ concentration of $CaCO_3$ in carbon units (`p_caco3(i,j,k,tau)`, [mol C kg<sup>-1</sup>]) <br>

For $D_{CaCO_3}^{\Omega_{cal}}$ and $D_{CaCO_3}^{\Omega_{ara}}$, dissolution is activated only under undersaturated conditions ($\Omega_{cal} < 1$; $\Omega_{ara} < 1$) and increases nonlinearly with increasing undersaturation. In contrast, $D_{CaCO_3}^{\Gamma_{det}^{\rightarrow C}}$ represents shallow water dissolution due to reducing microenvironments. In this scenario, $\Omega_{cal}$ and $\Omega_{ara}$ tend to be > 1 ([Sulpis et al., 2021](https://doi.org/10.1038/s41561-021-00743-y)) but dissolution nonetheless occurs in microenvironments enriched in $CO_{2}^{*}$ due to heterotrophic activity ([Borer et al., 2026](https://doi.org/10.1073/pnas.2510025123)).

The fourth term (`zoodiss(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]), which represents dissolution of $CaCO_3$ during zooplankton consumption of particulates and is calculated as:

$$
\begin{align}
(4) & \qquad D_{CaCO_3}^{g_{zoo}^{\leftarrow B_{det}^{\leftarrow C}}} = \quad g_{zoo}^{\leftarrow B_{det}^{C}} F_{gut} \dfrac{B_{CaCO_3}^{C}}{B_{det}^{C}} \\
\end{align}
$$

_where_ <br>
- $g_{zoo}^{\leftarrow B_{det}^{C}}$ is the grazing rate of particulate detritus by zooplankton (`zoograzdet(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $F_{gut}$ is the fraction of $CaCO_3$ that is dissolved within zooplankton guts (`fgutdiss`, [mol C (mol C)<sup>-1</sup>]) <br>
- $\dfrac{B_{CaCO_3}^{C}}{B_{det}^{C}}$ is the in situ ratio of $CaCO_3$ to organic carbon detritus (`caco3_mmolm3/det_mmolm3`, [mol C (mol C)<sup>-1</sup>]) <br>

Here we note that the processing of $CaCO_3$ by zooplankton grazing is treated differently to processing of organic carbon. For organic carbon, we route the biomass between zooplankton biomass (assimilation), inorganic nutrients (excretion) and particulate detritus (egestion). For $CaCO_3$ consumption by zooplankton the $CaCO_3$ is not assimilated since it does not contain nitrogen or other key elements for biosynthesis, and so is only routed between excretion to DIC and alkalinity or goes undissolved and remains $CaCO_3$ that sinks through the water column. This is supported by the fact that micro- and meso-zooplankton may dissolve 92±7% and 38-73% of coccolithophore calcite during feeding, respectively ([Smith et al., 2024](https://www.science.org/doi/10.1126/sciadv.adr5453); [White et al., 2018](https://www.nature.com/articles/s41598-018-28073-x); [Harris 1994](https://link.springer.com/article/10.1007/BF00347540)), and that the remainder is excreted and not assimilated ([Mayers et al., 2020](https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2020.569896/full)).


**Static $CaCO_3$ production and dissolution**

When $CaCO_3$ dynamics are disabled (`do_caco3_dynamics = .false.`), the model uses a static PIC:POC ratio (`f_inorg + 0.025`, [mol C (mol C)<sup>-1</sup>]) and a constant linear $CaCO_3$ dissolution rate (`caco3lrem`, $d_{cal}$, [s<sup>-1</sup>]). These are set as input parameters to the model.

---


### 12. Implicit nitrogen fixation.

Because we do not consider diazotrophs as an explicit phytoplankton functional type, we represent the fixation of nitrogen implicitly using a simple parameterization dependent on temperature, nutrient and light availability when `do_nitrogen_fixation == .true.`. The equation for new nitrogen (specifically NO<sub>3</sub>) added via diazotrophy is:

$$
\begin{align}
\mu_{diazo}^{\rightarrow NO_3} =& \quad \mu_{diazo}^{max} \left(1 - L_{phy}^{N} \right) \\
                                & \min\left(L_{diazo}^{Fe}, L_{diazo}^{PAR}\right) R_{diazo}^{N:C} \cdot 1 \times 10^{-6}
\end{align}
$$

_where_ <br>
- $\mu_{diazo}^{max}$ is the temperature-dependent maximum growth rate of diazotrophs (`tri_mumax(i,j,k)`, [s<sup>-1</sup>]) <br>
- $L_{phy}^{N}$ is the limitation term of nano-phytoplankton growth on nitrogen (`phy_lnit(i,j,k)`, [dimensionless])  <br>
- $L_{diazo}^{Fe}$ is the limitation term of diazotroph growth on iron (`tri_lfer(i,j,k)`, [dimensionless])  <br>
- $L_{diazo}^{PAR}$ is the limitation term of diazotroph growth on light (`tri_lpar(i,j,k)`, [dimensionless])  <br>
- $R_{diazo}^{N:C}$ is the ratio of N:C within diazotrophic biomass (`trin2c`, [mol N (mol C)<sup>-1</sup>]) <br>
- $1 \times 10^{-6}$ is a conversion factor to mol kg<sup>-1</sup> <br>

The temperature-dependent maximum growth rate ($\mu_{diazo}^{max}$) is taken directly from [Wrightson et al. (2022)]( https://doi.org/10.1111/gcb.16399) who based their formulation on the work of [Jiang et al. (2018)](https://doi.org/10.1038/s41558-018-021):

$$
\begin{align}
\mu_{diazo}^{max} =& \quad \left( -0.000399(T)^{3} + 0.02685(T)^{2} - 0.555T + 3.633 \right) \dfrac{1}{86400}
\end{align}
$$

_where_ <br>
- $T$ is in situ water temperature (`Temp(i,j,k)`, [ºC]) and we only consider $T > 15.8$ºC <br>
- $\dfrac{1}{86400}$ converts their formula from units of [day<sup>-1</sup>] to [s<sup>-1</sup>] <br>

The iron and light limitation terms are as follows:

$$
\begin{align}
L_{diazo}^{Fe} =& \quad \dfrac{dFe}{dFe + K_{diazo}^{Fe}} \\
L_{diazo}^{PAR} =& \quad 1 - e^{- \alpha_{diazo} PAR}
\end{align}
$$

_where_ <br>
- $dFe$ is the in situ concentration of dissolved iron (`fe_umolm3`, [nmol Fe kg<sup>-1</sup>]) <br>
- $K_{diazo}^{Fe}$ is the half-saturation coefficient for uptake of dissolved iron by diazotrophs (`trikf`, [nmol Fe kg<sup>-1</sup>]) <br>
- $\alpha_{diazo}$ is the chlorophyll-adjusted slope of the photosynthesis-irradience curve of diazotrophs (`alphabio_tri * trichlc`, [(W m<sup>-2</sup>)<sup>-1</sup>]) <br>
- $PAR$ is the downwelling photosynthetically available radiation (`radbio`, [W m<sup>-2</sup>]) <br>

---

### 13. Tracer tendencies.

**Nitrate** (`p_no3(i,j,k,tau)`, $NO_3$, [mol N kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta NO_3}{\Delta t} =& \quad \bigg(\Gamma_{det}^{\rightarrow C} + \gamma_{zoo}^{\rightarrow C} \\
                               & \quad + \gamma_{phy}^{\rightarrow C} + X_{zoo}^{\leftarrow C} \\
                               & \quad - \mu_{phy}^{\leftarrow C} \bigg) \cdot \dfrac{16}{122} + \mu_{diazo}^{\rightarrow NO_{3}}
\end{align}
$$

**Oxygen** (`p_o2(i,j,k,tau)`, $O_2$, [mol O<sub>2</sub> kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta O_2}{\Delta t} =& \quad \bigg( \mu_{phy}^{\leftarrow C} \\
                              & \quad - \left( \Gamma_{det}^{\rightarrow C} + \gamma_{zoo}^{\rightarrow C} + \gamma_{phy}^{\rightarrow C} + X_{zoo}^{\leftarrow C}\right) \bigg) \cdot \dfrac{172}{122}
\end{align}
$$

**Dissolved iron** (`p_fe(i,j,k,tau)`, $dFe$, [mol Fe kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta dFe}{\Delta t} =& \quad \Gamma_{det}^{\rightarrow C} Q_{det}^{Fe:C} 
                                    + \gamma_{phy}^{\rightarrow C} Q_{phy}^{Fe:C}
                                    + \gamma_{zoo}^{\rightarrow C} Q_{zoo}^{Fe:C} \\
                              & \quad + X_{zoo}^{\leftarrow Fe} - \mu_{phy}^{\leftarrow dFe} \\
                              & \quad - \left( Pr_{dFe}^{\rightarrow} - Sc_{dFe}^{\rightarrow} - Co_{dFe}^{\rightarrow B_{det}^{Fe}} \right)
\end{align}
$$

**Phytoplankton** (`p_phy(i,j,k,tau)`, $B_{phy}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{phy}^{C}}{\Delta t} =& \quad \mu_{phy}^{\leftarrow C} \\
                                        & - \left( \Gamma_{phy}^{\rightarrow C} + \gamma_{phy}^{\rightarrow C} + g_{zoo}^{\leftarrow B_{phy}^{C}} \right)
\end{align}
$$

**Phytoplankton chlorophyll** (`p_pchl(i,j,k,tau)`, $B_{phy}^{Chl}$, [mol Chl kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{phy}^{Chl}}{\Delta t} =& \quad \mu_{phy}^{\leftarrow Chl} \\
                                        & - \left( \Gamma_{phy}^{\rightarrow C} + \gamma_{phy}^{\rightarrow C} + g_{zoo}^{\leftarrow B_{phy}^{C}} \right) \cdot Q_{phy}^{Chl:C}
\end{align}
$$

**Phytoplankton iron** (`p_phyfe(i,j,k,tau)`, $B_{phy}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{phy}^{Fe}}{\Delta t} =& \quad \mu_{phy}^{\leftarrow dFe} \\
                                       & - \left( \Gamma_{phy}^{\rightarrow C} + \gamma_{phy}^{\rightarrow C} + g_{zoo}^{\leftarrow B_{phy}^{C}} \right) \cdot Q_{phy}^{Fe:C}
\end{align}
$$

**Zoooplankton** (`p_zoo(i,j,k,tau)`, $B_{zoo}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{zoo}^{C}}{\Delta t} =& \quad A_{zoo}^{\leftarrow C} - \Gamma_{zoo}^{\rightarrow C} - \gamma_{zoo}^{\rightarrow C}
\end{align}
$$

**Zoooplankton iron** (`p_zoofe(i,j,k,tau)`, $B_{zoo}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{zoo}^{Fe}}{\Delta t} =& \quad A_{zoo}^{\leftarrow Fe} 
                                       - \left( \Gamma_{zoo}^{\rightarrow C} + \gamma_{zoo}^{\rightarrow C} \right) \cdot Q_{zoo}^{Fe:C}
\end{align}
$$

**Detritus** (`p_det(i,j,k,tau)`, $B_{det}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{det}^{C}}{\Delta t} =& \quad E_{zoo}^{\leftarrow C} + \Gamma_{phy}^{\rightarrow C} + \Gamma_{zoo}^{\rightarrow C} \\
                                      & \quad - \left( \Gamma_{det}^{\rightarrow C} + g_{zoo}^{\leftarrow B_{det}^{C}} \right)
\end{align}
$$

**Detritus iron** (`p_detfe(i,j,k,tau)`, $B_{det}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{det}^{Fe}}{\Delta t} =& \quad E_{zoo}^{\leftarrow Fe} 
                                             + \Gamma_{phy}^{\rightarrow C} Q_{phy}^{Fe:C} 
                                             + \Gamma_{zoo}^{\rightarrow C} Q_{zoo}^{Fe:C} \\
                                       & \quad - \left( \Gamma_{det}^{\rightarrow C} + g_{zoo}^{\leftarrow B_{det}^{C}} \right) Q_{det}^{Fe:C} \\
                                       & \quad + Sc_{dFe}^{\rightarrow B_{det}^{Fe}} + Co_{dFe}^{\rightarrow B_{det}^{Fe}}
\end{align}
$$

**Calcium Carbonate** (`p_caco3(i,j,k,tau)`, $B_{CaCO_3}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{CaCO_3}^{C}}{\Delta t} =& \quad P_{CaCO_3}^{\Gamma_{phy}^{\rightarrow C}} 
                                               + P_{CaCO_3}^{\Gamma_{zoo}^{\rightarrow C}} 
                                               + P_{CaCO_3}^{g_{zoo}^{\leftarrow B_{phy}^{C}}} \\
                                       & \quad - D_{CaCO_3}^{\Omega_{cal}}
                                               - D_{CaCO_3}^{\Omega_{ara}}
                                               - D_{CaCO_3}^{\Gamma_{det}^{\rightarrow C}}
                                               - D_{CaCO_3}^{g_{zoo}^{\leftarrow B_{det}^{C}}}

\end{align}
$$

**Dissolved Inorganic Carbon** (`p_dic(i,j,k,tau)`, $DIC$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta DIC}{\Delta t} =& \quad \Gamma_{det}^{\rightarrow C} + \gamma_{zoo}^{\rightarrow C} + \gamma_{phy}^{\rightarrow C} + X_{zoo}^{\leftarrow C} \\
                              & \quad - \mu_{phy}^{\leftarrow C} - \dfrac{\Delta CaCO_3}{\Delta t}
\end{align}
$$

**Alkalinity** (`p_alk(i,j,k,tau)`, $Alk$, [mol Eq kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta Alk}{\Delta t} =& \quad - \dfrac{\Delta NO_3}{\Delta t}
                                      - 2 \cdot \dfrac{\Delta CaCO_3}{\Delta t}
\end{align}
$$


---


### 14. Check for conservation of mass.

When checks for the conservation of mass is enabled (`do_check_n_conserve = .true.` or `do_check_c_conserve = .true.`), the model will calculate the budget of nitrogen or carbon both before and after the ecosystem equations have completed. This checks that the ecosystem equations detailed above have indeed conserved the mass of both nitrogen and carbon within the ocean. In WOMBATlite, both $NO_3$ and $DIC$ should be perfectly conserved during ecosystem cycling. However, the exception to this is for nitrogen, where if either of `do_nitrogen_fixation = .true.` or `do_benthic_denitrification = .true.` then the model does not and should not be expected to conserve nitrogen.

---

### 15. Additional operations on tracers.

**First**, dissolved iron concentrations are set to equal 1 nM everywhere where the depth of the water column is less than 200 metres deep. WOMBAT-lite is not considered to be a model of the coastal ocean, but rather a model of the global pelagic ocean. Given that coastal waters are not limited in dissolved iron due to substantial interactions with sediments and exchange with the land, we universally set the dissolved iron concentration in these waters to 1 nM.

**Second**, if dissolved iron concentrations dip below that measureable by operational detection limits in waters deeper than 200 m, we reset these concentrations to this minimum  (`dfefloor`, $[dFe]^{min}$, [nmol Fe kg<sup>-1</sup>]). $[dFe]^{min}$ is set in the parameter list and is configurable at run time.

---


### 16. Sinking rate of particulates.

WOMBAT-lite functions with a spatially variable sinking rate of both organic and inorganic (i.e., $CaCO_3$) particulate matter. The variable sinking rates of detritus and $CaCO_3$ are computed within the `generic_WOMBATlite_update_from_source` subroutine. Values are positive downward. Once computed, these 3D arrays are transfered to pointer arrays (`p_wdet(i,j,k)` and `p_wcaco3(i,j,k)`) that interface with a tridiagonal solver within the `g_tracer_vertdiff_G` subroutine.

**Organic detritus**

We consider that sinking rates of detritus are varied as a function of phytoplankton concentration (in a similar fashion to half-saturation coefficients described earlier), as well as the fraction of particulate matter that is $CaCO_3$. This approach is taken to emulate observations of varying sinking speeds [(Riley et al., 2012)](https://doi.org/10.1029/2011GB004085) and because such variations may be strongly dependent on phytoplankton community composition [(Bach et al., 2016)](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016GB005372).

In accordance with a more general Navier–Stokes drag equation and using a compilation of particle sinking speeds, [Cael et al. (2021)](https://doi.org/10.1029/2020GL091771) identified that the sinking velocity of detrital particles ($\omega_{det}$) in [m s<sup>-1</sup>]) is proportional to their diameter raised to the power of roughly 0.63, such that

$$
\begin{align}
\omega_{det} \propto& \quad d^{0.63}
\end{align}
$$

Knowing that 

$$
\begin{align}
d =& \quad 6 V \pi \dfrac{1}{3}
\end{align}
$$

and given that the average volume of phytoplankton cells can be approximated by $V = (B_{phy}^{C})^{0.65}$ [(Wickman et al., 2024)](https://doi.org/10.1126/science.adk6901), we can relate $\omega_{det}$ to the biomass concentration of phytoplankton multiplied by the scaler $\omega_{det}^{0}$

$$
\begin{align}
\omega_{det} =& \quad \omega_{det}^{0} \cdot \left( B_{phy}^{C,k=1} \right)^{0.21}
\end{align}
$$

_where_ <br>
- $\omega_{det}^{0}$ is the sinking speed of organic detritus produced by a surface phytoplankton concentration of 1 mmol C m<sup>-3</sup> (`wdetbio`, [m s<sup>-1</sup>]) <br>
- $B_{phy}^{C}$ is the concentration of phytoplankton biomass (`phy_mmolm3`, [mmol C m<sup>-3</sup>]) <br>

This formula is identical to that presented by [Cael et al. (2021)](https://doi.org/10.1029/2020GL091771) in their Eq. (3), with the exception that we have related sinking rates to the biomass concentration of phytoplankton ($B_{phy}^{C}$) by assuming that $V = (B_{phy}^{C})^{0.65}$ based on marine phytoplankton data [(Wickman et al., 2024)](https://doi.org/10.1126/science.adk6901).

As phytoplankton concentrations are negligible beneath the euphotic zone we use $B_{phy}^{C}$ only in the uppermost grid cell (k = 1). This assumes that the sinking velocities of marine aggregates can be related to phytoplankton community composition at the surface ([Bach et al., 2016](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016GB005372); [Iversen and Lampitt, 2020](https://doi.org/10.1016/j.pocean.2020.102445)), which varies more horizontally across the ocean than vertically. Moreover, because we do not include dissolved/suspended organic matter as a tracer in WOMBAT-lite, we must also account for the large fraction of organics that are suspended and thus neutrally buoyant in the gyres. As such, we include a phytoplankton biomass threshold (`phybiot`, $B_{phy}^{thresh}$, [mmol C m<sup>-3</sup>]) above which sinking accelerates and beneath which any produced detritus emulates dissolved (neutrally buoyant) organic matter:

$$
\begin{align}
\omega_{det} =& \quad \omega_{det}^{0} \cdot \max \left( 0.0, B_{phy}^{C,k=1} - B_{phy}^{thresh} \right)^{0.21}
\end{align}
$$

Sinking speeds are then accelerated depending on the fraction of $CaCO_3$ within the particulate mass by

$$
\begin{align}
\omega_{det} =& \quad \omega_{det} + \dfrac{10}{86400} \cdot \min \left( 1, \dfrac{B_{CaCO_3}^{C}}{B_{CaCO_3}^{C} + B_{det}^{C}} \right)
\end{align}
$$

_where_ <br>
- $B_{phy}^{C,k=1}$ is the concentration of phytoplankton biomass at the surface (`phy_mmolm3`, [mmol C m<sup>-3</sup>]) <br>
- $B_{phy}^{thresh}$ is the phytoplankton biomass threshold above which the community cell size begins to increase (`phybiot`, [mmol C m<sup>-3</sup>]) <br>
- $B_{CaCO_3}^{C}$ is the in situ concentration of $CaCO_3$ biomass (`p_caco3(i,j,k,tau)`, [mol C kg<sup>-1</sup>]) <br>
- $B_{det}^{C}$ is the in situ concentration of organic detrital biomass (`p_det(i,j,k,tau)`, [mol C kg<sup>-1</sup>]) <br>


and where we chose a maximum increase of 10 metres per day based on the more modest effect observed in mesocosm experiments ([Bach et al., 2016](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016GB005372).

Finally, we apply a linear increase to sinking speeds with depth to ensure that the trend in the concentration of detritus with depth exhibits a power law behavior, which is widely observed ([Berelson, 2001](https://doi.org/10.1016/S0967-0645(01)00102-3); [Martin et al., 1987](https://doi.org/10.1016/0198-0149(87)90086-0)), thought to be associated with a greater attenuation of more slowly sinking particles, and shows better performance than a constant sinking rate in models ([Tjiputra et al., 2020](https://gmd.copernicus.org/articles/13/2393/2020/)). This is applied after the previous equation as:

$$
\begin{align}
\omega_{det} =& \quad \omega_{det} + \max\left(0.0, \dfrac{z}{5000} \cdot \left(\omega_{det}^{max} - \omega_{det} \right)\right)
\end{align}
$$

_where_ <br>
- $z$ is the in situ depth (`zw(i,j,k)`, [m]) <br>
- $\omega_{det}^{max}$ is the terminal sinking rate of organic detritus (`wdetmax`, [m s<sup>-1</sup>]) <br>


**Calcium carbonate**

The sinking rate of $CaCO_3$ is considered to always be a fraction of the sinking rate of organic detritus, such that

$$
\begin{align}
\omega_{CaCO_3} =& \quad \omega_{det} \cdot \dfrac{\omega_{CaCO_3}^{0}}{\omega_{det}^{0}}
\end{align}
$$

_where_ <br>
- $\omega_{det}^{0}$ is the sinking speed of organic detritus produced by a surface phytoplankton concentration of 1 mmol C m<sup>-3</sup> (`wdetbio`, [m s<sup>-1</sup>]) <br>
- $\omega_{CaCO_3}^{0}$ is the sinking speed of $CaCO_3$ produced by a surface phytoplankton concentration of 1 mmol C m<sup>-3</sup> (`wcaco3`, [m s<sup>-1</sup>]) <br>

Although more dense than organic matter, $CaCO_3$ particles tend to be smaller than organic aggregates and sink at a slower rate ([De La Rocha & Passow, 2007](https://doi.org/10.1016/j.dsr2.2007.01.004); [Zhang et al., 2018](https://doi.org/10.5194/bg-15-4759-2018)). Furthermore, the shedding of coccoliths by coccolithophores, which are near-neutrally bouyant, also contributes to a slower mean sinking speed of $CaCO_3$ ([Balch et al., 2009](https://doi.org/10.1029/2008JC004902)). Hence, $\omega_{CaCO_3}^{0} should be < \omega_{det}^{0}$.

The degree to which $CaCO_3$ particles sink more slowly than organic detritus is controlled by the user through altering the ratio of $\omega_{CaCO_3}^{0}$ to $\omega_{det}^{0}$.


**Detrital iron**

WOMBAT-lite sinks organic detrital iron at the same rate as organic detrital carbon.

---


### 17. Sedimentary processes.

WOMBAT-lite tracks the accumulation of organic detrital carbon (`p_det_sediment(i,j)`, $B_{det,sed}^{C}$, [mol m<sup>-2</sup>]), organic detrital iron (`p_detfe_sediment(i,j)`, $B_{det,sed}^{Fe}$, [mol m<sup>-2</sup>]) and $CaCO_3$ (`p_caco3_sediment(i,j)`, $B_{CaCO_3,sed}^{C}$, [mol m<sup>-2</sup>]) within sedimentary pools. The organic pools contribute to bottom fluxes of dissolved inorganic carbon (DIC), nitrate ($NO_3$), dissolved iron (dFe), oxygen ($O_2$) and alkalinity (Alk). Remineralisation of organic carbon ($\gamma_{det,sed}^{C}$) produces DIC and $NO_3$, but removes $O_2$ and Alk. Ratios of nitrogen to carbon and oxygen to carbon are static at 16:122 and -172:122. Remineralisation of organic iron produces dFe.

**Organics**

$$
\begin{align}
\gamma_{det,sed}^{C} =& \quad \gamma_{det,sed}^{0^{\circ}C} (β_{hete})^{T} B_{det,sed}^{C} \\
\gamma_{det,sed}^{N} =& \quad \gamma_{det,sed}^{C} R^{N:C} - \gamma_{det,sed}^{denit} \\
\gamma_{det,sed}^{O_2} =& \quad \gamma_{det,sed}^{C} R^{O_2:C} \cdot \left(1 - f_{denit} \right) \\
\gamma_{det,sed}^{Alk} =& \quad -\gamma_{det,sed}^{N} \\
\gamma_{det,sed}^{Fe} =& \quad \gamma_{det,sed}^{0^{\circ}C} (β_{hete})^{T} B_{det,sed}^{Fe}
\end{align}
$$

_where_ <br>
- $\gamma_{det,sed}^{0^{\circ}C}$ is the base linear rate of remineralisation of sedimentary organic detritus (`detlrem_sed`, [s<sup>-1</sup>]) <br>
- $\gamma_{det,sed}^{C}$ is the remineralisation rate of organic detrital carbon in the sediment pool (`det_sed_remin(i,j)`, [mol C m<sup>-2</sup> s<sup>-1</sup>]) <br>
- $\gamma_{det,sed}^{N}$ is the production of nitrate due to organic detrital remineralisation in the sediment pool ([mol N m<sup>-2</sup> s<sup>-1</sup>]) <br>
- $\gamma_{det,sed}^{denit}$ is the consumption rate of NO<sub>3</sub> due to benthic denitrification (`det_sed_denit(i,j)`, [mol N m<sup>-2</sup> s<sup>-1</sup>]) <br>
- $\gamma_{det,sed}^{O_2}$ is the consumption rate of oxygen during remineralisation of organics in the sediment pool ([mol $O_2$ m<sup>-2</sup> s<sup>-1</sup>]) <br>
- $\gamma_{det,sed}^{Alk}$ is the production of alkalinity due to organic detrital remineralisation in the sediment pool ([mol Alk m<sup>-2</sup> s<sup>-1</sup>]) <br>
- $\gamma_{det,sed}^{Fe}$ is the production of dissolved iron due to organic detrital remineralisation in the sediment pool ([mol Fe m<sup>-2</sup> s<sup>-1</sup>]) <br>
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless]) <br>
- $T$ is the in situ temperature at the sediment-water interface (`sedtemp(i,j)`, [ºC]) <br>
- $f_{denit}$ is the fraction of organic matter being remineralised by anaerobic metabolism, i.e., denitrification (`fdenit(i,j)`, [dimensionless]) <br>
- $B_{det,sed}^{C}$ is the total concentration of organic detrital carbon biomass (`p_det_sediment(i,j,1)`, [mol C m<sup>-2</sup>]) <br>
- $B_{det,sed}^{Fe}$ is the total concentration of organic detrital iron biomass (`p_detfe_sediment(i,j,1)`, [mol Fe m<sup>-2</sup>]) <br>
- $R^{N:C}$ is the static ratio of nitrogen to carbon in organic matter and is equal to 16:122 <br>
- $R^{O_2:C}$ is the static ratio of dissolved oxygen utilisation to carbon content in organic matter and is equal to -172:122 <br>

**Dissolution of $CaCO_3$**

Alk and DIC are in turn affected by dissolution of the sedimenary $CaCO_3$ pool, which is considered as entirely calcite. Sedimentary dissolution is controlled by bottom water temperature and an estimate of the pore-water calcite saturation state ($\Omega_{cal}^{sed}$):

$$
\begin{align}
D_{cal,sed} =& \quad d_{cal,sed} (β_{hete})^{T} \max\left(1 - \Omega_{cal,+sed},  1 - \Omega_{cal,sed}\right)^{4.5}
\end{align}
$$

_where_ <br>
- $d_{cal,sed}$ is a base rate of dissolution (`caco3lrem_sed`, [s<sup>-1</sup>]) <br>
- $\Omega_{cal,+sed}$ is the maximum possible saturation state within sediment pore water (`omegamax_sed`, [dimensionless]) <br>

The $\Omega_{cal,sed}$ is calculated using the [`mocsy`](https://github.com/ACCESS-NRI/mocsy) package for solving carbonate chemistry of seawater ([Orr & Epitalon, 2015](https://doi.org/10.5194/gmd-8-485-2015)). These routines require Alk and DIC as inputs, along with nutrient concentrations and temperature and salinity of bottom waters. For DIC, we chose to sum the water column concentration of DIC and the organic carbon content of the sediment to approximate the interstitial (i.e., porewater) DIC concentration. We assume that the organic carbon content of the sediment (`p_det_sediment` [mol m<sup>-2</sup>]) is relevant over one meter, and therefore can be automatically converted to [mol m<sup>-3</sup>]. With this assumption these arrays can be added together and approximates the reducing conditions of organic-rich sediments, which have lower $\Omega_{cal,sed}$. This ensures a greater rate of $CaCO_3$ dissolution within the sediment as organic matter accumulates.

**Benthic denitrification**

We also consider the consumption of NO<sub>3</sub> via benthic denitrification. When `do_benthic_denitrification = .true.`, a portion of the particulate organic matter within the sediments that remineralised is performed anaerobically (i.e., using NO<sub>3</sub> as the electron acceptor). Unlike this process in the water column, which is performed by bacterial metabolism, we estimate this process using an empirical parameterization from [Bohlen et al. (2012)](https://doi.org/10.1029/2011GB004198):

$$
\begin{align}
\gamma_{det,sed}^{denit} =& \quad \gamma_{det,sed}^{C} \min\left(0.9 \dfrac{94}{122}, \left(0.083 + 0.21 \cdot 0.98^{O_2 - NO_3} \right) \right)
\end{align}
$$

_where_ <br>
- $0.9$ is a hard upper limit stating that 90% of organic matter hydrolysation can potentially be performed anaerobically via denitrification <br>
- $\dfrac{94}{122}$ is the stoichiometry of nitrate demand per mol of organic carbon hydrolysed ([Paulmier et al., 2009](https://doi.org/10.5194/bg-6-923-2009)) <br>
- O<sub>2</sub> is the bottom water concentration of dissolved oxygen (mmol m<sup>-3</sup>) <br>
- NO<sub>3</sub> is the bottom water concentration of nitrate (mmol m<sup>-3</sup>) <br>

and where the fraction of organic matter that is remineralised via denitrification is equal to:

$$
\begin{align}
f_{denit} =& \quad \gamma_{det,sed}^{denit} \dfrac{122/94}{\gamma_{det,sed}^{C}}
\end{align}
$$

Overall bottom fluxes of tracers are:

$$
\begin{align}
\dfrac{\Delta NO_3}{\Delta t} =& \quad \gamma_{det,sed}^{N} \\
\dfrac{\Delta O_2}{\Delta t} =& \quad \gamma_{det,sed}^{O_2} \\
\dfrac{\Delta dFe}{\Delta t} =& \quad \gamma_{det,sed}^{Fe} \\
\dfrac{\Delta DIC}{\Delta t} =& \quad \gamma_{det,sed}^{C} + D_{cal,sed} \\
\dfrac{\Delta Alk}{\Delta t} =& \quad \gamma_{det,sed}^{Alk} + 2 \cdot D_{cal,sed} 
\end{align}
$$

_where_ <br>
- $\dfrac{\Delta NO_3}{\Delta t}$ is the total flux of nitrate into the water column (`b_no3(i,j)`, [mol N m<sup>-2</sup> s<sup>-1</sup>]) <br>
- $\dfrac{\Delta O_2}{\Delta t}$ is the total flux of oxygen into the water column (`b_o2(i,j)`, [mol $O_2$ m<sup>-2</sup> s<sup>-1</sup>]) <br>
- $\dfrac{\Delta dFe}{\Delta t}$ is the total flux of dissolved iron into the water column (`b_fe(i,j)`, [mol Fe m<sup>-2</sup> s<sup>-1</sup>]) <br>
- $\dfrac{\Delta DIC}{\Delta t}$ is the total flux of dissolved inorganic carbon into the water column (`b_dic(i,j)`, [mol C m<sup>-2</sup>/s<sup>-1</sup>]) <br>
- $\dfrac{\Delta Alk}{\Delta t}$ is the total flux of alkalinity into the water column (`b_alk(i,j)`, [mol Alk m<sup>-2</sup> s<sup>-1</sup>]) <br>

---

## Subroutine - "update_from_bottom"

The subroutine `generic_WOMBATlite_update_from_bottom` moves sinking organic material from the water column into the sediment pools.
It is at this point that the model performs permanent burial of sinking organic matter if desired.

---

### Permanent burial of particulates.

If `do_burial = .true.`, we compute the fraction of incident sinking organic matter and $CaCO_3$ that is permanently buried in the sediments. This permanently buried fraction is effectively removed from the model and therefore is not accumulated within the sedimentary pools.

The fraction buried is calculated according to Equation 3 of [Dunne et al. (2007)](https://doi.org/10.1029/2006GB002907):

$$
\begin{align}
F_{bury} =& \quad 0.013 \cdot 0.53 \dfrac{\left(f_{org}\right)^{2}}{\left(7 + f_{org}\right)^{2}}
\end{align}
$$

where $f_{org}$ is the rain rate of organic carbon detritus on the seafloor in [mmol C m<sup>-2</sup> s<sup>-1</sup>].

As organic matter rains down at a more rapid rate, the fraction of incident organic carbon, organic iron and $CACO_3$ that is buried increases.

If `do_conserve_tracers = .true.`, then we capture the total loss of both Alk and $NO_3$ via burial or organic detritus and $CaCO_3$ and redistribute the Alk and $NO_3$ amount back at the ocean surface. This amount of each tracer is redistributed uniformly to avoid strong gradients.

---
