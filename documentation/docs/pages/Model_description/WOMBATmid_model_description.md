# Description of the WOMBATmid ocean biogeochemical model

```
        (\___/)  .-.   .-. .--. .-..-..---.  .--. .-----.
        / o o \  : :.-.: :: ,. :: '' :: .; :: .; :'-. .-'
       (   "   ) : :: :: :: :: :: .. ::   .':    :  : :
        \__ __/  : '' '' ;: :; :: :; :: .; :: :: :  : :
                  '.,'.,' '.__.':_;:_;:___.':_;:_;  :_;

 World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT)
```

_Contact Pearse J. Buchanan and/or Dougal Squire for any questions_

_Pearse.Buchanan@csiro.au_
_Dougie.Squire@anu.edu.au_

---

![schematic](/assets/schematic-mid.png)

## Tracers

The following are the active tracers in WOMBAT-mid

| #    | Tracer                        | Code name  | Description                                         | Units                     | Default on? |
|------|-------------------------------|------------|-----------------------------------------------------|---------------------------|-------------|
| 1    | O<sub>2</sub>                 | `p_o2`     | Dissolved oxygen                                    | mol O2 kg<sup>-1</sup>    | Yes         |
| 2    | NH<sub>4</sub>                | `p_nh4`    | Ammonium                                            | mol N kg<sup>-1</sup>     | Yes         |
| 3    | NO<sub>3</sub>                | `p_no3`    | Nitrate                                             | mol N kg<sup>-1</sup>     | Yes         |
| 4    | Si(OH)<sub>4</sub>            | `p_sil`    | Silicic acid                                        | mol Si kg<sup>-1</sup>    | Yes         |
| 5    | N<sub>2</sub>O                | `p_n2o`    | Nitrous oxide                                       | mol N kg<sup>-1</sup>     | Yes         |
| 6    | dFe                           | `p_fe`     | Dissolved iron                                      | mol Fe kg<sup>-1</sup>    | Yes         |
| 7    | Fe<sub>sA</sub>               | `p_afe`    | Small sinking authigenic iron                       | mol Fe kg<sup>-1</sup>    | Yes         |
| 8    | Fe<sub>lA</sub>               | `p_bafe`   | Large sinking authigenic iron                       | mol Fe kg<sup>-1</sup>    | Yes         |
| 9    | B<sub>np</sub><sup>C</sup>    | `p_phy`    | Nano-phytoplankton                                  | mol C kg<sup>-1</sup>     | Yes         |
| 10   | B<sub>mp</sub><sup>C</sup>    | `p_dia`    | Micro-phytoplankton                                 | mol C kg<sup>-1</sup>     | Yes         |
| 11   | B<sub>mz</sub><sup>C</sup>    | `p_zoo`    | Micro-zooplankton                                   | mol C kg<sup>-1</sup>     | Yes         |
| 12   | B<sub>Mz</sub><sup>C</sup>    | `p_mes`    | Meso-zooplankton                                    | mol C kg<sup>-1</sup>     | Yes         |
| 13   | B<sub>sd</sub><sup>C</sup>    | `p_det`    | Small sinking detritus                              | mol C kg<sup>-1</sup>     | Yes         |
| 14   | B<sub>ld</sub><sup>C</sup>    | `p_bdet`   | Large sinking detritus                              | mol C kg<sup>-1</sup>     | Yes         |
| 15   | B<sub>np</sub><sup>Chl</sup>  | `p_pchl`   | Nano-phytoplankton chlorophyll content              | mol C kg<sup>-1</sup>     | Yes         |
| 16   | B<sub>mp</sub><sup>Chl</sup>  | `p_dchl`   | Micro-phytoplankton chlorophyll content             | mol C kg<sup>-1</sup>     | Yes         |
| 17   | B<sub>np</sub><sup>Fe</sup>   | `p_phyfe`  | Nano-phytoplankton iron content                     | mol Fe kg<sup>-1</sup>    | Yes         |
| 18   | B<sub>mp</sub><sup>Fe</sup>   | `p_diafe`  | Micro-phytoplankton iron content                    | mol Fe kg<sup>-1</sup>    | Yes         |
| 19   | B<sub>mp</sub><sup>Si</sup>   | `p_diasi`  | Micro-phytoplankton silicon content                 | mol Si kg<sup>-1</sup>    | Yes         |
| 20   | B<sub>mz</sub><sup>Fe</sup>   | `p_zoofe`  | Micro-zooplankton iron content                      | mol Fe kg<sup>-1</sup>    | Yes         |
| 21   | B<sub>Mz</sub><sup>Fe</sup>   | `p_mesfe`  | Meso-zooplankton iron content                       | mol Fe kg<sup>-1</sup>    | Yes         |
| 22   | B<sub>sd</sub><sup>Fe</sup>   | `p_detfe`  | Small sinking detritus iron content                 | mol Fe kg<sup>-1</sup>    | Yes         |
| 23   | B<sub>ld</sub><sup>Fe</sup>   | `p_bdetfe` | Large sinking detritus iron content                 | mol Fe kg<sup>-1</sup>    | Yes         |
| 24   | B<sub>ld</sub><sup>Si</sup>   | `p_bdetsi` | Large sinking detritus silicon content              | mol Si kg<sup>-1</sup>    | Yes         |
| 25   | B<sub>DOM</sub><sup>C</sup>   | `p_doc`    | Dissolved organic carbon                            | mol C kg<sup>-1</sup>     | Yes         |
| 26   | B<sub>DOM</sub><sup>H</sup>   | `p_doh`    | Dissolved organic hydrogen                          | mol N kg<sup>-1</sup>     | Yes         |
| 27   | B<sub>DOM</sub><sup>O</sup>   | `p_doo`    | Dissolved organic oxygen                            | mol N kg<sup>-1</sup>     | Yes         |
| 28   | B<sub>DOM</sub><sup>N</sup>   | `p_don`    | Dissolved organic nitrogen                          | mol N kg<sup>-1</sup>     | Yes         |
| 29   | B<sub>aoa</sub><sup>C</sup>   | `p_aoa`    | Ammonia oxidizing archaea                           | mol C kg<sup>-1</sup>     | Yes         |
| 30   | B<sub>b-p</sub><sup>C</sup>   | `p_bacp`   | Particle-associated heterotrophic bacteria          | mol C kg<sup>-1</sup>     | Yes         |
| 31   | B<sub>b-f1</sub><sup>C</sup>  | `p_bacf1`  | Free-living heterotrophic bacterial type #1         | mol C kg<sup>-1</sup>     | Yes         |
| 32   | B<sub>b-f2</sub><sup>C</sup>  | `p_bacf2`  | Free-living heterotrophic bacterial type #2         | mol C kg<sup>-1</sup>     | Yes         |
| 33   | DIC                           | `p_dic`    | Dissolved inorganic carbon                          | mol C kg<sup>-1</sup>     | Yes         |
| 34   | Alk                           | `p_alk`    | Dissolved alkalinity                                | mol Eq kg<sup>-1</sup>    | Yes         |
| 35   | CaCO<sub>3</sub>              | `p_caco3`  | Calcium carbonate                                   | mol C kg<sup>-1</sup>     | Yes         |
| 36   | DICp                          | -          | Preformed dissolved inorganic carbon                | mol C kg<sup>-1</sup>     | No          |
| 37   | DICr                          | `p_dicr`   | Remineralised dissolved inorganic carbon            | mol C kg<sup>-1</sup>     | No          |

---

## Logical controls

The following are logical statements within the `input.nml` namelist file that can be switched to TRUE or FALSE at runtime. 

| Logical                      | Description                                                                                | Default   |
|------------------------------|--------------------------------------------------------------------------------------------|-----------|
| `do_caco3_dynamics`          | Production and dissolution of CaCO3 depends on carbon system state                         | .true.    |
| `do_colloidal_shunt`         | Fraction of dissolved iron is colloids that coagulate onto sinking material                | .true.    |
| `do_two_ligands`             | Complex soluble iron using two ligands (weak + strong) rather than one                     | .false.   |
| `do_burial`                  | Permanently bury a fraction of sinking detrital material into the sediments                | .false.   |
| `do_nitrogen_fixation`       | Do implicit nitrogen fixation                                                              | .true.    |
| `do_anammox`                 | Do implicit anaerobic ammonium oxidation                                                   | .true.    |
| `do_wc_denitrification`      | Do anaerobic heterotrophic metabolism of bacteria using NO<sub>3</sub> and N<sub>2</sub>O  | .true.    |
| `do_benthic_denitrification` | Do implicit reduction of NO<sub>3</sub> in the sediment                                    | .true.    |
| `do_tracer_dicp`             | Carry preformed dissolved inorganic carbon (dicp) as a tracer                              | .false.   |
| `do_tracer_dicr`             | Carry remineralised dissolved inorganic carbon (dicr) as a tracer                          | .false.   |
| `do_viscous_sinking`         | Rubey's formula uses a non-constant dynamic viscosity of seawater                          | .true.    |
| `do_check_n_conserve`        | Checks that the ecosystem calculations are conserving the mass of nitrogen                 | .false.   |
| `do_check_c_conserve`        | Checks that the ecosystem calculations are conserving the mass of carbon                   | .false.   |
| `do_check_si_conserve`       | Checks that the ecosystem calculations are conserving the mass of silicon                  | .false.   |
| `do_check_fe_conserve`       | Checks that the ecosystem calculations are conserving the mass of iron                     | .false.   |

We note that when `do_two_ligands` is set to `.true.`, the `ligK` diagnostic variable reflects the binding strength of the strong ligand. However, when `do_two_ligands` is set to `.false.`, this diagnostic (`ligK`) reflects the binding strength of the bulk ligand pool.

---

## Diagnostic outputs

The following are all **2D** diagnostic output variables from WOMBAT-mid.

| Diagnostic        | Description                                                                                          | Units                                |
| ----------------- | ---------------------------------------------------------------------------------------------------- | ------------------------------------ |
| `pco2`            | Surface aqueous partial pressure of CO₂                                                              | µatm                                 |
| `npp2d`           | Vertically integrated net primary production                                                         | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `rpp2d`           | Vertically integrated regenerated primary production                                                 | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `zsp2d`           | Vertically integrated zooplankton secondary production                                               | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `det_radius`      | Mean radius of small detrital particles                                                              | m                                    |
| `bdet_radius`     | Mean radius of large detrital particles                                                              | m                                    |
| `det_sed_remin`   | Rate of remineralisation of detritus in accumulated sediment                                         | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `det_sed_depst`   | Rate of deposition of detritus to sediment at base of water column                                   | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `det_sed_denit`   | Rate of benthic denitrification (removal of NO<sub>3</sub>) in accumulated sediment                  | mol N m<sup>-2</sup> s<sup>-1</sup>  |
| `fbury`           | Fraction of deposited detritus permanently buried beneath sediment                                   | dimensionless                        |
| `fdenit`          | Fraction of sedimentary detritus remineralised via denitrification                                   | dimensionless                        |
| `detfe_sed_remin` | Rate of remineralisation of detrital iron in accumulated sediment                                    | mol Fe m<sup>-2</sup> s<sup>-1</sup> |
| `detfe_sed_depst` | Rate of deposition of detrital iron to sediment at base of water column                              | mol Fe m<sup>-2</sup> s<sup>-1</sup> |
| `detsi_sed_remin` | Rate of remineralisation of detrital silicon in accumulated sediment                                 | mol Si m<sup>-2</sup> s<sup>-1</sup> |
| `detsi_sed_depst` | Rate of deposition of detrital silicon to sediment at base of water column                           | mol Si m<sup>-2</sup> s<sup>-1</sup> |
| `caco3_sed_remin` | Rate of remineralisation of CaCO₃ in accumulated sediment                                            | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `caco3_sed_depst` | Rate of deposition of CaCO₃ to sediment at base of water column                                      | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `zeuphot`         | Depth of the euphotic zone (1% incident light)                                                       | m                                    |
| `seddep`          | Depth of the bottom layer                                                                            | m                                    |
| `sedmask`         | Mask of active sediment points                                                                       | dimensionless                        |
| `sedtemp`         | Temperature in the bottom layer                                                                      | °C                                   |
| `sedsalt`         | Salinity in the bottom layer                                                                         | psu                                  |
| `sedno3`          | Nitrate concentration in the bottom layer                                                            | mol N kg<sup>-1</sup>                |
| `sednh4`          | Ammonium concentration in the bottom layer                                                           | mol N kg<sup>-1</sup>                |
| `sedsil`          | Silicic acid concentration in the bottom layer                                                       | mol Si kg<sup>-1</sup>               |
| `seddic`          | Dissolved inorganic carbon concentration in the bottom layer                                         | mol C kg<sup>-1</sup>                |
| `sedalk`          | Alkalinity concentration in the bottom layer                                                         | mol Eq kg<sup>-1</sup>               |
| `sedhtotal`       | H<sup>+</sup> ion concentration in the bottom layer                                                  | mol H<sup>+</sup> kg<sup>-1</sup>    |
| `sedco3`          | CO₃<sup>2−</sup> ion concentration in the bottom layer                                               | mol C kg<sup>-1</sup>                |
| `sedomega_cal`    | Calcite saturation state in the bottom layer                                                         | dimensionless                        |
| `o2_stf`          | Surface flux of dissolved oxygen into ocean                                                          | mol O2 m<sup>-2</sup> s<sup>-1</sup> |
| `n2o_stf`         | Surface flux of nitrous oxide into ocean                                                             | mol N m<sup>-2</sup> s<sup>-1</sup>  |
| `nh4_stf`         | Surface flux of ammonium into ocean                                                                  | mol N m<sup>-2</sup> s<sup>-1</sup>  |
| `no3_stf`         | Surface flux of nitrate into ocean                                                                   | mol N m<sup>-2</sup> s<sup>-1</sup>  |
| `sil_stf`         | Surface flux of silicic acid into ocean                                                              | mol Si m<sup>-2</sup> s<sup>-1</sup> |
| `fe_stf`          | Surface flux of dissolved iron into ocean                                                            | mol Fe m<sup>-2</sup> s<sup>-1</sup> |
| `det_stf`         | Surface flux of small sinking detritus into ocean                                                    | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `bdet_stf`        | Surface flux of large sinking detritus into ocean                                                    | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `doc_stf`         | Surface flux of dissolved organic carbon into ocean                                                  | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `doh_stf`         | Surface flux of dissolved organic hydrogen into ocean                                                | mol H m<sup>-2</sup> s<sup>-1</sup>  |
| `doo_stf`         | Surface flux of dissolved organic oxygen into ocean                                                  | mol O m<sup>-2</sup> s<sup>-1</sup>  |
| `don_stf`         | Surface flux of dissolved organic nitrogen into ocean                                                | mol N m<sup>-2</sup> s<sup>-1</sup>  |
| `dic_stf`         | Surface flux of dissolved inorganic carbon into ocean                                                | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `dicp_stf`        | Surface flux of preformed dissolved inorganic carbon into ocean                                      | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `alk_stf`         | Surface flux of alkalinity into ocean                                                                | mol Eq m<sup>-2</sup> s<sup>-1</sup> |
| `no3_vstf`        | Virtual flux of nitrate into ocean due to salinity restoring/correction                              | mol N m<sup>-2</sup> s<sup>-1</sup>  |
| `nh4_vstf`        | Virtual flux of ammonium into ocean due to salinity restoring/correction                             | mol N m<sup>-2</sup> s<sup>-1</sup>  |
| `dic_vstf`        | Virtual flux of dissolved inorganic carbon into ocean due to salinity restoring/correction           | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `dicp_vstf`       | Virtual flux of preformed dissolved inorganic carbon into ocean due to salinity restoring/correction | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `alk_vstf`        | Virtual flux of alkalinity into ocean due to salinity restoring/correction                           | mol Eq m<sup>-2</sup> s<sup>-1</sup> |
| `o2_btf`          | Bottom flux of dissolved oxygen into ocean                                                           | mol O2 m<sup>-2</sup> s<sup>-1</sup> |
| `no3_btf`         | Bottom flux of nitrate into ocean                                                                    | mol N m<sup>-2</sup> s<sup>-1</sup>  |
| `sil_btf`         | Bottom flux of silicic acid into ocean                                                               | mol Si m<sup>-2</sup> s<sup>-1</sup> |
| `doc_btf`         | Bottom flux of dissolved organic carbon into ocean                                                   | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `doh_btf`         | Bottom flux of dissolved organic hydrogen into ocean                                                 | mol H m<sup>-2</sup> s<sup>-1</sup>  |
| `doo_btf`         | Bottom flux of dissolved organic oxygen into ocean                                                   | mol O m<sup>-2</sup> s<sup>-1</sup>  |
| `don_btf`         | Bottom flux of dissolved organic nitrogen into ocean                                                 | mol N m<sup>-2</sup> s<sup>-1</sup>  |
| `fe_btf`          | Bottom flux of dissolved iron into ocean                                                             | mol Fe m<sup>-2</sup> s<sup>-1</sup> |
| `dic_btf`         | Bottom flux of dissolved inorganic carbon into ocean                                                 | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `dicr_btf`        | Bottom flux of preformed dissolved inorganic carbon into ocean                                       | mol C m<sup>-2</sup> s<sup>-1</sup>  |
| `alk_btf`         | Bottom flux of alkalinity into ocean                                                                 | mol Eq m<sup>-2</sup> s<sup>-1</sup> |


The following are all **3D** diagnostic output variables from WOMBAT-mid.

| Diagnostic         | Description                                                                        | Units                                           |
| ------------------ | ---------------------------------------------------------------------------------- | ----------------------------------------------- |
| `htotal`           | Concentration of H<sup>+</sup> ion                                                 | mol H<sup>+</sup> kg<sup>-1</sup>               |
| `omega_ara`        | Saturation state of aragonite                                                      | dimensionless                                   |
| `omega_cal`        | Saturation state of calcite                                                        | dimensionless                                   |
| `co3`              | Carbonate ion concentration                                                        | mol C kg<sup>-1</sup>                           |
| `co2_star`         | CO2* (CO2(g) + H2CO3) concentration                                                | mol C kg<sup>-1</sup>                           |
| `dynvis_sw`        | Seawater dynamic viscosity                                                         | kg m<sup>-1</sup> s<sup>-1</sup>                |
| `radbio`           | Photosynthetically active radiation available for phytoplankton growth             | W m<sup>-2</sup>                                |
| `radmid`           | Photosynthetically active radiation at centre point of grid cell                   | W m<sup>-2</sup>                                |
| `radmld`           | Photosynthetically active radiation averaged in mixed layer                        | W m<sup>-2</sup>                                |
| `npp3d`            | Net primary productivity                                                           | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `rpp3d`            | Regenerated primary productivity                                                   | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zsp3d`            | Zooplankton secondary productivity                                                 | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `phy_mumax`        | Maximum growth rate of nano-phytoplankton                                          | s<sup>-1</sup>                                  |
| `phy_mu`           | Realised growth rate of nano-phytoplankton                                         | s<sup>-1</sup>                                  |
| `pchl_mu`          | Realised growth rate of nano-phytoplankton chlorophyll                             | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `phy_lpar`         | Limitation of nano-phytoplankton by light                                          | dimensionless                                   |
| `phy_kni`          | Half-saturation coefficient of nitrogen uptake by nano-phytoplankton               | mmol N m<sup>-3</sup>                           |
| `phy_kfe`          | Half-saturation coefficient of iron uptake by nano-phytoplankton                   | µmol Fe m<sup>-3</sup>                          |
| `phy_lnit`         | Limitation of nano-phytoplankton by nitrogen                                       | dimensionless                                   |
| `phy_lnh4`         | Limitation of nano-phytoplankton by ammonium                                       | dimensionless                                   |
| `phy_lno3`         | Limitation of nano-phytoplankton by nitrate                                        | dimensionless                                   |
| `phy_lfer`         | Limitation of nano-phytoplankton by iron                                           | dimensionless                                   |
| `phy_dfeupt`       | Uptake of dFe by nano-phytoplankton                                                | mol Fe kg<sup>-1</sup> s<sup>-1</sup>           |
| `phy_feupreg`      | Factor up regulation of dFe uptake by nano-phytoplankton                           | dimensionless                                   |
| `phy_fedoreg`      | Factor down regulation of dFe uptake by nano-phytoplankton                         | dimensionless                                   |
| `phygrow`          | Growth of nano-phytoplankton                                                       | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `phydoc`           | Overflow exudation of DOC by nano-phytoplankton                                    | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `phymorl`          | Linear mortality of nano-phytoplankton                                             | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `phymorq`          | Quadratic mortality of nano-phytoplankton                                          | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `dia_mumax`        | Maximum growth rate of micro-phytoplankton                                         | s<sup>-1</sup>                                  |
| `dia_mu`           | Realised growth rate of micro-phytoplankton                                        | s<sup>-1</sup>                                  |
| `dchl_mu`          | Realised growth rate of micro-phytoplankton chlorophyll                            | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `dia_lpar`         | Limitation of micro-phytoplankton by light                                         | dimensionless                                   |
| `dia_kni`          | Half-saturation coefficient of nitrogen uptake by micro-phytoplankton              | mmol N m<sup>-3</sup>                           |
| `dia_kfe`          | Half-saturation coefficient of iron uptake by micro-phytoplankton                  | µmol Fe m<sup>-3</sup>                          |
| `dia_ksi`          | Half-saturation coefficient of silicic acid uptake by micro-phytoplankton          | µmol Si m<sup>-3</sup>                          |
| `dia_lnit`         | Limitation of micro-phytoplankton by nitrogen                                      | dimensionless                                   |
| `dia_lnh4`         | Limitation of micro-phytoplankton by ammonium                                      | dimensionless                                   |
| `dia_lno3`         | Limitation of micro-phytoplankton by nitrate                                       | dimensionless                                   |
| `dia_lfer`         | Limitation of micro-phytoplankton by iron                                          | dimensionless                                   |
| `dia_lsil`         | Limitation of micro-phytoplankton by silicic acid                                  | dimensionless                                   |
| `dia_dfeupt`       | Uptake of dFe by micro-phytoplankton                                               | mol Fe kg<sup>-1</sup> s<sup>-1</sup>           |
| `dia_feupreg`      | Factor up regulation of dFe uptake by micro-phytoplankton                          | dimensionless                                   |
| `dia_fedoreg`      | Factor down regulation of dFe uptake by micro-phytoplankton                        | dimensionless                                   |
| `dia_silupt`       | Uptake of silicic acid by micro-phytoplankton                                      | mol Si kg<sup>-1</sup> s<sup>-1</sup>           |
| `dia_sidoreg`      | Factor down regulation of silicic acid uptake by micro-phytoplankton               | dimensionless                                   |
| `diagrow`          | Growth of micro-phytoplankton                                                      | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `diadoc`           | Overflow exudation of DOC by micro-phytoplankton                                   | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `diamorl`          | Linear mortality of micro-phytoplankton                                            | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `diamorq`          | Quadratic (density-dependent) mortality of micro-phytoplankton                     | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `nitrfix`          | Nitrogen fixation rate (NH<sub>4</sub> production)                                 | mol N kg<sup>-1</sup> s<sup>-1</sup>            |
| `tri_lpar`         | Limitation of implicit trichodesmium by light                                      | dimensionless                                   |
| `tri_lfer`         | Limitation of implicit trichodesmium by iron                                       | dimensionless                                   |
| `trimumax`         | Maximum growth rate of implicit trichodesmium                                      | s<sup>-1</sup>                                  |
| `sileqc`           | Equilibrium concentration of silicic acid                                          | mol Si kg<sup>-1</sup>                          |
| `disssi`           | Dissolution rate of biogenic silica                                                | s<sup>-1</sup>                                  |
| `bsidiss`          | Dissolution of biogenic silica                                                     | mol Si kg<sup>-1</sup> s<sup>-1</sup>           |
| `feIII`            | Free iron (Fe<sup>3+</sup>)                                                        | mol Fe kg<sup>-1</sup>                          |
| `ligK`             | Ligand stability constant                                                          | L mol<sup>-1</sup>                              |
| `felig`            | Ligand-bound dissolved iron                                                        | mol Fe kg<sup>-1</sup>                          |
| `fecol`            | Colloidal dissolved iron                                                           | mol Fe kg<sup>-1</sup>                          |
| `fescaafe`         | Scavenging of free Fe onto authigenic particles due to smaller organics            | mol Fe kg<sup>-1</sup> s<sup>-1</sup>           |
| `fescabafe`        | Scavenging of free Fe onto authigenic particles due to larger organics             | mol Fe kg<sup>-1</sup> s<sup>-1</sup>           |
| `fecaog2afe`       | Coagulation of colloidal dFe onto small authigenic particles                       | mol Fe kg<sup>-1</sup> s<sup>-1</sup>           |
| `fecoag2bafe`      | Coagulation of colloidal dFe onto large authigenic particles                       | mol Fe kg<sup>-1</sup> s<sup>-1</sup>           |
| `afediss`          | Dissolution of small colloidal authigenic Fe particles                             | mol Fe kg<sup>-1</sup> s<sup>-1</sup>           |
| `bafediss`         | Dissolution of large colloidal authigenic Fe particles                             | mol Fe kg<sup>-1</sup> s<sup>-1</sup>           |
| `fesources`        | Total source of dFe in water column                                                | mol Fe kg<sup>-1</sup> s<sup>-1</sup>           |
| `fesinks`          | Total sink of dFe in water column                                                  | mol Fe kg<sup>-1</sup> s<sup>-1</sup>           |
| `zooeps`           | Micro-zooplankton community-wide prey capture rate coefficient                     | m<sup>6</sup> mmolC<sup>-2</sup> s<sup>-1</sup> |
| `zooprefbacp`      | Grazing dietary fraction of micro-zooplankton on bacteria 1                        | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooprefbacf1`      | Grazing dietary fraction of micro-zooplankton on bacteria 2                        | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooprefaoa`       | Grazing dietary fraction of micro-zooplankton on ammonia oxidizing archaea         | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooprefphy`       | Grazing dietary fraction of micro-zooplankton on nano-phytoplankton                | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooprefdia`       | Grazing dietary fraction of micro-zooplankton on micro-phytoplankton               | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooprefdet`       | Grazing dietary fraction of micro-zooplankton on small detritus                    | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zoograzbacp`      | Grazing rate of micro-zooplankton on bacteria 1                                    | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zoograzbacf1`      | Grazing rate of micro-zooplankton on bacteria 2                                    | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zoograzaoa`       | Grazing rate of micro-zooplankton on ammonia oxidizing archaea                     | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zoograzphy`       | Grazing rate of micro-zooplankton on nano-phytoplankton                            | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zoograzdia`       | Grazing rate of micro-zooplankton on micro-phytoplankton                           | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zoograzdet`       | Grazing rate of micro-zooplankton on small detritus                                | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zoomorl`          | Linear mortality of micro-zooplankton                                              | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zoomorq`          | Quadratic (density-dependent) mortality of micro-zooplankton                       | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooexcrbacp`      | Excretion rate of micro-zooplankton eating bacteria 1                              | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooexcrbacf1`      | Excretion rate of micro-zooplankton eating bacteria 2                              | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooexcraoa`       | Excretion rate of micro-zooplankton eating ammonia oxidizing archaea               | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooexcrphy`       | Excretion rate of micro-zooplankton eating nano-phytoplankton                      | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooexcrdia`       | Excretion rate of micro-zooplankton eating micro-phytoplankton                     | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooexcrdet`       | Excretion rate of micro-zooplankton eating small detritus                          | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooegesbacp`      | Egestion rate of micro-zooplankton on bacteria 1                                   | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooegesbacf1`      | Egestion rate of micro-zooplankton on bacteria 2                                   | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooegesaoa`       | Egestion rate of micro-zooplankton on ammonia oxidizing archaea                    | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooegesphy`       | Egestion rate of micro-zooplankton on nano-phytoplankton                           | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooegesdia`       | Egestion rate of micro-zooplankton on micro-phytoplankton                          | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `zooegesdet`       | Egestion rate of micro-zooplankton on small detritus                               | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `meseps`           | Meso-zooplankton community-wide prey capture rate coefficient                      | m<sup>6</sup> mmolC<sup>-2</sup> s<sup>-1</sup> |
| `mesprefbacp`      | Grazing dietary fraction of meso-zooplankton on bacteria 1                         | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesprefbacf1`      | Grazing dietary fraction of meso-zooplankton on bacteria 2                         | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesprefaoa`       | Grazing dietary fraction of meso-zooplankton on ammonia oxidizing archaea          | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesprefphy`       | Grazing dietary fraction of meso-zooplankton on nano-phytoplankton                 | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesprefdia`       | Grazing dietary fraction of meso-zooplankton on micro-phytoplankton                | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesprefdet`       | Grazing dietary fraction of meso-zooplankton on small detritus                     | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesprefbdet`      | Grazing dietary fraction of meso-zooplankton on large detritus                     | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesprefzoo`       | Grazing dietary fraction of meso-zooplankton on micro-zooplankton                  | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesgrazbacp`      | Grazing rate of meso-zooplankton on bacteria 1                                     | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesgrazbacf1`      | Grazing rate of meso-zooplankton on bacteria 2                                     | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesgrazaoa`       | Grazing rate of meso-zooplankton on ammonia oxidizing archaea                      | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesgrazphy`       | Grazing rate of meso-zooplankton on nano-phytoplankton                             | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesgrazdia`       | Grazing rate of meso-zooplankton on micro-phytoplankton                            | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesgrazdet`       | Grazing rate of meso-zooplankton on small detritus                                 | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesgrazbdet`      | Grazing rate of meso-zooplankton on large detritus                                 | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesgrazzoo`       | Grazing rate of meso-zooplankton on micro-zooplankton                              | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesmorl`          | Linear mortality of meso-zooplankton                                               | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesmorq`          | Quadratic (density-dependent) mortality of meso-zooplankton                        | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesexcrbacp`      | Excretion rate of meso-zooplankton eating bacteria 1                               | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesexcrbacf1`      | Excretion rate of meso-zooplankton eating bacteria 2                               | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesexcraoa`       | Excretion rate of meso-zooplankton eating ammonia oxidizing archaea                | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesexcrphy`       | Excretion rate of meso-zooplankton eating nano-phytoplankton                       | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesexcrdia`       | Excretion rate of meso-zooplankton eating micro-phytoplankton                      | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesexcrdet`       | Excretion rate of meso-zooplankton eating small detritus                           | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesexcrbdet`      | Excretion rate of meso-zooplankton eating large detritus                           | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesexcrzoo`       | Excretion rate of meso-zooplankton eating micro-zooplankton                        | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesegesbacp`      | Egestion rate of meso-zooplankton on bacteria 1                                    | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesegesbacf1`      | Egestion rate of meso-zooplankton on bacteria 2                                    | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesegesaoa`       | Egestion rate of meso-zooplankton on ammonia oxidizing archaea                     | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesegesphy`       | Egestion rate of meso-zooplankton on nano-phytoplankton                            | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesegesdia`       | Egestion rate of meso-zooplankton on micro-phytoplankton                           | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesegesdet`       | Egestion rate of meso-zooplankton on small detritus                                | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesegesbdet`      | Egestion rate of meso-zooplankton on large detritus                                | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `mesegeszoo`       | Egestion rate of meso-zooplankton on micro-zooplankton                             | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `reminr`           | Rate of remineralisation                                                           | s<sup>-1</sup>                                  |
| `detremi`          | Hydrolysation of small sinking detritus                                            | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `bdetremi`         | Hydrolysation of large sinking detritus                                            | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `ammox`            | Ammonia oxidation rate (NH<sub>4</sub> consumption)                                | mol N kg<sup>-1</sup> s<sup>-1</sup>            |
| `aoa_loxy`         | Limitation of ammonia oxidizing archaea by oxygen                                  | dimensionless                                   |
| `aoa_lnh4`         | Limitation of ammonia oxidizing archaea by ammonium                                | dimensionless                                   |
| `aoa_en2o`         | Excretion of N<sub>2</sub>O produced by ammonia oxidizing archaea during oxidation | mol N (mol C Biomass)<sup>-1</sup>              |
| `aoa_eno3`         | Excretion of NO<sub>3</sub> produced by ammonia oxidizing archaea during oxidation | mol N (mol C Biomass)<sup>-1</sup>              |
| `aoa_mumax`        | Maximum growth rate of ammonia oxidizing archaea                                   | s<sup>-1</sup>                                  |
| `aoa_mu`           | Realized growth rate of ammonia oxidizing archaea                                  | s<sup>-1</sup>                                  |
| `aoagrow`          | Growth of ammonia oxidizing archaea                                                | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `aoaresp`          | Oxygen consumption of ammonia oxidizing archaea                                    | mol O2 kg<sup>-1</sup> s<sup>-1</sup>           |
| `aoamorl`          | Linear mortality of ammonia oxidizing archaea                                      | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `aoamorq`          | Quadratic mortality of ammonia oxidizing archaea                                   | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `anammox`          | Anammox rate (NH<sub>4</sub> consumption)                                          | mol kg<sup>-1</sup> s<sup>-1</sup>              |
| `aox_lnh4`         | Limitation of anammox bacteria by ammonium                                         | dimensionless                                   |
| `aox_mu`           | Realized growth rate of anammox bacteria                                           | s<sup>-1</sup>                                  |
| `pic2poc`          | Inorganic (CaCO3) to organic carbon ratio                                          | dimensionless                                   |
| `dissratcal`       | Dissolution rate of calcite CaCO3                                                  | s<sup>-1</sup>                                  |
| `dissratara`       | Dissolution rate of aragonite CaCO3                                                | s<sup>-1</sup>                                  |
| `dissratpoc`       | Dissolution rate of CaCO3 due to POC (detritus) remineralization                   | s<sup>-1</sup>                                  |
| `zoodiss`          | Dissolution of CaCO3 due to micro-zooplankton grazing                              | mol CaCO3 kg<sup>-1</sup> s<sup>-1</sup>        |
| `mesdiss`          | Dissolution of CaCO3 due to meso-zooplankton grazing                               | mol CaCO3 kg<sup>-1</sup> s<sup>-1</sup>        |
| `caldiss`          | Dissolution of calcite CaCO3                                                       | mol CaCO3 kg<sup>-1</sup> s<sup>-1</sup>        |
| `aradiss`          | Dissolution of aragonite CaCO3                                                     | mol CaCO3 kg<sup>-1</sup> s<sup>-1</sup>        |
| `pocdiss`          | Dissolution of CaCO3 due to POC remin                                              | mol CaCO3 kg<sup>-1</sup> s<sup>-1</sup>        |
| `poc1remi`         | Remineralisation of particulate organic carbon by particle-associated bacteria     | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `doc2remi`         | Remineralisation of dissolved organic carbon by free-living bacteria #1            | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `doc3remi`         | Remineralisation of dissolved organic nitrogen by free-living bacteria #2          | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `doc1prod`         | Production of dissolved organic carbon by particle-associated bacteria             | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `doc2prod`         | Production of dissolved organic carbon by free-living bacteria #1                     | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `doc3prod`         | Production of dissolved organic nitrogen by free-living bacteria #2                   | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `bacp_ypoc`        | Biomass yield of particle-associated bacteria (mol DOC per mol C biomass)          | mol DOC (mol C biomass)<sup>-1</sup>            |
| `bacf1_ydoc`        | Biomass yield of free-living bacteria (mol DOC per mol C biomass) #1                  | mol DOC (mol C biomass)<sup>-1</sup>            |
| `bacf2_ydoc`        | Biomass yield of free-living bacteria (mol DOC per mol C biomass) #2                  | mol DOC (mol C biomass)<sup>-1</sup>            |
| `bacpgrow`         | Growth of particle-associated heterotrophic bacteria                               | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `bacf1grow`         | Growth of free-living heterotrophic bacteria #1                                       | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `bacf2grow`         | Growth of free-living heterotrophic bacteria #2                                       | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `bacpresp`         | Oxygen consumption of particle-associated heterotrophic bacteria                   | mol O2 kg<sup>-1</sup> s<sup>-1</sup>           |
| `bacf1resp`         | Oxygen consumption of free-living heterotrophic bacteria #1                           | mol O2 kg<sup>-1</sup> s<sup>-1</sup>           |
| `bacf2resp`         | Oxygen consumption of free-living heterotrophic bacteria #2                           | mol O2 kg<sup>-1</sup> s<sup>-1</sup>           |
| `bacppco2`         | CO2 production by particle-associated heterotrophic bacteria                       | mol CO2 kg<sup>-1</sup> s<sup>-1</sup>          |
| `bacf1pco2`         | CO2 production by free-living heterotrophic bacteria #1                               | mol CO2 kg<sup>-1</sup> s<sup>-1</sup>          |
| `bacf2pco2`         | CO2 production by free-living heterotrophic bacteria #2                               | mol CO2 kg<sup>-1</sup> s<sup>-1</sup>          |
| `bacppnh4`         | NH4 production by particle-associated heterotrophic bacteria                       | mol NH4 kg<sup>-1</sup> s<sup>-1</sup>          |
| `bacf1pnh4`         | NH4 production by free-living heterotrophic bacteria #1                               | mol NH4 kg<sup>-1</sup> s<sup>-1</sup>          |
| `bacf2pnh4`         | NH4 production by free-living heterotrophic bacteria #2                               | mol NH4 kg<sup>-1</sup> s<sup>-1</sup>          |
| `bacpufer`         | Uptake of dFe by particle-associated heterotrophic bacteria                        | mol dFe kg<sup>-1</sup> s<sup>-1</sup>          |
| `bacf1ufer`         | Uptake of dFe by free-living heterotrophic bacteria #1                                | mol dFe kg<sup>-1</sup> s<sup>-1</sup>          |
| `bacf2ufer`         | Uptake of dFe by free-living heterotrophic bacteria #2                               | mol dFe kg<sup>-1</sup> s<sup>-1</sup>          |
| `bacp_mu`          | Realized growth rate of particle-associated heterotrophic bacteria                | s<sup>-1</sup>                                  |
| `bacf1_mu`          | Realized growth rate of free-living heterotrophic bacteria #1                         | s<sup>-1</sup>                                  |
| `bacf2_mu`          | Realized growth rate of free-living heterotrophic bacteria #2                         | s<sup>-1</sup>                                  |
| `bacp_fanaer`      | Fraction of particle-associated bacteria growth supported by anaeroby              | dimensionless                                   |
| `bacf1_fanaer`      | Fraction of free-living bacteria #1 growth supported by anaeroby                      | dimensionless                                   |
| `bacf2_fanaer`      | Fraction of free-living bacteria #2 growth supported by anaeroby                      | dimensionless                                   |
| `bacp_ffelim`      | Particle-associated bacteria growth limited by iron?                               | dimensionless                                   |
| `bacf1_ffelim`      | Free-living bacteria #1 growth limited by iron?                                       | dimensionless                                   |
| `bacf2_ffelim`      | Free-living bacteria #2 growth limited by iron?                                       | dimensionless                                   |
| `bacpmorl`         | Linear mortality of particle-associated heterotrophic bacteria                     | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `bacf1morl`         | Linear mortality of free-living heterotrophic bacteria #1                             | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `bacf2morl`         | Linear mortality of free-living heterotrophic bacteria #2                             | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `bacpmorq`         | Quadratic mortality of particle-associated heterotrophic bacteria                  | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `bacf1morq`         | Quadratic mortality of free-living heterotrophic bacteria #1                          | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `bacf2morq`         | Quadratic mortality of free-living heterotrophic bacteria #2                          | mol C kg<sup>-1</sup> s<sup>-1</sup>            |
| `bacpdeni`         | Particle-associated bacterial denitrification rate (NO<sub>3</sub> -> N<sub>2</sub>)                   | mol N kg<sup>-1</sup> s<sup>-1</sup>            |
| `bacf1deni`         | Free-living bacteria #1 denitrification rate (NO<sub>3</sub> -> N<sub>2</sub>O)                  | mol N kg<sup>-1</sup> s<sup>-1</sup>            |
| `bacf2deni`         | Free-living bacteria #2 denitrification rate (N<sub>2</sub>O consumption)                        | mol N2O kg<sup>-1</sup> s<sup>-1</sup>          |
| `bacp_rq`          | Respiratory quotient of particle-associated bacteria                               | mol CO<sub>2</sub> (mol O<sub>2</sub>)<sup>-1</sup> |
| `bacf1_rq`          | Respiratory quotient of free-living bacteria #1                                       | mol CO<sub>2</sub> (mol O<sub>2</sub>)<sup>-1</sup> |
| `bacf2_rq`          | Respiratory quotient of free-living bacteria #2                                       | mol CO<sub>2</sub> (mol O<sub>2</sub>)<sup>-1</sup> |
| `det_density`      | Mean density of small detrital particles                                           | kg m<sup>-3</sup>                               |
| `bdet_density`     | Mean density of large detrital particles                                           | kg m<sup>-3</sup>                               |
| `det_vmove`        | Sinking rate of small detritus                                                     | m s<sup>-1</sup>                                |
| `detfe_vmove`      | Sinking rate of small detrital iron                                                | m s<sup>-1</sup>                                |
| `bdet_vmove`       | Sinking rate of large detritus                                                     | m s<sup>-1</sup>                                |
| `bdetfe_vmove`     | Sinking rate of large detrital iron                                                | m s<sup>-1</sup>                                |
| `bdetsi_vmove`     | Sinking rate of large detrital silicon                                             | m s<sup>-1</sup>                                |
| `caco3_vmove`      | Sinking rate of CaCO3                                                              | m s<sup>-1</sup>                                |

---


## Subroutine - "update_from_source"

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
7. Synthesis of chlorophyll.
8. Phytoplankton uptake of iron.
9. Phytoplankton uptake of silicic acid.
10. Iron chemistry (scavenging, coagulation, dissolution).
11. Biogenic silica dissolution. 
12. Mortality terms.
13. Zooplankton grazing, egestion, excretion and assimilation.
14. Implicit nitrogen fixation.
15. Facultative bacterial heterotrophy.
16. Calcium carbonate production and dissolution.
17. Chemoautotrophy.
18. Tracer tendencies.
19. Check for conservation of mass.
20. Additional operations on tracers.
21. Sinking rate of particulates.
22. Sedimentary processes.

Below is a step‑by‑step explanation of each section together with the key equations. Variable names in `grey` follow the Fortran code, while 
variable names in $math font$ are pointers to the equations; `i,j,k` refer to horizontal and vertical indices; [square brackets] denote units. 
If a variable is without i,j,k dimensions, this variable is held as a scalar and not an array.

The model carries tracers in [mol kg-1]. That is, moles of solute/tracer per kilogram of seawater (i.e., molality). Some calculations herein are performed by converting tracers to units of [mmol m-3] or in the case of dissolved iron [µmol m-3]. However, we stress that all tracer tendency terms are converted back to [mol kg-1 s-1] when sources and sinks are applied.

---

### Parameter set and default values

| Parameter          | Description                                                                 | Value         | Units                                                                |
| ------------------ | --------------------------------------------------------------------------- | ------------- | -------------------------------------------------------------------- |
| `alphabio_phy`     | Initial slope of P–I curve (nano-phytoplankton)                             | 1.5           | mol C (mol Chl)<sup>-1</sup> (W m<sup>-2</sup>)<sup>-1</sup>         |
| `abioa_phy`        | Max growth rate parameter a (nano-phytoplankton)                            | 0.7/86400.0   | s<sup>-1</sup>                                                       |
| `bbioa_phy`        | Max growth rate parameter b (nano-phytoplankton) (Q10 = b^(10))             | 1.060         | dimensionless                                                        |
| `phyprefnh4`       | NH<sub>4</sub> preference over NO<sub>3</sub> (nano-phytoplankton)          | 5.0           | dimensionless                                                        |
| `phykn`            | Half-saturation coefficient N uptake (nano-phytoplankton)                   | 1.0           | mmol N m<sup>-3</sup>                                                |
| `phykf`            | Half-saturation coefficient Fe uptake (nano-phytoplankton)                  | 1.0           | µmol Fe m<sup>-3</sup>                                               |
| `phyminqc`         | Min Chl:C (nano-phytoplankton)                                              | 0.008         | mol Chl (mol C)<sup>-1</sup>                                         |
| `phymaxqc`         | Max Chl:C (nano-phytoplankton)                                              | 0.065         | mol Chl (mol C)<sup>-1</sup>                                         |
| `phyoptqf`         | Optimal Fe:C (nano-phytoplankton)                                           | 10e-6         | mol Fe (mol C)<sup>-1</sup>                                          |
| `phymaxqf`         | Max Fe:C (nano-phytoplankton)                                               | 50e-6         | mol Fe (mol C)<sup>-1</sup>                                          |
| `phylmor`          | Linear mortality rate (nano-phytoplankton)                                  | 0.001/86400.0 | s<sup>-1</sup>                                                       |
| `phyqmor`          | Quadratic mortality rate (nano-phytoplankton)                               | 0.05/86400.0  | (mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>                  |
| `phybiot`          | Biomass threshold (nano-phytoplankton)                                      | 1.0           | mmol C m<sup>-3</sup>                                                |
| `alphabio_dia`     | Initial slope of P–I curve (micro-phytoplankton)                            | 2.5           | mol C (mol Chl)<sup>-1</sup> (W m<sup>-2</sup>)<sup>-1</sup>         |
| `abioa_dia`        | Max growth rate parameter a (micro-phytoplankton)                           | 1.0/86400.0   | s<sup>-1</sup>                                                       |
| `bbioa_dia`        | Max growth rate parameter b (micro-phytoplankton) (Q10 = b^(10))            | 1.070         | dimensionless                                                        |
| `diaprefnh4`       | NH<sub>4</sub> preference over NO<sub>3</sub> (micro-phytoplankton)         | 5.0           | dimensionless                                                        |
| `diakn`            | Half-saturation coefficient N uptake (micro-phytoplankton)                  | 2.4           | mmol N m<sup>-3</sup>                                                |
| `diakf`            | Half-saturation coefficient Fe uptake (micro-phytoplankton)                 | 2.7           | µmol Fe m<sup>-3</sup>                                               |
| `diaks`            | Half-saturation coefficient Si uptake (micro-phytoplankton)                 | 6.7           | mmol Si m<sup>-3</sup>                                               |
| `diaminqc`         | Min Chl:C (micro-phytoplankton)                                             | 0.004         | mol Chl (mol C)<sup>-1</sup>                                         |
| `diamaxqc`         | Max Chl:C (micro-phytoplankton)                                             | 0.060         | mol Chl (mol C)<sup>-1</sup>                                         |
| `diaoptqf`         | Optimal Fe:C (micro-phytoplankton)                                          | 10e-6         | mol Fe (mol C)<sup>-1</sup>                                          |
| `diamaxqf`         | Max Fe:C (micro-phytoplankton)                                              | 65e-6         | mol Fe (mol C)<sup>-1</sup>                                          |
| `diaminqs`         | Min Si:C (micro-phytoplankton)                                              | 0.04          | mol Si (mol C)<sup>-1</sup>                                          |
| `diaoptqs`         | Optimal Si:C (micro-phytoplankton)                                          | 0.13          | mol Si (mol C)<sup>-1</sup>                                          |
| `diamaxqs`         | Max Si:C (micro-phytoplankton)                                              | 0.60          | mol Si (mol C)<sup>-1</sup>                                          |
| `diaVmaxs`         | Max Si uptake (micro-phytoplankton)                                         | 0.1/86400.0   | mol Si (mol C)<sup>-1</sup> s<sup>-1</sup>                           |
| `dialmor`          | Linear mortality rate (micro-phytoplankton)                                 | 0.001/86400.0 | s<sup>-1</sup>                                                       |
| `diaqmor`          | Quadratic mortality rate (micro-phytoplankton)                              | 0.05/86400.0  | (mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>                  |
| `diabiot`          | Biomass threshold (micro-phytoplankton)                                     | 0.5           | mmol C m<sup>-3</sup>                                                |
| `alphabio_tri`     | Initial slope of P–I curve (trichodesmium)                                  | 1.8           | mol C (mol Chl)<sup>-1</sup> (W m<sup>-2</sup>)<sup>-1</sup>         |
| `trikf`            | Fe half-saturation coefficient (trichodesmium)                              | 0.125         | µmol Fe m<sup>-3</sup>                                               |
| `trichlc`          | Chl:C (trichodesmium)                                                       | 0.01          | mol Chl (mol C)<sup>-1</sup>                                         |
| `trin2c`           | N:C (trichodesmium)                                                         | 50/300        | mol N (mol C)<sup>-1</sup>                                           |
| `chltau`           | Chlorophyll adjustment timescale                                            | 86400         | s                                                                    |
| `overflow`         | Max DOC exudation fraction by phytoplankton                                 | 0.75          | dimensionless                                                        |
| `bbioh`            | Heterotrophic growth scaling parameter b (Q10 = b^(10))                     | 1.072         | dimensionless                                                        |
| `zooCingest`       | Micro-zooplankton C ingestion efficiency                                    | 0.70          | mol C (mol C)<sup>-1</sup>                                           |
| `zooCassim`        | Micro-zooplankton C assimilation efficiency                                 | 0.40          | mol C (mol C)<sup>-1</sup>                                           |
| `zooFeingest`      | Micro-zooplankton Fe ingestion efficiency                                   | 0.06          | mol Fe (mol Fe)<sup>-1</sup>                                         |
| `zooFeassim`       | Micro-zooplankton Fe assimilation efficiency                                | 0.60          | mol Fe (mol Fe)<sup>-1</sup>                                         |
| `zooexcrdom`       | Micro-zooplankton excretion fraction routed to DOM                          | 0.70          | dimensionless                                                        |
| `zookz`            | Micro-zooplankton mortality half-saturation coefficient                     | 0.25          | mmol C m<sup>-3</sup>                                                |
| `zoogmax`          | Micro-zooplankton max grazing rate                                          | 3.3/86400.0   | s<sup>-1</sup>                                                       |
| `zooepsbacp`       | Micro-zooplankton prey capture efficiency (bacp)                            | 0.10/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `zooepsbacf1`      | Micro-zooplankton prey capture efficiency (bacf1)                           | 0.10/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `zooepsbacf2`      | Micro-zooplankton prey capture efficiency (bacf2)                           | 0.10/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `zooepsaoa`        | Micro-zooplankton prey capture efficiency (AOA)                             | 0.25/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `zooepsphy`        | Micro-zooplankton prey capture efficiency (nano-phytoplankton)              | 0.40/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `zooepsdia`        | Micro-zooplankton prey capture efficiency (micro-phytoplankton)             | 0.40/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `zooepsdet`        | Micro-zooplankton prey capture efficiency (small detritus)                  | 0.25/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `zprefbacp`        | Micro-zooplankton preference (bacp)                                         | 0.25          | dimensionless                                                        |
| `zprefbacf1`       | Micro-zooplankton preference (bacf1)                                        | 0.25          | dimensionless                                                        |
| `zprefbacf2`       | Micro-zooplankton preference (bacf2)                                        | 0.25          | dimensionless                                                        |
| `zprefaoa`         | Micro-zooplankton preference (AOA)                                          | 0.40          | dimensionless                                                        |
| `zprefphy`         | Micro-zooplankton preference (nano-phytoplankton)                           | 1.0           | dimensionless                                                        |
| `zprefdia`         | Micro-zooplankton preference (micro-phytoplankton)                          | 0.25          | dimensionless                                                        |
| `zprefdet`         | Micro-zooplankton preference (small detritus)                               | 0.80          | dimensionless                                                        |
| `zoolmor`          | Micro-zooplankton linear mortality rate                                     | 0.002/86400.0 | s<sup>-1</sup>                                                       |
| `zooqmor`          | Micro-zooplankton quadratic mortality rate                                  | 0.05/86400.0  | (mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>                  |
| `zoopreyswitch`    | Micro-zooplankton prey switching exponent                                   | 1.8           | dimensionless                                                        |
| `mesCingest`       | Meso-zooplankton C ingestion                                                | 0.75          | mol C (mol C)<sup>-1</sup>                                           |
| `mesCassim`        | Meso-zooplankton C assimilation                                             | 0.30          | mol C (mol C)<sup>-1</sup>                                           |
| `mesFeingest`      | Meso-zooplankton Fe ingestion                                               | 0.43          | mol Fe (mol Fe)<sup>-1</sup>                                         |
| `mesFeassim`       | Meso-zooplankton Fe assimilation                                            | 0.75          | mol Fe (mol Fe)<sup>-1</sup>                                         |
| `mesexcrdom`       | Meso-zooplankton excretion fraction routed to DOM                           | 0.35          | dimensionless                                                        |
| `mesgmax`          | Meso-zooplankton maximum grazing rate                                       | 1.0/86400.0   | s<sup>-1</sup>                                                       |
| `mesepsbacp`       | Meso-zooplankton prey capture efficiency (bacp)                             | 0.11/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `mesepsbacf1`      | Meso-zooplankton prey capture efficiency (bacf1)                            | 0.11/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `mesepsbacf2`      | Meso-zooplankton prey capture efficiency (bacf2)                            | 0.11/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `mesepsaoa`        | Meso-zooplankton prey capture efficiency (AOA)                              | 0.11/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `mesepsphy`        | Meso-zooplankton prey capture efficiency (nano-phytoplankton)               | 0.11/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `mesepsdia`        | Meso-zooplankton prey capture efficiency (micro-phytoplankton)              | 0.20/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `mesepsdet`        | Meso-zooplankton prey capture efficiency (small detritus)                   | 0.05/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `mesepsbdet`       | Meso-zooplankton prey capture efficiency (large detritus)                   | 0.10/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `mesepszoo`        | Meso-zooplankton prey capture efficiency (micro-zooplankton)                | 0.10/86400.0  | m<sup>6</sup> mmol<sup>-2</sup> s<sup>-1</sup>                       |
| `mprefbacp`        | Meso-zooplankton preference (bacp)                                          | 0.25          | dimensionless                                                        |
| `mprefbacf1`       | Meso-zooplankton preference (bacf1)                                         | 0.25          | dimensionless                                                        |
| `mprefbacf2`       | Meso-zooplankton preference (bacf2)                                         | 0.25          | dimensionless                                                        |
| `mprefaoa`         | Meso-zooplankton preference (AOA)                                           | 0.4           | dimensionless                                                        |
| `mprefphy`         | Meso-zooplankton preference (nano-phytoplankton)                            | 0.1           | dimensionless                                                        |
| `mprefdia`         | Meso-zooplankton preference (micro-phytoplankton)                           | 0.85          | dimensionless                                                        |
| `mprefdet`         | Meso-zooplankton preference (small detritus)                                | 0.80          | dimensionless                                                        |
| `mprefbdet`        | Meso-zooplankton preference (large detritus)                                | 0.80          | dimensionless                                                        |
| `mprefzoo`         | Meso-zooplankton preference (micro-zooplankton)                             | 0.85          | dimensionless                                                        |
| `meslmor`          | Meso-zooplankton linear mortality rate                                      | 0.002/86400.0 | s<sup>-1</sup>                                                       |
| `mesqmor`          | Meso-zooplankton quadratic mortality rate                                   | 0.75/86400.0  | (mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>                  |
| `mespreyswitch`    | Meso-zooplankton prey switching exponent                                    | 1.8           | dimensionless                                                        |
| `detlrem`          | Detritus hydrolysation rate                                                 | 0.7/86400.0   | (mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>                  |
| `detlrem_sed`      | Sediment detritus hydrolysation rate                                        | 0.005/86400.0 | s<sup>-1</sup>                                                       |
| `detphi`           | Porosity (small detritus)                                                   | 0.25          | dimensionless                                                        |
| `bdetphi`          | Porosity (large detritus)                                                   | 0.75          | dimensionless                                                        |
| `detrho`           | Detritus density                                                            | 1375          | kg m<sup>-3</sup>                                                    |
| `caco3rho`         | CaCO3 density                                                               | 2710          | kg m<sup>-3</sup>                                                    |
| `bsirho`           | Opal density                                                                | 2000          | kg m<sup>-3</sup>                                                    |
| `phyrad0`          | Nano-phytoplankton mean radius                                              | 10            | µm                                                                   |
| `diarad0`          | Micro-phytoplankton mean radius                                             | 50            | µm                                                                   |
| `zoorad0`          | Micro-zooplankton mean radius                                               | 30            | µm                                                                   |
| `mesrad0`          | Meso-zooplankton mean radius                                                | 1000          | µm                                                                   |
| `caco3lrem`        | CaCO<sub>3</sub> dissolution rate                                           | 0.01/86400.0  | s<sup>-1</sup>                                                       |
| `caco3lrem_sed`    | Sediment CaCO<sub>3</sub> dissolution rate                                  | 0.01/86400.0  | s<sup>-1</sup>                                                       |
| `f_inorg`          | Base inorganic fraction (PIC:POC ratio)                                     | 0.04          | mol CaCO3 (mol C)<sup>-1</sup>                                       |
| `disscal`          | Calcite dissolution rate                                                    | 0.10/86400.0  | s<sup>-1</sup>                                                       |
| `dissara`          | Aragonite dissolution rate                                                  | 0.10/86400.0  | s<sup>-1</sup>                                                       |
| `dissdet`          | Fraction CaCO3 dissolved per detritus hydrolyzed                            | 0.20          | mol CaCO3 (mol C)<sup>-1</sup>                                       |
| `fgutdiss`         | Zooplankton gut CaCO3 dissolution efficiency                                | 0.80          | dimensionless                                                        |
| `ligW`             | Weak ligand concentration                                                   | 1.7           | µmol m<sup>-3</sup>                                                  |
| `ligS`             | Strong ligand concentration                                                 | 0.4           | µmol m<sup>-3</sup>                                                  |
| `dfefloor`         | Minimum open water concentration of dissolved iron (detection limit)        | 0.001         | µmol Fe m<sup>-3</sup>                                               |
| `kscav_dfe`        | Free dissolved iron scavenging rate                                         | 0.01/86400.0  | (mmol mass of particle m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>   |
| `kcoag_dfe`        | Colloidal dissolved iron coagulation rate                                   | 1e-6/86400.0  | (mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>                  |
| `kagg_col`         | Colloidal dissolved iron aggregation rate                                   | 0.1/86400.0   | s<sup>-1</sup>                                                       |
| `kagg_kcol`        | Half-saturation coefficient for colloidal iron aggregation                  | 2.0           | µmol Fe m<sup>-3</sup>                                               |
| `kafe_dfe`         | Authigenic iron dissolution rate (small)                                    | 1e-4/86400    | s<sup>-1</sup>                                                       |
| `kbafe_dfe`        | Authigenic iron dissolution rate (large)                                    | 1e-4/86400    | s<sup>-1</sup>                                                       |
| `wafe`             | Authigenic iron sinking rate (small)                                        | 0.5/86400     | m s<sup>-1</sup>                                                     |
| `wbafe`            | Authigenic iron sinking rate (large)                                        | 5.0/86400     | m s<sup>-1</sup>                                                     |
| `bsi_fbac`         | Bacterial enhancement factor for silica dissolution                         | 20            | dimensionless                                                        |
| `bsi_kbac`         | Half-saturation coefficient for bacterial enhancement of silica dissolution | 0.5           | mmol C m<sup>-3</sup>                                                |
| `bsilrem_sed`      | Base sediment biogenic silica dissolution rate                              | 2.8e-8        | s<sup>-1</sup>                                                       |
| `aoa_knh4`         | AOA NH<sub>4</sub> half-saturation coefficient                              | 0.1           | mmol N m<sup>-3</sup>                                                |
| `aoa_poxy`         | AOA O<sub>2</sub> diffusive uptake limit                                    | 275/86400     | (mmol C biomass m<sup>3</sup>)<sup>-1</sup> s<sup>-1</sup>           |
| `aoa_ynh4`         | AOA NH<sub>4</sub> growth requirement                                       | 11            | mol NH<sub>4</sub> (mol C biomass)<sup>-1</sup>                      |
| `aoa_yoxy`         | AOA O<sub>2</sub> growth requirement                                        | 15.5          | mol O2 (mol C biomass)<sup>-1</sup>                                  |
| `aoa_en2omin`      | AOA minimum N<sub>2</sub>O yield per mol NH<sub>4</sub>                     | 0.0008        | mol N<sub>2</sub>O (mol NH<sub>4</sub>)<sup>-1</sup>                 |
| `aoa_C2N`          | AOA C:N ratio                                                               | 5             | mol C (mol N)<sup>-1</sup>                                           |
| `aoa_C2Fe`         | AOA C:Fe ratio                                                              | 1/(20e-6)     | mol C (mol Fe)<sup>-1</sup>                                          |
| `aoalmor`          | AOA linear mortality rate                                                   | 0.005/86400   | s<sup>-1</sup>                                                       |
| `aoaqmor`          | AOA quadratic mortality rate                                                | 0.001/86400   | (mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>                  |
| `bacp_Vmax_poc`    | Particle-associated bacteria POC uptake maximum                             | 6.7/86400     | mmol C m<sup>-3</sup> s<sup>-1</sup>                                 |
| `bacp_Vmax_no3`    | Particle-associated bacteria NO<sub>3</sub> uptake maximum                  | 7.2/86400     | mmol N m<sup>-3</sup> s<sup>-1</sup>                                 |
| `bacp_Vmax_dFe`    | Particle-associated bacteria dFe uptake maximum                             | 0.1/86400     | µmol Fe m<sup>-3</sup> s<sup>-1</sup>                                |
| `bacp_poxy`        | Particle-associated bacteria O<sub>2</sub> diffusive uptake limit           | 450/86400     | (mmol C biomass m<sup>3</sup>)<sup>-1</sup> s<sup>-1</sup>           |
| `bacp_kno3`        | Particle-associated bacteria NO<sub>3</sub> half-saturation coefficient     | 15            | mmol N m<sup>-3</sup>                                                |
| `bacp_kpoc`        | Particle-associated bacteria POC half-saturation coefficient                | 5             | mmol C m<sup>-3</sup>                                                |
| `bacp_kfer`        | Particle-associated bacteria dFe half-saturation coefficient                | 0.35          | µmol Fe m<sup>-3</sup>                                               |
| `bacp_alpha`       | Particle-associated bacteria degree of partial oxidation                    | 0.25          | mol C (mol C)<sup>-1</sup>                                           |
| `bacp_fele`        | Particle-associated bacteria fraction of electrons to biosynthesis          | 0.15          | e (e)<sup>-1</sup>                                                   |
| `bacf1_Vmax_doc`   | Free-living bacteria #1 DOC uptake maximum                                  | 6.7/86400     | mmol C m<sup>-3</sup> s<sup>-1</sup>                                 |
| `bacf1_Vmax_no3`   | Free-living bacteria #1 NO<sub>3</sub> uptake maximum                       | 7.2/86400     | mmol N m<sup>-3</sup> s<sup>-1</sup>                                 |
| `bacf1_Vmax_dFe`   | Free-living bacteria #1 dFe uptake maximum                                  | 0.1/86400     | µmol Fe m<sup>-3</sup> s<sup>-1</sup>                                |
| `bacf1_poxy`       | Free-living bacteria #1 O<sub>2</sub> diffusive uptake limit                | 450/86400     | (mmol C biomass m<sup>3</sup>)<sup>-1</sup> s<sup>-1</sup>           |
| `bacf1_kdoc`       | Free-living bacteria #1 DOC half-saturation coefficient                     | 60            | mmol C m<sup>-3</sup>                                                |
| `bacf1_kfer`       | Free-living bacteria #1 dFe half-saturation coefficient                     | 0.35          | µmol Fe m<sup>-3</sup>                                               |
| `bacf1_alpha`      | Free-living bacteria #1 degree of partial oxidation                         | 0.50          | mol C (mol C)<sup>-1</sup>                                           |
| `bacf1_fele`       | Free-living bacteria #1 fraction of electrons to biosynthesis               | 0.30          | e (e)<sup>-1</sup>                                                   |
| `bacf1_nosc_opt`   | Free-living bacteria #1 target oxidation state of organic carbon            | -0.5          | dimensionless                                                        |
| `bacf1_nosc_sig`   | Free-living bacteria #1 standard deviation of target oxidation states       | 1.0           | dimensionless                                                        |
| `bacf2_Vmax_doc`   | Free-living bacteria #2 DOC uptake maximum                                  | 6.7/86400     | mmol C m<sup>-3</sup> s<sup>-1</sup>                                 |
| `bacf2_Vmax_dFe`   | Free-living bacteria #2 dFe uptake maximum                                  | 0.1/86400     | µmol Fe m<sup>-3</sup> s<sup>-1</sup>                                |
| `bacf2_poxy`       | Free-living bacteria #2 O<sub>2</sub> diffusive uptake limit                | 450/86400     | (mmol C biomass m<sup>3</sup>)<sup>-1</sup> s<sup>-1</sup>           |
| `bacf2_pn2o`       | Free-living bacteria #2 N<sub>2</sub>O diffusive uptake limit               | 452/86400     | (mmol C biomass m<sup>3</sup>)<sup>-1</sup> s<sup>-1</sup>           |
| `bacf2_kdoc`       | Free-living bacteria #2 DOC half-saturation coefficient                     | 60            | mmol C m<sup>-3</sup>                                                |
| `bacf2_kfer`       | Free-living bacteria #2 dFe half-saturation coefficient                     | 0.35          | µmol Fe m<sup>-3</sup>                                               |
| `bacf2_alpha`      | Free-living bacteria #2 degree of partial oxidation                         | 0.25          | mol C (mol C)<sup>-1</sup>                                           |
| `bacf2_fele`       | Free-living bacteria #2 fraction of electrons to biosynthesis               | 0.15          | e (e)<sup>-1</sup>                                                   |
| `bacf2_nosc_opt`   | Free-living bacteria #2 target oxidation state of organic carbon            | 1.0           | dimensionless                                                        |
| `bacf2_nosc_sig`   | Free-living bacteria #2 standard deviation of target oxidation states       | 1.0           | dimensionless                                                        |
| `bac_C2N`          | Heterotrophic bacteria C:N                                                  | 5             | mol C (mol N)<sup>-1</sup>                                           |
| `bac_C2Fe`         | Heterotrophic bacteria C:Fe                                                 | 1/(40e-6)     | mol C (mol Fe)<sup>-1</sup>                                          |
| `baclmor`          | Heterotrophic bacteria linear mortality rate                                | 0.005/86400   | s<sup>-1</sup>                                                       |
| `bacqmor`          | Heterotrophic bacteria quadratic mortality rate                             | 0.05/86400    | (mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>                  |
| `aoxkn`            | Anammox NH<sub>4</sub> half-saturation coefficient                          | 0.5           | mmol N m<sup>-3</sup>                                                |
| `aoxmumax`         | Anammox maximum growth rate                                                 | 0.0025/86400  | s<sup>-1</sup>                                                       |
| `bottom_thickness` | Bottom layer thickness                                                      | 1.0           | m                                                                    |


---


### 1. Light attenuation through the water column.

Photosynthetically available radiation (PAR) is split into blue, green and red wavelengths. The incoming visible (photosynthetically available) short wave radiation flux (PAR, [W m<sup>-2</sup>]) is received from the physical model, and is then split evenly into each of blue, green and red light bands.

At the top (`par_bgr_top(k,b)`, $PAR^{top}$) and mid‑point (`par_bgr_mid(k,b)`, $PAR^{mid}$) of each layer `k` we calculate the downward irradiance by exponential decay of each band `b` through the layer thickness (`dzt(i,j,k)`, $\Delta z$, [m]) using band‑specific attenuation coefficients. These attenuation coefficients are related to the concentration of chlorophyll (`chl`, [mg m<sup>-3</sup>]), organic detritus (`ndet`, [mg N m<sup>-3</sup>]) and calcium carbonate (`carb`, [kg m<sup>-3</sup>]) in the water column. 

For chlorophyll, attenuation coefficients for each of blue, green and red light (`zbgr(ichl,b)`, [m<sup>-1</sup>]) are retrieved from the look-up table of [Morel & Maritorena (2001)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2000JC000319) (their Table 2) that explicitly relates chlorophyll concentration to attenuation rates and accounts for the packaging effect of chlorophyll in larger cells. Within `zbgr(ichl,b)`, `ichl` is an integer that corresponds to a particular band of chlorophyll concentration, with increasing chlorophyll concentrations associated with increasing attenuation.  

For organic detritus, attenuation coefficients for blue, green and red light (`dbgr(b)`, [(mg N m<sup>-3</sup>)<sup>-1</sup> m<sup>-1</sup>]) are taken from [Dutkiewicz et al. (2015)](https://bg.copernicus.org/articles/12/4447/2015/bg-12-4447-2015.html) (their Fig. 1b), while for calcium carbonate (`cbgr(b)`, [(kg CaCO<sub>3</sub> m<sup>-3</sup>)<sup>-1</sup>m<sup>-1</sup>]) we take the coefficients defined in [Soja-Wozniak et al. (2019)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JC014998). For both detritus and calcium carbonate, these studies provide concentration-normalized attenuation coefficients, which must be multiplied against concentrations to retrieve the correct units of [m<sup>-1</sup>].

Because WOMBAT-mid has two forms of phytoplankton (nanophytoplankton and microphytoplankton) with their own chlorophyll quotas and two forms of particulate detritus (small and large), we sum both chlorophyll pools and particulate detritus pools to return the total chlorophyll and the total particulate detritus. 

As an example, the PAR in the blue band (`b=1`) at the top of level k is computed as

$$
\begin{align}
PAR^{top}(k,1) =& \quad PAR^{top}(k-1,1) * e^{(-ex_{bgr}(k-1,1) * \Delta z(k-1))}
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
PAR^{mid}(k,3) =& \quad PAR^{mid}(k-1,3) * e^{(-0.5*(ex_{bgr}(k-1,3) * \Delta z(k-1) + ex_{bgr}(k,3) * \Delta z(k)))}
\end{align}
$$

_where_ <br>
- $PAR^{mid}(k-1,3)$ is the red light (`b=3`) at the mid-point of the overlying grid cell (`par_bgr_mid(k-1,3)`, [W m<sup>-2</sup>]) <br>
- $ex_{bgr}(k-1,3)$ is the total attenuation of red light (`b=3`) in the overlying grid cell (`ek_bgr(k-1,3)`, [m<sup>-1</sup>]) <br>
- $ex_{bgr}(k,3)$ is the total attenuation of red light (`b=3`) in the current grid cell (`ek_bgr(k,3)`, [m<sup>-1</sup>]) <br>
- $\Delta z(k-1)$ and $\Delta z(k)$ are the grid cell thicknesses of the overlying and current grid cells (`dzt(i,j,k)`, [m]) <br>

The total PAR available to phytoplantkon is assumed to be the sum of the blue, green and red bands. Because we assume that phytoplankton are homogenously distributed within a layer `k`, but we do not assume that light is homogenously distributed within that layer, we solve for the PAR that is seen by the average phytoplankton within that cell (`radbio`, $PAR$, [W m<sup>-2</sup>])

$$
\begin{align}
PAR(k) =& \quad \sum_{b=1}^3 \dfrac{(PAR^{top}(k,b) - PAR^{top}(k+1,b))}{ex_{bgr}(k,b) * \Delta z(k)}
\end{align}
$$

_where_ <br>
- $PAR^{top}(k,b)$ is the incoming photosynthetically active radiation at the top of grid cell `k` and light band `b` (`par_bgr_top(k,b)`, [W m<sup>-2</sup>]) <br>
- $ex_{bgr}(k,b)$ is the attenuation rate of light band `b` in grid cell `k` (`ek_bgr(k,b)`, [m<sup>-1</sup>]) <br>
- $\Delta z(k)$ is the grid cell thickness of grid cell `k` (`dzt(i,j,k)`, [m]) <br>

This ensures phytoplankton growth in the model responds to the mean light they experience in the cell, not just light at one point. See Eq. 19 from [Baird et al. (2020)](https://gmd.copernicus.org/articles/13/4503/2020/).

The euphotic depth (`zeuphot(i,j)`, [m]) is defined as the depth where `radbio` falls below the 1% threshold of incidient shortwave radiation or below 0.01 W m<sup>-2</sup>, whichever is shallower.

---

### 2. Nutrient limitation of phytoplankton.

At the start of each vertical loop `k=1` through `k=kmax` the code computes the biomass of nano-phytoplankton (`biophy`, $B_{np}$, [mmol C m<sup>-3</sup>]) and micro-phytoplankton (`biodia`, $B_{mp}$, [mmol C m<sup>-3</sup>]). Phytoplankton biomass is used to scale how nitrogen in the form of nitrate (`biono3`, NO<sub>3</sub>, [mmol N m<sup>-3</sup>]) and ammonium (`bionh4`, NH<sub>4</sub>, [mmol N m<sup>-3</sup>]), dissolved iron (`biofer`, $dFe$, [µmol dFe m<sup>-3</sup>]) and silicic acid in the case of micro-phytoplankton (`biosil`, H<sub>4</sub>SiO<sub>4</sub>, [mmol S m<sup>-3</sup>]) affect the growth of phytoplankton. Using compilations of marine phytoplankton and zooplankton communities, [Wickman et al. (2024)](https://www.science.org/doi/10.1126/science.adk6901) show that the nutrient affinity, $aff$, of a phytoplankton cell is related to its volume, $V$, via

$$
\begin{align}
aff =& \quad V^{-0.57}
\end{align}
$$

Additionally, the authors demonstrate that the volume of the average phytoplankton cell is related to the density (i.e., concentration) of phytoplankton via

$$
\begin{align}
V =& \quad (B_{phy})^{0.65}
\end{align}
$$
 
when combining panels c and f of their Figure 1. This then relates the affinity of an average cell to the concentration of phytoplankton biomass as

$$
\begin{align}
aff =& \quad (B_{phy})^{-0.37}
\end{align}
$$

With this information, we allow the half-saturation terms for nitrogen (`phy_kni(i,j,k)`, $K_{np}^{N}$, [mmol N m<sup>-3</sup>]; `dia_kni(i,j,k)`, $K_{mp}^{N}$, [mmol N m<sup>-3</sup>]), dissolved iron  (`phy_kfe(i,j,k)`, $K_{np}^{Fe}$, [µmol dFe m<sup>-3</sup>]; `dia_kfe(i,j,k)`, $K_{mp}^{Fe}$, [µmol dFe m<sup>-3</sup>]) and silicic acid (`dia_ksi(i,j,k)`, $K_{mp}^{Si}$, [mmol Si m<sup>-3</sup>]) uptake to vary as a function of phytoplankton biomass concentration. We set reference values for the half-saturation coefficient of nitrogen (`phykn`, $K_{np}^{N,0}$, [mmol N m<sup>-3</sup>]; `diakn`, $K_{mp}^{N,0}$, [mmol N m<sup>-3</sup>]), dissolved iron (`phykf`, $K_{np}^{Fe,0}$, [µmol dFe m<sup>-3</sup>]; `diakf`, $K_{mp}^{Fe,0}$, [µmol dFe m<sup>-3</sup>]) and silicic acid (`diaks`, $K_{mp}^{Si,0}$, [mmol Si m<sup>-3</sup>]) as input parameters to the model, and also set thresholds of nano-phytoplankton concentration (`phybiot`, $B_{np}^{thresh}$, [mmol C m<sup>-3</sup>]) and micro-phytoplankton concentration (`diabiot`, $B_{mp}^{thresh}$, [mmol C m<sup>-3</sup>]) beneath which cell size cannot decrease and affinity can no longer increase. At this minimum, where affinity is maximised, the half-saturation coefficients are bounded to be 10% of their reference values.

$$
\begin{align}
K_{np}^{N} =& \quad K_{np}^{N,0} * \max(0.1, \max(0.0, (B_{np}-B_{np}^{thresh}))^{0.37} ) \\
K_{np}^{Fe} =& \quad K_{np}^{Fe,0} * \max(0.1, \max(0.0, (B_{np}-B_{np}^{thresh}))^{0.37} ) \\
K_{mp}^{N} =& \quad K_{mp}^{N,0} * \max(0.1, \max(0.0, (B_{mp}-B_{mp}^{thresh}))^{0.37} ) \\
K_{mp}^{Fe} =& \quad K_{mp}^{Fe,0} * \max(0.1, \max(0.0, (B_{mp}-B_{mp}^{thresh}))^{0.37} ) \\
K_{mp}^{Si} =& \quad K_{mp}^{Si,0} * \max(0.1, \max(0.0, (B_{mp}-B_{mp}^{thresh}))^{0.37} )
\end{align}
$$

_where_ <br>
- $K_{np}^{N}$ and $K_{mp}^{N}$ are the half-saturation coefficients for nitrogen uptake by nano- and micro-phytoplankton (`phy_kni(i,j,k)` and `dia_kni(i,j,k)`, [mmol N m<sup>-3</sup>]) <br>
- $K_{np}^{Fe}$ and $K_{mp}^{Fe}$ are the half-saturation coefficients for iron uptake by nano- and micro-phytoplankton (`phy_kfe(i,j,k)` and `dia_kfe(i,j,k)`, [µmol Fe m<sup>-3</sup>]) <br>
- $K_{mp}^{Si}$ is the half-saturation coefficient of silicic acid uptake by micro-phytoplankton (`dia_ksi(i,j,k)`, [mmol Si m<sup>-3</sup>]) <br>


**Limitation of phytoplankton growth by nitrogen** (`phy_lnit(i,j,k)`, $L_{np}^{N}$), [dimensionless]; `dia_lnit(i,j,k)`, $L_{mp}^{N}$), [dimensionless]) is split between ammonium (`phy_lnh4(i,j,k)`, $L_{np}^{NH_4}$), [dimensionless]; `dia_lnh4(i,j,k)`, $L_{mp}^{NH_4}$), [dimensionless]) and nitrate (`phy_lno3(i,j,k)`, $L_{np}^{NO_3}$), [dimensionless]; `dia_lno3(i,j,k)`, $L_{mp}^{NO_3}$), [dimensionless]). Phytoplankton preferentially consume and grow on ammonium because it is most efficiently converted to glutamate for biomass synthesis, while nitrate must be first reduced within the cell ([Dortch, 1990](https://www.jstor.org/stable/24842258)). To represent this preference, we follow [Buchanan et al., 2025](https://bg.copernicus.org/articles/22/4865/2025/) who assert a 5-fold preference of phytoplankton for ammonium over nitrate and show that this reproduces preferences of ammonium-fueled growth in ocean field data.

$$
\begin{align}
l_{np}^{NH_4} =& \quad \dfrac{NH_4}{NH_4 + K_{np}^{N}} \\
l_{np}^{NO_3} =& \quad \dfrac{NO_3}{NO_3 + K_{np}^{N}} \\
l_{np}^{N} =& \quad \dfrac{NH_4 + NO_3}{NH_4 + NO_3 + K_{np}^{N}} \\
L_{np}^{NH_4} =& \quad \dfrac{5 \cdot l_{np}^{N} l_{np}^{NH_4}}{l_{np}^{NO_3} + 5 \cdot l_{np}^{NH_4}} \\
L_{np}^{NO_3} =& \quad \dfrac{l_{np}^{N} l_{np}^{NO_3}}{l_{np}^{NO_3} + 5 \cdot l_{np}^{NH_4}} \\
L_{np}^{N} =& \quad L_{np}^{NH_4} + L_{np}^{NO_3}
\end{align}
$$

_where_ <br>
- NH<sub>4</sub> is the in situ concentration of ammonium (`bionh4`, [mmol N m<sup>-3</sup>])  <br>
- NO<sub>3</sub> is the in situ concentration of nitrate (`biono3`, [mmol N m<sup>-3</sup>])  <br>
- $l_{np}^{NH_4}$ is the limitation term of nano-phytoplankton growth on ammonium before preferencing (`phy_limnh4`, [dimensionless])  <br>
- $l_{np}^{NO_3}$ is the limitation term of nano-phytoplankton growth on nitrate before preferencing (`phy_limno3`, [dimensionless])  <br>
- $l_{np}^{N}$ is the limitation term of nano-phytoplankton growth on nitrogen before preferencing (`phy_limdin`, [dimensionless])  <br>
- $L_{np}^{NH_4}$ is the limitation term of nano-phytoplankton growth on ammonium (`phy_lnh4(i,j,k)`, [dimensionless])  <br>
- $L_{np}^{NO_3}$ is the limitation term of nano-phytoplankton growth on nitrate (`phy_lno3(i,j,k)`, [dimensionless])  <br>
- $L_{np}^{N}$ is the limitation term of nano-phytoplankton growth on nitrogen (`phy_lnit(i,j,k)`, [dimensionless])  <br>


The same set of equations are applied to micro-phytoplankton:

$$
\begin{align}
l_{mp}^{NH_4} =& \quad \dfrac{NH_4}{NH_4 + K_{mp}^{N}} \\
l_{mp}^{NO_3} =& \quad \dfrac{NO_3}{NO_3 + K_{mp}^{N}} \\
l_{mp}^{N} =& \quad \dfrac{NH_4 + NO_3}{NH_4 + NO_3 + K_{mp}^{N}} \\
L_{mp}^{NH_4} =& \quad \dfrac{5 \cdot l_{mp}^{N} l_{mp}^{NH_4}}{l_{mp}^{NO_3} + 5 \cdot l_{mp}^{NH_4}} \\
L_{mp}^{NO_3} =& \quad \dfrac{l_{mp}^{N} l_{mp}^{NO_3}}{l_{mp}^{NO_3} + 5 \cdot l_{mp}^{NH_4}} \\
L_{mp}^{N} =& \quad L_{mp}^{NH_4} + L_{mp}^{NO_3}
\end{align}
$$

_where_ <br>
- NH<sub>4</sub> is the in situ concentration of ammonium (`bionh4`, [mmol N m<sup>-3</sup>])  <br>
- NO<sub>3</sub> is the in situ concentration of nitrate (`biono3`, [mmol N m<sup>-3</sup>])  <br>
- $l_{mp}^{NH_4}$ is the limitation term of micro-phytoplankton growth on ammonium before preferencing (`dia_limnh4`, [dimensionless])  <br>
- $l_{mp}^{NO_3}$ is the limitation term of micro-phytoplankton growth on nitrate before preferencing (`dia_limno3`, [dimensionless])  <br>
- $l_{mp}^{N}$ is the limitation term of micro-phytoplankton growth on nitrogen before preferencing (`dia_limdin`, [dimensionless])  <br>
- $L_{mp}^{NH_4}$ is the limitation term of micro-phytoplankton growth on ammonium (`dia_lnh4(i,j,k)`, [dimensionless])  <br>
- $L_{mp}^{NO_3}$ is the limitation term of micro-phytoplankton growth on nitrate (`dia_lno3(i,j,k)`, [dimensionless])  <br>
- $L_{mp}^{N}$ is the limitation term of micro-phytoplankton growth on nitrogen (`dia_lnit(i,j,k)`, [dimensionless])  <br>


Note that although phytoplankton prefer NH<sub>4</sub> over NO<sub>3</sub>, as NO<sub>3</sub> becomes more abundant than NH<sub>4</sub> the $L_{mp}^{NO_3}$ term begins to exceed the $L_{mp}^{NH_4}$ term such that phytoplankton switch from regenerated production (NH<sub>4</sub>-based) to new production (NO<sub>3</sub>-based). This reproduces the known switch of phytoplankton from regenerated to new production that is observed in the real ocean ([Dugdale & Goering, 1967](https://doi.org/10.4319/lo.1967.12.2.0196), [Buchanan et al., 2025](https://bg.copernicus.org/articles/22/4865/2025/)). Furthermore, if $K_{mp}^{N}$ > $K_{np}^{N}$, this ensures that (i) micro-phytoplankton are less competitive for NH<sub>4</sub> than nano-phytoplankton at any concentration and (ii) micro-phytoplankton growth is greater than nano-phytoplankton under abundant NO<sub>3</sub>, which is consistent with theory and observations ([Fawcett et al., 2011](https://doi.org/10.1038/ngeo1265), [Glibert et al., 2016](https://doi.org/10.1002/lno.10203))

**Limitation of phytoplankton growth by iron** follows an internal quota approach ([Droop, 1983](https://www.degruyterbrill.com/document/doi/10.1515/botm.1983.26.3.99/html)). Phytoplankton have a minimum iron quota (`phy_minqfe`, $Q_{np}^{-Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]; `dia_minqfe`, $Q_{mp}^{-Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]) and an optimal quota for growth (`phyoptqf`, $Q_{np}^{*Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]; `diaoptqf`, $Q_{mp}^{*Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]). The minimum iron quota, $Q_{np}^{-Fe:C}$ and $Q_{mp}^{-Fe:C}$, is dependent on three terms that each correspond to the iron required by photosystems, respiration and nitrate reduction ([Flynn & Hipkin, 1999](https://onlinelibrary.wiley.com/doi/10.1046/j.1529-8817.1999.3561171.x)):

$$
\begin{align}
Q_{np}^{-Fe:C} =& \quad 0.00167 / 55.85 \cdot \max\left( Q_{np}^{Chl:C}, Q_{np}^{-Chl:C} \right) \cdot 12 \\
                & + 1.21 \times 10^{-5} \cdot 14.0 / 55.85 / 7.625 \cdot 0.5 \cdot 1.5 \cdot L_{np}^{N} \\
                & + 1.15 \times 10^{-4} \cdot 14.0 / 55.85 / 7.625 \cdot 0.5 \cdot L_{np}^{NO_3}
\end{align}
$$

$$
\begin{align}
Q_{mp}^{-Fe:C} =& \quad 0.00167 / 55.85 \cdot \max\left( Q_{mp}^{Chl:C}, Q_{mp}^{-Chl:C} \right) \cdot 12 \\
                & + 1.21 \times 10^{-5} \cdot 14.0 / 55.85 / 7.625 \cdot 0.5 \cdot 1.5 \cdot L_{mp}^{N} \\
                & + 1.15 \times 10^{-4} \cdot 14.0 / 55.85 / 7.625 \cdot 0.5 \cdot L_{mp}^{NO_3}
\end{align}
$$

The first term reflects the amount of iron required for photosystems I and II. `0.00167/55.85` is equivalent to the grams of Fe per gram of chlorophyll divided by the grams of Fe per mol Fe, giving mol Fe per gram chlorophyll. This term is multipled by the chlorophyll to carbon ratio of the phytoplantkon cell (`phy_chlc`, $Q_{np}^{Chl:C}$, [mol C (mol C)<sup>-1</sup>]; `dia_chlc`, $Q_{mp}^{Chl:C}$, [mol C (mol C)<sup>-1</sup>]) and grams of C per mol C, returning mol Fe per mol C. At a healthy chlorophyll:C ratio of 0.03, this term returns an Fe:C ratio of roughly 10 µmol:mol, which reproduces well known requirements of phytoplankton cells ([Morel, Rueter & Price, 1991](https://www.jstor.org/stable/43924569)). The second term, representing the respiratory iron requirement, is derived from [Flynn & Hipkin (1999)](https://onlinelibrary.wiley.com/doi/10.1046/j.1529-8817.1999.3561171.x) who estimated 1.21 $\times 10^{-5}$ grams Fe per gram N assimilated into the cell, which is converted to mol Fe per mol C with 14 g N per mol N divided by 55.85 g Fe per mol Fe $\times$ 7.625 mol C per mol N. This second term assumes that respiration is reduced as growth becomes more limited by available nitrogen (`phy_lnit(i,j,k)`, $L_{np}^{N}$, [dimensionless]; `dia_lnit(i,j,k)`, $L_{mp}^{N}$, [dimensionless]). Finally, the third term represents the iron required by nitrate/nitrite reduction. Nitrate assimilation requires roughly 1.8-fold more iron than ammonia assimilation ([Raven, 1988](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/j.1469-8137.1988.tb04196.x)). [Flynn & Hipkin (1999)](https://onlinelibrary.wiley.com/doi/10.1046/j.1529-8817.1999.3561171.x) estimated a demand of 1.15 $\times 10^{-4}$ g Fe per mol NO$_3$ reduced, which is accounted for by the nitrate limitation term (`phy_lno3(i,j,k)`, $L_{np}^{NO_3}$, [dimensionless]; `dia_lno3(i,j,k)`, $L_{mp}^{NO_3}$, [dimensionless])). Note that the `1.5` is designed to account for dark respiration (i.e., respiration when the cells are not growing) and the `0.5` refers to the fact that during cell division the cell must reinstate half of its Fe reserves. 

The Fe limitation factor (`phy_lfer(i,j,k)`, $L_{np}^{Fe}$, [dimensionless]; `dia_lfer(i,j,k)`, $L_{mp}^{Fe}$, [dimensionless]) is then computed from the present Fe:C quota of the phytoplankton cells (`phy_Fe2C`, $Q_{np}^{Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]; `dia_Fe2C`, $Q_{mp}^{Fe:C}$, [mol Fe (mol C)<sup>-1</sup>]) relative to the minimum and optimal quotas.

$$
\begin{align}
L_{np}^{Fe} =& \quad \max\left(0.0, \min\left(1.0, \dfrac{ Q_{np}^{Fe:C} - Q_{np}^{-Fe:C} }{Q_{np}^{*Fe:C}} \right)\right)
\end{align}
$$

_where_ <br>
- $Q_{np}^{-Fe:C}$ is the minimum Fe:C quota of the nano-phytoplankton cell (`phy_minqfe`, [mol Fe (mol C)<sup>-1</sup>]) <br>
- $Q_{np}^{*Fe:C}$ is the optimal Fe:C quota of the nano-phytoplankton cell (`phyoptqf`, [mol Fe (mol C)<sup>-1</sup>]) <br>
- $Q_{np}^{Fe:C}$ is the in situ Fe:C quota of the nano-phytoplankton cell (`phy_Fe2C`, [mol Fe (mol C)<sup>-1</sup>]) <br>

$$
\begin{align}
L_{mp}^{Fe} =& \quad \max\left(0.0, \min\left(1.0, \dfrac{ Q_{mp}^{Fe:C} - Q_{mp}^{-Fe:C} }{Q_{mp}^{*Fe:C}} \right)\right)
\end{align}
$$

_where_ <br>
- $Q_{mp}^{-Fe:C}$ is the minimum Fe:C quota of the micro-phytoplankton cell (`dia_minqfe`, [mol Fe (mol C)<sup>-1</sup>]) <br>
- $Q_{mp}^{*Fe:C}$ is the optimal Fe:C quota of the micro-phytoplankton cell (`diaoptqf`, [mol Fe (mol C)<sup>-1</sup>]) <br>
- $Q_{mp}^{Fe:C}$ is the in situ Fe:C quota of the micro-phytoplankton cell (`dia_Fe2C`, [mol Fe (mol C)<sup>-1</sup>]) <br>

If the cell is Fe‑replete with a quota that exceeds the minimum quota by as much as the optimal quota, then Fe does not limit growth ($L_{np}^{Fe}$ = 1; $L_{mp}^{Fe}$ = 1). If the cell is Fe‑deplete with a quota equal to or less than the minimum quota, then the growth rate is reduced to zero. The optimal quota ($Q_{np}^{*Fe:C}$; $Q_{mp}^{*Fe:C}$) is therefore a measure of how much excess Fe is required to allow unrestricted growth.

**Limitation of micro-phytoplankton growth by silicic acid** is computed as a gating constraint on division via:

$$
\begin{align}
L_{mp}^{Si} =& \quad \min\left( 1.0, \max\left( 0.0, \dfrac{ Q_{mp}^{Si:C} - Q_{mp}^{-Si:C} }{ Q_{mp}^{*Si:C} - Q_{mp}^{-Si:C} }  \right) \right)
\end{align}
$$

_where_ <br>
- $Q_{mp}^{-Si:C}$ is the minimum Si:C quota of the micro-phytoplankton cell (`diaminqs`, [mol Si (mol C)<sup>-1</sup>]) <br>
- $Q_{mp}^{*Si:C}$ is the optimal Si:C quota of the micro-phytoplankton cell (`diaoptqs`, [mol Si (mol C)<sup>-1</sup>]) <br>

This formulation treats silicification as linearly limiting to growth between the minimum and optimal quotas. Above the optimal quota silica limitation does not exist. This reflects evidence that diatoms division is structurally constrained by silica until a threshold reserve is reached, at which point division can proceed ([Martin-Jézéquel, Hildebrand & Brzezinski, 2003](https://onlinelibrary.wiley.com/doi/full/10.1046/j.1529-8817.2000.00019.x)). This treatment is also supported by weak or even negative relationships between Si:C quotas and growth rates of marine diatoms ([María Mejía et al., 2013](https://www.sciencedirect.com/science/article/pii/S001670371300344X?via%3Dihub)) and is consistent with the apparent increase in Si:C quotas under Fe-limited growth ([Hutchins & Bruland, 1998](https://www.nature.com/articles/31203), [Takeda, 1998](https://www.nature.com/articles/31674)), which suggests that Si:C quotas can be decoupled from growth. 

---


### 3. Temperature-dependent metabolism.

**Autotrophy**

The maximum potential growth rate for nano-phytoplankton (`phy_mumax(i,j,k)`, $\mu_{np}^{max}$, [day<sup>-1</sup>]) and micro-phytoplankton (`dia_mumax(i,j,k)`, $\mu_{mp}^{max}$, [day<sup>-1</sup>]) is prescribed by the temperature-dependent Eppley curve ([Eppley, 1972](https://spo.nmfs.noaa.gov/content/temperature-and-phytoplankton-growth-sea)). This formulation scales a reference growth rate at 0ºC via a power-law scaling with temperature (`Temp(i,j,k)`, $T$, [ºC]).

$$
\begin{align}
\mu_{np}^{max} =& \quad \mu_{np}^{0^{\circ}C} \cdot (β_{np})^{T} \\
\mu_{mp}^{max} =& \quad \mu_{mp}^{0^{\circ}C} \cdot (β_{mp})^{T}
\end{align}
$$

_where_ <br>
- $\mu_{np}^{0^{\circ}}C$ is the rate of nano-phytoplankton growth at 0ºC (`abioa_phy`, [s<sup>-1</sup>]) <br>
- $β_{np}$ is the base temperature-sensitivity coefficient for autotrophy by nano-phytoplankton (`bbioa_phy`, [dimenionless]) <br>
- $\mu_{mp}^{0^{\circ}}C$ is the rate of micro-phytoplankton growth at 0ºC (`abioa_dia`, [s<sup>-1</sup>]) <br>
- $β_{mp}$ is the base temperature-sensitivity coefficient for autotrophy by micro-phytoplankton (`bbioa_dia`, [dimenionless]) <br>
- $T$ is in situ water temperature (`Temp(i,j,k)`, [ºC])

In the above, $\mu_{np}^{0ºC}$, $\mu_{mp}^{0ºC}$, $β_{np}$ and $β_{mp}$ are reference values input to the model at run time. This allows the user to configure nano-phytoplankton and micro-phytoplankton with different maximum potential growth rates and different sensitivities to temperature ([Anderson et al., 2021](https://www.nature.com/articles/s41467-021-26651-8)).

**Heterotrophy**

Heterotrophic processes include mortality of ecosystem functional types, grazing rates of zooplankton, growth rates of heterotrophic bacteria consuming organic matter (both particulate and dissolved) and the hydrolysation rate of particulate detritus in the sediments. These processes are scaled similarly to autotrophy, where some reference rate at 0ºC ($\mu_{het}^{0ºC}$, [<sup>s-1</sup>]) is multiplied by a power-law with temperature ($β_{hete}$). Each heterotrophic process has a different $\mu_{het}^{0ºC}$ value and we expand on this later under the mortality, grazing and bacterial heterotrophy sections. However, the basic formulation for scaling heterotrophic metabolisms with temperature takes the form:

$$
\begin{align}
\mu_{het} =& \quad \mu_{het}^{0ºC} \left(β_{hete}\right)^{T}
\end{align}
$$

_where_ <br>
- $\mu_{het}^{0ºC}$ is the rate of some heterotrophic metabolism at 0ºC ([s<sup>-1</sup>]) <br>
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless]) <br>
- $T$ is the in situ temperature of seawater (`Temp(i,j,k)`, [ºC])  <br>

In the code, the combined term $\left(β_{hete}\right)^{T}$ is saved as `fbc`. See sections below for further details on heterotrophic metabolisms.

---


### 4. Light limitation of phytoplankton

Phytoplankton growth is limited by light through a photosynthesis–irradiance (P–I) relationship that links cellular chlorophyll content and photosynthetically available radiation (`radbio`, $PAR$, [W m<sup>-2</sup>]).

First, The initial slope of the P–I curve, (`phy_pisl`, $\alpha_{np}$, [(W m<sup>-2</sup>)<sup>-1</sup>]; `dia_pisl`, $\alpha_{mp}$, [(W m<sup>-2</sup>)<sup>-1</sup>]), determines how efficiently phytoplankton convert light into carbon fixation. It is scaled by the cellular chlorophyll-to-carbon ratio (`phy_chlc`, $Q_{np}^{Chl:C}$, [mol C (mol C)<sup>-1</sup>]; `dia_chlc`, $Q_{mp}^{Chl:C}$, [mol C (mol C)<sup>-1</sup>]).

$$
\begin{align}
\alpha_{np} =& \quad \max(\alpha_{np}^{Chl} \cdot Q_{np}^{Chl:C} \ , \ \alpha_{np}^{Chl} \cdot Q_{np}^{-Chl:C}) \\
\alpha_{mp} =& \quad \max(\alpha_{mp}^{Chl} \cdot Q_{mp}^{Chl:C} \ , \ \alpha_{mp}^{Chl} \cdot Q_{mp}^{-Chl:C})
\end{align}
$$

_where_ <br>
- $\alpha_{np}^{Chl}$ is the photosynthetic efficiency per unit chlorophyll in nano-phytoplankton (`alphabio_phy`, [(W m<sup>-2</sup>)<sup>-1</sup> (mol C (mol C)<sup>-1</sup>)<sup>-1</sup>]) <br>
- $\alpha_{mp}^{Chl}$ is the photosynthetic efficiency per unit chlorophyll in micro-phytoplankton (`alphabio_dia`, [(W m<sup>-2</sup>)<sup>-1</sup> (mol C (mol C)<sup>-1</sup>)<sup>-1</sup>]) <br>
- $Q_{np}^{-Chl:C}$ is the minimum chlorophyll to carbon ratio of nano-phytoplankton cells (`phyminqc`, [mol C (mol C)<sup>-1</sup>]) <br>
- $Q_{mp}^{-Chl:C}$ is the minimum chlorophyll to carbon ratio of micro-phytoplankton cells (`diaminqc`, [mol C (mol C)<sup>-1</sup>]) <br>
- $Q_{np}^{Chl:C}$ is the in situ chlorophyll to carbon ratio of nano-phytoplankton cells (`phy_chlc`, [mol C (mol C)<sup>-1</sup>]) <br>
- $Q_{mp}^{Chl:C}$ is the in situ chlorophyll to carbon ratio of micro-phytoplankton cells (`dia_chlc`, [mol C (mol C)<sup>-1</sup>]) <br>

This constraint prevents photosynthesis from collapsing unrealistically at low chlorophyll concentrations. These values are parameter inputs at run time and can differ between nano-phytoplankton and micro-phytoplankton ([Edwards et al., 2015](https://doi.org/10.1002/lno.10033), [Litchman 2022](https://link.springer.com/chapter/10.1007/978-3-030-92499-7_1)).

Second, light limitation (`phy_lpar(i,j,k)`, $L_{np}^{PAR}$), [dimensionless]; `dia_lpar(i,j,k)`, $L_{mp}^{PAR}$), [dimensionless]) is calculated using an exponential P–I formulation.

$$
\begin{align}
L_{np}^{PAR} =& \quad 1 - e^{- \alpha_{np} PAR } \\
L_{mp}^{PAR} =& \quad 1 - e^{- \alpha_{mp} PAR }
\end{align}
$$

_where_ <br>
- $PAR$ is the downwelling photosynthetically available radiation (`radbio`, [W m<sup>-2</sup>]) <br>

At low irradiance ($PAR$), growth increases approximately linearly with light, while at high irradiance photosynthesis asymptotically saturates. We do not account for photoinhibition at very high irradiances.

---


### 5. Realized growth rate of phytoplankton.

Realized growth of nano-phytoplankton (`phy_mu(i,j,k)`, $\mu_{np}$, [s<sup>-1</sup>]) and micro-phytoplankton (`dia_mu(i,j,k)`, $\mu_{mp}$, [s<sup>-1</sup>]) is calculated as:

$$
\begin{align}
\mu_{np} =& \quad \mu_{np}^{max} L_{np}^{PAR} \min(L_{np}^{N}, L_{np}^{Fe}) \\
\mu_{mp} =& \quad \mu_{mp}^{max} L_{mp}^{PAR} \min(L_{mp}^{N}, L_{mp}^{Fe}) L_{mp}^{Si}
\end{align}
$$

_where_ <br>
- $\mu_{np}^{max}$ is the maximum potential rate of carbon fixation by nano-phytoplankton (`phy_mumax`, [s<sup>-1</sup>]) <br>
- $L_{np}^{PAR}$ is the growth limiter by light of nano-phytoplankton (`phy_lpar(i,j,k)`, [dimensionless]) <br>
- $L_{np}^{N}$ is the growth limiter by nitrogen of nano-phytoplankton (`phy_lnit(i,j,k)`, [dimensionless]) <br>
- $L_{np}^{Fe}$ is the growth limiter by iron of nano-phytoplankton (`phy_lfer(i,j,k)`, [dimensionless]) <br>
- $\mu_{mp}^{max}$ is the maximum potential rate of carbon fixation by micro-phytoplankton (`dia_mumax`, [s<sup>-1</sup>]) <br>
- $L_{mp}^{PAR}$ is the growth limiter by light of micro-phytoplankton (`dia_lpar(i,j,k)`, [dimensionless]) <br>
- $L_{mp}^{N}$ is the growth limiter by nitrogen of micro-phytoplankton (`dia_lnit(i,j,k)`, [dimensionless]) <br>
- $L_{mp}^{Fe}$ is the growth limiter by iron of micro-phytoplankton (`dia_lfer(i,j,k)`, [dimensionless]) <br>
- $L_{mp}^{Si}$ is the growth limiter by silicic acid of micro-phytoplankton (`dia_lsil(i,j,k)`, [dimensionless]) <br>

Liebig's law of the minimum ([Liebig, 1840](https://archive.org/details/organicchemistry00liebrich/mode/2up), [Blackman, 1905](https://doi.org/10.1093/oxfordjournals.aob.a089000)) is applied to resources that are required for biomass synthesis (N and Fe). For micro-phytoplankton, their growth is additionally restricted by silica limitation applied outside of Liebig's law because we treat silica limitation (`dia_lsil(i,j,k)`, $L_{mp}^{Si}$, [dimensionless]) as a structural threshold, rather than as a metabolic throttle (see below).

Carbon fixation by phytoplankton is then calculated as:

$$
\begin{align}
\mu_{np}^{\leftarrow C} =& \quad \mu_{np} B_{np}^{C} \\
\mu_{mp}^{\leftarrow C} =& \quad \mu_{mp} B_{mp}^{C}
\end{align}
$$

_where_ <br>
- $\mu_{np}^{\leftarrow C}$ is the realized rate of carbon biomass growth by nano-phytoplankton (`phygrow(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\mu_{mp}^{\leftarrow C}$ is the realized rate of carbon biomass growth by micro-phytoplankton (`diagrow(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $B_{np}^{C}$ is the in situ concentration of nano-phytoplankton biomass (`f_phy(i,j,k)`, [mol C kg<sup>-1</sup>]) <br>
- $B_{mp}^{C}$ is the in situ concentration of micro-phytoplankton biomass (`f_dia(i,j,k)`, [mol C kg<sup>-1</sup>]) <br>

---


### 6. Dissolved organic carbon release by phytoplankton.

We implement the overflow hypothesis ([Fogg, 1983](https://doi.org/10.1515/botm.1983.26.1.3); [Hansell & Carlson, 2014](https://books.google.com.au/books?id=7iKOAwAAQBAJ&lpg=PP1&ots=kzkdHuHMF_&dq=Carlson%20Hansell%202014%20doi&lr&pg=PP1#v=onepage&q&f=false)), which posits that phytoplankton can exude their assimilated carbon as dissolved organic carbon (DOC) in high light, low nutrient conditions. We thus account for a phytoplankton-mediated creation of DOC from dissolved inorganic carbon (DIC) via:

$$
\begin{align}
\mu_{np}^{\rightarrow DOC} =& \quad \min\left( f_{overflow} \mu_{np}^{totalC}, \max\left(0.02 \cdot \mu_{np}^{totalC}, \mu_{np}^{totalC} - \mu_{np}^{\leftarrow C}\right) \right) \\
\mu_{mp}^{\rightarrow DOC} =& \quad \min\left( f_{overflow} \mu_{mp}^{totalC}, \max\left(0.02 \cdot \mu_{mp}^{totalC}, \mu_{mp}^{totalC} - \mu_{mp}^{\leftarrow C}\right) \right)
\end{align}
$$

_where_ <br>
- $\mu_{np}^{\rightarrow DOC}$ is the overflow production of DOC by nano-phytoplankton (`phydoc(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\mu_{mp}^{\rightarrow DOC}$ is the overflow production of DOC by micro-phytoplankton (`diadoc(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\mu_{np}^{totalC}$ is the total carbon fixation rate of nano-phytoplankton (`zval`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\mu_{mp}^{totalC}$ is the total carbon fixation rate rate of micro-phytoplankton (`zval`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\mu_{np}^{\leftarrow C}$ is the realized biomass growth rate of nano-phytoplankton (`phygrow(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\mu_{mp}^{\leftarrow C}$ is the realized biomass growth rate of micro-phytoplankton (`diagrow(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $f_{overflow}$ is the maximum fraction total carbon fixation that goes to DOC exudation (`overflow`, [dimenionless])

The total carbon fixation rate of phytoplankton type $p$ is 

$$
\begin{align}
\mu_{p}^{totalC} =& \quad \mu_{p}^{\rightarrow DOC} + \mu_{p}^{\leftarrow C} = \mu_{p}^{max} L_{p}^{PAR} B_{p}^{C}
\end{align}
$$

This formulation is derived from the idea that DOC exudation occurs as a result of the difference between carbon fixation capacity, which is bounded by light, and biosynthesis, which is bounded by light and nutrient resources. Since [Thornton (2014)](https://doi.org/10.1080/09670262.2013.875596) identified that as much as 50% of total phytoplankton carbon fixation can be routed to DOC exudation, we cap DOC exudation at $f_{overflow}$ of total carbon fixation, which is set as to a default of 0.75. We also set a hard bound that 2% of total carbon fixation must at minimum go to DOC production based on the findings of [Bjørnsen (1988)](https://doi.org/10.4319/lo.1988.33.1.0151) who identified that even the healthiest cells lose a small fraction of their assimilated carbon as DOC via passive diffusion across the cell membrane.

---


### 7. Synthesis of chlorophyll

This step diagnoses the **rate of chlorophyll synthesis** as a function of mixed-layer light, the phytoplankton growth rate and nutrient availability. The structure is consistent with the [Geider, MacIntyre & Kana (1997)](https://doi.org/10.3354/meps148187) formulation that relaxes the chlorophyll-to-carbon ratio towards an optimal value that supports photosynthetic growth under prevailing light and nutrient conditions.

We first solve for the optimal chlorophyll-to-carbon ratio (`phy_chlc`, $Q_{np}^{*Chl:C}$, [mol C (mol C)<sup>-1</sup>]; `dia_chlc`, $Q_{mp}^{*Chl:C}$, [mol C (mol C)<sup>-1</sup>]), which is diagnosed as the ratio required to support maximal photosynthetic carbon fixation under the ambient mean light level in the mixed layer, while accounting for nutrient limitation of biosynthesis:

$$
\begin{align}
Q_{np}^{*Chl:C} =& \quad \dfrac{Q_{np}^{+Chl:C}}{1 + \dfrac{\alpha_{np} PAR_{MLD} Q_{np}^{+Chl:C}}{2 \cdot 86400 \cdot \mu_{np}^{max} \min \left(L_{np}^{N}, L_{np}^{Fe} \right) }} \\
Q_{mp}^{*Chl:C} =& \quad \dfrac{Q_{mp}^{+Chl:C}}{1 + \dfrac{\alpha_{mp} PAR_{MLD} Q_{mp}^{+Chl:C}}{2 \cdot 86400 \cdot \mu_{mp}^{max} \min \left(L_{mp}^{N}, L_{mp}^{Fe} \right) }}
\end{align}
$$

_where_ <br>
- $Q_{np}^{+Chl:C}$ and $Q_{mp}^{+Chl:C}$ are the maximum allowable chlorophyll-to-carbon ratios (`phymaxqc`; `diamaxqc`, [mol C (mol C)<sup>-1</sup>]) <br>
- $\alpha_{np}$ and $\alpha_{mp}$ are the chlorophyll-specific initial slopes of the P–I curve (`alphabio_phy`; `alphabio_dia`, [(W m<sup>-2</sup>)<sup>-1</sup> (mol C (mol C)<sup>-1</sup>)<sup>-1</sup>]) <br>
- $PAR_{MLD}$ is mean photosynthetically available radiation over the mixed layer (`radmld(i,j,k)`, [W m<sup>-2</sup>]) <br>
- $\mu_{np}^{max}$ and $\mu_{mp}^{max}$ are the temperature-dependent maximum phytoplankton growth rates (`phy_mumax(i,j,k)`; `dia_mumax(i,j,k)`, [s<sup>-1</sup>] ) <br>
- $L_{np}^{N}$ and $L_{np}^{Fe}$ are the nano-phytoplankton limitation factors for growth on N and Fe (`phy_lnit(i,j,k)`; `phy_lfer(i,j,k)`, [dimensionless]) <br>
- $L_{mp}^{N}$ and $L_{mp}^{Fe}$ are the micro-phytoplankton limitation factors for growth on N and Fe (`dia_lnit(i,j,k)`; `dia_lfer(i,j,k)`, [dimensionless]) <br>

We set a floor for the minimum chlorophyll-to-carbon ratio of phytoplankton via:

$$
\begin{align}
Q_{np}^{*Chl:C} =& \quad \max \left( Q_{np}^{*Chl:C}, Q_{np}^{-Chl:C} \right) \\
Q_{mp}^{*Chl:C} =& \quad \max \left( Q_{mp}^{*Chl:C}, Q_{mp}^{-Chl:C} \right)
\end{align}
$$

_where_ <br>
- $Q_{np}^{-Chl:C}$ and $Q_{mp}^{-Chl:C}$ are the minimum allowable chlorophyll-to-carbon ratios (`phyminqc`; `diaminqc`, [mol C (mol C)<sup>-1</sup>]) <br>

Synthesis of chlorophyll by nano-phytoplankton and micro-phytoplankton (`pchl_mu(i,j,k)`; `dchl_mu(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) is then calculated as:

$$
\begin{align}
\mu_{np}^{\leftarrow Chl} =& \quad \mu_{np} B_{np}^{Chl} + \dfrac{ Q_{np}^{*Chl:C} - Q_{np}^{Chl:C} }{\tau^{Chl}} \cdot B_{np}^{C} \\
\mu_{mp}^{\leftarrow Chl} =& \quad \mu_{mp} B_{mp}^{Chl} + \dfrac{ Q_{mp}^{*Chl:C} - Q_{mp}^{Chl:C} }{\tau^{Chl}} \cdot B_{mp}^{C}
\end{align}
$$

_where_ <br>
- $Q_{np}^{Chl:C}$ and $Q_{mp}^{Chl:C}$ are the in-situ chlorophyll-to-carbon ratios (`phy_chlc`; `dia_chlc`, [mol C (mol C)<sup>-1</sup>]) <br>
- $B_{np}^{Chl}$ and $B_{mp}^{Chl}$ are the in-situ concentations of phytoplankton chlorophyll (`f_pchl(i,j,k)`; `f_dchl(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $\mu_{np}$ and $\mu_{mp}$ are the realized growth rates of phytoplankton (`phy_mu(i,j,k)`; `dia_mu(i,j,k)`, [s<sup>-1</sup>] ) <br>
- $\tau^{Chl}$ is the timescale over which chlorophyll synthesis occurs within the cell (`chltau`, [s]) <br>
- $B_{np}^{C}$ and $B_{mp}^{C}$ are the in-situ concentations of phytoplankton carbon (`f_phy(i,j,k)`; `f_dia(i,j,k)`, [mol kg<sup>-1</sup>]) <br>

This formulation elevates chlorophyll-to-carbon ratios in low light and supresses synthesis when nutrients are low. $\tau^{Chl}$ is an input parameter at run time and should ideally be less than the doubling time of phytplankton given that phytoplankton can internally regulate their chlorophyll stores at rates greater than their overall growth.

---


### 8. Phytoplankton uptake of iron

Like chlorophyll, the iron content of phytoplankton is explicitly tracked as a tracer in WOMBAT-mid. First, a maximum quota is found based on the maximum Fe:C ratio of the phytoplankton type:

$$
\begin{align}
B_{np}^{+Fe} =& \quad B_{np}^{C} Q_{np}^{+Fe:C} \\
B_{mp}^{+Fe} =& \quad B_{mp}^{C} Q_{mp}^{+Fe:C}
\end{align}
$$

_where_ <br>
- $B_{np}^{+Fe}$ and $B_{mp}^{+Fe}$ are the maximum Fe quotas of the nano-phytoplankton and micro-phytoplankton cells (`phy_maxqfe`; `dia_maxqfe`, [mmol Fe m-3]) <br>
- $B_{np}^{C}$ and $B_{mp}^{C}$ are the in situ concentrations of nano-phytoplankton and micro-phytoplankton (`biophy`; `biodia`, [mmol C m<sup>-3</sup>]) <br>
- $Q_{np}^{+Fe:C}$ and $Q_{mp}^{+Fe:C}$ are the maximum Fe:C ratios of nano-phytoplankton and micro-phytoplankton cells (`phymaxqf`; `diamaxqf`, [mol Fe (mol C)-1]) <br>

Following [Aumont et al. (2015)](https://gmd.copernicus.org/articles/8/2465/2015/), this rate is scaled by three terms relating to (i) michaelis-menten type affinity for dFe, (ii) up-regulation of dFe uptake representing investment in transporters when cell quotas are limiting to growth, and (iii) down regulation of dFe uptake associated with enriched cellular quotas.

$$
\begin{align}
i_{np} =& \quad \dfrac{dFe}{dFe + K_{np}^{Fe}} \\
i_{mp} =& \quad \dfrac{dFe}{dFe + K_{mp}^{Fe}} \\
ii_{np} =& \quad 4 - \dfrac{4.5 L_{np}^{Fe}}{0.5 + L_{np}^{Fe}} \\
ii_{mp} =& \quad 4 - \dfrac{4.5 L_{mp}^{Fe}}{0.5 + L_{mp}^{Fe}} \\
iii_{np} =& \quad \max\left(0, 1 - \dfrac{B_{np}^{Fe} / B_{np}^{+Fe}}{\left|1.05 - B_{np}^{Fe} / B_{np}^{+Fe}\right|} \right) \\
iii_{mp} =& \quad \max\left(0, 1 - \dfrac{B_{mp}^{Fe} / B_{mp}^{+Fe}}{\left|1.05 - B_{mp}^{Fe} / B_{mp}^{+Fe}\right|} \right)
\end{align}
$$

_where_ <br>
- $dFe$ is the in situ dissolved iron concentration (`biofer`, [µmol Fe m<sup>-3</sup>]) <br>
- $K_{np}^{Fe}$ and $K_{mp}^{Fe}$ are the half-saturation coefficients for dFe uptake by nano-phytoplankton and micro-phytoplankton (`phy_kfe(i,j,k)`; `dia_kfe(i,j,k)`, [µmol Fe <sup>m-3</sup>]) <br>
- $L_{np}^{Fe}$ and $L_{mp}^{Fe}$ are the growth limiters of nano-phytoplankton and micro-phytoplankton by iron (`phy_lfer(i,j,k)`; `dia_lfer(i,j,k)`, [dimensionless]) <br>
- $B_{np}^{Fe}$ and $B_{mp}^{Fe}$ are the in situ Fe quotas of nano-phytoplankton and micro-phytoplankton cells (`biophyfe`; `biodiafe`, [mmol Fe <sup>m-3</sup>]) <br>
- $B_{np}^{+Fe}$ and $B_{mp}^{+Fe}$ are the maximum Fe quotas of nano-phytoplankton and micro-phytoplankton cells (`phy_maxqfe`; `dia_maxqfe`, [mmol Fe <sup>m-3</sup>]) <br>

Note that we additionally include a fourth term that decreases the maximum dFe uptake of a cell under light limitation. This is informed by slower uptake of Fe by cells grown in darkness compared to those grown in light by roughly 10-fold ([Strzepek et al., 2025](https://doi.org/10.1093/ismejo/wraf015)), which may be due to physiological stimulation of Fe uptake machinery or photoreduction of ligand-bound iron complexes ([Kong et al., 2023](https://doi.org/10.1002/lno.12331); [Maldonado et al., 2005](https://doi.org/10.1029/2005GB002481)), or possibly a combination of both. To obtain a 10-fold relative increase in Fe uptake rates under light, we applied the following term:

$$
\begin{align}
iv_{np} =& \quad \max\left(0.01, L_{np}^{PAR}\right)^{0.5} \\
iv_{mp} =& \quad \max\left(0.01, L_{mp}^{PAR}\right)^{0.5}
\end{align}
$$

_where_ <br>
- $L_{np}^{PAR}$ and $L_{mp}^{PAR}$ are the growth limiters of nano-phytoplankton and micro-phytoplankton by light (`phy_lpar(i,j,k)`; `dia_lpar(i,j,k)`, [dimensionless]) <br>

Under very low light, this fourth term reduces maximum potential Fe uptake by 10-fold than what it otherwise would be. All four terms are dimensionless and are designed to scale dissolved iron uptake either up or down. Dissolved iron uptake by nano-phytoplankton and micro-phytoplankton (`phy_dfeupt(i,j,k)`; `dia_dfeupt(i,j,k)`, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]) is then calculated as:

$$
\begin{align}
\mu_{np}^{\leftarrow dFe} =& \quad \mu_{np}^{max} B_{np}^{+Fe} \cdot (i_{np}) \cdot (ii_{np}) \cdot (iii_{np}) \cdot (iv_{np}) \\
\mu_{mp}^{\leftarrow dFe} =& \quad \mu_{mp}^{max} B_{mp}^{+Fe} \cdot (i_{mp}) \cdot (ii_{mp}) \cdot (iii_{mp}) \cdot (iv_{mp})
\end{align}
$$

_where_ <br>
- $\mu_{np}^{\leftarrow dFe}$ and $\mu_{mp}^{\leftarrow dFe}$ are the realized uptake rate of dissolved iron by nano-phytoplankton and micro-phytoplankton (`phy_dfeupt(i,j,k)`; `dia_dfeupt(i,j,k)`, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\mu_{np}^{max}$ and $\mu_{mp}^{max}$ are the maximum potential growth rates of nano-phytoplankton and micro-phytoplankton (`phy_mumax(i,j,k)`; `dia_mumax(i,j,k)`, [s<sup>-1</sup>]) <br>
- $B_{np}^{+Fe}$ and $B_{mp}^{+Fe}$ are the maximum Fe quotas of nano-phytoplankton and micro-phytoplankton cells (`phy_maxqfe`; `dia_maxqfe`, [mol Fe kg<sup>-1</sup>]) <br>

---


### 9. Phytoplankton uptake of silicic acid.

Like chlorophyll and iron, the silicon content of micro-phytoplankton is explicitly tracked as a tracer in WOMBAT-mid. Uptake of silicic acid by micro-phytoplankton (`dia_silupt(i,j,k)`, $\mu_{mp}^{\leftarrow Si}$, [mol Si kg<sup>-1</sup> s<sup>-1</sup>]) is scaled by two terms relating to (i) michaelis-menten type affinity for H<sub>4</sub>SiO<sub>4</sub> and (ii) down regulation of H<sub>4</sub>SiO<sub>4</sub> uptake associated with enriched cellular quotas. 

$$
\begin{align}
(i) & \quad \dfrac{H_{4}SiO_{4}}{H_{4}SiO_{4} + K_{mp}^{Si}} \\
(ii) & \quad \left(1.0 - \max\left(0.0, \dfrac{Q_{mp}^{Si:C} - Q_{mp}^{-Si:C}}{Q_{mp}^{+Si:C} - Q_{mp}^{-Si:C}} \right)\right)^{0.5}
\end{align}
$$

_where_ <br>
- H<sub>4</sub>SiO<sub>4</sub> is the in situ silicic acid concentration (`biosil`, [mmol Si m<sup>-3</sup>]) <br>
- $K_{mp}^{Si}$ is the half-saturation coefficient for siliic acid uptake by micro-phytoplankton (`dia_ksi(i,j,k)`, [mmol Si m<sup>-3</sup>]) <br>
- $Q_{mp}^{Si:C}$ is the in situ Si:C ratios of micro-phytoplankton cells (`dia_Si2C`, [mol Si (mol C)-1]) <br>
- $Q_{mp}^{+Si:C}$ is the maximum Si:C ratios of micro-phytoplankton cells (`diamaxqs`, [mol Si (mol C)-1]) <br>
- $Q_{mp}^{-Si:C}$ is the minimum Si:C ratios of micro-phytoplankton cells (`diaminqs`, [mol Si (mol C)-1]) <br>

Uptake is then calculated as

$$
\begin{align}
\mu_{mp}^{\leftarrow Si} =& \quad \max \left(0, V_{mp}^{Si} \cdot B_{mp}^{C} \cdot (i) \cdot (ii) \right)
\end{align}
$$

_where_ <br>
- $V_{mp}^{Si}$ is the maximum uptake rate of silicon to carbon by a micro-phytoplankton cell (`diaVmaxs`, [mol Si (mol C)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $B_{mp}^{C}$ is the in situ concentration of micro-phytoplankton carbon biomass (`f_dia(i,j,k)`, [mol C kg<sup>-1</sup>]) <br>

Unlike iron uptake, we do not include upregulation terms for silicic acid uptake. This is on the basis that highly silicified diatoms are caused by slow growth rather than increased/luxury uptake. Both light-limited and iron-limited diatoms show increases in their Si:C content by roughly 3-fold and it is suggested that this due to decoupling of biogenic silica precipitation from slowing carbon fixation ([Liu et al., 2016](https://doi.org/10.3389/fmars.2016.00089); [Hutchins & Bruland 1998](https://doi.org/10.1038/31203); [Takeda 1998](https://doi.org/10.1038/31674)). To properly decouple silicification from biomass growth we therefore make $V_{mp}^{Si}$ temperature-independent, which ensures that polar diatoms have a tendency towards heavier silicification then tropical diatoms ([Baines et al., 2010](https://doi.org/10.1029/2010GB003856)). 

---


### 10. Iron chemistry (scavenging, coagulation, dissolution).

Treatment of dissolved iron (`f_fe(i,j,k)`, $dFe$, mol kg<sup>-1</sup>) follows a combination of [Aumont et al. (2015)](https://gmd.copernicus.org/articles/8/2465/2015/) and [Tagliabue et al. (2023)](https://www.nature.com/articles/s41586-023-06210-5). Our calculations involve:
1. Solving for the distinct pools of dissolved iron: free iron, ligand-bound iron and colloidal iron.
2. Computing scavenging of free iron to authigenic sinking phases.
3. Computing coagulation of colloidal iron to authigenic sinking phases.
4. Computing dissolution of authigenic sinking phases back to dissolved iron.

We first estimate the **solubility of free Fe from Fe<sup>3+</sup>** in solution using temperature, pH and salinity using the thermodynamic equilibrium equations of [Liu & Millero (2002)](https://www.sciencedirect.com/science/article/abs/pii/S0304420301000743). 

$$
\begin{align}
T_K =& \quad \max(5.0, T) + 273.15 \\
\left(T_K\right)^{-1} =& \quad \dfrac{1}{T_K} \\
I_{S} =& \quad \dfrac{19.924,S}{1000 - 1.005,S}
\end{align}
$$

Solubility constants:

$$
\begin{align}
Fe_{sol1} =& \quad 10^{\left(-13.486 - 0.1856\sqrt{I_S} + 0.3073 I_S + 5254,\left(T_K\right)^{-1}\right)} \\
Fe_{sol2} =& \quad 10^{\left(2.517 - 0.8885\sqrt{I_S} + 0.2139 I_S - 1320,\left(T_K\right)^{-1}\right)} \\
Fe_{sol3} =& \quad 10^{\left(0.4511 - 0.3305\sqrt{I_S} - 1996,\left(T_K\right)^{-1}\right)} \\
Fe_{sol4} =& \quad 10^{\left(-0.2965 - 0.7881\sqrt{I_S} - 4086,\left(T_K\right)^{-1}\right)} \\
Fe_{sol5} =& \quad 10^{\left(4.4466 - 0.8505\sqrt{I_S} - 7980,\left(T_K\right)^{-1}\right)}
\end{align}
$$

Final Fe(III) solubility:

$$
\begin{align}
dFe_{sol} =& \quad Fe_{sol1}\left([H^+]^3 + Fe_{sol2}[H^+]^2 + Fe_{sol3}[H^+] + Fe_{sol4} + \dfrac{Fe_{sol5}}{[H^+]}\right)\times10^{9}
\end{align}
$$

_where_ <br>
- $T_{K}$ is in situ water temperature (`ztemk`, [ºK]) <br>
- $I_{S}$ is a salinity coefficient (`zval`, [dimenionless]) <br>
- $[H^+]$ is in situ hydrogen ion concentration (`hp`, [mol L<sup>-1</sup>]) <br>
- $dFe_{sol}$ is the final estimated solubility of dissolved iron in seawater (`fe3sol`, [nmol Fe kg<sup>-1</sup>]) <br>

Next we **estimate the concentration of colloidal iron** in solution following [Tagliabue et al. 2023](https://www.nature.com/articles/s41586-023-06210-5) in the case that `do_colloidal_shunt == .true.`. If `do_colloidal_shunt == .false.` we consider no dissolved Fe to be in colloidal form. Colloidal dissolved Fe (`fecol(i,j,k)`, $dFe_{col}$, [mmol Fe m<sup>-3</sup>]) is whatever exceeds the inorganic solubility ceiling (`fe3sol`, $dFe_{sol}$, [mmol Fe m<sup>-3</sup>]), but we enforce a hard minimum that colloids are at least 10% of total dissolved Fe (`biofer`, $dFe$, [mmol Fe m<sup>-3</sup>]).

$$
\begin{align}
dFe_{col} =& \quad \max\left(0.1 dFe, \ dFe - dFe_{sol} \right)
\end{align}
$$

Following solving for colloidal Fe, we **partition the remaining dissolved Fe into ligand-bound and free iron**. To do so, we find the remaining dissolved iron not in colloidal form (`fe_sfe`, $dFe_{sFe}$, [mmol Fe m<sup>-3</sup>]), 

$$
\begin{align}
dFe_{sFe} =& \quad \max\left(0.0,\ dFe - dFe_{col} \right)
\end{align}
$$

Partitioning of iron between free and ligand-bound forms is done using one of two approaches.

When `do_two_ligands == .false.`, we use a single ligand class and solve for the equilibrium fractionation between ligand-bound and free iron using a standard quadratic form. When `do_two_ligands == .true.`, we assume complexation of iron by a weak and a strong ligand and therefore solve for the equilibrium fractionation between free iron, weakly ligand-bound iron and strongly ligand-bound iron via an iterative root solver.

In either case, we first determine the conditional stability constant(s) of the ligand(s). In the case of `do_two_ligands == .true.`, we solve for the stability constant of a strong ligand (`ligK(i,j,k)`, $Lig_{s}^{K}$, [kg mol<sup>-1</sup>]) and then consider the stability constant of a weak ligand to be a constant offset equal to -1.5 log<sub>10</sub> units based on [Gledhill & Buck (2012)](https://doi.org/10.3389/fmicb.2012.00069). In the case of `do_two_ligands == .false.`, we solve for the stability constant of the strong (`ligK(i,j,k)`) and weak ligands (`ligW_K`), but take the concentration-weighted average binding strength to get the bulk ligand binding stregnth.

The stability constant (`ligK(i,j,k)`, $Lig_{s}^{K}$, [kg mol<sup>-1</sup>]) is known to vary with the environmental conditions. In WOMBAT-mid, we consider the effect of temperature, light, pH and the concentration of labile DOC on the binding strength. The temperature dependency comes from [Volker & Tagliabue (2015)](https://doi.org/10.1016/j.marchem.2014.11.008) and warmer waters increase binding strength. The light-dependency accounts for the photoreduction of photoreactive ligands, which was identified to reduce the conditional stability constant of aquachelin by 0.7 log<sub>10</sub> units ([Barbeau et al., 2001](https://doi.org/10.1038/35096545); [Vraspir & Butler, 2009](https://doi.org/10.1146/annurev.marine.010908.163712)). The pH and DOC concentration dependency comes from [Ye et al. (2020)](https://doi.org/10.1029/2019GB006425) and increases binding strength at lower pH and higher concentrations of DOC.

$$
\begin{align}
Lig_{s}^{K} =& \quad 10^{-9} \cdot \bigg( 10^{ \left(17.27 - 1565.7 \left(T_K\right)^{-1} \right)} \\
             & \qquad  10^{\left(-0.7 \dfrac{PAR}{PAR + 10}\right)} \\
             & \qquad 10^{\left(-0.0002 \left(B_{DOM}^{C}\right)^{2} + 0.034 \cdot B_{DOM}^{C} - 1.67 \cdot pH + 24.36\right)} \bigg)
\end{align}
$$

_where_ <br>
- $T_K$ is in situ water temperature (`ztemk`, [ºK]) <br>
- $PAR$ is the total photosynthetically available radiation (`radbio`, [W m<sup>-2</sup>]) <br>
- pH is the in situ pH <br>
- $B_{DOM}^{C}$ is the in situ concentration of dissolved organic carbon (`biodoc`, [mmol m<sup>-3</sup>]) <br>


After finding $Lig_{s}^{K}$ we solve for the free dissolved Fe concentration (`feIII`, $dFe_{free}$, [nmol Fe kg<sup>-1</sup>]) via the analytic method when `do_two_ligands == .false.`:

$$
\begin{align}
z =& \quad 1.0 + [Ligand] \cdot Lig_{bulk}^{K} - dFe_{sFe}\cdot Lig_{bulk}^{K} \\
dFe_{free} =& \quad \dfrac{-z + \sqrt{z^2 + 4.0 Lig_{bulk}^{K} dFe_{sFe}}}{2 Lig_{bulk}^{K} + \varepsilon} \\
dFe_{free} =& \quad \max\left(0,\ \min(dFe_{free}, dFe_{sFe})\right)
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

Whatever soluble dissolved iron is not present as inorganic free iron is assigned to ligand-bound dissolved iron:

$$
\begin{align}
dFe_{lig} =& \quad dFe_{sFe} - dFe_{free}
\end{align}
$$

Now that we have separated the dissolved Fe pool into its subcomponents of free, ligand-bound and colloidal Fe, we solve for scavenging of free iron and coagulation of colloidal, both of which remove dissolved iron and transfer these to two sinking authigenic particles. These authigenic sinking particles include a small, slowly sinking type (`f_afe(i,j,k)`, $Fe_{sA}$, [mol Fe kg<sup>-1</sup>]) and a large, fast sinking type (`f_bafe(i,j,k)`, $Fe_{lA}$, [mol Fe kg<sup>-1</sup>]). Their sinking rates are controlled by the input parameters `wafe` and `wbafe`. Both scavenging and colloidal coagulation are the major sinks of dissolved iron outside of phytoplankton uptake and this dissolved iron is transferred to the sinking authigenic pools.

**Scavenging:**

Scavenging of dissolved iron specifically affects free iron, is accelerated by the presence of particles in the water column and routes this iron to two sinking authigenic phases. Total scavenging of dissolved iron (`fescaven(i,j,k)`, $Sc_{dFe}^{\rightarrow}$, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]) is calculated as

$$
\begin{align}
Sc_{dFe}^{\rightarrow} =& \quad dFe_{free} \left(10^{-7} + \gamma_{dFe}^{scav} \cdot B_{particles}^{M} \right)
\end{align}
$$

_where_ <br>
- $dFe_{free}$ is the in situ concentration of dissolved free iron (`feIII(i,j,k)`, [nmol Fe kg<sup>-1</sup>]) <br>
- $\gamma_{dFe}^{scav}$ is the rate constant of scavenging (`kscav_dfe`, [(mmol m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $B_{particles}^{M}$ is the in situ concentration of detrital particles in the water column (`partic`, [mmol m<sup>-3</sup>]) <br>

$$
\begin{align}
B_{particles}^{M} =& \quad 2 \cdot B_{sd}^{C} + 2 \cdot B_{ld}^{C} + 2 \cdot B_{ld}^{Si} + 8.3 \cdot B_{CaCO_3}^{C}
\end{align}
$$

_where_ <br>
- $B_{sd}^{C}$ is the in situ concentration of small organic carbon detritus (`biodet`, [mmol C m<sup>-3</sup>]) <br>
- $B_{ld}^{C}$ is the in situ concentration of large organic carbon detritus (`biobdet`, [mmol C m<sup>-3</sup>]) <br>
- $B_{ld}^{Si}$ is the in situ concentration of biogenic silica detritus (`biobdetsi`, [mmol Si m<sup>-3</sup>]) <br>
- $B_{CaCO_3}^{C}$ is the in situ concentration of calcium carbonate detritus (`biocaco3`, [mmol C m<sup>-3</sup>]) <br>

Organic carbon-based particle types $B_{sd}^{C}$ and $B_{ld}^{C}$ are multipled by 2 assuming that carbon represents half the mass of the particle, $B_{ld}^{Si}$ is multipled by 2 assuming that it represents biogenic silica with a molecular mass of 60 g mol<sup>-1</sup>, and inorganic carbon-based particles $B_{CaCO_3}^{C}$ is multipled by 8.3 since the molecular weight of calcium carbonate is 100 g mol<sup>-1</sup>. 

Total scavenging ($Sc_{dFe}^{\rightarrow}$) of free iron is then broken into two parts: scavenging to small authigenic particles (`fescaafe(i,j,k)`, $Sc_{dFe}^{\rightarrow Fe_{sA}}$, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]) and scavenging to large authigenic particles (`fescabafe(i,j,k)`, $Sc_{dFe}^{\rightarrow Fe_{lA}}$, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]). 

$$
\begin{align}
Sc_{dFe}^{\rightarrow Fe_{sA}} =& \quad Sc_{dFe}^{\rightarrow} \cdot \dfrac{ 2 \cdot B_{sd}^{C} + 8.3 \cdot B_{CaCO_3}^{C} }{ B_{particles}^{M} } \\
Sc_{dFe}^{\rightarrow Fe_{lA}} =& \quad Sc_{dFe}^{\rightarrow} \cdot \dfrac{ 2 \cdot B_{ld}^{C} + 2 \cdot B_{ld}^{Si} }{ B_{particles}^{M} }
\end{align}
$$


**Coagulation:**

Similarly to scavenging of free iron, coagulation routes dissolved iron to two sinking authigenic phases. However, coagulation acts on the colloidal fraction of dissolved iron ([Tagliabue et al., 2023](https://www.nature.com/articles/s41586-023-06210-5)). Rates of coagulation of colloidal iron to small, slowly sinking authigenic iron (`fecoag2afe(i,j,k)`, $Co_{dFe}^{\rightarrow Fe_{sA}}$, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]) and large, fast sinking authigenic iron (`fecoag2bafe(i,j,k)`, $Co_{dFe}^{\rightarrow Fe_{lA}}$, [mol Fe kg<sup>-1</sup> s<sup>-1</sup>]) follow the form:

$$
\begin{align}
Co_{dFe}^{\rightarrow Fe_{sA}} =& \quad dFe_{col} \gamma_{dFe}^{coag} \cdot S_{coag}^{sA} \\
Co_{dFe}^{\rightarrow Fe_{lA}} =& \quad dFe_{col} \gamma_{dFe}^{coag} \cdot S_{coag}^{lA}
\end{align}
$$

_where_ <br>
- $dFe_{col}$ is the in situ concentration of dissolved colloidal iron (`fecol(i,j,k)`, [mol Fe kg<sup>-1</sup>]) <br>
- $\gamma_{dFe}^{coag}$ is the iron coagulation rate constant (`kcoag_dfe`, [(mmol m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $S_{coag}^{sA}$ and $S_{coag}^{lA}$ are scaling coefficients to decelerate or accelerate coagulation of small and large particles (`zval`, [mmol C m<sup>-3</sup>]) <br>

The coagulation scaling coefficients are themselves dependent on the concentrations of dissolved organic carbon, particulate organic carbon, phytoplankton biomass and the rate of mixing. For small particle coagulation:

$$
\begin{align}
S_{coag}^{sA} =& \quad H_{mix} \left(10.8 \cdot F_{coag} \left(B_{DOM}^{C} + 40\right) + 9.05 \cdot B_{sd}^{C}\right) \\
& \quad + 2.49 \cdot B_{sd}^{C} + 115.02 \cdot F_{coag} \left(B_{DOC}^{C} + 40\right) + 725.7 \cdot B_{sd}^{C} \\
& \quad + \gamma_{dFe}^{agg} \cdot \dfrac{\left(dFe_{col}\right)^{4}}{\left(dFe_{col}\right)^{4} + \left(K_{dFe}^{agg}\right)^{4}} \\ 
F_{coag} =& \quad \dfrac{B_{np}^{C} + B_{mp}^{C}}{B_{np}^{C} + B_{mp}^{C} + 0.03}
\end{align}
$$

_where_ <br>
- $H_{mix}$ is a Heaviside step function that is equalt to 1 in the mixed layer and 0.01 beneath the mixed layer (`shear`, [dimensionless]) <br>
- $F_{coag}$ is a phytoplankton concentration dependent coagulation factor (`biof`, [dimensionless])  <br>
- $B_{np}^{C}$ and $B_{mp}^{C}$ are the concentrations of nano- and micro-phytoplankton biomass (`biophy`; `biodia`, [mmol C m<sup>-3</sup>])  <br>
- $B_{DOM}^{C}$ is the concentration of dissolved organic matter in carbon (`biodoc`, [mmol C m<sup>-3</sup>]) <br>
- $B_{sd}^{C}$ is the concentration of small organic detrital particles (`biodet`, [mmol C m<sup>-3</sup>]) <br>
- $\gamma_{dFe}^{agg}$ is the colloidal iron aggregation rate constant (`kagg_col`, [s<sup>-1</sup>]) <br>
- $K_{dFe}^{agg}$ is the half-saturation coefficient for colloidal iron aggregation (`kagg_kcol`, [µmol m<sup>-3</sup>]) <br>

For large particle coagulation:

$$
\begin{align}
S_{coag}^{lA} =& \quad \left(2 \cdot H_{mix} + 1.37 \right) B_{ld}^{C} + 1.94 \cdot B_{ld}^{C}
\end{align}
$$

_where_ <br>
- $H_{mix}$ is a Heaviside step function that is equalt to 1 in the mixed layer and 0.01 beneath the mixed layer (`shear`, [dimensionless]) <br>
- $B_{ld}^{C}$ is the concentration of large organic detrital particles (`biobdet`, [mmol C m<sup>-3</sup>]) <br>

Together, these terms implement a biologically mediated coagulation pathway in which iron removal from the dissolved pool is tightly coupled to ecosystem state. The formulation reflects the central conclusion of [Tagliabue et al. (2023)](https://www.nature.com/articles/s41586-023-06210-5): that iron cycling is not governed solely by inorganic chemistry, but is strongly regulated by biological activity, organic matter dynamics, and particle ecology across the upper ocean. 


**Dissolution:**

Small, slow sinking authigenic (`f_afe(i,j,k)`, $Fe_{sA}$, [mol Fe kg<sup>-1</sup>]) and a large, fast sinking authigenic iron (`f_bafe(i,j,k)`, $Fe_{lA}$, [mol Fe kg<sup>-1</sup>]) are returned back to the dissolved iron phase through reductive processes or complexation with ligands ([Tagliabue et al., 2023](https://www.nature.com/articles/s41586-023-06210-5)). We represent this process simply via dissolution rate cofficients:

$$
\begin{align}
D_{sA}^{\rightarrow dFe} =& \quad Fe_{sA} \gamma_{sA}^{diss} \\
D_{lA}^{\rightarrow dFe} =& \quad Fe_{lA} \gamma_{lA}^{diss}
\end{align}
$$

_where_ <br>
- $\gamma_{sA}^{diss}$ is the constant dissolution rate of the small sinking authigenic iron (`kafe_dfe`, [s<sup>-1</sup>]) <br>
- $\gamma_{lA}^{diss}$ is the constant dissolution rate of the large sinking authigenic iron (`kbafe_dfe`, [s<sup>-1</sup>]) <br>

---


### 11. Biogenic silica dissolution. 

**Silicic acid equilibrium concentration**

To determine the rate of biogenic silica dissolution we must first determine the equilibrium concentration of silicic acid (H<sub>4</sub>SiO<sub>4</sub>) in seawater. To do so, we solve for this equilibrium concentration via thermodynamic first-principles:

$$
\begin{align}
K_{H_{4}SiO_{4}}(T,P) =& \quad \dfrac{\gamma_{H_{4}SiO_{4}^{0}} \cdot [H_{4}SiO_{4}]^{eq}}{(a_{H_{2}O})^{2}}
\end{align}
$$

_where_ <br>
- $K_{H_{4}SiO_{4}}(T,P)$ is the thermodynamic equilibrium constant in seawater at a given temperature and pressure (`K_am_silica`, [mol Si kg<sup>-1</sup>]) <br>
- $\gamma_{H_{4}SiO_{4}^{0}}$ is the activity ratio of H<sub>4</sub>SiO<sub>4</sub> in seawater (`gamma0`, [dimensionless]) <br>
- $[H_{4}SiO_{4}]^{eq}$ is the equilibrium concentration of H<sub>4</sub>SiO<sub>4</sub> (`sileqc(i,j,k)`, [mol Si kg<sup>-1</sup>]) <br>
- $a_{H_{2}O}$ is the activity of seawater (`alphaH2O`, [dimensionless]) <br>

The equation is rearranged such that:

$$
\begin{align}
[H_{4}SiO_{4}]^{eq} =& \quad \dfrac{K_{H_{4}SiO_{4}}(T,P) \cdot (a_{H_{2}O})^{2}}{\gamma_{H_{4}SiO_{4}^{0}}}
\end{align}
$$

The activity of seawater is slightly less than 1 due to dissolved salts lowering its chemical potential and so we set $a_{H_{2}O}$ equal to 0.999 ([IOC, SCOR & IAPSO, 2010](https://www.teos-10.org/pubs/TEOS-10_Manual.pdf)). For $\gamma_{H_{4}SiO_{4}^{0}}$ we follow [Savenko 2014](https://doi.org/10.1134/S0001437014020222) who demonstrated that the solubility of H<sub>4</sub>SiO<sub>4</sub> decreases predictably with salinity according to

$$
\begin{align}
\gamma_{H_{4}SiO_{4}^{0}} =& \quad 1 + 0.0053 \cdot S - 0.000034 \cdot S^{2}
\end{align}
$$

_where_ <br>
- $S$ is the in situ salinity of seawater (`Salt(i,j,k)`, [psu]) <br>

For $K_{H_{4}SiO_{4}}(T,P)$ we follow the derivation of [Gunnarsson & Arnórsson (2000)](https://doi.org/10.1016/S0016-7037(99)00426-3) who relate the thermodynamic equilibrium constant of H<sub>4</sub>SiO<sub>4</sub> to variations in temperature at a constant pressure of 1 bar ($P^{1}$):

$$
\begin{align}
K(T,P^{1}) =& \quad 10^{ \left( -8.476 - \frac{485.24}{T_{K}} - 2.268 \times 10^{-6} \cdot (T_{K})^{2} + 3.068 \cdot log_{10}(T_{K}) \right)}
\end{align}
$$

_where_ <br>
- $T_{K}$ is the in situ temperature of seawater ([`zval`, [ºK]]) <br>

We add a classic pressure correction to $K(T,P^{1})$ to retrieve $K(T,P)$ of the form:

$$
\begin{align}
K(T,P) =& \quad K(T,P^{1}) \cdot e^{- \dfrac{\Delta V^{0}}{RT_{K}}P}
\end{align}
$$

_where_ <br>
- $\Delta V^{0}$ is the partial molal volume change (`deltaV0`, [m<sup>3</sup> mol<sup>-1</sup>]) <br>
- $R$ is the universal gas constant (`Rgas`, [J ºK<sup>-1</sup> mol<sup>-1</sup>]) <br>
- $T_{K}$ is the in situ temperature of seawater ([`zval`, [ºK]]) <br>
- $P$ is the in situ pressure [`zm(i,j,k) * 1.0e4`, [Pa]] <br>

We obtain $\Delta V^{0}$ from [Willey, 1982](https://doi.org/10.1016/0016-7037(82)90015-1) and [Loucaides et al., 2012](https://doi.org/10.1016/j.marchem.2012.04.003) of roughly -9.0 [cm<sup>3</sup> mol<sup>-1</sup>], which we convert to [m<sup>3</sup> mol<sup>-1</sup>] by multiplying by $10^{-6}$. The negative value of $\Delta V^{0}$ implies an increase in dissolution of silica at higher pressures. These value return equilibrium concentrations of silicic acid on the order of 1000 to 1800 mmol m<sup>-3</sup>. Temperature increases are the largest control, while pressure increases from the surface to the ocean bottom increase solubility by 15-20%. 


**Biogenic silica dissolution**

Biogenic silica is only considered associated with the large type of sinking particulate organic matter and dissolution (`bsidiss(i,j,k)`, $D_{B_{ld}^{Si}}^{\rightarrow Si}$, [mol Si kg<sup>-1</sup> s<sup>-1</sup>]) occurs via

$$
\begin{align}
D_{B_{ld}^{Si}}^{\rightarrow Si} =& \quad d_{B_{ld}^{Si}} \cdot B_{ld}^{Si}
\end{align}
$$

_where_ <br>
- $d_{B_{ld}^{Si}}$ is the rate of biogenic silica dissolution (`disssi(i,j,k)`, [s<sup>-1</sup>]) <br>
- $B_{ld}^{Si}$ is the in situ concentration of biogenic silica (`f_bdetsi(i,j,k)`, [mol Si kg<sup>-1</sup>]) <br>

We treat the dissolution rate of biogenic silica ($d_{B_{ld}^{Si}}$) as dependent on three conditions: the degree of undersaturation ([Rickert, 2000](https://epic.awi.de/id/eprint/26530/1/BerPolarforsch2000351.pdf); [Van Cappellen & Qiu, 1997](https://doi.org/10.1016/S0967-0645(96)00112-9); [Van Cappellen et al., 2002](https://doi.org/10.1029/2001GB001431)), in situ temperature ([Kamatani, 1982](https://doi.org/10.1007/BF00393146); [Greenwood et al., 2005](https://doi.org/10.1007/s10498-004-9515-y)) and the activity of heterotrophic microbes ([Bidle & Azam, 1999](https://doi.org/10.1038/17351); [Bidle et al., 2003](https://doi.org/10.4319/lo.2003.48.5.1855)). To account for these conditions, we formulate the rate of dissolution of biogenic silica (`disssi(i,j,k)`, $d_{Si}$, [s<sup>-1</sup>]) as the product of three terms:

$$
\begin{align}
d_{B_{ld}^{Si}} =& \quad d_{B_{ld}^{Si}}^{T} \cdot S_{B_{ld}^{Si}}^{Sat} \cdot S_{B_{ld}^{Si}}^{bio}
\end{align}
$$

_where_ <br>
- $d_{B_{ld}^{Si}}^{T}$ is the temperature-dependent rate of dissolution (`disssi_temp`, [s<sup>-1</sup>]) <br>
- $S_{B_{ld}^{Si}}^{Sat}$ is a scaling factor that decelerates dissolution as the in situ concentration approachs the equilibrium concentration (`disssi_usat`, [dimenionless]) <br>
- $S_{B_{ld}^{Si}}^{bio}$ is a scaling factor that accelerates dissolution in the presence of heterotrophic bacterial biomass (`disssi_bact`, [dimenionless]) <br>

First, we solve for $d_{B_{ld}^{Si}}^{T}$. [Kamatani, 1982](https://doi.org/10.1007/BF00393146) measured dissolution rates of biogenic silica collected in Tokyo Bay between 8ºC and 30ºC and identified that these roughly obeyed the equation:

$$
\begin{align}
d_{B_{ld}^{Si}}^{T} =& \quad \dfrac{e^{\left(\alpha + β T\right)}}{3600}
\end{align}
$$

_where_ <br>
- $\alpha$ is a species-dependent dissolution intercept that ranges between -7.35 and -10.38. We set $\alpha$ = -8.0 <br>
- $β$ is the slope common to all species and is equal to 0.0833 <br>
- $T$ is the in situ temperature of seawater ([`Temp(i,j,k)`, [ºC]]) <br>
- $3600$ converts the rate from [hour<sup>-1</sup>] to [s<sup>-1</sup>] <br>

Next, we apply scaling terms that either decelerate or accelerate dissolution. Given that equilibrium concentrations of H<sub>4</sub>SiO<sub>4</sub> vary between 1000 to 1800 mmol m<sup>-3</sup> in the ocean, while actual in situ concentrations rarely exceed 200 mmol m<sup>-3</sup>, H<sub>4</sub>SiO<sub>4</sub> is always undersaturated. We therefore assume that H<sub>4</sub>SiO<sub>4</sub> is highly undersaturated everywhere in the ocean. According to [Van Cappellen et al., 2002](https://doi.org/10.1029/2001GB001431) "Detailed kinetic studies of biogenic silica dissolution conducted in flow-through reactors demonstrate that at very high degrees of undersaturation the dissolution kinetics switch from a linear dependence on the degree of undersaturation to an exponential one". Hence, we apply equation 2.13 from [Rickert, 2000](https://epic.awi.de/id/eprint/26530/1/BerPolarforsch2000351.pdf):

$$
\begin{align}
S_{B_{ld}^{Si}}^{Sat} =& \quad \left(1 - \dfrac{[H_{4}SiO_{4}]}{[H_{4}SiO_{4}]^{eq}}\right)^{2}
\end{align}
$$

_where_ <br>
- $[H_{4}SiO_{4}]$ is the in situ concentration of H<sub>4</sub>SiO<sub>4</sub> (`f_sil(i,j,k)`, [mol Si kg<sup>-1</sup>]) <br>
- $[H_{4}SiO_{4}]^{eq}$ is the equilibrium concentration of H<sub>4</sub>SiO<sub>4</sub> (`sileqc(i,j,k)`, [mol Si kg<sup>-1</sup>]) <br>
- we use an exponent of $2$ informed by the organic carbon-rich biogenic silica dissolution kinetics reported in Table 3.4 of [Rickert, 2000](https://epic.awi.de/id/eprint/26530/1/BerPolarforsch2000351.pdf) <br>

The scaling term associated with activity of heterotrophic bacteria is informed by substantial evidence. According to [Rickert et al., 2002](https://doi.org/10.1016/S0016-7037(01)00757-8) "The removal of organic or inorganic coatings enhance the reactivity by at least an order of magnitude". Order of magnitude increases in silica dissolution have been reported for diatom frustules in contact with bacteria ([Bidle & Azam, 1999](https://doi.org/10.1038/17351)), while anti-biotic treatments to mesocosms off Monterey Bay caused silica dissolution to reduced by 50% ([Bidle et al., 2003](https://doi.org/10.4319/lo.2003.48.5.1855)). We represent this bacterially-induced stimulation of dissolution with

$$
\begin{align}
S_{B_{ld}^{Si}}^{bio} =& \quad 1 + F_{B_{ld}^{Si}}^{bac} \cdot \dfrac{B_{b-p}^{C}}{B_{b-p}^{C} + K_{B_{ld}^{Si}}^{bac}}
\end{align}
$$

_where_ <br>
- $F_{B_{ld}^{Si}}^{bac}$ is the factor increase in dissolution caused by peak particle-associated bacterial biomass (`bsi_fbac`, [dimenionless]) <br>
- $K_{B_{ld}^{Si}}^{bac}$ is the half-saturation coefficient for stimulation of silica dissolution in the presence of particle-associated bacterial biomass (`bsi_kbac`, [mmol C m<sup>-3</sup>]) <br>
- $B_{b-p}^{C}$ is the in situ concentration of particle-associated bacterial biomass (`biobacp`, [mmol C m<sup>-3</sup>]) <br>

---


### 12. Mortality terms

Mortality of ecological functional types are affected by both linear ($\gamma$) and quadratic ($\Gamma$) terms. Linear terms are per-capita losses associated with the costs of basal metabolism. Quadratic, and thus density-dependent losses, are associated with disease, aggregation and coagulation, viruses, infection and cannibalism. None of these processes are represented explicitly within the model, so we represent them implicitly.

**Linear losses** of nano-phytoplankton (<sub>np</sub>), micro-phytoplankton (<sub>mp</sub>), micro-zooplankton (<sub>mz</sub>), meso-zooplankton (<sub>Mz</sub>), facultative heterotrophic bacterial types (<sub>b-p</sub>, <sub>b-f1</sub>, <sub>b-f2</sub>), and ammonia oxidizing archaea (<sub>aoa</sub>) in [mol kg<sup>-1</sup> s<sup>-1</sup>] are modelled as

$$
\begin{align}
\gamma_{np}^{\rightarrow C} =& \quad \gamma_{np}^{0ºC} (β_{hete})^{T} B_{np}^{C} \\
\gamma_{mp}^{\rightarrow C} =& \quad \gamma_{mp}^{0ºC} (β_{hete})^{T} B_{mp}^{C} \\
\gamma_{mz}^{\rightarrow C} =& \quad \gamma_{mz}^{0ºC} (β_{hete})^{T} B_{mz}^{C} \\
\gamma_{Mz}^{\rightarrow C} =& \quad \gamma_{Mz}^{0ºC} (β_{hete})^{T} B_{Mz}^{C} \\
\gamma_{b-p}^{\rightarrow C} =& \quad \gamma_{b}^{0ºC} (β_{hete})^{T} B_{b-p}^{C} \\
\gamma_{b-f1}^{\rightarrow C} =& \quad \gamma_{b}^{0ºC} (β_{hete})^{T} B_{b-f1}^{C} \\
\gamma_{b-f2}^{\rightarrow C} =& \quad \gamma_{b}^{0ºC} (β_{hete})^{T} B_{b-f2}^{C} \\
\gamma_{aoa}^{\rightarrow C} =& \quad \gamma_{aoa}^{0ºC} (β_{hete})^{T} B_{aoa}^{C}
\end{align}
$$

_where_ <br>
- $\gamma_{np}^{0ºC}$ is the rate of linear mortality of nano-phytoplankton at 0ºC (`phylmor`, [s<sup>-1</sup>]) <br>
- $\gamma_{mp}^{0ºC}$ is the rate of linear mortality of micro-phytoplankton at 0ºC (`dialmor`, [s<sup>-1</sup>]) <br>
- $\gamma_{mz}^{0ºC}$ is the rate of linear mortality of micro-zooplankton at 0ºC (`zoolmor`, [s<sup>-1</sup>]) <br>
- $\gamma_{Mz}^{0ºC}$ is the rate of linear mortality of meso-zooplankton at 0ºC (`meslmor`, [s<sup>-1</sup>]) <br>
- $\gamma_{b}^{0ºC}$ is the rate of linear mortality of all heterotrophic bacterial types at 0ºC (`baclmor`, [s<sup>-1</sup>]) <br>
- $\gamma_{aoa}^{0ºC}$ is the rate of linear mortality of ammonia oxidizing archaea at 0ºC (`aoalmor`, [s<sup>-1</sup>]) <br>
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless]) <br>
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC]) <br>
- $B_{np}^{C}$ is the concentration of nano-phytoplankton carbon biomass (`f_phy(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{mp}^{C}$ is the concentration of micro-phytoplankton carbon biomass (`f_dia(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{mz}^{C}$ is the concentration of micro-zooplankton carbon biomass (`f_zoo(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{Mz}^{C}$ is the concentration of meso-zooplankton carbon biomass (`f_mes(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{b-p}^{C}$ is the concentration of facultative NO<sub>3</sub>-reducing particle-associated bacteria carbon biomass (`f_bacp(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{b-f1}^{C}$ is the concentration of facultative NO<sub>3</sub>-reducing free-living bacteria carbon biomass (`f_bacf1(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{b-f2}^{C}$ is the concentration of facultative N<sub>2</sub>O-reducing free-living bacteria carbon biomass (`f_bacf2(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{aoa}^{C}$ is the concentration of ammonia oxidizing archaea carbon biomass (`f_aoa(i,j,k)`, [mol kg<sup>-1</sup>]) <br>


**Quadratic losses** of nano-phytoplankton (<sub>np</sub>), micro-phytoplankton (<sub>mp</sub>), micro-zooplankton (<sub>mz</sub>), meso-zooplankton (<sub>Mz</sub>), ffacultative heterotrophic bacterial types (<sub>b-p</sub>, <sub>b-f1</sub>, <sub>b-f2</sub>), and ammonia oxidizing archaea (<sub>aoa</sub>) in [mol kg<sup>-1</sup> s<sup>-1</sup>] are modelled as

$$
\begin{align}
\Gamma_{np}^{\rightarrow C} =& \quad \Gamma_{np}^{0ºC} (β_{hete})^{T} \left(B_{np}^{C}\right)^{2} \\
\Gamma_{mp}^{\rightarrow C} =& \quad \Gamma_{mp}^{0ºC} (β_{hete})^{T} \left(B_{mp}^{C}\right)^{2} \\
\Gamma_{mz}^{\rightarrow C} =& \quad \Gamma_{mz}^{0ºC} (β_{hete})^{T} \left(B_{mz}^{C}\right)^{2} \\
\Gamma_{Mz}^{\rightarrow C} =& \quad \Gamma_{Mz}^{0ºC} (β_{hete})^{T} \left(B_{Mz}^{C}\right)^{2} \\
\Gamma_{b-p}^{\rightarrow C} =& \quad \Gamma_{b}^{0ºC} (β_{hete})^{T} \left(B_{b-p}^{C}\right)^{2} \\
\Gamma_{b-f1}^{\rightarrow C} =& \quad \Gamma_{b}^{0ºC} (β_{hete})^{T} \left(B_{b-f1}^{C}\right)^{2} \\
\Gamma_{b-f2}^{\rightarrow C} =& \quad \Gamma_{b}^{0ºC} (β_{hete})^{T} \left(B_{b-f2}^{C}\right)^{2} \\
\Gamma_{aoa}^{\rightarrow C} =& \quad \Gamma_{aoa}^{0ºC} (β_{hete})^{T} \left(B_{aoa}^{C}\right)^{2}
\end{align}
$$

_where_ <br>
- $\Gamma_{np}^{0ºC}$ is the rate of quadratic mortality of nano-phytoplankton at 0ºC (`phyqmor`, [(mol C kg<sup>-1</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\Gamma_{mp}^{0ºC}$ is the rate of quadratic mortality of micro-phytoplankton at 0ºC (`diaqmor`, [(mol C kg<sup>-1</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\Gamma_{mz}^{0ºC}$ is the rate of quadratic mortality of micro-zooplankton at 0ºC (`zooqmor`, [(mol C kg<sup>-1</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\Gamma_{Mz}^{0ºC}$ is the rate of quadratic mortality of meso-zooplankton at 0ºC (`mesqmor`, [(mol C kg<sup>-1</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\Gamma_{b}^{0ºC}$ is the rate of quadratic mortality of all heterotrophic bacterial types at 0ºC (`bacqmor`, [(mol C kg<sup>-1</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\Gamma_{aoa}^{0ºC}$ is the rate of quadratic mortality of ammonia oxidizing archaea at 0ºC (`aoaqmor`, [(mol C kg<sup>-1</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless]) <br>
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC]) <br>
- $B_{np}^{C}$ is the concentration of nano-phytoplankton carbon biomass (`f_phy(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{mp}^{C}$ is the concentration of micro-phytoplankton carbon biomass (`f_dia(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{mz}^{C}$ is the concentration of micro-zooplankton carbon biomass (`f_zoo(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{Mz}^{C}$ is the concentration of meso-zooplankton carbon biomass (`f_mes(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{b-p}^{C}$ is the concentration of facultative NO<sub>3</sub>-reducing particle-associated bacteria carbon biomass (`f_bacp(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{b-f1}^{C}$ is the concentration of facultative NO<sub>3</sub>-reducing free-living bacteria carbon biomass (`f_bacf1(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{b-f2}^{C}$ is the concentration of facultative N<sub>2</sub>O-reducing free-living bacteria carbon biomass (`f_bacf2(i,j,k)`, [mol kg<sup>-1</sup>]) <br>
- $B_{aoa}^{C}$ is the concentration of ammonia oxidizing archaea carbon biomass (`f_aoa(i,j,k)`, [mol kg<sup>-1</sup>]) <br>

---


### 13. Zooplankton grazing, egestion, excretion and assimilation.

**Grazing** by micro-zooplankton (`g_zoo`, $g_{mz}$, [s<sup>-1</sup>]) and meso-zooplankton (`g_mes`, $g_{Mz}$, [s<sup>-1</sup>]) is computed using a Holling Type III functional response [Holling, 1959](https://doi.org/10.4039/Ent91385-7), where:

$$
\begin{align}
g_{mz} =& \quad \dfrac{\mu_{mz}^{max} (β_{hete})^{T} \sum_{i} \left(\varepsilon_{mz}^{i} \left(\phi_{mz}^{i} B_{i}^{C}\right)^{2}\right)}{\mu_{mz}^{max} (β_{hete})^{T} + \sum_{i} \left( \varepsilon_{mz}^{i} \left(\phi_{mz}^{i} B_{i}^{C}\right)^{2}\right)} \\
g_{Mz} =& \quad \dfrac{\mu_{Mz}^{max} (β_{hete})^{T} \sum_{i} \left(\varepsilon_{Mz}^{i} \left(\phi_{Mz}^{i} B_{i}^{C}\right)^{2}\right)}{\mu_{Mz}^{max} (β_{hete})^{T} + \sum_{i} \left( \varepsilon_{Mz}^{i} \left(\phi_{Mz}^{i} B_{i}^{C}\right)^{2}\right)}
\end{align}
$$

_where_ <br>
- $\mu_{mz}^{max}$ is the maximum rate of micro-zooplankton grazing at 0ºC (`zoogmax`, [s<sup>-1</sup>]) <br>
- $\mu_{Mz}^{max}$ is the maximum rate of meso-zooplankton grazing at 0ºC (`mesgmax`, [s<sup>-1</sup>]) <br>
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless]) <br>
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC]) <br>
- $B_{i}^{C}$ is the concentration of prey type $i$ carbon biomass ([mmol C m<sup>-3</sup>]) <br>
- $\phi_{mz}^{i}$ is the relative prey preference of micro-zooplankton for prey type $i$ ([dimenionless]) <br>
- $\phi_{Mz}^{i}$ is the relative prey preference of meso-zooplankton for prey type $i$ ([dimenionless]) <br>
- $\varepsilon_{mz}^{i}$ is the prey capture rate coefficient of micro-zooplankton for prey type $i$ ([(mmol C m<sup>-3</sup>)<sup>-2</sup>]) <br>
- $\varepsilon_{Mz}^{i}$ is the prey capture rate coefficient of meso-zooplankton for prey type $i$ ([(mmol C m<sup>-3</sup>)<sup>-2</sup>]) <br>

and where

$$
\begin{align}
\sum_{i} \phi_{mz}^{i} =& \quad \sum_{i} \phi_{Mz}^{i} = 1.0
\end{align}
$$

This formulation suppresses grazing at low prey biomass ($B_{i}^{C}$) due to reduced encounter and clearance rates, accelerates grazing at intermediate prey biomass as zooplankton effectively learn and switch to available prey, and saturates at high prey biomass due to handling-time limitation ([Gentleman and Neuheimer, 2008](https://doi.org/10.1093/plankt/fbn078); Rohr et al., [2022](https://doi.org/10.1016/j.pocean.2022.102878), [2024](https://doi.org/10.1029/2023GL107732)). This choice increases ecosystem stability and prolongs phytoplankton blooms relative to a Type II formulation.

The application of the temperature-dependent maximum growth rate in both the numerator and denominator makes this grazing formula unique [(Rohr et al., 2023)](https://www.nature.com/articles/s43247-023-00871-w) and equivalent to a disk formulation, rather than a Michaelis–Menten formulation [(Rohr et al., 2022)](https://doi.org/10.1016/j.pocean.2022.102878). Practically, this amplifies grazing in warmer climes, but to a lesser extent than other formulations that apply the temperature amplification ($(β_{hete})^{T}$) only in the numerator [(Rohr et al., 2023)](https://www.nature.com/articles/s43247-023-00871-w). This dampens the effect that variations in temperature have on grazing activity, amplifying the effect of $\varepsilon^{i}$ and aligning with observations that the ratio of grazing to phytoplankton growth varies little between tropical and polar climes [(Calbet and Landry, 2004)](https://doi.org/10.4319/lo.2004.49.1.0051). Theoretically, this assumes some evolutionary adaptation to account for the physiological effects of temperature across environmental niches, such that the efficiency of prey capture and handling becomes more important to grazers than metabolic constraints due to temperature.

The normalized prey preferences (i.e., dietary fractions) are further modified by prey switching prior to computation of total prey biomass ([Gentleman et al., 2003](https://doi.org/10.1016/j.dsr2.2003.07.001)) such that

$$
\begin{align}
\phi_{z}^{i} =& \quad \left( \phi_{z}^{i} B_{i}^{C} \right)^{s_{z}}
\end{align}
$$

_where_ <br>
- $\phi_{z}^{i}$ is the relative prey preference of zooplankton type $z$ for prey type $i$ <br>
- $B_{i}^{C}$ is the concentration of prey type $i$ in carbon biomass <br>
- $s_{z}$ is the prey-switching exponent of zooplankton type $z$ (`zoopreyswitch`; `mespreyswitch`) <br>

When $s_{z} < 1$, zooplankton feed equally across all prey items irrespective of availability  <br>
When $s_{z} = 1$, zooplankton feed according to pre-defined dietary fractions  <br>
When $s_{z} > 1$, zooplankton exhibit prey-switching and feed disproportionately on most abundant prey  <br>

Again, prey preferences are normalized to ensure 

$$
\begin{align}
\sum_{i} \phi_{mz}^{i} =& \quad \sum_{i} \phi_{Mz}^{i} = 1.0
\end{align}
$$

The community average prey capture rate coefficients of micro-zooplankton (`zooeps(i,j,k)`, $\varepsilon_{mz}$, [(mmol C m<sup>-3</sup>)<sup>-2</sup>]) and meso-zooplankton (`meseps(i,j,k)`, $\varepsilon_{Mz}$, [(mmol C m<sup>-3</sup>)<sup>-2</sup>]) vary as a function of the prey biomasses and the consequential variations in prey preferences associated with prey-switching, which is consistent with the prey-dependent behaviour described by [Rohr et al. (2024)](doi.org/10.1029/2023GL107732). 

Total grazing of biomass by micro-zooplankton ([mol C kg<sup>-1</sup> day<sup>-1</sup>]) is therefore

$$
\begin{align}
g_{mz}^{\leftarrow C} =& \quad g_{mz} B_{mz}^{C} \\
g_{Mz}^{\leftarrow C} =& \quad g_{Mz} B_{Mz}^{C}
\end{align}
$$

_where_ <br>
- $g_{mz}$ is the total specific rate of grazing of micro-zooplankton (`g_zoo`, [s<sup>-1</sup>]) <br>
- $g_{Mz}$ is the total specific rate of grazing of meso-zooplankton (`g_mes`, [s<sup>-1</sup>]) <br>
- $B_{mz}^{C}$ is the in situ concentration of micro-zooplankton carbon biomass (`f_zoo(i,j,k)`, [mol C kg<sup>-1</sup>]) <br>
- $B_{Mz}^{C}$ is the in situ concentration of meso-zooplankton carbon biomass (`f_mes(i,j,k)`, [mol C kg<sup>-1</sup>]) <br>

Total grazing of prey can also be expressed as the sum of individual prey type consumption:

$$
\begin{align}
g_{mz}^{\leftarrow C} =& \quad g_{mz}^{\leftarrow B_{np}^{C}} + g_{mz}^{\leftarrow B_{mp}^{C}} + g_{mz}^{\leftarrow B_{sd}^{C}} + g_{mz}^{\leftarrow B_{b-p}^{C}} + g_{mz}^{\leftarrow B_{b-f1}^{C}} + g_{mz}^{\leftarrow B_{b-f2}^{C}} + g_{mz}^{\leftarrow B_{aoa}^{C}} \\
g_{Mz}^{\leftarrow C} =& \quad g_{Mz}^{\leftarrow B_{np}^{C}} + g_{Mz}^{\leftarrow B_{mp}^{C}} + g_{Mz}^{\leftarrow B_{sd}^{C}} + g_{Mz}^{\leftarrow B_{ld}^{C}} + g_{Mz}^{\leftarrow B_{b-p}^{C}} + g_{Mz}^{\leftarrow B_{b-f1}^{C}} + g_{Mz}^{\leftarrow B_{b-f2}^{C}} + g_{Mz}^{\leftarrow B_{aoa}^{C}} + g_{Mz}^{\leftarrow B_{mz}^{C}}
\end{align}
$$

In this formulation, consumption of each prey item $i$ in [mol C kg<sup>-1</sup>] is equal to:

$$
\begin{align}
g_{z}^{\leftarrow B_{i}^{C}} =& \quad g_{z} B_{z}^{C} \cdot \dfrac{\varepsilon_{mz}^{i} \left(\phi_{mz}^{i} B_{i}^{C} \right)^{2}}{\sum_{i} \varepsilon_{mz}^{i} \left(\phi_{mz}^{i} B_{i}^{C} \right)^{2}}
\end{align}
$$

Thus: <br>
- $g_{mz}^{\leftarrow B_{np}^{C}}$ is the grazing rate of nano-phytoplankton by micro-zooplankton (`zoograzphy(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{mz}^{\leftarrow B_{mp}^{C}}$ is the grazing rate of micro-phytoplankton by micro-zooplankton (`zoograzdia(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{mz}^{\leftarrow B_{sd}^{C}}$ is the grazing rate of small particulate detritus by micro-zooplankton (`zoograzdet(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{mz}^{\leftarrow B_{b-p}^{C}}$ is the grazing rate of facultative NO<sub>3</sub>-reducing particle-associated bacteria by micro-zooplankton (`zoograzbacp(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{mz}^{\leftarrow B_{b-f1}^{C}}$ is the grazing rate of facultative NO<sub>3</sub>-reducing free-living bacteria by micro-zooplankton (`zoograzbacf1(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{mz}^{\leftarrow B_{b-f2}^{C}}$ is the grazing rate of facultative N<sub>2</sub>O-reducing free-living bacteria by micro-zooplankton (`zoograzbacf2(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{mz}^{\leftarrow B_{aoa}^{C}}$ is the grazing rate of ammonia oxidizing archaea by micro-zooplankton (`zoograzaoa(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{Mz}^{\leftarrow B_{np}^{C}}$ is the grazing rate of nano-phytoplankton by meso-zooplankton (`mesgrazphy(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{Mz}^{\leftarrow B_{mp}^{C}}$ is the grazing rate of micro-phytoplankton by meso-zooplankton (`mesgrazdia(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{Mz}^{\leftarrow B_{sd}^{C}}$ is the grazing rate of small particulate detritus by meso-zooplankton (`mesgrazdet(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{Mz}^{\leftarrow B_{ld}^{C}}$ is the grazing rate of large particulate detritus by meso-zooplankton (`mesgrazbdet(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{Mz}^{\leftarrow B_{b-p}^{C}}$ is the grazing rate of facultative NO<sub>3</sub>-reducing particle-associated bacteria by meso-zooplankton (`mesgrazbacp(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{Mz}^{\leftarrow B_{b-f1}^{C}}$ is the grazing rate of facultative NO<sub>3</sub>-reducing free-living bacteria by meso-zooplankton (`mesgrazbacf1(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{Mz}^{\leftarrow B_{b-f2}^{C}}$ is the grazing rate of facultative N<sub>2</sub>O-reducing free-living bacteria by meso-zooplankton (`mesgrazbacf2(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{Mz}^{\leftarrow B_{aoa}^{C}}$ is the grazing rate of ammonia oxidizing archaea by meso-zooplankton (`mesgrazaoa(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{Mz}^{\leftarrow B_{mz}^{C}}$ is the grazing rate of micro-zooplankton by meso-zooplankton (`mesgrazzoo(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>


**Zooplankton egestion, excretion and assimilation** are then calculated assuming static assimilation coefficients. Grazed biomass is routed to either egestion or ingestion via an ingestion coefficient ($\lambda^{C}$, [mol C (mol C)<sup>-1</sup>]), with the egested fraction being equal to $1.0 - \lambda^{C}$. The biomass that is ingested is then split between assimilation and excretion based on an assimilation coefficient ($\eta^{C}$, [mol C (mol C)<sup>-1</sup>]) with the excreted fraction being equal to $1.0 - \eta^{C}$. Egestion ($E$), excretion ($X$) and assimilation ($A$) of organic carbon due to grazing of prey type $i$ by zooplankton type $z$ are:

$$
\begin{align}
E_{z}^{\leftarrow B_{i}^{C}} =& \quad g_{z}^{\leftarrow B_{i}^{C}} \left(1 - \lambda_{z}^{C} \right) \\
X_{z}^{\leftarrow B_{i}^{C}} =& \quad g_{z}^{\leftarrow B_{i}^{C}} \lambda_{z}^{C} \left(1 - \eta_{z}^{C} \right) \\
A_{z}^{\leftarrow B_{i}^{C}} =& \quad g_{z}^{\leftarrow B_{i}^{C}} \lambda_{z}^{C} \eta_{z}^{C}
\end{align}
$$

_where_ <br>
- $E_{z}^{\leftarrow B_{i}^{C}}$ is the rate of egestion of carbon biomass by zooplankton type $z$ feeding on prey type $i$ ([mol C kg<sup>-1</sup>]) <br>
- $X_{z}^{\leftarrow B_{i}^{C}}$ is the rate of excretion of carbon biomass by zooplankton type $z$ feeding on prey type $i$ ([mol C kg<sup>-1</sup>]) <br>
- $A_{z}^{\leftarrow B_{i}^{C}}$ is the rate of assimilation of carbon biomass by zooplankton type $z$ feeding on prey type $i$ ([mol C kg<sup>-1</sup>]) <br>
- $g_{z}^{\leftarrow B_{i}^{C}}$ is the grazing rate of zooplankton type $z$ on prey type $i$ ([mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\lambda_{z}^{C}$ is the fraction of prey carbon biomass that is ingested by zooplankton type $z$ (`zooCingest`; `mesCingest`, [mol C (mol C)<sup>-1</sup>]) <br>
- $\eta_{z}^{C}$ is the fraction of ingested prey carbon biomass that is assimilated by zooplankton type $z$ (`zooCassim`; `mesCassim`, [mol C (mol C)<sup>-1</sup>]) <br>

Total egestion, excretion and assimilation or carbon are therefore:

$$
\begin{align}
E_{z}^{\leftarrow C} =& \quad g_{z}^{\leftarrow C} \left(1 - \lambda_{z}^{C} \right) \\
X_{z}^{\leftarrow C} =& \quad g_{z}^{\leftarrow C} \lambda_{z}^{C} \left(1 - \eta_{z}^{C} \right) \\
A_{z}^{\leftarrow C} =& \quad g_{z}^{\leftarrow C} \lambda_{z}^{C} \eta_{z}^{C}
\end{align}
$$

Because we track both carbon and iron through the ecosystem components, we assign unique ingestion and assimilation coefficients to carbon and iron. This separation of ingestion and assimilation coefficients for iron and carbon follows [Le Mézo & Galbraith (2021)](https://doi.org/10.1002/lno.11597). For iron, we apply unique ingestion ($\lambda^{Fe}$, [mol Fe (mol Fe)<sup>-1</sup>]) and assimilation coefficients ($\eta^{Fe}$, [mol Fe (mol Fe)<sup>-1</sup>]). [Le Mézo & Galbraith (2021)](https://doi.org/10.1002/lno.11597) show that if $\lambda^{Fe} << \lambda^{C}$ then egestion is enriched in Fe:C, and it follows that $\eta^{Fe} >> \eta^{C}$ so that zooplankton can absorb sufficient iron from their prey. Consequently:

$$
\begin{align}
E_{z}^{\leftarrow B_{i}^{Fe}} =& \quad g_{z}^{\leftarrow B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \left(1 - \lambda_{z}^{Fe} \right) \\
X_{z}^{\leftarrow B_{i}^{Fe}} =& \quad g_{z}^{\leftarrow B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \lambda_{z}^{Fe} \left(1 - \eta_{z}^{Fe} \right) \\
A_{z}^{\leftarrow B_{i}^{Fe}} =& \quad g_{z}^{\leftarrow B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \lambda_{z}^{Fe} \eta_{z}^{Fe}
\end{align}
$$

_where_ <br>
- $E_{z}^{\leftarrow B_{i}^{Fe}}$ is the rate of egestion of iron biomass by zooplankton type $z$ feeding on prey type $i$ ([mol Fe kg<sup>-1</sup>]) <br>
- $X_{z}^{\leftarrow B_{i}^{Fe}}$ is the rate of excretion of iron biomass by zooplankton type $z$ feeding on prey type $i$ ([mol Fe kg<sup>-1</sup>]) <br>
- $A_{z}^{\leftarrow B_{i}^{Fe}}$ is the rate of assimilation of iron biomass by zooplankton type $z$ feeding on prey type $i$ ([mol Fe kg<sup>-1</sup>]) <br>
- $g_{z}^{\leftarrow B_{i}^{C}}$ is the grazing rate of zooplankton type $z$ on prey type $i$ ([mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\dfrac{B_{i}^{Fe}}{B_{i}^{C}}$ is the Fe:C ratio of prey type $i$ ([mol Fe (mol C)<sup>-1</sup>]) <br>
- $\lambda_{z}^{Fe}$ is the fraction of prey iron biomass that is ingested by zooplankton type $z$ (`zooFeingest`; `mesFeingest`, [mol Fe (mol Fe)<sup>-1</sup>]) <br>
- $\eta_{z}^{Fe}$ is the fraction of ingested prey iron biomass that is assimilated by zooplankton type $z$ (`zooFeassim`; `mesFeassim`, [mol Fe (mol Fe)<sup>-1</sup>]) <br>

Total egestion, excretion and assimilation or iron are therefore:

$$
\begin{align}
E_{z}^{\leftarrow Fe} =& \quad \sum_{i} \left( g_{z}^{\leftarrow B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \right) \cdot \left(1 - \lambda_{z}^{Fe} \right) \\
X_{z}^{\leftarrow Fe} =& \quad \sum_{i} \left( g_{z}^{\leftarrow B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \right) \cdot \lambda_{z}^{Fe} \left(1 - \eta_{z}^{Fe} \right) \\
A_{z}^{\leftarrow Fe} =& \quad \sum_{i} \left( g_{z}^{\leftarrow B_{i}^{C}} \dfrac{B_{i}^{Fe}}{B_{i}^{C}} \right) \cdot \lambda_{z}^{Fe} \eta_{z}^{Fe}
\end{align}
$$

**Excretion of nitrogen**

For zooplankton preying on phytoplankton, other zooplankton and detritus, their excretion of nitrogen will be equal to:

$$
\begin{align}
X_{z}^{\leftarrow B_{i}^{N}} =& \quad X_{z}^{\leftarrow B_{i}^{C}} \dfrac{16}{122}
\end{align}
$$

However, since both micro-zooplankton and meso-zooplankton consume heterotrophic bacteria and ammonia oxidizing arcaheal types, which have different C:N ratios to other ecosystem biomass components, we must also compute the specific excretion of NH<sub>4</sub> and $B_{DOM}^{N}$ by zooplankton when feeding on these types. Since these types are richer in N than the other prey types, zooplankton excrete more NH<sub>4</sub> and $B_{DOM}^{N}$ when bacteria and archaea represent a greater proportion of their diet ([Sterner & Elser, 2002](https://press.princeton.edu/books/ebook/9781400885695/ecological-stoichiometry-pdf)). Total excretion of nitrogen from bacterial/archaeal type $i$ by zooplankton type $z$ is as follows:

$$
\begin{align}
X_{z}^{\leftarrow B_{i}^{N}} =& \quad g_{z}^{\leftarrow B_{i}^{C}} \dfrac{1}{R_{i}^{C:N}} - \dfrac{A_{z}^{\leftarrow B_{i}^{C}} + E_{z}^{\leftarrow B_{i}^{C}}}{R_{z}^{C:N}}
\end{align}
$$

---


### 14. Implicit nitrogen fixation.

Because we do not consider diazotrophs as an explicit phytoplankton functional type, we represent the fixation of nitrogen implicitly using a simple parameterization dependent on temperature, nutrient and light availability when `do_nitrogen_fixation == .true.`. The equation for new nitrogen (specifically NH<sub>4</sub>) added via diazotrophy is:

$$
\begin{align}
\mu_{diazo}^{\rightarrow NH_4} =& \quad \mu_{diazo}^{max} \left(1 - L_{np}^{N} \right) \\
                                & \min\left(L_{diazo}^{Fe}, L_{diazo}^{PAR}\right) R_{diazo}^{N:C} \cdot 1 \times 10^{-6}
\end{align}
$$

_where_ <br>
- $\mu_{diazo}^{max}$ is the temperature-dependent maximum growth rate of diazotrophs (`trimumax(i,j,k)`, [s<sup>-1</sup>]) <br>
- $L_{np}^{N}$ is the limitation term of nano-phytoplankton growth on nitrogen (`phy_lnit(i,j,k)`, [dimensionless])  <br>
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
- $dFe$ is the in situ concentration of dissolved iron (`biofer`, [nmol Fe kg<sup>-1</sup>]) <br>
- $K_{diazo}^{Fe}$ is the half-saturation coefficient for uptake of dissolved iron by diazotrophs (`trikf`, [nmol Fe kg<sup>-1</sup>]) <br>
- $\alpha_{diazo}$ is the chlorophyll-adjusted slope of the photosynthesis-irradience curve of diazotrophs (`alphabio_tri * trichlc`, [(W m<sup>-2</sup>)<sup>-1</sup>]) <br>
- $PAR$ is the downwelling photosynthetically available radiation (`radbio`, [W m<sup>-2</sup>]) <br>

---

### 15. Facultative bacterial heterotrophy.

We remineralise organic matter via the activity of three facultative bacterial heterotrophs: one particle-associated (`f_bacp(i,j,k)`, $B_{b-p}$, mol C kg<sup>-1</sup>) and two free-living (`f_bacf1(i,j,k)`, $B_{b-f1}$, mol C kg<sup>-1</sup>; `f_bacf2(i,j,k)`, $B_{b-f2}$, mol C kg<sup>-1</sup>). The particle-associated bacteria oxidise particulate organic matter ($B_{sd}^{C}$ and $B_{ld}^{C}$) while free-living types oxidise dissolved organic matter ($B_{DOM}^{C}$). All types reduce dissolved oxygen (O<sub>2</sub>), but they also perform different steps of denitrification when oxygen is limiting [(Sun et al., 2024)](https://doi.org/10.1073/pnas.2417421121) by being facultatively anaerobic ([Zakem et al., 2020](https://doi.org/10.1038/s41396-019-0523-8)).

**Growth requirements and products during full oxidation**

We compute the resource requirements and products of heterotrophic bacteria performing full oxidation following [Zakem et al. (2020)](https://doi.org/10.1038/s41396-019-0523-8) and [Rittman & McCarty (2001)](https://books.google.com.au/books/about/Environmental_Biotechnology_Principles_a.html?id=1PMeAQAAIAAJ&redir_esc=y). For a heterotrophic bacteria performing aerobic metabolism, that is consuming an organic matter substrate (S) and oxygen (O<sub>2</sub>) to produce biomass (B), CO<sub>2</sub> and inorganic nutrient, we split their metabolism into three half reactions and normalize by carbon:

$$
\begin{align}
CH_{h_S}O_{o_S}N_{n_S} + (2 - o_{S} + n_{S})H_{2}O \rightarrow& n_{S} NH_{4}^{+} + (1-n_{S})CO_{2} + n_{S}HCO_{3}^{-} + d_{S}H^{+} + d_{S}e^{-} \\
(1-f)d_{S} [ \dfrac{1}{4}O_{2} + H^{+} + e^{-} \rightarrow& \dfrac{1}{2} H_{2}O ] \\
\dfrac{f\cdot d_{S}}{d_{B}} [ n_{B} NH_{4}^{+} + (1-n_{B})CO_{2} + n_{B}HCO_{3}^{-} + d_{B} H^{+} + d_{B} e^{-} \rightarrow& CH_{h_{B}}O_{o_{B}}N_{n_{B}} + (2 - o_{B} + n_{B})H_{2}O 
\end{align}
$$

_where_ <br>
- $CH_{h_S}O_{o_S}N_{n_S}$ is the carbon, hydrogen, oxygen and nitrogen stoichiometry of the organic substrate (S), normalized by carbon  <br>
- $CH_{h_B}O_{o_B}N_{n_B}$ is the carbon, hydrogen, oxygen and nitrogen stoichiometry of the bacterial biomass (B), normalized by carbon  <br>
- $f$ is the fraction of electrons that are routed to biomass synthesis, rather than oxygen reduction via the respiratory electron chain (`bacp_fele`; `bacf1_fele`; `bacf2_fele`, [e (e)<sup>-1</sup>]) <br>
- $d_{S}$ is the number of electrons per carbon atom within the substrate (S) (`e_pom`; `e_dom`, [e]) <br>
- $d_{B}$ is the number of electrons per carbon atom within the bacterial biomass (B) (`e_bac`, [e]) <br>

The first half reaction is the electron donor oxidation, where organic matter substrate is fully oxidized to its inorganic constituents and electrons are released. The second and third constitute the electron acceptor reduction, which generates ATP, and biomass synthesis, which uses ATP. $f$ in this case represents the fraction of electrons that are routed from oxidation of organic matter substrate to biomass synthesis and can be thought of as the underlying growth efficiency of the cell. All reactions are scaled to one mole carbon in the substrate.

If we sum these three equations and treat S as $CH_{h_S}O_{o_S}N_{n_S}$ and B as $CH_{h_{B}}O_{o_{B}}N_{n_{B}}$:

$$
\begin{align}
S + \dfrac{(1 - f)d_{S}}{4} O_{2} \rightarrow& \dfrac{f\cdot d_{S}}{d_{B}} B + (1 - \dfrac{f \cdot d_{S}}{d_{B}})CO_{2} + (n_{S} - \dfrac{n_{B} \cdot f \cdot d_{S}}{d_{B}}) NH_{4}^{+} 
\end{align}
$$

We calculate $d_{S}$ and $d_{B}$ as equal to:

$$
\begin{align}
d_{S} =& 4 + h_{S} - 2 o_{S} - 3 n_{S} \\
d_{B} =& 4 + h_{B} - 2 o_{B} - 3 n_{B} 
\end{align}
$$

_where_ <br>
- $h_{S}$, $o_{S}$, and $n_{S}$ are the ratios of hydrogen, oxygen and nitrogen to carbon within the substrate <br>
- $h_{B}$, $o_{B}$, and $n_{B}$ are the ratios of hydrogen, oxygen and nitrogen to carbon within the bacterial biomass <br>

A high $d$ value indicates an organic molecule that is highly reduced with a high energy content (i.e., many electrons per carbon atom), whereas a low $d$ value means the molecule is more oxidized with less electrons and less potential energy [(Rittman & McCarty, 2001)](https://books.google.com.au/books/about/Environmental_Biotechnology_Principles_a.html?id=1PMeAQAAIAAJ&redir_esc=y). For the biomass of heterotrophic marine bacteria we assume a stoichiometry of $CH_{1.4}O_{0.4}N_{0.2}$ ([White et al., 2019](https://doi.org/10.1002/lol2.10103); [Zimmerman et al., 2014](https://doi.org/10.1111/1462-2920.12329)), which returns $d_{B}$ = 4. 

For the particle-associated bacteria, we assume labile marine organic matter with a typical stoichiometry of CH<sub>1.65</sub>O<sub>0.4</sub>N<sub>0.131</sub> ([Anderson, 1995](https://doi.org/10.1016/0967-0637(95)00072-E)), which gives a $d_{S}$ = 4.45. This returns an overall stoichiometry of heterotrophic bacteria performing complete and aerobic oxidation of organic matter of:

$$
\begin{align}
S + \dfrac{(1 - f)4.45}{4} O_{2} \rightarrow& \dfrac{f\cdot 4.45}{4} B + (1 - \dfrac{f \cdot 4.45}{4})CO_{2} + (0.131 - \dfrac{0.2 \cdot 4.45 f }{4}) NH_{4}^{+} 
\end{align}
$$

In this case, with full oxidation of organic matter, the biomass yield of the bacteria, $y_{B}$ is equal to $\min \left(1 - \alpha, \dfrac{f D}{d_{B}} \right)$, which is equivalent to $f \cdot \dfrac{4.45}{4}$ with our assumptions of substrate and bacterial stoichtry.

For the free-living bacteria that consume the dissolved organic matter substrate, the $d_{S}$ varies dynamically because we allow the hydrogen, oxygen and nitrogen content of dissolved organic matter to vary, meaning that $h_{S}$, $o_{S}$ and $n_{S}$ ratios change in space and time. If free-living bacteria consume a more reduced, labile organic substrate, this increases their growth yields.

**Growth requirements and products during partial oxidation**

Heterotrophic bacteria do not always completely oxidize the organic substrates they feed on to CO<sub>2</sub>. It is well appreciated that cross-feeding across different bacterial types occurs in nature, where one bacteria will excrete partially oxidized material and this will be used by another ([Amarnath et al., 2023](https://doi.org/10.1038/s41467-023-38913-8); [Braakman et al., 2025](https://www.science.org/doi/full/10.1126/sciadv.adp1949); [Pontrelli et al., 2022](https://www.science.org/doi/10.1126/sciadv.abk3076); [Reintjes et al., 2019](https://doi.org/10.1038/s41396-018-0326-3)).

We build from the above equations to include partial oxidation by considering an additional dissolved organic matter product (P) on the right-hand-side of the half reactions:

$$
\begin{align}
CH_{h_S}O_{o_S}N_{n_S} + x H_{2}O \rightarrow& \alpha CH_{h_P}O_{o_P}N_{n_P} + y CO_{2} + z NH_{4}^{+} + z HCO_{3}^{-} + (d_{S} - \alpha d_{P})H^{+} + (d_{S} - \alpha d_{P})e^{-} \\
(1-f)(d_{S} - \alpha d_{P}) [ \dfrac{1}{4}O_{2} + H^{+} + e^{-} \rightarrow& \dfrac{1}{2} H_{2}O ] \\
\dfrac{f\cdot(d_{S} - \alpha d_{P})}{d_{B}} [ n_{B} NH_{4}^{+} + (1-n_{B})CO_{2} + n_{B}HCO_{3}^{-} + d_{B} H^{+} + d_{B} e^{-} \rightarrow& CH_{h_{B}}O_{o_{B}}N_{n_{B}} + (2 - o_{B} + n_{B})H_{2}O \\
\end{align}
$$

_where_ <br>
- $\alpha$ is the fraction of organic substrate that undergoes partial oxidation (`bacp_alpha`; `bacf1_alpha`; `bacf2_alpha`, [dimensionless])
- $CH_{h_S}O_{o_S}N_{n_S}$ is the carbon, hydrogen, oxygen and nitrogen stoichiometry of the organic substrate (S), normalized by carbon  <br>
- $CH_{h_B}O_{o_B}N_{n_B}$ is the carbon, hydrogen, oxygen and nitrogen stoichiometry of the bacterial biomass (B), normalized by carbon  <br>
- $CH_{h_P}O_{o_P}N_{n_P}$ is the carbon, hydrogen, oxygen and nitrogen stoichiometry of the organic product (P), normalized by carbon  <br>
- $f$ is the fraction of electrons that are routed to biomass synthesis, rather than oxygen reduction via the respiratory electron chain (`bacp_fele`; `bacf1_fele`; `bacf2_fele`, [e (e)<sup>-1</sup>]) <br>
- $d_{S}$ is the number of electrons per carbon atom within the substrate (S) (`e_pom`; `e_dom`, [e]) <br>
- $d_{B}$ is the number of electrons per carbon atom within the bacterial biomass (B) (`e_bac`, [e]) <br>
- $d_{P}$ is the number of electrons per carbon atom within the product (P) (`e_domp`, [e]) <br>

The first half reaction is the electron donor oxidation, where organic matter is partially oxidized to its inorganic constituents, an organic matter product and some electrons are released. The second and third constitute the electron acceptor reduction, which generates ATP, and biomass synthesis, which uses ATP. Again, $f$ represents the fraction of electrons that are routed from oxidation of organic matter to biomass synthesis and is an overall cellular efficiency for growth.

$$
\begin{align}
x =& 2 - o_{S} + n_{S} + \alpha (o_{P} - 2 - n_{P}) \\
y =& (1 - \alpha - n_{S} + \alpha n_{P}) \\
z =& (n_{S} - \alpha n_{P})   
\end{align}
$$

If we sum the half reactions with partial oxidation (ignoring H<sub>2</sub>O), treating organic matter substrate S as $CH_{h_S}O_{o_S}N_{n_S}$, B as $CH_{h_{B}}O_{o_{B}}N_{n_{B}}$, and P as $CH_{h_{P}}O_{o_{P}}N_{n_{P}}$ we retrieve:

$$
\begin{align}
S + \dfrac{(1 - f)D}{4} O_{2} \rightarrow& \dfrac{f D}{d_{B}} B + \alpha P + (1 - \alpha - \dfrac{f D}{d_{B}})CO_{2} + (n_{S} - \alpha n_{P} - n_{B} \cdot \dfrac{f D}{d_{B}}) NH_{4}^{+} 
\end{align}
$$

_where_ <br>
- $D$ is the number of electrons released per carbon atom during oxidation (`e1_res`; `e2_res`; `e3_res`, [e]) <br>

and is equal to

$$
\begin{align}
D =& (d_{S} - \alpha d_{P})
\end{align}
$$

The biomass yield of the bacteria, $y_{B}$, performing partial oxidation is dependent on $D$ and is equal to

$$
\begin{align}
y_{B} =& \min \left(1 - \alpha, \dfrac{f D}{d_{B}} \right) \\
\end{align}
$$

Both $d_{S}$ and $d_{P}$ vary dynamically because we allow the hydrogen, oxygen and nitrogen content of dissolved organic matter to vary, meaning that $h_{S}$ and $h_{P}$, $o_{S}$ and $o_{P}$, and $n_{S}$ and $n_{P}$ ratios change in space and time. A higher $D$ increases the potential biomass yield. Increasing $\alpha$ generally lowers $D$, provided $d_{P} > 0$, because more substrate carbon is retained in the partially oxidized organic product rather than oxidized to release electrons.


**Organic products of partial oxidation**

We assume for simplicity that the dissolved organic matter coming from the breakdown of phytoplankton, zooplankton and particulate sinking detritus has a shared and constant stoichiometry equal to CH<sub>1.65</sub>O<sub>0.4</sub>N<sub>0.131</sub>, which reflects labile reduced organic matter ([Anderson et al., 1995](https://doi.org/10.1016/0967-0637(95)00072-E)). Phytoplankton overflow production of DOC (see Step 6), has a unique but constant stoichiometry of CH<sub>2.0</sub>O to represent exudation of carbohydrates ([Hansell & Carlson, 2014](https://books.google.com.au/books?id=7iKOAwAAQBAJ&lpg=PP1&ots=kzkdHuHMF_&dq=Carlson%20Hansell%202014%20doi&lr&pg=PP1#v=onepage&q&f=false)). 

These sources of dissolved organic matter are then reworked by heterotrophic bacteria, who may partially oxidize it. The partial oxidation of the organic matter substrate, whether particulate or dissolved, prefentially breaks C-H and C-N bonds while creating C-O bonds. We apply oxidation factors of 0.5 to C-H bonds (`Hox_fac`, $H_{ox}$), 1.5 to C-O bonds (`Oox_fac`, $O_{ox}$) and 0.6 to C-N bonds (`Nox_fac`, $N_{ox}$), which emulates the oxidation of glucose (CH<sub>2</sub>O) to glyoxylate-like (CHO<sub>1.5</sub>) and also reflects the preferential remineralisation of nitrogen ([]()). Therefore, for labile, newly produced organic matter of the form CH<sub>1.65</sub>O<sub>0.4</sub>N<sub>0.131</sub>, the partially oxidized product is of the form CH<sub>0.825</sub>O<sub>0.6</sub>N<sub>0.066</sub>. These factors can be altered at model run time, and altering them will affect the amount of energy that is released during partial oxidation via the $D$ value, as well as bacterial growth yield $y_{B}$, since

$$
\begin{align}
y_{B} =& \quad \min \left(1 - \alpha, \dfrac{f D}{d_{B}} \right) \\
D =& \quad (d_{S} - \alpha d_{P})
\end{align}
$$

and

$$
\begin{align}
d_{P} =& \quad 4 + h_{P} - 2 o_{P} - 3 n_{P}
\end{align}
$$

_where_ <br>
- $h_{P} = \quad h_{S} H_{ox}$ <br>
- $o_{P} = \quad o_{S} O_{ox}$ <br>
- $n_{P} = \quad n_{S} N_{ox}$ <br>

As a result, both particle-associated and free-living types produce more oxidized dissolved organic matter and preferentially release nitrogen to ammonium when performing partial oxidation.

**Competition and cross-feeding between free-living types**

The two free-living bacterial types compete for the same resource: dissolved organic matter. However, we control this competition by assigning each free-living bacterial type an optimal nominal oxidation state of carbon (`bacf1_nosc_opt`; `bacf2_nosc_opt`, $\rho_{b}^{NOSC}$) and a range around this optimum defined by a gaussian distribution (`bacf1_nosc_sig `; `bacf2_nosc_sig`, $\sigma_{b}^{NOSC}$) within which they feed. 

The fractional availability of dissolved organic carbon, $f_{b}^{[DOC]}$, that is "seen" by a bacterial type $b$ is calculated as

$$
\begin{align}
f_{b}^{[DOC]} =& \quad \dfrac{1}{\sigma_{b}^{NOSC} \sqrt{2 \pi}} e^{ -\dfrac{1}{2} \left(\dfrac{NOSC - \rho_{b}^{NOSC}}{\sigma_{b}^{NOSC}}\right)^{2}} 
\end{align}
$$

_where_ <br>
- $NOSC$ is the nominal oxidation state of dissolved organic carbon (`nosc`, [dimensionless]) <br>
- $\rho_{b}^{NOSC}$ is the target NOSC of bacterial type $b$ (`bacf1_nosc_opt`; `bacf2_nosc_opt`, [dimensionless]) <br>
- $\sigma_{b}^{NOSC}$ is the standard deviation around the optimal NOSC of bacterial type $b$ (`bacf1_nosc_sig`; `bacf2_nosc_sig`, [dimensionless]) <br>

Here, the nominal oxidation state of carbon ([La Rowe & Van Cappellen, 2011](https://www.sciencedirect.com/science/article/pii/S0016703711000378)) in organic matter substrate $S$ is equal to 

$$
\begin{align}
NOSC =& \quad \max \left(-4, \min \left(4, -\left(h_{S} - 2 o_{S} - 3 n_{S} \right) \right) \right)
\end{align}
$$

_where_ <br>
- $h_{S}$ is the H:C ratio of the organic matter substrate (`dom_H2C`, [mol H (mol C)<sup>-1</sup>]) <br>
- $o_{S}$ is the O:C ratio of the organic matter substrate (`dom_O2C`, [mol O (mol C)<sup>-1</sup>]) <br>
- $n_{S}$ is the N:C ratio of the organic matter substrate (`dom_N2C`, [mol N (mol C)<sup>-1</sup>]) <br>

The relative availability of DOC to each free-living bacterial type $b$ is then

$$
\begin{align}
[DOC]_{b} =& \quad \dfrac{f_{b}^{[DOC]}}{\sum_{b=1}^{2} f_{b}^{[DOC]} } [DOC] \\
\end{align}
$$

_where_ <br>
- $[DOC]$ is the ambient concentration of dissolved organic carbon (`biodoc`, [mol C kg<sup>-1</sup>])

Assigning each bacterial type the same target NOSC and distribution would result in direct and maximal competition. However, competition can be reduced by setting different NOSC optimums with non-intersecting, or at least partially intersecting, distributions. Further, a lack of resource niche differentiation can be turned into cross-feeding, and thus commensualism, if non-overlapping resource niches are combined with partial oxidation that converts reduced DOM with a low NOSC to more oxidized DOM with a higher NOSC. The targeting of different forms of DOM by different bacterial types is apparent in nature and may be a major contributor to diversity in bacterial functional types ([Reynolds et al., 2026](https://www.science.org/doi/full/10.1126/sciadv.adz0537)).

**Uptake of inorganic nutrients**

All heterotrophic bacteria assimilate dissolved iron ($dFe$) to support biosynthesis and will assimilate ammonium (NH<sub>4</sub>) if limited by nitrogen (very low N:C ratios of the organic substrate). Ammonium production is computed from stoichiometric balance, and it will become negative if biomass N demand exceeds the N supplied by the substrate. By taking up NH<sub>4</sub> and $dFe$, bacteria may compete directly with phytoplankton, consistent with prior observations ([Kirchman, 1994](https://www.jstor.org/stable/4251383); [Tortell et al., 1996](https://doi.org/10.1038/383330a0); [Kirchman & Wheeler, 1998](https://doi.org/10.1016/S0967-0637(97)00075-7); [Fourquez et al., 2015](https://doi.org/10.5194/bg-12-1893-2015); [Deng et al., 2021](https://doi.org/10.1002/lno.11883); [Strzepek et al., 2025](https://doi.org/10.1093/ismejo/wraf015)).

**Anaerobic growth**

We consider these bacterial types to be facultatively anaerobic to reflect the presence of denitrifying genes in ubiquitous SAR11 bacteria ([Zumft, 1997](https://doi.org/10.1128/mmbr.61.4.533-616.1997); [Tsementzi et al., 2016](https://doi.org/10.1038/nature19068)). This means that they can shift their metabolism to using either nitrate (NO<sub>3</sub>) or nitrous oxide (N<sub>2</sub>O) as an alternative electron acceptor when O<sub>2</sub> is limiting. Following [Sun et al. (2024)](https://doi.org/10.1073/pnas.2417421121), we consider the particle-associated type to perform complete denitrification from NO<sub>3</sub> to N<sub>2</sub> due to the abundance of DOC in particle pore water, while we consider the two free-living types to perform NO<sub>3</sub> reduction to N<sub>2</sub>O and N<sub>2</sub>O reduction to N<sub>2</sub>, respectively.

**Overall growth of bacterial types**

The realized biomass growth rate (integration of carbon into biomass) of bacterial functional type $b$ (`bacpgrow(i,j,k)`; `bacf1grow(i,j,k)`; `bacf2grow(i,j,k)`, $\mu_{b}^{\leftarrow C}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) is defined by:

$$
\begin{align}
\mu_{b}^{\leftarrow C} =& \quad \max\left(\mu_{b}^{aer}, \mu_{b}^{ana} \right) B_{b}^{C}
\end{align}
$$

_where_ <br>
- $\mu_{b}^{aer}$ is the realized growth rate due to aerobic metabolism (`bac_muaer`, [s<sup>-1</sup>]) <br>
- $\mu_{b}^{ana}$ is the realized growth rate due to anaerobic metabolism (`bac_muana`, [s<sup>-1</sup>]) <br>
- $B_{b}^{C}$ is the in situ concentration of bacterial functional type $b$ (`f_bacp(i,j,k)`; `f_bacf1(i,j,k)`; `f_bacf2(i,j,k)`, [mol C kg<sup>-1</sup>]) <br>

Thus, when `do_wc_denitrification == .true.` bacteria use whichever of aerobic and anaerobic metabolism offers the greatest growth rate. 

Both aerobic ($\mu_{b}^{aer}$) and anaerobic ($\mu_{b}^{ana}$) growth rates are calculated as the minimum of three resource-specific rates: growth supported by organic carbon substrate ($OC$), growth supported by dissolved iron ($dFe$), and growth supported by the electron acceptor ($EA$). For particle-associated bacteria, the organic carbon substrate is particulate organic carbon:

$$
\begin{align}
OC_{b-p} =& \quad B_{sd}^{C} + B_{ld}^{C}
\end{align}
$$

For free-living bacteria, the organic carbon substrate is dissolved organic carbon, $B_{DOM}^{C}$, partitioned between the two free-living types according to their NOSC-dependent availability functions.

For aerobic growth:

$$
\begin{align}
\mu_{b}^{aer} =& \quad \min \left(\mu_{b}^{aer(OC)}, \mu_{b}^{aer(dFe)}, \mu_{b}^{aer(EA)} \right) (β_{hete})^{T}
\end{align}
$$

For anaerobic growth:

$$
\begin{align}
\mu_{b}^{ana} =& \quad \min \left(\mu_{b}^{ana(OC)}, \mu_{b}^{ana(dFe)}, \mu_{b}^{ana(EA)} \right) (β_{hete})^{T}
\end{align}
$$

_where_ <br>
- $(β_{hete})^{T}$ is the temperature-dependent scaling on heterotrophic metabolism (`fbc`, [dimensionless]) <br>
- $\mu_{b}^{aer(OC)}$ and $\mu_{b}^{ana(OC)}$ are the potential specific growth rates supported by organic carbon substrate uptake ([s<sup>-1</sup>]) <br>
- $\mu_{b}^{aer(dFe)}$ and $\mu_{b}^{ana(dFe)}$ are the potential specific growth rates supported by dissolved iron uptake ([s<sup>-1</sup>]) <br>
- $\mu_{b}^{aer(EA)}$ and $\mu_{b}^{ana(EA)}$ are the potential specific growth rates supported by electron-acceptor uptake ([s<sup>-1</sup>]) <br>

When `do_wc_denitrification == .true.`, $\mu_{b}^{ana}$ ≥ 0.0. However, when `do_wc_denitrification == .false.`, $\mu_{b}^{ana}$ = 0.0. 

For aerobic growth, the resource-specific growth rates are:

$$
\begin{align}
\mu_{b}^{aer(OC)} =& \quad V_{b}^{OC} y_{b}^{aer(OC)} \\
\mu_{b}^{aer(dFe)} =& \quad V_{b}^{dFe} R_{b}^{C:Fe} \\
\mu_{b}^{aer(EA)} =& \quad \dfrac{V_{b}^{O_{2}}}{c_{b}^{O_{2}}}
\end{align}
$$

_where_ <br>
- $V_{b}^{OC}$ is the potential uptake rate of organic carbon substrate by bacterial type $b$ (`bac_Voc`, [mol C substrate (mol C biomass)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $V_{b}^{dFe}$ is the potential uptake rate of dissolved iron by bacterial type $b$ (`bac_VdFe`, [mol Fe (mol C biomass)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $V_{b}^{O_2}$ is the potential uptake rate of dissolved oxygen by bacterial type $b$ (`bac_Voxy`, [mol O<sub>2</sub> (mol C biomass)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $y_{b}^{aer(OC)}$ is the aerobic biomass yield on organic carbon substrate (`bacp_ypoc(i,j,k)`; `bacf1_ydoc(i,j,k)`; `bacf2_ydoc(i,j,k)`, [mol C biomass (mol C substrate)<sup>-1</sup>]) <br>
- $R_{b}^{C:Fe}$ is the bacterial carbon-to-iron ratio (`bac_C2Fe`, [mol C (mol Fe)<sup>-1</sup>]) <br>
- $c_{b}^{O_{2}}$ is the oxygen requirement per unit bacterial biomass produced (`bacp_coxy`; `bacf1_coxy`; `bacf2_coxy`, [mol O<sub>2</sub> (mol C biomass)<sup>-1</sup>]) <br>

For anaerobic growth, the organic carbon and iron terms are calculated in the same way, but the electron acceptor is nitrate or nitrous oxide rather than oxygen. Particle-associated bacteria reduce NO<sub>3</sub> to N<sub>2</sub>, the first free-living type reduces NO<sub>3</sub> to N<sub>2</sub>O, and the second free-living type reduces N<sub>2</sub>O to N<sub>2</sub>:

$$
\begin{align}
\mu_{b}^{ana(OC)} =& \quad V_{b}^{OC} y_{b}^{ana(OC)} \\
\mu_{b}^{ana(dFe)} =& \quad V_{b}^{dFe} R_{b}^{C:Fe} \\
\mu_{b-p}^{ana(EA)} =& \quad \dfrac{V_{b-p}^{NO_{3}}}{c_{b-p}^{NO_{3} \rightarrow N_{2}}} \\
\mu_{b-f1}^{ana(EA)} =& \quad \dfrac{V_{b-f1}^{NO_{3}}}{c_{b-f1}^{NO_{3} \rightarrow N_{2}O}} \\
\mu_{b-f2}^{ana(EA)} =& \quad \dfrac{V_{b-f2}^{N_{2}O}}{c_{b-f2}^{N_{2}O \rightarrow N_{2}}}
\end{align}
$$

_where_ <br>
- $V_{b-p}^{NO_3}$ is the potential nitrate uptake rate by particle-associated bacteria (`bac_Vno3`, [mol NO<sub>3</sub> (mol C biomass)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $V_{b-f1}^{NO_3}$ is the potential nitrate uptake rate by the first free-living bacterial type (`bac_Vno3`, [mol NO<sub>3</sub> (mol C biomass)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $V_{b-f2}^{N_{2}O}$ is the potential nitrous oxide uptake rate by the second free-living bacterial type (`bac_Vn2o`, [mol N<sub>2</sub>O (mol C biomass)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $y_{b}^{ana(OC)}$ is the anaerobic biomass yield on organic carbon substrate (`bacp_ypoc_ana`; `bacf1_ydoc_ana`; `bacf2_ydoc_ana`, [mol C biomass (mol C substrate)<sup>-1</sup>]) <br>
- $c_{b-p}^{NO_{3} \rightarrow N_{2}}$ is the nitrate requirement for complete denitrification by particle-associated bacteria (`bacp_cno3_ana`, [mol NO<sub>3</sub> (mol C biomass)<sup>-1</sup>]) <br>
- $c_{b-f1}^{NO_{3} \rightarrow N_{2}O}$ is the nitrate requirement for N<sub>2</sub>O production by the first free-living bacterial type (`bacf1_cno3_ana`, [mol NO<sub>3</sub> (mol C biomass)<sup>-1</sup>]) <br>
- $c_{b-f2}^{N_{2}O \rightarrow N_{2}}$ is the nitrous oxide requirement for N<sub>2</sub>O reduction by the second free-living bacterial type (`bacf2_cn2o_ana`, [mol N<sub>2</sub>O (mol C biomass)<sup>-1</sup>]) <br>

Whether aerobic or anaerobic metabolism results in higher bacterial growth therefore depends on both substrate uptake rates and the stoichiometric cost of converting those substrates into biomass. The model does not mix aerobic and anaerobic metabolism fractionally. Instead, it uses a binary anaerobic pathway selector:

$$
\begin{align}
I_{b}^{ana} =
\begin{cases}
1, & \mu_{b}^{ana} > \mu_{b}^{aer} \\
0, & \mu_{b}^{ana} \leq \mu_{b}^{aer}
\end{cases}
\end{align}
$$

This selector is then used to apply either the aerobic or anaerobic source-sink stoichiometry.

**Uptake rates** 

Potential uptake rates of organic carbon substrate, dissolved iron and electron acceptors are calculated as:

$$
\begin{align}
V_{b-p}^{POC} =& \quad V_{b-p}^{max,POC} \cdot \dfrac{B_{sd}^{C}+B_{ld}^{C}}{B_{sd}^{C}+B_{ld}^{C} + K_{b-p}^{POC}} \\
V_{b-f1}^{DOC} =& \quad V_{b-f1}^{max,DOC} \cdot \dfrac{B_{DOM,b-f1}^{C}}{B_{DOM,b-f1}^{C} + K_{b-f1}^{DOC}} \\
V_{b-f2}^{DOC} =& \quad V_{b-f2}^{max,DOC} \cdot \dfrac{B_{DOM,b-f2}^{C}}{B_{DOM,b-f2}^{C} + K_{b-f2}^{DOC}} \\
V_{b}^{dFe} =& \quad V_{b}^{max,dFe} \cdot \dfrac{dFe}{dFe + K_{b}^{dFe}} \\
V_{b}^{O_{2}} =& \quad \rho_{b}^{O_2} \cdot O_{2} \\
V_{b}^{NO_{3}} =& \quad V_{b}^{max,NO_{3}} \cdot \dfrac{NO_{3}}{NO_{3} + K_{b}^{NO_{3}}} \\
V_{b-f2}^{N_{2}O} =& \quad \rho_{b-f2}^{N_{2}O} \cdot N_{2}O
\end{align}
$$

_where_ <br>
- $V_{b-p}^{max,POC}$ is the maximum uptake rate of particulate organic carbon by particle-associated bacteria (`bacp_Vmax_poc`, [mol C substrate (mol C biomass)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $V_{b-f1}^{max,DOC}$ and $V_{b-f2}^{max,DOC}$ are the maximum uptake rates of DOC by the two free-living bacterial types (`bacf1_Vmax_doc`; `bacf2_Vmax_doc`, [mol C substrate (mol C biomass)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $V_{b}^{max,dFe}$ is the maximum uptake rate of $dFe$ by bacterial functional type $b$ (`bacp_Vmax_dfe`; `bacf1_Vmax_dfe`; `bacf2_Vmax_dfe`, [mol Fe (mol C biomass)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $V_{b}^{max,NO_{3}}$ is the maximum uptake rate of NO<sub>3</sub> by nitrate-reducing bacteria (`bacp_Vmax_no3`; `bacf1_Vmax_no3`, [mol NO<sub>3</sub> (mol C biomass)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\rho_{b}^{O_{2}}$ is the diffusive uptake coefficient for O<sub>2</sub> by bacterial functional type $b$ (`bacp_poxy`; `bacf1_poxy`; `bacf2_poxy`, [(mmol O<sub>2</sub> m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\rho_{b-f2}^{N_{2}O}$ is the diffusive uptake coefficient for N<sub>2</sub>O by the second free-living bacterial type (`bacf2_pn2o`, [(mmol N<sub>2</sub>O m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $K_{b-p}^{POC}$ is the half-saturation coefficient for POC uptake by particle-associated bacteria (`bacp_kpoc`, [mmol C m<sup>-3</sup>]) <br>
- $K_{b-f1}^{DOC}$ and $K_{b-f2}^{DOC}$ are the half-saturation coefficients for DOC uptake by free-living bacteria (`bacf1_kdoc`; `bacf2_kdoc`, [mmol C m<sup>-3</sup>]) <br>
- $K_{b}^{dFe}$ is the half-saturation coefficient for $dFe$ uptake by bacterial functional type $b$ (`bacp_kfer`; `bacf1_kfer`; `bacf2_kfer`, [µmol Fe m<sup>-3</sup>]) <br>
- $K_{b}^{NO_{3}}$ is the half-saturation coefficient for NO<sub>3</sub> uptake by nitrate-reducing bacteria (`bacp_kno3`; `bacf1_kno3`, [mmol N m<sup>-3</sup>]) <br>
- $B_{sd}^{C}$ and $B_{ld}^{C}$ are small and large detrital particulate organic carbon (`biodet`; `biobdet`, [mmol C m<sup>-3</sup>]) <br>
- $B_{DOM,b-f1}^{C}$ and $B_{DOM,b-f2}^{C}$ are the NOSC-partitioned DOC concentrations available to the two free-living bacterial types (`m2_doc`; `m3_doc`, [mmol C m<sup>-3</sup>]) <br>
- $dFe$ is the in situ concentration of dissolved iron (`biofer`, [µmol Fe m<sup>-3</sup>]) <br>
- $NO_3$ is the in situ concentration of nitrate (`biono3`, [mmol N m<sup>-3</sup>]) <br>
- $O_2$ is the in situ concentration of oxygen (`biooxy`, [mmol O<sub>2</sub> m<sup>-3</sup>]) <br>
- $N_2O$ is the in situ concentration of nitrous oxide (`bion2o`, [mmol N<sub>2</sub>O m<sup>-3</sup>]) <br>

NH<sub>4</sub> is not treated as an independent uptake-rate limitation in this bacterial growth calculation. Instead, net NH<sub>4</sub> production or uptake is diagnosed from the nitrogen balance of substrate consumption, partially oxidized DOC production, and bacterial biomass synthesis.

---


### 16. Calcium carbonate production and dissolution.

**Dynamic $CaCO_3$ production and dissolution**

When $CaCO_3$ dynamics are enabled (`do_caco3_dynamics = .true.`), the model computes both particulate inorganic carbon production (via the PIC:POC ratio) and $CaCO_3$ dissolution rates as functions of carbonate chemistry, temperature, and organic matter availability.

**Production** of $CaCO_3$ in WOMBAT-mid comes from five sources: (1) density-dependent mortality of nano-phytoplankton (i.e., coccolithophorids), (2) density-dependent mortality of micro-zooplankton (i.e., foraminifera), (3) micro-zooplankton egestion of grazed nano-phytoplankton, (4) meso-zooplankton egestion of grazed nano-phytoplankton, and (5) meso-zooplankton egestion of grazed micro-zooplankton. Each term is multiplied by the particulate inorganic to organic carbon production ratio (`pic2poc`, $PIC:POC$, [mol/mol]) to return a rate of $CaCO_3$ production in mol C kg<sup>-1</sup> s<sup>-1</sup>.

$$
\begin{align}
(1) & P_{CaCO_3}^{\Gamma_{np}^{C}} = \quad \Gamma_{np}^{\rightarrow C} \cdot PIC:POC \\
(2) & P_{CaCO_3}^{\Gamma_{mz}^{C}} = \quad \Gamma_{mz}^{\rightarrow C} \cdot PIC:POC \\
(3) & P_{CaCO_3}^{g_{mz}^{\leftarrow B_{np}^{C}}} = \quad g_{mz}^{\leftarrow B_{np}^{C}} \cdot PIC:POC \left(1 - F_{gut}\right) \\
(4) & P_{CaCO_3}^{g_{Mz}^{\leftarrow B_{np}^{C}}} = \quad g_{Mz}^{\leftarrow B_{np}^{C}} \cdot PIC:POC \left(1 - F_{gut}\right) \\
(5) & P_{CaCO_3}^{g_{Mz}^{\leftarrow B_{mz}^{C}}} = \quad g_{Mz}^{\leftarrow B_{mz}^{C}} \cdot PIC:POC \left(1 - F_{gut}\right)
\end{align}
$$

_where_ <br>
- $\Gamma_{np}^{\rightarrow C}$ is the quadratic (density-dependent) loss rate of nano-phytoplankton biomass (`phymorq`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\Gamma_{mz}^{\rightarrow C}$ is the quadratic (density-dependent) loss rate of micro-zooplankton biomass (`zoomorq`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{mz}^{\leftarrow B_{np}^{C}}$ is the grazing rate of nano-phytoplankton by micro-zooplankton (`zoograzphy(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{Mz}^{\leftarrow B_{np}^{C}}$ is the grazing rate of nano-phytoplankton by meso-zooplankton (`mesgrazphy(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{Mz}^{\leftarrow B_{mz}^{C}}$ is the grazing rate of micro-zooplankton by meso-zooplankton (`mesgrazzoo(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $F_{gut}$ is the fraction of $CaCO_3$ that is dissolved within zooplankton guts (`fgutdiss`, [mol C (mol C)<sup>-1</sup>]) <br>

In the above, the $PIC:POC$ ratio is formulated as 

$$
\begin{align}
PIC:POC =& \quad \min \left( 0.3,  \left( f_{inorg} + 10^{-3 + 4.31 \times 10^{-6} \left( \dfrac{[HCO_3^-]}{[H^+]} \right)} \right) F_T \right)
\end{align}
$$

_where_ <br>
- $f_{inorg}$ is the background PIC:POC ratio (f_inorg, [mol C (mol C)-1]) <br>
- $[HCO_{3}^{-}]$ is the concentration of bicarbonate ions (hco3, [mol kg-1]) <br>
- $[H^{+}]$ is concentration of free hydrogen ions (htotal(i,j,k), [µmol kg-1]) <br>
- $F_{T}$ is a temperature-dependent suppression term and if defined by $F_{T} = 0.55 + 0.45 \cdot \tanh\left(T - 4 \right)$  <br>

This formulation of $PIC:POC$ is therefore a function of the substrate–inhibitor ratio between bicarbonate and free hydrogen ions (`hco3 / htotal(i,j,k)`, $\dfrac{[HCO_{3}^{-}]}{[H^{+}]}$, [mol µmol<sup>-1</sup>]), following [Lehmann & Bach (2025)](https://www.nature.com/articles/s41561-025-01644-0). This reflects the sensitivity of calcification to carbonate system speciation, which is nonlinearly enhanced with an increasing substrate-inhibitor ratio. Moreover, the $F_{T}$ term strongly reduces production in cold waters, enforcing near-zero calcification below approximately 3°C consistent with observations of _Emiliania huxleyi_ growth limits in polar environments ([Fielding, 2013](https://doi.org/10.4319/lo.2013.58.2.0663)). Finally, we also cap the $PIC:POC$ ratio at an upper bound of 0.3 to prevent unrealistically high inorganic carbon production and accord with the highest measured ratios in the ocean.

**Dissolution** of $CaCO_3$ is computed as the sum of five contributions: 

(1) undersaturation-driven dissolution of calcite (`caldiss(i,j,k)`, $D_{CaCO_3}^{\Omega_{cal}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
(2) undersaturation-driven dissolution of aragonite (`aradiss(i,j,k)`, $D_{CaCO_3}^{\Omega_{ara}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
(3) biologically-mediated dissolution associated with degredation of small detrital organic matter (`pocdiss(i,j,k)`, $D_{CaCO_3}^{\Gamma_{sd}^{\rightarrow C}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
(4) dissolution within micro-zooplankton during their digestion of detrital aggregates (`zoodiss(i,j,k)`, $D_{CaCO_3}^{g_{mz}^{\leftarrow B_{sd}^{C}}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
(5) dissolution within meso-zooplankton during their digestion of detrital aggregates (`mesdiss(i,j,k)`, $D_{CaCO_3}^{g_{Mz}^{\leftarrow B_{sd}^{C}}}$, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>

Total $CaCO_3$ dissolution is:

$$
\begin{align}
D_{CaCO_3} =& \quad D_{CaCO_3}^{\Omega_{cal}} + D_{CaCO_3}^{\Omega_{ara}} + D_{CaCO_3}^{\Gamma_{sd}^{\rightarrow C}} + D_{CaCO_3}^{g_{mz}^{\leftarrow B_{sd}^{C}}} + D_{CaCO_3}^{g_{Mz}^{\leftarrow B_{sd}^{C}}}
\end{align}
$$

The first three terms follow [Kwon et al. (2024)](https://www.science.org/doi/full/10.1126/sciadv.adl0779):

$$
\begin{align}
(1) & D_{CaCO_3}^{\Omega_{cal}} = \quad d_{CaCO_3}^{\Omega_{cal}} \max\left(0,  1 - \Omega_{cal}\right)^{2.2} B_{CaCO_3}^{C} \\ 
(2) & D_{CaCO_3}^{\Omega_{ara}} = \quad d_{CaCO_3}^{\Omega_{ara}} \max\left(0,  1 - \Omega_{ara}\right)^{1.5} B_{CaCO_3}^{C} \\
(3) & D_{CaCO_3}^{\Gamma_{sd}^{\rightarrow C}} = \quad d_{CaCO_3}^{\Gamma_{sd}} \Gamma_{sd}^{\rightarrow C} B_{CaCO_3}^{C}
\end{align}
$$

_where_ <br>
- $\Omega_{cal}$ is the saturation state of calcite (`omega_cal(i,j,k)`, [dimenionless]) <br>
- $\Omega_{ara}$ is the saturation state of aragonite (`omega_ara(i,j,k)`, [dimenionless]) <br>
- $d_{CaCO_3}^{\Omega_{cal}}$ is the reference dissolution rate constant for calcite (`disscal`, [s<sup>-1</sup>])  <br>
- $d_{CaCO_3}^{\Omega_{ara}}$ is the reference dissolution rate constant for aragonite (`dissara`, [s<sup>-1</sup>])  <br>
- $d_{CaCO_3}^{\Gamma_{sd}}$ is the reference dissolution rate constant per unit of small detrital organic carbon remineralised (`dissdet`, [(mmol C m<sup>-3</sup>)<sup>-1</sup>])  <br>
- $\Gamma_{sd}^{\rightarrow C}$ is the in situ remineralisation rate of small detrital organic carbon (`detremi(i,j,k)`, [mmol C m<sup>-3</sup> s<sup>-1</sup>]) <br>
- $B_{CaCO_3}^{C}$ is the in situ concentration of $CaCO_3$ in carbon units (`f_caco3(i,j,k)`, [mol C kg<sup>-1</sup>]) <br>

For $D_{CaCO_3}^{\Omega_{cal}}$ and $D_{CaCO_3}^{\Omega_{ara}}$, dissolution is activated only under undersaturated conditions ($\Omega_{cal} < 1$; $\Omega_{ara} < 1$) and increases nonlinearly with increasing undersaturation. In contrast, $D_{CaCO_3}^{\Gamma_{sd}^{\rightarrow C}}$ represents shallow water dissolution due to reducing microenvironments. In this scenario, $\Omega_{cal}$ and $\Omega_{ara}$ tend to be > 1 ([Sulpis et al., 2021](https://doi.org/10.1038/s41561-021-00743-y)) but dissolution nonetheless occurs in microenvironments enriched in $CO_{2}^{*}$ due to heterotrophic activity ([Borer et al., 2026](https://doi.org/10.1073/pnas.2510025123)).

The fourth and fifth terms (`zoodiss(i,j,k)`; `mesdiss(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) represent dissolution of $CaCO_3$ during zooplankton digestion of detrital particulates. 

$$
\begin{align}
(4) & D_{CaCO_3}^{g_{mz}^{\leftarrow B_{sd}^{\leftarrow C}}} = \quad g_{mz}^{\leftarrow B_{sd}^{C}} F_{gut} \dfrac{B_{CaCO_3}^{C}}{B_{sd}^{C}} \\
(5) & D_{CaCO_3}^{g_{Mz}^{\leftarrow B_{sd}^{\leftarrow C}}} = \quad g_{Mz}^{\leftarrow B_{sd}^{C}} F_{gut} \dfrac{B_{CaCO_3}^{C}}{B_{sd}^{C}}
\end{align}
$$

_where_ <br>
- $g_{mz}^{\leftarrow B_{sd}^{C}}$ is the grazing rate of small particulate detritus by micro-zooplankton (`zoograzdet(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $g_{Mz}^{\leftarrow B_{sd}^{C}}$ is the grazing rate of small particulate detritus by meso-zooplankton (`mesgrazdet(i,j,k)`, [mol C kg<sup>-1</sup> s<sup>-1</sup>]) <br>
- $F_{gut}$ is the fraction of $CaCO_3$ that is dissolved within zooplankton guts (`fgutdiss`, [mol C (mol C)<sup>-1</sup>]) <br>
- $\dfrac{B_{CaCO_3}^{C}}{B_{sd}^{C}}$ is the in situ ratio of $CaCO_3$ to small organic carbon detritus (`biocaco3/biodet`, [mol C (mol C)<sup>-1</sup>]) <br>

Here we note that the processing of $CaCO_3$ by zooplankton grazing is treated differently to processing of organic carbon. For organic carbon, we route the biomass between zooplankton biomass (assimilation), inorganic nutrients (excretion) and particulate detritus (egestion). For $CaCO_3$ consumption by both micro-zooplankton and meso-zooplankton the $CaCO_3$ is not assimilated since it does not contain nitrogen or other key elements for biosynthesis, and so is only routed between excretion to DIC and alkalinity or goes undissolved and remains $CaCO_3$ that sinks through the water column. This is supported by the fact that micro- and meso-zooplankton may dissolve 92±7% and 38-73% of coccolithophore calcite during feeding, respectively ([Smith et al., 2024](https://doi.org/10.1126/sciadv.adr5453); [White et al., 2018](https://doi.org/10.1038/s41598-018-28073-x); [Harris, 1994](https://doi.org/10.1007/BF00347540)), and that the remainder is excreted and not assimilated ([Mayers et al., 2020](https://doi.org/10.3389/fmars.2020.569896)).

**Static $CaCO_{3}$ production and dissolution**

When $CaCO_3$ dynamics are disabled (`do_caco3_dynamics = .false.`), the model uses a static PIC:POC ratio (`f_inorg + 0.025`, [mol C (mol C)<sup>-1</sup>]) and $CaCO_3$ dissolution rate (`caco3lrem`, [s<sup>-1</sup>]). These are set as input parameters to the model.

---


### 17. Chemoautotrophy.

We consider two forms of chemoautotrophy carried out by two distinct forms of microbes: ammonia oxidizing archaea and anaerobic ammonia oxidizing (anammox) bacteria. Ammonia oxidizing archaea are considered explicitly within WOMBAT-mid (`f_aoa(i,j,k)`, [mol C kg<sup>-1</sup>]), while anammox bacteria are considered implicitly (when `do_anammox == .true.`) and therefore do not have varying biomasses (i.e., we only compute rates of anammox).

**Ammonia oxidizing archaea**

Growth of ammonia oxidizing archaea (`aoagrow(i,j,k)`, $\mu_{aoa}^{C}$, [mol C kg<sup>-1</sup>]) is defined similarly to other microbes:

$$
\begin{align}
\mu_{aoa}^{C} =& \quad \mu_{aoa} B_{aoa}^{C}
\end{align}
$$

_where_ <br>
- $\mu_{aoa}$ is the realized growth rate of ammonia oxidizing archaea (`aoa_mu(i,j,k)`, [s<sup>-1</sup>]) <br>
- $B_{aoa}^{C}$ is the in situ concentration of carbon biomass of ammonia oxidizing archaea (`f_aoa(i,j,k)`, [mol C kg<sup>-1</sup>]) <br>

The realized growth rate, $\mu_{aoa}$, is the minimum growth achievable on oxygen and ammonium:

$$
\begin{align}
\mu_{aoa} =& \quad \min\left(\mu_{aoa}^{NH_4}, \mu_{aoa}^{O_2}\right) \\
\mu_{aoa}^{NH_4} =& \quad \mu_{aoa}^{max} \dfrac{NH_4}{NH_4 + K_{aoa}^{NH_4}} \\
\mu_{aoa}^{O_2} =& \quad \dfrac{\rho_{aoa}^{O_2} O_2}{y_{aoa}^{O_2}}
\end{align}
$$

_where_ <br>
- $\mu_{aoa}^{max}$ is a temperature-dependent maximum growh rate of ammonia oxidizing archaea (`aoa_mumax(i,j,k)`, [s<sup>-1</sup>]) <br>
- $K_{aoa}^{NH_4}$ is the half-saturation coefficient for uptake of NH<sub>4</sub> by ammonia oxidizing archaea (`aoa_knh4`, [mmol N m<sup>-3</sup>]) <br>
- $\rho_{aoa}^{O_2}$ is the diffusive uptake limit of O<sub>2</sub> by ammonia oxidizing archaea (`aoa_poxy`, [(mmol C m<sup>-3</sup>)<sup>-1</sup> s<sup>-1</sup>]) <br>
- $y_{aoa}^{O_2}$ is the aerobic growth yield of ammonia oxidizing archaea on O<sub>2</sub> (`aoa_yoxy`, [mol C biomass (mol O<sub>2</sub>)<sup>-1</sup>]) <br>
- NH<sub>4</sub> is the in situ concentration of NH<sub>4</sub> (`bionh4`, [mmol N m<sup>-3</sup>]) <br>
- O<sub>2</sub> is the in situ concentration of O<sub>2</sub>$ (`biooxy`, [mmol O<sub>2</sub> m<sup>-3</sup>]) <br>

The temperature-dependent maximum growth rate of ammonia oxidizing archaea is informed by the cultures of [Qin et al. (2015)](https://doi.org/10.1073/pnas.1501568112):

$$
\begin{align}
\mu_{aoa}^{max} =& \quad \dfrac{\max\left(0.2, 0.029 \cdot T - 0.147 \right)}{86400}
\end{align}
$$

_where_ <br>
- $T$ is the in situ temperature of seawater (`Temp(i,j,k)`, [ºC]) <br>

In reality, ammonia oxidizing archaea perform the first step of the nitrification process by oxidizing ammonia through to nitrite. However, in WOMBAT-mid we do not consider nitrite oxidizing bacteria that then complete the second step of the nitrification process to produce nitrate. Hence, in this version of WOMBAT-mid we consider ammonia oxidizing archaea to perform full nitrification and oxidize NH<sub>4</sub> direclty to NO<sub>3</sub>. Consumption of NH<sub>4</sub> (`ammox(i,j,k)`, [mol N kg<sup>-1</sup> s<sup>-1</sup>]) and O<sub>2</sub> (`aoaresp(i,j,k)`, [mol O<sub>2</sub> kg<sup>-1</sup> s<sup>-1</sup>]) are calculated as:

$$
\begin{align}
\mu_{aoa}^{\leftarrow NH_4} =& \quad \dfrac{\mu_{aoa}^{C}}{y_{aoa}^{NH_4}} \\
\mu_{aoa}^{\leftarrow O_2} =& \quad \dfrac{\mu_{aoa}^{C}}{y_{aoa}^{O_2}}
\end{align}
$$

_where_ <br>
- $y_{aoa}^{NH_4}$ is the aerobic growth yield of ammonia oxidizing archaea on NH<sub>4</sub> (`aoa_ynh4`, [mol C biomass (mol NH<sub>4</sub>)<sup>-1</sup>]) <br>
- $y_{aoa}^{O_2}$ is the aerobic growth yield of ammonia oxidizing archaea on O<sub>2</sub> (`aoa_yoxy`, [mol C biomass (mol O<sub>2</sub>)<sup>-1</sup>]) <br>

Ammonia ozidizing archaea produce both NO<sub>3</sub> (in our model, while in reality they produce nitrite (NO<sub>2</sub>)) and a small amount of N<sub>2</sub>O as they oxidize NH<sub>4</sub>. The production of N<sub>2</sub>O is informed by numerous studies that show a positive relationship between the production of N<sub>2</sub>O and declining ambient dissolved oxygen concentrations ([Goreau et al., 1980](https://doi.org/10.1128/aem.40.3.526-532.1980); [Santoro et al., 2011](https://www.science.org/doi/full/10.1126/science.1208239); [Qin et al., 2017](https://doi.org/10.1111/1758-2229.12525); [Ji et al., 2018](https://doi.org/10.1029/2018GB005887); [Frey et al., 2023](https://doi.org/10.1002/lno.12283); [Kelly et al., 2024](https://doi.org/10.5194/bg-21-3215-2024)). We parameterize N<sub>2</sub>O production by ammonia oxidation via the relationships reported in [Kelly et al. (2024)](https://doi.org/10.5194/bg-21-3215-2024). These authors determined the fraction of N<sub>2</sub>O produced per unit of NO<sub>2</sub> via two pathways:

$$
\begin{align}
f_{aoa ( NH_4 + NH_4 )}^{\rightarrow N_{2}O} =& \quad 0.022 \cdot e^{\left( -1.50 \cdot O_2 \right)} + f_{min}^{\rightarrow N_{2}O} \\
f_{aoa ( NH_4 + NO_2 )}^{\rightarrow N_{2}O} =& \quad 0.204 \cdot e^{\left( -0.58 \cdot O_2 \right)}
\end{align}
$$

_where_ <br>
- $f_{aoa ( NH_4 + NH_4 )}^{\rightarrow N_{2}O}$ is the fraction of NH<sub>4</sub> diverted to N<sub>2</sub>O via the combination of two NH<sub>4</sub> (`aoa_en2o_nh4`, [mol N (mol N)<sup>-1</sup>]) <br>
- $f_{aoa ( NH_4 + NO_2 )}^{\rightarrow N_{2}O}$ is the fraction of NH<sub>4</sub> diverted to N<sub>2</sub>O via the combination of one NH<sub>4</sub> and one NO<sub>2</sub> (`aoa_en2o_hyb`, [mol N (mol N)<sup>-1</sup>]) <br>
- $O_2$ is the in situ concentration of dissolved oxygen (`biooxy`, [mmol $O_2$ m<sup>-3</sup>]) <br>
- $f_{min}^{\rightarrow N_{2}O}$ is the minimum fraction of NH<sub>4</sub> diverted to N<sub>2</sub>O in fully oxic conditions (`aoa_en2omin`, [mol N (mol N)<sup>-1</sup>]) <br>

Both pathways of N<sub>2</sub>O production show an exponential decline as dissolved oxygen increases. We then calculate the production of N<sub>2</sub>O per unit of biomass growth of ammonia oxidizing archaea (`aoa_en2o(i,j,k)`, $e_{aoa}^{\rightarrow N_{2}O}$, [mol N<sub>2</sub>O (mol C biomass)<sup>-1</sup>]) with:

$$
\begin{align}
e_{aoa}^{N_{2}O} =& \quad \left( \dfrac{1}{y_{aoa}^{NH_4} - R_{aoa}^{C:N}} \right) \cdot \left(0.5 \cdot f_{aoa ( NH_4 + NH_4 )}^{\rightarrow N_{2}O} + f_{aoa ( NH_4 + NO_2 )}^{\rightarrow N_{2}O} \right)
\end{align}
$$

_where_ <br>
- $y_{aoa}^{NH_4}$ is the aerobic growth yield of ammonia oxidizing archaea on NH<sub>4</sub> (`aoa_ynh4`, [mol C biomass (mol NH<sub>4</sub>)<sup>-1</sup>]) <br>
- $R_{aoa}^{C:N}$ is the ratio of carbon to nitrogen within the biomass of ammonia oxidizing archaea (`aoa_C2N`, [mol C (mol N)<sup>-1</sup>]) <br>
- and $0.5$ is applied against $f_{aoa ( NH_4 + NH_4 )}^{\rightarrow N_{2}O}$ because two molecules of NH<sub>4</sub> are combined to one molecule of N<sub>2</sub> <br>

We then finally compute the amount of NO<sub>3</sub> that is excreted by the ammonia oxidizing archaea as the residual of the remaining nitrogen atoms.

$$
\begin{align}
e_{aoa}^{NO_{3}} =& \quad \left( \dfrac{1}{y_{aoa}^{NH_4} - R_{aoa}^{C:N}} \right) \cdot \left(1 - f_{aoa ( NH_4 + NH_4 )}^{\rightarrow N_{2}O} - 2\cdot f_{aoa ( NH_4 + NO_2 )}^{\rightarrow N_{2}O} \right)
\end{align}
$$

_where_ <br>
- the $2$ is applied against $f_{aoa ( NH_4 + NO_2 )}^{\rightarrow N_{2}O}$ because one molecule of NO<sub>2</sub> is consumed during hybrid N<sub>2</sub>O production (from NH<sub>4</sub> and NO<sub>2</sub>). Because we don't represent NO<sub>2</sub>, this instead removes NO<sub>3</sub> <br>

Thus, the production of N<sub>2</sub>O and NO<sub>3</sub> are:

$$
\begin{align}
\mu_{aoa}^{\rightarrow N_{2}O} =& \quad \mu_{aoa}^{C} e_{aoa}^{N_{2}O} \\
\mu_{aoa}^{\rightarrow NO_3} =& \quad \mu_{aoa}^{C} e_{aoa}^{NO_{3}}
\end{align}
$$


**Anaerobic ammonia oxidizing (Anammox) bacteria**

Anammox bacteria are considered to be an implicit population within WOMBAT-mid when `do_anammox == .true.` and we do not track variations in their biomass. Rather then computing growth of anammox bacteria we therefore compute rates of anammox, which convert NH<sub>4</sub> to $N_2$. As with N<sub>2</sub>O-reducing heterotrophic bacteria, this nitrogen is then permanently lost from the ocean. We perform this metabolism as:

$$
\begin{align}
\mu_{aox}^{NH_4 \rightarrow N_2} =& \quad \mu_{aox}^{max} \left(β_{hete}\right)^{T} f_{ana} L_{aox}^{NH_4}
\end{align}
$$

_where_ <br>
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless]) <br>
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC]) <br>
- $f_{ana}$ is the binary anaerobic pathway selector for the free-living bacteria (`bacf1_fanaer(i,j,k)`, [dimenionless]) <br>
- $L_{aox}^{NH_4}$ is the growth limiter of anammox associated with NH<sub>4</sub> availability (`aox_lnh4(i,j,k)`, [dimensionless]) <br>

Note that anammox is considered to be present only when anaerobic metabolisms are ocurring. While anammox bacteria can perform anammox in oxygenated and deoxygenated environments, this metabolism is only appreciably measured in deoxygenated environments due to reduced competition with ammonia oxidizing archaea for a limited supply of NH<sub>4</sub>. Because we do not resolve this competition explicitly, we apply $f_{ana}$ here. The growth limiter due to ammonium availability is a simple michealis-menten limitation function:

$$
\begin{align}
L_{aox}^{NH_4} =& \quad \dfrac{NH_4}{NH_4 + K_{aox}^{NH_4}}
\end{align}
$$

_where_ <br>
- $K_{aoa}^{NH_4}$ is the half-saturation coefficient for uptake of NH<sub>4</sub> by anammox bacteria (`aox_knh4`, [mmol N m<sup>-3</sup>]) <br>
- NH<sub>4</sub> is the in situ concentration of NH<sub>4</sub> (`bionh4`, [mmol N m<sup>-3</sup>]) <br>

---


### 18. Tracer tendencies

**Nitrate** (`f_no3(i,j,k)`, NO<sub>3</sub>, [mol N kg<sup>-1</sup>])

$$
\begin{aligned}
\dfrac{\Delta NO_3}{\Delta t} =& \quad \mu_{aoa}^{\rightarrow NO_3} - \mu_{b-p}^{\leftarrow NO_3} - \mu_{b-f1}^{\leftarrow NO_3} \\
                               &  - \left( \mu_{np}^{\leftarrow C} \dfrac{L_{np}^{NO_3}}{L_{np}^{N}} + \mu_{mp}^{\leftarrow C} \dfrac{L_{mp}^{NO_3}}{L_{mp}^{N}} \right) \cdot \dfrac{16}{122}
\end{aligned}
$$

**Ammonium** (`f_nh4(i,j,k)`, NH<sub>4</sub>, [mol N kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta NH_4}{\Delta t} =& \quad \bigg( X_{mz}^{\leftarrow B_{np}^{N}} 
                                        + X_{mz}^{\leftarrow B_{mp}^{N}} 
                                        + X_{mz}^{\leftarrow B_{sd}^{N}} \\
                               & \qquad + X_{mz}^{\leftarrow B_{b-p}^{N}} 
                                        + X_{mz}^{\leftarrow B_{b-f1}^{N}} 
                                        + X_{mz}^{\leftarrow B_{b-f2}^{N}} 
                                        + X_{mz}^{\leftarrow B_{aoa}^{N}} \bigg) \left(1 - f_{mz}^{X \rightarrow DOM} \right) \\
                               & + \bigg( X_{Mz}^{\leftarrow B_{np}^{N}} 
                                        + X_{Mz}^{\leftarrow B_{mp}^{N}} 
                                        + X_{Mz}^{\leftarrow B_{sd}^{N}} \\
                               & \qquad + X_{Mz}^{\leftarrow B_{ld}^{N}} 
                                        + X_{Mz}^{\leftarrow B_{mz}^{N}} 
                                        + X_{Mz}^{\leftarrow B_{b-p}^{N}} \\
                               & \qquad + X_{Mz}^{\leftarrow B_{b-f1}^{N}} 
                                        + X_{Mz}^{\leftarrow B_{b-f2}^{N}} 
                                        + X_{Mz}^{\leftarrow B_{aoa}^{N}} \bigg) \cdot \left(1 - f_{Mz}^{X \rightarrow DOM} \right) \\
                               & + \bigg( \gamma_{mz}^{\rightarrow C}  
                                        + \gamma_{Mz}^{\rightarrow C} \bigg) \cdot \dfrac{16}{122} \\
                               & + \mu_{diazo}^{\rightarrow NH_4} 
                                 + \mu_{b-p}^{\rightarrow NH_4} 
                                 + \mu_{b-f1}^{\rightarrow NH_4} 
                                 + \mu_{b-f2}^{\rightarrow NH_4} \\
                               & - \mu_{aox}^{NH_4 \rightarrow N_2} 
                                 - \mu_{aoa}^{\leftarrow NH_4} \\
                               & - \bigg( \mu_{np}^{\leftarrow C} \dfrac{L_{np}^{NH_4}}{L_{np}^{N}} 
                                        + \mu_{mp}^{\leftarrow C} \dfrac{L_{mp}^{NH_4}}{L_{mp}^{N}} \bigg) \cdot \dfrac{16}{122}
\end{align}
$$

**Silicic acid** (`f_sil(i,j,k)`, $H_{4}SiO_{4}$, [mol Si kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta H_{4}SiO_{4}}{\Delta t} =& \quad \left( \gamma_{mp}^{\rightarrow C} + g_{mz}^{\leftarrow B_{mp}^{C}} \right) \cdot Q_{mp}^{Si:C} \\
                                       & + D_{B_{ld}^{Si}}^{\rightarrow Si} - \mu_{mp}^{\leftarrow Si}
\end{align}
$$

**Nitrous oxide** (`f_n2o(i,j,k)`, N<sub>2</sub>O, [mol N<sub>2</sub>O kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta N_{2}O}{\Delta t} =& \quad \mu_{aoa}^{\rightarrow N_{2}O} + \mu_{b-f1}^{\rightarrow N_{2}O} - \mu_{b-f2}^{\leftarrow N_{2}O}
\end{align}
$$

**Oxygen** (`f_o2(i,j,k)`, O<sub>2</sub>, [mol O<sub>2</sub> kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta O_2}{\Delta t} =& \quad \bigg( X_{mz}^{\leftarrow C} \left(1 - f_{mz}^{X \rightarrow DOM} \right)
                                            + X_{Mz}^{\leftarrow C} \left(1 - f_{Mz}^{X \rightarrow DOM} \right) \\
                              &             + \gamma_{mz}^{\rightarrow C}
                                            + \gamma_{Mz}^{\rightarrow C} \\
                              &             - \mu_{np}^{\leftarrow C}
                                            - \mu_{mp}^{\leftarrow C} \bigg) \dfrac{-132}{122} \\
                              & - \mu_{b-p}^{\leftarrow O_2}
                                - \mu_{b-f1}^{\leftarrow O_2}
                                - \mu_{b-f2}^{\leftarrow O_2}
                                - \mu_{aoa}^{\leftarrow O_2}
\end{align}
$$
 
**Dissolved iron** (`f_fe(i,j,k)`, $dFe$, [mol Fe kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta dFe}{\Delta t} =& \quad \Gamma_{sd}^{\rightarrow C} Q_{sd}^{Fe:C} 
                                    + \Gamma_{ld}^{\rightarrow C} Q_{ld}^{Fe:C} \\
                              &     + \gamma_{np}^{\rightarrow C} Q_{np}^{Fe:C} 
                                    + \gamma_{mp}^{\rightarrow C} Q_{mp}^{Fe:C} \\
                              &     + \gamma_{mz}^{\rightarrow C} Q_{mz}^{Fe:C} 
                                    + \gamma_{Mz}^{\rightarrow C} Q_{Mz}^{Fe:C} \\
                              &     + \dfrac{\gamma_{b-p}^{\rightarrow C} + \Gamma_{b-p}^{\rightarrow C}}{R_{b}^{C:Fe}}
                                    + \dfrac{\gamma_{b-f1}^{\rightarrow C} + \Gamma_{b-f1}^{\rightarrow C}}{R_{b}^{C:Fe}}
                                    + \dfrac{\gamma_{b-f2}^{\rightarrow C} + \Gamma_{b-f2}^{\rightarrow C}}{R_{b}^{C:Fe}}
                                    + \dfrac{\gamma_{aoa}^{\rightarrow C} + \Gamma_{aoa}^{\rightarrow C}}{R_{aoa}^{C:Fe}} \\
                              &     + X_{mz}^{\leftarrow Fe} 
                                    + X_{Mz}^{\leftarrow Fe} \\
                              &     + D_{Fe_{sA}}^{\rightarrow dFe} 
                                    + D_{Fe_{lA}}^{\rightarrow dFe} \\
                              &     - \mu_{np}^{\leftarrow dFe} 
                                    - \mu_{mp}^{\leftarrow dFe}
                                    - \mu_{b-p}^{\leftarrow dFe}
                                    - \mu_{b-f1}^{\leftarrow dFe}
                                    - \mu_{b-f2}^{\leftarrow dFe}
                                    - \mu_{aoa}^{\leftarrow dFe} \\
                              &     - Sc_{dFe}^{\rightarrow Fe_{sA}} 
                                    - Sc_{dFe}^{\rightarrow Fe_{lA}}
                                    - Co_{dFe}^{\rightarrow Fe_{sA}} 
                                    - Co_{dFe}^{\rightarrow Fe_{lA}}
\end{align}
$$

**Small authigenic iron** (`f_afe(i,j,k)`, $Fe_{sA}$, [mol Fe kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta Fe_{sA}}{\Delta t} =& \quad Sc_{dFe}^{\rightarrow Fe_{sA}}
                                        + Co_{dFe}^{\rightarrow Fe_{sA}}
                                        - D_{Fe_{sA}}^{\rightarrow dFe}
\end{align}
$$

**Large authigenic iron** (`f_bafe(i,j,k)`, $Fe_{lA}$, [mol Fe kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta Fe_{lA}}{\Delta t} =& \quad Sc_{dFe}^{\rightarrow Fe_{lA}}
                                        + Co_{dFe}^{\rightarrow Fe_{lA}}
                                        - D_{Fe_{lA}}^{\rightarrow dFe}
\end{align}
$$

**Nano-phytoplankton** (`f_phy(i,j,k)`, $B_{np}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{np}^{C}}{\Delta t} =& \quad \mu_{np}^{\leftarrow C} 
                                           - \Gamma_{np}^{\rightarrow C} 
                                           - \gamma_{np}^{\rightarrow C}
                                           - g_{mz}^{\leftarrow B_{np}^{C}} 
                                           - g_{Mz}^{\leftarrow B_{np}^{C}}
\end{align}
$$

**Nano-phytoplankton chlorophyll** (`f_pchl(i,j,k)`, $B_{np}^{chl}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{np}^{Chl}}{\Delta t} =& \quad \mu_{np}^{\leftarrow Chl}
                                       - \bigg( \Gamma_{np}^{\rightarrow C} + \gamma_{np}^{\rightarrow C} 
                                               + g_{mz}^{\leftarrow B_{np}^{C}} + g_{Mz}^{\leftarrow B_{np}^{C}} \bigg) \cdot Q_{np}^{Chl:C}
\end{align}
$$

**Nano-phytoplankton iron** (`f_phyfe(i,j,k)`, $B_{np}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{np}^{Fe}}{\Delta t} =& \quad \mu_{np}^{\leftarrow dFe} - \bigg( \Gamma_{np}^{\rightarrow C} + \gamma_{np}^{\rightarrow C} 
                                              + g_{mz}^{\leftarrow B_{np}^{C}} + g_{Mz}^{\leftarrow B_{np}^{C}} \bigg) \cdot Q_{np}^{Fe:C}
\end{align}
$$

**Micro-phytoplankton** (`f_dia(i,j,k)`, $B_{mp}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{mp}^{C}}{\Delta t} =& \quad \mu_{mp}^{\leftarrow C} 
                                           - \Gamma_{mp}^{\rightarrow C} 
                                           - \gamma_{mp}^{\rightarrow C} 
                                           - g_{mz}^{\leftarrow B_{mp}^{C}} 
                                           - g_{Mz}^{\leftarrow B_{mp}^{C}} 
\end{align}
$$

**Micro-phytoplankton chlorophyll** (`f_dchl(i,j,k)`, $B_{mp}^{chl}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{mp}^{Chl}}{\Delta t} =& \quad \mu_{mp}^{\leftarrow Chl} 
                                       - \bigg( \Gamma_{mp}^{\rightarrow C} 
                                              + \gamma_{mp}^{\rightarrow C} 
                                              + g_{mz}^{\leftarrow B_{mp}^{C}} 
                                              + g_{Mz}^{\leftarrow B_{mp}^{C}} \bigg) \cdot Q_{mp}^{Chl:C}
\end{align}
$$

**Micro-phytoplankton iron** (`f_diafe(i,j,k)`, $B_{mp}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{mp}^{Fe}}{\Delta t} =& \quad \mu_{mp}^{\leftarrow dFe} 
                                      - \bigg( \Gamma_{mp}^{\rightarrow C} 
                                             + \gamma_{mp}^{\rightarrow C} 
                                             + g_{mz}^{\leftarrow B_{mp}^{C}} 
                                             + g_{Mz}^{\leftarrow B_{mp}^{C}} \bigg) \cdot Q_{mp}^{Fe:C}
\end{align}
$$

**Micro-phytoplankton silica** (`f_diasi(i,j,k)`, $B_{mp}^{Si}$, [mol Si kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{mp}^{Si}}{\Delta t} =& \quad \mu_{mp}^{\leftarrow Si}
                                      - \bigg( \Gamma_{mp}^{\rightarrow C} 
                                             + \gamma_{mp}^{\rightarrow C} 
                                             + g_{mz}^{\leftarrow B_{mp}^{C}} 
                                             + g_{Mz}^{\leftarrow B_{mp}^{C}} \bigg) \cdot Q_{mp}^{Si:C}
\end{align}
$$

**Micro-zooplankton** (`f_zoo(i,j,k)`, $B_{mz}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{mz}^{C}}{\Delta t} =& \quad A_{mz}^{\leftarrow C} 
                                           - \Gamma_{mz}^{\rightarrow C} 
                                           - \gamma_{mz}^{\rightarrow C}
                                           - g_{Mz}^{\leftarrow B_{mz}^{C}}
\end{align}
$$

**Micro-zooplankton iron** (`f_zoofe(i,j,k)`, $B_{mz}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{mz}^{Fe}}{\Delta t} =& \quad A_{mz}^{\leftarrow Fe}
                                      - \bigg( \Gamma_{mz}^{\rightarrow C} 
                                             + \gamma_{mz}^{\rightarrow C} 
                                             + g_{Mz}^{\leftarrow B_{mz}^{C}} \bigg) \cdot Q_{mz}^{Fe:C}
\end{align}
$$

**Meso-zooplankton** (`f_mes(i,j,k)`, $B_{Mz}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{Mz}^{C}}{\Delta t} =& \quad A_{Mz}^{\leftarrow C}
                                           - \Gamma_{Mz}^{\rightarrow C} 
                                           - \gamma_{Mz}^{\rightarrow C}
\end{align}
$$

**Meso-zooplankton iron** (`f_mesfe(i,j,k)`, $B_{Mz}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{Mz}^{Fe}}{\Delta t} =& \quad A_{Mz}^{\leftarrow Fe}
                                      - \bigg( \Gamma_{Mz}^{\rightarrow C} 
                                             + \gamma_{Mz}^{\rightarrow C} \bigg) \cdot Q_{Mz}^{Fe:C}
\end{align}
$$

**Small detritus** (`f_det(i,j,k)`, $B_{sd}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{sd}^{C}}{\Delta t} =& \quad E_{mz}^{\leftarrow C}
                                           + \Gamma_{np}^{\rightarrow C} 
                                           + \Gamma_{mz}^{\rightarrow C}
                                           - g_{mz}^{\leftarrow B_{sd}^{C}}
                                           - g_{Mz}^{\leftarrow B_{sd}^{C}}
                                           - \Gamma_{sd}^{\rightarrow C}
\end{align}
$$

**Small detritus iron** (`f_detfe(i,j,k)`, $B_{sd}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{sd}^{Fe}}{\Delta t} = \quad E_{mz}^{\leftarrow Fe}
                                      + \Gamma_{np}^{\rightarrow C} Q_{np}^{Fe:C} 
                                      + \Gamma_{mz}^{\rightarrow C} Q_{mz}^{Fe:C}
                                      - \left( g_{mz}^{\leftarrow B_{sd}^{C}}
                                             + g_{Mz}^{\leftarrow B_{sd}^{C}}
                                             + \Gamma_{sd}^{\rightarrow C} \right) Q_{sd}^{Fe:C}
\end{align}
$$

**Large detritus** (`f_bdet(i,j,k)`, $B_{ld}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{ld}^{C}}{\Delta t} =& \quad E_{Mz}^{\leftarrow C}
                                           + \Gamma_{mp}^{\rightarrow C} 
                                           + \Gamma_{Mz}^{\rightarrow C}
                                           - g_{Mz}^{\leftarrow B_{ld}^{C}}
                                           - \Gamma_{ld}^{\rightarrow C}
\end{align}
$$

**Large detritus iron** (`f_bdetfe(i,j,k)`, $B_{ld}^{Fe}$, [mol Fe kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{ld}^{Fe}}{\Delta t} =& \quad E_{Mz}^{\leftarrow Fe}
                                            + \Gamma_{mp}^{\rightarrow C} Q_{mp}^{Fe:C} 
                                            + \Gamma_{Mz}^{\rightarrow C} Q_{Mz}^{Fe:C}
                                            - \bigg( g_{Mz}^{\leftarrow B_{ld}^{C}} + \Gamma_{ld}^{\rightarrow C} \bigg) Q_{ld}^{Fe:C}
\end{align}
$$

**Large detritus silicon** (`f_bdetsi(i,j,k)`, $B_{ld}^{Si}$, [mol Si kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{ld}^{Si}}{\Delta t} =& \quad \bigg( \Gamma_{mp}^{\rightarrow C} 
                                                    + g_{Mz}^{\leftarrow B_{mp}^{C}} \bigg) Q_{mp}^{Si:C}  
                                      - D_{B_{ld}^{Si}}^{\rightarrow Si} 
\end{align}
$$

**Faculative NO<sub>3</sub>-reducing, particle-associated heterotrophic bacteria** (`f_bacp(i,j,k)`, $B_{b-p}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{b-p}^{C}}{\Delta t} =& \quad \mu_{b-p}^{\leftarrow C} 
                                            - g_{mz}^{\leftarrow B_{b-p}^{C}}
                                            - g_{Mz}^{\leftarrow B_{b-p}^{C}}
                                            - \gamma_{b-p}^{\rightarrow C} 
                                            - \Gamma_{b-p}^{\rightarrow C}
\end{align}
$$

**Faculative NO<sub>3</sub>-reducing, free-living heterotrophic bacteria** (`f_bacf1(i,j,k)`, $B_{b-f1}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{b-f1}^{C}}{\Delta t} =& \quad \mu_{b-f1}^{\leftarrow C} 
                                           - g_{mz}^{\leftarrow B_{b-f1}^{C}}
                                           - g_{Mz}^{\leftarrow B_{b-f1}^{C}}
                                           - \gamma_{b-f1}^{\rightarrow C} 
                                           - \Gamma_{b-f1}^{\rightarrow C}
\end{align}
$$

**Faculative N<sub>2</sub>O-reducing, free-living heterotrophic bacteria** (`f_bacf2(i,j,k)`, $B_{b-f2}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{b-f2}^{C}}{\Delta t} =& \quad \mu_{b-f2}^{\leftarrow C} 
                                           - g_{mz}^{\leftarrow B_{b-f2}^{C}}
                                           - g_{Mz}^{\leftarrow B_{b-f2}^{C}}
                                           - \gamma_{b-f2}^{\rightarrow C} 
                                           - \Gamma_{b-f2}^{\rightarrow C}
\end{align}
$$

**Ammonia oxidizing archaea** (`f_aoa(i,j,k)`, $B_{aoa}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{aoa}^{C}}{\Delta t} =& \quad \mu_{aoa}^{\leftarrow C} 
                                            - g_{mz}^{\leftarrow B_{aoa}^{C}}
                                            - g_{Mz}^{\leftarrow B_{aoa}^{C}}
                                            - \gamma_{aoa}^{\rightarrow C} 
                                            - \Gamma_{aoa}^{\rightarrow C}
\end{align}
$$

**Dissolved organic carbon** (`f_doc(i,j,k)`, $B_{DOM}^{C}$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{DOM}^{C}}{\Delta t} =& \quad \mu_{np}^{\rightarrow DOC} 
                                            + \mu_{mp}^{\rightarrow DOC} \\
                                      &     + \gamma_{np}^{\rightarrow C}
                                            + \gamma_{mp}^{\rightarrow C}
                                            + \gamma_{b-p}^{\rightarrow C}
                                            + \gamma_{b-f1}^{\rightarrow C}
                                            + \gamma_{b-f2}^{\rightarrow C} \\
                                      &     + \Gamma_{b-p}^{\rightarrow C}
                                            + \Gamma_{b-f1}^{\rightarrow C}
                                            + \Gamma_{b-f2}^{\rightarrow C}
                                            + \Gamma_{aoa}^{\rightarrow C}
                                            + \gamma_{aoa}^{\rightarrow C} \\
                                      &     + X_{mz}^{\leftarrow C} f_{mz}^{X \rightarrow DOM} 
                                            + X_{Mz}^{\leftarrow C} f_{Mz}^{X \rightarrow DOM} \\
                                      &     + \mu_{b-p}^{\rightarrow B_{DOM}^{C}}
                                            + \mu_{b-f1}^{\rightarrow B_{DOM}^{C}}
                                            + \mu_{b-f2}^{\rightarrow B_{DOM}^{C}} \\
                                      &     - \mu_{b-f1}^{\leftarrow B_{DOM}^{C}}
                                            - \mu_{b-f2}^{\leftarrow B_{DOM}^{C}}
\end{align}
$$

**Dissolved organic hydrogen** (`f_doh(i,j,k)`, $B_{DOM}^{H}$, [mol H kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{DOM}^{H}}{\Delta t} =& \quad \bigg( \mu_{np}^{\rightarrow DOC}
                                                   + \mu_{mp}^{\rightarrow DOC} \bigg) \cdot 2.0 \\
                                      &     + \bigg( \gamma_{np}^{\rightarrow C}
                                                   + \gamma_{mp}^{\rightarrow C} \\
                                      &            + X_{mz}^{\leftarrow C} f_{mz}^{X \rightarrow DOM} 
                                                   + X_{Mz}^{\leftarrow C} f_{Mz}^{X \rightarrow DOM} \bigg) \cdot 1.65 \\
                                      &     + \bigg( \gamma_{b-p}^{\rightarrow C}
                                                   + \gamma_{b-f1}^{\rightarrow C}
                                                   + \gamma_{b-f2}^{\rightarrow C} \\
                                      &            + \Gamma_{b-p}^{\rightarrow C}
                                                   + \Gamma_{b-f1}^{\rightarrow C}
                                                   + \Gamma_{b-f2}^{\rightarrow C} \\
                                      &            + \Gamma_{aoa}^{\rightarrow C}
                                                   + \gamma_{aoa}^{\rightarrow C} \bigg) \cdot 1.40 \\
                                      &     + \bigg( \mu_{b-p}^{\rightarrow B_{DOM}^{C}} \cdot 1.65 
                                                   + \mu_{b-f1}^{\rightarrow B_{DOM}^{C}} \cdot \dfrac{B_{DOM}^{H}}{B_{DOM}^{C}}
                                                   + \mu_{b-f2}^{\rightarrow B_{DOM}^{C}} \cdot \dfrac{B_{DOM}^{H}}{B_{DOM}^{C}} \bigg) \cdot H_{ox}\\
                                      &     - \bigg( \mu_{b-f1}^{\leftarrow B_{DOM}^{C}} 
                                                   + \mu_{b-f2}^{\leftarrow B_{DOM}^{C}} \bigg) \dfrac{B_{DOM}^{H}}{B_{DOM}^{C}}
 
\end{align}
$$

**Dissolved organic oxygen** (`f_doo(i,j,k)`, $B_{DOM}^{O}$, [mol O kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{DOM}^{O}}{\Delta t} =& \quad \bigg( \mu_{np}^{\rightarrow DOC}
                                                   + \mu_{mp}^{\rightarrow DOC} \bigg) \cdot 1.0 \\
                                      &     + \bigg( \gamma_{np}^{\rightarrow C}
                                                   + \gamma_{mp}^{\rightarrow C} \\
                                      &            + X_{mz}^{\leftarrow C} f_{mz}^{X \rightarrow DOM} 
                                                   + X_{Mz}^{\leftarrow C} f_{Mz}^{X \rightarrow DOM} \bigg) \cdot 0.40 \\
                                      &     + \bigg( \gamma_{b-p}^{\rightarrow C}
                                                   + \gamma_{b-f1}^{\rightarrow C}
                                                   + \gamma_{b-f2}^{\rightarrow C} \\
                                      &            + \Gamma_{b-p}^{\rightarrow C}
                                                   + \Gamma_{b-f1}^{\rightarrow C}
                                                   + \Gamma_{b-f2}^{\rightarrow C} \\
                                      &            + \Gamma_{aoa}^{\rightarrow C}
                                                   + \gamma_{aoa}^{\rightarrow C} \bigg) \cdot 0.40 \\
                                      &     + \bigg( \mu_{b-p}^{\rightarrow B_{DOM}^{C}} \cdot 0.4 
                                                   + \mu_{b-f1}^{\rightarrow B_{DOM}^{C}} \cdot \dfrac{B_{DOM}^{O}}{B_{DOM}^{C}}
                                                   + \mu_{b-f2}^{\rightarrow B_{DOM}^{C}} \cdot \dfrac{B_{DOM}^{O}}{B_{DOM}^{C}} \bigg) \cdot O_{ox}\\
                                      &     - \bigg( \mu_{b-f1}^{\leftarrow B_{DOM}^{C}} 
                                                   + \mu_{b-f2}^{\leftarrow B_{DOM}^{C}} \bigg) \dfrac{B_{DOM}^{O}}{B_{DOM}^{C}}
 
\end{align}
$$

**Dissolved organic nitrogen** (`f_don(i,j,k)`, $B_{DOM}^{N}$, [mol N kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta B_{DOM}^{N}}{\Delta t} =& \quad \bigg( \Gamma_{b-p}^{\rightarrow C} + \gamma_{b-p}^{\rightarrow C}
                                                   + \Gamma_{b-f1}^{\rightarrow C} + \gamma_{b-f1}^{\rightarrow C}
                                                   + \Gamma_{b-f2}^{\rightarrow C} + \gamma_{b-f2}^{\rightarrow C}
                                                   + \Gamma_{aoa}^{\rightarrow C} + \gamma_{aoa}^{\rightarrow C} \\
                                      &            + X_{mz}^{\leftarrow B_{b-p}^{C}} + X_{mz}^{\leftarrow B_{b-f1}^{C}} + X_{mz}^{\leftarrow B_{b-f2}^{C}}
                                                   + X_{mz}^{\leftarrow B_{aoa}^{C}} + X_{Mz}^{\leftarrow B_{b-p}^{C}} + X_{Mz}^{\leftarrow B_{b-f1}^{C}}
                                                   + X_{Mz}^{\leftarrow B_{b-f2}^{C}} + X_{Mz}^{\leftarrow B_{aoa}^{C}} \bigg) \cdot \dfrac{1}{R_{b}^{C:N}} \\
                                      &  \quad + \Bigg( \bigg( X_{mz}^{\leftarrow B_{np}^{C}} + X_{mz}^{\leftarrow B_{mp}^{C}}
                                                             + X_{mz}^{\leftarrow B_{sd}^{C}} \bigg) f_{mz}^{X \rightarrow DOM} \\
                                      &  \qquad       + \bigg( X_{Mz}^{\leftarrow B_{np}^{C}} + X_{Mz}^{\leftarrow B_{mp}^{C}}
                                                             + X_{Mz}^{\leftarrow B_{sd}^{C}} + X_{Mz}^{\leftarrow B_{ld}^{C}}
                                                             + X_{Mz}^{\leftarrow B_{mz}^{C}} \bigg) f_{Mz}^{X \rightarrow DOM} \Bigg) \cdot \dfrac{16}{122} \\
                                      &     + \bigg( \mu_{b-p}^{\rightarrow B_{DOM}^{C}} \cdot \dfrac{16}{122} 
                                                   + \mu_{b-f1}^{\rightarrow B_{DOM}^{C}} \cdot \dfrac{B_{DOM}^{N}}{B_{DOM}^{C}}
                                                   + \mu_{b-f2}^{\rightarrow B_{DOM}^{C}} \cdot \dfrac{B_{DOM}^{N}}{B_{DOM}^{C}} \bigg) \cdot N_{ox}\\
                                      &     - \bigg( \mu_{b-f1}^{\leftarrow B_{DOM}^{C}} 
                                                   + \mu_{b-f2}^{\leftarrow B_{DOM}^{C}} \bigg) \dfrac{B_{DOM}^{N}}{B_{DOM}^{C}}
 
\end{align}
$$

**Calcium Carbonate** (`f_caco3(i,j,k)`, $CaCO_3$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta CaCO_3}{\Delta t} =& \quad P_{CaCO_3}^{\Gamma_{np}^{\rightarrow C}} 
                                       + P_{CaCO_3}^{\Gamma_{mz}^{\rightarrow C}}
                                       + P_{CaCO_3}^{g_{mz}^{\leftarrow B_{np}^{C}}}
                                       + P_{CaCO_3}^{g_{Mz}^{\leftarrow B_{np}^{C}}}
                                       + P_{CaCO_3}^{g_{Mz}^{\leftarrow B_{mz}^{C}}} \\
                                 &     - D_{CaCO_3}^{\Omega_{cal}}
                                       - D_{CaCO_3}^{\Omega_{ara}}
                                       - D_{CaCO_3}^{\Gamma_{sd}^{\rightarrow C}}
                                       - D_{CaCO_3}^{g_{mz}^{\leftarrow B_{sd}^{C}}}
                                       - D_{CaCO_3}^{g_{Mz}^{\leftarrow B_{sd}^{C}}}
\end{align}
$$

**Dissolved Inorganic Carbon** (`f_dic(i,j,k)`, $DIC$, [mol C kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta DIC}{\Delta t} =& \quad X_{mz}^{\leftarrow C} \left(1 - f_{mz}^{X \rightarrow DOM} \right)
                                    + X_{Mz}^{\leftarrow C} \left(1 - f_{Mz}^{X \rightarrow DOM} \right) \\
                              &     + \gamma_{mz}^{\rightarrow C}
                                    + \gamma_{Mz}^{\rightarrow C}
                                    + \mu_{b-p}^{\rightarrow DIC}
                                    + \mu_{b-f1}^{\rightarrow DIC}
                                    + \mu_{b-f2}^{\rightarrow DIC} \\
                              &     + D_{CaCO_3}^{\Omega_{cal}}
                                    + D_{CaCO_3}^{\Omega_{ara}}
                                    + D_{CaCO_3}^{\Gamma_{sd}^{\rightarrow C}}
                                    + D_{CaCO_3}^{g_{mz}^{\leftarrow B_{sd}^{C}}}
                                    + D_{CaCO_3}^{g_{Mz}^{\leftarrow B_{sd}^{C}}} \\
                              &     - \mu_{np}^{\leftarrow C}
                                    - \mu_{mp}^{\leftarrow C}
                                    - \mu_{aoa}^{\leftarrow C}
                                    - \mu_{np}^{\rightarrow DOC}
                                    - \mu_{mp}^{\rightarrow DOC} \\
                              &     - P_{CaCO_3}^{\Gamma_{np}^{\rightarrow C}}
                                    - P_{CaCO_3}^{\Gamma_{mz}^{\rightarrow C}}
                                    - P_{CaCO_3}^{g_{mz}^{\leftarrow B_{np}^{C}}}
                                    - P_{CaCO_3}^{g_{Mz}^{\leftarrow B_{np}^{C}}}
                                    - P_{CaCO_3}^{g_{Mz}^{\leftarrow B_{mz}^{C}}}
\end{align}
$$


**Alkalinity** (`f_alk(i,j,k)`, $Alk$, [mol Eq kg<sup>-1</sup>])

$$
\begin{align}
\dfrac{\Delta Alk}{\Delta t} =& \quad \bigg( \mu_{np}^{\leftarrow C} \cdot \dfrac{L_{np}^{NO_3}}{L_{np}^{N}} 
                                           + \mu_{mp}^{\leftarrow C} \cdot \dfrac{L_{mp}^{NO_3}}{L_{mp}^{N}} 
                                           - \mu_{np}^{\leftarrow C} \cdot \dfrac{L_{np}^{NH_4}}{L_{np}^{N}} \\
                              &            - \mu_{mp}^{\leftarrow C} \cdot \dfrac{L_{mp}^{NH_4}}{L_{mp}^{N}} 
                                           + \gamma_{mz}^{\rightarrow C} 
                                           + \gamma_{Mz}^{\rightarrow C} \bigg) \cdot \dfrac{16}{122} \\
                              &     + \bigg( X_{mz}^{\leftarrow B_{np}^{N}}
                                           + X_{mz}^{\leftarrow B_{mp}^{N}}
                                           + X_{mz}^{\leftarrow B_{sd}^{N}} \\
                              &            + X_{mz}^{\leftarrow B_{b-p}^{N}}
                                           + X_{mz}^{\leftarrow B_{b-f1}^{N}}
                                           + X_{mz}^{\leftarrow B_{b-f2}^{N}}
                                           + X_{mz}^{\leftarrow B_{aoa}^{N}} \bigg) \left(1 - f_{mz}^{X \rightarrow DOM} \right) \\
                              &     + \bigg( X_{Mz}^{\leftarrow B_{np}^{N}}
                                           + X_{Mz}^{\leftarrow B_{mp}^{N}}
                                           + X_{Mz}^{\leftarrow B_{sd}^{N}} \\
                              &            + X_{Mz}^{\leftarrow B_{ld}^{N}}
                                           + X_{Mz}^{\leftarrow B_{mz}^{N}}
                                           + X_{Mz}^{\leftarrow B_{b-p}^{N}} \\
                              &            + X_{Mz}^{\leftarrow B_{b-f1}^{N}}
                                           + X_{Mz}^{\leftarrow B_{b-f2}^{N}}
                                           + X_{Mz}^{\leftarrow B_{aoa}^{N}} \bigg) \cdot \left(1 - f_{Mz}^{X \rightarrow DOM} \right) \\
                              &            + \mu_{b-p}^{\rightarrow NH_4} 
                                           + \mu_{b-f1}^{\rightarrow NH_4}
                                           + \mu_{b-f2}^{\rightarrow NH_4}
                                           + \mu_{b-p}^{\leftarrow NO_{3}}
                                           + \mu_{b-f1}^{\leftarrow NO_{3}}
                                           - \mu_{aoa}^{\leftarrow NH_{4}}
                                           - \mu_{aoa}^{\rightarrow NO_{3}}
                                           - \mu_{aox}^{NH_4 \rightarrow N_2} \\
                              &            + \bigg( D_{CaCO_3}^{\Omega_{cal}}
                                           + D_{CaCO_3}^{\Omega_{ara}}
                                           + D_{CaCO_3}^{\Gamma_{sd}^{\rightarrow C}}
                                           + D_{CaCO_3}^{g_{mz}^{\leftarrow B_{sd}^{C}}}
                                           + D_{CaCO_3}^{g_{Mz}^{\leftarrow B_{sd}^{C}}} \\
                              &            - P_{CaCO_3}^{\Gamma_{np}^{\rightarrow C}}
                                           - P_{CaCO_3}^{\Gamma_{mz}^{\rightarrow C}}
                                           - P_{CaCO_3}^{g_{mz}^{\leftarrow B_{np}^{C}}}
                                           - P_{CaCO_3}^{g_{Mz}^{\leftarrow B_{np}^{C}}}
                                           - P_{CaCO_3}^{g_{Mz}^{\leftarrow B_{mz}^{C}}} \bigg) \cdot 2
\end{align}
$$


---


### 19. Check for conservation of mass

When checks for the conservation of mass is enabled (`do_check_n_conserve = .true.` or `do_check_c_conserve = .true.` or `do_check_si_conserve = .true.`), the model will calculate the budget of nitrogen or carbon or silicon before and after the ecosystem equations have completed. This checks that the ecosystem equations detailed above have indeed conserved the mass of both nitrogen and carbon within the ocean. In WOMBAT-mid, both nitrogen and carbon and silicon should be perfectly conserved during ecosystem cycling. The only exception to this is for nitrogen, where if any of `do_nitrogen_fixation = .true.`, `do_anammox = .true.`, `do_wc_denitrification = .true.` or `do_benthic_denitrification = .true.` then the model does not and should not be expected to conserve nitrogen.

---

### 20. Additional operations on tracers

**First**, dissolved iron concentrations are set to equal 1 nM everywhere where the depth of the water column is less than 200 metres deep. WOMBAT-mid is not considered to be a model of the coastal ocean, but rather a model of the global pelagic ocean. Given that coastal waters are not limited in dissolved iron due to substantial interactions with sediments and exchange with the land, we universally set the dissolved iron concentration in these waters to 1 nM.

**Second**, if dissolved iron concentrations dip below that measureable by operational detection limits considered to be roughlly 50 pM ([Worsford et al., 2014](https://doi.org/10.1016/j.marchem.2014.08.009)), we reset these concentrations to this minimum (`zfermin`, $[dFe]^{min}$, [µmol m<sup>-3</sup>]):

$$
\begin{align}
[dFe]^{min} =& \quad 0.05
\end{align}
$$

This resetting of minimum dFe concentration essentially copies what is done in the PISCES ocean model and functions as a constant source of dFe to the ocean when surface concentrations are drawn down to near zero values. Ideally, complexation by ligands would function to maintain iron in dissolved, biologically available form without the need for an addition source at low concentrations.

---


### 21. Sinking rate of particulates.

WOMBAT-mid functions with a spatially variable sinking rate of organic detritus (`f_det(i,j,k)`; `f_bdet(i,j,k)`), calcium carbonate (`f_caco3(i,j,k)`) and biogenic silica (`f_bdetsi(i,j,k)`). Sinking of organic iron (`f_detfe(i,j,k)`; `f_bdetfe(i,j,k)`)occurs at the same rate as their respective organic particulate carbon types, while small and large authigenic iron particles (`f_afe(i,j,k)`; `f_bafe(i,j,k)`) sink at their own unique rates. The algorithm to compute sinking rates functions by computing:

1. the average radii of particles in the community;
2. the seawater dynamic viscosity (if `do_viscous_sinking =.true.`);
3. the effect of mineral ballasting, specifically $CaCO_3$ and biogenic silica, on particulate excess density;
4. Rubey's equation for sinking rates of particles.

This approach is inspired by the study of [Dinauer et al. (2022)](https://doi.org/10.1029/2021GB007131). We deal with each of these steps below.

**Average radii of particulates**

We first estimate the average radius of sinking particles belonging to small and large detritus. Nano-phytoplankton and micro-zooplankton contribute to the small detritus pool, and as such variations in the mean size of these plankton types determine the mean radius of small particles. Similarly, micro-phytoplankton and meso-zooplankton contribute to the large detritus pool, and their sizes determine the mean radius of large particles. According to [Wickman et al. (2024)]((https://doi.org/10.1126/science.adk6901)), the average volume, $V$, of phytoplankton, $p$, in the marine community scales with the biomass density according to:

$$
\begin{align}
V =& \quad \left(B_{p}^{C}\right)^{0.65}
\end{align}
$$

We can relate the radius $r$ in units of µm to the volume of phytoplankton cells via:

$$
\begin{align}
r =& \quad 0.5 \cdot \left(\dfrac{6}{\pi} \left(B_{p}^{C}\right)^{0.65} \right)^{\dfrac{1}{3}}
\end{align}
$$

Therefore, the average radius of nano-phytoplankton (`rad_phy`, [m]) and micro-phytoplankton (`rad_dia`, [m]) is equal to:

$$
\begin{align}
r_{np} =& \quad r_{np}^{-} \cdot 0.5 \cdot \left(\dfrac{6}{\pi} \left(B_{np}^{C}\right)^{0.65} \right)^{\dfrac{1}{3}} \cdot 1 \times 10^{-6} \\
r_{mp} =& \quad r_{mp}^{-} \cdot 0.5 \cdot \left(\dfrac{6}{\pi} \left(B_{mp}^{C}\right)^{0.65} \right)^{\dfrac{1}{3}} \cdot 1 \times 10^{-6}
\end{align}
$$

which simplifies to:

$$
\begin{align}
r_{np} =& \quad r_{np}^{-} \left(B_{np}^{C}\right)^{0.217} \cdot 0.62 \times 10^{-6} \\
r_{mp} =& \quad r_{mp}^{-} \left(B_{mp}^{C}\right)^{0.217} \cdot 0.62 \times 10^{-6}
\end{align}
$$

For micro-zooplankton, we use a relationship presented by [Menden-Deuer & Lessard (2000)](https://doi.org/10.4319/lo.2000.45.3.0569) who identified that the carbon concentration of diverse protists, including heterotrophic dinoflagellates and other micro-zooplankton, was related to their cell volume to the power of 0.939. Hence, we estimate the radius of micro-zooplankton (`rad_zoo`, [m]) from their carbon biomass concentration by inverting this exponent:

$$
\begin{align}
r_{mz} =& \quad r_{mz}^{-} \cdot 0.5 \cdot \left(\dfrac{6}{\pi} \left(B_{mz}^{C}\right)^{1.065} \right)^{\dfrac{1}{3}} \cdot 1 \times 10^{-6}
\end{align}
$$

which simplifies to:

$$
\begin{align}
r_{mz} =& \quad r_{mz}^{-} \left(B_{mz}^{C}\right)^{0.355} \cdot 0.62 \times 10^{-6}
\end{align}
$$

For meso-zooplankton, we assume that the dry carbon biomass scales with the body length to the power of 3, such that $B_{Mz}^{C} \propto L^{3}$ ([Uye, 1982](https://doi.org/10.1007/BF02110286); [Lehette & Hernandez-Leon, 2009](https://doi.org/10.4319/lom.2009.7.304)). Thus, $L \propto \left(B_{Mz}^{C}\right)^{\dfrac{1}{3}}$, and:

$$
\begin{align}
r_{Mz} =& \quad r_{Mz}^{-} \cdot 0.5 \cdot \left(B_{Mz}^{C}\right)^{\dfrac{1}{3}} \cdot 1 \times 10^{-6}
\end{align}
$$

We determine the mean radius of small (`rad_det`, [m]) and large particulate detritus (`rad_bdet`, [m]) by taking the biomass-weighted means of each plankton functional type:

$$
\begin{align}
r_{s} =& \quad \dfrac{B_{np}^{C} r_{np} + B_{mz}^{C} r_{mz}}{B_{np}^{C} + B_{mz}^{C}} \\
r_{l} =& \quad \dfrac{B_{mp}^{C} r_{mp} + B_{Mz}^{C} r_{Mz}}{B_{mp}^{C} + B_{Mz}^{C}}
\end{align}
$$

_where_ <br>
- $B_{np}^{C}$ is the concentration of nano-phytoplankton biomass at the surface of the water column (`biophy`, [mmol C m<sup>-3</sup>]) <br>
- $B_{mp}^{C}$ is the concentration of micro-phytoplankton biomass at the surface of the water column (`biodia`, [mmol C m<sup>-3</sup>]) <br>
- $B_{mz}^{C}$ is the concentration of micro-zooplankton biomass at the surface of the water column (`biozoo`, [mmol C m<sup>-3</sup>]) <br>
- $B_{Mz}^{C}$ is the concentration of meso-zooplankton biomass at the surface of the water column (`biomes`, [mmol C m<sup>-3</sup>]) <br>


**Seawater dynamic viscosity**

If `do_viscous_sinking = .true.`, we calculate the dynamic viscosity of the in situ seawater (`dynvis_sw(i,j,k)`, $\eta_{sw}$, [kg m<sup>-1</sup> s<sup>-1</sup>]) that particulates are sinking through. This involves three steps and is dependent on temperature, salinity and pressure.

The dynamic viscosity of seawater at atmospheric pressure ($\eta_{sw}^{1atm}$) is described by equations 22 and 23 from [Sharqawy et al. (2010)](https://doi.org/10.5004/dwt.2010.1079), which are based on [Isdale et al. (1972)](https://doi.org/10.1016/S0011-9164(00)80002-8):

$$
\begin{align}
\eta_{sw}^{1atm} =& \quad \eta_{w}^{1atm} \left(1 + A \dfrac{S}{1000} + B \left(\dfrac{S}{1000}\right)^{2} \right)
\end{align}
$$

_where_ <br>

$$
\begin{align}
A =& \quad 1.541 + 1.998 \times 10^{-2} \cdot T - 9.52 \times 10^{-5} \cdot T^{2} \\
B =& \quad 7.974 - 7.561 \times 10^{-2} \cdot T + 4.724 \times 10^{-4} \cdot T^{2} \\
\eta_{w}^{1atm} =& \quad 4.2844 \times 10^{-5} + \left(0.157 \left( T + 64.993 \right)^{2} - 91.296 \right)^{-1}
\end{align}
$$

_where_ <br>
- T is in situ seawater temperature (`Temp(i,j,k)`, [ºC]) <br>
- S is in situ seawater salinity (`Salt(i,j,k)`, [g kg<sup>-1</sup>]) <br>

After calculating $\eta_{sw}^{1atm}$, we must correct for pressure changes in the water column. This is done by calculating the effect of pressure and temperature changes on the dynamic viscosity of pure water ($\eta_{w}$), and then applying this correction to our estimate of seawater dynamic viscosity at atmospheric pressure, such that:

$$
\begin{align}
\eta_{sw} =& \quad \eta_{sw}^{1atm} \dfrac{\eta_{w}}{\eta_{w}^{1atm}}
\end{align}
$$

We solve for $\eta_{w}$, the dynamic viscosity of pure water corrected for pressure effects, by following the [IAPWS (2008)](https://iapws.org/documents/release/viscosity). This formulation requires multiple steps. First, we estimate the density of pure water, $\rho_{w}$, at a given temperature $T$ and pressure $P$ using equation 14 from the UNESCO EOS-80:

$$
\begin{align}
P_{MPa} =& \quad \left(101325 + 9.81 * 1025.0 * z\right) 1 \times 10^{-6} \\
\rho_{0} =& \quad 999.842594 + 6.793952 \times 10^{-2} \cdot T \\
          &     - 9.095290 \times 10^{-3} \cdot T^{2} \\
          &     + 1.001685 \times 10^{-4} \cdot T^{3} \\
          &     - 1.120083 \times 10^{-6} \cdot T^{4} \\
          &     + 6.536332 \times 10^{-9} \cdot T^{5} \\
\rho_{w} =& \quad \dfrac{\rho_{0}}{1 - \dfrac{P_{MPa} - 0.101325}{2.2 \times 10^{3}}}
\end{align}
$$

_where_ <br>
- $P_{MPa}$ is the in situ pressure at a given depth [`P_MPa`, [MPa]] <br>
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC]) <br>

Next, we solve for the dynamic viscosity at the dilute gas-limit, $\hat{\eta_{0}}$, detailed in equation 11 in the [IAPWS (2008)](https://iapws.org/documents/release/viscosity) and with $H_{i}$ coefficients detailed in their Table 1.

$$
\begin{align}
\hat{\eta_{0}} =& \quad \dfrac{100 \sqrt{\hat{T}}}{\sum_{i=0}^{3} \dfrac{H_{i}}{\left(\hat{T}\right)^{i}}}
\end{align}
$$

_where_ <br>
- $\hat{T} = \dfrac{T + 273.15}{647.096}$ (`T_hat`, [dimenionless]) <br>

Next, we solve for the contribution of finite density to dynamic viscosity, $\hat{\eta_{1}}$, detailed in equation 12 in the [IAPWS (2008)](https://iapws.org/documents/release/viscosity) and with $H_{ij}$ coefficients detailed in their Table 2.

$$
\begin{align}
\hat{\eta_{1}} =& \quad e^{\left[ \hat{\rho} \sum_{i=0}^{5}\left( \dfrac{1}{\hat{T} - 1} \right)^{i} \sum_{j=0}^{6} H_{ij}\left(\hat{\rho} - 1\right)^{j} \right]}
\end{align}
$$
 
_where_ <br>
- $\hat{\rho} = \dfrac{\rho_{w}}{322}$ (`rho_hat`, [dimensionless]) <br>

Finally, we compute the density-corrected pure water dynamic viscosity:

$$
\begin{align}
\eta_{w} =& \quad \left(\hat{\eta_{0}} \cdot \hat{\eta_{1}} \right) \eta^{*}
\end{align}
$$

_where_ <br>
- $\eta^{*} = 1 \times 10^{-6}$ (`mu_star`, [kg m<sup>-1</sup> s<sup>-1</sup>]) <br>

which we apply above to calculate the dynamic viscosity of seawater ($\eta_{sw}$) for a given temperature, salinity and pressure. We note that it is expected that the dynamic viscosity of water actually decreases with increasing pressure at low temperatures ([Percy W. Bridgman (1925)](https://api.semanticscholar.org/CorpusID:27500909)). 


**Mineral ballasting and excess density**

WOMBAT-mid explicitly considers small organic carbon, large aggregates of organic carbon, $CaCO_3$ and biogenic silica. Each of these particulate types have unique densities. We compute the mass of each particulate type in [kg m</sup>-3</sup>]:

$$
\begin{align}
M_{sd} =& \quad B_{sd}^{C} \cdot 1 \times 10^{-6} \cdot \dfrac{12}{0.4} \\
M_{ld} =& \quad B_{ld}^{C} \cdot 1 \times 10^{-6} \cdot \dfrac{12}{0.4} \\
M_{CaCO_3} =& \quad B_{CaCO_3}^{C} \cdot 1 \times 10^{-6} \cdot 100 \\
M_{BSi} =& \quad B_{ld}^{Si} \cdot 1 \times 10^{-6} \cdot 60
\end{align}
$$

_where_ <br>
- $B_{sd}^{C}$ is the in situ concentration of small particulate organic carbon (`biodet`, [mmol C m<sup>-3</sup>]) <br>
- $B_{ld}^{C}$ is the in situ concentration of large particulate organic carbon (`biobdet`, [mmol C m<sup>-3</sup>]) <br>
- $B_{CaCO_3}^{C}$ is the in situ concentration of calcium carbonate carbon (`biocaco3`, [mmol C m<sup>-3</sup>]) <br>
- $B_{ld}^{Si}$ is the in situ concentration of biogenic silica (`biobdetsi`, [mmol Si m<sup>-3</sup>]) <br>
- $\dfrac{12}{0.4}$ is the g (mol C)<sup>-1</sup> and assuming that 40% of the total biomass of particulate organics is carbon. <br>
- $100$ is the g (mol C)<sup>-1</sup> of calcium carbonate. <br>
- $60$ is the g (mol Si)<sup>-1</sup> of biogenic silica. <br>

We consider $CaCO_3$ to be part of the small sinking particulates because, although more dense than organic matter, $CaCO_3$ particles tend to be smaller than organic aggregates and sink at a slower rate ([De La Rocha & Passow, 2007](https://doi.org/10.1016/j.dsr2.2007.01.004); [Zhang et al., 2018](https://doi.org/10.5194/bg-15-4759-2018)). Furthermore, the shedding of coccoliths by coccolithophores, which are near-neutrally bouyant, also contributes to a slower mean sinking speed of $CaCO_3$ ([Balch et al., 2009](https://doi.org/10.1029/2008JC004902)). In contrast, we consider biogenic silica to be part of the large sinking particulates:

$$
\begin{align}
M_{s} =& \quad M_{sd} + M_{CaCO_3} \\
M_{l} =& \quad M_{ld} + M_{BSi}
\end{align}
$$

And we compute the harmonic mean density of the small particulates (`rho_small`, $\rho_{s}$, [kg m<sup>-3</sup>]) and large particulates (`rho_large`, $\rho_{l}$, [kg m<sup>-3</sup>]) weighted by mass fractions, which accounts for the fact that less dense mass fractions account for greater volume within aggregates:

$$
\begin{align}
\rho_{s} =& \quad \dfrac{1}{\dfrac{M_{sd} / M_{s}}{\rho_{orgC}} + \dfrac{M_{CaCO_3} / M_{s}}{\rho_{CaCO_3}} } \\
\rho_{l} =& \quad \dfrac{1}{\dfrac{M_{ld} / M_{l}}{\rho_{orgC}} + \dfrac{M_{Bsi} / M_{l}}{\rho_{Bsi}} }
\end{align}
$$

_where_ <br>
- $\rho_{orgC}$ is the density of organic carbon particles (`detrho`, [kg m<sup>-3</sup>]) <br>
- $\rho_{CaCO_3}$ is the density of calcium carbonate particles (`caco3rho`, [kg m<sup>-3</sup>]) <br>
- $\rho_{BSi}$ is the density of biogenic silica particles (`bsirho`, [kg m<sup>-3</sup>]) <br>

Finally, we incorporate estimates of particle porosity to their density:

$$
\begin{align}
\rho_{s} =& \quad \left(1 - p_{s}\right) \rho_{s} + p_{s} \rho_{sw} \\
\rho_{l} =& \quad \left(1 - p_{l}\right) \rho_{l} + p_{l} \rho_{sw}
\end{align}
$$

_where_ <br>
- $p_{s}$ is the porosity of small particles (`detphi`, [dimensionless]) <br>
- $p_{l}$ is the porosity of large particles (`bdetphi`, [dimensionless]) <br>
- $\rho_{sw}$ is the density of seawater, which we set here to a constant 1025 (kg m<sup>-3</sup>) <br>


**Rubey's equation**

Rubey's equation ([Rubey, 1933](https://doi.org/10.2475/ajs.s5-25.148.325)) combines the radius of particles, their excess density relative to fluid and the dynamic viscosity of that fluid to compute the sinking rate of particles. We find the sinking rate of small (`wsink1(k)`, [m s<sup>-1</sup>]) and large particles (`wsink2(k)`, [m s<sup>-1</sup>]) using Rubey's equation:

$$
\begin{align}
\omega_{s} =& \quad \dfrac{ \sqrt{\dfrac{4}{3} 9.8 \cdot \rho_{sw} \left(\rho_{s} - \rho_{sw} \right) \left(r_{s}\right)^{3} + 9 \left(\eta_{sw}\right)^{2}} - 3 \eta_{sw} }{\rho_{sw} \cdot r_{s}}
\end{align}
$$

$$
\begin{align}
\omega_{l} =& \quad \dfrac{ \sqrt{\dfrac{4}{3} 9.8 \cdot \rho_{sw} \left(\rho_{l} - \rho_{sw} \right) \left(r_{l}\right)^{3} + 9 \left(\eta_{sw}\right)^{2}} - 3 \eta_{sw} }{\rho_{sw} \cdot r_{l}}
\end{align}
$$

_where_ <br>
- $r_{s}$ and $r_{l}$ are the mean radii of small and large particles (`rad_det`; `rad_bdet`; [m]) <br>
- $\eta_{sw}$ is the dynamic viscosity of seawater at in situ temperature, salinity and pressure (`dynvis_sw(i,j,k)`, [kg m<sup>-1</sup> s<sup>-1</sup>]) <br>
- $\rho_{s}$ and $\rho_{l}$ are the harmonic mean densities of small and large particles (`rho_small`; `rho_large`, [kg m<sup>-3</sup>]) <br>
- $\rho_{sw}$ is the density of seawater, which we set here to a constant 1025 (kg m<sup>-3</sup>) <br>

Our approach therefore considers mineral ballasting on particle excess density, particle size and the viscosity of fluid in determining sinking rates. This allows for "an environmentally dependent, space-varying $\omega_{s}$ and $\omega_{l}$" ([Dinauer et al., 2022](https://doi.org/10.1029/2021GB007131)).


---


### 22. Sedimentary processes.

WOMBAT-mid tracks the accumulation of organic detrital carbon (`p_det_sediment(i,j)`, $B_{det,sed}^{C}$, [mol C m<sup>-2</sup>]), organic detrital iron (`p_detfe_sediment(i,j)`, $B_{det,sed}^{Fe}$, [mol Fe m<sup>-2</sup>]), organic detrital silica (`p_detsi_sediment(i,j)`, $B_{det,sed}^{Si}$, [mol Si m<sup>-2</sup>]) and $CaCO_3$ (`p_caco3_sediment(i,j)`, $B_{CaCO_3,sed}^{C}$, [mol C m<sup>-2</sup>]) within sedimentary pools. The organic pools contribute to bottom fluxes of dissolved organic carbon (DOC), dissolved organic nitrogen (DON), dissolved inorganic carbon (DIC), dissolved iron (dFe), silicic acid (H<sub>4</sub>SiO<sub>4</sub>), oxygen (O<sub>2</sub>) and alkalinity (Alk). 


**Organics**

Remineralisation of organic carbon ($\gamma_{sed}^{\rightarrow C}$) produces DOC and DON and removes O<sub>2</sub>. Remineralisation of organic iron produces dFe and remineralisation of organic silica produces silicic acid. Ratios of nitrogen to carbon and oxygen to carbon are static at 16:122 and 132:122.

$$
\begin{align}
\gamma_{sed}^{\rightarrow DOC} =& \quad \gamma_{sed}^{0^{\circ}C} (β_{hete})^{T} B_{sed}^{C} \\
\gamma_{sed}^{\rightarrow DON} =& \quad \gamma_{sed}^{\rightarrow DOC} R^{N:C} \\
\gamma_{sed}^{\leftarrow O_2} =& \quad \gamma_{sed}^{\rightarrow DOC} R^{O_2:C} \\
\gamma_{sed}^{\rightarrow dFe} =& \quad \gamma_{sed}^{0^{\circ}C} (β_{hete})^{T} B_{sed}^{Fe}
\end{align}
$$

_where_ <br>
- $\gamma_{sed}^{0^{\circ}C}$ is a base rate of organic matter hydrolysation at 0ºC in the sediments (`detlrem_sed`, [s<sup>-1</sup>]) <br>
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless]) <br>
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC]) <br>
- $B_{sed}^{C}$ is the concentration of organic carbon in the sediment pool (`p_det_sediment(i,j,1)`, [mol C m<sup>-2</sup>]) <br>
- $B_{sed}^{Fe}$ is the concentration of organic iron in the sediment pool (`p_detfe_sediment(i,j,1)`, [mol Fe m<sup>-2</sup>]) <br>
- $R^{N:C}$ is the static Redfield ratio of nitrogen to carbon in the organic matter (`16/122)` [mol N (mol C)<sup>-1</sup>]) <br>
- $R^{O_2:C}$ is the static Redfield ratio of dissolved oxygen to carbon required to hydrolyse organic matter (`132/122)` [mol O<sub>2</sub> (mol C)<sup>-1</sup>]) <br>


**Dissolution of biogenic silica**

With regard to the dissolution of sedimentary biogenic silica, we compute it in the same way as how it is computed in the water column:

$$
\begin{align}
\gamma_{sed}^{\rightarrow Si} =& \quad \left(d_{sed^{Si}}^{T} \cdot S_{sed^{Si}}^{Sat} \cdot S_{sed^{Si}}^{bio}\right) B_{sed}^{Si}
\end{align}
$$

_where_ <br>
- $B_{sed}^{Si}$ is the concentration of organic silicon in the sediment pool (`p_detsi_sediment(i,j,1)`, [mol Si m<sup>-2</sup>]) 
- $d_{sed^{Si}}^{T}$ is the temperature-dependent rate of dissolution (`disssi_temp`, [s<sup>-1</sup>]) <br>
- $S_{sed^{Si}}^{Sat}$ is a scaling factor that decelerates dissolution as the in situ concentration approachs the equilibrium concentration (`disssi_usat`, [dimenionless]) <br>
- $S_{sed^{Si}}^{bio}$ is a scaling factor that accelerates dissolution in the presence of heterotrophic bacterial biomass (`disssi_bact`, [dimenionless]) <br>

Please refer to the description above for the equations that describe these terms.


**Dissolution of $CaCO_3$**

Dissolution of $CaCO_3$ produces DIC and alkalinity. The sedimenary $CaCO_3$ pool is considered as entirely calcite. If `do_caco3_dynamics = .true.`, then sedimentary dissolution is controlled by bottom water temperature and an estimate of the pore-water calcite saturation state ($\Omega_{cal,sed}$):

$$
\begin{align}
D_{CaCO_{3},sed} =& \quad d_{CaCO_{3},sed} \left(β_{hete}\right)^{T} \left(1 - \Omega_{cal,sed}\right)^{4.5}
\end{align}
$$

_where_ <br>
- $d_{CaCO_{3},sed}$ is a base rate of dissolution in units of [day<sup>-1</sup>] <br>
- $β_{hete}$ is the base temperature-sensitivity coefficient for heterotrophy (`bbioh`, [dimenionless]) <br>
- $T$ is the in situ temperature (`Temp(i,j,k)`, [ºC]) <br>
- $\Omega_{cal,sed}$ is the calcite saturation state within sedimentary pore waters (`sedomega_cal(i,j)`, [dimensionless]) <br>

The $\Omega_{cal,sed}$ is calculated using the `mocsy` package for solving carbonate chemistry of seawater ([Orr & Epitalon, 2015](https://doi.org/10.5194/gmd-8-485-2015)). These routines require Alk and DIC as inputs, along with nutrient concentrations and temperature and salinity of bottom waters. For DIC, we chose to sum the water column concentration of DIC and the organic carbon content of the sediment to approximate the interstitial (i.e., porewater) DIC concentration. We assume that the organic carbon content of the sediment (`p_det_sediment`), which is held in units of in [mol m<sup>-2</sup>] is relevant over 10 centimeters, and therefore can be automatically converted to [mol m<sup>-3</sup>] via division by 0.1. With this assumption these arrays can be added together and approximates the reducing conditions of organic-rich sediments, which have lower $\Omega_{cal,sed}$. This ensures a greater rate of $CaCO_3$ dissolution within the sediment as organic matter accumulates.

However, if `do_caco3_dynamics = .false.`, then dissolution of $CaCO_3$ in the sediments proceeds according to a constant assumed $\Omega_{cal,sed}$ of 0.2. We are aware that such a low $\Omega_{cal,sed}$ would not occur in real sediments given the buffering effect of dissolving $CaCO_3$ and the subsequent release of alkalinity. However, in the absence of feedbacks between organic carbon remineralisation and $CaCO_3$ dissolution, we assert a low $\Omega_{cal,sed}$ to ensure that sufficient $CaCO_3$ is dissolved back into the water column.


**Benthic denitrification**

We also consider the consumption of NO<sub>3</sub> via benthic denitrification. When `do_benthic_denitrification = .true.`, a portion of the particulate organic matter within the sediments that is hydrolysed to DOC and DON is performed anaerobically (i.e., using NO<sub>3</sub> as the electron acceptor). Unlike this process in the water column, which is performed by bacterial metabolism, we estimate this process using an empirical parameterization from [Bohlen et al. (2012)](https://doi.org/10.1029/2011GB004198):

$$
\begin{align}
\gamma_{sed}^{\leftarrow NO_3} =& \quad \gamma_{sed}^{\rightarrow DOC} \min\left(0.9 \dfrac{94}{122}, \left(0.083 + 0.21 \cdot 0.98^{O_2 - NO_3} \right) \right)
\end{align}
$$

_where_ <br>
- $0.9$ is a hard upper limit stating that 90% of organic matter hydrolysation can potentially be performed anaerobically via denitrification <br>
- $\dfrac{94}{122}$ is the stoichiometry of nitrate demand per mol of organic carbon hydrolysed ([Paulmier et al., 2009](https://doi.org/10.5194/bg-6-923-2009)) <br>
- O<sub>2</sub> is the bottom water concentration of dissolved oxygen (mmol m<sup>-3</sup>) <br>
- NO<sub>3</sub> is the bottom water concentration of nitrate (mmol m<sup>-3</sup>) <br>

and where the fraction of organic matter that is hydrolysed via denitrification is equal to:

$$
\begin{align}
f_{sed}^{denit} =& \quad \gamma_{sed}^{\leftarrow NO_3} \dfrac{122/94}{\gamma_{sed}^{\rightarrow DOC}}
\end{align}
$$

**Tendencies from sediment processes**

Overall bottom fluxes of tracers are:

$$
\begin{align}
\dfrac{\Delta DOC}{\Delta t} =& \quad \gamma_{det,sed}^{\rightarrow DOC} \\
\dfrac{\Delta DON}{\Delta t} =& \quad \gamma_{det,sed}^{\rightarrow DON} \\
\dfrac{\Delta NO_3}{\Delta t} =& \quad \gamma_{det,sed}^{\leftarrow NO_3} \\
\dfrac{\Delta O_2}{\Delta t} =& \quad \gamma_{det,sed}^{\leftarrow O_2} \left(1 - f_{sed}^{denit}\right) \\
\dfrac{\Delta Si}{\Delta t} =& \quad \gamma_{det,sed}^{\rightarrow Si} \\
\dfrac{\Delta dFe}{\Delta t} =& \quad \gamma_{det,sed}^{\rightarrow dFe} \\
\dfrac{\Delta DIC}{\Delta t} =& \quad D_{CaCO_{3},sed} \\
\dfrac{\Delta Alk}{\Delta t} =& \quad 2 \cdot D_{CaCO_{3},sed} - \dfrac{\Delta NO_3}{\Delta t} \\
\end{align}
$$

---

## Subroutine - "update_from_bottom"

The subroutine `generic_WOMBATmid_update_from_bottom` moves sinking organic material from the water column into the sediment pools.
It is at this point that the model performs permanent burial of sinking organic matter if desired.

---

### Permanent burial of particulates.

If `do_burial = .true.`, we compute the fraction of incident sinking particualte carbon, iron, silicon and $CaCO_3$ that is permanently buried in the sediments. This permanently buried fraction is effectively removed from the model and therefore is not accumulated within the sedimentary pools.

The fraction buried is calculated according to Equation 3 of [Dunne et al. (2007)](https://doi.org/10.1029/2006GB002907):

$$
\begin{align}
F_{bury} =& \quad 0.013 + 0.53 \dfrac{(f_{org})^{2}}{\left(7 + f_{org}\right)^{2}}
\end{align}
$$

where $f_{org}$ is the rain rate of organic carbon detritus on the seafloor in [mmol C m<sup>-2</sup> day<sup>-1</sup>].

As organic matter rains down at a more rapid rate, the fraction of incident organic carbon, organic iron, organic silicon and $CaCO_3$ that is buried increases.

---
