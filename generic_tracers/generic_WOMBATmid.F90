!-----------------------------------------------------------------------
!
!                       PRIMARY DEVELOPERS
!
! <CONTACT EMAIL="Pearse.Buchanan@csiro.au"> Pearse Buchanan
! </CONTACT>
!
! <CONTACT EMAIL="Dougal.Squire@anu.edu.au"> Dougie Squire
! </CONTACT>
!
!
!                        OTHER CONTACTS
!
! <CONTACT EMAIL="Richard.Matear@csiro.au"> Richard Matear
! </CONTACT>
!
! <CONTACT EMAIL="Matthew.Chamberlain@csiro.au"> Matt Chamberlain
! </CONTACT>
!
!
!
! <OVERVIEW>
!  This module contains the generic version of WOMBATmid.
!  It is designed so that both GFDL Ocean models, GOLD and MOM, can use
!  it.
! </OVERVIEW>
!
! <DESCRIPTION>
!
!
!     (\___/)  .-.   .-. .--. .-..-..---.  .--. .-----.
!     / o o \  : :.-.: :: ,. :: `' :: .; :: .; :`-. .-'
!    (   "   ) : :: :: :: :: :: .. ::   .':    :  : :
!     \__ __/  : `' `' ;: :; :: :; :: .; :: :: :  : :
!               `.,`.,' `.__.':_;:_;:___.':_;:_;  :_;
!
!  World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT)
!
!  World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT) is
!  based on a NPZD (nutrient–phytoplankton–zooplankton–detritus) model.
!  This is the "mid" version of WOMBAT which includes two classes each of
!  phytoplankton, zooplankton, particulate organics, dissolved organics,
!  and heterotrophic bacteria. Nutrients are nitrate (NO3), ammonium (NH4),
!  dissolved iron (Fe), and silicic acid (SIL). Dissolved organic matter is
!  split into large and small molecules of carbon (lDOC, sDOC) and nitrogen
!  (lDON, sDON). Other tracers are dissolved inorganic carbon (DIC), calcium
!  carbonate (CaCO3), alkalinity (ALK), and oxygen (O2). Varying Fe quotas
!  exist in phytoplankton, zooplankton and particulate organic matter. Fe is
!  additionally routed to small and large authigenic particle pools (sAFe
!  and lAFe) via the "colloidal shunt". C:N ratios are fixed in all organic
!  pools except for dissolved organics. Si is carried through micro-
!  phytoplankton, large particulate organic matter and, like for carbon and
!  Fe, is deposited into a sediment pool.
!  Gas exchange follows MOCSY protocols.
! </DESCRIPTION>
!
! <INFO>
!  <REFERENCE>
!   This model is available for public use. Note that different tracers
!   may exist in different versions of the module.
!  </REFERENCE>
!
!  <DEVELOPER_NOTES>
!   This code was originally ported from WOMBAT v3 here:
!   https://github.com/mom-ocean/MOM5/tree/d7ba13a3f364ce130b6ad0ba813f01832cada7a2/src/mom5/ocean_csiro_bgc
!   using generic_BLING.F90 as a template.
!  </DEVELOPER_NOTES>
! </INFO>
!
! <NAMELIST NAME="generic_wombatmid_nml">
!  <DATA NAME="co2_calc" TYPE="character">
!   Defines the carbon equiliabration method.  Default is 'ocmip2' which
!   uses the FMS_ocmip2_co2calc routine.  The other option is 'mocsy',
!   which uses the set of routines authored by J. Orr. See reference at:
!   http://ocmip5.ipsl.jussieu.fr/mocsy/index.html
!  </DATA>
!
!  <DATA NAME="do_caco3_dynamics" TYPE="logical">
!   If true, do dynamic CaCO3 precipitation, dissolution and ballasting
!  </DATA>
!
!  <DATA NAME="do_colloidal_shunt" TYPE="logical">
!   If true, do colloidal shunt and coagulation to authigenic pools
!  </DATA>
!
!  <DATA NAME="do_two_ligands" TYPE="logical">
!   If true, do two ligands (one strong, one weak) for iron complexation.
!  </DATA>
!
!  <DATA NAME="do_burial" TYPE="logical">
!   If true, permanently bury organics and CaCO3 in sediments
!  </DATA>
!
!  <DATA NAME="do_nitrogen_fixation" TYPE="logical">
!   If true, do nitrogen fixation
!  </DATA>
!
!  <DATA NAME="do_anammox" TYPE="logical">
!   If true, do anammox
!  </DATA>
!
!  <DATA NAME="do_wc_denitrification" TYPE="logical">
!   If true, do water column denitrification
!  </DATA>
!
!  <DATA NAME="do_benthic_denitrification" TYPE="logical">
!   If true, do benthic denitrification
!  </DATA>
!
!  <DATA NAME="do_tracer_dicp" TYPE="logical">
!   If true, do carry preformed dissolved inorganic carbon (dicp) as a tracer
!  </DATA>
!
!  <DATA NAME="do_tracer_dicr" TYPE="logical">
!   If true, do carry remineralised dissolved inorganic carbon (dicr) as a tracer
!  </DATA>
!
!  <DATA NAME="do_viscous_sinking" TYPE="logical">
!   If true, do computation of dynamic viscosity
!  </DATA>
!
!  <DATA NAME="do_check_n_conserve" TYPE="logical">
!   If true, check that the ecosystem model conserves nitrogen. NOTE:
!   not appropriate if dentirification, anammox and nitrogen fixation are on.
!  </DATA>
!
!  <DATA NAME="do_check_c_conserve" TYPE="logical">
!   If true, check that the ecosystem model conserves carbon
!  </DATA>
!
!  <DATA NAME="do_check_si_conserve" TYPE="logical">
!   If true, check that the ecosystem model conserves silicon
!  </DATA>
!
!  <DATA NAME="do_check_fe_conserve" TYPE="logical">
!   If true, check that the ecosystem model conserves iron
!  </DATA>
!
! </NAMELIST>
!
!-----------------------------------------------------------------------

module generic_WOMBATmid

  use field_manager_mod, only: fm_string_len
  use mpp_mod,           only: input_nml_file, mpp_error, FATAL, WARNING
  use fms_mod,           only: write_version_number, check_nml_error, stdout, stdlog
  use time_manager_mod,  only: time_type
  use constants_mod,     only: WTMCO2, WTMO2

  use g_tracer_utils, only : g_diag_type, g_tracer_type
  use g_tracer_utils, only : g_tracer_start_param_list, g_tracer_end_param_list
  use g_tracer_utils, only : g_tracer_add, g_tracer_add_param, g_tracer_set_files
  use g_tracer_utils, only : g_tracer_get_common, g_tracer_get_pointer
  use g_tracer_utils, only : g_tracer_get_values, g_tracer_set_values
  use g_tracer_utils, only : register_diag_field=>g_register_diag_field
  use g_tracer_utils, only : g_send_data

  use FMS_ocmip2_co2calc_mod, only : FMS_ocmip2_co2calc, CO2_dope_vector

  implicit none ; private

  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

  character(len=fm_string_len), parameter :: mod_name     = 'generic_WOMBATmid'
  character(len=fm_string_len), parameter :: package_name = 'generic_wombatmid'

  public do_generic_WOMBATmid
  public generic_WOMBATmid_register
  public generic_WOMBATmid_init
  public generic_WOMBATmid_register_diag
  public generic_WOMBATmid_update_from_coupler
  public generic_WOMBATmid_update_from_source
  public generic_WOMBATmid_update_from_bottom
  public generic_WOMBATmid_set_boundary_values
  public generic_WOMBATmid_end

  ! The following variable for using this module is overwritten by
  ! generic_tracer_nml namelist
  logical, save :: do_generic_WOMBATmid = .false.

  real, parameter :: missing_value1 = -1.0e+10

  !=======================================================================
  ! Namelist Options
  !=======================================================================
  character(len=10) :: co2_calc  = 'mocsy' ! other option is 'ocmip2'
  logical :: do_caco3_dynamics          = .true.  ! do dynamic CaCO3 precipitation, dissolution and ballasting?
  logical :: do_colloidal_shunt         = .true.  ! do colloidal shunt and coagulation to authigenic pools?
  logical :: do_two_ligands             = .false. ! do two ligands (one strong, one weak) for iron complexation?
  logical :: do_burial                  = .false. ! permanently bury organics and CaCO3 in sediments?
  logical :: do_nitrogen_fixation       = .true.  ! N cycle has nitrogen fixation?
  logical :: do_anammox                 = .true.  ! N cycle has anammox?
  logical :: do_wc_denitrification      = .true.  ! N cycle has water column denitrification?
  logical :: do_benthic_denitrification = .true.  ! N cycle has N loss in sediments?
  logical :: do_tracer_dicp             = .false. ! Enable preformed dissolved inorganic carbon tracer, dicp?
  logical :: do_tracer_dicr             = .false. ! Enable remineralized dissolved inorganic carbon tracer, dicr?
  logical :: do_viscous_sinking         = .true.  ! Rubey's formula uses a non-constant dynamic viscosity?
  logical :: do_check_n_conserve        = .false. ! check that the N fluxes balance in the ecosystem
  logical :: do_check_c_conserve        = .false. ! check that the C fluxes balance in the ecosystem
  logical :: do_check_si_conserve       = .false. ! check that the Si fluxes balance in the ecosystem
  logical :: do_check_fe_conserve       = .false. ! check that the Fe fluxes balance in the ecosystem

  namelist /generic_wombatmid_nml/ co2_calc, do_caco3_dynamics, do_colloidal_shunt, do_two_ligands, do_burial, &
                                   do_nitrogen_fixation, do_anammox, do_wc_denitrification, do_benthic_denitrification, &
                                   do_tracer_dicp, do_tracer_dicr, do_viscous_sinking, &
                                   do_check_n_conserve, do_check_c_conserve, do_check_si_conserve, do_check_fe_conserve

  !=======================================================================
  ! This type contains all the parameters and arrays used in this module
  !=======================================================================
  type generic_WOMBATmid_type
    !-----------------------------------------------------------------------
    ! Configurable parameters
    !-----------------------------------------------------------------------
    ! See user_add_params for descriptions of each parameter
    logical :: &
        init, &
        force_update_fluxes ! Set in generic_tracer_nml

    real :: &
        alphabio_phy, &
        abioa_phy, &
        bbioa_phy, &
        alphabio_dia, &
        abioa_dia, &
        bbioa_dia, &
        alphabio_tri, &
        bbioh, &
        phykn, &
        phykf, &
        phyprefnh4, &
        phyminqc, &
        phymaxqc, &
        phyoptqf, &
        phymaxqf, &
        phylmor, &
        phyqmor, &
        diakn, &
        diakf, &
        diaks, &
        diaprefnh4, &
        diaminqc, &
        diamaxqc, &
        diaoptqf, &
        diamaxqf, &
        diaminqs, &
        diaoptqs, &
        diamaxqs, &
        diaVmaxs, &
        dialmor, &
        diaqmor, &
        chltau, &
        overflow, &
        trikf, &
        trichlc, &
        trin2c, &
        zooCingest, &
        zooCassim, &
        zooFeingest, &
        zooFeassim, &
        zooexcrdom, &
        zoogmax, &
        zooepslbac, &
        zooepssbac, &
        zooepsaoa, &
        zooepsphy, &
        zooepsdia, &
        zooepssdet, &
        zpreflbac, &
        zprefsbac, &
        zprefaoa, &
        zprefphy, &
        zprefdia, &
        zprefsdet, &
        zoolmor, &
        zooqmor, &
        mesCingest, &
        mesCassim, &
        mesFeingest, &
        mesFeassim, &
        mesexcrdom, &
        fgutdiss, &
        mesgmax, &
        mesepslbac, &
        mesepssbac, &
        mesepsaoa, &
        mesepsphy, &
        mesepsdia, &
        mesepssdet, &
        mesepsldet, &
        mesepszoo, &
        mpreflbac, &
        mprefsbac, &
        mprefaoa, &
        mprefphy, &
        mprefdia, &
        mprefsdet, &
        mprefldet, &
        mprefzoo, &
        meslmor, &
        mesqmor, &
        zoopreyswitch, &
        mespreyswitch, &
        sdetlrem, &
        bottom_thickness, &
        detlrem_sed, &
        sdetphi, &
        ldetphi, &
        phyrad0, &
        diarad0, &
        zoorad0, &
        mesrad0, &
        sdetrho, &
        caco3rho, &
        bsirho, &
        phybiot, &
        diabiot, &
        caco3lrem, &
        caco3lrem_sed, &
        f_inorg, &
        disscal, &
        dissara, &
        dissdet, &
        ligW, &
        ligS, &
        dfefloor, &
        kscav_dfe, &
        kcoag_dfe, &
        kagg_col, &
        kagg_kcol, &
        ksafe_dfe, &
        klafe_dfe, &
        wsafe, &
        wlafe, &
        bsi_alpha, &
        bsi_fbac, &
        bsi_kbac, &
        bsilrem_sed, &
        aoa_knh4, &
        aoa_poxy, &
        aoa_ynh4, &
        aoa_yoxy, &
        aoa_C2N, &
        aoa_C2Fe, &
        aoalmor, &
        aoaqmor, &
        lbac_Vmax_doc, &
        lbac_Vmax_dfe, &
        lbac_Vmax_no3, &
        lbac_poxy, &
        lbac_kno3, &
        lbac_kdoc, &
        lbac_kfer, &
        lbac_alpha, &
        lbac_fele, &
        sbac_Vmax_doc, &
        sbac_Vmax_dfe, &
        sbac_Vmax_no3, &
        sbac_poxy, &
        sbac_kno3, &
        sbac_kdoc, &
        sbac_kfer, &
        sbac_fele, &
        bac_C2N, &
        bac_C2Fe, &
        baclmor, &
        bacqmor, &
        aoxkn, &
        aoxmumax, &
        dt_npzd, &
        sal_global, &
        dic_global, &
        alk_global, &
        no3_global, &
        nh4_global, &
        htotal_scale_lo, &
        htotal_scale_hi, &
        htotal_in, &
        Rho_0, &
        a_0, a_1, a_2, a_3, a_4, a_5, &
        b_0, b_1, b_2, b_3, c_0, &
        a1_co2, a2_co2, a3_co2, a4_co2, a5_co2, &
        a1_o2, a2_o2, a3_o2, a4_o2, a5_o2

    character(len=fm_string_len) :: ice_restart_file
    character(len=fm_string_len) :: ocean_restart_file
    character(len=fm_string_len) :: IC_file

    !-----------------------------------------------------------------------
    ! Arrays for surface gas fluxes
    !-----------------------------------------------------------------------
    real, dimension(:,:), allocatable :: &
        htotallo, htotalhi,  &
        co2_csurf, co2_alpha, co2_sc_no, pco2_csurf, &
        o2_csurf, o2_alpha, o2_sc_no, &
        no3_vstf, nh4_vstf, dic_vstf, alk_vstf


    real, dimension(:,:,:), allocatable :: &
        htotal, &
        omega_ara, omega_cal, &
        co3, co2_star

    !-----------------------------------------------------------------------
    ! Arrays for tracer fields and source terms
    !-----------------------------------------------------------------------
    ! The prefixes "f_" refers to a "field", "j" to a volumetric rate, "b_"
    ! to a bottom flux and "p_" to a "pointer".
    real, dimension(:,:), allocatable :: &
        b_dic, &
        b_dicr, &
        b_alk, &
        b_sdoc, &
        b_ldoc, &
        b_ldon, &
        b_no3, &
        b_sil, &
        b_o2, &
        b_fe, &
        npp2d, &
        rpp2d, &
        zsp2d, &
        sdet_btm, &
        sdetfe_btm, &
        ldet_btm, &
        ldetfe_btm, &
        ldetsi_btm, &
        caco3_btm, &
        safe_btm, &
        lafe_btm, &
        det_sed_remin, &
        detfe_sed_remin, &
        detsi_sed_remin, &
        caco3_sed_remin, &
        det_sed_denit, &
        fbury, &
        fdenit, &
        zeuphot, &
        sdet_radius, &
        ldet_radius, &
        seddep, &
        sedmask, &
        sedtemp, &
        sedsalt, &
        sedno3, &
        sednh4, &
        sedsil, &
        sedo2, &
        seddic, &
        sedalk, &
        sedhtotal, &
        sedco3, &
        sedomega_cal

    real, dimension(:,:,:), allocatable :: &
        f_dic, &
        f_dicr, &
        f_alk, &
        f_no3, &
        f_nh4, &
        f_sil, &
        f_phy, &
        f_dia, &
        f_pchl, &
        f_dchl, &
        f_phyfe, &
        f_diafe, &
        f_diasi, &
        f_zoo, &
        f_zoofe, &
        f_mes, &
        f_mesfe, &
        f_sdet, &
        f_sdetfe, &
        f_ldet, &
        f_ldetfe, &
        f_ldetsi, &
        f_ldoc, &
        f_ldon, &
        f_sdoc, &
        f_lbac, &
        f_sbac, &
        f_aoa, &
        f_o2, &
        f_caco3, &
        f_fe, &
        f_safe, &
        f_lafe, &
        radbio, &
        radmid, &
        radmld, &
        npp3d, &
        rpp3d, &
        zsp3d, &
        phy_mumax, &
        phy_mu, &
        pchl_mu, &
        phy_kni, &
        phy_kfe, &
        phy_lpar, &
        phy_lnit, &
        phy_lnh4, &
        phy_lno3, &
        phy_lfer, &
        phy_dfeupt, &
        dia_mumax, &
        dia_mu, &
        dchl_mu, &
        dia_kni, &
        dia_kfe, &
        dia_ksi, &
        dia_lpar, &
        dia_lnit, &
        dia_lnh4, &
        dia_lno3, &
        dia_lfer, &
        dia_lsil, &
        dia_dfeupt, &
        dia_silupt, &
        trimumax, &
        tri_lfer, &
        tri_lpar, &
        sileqc, &
        disssi, &
        bsidiss, &
        feIII, &
        felig, &
        ligK, &
        fecol, &
        fescasafe, &
        fescalafe, &
        fecoag2safe, &
        fecoag2lafe, &
        safediss, &
        lafediss, &
        fesources, &
        fesinks, &
        phy_feupreg, &
        phy_fedoreg, &
        phygrow, &
        phydoc, &
        phymorl, &
        phymorq, &
        dia_feupreg, &
        dia_fedoreg, &
        dia_sidoreg, &
        diagrow, &
        diadoc, &
        diamorl, &
        diamorq, &
        zooeps, &
        zoopreflbac, &
        zooprefsbac, &
        zooprefaoa, &
        zooprefphy, &
        zooprefdia, &
        zooprefsdet, &
        zoograzlbac, &
        zoograzsbac, &
        zoograzaoa, &
        zoograzphy, &
        zoograzdia, &
        zoograzsdet, &
        zoomorl, &
        zoomorq, &
        zooexcrlbac, &
        zooexcrsbac, &
        zooexcraoa, &
        zooexcrphy, &
        zooexcrdia, &
        zooexcrsdet, &
        zooegeslbac, &
        zooegessbac, &
        zooegesaoa, &
        zooegesphy, &
        zooegesdia, &
        zooegessdet, &
        meseps, &
        mespreflbac, &
        mesprefsbac, &
        mesprefaoa, &
        mesprefphy, &
        mesprefdia, &
        mesprefsdet, &
        mesprefldet, &
        mesprefzoo, &
        mesgrazlbac, &
        mesgrazsbac, &
        mesgrazaoa, &
        mesgrazphy, &
        mesgrazdia, &
        mesgrazsdet, &
        mesgrazldet, &
        mesgrazzoo, &
        mesmorl, &
        mesmorq, &
        mesexcrlbac, &
        mesexcrsbac, &
        mesexcraoa, &
        mesexcrphy, &
        mesexcrdia, &
        mesexcrsdet, &
        mesexcrldet, &
        mesexcrzoo, &
        mesegeslbac, &
        mesegessbac, &
        mesegesaoa, &
        mesegesphy, &
        mesegesdia, &
        mesegessdet, &
        mesegesldet, &
        mesegeszoo, &
        reminr, &
        ldetremi, &
        sdetremi, &
        ldocremi, &
        sdocprod, &
        sdocremi, &
        pic2poc, &
        sdet_density, &
        ldet_density, &
        dissratcal, &
        dissratara, &
        dissratpoc, &
        caldiss, &
        aradiss, &
        pocdiss, &
        zoodiss, &
        mesdiss, &
        aoa_loxy, &
        aoa_lnh4, &
        aoa_eno3, &
        aoa_mumax, &
        aoa_mu, &
        aoagrow, &
        aoaresp, &
        aoamorl, &
        aoamorq, &
        lbacydoc, &
        lbacgrow, &
        lbacresp, &
        lbacpnh4, &
        lbacpco2, &
        lbacufer, &
        lbac_mu, &
        lbac_fanaer, &
        lbac_ffelim, &
        lbacmorl, &
        lbacmorq, &
        lbacdeni, &
        sbacydoc, &
        sbacgrow, &
        sbacresp, &
        sbacpnh4, &
        sbacpco2, &
        sbacufer, &
        sbac_mu, &
        sbac_fanaer, &
        sbac_ffelim, &
        sbacmorl, &
        sbacmorq, &
        sbacdeni, &
        aox_lnh4, &
        aox_mu, &
        nitrfix, &
        ammox, &
        anammox, &
        dynvis_sw, &
        zw, &
        zm

    real, dimension(:,:,:,:), pointer :: &
        p_o2

    real, dimension(:,:,:), pointer :: &
        p_det_sediment, &
        p_detfe_sediment, &
        p_detsi_sediment, &
        p_caco3_sediment

    real, dimension(:,:,:), pointer :: &
        p_wsdet, &
        p_wsdetfe, &
        p_wldet, &
        p_wldetfe, &
        p_wldetsi, &
        p_wcaco3, &
        p_wsafe, &
        p_wlafe

    real, dimension(:,:), pointer :: &
        p_no3_stf, &
        p_nh4_stf, &
        p_dic_stf, &
        p_alk_stf

    !-----------------------------------------------------------------------
    ! IDs for diagnostics
    !-----------------------------------------------------------------------
    ! See register_diagnostics for descriptions of each diagnostic
    integer :: &
        id_pco2 = -1, &
        id_htotal = -1, &
        id_omega_ara = -1, &
        id_omega_cal = -1, &
        id_co3 = -1, &
        id_co2_star = -1, &
        id_no3_vstf = -1, &
        id_nh4_vstf = -1, &
        id_dic_vstf = -1, &
        id_dicp_vstf = -1, &
        id_alk_vstf = -1, &
        id_dynvis_sw = -1, &
        id_radbio = -1, &
        id_radmid = -1, &
        id_radmld = -1, &
        id_phy_kni = -1, &
        id_phy_kfe = -1, &
        id_phy_lpar = -1, &
        id_phy_lnit = -1, &
        id_phy_lnh4 = -1, &
        id_phy_lno3 = -1, &
        id_phy_lfer = -1, &
        id_phy_dfeupt = -1, &
        id_dia_kni = -1, &
        id_dia_kfe = -1, &
        id_dia_ksi = -1, &
        id_dia_lpar = -1, &
        id_dia_lnit = -1, &
        id_dia_lnh4 = -1, &
        id_dia_lno3 = -1, &
        id_dia_lfer = -1, &
        id_dia_lsil = -1, &
        id_dia_dfeupt = -1, &
        id_dia_silupt = -1, &
        id_trimumax = -1, &
        id_tri_lfer = -1, &
        id_tri_lpar = -1, &
        id_sileqc = -1, &
        id_disssi = -1, &
        id_bsidiss = -1, &
        id_feIII = -1, &
        id_felig = -1, &
        id_ligK = -1, &
        id_fecol = -1, &
        id_fescasafe = -1, &
        id_fescalafe = -1, &
        id_fecoag2safe = -1, &
        id_fecoag2lafe = -1, &
        id_safediss = -1, &
        id_lafediss = -1, &
        id_fesources = -1, &
        id_fesinks = -1, &
        id_phy_feupreg = -1, &
        id_phy_fedoreg = -1, &
        id_phygrow = -1, &
        id_phydoc = -1, &
        id_phymorl = -1, &
        id_phymorq = -1, &
        id_dia_feupreg = -1, &
        id_dia_fedoreg = -1, &
        id_dia_sidoreg = -1, &
        id_diagrow = -1, &
        id_diadoc = -1, &
        id_diamorl = -1, &
        id_diamorq = -1, &
        id_zooeps = -1, &
        id_zoopreflbac = -1, &
        id_zooprefsbac = -1, &
        id_zooprefaoa = -1, &
        id_zooprefphy = -1, &
        id_zooprefdia = -1, &
        id_zooprefsdet = -1, &
        id_zoograzlbac = -1, &
        id_zoograzsbac = -1, &
        id_zoograzaoa = -1, &
        id_zoograzphy = -1, &
        id_zoograzdia = -1, &
        id_zoograzsdet = -1, &
        id_zoomorl = -1, &
        id_zoomorq = -1, &
        id_zooexcrlbac = -1, &
        id_zooexcrsbac = -1, &
        id_zooexcraoa = -1, &
        id_zooexcrphy = -1, &
        id_zooexcrdia = -1, &
        id_zooexcrsdet = -1, &
        id_zooegeslbac = -1, &
        id_zooegessbac = -1, &
        id_zooegesaoa = -1, &
        id_zooegesphy = -1, &
        id_zooegesdia = -1, &
        id_zooegessdet = -1, &
        id_meseps = -1, &
        id_mespreflbac = -1, &
        id_mesprefsbac = -1, &
        id_mesprefaoa = -1, &
        id_mesprefphy = -1, &
        id_mesprefdia = -1, &
        id_mesprefsdet = -1, &
        id_mesprefldet = -1, &
        id_mesprefzoo = -1, &
        id_mesgrazlbac = -1, &
        id_mesgrazsbac = -1, &
        id_mesgrazaoa = -1, &
        id_mesgrazphy = -1, &
        id_mesgrazdia = -1, &
        id_mesgrazsdet = -1, &
        id_mesgrazldet = -1, &
        id_mesgrazzoo = -1, &
        id_mesmorl = -1, &
        id_mesmorq = -1, &
        id_mesexcrlbac = -1, &
        id_mesexcrsbac = -1, &
        id_mesexcraoa = -1, &
        id_mesexcrphy = -1, &
        id_mesexcrdia = -1, &
        id_mesexcrsdet = -1, &
        id_mesexcrldet = -1, &
        id_mesexcrzoo = -1, &
        id_mesegeslbac = -1, &
        id_mesegessbac = -1, &
        id_mesegesaoa = -1, &
        id_mesegesphy = -1, &
        id_mesegesdia = -1, &
        id_mesegessdet = -1, &
        id_mesegesldet = -1, &
        id_mesegeszoo = -1, &
        id_reminr = -1, &
        id_ldetremi = -1, &
        id_sdetremi = -1, &
        id_ldocremi = -1, &
        id_sdocprod = -1, &
        id_sdocremi = -1, &
        id_pic2poc = -1, &
        id_dissratcal = -1, &
        id_dissratara = -1, &
        id_dissratpoc = -1, &
        id_caldiss = -1, &
        id_aradiss = -1, &
        id_pocdiss = -1, &
        id_zoodiss = -1, &
        id_mesdiss = -1, &
        id_aoa_loxy = -1, &
        id_aoa_lnh4 = -1, &
        id_aoa_eno3 = -1, &
        id_aoa_mumax = -1, &
        id_aoa_mu = -1, &
        id_aoagrow = -1, &
        id_aoaresp = -1, &
        id_aoamorl = -1, &
        id_aoamorq = -1, &
        id_lbacydoc = -1, &
        id_lbacgrow = -1, &
        id_lbacresp = -1, &
        id_lbacpnh4 = -1, &
        id_lbacpco2 = -1, &
        id_lbacufer = -1, &
        id_lbac_mu = -1, &
        id_lbac_fanaer = -1, &
        id_lbac_ffelim = -1, &
        id_lbacmorl = -1, &
        id_lbacmorq = -1, &
        id_lbacdeni = -1, &
        id_sbacydoc = -1, &
        id_sbacgrow = -1, &
        id_sbacresp = -1, &
        id_sbacpnh4 = -1, &
        id_sbacpco2 = -1, &
        id_sbacufer = -1, &
        id_sbac_mu = -1, &
        id_sbac_fanaer = -1, &
        id_sbac_ffelim = -1, &
        id_sbacmorl = -1, &
        id_sbacmorq = -1, &
        id_sbacdeni = -1, &
        id_aox_lnh4 = -1, &
        id_aox_mu = -1, &
        id_nitrfix = -1, &
        id_ammox = -1, &
        id_anammox = -1, &
        id_phy_mumax = -1, &
        id_phy_mu = -1, &
        id_pchl_mu = -1, &
        id_dia_mumax = -1, &
        id_dia_mu = -1, &
        id_dchl_mu = -1, &
        id_npp3d = -1, &
        id_rpp3d = -1, &
        id_zsp3d = -1, &
        id_npp2d = -1, &
        id_rpp2d = -1, &
        id_zsp2d = -1, &
        id_det_sed_remin = -1, &
        id_det_sed_depst = -1, &
        id_det_sed_denit = -1, &
        id_fbury = -1, &
        id_fdenit = -1, &
        id_detfe_sed_remin = -1, &
        id_detfe_sed_depst = -1, &
        id_detsi_sed_remin = -1, &
        id_detsi_sed_depst = -1, &
        id_caco3_sed_remin = -1, &
        id_caco3_sed_depst = -1, &
        id_zeuphot = -1, &
        id_sdet_radius = -1, &
        id_ldet_radius = -1, &
        id_sdet_density = -1, &
        id_ldet_density = -1, &
        id_seddep = -1, &
        id_sedmask = -1, &
        id_sedtemp = -1, &
        id_sedsalt = -1, &
        id_sedno3 = -1, &
        id_sednh4 = -1, &
        id_sedsil = -1, &
        id_sedo2 = -1, &
        id_seddic = -1, &
        id_sedalk = -1, &
        id_sedhtotal = -1, &
        id_sedco3 = -1, &
        id_sedomega_cal = -1

  end type generic_WOMBATmid_type

  type(generic_WOMBATmid_type), save :: wombat

  ! An auxiliary type for storing varible names
  type, public :: vardesc
    character(len=fm_string_len) :: name     ! The variable name in a NetCDF file.
    character(len=fm_string_len) :: longname ! The long name of that variable.
    character(len=1)             :: hor_grid ! The hor. grid:  u, v, h, q, or 1.
    character(len=1)             :: z_grid   ! The vert. grid:  L, i, or 1.
    character(len=1)             :: t_grid   ! The time description: s, a, m, or 1.
    character(len=fm_string_len) :: units    ! The dimensions of the variable.
    character(len=1)             :: mem_size ! The size in memory: d or f.
  end type vardesc

  type(CO2_dope_vector) :: CO2_dope_vec

  contains

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATmid_register">
  !  <OVERVIEW>
  !   Register the generic WOMBATmid module
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   This subroutine reads and checks the WOMBATmid namelist and adds all
  !   WOMBATmid tracers via subroutine user_add_tracers()
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATmid_register(tracer_list, force_update_fluxes)
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATmid_register(tracer_list)
    type(g_tracer_type), pointer, intent(in) :: tracer_list

    integer                                 :: ierr
    integer                                 :: io_status
    integer                                 :: stdoutunit, stdlogunit
    character(len=fm_string_len), parameter :: sub_name = 'generic_WOMBATmid_register'
    character(len=256), parameter           :: error_header = &
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
    character(len=256), parameter           :: warn_header =  &
        '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
    character(len=256), parameter           :: note_header =  &
        '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '

    ! Provide for namelist over-ride
    ! This needs to go before user_add_tracers in order to allow the namelist
    ! settings to switch tracers on and off.
    stdoutunit = stdout(); stdlogunit = stdlog()

    read (input_nml_file, nml=generic_wombatmid_nml, iostat=io_status)
    ierr = check_nml_error(io_status, 'generic_wombatmid_nml')

    write (stdoutunit,'(/)')
    write (stdoutunit, generic_wombatmid_nml)
    write (stdlogunit, generic_wombatmid_nml)

    if (trim(co2_calc) == 'ocmip2') then
      write (stdoutunit,*) trim(note_header), 'Using FMS OCMIP2 CO2 routine'
    else if (trim(co2_calc) == 'mocsy') then
      write (stdoutunit,*) trim(note_header), 'Using Mocsy CO2 routine'
    else
      call mpp_error(FATAL,"Unknown co2_calc option specified in generic_wombatmid_nml")
    endif

    if (do_caco3_dynamics) then
      write (stdoutunit,*) trim(note_header), &
          'Doing dynamic CaCO3 precipitation and dissolution'
    endif

    if (do_colloidal_shunt) then
      write (stdoutunit,*) trim(note_header), &
          'Doing colloidal shunt and coagulation to authigenic pools'
    endif

    if (do_two_ligands) then
      write (stdoutunit,*) trim(note_header), &
          'Using two ligands (one strong, one weak) for iron complexation'
    endif

    if (do_burial) then
      write (stdoutunit,*) trim(note_header), &
          'Permanently burying organics and CaCO3 in sediments'
    endif

    if (do_nitrogen_fixation .or. do_anammox .or. do_wc_denitrification .or. do_benthic_denitrification) then
      write (stdoutunit,*) trim(note_header), &
          'Nitrogen cycle has one or more of the following: nitrogen fixation, anammox, water column denitrification, benthic denitrification'
      if (do_check_n_conserve) then
        call mpp_error(FATAL, "do_check_n_conserve = .true. is going to fail because the N cycle is open")
      endif
    endif

    if (do_nitrogen_fixation) then
      write (stdoutunit,*) trim(note_header), &
          'Doing nitrogen fixation'
    endif

    if (do_anammox) then
      write (stdoutunit,*) trim(note_header), &
          'Doing anammox'
    endif

    if (do_wc_denitrification) then
      write (stdoutunit,*) trim(note_header), &
          'Doing water column denitrification'
    endif

    if (do_benthic_denitrification) then
      write (stdoutunit,*) trim(note_header), &
          'Doing benthic denitrification'
    endif

    if (do_tracer_dicp) then
      write (stdoutunit,*) trim(note_header), &
          'Including preformed dissolved inorganic carbon tracer, dicp'
    endif

    if (do_tracer_dicr) then
      write (stdoutunit,*) trim(note_header), &
          'Including remineralised dissolved inorganic carbon tracer, dicr'
    endif

    if (do_viscous_sinking) then
      write (stdoutunit,*) trim(note_header), &
          'Doing dynamic viscosity calculation for input to sinking scheme'
    endif

    if (do_check_n_conserve) then
      write (stdoutunit,*) trim(note_header), &
          'Checking that the ecosystem model conserves nitrogen'
    endif

    if (do_check_c_conserve) then
      write (stdoutunit,*) trim(note_header), &
          'Checking that the ecosystem model conserves carbon'
    endif

    if (do_check_si_conserve) then
      write (stdoutunit,*) trim(note_header), &
          'Checking that the ecosystem model conserves silicon'
    endif

    if (do_check_fe_conserve) then
      write (stdoutunit,*) trim(note_header), &
          'Checking that the ecosystem model conserves iron'
    endif

    ! Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)

  end subroutine generic_WOMBATmid_register

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATmid_init">
  !  <OVERVIEW>
  !   Initialize the generic WOMBATmid module
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   This subroutine: adds all the WOMBATmid tracers to the list of
  !   generic tracers passed to it via utility subroutine g_tracer_add();
  !   adds all the parameters used by this module via utility subroutine
  !   g_tracer_add_param(); and allocates all work arrays used in the
  !   module.
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATmid_init(tracer_list, force_update_fluxes)
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !
  !  <IN NAME="force_update_fluxes" TYPE="logical">
  !   Flag to force update the fluxes every timestep. This maybe be
  !   necessary in situations where the column_physics (update_from_source)
  !   is not called every timestep such as when MOM6
  !   THERMO_SPANS_COUPLING=True
  !  </IN>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATmid_init(tracer_list, force_update_fluxes)
    type(g_tracer_type), pointer, intent(in) :: tracer_list
    logical, intent(in)                      :: force_update_fluxes

    character(len=fm_string_len), parameter :: sub_name = 'generic_WOMBATmid_init'

    wombat%force_update_fluxes = force_update_fluxes

    call write_version_number( version, tagname )

    ! Specify and initialize all parameters used by this package
    call user_add_params

    ! Allocate all the private work arrays used by this module.
    call user_allocate_arrays

  end subroutine generic_WOMBATmid_init

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATmid_register_diag">
  !  <OVERVIEW>
  !   Register diagnostic fields to be used in this module.
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   Register diagnostic fields to be used in this module. Note that the
  !   tracer fields are automatically registered in user_add_tracers. User
  !   adds only diagnostics for fields that are not a member of
  !   g_tracer_type
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATmid_register_diag(diag_list)
  !  </TEMPLATE>
  !
  !  <IN NAME="g_diag_type" TYPE="type(g_diag_type), pointer">
  !   Pointer to the head of generic diag list. Currently, this is not
  !   actually used.
  !  </IN>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATmid_register_diag(diag_list)
    type(g_diag_type), pointer, intent(in) :: diag_list ! dts: this is not actually used

    type(vardesc)   :: vardesc_temp
    integer         :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, axes(3)
    type(time_type) :: init_time

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, &
        axes=axes, init_time=init_time)

    !=======================================================================
    ! Register all diagnostics in this module
    !=======================================================================
    !
    ! The following vardesc types contain a package of metadata about each
    ! tracer, including, in order, the following elements: name; longname;
    ! horizontal staggering ('h') for collocation with thickness points;
    ! vertical staggering ('L') for a layer variable ; temporal staggering
    ! ('s' for snapshot) ; units; and precision in non-restart output files
    ! ('f' for 32-bit float or 'd' for 64-bit doubles). For most tracers, only
    ! the name, longname and units should be changed.
    !
    ! Niki: The register_diag_field interface needs to be extended to take the
    ! MOM6 axes_grp as argument instead of this integer array axes_grp%handle.
    ! Currently the actual MOM6 diag axes is chosen to be T or Tl based on the
    ! size of the axes argument, 2 or 3. The actual values of these axes
    ! argument are not used, only their size is checked to determine the diag
    ! axes! This is not correct since axesTi and axesTl are both of size 3,
    ! likewise there are many axes of size 2. To accomodate axesTi with the
    ! least amount of code modification we can set and check for an input
    ! array of size 1.

    !=======================================================================
    ! Surface flux diagnostics
    !=======================================================================
    !
    ! dts: other gas exchange diagnostics are available via the "ocean_flux",
    ! "ocean_sfc", "atmos_sfc" diagnostic module_names
    vardesc_temp = vardesc( &
        'pco2', 'Surface aqueous partial pressure of CO2', 'h', '1', 's', 'uatm', 'f')
    wombat%id_pco2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'htotal', 'H+ concentration', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_htotal = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'omega_ara', 'Saturation state of aragonite', 'h', 'L', 's', ' ', 'f')
    wombat%id_omega_ara = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'omega_cal', 'Saturation state of calcite', 'h', 'L', 's', ' ', 'f')
    wombat%id_omega_cal = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'co3', 'Carbonate ion concentration', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_co3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'co2_star', 'CO2* (CO2(g) + H2CO3)) concentration', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_co2_star = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'no3_vstf', 'Virtual flux of nitrate into ocean due to salinity restoring/correction', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_no3_vstf = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'nh4_vstf', 'Virtual flux of ammonium into ocean due to salinity restoring/correction', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_nh4_vstf = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dic_vstf', 'Virtual flux of dissolved inorganic carbon into ocean due to '// &
        'salinity restoring/correction', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_dic_vstf = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    if (do_tracer_dicp) then
      vardesc_temp = vardesc( &
          'dicp_vstf', 'Virtual flux of preformed dissolved inorganic carbon into ocean due to '// &
          'salinity restoring/correction', 'h', '1', 's', 'mol/m^2/s', 'f')
      wombat%id_dicp_vstf = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
          init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    endif

    vardesc_temp = vardesc( &
        'alk_vstf', 'Virtual flux of alkalinity into ocean due to salinity restoring/correction', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_alk_vstf = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    !=======================================================================
    ! Tracer and source term diagnostics
    !=======================================================================

    vardesc_temp = vardesc( &
        'dynvis_sw', 'Seawater dynamic viscosity',  'h', 'L', 's', 'kg/m/s', 'f')
    wombat%id_dynvis_sw = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'radbio', 'Photosynthetically active radiation available for phytoplankton growth', &
        'h', 'L', 's', 'W m-2', 'f')
    wombat%id_radbio = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'radmid', 'Photosynthetically active radiation at centre point of grid cell', &
        'h', 'L', 's', 'W m-2', 'f')
    wombat%id_radmid = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'radmld', 'Photosynthetically active radiation averaged in mixed layer', &
        'h', 'L', 's', 'W m-2', 'f')
    wombat%id_radmld = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'det_sed_remin', 'Rate of remineralisation of detritus in accumulated sediment', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_det_sed_remin = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'det_sed_depst', 'Rate of deposition of detritus to sediment at base of water column', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_det_sed_depst = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'det_sed_denit', 'Rate of denitrification (consuming NO3) in accumulated sediment', &
        'h', '1', 's', 'molNO3/m^2/s', 'f')
    wombat%id_det_sed_denit = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fbury', 'Fraction of deposited detritus permanently buried beneath sediment', &
        'h', '1', 's', '[0-1]', 'f')
    wombat%id_fbury = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fdenit', 'Fraction of detritus remineralised via denitrification', &
        'h', '1', 's', '[0-1]', 'f')
    wombat%id_fdenit = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'detfe_sed_remin', 'Rate of remineralisation of detrital iron in accumulated sediment', &
        'h', '1', 's', 'molFe/m^2/s', 'f')
    wombat%id_detfe_sed_remin = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'detfe_sed_depst', 'Rate of deposition of detrital iron to sediment at base of water column', &
        'h', '1', 's', 'molFe/m^2/s', 'f')
    wombat%id_detfe_sed_depst = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'detsi_sed_remin', 'Rate of remineralisation of detrital silicon in accumulated sediment', &
        'h', '1', 's', 'molSi/m^2/s', 'f')
    wombat%id_detsi_sed_remin = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'detsi_sed_depst', 'Rate of deposition of detrital silicon to sediment at base of water column', &
        'h', '1', 's', 'molSi/m^2/s', 'f')
    wombat%id_detsi_sed_depst = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'caco3_sed_remin', 'Rate of remineralisation of CaCO3 in accumulated sediment', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_caco3_sed_remin = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'caco3_sed_depst', 'Rate of deposition of CaCO3 to sediment at base of water column', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_caco3_sed_depst = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'npp3d', 'Net primary productivity', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_npp3d = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'rpp3d', 'Regenerated phytoplankton production', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_rpp3d = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zsp3d', 'Gross zooplankton production', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_zsp3d = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'npp2d', 'Vertically integrated net primary productivity', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_npp2d = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'rpp2d', 'Vertically integrated regenerated primary productivity', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_rpp2d = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zsp2d', 'Vertically integrated zooplankton secondary productivity', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_zsp2d = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_mumax', 'Maximum growth rate of phytoplankton', 'h', 'L', 's', '/s', 'f')
    wombat%id_phy_mumax = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_mu', 'Realised growth rate of phytoplankton', 'h', 'L', 's', '/s', 'f')
    wombat%id_phy_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'pchl_mu', 'Realised growth rate of phytoplankton chlorophyll', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_pchl_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_lpar', 'Limitation of phytoplankton by light', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_phy_lpar = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_kni', 'Half-saturation coefficient of nitrogen uptake by phytoplankton', 'h', 'L', 's', 'mmol/m3', 'f')
    wombat%id_phy_kni = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_kfe', 'Half-saturation coefficient of iron uptake by phytoplankton', 'h', 'L', 's', 'umol/m3', 'f')
    wombat%id_phy_kfe = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_lnit', 'Limitation of phytoplankton by nitrogen', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_phy_lnit = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_lnh4', 'Limitation of phytoplankton by ammonium', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_phy_lnh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_lno3', 'Limitation of phytoplankton by nitrate', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_phy_lno3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_lfer', 'Limitation of phytoplankton by iron', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_phy_lfer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_dfeupt', 'Uptake of dFe by phytoplankton', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_phy_dfeupt = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_mumax', 'Maximum growth rate of microphytoplankton', 'h', 'L', 's', '/s', 'f')
    wombat%id_dia_mumax = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_mu', 'Realised growth rate of microphytoplankton', 'h', 'L', 's', '/s', 'f')
    wombat%id_dia_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dchl_mu', 'Realised growth rate of microphytoplankton chlorophyll', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_dchl_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_lpar', 'Limitation of microphytoplankton by light', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_dia_lpar = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_kni', 'Half-saturation coefficient of nitrogen uptake by microphytoplankton', 'h', 'L', 's', 'mmol/m3', 'f')
    wombat%id_dia_kni = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_kfe', 'Half-saturation coefficient of iron uptake by microphytoplankton', 'h', 'L', 's', 'umol/m3', 'f')
    wombat%id_dia_kfe = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_ksi', 'Half-saturation coefficient of silicic acid uptake by microphytoplankton', 'h', 'L', 's', 'mmol/m3', 'f')
    wombat%id_dia_ksi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_lnit', 'Limitation of microphytoplankton by nitrogen', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_dia_lnit = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_lnh4', 'Limitation of microphytoplankton by ammonium', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_dia_lnh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_lno3', 'Limitation of microphytoplankton by nitrate', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_dia_lno3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_lfer', 'Limitation of microphytoplankton by iron', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_dia_lfer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_lsil', 'Limitation of microphytoplankton by silicic acid', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_dia_lsil = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_dfeupt', 'Uptake of dFe by microphytoplankton', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_dia_dfeupt = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_silupt', 'Uptake of silicic acid by microphytoplankton', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_dia_silupt = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'tri_lfer', 'Limitation of trichodesmium by iron', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_tri_lfer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'tri_lpar', 'Limitation of trichodesmium by light', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_tri_lpar = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'trimumax', 'Trichodesmium temperature-dependent maximum growth rate', 'h', 'L', 's', '/s', 'f')
    wombat%id_trimumax = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sileqc', 'equilibrium concentration of silicic acid', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_sileqc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'disssi', 'dissolution rate of biogenic silica', 'h', 'L', 's', '/s', 'f')
    wombat%id_disssi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bsidiss', 'dissolution of biogenic silica', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_bsidiss = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'feIII', 'free iron (Fe3+)', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_feIII = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'felig', 'ligand-bound dissolved iron', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_felig = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'ligK', 'Ligand stability constant', 'h', 'L', 's', 'L/mol', 'f')
    wombat%id_ligK = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fecol', 'Colloidal dissolved iron', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_fecol = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fescasafe', 'Scavenging of free Fe onto authigenic particles due to smaller organics', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fescasafe = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fescalafe', 'Scavenging of free Fe onto authigenic particles due to larger organics', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fescalafe = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fecoag2safe', 'Coagulation of colloidal dFe onto small authigenic particles', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fecoag2safe = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fecoag2lafe', 'Coagulation of colloidal dFe onto large authigenic particles', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fecoag2lafe = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'safediss', 'Dissolution of small colloidal authigenic Fe particles', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_safediss = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'lafediss', 'Dissolution of large colloidal authigenic Fe particles', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_lafediss = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fesources', 'Total source of dFe in water column', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fesources = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fesinks', 'Total sink of dFe in water column', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fesinks = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_feupreg', 'Factor up regulation of dFe uptake by phytoplankton', 'h', 'L', 's', ' ', 'f')
    wombat%id_phy_feupreg = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
      'phy_fedoreg', 'Factor down regulation of dFe uptake by phytoplankton', 'h', 'L', 's', ' ', 'f')
    wombat%id_phy_fedoreg = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phygrow', 'Growth of phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_phygrow = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phydoc', 'Overflow exudation of DOC by phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_phydoc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phymorl', 'Linear mortality of phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_phymorl = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phymorq', 'Quadratic (density-dependent) mortality of phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_phymorq = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dia_feupreg', 'Factor up regulation of dFe uptake by microphytoplankton', 'h', 'L', 's', ' ', 'f')
    wombat%id_dia_feupreg = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
      'dia_fedoreg', 'Factor down regulation of dFe uptake by microphytoplankton', 'h', 'L', 's', ' ', 'f')
    wombat%id_dia_fedoreg = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
      'dia_sidoreg', 'Factor down regulation of silicic acid uptake by microphytoplankton', 'h', 'L', 's', ' ', 'f')
    wombat%id_dia_sidoreg = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'diagrow', 'Growth of microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_diagrow = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'diadoc', 'Overflow exudation of DOC by microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_diadoc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'diamorl', 'Linear mortality of microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_diamorl = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'diamorq', 'Quadratic (density-dependent) mortality of microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_diamorq = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooeps', 'Zooplankton community-wide prey capture rate coefficient', 'h', 'L', 's', 'm^6/mmolC^2/s', 'f')
    wombat%id_zooeps = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoopreflbac', 'Grazing dietary fraction of zooplankton on large, sharing bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoopreflbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooprefsbac', 'Grazing dietary fraction of zooplankton on small, selfish bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooprefsbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooprefaoa', 'Grazing dietary fraction of zooplankton on ammonia oxidizing archaea', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooprefaoa = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooprefphy', 'Grazing dietary fraction of zooplankton on phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooprefphy = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooprefdia', 'Grazing dietary fraction of zooplankton on microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooprefdia = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooprefsdet', 'Grazing dietary fraction of zooplankton on small detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooprefsdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoograzlbac', 'Grazing rate of zooplankton on large bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoograzlbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoograzsbac', 'Grazing rate of zooplankton on small bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoograzsbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoograzaoa', 'Grazing rate of zooplankton on ammonia oxidizing archaea', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoograzaoa = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoograzphy', 'Grazing rate of zooplankton on phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoograzphy = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoograzdia', 'Grazing rate of zooplankton on microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoograzdia = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoograzsdet', 'Grazing rate of zooplankton on small detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoograzsdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoomorl', 'Linear mortality of zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoomorl = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoomorq', 'Quadratic (density-dependent) mortality of zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoomorq = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooexcrlbac', 'Excretion rate of zooplankton eating large bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooexcrlbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooexcrsbac', 'Excretion rate of zooplankton eating small bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooexcrsbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooexcraoa', 'Excretion rate of zooplankton eating ammonia oxidizing archaea', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooexcraoa = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooexcrphy', 'Excretion rate of zooplankton eating phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooexcrphy = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooexcrdia', 'Excretion rate of zooplankton eating microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooexcrdia = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooexcrsdet', 'Excretion rate of zooplankton eating small detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooexcrsdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooegeslbac', 'Egestion rate of zooplankton on large bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooegeslbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooegessbac', 'Egestion rate of zooplankton on small bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooegessbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooegesaoa', 'Egestion rate of zooplankton on ammonia oxidizing archaea', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooegesaoa = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooegesphy', 'Egestion rate of zooplankton on phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooegesphy = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooegesdia', 'Egestion rate of zooplankton on microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooegesdia = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooegessdet', 'Egestion rate of zooplankton on small detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooegessdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'meseps', 'mesozooplankton community-wide prey capture rate coefficient', 'h', 'L', 's', 'm^6/mmolC^2/s', 'f')
    wombat%id_meseps = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mespreflbac', 'Grazing dietary fraction of mesozooplankton on large bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mespreflbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesprefsbac', 'Grazing dietary fraction of mesozooplankton on small bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesprefsbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesprefaoa', 'Grazing dietary fraction of mesozooplankton on ammonia oxidizing archaea', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesprefaoa = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesprefphy', 'Grazing dietary fraction of mesozooplankton on phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesprefphy = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesprefdia', 'Grazing dietary fraction of mesozooplankton on microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesprefdia = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesprefsdet', 'Grazing dietary fraction of mesozooplankton on small detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesprefsdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesprefldet', 'Grazing dietary fraction of mesozooplankton on large detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesprefldet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesprefzoo', 'Grazing dietary fraction of mesozooplankton on zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesprefzoo = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesgrazlbac', 'Grazing rate of mesozooplankton on large bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesgrazlbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesgrazsbac', 'Grazing rate of mesozooplankton on small bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesgrazsbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesgrazaoa', 'Grazing rate of mesozooplankton on ammonia oxidizing archaea', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesgrazaoa = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesgrazphy', 'Grazing rate of mesozooplankton on phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesgrazphy = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesgrazdia', 'Grazing rate of mesozooplankton on microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesgrazdia = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesgrazsdet', 'Grazing rate of mesozooplankton on small detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesgrazsdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesgrazldet', 'Grazing rate of mesozooplankton on large detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesgrazldet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesgrazzoo', 'Grazing rate of mesozooplankton on zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesgrazzoo = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesmorl', 'Linear mortality of mesozooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesmorl = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesmorq', 'Quadratic (density-dependent) mortality of mesozooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesmorq = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesexcrlbac', 'Excretion rate of mesozooplankton eating large bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesexcrlbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesexcrsbac', 'Excretion rate of mesozooplankton eating small bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesexcrsbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesexcraoa', 'Excretion rate of mesozooplankton eating ammonia oxidizing archaea', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesexcraoa = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesexcrphy', 'Excretion rate of mesozooplankton eating phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesexcrphy = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesexcrdia', 'Excretion rate of mesozooplankton eating microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesexcrdia = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesexcrsdet', 'Excretion rate of mesozooplankton eating small detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesexcrsdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesexcrldet', 'Excretion rate of mesozooplankton eating large detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesexcrldet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesexcrzoo', 'Excretion rate of mesozooplankton eating zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesexcrzoo = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesegeslbac', 'Egestion rate of mesozooplankton on large bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesegeslbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesegessbac', 'Egestion rate of mesozooplankton on small bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesegessbac = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesegesaoa', 'Egestion rate of mesozooplankton on ammonia oxidizing archaea', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesegesaoa = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesegesphy', 'Egestion rate of mesozooplankton on phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesegesphy = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesegesdia', 'Egestion rate of mesozooplankton on microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesegesdia = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesegessdet', 'Egestion rate of mesozooplankton on small detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesegessdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesegesldet', 'Egestion rate of mesozooplankton on large detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesegesldet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesegeszoo', 'Egestion rate of mesozooplankton on zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesegeszoo = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'reminr', 'Rate of remineralisation', 'h', 'L', 's', '/s', 'f')
    wombat%id_reminr = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'ldocremi', 'Remineralisation of large dissolved organic carbon by bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_ldocremi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sdocremi', 'Remineralisation of small dissolved organic carbon by bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_sdocremi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sdocprod', 'Production of small dissolved organic carbon by bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_sdocprod = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sdetremi', 'Hydrolysation of small detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_sdetremi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'ldetremi', 'Hydrolysation of large detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_ldetremi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'pic2poc', 'Inorganic (CaCO3) to organic carbon ratio', 'h', 'L', 's', ' ', 'f')
    wombat%id_pic2poc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dissratcal', 'Dissolution rate of Calcite CaCO3', 'h', 'L', 's', '/s', 'f')
    wombat%id_dissratcal = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dissratara', 'Dissolution rate of Aragonite CaCO3', 'h', 'L', 's', '/s', 'f')
    wombat%id_dissratara = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dissratpoc', 'Dissolution rate of CaCO3 due to POC remin', 'h', 'L', 's', '/s', 'f')
    wombat%id_dissratpoc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoodiss', 'Dissolution of CaCO3 due to microzooplankton grazing', 'h', 'L', 's', 'molCaCO3/kg/s', 'f')
    wombat%id_zoodiss = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesdiss', 'Dissolution of CaCO3 due to mesozooplankton grazing', 'h', 'L', 's', 'molCaCO3/kg/s', 'f')
    wombat%id_mesdiss = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'caldiss', 'Dissolution of Calcite CaCO3', 'h', 'L', 's', 'molCaCO3/kg/s', 'f')
    wombat%id_caldiss = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'aradiss', 'Dissolution of Aragonite CaCO3', 'h', 'L', 's', 'molCaCO3/kg/s', 'f')
    wombat%id_aradiss = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'pocdiss', 'Dissolution of CaCO3 due to POC remin', 'h', 'L', 's', 'molCaCO3/kg/s', 'f')
    wombat%id_pocdiss = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'aoa_loxy', 'Limitation of Ammonia Oxidizing Archaea by oxygen', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_aoa_loxy = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'aoa_lnh4', 'Limitation of Ammonia Oxidizing Archaea by ammonium', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_aoa_lnh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'aoa_eno3', 'Excretion of NO3 produced by Ammonia Oxidizing Archaea during oxidation', 'h', 'L', 's', &
        'mol N / mol Biomass', 'f')
    wombat%id_aoa_eno3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'aoa_mumax', 'Maximum growth rate of Ammonia Oxidizing Archaea', 'h', 'L', 's', '/s', 'f')
    wombat%id_aoa_mumax = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'aoa_mu', 'Realized growth rate of Ammonia Oxidizing Archaea', 'h', 'L', 's', '/s', 'f')
    wombat%id_aoa_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'aoagrow', 'Growth of Ammonia Oxidizing Archaea', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_aoagrow = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'aoaresp', 'Oxygen consumption of Ammonia Oxidizing Archaea', 'h', 'L', 's', 'molO2/kg/s', 'f')
    wombat%id_aoaresp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'aoamorl', 'Linear mortality of Ammonia Oxidizing Archaea', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_aoamorl = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'aoamorq', 'Quadratic mortality of Ammonia Oxidizing Archaea', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_aoamorq = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'lbacydoc', 'Biomass yield of large bacteria (mol DOC per mol biomass grown)', 'h', 'L', 's', 'molDOC/molB', 'f')
    wombat%id_lbacydoc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'lbacgrow', 'Growth of large bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_lbacgrow = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'lbacresp', 'Oxygen consumption by large bacteria', 'h', 'L', 's', 'molO2/kg/s', 'f')
    wombat%id_lbacresp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'lbacpnh4', 'Production of NH4 by large bacteria', 'h', 'L', 's', 'molNH4/kg/s', 'f')
    wombat%id_lbacpnh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'lbacpco2', 'Production of CO2 by large bacteria', 'h', 'L', 's', 'molCO2/kg/s', 'f')
    wombat%id_lbacpco2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'lbacufer', 'Uptake of dFe of large bacteria', 'h', 'L', 's', 'moldFe/kg/s', 'f')
    wombat%id_lbacufer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'lbac_mu', 'Realized growth rate of large bacteria', 'h', 'L', 's', '/s', 'f')
    wombat%id_lbac_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'lbac_fanaer', 'Fraction of growth supported by anaerobic metabolism', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_lbac_fanaer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'lbac_ffelim', 'Bacteria growth limited by iron?', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_lbac_ffelim = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'lbacmorl', 'Linear mortality of large bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_lbacmorl = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'lbacmorq', 'Quadratic mortality of large bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_lbacmorq = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'lbacdeni', 'bacterial denitrification rate (NO3 consumption)', 'h', 'L', 's', '[molN/kg/s]', 'f')
    wombat%id_lbacdeni = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sbacydoc', 'Biomass yield of small bacteria (mol DOC per mol biomass grown)', 'h', 'L', 's', 'molDOC/molB', 'f')
    wombat%id_sbacydoc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sbacgrow', 'Growth of small bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_sbacgrow = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sbacresp', 'Oxygen consumption by small bacteria', 'h', 'L', 's', 'molO2/kg/s', 'f')
    wombat%id_sbacresp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sbacpnh4', 'Production of NH4 by small bacteria', 'h', 'L', 's', 'molNH4/kg/s', 'f')
    wombat%id_sbacpnh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sbacpco2', 'Production of CO2 by small bacteria', 'h', 'L', 's', 'molCO2/kg/s', 'f')
    wombat%id_sbacpco2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sbacufer', 'Uptake of dFe of small bacteria', 'h', 'L', 's', 'moldFe/kg/s', 'f')
    wombat%id_sbacufer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sbac_mu', 'Realized growth rate of small bacteria', 'h', 'L', 's', '/s', 'f')
    wombat%id_sbac_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sbac_fanaer', 'Fraction of growth supported by anaerobic metabolism', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_sbac_fanaer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sbac_ffelim', 'Bacteria growth limited by iron?', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_sbac_ffelim = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sbacmorl', 'Linear mortality of small bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_sbacmorl = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sbacmorq', 'Quadratic mortality of small bacteria', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_sbacmorq = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sbacdeni', 'bacterial denitrification rate (NO3 consumption)', 'h', 'L', 's', '[molN/kg/s]', 'f')
    wombat%id_sbacdeni = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'aox_lnh4', 'Limitation of Anammox bacteria by ammonium', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_aox_lnh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'aox_mu', 'Realized growth rate of Anammox bacteria', 'h', 'L', 's', '/s', 'f')
    wombat%id_aox_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'nitrfix', 'Nitrogen fixation rate (NH4 production)', 'h', 'L', 's', '[mol/kg/s]', 'f')
    wombat%id_nitrfix = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

        vardesc_temp = vardesc( &
        'ammox', 'Ammonia Oxidation rate (NH4 consumption)', 'h', 'L', 's', '[mol/kg/s]', 'f')
    wombat%id_ammox = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'anammox', 'Anammox rate (NH4 consumption)', 'h', 'L', 's', '[mol/kg/s]', 'f')
    wombat%id_anammox = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sdet_density', 'Mean density of small detrital particles', 'h', 'L', 's', 'kg/m3', 'f')
    wombat%id_sdet_density = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'ldet_density', 'Mean density of large detrital particles', 'h', 'L', 's', 'kg/m3', 'f')
    wombat%id_ldet_density = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zeuphot', 'Depth of the euphotic zone (%1 incident light)', 'h', '1', 's', 'm', 'f')
    wombat%id_zeuphot = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sdet_radius', 'Mean radius of small detrital particles', 'h', '1', 's', 'm', 'f')
    wombat%id_sdet_radius = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'ldet_radius', 'Mean radius of large detrital particles', 'h', '1', 's', 'm', 'f')
    wombat%id_ldet_radius = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'seddep', 'Depth of the bottom layer', 'h', '1', 's', 'm', 'f')
    wombat%id_seddep = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedmask', 'Mask of active sediment points', 'h', '1', 's', ' ', 'f')
    wombat%id_sedmask = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedtemp', 'Temperature in the bottom layer', 'h', '1', 's', 'deg C', 'f')
    wombat%id_sedtemp = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedsalt', 'Salinity in the bottom layer', 'h', '1', 's', 'psu', 'f')
    wombat%id_sedsalt = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedno3', 'Nitrate concentration in the bottom layer', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_sedno3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sednh4', 'Ammonium concentration in the bottom layer', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_sednh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedsil', 'Silicic acid  in the bottom layer', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_sedsil = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedo2', 'Oxygen concentration in the bottom layer', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_sedo2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'seddic', 'Dissolved inorganic carbon concentration in the bottom layer', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_seddic = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedalk', 'Alkalinity concentration in the bottom layer', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_sedalk = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedhtotal', 'H+ ion concentration in the bottom layer', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_sedhtotal = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedco3', 'CO3 ion concentration in the bottom layer', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_sedco3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedomega_cal', 'Calcite saturation state in the bottom layer', 'h', '1', 's', ' ', 'f')
    wombat%id_sedomega_cal = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

  end subroutine generic_WOMBATmid_register_diag

  !#######################################################################
  !
  ! This is an internal sub, not a public interface.
  ! Add all the parameters to be used in this module.
  !
  subroutine user_add_params

    !=======================================================================
    ! Specify all parameters used in this modules.
    !=======================================================================
    !
    ! Add the known experimental parameters used for calculations in this
    ! module. All the g_tracer_add_param calls must happen between
    ! g_tracer_start_param_list and g_tracer_end_param_list calls. This
    ! implementation enables runtime overwrite via field_table.
    ! dts: Note, some parameters are required by the user_add_tracers routine
    ! which is run _before_ this one. Those parameters are added in
    ! user_add_tracers.

    ! User adds one call for each parameter below with the template
    ! g_tracer_add_param(name, variable,  default_value)
    call g_tracer_start_param_list(package_name)

    !=======================================================================
    ! General parameters
    !=======================================================================
    !
    ! dts: This was copied from BLING to enable calculation of surface flux
    ! terms when the update_from_source routine is commented out for
    ! debugging.
    call g_tracer_add_param('init', wombat%init, .false. )

    ! Average density of sea water [kg/m^3]
    !-----------------------------------------------------------------------
    ! Rho_0 is used in the Boussinesq approximation to calculations of
    ! pressure and pressure gradients, in units of kg m-3.
    call g_tracer_add_param('Rho_0', wombat%Rho_0, 1035.0)

    !=======================================================================
    ! Surface gas flux parameters
    !=======================================================================

    ! Coefficients for O2 saturation [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('a_0', wombat%a_0, 2.00907)
    call g_tracer_add_param('a_1', wombat%a_1, 3.22014)
    call g_tracer_add_param('a_2', wombat%a_2, 4.05010)
    call g_tracer_add_param('a_3', wombat%a_3, 4.94457)
    call g_tracer_add_param('a_4', wombat%a_4, -2.56847e-01)
    call g_tracer_add_param('a_5', wombat%a_5, 3.88767)
    call g_tracer_add_param('b_0', wombat%b_0, -6.24523e-03)
    call g_tracer_add_param('b_1', wombat%b_1, -7.37614e-03)
    call g_tracer_add_param('b_2', wombat%b_2, -1.03410e-02)
    call g_tracer_add_param('b_3', wombat%b_3, -8.17083e-03)
    call g_tracer_add_param('c_0', wombat%c_0, -4.88682e-07)

    ! Schmidt number coefficients [1]
    !-----------------------------------------------------------------------
    ! Compute the Schmidt number of CO2 in seawater using the
    ! formulation proposed by Wanninkhof (2014, Limnology and Oceanography: Methods, 12, 351-362)
    call g_tracer_add_param('a1_co2', wombat%a1_co2,  2116.8)
    call g_tracer_add_param('a2_co2', wombat%a2_co2, -136.25)
    call g_tracer_add_param('a3_co2', wombat%a3_co2,  4.7353)
    call g_tracer_add_param('a4_co2', wombat%a4_co2, -0.092307)
    call g_tracer_add_param('a5_co2', wombat%a5_co2,  0.0007555)

    ! Compute the Schmidt number of O2 in seawater using the
    ! formulation proposed by Wanninkhof (2014, Limnology and Oceanography: Methods, 12, 351-362)
    call g_tracer_add_param('a1_o2', wombat%a1_o2, 1920.4)
    call g_tracer_add_param('a2_o2', wombat%a2_o2, -135.6)
    call g_tracer_add_param('a3_o2', wombat%a3_o2, 5.2122)
    call g_tracer_add_param('a4_o2', wombat%a4_o2, -0.10939)
    call g_tracer_add_param('a5_o2', wombat%a5_o2, 0.00093777)

    ! Initial H+ concentration [mol/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('htotal_in', wombat%htotal_in, 1.e-8) ! dts: default conc from COBALT

    ! Scale factor to set lower limit of htotal range [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('htotal_scale_lo', wombat%htotal_scale_lo, 0.1)

    ! Scale factor to set upper limit of htotal range [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('htotal_scale_hi', wombat%htotal_scale_hi, 100.0)

    !=======================================================================
    ! NPZD parameters
    !=======================================================================
    ! dts: note the parameter units and default values are as used by WOMBAT
    ! v3 in ACCESS-OM2 and ACCESS-ESM1.5. Unit conversions are done
    ! internally to account for the different units carried in this generic
    ! version of WOMBATmid.

    ! Initial slope of P-I curve for phytoplankton [mol Chl (mol C)-1 (W m-2)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('alphabio_phy', wombat%alphabio_phy, 1.5)

    ! Autotrophy maximum growth rate parameter a for phytoplankton [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('abioa_phy', wombat%abioa_phy, 0.7/86400.0)

    ! Autotrophy maximum growth rate parameter b for phytoplankton [dimensionless]
    ! Q10 = b^(10)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bbioa_phy', wombat%bbioa_phy, 1.055)

    ! Initial slope of P-I curve for microphytoplankton [mol Chl (mol C)-1 (W m-2)-1]
    ! "When diatoms are compared to other groups of phytoplankton, they tend to differ
    !  in their light-related traits and appear to have high α and high μmax (Edwards et al., 2015)."
    !  [quote from Litchman (2022): Trait-based diatom ecology, in "The Molecular Life of Diatoms"]
    !  [Edwards et al. 2015 Limnol Oceanogr 60:540–552]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('alphabio_dia', wombat%alphabio_dia, 2.5)

    ! Initial slope of P-I curve for trichodesmium [mol Chl (mol C)-1 (W m-2)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('alphabio_tri', wombat%alphabio_tri, 1.8)

    ! Autotrophy maximum growth rate parameter a for microphytoplankton [s-1]
    ! [Anderson et al., 2021 Nat Communications]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('abioa_dia', wombat%abioa_dia, 1.0/86400.0)

    ! Autotrophy maximum growth rate parameter b for microphytoplankton [dimensionless]
    ! Q10 = b^(10)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bbioa_dia', wombat%bbioa_dia, 1.070)

    ! Heterotrophy maximum growth rate parameter b [dimensionless]
    ! Q10 = b^(10)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bbioh', wombat%bbioh, 1.072)

    ! Phytoplankton half saturation constant for nitrogen uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phykn', wombat%phykn, 1.0)

    ! Phytoplankton half saturation constant for iron uptake [umolFe/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phykf', wombat%phykf, 1.0)

    ! Nano-phytoplankton preference for ammonium over nitrate [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phyprefnh4', wombat%phyprefnh4, 5.0)

    ! Phytoplankton minimum quota of chlorophyll to carbon [mol Chl (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phyminqc', wombat%phyminqc, 0.008)

    ! Phytoplankton maximum quota of chlorophyll to carbon [mol Chl (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phymaxqc', wombat%phymaxqc, 0.065)

    ! Phytoplankton optimal quota of iron to carbon [mol Fe (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phyoptqf', wombat%phyoptqf, 10e-6)

    ! Phytoplankton maximum quota of iron to carbon [mol Fe (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phymaxqf', wombat%phymaxqf, 50e-6)

    ! Phytoplankton linear mortality rate constant [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phylmor', wombat%phylmor, 0.001/86400.0)

    ! Phytoplankton quadratic mortality rate constant [(mmol C m-3)-1 s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phyqmor', wombat%phyqmor, 0.05/86400.0)

    ! microphytoplankton half saturation constant for nitrogen uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diakn', wombat%diakn, 2.4)

    ! microphytoplankton half saturation constant for iron uptake [umolFe/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diakf', wombat%diakf, 2.7)

    ! microphytoplankton half saturation constant for silicic acid uptake [mmolSi/m3]
    !-----------------------------------------------------------------------
    !  - minimal Ks in natural assemblages of 0.5 - 0.9 mmolSi/m3
    !    [Nelson & Brzezinski, 1990, Marine Ecology Progress Series, 62, 283-292]
    !  - We set 5.0 here as default due to recalculation of variable Ksi below
    call g_tracer_add_param('diaks', wombat%diaks, 6.7)

    ! Micro-phytoplankton preference for ammonium over nitrate [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diaprefnh4', wombat%diaprefnh4, 5.0)

    ! microphytoplankton minimum quota of chlorophyll to carbon [mol Chl (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diaminqc', wombat%diaminqc, 0.004)

    ! microphytoplankton maximum quota of chlorophyll to carbon [mol Chl (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diamaxqc', wombat%diamaxqc, 0.060)

    ! microphytoplankton optimal quota of iron to carbon [mol Fe (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diaoptqf', wombat%diaoptqf, 10e-6)

    ! microphytoplankton maximum quota of iron to carbon [mol Fe (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diamaxqf', wombat%diamaxqf, 65e-6)

    ! microphytoplankton minimal quota of silicon to carbon to build a valve [mol Si (mol C)-1]
    !   Brzezinksi (1985) J. Phycology, 21, 347-357
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diaminqs', wombat%diaminqs, 0.04)

    ! microphytoplankton optimal quota of silicon to carbon [mol Si (mol C)-1]
    !   Brzezinksi (1985) J. Phycology, 21, 347-357
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diaoptqs', wombat%diaoptqs, 0.13)

    ! microphytoplankton maximum quota of silicon to carbon [mol Si (mol C)-1]
    !  Brzezinski et al., (2003) Deep-Sea Research II, 50, 619-633
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diamaxqs', wombat%diamaxqs, 0.60)

    ! microphytoplankton maximum uptake rate of silicon to carbon [mol Si (mol C)-1 s-1]
    !  1.2 - 950 fmol Si cell h-1 [Kolbe & Brunner 2022, in "Molecular Life of Diatoms"]
    !  We take a constrained range of 10-100 fmol Si cell h-1
    !  for a 100 pg C cell (8.3 pmol), this maps to roughly 0.0012 - 0.012 mol Si/mol C h-1
    !  which is roughly 0.03 to 0.3 mol Si/mol C day-1
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diaVmaxs', wombat%diaVmaxs, 0.1/86400.0)

    ! microphytoplankton linear mortality rate constant [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('dialmor', wombat%dialmor, 0.001/86400.0)

    ! microphytoplankton quadratic mortality rate constant [(mmol C m-3)-1 s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diaqmor', wombat%diaqmor, 0.05/86400.0)

    ! Timescale of chlorophyll synthesis by phytoplankton [s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('chltau', wombat%chltau, 86400.0)

    ! Maximum fraction of NPP that can be routed to DOC exudation by phytoplankton [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('overflow', wombat%overflow, 0.50)

    ! Trichodesmium half saturation constant for iron uptake [µmolFe/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('trikf', wombat%trikf, 0.125)

    ! Trichodesmium typical chlorophyll to carbon ratio [mol Chl (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('trichlc', wombat%trichlc, 0.01)

    ! Trichodesmium typical nitrogen to carbon ratio [mol N (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('trin2c', wombat%trin2c, 50.0/300.0)

    ! Zooplankton carbon bulk ingestion efficiency [mol C (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooCingest', wombat%zooCingest, 0.70)

    ! Zooplankton carbon assimilation efficiency [mol C (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooCassim', wombat%zooCassim, 0.40)

    ! Zooplankton iron bulk ingestion efficiency [mol Fe (mol Fe)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooFeingest', wombat%zooFeingest, 0.06)

    ! Zooplankton iron assimilation efficiency [mol Fe (mol Fe)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooFeassim', wombat%zooFeassim, 0.60)

    ! Zooplankton fraction of excretion to DOM [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooexcrdom', wombat%zooexcrdom, 0.70)

    ! Zooplankton maximum grazing rate constant [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zoogmax', wombat%zoogmax, 3.3/86400.0)

    ! Zooplankton prey capture rate constant for large bacteria [(mmol C m-3)-2 s-1]
    !  - e.g., protozoans feeding on large bacteria
    !  - aim for half-saturation coefficent B1/2 = 5.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepslbac', wombat%zooepslbac, 0.10/86400.0)

    ! Zooplankton prey capture rate constant for small bacteria [(mmol C m-3)-2 s-1]
    !  - e.g., protozoans feeding on large bacteria
    !  - aim for half-saturation coefficent B1/2 = 5.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepssbac', wombat%zooepssbac, 0.10/86400.0)

    ! Zooplankton prey capture rate constant for ammonia oxidizing archaea [(mmol C m-3)-2 s-1]
    !  - e.g., ciliates feeding on ammonia oxidizing archaea (similar size or larger than pico)
    !  - aim for half-saturation coefficent B1/2 = 5.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepsaoa', wombat%zooepsaoa, 0.25/86400.0)

    ! Zooplankton prey capture rate constant for nanophytoplankton [(mmol C m-3)-2 s-1]
    !  - e.g., ciliates feeding on small (nano/pico)phytoplankton
    !  - aim for half-saturation coefficent B1/2 = 2.5 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepsphy', wombat%zooepsphy, 0.40/86400.0)

    ! Zooplankton prey capture rate constant for microphytoplankton [(mmol C m-3)-2 s-1]
    !  - e.g., larger ciliates feeding on smaller diatoms and other microphytoplankton
    !  - aim for half-saturation coefficent B1/2 = 3.5 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepsdia', wombat%zooepsdia, 0.40/86400.0)

    ! Zooplankton prey capture rate constant for small detritus [(mmol C m-3)-2 s-1]
    !  - e.g., protozoa grazing on slowly sinking detrital particles
    !  - aim for half-saturation coefficent B1/2 = 5.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepssdet', wombat%zooepssdet, 0.25/86400.0)

    ! Zooplankton preference for large bacteria [dimensionless]
    ! Landry (2025) J. Plankton Res. --> find that ~100 mg C m-2 day-1 of ~500 mg C m-2 d-1
    !  of microzooplankton grazing/biomass gain comes from large bacterivory
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zpreflbac', wombat%zpreflbac, 0.50)

    ! Zooplankton preference for small bacteria [dimensionless]
    ! Landry (2025) J. Plankton Res. --> find that ~100 mg C m-2 day-1 of ~500 mg C m-2 d-1
    !  of microzooplankton grazing/biomass gain comes from large bacterivory
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefsbac', wombat%zprefsbac, 0.50)

    ! Zooplankton preference for ammonia oxidizing archaea [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefaoa', wombat%zprefaoa, 0.50)

    ! Zooplankton preference for nano-phytoplankton [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefphy', wombat%zprefphy, 1.0)

    ! Zooplankton preference for microphytoplankton [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefdia', wombat%zprefdia, 0.25)

    ! Zooplankton preference for small detritus [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefsdet', wombat%zprefsdet, 1.0)

    ! Zooplankton respiration rate constant [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zoolmor', wombat%zoolmor, 0.002/86400.0)

    ! Zooplankton quadratic mortality rate constant [(mmol C m-3)-1 s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooqmor', wombat%zooqmor, 0.05/86400.0)

    ! Mesozooplankton carbon bulk ingestion efficiency [mol C (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesCingest', wombat%mesCingest, 0.75)

    ! Mesozooplankton carbon assimilation efficiency [mol C (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesCassim', wombat%mesCassim, 0.50)

    ! Mesozooplankton iron bulk ingestion efficiency [mol Fe (mol Fe)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesFeingest', wombat%mesFeingest, 0.43)

    ! Mesozooplankton iron assimilation efficiency [mol Fe (mol Fe)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesFeassim', wombat%mesFeassim, 0.75)

    ! Mesozooplankton fraction of excretion to DOM [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesexcrdom', wombat%mesexcrdom, 0.35)

    ! Zooplankton dissolution efficiency of CaCO3 within guts [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('fgutdiss', wombat%fgutdiss, 0.80)

    ! Mesozooplankton maximum grazing rate constant [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesgmax', wombat%mesgmax, 0.30/86400.0)

    ! Mesozooplankton prey capture rate constant for large bacteria [(mmol C m-3)-2 s-1]
    !  - e.g., appendicularians filter feeding on large bacteria
    !  - aim for half-saturation coefficent B1/2 = 3.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepslbac', wombat%mesepslbac, 0.11/86400.0)

    ! Mesozooplankton prey capture rate constant for small bacteria [(mmol C m-3)-2 s-1]
    !  - e.g., appendicularians filter feeding on small bacteria
    !  - aim for half-saturation coefficent B1/2 = 5.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepssbac', wombat%mesepssbac, 0.10/86400.0)

    ! Mesozooplankton prey capture rate constant for ammonia oxidizing archaea [(mmol C m-3)-2 s-1]
    !  - e.g., appendicularians filter feeding on ammonia oxidizing archaea
    !  - aim for half-saturation coefficent B1/2 = 3.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsaoa', wombat%mesepsaoa, 0.11/86400.0)

    ! Mesozooplankton prey capture rate constant for nanophytoplankton [(mmol C m-3)-2 s-1]
    !  - e.g., appendicularians filter feeding on small phytoplankton
    !  - aim for half-saturation coefficent B1/2 = 3.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsphy', wombat%mesepsphy, 0.11/86400.0)

    ! Mesozooplankton prey capture rate constant for microphytoplankton [(mmol C m-3)-2 s-1]
    !  - e.g., copepods preying on diatoms and other microphytoplankton
    !  - aim for half-saturation coefficent B1/2 = 2.5 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsdia', wombat%mesepsdia, 0.20/86400.0)

    ! Mesozooplankton prey capture rate constant for small detritus [(mmol C m-3)-2 s-1]
    !  - e.g., appendicularians filter feeding on fine detritus
    !  - aim for half-saturation coefficent B1/2 = 3.5 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepssdet', wombat%mesepssdet, 0.05/86400.0)

    ! Mesozooplankton prey capture rate constant for large detritus [(mmol C m-3)-2 s-1]
    !  - e.g., copepods consuming sinking aggregates of marine snow
    !  - aim for half-saturation coefficent B1/2 = 10.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsldet', wombat%mesepsldet, 0.10/86400.0)

    ! Mesozooplankton prey capture rate constant for microzooplankton [(mmol C m-3)-2 s-1]
    !  - e.g., chaetognaths preying on copepods; copepods consuming ciliates
    !  - aim for half-saturation coefficent B1/2 = 5.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepszoo', wombat%mesepszoo, 0.10/86400.0)

    ! Mesozooplankton preference for large bacteria [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mpreflbac', wombat%mpreflbac, 0.0)

    ! Mesozooplankton preference for small bacteria [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefsbac', wombat%mprefsbac, 0.0)

    ! Mesozooplankton preference for ammonia oxidizing archaea [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefaoa', wombat%mprefaoa, 0.0)

    ! Mesozooplankton preference for phytoplankton [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefphy', wombat%mprefphy, 0.5)

    ! Mesozooplankton preference for microphytoplankton [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefdia', wombat%mprefdia, 1.0)

    ! Mesozooplankton preference for small detritus [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefsdet', wombat%mprefsdet, 1.0)

    ! Mesozooplankton preference for large detritus (aggregates) [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefldet', wombat%mprefldet, 1.0)

    ! Mesozooplankton preference for micro-zooplankton [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefzoo', wombat%mprefzoo, 1.0)

    ! Mesozooplankton respiration rate constant [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('meslmor', wombat%meslmor, 0.002/86400.0)

    ! Mesozooplankton quadratic mortality rate constant [(mmol C m-3)-1 s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesqmor', wombat%mesqmor, 0.75/86400.0)

    ! Prey switching exponent for microzooplantkon [van Leeuwen et al. (2013), J. Theor. Biol.]
    ! when <1, more even feeding across prey items
    ! when =1, grazing proportional to prey biomasses
    ! when >1, overweighting abundant prey and downweighting scarce prey
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zoopreyswitch', wombat%zoopreyswitch, 1.8)

    ! Prey switching exponent for mesozooplantkon [van Leeuwen et al. (2013), J. Theor. Biol.]
    ! when <1, more even feeding across prey items
    ! when =1, grazing proportional to prey biomasses
    ! when >1, overweighting abundant prey and downweighting scarce prey
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mespreyswitch', wombat%mespreyswitch, 1.8)

    ! Small detritus hydrolyzation rate constant [(mmol C m-3)-1 s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('sdetlrem', wombat%sdetlrem, 0.7/86400.0)

    ! Detritus hydrolyzation rate constant in sediments [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('detlrem_sed', wombat%detlrem_sed, 0.005/86400.0)

    ! Porosity of sinking small detritus [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('sdetphi', wombat%sdetphi, 0.25)

    ! Porosity of sinking large aggregated detritus [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('ldetphi', wombat%ldetphi, 0.75)

    ! Base radius of nanophytoplankton [µm]
    !  2–20 µm in diameter
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phyrad0', wombat%phyrad0, 10.0)

    ! Base radius of microphytoplankton [µm]
    !  20–200 µm in diameter
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diarad0', wombat%diarad0, 50.0)

    ! Base radius of microzooplankton [µm]
    !  takes into account fecal pellets / waste produced
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zoorad0', wombat%zoorad0, 30.0)

    ! Base radius of mesoozooplankton [µm]
    !  takes into account fecal pellets / waste produced
    !  [Hatton et al. 2021 Sci Adv - near 1mm]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesrad0', wombat%mesrad0, 1000.0)

    ! Density of organic small detritus [kg/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('sdetrho', wombat%sdetrho, 1375.0)

    ! Density of calcium carbonate [kg/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('caco3rho', wombat%caco3rho, 2710.0)

    ! Density of biogenic opal [kg/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bsirho', wombat%bsirho, 2000.0)

    ! Phytoplankton biomass threshold to scale recycling [mmolC/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phybiot', wombat%phybiot, 1.0)

    ! Microphytoplankton biomass threshold to scale recycling [mmolC/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diabiot', wombat%diabiot, 0.5)

    ! CaCO3 dissolution rate constant (base rate) [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('caco3lrem', wombat%caco3lrem, 0.01/86400.0)

    ! CaCO3 dissolution rate constant in sediments [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('caco3lrem_sed', wombat%caco3lrem_sed, 0.01/86400.0)

    ! CaCO3 inorganic fraction (PIC:POC) [mol C (mol C)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('f_inorg', wombat%f_inorg, 0.04)

    ! CaCO3 dissolution factor due to calcite undersaturation (s-1)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('disscal', wombat%disscal, 0.10/86400.0)

    ! CaCO3 dissolution factor due to aragonite undersaturation (s-1)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('dissara', wombat%dissara, 0.10/86400.0)

    ! CaCO3 dissolution factor due to detritus remineralisation creating
    !  anoxic microenvironment (mol C in CaCO3 per mol C remineralised)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('dissdet', wombat%dissdet, 0.20)

    ! Background concentration of weak iron-binding ligand [µmol/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('ligW', wombat%ligW, 1.7)

    ! Background concentration of strong iron-binding ligand [µmol/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('ligS', wombat%ligS, 0.4)

    ! Set floor for dissolved iron concentration based on the measurement detection limit [µmol/m3]
    ! Worsford et al., 2014 Mar. Chem. says anywhere between 10 - 50 pM
    !-----------------------------------------------------------------------
    call g_tracer_add_param('dfefloor', wombat%dfefloor, 0.025)

    ! Scavenging of Fe` onto biogenic particles [(mmolC/m3)-1 s-1]
    !-----------------------------------------------------------------------
    ! Ye et al., 2011 (Biogeosciences) find scavenging rates of 30 - 750
    ! (kg/m3)-1 day-1 in mesocosm experiments. Assuming that there are
    ! 40,000 mmol C kg-1 (1 kg of pure carbon contains 83 mol and assuming
    ! that half of marine organic particles are pure carbon by mass means
    ! that roughly 40,000 mmol C kg-1), this translates to scavenging rates
    ! of 0.001 to 0.02 (mmol mass particles / m3)-1 day-1.
    ! NOTE that scavenging becomes less important when colloids dominate.
    call g_tracer_add_param('kscav_dfe', wombat%kscav_dfe, 0.01/86400.0)

    ! Coagulation of dFe onto organic particles [(mmolC/m3)-1 s-1]
    !-----------------------------------------------------------------------
    ! Colloigal coagulation rates are the principal way to remove dFe at high
    ! concentrations. Effectively, when `do_colloidal_shunt == .true.`, the
    ! `kcoag_dfe` parameter sets the maximum dFe concentration via:
    !    [dFe remaining in solution] = dFe_sources / kcoag_dfe
    !  1e-5 ---> coagulation at roughly 0.1 per day in productive surface waters
    !            and 1/1000 per day in deep ocean
    !  1e-6 ---> coagulation at roughly 0.01 per day in productive surface waters
    !            and 1/10000 per day in deep ocean
    call g_tracer_add_param('kcoag_dfe', wombat%kcoag_dfe, 1e-5/86400.0)

    ! Rate of aggregation of colloidal iron into authigenic Fe particles [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('kagg_col', wombat%kagg_col, 0.1/86400.0)

    ! Half-saturation coefficient modulating aggregation of colloidal iron [µmolFe/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('kagg_kcol', wombat%kagg_kcol, 2.0)

    ! Rate of dissolution of authigenic iron into dissolved Fe [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('ksafe_dfe', wombat%ksafe_dfe, 1e-4/86400)

    ! Rate of dissolution of larger authigenic iron into dissolved Fe [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('klafe_dfe', wombat%klafe_dfe, 1e-4/86400)

    ! Sinking speed of small authigenic iron (oxyhydroxide) [m/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('wsafe', wombat%wsafe, 0.5/86400.0)

    ! Sinking speed of large authigenic iron (oxyhydroxide) [m/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('wlafe', wombat%wlafe, 5.0/86400.0)

    ! Intercept of temperature-dependency of biogenic silica dissolution [Kamatani 1982 Marine Biology]
    !  According to Kamatani (1982), this intercept varies between -7.35 to -10.38 depending on
    !  the diatom species and possibly the degree roughness and biological film strength
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bsi_alpha', wombat%bsi_alpha, -10.0)

    ! Factor increase in biogenic silica dissolution caused by bacterial activity [dimensionless]
    !  Bidle & Azam (1999) found 10x increase in dissolution with bacterial activity
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bsi_fbac', wombat%bsi_fbac, 10.0)

    ! Half-saturation coefficient modulating increase in biogenic silica dissolution due to bacterial activity [mmolC/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bsi_kbac', wombat%bsi_kbac, 0.5)

    ! Rate of biogenic silica dissolution in the sediments when at near-total undersaturation [s-1]
    !-----------------------------------------------------------------------
    !  Maximum rate of ~1.0 nmol/g biogenic silica s-1 at 0 degrees and 0 cm sediment depth at near-total
    !  undersaturation, equivalent to 2.8e-8 /s dissolution rate [Fig. 11 in Van Cappellen & Qiu, 1997]
    !  - this very low rate accounts for the factors in sediments that retard dissolution [Van Cappellen et al., 2002 GBC]
    call g_tracer_add_param('bsilrem_sed', wombat%bsilrem_sed, 2.8e-8)

    ! Ammonia Oxidizing Archaea half saturation constant for NH4 uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoa_knh4', wombat%aoa_knh4, 0.1)

    ! Ammonia Oxidizing Archaea diffusive uptake limit for oxygen [(mmol C biomass m-3)-1 s-1)]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoa_poxy', wombat%aoa_poxy, 275.0/86400.0)

    ! Ammonia Oxidizing Archaea biomass yield per NH4 [mol NH4 (mol C biomass)-1]
    !-----------------------------------------------------------------------
    ! NOTE: in the WOMBAT-mid documentation, aoa_ynh4 is written as the inverse
    !       of what it is here in the code, in units of mol Biomass per mol NH4
    call g_tracer_add_param('aoa_ynh4', wombat%aoa_ynh4, 11.0)

    ! Ammonia Oxidizing Archaea biomass yield per O2 [mol O2 (mol C biomass)-1]
    !-----------------------------------------------------------------------
    ! NOTE: in the WOMBAT-mid documentation, aoa_yoxy is written as the inverse
    !       of what it is here in the code, in units of mol Biomass per mol O2
    call g_tracer_add_param('aoa_yoxy', wombat%aoa_yoxy, 15.5)

    ! Ammonia Oxidizing Archaea biomass carbon to nitrogen ratio [mol C (mol N)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoa_C2N', wombat%aoa_C2N, 5.0)

    ! Ammonia Oxidizing Archaea biomass carbon to iron ratio [mol C (mol Fe)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoa_C2Fe', wombat%aoa_C2Fe, 1.0/20e-6)

    ! Ammonia Oxidizing Archaea linear mortality rate constant [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoalmor', wombat%aoalmor, 0.005/86400.0)

    ! Ammonia Oxidizing Archaea quadratic mortality rate constant [(mmol C m-3)-1 s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoaqmor', wombat%aoaqmor, 0.001/86400.0)

    ! Large heterotrophic bacteria maximum rate of uptake of DOC [mmol C m-3 s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('lbac_Vmax_doc', wombat%lbac_Vmax_doc, 6.7/86400.0)

    ! Large heterotrophic bacteria maximum rate of uptake of NO3 [mmol N m-3 s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('lbac_Vmax_no3', wombat%lbac_Vmax_no3, 7.2/86400.0)

    ! Large heterotrophic bacteria maximum rate of uptake of dFe [mmol Fe m-3 s-1]
    !-----------------------------------------------------------------------
    ! From Fourquez et al., 2020 Frontiers in Marine Science: Heterotrophic bacteria
    ! took up dFe at a rate of 100 pmol L-1 day-1 --> 0.00010 mmol m-3 day-1
    ! in unfiltered seawater when they added Fe+C
    call g_tracer_add_param('lbac_Vmax_dFe', wombat%lbac_Vmax_dFe, 0.00010/86400.0)

    ! Large heterotrophic bacteria diffusive uptake limit of O2 [(mmol C biomass m-3)-1 s-1)]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('lbac_poxy', wombat%lbac_poxy, 450.0/86400.0)

    ! Large heterotrophic bacteria half saturation constant for nitrate uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('lbac_kno3', wombat%lbac_kno3, 15.0)

    ! Large heterotrophic bacteria half saturation constant for DOC uptake [mmolC/m3]
    !-----------------------------------------------------------------------
    ! DON is preferentially targeted by heterotrophs for remineralisation over DOC
    ! (Letscher & Moore, 2015 GBC; Hach et al., 2020 Sci. Rep; Zakem et al., 2019 GBC)
    call g_tracer_add_param('lbac_kdoc', wombat%lbac_kdoc, 60.0)

    ! Large heterotrophic bacteria half saturation constant for dissolved iron uptake [µmolFe/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('lbac_kfer', wombat%lbac_kfer, 0.35)

    ! Large heterotrophic bacteria degree of partial oxidation [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('lbac_alpha', wombat%lbac_alpha, 0.25)

    ! Large heterotrophic bacteria fraction of electrons to biosynthesis [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('lbac_fele', wombat%lbac_fele, 0.15)

    ! Small heterotrophic bacteria maximum rate of uptake of DOC [mmol C m-3 s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('sbac_Vmax_doc', wombat%sbac_Vmax_doc, 6.7/86400.0)

    ! Small heterotrophic bacteria maximum rate of uptake of NO3 [mmol N m-3 s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('sbac_Vmax_no3', wombat%sbac_Vmax_no3, 7.2/86400.0)

    ! Small heterotrophic bacteria maximum rate of uptake of dFe [mmol Fe m-3 s-1]
    !-----------------------------------------------------------------------
    ! From Fourquez et al., 2020 Frontiers in Marine Science: Heterotrophic bacteria
    ! took up dFe at a rate of 100 pmol L-1 day-1 --> 0.00010 mmol m-3 day-1
    ! in unfiltered seawater when they added Fe+C
    call g_tracer_add_param('sbac_Vmax_dFe', wombat%sbac_Vmax_dFe, 0.00010/86400.0)

    ! Small heterotrophic bacteria diffusive uptake limit of O2 [(mmol C biomass m-3)-1 s-1)]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('sbac_poxy', wombat%sbac_poxy, 450.0/86400.0)

    ! Small heterotrophic bacteria half saturation constant for nitrate uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('sbac_kno3', wombat%sbac_kno3, 15.0)

    ! Small heterotrophic bacteria half saturation constant for DOC uptake [mmolC/m3]
    !-----------------------------------------------------------------------
    ! DON is preferentially targeted by heterotrophs for remineralisation over DOC
    ! (Letscher & Moore, 2015 GBC; Hach et al., 2020 Sci. Rep; Zakem et al., 2019 GBC)
    call g_tracer_add_param('sbac_kdoc', wombat%sbac_kdoc, 60.0)

    ! Small heterotrophic bacteria half saturation constant for dissolved iron uptake [µmolFe/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('sbac_kfer', wombat%sbac_kfer, 0.35)

    ! Small heterotrophic bacteria fraction of electrons to biosynthesis [dimensionless]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('sbac_fele', wombat%sbac_fele, 0.15)

    ! Heterotrophic bacteria biomass carbon to nitrogen ratio [mol C (mol N)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac_C2N', wombat%bac_C2N, 5.0)

    ! Heterotrophic bacteria biomass carbon to iron ratio [mol C (mol Fe)-1]
    !-----------------------------------------------------------------------
    ! From Fourquez et al., 2020 Frontiers in Marine Science: Heterotrophic bacteria
    ! measured 20 - 50 µmol Fe per mol C in unfiltered seawater samples
    call g_tracer_add_param('bac_C2Fe', wombat%bac_C2Fe, 1.0/40e-6)

    ! Heterotrophic bacteria linear mortality rate constant [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('baclmor', wombat%baclmor, 0.005/86400.0)

    ! Heterotrophic bacteria quadratic mortality rate constant [(mmol C m-3)-1 s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bacqmor', wombat%bacqmor, 0.05/86400.0)

    ! Anammox bacteria half saturation constant for ammonium uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoxkn', wombat%aoxkn, 0.5)

    ! Anammox bacteria maximum growth * biomass rate [s-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoxmumax', wombat%aoxmumax, 0.0025/86400.0)

    ! Bottom thickness [m]
    !-----------------------------------------------------------------------
    ! Thickness over which tracer values are integrated to define the bottom layer
    call g_tracer_add_param('bottom_thickness', wombat%bottom_thickness, 0.1)

    ! Nested timestep for the ecosystem model [s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('dt_npzd', wombat%dt_npzd, 900.)

    ! Global average surface salinity used for virtual flux correction [g/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('sal_global', wombat%sal_global, 34.6)

    ! Global average surface dic used for virtual flux correction [mol/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('dic_global', wombat%dic_global, 1.90e-3)

    ! Global average surface alk used for virtual flux correction [mol/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('alk_global', wombat%alk_global, 2.225e-3)

    ! Global average surface no3 used for virtual flux correction [mol/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('no3_global', wombat%no3_global, 4.512e-06)

    ! Global average surface no3 used for virtual flux correction [mol/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('nh4_global', wombat%nh4_global, 0.05e-06)

    call g_tracer_end_param_list(package_name)

  end subroutine user_add_params

  !#######################################################################
  !
  ! This is an internal sub, not a public interface.
  ! Add all the tracers to be used in this module.
  !
  subroutine user_add_tracers(tracer_list)
    type(g_tracer_type), pointer, intent(in) :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'user_add_tracers'
    real                                    :: as_coeff_wombatmid

    !=======================================================================
    ! Parameters
    !=======================================================================
    ! Add here only the parameters that are required at the time of registeration
    ! (to make flux exchanging ocean tracers known for all PE's)

    ! Air-sea gas exchange coefficient presented in OCMIP2 protocol.
    !-----------------------------------------------------------------------
    ! From Wanninkhof 1992 for steady wind speed (in m/s)
    as_coeff_wombatmid = 0.31 / 3.6e5

    call g_tracer_start_param_list(package_name)

    call g_tracer_add_param('ice_restart_file', wombat%ice_restart_file, 'ice_wombatmid.res.nc')
    call g_tracer_add_param('ocean_restart_file', wombat%ocean_restart_file, 'ocean_wombatmid.res.nc')
    call g_tracer_add_param('IC_file', wombat%IC_file, '')

    call g_tracer_end_param_list(package_name)

    ! Set Restart files
    call g_tracer_set_files(ice_restart_file = wombat%ice_restart_file, &
        ocean_restart_file = wombat%ocean_restart_file )

    !=======================================================================
    ! Specify all tracers of this module
    !=======================================================================
    !
    ! User adds one call for each tracer below!
    ! User should specify if fluxes must be extracted from boundary by passing
    ! one or more of the following methods as .true. and provide the corresponding
    ! parameters array methods: flux_gas, flux_runoff, flux_wetdep, flux_drydep.
    ! Pass an init_value arg if the tracers should be initialized to a nonzero
    ! value everywhere otherwise they will be initialized to zero.
    !
    ! dts: diagnostic (prog = .false.) tracers added here are automatically
    ! registered for restart but not for horizontal advection and diffusion. All
    ! tracer fields are registered for diag output.

    !=======================================================================
    ! Prognostic Tracers
    !=======================================================================

    ! Nitrate
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Nitrate
    call g_tracer_add(tracer_list, package_name, &
        name = 'no3', &
        longname = 'Nitrate', &
        units = 'mol/kg', &
        prog = .true., &
        flux_bottom = .true., &
        flux_virtual = .true.)

    ! Ammonium
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Ammonium
    call g_tracer_add(tracer_list, package_name, &
        name = 'nh4', &
        longname = 'Ammonium', &
        units = 'mol/kg', &
        prog = .true., &
        flux_virtual = .true.)

    ! Silicic acid (H4SiO4)
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'sil', &
        longname = 'Silicic acid', &
        units = 'mol/kg', &
        prog = .true., &
        flux_bottom = .true., &
        flux_virtual = .true.)

    ! Phytoplankton
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Phytoplankton
    call g_tracer_add(tracer_list, package_name, &
        name = 'phy', &
        longname = 'Phytoplankton', &
        units = 'mol/kg', &
        prog = .true.)

    ! Phytoplankton Chlorophyll
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Phytoplankton
    call g_tracer_add(tracer_list, package_name, &
        name = 'pchl', &
        longname = 'Phytoplankton chlorophyll', &
        units = 'mol/kg', &
        prog = .true.)

    ! Phytoplankton Iron content
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Phytoplankton
    call g_tracer_add(tracer_list, package_name, &
        name = 'phyfe', &
        longname = 'Phytoplankton iron content', &
        units = 'mol/kg', &
        prog = .true.)

    ! Microphytoplankton (diatoms)
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Microphytoplankton (diatoms)
    call g_tracer_add(tracer_list, package_name, &
        name = 'dia', &
        longname = 'Microphytoplankton (diatoms)', &
        units = 'mol/kg', &
        prog = .true.)

    ! Microphytoplankton (diatoms) Chlorophyll
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Microphytoplankton (diatoms)
    call g_tracer_add(tracer_list, package_name, &
        name = 'dchl', &
        longname = 'Microphytoplankton (diatoms) chlorophyll', &
        units = 'mol/kg', &
        prog = .true.)

    ! Microphytoplankton (diatoms) Iron content
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Microphytoplankton (diatoms)
    call g_tracer_add(tracer_list, package_name, &
        name = 'diafe', &
        longname = 'Microphytoplankton (diatoms) iron content', &
        units = 'mol/kg', &
        prog = .true.)

    ! Microphytoplankton (diatoms) Silicon content
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Microphytoplankton (diatoms)
    call g_tracer_add(tracer_list, package_name, &
        name = 'diasi', &
        longname = 'Microphytoplankton (diatoms) silicon content', &
        units = 'mol/kg', &
        prog = .true.)

    ! Oxygen
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'o2', &
        longname = 'Oxygen', &
        units = 'mol/kg', &
        prog = .true., &
        flux_gas = .true.,  &
        flux_bottom = .true., &
        flux_gas_name = 'o2_flux', &
        flux_gas_type = 'air_sea_gas_flux_generic', &
        flux_gas_molwt = WTMO2, &
        flux_gas_param = [ as_coeff_wombatmid, 9.7561e-06 ], & ! dts: param(2) converts Pa -> atm
        flux_gas_restart_file = 'ocean_wombatmid_airsea_flux.res.nc')

    ! Zooplankton
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'zoo', &
        longname = 'Zooplankton', &
        units = 'mol/kg', &
        prog = .true.)

    ! Zooplankton iron content
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'zoofe', &
        longname = 'Zooplankton iron content', &
        units = 'mol/kg', &
        prog = .true.)

    ! Mesozooplankton
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'mes', &
        longname = 'Mesozooplankton', &
        units = 'mol/kg', &
        prog = .true.)

    ! Mesozooplankton iron content
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'mesfe', &
        longname = 'Mesozooplankton iron content', &
        units = 'mol/kg', &
        prog = .true.)

    ! Small detritus
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'sdet', &
        longname = 'Small detritus', &
        units = 'mol/kg', &
        prog = .true., &
        move_vertical = .true., &
        btm_reservoir = .true.)

    ! Small detrital iron content
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'sdetfe', &
        longname = 'Small detrital iron content', &
        units = 'mol/kg', &
        prog = .true., &
        move_vertical = .true., &
        btm_reservoir = .true.)

    ! Large detritus
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'ldet', &
        longname = 'Large detritus', &
        units = 'mol/kg', &
        prog = .true., &
        move_vertical = .true., &
        btm_reservoir = .true.)

    ! Large detrital iron content
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'ldetfe', &
        longname = 'Large detrital iron content', &
        units = 'mol/kg', &
        prog = .true., &
        move_vertical = .true., &
        btm_reservoir = .true.)

    ! Large detrital silicon content
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'ldetsi', &
        longname = 'Large detrital silicon content', &
        units = 'mol/kg', &
        prog = .true., &
        move_vertical = .true., &
        btm_reservoir = .true.)

    ! Small molecules (LMW) of dissolved organic carbon
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'sdoc', &
        longname = 'Small dissolved organic carbon', &
        units = 'mol/kg', &
        flux_bottom = .true., &
        prog = .true.)

    ! Large molecules (HMW) of dissolved organic carbon
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'ldoc', &
        longname = 'Large dissolved organic carbon', &
        units = 'mol/kg', &
        flux_bottom = .true., &
        prog = .true.)

    ! Large molecules (HMW) of dissolved organic nitrogen
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'ldon', &
        longname = 'Large dissolved organic nitrogen', &
        units = 'mol/kg', &
        flux_bottom = .true., &
        prog = .true.)

    ! Large, sharing heterotrophic bacteria
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'lbac', &
        longname = 'Large, sharing heterotrophic bacteria', &
        units = 'mol/kg', &
        prog = .true.)

    ! Small, selfish heterotrophic bacteria
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'sbac', &
        longname = 'Small, selfish heterotrophic bacteria', &
        units = 'mol/kg', &
        prog = .true.)

    ! Ammonia oxidizing archaea
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'aoa', &
        longname = 'Ammonia oxidizing archaea', &
        units = 'mol/kg', &
        prog = .true.)

    ! CaCO3
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'caco3', &
        longname = 'CaCO3', &
        units = 'mol/kg', &
        prog = .true., &
        move_vertical = .true., &
        btm_reservoir = .true.)

    ! DIC (Dissolved inorganic carbon)
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'dic', &
        longname = 'Dissolved Inorganic Carbon', &
        units = 'mol/kg', &
        prog = .true., &
        flux_gas = .true., &
        flux_bottom = .true., &
        flux_gas_name = 'co2_flux', &
        flux_gas_type = 'air_sea_gas_flux_generic', &
        flux_gas_molwt = WTMCO2, &
        flux_gas_param = [ as_coeff_wombatmid, 9.7561e-06 ], & ! dts: param(2) converts Pa -> atm
        flux_gas_restart_file = 'ocean_wombatmid_airsea_flux.res.nc', &
        flux_virtual = .true.)

    ! DICp (preformed Dissolved inorganic carbon)
    ! dts: Note, we use flux_virtual=.true. only to ensure that an stf array is allocated for dicp.
    ! The dicp stf is set to equal the dic stf in update_from_coupler.
    !-----------------------------------------------------------------------
    if (do_tracer_dicp) then
      call g_tracer_add(tracer_list, package_name, &
          name = 'dicp', &
          longname = 'preformed Dissolved Inorganic Carbon', &
          units = 'mol/kg', &
          prog = .true., &
          flux_virtual = .true.)
    endif

    ! DICr (remineralised dissolved inorganic carbon)
    !-----------------------------------------------------------------------
    if (do_tracer_dicr) then
      call g_tracer_add(tracer_list, package_name, &
          name = 'dicr', &
          longname = 'remineralised Dissolved Inorganic Carbon', &
          units = 'mol/kg', &
          prog = .true., &
          flux_bottom = .true.)
    endif

    ! Alk (Total carbonate alkalinity)
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'alk', &
        longname = 'Alkalinity', &
        units = 'mol/kg', &
        prog = .true., &
        flux_bottom = .true., &
        flux_virtual = .true.)

    ! Dissolved Iron
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'fe', &
        longname = 'Dissolved Iron', &
        units = 'mol/kg', &
        prog = .true., &
        flux_drydep = .true., &
        flux_param  = [ 1.0 ], & ! dts: conversion to mol/m2/s done in data_table
        flux_bottom = .true.)

    ! Small authigenic Fe
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'safe', &
        longname = 'Small authigenic iron', &
        units = 'mol/kg', &
        prog = .true., &
        move_vertical = .true., &
        btm_reservoir = .true.)

    ! Large authigenic Fe (from big aggregates)
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'lafe', &
        longname = 'Large authigenic iron', &
        units = 'mol/kg', &
        prog = .true., &
        move_vertical = .true., &
        btm_reservoir = .true.)


    !=======================================================================
    ! Diagnostic Tracers
    !=======================================================================

    ! Detritus sitting at base of column as sediment
    !-----------------------------------------------------------------------
    ! dts: included here so included in restart
    call g_tracer_add(tracer_list, package_name, &
        name = 'det_sediment', &
        longname = 'Detritus at base of column as sediment', &
        units = 'mol m-2', &
        prog = .false.)

    ! Detrital iron sitting at base of column as sediment
    !-----------------------------------------------------------------------
    ! dts: included here so included in restart
    call g_tracer_add(tracer_list, package_name, &
        name = 'detfe_sediment', &
        longname = 'Detrital iron at base of column as sediment', &
        units = 'mol m-2', &
        prog = .false.)

    ! Detrital silicon sitting at base of column as sediment
    !-----------------------------------------------------------------------
    ! dts: included here so included in restart
    call g_tracer_add(tracer_list, package_name, &
        name = 'detsi_sediment', &
        longname = 'Detrital silicon at base of column as sediment', &
        units = 'mol m-2', &
        prog = .false.)

    ! CaCO3 sitting at base of column as sediment
    !-----------------------------------------------------------------------
    ! dts: included here so included in restart
    call g_tracer_add(tracer_list, package_name, &
        name  = 'caco3_sediment', &
        longname = 'CaCO3 at base of column as sediment', &
        units = 'mol m-2', &
        prog = .false.)

  end subroutine user_add_tracers

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATmid_update_from_coupler">
  !  <OVERVIEW>
  !     Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !    Some tracer fields could be modified after values are obtained from
  !    the coupler. This subroutine is the place for specific tracer
  !    manipulations. In WOMBATmid we apply virtual flux corrections due
  !    to salt flux restoring/correction here.
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATmid_update_from_coupler(tracer_list)
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !
  !  <IN NAME="salt_flux_added" TYPE="real, dimension(ilb:,jlb:), optional">
  !   Surface salt flux into ocean from restoring or flux adjustment [g/m2/s]
  !  </IN>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATmid_update_from_coupler(tracer_list, ilb, jlb, salt_flux_added)
    type(g_tracer_type), pointer, intent(in) :: tracer_list
    integer, intent(in)                      :: ilb, jlb
    real, dimension(ilb:,jlb:), intent(in)   :: salt_flux_added

    integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau)

    ! Account for virtual fluxes due to salt flux restoring/correction
    !-----------------------------------------------------------------------
    call g_tracer_get_pointer(tracer_list, 'no3', 'stf', wombat%p_no3_stf)
    wombat%no3_vstf(:,:) = (wombat%no3_global / wombat%sal_global) * salt_flux_added(:,:) ! [mol/m2/s]
    wombat%p_no3_stf(:,:) = wombat%p_no3_stf(:,:) + wombat%no3_vstf(:,:) ! [mol/m2/s]

    call g_tracer_get_pointer(tracer_list, 'nh4', 'stf', wombat%p_nh4_stf)
    wombat%nh4_vstf(:,:) = (wombat%nh4_global / wombat%sal_global) * salt_flux_added(:,:) ! [mol/m2/s]
    wombat%p_nh4_stf(:,:) = wombat%p_nh4_stf(:,:) + wombat%nh4_vstf(:,:) ! [mol/m2/s]

    call g_tracer_get_pointer(tracer_list, 'dic', 'stf', wombat%p_dic_stf)
    wombat%dic_vstf(:,:) = (wombat%dic_global / wombat%sal_global) * salt_flux_added(:,:) ! [mol/m2/s]
    wombat%p_dic_stf(:,:) = wombat%p_dic_stf(:,:) + wombat%dic_vstf(:,:) ! [mol/m2/s]

    call g_tracer_get_pointer(tracer_list, 'alk', 'stf', wombat%p_alk_stf)
    wombat%alk_vstf(:,:) = (wombat%alk_global / wombat%sal_global) * salt_flux_added(:,:) ! [mol/m2/s]
    wombat%p_alk_stf(:,:) = wombat%p_alk_stf(:,:) + wombat%alk_vstf(:,:) ! [mol/m2/s]

    if (do_tracer_dicp) then
      ! Set dicp stf equal to dic stf
      call g_tracer_set_values(tracer_list, 'dicp', 'stf', wombat%p_dic_stf, isd, jsd)
    endif

  end subroutine generic_WOMBATmid_update_from_coupler

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATmid_update_from_bottom">
  !  <OVERVIEW>
  !   Set values of bottom fluxes and reservoirs
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   Some tracers could have bottom fluxes and reservoirs. This subroutine
  !   is the place for specific tracer manipulations.
  !   In WOMBATmid, remineralization from the sediment tracers (which
  !   requires temperature) is done in update_from_source. Deposition from
  !   sinking is handled here.
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATmid_update_from_bottom(tracer_list, dt, tau)
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !
  !  <IN NAME="dt" TYPE="real">
  !   Time step increment
  !  </IN>
  !
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index to be used for %field
  !  </IN>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATmid_update_from_bottom(tracer_list, dt, tau, model_time)
    type(g_tracer_type), pointer, intent(in) :: tracer_list
    real, intent(in)                         :: dt
    integer, intent(in)                      :: tau
    type(time_type), intent(in)              :: model_time

    integer                         :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, i, j
    real, dimension(:,:,:), pointer :: grid_tmask
    real                            :: orgflux
    logical                         :: used

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, &
        grid_tmask=grid_tmask)

    ! Move bottom reservoirs to sediment tracers
    !-----------------------------------------------------------------------
    call g_tracer_get_values(tracer_list, 'sdet', 'btm_reservoir', wombat%sdet_btm, isd, jsd)
    call g_tracer_get_values(tracer_list, 'sdetfe', 'btm_reservoir', wombat%sdetfe_btm, isd, jsd)
    call g_tracer_get_values(tracer_list, 'ldet', 'btm_reservoir', wombat%ldet_btm, isd, jsd)
    call g_tracer_get_values(tracer_list, 'ldetfe', 'btm_reservoir', wombat%ldetfe_btm, isd, jsd)
    call g_tracer_get_values(tracer_list, 'ldetsi', 'btm_reservoir', wombat%ldetsi_btm, isd, jsd)
    call g_tracer_get_values(tracer_list, 'caco3', 'btm_reservoir', wombat%caco3_btm, isd, jsd)
    call g_tracer_get_values(tracer_list, 'safe', 'btm_reservoir', wombat%safe_btm, isd, jsd)
    call g_tracer_get_values(tracer_list, 'lafe', 'btm_reservoir', wombat%lafe_btm, isd, jsd)

    ! Calculate burial of deposited detritus (Dunne et al., 2007)
    wombat%fbury(:,:) = 0.0
    if (do_burial) then
      do i = isc, iec
        do j = jsc, jec
          orgflux = (wombat%sdet_btm(i,j) + wombat%ldet_btm(i,j)) / dt * 86400 * 1e3 ! mmol C m-2 day-1
          wombat%fbury(i,j) = 0.013 + 0.53 * orgflux**2.0 / (7.0 + orgflux)**2.0  ! Eq. 3 Dunne et al. 2007
        enddo
      enddo
    endif

    call g_tracer_get_pointer(tracer_list, 'det_sediment', 'field', wombat%p_det_sediment)
    wombat%p_det_sediment(:,:,1) = wombat%p_det_sediment(:,:,1) + (wombat%sdet_btm(:,:) + wombat%ldet_btm(:,:)) * (1.0-wombat%fbury(:,:)) ! [mol/m2]
    call g_tracer_set_values(tracer_list, 'sdet', 'btm_reservoir', 0.0)
    call g_tracer_set_values(tracer_list, 'ldet', 'btm_reservoir', 0.0)

    call g_tracer_get_pointer(tracer_list, 'detfe_sediment', 'field', wombat%p_detfe_sediment)
    wombat%p_detfe_sediment(:,:,1) = wombat%p_detfe_sediment(:,:,1) + (wombat%sdetfe_btm(:,:) + wombat%ldetfe_btm(:,:)) * (1.0-wombat%fbury(:,:)) ! [mol/m2]
    call g_tracer_set_values(tracer_list, 'sdetfe', 'btm_reservoir', 0.0)
    call g_tracer_set_values(tracer_list, 'ldetfe', 'btm_reservoir', 0.0)

    call g_tracer_get_pointer(tracer_list, 'detsi_sediment', 'field', wombat%p_detsi_sediment)
    wombat%p_detsi_sediment(:,:,1) = wombat%p_detsi_sediment(:,:,1) + wombat%ldetsi_btm(:,:) * (1.0-wombat%fbury(:,:)) ! [mol/m2]
    call g_tracer_set_values(tracer_list, 'ldetsi', 'btm_reservoir', 0.0)

    call g_tracer_get_pointer(tracer_list, 'caco3_sediment', 'field', wombat%p_caco3_sediment)
    wombat%p_caco3_sediment(:,:,1) =  wombat%p_caco3_sediment(:,:,1) + wombat%caco3_btm(:,:) * (1.0-wombat%fbury(:,:)) ! [mol/m2]
    call g_tracer_set_values(tracer_list, 'caco3', 'btm_reservoir', 0.0)

    ! Bury all authigenic Fe in sediment
    call g_tracer_set_values(tracer_list, 'safe', 'btm_reservoir', 0.0)
    call g_tracer_set_values(tracer_list, 'lafe', 'btm_reservoir', 0.0)

    ! Send diagnostics
    !-----------------------------------------------------------------------
    if (wombat%id_det_sed_depst > 0) &
      used = g_send_data(wombat%id_det_sed_depst, (wombat%sdet_btm + wombat%ldet_btm) / dt, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_detfe_sed_depst > 0) &
      used = g_send_data(wombat%id_detfe_sed_depst, (wombat%sdetfe_btm + wombat%ldetfe_btm) / dt, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_detsi_sed_depst > 0) &
      used = g_send_data(wombat%id_detsi_sed_depst, wombat%ldetsi_btm / dt, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_caco3_sed_depst > 0) &
      used = g_send_data(wombat%id_caco3_sed_depst, wombat%caco3_btm / dt, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_fbury > 0) &
      used = g_send_data(wombat%id_fbury, wombat%fbury, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  end subroutine generic_WOMBATmid_update_from_bottom

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATmid_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   This is the subroutine to contain most of the biogeochemistry for
  !   calculating the interaction of tracers with each other and with outside
  !   forcings.
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATmid_update_from_source(tracer_list, Temp, Salt, &
  !     dzt, hblt_depth, ilb, jlb, tau, dt, grid_dat, sw_pen, opacity)
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !
  !  <IN NAME="Temp" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean temperature
  !  </IN>
  !
  !  <IN NAME="Salt" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean salinity
  !  </IN>
  !
  !  <IN NAME="dzt" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean layer thickness (meters)
  !  </IN>
  !
  !  <IN NAME="opacity" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean opacity
  !  </IN>
  !
  !  <IN NAME="sw_pen" TYPE="real, dimension(ilb:,jlb:)">
  !   Shortwave peneteration
  !  </IN>
  !
  !  <IN NAME="hblt_depth" TYPE="real, dimension(ilb:,jlb:)">
  !   Depth of actively mixing layer
  !  </IN>
  !
  !  <IN NAME="grid_dat" TYPE="real, dimension(ilb:,jlb:)">
  !   Grid area
  !  </IN>
  !
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index of %field
  !  </IN>
  !
  !  <IN NAME="dt" TYPE="real">
  !   Time step increment
  !  </IN>
  !
  ! </SUBROUTINE>
  subroutine generic_WOMBATmid_update_from_source(tracer_list, Temp, Salt,  &
      rho_dzt, dzt, hblt_depth, ilb, jlb, tau, dt, grid_dat, model_time, nbands, &
      max_wavelength_band, sw_pen_band, opacity_band)
    type(g_tracer_type), pointer, intent(in)   :: tracer_list
    real, dimension(ilb:,jlb:,:), intent(in)   :: Temp, Salt, rho_dzt, dzt
    real, dimension(ilb:,jlb:), intent(in)     :: hblt_depth
    integer, intent(in)                        :: ilb, jlb, tau
    real, intent(in)                           :: dt
    real, dimension(ilb:,jlb:), intent(in)     :: grid_dat
    type(time_type), intent(in)                :: model_time
    integer, intent(in)                        :: nbands
    real, dimension(:), intent(in)             :: max_wavelength_band
    real, dimension(:,ilb:,jlb:), intent(in)   :: sw_pen_band
    real, dimension(:,ilb:,jlb:,:), intent(in) :: opacity_band

    integer                                 :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, tn
    integer                                 :: i, j, k, n, nz, k_bot
    real, dimension(:,:,:), pointer         :: grid_tmask
    integer, dimension(:,:), pointer        :: grid_kmt
    integer, dimension(:,:), allocatable    :: kmeuph ! deepest level of euphotic zone
    integer, dimension(:,:), allocatable    :: k100 ! deepest level less than 100 m
    real                                    :: mmol_m3_to_mol_kg, umol_m3_to_mol_kg
    integer                                 :: ts_npzd ! number of time steps within NPZD model
    real                                    :: dtsb ! number of seconds per NPZD timestep
    real                                    :: rdtts ! 1 / dt
    real, dimension(nbands)                 :: sw_pen
    real                                    :: swpar
    real                                    :: g_zoo, g_mes, Xzoo, I_Xzoo, Xmes, I_Xmes
    real                                    :: no3_p, nh4_p, oxy_p, fe_p, sil_p, sdoc_p, ldoc_p, ldon_p, caco3_p, safe_p
    real                                    :: phy_p, dia_p, zoo_p, mes_p, sdet_p, ldet_p, ldetsi_p, lbac_p, sbac_p, aoa_p, lafe_p
    real                                    :: phyfe_p, diafe_p, diasi_p, pchl_p, dchl_p, zoofe_p, mesfe_p, sdetfe_p, ldetfe_p
    real                                    :: no3_mmolm3, nh4_mmolm3, oxy_mmolm3, fe_umolm3, sil_mmolm3
    real                                    :: sdoc_mmolm3, ldoc_mmolm3, ldon_mmolm3, caco3_mmolm3
    real                                    :: phy_mmolm3, dia_mmolm3, zoo_mmolm3, mes_mmolm3, sdet_mmolm3, ldet_mmolm3
    real                                    :: ldetsi_mmolm3, lbac_mmolm3, sbac_mmolm3, aoa_mmolm3, phyfe_mmolm3, diafe_mmolm3
    real                                    :: I_denom, wzlbac, wzsbac, wzaoa, wzphy, wzdia, wzsdet, wzldet, wzzoo, I_wzsum
    real                                    :: fbc
    real, parameter                         :: epsi = 1.0e-30
    real, parameter                         :: pi = 3.14159265358979
    real, parameter                         :: Rgas = 8.314462168 ! J/(K mol)
    real, parameter                         :: T_star = 647.096 ! K
    real, parameter                         :: rho_star = 322.0 ! kg/m3
    real, parameter                         :: mu_star = 1e-6 ! Pa s
    real, dimension(4)                      :: mu0_H
    integer, dimension(21)                  :: mu1_i, mu1_j
    real, dimension(21)                     :: mu1_H
    real                                    :: P_MPa, rho0, rho_w, at, bt, mu_w, mu_w_iapws
    real                                    :: T_hat, rho_hat, invT_hat, mu0, mu1, poly
    real                                    :: invT_hat2, invT_hat3, temp2, temp3, temp4, temp5
    integer                                 :: ichl, iter
    real                                    :: rad_phy, rad_dia, rad_zoo, rad_mes, rad_sdet, rad_ldet
    real                                    :: mass_sdet, mass_ldet, mass_caco3, mass_bsi, mass_small, mass_large
    real                                    :: w1, w2, rho_small, rho_large
    real                                    :: par_phy_mldsum, par_z_mldsum
    real                                    :: chl, ndet, carb, zchl, zval, sqrt_zval, phy_chlc, dia_chlc
    real                                    :: phy_limnh4, phy_limno3, phy_limdin
    real                                    :: dia_limnh4, dia_limno3, dia_limdin
    real                                    :: phy_pisl, dia_pisl
    real                                    :: zooegeslbacfe, zooegessbacfe, zooegesaoafe, zooegesphyfe, zooegesdiafe, zooegessdetfe
    real                                    :: zooassilbacfe, zooassisbacfe, zooassiaoafe, zooassiphyfe, zooassidiafe, zooassisdetfe
    real                                    :: zooexcrlbacfe, zooexcrsbacfe, zooexcraoafe, zooexcrphyfe, zooexcrdiafe, zooexcrsdetfe
    real                                    :: mesegeslbacfe, mesegessbacfe, mesegesaoafe, mesegesphyfe, mesegesdiafe
    real                                    :: mesegessdetfe, mesegesldetfe, mesegeszoofe
    real                                    :: mesassilbacfe, mesassisbacfe, mesassiaoafe, mesassiphyfe, mesassidiafe
    real                                    :: mesassisdetfe, mesassildetfe, mesassizoofe
    real                                    :: mesexcrlbacfe, mesexcrsbacfe, mesexcraoafe, mesexcrphyfe, mesexcrdiafe
    real                                    :: mesexcrsdetfe, mesexcrldetfe, mesexcrzoofe
    real                                    :: zooexcrlbacn, zooexcrsbacn, mesexcrlbacn, mesexcrsbacn, zooexcraoan, mesexcraoan
    real, dimension(:,:), allocatable       :: ek_bgr, par_bgr_mid, par_bgr_top
    real, dimension(:), allocatable         :: wsink1, wsink2
    real, dimension(4,61)                   :: zbgr
    real, dimension(3)                      :: dbgr, cbgr
    real                                    :: ztemk, I_ztemk, fe_sfe, partic, fescaven, inv1, inv2, fe3
    real                                    :: ligW_K, FeL1_mid, FeL2_mid
    real                                    :: fesol1, fesol2, fesol3, fesol4, fesol5, hp, fe3sol
    real                                    :: feagg1, feagg2, feagg3, feagg4, feagg5
    real                                    :: biof, shear
    real                                    :: phy_Fe2C, dia_Fe2C, zoo_Fe2C, mes_Fe2C, sdet_Fe2C, ldet_Fe2C
    real                                    :: sdom_N2C, ldom_N2C, dia_Si2C
    real                                    :: theta_opt
    real                                    :: phy_minqfe, phy_maxqfe
    real                                    :: dia_minqfe, dia_maxqfe
    real                                    :: hco3
    real                                    :: dzt_bot, dzt_bot_os
    real                                    :: e_pom, e_ldom, e_sdom, e_bac, e_lres
    real                                    :: lbac_cdoc, lbac_coxy, lbac_pdoc, lbac_pco2, lbac_pnh4, lbacydoc_ana
    real                                    :: lbac_cdoc_ana, lbac_cno3_ana, lbac_pdoc_ana, lbac_pco2_ana, lbac_pnh4_ana
    real                                    :: sbac_cdoc, sbac_coxy, sbac_pco2, sbac_pnh4, sbacydoc_ana
    real                                    :: sbac_cdoc_ana, sbac_cno3_ana, sbac_pco2_ana, sbac_pnh4_ana
    real                                    :: bac_Voc, bac_VdFe, bac_Voxy, bac_Vno3
    real                                    :: bac_gC, bac_gFe, bac_gEA
    real                                    :: bac_muana, bac_muaer
    real                                    :: aoa_Vnh4, aoa_Voxy
    real                                    :: K_am_silica, gamma0, alphaH2O, deltaV0, spmvcorrect
    real                                    :: disssi_temp, disssi_usat, disssi_bact
    real, dimension(:,:,:,:), allocatable   :: n_pools, c_pools, si_pools, fe_pools
    logical                                 :: used

    character(len=fm_string_len), parameter :: sub_name = 'generic_WOMBATmid_update_from_source'
    character(len=256), parameter           :: error_header = &
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, &
        grid_tmask=grid_tmask, grid_kmt=grid_kmt)

    ! dts: Note, other generic_tracer modules call this zt. However here we
    ! call this zw to be consistent with WOMBAT v3, which uses MOM5 terminology.
    ! zm here is halfway between interfaces
    wombat%zw(:,:,:) = 0.0
    wombat%zm(:,:,:) = 0.0
    do j = jsc,jec; do i = isc,iec
      wombat%zw(i,j,1) = dzt(i,j,1)
      wombat%zm(i,j,1) = 0.5 * dzt(i,j,1)
    enddo; enddo
    do k = 2,nk; do j = jsc,jec ; do i = isc,iec
      wombat%zw(i,j,k) = wombat%zw(i,j,k-1) + dzt(i,j,k)
      wombat%zm(i,j,k) = wombat%zw(i,j,k-1) + 0.5 * dzt(i,j,k)
    enddo; enddo ; enddo

    ! Some unit conversion factors
    mmol_m3_to_mol_kg = 1.e-3 / wombat%Rho_0
    umol_m3_to_mol_kg = 1.e-3 * mmol_m3_to_mol_kg

    !==========================================================================
    ! Constants for calculating water viscosity (IAPWS 2008)
    !==========================================================================
    mu0_H(1) = 1.67752; mu0_H(2) = 2.20462; mu0_H(3) = 0.6366564; mu0_H(4) = -0.241605
    mu1_i = [ &
      0, 1, 2, 3,       & ! j=0
      0, 1, 2, 3,    5, & ! j=1
      0, 1, 2, 3, 4,    & ! j=2
      0, 1,             & ! j=3
      0,       3,       & ! j=4
                  4,    & ! j=5
               3,    5 &
            ] ! j=6
    mu1_j = [ &
      0, 0, 0, 0,       & ! i=0
      1, 1, 1, 1, 1,    & ! i=1
      2, 2, 2, 2, 2,    & ! i=2
      3, 3,             & ! i=3
      4,       4,       & ! i=4
                  5,    & ! i=5
               6,    6  &
            ] ! i=6
    mu1_H = [ &
        0.520094,        & ! i=0, j=0
        0.0850895,       & ! i=1, j=0
       -1.08374,        & ! i=2, j=0
       -0.289555,       & ! i=3, j=0
        0.222531,       & ! i=0, j=1
        0.999115,       & ! i=1, j=1
        1.88797,        & ! i=2, j=1
        1.26613,        & ! i=3, j=1
        0.120573,       & ! i=5, j=1
       -0.281378,       & ! i=0, j=2
       -0.906851,       & ! i=1, j=2
       -0.772479,       & ! i=2, j=2
       -0.489837,       & ! i=3, j=2
       -0.257040,       & ! i=4, j=2
        0.161913,       & ! i=0, j=3
        0.257399,       & ! i=1, j=3
       -0.0325372,      & ! i=0, j=4
        0.0698452,      & ! i=3, j=4
        0.00872102,     & ! i=4, j=5
        -0.00435673,    & ! i=3, j=6
        -0.000593264     & ! i=5, j=6
            ]

    !==========================================================================
    ! Attenuation coefficients for blue, green and red light due to chlorophyll
    !==========================================================================
    ! Chlorophyll      ! Blue attenuation    ! Green attenuation   ! Red attenuation
    zbgr(1, 1) =  0.010; zbgr(2, 1) = 0.01618; zbgr(3, 1) = 0.07464; zbgr(4, 1) = 0.3780
    zbgr(1, 2) =  0.011; zbgr(2, 2) = 0.01654; zbgr(3, 2) = 0.07480; zbgr(4, 2) = 0.37823
    zbgr(1, 3) =  0.013; zbgr(2, 3) = 0.01693; zbgr(3, 3) = 0.07499; zbgr(4, 3) = 0.37840
    zbgr(1, 4) =  0.014; zbgr(2, 4) = 0.01736; zbgr(3, 4) = 0.07518; zbgr(4, 4) = 0.37859
    zbgr(1, 5) =  0.016; zbgr(2, 5) = 0.01782; zbgr(3, 5) = 0.07539; zbgr(4, 5) = 0.37879
    zbgr(1, 6) =  0.018; zbgr(2, 6) = 0.01831; zbgr(3, 6) = 0.07562; zbgr(4, 6) = 0.37900
    zbgr(1, 7) =  0.020; zbgr(2, 7) = 0.01885; zbgr(3, 7) = 0.07586; zbgr(4, 7) = 0.37923
    zbgr(1, 8) =  0.022; zbgr(2, 8) = 0.01943; zbgr(3, 8) = 0.07613; zbgr(4, 8) = 0.37948
    zbgr(1, 9) =  0.025; zbgr(2, 9) = 0.02005; zbgr(3, 9) = 0.07641; zbgr(4, 9) = 0.37976
    zbgr(1,10) =  0.028; zbgr(2,10) = 0.02073; zbgr(3,10) = 0.07672; zbgr(4,10) = 0.38005
    zbgr(1,11) =  0.032; zbgr(2,11) = 0.02146; zbgr(3,11) = 0.07705; zbgr(4,11) = 0.38036
    zbgr(1,12) =  0.035; zbgr(2,12) = 0.02224; zbgr(3,12) = 0.07741; zbgr(4,12) = 0.38070
    zbgr(1,13) =  0.040; zbgr(2,13) = 0.02310; zbgr(3,13) = 0.07780; zbgr(4,13) = 0.38107
    zbgr(1,14) =  0.045; zbgr(2,14) = 0.02402; zbgr(3,14) = 0.07821; zbgr(4,14) = 0.38146
    zbgr(1,15) =  0.050; zbgr(2,15) = 0.02501; zbgr(3,15) = 0.07866; zbgr(4,15) = 0.38189
    zbgr(1,16) =  0.056; zbgr(2,16) = 0.02608; zbgr(3,16) = 0.07914; zbgr(4,16) = 0.38235
    zbgr(1,17) =  0.063; zbgr(2,17) = 0.02724; zbgr(3,17) = 0.07967; zbgr(4,17) = 0.38285
    zbgr(1,18) =  0.071; zbgr(2,18) = 0.02849; zbgr(3,18) = 0.08023; zbgr(4,18) = 0.38338
    zbgr(1,19) =  0.079; zbgr(2,19) = 0.02984; zbgr(3,19) = 0.08083; zbgr(4,19) = 0.38396
    zbgr(1,20) =  0.089; zbgr(2,20) = 0.03131; zbgr(3,20) = 0.08149; zbgr(4,20) = 0.38458
    zbgr(1,21) =  0.100; zbgr(2,21) = 0.03288; zbgr(3,21) = 0.08219; zbgr(4,21) = 0.38526
    zbgr(1,22) =  0.112; zbgr(2,22) = 0.03459; zbgr(3,22) = 0.08295; zbgr(4,22) = 0.38598
    zbgr(1,23) =  0.126; zbgr(2,23) = 0.03643; zbgr(3,23) = 0.08377; zbgr(4,23) = 0.38676
    zbgr(1,24) =  0.141; zbgr(2,24) = 0.03842; zbgr(3,24) = 0.08466; zbgr(4,24) = 0.38761
    zbgr(1,25) =  0.158; zbgr(2,25) = 0.04057; zbgr(3,25) = 0.08561; zbgr(4,25) = 0.38852
    zbgr(1,26) =  0.178; zbgr(2,26) = 0.04289; zbgr(3,26) = 0.08664; zbgr(4,26) = 0.38950
    zbgr(1,27) =  0.200; zbgr(2,27) = 0.04540; zbgr(3,27) = 0.08775; zbgr(4,27) = 0.39056
    zbgr(1,28) =  0.224; zbgr(2,28) = 0.04811; zbgr(3,28) = 0.08894; zbgr(4,28) = 0.39171
    zbgr(1,29) =  0.251; zbgr(2,29) = 0.05103; zbgr(3,29) = 0.09023; zbgr(4,29) = 0.39294
    zbgr(1,30) =  0.282; zbgr(2,30) = 0.05420; zbgr(3,30) = 0.09162; zbgr(4,30) = 0.39428
    zbgr(1,31) =  0.316; zbgr(2,31) = 0.05761; zbgr(3,31) = 0.09312; zbgr(4,31) = 0.39572
    zbgr(1,32) =  0.355; zbgr(2,32) = 0.06130; zbgr(3,32) = 0.09474; zbgr(4,32) = 0.39727
    zbgr(1,33) =  0.398; zbgr(2,33) = 0.06529; zbgr(3,33) = 0.09649; zbgr(4,33) = 0.39894
    zbgr(1,34) =  0.447; zbgr(2,34) = 0.06959; zbgr(3,34) = 0.09837; zbgr(4,34) = 0.40075
    zbgr(1,35) =  0.501; zbgr(2,35) = 0.07424; zbgr(3,35) = 0.10040; zbgr(4,35) = 0.40270
    zbgr(1,36) =  0.562; zbgr(2,36) = 0.07927; zbgr(3,36) = 0.10259; zbgr(4,36) = 0.40480
    zbgr(1,37) =  0.631; zbgr(2,37) = 0.08470; zbgr(3,37) = 0.10495; zbgr(4,37) = 0.40707
    zbgr(1,38) =  0.708; zbgr(2,38) = 0.09056; zbgr(3,38) = 0.10749; zbgr(4,38) = 0.40952
    zbgr(1,39) =  0.794; zbgr(2,39) = 0.09690; zbgr(3,39) = 0.11024; zbgr(4,39) = 0.41216
    zbgr(1,40) =  0.891; zbgr(2,40) = 0.10374; zbgr(3,40) = 0.11320; zbgr(4,40) = 0.41502
    zbgr(1,41) =  1.000; zbgr(2,41) = 0.11114; zbgr(3,41) = 0.11639; zbgr(4,41) = 0.41809
    zbgr(1,42) =  1.122; zbgr(2,42) = 0.11912; zbgr(3,42) = 0.11984; zbgr(4,42) = 0.42142
    zbgr(1,43) =  1.259; zbgr(2,43) = 0.12775; zbgr(3,43) = 0.12356; zbgr(4,43) = 0.42500
    zbgr(1,44) =  1.413; zbgr(2,44) = 0.13707; zbgr(3,44) = 0.12757; zbgr(4,44) = 0.42887
    zbgr(1,45) =  1.585; zbgr(2,45) = 0.14715; zbgr(3,45) = 0.13189; zbgr(4,45) = 0.43304
    zbgr(1,46) =  1.778; zbgr(2,46) = 0.15803; zbgr(3,46) = 0.13655; zbgr(4,46) = 0.43754
    zbgr(1,47) =  1.995; zbgr(2,47) = 0.16978; zbgr(3,47) = 0.14158; zbgr(4,47) = 0.44240
    zbgr(1,48) =  2.239; zbgr(2,48) = 0.18248; zbgr(3,48) = 0.14701; zbgr(4,48) = 0.44765
    zbgr(1,49) =  2.512; zbgr(2,49) = 0.19620; zbgr(3,49) = 0.15286; zbgr(4,49) = 0.45331
    zbgr(1,50) =  2.818; zbgr(2,50) = 0.21102; zbgr(3,50) = 0.15918; zbgr(4,50) = 0.45942
    zbgr(1,51) =  3.162; zbgr(2,51) = 0.22703; zbgr(3,51) = 0.16599; zbgr(4,51) = 0.46601
    zbgr(1,52) =  3.548; zbgr(2,52) = 0.24433; zbgr(3,52) = 0.17334; zbgr(4,52) = 0.47313
    zbgr(1,53) =  3.981; zbgr(2,53) = 0.26301; zbgr(3,53) = 0.18126; zbgr(4,53) = 0.48080
    zbgr(1,54) =  4.467; zbgr(2,54) = 0.28320; zbgr(3,54) = 0.18981; zbgr(4,54) = 0.48909
    zbgr(1,55) =  5.012; zbgr(2,55) = 0.30502; zbgr(3,55) = 0.19903; zbgr(4,55) = 0.49803
    zbgr(1,56) =  5.623; zbgr(2,56) = 0.32858; zbgr(3,56) = 0.20898; zbgr(4,56) = 0.50768
    zbgr(1,57) =  6.310; zbgr(2,57) = 0.35404; zbgr(3,57) = 0.21971; zbgr(4,57) = 0.51810
    zbgr(1,58) =  7.079; zbgr(2,58) = 0.38154; zbgr(3,58) = 0.23129; zbgr(4,58) = 0.52934
    zbgr(1,59) =  7.943; zbgr(2,59) = 0.41125; zbgr(3,59) = 0.24378; zbgr(4,59) = 0.54147
    zbgr(1,60) =  8.912; zbgr(2,60) = 0.44336; zbgr(3,60) = 0.25725; zbgr(4,60) = 0.55457
    zbgr(1,61) = 10.000; zbgr(2,61) = 0.47804; zbgr(3,61) = 0.27178; zbgr(4,61) = 0.56870

    !===================================================================================
    ! Attenuation coefficients for blue, green and red light due to detritus (m2 / mg N)
    !  Source: Dutkiewicz et al.(2015) Biogeosciences 12, 4447-4481, Fig. 1b
    !          collated into NetCDF file by Mark Baird for EMS model
    !           - csiro_mass_specific_iops_library.nc
    !          assume blue (450-495 nm), green (495-570 nm) and red (620-750 nm)
    !          to create values, we average absorption within these wavelengths
    !===================================================================================
    ! Blue attenuation    ! Green attenuation   ! Red attenuation
    dbgr(1) = 0.01006;    dbgr(2) = 0.009007;   dbgr(3) = 0.007264

    !===================================================================================
    ! Attenuation coefficients for blue, green and red light due to CaCO3 (m2 / kg CaCO3)
    !  Source: Soja-Wozniak et al., 2019 J. Geophys. Res. (Oceans) 124 https://doi.org/10.1029/2019JC014998
    !          collated into NetCDF file by Mark Baird for EMS model
    !           - csiro_mass_specific_iops_library.nc
    !          assume blue (450-495 nm), green (495-570 nm) and red (620-750 nm)
    !          to create values, we average absorption within these wavelengths
    !===================================================================================
    ! Blue attenuation    ! Green attenuation   ! Red attenuation
    cbgr(1) = 1.55641;    cbgr(2) = 3.200139;   cbgr(3) = 20.068027


    !=======================================================================
    ! Surface gas fluxes
    !=======================================================================
    !
    ! Calculate the surface gas fluxes for the next round of exchange. This
    ! is done here to align with other generic_tracer modules (e.g. BLING).
    ! dts: I think this done here in other modules because they calculate
    ! 3D Carbonate ion concentration (co3_ion) here using the FMS_ocmip2_co2calc
    ! routine. The FMS_ocmip2_co2calc routine also calculates co2star, alpha
    ! and pco2surf, so it makes sense to set these values here rather than
    ! recalculating them in set_boundary_values.

    call g_tracer_get_values(tracer_list, 'dic', 'field', wombat%f_dic, isd, jsd, ntau=tau)
    call g_tracer_get_values(tracer_list, 'no3', 'field', wombat%f_no3, isd, jsd, ntau=tau)
    call g_tracer_get_values(tracer_list, 'sil', 'field', wombat%f_sil, isd, jsd, ntau=tau)
    call g_tracer_get_values(tracer_list, 'alk', 'field', wombat%f_alk, isd, jsd, ntau=tau)

    do k = 1,nk !{
     do j = jsc,jec; do i = isc,iec
       wombat%htotallo(i,j) = wombat%htotal_scale_lo * wombat%htotal(i,j,k)
       wombat%htotalhi(i,j) = wombat%htotal_scale_hi * wombat%htotal(i,j,k)
     enddo; enddo

     if (k==1) then !{
       call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,k), &
           Temp(:,:,k), max(1.0, Salt(:,:,k)), &
           max(0.0, wombat%f_dic(:,:,k)), &
           max(0.0, wombat%f_no3(:,:,k) / 16.), &
           max(0.0, wombat%f_sil(:,:,k)), &
           max(0.0, wombat%f_alk(:,:,k)), &
           wombat%htotallo(:,:), wombat%htotalhi(:,:), &
           wombat%htotal(:,:,k), &
           co2_calc=trim(co2_calc), &
           zt=wombat%zw(:,:,k), &
           co2star=wombat%co2_csurf(:,:), alpha=wombat%co2_alpha(:,:), co3_ion=wombat%co3(:,:,k), &
           pCO2surf=wombat%pco2_csurf(:,:), omega_arag=wombat%omega_ara(:,:,k), omega_calc=wombat%omega_cal(:,:,k))

       call g_tracer_set_values(tracer_list, 'dic', 'alpha', wombat%co2_alpha, isd, jsd)
       call g_tracer_set_values(tracer_list, 'dic', 'csurf', wombat%co2_csurf, isd, jsd)

       wombat%co2_star(:,:,1) = wombat%co2_csurf(:,:)

     else

       call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,k), &
           Temp(:,:,k), max(1.0, Salt(:,:,k)), &
           max(0.0, wombat%f_dic(:,:,k)), &
           max(0.0, wombat%f_no3(:,:,k) / 16.), &
           max(0.0, wombat%f_sil(:,:,k)), &
           max(0.0, wombat%f_alk(:,:,k)), &
           wombat%htotallo(:,:), wombat%htotalhi(:,:), &
           wombat%htotal(:,:,k), &
           co2_calc=trim(co2_calc), &
           zt=wombat%zw(:,:,k), &
           co2star=wombat%co2_star(:,:,k), co3_ion=wombat%co3(:,:,k), &
           omega_arag=wombat%omega_ara(:,:,k), omega_calc=wombat%omega_cal(:,:,k))

     endif !} if k.eq.1
    enddo !} do k = 1,nk

    !=======================================================================
    ! Calculate the source terms
    !=======================================================================

    wombat%dynvis_sw(:,:,:) = 1e-3
    wombat%radbio(:,:,:) = 0.0
    wombat%radmid(:,:,:) = 0.0
    wombat%radmld(:,:,:) = 0.0
    wombat%npp3d(:,:,:) = 0.0
    wombat%rpp3d(:,:,:) = 0.0
    wombat%zsp3d(:,:,:) = 0.0
    wombat%phy_mumax(:,:,:) = 0.0
    wombat%phy_mu(:,:,:) = 0.0
    wombat%pchl_mu(:,:,:) = 0.0
    wombat%phy_kni(:,:,:) = 0.0
    wombat%phy_kfe(:,:,:) = 0.0
    wombat%phy_lpar(:,:,:) = 0.0
    wombat%phy_lnit(:,:,:) = 0.0
    wombat%phy_lnh4(:,:,:) = 0.0
    wombat%phy_lno3(:,:,:) = 0.0
    wombat%phy_lfer(:,:,:) = 0.0
    wombat%phy_dfeupt(:,:,:) = 0.0
    wombat%dia_mumax(:,:,:) = 0.0
    wombat%dia_mu(:,:,:) = 0.0
    wombat%dchl_mu(:,:,:) = 0.0
    wombat%dia_kni(:,:,:) = 0.0
    wombat%dia_kfe(:,:,:) = 0.0
    wombat%dia_ksi(:,:,:) = 0.0
    wombat%dia_lpar(:,:,:) = 0.0
    wombat%dia_lnit(:,:,:) = 0.0
    wombat%dia_lnh4(:,:,:) = 0.0
    wombat%dia_lno3(:,:,:) = 0.0
    wombat%dia_lfer(:,:,:) = 0.0
    wombat%dia_lsil(:,:,:) = 0.0
    wombat%dia_dfeupt(:,:,:) = 0.0
    wombat%dia_silupt(:,:,:) = 0.0
    wombat%sileqc(:,:,:) = 0.0
    wombat%disssi(:,:,:) = 0.0
    wombat%bsidiss(:,:,:) = 0.0
    wombat%feIII(:,:,:) = 0.0
    wombat%felig(:,:,:) = 0.0
    wombat%ligK(:,:,:) = 0.0
    wombat%fecol(:,:,:) = 0.0
    wombat%fescasafe(:,:,:) = 0.0
    wombat%fescalafe(:,:,:) = 0.0
    wombat%fesources(:,:,:) = 0.0
    wombat%fesinks(:,:,:) = 0.0
    wombat%fecoag2safe(:,:,:) = 0.0
    wombat%fecoag2lafe(:,:,:) = 0.0
    wombat%safediss(:,:,:) = 0.0
    wombat%lafediss(:,:,:) = 0.0
    wombat%phy_feupreg(:,:,:) = 0.0
    wombat%phy_fedoreg(:,:,:) = 0.0
    wombat%phygrow(:,:,:) = 0.0
    wombat%phydoc(:,:,:) = 0.0
    wombat%phymorl(:,:,:) = 0.0
    wombat%phymorq(:,:,:) = 0.0
    wombat%dia_feupreg(:,:,:) = 0.0
    wombat%dia_fedoreg(:,:,:) = 0.0
    wombat%dia_sidoreg(:,:,:) = 0.0
    wombat%diagrow(:,:,:) = 0.0
    wombat%diadoc(:,:,:) = 0.0
    wombat%diamorl(:,:,:) = 0.0
    wombat%diamorq(:,:,:) = 0.0
    wombat%zooeps(:,:,:) = 0.0
    wombat%zoopreflbac(:,:,:) = 0.0
    wombat%zooprefsbac(:,:,:) = 0.0
    wombat%zooprefaoa(:,:,:) = 0.0
    wombat%zooprefphy(:,:,:) = 0.0
    wombat%zooprefdia(:,:,:) = 0.0
    wombat%zooprefsdet(:,:,:) = 0.0
    wombat%zoograzlbac(:,:,:) = 0.0
    wombat%zoograzsbac(:,:,:) = 0.0
    wombat%zoograzaoa(:,:,:) = 0.0
    wombat%zoograzphy(:,:,:) = 0.0
    wombat%zoograzdia(:,:,:) = 0.0
    wombat%zoograzsdet(:,:,:) = 0.0
    wombat%zoomorl(:,:,:) = 0.0
    wombat%zoomorq(:,:,:) = 0.0
    wombat%zooexcrlbac(:,:,:) = 0.0
    wombat%zooexcrsbac(:,:,:) = 0.0
    wombat%zooexcraoa(:,:,:) = 0.0
    wombat%zooexcrphy(:,:,:) = 0.0
    wombat%zooexcrdia(:,:,:) = 0.0
    wombat%zooexcrsdet(:,:,:) = 0.0
    wombat%zooegeslbac(:,:,:) = 0.0
    wombat%zooegessbac(:,:,:) = 0.0
    wombat%zooegesaoa(:,:,:) = 0.0
    wombat%zooegesphy(:,:,:) = 0.0
    wombat%zooegesdia(:,:,:) = 0.0
    wombat%zooegessdet(:,:,:) = 0.0
    wombat%meseps(:,:,:) = 0.0
    wombat%mespreflbac(:,:,:) = 0.0
    wombat%mesprefsbac(:,:,:) = 0.0
    wombat%mesprefaoa(:,:,:) = 0.0
    wombat%mesprefphy(:,:,:) = 0.0
    wombat%mesprefdia(:,:,:) = 0.0
    wombat%mesprefsdet(:,:,:) = 0.0
    wombat%mesprefldet(:,:,:) = 0.0
    wombat%mesprefzoo(:,:,:) = 0.0
    wombat%mesgrazlbac(:,:,:) = 0.0
    wombat%mesgrazsbac(:,:,:) = 0.0
    wombat%mesgrazaoa(:,:,:) = 0.0
    wombat%mesgrazphy(:,:,:) = 0.0
    wombat%mesgrazdia(:,:,:) = 0.0
    wombat%mesgrazsdet(:,:,:) = 0.0
    wombat%mesgrazldet(:,:,:) = 0.0
    wombat%mesgrazzoo(:,:,:) = 0.0
    wombat%mesmorl(:,:,:) = 0.0
    wombat%mesmorq(:,:,:) = 0.0
    wombat%mesexcrlbac(:,:,:) = 0.0
    wombat%mesexcrsbac(:,:,:) = 0.0
    wombat%mesexcraoa(:,:,:) = 0.0
    wombat%mesexcrphy(:,:,:) = 0.0
    wombat%mesexcrdia(:,:,:) = 0.0
    wombat%mesexcrsdet(:,:,:) = 0.0
    wombat%mesexcrldet(:,:,:) = 0.0
    wombat%mesexcrzoo(:,:,:) = 0.0
    wombat%mesegeslbac(:,:,:) = 0.0
    wombat%mesegessbac(:,:,:) = 0.0
    wombat%mesegesaoa(:,:,:) = 0.0
    wombat%mesegesphy(:,:,:) = 0.0
    wombat%mesegesdia(:,:,:) = 0.0
    wombat%mesegessdet(:,:,:) = 0.0
    wombat%mesegesldet(:,:,:) = 0.0
    wombat%mesegeszoo(:,:,:) = 0.0
    wombat%reminr(:,:,:) = 0.0
    wombat%ldetremi(:,:,:) = 0.0
    wombat%sdetremi(:,:,:) = 0.0
    wombat%ldocremi(:,:,:) = 0.0
    wombat%sdocprod(:,:,:) = 0.0
    wombat%sdocremi(:,:,:) = 0.0
    wombat%pic2poc(:,:,:) = 0.0
    wombat%dissratcal(:,:,:) = 0.0
    wombat%dissratara(:,:,:) = 0.0
    wombat%dissratpoc(:,:,:) = 0.0
    wombat%caldiss(:,:,:) = 0.0
    wombat%aradiss(:,:,:) = 0.0
    wombat%pocdiss(:,:,:) = 0.0
    wombat%zoodiss(:,:,:) = 0.0
    wombat%mesdiss(:,:,:) = 0.0
    wombat%aoa_loxy(:,:,:) = 0.0
    wombat%aoa_lnh4(:,:,:) = 0.0
    wombat%aoa_eno3(:,:,:) = 0.0
    wombat%aoa_mumax(:,:,:) = 0.0
    wombat%aoa_mu(:,:,:) = 0.0
    wombat%aoagrow(:,:,:) = 0.0
    wombat%aoaresp(:,:,:) = 0.0
    wombat%aoamorl(:,:,:) = 0.0
    wombat%aoamorq(:,:,:) = 0.0
    wombat%lbacydoc(:,:,:) = 1.0
    wombat%lbacgrow(:,:,:) = 0.0
    wombat%lbacresp(:,:,:) = 0.0
    wombat%lbacpnh4(:,:,:) = 0.0
    wombat%lbacpco2(:,:,:) = 0.0
    wombat%lbacufer(:,:,:) = 0.0
    wombat%lbac_mu(:,:,:) = 0.0
    wombat%lbac_fanaer(:,:,:) = 0.0
    wombat%lbac_ffelim(:,:,:) = 0.0
    wombat%lbacmorl(:,:,:) = 0.0
    wombat%lbacmorq(:,:,:) = 0.0
    wombat%lbacdeni(:,:,:) = 0.0
    wombat%sbacydoc(:,:,:) = 1.0
    wombat%sbacgrow(:,:,:) = 0.0
    wombat%sbacresp(:,:,:) = 0.0
    wombat%sbacpnh4(:,:,:) = 0.0
    wombat%sbacpco2(:,:,:) = 0.0
    wombat%sbacufer(:,:,:) = 0.0
    wombat%sbac_mu(:,:,:) = 0.0
    wombat%sbac_fanaer(:,:,:) = 0.0
    wombat%sbac_ffelim(:,:,:) = 0.0
    wombat%sbacmorl(:,:,:) = 0.0
    wombat%sbacmorq(:,:,:) = 0.0
    wombat%sbacdeni(:,:,:) = 0.0
    wombat%aox_lnh4(:,:,:) = 0.0
    wombat%aox_mu(:,:,:) = 0.0
    wombat%ammox(:,:,:) = 0.0
    wombat%anammox(:,:,:) = 0.0
    wombat%zeuphot(:,:) = 0.0
    wombat%sdet_radius(:,:) = 0.0
    wombat%ldet_radius(:,:) = 0.0
    wombat%sdet_density(:,:,:) = 0.0
    wombat%ldet_density(:,:,:) = 0.0
    wombat%fbury(:,:) = 0.0
    wombat%seddep(:,:) = 0.0
    wombat%sedmask(:,:) = 0.0
    wombat%sedtemp(:,:) = 0.0
    wombat%sedsalt(:,:) = 0.0
    wombat%sedno3(:,:) = 0.0
    wombat%sednh4(:,:) = 0.0
    wombat%sedsil(:,:) = 0.0
    wombat%sedo2(:,:) = 0.0
    wombat%seddic(:,:) = 0.0
    wombat%sedalk(:,:) = 0.0
    wombat%sedhtotal(:,:) = 0.0
    wombat%sedco3(:,:) = 0.0
    wombat%sedomega_cal(:,:) = 0.0
    if (do_nitrogen_fixation) then
      wombat%trimumax(:,:,:) = 0.0
      wombat%tri_lpar(:,:,:) = 0.0
      wombat%tri_lfer(:,:,:) = 0.0
      wombat%nitrfix(:,:,:) = 0.0
    endif

    ! Allocate and initialise some multi-dimensional variables
    allocate(wsink1(nk)); wsink1(:)=0.0
    allocate(wsink2(nk)); wsink2(:)=0.0
    allocate(ek_bgr(nk,3)); ek_bgr(:,:)=0.0
    allocate(par_bgr_mid(nk,3)); par_bgr_mid(:,:)=0.0
    allocate(par_bgr_top(nk,3)); par_bgr_top(:,:)=0.0
    allocate(n_pools(isc:iec,jsc:jec,nk,2)); n_pools(:,:,:,:)=0.0
    allocate(c_pools(isc:iec,jsc:jec,nk,2)); c_pools(:,:,:,:)=0.0
    allocate(si_pools(isc:iec,jsc:jec,nk,2)); si_pools(:,:,:,:)=0.0
    allocate(fe_pools(isc:iec,jsc:jec,nk,2)); fe_pools(:,:,:,:)=0.0

    ! Set the maximum index for euphotic depth
    ! dts: in WOMBAT v3, kmeuph and k100 are integers but here they are arrays since zw
    ! may vary spatially
    allocate(kmeuph(isc:iec, jsc:jec)); kmeuph(:,:)=1
    allocate(k100(isc:iec, jsc:jec)); k100(:,:)=1
    do j = jsc,jec; do i = isc,iec;
      nz = grid_kmt(i,j)
      do k = 1,nz
        if (wombat%zw(i,j,k) <= 400) kmeuph(i,j)=k
        if (wombat%zw(i,j,k) <= 100) k100(i,j)=k
      enddo
    enddo; enddo

    ! Get the timestep for the ecosystem model
    ts_npzd = max(1, nint(dt / wombat%dt_npzd)) ! number of ecosystem timesteps per model timestep
    rdtts = 1 / dt
    dtsb = dt / float(ts_npzd) ! number of seconds per nested ecosystem timestep

    ! Get the prognostic tracer values
    ! dts attn: we should probably update prognostic tracers via pointers to avoid
    ! having to allocate all these field arrays
    ! dts attn: do we really want/need to force these to be positive?
    call g_tracer_get_values(tracer_list, 'no3', 'field', wombat%f_no3, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'nh4', 'field', wombat%f_nh4, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'sil', 'field', wombat%f_sil, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'phy', 'field', wombat%f_phy, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'pchl', 'field', wombat%f_pchl, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'phyfe', 'field', wombat%f_phyfe, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'dia', 'field', wombat%f_dia, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'dchl', 'field', wombat%f_dchl, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'diafe', 'field', wombat%f_diafe, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'diasi', 'field', wombat%f_diasi, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'zoo', 'field', wombat%f_zoo, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'zoofe', 'field', wombat%f_zoofe, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'mes', 'field', wombat%f_mes, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'mesfe', 'field', wombat%f_mesfe, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'sdet', 'field', wombat%f_sdet, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'sdetfe', 'field', wombat%f_sdetfe, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'ldet', 'field', wombat%f_ldet, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'ldetfe', 'field', wombat%f_ldetfe, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'ldetsi', 'field', wombat%f_ldetsi, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'sdoc', 'field', wombat%f_sdoc, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'ldoc', 'field', wombat%f_ldoc, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'ldon', 'field', wombat%f_ldon, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'lbac', 'field', wombat%f_lbac, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'sbac', 'field', wombat%f_sbac, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'aoa', 'field', wombat%f_aoa, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'o2', 'field', wombat%f_o2, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'caco3', 'field', wombat%f_caco3, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'fe', 'field', wombat%f_fe, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'safe', 'field', wombat%f_safe, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'lafe', 'field', wombat%f_lafe, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'dic', 'field', wombat%f_dic, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'alk', 'field', wombat%f_alk, isd, jsd, ntau=tau) ! [mol/kg]
    if (do_tracer_dicr) call g_tracer_get_values(tracer_list, 'dicr', 'field', wombat%f_dicr, isd, jsd, ntau=tau) ! [mol/kg]


    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !     (\___/)  .-.   .-. .--. .-..-..---.  .--. .-----.                 !
    !     / o o \  : :.-.: :: ,. :: `' :: .; :: .; :`-. .-'                 !
    !    (   "   ) : :: :: :: :: :: .. ::   .':    :  : :                   !
    !     \__ __/  : `' `' ;: :; :: :; :: .; :: :: :  : :                   !
    !               `.,`.,' `.__.':_;:_;:___.':_;:_;  :_;                   !
    !                                                                       !
    ! World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT)    !
    !                                                                       !
    !  Steps:                                                               !
    !    1.  Light attenuation through water column                         !
    !    2.  Nutrient limitation of phytoplankton                           !
    !    3.  Temperature-dependent heterotrophy and POM-->DOM               !
    !    4.  Light limitation of phytoplankton                              !
    !    5.  Realized growth rate of phytoplankton                          !
    !    6.  Dissolved organic carbon release by phytoplankton              !
    !    7.  Synthesis of chlorophyll                                       !
    !    8.  Phytoplankton uptake of iron                                   !
    !    9.  Phytoplankton uptake of silicic acid                           !
    !    10. Iron chemistry (scavenging, coagulation, dissolution)          !
    !    11. Biogenic silica dissolution                                    !
    !    12. Mortality terms                                                !
    !    13. Zooplankton grazing, egestion, excretion, assimilation         !
    !    14. Implicit nitrogen fixation                                     !
    !    15. Facultative bacterial heterotrophy                             !
    !    16. Calcium carbonate production and dissolution                   !
    !    17. Chemoautotrophy                                                !
    !    18. Tracer tendencies                                              !
    !    19. Check for conservation of mass                                 !
    !    20. Additional operations on tracers                               !
    !    21. Sinking rates of particulates                                  !
    !    22. Sedimentary processes                                          !
    !                                                                       !
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!


    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !  [Step 1] Light attenutation through water column                     !
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    do j = jsc,jec; do i = isc,iec
      !--- Retrieve incident light ---!
      ! dts: keep only shortwave wavelengths < 710 nm (taken from COBALT/BLING)
      swpar = 0.0
      do n = 1,nbands;
        if (max_wavelength_band(n) < 710) then
          sw_pen(n) = max(0.0, sw_pen_band(n,i,j))
        else
          sw_pen(n) = 0.0
        endif
        swpar = swpar + sw_pen(n)
      enddo

      !--- Light field ---!

      par_bgr_top(:,:) = 0.0
      par_bgr_mid(:,:) = 0.0

      do k = 1,grid_kmt(i,j)  !{

        ! chlorophyll concentration conversion from mol/kg --> mg/m3 for look-up table
        chl = max(0.0, (wombat%f_pchl(i,j,k) + wombat%f_dchl(i,j,k))) * 12.0 / mmol_m3_to_mol_kg
        ! detritus concentration conversion from mol/kg --> mgN/m3 for look-up table
        ndet = max(0.0, (wombat%f_sdet(i,j,k) + wombat%f_ldet(i,j,k))) * 16.0/122.0 * 14.0 / mmol_m3_to_mol_kg
        ! CaCO3 concentration conversion from mol/kg --> kg/m3 for look-up table
        carb = max(0.0,wombat%f_caco3(i,j,k)) / mmol_m3_to_mol_kg * 100.09 * 1e-3 * 1e-3 ! convert to kg/m3

        ! Attenuation coefficients given chlorophyll concentration
        zchl = max(0.01, min(10.0, chl) )
        ichl = nint( 41 + 20.0*log10(zchl) + epsi )
        ek_bgr(k,1) = (zbgr(2,ichl) + ndet * dbgr(1) + carb * cbgr(1)) * dzt(i,j,k) ! [/m * m]
        ek_bgr(k,2) = (zbgr(3,ichl) + ndet * dbgr(2) + carb * cbgr(2)) * dzt(i,j,k) ! [/m * m]
        ek_bgr(k,3) = (zbgr(4,ichl) + ndet * dbgr(3) + carb * cbgr(3)) * dzt(i,j,k) ! [/m * m]

        ! BGR light available in the water column
        if (swpar>0.0) then
          if (k==1) then
            par_bgr_top(k,1) = swpar * 1.0/3.0  ! Assumes 1/3 of PAR is blue (coming through top of cell)
            par_bgr_top(k,2) = swpar * 1.0/3.0  ! Assumes 1/3 of PAR is green (coming through top of cell)
            par_bgr_top(k,3) = swpar * 1.0/3.0  ! Assumes 1/3 of PAR is red (coming through top of cell)
            par_bgr_mid(k,1) = swpar * 1.0/3.0 * exp(-0.5 * ek_bgr(k,1)) ! Assumes 1/3 of PAR is blue (at centre point)
            par_bgr_mid(k,2) = swpar * 1.0/3.0 * exp(-0.5 * ek_bgr(k,2)) ! Assumes 1/3 of PAR is green (at centre point)
            par_bgr_mid(k,3) = swpar * 1.0/3.0 * exp(-0.5 * ek_bgr(k,3)) ! Assumes 1/3 of PAR is red (at centre point)
          else
            par_bgr_top(k,1) = par_bgr_top(k-1,1) * exp(-ek_bgr(k-1,1))
            par_bgr_top(k,2) = par_bgr_top(k-1,2) * exp(-ek_bgr(k-1,2))
            par_bgr_top(k,3) = par_bgr_top(k-1,3) * exp(-ek_bgr(k-1,3))
            par_bgr_mid(k,1) = par_bgr_mid(k-1,1) * exp(-0.5 * (ek_bgr(k-1,1)+ek_bgr(k,1)))
            par_bgr_mid(k,2) = par_bgr_mid(k-1,2) * exp(-0.5 * (ek_bgr(k-1,2)+ek_bgr(k,2)))
            par_bgr_mid(k,3) = par_bgr_mid(k-1,3) * exp(-0.5 * (ek_bgr(k-1,3)+ek_bgr(k,3)))
          endif
        endif

        ! Light attenuation at mid point of the cells
        wombat%radmid(i,j,k) = par_bgr_mid(k,1) + par_bgr_mid(k,2) + par_bgr_mid(k,3)

      enddo !} k

      ! initialise some variables
      par_phy_mldsum = 0.0
      par_z_mldsum = 0.0
      wombat%zeuphot(i,j) = 0.0  ! Euphotic zone depth set at 0.0 meters when there is no light

      !--- Collect euphotic zone depth, light at mid points, and mean light ---!
      do k = 1,grid_kmt(i,j)  !{

        if (swpar>0.0) then
          ! Euphotic zone
          if (wombat%radmid(i,j,k)>(swpar*0.01) .and. wombat%radmid(i,j,k)>0.01) then
            wombat%zeuphot(i,j) = wombat%zw(i,j,k)
          endif
          ! Light attenuation mean over the grid cells (Eq. 19 of Baird et al., 2020 GMD)
          if (k<grid_kmt(i,j)) then
            wombat%radbio(i,j,k) = ((par_bgr_top(k,1) - par_bgr_top(k+1,1)) / ek_bgr(k,1)) + &
                                   ((par_bgr_top(k,2) - par_bgr_top(k+1,2)) / ek_bgr(k,2)) + &
                                   ((par_bgr_top(k,3) - par_bgr_top(k+1,3)) / ek_bgr(k,3))
          else
            wombat%radbio(i,j,k) = ((par_bgr_top(k,1) - par_bgr_top(k,1)*exp(-ek_bgr(k,1))) / ek_bgr(k,1)) + &
                                   ((par_bgr_top(k,2) - par_bgr_top(k,2)*exp(-ek_bgr(k,2))) / ek_bgr(k,2)) + &
                                   ((par_bgr_top(k,3) - par_bgr_top(k,3)*exp(-ek_bgr(k,3))) / ek_bgr(k,3))
          endif
        else
          wombat%radbio(i,j,k) = 0.0
        endif

        ! Integrated light field in the mixed layer
        if (wombat%zw(i,j,k)<=hblt_depth(i,j)) then
          par_phy_mldsum = par_phy_mldsum + wombat%radbio(i,j,k) * dzt(i,j,k)
          par_z_mldsum = par_z_mldsum + dzt(i,j,k)
        endif

      enddo !}

      !--- Aggregate light in mixed layer and calculate maximum growth rates ---!
      do k = 1,grid_kmt(i,j)  !{

        ! Calculate average light level in the mixed layer
        if (wombat%zw(i,j,k)<=hblt_depth(i,j)) then
          zval = 1.0/(par_z_mldsum + epsi)
          wombat%radmld(i,j,k) = par_phy_mldsum * zval
        else
          wombat%radmld(i,j,k) = wombat%radbio(i,j,k)
        endif

        ! Temperature-dependent maximum growth rate (Eppley curve)
        wombat%phy_mumax(i,j,k) = wombat%abioa_phy * wombat%bbioa_phy ** (Temp(i,j,k)) ! [1/s]
        wombat%dia_mumax(i,j,k) = wombat%abioa_dia * wombat%bbioa_dia ** (Temp(i,j,k)) ! [1/s]

      enddo  !} k

    enddo; enddo

    !! pjb - these are no longer used but it is a useful place to save this information
    !wombat%no3_prev(:,:,:) = wombat%f_no3(:,:,:)
    !wombat%nh4_prev(:,:,:) = wombat%f_nh4(:,:,:)
    !wombat%caco3_prev(:,:,:) = wombat%f_caco3(:,:,:)

    ! Arrays for assessing conservation of mass within ecosystem component
    n_pools(:,:,:,:) = 0.0
    c_pools(:,:,:,:) = 0.0
    si_pools(:,:,:,:) = 0.0
    fe_pools(:,:,:,:) = 0.0

    do tn = 1,ts_npzd  !{

      n_pools(:,:,:,1) = n_pools(:,:,:,2)
      c_pools(:,:,:,1) = c_pools(:,:,:,2)
      si_pools(:,:,:,1) = si_pools(:,:,:,2)
      fe_pools(:,:,:,1) = fe_pools(:,:,:,2)
      do k = 1,nk; do j = jsc,jec; do i = isc,iec;

      ! Initialise some floored values and ratios (put into nicer units than mol/kg)
      ! Loss/source coefficients are calculated using the floored values but are applied to the
      ! un-floored tracer values
      phy_p      = max(0.0, wombat%f_phy(i,j,k) )
      dia_p      = max(0.0, wombat%f_dia(i,j,k) )
      pchl_p     = max(0.0, wombat%f_pchl(i,j,k) )
      dchl_p     = max(0.0, wombat%f_dchl(i,j,k) )
      phyfe_p    = max(0.0, wombat%f_phyfe(i,j,k))
      diafe_p    = max(0.0, wombat%f_diafe(i,j,k))
      diasi_p    = max(0.0, wombat%f_diasi(i,j,k))
      zoo_p      = max(0.0, wombat%f_zoo(i,j,k) )
      mes_p      = max(0.0, wombat%f_mes(i,j,k) )
      sdet_p     = max(0.0, wombat%f_sdet(i,j,k) )
      ldet_p     = max(0.0, wombat%f_ldet(i,j,k) )
      ldetsi_p   = max(0.0, wombat%f_ldetsi(i,j,k) )
      ldoc_p     = max(0.0, wombat%f_ldoc(i,j,k) )
      ldon_p     = max(0.0, wombat%f_ldon(i,j,k) )
      sdoc_p     = max(0.0, wombat%f_sdoc(i,j,k) )
      lbac_p     = max(0.0, wombat%f_lbac(i,j,k) )
      sbac_p     = max(0.0, wombat%f_sbac(i,j,k) )
      aoa_p      = max(0.0, wombat%f_aoa(i,j,k) )
      no3_p      = max(0.0, wombat%f_no3(i,j,k) )
      nh4_p      = max(0.0, wombat%f_nh4(i,j,k) )
      oxy_p      = max(0.0, wombat%f_o2(i,j,k) )
      fe_p       = max(0.0, wombat%f_fe(i,j,k)  )
      sil_p      = max(0.0, wombat%f_sil(i,j,k) )
      caco3_p    = max(0.0, wombat%f_caco3(i,j,k))
      safe_p     = max(0.0, wombat%f_safe(i,j,k) )
      lafe_p     = max(0.0, wombat%f_lafe(i,j,k) )

      phy_mmolm3 = phy_p / mmol_m3_to_mol_kg  ![mmol/m3]
      dia_mmolm3 = dia_p / mmol_m3_to_mol_kg  ![mmol/m3]
      phyfe_mmolm3 = phyfe_p / mmol_m3_to_mol_kg  ![mmol/m3]
      diafe_mmolm3 = diafe_p / mmol_m3_to_mol_kg  ![mmol/m3]
      zoo_mmolm3 = zoo_p / mmol_m3_to_mol_kg  ![mmol/m3]
      mes_mmolm3 = mes_p / mmol_m3_to_mol_kg  ![mmol/m3]
      sdet_mmolm3 = sdet_p / mmol_m3_to_mol_kg  ![mmol/m3]
      ldet_mmolm3 = ldet_p / mmol_m3_to_mol_kg  ![mmol/m3]
      ldetsi_mmolm3 = ldetsi_p / mmol_m3_to_mol_kg  ![mmol/m3]
      ldoc_mmolm3 = ldoc_p / mmol_m3_to_mol_kg  ![mmol/m3]
      ldon_mmolm3 = ldon_p / mmol_m3_to_mol_kg  ![mmol/m3]
      sdoc_mmolm3 = sdoc_p / mmol_m3_to_mol_kg  ![mmol/m3]
      lbac_mmolm3 = lbac_p / mmol_m3_to_mol_kg  ![mmol/m3]
      sbac_mmolm3 = sbac_p / mmol_m3_to_mol_kg  ![mmol/m3]
      aoa_mmolm3 = aoa_p / mmol_m3_to_mol_kg  ![mmol/m3]
      no3_mmolm3 = no3_p / mmol_m3_to_mol_kg  ![mmol/m3]
      nh4_mmolm3 = nh4_p / mmol_m3_to_mol_kg  ![mmol/m3]
      oxy_mmolm3 = oxy_p / mmol_m3_to_mol_kg  ![mmol/m3]
      fe_umolm3 = fe_p / umol_m3_to_mol_kg  ![umol/m3]
      sil_mmolm3 = sil_p / mmol_m3_to_mol_kg  ![mmol/m3]
      caco3_mmolm3 = caco3_p / mmol_m3_to_mol_kg  ![mmol/m3]

      ! Initialise ratios
      phy_chlc = 0.0
      dia_chlc = 0.0
      phy_Fe2C = 0.0
      dia_Fe2C = 0.0
      zoo_Fe2C = 0.0
      mes_Fe2C = 0.0
      sdet_Fe2C = 0.0
      ldet_Fe2C = 0.0
      ldom_N2C = 0.0
      sdom_N2C = 0.0
      dia_Si2C = 0.0

      ! Only calculate these ratios when the denominator is > 0.0
      if (phy_p > 0.0) then
        phy_chlc = pchl_p / phy_p
        phy_Fe2C = phyfe_p / phy_p
      endif
      if (dia_p > 0.0) then
        dia_chlc = dchl_p / dia_p
        dia_Fe2C = diafe_p / dia_p
        dia_Si2C = diasi_p / dia_p
      endif
      if (zoo_p > 0.0) then
        zoo_Fe2C = zoofe_p / zoo_p
      end if
      if (mes_p > 0.0) then
        mes_Fe2C = mesfe_p / mes_p
      end if
      if (sdet_p > 0.0) then
        sdet_Fe2C = sdetfe_p / sdet_p
      end if
      if (ldet_p > 0.0) then
        ldet_Fe2C = ldetfe_p / ldet_p
      end if
      if (ldoc_p > 0.0) then
        ldom_N2C  = ldon_p / ldoc_p
      end if
      !if (sdoc_p > 0.0) then
      !  sdom_N2C  = sdon_p / sdoc_p
      !end if


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 2] Nutrient limitation of phytoplankton                        !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! 1. Allometric scaling of 0.37 (Wickman et al., 2024; Science)
      ! 2. Nitrogen limitation, with NH4 preferred for supporting growth
      ! 3. Iron limiation, with distinct costs associated with chlorophyll, respiration and nitrate reduction

      !!!~~~ Phytoplankton ~~~!!!
      wombat%phy_kni(i,j,k) = wombat%phykn * max(0.1, max(0.0, (phy_mmolm3-wombat%phybiot))**0.37)
      wombat%phy_kfe(i,j,k) = wombat%phykf * max(0.1, max(0.0, (phy_mmolm3-wombat%phybiot))**0.37)
      ! Nitrogen limitation (split into NH4 and NO3 uptake limitation terms, with NH4 uptake being 5x more preferred)
      !   See Buchanan et al., (2025) Biogeosciences
      phy_limnh4 = nh4_mmolm3 / (nh4_mmolm3 + wombat%phy_kni(i,j,k) + epsi)
      phy_limno3 = no3_mmolm3 / (no3_mmolm3 + wombat%phy_kni(i,j,k) + epsi)
      phy_limdin = (no3_mmolm3 + nh4_mmolm3) / (no3_mmolm3 + nh4_mmolm3 + wombat%phy_kni(i,j,k) + epsi)
      wombat%phy_lnh4(i,j,k) = wombat%phyprefnh4 * phy_limdin * phy_limnh4 / (phy_limno3 + wombat%phyprefnh4 * phy_limnh4 + epsi)
      wombat%phy_lno3(i,j,k) = phy_limdin * phy_limno3 / (phy_limno3 + wombat%phyprefnh4 * phy_limnh4 + epsi)
      wombat%phy_lnit(i,j,k) = wombat%phy_lno3(i,j,k) + wombat%phy_lnh4(i,j,k)
      ! Iron limitation (Quota model, constants from Flynn & Hipkin 1999)
      phy_minqfe = 0.00167 / 55.85 * max(wombat%phyminqc, phy_chlc)*12 + &
                   1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * wombat%phy_lnit(i,j,k) + &
                   1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * wombat%phy_lno3(i,j,k)
      wombat%phy_lfer(i,j,k) = min(1.0, max(0.0, (phy_Fe2C - phy_minqfe) / wombat%phyoptqf ))

      !!!~~~ Microphytoplankton ~~~!!!
      wombat%dia_kni(i,j,k) = wombat%diakn * max(0.1, max(0.0, (dia_mmolm3-wombat%diabiot))**0.37)
      wombat%dia_kfe(i,j,k) = wombat%diakf * max(0.1, max(0.0, (dia_mmolm3-wombat%diabiot))**0.37)
      wombat%dia_ksi(i,j,k) = wombat%diaks * max(0.1, max(0.0, (dia_mmolm3-wombat%diabiot))**0.37)
      ! Nitrogen limitation (split into NH4 and NO3 uptake limitation terms, with NH4 uptake being 5x more preferred)
      !   See Buchanan et al., (2025) Biogeosciences
      dia_limnh4 = nh4_mmolm3 / (nh4_mmolm3 + wombat%dia_kni(i,j,k) + epsi)
      dia_limno3 = no3_mmolm3 / (no3_mmolm3 + wombat%dia_kni(i,j,k) + epsi)
      dia_limdin = (no3_mmolm3 + nh4_mmolm3) / (no3_mmolm3 + nh4_mmolm3 + wombat%dia_kni(i,j,k) + epsi)
      wombat%dia_lnh4(i,j,k) = wombat%diaprefnh4 * dia_limdin * dia_limnh4 / (dia_limno3 + wombat%diaprefnh4 * dia_limnh4 + epsi)
      wombat%dia_lno3(i,j,k) = dia_limdin * dia_limno3 / (dia_limno3 + wombat%diaprefnh4 * dia_limnh4 + epsi)
      wombat%dia_lnit(i,j,k) = wombat%dia_lno3(i,j,k) + wombat%dia_lnh4(i,j,k)
      ! Iron limitation (Quota model, constants from Flynn & Hipkin 1999)
      dia_minqfe = 0.00167 / 55.85 * max(wombat%diaminqc, dia_chlc)*12 + &
                   1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * wombat%dia_lnit(i,j,k) + &
                   1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * wombat%dia_lno3(i,j,k)
      wombat%dia_lfer(i,j,k) = min(1.0, max(0.0, (dia_Fe2C - dia_minqfe) / wombat%diaoptqf ))
      ! Silicic acid limitation (gating constraint on division)
      wombat%dia_lsil(i,j,k) = min(1.0, max(0.0, (dia_Si2C - wombat%diaminqs) / (wombat%diaoptqs - wombat%diaminqs) ))


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 3] Temperature-dependent heterotrophy and POM-->DOM            !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Temperature dependance of heterotrophy (applies to bact and zoo)
      fbc = wombat%bbioh ** (Temp(i,j,k))

      ! Variable rates of remineralisation
      wombat%reminr(i,j,k) = wombat%sdetlrem * fbc

      ! remineralisation of POM --> lDOM
      wombat%sdetremi(i,j,k) = wombat%reminr(i,j,k) / mmol_m3_to_mol_kg * sdet_p * sdet_p ! [molC/kg/s]
      wombat%ldetremi(i,j,k) = wombat%reminr(i,j,k) / mmol_m3_to_mol_kg * ldet_p * ldet_p ! [molC/kg/s]


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 4] Light limitation of phytoplankton                           !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! 1. initial slope of Photosynthesis-Irradiance curve
      ! 2. Light limitation

      !!!~~~ Phytoplankton ~~~!!!
      phy_pisl  = max(phy_chlc, wombat%phyminqc) * wombat%alphabio_phy
      wombat%phy_lpar(i,j,k) = (1. - exp(-phy_pisl * wombat%radbio(i,j,k)))

      !!!~~~ Microphytoplankton ~~~!!!
      dia_pisl  = max(dia_chlc, wombat%diaminqc) * wombat%alphabio_dia
      wombat%dia_lpar(i,j,k) = (1. - exp(-dia_pisl * wombat%radbio(i,j,k)))


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 5] Realized growth rate of phytoplankton                       !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! 1. Apply light and nutrient limitations to maximum growth rate
      wombat%phy_mu(i,j,k) = wombat%phy_mumax(i,j,k) * wombat%phy_lpar(i,j,k) * &
                             min(wombat%phy_lnit(i,j,k), wombat%phy_lfer(i,j,k))
      wombat%dia_mu(i,j,k) = wombat%dia_mumax(i,j,k) * wombat%dia_lpar(i,j,k) * &
                             min(wombat%dia_lnit(i,j,k), wombat%dia_lfer(i,j,k)) * &
                             wombat%dia_lsil(i,j,k)

      if (nh4_p + no3_p > epsi) then
        wombat%phygrow(i,j,k) = wombat%phy_mu(i,j,k) * phy_p ! [molC/kg/s]
        wombat%diagrow(i,j,k) = wombat%dia_mu(i,j,k) * dia_p ! [molC/kg/s]
      else
        wombat%phygrow(i,j,k) = 0.0
        wombat%diagrow(i,j,k) = 0.0
      endif


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 6] Dissolved organic carbon release by phytoplankton           !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Excess DOC exudation (active exudation via overflow hypothesis; Fogg 1966, 1983; Williams 1990; Carlson & Hansell 2014)
        ! Up to 50% (set by `overflow`) of assimilated carbon can be exuded by phytoplankton as DOC in high light, low nutrient
        ! conditions (Thornton 2014)
        ! Some small amount of DOC is exuded via passive diffusion even in the healthiest phytoplankton (Bjornsen 1988)
        ! If too much DOC is exuded, bacterial competition for nutrients can limit phytoplankton growth (Bratbak & Thingstad, 1985;
        ! Ratnarajah et al. 2021)
        ! However, active release of DOM by mixotrophic phytoplankton can "farm" heterotrophic bacteria (Mitra et al. 2013) (NOT
        ! YET IMPLEMENTED)
      if (phy_p > epsi) then
        zval = wombat%phy_mumax(i,j,k) * wombat%phy_lpar(i,j,k) * phy_p ! Gross carbon fixation
        wombat%phydoc(i,j,k) = min( wombat%overflow * zval, &
                                    max( 0.02 * zval, max(0.0, zval - wombat%phygrow(i,j,k)) ) ) !
      else
        wombat%phydoc(i,j,k) = 0.0
      endif
      if (dia_p > epsi) then
        zval = wombat%dia_mumax(i,j,k) * wombat%dia_lpar(i,j,k) * dia_p ! Gross carbon fixation
        wombat%diadoc(i,j,k) = min( wombat%overflow * zval, &
                                    max( 0.02 * zval, max(0.0, zval - wombat%diagrow(i,j,k)) ) ) !
      else
        wombat%diadoc(i,j,k) = 0.0
      endif


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 7] Synthesis of chlorophyll                                    !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! 1. Estimate the target Chl:C ratio required to support maximum growth
      !    This is a direct approximation of Geider, MacIntyre & Kana (1997):
      !     - Chl increased in response to low light (we use the mean of the MLD)
      !     - Chl decreased in response to N and Fe limitation
      ! 2. Ensure that a minimum chl:c ratio must be maintained by the cell
      ! 3. Estimate the rate that chlorophyll is synthesized towards this optimal
      !    Rates of chlorophyll synthesis are not instantaneous, and take hours to days
      !    Here, we make chlorophyll synthesis respond on a timescale of chltau

      !!!~~~ Nanophytoplankton ~~~!!!
      theta_opt = wombat%phymaxqc / (1.0 + &
                  ( wombat%alphabio_phy * wombat%radmld(i,j,k) * wombat%phymaxqc ) &
                 /( epsi + 2.0 * wombat%phy_mumax(i,j,k) * 86400.0 &
                    * max(0.01, min(wombat%phy_lnit(i,j,k), wombat%phy_lfer(i,j,k))) ) )
      theta_opt = max(wombat%phyminqc, theta_opt)
      wombat%pchl_mu(i,j,k) = wombat%phy_mu(i,j,k) * pchl_p &
                            + (theta_opt - phy_chlc) / wombat%chltau * phy_p

      !!!~~~ Microphytoplankton ~~~!!!
      theta_opt = wombat%diamaxqc / (1.0 + &
                  ( wombat%alphabio_dia * wombat%radmld(i,j,k) * wombat%diamaxqc ) &
                 /( epsi + 2.0 * wombat%dia_mumax(i,j,k) * 86400.0 &
                    * max(0.01, min(wombat%dia_lnit(i,j,k), wombat%dia_lfer(i,j,k))) ) )
      theta_opt = max(wombat%diaminqc, theta_opt)
      wombat%dchl_mu(i,j,k) = wombat%dia_mu(i,j,k) * dchl_p &
                            + (theta_opt - dia_chlc) / wombat%chltau * dia_p


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 8] Phytoplankton uptake of iron                                !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! 1. Maximum iron content of phytoplankton cell
      ! 2. Ensure that dFe uptake increases or decreases in response to cell quota
      ! 3. Iron uptake of phytoplankton (reduced 10-fold in darkness)

      !!!~~~ Nano-phytoplankton ~~~!!!
      if (phy_p > epsi) then
        phy_maxqfe = phy_mmolm3 * wombat%phymaxqf  !mmol Fe / m3
        wombat%phy_feupreg(i,j,k) = (4.0 - 4.5 * wombat%phy_lfer(i,j,k) / &
                                    (wombat%phy_lfer(i,j,k) + 0.5) )
        wombat%phy_fedoreg(i,j,k) = max(0.0, (1.0 - phyfe_mmolm3/phy_maxqfe) / &
                                    abs(1.05 - phyfe_mmolm3/phy_maxqfe) )
        wombat%phy_dfeupt(i,j,k) = (wombat%phy_mumax(i,j,k) * phy_maxqfe * &
                                    max(0.01, wombat%phy_lpar(i,j,k))**0.5 * &
                                    fe_umolm3 / (fe_umolm3 + wombat%phy_kfe(i,j,k)) * &
                                    wombat%phy_feupreg(i,j,k) * &
                                    wombat%phy_fedoreg(i,j,k) ) * mmol_m3_to_mol_kg
      endif

      !!!~~~ Microphytoplankton ~~~!!!
      if (dia_p > epsi) then
        dia_maxqfe = dia_mmolm3 * wombat%diamaxqf  !mmol Fe / m3
        wombat%dia_feupreg(i,j,k) = (4.0 - 4.5 * wombat%dia_lfer(i,j,k) / &
                                    (wombat%dia_lfer(i,j,k) + 0.5) )
        wombat%dia_fedoreg(i,j,k) = max(0.0, (1.0 - diafe_mmolm3/dia_maxqfe) / &
                                    abs(1.05 - diafe_mmolm3/dia_maxqfe) )
        wombat%dia_dfeupt(i,j,k) = max(0.0, ( wombat%dia_mumax(i,j,k) * dia_maxqfe &
                                       * max(0.01, wombat%dia_lpar(i,j,k))**0.5 &
                                       * fe_umolm3 / (fe_umolm3 + wombat%dia_kfe(i,j,k)) &
                                       * wombat%dia_feupreg(i,j,k) &
                                       * wombat%dia_fedoreg(i,j,k))) * mmol_m3_to_mol_kg
      endif


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 9] Phytoplankton uptake of silicic acid                        !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! 1. Ensure that silicic acid uptake can decrease in response to high cell quotas
      ! 2. We purposefully do not up-regulate silicic acid uptake because the mechanism
      !    of highly silicified diatoms is associated with slow growth (and vice versa)
      !    and therefore less dilution of the cell quota, not increased silicic acid uptake
      !    - light limitation of growth increases Si:C by 3-fold [Liu et al., 2016 Frontiers in Marine Science]
      !    - iron limitation of growth increases Si:C by 3-fold [Hutchins & Bruland 1998; Takeda 1998]
      !    These studies report such effects because growth rates slow, but Si uptake does not change
      !    Also, note we use a T-independent maximum uptake rate for silicic acid "diaVmaxs",
      !    which ensures that slower growing diatoms in the cold polar regions take up more Si
      !    per C than faster growing diatoms in the warm tropics [Baines et al., 2010 GBC]

      !!!~~~ Microphytoplankton ~~~!!!
      !  NOTE: "zval" below is different from "dia_lsil" above due to use of "diamaxqs" instead of "diaoptqs"
      zval = min(1.0, max(0.0, (dia_Si2C - wombat%diaminqs) / (wombat%diamaxqs - wombat%diaminqs) ))
      wombat%dia_sidoreg(i,j,k) = max(0.0, (1.0 - zval)**0.5)
      wombat%dia_silupt(i,j,k) = max(0.0, (wombat%diaVmaxs * dia_mmolm3 &
                                     * sil_mmolm3 / (sil_mmolm3 + wombat%dia_ksi(i,j,k)) &
                                     * wombat%dia_sidoreg(i,j,k) )) * mmol_m3_to_mol_kg ! [molSi/kg/s]


      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !  [Step 10] Iron chemistry (scavenging, coagulation, dissolution)       !
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!

      ! Estimate solubility of Fe3+ (free Fe) in solution using temperature,
      ! pH and salinity using the equations of Liu & Millero (2002)
      ztemk = max(5.0, Temp(i,j,k)) + 273.15    ! temperature in kelvin
      I_ztemk = 1.0 / ztemk
      zval = 19.924 * Salt(i,j,k) / ( 1000. - 1.005 * Salt(i,j,k))
      sqrt_zval = sqrt(zval)
      fesol1 = 10.0**(-13.486 - 0.1856*sqrt_zval + 0.3073*zval + 5254.0*I_ztemk)
      fesol2 = 10.0**(2.517 - 0.8885*sqrt_zval + 0.2139*zval - 1320.0*I_ztemk)
      fesol3 = 10.0**(0.4511 - 0.3305*sqrt_zval - 1996.0*I_ztemk)
      fesol4 = 10.0**(-0.2965 - 0.7881*sqrt_zval - 4086.0*I_ztemk)
      fesol5 = 10.0**(4.4466 - 0.8505*sqrt_zval - 7980.0*I_ztemk)
      if (wombat%htotal(i,j,k)>0.0) then
        hp = wombat%htotal(i,j,k)
      else
        hp = 1.25893e-08 ! dts: =10.0**(-7.9)
      endif
      fe3sol = fesol1 * ( hp*hp*hp + fesol2*hp*hp + fesol3*hp + fesol4 + fesol5/hp ) *1e9

      ! Estimate total colloidal iron following Tagliabue et al. (2023).
      ! Colloidal dFe is considered to be whatever exceeds the inorganic solubility
      ! ceiling, although there is always a hard lower limit of 10% of total dFe.
      if (do_colloidal_shunt) then
        wombat%fecol(i,j,k) = max(0.1 * fe_umolm3, fe_umolm3 - fe3sol)
      else
        wombat%fecol(i,j,k) = 0.0
      endif

      ! Determine equilibriuim fractionation of the remaining dFe (non-colloidal)
      ! between Fe' and ligand-bound iron (L-Fe). Below, temperature increases Keq
      ! (reducing free Fe), light decreases Keq (increasing free Fe), pH decreases Keq
      ! and DOC increases Keq. The temperature-dependency comes from Volker & Tagliabue
      ! (2015), while the light dependency is informed by Barbeau et al. (2001) who saw
      ! a 0.7 log10 unit decrease in K in high light. The pH and DOC dependency (3rd term)
      ! comes from Ye et al. (2020) and increases binding strength at lower pH and higher
      ! concentrations of DOC.
      fe_sfe = max(0.0, fe_umolm3 - wombat%fecol(i,j,k))
      wombat%ligK(i,j,k) = 1e-9 * ( 10.0**( (17.27 - 1565.7 * I_ztemk ) &
                               - 0.7 * wombat%radbio(i,j,k) / (wombat%radbio(i,j,k) + 10.0) ) &
                               + 10.0**( (-2e-4*(sdoc_mmolm3 + ldoc_mmolm3) + 0.034) &
                                        *(sdoc_mmolm3 + ldoc_mmolm3)  - 1.67*(-log10(hp)) + 24.36 ) )
      ligW_K = wombat%ligK(i,j,k) * 10.0**(-1.5) ! ligand binding constant for weak ligands (assumed to be -1.5 log10 units weaker)
      if (do_two_ligands) then
        ! Newton-Raphson to solve: x + K1*L1*x/(1+K1*x) + K2*L2*x/(1+K2*x) = fe_sfe
        ! Initial guess is the max of two limiting-case approximations:
        !   - low-Fe (ligand-unsaturated, K*x<<1): x ~ fe_sfe / (1 + K1*L1 + K2*L2)
        !   - high-Fe (ligand-saturated, K*x>>1): x ~ fe_sfe - L1 - L2
        ! Both underestimate x*, so max() picks whichever regime applies.
        fe3 = max( fe_sfe / (1.0 + wombat%ligK(i,j,k)*wombat%ligS + ligW_K*wombat%ligW + epsi), &
                    fe_sfe - wombat%ligS - wombat%ligW )
        fe3 = max(0.0, min(fe3, fe_sfe))
        do iter = 1, 8
          inv1 = 1.0 / (1.0 + wombat%ligK(i,j,k) * fe3) ! 1/(1 + K1*x)
          inv2 = 1.0 / (1.0 + ligW_K * fe3)               ! 1/(1 + K2*x)
          FeL1_mid = wombat%ligK(i,j,k) * wombat%ligS * fe3 * inv1
          FeL2_mid = ligW_K * wombat%ligW * fe3 * inv2
          zval = fe3 + FeL1_mid + FeL2_mid - fe_sfe       ! f(x)
          fe3 = fe3 - zval / &                            ! x - f(x)/f'(x)
              (1.0 + wombat%ligK(i,j,k)*wombat%ligS*inv1*inv1 + ligW_K*wombat%ligW*inv2*inv2)
          fe3 = max(0.0, min(fe3, fe_sfe))
        enddo
        wombat%feIII(i,j,k) = fe3
      else
        wombat%ligK(i,j,k) = (wombat%ligK(i,j,k)*wombat%ligS + ligW_K*wombat%ligW) / (wombat%ligS + wombat%ligW)
        zval = 1.0 + (wombat%ligS + wombat%ligW) * wombat%ligK(i,j,k) - fe_sfe * wombat%ligK(i,j,k)
        wombat%feIII(i,j,k) = ( -zval + SQRT( zval*zval + 4.0*wombat%ligK(i,j,k)*fe_sfe ) ) &
                              / ( 2.*wombat%ligK(i,j,k) + epsi )
        wombat%feIII(i,j,k) = max(0.0, min(wombat%feIII(i,j,k), fe_sfe) )
      endif
      wombat%felig(i,j,k) = max(0.0, fe_sfe - wombat%feIII(i,j,k))

      ! Scavenging of Fe` onto biogenic particles
      partic = (sdet_mmolm3*2 + ldet_mmolm3*2 + ldetsi_mmolm3*2 + caco3_mmolm3*8.3) ! total particle concentration [mmol/m3]
      fescaven = wombat%feIII(i,j,k) * (1e-7/86400.0 + wombat%kscav_dfe * partic)
      wombat%fescasafe(i,j,k) = fescaven * (sdet_mmolm3*2 + caco3_mmolm3*8.3) / (partic+epsi)
      wombat%fescalafe(i,j,k) = fescaven * (ldet_mmolm3*2 + ldetsi_mmolm3*2) / (partic+epsi)

      ! Coagulation of colloidal Fe (umol/m3) to form sinking particles (mmol/m3)
      ! Following Tagliabue et al. (2023), make coagulation rate dependent on DOC and Phytoplankton biomass
      biof = (phy_mmolm3 + dia_mmolm3) / (phy_mmolm3 + dia_mmolm3 + 0.03)
      shear = merge(1.0, 0.01, wombat%zw(i,j,k) <= hblt_depth(i,j))
      ! Colloidal shunt associated with small particles and DOC (Tagliabue et al., 2023)
      ! NOTE: we recycle "fesol3" because they are already defined as real variables
      feagg1 = 10.8 * biof ! 12 * 3 * 0.3 (Tagliabue et al., 2023; *3 (DOC effect) *0.3 (phytoplankton effect))
      feagg2 = 9.05
      feagg3 = 2.49
      feagg4 = 115.02 * biof ! 127.8 * 3 * 0.3 (Tagliabue et al., 2023; *3 (DOC effect) *0.3 (phytoplankton effect))
      feagg5 = 725.7
      zval = ( shear*(feagg1*(sdoc_mmolm3 + ldoc_mmolm3) + feagg2*sdet_mmolm3) + feagg3*sdet_mmolm3 &
               + feagg4*(sdoc_mmolm3 + ldoc_mmolm3) + feagg5*sdet_mmolm3 ) * wombat%kcoag_dfe
      wombat%fecoag2safe(i,j,k) = wombat%fecol(i,j,k) * zval
      ! Include an aggregation of colloidal authigenic Fe when concentration of colloidal Fe is high
      wombat%fecoag2safe(i,j,k) = wombat%fecoag2safe(i,j,k) + wombat%kagg_col &
                                 * wombat%fecol(i,j,k)**4 / (wombat%fecol(i,j,k)**4 + wombat%kagg_kcol**4)
      ! Colloidal shunt associated with big particles (Tagliabue et al., 2023)
      feagg1 = 1.37
      feagg2 = 1.94
      zval = (( shear*2.0 + feagg1)*ldet_mmolm3 + feagg2*ldet_mmolm3 ) * wombat%kcoag_dfe
      wombat%fecoag2lafe(i,j,k) = wombat%fecol(i,j,k) * zval

      ! dissolution of Fe from authigenic particles back to dissolved phase
      wombat%safediss(i,j,k) = wombat%ksafe_dfe * safe_p
      wombat%lafediss(i,j,k) = wombat%klafe_dfe * lafe_p

      ! Convert the terms back to mol/kg
      wombat%fescasafe(i,j,k) = wombat%fescasafe(i,j,k) * umol_m3_to_mol_kg
      wombat%fescalafe(i,j,k) = wombat%fescalafe(i,j,k) * umol_m3_to_mol_kg
      wombat%fecoag2safe(i,j,k) = wombat%fecoag2safe(i,j,k) * umol_m3_to_mol_kg
      wombat%fecoag2lafe(i,j,k) = wombat%fecoag2lafe(i,j,k) * umol_m3_to_mol_kg
      wombat%feIII(i,j,k) = wombat%feIII(i,j,k) * umol_m3_to_mol_kg
      wombat%felig(i,j,k) = wombat%felig(i,j,k) * umol_m3_to_mol_kg
      wombat%ligK(i,j,k) = wombat%ligK(i,j,k) * 1e9 ! PJB: convert back to L/mol
      wombat%fecol(i,j,k) = wombat%fecol(i,j,k) * umol_m3_to_mol_kg


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 11] Silicic acid chemistry and biogenic silica dissolution     !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Solve for the equilibrium of silicic acid (Si(OH)4) in seawater
      ! Do so via first-principles thermodynamics, where
      ! K(T,P) = y * m(eq) / ( a(H2O)**2 )
        ! 1. K(T,P0) is the thermodynamic equilibrium constant
        ! 2. gamma0 is the activity coefficient of neutral H4Si(OH)4 in seawater
        ! 3. alphaH2O is the activity of water in seawater and ~ 0.999
        ! 4. deltaV0 is the change in standard partial molar volume of adding
        !    Si(OH)4 to seawater from amorphorus silica
        ! 5. spmvcorrect is the pressure correction. The P-dependence follows the
        !    HKF/SUPCRT approach using standard partial molal volumes
        ! 6. sileqc is the equilibrium molality of Si(OH)4 in seawater (mol/kg)
        !    accounting for temperature, salinity and pressure effects
      zval = max(273.15, Temp(i,j,k) + 273.15)  ! temperature in Kelvin
      K_am_silica = 10**( -8.476 - 485.24/zval - 2.268e-6*zval*zval + 3.068 * log10(zval) ) ! Gunnarsson & Arnorsson (2000)
      gamma0 = 1.0 + 0.0053 * Salt(i,j,k) - 0.000034 * Salt(i,j,k)**2 ! Savenko 2014 Oceanology
      alphaH2O = 0.999  ! activity of water in seawater (from TEOS-10)
      deltaV0 = -9.0 * 1e-6 ! m3/mol, (Willey 1982 Geochim. et Cosmochim. Acta) (also see Loucaides et al., 2012 Mar. Chem.)
      spmvcorrect = exp( -deltaV0/(Rgas * zval) * wombat%zm(i,j,k) * 1.0e4 ) ! 1.0e4 converts dbar to Pa
      wombat%sileqc(i,j,k) = (K_am_silica * spmvcorrect) * alphaH2O**2.0 / gamma0 ! mol/kg

      ! Dissolution of biogenic silica into silicic acid
        ! 1. Temperature effect
        !    - activation energy of 58 kJ/mol (Kamatani 1982 Marine Biology)
        !    - Kamatani (1982) measure dissolution rates (K, hr-1) dependent on temperature (T, ºC)
        !      according to ln(K) = alpha + 0.0833*T, where alpha is ~-8 to -10.
        !    - Greenwood et al. (2005) finds dissolution rates [1/s] at a given temperature
        !      equal to exp(20 - 10050/T), with T in Kelvin (Figure 4). But these T are too high.
        !    We use Kamatani (1982):
        disssi_temp = exp(wombat%bsi_alpha + 0.0833*Temp(i,j,k)) / 3600.0 ! [1/s]
        ! 2. Undersaturation term
        !    - see Eq. 2.13 and fits of this equation to ocean data in Figures 3.20 and 3.21 in
        !      Rickert, D., Dissolution kinetics of biogenic silica in marine environments, Ber. Polarforsch., 351, 2000.
        !    - From Van Cappellen et al., (2002) Global Biogeochemical Cycles:
        !      "Detailed kinetic studies of biogenic silica dissolution conducted in flow-through reactors
        !       demonstrate that at very high degrees of undersaturation the dissolution kinetics switch
        !       from a linear dependence on the degree of undersaturation to an exponential one [Van Cappellen
        !       and Qiu, 1997b; Rickert, 2000]."
        !    - However, we apply an exponent of 1 because on acid-cleaned silica, the dissolution proceeds with
        !      undersaturation linearly (Rickert 2000 - Table 3.4). This is thermodynamically defendable.
      disssi_usat = (1 - min(1.0, sil_p / wombat%sileqc(i,j,k)) )
      ! 3. Bio-interference term?
        !    - "The removal of organic or inorganic coatings enhance the reactivity by at least an order of magnitude."
        !      Ricket et al., 2002 Geochim. et Cosmochim. Acta
        !    - Diatom frustule dissolution increased by order of magnitude with bacteria (Bidle & Azam 1999 Nature)
        !    - During a bloom off Monterey Bay, anti-biotics decreased dissolution by ~50% (Bidle et al., 2003 Limnol. Oceanogr.)
      disssi_bact = 1.0 + wombat%bsi_fbac * (lbac_mmolm3 / ( lbac_mmolm3 + wombat%bsi_kbac ))
        ! 4. Dissolution rate of biogenic silica (/s) composed of the above terms
      wombat%disssi(i,j,k) = disssi_temp * disssi_usat * disssi_bact

      ! Biogenic silica dissolution
      wombat%bsidiss(i,j,k) = wombat%disssi(i,j,k) * ldetsi_p ! [mol/kg/s]


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 12] Mortality terms (will need to make implicit later)         !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Mortality terms
      if (phy_mmolm3>1e-3) then
        wombat%phymorl(i,j,k) = wombat%phylmor * fbc * phy_p ! [molC/kg/s]
        wombat%phymorq(i,j,k) = wombat%phyqmor / mmol_m3_to_mol_kg * phy_p * phy_p ! [molC/kg/s]
      else
        wombat%phymorl(i,j,k) = 0.0
        wombat%phymorq(i,j,k) = 0.0
      endif
      if (dia_mmolm3>1e-3) then
        wombat%diamorl(i,j,k) = wombat%dialmor * fbc * dia_p ! [molC/kg/s]
        wombat%diamorq(i,j,k) = wombat%diaqmor / mmol_m3_to_mol_kg * dia_p * dia_p ! [molC/kg/s]
      else
        wombat%diamorl(i,j,k) = 0.0
        wombat%diamorq(i,j,k) = 0.0
      endif
      if (zoo_mmolm3>1e-3) then
        wombat%zoomorl(i,j,k) = wombat%zoolmor * fbc * zoo_p ! [molC/kg/s]
        wombat%zoomorq(i,j,k) = wombat%zooqmor / mmol_m3_to_mol_kg * zoo_p * zoo_p ! [molC/kg/s]
      else
        wombat%zoomorl(i,j,k) = 0.0
        wombat%zoomorq(i,j,k) = 0.0
      endif
      if (mes_mmolm3>1e-3) then
        wombat%mesmorl(i,j,k) = wombat%meslmor * fbc * mes_p ! [molC/kg/s]
        wombat%mesmorq(i,j,k) = wombat%mesqmor / mmol_m3_to_mol_kg * mes_p * mes_p ! [molC/kg/s]
      else
        wombat%mesmorl(i,j,k) = 0.0
        wombat%mesmorq(i,j,k) = 0.0
      endif
      if (lbac_mmolm3>1e-3) then
        wombat%lbacmorl(i,j,k) = wombat%baclmor * fbc * lbac_p ! [molC/kg/s]
        wombat%lbacmorq(i,j,k) = wombat%bacqmor / mmol_m3_to_mol_kg * lbac_p * lbac_p ! [molC/kg/s]
      else
        wombat%lbacmorl(i,j,k) = 0.0
        wombat%lbacmorq(i,j,k) = 0.0
      endif
      if (sbac_mmolm3>1e-3) then
        wombat%sbacmorl(i,j,k) = wombat%baclmor * fbc * sbac_p ! [molC/kg/s]
        wombat%sbacmorq(i,j,k) = wombat%bacqmor / mmol_m3_to_mol_kg * sbac_p * sbac_p ! [molC/kg/s]
      else
        wombat%sbacmorl(i,j,k) = 0.0
        wombat%sbacmorq(i,j,k) = 0.0
      endif
      if (aoa_mmolm3>1e-3) then
        wombat%aoamorl(i,j,k) = wombat%aoalmor * fbc * aoa_p ! [molC/kg/s]
        wombat%aoamorq(i,j,k) = wombat%aoaqmor / mmol_m3_to_mol_kg * aoa_p * aoa_p ! [molC/kg/s]
      else
        wombat%aoamorl(i,j,k) = 0.0
        wombat%aoamorq(i,j,k) = 0.0
      endif


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 13] Zooplankton grazing                                        !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      !!!~~~ Zooplankton ~~~!!!
      ! Grazing function ! [1/s]
      ! normalize the prey preference kernal to reflect dietary fractions (Gentleman et al., (2003) DSRII)
      I_denom = 1.0 / ( wombat%zpreflbac + wombat%zprefsbac +wombat%zprefaoa &
                      + wombat%zprefphy + wombat%zprefdia + wombat%zprefsdet )
      wombat%zoopreflbac(i,j,k) = wombat%zpreflbac * I_denom
      wombat%zooprefsbac(i,j,k) = wombat%zprefsbac * I_denom
      wombat%zooprefaoa(i,j,k) = wombat%zprefaoa * I_denom
      wombat%zooprefphy(i,j,k) = wombat%zprefphy * I_denom
      wombat%zooprefdia(i,j,k) = wombat%zprefdia * I_denom
      wombat%zooprefsdet(i,j,k) = wombat%zprefsdet * I_denom

      ! Gentleman et al. (2003) DSRII
      !   - add a switching component designed to weight the diet towards abundant prey
      !   - see their Eq. 19
      ! Emulates empirical basis of selective feeding on more abundant prey (Kiorboe et al., 2017; L&O)
      !   ... if denominator is zero, then set all preferences to 1/3 (this is a failsafe, but it should not happen)
      zval = wombat%zoopreflbac(i,j,k) + wombat%zooprefsbac(i,j,k) + wombat%zooprefaoa(i,j,k) &
              + wombat%zooprefphy(i,j,k) + wombat%zooprefdia(i,j,k) + wombat%zooprefsdet(i,j,k)
      if (zval < epsi) then
        wombat%zoopreflbac(i,j,k) = 1.0/6.0; wombat%zooprefsbac(i,j,k) = 1.0/6.0; wombat%zooprefaoa(i,j,k) = 1.0/6.0
        wombat%zooprefphy(i,j,k) = 1.0/6.0; wombat%zooprefdia(i,j,k) = 1.0/6.0; wombat%zooprefsdet(i,j,k) = 1.0/6.0
      else
        wzlbac = (wombat%zoopreflbac(i,j,k) * lbac_mmolm3)**wombat%zoopreyswitch
        wzsbac = (wombat%zooprefsbac(i,j,k) * sbac_mmolm3)**wombat%zoopreyswitch
        wzaoa = (wombat%zooprefaoa(i,j,k) * aoa_mmolm3)**wombat%zoopreyswitch
        wzphy = (wombat%zooprefphy(i,j,k) * phy_mmolm3)**wombat%zoopreyswitch
        wzdia = (wombat%zooprefdia(i,j,k) * dia_mmolm3)**wombat%zoopreyswitch
        wzsdet = (wombat%zooprefsdet(i,j,k) * sdet_mmolm3)**wombat%zoopreyswitch
        I_wzsum = 1.0 / (wzlbac + wzsbac + wzaoa + wzphy + wzdia + wzsdet + epsi)
        wombat%zoopreflbac(i,j,k) = wzlbac * I_wzsum
        wombat%zooprefsbac(i,j,k) = wzsbac * I_wzsum
        wombat%zooprefaoa(i,j,k) = wzaoa * I_wzsum
        wombat%zooprefphy(i,j,k) = wzphy * I_wzsum
        wombat%zooprefdia(i,j,k) = wzdia * I_wzsum
        wombat%zooprefsdet(i,j,k) = wzsdet * I_wzsum
      endif
      ! Compute sum of prey-specific Type-III terms to obtain grazing rate [1/s]
      !  - this avoids "perfect substitution" of prey types and aligns with reccommendations of Gentleman et al. (2003)
      Xzoo = (  wombat%zooepslbac * (wombat%zoopreflbac(i,j,k) * lbac_mmolm3)**2 &
              + wombat%zooepssbac * (wombat%zooprefsbac(i,j,k) * sbac_mmolm3)**2 &
              + wombat%zooepsaoa * (wombat%zooprefaoa(i,j,k) * aoa_mmolm3)**2 &
              + wombat%zooepsphy * (wombat%zooprefphy(i,j,k) * phy_mmolm3)**2 &
              + wombat%zooepsdia * (wombat%zooprefdia(i,j,k) * dia_mmolm3)**2 &
              + wombat%zooepssdet * (wombat%zooprefsdet(i,j,k) * sdet_mmolm3)**2)
      g_zoo = wombat%zoogmax * fbc * Xzoo / (wombat%zoogmax * fbc + Xzoo)

      ! Grazing, egestion, excretion and assimilation
      if (Xzoo>epsi) then
        I_Xzoo = 1.0 / Xzoo
        ! find "apparent" community epsilon (prey capture rate coefficient)
        wombat%zooeps(i,j,k) = Xzoo / ( (wombat%zoopreflbac(i,j,k) * lbac_mmolm3)**2 &
                                      + (wombat%zooprefsbac(i,j,k) * sbac_mmolm3)**2 &
                                      + (wombat%zooprefaoa(i,j,k) * aoa_mmolm3)**2 &
                                      + (wombat%zooprefphy(i,j,k) * phy_mmolm3)**2 &
                                      + (wombat%zooprefdia(i,j,k) * dia_mmolm3)**2 &
                                      + (wombat%zooprefsdet(i,j,k) * sdet_mmolm3)**2 )
        wombat%zoograzlbac(i,j,k) = g_zoo * zoo_p * wombat%zooepslbac*(wombat%zoopreflbac(i,j,k)*lbac_mmolm3)**2 * I_Xzoo ! [molC/kg/s]
        wombat%zoograzsbac(i,j,k) = g_zoo * zoo_p * wombat%zooepssbac*(wombat%zooprefsbac(i,j,k)*sbac_mmolm3)**2 * I_Xzoo ! [molC/kg/s]
        wombat%zoograzaoa(i,j,k) = g_zoo * zoo_p * wombat%zooepsaoa*(wombat%zooprefaoa(i,j,k)*aoa_mmolm3)**2 * I_Xzoo ! [molC/kg/s]
        wombat%zoograzphy(i,j,k) = g_zoo * zoo_p * wombat%zooepsphy*(wombat%zooprefphy(i,j,k)*phy_mmolm3)**2 * I_Xzoo ! [molC/kg/s]
        wombat%zoograzdia(i,j,k) = g_zoo * zoo_p * wombat%zooepsdia*(wombat%zooprefdia(i,j,k)*dia_mmolm3)**2 * I_Xzoo ! [molC/kg/s]
        wombat%zoograzsdet(i,j,k) = g_zoo * zoo_p * wombat%zooepssdet*(wombat%zooprefsdet(i,j,k)*sdet_mmolm3)**2 * I_Xzoo ! [molC/kg/s]
      else
        wombat%zoograzlbac(i,j,k) = 0.0
        wombat%zoograzsbac(i,j,k) = 0.0
        wombat%zoograzaoa(i,j,k) = 0.0
        wombat%zoograzphy(i,j,k) = 0.0
        wombat%zoograzdia(i,j,k) = 0.0
        wombat%zoograzsdet(i,j,k) = 0.0
      endif
      ! We follow Le Mezo & Galbraith (2021) L&O - The fecal iron pump: Global impact of animals on the iron stoichiometry...
      !  - ingestion, assimilation and excretion of carbon and iron by zooplankton are calculated separately
      !  - the idea is to enrich fecal pellets in iron compared to carbon
      wombat%zooexcrlbac(i,j,k) = wombat%zoograzlbac(i,j,k) * wombat%zooCingest*(1.0 - wombat%zooCassim)
      wombat%zooexcrsbac(i,j,k) = wombat%zoograzsbac(i,j,k) * wombat%zooCingest*(1.0 - wombat%zooCassim)
      wombat%zooexcraoa(i,j,k) = wombat%zoograzaoa(i,j,k) * wombat%zooCingest*(1.0 - wombat%zooCassim)
      wombat%zooexcrphy(i,j,k) = wombat%zoograzphy(i,j,k) * wombat%zooCingest*(1.0 - wombat%zooCassim)
      wombat%zooexcrdia(i,j,k) = wombat%zoograzdia(i,j,k) * wombat%zooCingest*(1.0 - wombat%zooCassim)
      wombat%zooexcrsdet(i,j,k) = wombat%zoograzsdet(i,j,k) * wombat%zooCingest*(1.0 - wombat%zooCassim)
      wombat%zooegeslbac(i,j,k) = wombat%zoograzlbac(i,j,k) * (1.0-wombat%zooCingest)
      wombat%zooegessbac(i,j,k) = wombat%zoograzsbac(i,j,k) * (1.0-wombat%zooCingest)
      wombat%zooegesaoa(i,j,k) = wombat%zoograzaoa(i,j,k) * (1.0-wombat%zooCingest)
      wombat%zooegesphy(i,j,k) = wombat%zoograzphy(i,j,k) * (1.0-wombat%zooCingest)
      wombat%zooegesdia(i,j,k) = wombat%zoograzdia(i,j,k) * (1.0-wombat%zooCingest)
      wombat%zooegessdet(i,j,k) = wombat%zoograzsdet(i,j,k) * (1.0-wombat%zooCingest)
      zooegeslbacfe = wombat%zoograzlbac(i,j,k) / wombat%bac_C2Fe * (1.0-wombat%zooFeingest)
      zooegessbacfe = wombat%zoograzsbac(i,j,k) / wombat%bac_C2Fe * (1.0-wombat%zooFeingest)
      zooegesaoafe = wombat%zoograzaoa(i,j,k) / wombat%aoa_C2Fe * (1.0-wombat%zooFeingest)
      zooegesphyfe = wombat%zoograzphy(i,j,k) * phy_Fe2C * (1.0-wombat%zooFeingest)
      zooegesdiafe = wombat%zoograzdia(i,j,k) * dia_Fe2C * (1.0-wombat%zooFeingest)
      zooegessdetfe = wombat%zoograzsdet(i,j,k) * sdet_Fe2C * (1.0-wombat%zooFeingest)
      zooassilbacfe = wombat%zoograzlbac(i,j,k) / wombat%bac_C2Fe * wombat%zooFeingest*wombat%zooFeassim
      zooassisbacfe = wombat%zoograzsbac(i,j,k) / wombat%bac_C2Fe * wombat%zooFeingest*wombat%zooFeassim
      zooassiaoafe = wombat%zoograzaoa(i,j,k) / wombat%aoa_C2Fe * wombat%zooFeingest*wombat%zooFeassim
      zooassiphyfe = wombat%zoograzphy(i,j,k) * phy_Fe2C * wombat%zooFeingest*wombat%zooFeassim
      zooassidiafe = wombat%zoograzdia(i,j,k) * dia_Fe2C * wombat%zooFeingest*wombat%zooFeassim
      zooassisdetfe = wombat%zoograzsdet(i,j,k) * sdet_Fe2C * wombat%zooFeingest*wombat%zooFeassim
      zooexcrlbacfe = wombat%zoograzlbac(i,j,k) / wombat%bac_C2Fe * wombat%zooFeingest*(1.0 - wombat%zooFeassim)
      zooexcrsbacfe = wombat%zoograzsbac(i,j,k) / wombat%bac_C2Fe * wombat%zooFeingest*(1.0 - wombat%zooFeassim)
      zooexcraoafe = wombat%zoograzaoa(i,j,k) / wombat%aoa_C2Fe * wombat%zooFeingest*(1.0 - wombat%zooFeassim)
      zooexcrphyfe = wombat%zoograzphy(i,j,k) * phy_Fe2C * wombat%zooFeingest*(1.0 - wombat%zooFeassim)
      zooexcrdiafe = wombat%zoograzdia(i,j,k) * dia_Fe2C * wombat%zooFeingest*(1.0 - wombat%zooFeassim)
      zooexcrsdetfe = wombat%zoograzsdet(i,j,k) * sdet_Fe2C * wombat%zooFeingest*(1.0 - wombat%zooFeassim)
      zooexcrlbacn  = wombat%zoograzlbac(i,j,k) / wombat%bac_C2N &
                     - (wombat%zoograzlbac(i,j,k) * wombat%zooCingest * wombat%zooCassim / (122.0/16.0)) &
                     - (wombat%zooegeslbac(i,j,k) / (122.0/16.0)) ! [molN/kg/s]
      zooexcrsbacn  = wombat%zoograzsbac(i,j,k) / wombat%bac_C2N &
                     - (wombat%zoograzsbac(i,j,k) * wombat%zooCingest * wombat%zooCassim / (122.0/16.0)) &
                     - (wombat%zooegessbac(i,j,k) / (122.0/16.0)) ! [molN/kg/s]
      zooexcraoan  = wombat%zoograzaoa(i,j,k) / wombat%aoa_C2N &
                     - (wombat%zoograzaoa(i,j,k) * wombat%zooCingest * wombat%zooCassim / (122.0/16.0)) &
                     - (wombat%zooegesaoa(i,j,k) / (122.0/16.0)) ! [molN/kg/s]


      !!!~~~ Mesozooplankton ~~~!!!
      ! Grazing function ! [1/s]
      ! normalize the prey preference kernal to reflect dietary fractions (Gentleman et al., (2003) DSRII)
      I_denom = 1.0 / ( wombat%mpreflbac + wombat%mprefsbac + wombat%mprefaoa &
                      + wombat%mprefphy + wombat%mprefdia &
                      + wombat%mprefsdet + wombat%mprefzoo )
      wombat%mespreflbac(i,j,k) = wombat%mpreflbac * I_denom
      wombat%mesprefsbac(i,j,k) = wombat%mprefsbac * I_denom
      wombat%mesprefaoa(i,j,k) = wombat%mprefaoa * I_denom
      wombat%mesprefphy(i,j,k) = wombat%mprefphy * I_denom
      wombat%mesprefdia(i,j,k) = wombat%mprefdia * I_denom
      wombat%mesprefsdet(i,j,k) = wombat%mprefsdet * I_denom
      wombat%mesprefldet(i,j,k) = wombat%mprefldet * I_denom
      wombat%mesprefzoo(i,j,k) = wombat%mprefzoo * I_denom
      zval = wombat%mespreflbac(i,j,k) + wombat%mesprefsbac(i,j,k) + wombat%mesprefaoa(i,j,k) + wombat%mesprefphy(i,j,k) &
           + wombat%mesprefdia(i,j,k) + wombat%mesprefsdet(i,j,k) + wombat%mesprefldet(i,j,k) + wombat%mesprefzoo(i,j,k)
      ! Gentleman et al. (2003) DSRII
      !   - add a switching component designed to weight the diet towards abundant prey
      !   - see their Eq. 19
      ! Emulates empirical basis of selective feeding on more abundant prey (Kiorboe et al., 2017; L&O)
      !   ... if denominator is zero, then set all preferences to 1/3 (this is a failsafe, but it should not happen)
      if (zval < 1e-20) then
        wombat%mespreflbac(i,j,k) = 1.0/8.0; wombat%mesprefsbac(i,j,k) = 1.0/8.0; wombat%mesprefaoa(i,j,k) = 1.0/8.0
        wombat%mesprefphy(i,j,k) = 1.0/8.0; wombat%mesprefdia(i,j,k) = 1.0/8.0; wombat%mesprefsdet(i,j,k) = 1.0/8.0
        wombat%mesprefldet(i,j,k) = 1.0/8.0; wombat%mesprefzoo(i,j,k) = 1.0/8.0
      else
        wzlbac = (wombat%mespreflbac(i,j,k) * lbac_mmolm3)**wombat%mespreyswitch
        wzsbac = (wombat%mesprefsbac(i,j,k) * sbac_mmolm3)**wombat%mespreyswitch
        wzaoa = (wombat%mesprefaoa(i,j,k) * aoa_mmolm3)**wombat%mespreyswitch
        wzphy = (wombat%mesprefphy(i,j,k) * phy_mmolm3)**wombat%mespreyswitch
        wzdia = (wombat%mesprefdia(i,j,k) * dia_mmolm3)**wombat%mespreyswitch
        wzsdet= (wombat%mesprefsdet(i,j,k) * sdet_mmolm3)**wombat%mespreyswitch
        wzldet= (wombat%mesprefldet(i,j,k) * ldet_mmolm3)**wombat%mespreyswitch
        wzzoo = (wombat%mesprefzoo(i,j,k) * zoo_mmolm3)**wombat%mespreyswitch
        I_wzsum = 1.0 / (wzlbac + wzsbac + wzaoa + wzphy + wzdia + wzsdet + wzldet + wzzoo + epsi)
        wombat%mespreflbac(i,j,k) = wzlbac * I_wzsum
        wombat%mesprefsbac(i,j,k) = wzsbac * I_wzsum
        wombat%mesprefaoa(i,j,k) = wzaoa * I_wzsum
        wombat%mesprefphy(i,j,k) = wzphy * I_wzsum
        wombat%mesprefdia(i,j,k) = wzdia * I_wzsum
        wombat%mesprefsdet(i,j,k) = wzsdet * I_wzsum
        wombat%mesprefldet(i,j,k) = wzldet * I_wzsum
        wombat%mesprefzoo(i,j,k) = wzzoo * I_wzsum
      endif
      ! Compute sum of prey-specific Type-III terms to obtain grazing rate [1/s]
      !  - this avoids "perfect substitution" of prey types and aligns with reccommendations of Gentleman et al. (2003)
      Xmes = (  wombat%mesepslbac * (wombat%mespreflbac(i,j,k) * lbac_mmolm3)**2 &
              + wombat%mesepssbac * (wombat%mesprefsbac(i,j,k) * sbac_mmolm3)**2 &
              + wombat%mesepsaoa * (wombat%mesprefaoa(i,j,k) * aoa_mmolm3)**2 &
              + wombat%mesepsphy * (wombat%mesprefphy(i,j,k) * phy_mmolm3)**2 &
              + wombat%mesepsdia * (wombat%mesprefdia(i,j,k) * dia_mmolm3)**2 &
              + wombat%mesepssdet * (wombat%mesprefsdet(i,j,k) * sdet_mmolm3)**2 &
              + wombat%mesepsldet * (wombat%mesprefldet(i,j,k) * ldet_mmolm3)**2 &
              + wombat%mesepszoo * (wombat%mesprefzoo(i,j,k) * zoo_mmolm3)**2 )
      g_mes = wombat%mesgmax * fbc * Xmes / (wombat%mesgmax * fbc + Xmes)

      ! Grazing, egestion, excretion and assimilation
      if (Xmes>epsi) then
        I_Xmes = 1.0 / Xmes
        ! find "apparent" community epsilon (prey capture rate coefficient)
        wombat%meseps(i,j,k) = Xmes / ( (wombat%mespreflbac(i,j,k) * lbac_mmolm3)**2 &
                                      + (wombat%mesprefsbac(i,j,k) * sbac_mmolm3)**2 &
                                      + (wombat%mesprefaoa(i,j,k) * aoa_mmolm3)**2 &
                                      + (wombat%mesprefphy(i,j,k) * phy_mmolm3)**2 &
                                      + (wombat%mesprefdia(i,j,k) * dia_mmolm3)**2 &
                                      + (wombat%mesprefsdet(i,j,k) * sdet_mmolm3)**2 &
                                      + (wombat%mesprefldet(i,j,k) * ldet_mmolm3)**2 &
                                      + (wombat%mesprefzoo(i,j,k) * zoo_mmolm3)**2 )
        wombat%mesgrazlbac(i,j,k) = g_mes * mes_p * wombat%mesepslbac*(wombat%mespreflbac(i,j,k)*lbac_mmolm3)**2 * I_Xmes ! [molC/kg/s]
        wombat%mesgrazsbac(i,j,k) = g_mes * mes_p * wombat%mesepssbac*(wombat%mesprefsbac(i,j,k)*sbac_mmolm3)**2 * I_Xmes ! [molC/kg/s]
        wombat%mesgrazaoa(i,j,k) = g_mes * mes_p * wombat%mesepsaoa*(wombat%mesprefaoa(i,j,k)*aoa_mmolm3)**2 * I_Xmes ! [molC/kg/s]
        wombat%mesgrazphy(i,j,k) = g_mes * mes_p * wombat%mesepsphy*(wombat%mesprefphy(i,j,k)*phy_mmolm3)**2 * I_Xmes ! [molC/kg/s]
        wombat%mesgrazdia(i,j,k) = g_mes * mes_p * wombat%mesepsdia*(wombat%mesprefdia(i,j,k)*dia_mmolm3)**2 * I_Xmes ! [molC/kg/s]
        wombat%mesgrazsdet(i,j,k) = g_mes * mes_p * wombat%mesepssdet*(wombat%mesprefsdet(i,j,k)*sdet_mmolm3)**2 * I_Xmes ! [molC/kg/s]
        wombat%mesgrazldet(i,j,k) = g_mes * mes_p * wombat%mesepsldet*(wombat%mesprefldet(i,j,k)*ldet_mmolm3)**2 * I_Xmes ! [molC/kg/s]
        wombat%mesgrazzoo(i,j,k) = g_mes * mes_p * wombat%mesepszoo*(wombat%mesprefzoo(i,j,k)*zoo_mmolm3)**2 * I_Xmes ! [molC/kg/s]
      else
        wombat%mesgrazlbac(i,j,k) = 0.0
        wombat%mesgrazsbac(i,j,k) = 0.0
        wombat%mesgrazaoa(i,j,k) = 0.0
        wombat%mesgrazphy(i,j,k) = 0.0
        wombat%mesgrazdia(i,j,k) = 0.0
        wombat%mesgrazsdet(i,j,k) = 0.0
        wombat%mesgrazldet(i,j,k) = 0.0
        wombat%mesgrazzoo(i,j,k) = 0.0
      endif
      ! We follow Le Mezo & Galbraith (2021) L&O - The fecal iron pump: Global impact of animals on the iron stoichiometry...
      !  - ingestion, assimilation and excretion of carbon and iron by zooplankton are calculated separately
      !  - the idea is to enrich fecal pellets in iron compared to carbon
      wombat%mesexcrlbac(i,j,k) = wombat%mesgrazlbac(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcrsbac(i,j,k) = wombat%mesgrazsbac(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcraoa(i,j,k) = wombat%mesgrazaoa(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcrphy(i,j,k) = wombat%mesgrazphy(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcrdia(i,j,k) = wombat%mesgrazdia(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcrsdet(i,j,k) = wombat%mesgrazsdet(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcrldet(i,j,k) = wombat%mesgrazldet(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcrzoo(i,j,k) = wombat%mesgrazzoo(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesegeslbac(i,j,k) = wombat%mesgrazlbac(i,j,k) * (1.0 - wombat%mesCingest)
      wombat%mesegessbac(i,j,k) = wombat%mesgrazsbac(i,j,k) * (1.0 - wombat%mesCingest)
      wombat%mesegesaoa(i,j,k) = wombat%mesgrazaoa(i,j,k) * (1.0 - wombat%mesCingest)
      wombat%mesegesphy(i,j,k) = wombat%mesgrazphy(i,j,k) * (1.0 - wombat%mesCingest)
      wombat%mesegesdia(i,j,k) = wombat%mesgrazdia(i,j,k) * (1.0 - wombat%mesCingest)
      wombat%mesegessdet(i,j,k) = wombat%mesgrazsdet(i,j,k) * (1.0 - wombat%mesCingest)
      wombat%mesegesldet(i,j,k) = wombat%mesgrazldet(i,j,k) * (1.0 - wombat%mesCingest)
      wombat%mesegeszoo(i,j,k) = wombat%mesgrazzoo(i,j,k) * (1.0 - wombat%mesCingest)
      mesegeslbacfe = wombat%mesgrazlbac(i,j,k) / wombat%bac_C2Fe * (1.0-wombat%mesFeingest)
      mesegessbacfe = wombat%mesgrazsbac(i,j,k) / wombat%bac_C2Fe * (1.0-wombat%mesFeingest)
      mesegesaoafe = wombat%mesgrazaoa(i,j,k) / wombat%aoa_C2Fe * (1.0-wombat%mesFeingest)
      mesegesphyfe = wombat%mesgrazphy(i,j,k) * phy_Fe2C * (1.0-wombat%mesFeingest)
      mesegesdiafe = wombat%mesgrazdia(i,j,k) * dia_Fe2C * (1.0-wombat%mesFeingest)
      mesegessdetfe = wombat%mesgrazsdet(i,j,k) * sdet_Fe2C * (1.0-wombat%mesFeingest)
      mesegesldetfe = wombat%mesgrazldet(i,j,k) * ldet_Fe2C * (1.0-wombat%mesFeingest)
      mesegeszoofe = wombat%mesgrazzoo(i,j,k) * zoo_Fe2C * (1.0-wombat%mesFeingest)
      mesassilbacfe = wombat%mesgrazlbac(i,j,k) / wombat%bac_C2Fe * wombat%mesFeingest*wombat%mesFeassim
      mesassisbacfe = wombat%mesgrazsbac(i,j,k) / wombat%bac_C2Fe * wombat%mesFeingest*wombat%mesFeassim
      mesassiaoafe = wombat%mesgrazaoa(i,j,k) / wombat%aoa_C2Fe * wombat%mesFeingest*wombat%mesFeassim
      mesassiphyfe = wombat%mesgrazphy(i,j,k) * phy_Fe2C * wombat%mesFeingest*wombat%mesFeassim
      mesassidiafe = wombat%mesgrazdia(i,j,k) * dia_Fe2C * wombat%mesFeingest*wombat%mesFeassim
      mesassisdetfe = wombat%mesgrazsdet(i,j,k) * sdet_Fe2C * wombat%mesFeingest*wombat%mesFeassim
      mesassildetfe = wombat%mesgrazldet(i,j,k) * ldet_Fe2C * wombat%mesFeingest*wombat%mesFeassim
      mesassizoofe = wombat%mesgrazzoo(i,j,k) * zoo_Fe2C * wombat%mesFeingest*wombat%mesFeassim
      mesexcrlbacfe = wombat%mesgrazlbac(i,j,k) / wombat%bac_C2Fe * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrsbacfe = wombat%mesgrazsbac(i,j,k) / wombat%bac_C2Fe * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcraoafe = wombat%mesgrazaoa(i,j,k) / wombat%aoa_C2Fe * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrphyfe = wombat%mesgrazphy(i,j,k) * phy_Fe2C * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrdiafe = wombat%mesgrazdia(i,j,k) * dia_Fe2C * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrsdetfe = wombat%mesgrazsdet(i,j,k) * sdet_Fe2C * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrldetfe = wombat%mesgrazldet(i,j,k) * ldet_Fe2C * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrzoofe = wombat%mesgrazzoo(i,j,k) * zoo_Fe2C * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrlbacn  = wombat%mesgrazlbac(i,j,k) / wombat%bac_C2N &
                     - (wombat%mesgrazlbac(i,j,k) * wombat%mesCingest * wombat%mesCassim / (122.0/16.0)) &
                     - (wombat%mesegeslbac(i,j,k) / (122.0/16.0)) ! [molN/kg/s]
      mesexcrsbacn  = wombat%mesgrazsbac(i,j,k) / wombat%bac_C2N &
                     - (wombat%mesgrazsbac(i,j,k) * wombat%mesCingest * wombat%mesCassim / (122.0/16.0)) &
                     - (wombat%mesegessbac(i,j,k) / (122.0/16.0)) ! [molN/kg/s]
      mesexcraoan  = wombat%mesgrazaoa(i,j,k) / wombat%aoa_C2N &
                     - (wombat%mesgrazaoa(i,j,k) * wombat%mesCingest * wombat%mesCassim / (122.0/16.0)) &
                     - (wombat%mesegesaoa(i,j,k) / (122.0/16.0)) ! [molN/kg/s]


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 14] Implicit nitrogen fixation                                 !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      if (do_nitrogen_fixation) then
        ! Temperature dependent maximum growth rate of Trichodesmium (Jiang et al., 2018)
        if (Temp(i,j,k)>15.8) then
          wombat%trimumax(i,j,k) = ( ( -3.99e-4 * Temp(i,j,k)**3.0 ) + &
                                     (  0.02685 * Temp(i,j,k)**2.0 ) + &
                                     ( -0.555 * Temp(i,j,k) ) + 3.633 ) / 86400.0
        endif
        ! Nutrient and light limitation terms
        wombat%tri_lfer(i,j,k) = fe_umolm3 / (fe_umolm3 + wombat%trikf)
        wombat%tri_lpar(i,j,k) = (1. - exp(-wombat%alphabio_tri * wombat%trichlc * wombat%radbio(i,j,k)))
        ! Nitrogen fixation rate of Trichodesmium
        wombat%nitrfix(i,j,k) = wombat%trimumax(i,j,k) * (1.0 - wombat%phy_lnit(i,j,k)) &
                                * min(wombat%tri_lfer(i,j,k), wombat%tri_lpar(i,j,k)) &
                                * wombat%trin2c * 1e-6  ! 1e-6 scaler to account for biomass in mol/kg
      endif


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 15] Facultative bacterial heterotrophy                         !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! From the base biomass yield on N, compute yields for O2 and anaerobic growth on alternative electron acceptors and DOC
      !  [ Zakem et al., 2020 ISME; Buchanan et al., 2025 Science]
      !  Find electron potential of the bacterial biomass and DOM
      e_pom = max(epsi, 4.0 + 1.65 - 2.0*0.4 - 3.0*16./122.)  ! [Anderson et al., 1995]
      e_ldom = max(epsi, 4.0 + 1.65 - 2.0*0.4 - 3.0*ldom_N2C)  ! [Anderson et al., 1995]
      e_sdom = max(epsi, 4.0 + 1.0 - 2.0*0.6 - 3.0*ldom_N2C*0.0)  ! [Lechtenfeld et al., 2014] - CRAM has H:C of ~1.0, O:C of ~0.6
      e_bac = max(epsi, 4.0 + 1.4 - 2.0*0.4 - 3.0/wombat%bac_C2N) ! [Zimmerman et al., 2014]
      e_lres = max(epsi, e_ldom - wombat%lbac_alpha * e_sdom)
      sdom_N2C = 0.0

      ! Determine the biomass yield in terms of carbon (mol C-biomass per mol DOC consumed)
      wombat%lbacydoc(i,j,k) = min(1.0 - wombat%lbac_alpha, wombat%lbac_fele * e_lres/e_bac)
      wombat%sbacydoc(i,j,k) = wombat%sbac_fele * e_sdom/e_bac

      lbac_cdoc = 0.0; !sbac_cdoc = 0.0 ! reinitialise for safety
      lbac_cdoc_ana = 0.0; !sbac_cdoc_ana = 0.0 ! reinitialise for safety

      ! Determine other resource yields during bacterial growth
      !!!~~~ Large, sharing bacteria ~~~!!!
      ! aerobic conditions
      if (wombat%lbacydoc(i,j,k) > 0.0) lbac_cdoc = 1.0 / wombat%lbacydoc(i,j,k) ! The amount of DOC consumed per mol of large bacterial C-biomass produced
      lbac_coxy = max(0.0, e_lres - wombat%lbacydoc(i,j,k)*e_bac)/4.0 * lbac_cdoc ! The amount of oxygen consumed per mol of C-biomass produced
      lbac_pdoc = wombat%lbac_alpha * lbac_cdoc ! The amount of partially oxidized product produced per mol of C-biomass produced
      lbac_pco2 = (1.0 - wombat%lbac_alpha - wombat%lbacydoc(i,j,k)) * lbac_cdoc ! The amount of CO2 produced per mol of C-biomass produced
      lbac_pnh4 = (ldom_N2C - wombat%lbac_alpha * ldom_N2C*0.0 &
                 - 1.0/wombat%bac_C2N * wombat%lbacydoc(i,j,k)) * lbac_cdoc ! The amount of NH4 produced per mol of C-biomass produced
      ! anaerobic conditions
      lbacydoc_ana = min(1.0 - wombat%lbac_alpha, wombat%lbac_fele * 0.9 * e_lres/e_bac)
      if (lbacydoc_ana > 0.0) lbac_cdoc_ana = 1.0 / lbacydoc_ana ! The amount of DOC consumed per mol of lbacterial C-biomass produced
      lbac_cno3_ana = max(0.0, e_lres * lbacydoc_ana*e_bac)/4.0 * lbac_cdoc_ana ! The amount of N (NO3 --> N2O) molecules consumed per mol of C-biomass produced
      lbac_pdoc_ana = wombat%lbac_alpha * lbac_cdoc_ana ! The amount of partially oxidized product produced per mol of C-biomass produced
      lbac_pco2_ana = (1.0 - wombat%lbac_alpha - lbacydoc_ana) * lbac_cdoc_ana ! The amount of CO2 produced per mol of C-biomass produced
      lbac_pnh4_ana = (ldom_N2C - wombat%lbac_alpha * ldom_N2C*0.0 &
                     - 1.0/wombat%bac_C2N * lbacydoc_ana) * lbac_cdoc_ana ! The amount of NH4 produced per mol of C-biomass produced
      !!!~~~ Small, selfish bacteria ~~~!!!
      ! aerobic conditions
      if (wombat%sbacydoc(i,j,k) > 0.0) sbac_cdoc = 1.0 / wombat%sbacydoc(i,j,k) ! The amount of DOC consumed per mol of small bacterial C-biomass produced
      sbac_coxy = max(0.0, e_sdom - wombat%sbacydoc(i,j,k)*e_bac)/4.0 * sbac_cdoc ! The amount of oxygen consumed per mol of C-biomass produced
      sbac_pco2 = (1.0 - wombat%sbacydoc(i,j,k)) * sbac_cdoc ! The amount of CO2 produced per mol of C-biomass produced
      sbac_pnh4 = (sdom_N2C - 1.0/wombat%bac_C2N * wombat%sbacydoc(i,j,k)) * sbac_cdoc ! The amount of NH4 produced per mol of C-biomass produced
      ! anaerobic conditions
      sbacydoc_ana = wombat%sbac_fele * 0.9 * e_sdom/e_bac
      if (sbacydoc_ana > 0.0) sbac_cdoc_ana = 1.0 / sbacydoc_ana ! The amount of DOC consumed per mol of lbacterial C-biomass produced
      sbac_cno3_ana = max(0.0, e_sdom * sbacydoc_ana*e_bac)/4.0 * sbac_cdoc_ana ! The amount of N (NO3 --> N2O) molecules consumed per mol of C-biomass produced
      sbac_pco2_ana = (1.0 - sbacydoc_ana) * sbac_cdoc_ana ! The amount of CO2 produced per mol of C-biomass produced
      sbac_pnh4_ana = (sdom_N2C - 1.0/wombat%bac_C2N * sbacydoc_ana) * sbac_cdoc_ana ! The amount of NH4 produced per mol of C-biomass produced

      ! Determine growth rates
      !!!~~~ Large, sharing bacteria ~~~!!!
      ! Aerobic
      bac_Voxy = oxy_mmolm3 * wombat%lbac_poxy ! Uptake of O2 (i.e., O2-limited growth)
      bac_VdFe = wombat%lbac_Vmax_dfe * fe_umolm3 / (fe_umolm3 + wombat%lbac_kfer + epsi) ! Uptake of dFe (i.e., Fe-limited growth)
      bac_Voc = wombat%lbac_Vmax_doc * ldoc_mmolm3 / (ldoc_mmolm3 + wombat%lbac_kdoc + epsi) ! Uptake of DOC (i.e., DOC-limited growth)
      bac_gEA = bac_Voxy / (lbac_coxy + epsi) ! Growth of C biomass due to electron acceptor (O2) uptake
      bac_gFe = bac_VdFe * wombat%bac_C2Fe ! Growth of C biomass due to Fe uptake
      bac_gC = bac_Voc * wombat%lbacydoc(i,j,k) ! Growth of C biomass due to DOC uptake
      bac_muaer = max(0.0, min( bac_gC, bac_gFe, bac_gEA ) ) * fbc
      if (bac_gFe<min(bac_gC,bac_gEA)) wombat%lbac_ffelim(i,j,k) = 1.0
      ! Anaerobic
      bac_Vno3 = wombat%lbac_Vmax_no3 * no3_mmolm3 / (no3_mmolm3 + wombat%lbac_kno3 + epsi) ! Uptake of NO3 (i.e., NO3-limited growth)
      bac_gEA = bac_Vno3 / (lbac_cno3_ana + epsi) ! Growth of C biomass due to electron acceptor (NO3) uptake
      bac_gFe = bac_VdFe * wombat%bac_C2Fe ! Growth of C biomass due to Fe uptake
      bac_gC = bac_Voc * lbacydoc_ana ! Growth of C biomass due to DOC uptake
      bac_muana = max(0.0, min( bac_gC, bac_gFe, bac_gEA ) ) * fbc
      if (.not.do_wc_denitrification) bac_muana = 0.0 ! If no denitrification, anaerobic growth is zero
      ! Save occurance of anaerobic growth to array
      if (bac_muana>bac_muaer) wombat%lbac_fanaer(i,j,k) = 1.0
      ! Take the maximum growth rate as the realised growth rate
      wombat%lbac_mu(i,j,k) = max(bac_muaer, bac_muana)

      !!!~~~ Small, selfish bacteria ~~~!!!
      ! Aerobic
      bac_Voxy = oxy_mmolm3 * wombat%sbac_poxy ! Uptake of O2 (i.e., O2-limited growth)
      bac_VdFe = wombat%sbac_Vmax_dfe * fe_umolm3 / (fe_umolm3 + wombat%sbac_kfer + epsi) ! Uptake of dFe (i.e., Fe-limited growth)
      bac_Voc = wombat%sbac_Vmax_doc * sdoc_mmolm3 / (sdoc_mmolm3 + wombat%sbac_kdoc + epsi) ! Uptake of DOC (i.e., DOC-limited growth)
      bac_gEA = bac_Voxy / (sbac_coxy + epsi) ! Growth of C biomass due to electron acceptor (O2) uptake
      bac_gFe = bac_VdFe * wombat%bac_C2Fe ! Growth of C biomass due to Fe uptake
      bac_gC = bac_Voc * wombat%sbacydoc(i,j,k) ! Growth of C biomass due to DOC uptake
      bac_muaer = max(0.0, min( bac_gC, bac_gFe, bac_gEA ) ) * fbc
      if (bac_gFe<min(bac_gC,bac_gEA)) wombat%sbac_ffelim(i,j,k) = 1.0
      ! Anaerobic
      bac_Vno3 = wombat%sbac_Vmax_no3 * no3_mmolm3 / (no3_mmolm3 + wombat%sbac_kno3 + epsi) ! Uptake of NO3 (i.e., NO3-limited growth)
      bac_gEA = bac_Vno3 / (sbac_cno3_ana + epsi) ! Growth of C biomass due to electron acceptor (NO3) uptake
      bac_gFe = bac_VdFe * wombat%bac_C2Fe ! Growth of C biomass due to Fe uptake
      bac_gC = bac_Voc * sbacydoc_ana ! Growth of C biomass due to DOC uptake
      bac_muana = max(0.0, min( bac_gC, bac_gFe, bac_gEA ) ) * fbc
      if (.not.do_wc_denitrification) bac_muana = 0.0 ! If no denitrification, anaerobic growth is zero
      ! Save occurance of anaerobic growth to array
      if (bac_muana>bac_muaer) wombat%sbac_fanaer(i,j,k) = 1.0
      ! Take the maximum growth rate as the realised growth rate
      wombat%sbac_mu(i,j,k) = max(bac_muaer, bac_muana)

      ! Sources and sinks due to heterotrophic bacterial activity
      wombat%lbacgrow(i,j,k) = wombat%lbac_mu(i,j,k) * lbac_p ! [molC/kg/s]
      wombat%ldocremi(i,j,k) = wombat%lbacgrow(i,j,k) * lbac_cdoc * (1. - wombat%lbac_fanaer(i,j,k)) &
                             + wombat%lbacgrow(i,j,k) * lbac_cdoc_ana * wombat%lbac_fanaer(i,j,k) ! [molC/kg/s]
      wombat%sdocprod(i,j,k) = wombat%lbacgrow(i,j,k) * lbac_pdoc * (1. - wombat%lbac_fanaer(i,j,k)) &
                             + wombat%lbacgrow(i,j,k) * lbac_pdoc_ana * wombat%lbac_fanaer(i,j,k) ! [molC/kg/s]
      wombat%lbacresp(i,j,k) = wombat%lbacgrow(i,j,k) * lbac_coxy * (1. - wombat%lbac_fanaer(i,j,k)) ! [molO2/kg/s]
      wombat%lbacpco2(i,j,k) = wombat%lbacgrow(i,j,k) * lbac_pco2 * (1. - wombat%lbac_fanaer(i,j,k)) &
                            + wombat%lbacgrow(i,j,k) * lbac_pco2_ana * wombat%lbac_fanaer(i,j,k) ! [molCO2/kg/s]
      wombat%lbacdeni(i,j,k) = wombat%lbacgrow(i,j,k) * lbac_cno3_ana * wombat%lbac_fanaer(i,j,k) ! [molNO3/kg/s]
      wombat%lbacufer(i,j,k) = wombat%lbacgrow(i,j,k) / wombat%bac_C2Fe ! [molFe/kg/s]
      wombat%lbacpnh4(i,j,k) = wombat%lbacgrow(i,j,k) * lbac_pnh4 * (1. - wombat%lbac_fanaer(i,j,k)) &
                            + wombat%lbacgrow(i,j,k) * lbac_pnh4_ana * wombat%lbac_fanaer(i,j,k) ! [molN/kg/s]

      wombat%sbacgrow(i,j,k) = wombat%sbac_mu(i,j,k) * sbac_p ! [molC/kg/s]
      wombat%sdocremi(i,j,k) = wombat%sbacgrow(i,j,k) * sbac_cdoc * (1. - wombat%sbac_fanaer(i,j,k)) &
                             + wombat%sbacgrow(i,j,k) * sbac_cdoc_ana * wombat%sbac_fanaer(i,j,k) ! [molC/kg/s]
      wombat%sbacresp(i,j,k) = wombat%sbacgrow(i,j,k) * sbac_coxy * (1. - wombat%sbac_fanaer(i,j,k)) ! [molO2/kg/s]
      wombat%sbacpco2(i,j,k) = wombat%sbacgrow(i,j,k) * sbac_pco2 * (1. - wombat%sbac_fanaer(i,j,k)) &
                            + wombat%sbacgrow(i,j,k) * sbac_pco2_ana * wombat%sbac_fanaer(i,j,k) ! [molCO2/kg/s]
      wombat%sbacdeni(i,j,k) = wombat%sbacgrow(i,j,k) * sbac_cno3_ana * wombat%sbac_fanaer(i,j,k) ! [molNO3/kg/s] !pjb - change this to N2O
      wombat%sbacufer(i,j,k) = wombat%sbacgrow(i,j,k) / wombat%bac_C2Fe ! [molFe/kg/s]
      wombat%sbacpnh4(i,j,k) = wombat%sbacgrow(i,j,k) * sbac_pnh4 * (1. - wombat%sbac_fanaer(i,j,k)) &
                            + wombat%sbacgrow(i,j,k) * sbac_pnh4_ana * wombat%sbac_fanaer(i,j,k) ! [molN/kg/s]



      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 16] Calcium carbonate production and dissolution               !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      if (do_caco3_dynamics) then
        ! PIC:POC ratio is a function of the substrate:inhibitor ratio, which is the
        !  HCO3- to free H+ ions ratio (mol/umol), following Lehmann & Bach (2024).
        !  We also add a T-dependent function to scale down CaCO3 production in waters colder
        !  than 3 degrees C based off the observation of no E hux growth beneath this (Fielding 2013; L&O)
        hco3 = wombat%f_dic(i,j,k) - wombat%co3(i,j,k) - wombat%co2_star(i,j,k)
        wombat%pic2poc(i,j,k) = min(0.3, (wombat%f_inorg + 10.0**(min(2.0, -3.0 + 4.31e-6 * &
                                          hco3 / wombat%htotal(i,j,k)))) * &
                                         (0.55 + 0.45 * tanh(Temp(i,j,k) - 4.0)) )

        ! The dissolution rate is a function of omegas for calcite and aragonite, as well the
        !  concentration of POC, following Kwon et al., 2024, Science Advances; Table S1, and
        !  we account for the dissolution due to zooplankton grazing on particulates
        wombat%dissratcal(i,j,k) = (wombat%disscal * max(0.0, 1.0 - wombat%omega_cal(i,j,k))**2.2)
        wombat%dissratara(i,j,k) = (wombat%dissara * max(0.0, 1.0 - wombat%omega_ara(i,j,k))**1.5)
        wombat%dissratpoc(i,j,k) = (wombat%dissdet * wombat%reminr(i,j,k) * sdet_mmolm3**2.0)
      else
        wombat%pic2poc(i,j,k) = wombat%f_inorg + 0.025
        wombat%dissratcal(i,j,k) = wombat%caco3lrem
        wombat%dissratara(i,j,k) = 0.0
        wombat%dissratpoc(i,j,k) = 0.0
      endif
      if (sdet_p > epsi) then
        wombat%zoodiss(i,j,k) = wombat%zoograzsdet(i,j,k) * wombat%fgutdiss * caco3_mmolm3/sdet_mmolm3
        wombat%mesdiss(i,j,k) = wombat%mesgrazsdet(i,j,k) * wombat%fgutdiss * caco3_mmolm3/sdet_mmolm3
      endif
      wombat%caldiss(i,j,k) = wombat%dissratcal(i,j,k) * caco3_p ! [mol/kg/s]
      wombat%aradiss(i,j,k) = wombat%dissratara(i,j,k) * caco3_p ! [mol/kg/s]
      wombat%pocdiss(i,j,k) = wombat%dissratpoc(i,j,k) * caco3_p ! [mol/kg/s]


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 17] Chemoautotrophy                                            !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! 1. find max growth rate of AOA dependent on temperature (Qin et al., 2015 PNAS)
      wombat%aoa_mumax(i,j,k) = max(0.2, 0.029 * Temp(i,j,k) - 0.147) / 86400.0
      ! 2. Limitation terms of oxygen and ammonium substrate affecting uptake
      wombat%aoa_loxy(i,j,k) = min(1.0, oxy_mmolm3 * wombat%aoa_poxy)
      wombat%aoa_lnh4(i,j,k) = nh4_mmolm3 / (nh4_mmolm3 + wombat%aoa_knh4)
      aoa_Voxy = oxy_mmolm3 * wombat%aoa_poxy
      aoa_Vnh4 = wombat%aoa_ynh4 * wombat%aoa_mumax(i,j,k) * wombat%aoa_lnh4(i,j,k)  ! Note: yield * max growth rate = Vmax
      ! 3. Redefine growth rate based on these limitations
      wombat%aoa_mu(i,j,k) = min( (aoa_Voxy/wombat%aoa_yoxy), (aoa_Vnh4/wombat%aoa_ynh4) )
      wombat%aoa_eno3(i,j,k) = wombat%aoa_ynh4 - 1.0/wombat%aoa_C2N

      if (do_anammox) then
        ! Anaerobic ammonium oxidation (anammox)
        wombat%aox_lnh4(i,j,k) = nh4_mmolm3 / (nh4_mmolm3 + wombat%aoxkn)
        wombat%aox_mu(i,j,k) = wombat%aoxmumax * wombat%bbioh**(Temp(i,j,k)) &
                               * wombat%lbac_fanaer(i,j,k) * wombat%aox_lnh4(i,j,k)
      endif

      ! Chemoautotrophy
      if (aoa_p > epsi) then
        wombat%aoagrow(i,j,k) = wombat%aoa_mu(i,j,k) * aoa_p ! [molC/kg/s]
        wombat%ammox(i,j,k) = wombat%aoagrow(i,j,k) * wombat%aoa_ynh4 ! [molNH4/kg/s]
        wombat%aoaresp(i,j,k) = wombat%aoagrow(i,j,k) * wombat%aoa_yoxy ! [molO2/kg/s]
      else
        wombat%aoagrow(i,j,k) = 0.0
        wombat%ammox(i,j,k) = 0.0
        wombat%aoaresp(i,j,k) = 0.0
      endif

      if (nh4_p > epsi) then
        wombat%anammox(i,j,k) = wombat%aox_mu(i,j,k) * nh4_p ! [molN/kg/s]
      else
        wombat%anammox(i,j,k) = 0.0
      endif


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 18] Tracer tendencies                                          !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Nitrate equation ! [molN/kg]
      !----------------------------------------------------------------------
      wombat%f_no3(i,j,k) = wombat%f_no3(i,j,k) + dtsb * ( &
                            wombat%aoagrow(i,j,k) * wombat%aoa_eno3(i,j,k) &
                          - wombat%lbacdeni(i,j,k) - wombat%sbacdeni(i,j,k) ) &
                          - dtsb * 16./122. * ( &
                            wombat%phygrow(i,j,k) * wombat%phy_lno3(i,j,k) / ( wombat%phy_lnit(i,j,k) + epsi ) &
                          + wombat%diagrow(i,j,k) * wombat%dia_lno3(i,j,k) / ( wombat%dia_lnit(i,j,k) + epsi ) )

      ! Ammonium equation ! [molN/kg]
      !----------------------------------------------------------------------
      wombat%f_nh4(i,j,k) = wombat%f_nh4(i,j,k) + dtsb * ( &
                            ( zooexcrlbacn &
                            + zooexcrsbacn &
                            + zooexcraoan ) * (1.0-wombat%zooexcrdom) &
                          + ( mesexcrlbacn &
                            + mesexcrsbacn &
                            + mesexcraoan ) * (1.0-wombat%mesexcrdom) &
                          + wombat%lbacpnh4(i,j,k) &
                          + wombat%sbacpnh4(i,j,k) &
                          - wombat%ammox(i,j,k) &
                          - wombat%anammox(i,j,k) ) + dtsb * 16./122. * ( &
                            wombat%zoomorl(i,j,k) &
                          + wombat%mesmorl(i,j,k) &
                          + ( wombat%zooexcrphy(i,j,k) &
                            + wombat%zooexcrdia(i,j,k) &
                            + wombat%zooexcrsdet(i,j,k) ) * (1.0-wombat%zooexcrdom) &
                          + ( wombat%mesexcrphy(i,j,k) &
                            + wombat%mesexcrdia(i,j,k) &
                            + wombat%mesexcrsdet(i,j,k) &
                            + wombat%mesexcrldet(i,j,k) &
                            + wombat%mesexcrzoo(i,j,k) ) * (1.0-wombat%mesexcrdom) &
                          - wombat%phygrow(i,j,k) * wombat%phy_lnh4(i,j,k) / ( wombat%phy_lnit(i,j,k) + epsi ) &
                          - wombat%diagrow(i,j,k) * wombat%dia_lnh4(i,j,k) / ( wombat%dia_lnit(i,j,k) + epsi ) )
      if (do_nitrogen_fixation) &
        wombat%f_nh4(i,j,k) = wombat%f_nh4(i,j,k) + dtsb * wombat%nitrfix(i,j,k)

      ! Silicic acid equation ! [molSi/kg]
      !   Microzooplankton grazing on diatoms produces clean, small, largely suspended
      !   pieces of biogenic silica prone to rapid dissolution [Krause et al., 2010 L&O]
      !----------------------------------------------------------------------
      wombat%f_sil(i,j,k) = wombat%f_sil(i,j,k) + dtsb * ( &
                            wombat%diamorl(i,j,k) * dia_Si2C &
                          - wombat%dia_silupt(i,j,k) &
                          + wombat%zoograzdia(i,j,k) * dia_Si2C &
                          + wombat%bsidiss(i,j,k) )

      ! Phytoplankton equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_phy(i,j,k) = wombat%f_phy(i,j,k) + dtsb * ( &
                            wombat%phygrow(i,j,k) &
                          - wombat%phymorl(i,j,k) &
                          - wombat%phymorq(i,j,k) &
                          - wombat%zoograzphy(i,j,k) &
                          - wombat%mesgrazphy(i,j,k) )

      ! Phytoplankton chlorophyll equation ! [molChl/kg]
      !-----------------------------------------------------------------------
      wombat%f_pchl(i,j,k) = wombat%f_pchl(i,j,k) + dtsb * ( &
                             wombat%pchl_mu(i,j,k) &
                           - ( wombat%phymorl(i,j,k) &
                             + wombat%phymorq(i,j,k) &
                             + wombat%zoograzphy(i,j,k) &
                             + wombat%mesgrazphy(i,j,k) ) * phy_chlc )

      ! Phytoplankton iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_phyfe(i,j,k) = wombat%f_phyfe(i,j,k) + dtsb * ( &
                              wombat%phy_dfeupt(i,j,k) &
                            - ( wombat%phymorl(i,j,k) &
                              + wombat%phymorq(i,j,k) &
                              + wombat%zoograzphy(i,j,k) &
                              + wombat%mesgrazphy(i,j,k) ) * phy_Fe2C )

      ! Microphytoplankton equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_dia(i,j,k)  = wombat%f_dia(i,j,k) + dtsb * ( &
                             wombat%diagrow(i,j,k) &
                           - wombat%diamorl(i,j,k) &
                           - wombat%diamorq(i,j,k) &
                           - wombat%mesgrazdia(i,j,k) &
                           - wombat%zoograzdia(i,j,k) )

      ! Microphytoplankton chlorophyll equation ! [molChl/kg]
      !-----------------------------------------------------------------------
      wombat%f_dchl(i,j,k) = wombat%f_dchl(i,j,k) + dtsb * ( &
                             wombat%dchl_mu(i,j,k) &
                           - ( wombat%diamorl(i,j,k) &
                             + wombat%diamorq(i,j,k) &
                             + wombat%zoograzdia(i,j,k) &
                             + wombat%mesgrazdia(i,j,k) ) * dia_chlc )

      ! Microphytoplankton iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_diafe(i,j,k) = wombat%f_diafe(i,j,k) + dtsb * ( &
                              wombat%dia_dfeupt(i,j,k) &
                            - ( wombat%diamorl(i,j,k) &
                              + wombat%diamorq(i,j,k) &
                              + wombat%zoograzdia(i,j,k) &
                              + wombat%mesgrazdia(i,j,k) ) * dia_Fe2C )

      ! Microphytoplankton silicon equation ! [molSi/kg]
      !----------------------------------------------------------------------
      wombat%f_diasi(i,j,k) = wombat%f_diasi(i,j,k) + dtsb * ( &
                              wombat%dia_silupt(i,j,k) &
                            - ( wombat%diamorl(i,j,k) &
                              + wombat%diamorq(i,j,k) &
                              + wombat%mesgrazdia(i,j,k) &
                              + wombat%zoograzdia(i,j,k) ) * dia_Si2C )

      ! Estimate primary productivity from phytoplankton growth ! [molC/kg/s]
      wombat%rpp3d(i,j,k) = wombat%rpp3d(i,j,k) + dtsb * ( &
                              wombat%phygrow(i,j,k) * wombat%phy_lnh4(i,j,k) / ( wombat%phy_lnit(i,j,k) + epsi ) &
                            + wombat%diagrow(i,j,k) * wombat%dia_lnh4(i,j,k) / ( wombat%dia_lnit(i,j,k) + epsi ) )

      ! Net primary productivity (gross PP minus linear mortality) ! [molC/kg/s]
      wombat%npp3d(i,j,k) = wombat%npp3d(i,j,k) + dtsb * ( &
                            wombat%phygrow(i,j,k) + wombat%diagrow(i,j,k) )

      ! Zooplankton equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_zoo(i,j,k) = wombat%f_zoo(i,j,k) + dtsb * ( &
                            ( wombat%zoograzlbac(i,j,k) &
                            + wombat%zoograzsbac(i,j,k) &
                            + wombat%zoograzaoa(i,j,k) &
                            + wombat%zoograzphy(i,j,k) &
                            + wombat%zoograzdia(i,j,k) &
                            + wombat%zoograzsdet(i,j,k) ) * wombat%zooCingest*wombat%zooCassim &
                          - wombat%mesgrazzoo(i,j,k) &
                          - wombat%zoomorl(i,j,k) &
                          - wombat%zoomorq(i,j,k) )

      ! Zooplankton iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_zoofe(i,j,k) = wombat%f_zoofe(i,j,k) + dtsb * ( &
                              zooassilbacfe &
                            + zooassisbacfe &
                            + zooassiaoafe &
                            + zooassiphyfe &
                            + zooassidiafe &
                            + zooassisdetfe &
                            - ( wombat%mesgrazzoo(i,j,k) &
                              + wombat%zoomorl(i,j,k) &
                              + wombat%zoomorq(i,j,k) ) * zoo_Fe2C )

      ! Mesozooplankton equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_mes(i,j,k) = wombat%f_mes(i,j,k) + dtsb * ( &
                            ( wombat%mesgrazlbac(i,j,k) &
                            + wombat%mesgrazsbac(i,j,k) &
                            + wombat%mesgrazaoa(i,j,k) &
                            + wombat%mesgrazphy(i,j,k) &
                            + wombat%mesgrazdia(i,j,k) &
                            + wombat%mesgrazsdet(i,j,k) &
                            + wombat%mesgrazldet(i,j,k) &
                            + wombat%mesgrazzoo(i,j,k) ) * wombat%mesCingest*wombat%mesCassim &
                           - wombat%mesmorl(i,j,k) &
                           - wombat%mesmorq(i,j,k) )

      ! Mesozooplankton iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_mesfe(i,j,k) = wombat%f_mesfe(i,j,k) + dtsb * ( &
                              mesassilbacfe &
                            + mesassisbacfe &
                            + mesassiaoafe &
                            + mesassiphyfe &
                            + mesassidiafe &
                            + mesassisdetfe &
                            + mesassildetfe &
                            + mesassizoofe &
                            - ( wombat%mesmorl(i,j,k) &
                              + wombat%mesmorq(i,j,k) ) * mes_Fe2C )

      ! Estimate secondary productivity from zooplankton growth ! [molC/kg/s]
      wombat%zsp3d(i,j,k) = wombat%zsp3d(i,j,k) + dtsb * ( &
                            ( wombat%zoograzlbac(i,j,k) &
                            + wombat%zoograzsbac(i,j,k) &
                            + wombat%zoograzaoa(i,j,k) &
                            + wombat%zoograzphy(i,j,k) &
                            + wombat%zoograzdia(i,j,k) &
                            + wombat%zoograzsdet(i,j,k) ) * wombat%zooCingest*wombat%zooCassim &
                          + ( wombat%mesgrazlbac(i,j,k) &
                            + wombat%mesgrazsbac(i,j,k) &
                            + wombat%mesgrazaoa(i,j,k) &
                            + wombat%mesgrazphy(i,j,k) &
                            + wombat%mesgrazdia(i,j,k) &
                            + wombat%mesgrazsdet(i,j,k) &
                            + wombat%mesgrazldet(i,j,k) &
                            + wombat%mesgrazzoo(i,j,k) ) * wombat%mesCingest*wombat%mesCassim )

      ! Small detritus equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_sdet(i,j,k) = wombat%f_sdet(i,j,k) + dtsb * ( &
                            wombat%zooegeslbac(i,j,k) &
                          + wombat%zooegessbac(i,j,k) &
                          + wombat%zooegesaoa(i,j,k) &
                          + wombat%zooegesphy(i,j,k) &
                          + wombat%zooegesdia(i,j,k) &
                          + wombat%zooegessdet(i,j,k) &
                          + wombat%phymorq(i,j,k) &
                          + wombat%zoomorq(i,j,k) &
                          - wombat%zoograzsdet(i,j,k) &
                          - wombat%mesgrazsdet(i,j,k) &
                          - wombat%sdetremi(i,j,k) )

      ! Small detrital iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_sdetfe(i,j,k) = wombat%f_sdetfe(i,j,k) + dtsb * ( &
                              zooegeslbacfe &
                            + zooegessbacfe &
                            + zooegesaoafe &
                            + zooegesphyfe &
                            + zooegesdiafe &
                            + zooegessdetfe &
                            + wombat%phymorq(i,j,k) * phy_Fe2C &
                            + wombat%zoomorq(i,j,k) * zoo_Fe2C &
                            - ( wombat%zoograzsdet(i,j,k) &
                              + wombat%mesgrazsdet(i,j,k) &
                              + wombat%sdetremi(i,j,k) ) * sdet_Fe2C )

      ! Large detritus equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_ldet(i,j,k) = wombat%f_ldet(i,j,k) + dtsb * ( &
                             wombat%mesegeslbac(i,j,k) &
                           + wombat%mesegessbac(i,j,k) &
                           + wombat%mesegesaoa(i,j,k) &
                           + wombat%mesegesphy(i,j,k) &
                           + wombat%mesegesdia(i,j,k) &
                           + wombat%mesegessdet(i,j,k) &
                           + wombat%mesegesldet(i,j,k) &
                           + wombat%mesegeszoo(i,j,k) &
                           + wombat%diamorq(i,j,k) &
                           + wombat%mesmorq(i,j,k) &
                           - wombat%mesgrazldet(i,j,k) &
                           - wombat%ldetremi(i,j,k) )

      ! Compact, fast sinking detrital iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_ldetfe(i,j,k) = wombat%f_ldetfe(i,j,k) + dtsb * ( &
                               mesegeslbacfe &
                             + mesegessbacfe &
                             + mesegesaoafe &
                             + mesegesphyfe &
                             + mesegesdiafe &
                             + mesegessdetfe &
                             + mesegesldetfe &
                             + mesegeszoofe &
                             + wombat%diamorq(i,j,k) * dia_Fe2C &
                             + wombat%mesmorq(i,j,k) * mes_Fe2C &
                             - ( wombat%mesgrazldet(i,j,k) &
                               + wombat%ldetremi(i,j,k) ) * ldet_Fe2C )

      ! Compact, fast sinking detrital silicon equation ! [molSi/kg]
      !   Copepod egestion (fecal pellets) represented 42-107% of biogenic silica export at
      !   100 metres in the spring bloom at the Antarctic Polar Front [Dagg et al., 2003 DSRII]
      !   All mesozooplankton grazing on diatoms goes to egestion, no dissolution and no assimilation
      !----------------------------------------------------------------------
      wombat%f_ldetsi(i,j,k) = wombat%f_ldetsi(i,j,k) + dtsb * ( &
                               ( wombat%diamorq(i,j,k) &
                               + wombat%mesgrazdia(i,j,k) ) * dia_Si2C &
                             - wombat%bsidiss(i,j,k)  )

      ! Large dissolved organic carbon equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_ldoc(i,j,k) = wombat%f_ldoc(i,j,k) + dtsb * ( &
                             wombat%sdetremi(i,j,k) &
                           + wombat%ldetremi(i,j,k) &
                           + wombat%phymorl(i,j,k) &
                           + wombat%diamorl(i,j,k) &
                           + wombat%lbacmorl(i,j,k) &
                           + wombat%lbacmorq(i,j,k) &
                           + wombat%sbacmorl(i,j,k) &
                           + wombat%sbacmorq(i,j,k) &
                           + wombat%aoamorl(i,j,k) &
                           + wombat%aoamorq(i,j,k) &
                           + ( wombat%zooexcrlbac(i,j,k) &
                             + wombat%zooexcrsbac(i,j,k) &
                             + wombat%zooexcraoa(i,j,k) &
                             + wombat%zooexcrphy(i,j,k) &
                             + wombat%zooexcrdia(i,j,k) &
                             + wombat%zooexcrsdet(i,j,k) ) * wombat%zooexcrdom &
                           + ( wombat%mesexcrlbac(i,j,k) &
                             + wombat%mesexcrsbac(i,j,k) &
                             + wombat%mesexcraoa(i,j,k) &
                             + wombat%mesexcrphy(i,j,k) &
                             + wombat%mesexcrdia(i,j,k) &
                             + wombat%mesexcrsdet(i,j,k) &
                             + wombat%mesexcrldet(i,j,k) &
                             + wombat%mesexcrzoo(i,j,k) ) * wombat%mesexcrdom &
                           - wombat%ldocremi(i,j,k) )

      ! Large dissolved organic nitrogen equation ! [molN/kg]
      !-----------------------------------------------------------------------
      wombat%f_ldon(i,j,k) = wombat%f_ldon(i,j,k) + dtsb * ( &
                            ( wombat%lbacmorl(i,j,k) &
                            + wombat%lbacmorq(i,j,k) &
                            + wombat%sbacmorl(i,j,k) &
                            + wombat%sbacmorq(i,j,k) ) / wombat%bac_C2N &
                          + ( wombat%aoamorl(i,j,k) &
                            + wombat%aoamorq(i,j,k) ) / wombat%aoa_C2N &
                          + ( zooexcrlbacn &
                            + zooexcrsbacn &
                            + zooexcraoan ) * wombat%zooexcrdom &
                          + ( mesexcrlbacn &
                            + mesexcrsbacn &
                            + mesexcraoan ) * wombat%mesexcrdom &
                          - wombat%ldocremi(i,j,k) * ldom_N2C ) &
                          + dtsb * 16./122.0 * ( &
                            wombat%sdetremi(i,j,k) &
                          + wombat%ldetremi(i,j,k) &
                          + wombat%phymorl(i,j,k) &
                          + wombat%diamorl(i,j,k) &
                          + ( wombat%zooexcrphy(i,j,k) &
                            + wombat%zooexcrdia(i,j,k) &
                            + wombat%zooexcrsdet(i,j,k) ) *wombat%zooexcrdom &
                          + ( wombat%mesexcrphy(i,j,k) &
                            + wombat%mesexcrdia(i,j,k) &
                            + wombat%mesexcrsdet(i,j,k) &
                            + wombat%mesexcrldet(i,j,k) &
                            + wombat%mesexcrzoo(i,j,k) ) * wombat%mesexcrdom )

      ! Small dissolved organic carbon equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_sdoc(i,j,k) = wombat%f_sdoc(i,j,k) + dtsb * ( &
                             wombat%phydoc(i,j,k) &
                           + wombat%diadoc(i,j,k) &
                           + wombat%sdocprod(i,j,k) &
                           - wombat%sdocremi(i,j,k) )

      ! Large, sharing heterotrophic bacteria ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_lbac(i,j,k) = wombat%f_lbac(i,j,k) + dtsb * ( &
                             wombat%lbacgrow(i,j,k) &
                           - wombat%zoograzlbac(i,j,k) &
                           - wombat%mesgrazlbac(i,j,k) &
                           - wombat%lbacmorl(i,j,k) &
                           - wombat%lbacmorq(i,j,k) )

      ! Small, selfish heterotrophic bacteria ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_sbac(i,j,k) = wombat%f_sbac(i,j,k) + dtsb * ( &
                             wombat%sbacgrow(i,j,k) &
                           - wombat%zoograzsbac(i,j,k) &
                           - wombat%mesgrazsbac(i,j,k) &
                           - wombat%sbacmorl(i,j,k) &
                           - wombat%sbacmorq(i,j,k) )

      ! AOA ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_aoa(i,j,k) = wombat%f_aoa(i,j,k) + dtsb * ( &
                            wombat%aoagrow(i,j,k) &
                          - wombat%zoograzaoa(i,j,k) &
                          - wombat%mesgrazaoa(i,j,k) &
                          - wombat%aoamorl(i,j,k) &
                          - wombat%aoamorq(i,j,k) )

      ! Oxygen equation ! [molO2/kg]
      !-----------------------------------------------------------------------
      if (wombat%f_o2(i,j,k) > epsi) &
        wombat%f_o2(i,j,k) = wombat%f_o2(i,j,k) - 132./122. * dtsb * ( &
                             ( wombat%zooexcrlbac(i,j,k) &
                             + wombat%zooexcrsbac(i,j,k) &
                             + wombat%zooexcraoa(i,j,k) &
                             + wombat%zooexcrphy(i,j,k) &
                             + wombat%zooexcrdia(i,j,k) &
                             + wombat%zooexcrsdet(i,j,k) ) * (1.0-wombat%zooexcrdom) &
                           + ( wombat%mesexcrlbac(i,j,k) &
                             + wombat%mesexcrsbac(i,j,k) &
                             + wombat%mesexcraoa(i,j,k) &
                             + wombat%mesexcrphy(i,j,k) &
                             + wombat%mesexcrdia(i,j,k) &
                             + wombat%mesexcrsdet(i,j,k) &
                             + wombat%mesexcrldet(i,j,k) &
                             + wombat%mesexcrzoo(i,j,k) ) * (1.0-wombat%mesexcrdom) &
                           + wombat%zoomorl(i,j,k) &
                           + wombat%mesmorl(i,j,k) &
                           - wombat%phygrow(i,j,k) &
                           - wombat%diagrow(i,j,k) ) &
                           - dtsb * ( &
                             wombat%lbacresp(i,j,k) &
                           + wombat%sbacresp(i,j,k) &
                           + wombat%aoaresp(i,j,k) )


      ! Equation for CaCO3 ! [molCaCO3/kg]
      !-----------------------------------------------------------------------
      wombat%f_caco3(i,j,k) = wombat%f_caco3(i,j,k) + dtsb * ( &
                              ( ( wombat%zoograzphy(i,j,k) &
                                + wombat%mesgrazphy(i,j,k) &
                                + wombat%mesgrazzoo(i,j,k) ) * (1.0-wombat%fgutdiss) &
                              + wombat%phymorq(i,j,k)  &
                              + wombat%zoomorq(i,j,k) ) * wombat%pic2poc(i,j,k) &
                            - wombat%zoodiss(i,j,k) &
                            - wombat%mesdiss(i,j,k) &
                            - wombat%caldiss(i,j,k) &
                            - wombat%aradiss(i,j,k) &
                            - wombat%pocdiss(i,j,k) )

      ! Equation for DIC ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_dic(i,j,k) = wombat%f_dic(i,j,k) + dtsb * ( &
                            ( wombat%zooexcrlbac(i,j,k) &
                            + wombat%zooexcrsbac(i,j,k) &
                            + wombat%zooexcraoa(i,j,k) &
                            + wombat%zooexcrphy(i,j,k) &
                            + wombat%zooexcrdia(i,j,k) &
                            + wombat%zooexcrsdet(i,j,k) ) * (1.0-wombat%zooexcrdom) &
                          + ( wombat%mesexcrlbac(i,j,k) &
                            + wombat%mesexcrsbac(i,j,k) &
                            + wombat%mesexcraoa(i,j,k) &
                            + wombat%mesexcrphy(i,j,k) &
                            + wombat%mesexcrdia(i,j,k) &
                            + wombat%mesexcrsdet(i,j,k) &
                            + wombat%mesexcrldet(i,j,k) &
                            + wombat%mesexcrzoo(i,j,k) ) * (1.0-wombat%mesexcrdom) &
                          + wombat%zoomorl(i,j,k) &
                          + wombat%mesmorl(i,j,k) &
                          + wombat%lbacpco2(i,j,k) &
                          + wombat%sbacpco2(i,j,k) &
                          + wombat%zoodiss(i,j,k) &
                          + wombat%mesdiss(i,j,k) &
                          + wombat%caldiss(i,j,k) &
                          + wombat%aradiss(i,j,k) &
                          + wombat%pocdiss(i,j,k) &
                          - wombat%phygrow(i,j,k) &
                          - wombat%diagrow(i,j,k) &
                          - wombat%aoagrow(i,j,k) &
                          - wombat%phydoc(i,j,k) &
                          - wombat%diadoc(i,j,k) &
                          - ( ( wombat%zoograzphy(i,j,k) &
                              + wombat%mesgrazphy(i,j,k) &
                              + wombat%mesgrazzoo(i,j,k) ) * (1.0-wombat%fgutdiss) &
                              + wombat%phymorq(i,j,k) &
                              + wombat%zoomorq(i,j,k) ) * wombat%pic2poc(i,j,k) )

      ! Equation for DICr ! [molC/kg]
      !-----------------------------------------------------------------------
      if (do_tracer_dicr) then
        wombat%f_dicr(i,j,k) = wombat%f_dicr(i,j,k) + dtsb * ( &
                                ( wombat%zooexcrlbac(i,j,k) &
                                + wombat%zooexcrsbac(i,j,k) &
                                + wombat%zooexcraoa(i,j,k) &
                                + wombat%zooexcrphy(i,j,k) &
                                + wombat%zooexcrdia(i,j,k) &
                                + wombat%zooexcrsdet(i,j,k) ) * (1.0-wombat%zooexcrdom) &
                            + ( wombat%mesexcrlbac(i,j,k) &
                              + wombat%mesexcrsbac(i,j,k) &
                              + wombat%mesexcraoa(i,j,k) &
                              + wombat%mesexcrphy(i,j,k) &
                              + wombat%mesexcrdia(i,j,k) &
                              + wombat%mesexcrsdet(i,j,k) &
                              + wombat%mesexcrldet(i,j,k) &
                              + wombat%mesexcrzoo(i,j,k) ) * (1.0-wombat%mesexcrdom) &
                            + wombat%zoomorl(i,j,k) &
                            + wombat%mesmorl(i,j,k) &
                            + wombat%lbacpco2(i,j,k) &
                            + wombat%sbacpco2(i,j,k) &
                            + wombat%zoodiss(i,j,k) &
                            + wombat%mesdiss(i,j,k) &
                            + wombat%caldiss(i,j,k) &
                            + wombat%aradiss(i,j,k) &
                            + wombat%pocdiss(i,j,k) &
                            - wombat%phygrow(i,j,k) &
                            - wombat%diagrow(i,j,k) &
                            - wombat%aoagrow(i,j,k) &
                            - wombat%phydoc(i,j,k) &
                            - wombat%diadoc(i,j,k) &
                            - ( ( wombat%zoograzphy(i,j,k) &
                                + wombat%mesgrazphy(i,j,k) &
                                + wombat%mesgrazzoo(i,j,k) ) * (1.0-wombat%fgutdiss) &
                                + wombat%phymorq(i,j,k) &
                                + wombat%zoomorq(i,j,k) ) * wombat%pic2poc(i,j,k) )
      endif


      ! Equation for ALK ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_alk(i,j,k) = wombat%f_alk(i,j,k) + dtsb * 16.0/122.0 * ( &
                            wombat%phygrow(i,j,k) * wombat%phy_lno3(i,j,k) / ( wombat%phy_lnit(i,j,k) + epsi ) &
                          + wombat%diagrow(i,j,k) * wombat%dia_lno3(i,j,k) / ( wombat%dia_lnit(i,j,k) + epsi ) &
                          - wombat%phygrow(i,j,k) * wombat%phy_lnh4(i,j,k) / ( wombat%phy_lnit(i,j,k) + epsi ) &
                          - wombat%diagrow(i,j,k) * wombat%dia_lnh4(i,j,k) / ( wombat%dia_lnit(i,j,k) + epsi ) &
                          + wombat%zoomorl(i,j,k) &
                          + wombat%mesmorl(i,j,k) &
                          + ( wombat%zooexcrphy(i,j,k) &
                            + wombat%zooexcrdia(i,j,k) &
                            + wombat%zooexcrsdet(i,j,k) ) * (1.0-wombat%zooexcrdom) &
                          + ( wombat%mesexcrphy(i,j,k) &
                            + wombat%mesexcrdia(i,j,k) &
                            + wombat%mesexcrsdet(i,j,k) &
                            + wombat%mesexcrldet(i,j,k) &
                            + wombat%mesexcrzoo(i,j,k) ) * (1.0-wombat%mesexcrdom) ) &
                          + dtsb * ( &
                            wombat%lbacpnh4(i,j,k) &
                          + wombat%sbacpnh4(i,j,k) &
                          + ( zooexcrlbacn &
                            + zooexcrsbacn &
                            + zooexcraoan ) * (1.0-wombat%zooexcrdom) &
                          + ( mesexcrlbacn &
                            + mesexcrsbacn &
                            + mesexcraoan ) * (1.0-wombat%mesexcrdom) &
                          + wombat%lbacdeni(i,j,k) &
                          + wombat%sbacdeni(i,j,k) &
                          - 2.0 * wombat%ammox(i,j,k) + wombat%aoagrow(i,j,k)/wombat%aoa_C2N &
                          - wombat%anammox(i,j,k) ) &
                          + dtsb * 2.0 * ( &
                            wombat%zoodiss(i,j,k) &
                          + wombat%mesdiss(i,j,k) &
                          + wombat%caldiss(i,j,k) &
                          + wombat%aradiss(i,j,k) &
                          + wombat%pocdiss(i,j,k) &
                          - ( ( wombat%zoograzphy(i,j,k) &
                              + wombat%mesgrazphy(i,j,k) &
                              + wombat%mesgrazzoo(i,j,k) ) * (1.0-wombat%fgutdiss) &
                            + wombat%phymorq(i,j,k) &
                            + wombat%zoomorq(i,j,k) ) * wombat%pic2poc(i,j,k) )


      ! Equation for dissolved iron ! [molFe/kg]
      !----------------------------------------------------------------------
      wombat%f_fe(i,j,k) = wombat%f_fe(i,j,k) + dtsb * ( &
                           wombat%sdetremi(i,j,k) * sdet_Fe2C &
                         + wombat%ldetremi(i,j,k) * ldet_Fe2C &
                         + wombat%phymorl(i,j,k) * phy_Fe2C &
                         + wombat%diamorl(i,j,k) * dia_Fe2C &
                         + wombat%zoomorl(i,j,k) * zoo_Fe2C &
                         + wombat%mesmorl(i,j,k) * mes_Fe2C &
                         + ( wombat%lbacmorl(i,j,k) &
                           + wombat%lbacmorq(i,j,k) &
                           + wombat%sbacmorl(i,j,k) &
                           + wombat%sbacmorq(i,j,k) ) / wombat%bac_C2Fe &
                         + ( wombat%aoamorl(i,j,k) &
                           + wombat%aoamorq(i,j,k) ) / wombat%aoa_C2Fe &
                         + zooexcrlbacfe &
                         + zooexcrsbacfe &
                         + zooexcraoafe &
                         + zooexcrphyfe &
                         + zooexcrdiafe &
                         + zooexcrsdetfe &
                         + mesexcrlbacfe &
                         + mesexcrsbacfe &
                         + mesexcraoafe &
                         + mesexcrphyfe &
                         + mesexcrdiafe &
                         + mesexcrsdetfe &
                         + mesexcrldetfe &
                         + mesexcrzoofe &
                         + wombat%safediss(i,j,k) &
                         + wombat%lafediss(i,j,k) &
                         - wombat%phy_dfeupt(i,j,k) &
                         - wombat%dia_dfeupt(i,j,k) &
                         - wombat%lbacufer(i,j,k) &
                         - wombat%sbacufer(i,j,k) &
                         - wombat%aoagrow(i,j,k) / wombat%aoa_C2Fe &
                         - wombat%fescasafe(i,j,k) &
                         - wombat%fescalafe(i,j,k) &
                         - wombat%fecoag2safe(i,j,k) &
                         - wombat%fecoag2lafe(i,j,k) )

      ! Collect dFe sources and sinks for diagnostic output
      wombat%fesources(i,j,k) = wombat%fesources(i,j,k) + dtsb * ( &
                                  wombat%sdetremi(i,j,k) * sdet_Fe2C &
                                + wombat%ldetremi(i,j,k) * ldet_Fe2C &
                                + wombat%zoomorl(i,j,k) * zoo_Fe2C &
                                + wombat%mesmorl(i,j,k) * mes_Fe2C &
                                + wombat%lbacmorl(i,j,k) / wombat%bac_C2Fe &
                                + wombat%lbacmorq(i,j,k) / wombat%bac_C2Fe &
                                + wombat%sbacmorl(i,j,k) / wombat%bac_C2Fe &
                                + wombat%sbacmorq(i,j,k) / wombat%bac_C2Fe &
                                + wombat%aoamorl(i,j,k) / wombat%aoa_C2Fe &
                                + wombat%aoamorq(i,j,k) / wombat%aoa_C2Fe &
                                + zooexcrlbacfe &
                                + zooexcrsbacfe &
                                + zooexcraoafe &
                                + zooexcrphyfe &
                                + zooexcrdiafe &
                                + zooexcrsdetfe &
                                + mesexcrlbacfe &
                                + mesexcrsbacfe &
                                + mesexcraoafe &
                                + mesexcrphyfe &
                                + mesexcrdiafe &
                                + mesexcrsdetfe &
                                + mesexcrldetfe &
                                + mesexcrzoofe &
                                + wombat%phymorl(i,j,k) * phy_Fe2C &
                                + wombat%diamorl(i,j,k) * dia_Fe2C &
                                + wombat%safediss(i,j,k) &
                                + wombat%lafediss(i,j,k))
      wombat%fesinks(i,j,k) = wombat%fesinks(i,j,k) + dtsb * ( &
                                wombat%phy_dfeupt(i,j,k) &
                              + wombat%dia_dfeupt(i,j,k) &
                              + wombat%lbacufer(i,j,k) &
                              + wombat%sbacufer(i,j,k) &
                              + wombat%aoagrow(i,j,k) / wombat%aoa_C2Fe &
                              + wombat%fescasafe(i,j,k) &
                              + wombat%fescalafe(i,j,k) &
                              + wombat%fecoag2safe(i,j,k) &
                              + wombat%fecoag2lafe(i,j,k))

      ! Equation for small authigenic iron (oxyhydroxide) ! [molFe/kg]
      !----------------------------------------------------------------------
      wombat%f_safe(i,j,k) = wombat%f_safe(i,j,k) + dtsb * ( &
                              wombat%fecoag2safe(i,j,k) &
                            + wombat%fescasafe(i,j,k) &
                            - wombat%safediss(i,j,k) )

      ! Equation for large authigenic iron (oxyhydroxide) ! [molFe/kg]
      !----------------------------------------------------------------------
      wombat%f_lafe(i,j,k) = wombat%f_lafe(i,j,k) + dtsb * ( &
                               wombat%fecoag2lafe(i,j,k) &
                             + wombat%fescalafe(i,j,k) &
                             - wombat%lafediss(i,j,k) )


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 19] Check for conservation of mass by ecosystem component      !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      n_pools(i,j,k,2) = wombat%f_no3(i,j,k) + wombat%f_nh4(i,j,k) + wombat%f_ldon(i,j,k) + &
                         ( wombat%f_phy(i,j,k) + wombat%f_sdet(i,j,k) + wombat%f_ldet(i,j,k) + &
                         wombat%f_zoo(i,j,k) + wombat%f_mes(i,j,k) + wombat%f_dia(i,j,k) ) * 16/122.0 + &
                         ( wombat%f_lbac(i,j,k) / wombat%bac_C2N + wombat%f_sbac(i,j,k) / wombat%bac_C2N + &
                           wombat%f_aoa(i,j,k) / wombat%aoa_C2N )
      c_pools(i,j,k,2) = wombat%f_dic(i,j,k) + wombat%f_phy(i,j,k) + wombat%f_sdet(i,j,k) + wombat%f_ldet(i,j,k) + &
                         wombat%f_zoo(i,j,k) + wombat%f_mes(i,j,k) + wombat%f_caco3(i,j,k) + wombat%f_dia(i,j,k) + &
                         wombat%f_ldoc(i,j,k) + wombat%f_sdoc(i,j,k) + wombat%f_lbac(i,j,k) + wombat%f_sbac(i,j,k) + &
                         wombat%f_aoa(i,j,k)
      si_pools(i,j,k,2) = wombat%f_sil(i,j,k) + wombat%f_diasi(i,j,k) + wombat%f_ldetsi(i,j,k)
      fe_pools(i,j,k,2) = wombat%f_fe(i,j,k) + wombat%f_safe(i,j,k) + wombat%f_lafe(i,j,k) + &
                          wombat%f_phyfe(i,j,k) + wombat%f_diafe(i,j,k) + wombat%f_zoofe(i,j,k) + wombat%f_mesfe(i,j,k) + &
                          wombat%f_sdetfe(i,j,k) + wombat%f_ldetfe(i,j,k) + wombat%f_lbac(i,j,k) / wombat%bac_C2Fe + &
                          wombat%f_sbac(i,j,k) / wombat%bac_C2Fe + wombat%f_aoa(i,j,k) / wombat%aoa_C2Fe


      if (tn>1) then
        if (do_check_n_conserve) then
          if (abs(n_pools(i,j,k,2) - n_pools(i,j,k,1))>1e-16) then
            print *, "--------------------------------------------"
            print *, trim(error_header) // " Ecosystem model is not conserving nitrogen"
            print *, "       Longitude index =", i
            print *, "       Latitude index =", j
            print *, "       Depth index and value =", k, wombat%zm(i,j,k)
            print *, "       Nested timestep number =", tn
            print *, " "
            print *, "       Biological N budget (molN/kg) at two timesteps =", n_pools(i,j,k,1), n_pools(i,j,k,2)
            print *, "       Difference in budget between timesteps =", n_pools(i,j,k,2) - n_pools(i,j,k,1)
            print *, " "
            print *, "       NO3 (molNO3/kg) =", wombat%f_no3(i,j,k)
            print *, "       NH4 (molNH4/kg) =", wombat%f_nh4(i,j,k)
            print *, "       PHY (molN/kg) =", wombat%f_phy(i,j,k) * 16.0 / 122.0
            print *, "       DIA (molN/kg) =", wombat%f_dia(i,j,k) * 16.0 / 122.0
            print *, "       ZOO (molN/kg) =", wombat%f_zoo(i,j,k) * 16.0 / 122.0
            print *, "       MES (molN/kg) =", wombat%f_mes(i,j,k) * 16.0 / 122.0
            print *, "       sDET (molN/kg) =", wombat%f_sdet(i,j,k) * 16.0 / 122.0
            print *, "       lDET (molN/kg) =", wombat%f_ldet(i,j,k) * 16.0 / 122.0
            print *, "       lDON (molN/kg) =", wombat%f_ldon(i,j,k)
            print *, "       lBAC (molN/kg) =", wombat%f_lbac(i,j,k) / wombat%bac_C2N
            print *, "       sBAC (molN/kg) =", wombat%f_sbac(i,j,k) / wombat%bac_C2N
            print *, "       AOA (molN/kg) =", wombat%f_aoa(i,j,k) / wombat%aoa_C2N
            print *, " "
            print *, "--------------------------------------------"
            call mpp_error(FATAL, trim(error_header) // " Terminating run due to non-conservation of tracer")
          endif
        endif
        if (do_check_c_conserve) then
          if (abs(c_pools(i,j,k,2) - c_pools(i,j,k,1))>1e-16) then
            print *, "--------------------------------------------"
            print *, trim(error_header) // " Ecosystem model is not conserving carbon"
            print *, "       Longitude index =", i
            print *, "       Latitude index =", j
            print *, "       Depth index and value =", k, wombat%zm(i,j,k)
            print *, "       Nested timestep number =", tn
            print *, " "
            print *, "       Biological C budget (molC/kg) at two timesteps =", c_pools(i,j,k,1), c_pools(i,j,k,2)
            print *, "       Difference in budget between timesteps =", c_pools(i,j,k,2) - c_pools(i,j,k,1)
            print *, " "
            print *, "       DIC (molC/kg) =", wombat%f_dic(i,j,k)
            print *, "       ALK (molC/kg) =", wombat%f_alk(i,j,k)
            print *, "       PHY (molC/kg) =", wombat%f_phy(i,j,k)
            print *, "       DIA (molC/kg) =", wombat%f_dia(i,j,k)
            print *, "       ZOO (molC/kg) =", wombat%f_zoo(i,j,k)
            print *, "       MES (molC/kg) =", wombat%f_mes(i,j,k)
            print *, "       sDET (molC/kg) =", wombat%f_sdet(i,j,k)
            print *, "       lDET (molC/kg) =", wombat%f_ldet(i,j,k)
            print *, "       lBAC (molC/kg) =", wombat%f_lbac(i,j,k)
            print *, "       sBAC (molC/kg) =", wombat%f_sbac(i,j,k)
            print *, "       AOA (molC/kg) =", wombat%f_aoa(i,j,k)
            print *, "       lDOC (molC/kg) =", wombat%f_ldoc(i,j,k)
            print *, "       sDOC (molC/kg) =", wombat%f_sdoc(i,j,k)
            print *, "       CaCO3 (molC/kg) =", wombat%f_caco3(i,j,k)
            print *, "       Temp =", Temp(i,j,k)
            print *, "       Salt =", Salt(i,j,k)
            print *, "       surface pCO2 =", wombat%pco2_csurf(i,j)
            print *, "       htotal =", wombat%htotal(i,j,k)
            print *, " "
            print *, "--------------------------------------------"
            call mpp_error(FATAL, trim(error_header) // " Terminating run due to non-conservation of tracer")
          endif
        endif
        if (do_check_si_conserve) then
          if (abs(si_pools(i,j,k,2) - si_pools(i,j,k,1))>1e-16) then
            print *, "--------------------------------------------"
            print *, trim(error_header) // " Ecosystem model is not conserving silicon"
            print *, "       Longitude index =", i
            print *, "       Latitude index =", j
            print *, "       Depth index and value =", k, wombat%zm(i,j,k)
            print *, "       Nested timestep number =", tn
            print *, " "
            print *, "       Biological Si budget (molSi/kg) at two timesteps =", si_pools(i,j,k,1), si_pools(i,j,k,2)
            print *, "       Difference in budget between timesteps =", si_pools(i,j,k,2) - si_pools(i,j,k,1)
            print *, " "
            print *, "       SIL (molSi/kg) =", wombat%f_sil(i,j,k)
            print *, "       DIASI (molSi/kg) =", wombat%f_diasi(i,j,k)
            print *, "       ldetSI (molSi/kg) =", wombat%f_ldetsi(i,j,k)
            print *, "       Temp =", Temp(i,j,k)
            print *, "       Salt =", Salt(i,j,k)
            print *, " "
            print *, "--------------------------------------------"
            call mpp_error(FATAL, trim(error_header) // " Terminating run due to non-conservation of tracer")
          endif
        endif
        if (do_check_fe_conserve) then
          if (abs(fe_pools(i,j,k,2) - fe_pools(i,j,k,1))>1e-16) then
            print *, "--------------------------------------------"
            print *, trim(error_header) // " Ecosystem model is not conserving iron"
            print *, "       Longitude index =", i
            print *, "       Latitude index =", j
            print *, "       Depth index and value =", k, wombat%zm(i,j,k)
            print *, "       Nested timestep number =", tn
            print *, " "
            print *, "       Fe budget (molFe/kg) at two timesteps =", fe_pools(i,j,k,1), fe_pools(i,j,k,2)
            print *, "       Difference in budget between timesteps =", fe_pools(i,j,k,2) - fe_pools(i,j,k,1)
            print *, " "
            print *, "       dFE (molFe/kg) =", wombat%f_fe(i,j,k)
            print *, "       sAFe (molFe/kg) =", wombat%f_safe(i,j,k)
            print *, "       lAFe (molFe/kg) =", wombat%f_lafe(i,j,k)
            print *, "       PHY-Fe (molFe/kg) =", wombat%f_phyfe(i,j,k)
            print *, "       DIA-Fe (molFe/kg) =", wombat%f_diafe(i,j,k)
            print *, "       ZOO-Fe (molFe/kg) =", wombat%f_zoofe(i,j,k)
            print *, "       MES-Fe (molFe/kg) =", wombat%f_mesfe(i,j,k)
            print *, "       sDET-Fe (molFe/kg) =", wombat%f_sdetfe(i,j,k)
            print *, "       ldet-Fe (molFe/kg) =", wombat%f_ldetfe(i,j,k)
            print *, "       lBAC (molFe/kg) =", wombat%f_lbac(i,j,k) / wombat%bac_C2Fe
            print *, "       sBAC (molFe/kg) =", wombat%f_sbac(i,j,k) / wombat%bac_C2Fe
            print *, "       AOA (molFe/kg) =", wombat%f_aoa(i,j,k) / wombat%aoa_C2Fe
            print *, " "
            print *, "--------------------------------------------"
            call mpp_error(FATAL, trim(error_header) // " Terminating run due to non-conservation of tracer")
          endif
        endif
      endif

      enddo; enddo; enddo
    enddo !} nested timestep (ts_npzd)

    ! Add biotically induced tendency to biotracers
    !-----------------------------------------------------------------------

    do k = 1,nk; do j = jsc,jec; do i = isc,iec;
      ! PJB: it's possible that to reduce costs we could use the total change in NO3, NH4 and CaCO3
      ! tracers to solve for the change in alkalinity without needing the large number of terms above

      wombat%npp3d(i,j,k) = rdtts * wombat%npp3d(i,j,k) * grid_tmask(i,j,k) ! [mol/kg/s]
      wombat%rpp3d(i,j,k) = rdtts * wombat%rpp3d(i,j,k) * grid_tmask(i,j,k) ! [mol/kg/s]
      wombat%zsp3d(i,j,k) = rdtts * wombat%zsp3d(i,j,k) * grid_tmask(i,j,k) ! [mol/kg/s]
      wombat%fesources(i,j,k) = rdtts * wombat%fesources(i,j,k) * grid_tmask(i,j,k) ! [mol/kg/s]
      wombat%fesinks(i,j,k) = rdtts * wombat%fesinks(i,j,k) * grid_tmask(i,j,k) ! [mol/kg/s]

    enddo; enddo; enddo


    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !  [Step 20] Additional operations on tracers                           !
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!

    do j = jsc,jec; do i = isc,iec;
      ! mac: bottom dFe fix to 1 nM when the water is <= 200 m deep.
      if (grid_kmt(i,j) > 0) then
        k = grid_kmt(i,j)
        if (wombat%zw(i,j,k) <= 200) wombat%f_fe(i,j,k)= umol_m3_to_mol_kg * 0.999 ! [mol/kg]
      endif
      do k = 1,nk
        ! pjb: tune minimum dissolved iron concentration to detection limit...
        !       this is essential for ensuring dFe is replenished in upper ocean and actually
        !       looks to be the secret of PISCES ability to replicate dFe limitation in the right places
        wombat%f_fe(i,j,k) = max(wombat%dfefloor * umol_m3_to_mol_kg, wombat%f_fe(i,j,k)) * grid_tmask(i,j,k)
      enddo
    enddo; enddo


    ! Set tracers values
    call g_tracer_set_values(tracer_list, 'no3', 'field', wombat%f_no3, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'nh4', 'field', wombat%f_nh4, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'sil', 'field', wombat%f_sil, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'phy', 'field', wombat%f_phy, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'pchl', 'field', wombat%f_pchl, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'phyfe', 'field', wombat%f_phyfe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'dia', 'field', wombat%f_dia, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'dchl', 'field', wombat%f_dchl, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'diafe', 'field', wombat%f_diafe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'diasi', 'field', wombat%f_diasi, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'zoo', 'field', wombat%f_zoo, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'zoofe', 'field', wombat%f_zoofe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'mes', 'field', wombat%f_mes, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'mesfe', 'field', wombat%f_mesfe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'sdet', 'field', wombat%f_sdet, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'sdetfe', 'field', wombat%f_sdetfe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'ldet', 'field', wombat%f_ldet, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'ldetfe', 'field', wombat%f_ldetfe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'ldetsi', 'field', wombat%f_ldetsi, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'sdoc', 'field', wombat%f_sdoc, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'ldoc', 'field', wombat%f_ldoc, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'ldon', 'field', wombat%f_ldon, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'lbac', 'field', wombat%f_lbac, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'sbac', 'field', wombat%f_sbac, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'aoa', 'field', wombat%f_aoa, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'o2', 'field', wombat%f_o2, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'caco3', 'field', wombat%f_caco3, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'fe', 'field', wombat%f_fe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'safe', 'field', wombat%f_safe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'lafe', 'field', wombat%f_lafe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'dic', 'field', wombat%f_dic, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'alk', 'field', wombat%f_alk, isd, jsd, ntau=tau)
    if (do_tracer_dicr) &
      call g_tracer_set_values(tracer_list, 'dicr', 'field', wombat%f_dicr, isd, jsd, ntau=tau)


    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !  [Step 21] Compute sinking rates of detrital pools                    !
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!

    call g_tracer_get_pointer(tracer_list, 'sdet', 'vmove', wombat%p_wsdet) ! [m/s]
    call g_tracer_get_pointer(tracer_list, 'sdetfe', 'vmove', wombat%p_wsdetfe) ! [m/s]
    call g_tracer_get_pointer(tracer_list, 'ldet', 'vmove', wombat%p_wldet) ! [m/s]
    call g_tracer_get_pointer(tracer_list, 'ldetfe', 'vmove', wombat%p_wldetfe) ! [m/s]
    call g_tracer_get_pointer(tracer_list, 'ldetsi', 'vmove', wombat%p_wldetsi) ! [m/s]
    call g_tracer_get_pointer(tracer_list, 'caco3', 'vmove', wombat%p_wcaco3) ! [m/s]
    call g_tracer_get_pointer(tracer_list, 'safe', 'vmove', wombat%p_wsafe) ! [m/s]
    call g_tracer_get_pointer(tracer_list, 'lafe', 'vmove', wombat%p_wlafe) ! [m/s]

    ! Variable sinking rates of organic detritus (positive for sinking when GOLDtridiag == .true.)
    !                                            (negative for sinking when IOWtridiag ==.true.)
    ! Note: sinking distances are limited in the vertdiff solver to prevent characteristics crossing within a timestep

    do j = jsc,jec; do i = isc,iec
      if (grid_kmt(i,j)>0) then

        ! 1. Find the average radius of particles in the community
        !    For phytoplankton, we use Wickman et al. (2024):
        !    - Volume = Biomass ** 0.65
        !    - Diameter = ( (6/pi) * Volume ) ** (1/3)
        !    - therefore, diameter = ( (6/pi) * Biomass ** 0.65 ) ** (1/3)
        phy_mmolm3 = max(epsi, wombat%f_phy(i,j,1) ) / mmol_m3_to_mol_kg  ![mmol/m3]
        dia_mmolm3 = max(epsi, wombat%f_dia(i,j,1) ) / mmol_m3_to_mol_kg  ![mmol/m3]
        zoo_mmolm3 = max(epsi, wombat%f_zoo(i,j,1) ) / mmol_m3_to_mol_kg  ![mmol/m3]
        mes_mmolm3 = max(epsi, wombat%f_mes(i,j,1) ) / mmol_m3_to_mol_kg  ![mmol/m3]
        !rad_phy = wombat%phyrad0 * 0.5 * (6.0/pi * phy_mmolm3**0.65)**0.33 * 1e-6 ! [m]
        !rad_dia = wombat%diarad0 * 0.5 * (6.0/pi * dia_mmolm3**0.65)**0.33 * 1e-6 ! [m]
        ! ... simplifies to:
        rad_phy = wombat%phyrad0 * phy_mmolm3**0.217 * 0.62e-6 ! [m]
        rad_dia = wombat%diarad0 * dia_mmolm3**0.217 * 0.62e-6 ! [m]
        !    For microzoo (protists) we use Menden-Deuer & Lessard (2000) exponent
        !      of 0.939 that relates carbon biomass to volume
        !rad_zoo = wombat%zoorad0 * 0.5 * (6.0/pi * zoo_mmolm3**1.065)**0.33 * 1e-6 ! [m]
        ! simplifies to:
        rad_zoo = wombat%zoorad0 * zoo_mmolm3**0.355 * 0.62e-6 ! [m]
        !    For mesozoo we use the length-weight relationships for copepods from
        !      Uye (1982) and Lehette & Hernandez-Leon (2009)) of 1/3
        !rad_mes = wombat%mesrad0 * 0.5 * mes_mmolm3**0.33 * 1e-6 ! [m]
        rad_mes = wombat%mesrad0 * mes_mmolm3**0.33 * 0.5e-6 ! [m]
        if (phy_mmolm3 + zoo_mmolm3 > epsi) then
          rad_sdet = min(1e-3, max(1e-6, (rad_phy*phy_mmolm3 + rad_zoo*zoo_mmolm3) / (phy_mmolm3 + zoo_mmolm3)))
        else
          rad_sdet = min(1e-3, max(1e-6, 0.5 * (rad_phy + rad_zoo)))
        endif
        if (dia_mmolm3 + mes_mmolm3 > epsi) then
          rad_ldet = min(1e-3, max(1e-6, (rad_dia*dia_mmolm3 + rad_mes*mes_mmolm3) / (dia_mmolm3 + mes_mmolm3)))
        else
          rad_ldet = min(1e-3, max(1e-6, 0.5 * (rad_dia + rad_mes)))
        endif

        ! Run through the water column and compute the sinking rates via Rubey's equation:
        !  - radius,
        !  - density and
        !  - dynamic viscosity of seawater (that depends on T, S, and P)
        wsink1(:) = 0.0
        wsink2(:) = 0.0
        do k = 1,nk

          if (do_viscous_sinking) then
            temp2 = Temp(i,j,k) * Temp(i,j,k)
            temp3 = temp2 * Temp(i,j,k)
            temp4 = temp3 * Temp(i,j,k)
            temp5 = temp4 * Temp(i,j,k)
            ! 2. Compute seawater dynamic viscosity
            !  2a. Compute seawater dynamic viscosity at atmospheric pressure
            !    - Equations 22 and 23 from Sharqawy et al., 2010, https://doi.org/10.5004/dwt.2010.1079
            !      based off Isdale et al. (1972) "Physical properties of sea water solutions: viscosity, Desalination"
            !      that relate temperature (0 - 180 degC) and salinity (0 - 150 psu) to dynamic viscosity (kg/m/s)
            at = 1.541 + 1.998e-2*Temp(i,j,k) - 9.52e-5*temp2
            bt = 7.974 - 7.561e-2*Temp(i,j,k) + 4.724e-4*temp2
            mu_w = 4.2844e-5 + 1.0/(0.157 * (Temp(i,j,k) + 64.993)**2 - 91.296)
            wombat%dynvis_sw(i,j,k) = mu_w * (1 + at*Salt(i,j,k)*1e-3 + bt*(Salt(i,j,k)*1e-3)**2)
            !  2b. Find the effect of T and density changes on pure water using the IAPWS formulation
            ! ---- rho_w(z) : density estimate for pure water with pressure changes
            !        - rho0(T) comes from UNESCO/EOS-80 Eq. 14
            !        - rho_w(T,z) is the linearized form of the first-order compressibility correction
            ! ---- mu0(T^) : dilute-gas term (Eq. 11 in IAPWS 2008)
            ! ---- mu1(T^,rho^) : contribution to viscosity by finite density (Eq. 12 in IAPWS 2008)
            P_MPa  = (101325.0 + 9.81 * 1025.0 * wombat%zm(i,j,k)) * 1e-6  ! [MPa]
            rho0 = 999.842594 + 6.793952e-2*Temp(i,j,k) - 9.095290e-3*temp2 &
                   + 1.001685e-4*temp3 - 1.120083e-6*temp4 + 6.536332e-9*temp5  ! kg/m^3
            rho_w = rho0 / (1.0 - (P_MPa - 0.101325) / 2.2e3)
            T_hat   = (Temp(i,j,k)+273.15) / T_star
            rho_hat = rho_w / rho_star
            invT_hat = 1.0 / T_hat
            invT_hat2 = invT_hat * invT_hat
            invT_hat3 = invT_hat2 * invT_hat
            ! Find the dynamic viscosity at the dilute gas-limit (mu0)
            mu0 = 100*sqrt(T_hat) / ( mu0_H(1) + mu0_H(2)*invT_hat + mu0_H(3)*invT_hat2 + mu0_H(4)*invT_hat3 )
            ! Find the dynamic viscosity at the dilute gas-limit (mu1)
            at = 1.0 / T_hat - 1.0
            bt = rho_hat - 1.0
            poly = 0.0
            do iter = 1,21
              zval = mu1_H(iter) * at**mu1_i(iter) * bt**mu1_j(iter)
              poly = poly + zval
            enddo
            mu1 = exp(rho_hat * poly)
            mu_w_iapws = (mu0 * mu1) * mu_star  ! [Pa·s = kg/m/s]
            !  2c. Salinity corrected viscosity * pressure factor (i.e., mu_w_iapws / mu_w)
            !      NOTE: the dynamic viscosity of water actually decreases with pressure at
            !            low temperatures (Percy W. Bridgman 1925 "The Viscosity of Liquids under Pressure")
            wombat%dynvis_sw(i,j,k) = wombat%dynvis_sw(i,j,k) * mu_w_iapws / mu_w
          endif

          ! 3. Apply mineral ballasting to size classes (CaCO3 --> small, BSi --> large) and determine mass
          sdet_mmolm3 = max(epsi, wombat%f_sdet(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
          ldet_mmolm3 = max(epsi, wombat%f_ldet(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
          caco3_mmolm3 = max(epsi, wombat%f_caco3(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
          ldetsi_mmolm3 = max(epsi, wombat%f_ldetsi(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
          mass_sdet = sdet_mmolm3 * 1e-3 * 12.0 * 1e-3 * 2.5    ! [kg/m3] (assuming C is 40% of total mass)
          mass_ldet = ldet_mmolm3 * 1e-3 * 12.0 * 1e-3 * 2.5  ! [kg/m3] (assuming C is 40% of total mass)
          mass_caco3 = caco3_mmolm3 * 1e-3 * 100.0 * 1e-3     ! [kg/m3] (assuming CaCO3 is 100% of total mass)
          mass_bsi = ldetsi_mmolm3 * 1e-3 * 60.0 * 1e-3       ! [kg/m3] (assuming SiO2 is 100% of total mass)
          mass_small = mass_sdet + mass_caco3
          mass_large = mass_ldet + mass_bsi
          if (mass_small > epsi) then
            w1 = mass_sdet / mass_small
            w2 = mass_caco3 / mass_small
            rho_small = 1.0 / (max( w1/wombat%sdetrho + w2/wombat%caco3rho, epsi)) ! [kg/m3]
          else
            rho_small = wombat%sdetrho
          endif
          if (mass_large > epsi) then
            w1 = mass_ldet / mass_large
            w2 = mass_bsi / mass_large
            rho_large = 1.0 / (max( w1/wombat%sdetrho + w2/wombat%bsirho, epsi)) ! [kg/m3]
          else
            rho_large = wombat%sdetrho
          endif

          ! 4. Compute excess density given the mass of particles and their porosity
          wombat%sdet_density(i,j,k) = (1.0 - wombat%sdetphi) * rho_small + wombat%sdetphi * 1025.0 ! [kg/m3]
          wombat%ldet_density(i,j,k) = (1.0 - wombat%ldetphi) * rho_large + wombat%ldetphi * 1025.0 ! [kg/m3]

          ! 4. Compute the Rubey equation for sinking rates of particles (see Rubey, 1933)
          wsink1(k) = ( sqrt( 4.0/3.0 * 9.8 * 1025.0 * (wombat%sdet_density(i,j,k) - 1025.0) &
                        * rad_sdet**3.0 + 9.0*wombat%dynvis_sw(i,j,k)**2.0 ) &
                        - 3.0*wombat%dynvis_sw(i,j,k) ) / (1025.0 * rad_sdet) ! [m/s]
          wsink2(k) = ( sqrt( 4.0/3.0 * 9.8 * 1025.0 * (wombat%ldet_density(i,j,k) - 1025.0) &
                        * rad_ldet**3.0 + 9.0*wombat%dynvis_sw(i,j,k)**2.0 ) &
                        - 3.0*wombat%dynvis_sw(i,j,k) ) / (1025.0 * rad_ldet) ! [m/s]

        enddo
        wombat%p_wsdet(i,j,:) = wsink1(:)
        wombat%p_wsdetfe(i,j,:) = wsink1(:)
        wombat%p_wldet(i,j,:) = wsink2(:)
        wombat%p_wldetfe(i,j,:) = wsink2(:)
        wombat%p_wldetsi(i,j,:) = wsink2(:)
        wombat%p_wcaco3(i,j,:) = wsink1(:)
        wombat%p_wsafe(i,j,:) = wombat%wsafe
        wombat%p_wlafe(i,j,:) = wombat%wlafe
        wombat%sdet_radius(i,j) = rad_sdet
        wombat%ldet_radius(i,j) = rad_ldet
      else
        wombat%p_wsdet(i,j,:) = 0.0
        wombat%p_wsdetfe(i,j,:) = 0.0
        wombat%p_wldet(i,j,:) = 0.0
        wombat%p_wldetfe(i,j,:) = 0.0
        wombat%p_wldetsi(i,j,:) = 0.0
        wombat%p_wcaco3(i,j,:) = 0.0
        wombat%p_wsafe(i,j,:) = 0.0
        wombat%p_wlafe(i,j,:) = 0.0
      endif
    enddo; enddo


    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !  [Step 22] Sedimentary processes                                      !
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------

    call g_tracer_get_pointer(tracer_list, 'det_sediment', 'field', wombat%p_det_sediment) ! [mol/m2]
    call g_tracer_get_pointer(tracer_list, 'detfe_sediment', 'field', wombat%p_detfe_sediment) ! [mol/m2]
    call g_tracer_get_pointer(tracer_list, 'detsi_sediment', 'field', wombat%p_detsi_sediment) ! [mol/m2]
    call g_tracer_get_pointer(tracer_list, 'caco3_sediment', 'field', wombat%p_caco3_sediment) ! [mol/m2]

    ! Get bottom conditions, including those that influence bottom fluxes. Bottom conditions are
    ! calculated over a layer defined by wombat%bottom_thickness (default 0.1 m). This is done because
    ! the bottom layers in MOM6 are usually "vanished" layers. This approach is based on what is done
    ! in COBALT v3.
    do j = jsc,jec; do i = isc,iec;
      if (grid_kmt(i,j)>0) then
        k_bot = 0
        dzt_bot = 0.0
        do k = grid_kmt(i,j),1,-1
            if (dzt_bot < wombat%bottom_thickness) then
            k_bot = k
            dzt_bot = dzt_bot + dzt(i,j,k) ! [m]
            wombat%sedtemp(i,j) = wombat%sedtemp(i,j) + Temp(i,j,k) * dzt(i,j,k) ! [m*degC]
            wombat%sedsalt(i,j) = wombat%sedsalt(i,j) + Salt(i,j,k) * dzt(i,j,k) ! [m*psu]
            wombat%sedno3(i,j) = wombat%sedno3(i,j) + wombat%f_no3(i,j,k) * dzt(i,j,k) ! [m*mol/kg]
            wombat%sednh4(i,j) = wombat%sednh4(i,j) + wombat%f_nh4(i,j,k) * dzt(i,j,k) ! [m*mol/kg]
            wombat%sedsil(i,j) = wombat%sedsil(i,j) + wombat%f_sil(i,j,k) * dzt(i,j,k) ! [m*mol/kg]
            wombat%sedo2(i,j) = wombat%sedo2(i,j) + wombat%f_o2(i,j,k) * dzt(i,j,k) ! [m*mol/kg]
            wombat%seddic(i,j) = wombat%seddic(i,j) + wombat%f_dic(i,j,k) * dzt(i,j,k) ! [m*mol/kg]
            wombat%sedalk(i,j) = wombat%sedalk(i,j) + wombat%f_alk(i,j,k) * dzt(i,j,k) ! [m*mol/kg]
            wombat%sedhtotal(i,j) = wombat%sedhtotal(i,j) + wombat%htotal(i,j,k) * dzt(i,j,k) ! [m*mol/kg]
            endif
        enddo
        ! Subtract off overshoot
        dzt_bot_os = dzt_bot - wombat%bottom_thickness
        wombat%sedtemp(i,j) = wombat%sedtemp(i,j) - Temp(i,j,k_bot) * dzt_bot_os ! [m*degC]
        wombat%sedsalt(i,j) = wombat%sedsalt(i,j) - Salt(i,j,k_bot) * dzt_bot_os ! [m*psu]
        wombat%sedno3(i,j) = wombat%sedno3(i,j) - wombat%f_no3(i,j,k_bot) * dzt_bot_os ! [m*mol/kg]
        wombat%sednh4(i,j) = wombat%sednh4(i,j) - wombat%f_nh4(i,j,k_bot) * dzt_bot_os ! [m*mol/kg]
        wombat%sedsil(i,j) = wombat%sedsil(i,j) - wombat%f_sil(i,j,k_bot) * dzt_bot_os ! [m*mol/kg]
        wombat%sedo2(i,j) = wombat%sedo2(i,j) - wombat%f_o2(i,j,k_bot) * dzt_bot_os ! [m*mol/kg]
        wombat%seddic(i,j) = wombat%seddic(i,j) - wombat%f_dic(i,j,k_bot) * dzt_bot_os ! [m*mol/kg]
        wombat%sedalk(i,j) = wombat%sedalk(i,j) - wombat%f_alk(i,j,k_bot) * dzt_bot_os ! [m*mol/kg]
        wombat%sedhtotal(i,j) = wombat%sedhtotal(i,j) - wombat%htotal(i,j,k_bot) * dzt_bot_os ! [m*mol/kg]
        ! Convert to mol/kg
        wombat%sedtemp(i,j) = wombat%sedtemp(i,j) / wombat%bottom_thickness ! [degC]
        wombat%sedsalt(i,j) = wombat%sedsalt(i,j) / wombat%bottom_thickness ! [psu]
        wombat%sedno3(i,j) = wombat%sedno3(i,j) / wombat%bottom_thickness ! [mol/kg]
        wombat%sednh4(i,j) = wombat%sednh4(i,j) / wombat%bottom_thickness ! [mol/kg]
        wombat%sedsil(i,j) = wombat%sedsil(i,j) / wombat%bottom_thickness ! [mol/kg]
        wombat%sedo2(i,j) = wombat%sedo2(i,j) / wombat%bottom_thickness ! [mol/kg]
        wombat%seddic(i,j) = wombat%seddic(i,j) / wombat%bottom_thickness ! [mol/kg]
        wombat%sedalk(i,j) = wombat%sedalk(i,j) / wombat%bottom_thickness ! [mol/kg]
        wombat%sedhtotal(i,j) = wombat%sedhtotal(i,j) / wombat%bottom_thickness ! [mol/kg]

        ! Set seddep as full depth minus half the bottom thickness and sedmask from bottom layer
        k = grid_kmt(i,j)
        wombat%seddep(i,j) = max(0.0, wombat%zw(i,j,k) - (wombat%bottom_thickness / 2.0))
        wombat%sedmask(i,j) = grid_tmask(i,j,k)

        ! pjb: Sum the water column concentration of DIC and the organic carbon content of the
        ! sediment to approximate the interstitial (i.e., porewater) DIC concentration.
        ! We assume that the organic carbon content of the sediment (p_det_sediment) in mol/m2 is
        ! relevant over 10 cm, and therefore can be converted to mol/m3 via division by 0.1 and then
        ! subsequently converted through the mol/kg using Rho_0. With this assumption these arrays
        ! can be added together.
        ! We add these arrays together to simulate the reducing conditions of organic-rich sediments,
        ! and to calculate a lower omega for calcite, which ensures greater rates of dissolution of
        ! CaCO3 within the sediment as organic matter accumulates.
        wombat%seddic(i,j) = wombat%seddic(i,j) + wombat%p_det_sediment(i,j,1) / wombat%bottom_thickness / wombat%Rho_0
      endif
    enddo; enddo

    call FMS_ocmip2_co2calc(CO2_dope_vec, wombat%sedmask(:,:), &
        wombat%sedtemp(:,:), max(1.0, wombat%sedsalt(:,:)), &
        wombat%seddic(:,:), &
        max(wombat%sedno3(:,:) / 16., 1e-9), &
        wombat%sedsil(:,:), &
        wombat%sedalk(:,:), &
        wombat%sedhtotal(:,:)*wombat%htotal_scale_lo, &
        wombat%sedhtotal(:,:)*wombat%htotal_scale_hi, &
        wombat%sedhtotal(:,:), &
        co2_calc=trim(co2_calc), zt=wombat%seddep(:,:), &
        co3_ion=wombat%sedco3(:,:), &
        omega_calc=wombat%sedomega_cal(:,:))

    do j = jsc,jec; do i = isc,iec;
      k = grid_kmt(i,j)
      if (k > 0) then
        fbc = wombat%bbioh ** (wombat%sedtemp(i,j))

        !!!~~~ Organic carbon and iron ~~~!!!
        wombat%det_sed_remin(i,j) = wombat%detlrem_sed * fbc * wombat%p_det_sediment(i,j,1) ! [mol/m2/s]
        wombat%detfe_sed_remin(i,j) = wombat%detlrem_sed * fbc * wombat%p_detfe_sediment(i,j,1) ! [mol/m2/s]

        !!!~~~ Biogenic silica ~~~!!!
        zval = max(273.15, wombat%sedtemp(i,j) + 273.15)  ! temperature in Kelvin
        lbac_mmolm3 = max(epsi, wombat%f_lbac(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3] ! pjb - ideally this would be particle associated
        disssi_temp = exp(wombat%bsi_alpha + 0.0833*Temp(i,j,k)) / 3600.0 ! [1/s]
        disssi_usat = 1 - min(1.0, wombat%f_sil(i,j,k) / max(wombat%sileqc(i,j,k), 1e-3))
        disssi_bact = 1.0 + wombat%bsi_fbac * (lbac_mmolm3 / ( lbac_mmolm3 + wombat%bsi_kbac ))
        wombat%detsi_sed_remin(i,j) = wombat%p_detsi_sediment(i,j,1) * disssi_temp * disssi_usat * disssi_bact! [mol/m2/s]

        !!!~~~ CaCO3 dissolution ~~~!!!
        if (do_caco3_dynamics) then
            wombat%caco3_sed_remin(i,j) = wombat%caco3lrem_sed * fbc * wombat%p_caco3_sediment(i,j,1) &
                                        * (1.0 - min(1.0, max(0.01, wombat%sedomega_cal(i,j))))**(4.5)
        else
            wombat%caco3_sed_remin(i,j) = wombat%caco3lrem_sed * fbc * wombat%p_caco3_sediment(i,j,1) &
                                        * (1.0 - 0.20)**(4.5)
        endif

        !!!~~~ Benthic denitrification ~~~!!!
        if (do_benthic_denitrification) then
          ! sedimentary denitrification (Bohlen et al., 2012 Global Biogeochemical Cycles)
          !   Stoichiometry of 94 mol NO3 used per 122 mol organic carbon oxidised (Paulmier et al, 2009 BG)
          !   Hard limit where denitrification is at maximum 90% responsible for organic matter remin
          wombat%det_sed_denit(i,j) = wombat%det_sed_remin(i,j) * min(0.9 * 94.0/122.0, &
                                      (0.083 + 0.21 * 0.98**((wombat%sedo2(i,j) - wombat%sedno3(i,j))/mmol_m3_to_mol_kg)))
          wombat%fdenit(i,j) = wombat%det_sed_denit(i,j) * 122.0/94.0 / (wombat%det_sed_remin(i,j) + epsi)
        else
          wombat%det_sed_denit(i,j) = 0.0 ! [mol/m2/s]
          wombat%fdenit(i,j) = 0.0
        endif

        ! Remineralisation of sediments to supply nutrient fields.
        ! btf values are positive from the water column into the sediment.
        wombat%b_ldoc(i,j) = -wombat%det_sed_remin(i,j) * 0.5 ! [mol/m2/s]
        wombat%b_sdoc(i,j) = -wombat%det_sed_remin(i,j) * 0.5 ! [mol/m2/s]
        wombat%b_ldon(i,j) = -16./122. * wombat%det_sed_remin(i,j) ! [mol/m2/s]
        wombat%b_no3(i,j) = wombat%det_sed_denit(i,j) ! [molN/m2/s]
        wombat%b_o2(i,j) = -132./122. * (wombat%b_ldoc(i,j) + wombat%b_sdoc(i,j)) * (1.0 - wombat%fdenit(i,j)) ! [mol/m2/s]
        wombat%b_dic(i,j) = -wombat%caco3_sed_remin(i,j) ! [mol/m2/s]
        if (do_tracer_dicr) wombat%b_dicr(i,j) = wombat%b_dic(i,j) ! [mol/m2/s]
        wombat%b_fe(i,j) = -1.0 * wombat%detfe_sed_remin(i,j) ! [mol/m2/s]
        wombat%b_sil(i,j) = -1.0 * wombat%detsi_sed_remin(i,j) ! [mol/m2/s]
        wombat%b_alk(i,j) = -2.0 * wombat%caco3_sed_remin(i,j) - wombat%b_no3(i,j) ! [mol/m2/s]
      endif
    enddo; enddo


    ! Apply remineralisation rates to sediment tracers
    !-----------------------------------------------------------------------
    do j = jsc,jec; do i = isc,iec;

      if (grid_kmt(i,j) > 0) then
        wombat%p_det_sediment(i,j,1) = wombat%p_det_sediment(i,j,1) - dt * wombat%det_sed_remin(i,j) ! [mol/m2]
        wombat%p_detfe_sediment(i,j,1) = wombat%p_detfe_sediment(i,j,1) - dt * wombat%detfe_sed_remin(i,j) ! [mol/m2]
        wombat%p_detsi_sediment(i,j,1) = wombat%p_detsi_sediment(i,j,1) - dt * wombat%detsi_sed_remin(i,j) ! [mol/m2]
        wombat%p_caco3_sediment(i,j,1) = wombat%p_caco3_sediment(i,j,1) - dt * wombat%caco3_sed_remin(i,j) ! [mol/m2]
      endif
    enddo; enddo

    call g_tracer_set_values(tracer_list, 'ldoc', 'btf', wombat%b_ldoc, isd, jsd)
    call g_tracer_set_values(tracer_list, 'sdoc', 'btf', wombat%b_sdoc, isd, jsd)
    call g_tracer_set_values(tracer_list, 'ldon', 'btf', wombat%b_ldon, isd, jsd)
    call g_tracer_set_values(tracer_list, 'no3', 'btf', wombat%b_no3, isd, jsd)
    call g_tracer_set_values(tracer_list, 'o2', 'btf', wombat%b_o2, isd, jsd)
    call g_tracer_set_values(tracer_list, 'dic', 'btf', wombat%b_dic, isd, jsd)
    if (do_tracer_dicr) call g_tracer_set_values(tracer_list, 'dicr', 'btf', wombat%b_dicr, isd, jsd)
    call g_tracer_set_values(tracer_list, 'fe', 'btf', wombat%b_fe, isd, jsd)
    call g_tracer_set_values(tracer_list, 'sil', 'btf', wombat%b_sil, isd, jsd)
    call g_tracer_set_values(tracer_list, 'alk', 'btf', wombat%b_alk, isd, jsd)


    !=======================================================================
    ! Send diagnostics
    !=======================================================================

    if (wombat%id_pco2 > 0) &
      used = g_send_data(wombat%id_pco2, wombat%pco2_csurf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_htotal > 0) &
      used = g_send_data(wombat%id_htotal, wombat%htotal, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_omega_ara > 0) &
      used = g_send_data(wombat%id_omega_ara, wombat%omega_ara, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_omega_cal > 0) &
      used = g_send_data(wombat%id_omega_cal, wombat%omega_cal, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_co3 > 0) &
      used = g_send_data(wombat%id_co3, wombat%co3, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_co2_star > 0) &
      used = g_send_data(wombat%id_co2_star, wombat%co2_star, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_no3_vstf > 0) &
      used = g_send_data(wombat%id_no3_vstf, wombat%no3_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_nh4_vstf > 0) &
      used = g_send_data(wombat%id_nh4_vstf, wombat%nh4_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_dic_vstf > 0) &
      used = g_send_data(wombat%id_dic_vstf, wombat%dic_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_dicp_vstf > 0) &
      used = g_send_data(wombat%id_dicp_vstf, wombat%dic_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_alk_vstf > 0) &
      used = g_send_data(wombat%id_alk_vstf, wombat%alk_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_npp3d > 0) &
      used = g_send_data(wombat%id_npp3d, wombat%npp3d, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_rpp3d > 0) &
      used = g_send_data(wombat%id_rpp3d, wombat%rpp3d, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zsp3d > 0) &
      used = g_send_data(wombat%id_zsp3d, wombat%zsp3d, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_npp2d > 0) then
      wombat%npp2d = 0.0
      do k = 1,nk
        wombat%npp2d(isc:iec,jsc:jec) = wombat%npp2d(isc:iec,jsc:jec) + &
            wombat%npp3d(isc:iec,jsc:jec,k) * rho_dzt(isc:iec,jsc:jec,k) ! [mol/m2/s]
      enddo
      used = g_send_data(wombat%id_npp2d, wombat%npp2d, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif

    if (wombat%id_rpp2d > 0) then
      wombat%rpp2d = 0.0
      do k = 1,nk
        wombat%rpp2d(isc:iec,jsc:jec) = wombat%rpp2d(isc:iec,jsc:jec) + &
            wombat%rpp3d(isc:iec,jsc:jec,k) * rho_dzt(isc:iec,jsc:jec,k) ! [mol/m2/s]
      enddo
      used = g_send_data(wombat%id_rpp2d, wombat%rpp2d, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif

    if (wombat%id_zsp2d > 0) then
      wombat%zsp2d = 0.0
      do k = 1,nk
        wombat%zsp2d(isc:iec,jsc:jec) = wombat%zsp2d(isc:iec,jsc:jec) + &
            wombat%zsp3d(isc:iec,jsc:jec,k) * rho_dzt(isc:iec,jsc:jec,k) ! [mol/m2/s]
      enddo
      used = g_send_data(wombat%id_zsp2d, wombat%zsp2d, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif

    if (wombat%id_dynvis_sw > 0) &
      used = g_send_data(wombat%id_dynvis_sw, wombat%dynvis_sw, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_radbio > 0) &
      used = g_send_data(wombat%id_radbio, wombat%radbio, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_radmid > 0) &
      used = g_send_data(wombat%id_radmid, wombat%radmid, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_radmld > 0) &
      used = g_send_data(wombat%id_radmld, wombat%radmld, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_mumax > 0) &
      used = g_send_data(wombat%id_phy_mumax, wombat%phy_mumax, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_mu > 0) &
      used = g_send_data(wombat%id_phy_mu, wombat%phy_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_pchl_mu > 0) &
      used = g_send_data(wombat%id_pchl_mu, wombat%pchl_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_kni > 0) &
      used = g_send_data(wombat%id_phy_kni, wombat%phy_kni, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_kfe > 0) &
      used = g_send_data(wombat%id_phy_kfe, wombat%phy_kfe, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_lpar > 0) &
      used = g_send_data(wombat%id_phy_lpar, wombat%phy_lpar, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_lnit > 0) &
      used = g_send_data(wombat%id_phy_lnit, wombat%phy_lnit, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_lnh4 > 0) &
      used = g_send_data(wombat%id_phy_lnh4, wombat%phy_lnh4, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_lno3 > 0) &
      used = g_send_data(wombat%id_phy_lno3, wombat%phy_lno3, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_lfer > 0) &
      used = g_send_data(wombat%id_phy_lfer, wombat%phy_lfer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_dfeupt > 0) &
      used = g_send_data(wombat%id_phy_dfeupt, wombat%phy_dfeupt, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_mumax > 0) &
      used = g_send_data(wombat%id_dia_mumax, wombat%dia_mumax, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_mu > 0) &
      used = g_send_data(wombat%id_dia_mu, wombat%dia_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dchl_mu > 0) &
      used = g_send_data(wombat%id_dchl_mu, wombat%dchl_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_kni > 0) &
      used = g_send_data(wombat%id_dia_kni, wombat%dia_kni, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_kfe > 0) &
      used = g_send_data(wombat%id_dia_kfe, wombat%dia_kfe, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_ksi > 0) &
      used = g_send_data(wombat%id_dia_ksi, wombat%dia_ksi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_lpar > 0) &
      used = g_send_data(wombat%id_dia_lpar, wombat%dia_lpar, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_lnit > 0) &
      used = g_send_data(wombat%id_dia_lnit, wombat%dia_lnit, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_lnh4 > 0) &
      used = g_send_data(wombat%id_dia_lnh4, wombat%dia_lnh4, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_lno3 > 0) &
      used = g_send_data(wombat%id_dia_lno3, wombat%dia_lno3, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_lfer > 0) &
      used = g_send_data(wombat%id_dia_lfer, wombat%dia_lfer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_lsil > 0) &
      used = g_send_data(wombat%id_dia_lsil, wombat%dia_lsil, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_dfeupt > 0) &
      used = g_send_data(wombat%id_dia_dfeupt, wombat%dia_dfeupt, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_silupt > 0) &
      used = g_send_data(wombat%id_dia_silupt, wombat%dia_silupt, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_tri_lfer > 0) &
      used = g_send_data(wombat%id_tri_lfer, wombat%tri_lfer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_tri_lpar > 0) &
      used = g_send_data(wombat%id_tri_lpar, wombat%tri_lpar, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_trimumax > 0) &
      used = g_send_data(wombat%id_trimumax, wombat%trimumax, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_feIII > 0) &
      used = g_send_data(wombat%id_feIII, wombat%feIII, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sileqc > 0) &
      used = g_send_data(wombat%id_sileqc, wombat%sileqc, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_disssi > 0) &
      used = g_send_data(wombat%id_disssi, wombat%disssi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bsidiss > 0) &
      used = g_send_data(wombat%id_bsidiss, wombat%bsidiss, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_felig > 0) &
      used = g_send_data(wombat%id_felig, wombat%felig, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_ligK > 0) &
      used = g_send_data(wombat%id_ligK, wombat%ligK, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fecol > 0) &
      used = g_send_data(wombat%id_fecol, wombat%fecol, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fescasafe > 0) &
      used = g_send_data(wombat%id_fescasafe, wombat%fescasafe, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fescalafe > 0) &
      used = g_send_data(wombat%id_fescalafe, wombat%fescalafe, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fecoag2safe > 0) &
      used = g_send_data(wombat%id_fecoag2safe, wombat%fecoag2safe, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fecoag2lafe > 0) &
      used = g_send_data(wombat%id_fecoag2lafe, wombat%fecoag2lafe, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_safediss > 0) &
      used = g_send_data(wombat%id_safediss, wombat%safediss, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_lafediss > 0) &
      used = g_send_data(wombat%id_lafediss, wombat%lafediss, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fesources > 0) &
      used = g_send_data(wombat%id_fesources, wombat%fesources, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fesinks > 0) &
      used = g_send_data(wombat%id_fesinks, wombat%fesinks, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_feupreg > 0) &
      used = g_send_data(wombat%id_phy_feupreg, wombat%phy_feupreg, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_fedoreg > 0) &
      used = g_send_data(wombat%id_phy_fedoreg, wombat%phy_fedoreg, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phygrow > 0) &
      used = g_send_data(wombat%id_phygrow, wombat%phygrow, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phydoc > 0) &
      used = g_send_data(wombat%id_phydoc, wombat%phydoc, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phymorl > 0) &
      used = g_send_data(wombat%id_phymorl, wombat%phymorl, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phymorq > 0) &
      used = g_send_data(wombat%id_phymorq, wombat%phymorq, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_feupreg > 0) &
      used = g_send_data(wombat%id_dia_feupreg, wombat%dia_feupreg, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_fedoreg > 0) &
      used = g_send_data(wombat%id_dia_fedoreg, wombat%dia_fedoreg, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_sidoreg > 0) &
      used = g_send_data(wombat%id_dia_sidoreg, wombat%dia_sidoreg, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_diagrow > 0) &
      used = g_send_data(wombat%id_diagrow, wombat%diagrow, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_diadoc > 0) &
      used = g_send_data(wombat%id_diadoc, wombat%diadoc, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_diamorl > 0) &
      used = g_send_data(wombat%id_diamorl, wombat%diamorl, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_diamorq > 0) &
      used = g_send_data(wombat%id_diamorq, wombat%diamorq, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooeps > 0) &
      used = g_send_data(wombat%id_zooeps, wombat%zooeps, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoopreflbac > 0) &
      used = g_send_data(wombat%id_zoopreflbac, wombat%zoopreflbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooprefsbac > 0) &
      used = g_send_data(wombat%id_zooprefsbac, wombat%zooprefsbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooprefaoa > 0) &
      used = g_send_data(wombat%id_zooprefaoa, wombat%zooprefaoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooprefphy > 0) &
      used = g_send_data(wombat%id_zooprefphy, wombat%zooprefphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooprefdia > 0) &
      used = g_send_data(wombat%id_zooprefdia, wombat%zooprefdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooprefsdet > 0) &
      used = g_send_data(wombat%id_zooprefsdet, wombat%zooprefsdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoograzlbac > 0) &
      used = g_send_data(wombat%id_zoograzlbac, wombat%zoograzlbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoograzsbac > 0) &
      used = g_send_data(wombat%id_zoograzsbac, wombat%zoograzsbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoograzaoa > 0) &
      used = g_send_data(wombat%id_zoograzaoa, wombat%zoograzaoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoograzphy > 0) &
      used = g_send_data(wombat%id_zoograzphy, wombat%zoograzphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoograzdia > 0) &
      used = g_send_data(wombat%id_zoograzdia, wombat%zoograzdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoograzsdet > 0) &
      used = g_send_data(wombat%id_zoograzsdet, wombat%zoograzsdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoomorl > 0) &
      used = g_send_data(wombat%id_zoomorl, wombat%zoomorl, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoomorq > 0) &
      used = g_send_data(wombat%id_zoomorq, wombat%zoomorq, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcrlbac > 0) &
      used = g_send_data(wombat%id_zooexcrlbac, wombat%zooexcrlbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcrsbac > 0) &
      used = g_send_data(wombat%id_zooexcrsbac, wombat%zooexcrsbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcraoa > 0) &
      used = g_send_data(wombat%id_zooexcraoa, wombat%zooexcraoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcrphy > 0) &
      used = g_send_data(wombat%id_zooexcrphy, wombat%zooexcrphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcrdia > 0) &
      used = g_send_data(wombat%id_zooexcrdia, wombat%zooexcrdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcrsdet > 0) &
      used = g_send_data(wombat%id_zooexcrsdet, wombat%zooexcrsdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooegeslbac > 0) &
      used = g_send_data(wombat%id_zooegeslbac, wombat%zooegeslbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooegessbac > 0) &
      used = g_send_data(wombat%id_zooegessbac, wombat%zooegessbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooegesaoa > 0) &
      used = g_send_data(wombat%id_zooegesaoa, wombat%zooegesaoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooegesphy > 0) &
      used = g_send_data(wombat%id_zooegesphy, wombat%zooegesphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooegesdia > 0) &
      used = g_send_data(wombat%id_zooegesdia, wombat%zooegesdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooegessdet > 0) &
      used = g_send_data(wombat%id_zooegessdet, wombat%zooegessdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_meseps > 0) &
      used = g_send_data(wombat%id_meseps, wombat%meseps, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mespreflbac > 0) &
      used = g_send_data(wombat%id_mespreflbac, wombat%mespreflbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefsbac > 0) &
      used = g_send_data(wombat%id_mesprefsbac, wombat%mesprefsbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefaoa > 0) &
      used = g_send_data(wombat%id_mesprefaoa, wombat%mesprefaoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefphy > 0) &
      used = g_send_data(wombat%id_mesprefphy, wombat%mesprefphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefdia > 0) &
      used = g_send_data(wombat%id_mesprefdia, wombat%mesprefdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefsdet > 0) &
      used = g_send_data(wombat%id_mesprefsdet, wombat%mesprefsdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefldet > 0) &
      used = g_send_data(wombat%id_mesprefldet, wombat%mesprefldet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefzoo > 0) &
      used = g_send_data(wombat%id_mesprefzoo, wombat%mesprefzoo, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesgrazlbac > 0) &
      used = g_send_data(wombat%id_mesgrazlbac, wombat%mesgrazlbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesgrazsbac > 0) &
      used = g_send_data(wombat%id_mesgrazsbac, wombat%mesgrazsbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesgrazaoa > 0) &
      used = g_send_data(wombat%id_mesgrazaoa, wombat%mesgrazaoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesgrazphy > 0) &
      used = g_send_data(wombat%id_mesgrazphy, wombat%mesgrazphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesgrazdia > 0) &
      used = g_send_data(wombat%id_mesgrazdia, wombat%mesgrazdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesgrazsdet > 0) &
      used = g_send_data(wombat%id_mesgrazsdet, wombat%mesgrazsdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesgrazldet > 0) &
      used = g_send_data(wombat%id_mesgrazldet, wombat%mesgrazldet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesgrazzoo > 0) &
      used = g_send_data(wombat%id_mesgrazzoo, wombat%mesgrazzoo, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesmorl > 0) &
      used = g_send_data(wombat%id_mesmorl, wombat%mesmorl, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesmorq > 0) &
      used = g_send_data(wombat%id_mesmorq, wombat%mesmorq, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrlbac > 0) &
      used = g_send_data(wombat%id_mesexcrlbac, wombat%mesexcrlbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrsbac > 0) &
      used = g_send_data(wombat%id_mesexcrsbac, wombat%mesexcrsbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcraoa > 0) &
      used = g_send_data(wombat%id_mesexcraoa, wombat%mesexcraoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrphy > 0) &
      used = g_send_data(wombat%id_mesexcrphy, wombat%mesexcrphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrdia > 0) &
      used = g_send_data(wombat%id_mesexcrdia, wombat%mesexcrdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrsdet > 0) &
      used = g_send_data(wombat%id_mesexcrsdet, wombat%mesexcrsdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrldet > 0) &
      used = g_send_data(wombat%id_mesexcrldet, wombat%mesexcrldet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrzoo > 0) &
      used = g_send_data(wombat%id_mesexcrzoo, wombat%mesexcrzoo, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegeslbac > 0) &
      used = g_send_data(wombat%id_mesegeslbac, wombat%mesegeslbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegessbac > 0) &
      used = g_send_data(wombat%id_mesegessbac, wombat%mesegessbac, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegesaoa > 0) &
      used = g_send_data(wombat%id_mesegesaoa, wombat%mesegesaoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegesphy > 0) &
      used = g_send_data(wombat%id_mesegesphy, wombat%mesegesphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegesdia > 0) &
      used = g_send_data(wombat%id_mesegesdia, wombat%mesegesdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegessdet > 0) &
      used = g_send_data(wombat%id_mesegessdet, wombat%mesegessdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegesldet > 0) &
      used = g_send_data(wombat%id_mesegesldet, wombat%mesegesldet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegeszoo > 0) &
      used = g_send_data(wombat%id_mesegeszoo, wombat%mesegeszoo, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_reminr > 0) &
      used = g_send_data(wombat%id_reminr, wombat%reminr, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_ldocremi > 0) &
      used = g_send_data(wombat%id_ldocremi, wombat%ldocremi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sdocremi > 0) &
      used = g_send_data(wombat%id_sdocremi, wombat%sdocremi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sdocprod > 0) &
      used = g_send_data(wombat%id_sdocprod, wombat%sdocprod, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sdetremi > 0) &
      used = g_send_data(wombat%id_sdetremi, wombat%sdetremi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_ldetremi > 0) &
      used = g_send_data(wombat%id_ldetremi, wombat%ldetremi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_pic2poc > 0) &
      used = g_send_data(wombat%id_pic2poc, wombat%pic2poc, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dissratcal > 0) &
      used = g_send_data(wombat%id_dissratcal, wombat%dissratcal, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dissratara > 0) &
      used = g_send_data(wombat%id_dissratara, wombat%dissratara, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dissratpoc > 0) &
      used = g_send_data(wombat%id_dissratpoc, wombat%dissratpoc, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoodiss > 0) &
      used = g_send_data(wombat%id_zoodiss, wombat%zoodiss, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesdiss > 0) &
      used = g_send_data(wombat%id_mesdiss, wombat%mesdiss, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_caldiss > 0) &
      used = g_send_data(wombat%id_caldiss, wombat%caldiss, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aradiss > 0) &
      used = g_send_data(wombat%id_aradiss, wombat%aradiss, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_pocdiss > 0) &
      used = g_send_data(wombat%id_pocdiss, wombat%pocdiss, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoa_loxy > 0) &
      used = g_send_data(wombat%id_aoa_loxy, wombat%aoa_loxy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoa_lnh4 > 0) &
      used = g_send_data(wombat%id_aoa_lnh4, wombat%aoa_lnh4, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoa_eno3 > 0) &
      used = g_send_data(wombat%id_aoa_eno3, wombat%aoa_eno3, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoa_mumax > 0) &
      used = g_send_data(wombat%id_aoa_mumax, wombat%aoa_mumax, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoa_mu > 0) &
      used = g_send_data(wombat%id_aoa_mu, wombat%aoa_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoagrow > 0) &
      used = g_send_data(wombat%id_aoagrow, wombat%aoagrow, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoaresp > 0) &
      used = g_send_data(wombat%id_aoaresp, wombat%aoaresp, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoamorl > 0) &
      used = g_send_data(wombat%id_aoamorl, wombat%aoamorl, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoamorq > 0) &
      used = g_send_data(wombat%id_aoamorq, wombat%aoamorq, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_lbacydoc > 0) &
      used = g_send_data(wombat%id_lbacydoc, wombat%lbacydoc, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_lbacgrow > 0) &
      used = g_send_data(wombat%id_lbacgrow, wombat%lbacgrow, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_lbacresp > 0) &
      used = g_send_data(wombat%id_lbacresp, wombat%lbacresp, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_lbacpnh4 > 0) &
      used = g_send_data(wombat%id_lbacpnh4, wombat%lbacpnh4, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_lbacpco2 > 0) &
      used = g_send_data(wombat%id_lbacpco2, wombat%lbacpco2, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_lbacufer > 0) &
      used = g_send_data(wombat%id_lbacufer, wombat%lbacufer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_lbac_mu > 0) &
      used = g_send_data(wombat%id_lbac_mu, wombat%lbac_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_lbac_fanaer > 0) &
      used = g_send_data(wombat%id_lbac_fanaer, wombat%lbac_fanaer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_lbac_ffelim > 0) &
      used = g_send_data(wombat%id_lbac_ffelim, wombat%lbac_ffelim, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_lbacmorl > 0) &
      used = g_send_data(wombat%id_lbacmorl, wombat%lbacmorl, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_lbacmorq > 0) &
      used = g_send_data(wombat%id_lbacmorq, wombat%lbacmorq, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_lbacdeni > 0) &
      used = g_send_data(wombat%id_lbacdeni, wombat%lbacdeni, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sbacydoc > 0) &
      used = g_send_data(wombat%id_sbacydoc, wombat%sbacydoc, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sbacgrow > 0) &
      used = g_send_data(wombat%id_sbacgrow, wombat%sbacgrow, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sbacresp > 0) &
      used = g_send_data(wombat%id_sbacresp, wombat%sbacresp, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sbacpnh4 > 0) &
      used = g_send_data(wombat%id_sbacpnh4, wombat%sbacpnh4, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sbacpco2 > 0) &
      used = g_send_data(wombat%id_sbacpco2, wombat%sbacpco2, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sbacufer > 0) &
      used = g_send_data(wombat%id_sbacufer, wombat%sbacufer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sbac_mu > 0) &
      used = g_send_data(wombat%id_sbac_mu, wombat%sbac_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sbac_fanaer > 0) &
      used = g_send_data(wombat%id_sbac_fanaer, wombat%sbac_fanaer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sbac_ffelim > 0) &
      used = g_send_data(wombat%id_sbac_ffelim, wombat%sbac_ffelim, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sbacmorl > 0) &
      used = g_send_data(wombat%id_sbacmorl, wombat%sbacmorl, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sbacmorq > 0) &
      used = g_send_data(wombat%id_sbacmorq, wombat%sbacmorq, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sbacdeni > 0) &
      used = g_send_data(wombat%id_sbacdeni, wombat%sbacdeni, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aox_lnh4 > 0) &
      used = g_send_data(wombat%id_aox_lnh4, wombat%aox_lnh4, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aox_mu > 0) &
      used = g_send_data(wombat%id_aox_mu, wombat%aox_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_nitrfix > 0) &
      used = g_send_data(wombat%id_nitrfix, wombat%nitrfix, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_ammox > 0) &
      used = g_send_data(wombat%id_ammox, wombat%ammox, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_anammox > 0) &
      used = g_send_data(wombat%id_anammox, wombat%anammox, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sdet_density > 0) &
      used = g_send_data(wombat%id_sdet_density, wombat%sdet_density, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_ldet_density > 0) &
      used = g_send_data(wombat%id_ldet_density, wombat%ldet_density, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_det_sed_remin > 0) &
      used = g_send_data(wombat%id_det_sed_remin, wombat%det_sed_remin, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_detfe_sed_remin > 0) &
      used = g_send_data(wombat%id_detfe_sed_remin, wombat%detfe_sed_remin, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_detsi_sed_remin > 0) &
      used = g_send_data(wombat%id_detsi_sed_remin, wombat%detsi_sed_remin, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_det_sed_denit > 0) &
      used = g_send_data(wombat%id_det_sed_denit, wombat%det_sed_denit, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_fdenit > 0) &
      used = g_send_data(wombat%id_fdenit, wombat%fdenit, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_caco3_sed_remin > 0) &
      used = g_send_data(wombat%id_caco3_sed_remin, wombat%caco3_sed_remin, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_zeuphot > 0) &
      used = g_send_data(wombat%id_zeuphot, wombat%zeuphot, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sdet_radius > 0) &
      used = g_send_data(wombat%id_sdet_radius, wombat%sdet_radius, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_ldet_radius > 0) &
      used = g_send_data(wombat%id_ldet_radius, wombat%ldet_radius, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_seddep > 0) &
      used = g_send_data(wombat%id_seddep, wombat%seddep, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedmask > 0) &
      used = g_send_data(wombat%id_sedmask, wombat%sedmask, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedtemp > 0) &
      used = g_send_data(wombat%id_sedtemp, wombat%sedtemp, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedsalt > 0) &
      used = g_send_data(wombat%id_sedsalt, wombat%sedsalt, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedno3 > 0) &
      used = g_send_data(wombat%id_sedno3, wombat%sedno3, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sednh4 > 0) &
      used = g_send_data(wombat%id_sednh4, wombat%sednh4, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedsil > 0) &
      used = g_send_data(wombat%id_sedsil, wombat%sedsil, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedo2 > 0) &
      used = g_send_data(wombat%id_sedo2, wombat%sedo2, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_seddic > 0) &
      used = g_send_data(wombat%id_seddic, wombat%seddic, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedalk > 0) &
      used = g_send_data(wombat%id_sedalk, wombat%sedalk, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedhtotal > 0) &
      used = g_send_data(wombat%id_sedhtotal, wombat%sedhtotal, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedco3 > 0) &
      used = g_send_data(wombat%id_sedco3, wombat%sedco3, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedomega_cal > 0) &
      used = g_send_data(wombat%id_sedomega_cal, wombat%sedomega_cal, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    deallocate(kmeuph, k100)

  end subroutine generic_WOMBATmid_update_from_source

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATmid_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   Calculate and set coupler values at the surface / bottom of the ocean.
  !   User must provide the calculations for these boundary values.
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATmid_set_boundary_values(tracer_list, SST, SSS, rho, ilb, jlb, tau, dzt)
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !
  !  <IN NAME="SST" TYPE="real, dimension(ilb:,jlb:)">
  !   Sea Surface Temperature
  !  </IN>
  !
  !  <IN NAME="SSS" TYPE="real, dimension(ilb:,jlb:)">
  !   Sea Surface Salinity
  !  </IN>
  !
  !  <IN NAME="rho" TYPE="real, dimension(ilb:,jlb:,:,:)">
  !   Ocean density
  !  </IN>
  !
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index of %field
  !  </IN>
  !
  !  <IN NAME="dzt" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Layer thickness
  !  </IN>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATmid_set_boundary_values(tracer_list, SST, SSS, rho, ilb, jlb, tau, dzt)
    type(g_tracer_type), pointer, intent(in)           :: tracer_list
    real, dimension(ilb:,jlb:), intent(in)             :: SST, SSS
    real, dimension(ilb:,jlb:,:,:), intent(in)         :: rho
    integer, intent(in)                                :: ilb, jlb, tau
    real, dimension(ilb:,jlb:,:), optional, intent(in) :: dzt

    integer                                 :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, i, j
    real                                    :: sal, ST, o2_solubility
    real                                    :: tt, tk, tk100, ts, ts2, ts3, ts4, ts5
    real                                    :: mmol_m3_to_mol_kg
    real, dimension(:,:,:), pointer         :: grid_tmask
    character(len=fm_string_len), parameter :: sub_name = 'generic_WOMBATmid_set_boundary_values'

    ! Get the necessary properties
    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, &
        grid_tmask=grid_tmask)

    call g_tracer_get_pointer(tracer_list, 'o2', 'field', wombat%p_o2)

    ! Some unit conversion factors
    mmol_m3_to_mol_kg = 1.e-3 / wombat%Rho_0

    ! nnz: Since the generic_WOMBATmid_update_from_source() subroutine is called by this time
    ! the following if block is not really necessary (since this calculation is already done in
    ! source).
    ! It is only neccessary if source routine is commented out for debugging.
    ! Note: In order for this to work we should NOT zero out the coupler values for generic tracers
    ! This zeroing is done for non-generic TOPAZ by calling zero_ocean_sfc.
    ! Since the coupler values here are non-cumulative there is no need to zero them out anyway.
    if (wombat%init .OR. wombat%force_update_fluxes) then
      ! Get necessary fields
      call g_tracer_get_values(tracer_list, 'dic', 'field', wombat%f_dic, isd, jsd, ntau=1)
      call g_tracer_get_values(tracer_list, 'no3', 'field', wombat%f_no3, isd, jsd, ntau=1)
      call g_tracer_get_values(tracer_list, 'sil', 'field', wombat%f_sil, isd, jsd, ntau=1)
      call g_tracer_get_values(tracer_list, 'alk', 'field', wombat%f_alk, isd, jsd, ntau=1)

      do j = jsc, jec; do i = isc, iec
          wombat%htotallo(i,j) = wombat%htotal_scale_lo * wombat%htotal(i,j,1)
          wombat%htotalhi(i,j) = wombat%htotal_scale_hi * wombat%htotal(i,j,1)
      enddo; enddo

      if ((trim(co2_calc) == 'mocsy') .and. (.not. present(dzt))) then
        call mpp_error(FATAL,"mocsy method of co2_calc needs dzt to be passed to the "// &
            "FMS_ocmip2_co2calc subroutine.")
      endif

      call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,1), &
          SST(:,:), max(1.0, SSS(:,:)), &
          max(0.0, wombat%f_dic(:,:,1)), &
          max(0.0, wombat%f_no3(:,:,1) / 16.), &
          max(0.0, wombat%f_sil(:,:,1)), &
          max(0.0, wombat%f_alk(:,:,1)), &
          wombat%htotallo(:,:), wombat%htotalhi(:,:), &
          wombat%htotal(:,:,1), &
          co2_calc=trim(co2_calc), &
          zt=dzt(:,:,1), &
          co2star=wombat%co2_csurf(:,:), alpha=wombat%co2_alpha(:,:), co3_ion=wombat%co3(:,:,1), &
          pCO2surf=wombat%pco2_csurf(:,:), omega_arag=wombat%omega_ara(:,:,1), omega_calc=wombat%omega_cal(:,:,1))

      call g_tracer_set_values(tracer_list, 'dic', 'alpha', wombat%co2_alpha, isd, jsd)
      call g_tracer_set_values(tracer_list, 'dic', 'csurf', wombat%co2_csurf, isd, jsd)

      ! nnz: If source is called uncomment the following
      wombat%init = .false. !nnz: This is necessary since the above calls appear in source subroutine too.
    endif

    call g_tracer_get_values(tracer_list, 'dic', 'alpha', wombat%co2_alpha ,isd, jsd)
    call g_tracer_get_values(tracer_list, 'dic', 'csurf', wombat%co2_csurf ,isd, jsd)

    do j=jsc,jec ; do i=isc,iec
      !-----------------------------------------------------------------------
      !  Compute the Schmidt number of O2 in seawater using Wanninkhof (2014)
      !-----------------------------------------------------------------------
      ST = SST(i,j)
      wombat%co2_sc_no(i,j) = wombat%a1_co2 + ST * (wombat%a2_co2 + ST * (wombat%a3_co2 + ST &
          * (wombat%a4_co2 + ST * wombat%a5_co2))) * grid_tmask(i,j,1)

      wombat%co2_alpha(i,j) = wombat%co2_alpha(i,j) * wombat%Rho_0 !nnz: MOM has rho(i,j,1,tau)
      wombat%co2_csurf(i,j) = wombat%co2_csurf(i,j) * wombat%Rho_0 !nnz: MOM has rho(i,j,1,tau)
    enddo; enddo

    ! Set %csurf, %alpha and %sc_no for these tracers. This will mark them
    ! for sending fluxes to coupler
    call g_tracer_set_values(tracer_list, 'dic', 'alpha', wombat%co2_alpha, isd, jsd)
    call g_tracer_set_values(tracer_list, 'dic', 'csurf', wombat%co2_csurf, isd, jsd)
    call g_tracer_set_values(tracer_list, 'dic', 'sc_no', wombat%co2_sc_no, isd, jsd)

    call g_tracer_get_values(tracer_list, 'o2', 'alpha', wombat%o2_alpha, isd, jsd)
    call g_tracer_get_values(tracer_list, 'o2', 'csurf', wombat%o2_csurf ,isd, jsd)

    do j=jsc,jec ; do i=isc,iec
      !-----------------------------------------------------------------------
      ! Compute the oxygen saturation concentration at 1 atm total
      ! pressure in mol/kg given the temperature (t, in deg C) and
      ! the salinity (s, in permil)
      !
      ! From Garcia and Gordon (1992), Limnology and Oceonography.
      ! The formula used is from page 1310, eq (8).
      !
      ! *** Note: the "a3*ts^2" term (in the paper) is incorrect. ***
      ! *** It shouldn't be there.                                ***
      ! *** PJB & DTS : Also note that Garcia & Gordon (1992) solve
      !     their polynomial for the concentration of O2 in mol/kg,
      !     which already includes the atm-1 concentration (0.20946)
      !
      ! o2_solubility is defined between T(freezing) <= T <= 40 deg C and
      !                                   0 permil <= S <= 42 permil
      ! We impose these bounds here.
      !
      ! check value: T = 10 deg C, S = 35 permil,
      !              o2_solubility = 0.282015 / 0.20946 mol m-3
      !-----------------------------------------------------------------------
      sal = SSS(i,j) ; ST = SST(i,j)

      ! jgj 2015/05/14 impose temperature and salinity bounds for o2sat
      sal = min(42.0, max(0.0, sal))
      tt = 298.15 - min(40.0, max(0.0, ST))
      tk = 273.15 + min(40.0, max(0.0, ST))
      ts = log(tt / tk)
      ts2 = ts  * ts
      ts3 = ts2 * ts
      ts4 = ts3 * ts
      ts5 = ts4 * ts

      ! Note that we divide by 0.20946 to convert mol m-3 to mol m-3 atm-1
      o2_solubility = (1000.0/22391.6) * grid_tmask(i,j,1) *  & !convert from ml/l to mol m-3
          exp( wombat%a_0 + wombat%a_1*ts + wombat%a_2*ts2 + wombat%a_3*ts3 + wombat%a_4*ts4 + &
              wombat%a_5*ts5 + (wombat%b_0 + wombat%b_1*ts + wombat%b_2*ts2 + wombat%b_3*ts3 + &
              wombat%c_0*sal)*sal) / 0.20946

      !-----------------------------------------------------------------------
      !  Compute the Schmidt number of O2 in seawater using Wanninkhof (2014)
      !-----------------------------------------------------------------------
      wombat%o2_sc_no(i,j) = wombat%a1_o2 + ST * (wombat%a2_o2 + ST * (wombat%a3_o2 + ST &
          * (wombat%a4_o2 + ST * wombat%a5_o2))) * grid_tmask(i,j,1)

      wombat%o2_alpha(i,j) = o2_solubility ! already in mol/m3 atm-1 (see above)
      wombat%o2_csurf(i,j) = wombat%p_o2(i,j,1,tau) * wombat%Rho_0 !nnz: MOM has rho(i,j,1,tau)
    enddo; enddo

    ! Set %csurf, %alpha and %sc_no for these tracers. This will mark them
    ! for sending fluxes to coupler
    call g_tracer_set_values(tracer_list, 'o2', 'alpha', wombat%o2_alpha, isd, jsd)
    call g_tracer_set_values(tracer_list, 'o2', 'csurf', wombat%o2_csurf, isd, jsd)
    call g_tracer_set_values(tracer_list, 'o2', 'sc_no', wombat%o2_sc_no, isd, jsd)

  end subroutine generic_WOMBATmid_set_boundary_values

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATmid_end">
  !  <OVERVIEW>
  !   End the module.
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   Deallocate all work arrays
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATmid_end
  !  </TEMPLATE>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATmid_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_WOMBATmid_end'

    call user_deallocate_arrays

  end subroutine generic_WOMBATmid_end

  !#######################################################################
  !
  ! This is an internal sub, not a public interface.
  ! Allocate all the work arrays to be used in this module.
  !
  subroutine user_allocate_arrays
    integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau)

    ! Used in ocmip2_co2calc
    CO2_dope_vec%isc = isc; CO2_dope_vec%iec = iec
    CO2_dope_vec%jsc = jsc; CO2_dope_vec%jec = jec
    CO2_dope_vec%isd = isd; CO2_dope_vec%ied = ied
    CO2_dope_vec%jsd = jsd; CO2_dope_vec%jed = jed

    allocate(wombat%htotallo(isd:ied, jsd:jed)); wombat%htotallo(:,:)=0.0
    allocate(wombat%htotalhi(isd:ied, jsd:jed)); wombat%htotalhi(:,:)=0.0
    allocate(wombat%htotal(isd:ied, jsd:jed, 1:nk)); wombat%htotal(:,:,:)=wombat%htotal_in
    allocate(wombat%omega_ara(isd:ied, jsd:jed, 1:nk)); wombat%omega_ara(:,:,:)=1.0
    allocate(wombat%omega_cal(isd:ied, jsd:jed, 1:nk)); wombat%omega_cal(:,:,:)=1.0
    allocate(wombat%co3(isd:ied, jsd:jed, 1:nk)); wombat%co3(:,:,:)=0.0
    allocate(wombat%co2_star(isd:ied, jsd:jed, 1:nk)); wombat%co2_star(:,:,:)=0.0
    allocate(wombat%co2_csurf(isd:ied, jsd:jed)); wombat%co2_csurf(:,:)=0.0
    allocate(wombat%co2_alpha(isd:ied, jsd:jed)); wombat%co2_alpha(:,:)=0.0
    allocate(wombat%co2_sc_no(isd:ied, jsd:jed)); wombat%co2_sc_no(:,:)=0.0
    allocate(wombat%pco2_csurf(isd:ied, jsd:jed)); wombat%pco2_csurf(:,:)=0.0
    allocate(wombat%o2_csurf(isd:ied, jsd:jed)); wombat%o2_csurf(:,:)=0.0
    allocate(wombat%o2_alpha(isd:ied, jsd:jed)); wombat%o2_alpha(:,:)=0.0
    allocate(wombat%o2_sc_no(isd:ied, jsd:jed)); wombat%o2_sc_no(:,:)=0.0
    allocate(wombat%no3_vstf(isd:ied, jsd:jed)); wombat%no3_vstf(:,:)=0.0
    allocate(wombat%nh4_vstf(isd:ied, jsd:jed)); wombat%nh4_vstf(:,:)=0.0
    allocate(wombat%dic_vstf(isd:ied, jsd:jed)); wombat%dic_vstf(:,:)=0.0
    allocate(wombat%alk_vstf(isd:ied, jsd:jed)); wombat%alk_vstf(:,:)=0.0

    allocate(wombat%f_dic(isd:ied, jsd:jed, 1:nk)); wombat%f_dic(:,:,:)=0.0
    allocate(wombat%f_dicr(isd:ied, jsd:jed, 1:nk)); wombat%f_dicr(:,:,:)=0.0
    allocate(wombat%f_alk(isd:ied, jsd:jed, 1:nk)); wombat%f_alk(:,:,:)=0.0
    allocate(wombat%f_no3(isd:ied, jsd:jed, 1:nk)); wombat%f_no3(:,:,:)=0.0
    allocate(wombat%f_nh4(isd:ied, jsd:jed, 1:nk)); wombat%f_nh4(:,:,:)=0.0
    allocate(wombat%f_sil(isd:ied, jsd:jed, 1:nk)); wombat%f_sil(:,:,:)=0.0
    allocate(wombat%f_phy(isd:ied, jsd:jed, 1:nk)); wombat%f_phy(:,:,:)=0.0
    allocate(wombat%f_pchl(isd:ied, jsd:jed, 1:nk)); wombat%f_pchl(:,:,:)=0.0
    allocate(wombat%f_phyfe(isd:ied, jsd:jed, 1:nk)); wombat%f_phyfe(:,:,:)=0.0
    allocate(wombat%f_dia(isd:ied, jsd:jed, 1:nk)); wombat%f_dia(:,:,:)=0.0
    allocate(wombat%f_dchl(isd:ied, jsd:jed, 1:nk)); wombat%f_dchl(:,:,:)=0.0
    allocate(wombat%f_diafe(isd:ied, jsd:jed, 1:nk)); wombat%f_diafe(:,:,:)=0.0
    allocate(wombat%f_diasi(isd:ied, jsd:jed, 1:nk)); wombat%f_diasi(:,:,:)=0.0
    allocate(wombat%f_zoo(isd:ied, jsd:jed, 1:nk)); wombat%f_zoo(:,:,:)=0.0
    allocate(wombat%f_zoofe(isd:ied, jsd:jed, 1:nk)); wombat%f_zoofe(:,:,:)=0.0
    allocate(wombat%f_mes(isd:ied, jsd:jed, 1:nk)); wombat%f_mes(:,:,:)=0.0
    allocate(wombat%f_mesfe(isd:ied, jsd:jed, 1:nk)); wombat%f_mesfe(:,:,:)=0.0
    allocate(wombat%f_sdet(isd:ied, jsd:jed, 1:nk)); wombat%f_sdet(:,:,:)=0.0
    allocate(wombat%f_sdetfe(isd:ied, jsd:jed, 1:nk)); wombat%f_sdetfe(:,:,:)=0.0
    allocate(wombat%f_ldet(isd:ied, jsd:jed, 1:nk)); wombat%f_ldet(:,:,:)=0.0
    allocate(wombat%f_ldetfe(isd:ied, jsd:jed, 1:nk)); wombat%f_ldetfe(:,:,:)=0.0
    allocate(wombat%f_ldetsi(isd:ied, jsd:jed, 1:nk)); wombat%f_ldetsi(:,:,:)=0.0
    allocate(wombat%f_ldoc(isd:ied, jsd:jed, 1:nk)); wombat%f_ldoc(:,:,:)=0.0
    allocate(wombat%f_sdoc(isd:ied, jsd:jed, 1:nk)); wombat%f_sdoc(:,:,:)=0.0
    allocate(wombat%f_ldon(isd:ied, jsd:jed, 1:nk)); wombat%f_ldon(:,:,:)=0.0
    allocate(wombat%f_lbac(isd:ied, jsd:jed, 1:nk)); wombat%f_lbac(:,:,:)=0.0
    allocate(wombat%f_sbac(isd:ied, jsd:jed, 1:nk)); wombat%f_sbac(:,:,:)=0.0
    allocate(wombat%f_aoa(isd:ied, jsd:jed, 1:nk)); wombat%f_aoa(:,:,:)=0.0
    allocate(wombat%f_o2(isd:ied, jsd:jed, 1:nk)); wombat%f_o2(:,:,:)=0.0
    allocate(wombat%f_caco3(isd:ied, jsd:jed, 1:nk)); wombat%f_caco3(:,:,:)=0.0
    allocate(wombat%f_fe(isd:ied, jsd:jed, 1:nk)); wombat%f_fe(:,:,:)=0.0
    allocate(wombat%f_safe(isd:ied, jsd:jed, 1:nk)); wombat%f_safe(:,:,:)=0.0
    allocate(wombat%f_lafe(isd:ied, jsd:jed, 1:nk)); wombat%f_lafe(:,:,:)=0.0

    allocate(wombat%b_ldoc(isd:ied, jsd:jed)); wombat%b_ldoc(:,:)=0.0
    allocate(wombat%b_sdoc(isd:ied, jsd:jed)); wombat%b_sdoc(:,:)=0.0
    allocate(wombat%b_ldon(isd:ied, jsd:jed)); wombat%b_ldon(:,:)=0.0
    allocate(wombat%b_no3(isd:ied, jsd:jed)); wombat%b_no3(:,:)=0.0
    allocate(wombat%b_o2(isd:ied, jsd:jed)); wombat%b_o2(:,:)=0.0
    allocate(wombat%b_dic(isd:ied, jsd:jed)); wombat%b_dic(:,:)=0.0
    if (do_tracer_dicr) then
      allocate(wombat%b_dicr(isd:ied, jsd:jed)); wombat%b_dicr(:,:)=0.0
    endif
    allocate(wombat%b_fe(isd:ied, jsd:jed)); wombat%b_fe(:,:)=0.0
    allocate(wombat%b_sil(isd:ied, jsd:jed)); wombat%b_sil(:,:)=0.0
    allocate(wombat%b_alk(isd:ied, jsd:jed)); wombat%b_alk(:,:)=0.0

    allocate(wombat%dynvis_sw(isd:ied, jsd:jed, 1:nk)); wombat%dynvis_sw(:,:,:)=0.0
    allocate(wombat%radbio(isd:ied, jsd:jed, 1:nk)); wombat%radbio(:,:,:)=0.0
    allocate(wombat%radmid(isd:ied, jsd:jed, 1:nk)); wombat%radmid(:,:,:)=0.0
    allocate(wombat%radmld(isd:ied, jsd:jed, 1:nk)); wombat%radmld(:,:,:)=0.0
    allocate(wombat%npp3d(isd:ied, jsd:jed, 1:nk)); wombat%npp3d(:,:,:)=0.0
    allocate(wombat%rpp3d(isd:ied, jsd:jed, 1:nk)); wombat%rpp3d(:,:,:)=0.0
    allocate(wombat%zsp3d(isd:ied, jsd:jed, 1:nk)); wombat%zsp3d(:,:,:)=0.0
    allocate(wombat%npp2d(isd:ied, jsd:jed)); wombat%npp2d(:,:)=0.0
    allocate(wombat%rpp2d(isd:ied, jsd:jed)); wombat%rpp2d(:,:)=0.0
    allocate(wombat%zsp2d(isd:ied, jsd:jed)); wombat%zsp2d(:,:)=0.0
    allocate(wombat%phy_mumax(isd:ied, jsd:jed, 1:nk)); wombat%phy_mumax(:,:,:)=0.0
    allocate(wombat%phy_mu(isd:ied, jsd:jed, 1:nk)); wombat%phy_mu(:,:,:)=0.0
    allocate(wombat%pchl_mu(isd:ied, jsd:jed, 1:nk)); wombat%pchl_mu(:,:,:)=0.0
    allocate(wombat%phy_kni(isd:ied, jsd:jed, 1:nk)); wombat%phy_kni(:,:,:)=0.0
    allocate(wombat%phy_kfe(isd:ied, jsd:jed, 1:nk)); wombat%phy_kfe(:,:,:)=0.0
    allocate(wombat%phy_lpar(isd:ied, jsd:jed, 1:nk)); wombat%phy_lpar(:,:,:)=0.0
    allocate(wombat%phy_lnit(isd:ied, jsd:jed, 1:nk)); wombat%phy_lnit(:,:,:)=0.0
    allocate(wombat%phy_lnh4(isd:ied, jsd:jed, 1:nk)); wombat%phy_lnh4(:,:,:)=0.0
    allocate(wombat%phy_lno3(isd:ied, jsd:jed, 1:nk)); wombat%phy_lno3(:,:,:)=0.0
    allocate(wombat%phy_lfer(isd:ied, jsd:jed, 1:nk)); wombat%phy_lfer(:,:,:)=0.0
    allocate(wombat%phy_dfeupt(isd:ied, jsd:jed, 1:nk)); wombat%phy_dfeupt(:,:,:)=0.0
    allocate(wombat%dia_mumax(isd:ied, jsd:jed, 1:nk)); wombat%dia_mumax(:,:,:)=0.0
    allocate(wombat%dia_mu(isd:ied, jsd:jed, 1:nk)); wombat%dia_mu(:,:,:)=0.0
    allocate(wombat%dchl_mu(isd:ied, jsd:jed, 1:nk)); wombat%dchl_mu(:,:,:)=0.0
    allocate(wombat%dia_kni(isd:ied, jsd:jed, 1:nk)); wombat%dia_kni(:,:,:)=0.0
    allocate(wombat%dia_kfe(isd:ied, jsd:jed, 1:nk)); wombat%dia_kfe(:,:,:)=0.0
    allocate(wombat%dia_ksi(isd:ied, jsd:jed, 1:nk)); wombat%dia_ksi(:,:,:)=0.0
    allocate(wombat%dia_lpar(isd:ied, jsd:jed, 1:nk)); wombat%dia_lpar(:,:,:)=0.0
    allocate(wombat%dia_lnit(isd:ied, jsd:jed, 1:nk)); wombat%dia_lnit(:,:,:)=0.0
    allocate(wombat%dia_lnh4(isd:ied, jsd:jed, 1:nk)); wombat%dia_lnh4(:,:,:)=0.0
    allocate(wombat%dia_lno3(isd:ied, jsd:jed, 1:nk)); wombat%dia_lno3(:,:,:)=0.0
    allocate(wombat%dia_lfer(isd:ied, jsd:jed, 1:nk)); wombat%dia_lfer(:,:,:)=0.0
    allocate(wombat%dia_lsil(isd:ied, jsd:jed, 1:nk)); wombat%dia_lsil(:,:,:)=0.0
    allocate(wombat%dia_dfeupt(isd:ied, jsd:jed, 1:nk)); wombat%dia_dfeupt(:,:,:)=0.0
    allocate(wombat%dia_silupt(isd:ied, jsd:jed, 1:nk)); wombat%dia_silupt(:,:,:)=0.0
    allocate(wombat%tri_lfer(isd:ied, jsd:jed, 1:nk)); wombat%tri_lfer(:,:,:)=0.0
    allocate(wombat%tri_lpar(isd:ied, jsd:jed, 1:nk)); wombat%tri_lpar(:,:,:)=0.0
    allocate(wombat%trimumax(isd:ied, jsd:jed, 1:nk)); wombat%trimumax(:,:,:)=0.0
    allocate(wombat%feIII(isd:ied, jsd:jed, 1:nk)); wombat%feIII(:,:,:)=0.0
    allocate(wombat%sileqc(isd:ied, jsd:jed, 1:nk)); wombat%sileqc(:,:,:)=0.0
    allocate(wombat%disssi(isd:ied, jsd:jed, 1:nk)); wombat%disssi(:,:,:)=0.0
    allocate(wombat%bsidiss(isd:ied, jsd:jed, 1:nk)); wombat%bsidiss(:,:,:)=0.0
    allocate(wombat%felig(isd:ied, jsd:jed, 1:nk)); wombat%felig(:,:,:)=0.0
    allocate(wombat%ligK(isd:ied, jsd:jed, 1:nk)); wombat%ligK(:,:,:)=0.0
    allocate(wombat%fecol(isd:ied, jsd:jed, 1:nk)); wombat%fecol(:,:,:)=0.0
    allocate(wombat%fescasafe(isd:ied, jsd:jed, 1:nk)); wombat%fescasafe(:,:,:)=0.0
    allocate(wombat%fescalafe(isd:ied, jsd:jed, 1:nk)); wombat%fescalafe(:,:,:)=0.0
    allocate(wombat%fecoag2safe(isd:ied, jsd:jed, 1:nk)); wombat%fecoag2safe(:,:,:)=0.0
    allocate(wombat%fecoag2lafe(isd:ied, jsd:jed, 1:nk)); wombat%fecoag2lafe(:,:,:)=0.0
    allocate(wombat%safediss(isd:ied, jsd:jed, 1:nk)); wombat%safediss(:,:,:)=0.0
    allocate(wombat%lafediss(isd:ied, jsd:jed, 1:nk)); wombat%lafediss(:,:,:)=0.0
    allocate(wombat%fesources(isd:ied, jsd:jed, 1:nk)); wombat%fesources(:,:,:)=0.0
    allocate(wombat%fesinks(isd:ied, jsd:jed, 1:nk)); wombat%fesinks(:,:,:)=0.0
    allocate(wombat%phy_feupreg(isd:ied, jsd:jed, 1:nk)); wombat%phy_feupreg(:,:,:)=0.0
    allocate(wombat%phy_fedoreg(isd:ied, jsd:jed, 1:nk)); wombat%phy_fedoreg(:,:,:)=0.0
    allocate(wombat%phygrow(isd:ied, jsd:jed, 1:nk)); wombat%phygrow(:,:,:)=0.0
    allocate(wombat%phydoc(isd:ied, jsd:jed, 1:nk)); wombat%phydoc(:,:,:)=0.0
    allocate(wombat%phymorl(isd:ied, jsd:jed, 1:nk)); wombat%phymorl(:,:,:)=0.0
    allocate(wombat%phymorq(isd:ied, jsd:jed, 1:nk)); wombat%phymorq(:,:,:)=0.0
    allocate(wombat%dia_feupreg(isd:ied, jsd:jed, 1:nk)); wombat%dia_feupreg(:,:,:)=0.0
    allocate(wombat%dia_fedoreg(isd:ied, jsd:jed, 1:nk)); wombat%dia_fedoreg(:,:,:)=0.0
    allocate(wombat%dia_sidoreg(isd:ied, jsd:jed, 1:nk)); wombat%dia_sidoreg(:,:,:)=0.0
    allocate(wombat%diagrow(isd:ied, jsd:jed, 1:nk)); wombat%diagrow(:,:,:)=0.0
    allocate(wombat%diadoc(isd:ied, jsd:jed, 1:nk)); wombat%diadoc(:,:,:)=0.0
    allocate(wombat%diamorl(isd:ied, jsd:jed, 1:nk)); wombat%diamorl(:,:,:)=0.0
    allocate(wombat%diamorq(isd:ied, jsd:jed, 1:nk)); wombat%diamorq(:,:,:)=0.0
    allocate(wombat%zooeps(isd:ied, jsd:jed, 1:nk)); wombat%zooeps(:,:,:)=0.0
    allocate(wombat%zoopreflbac(isd:ied, jsd:jed, 1:nk)); wombat%zoopreflbac(:,:,:)=0.0
    allocate(wombat%zooprefsbac(isd:ied, jsd:jed, 1:nk)); wombat%zooprefsbac(:,:,:)=0.0
    allocate(wombat%zooprefaoa(isd:ied, jsd:jed, 1:nk)); wombat%zooprefaoa(:,:,:)=0.0
    allocate(wombat%zooprefphy(isd:ied, jsd:jed, 1:nk)); wombat%zooprefphy(:,:,:)=0.0
    allocate(wombat%zooprefdia(isd:ied, jsd:jed, 1:nk)); wombat%zooprefdia(:,:,:)=0.0
    allocate(wombat%zooprefsdet(isd:ied, jsd:jed, 1:nk)); wombat%zooprefsdet(:,:,:)=0.0
    allocate(wombat%zoograzlbac(isd:ied, jsd:jed, 1:nk)); wombat%zoograzlbac(:,:,:)=0.0
    allocate(wombat%zoograzsbac(isd:ied, jsd:jed, 1:nk)); wombat%zoograzsbac(:,:,:)=0.0
    allocate(wombat%zoograzaoa(isd:ied, jsd:jed, 1:nk)); wombat%zoograzaoa(:,:,:)=0.0
    allocate(wombat%zoograzphy(isd:ied, jsd:jed, 1:nk)); wombat%zoograzphy(:,:,:)=0.0
    allocate(wombat%zoograzdia(isd:ied, jsd:jed, 1:nk)); wombat%zoograzdia(:,:,:)=0.0
    allocate(wombat%zoograzsdet(isd:ied, jsd:jed, 1:nk)); wombat%zoograzsdet(:,:,:)=0.0
    allocate(wombat%zoomorl(isd:ied, jsd:jed, 1:nk)); wombat%zoomorl(:,:,:)=0.0
    allocate(wombat%zoomorq(isd:ied, jsd:jed, 1:nk)); wombat%zoomorq(:,:,:)=0.0
    allocate(wombat%zooexcrlbac(isd:ied, jsd:jed, 1:nk)); wombat%zooexcrlbac(:,:,:)=0.0
    allocate(wombat%zooexcrsbac(isd:ied, jsd:jed, 1:nk)); wombat%zooexcrsbac(:,:,:)=0.0
    allocate(wombat%zooexcraoa(isd:ied, jsd:jed, 1:nk)); wombat%zooexcraoa(:,:,:)=0.0
    allocate(wombat%zooexcrphy(isd:ied, jsd:jed, 1:nk)); wombat%zooexcrphy(:,:,:)=0.0
    allocate(wombat%zooexcrdia(isd:ied, jsd:jed, 1:nk)); wombat%zooexcrdia(:,:,:)=0.0
    allocate(wombat%zooexcrsdet(isd:ied, jsd:jed, 1:nk)); wombat%zooexcrsdet(:,:,:)=0.0
    allocate(wombat%zooegeslbac(isd:ied, jsd:jed, 1:nk)); wombat%zooegeslbac(:,:,:)=0.0
    allocate(wombat%zooegessbac(isd:ied, jsd:jed, 1:nk)); wombat%zooegessbac(:,:,:)=0.0
    allocate(wombat%zooegesaoa(isd:ied, jsd:jed, 1:nk)); wombat%zooegesaoa(:,:,:)=0.0
    allocate(wombat%zooegesphy(isd:ied, jsd:jed, 1:nk)); wombat%zooegesphy(:,:,:)=0.0
    allocate(wombat%zooegesdia(isd:ied, jsd:jed, 1:nk)); wombat%zooegesdia(:,:,:)=0.0
    allocate(wombat%zooegessdet(isd:ied, jsd:jed, 1:nk)); wombat%zooegessdet(:,:,:)=0.0
    allocate(wombat%meseps(isd:ied, jsd:jed, 1:nk)); wombat%meseps(:,:,:)=0.0
    allocate(wombat%mespreflbac(isd:ied, jsd:jed, 1:nk)); wombat%mespreflbac(:,:,:)=0.0
    allocate(wombat%mesprefsbac(isd:ied, jsd:jed, 1:nk)); wombat%mesprefsbac(:,:,:)=0.0
    allocate(wombat%mesprefaoa(isd:ied, jsd:jed, 1:nk)); wombat%mesprefaoa(:,:,:)=0.0
    allocate(wombat%mesprefphy(isd:ied, jsd:jed, 1:nk)); wombat%mesprefphy(:,:,:)=0.0
    allocate(wombat%mesprefdia(isd:ied, jsd:jed, 1:nk)); wombat%mesprefdia(:,:,:)=0.0
    allocate(wombat%mesprefsdet(isd:ied, jsd:jed, 1:nk)); wombat%mesprefsdet(:,:,:)=0.0
    allocate(wombat%mesprefldet(isd:ied, jsd:jed, 1:nk)); wombat%mesprefldet(:,:,:)=0.0
    allocate(wombat%mesprefzoo(isd:ied, jsd:jed, 1:nk)); wombat%mesprefzoo(:,:,:)=0.0
    allocate(wombat%mesgrazlbac(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazlbac(:,:,:)=0.0
    allocate(wombat%mesgrazsbac(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazsbac(:,:,:)=0.0
    allocate(wombat%mesgrazaoa(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazaoa(:,:,:)=0.0
    allocate(wombat%mesgrazphy(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazphy(:,:,:)=0.0
    allocate(wombat%mesgrazdia(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazdia(:,:,:)=0.0
    allocate(wombat%mesgrazsdet(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazsdet(:,:,:)=0.0
    allocate(wombat%mesgrazldet(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazldet(:,:,:)=0.0
    allocate(wombat%mesgrazzoo(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazzoo(:,:,:)=0.0
    allocate(wombat%mesmorl(isd:ied, jsd:jed, 1:nk)); wombat%mesmorl(:,:,:)=0.0
    allocate(wombat%mesmorq(isd:ied, jsd:jed, 1:nk)); wombat%mesmorq(:,:,:)=0.0
    allocate(wombat%mesexcrlbac(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrlbac(:,:,:)=0.0
    allocate(wombat%mesexcrsbac(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrsbac(:,:,:)=0.0
    allocate(wombat%mesexcraoa(isd:ied, jsd:jed, 1:nk)); wombat%mesexcraoa(:,:,:)=0.0
    allocate(wombat%mesexcrphy(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrphy(:,:,:)=0.0
    allocate(wombat%mesexcrdia(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrdia(:,:,:)=0.0
    allocate(wombat%mesexcrsdet(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrsdet(:,:,:)=0.0
    allocate(wombat%mesexcrldet(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrldet(:,:,:)=0.0
    allocate(wombat%mesexcrzoo(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrzoo(:,:,:)=0.0
    allocate(wombat%mesegeslbac(isd:ied, jsd:jed, 1:nk)); wombat%mesegeslbac(:,:,:)=0.0
    allocate(wombat%mesegessbac(isd:ied, jsd:jed, 1:nk)); wombat%mesegessbac(:,:,:)=0.0
    allocate(wombat%mesegesaoa(isd:ied, jsd:jed, 1:nk)); wombat%mesegesaoa(:,:,:)=0.0
    allocate(wombat%mesegesphy(isd:ied, jsd:jed, 1:nk)); wombat%mesegesphy(:,:,:)=0.0
    allocate(wombat%mesegesdia(isd:ied, jsd:jed, 1:nk)); wombat%mesegesdia(:,:,:)=0.0
    allocate(wombat%mesegessdet(isd:ied, jsd:jed, 1:nk)); wombat%mesegessdet(:,:,:)=0.0
    allocate(wombat%mesegesldet(isd:ied, jsd:jed, 1:nk)); wombat%mesegesldet(:,:,:)=0.0
    allocate(wombat%mesegeszoo(isd:ied, jsd:jed, 1:nk)); wombat%mesegeszoo(:,:,:)=0.0
    allocate(wombat%reminr(isd:ied, jsd:jed, 1:nk)); wombat%reminr(:,:,:)=0.0
    allocate(wombat%ldocremi(isd:ied, jsd:jed, 1:nk)); wombat%ldocremi(:,:,:)=0.0
    allocate(wombat%sdocremi(isd:ied, jsd:jed, 1:nk)); wombat%sdocremi(:,:,:)=0.0
    allocate(wombat%sdocprod(isd:ied, jsd:jed, 1:nk)); wombat%sdocprod(:,:,:)=0.0
    allocate(wombat%sdetremi(isd:ied, jsd:jed, 1:nk)); wombat%sdetremi(:,:,:)=0.0
    allocate(wombat%ldetremi(isd:ied, jsd:jed, 1:nk)); wombat%ldetremi(:,:,:)=0.0
    allocate(wombat%pic2poc(isd:ied, jsd:jed, 1:nk)); wombat%pic2poc(:,:,:)=0.0
    allocate(wombat%dissratcal(isd:ied, jsd:jed, 1:nk)); wombat%dissratcal(:,:,:)=0.0
    allocate(wombat%dissratara(isd:ied, jsd:jed, 1:nk)); wombat%dissratara(:,:,:)=0.0
    allocate(wombat%dissratpoc(isd:ied, jsd:jed, 1:nk)); wombat%dissratpoc(:,:,:)=0.0
    allocate(wombat%zoodiss(isd:ied, jsd:jed, 1:nk)); wombat%zoodiss(:,:,:)=0.0
    allocate(wombat%mesdiss(isd:ied, jsd:jed, 1:nk)); wombat%mesdiss(:,:,:)=0.0
    allocate(wombat%caldiss(isd:ied, jsd:jed, 1:nk)); wombat%caldiss(:,:,:)=0.0
    allocate(wombat%aradiss(isd:ied, jsd:jed, 1:nk)); wombat%aradiss(:,:,:)=0.0
    allocate(wombat%pocdiss(isd:ied, jsd:jed, 1:nk)); wombat%pocdiss(:,:,:)=0.0
    allocate(wombat%aoa_loxy(isd:ied, jsd:jed, 1:nk)); wombat%aoa_loxy(:,:,:)=0.0
    allocate(wombat%aoa_lnh4(isd:ied, jsd:jed, 1:nk)); wombat%aoa_lnh4(:,:,:)=0.0
    allocate(wombat%aoa_eno3(isd:ied, jsd:jed, 1:nk)); wombat%aoa_eno3(:,:,:)=0.0
    allocate(wombat%aoa_mumax(isd:ied, jsd:jed, 1:nk)); wombat%aoa_mumax(:,:,:)=0.0
    allocate(wombat%aoa_mu(isd:ied, jsd:jed, 1:nk)); wombat%aoa_mu(:,:,:)=0.0
    allocate(wombat%aoagrow(isd:ied, jsd:jed, 1:nk)); wombat%aoagrow(:,:,:)=0.0
    allocate(wombat%aoaresp(isd:ied, jsd:jed, 1:nk)); wombat%aoaresp(:,:,:)=0.0
    allocate(wombat%aoamorl(isd:ied, jsd:jed, 1:nk)); wombat%aoamorl(:,:,:)=0.0
    allocate(wombat%aoamorq(isd:ied, jsd:jed, 1:nk)); wombat%aoamorq(:,:,:)=0.0
    allocate(wombat%lbacydoc(isd:ied, jsd:jed, 1:nk)); wombat%lbacydoc(:,:,:)=0.0
    allocate(wombat%lbacgrow(isd:ied, jsd:jed, 1:nk)); wombat%lbacgrow(:,:,:)=0.0
    allocate(wombat%lbacresp(isd:ied, jsd:jed, 1:nk)); wombat%lbacresp(:,:,:)=0.0
    allocate(wombat%lbacpnh4(isd:ied, jsd:jed, 1:nk)); wombat%lbacpnh4(:,:,:)=0.0
    allocate(wombat%lbacpco2(isd:ied, jsd:jed, 1:nk)); wombat%lbacpco2(:,:,:)=0.0
    allocate(wombat%lbacufer(isd:ied, jsd:jed, 1:nk)); wombat%lbacufer(:,:,:)=0.0
    allocate(wombat%lbac_mu(isd:ied, jsd:jed, 1:nk)); wombat%lbac_mu(:,:,:)=0.0
    allocate(wombat%lbac_fanaer(isd:ied, jsd:jed, 1:nk)); wombat%lbac_fanaer(:,:,:)=0.0
    allocate(wombat%lbac_ffelim(isd:ied, jsd:jed, 1:nk)); wombat%lbac_ffelim(:,:,:)=0.0
    allocate(wombat%lbacmorl(isd:ied, jsd:jed, 1:nk)); wombat%lbacmorl(:,:,:)=0.0
    allocate(wombat%lbacmorq(isd:ied, jsd:jed, 1:nk)); wombat%lbacmorq(:,:,:)=0.0
    allocate(wombat%lbacdeni(isd:ied, jsd:jed, 1:nk)); wombat%lbacdeni(:,:,:)=0.0
    allocate(wombat%sbacydoc(isd:ied, jsd:jed, 1:nk)); wombat%sbacydoc(:,:,:)=0.0
    allocate(wombat%sbacgrow(isd:ied, jsd:jed, 1:nk)); wombat%sbacgrow(:,:,:)=0.0
    allocate(wombat%sbacresp(isd:ied, jsd:jed, 1:nk)); wombat%sbacresp(:,:,:)=0.0
    allocate(wombat%sbacpnh4(isd:ied, jsd:jed, 1:nk)); wombat%sbacpnh4(:,:,:)=0.0
    allocate(wombat%sbacpco2(isd:ied, jsd:jed, 1:nk)); wombat%sbacpco2(:,:,:)=0.0
    allocate(wombat%sbacufer(isd:ied, jsd:jed, 1:nk)); wombat%sbacufer(:,:,:)=0.0
    allocate(wombat%sbac_mu(isd:ied, jsd:jed, 1:nk)); wombat%sbac_mu(:,:,:)=0.0
    allocate(wombat%sbac_fanaer(isd:ied, jsd:jed, 1:nk)); wombat%sbac_fanaer(:,:,:)=0.0
    allocate(wombat%sbac_ffelim(isd:ied, jsd:jed, 1:nk)); wombat%sbac_ffelim(:,:,:)=0.0
    allocate(wombat%sbacmorl(isd:ied, jsd:jed, 1:nk)); wombat%sbacmorl(:,:,:)=0.0
    allocate(wombat%sbacmorq(isd:ied, jsd:jed, 1:nk)); wombat%sbacmorq(:,:,:)=0.0
    allocate(wombat%sbacdeni(isd:ied, jsd:jed, 1:nk)); wombat%sbacdeni(:,:,:)=0.0
    allocate(wombat%aox_lnh4(isd:ied, jsd:jed, 1:nk)); wombat%aox_lnh4(:,:,:)=0.0
    allocate(wombat%aox_mu(isd:ied, jsd:jed, 1:nk)); wombat%aox_mu(:,:,:)=0.0
    allocate(wombat%nitrfix(isd:ied, jsd:jed, 1:nk)); wombat%nitrfix(:,:,:)=0.0
    allocate(wombat%ammox(isd:ied, jsd:jed, 1:nk)); wombat%ammox(:,:,:)=0.0
    allocate(wombat%anammox(isd:ied, jsd:jed, 1:nk)); wombat%anammox(:,:,:)=0.0
    allocate(wombat%det_sed_remin(isd:ied, jsd:jed)); wombat%det_sed_remin(:,:)=0.0
    allocate(wombat%det_sed_denit(isd:ied, jsd:jed)); wombat%det_sed_denit(:,:)=0.0
    allocate(wombat%sdet_btm(isd:ied, jsd:jed)); wombat%sdet_btm(:,:)=0.0
    allocate(wombat%ldet_btm(isd:ied, jsd:jed)); wombat%ldet_btm(:,:)=0.0
    allocate(wombat%fbury(isd:ied, jsd:jed)); wombat%fbury(:,:)=0.0
    allocate(wombat%fdenit(isd:ied, jsd:jed)); wombat%fdenit(:,:)=0.0
    allocate(wombat%detfe_sed_remin(isd:ied, jsd:jed)); wombat%detfe_sed_remin(:,:)=0.0
    allocate(wombat%detsi_sed_remin(isd:ied, jsd:jed)); wombat%detsi_sed_remin(:,:)=0.0
    allocate(wombat%sdetfe_btm(isd:ied, jsd:jed)); wombat%sdetfe_btm(:,:)=0.0
    allocate(wombat%ldetfe_btm(isd:ied, jsd:jed)); wombat%ldetfe_btm(:,:)=0.0
    allocate(wombat%ldetsi_btm(isd:ied, jsd:jed)); wombat%ldetsi_btm(:,:)=0.0
    allocate(wombat%caco3_sed_remin(isd:ied, jsd:jed)); wombat%caco3_sed_remin(:,:)=0.0
    allocate(wombat%caco3_btm(isd:ied, jsd:jed)); wombat%caco3_btm(:,:)=0.0
    allocate(wombat%safe_btm(isd:ied, jsd:jed)); wombat%safe_btm(:,:)=0.0
    allocate(wombat%lafe_btm(isd:ied, jsd:jed)); wombat%lafe_btm(:,:)=0.0
    allocate(wombat%zw(isd:ied, jsd:jed, 1:nk)); wombat%zw(:,:,:)=0.0
    allocate(wombat%zm(isd:ied, jsd:jed, 1:nk)); wombat%zm(:,:,:)=0.0
    allocate(wombat%sdet_density(isd:ied, jsd:jed, 1:nk)); wombat%sdet_density(:,:,:)=0.0
    allocate(wombat%ldet_density(isd:ied, jsd:jed, 1:nk)); wombat%ldet_density(:,:,:)=0.0

    allocate(wombat%zeuphot(isd:ied, jsd:jed)); wombat%zeuphot(:,:)=0.0
    allocate(wombat%sdet_radius(isd:ied, jsd:jed)); wombat%sdet_radius(:,:)=0.0
    allocate(wombat%ldet_radius(isd:ied, jsd:jed)); wombat%ldet_radius(:,:)=0.0
    allocate(wombat%seddep(isd:ied, jsd:jed)); wombat%seddep(:,:)=0.0
    allocate(wombat%sedmask(isd:ied, jsd:jed)); wombat%sedmask(:,:)=0.0
    allocate(wombat%sedtemp(isd:ied, jsd:jed)); wombat%sedtemp(:,:)=0.0
    allocate(wombat%sedsalt(isd:ied, jsd:jed)); wombat%sedsalt(:,:)=0.0
    allocate(wombat%sedno3(isd:ied, jsd:jed)); wombat%sedno3(:,:)=0.0
    allocate(wombat%sednh4(isd:ied, jsd:jed)); wombat%sednh4(:,:)=0.0
    allocate(wombat%sedsil(isd:ied, jsd:jed)); wombat%sedsil(:,:)=0.0
    allocate(wombat%sedo2(isd:ied, jsd:jed)); wombat%sedo2(:,:)=0.0
    allocate(wombat%seddic(isd:ied, jsd:jed)); wombat%seddic(:,:)=0.0
    allocate(wombat%sedalk(isd:ied, jsd:jed)); wombat%sedalk(:,:)=0.0
    allocate(wombat%sedhtotal(isd:ied, jsd:jed)); wombat%sedhtotal(:,:)=0.0
    allocate(wombat%sedco3(isd:ied, jsd:jed)); wombat%sedco3(:,:)=0.0
    allocate(wombat%sedomega_cal(isd:ied, jsd:jed)); wombat%sedomega_cal(:,:)=0.0

  end subroutine user_allocate_arrays

  !#######################################################################
  !
  !   This is an internal sub, not a public interface.
  !   Deallocate all the work arrays allocated by user_allocate_arrays.
  !
  subroutine user_deallocate_arrays

    deallocate( &
        wombat%htotallo, &
        wombat%htotalhi, &
        wombat%htotal, &
        wombat%omega_ara, &
        wombat%omega_cal, &
        wombat%co3, &
        wombat%co2_star, &
        wombat%co2_csurf, &
        wombat%co2_alpha, &
        wombat%co2_sc_no, &
        wombat%pco2_csurf, &
        wombat%o2_csurf, &
        wombat%o2_alpha, &
        wombat%o2_sc_no, &
        wombat%no3_vstf, &
        wombat%nh4_vstf, &
        wombat%dic_vstf, &
        wombat%alk_vstf)

    deallocate( &
        wombat%f_dic, &
        wombat%f_dicr, &
        wombat%f_alk, &
        wombat%f_no3, &
        wombat%f_nh4, &
        wombat%f_sil, &
        wombat%f_phy, &
        wombat%f_pchl, &
        wombat%f_phyfe, &
        wombat%f_dia, &
        wombat%f_dchl, &
        wombat%f_diafe, &
        wombat%f_diasi, &
        wombat%f_zoo, &
        wombat%f_zoofe, &
        wombat%f_mes, &
        wombat%f_mesfe, &
        wombat%f_sdet, &
        wombat%f_sdetfe, &
        wombat%f_ldet, &
        wombat%f_ldetfe, &
        wombat%f_ldetsi, &
        wombat%f_ldoc, &
        wombat%f_sdoc, &
        wombat%f_ldon, &
        wombat%f_lbac, &
        wombat%f_sbac, &
        wombat%f_aoa, &
        wombat%f_o2, &
        wombat%f_caco3, &
        wombat%f_fe, &
        wombat%f_safe, &
        wombat%f_lafe)

    deallocate( &
        wombat%b_ldoc, &
        wombat%b_sdoc, &
        wombat%b_ldon, &
        wombat%b_no3, &
        wombat%b_o2, &
        wombat%b_dic, &
        wombat%b_fe, &
        wombat%b_sil, &
        wombat%b_alk)
    if (do_tracer_dicr) deallocate(wombat%b_dicr)

    deallocate( &
        wombat%radbio, &
        wombat%radmid, &
        wombat%radmld, &
        wombat%zsp3d, &
        wombat%rpp3d, &
        wombat%npp3d, &
        wombat%npp2d, &
        wombat%rpp2d, &
        wombat%zsp2d, &
        wombat%phy_mumax, &
        wombat%phy_mu, &
        wombat%pchl_mu, &
        wombat%phy_kni, &
        wombat%phy_kfe, &
        wombat%phy_lpar, &
        wombat%phy_lnit, &
        wombat%phy_lnh4, &
        wombat%phy_lno3, &
        wombat%phy_lfer, &
        wombat%phy_dfeupt, &
        wombat%dia_mumax, &
        wombat%dia_mu, &
        wombat%dchl_mu, &
        wombat%dia_kni, &
        wombat%dia_kfe, &
        wombat%dia_ksi, &
        wombat%dia_lpar, &
        wombat%dia_lnit, &
        wombat%dia_lnh4, &
        wombat%dia_lno3, &
        wombat%dia_lfer, &
        wombat%dia_lsil, &
        wombat%dia_dfeupt, &
        wombat%dia_silupt, &
        wombat%tri_lfer, &
        wombat%tri_lpar, &
        wombat%trimumax, &
        wombat%feIII, &
        wombat%sileqc, &
        wombat%disssi, &
        wombat%bsidiss, &
        wombat%felig, &
        wombat%ligK, &
        wombat%fecol, &
        wombat%fescasafe, &
        wombat%fescalafe, &
        wombat%fecoag2safe, &
        wombat%fecoag2lafe, &
        wombat%safediss, &
        wombat%lafediss, &
        wombat%fesources, &
        wombat%fesinks, &
        wombat%phy_feupreg, &
        wombat%phy_fedoreg, &
        wombat%phygrow, &
        wombat%phydoc, &
        wombat%phymorl, &
        wombat%phymorq, &
        wombat%dia_feupreg, &
        wombat%dia_fedoreg, &
        wombat%dia_sidoreg, &
        wombat%diagrow, &
        wombat%diadoc, &
        wombat%diamorl, &
        wombat%diamorq, &
        wombat%zoopreflbac, &
        wombat%zooprefsbac, &
        wombat%zooprefaoa, &
        wombat%zooprefphy, &
        wombat%zooprefdia, &
        wombat%zooprefsdet, &
        wombat%zoograzlbac, &
        wombat%zoograzsbac, &
        wombat%zoograzaoa, &
        wombat%zoograzphy, &
        wombat%zoograzdia, &
        wombat%zoograzsdet, &
        wombat%zoomorl, &
        wombat%zoomorq, &
        wombat%zooexcrlbac, &
        wombat%zooexcrsbac, &
        wombat%zooexcraoa, &
        wombat%zooexcrphy, &
        wombat%zooexcrdia, &
        wombat%zooexcrsdet, &
        wombat%zooegeslbac, &
        wombat%zooegessbac, &
        wombat%zooegesaoa, &
        wombat%zooegesphy, &
        wombat%zooegesdia, &
        wombat%zooegessdet, &
        wombat%mespreflbac, &
        wombat%mesprefsbac, &
        wombat%mesprefaoa, &
        wombat%mesprefphy, &
        wombat%mesprefdia, &
        wombat%mesprefsdet, &
        wombat%mesprefldet, &
        wombat%mesprefzoo, &
        wombat%mesgrazlbac, &
        wombat%mesgrazsbac, &
        wombat%mesgrazaoa, &
        wombat%mesgrazphy, &
        wombat%mesgrazdia, &
        wombat%mesgrazsdet, &
        wombat%mesgrazldet, &
        wombat%mesgrazzoo, &
        wombat%mesmorl, &
        wombat%mesmorq, &
        wombat%mesexcrlbac, &
        wombat%mesexcrsbac, &
        wombat%mesexcraoa, &
        wombat%mesexcrphy, &
        wombat%mesexcrdia, &
        wombat%mesexcrsdet, &
        wombat%mesexcrldet, &
        wombat%mesexcrzoo, &
        wombat%mesegeslbac, &
        wombat%mesegessbac, &
        wombat%mesegesaoa, &
        wombat%mesegesphy, &
        wombat%mesegesdia, &
        wombat%mesegessdet, &
        wombat%mesegesldet, &
        wombat%mesegeszoo, &
        wombat%reminr, &
        wombat%ldocremi, &
        wombat%sdocremi, &
        wombat%sdocprod, &
        wombat%sdetremi, &
        wombat%ldetremi, &
        wombat%pic2poc, &
        wombat%dissratcal, &
        wombat%dissratara, &
        wombat%dissratpoc, &
        wombat%zoodiss, &
        wombat%mesdiss, &
        wombat%caldiss, &
        wombat%aradiss, &
        wombat%pocdiss, &
        wombat%aoa_loxy, &
        wombat%aoa_lnh4, &
        wombat%aoa_eno3, &
        wombat%aoa_mumax, &
        wombat%aoa_mu, &
        wombat%aoagrow, &
        wombat%aoaresp, &
        wombat%aoamorl, &
        wombat%aoamorq, &
        wombat%lbacydoc, &
        wombat%lbacgrow, &
        wombat%lbacresp, &
        wombat%lbacpnh4, &
        wombat%lbacpco2, &
        wombat%lbacufer, &
        wombat%lbac_mu, &
        wombat%lbac_fanaer, &
        wombat%lbac_ffelim, &
        wombat%lbacmorl, &
        wombat%lbacmorq, &
        wombat%lbacdeni, &
        wombat%sbacydoc, &
        wombat%sbacgrow, &
        wombat%sbacresp, &
        wombat%sbacpnh4, &
        wombat%sbacpco2, &
        wombat%sbacufer, &
        wombat%sbac_mu, &
        wombat%sbac_fanaer, &
        wombat%sbac_ffelim, &
        wombat%sbacmorl, &
        wombat%sbacmorq, &
        wombat%sbacdeni, &
        wombat%aox_lnh4, &
        wombat%aox_mu, &
        wombat%nitrfix, &
        wombat%ammox, &
        wombat%anammox, &
        wombat%det_sed_remin, &
        wombat%det_sed_denit, &
        wombat%sdet_btm, &
        wombat%ldet_btm, &
        wombat%fbury, &
        wombat%fdenit, &
        wombat%detfe_sed_remin, &
        wombat%detsi_sed_remin, &
        wombat%sdetfe_btm, &
        wombat%ldetfe_btm, &
        wombat%ldetsi_btm, &
        wombat%caco3_sed_remin, &
        wombat%caco3_btm, &
        wombat%safe_btm, &
        wombat%lafe_btm, &
        wombat%dynvis_sw, &
        wombat%zw, &
        wombat%sdet_density, &
        wombat%ldet_density, &
        wombat%zm)

    deallocate( &
        wombat%zeuphot, &
        wombat%sdet_radius, &
        wombat%ldet_radius, &
        wombat%seddep, &
        wombat%sedmask, &
        wombat%sedtemp, &
        wombat%sedsalt, &
        wombat%sedno3, &
        wombat%sednh4, &
        wombat%sedsil, &
        wombat%sedo2, &
        wombat%seddic, &
        wombat%sedalk, &
        wombat%sedhtotal, &
        wombat%sedco3, &
        wombat%sedomega_cal)

  end subroutine user_deallocate_arrays

end module generic_WOMBATmid
