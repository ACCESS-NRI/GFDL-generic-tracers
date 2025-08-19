!-----------------------------------------------------------------------
!
! <CONTACT EMAIL="Richard.Matear@csiro.au"> Richard Matear
! </CONTACT>
!
! <CONTACT EMAIL="Matthew.Chamberlain@csiro.au"> Matt Chamberlain
! </CONTACT>
!
! <CONTACT EMAIL="Dougal.Squire@anu.edu.au"> Dougie Squire
! </CONTACT>
!
! <CONTACT EMAIL="Pearse.Buchanan@csiro.au"> Pearse Buchanan
! </CONTACT>
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
!  phytoplankton, zooplankton and sinking detritus, as well as nitrate 
!  (NO3), ammonium (NH4), nitrous oxide (N2O), dissolved iron (Fe), 
!  silicic acid (SIL), dissolved organic matter that is split into carbon 
!  (DOC) and nitrogen (DON), two explicit heterotrophic bacterial types 
!  (BAC1 & BAC2) and ammonia oxidizing archaea (AOA), dissolved inorganic 
!  carbon (DIC), calcium carbonate (CaCO3), alkalinity (ALK), and oxygen 
!  (O2). Fe is carried through all exosystem biomass pools except bacteria 
!  and AOA, who have constant C:N:Fe ratios. C:N ratios are fixed in all 
!  biomass pools except for dissolved organics, since we represent both DOC 
!  and DON. Si is carried through microphytoplankton, large detrtis and, 
!  like for carbon and Fe, is deposited into a sediment pool.
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
!  <DATA NAME="do_burial" TYPE="logical">
!   If true, permanently bury organics and CaCO3 in sediments
!  </DATA>
!
!  <DATA NAME="do_conserve_tracers" TYPE="logical">
!   If true and do_burial is true, add back the lost NO3 and Alk due to
!   burial to surface
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
  logical :: do_caco3_dynamics   = .true.  ! do dynamic CaCO3 precipitation, dissolution and ballasting?
  logical :: do_burial           = .false. ! permanently bury organics and CaCO3 in sediments?
  logical :: do_conserve_tracers = .false. ! add back the lost NO3 and Alk due to burial to surface?
  logical :: do_nitrogen_fixation= .true.  ! N cycle has nitrogen fixation?
  logical :: do_anammox          = .true.  ! N cycle has anammox?
  logical :: do_wc_denitrification = .true.  ! N cycle has water column denitrification?
  logical :: do_benthic_denitrification = .true.  ! N cycle has N loss in sediments?
  logical :: do_check_n_conserve = .false. ! check that the N fluxes balance in the ecosystem
  logical :: do_check_c_conserve = .true.  ! check that the C fluxes balance in the ecosystem
  logical :: do_check_si_conserve = .true.  ! check that the Si fluxes balance in the ecosystem

  namelist /generic_wombatmid_nml/ co2_calc, do_caco3_dynamics, do_burial, do_conserve_tracers, &
                                   do_nitrogen_fixation, do_anammox, do_wc_denitrification, do_benthic_denitrification, &
                                   do_check_n_conserve, do_check_c_conserve, do_check_si_conserve

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
        phyminqc, &
        phyoptqc, &
        phyoptqf, &
        phymaxqf, &
        phylmor, &
        phyqmor, &
        diakn, &
        diakf, &
        diaminqc, &
        diaoptqc, &
        diaoptqf, &
        diamaxqf, &
        dialmor, &
        diaqmor, &
        chlkWm2, &
        overflow, &
        trikf, &
        trichlc, &
        trin2c, &
        zooCingest, &
        zooCassim, &
        zooFeingest, &
        zooFeassim, &
        zooexcrdom, &
        zookz, &
        zoogmax, &
        zooepsbac1, &
        zooepsbac2, &
        zooepsaoa, &
        zooepsphy, &
        zooepsdia, &
        zooepsdet, &
        zprefbac1, &
        zprefbac2, &
        zprefaoa, &
        zprefphy, &
        zprefdia, &
        zprefdet, &
        zoolmor, &
        zooqmor, &
        mesCingest, &
        mesCassim, &
        mesFeingest, &
        mesFeassim, &
        mesexcrdom, &
        meskz, &
        mesgmax, &
        mesepsbac1, &
        mesepsbac2, &
        mesepsaoa, &
        mesepsphy, &
        mesepsdia, &
        mesepsdet, &
        mesepsbdet, &
        mesepszoo, &
        mprefbac1, &
        mprefbac2, &
        mprefaoa, &
        mprefphy, &
        mprefdia, &
        mprefdet, &
        mprefbdet, &
        mprefzoo, &
        meslmor, &
        mesqmor, &
        zoopreyswitch, &
        mespreyswitch, &
        detlrem, &
        bottom_thickness, &
        detlrem_sed, &
        wdetbio, &
        wbdetbio, &
        wdetmax, &
        phybiot, &
        diabiot, &
        wcaco3, &
        caco3lrem, &
        caco3lrem_sed, &
        omegamax_sed, &
        f_inorg, &
        disscal, &
        dissara, &
        dissdet, &
        ligand, &
        fcolloid, &
        knano_dfe, &
        kscav_dfe, &
        kcoag_dfe, &
        bsi_fbac, &
        bsi_kbac, &
        aoa_knh4, &
        aoa_poxy, &
        aoa_ynh4, &
        aoa_yoxy, &
        aoa_C2N, &
        aoa_C2Fe, &
        aoalmor, &
        aoaqmor, &
        bac1_Vmax_doc, &
        bac1_Vmax_no3, &
        bac1_poxy, &
        bac1_kno3, &
        bac1_kdoc_min, &
        bac1_kdoc_max, &
        bac1_knh4, &
        bac1_kfer, &
        bac1_yoxy, &
        bac1_yaerC, &
        bac1_yno3, &
        bac1_yanaC, &
        bac1_C2N, &
        bac1_C2Fe, &
        bac1lmor, &
        bac1qmor, &
        bac2_Vmax_doc, &
        bac2_poxy, &
        bac2_pn2o, &
        bac2_kdoc_min, &
        bac2_kdoc_max, &
        bac2_knh4, &
        bac2_kfer, &
        bac2_yoxy, &
        bac2_yaerC, &
        bac2_yn2o, &
        bac2_yanaC, &
        bac2_C2N, &
        bac2_C2Fe, &
        bac2lmor, &
        bac2qmor, &
        aoxkn, &
        aoxmumax, &
        dt_npzd, &
        sal_global, &
        dic_global, &
        alk_global, &
        no3_global, &
        nh4_global, &
        sio2_surf, &
        dic_min, &
        dic_max, &
        alk_min, &
        alk_max, &
        htotal_scale_lo, &
        htotal_scale_hi, &
        htotal_in, &
        Rho_0, &
        a_0, a_1, a_2, a_3, a_4, a_5, &
        b_0, b_1, b_2, b_3, c_0, &
        a1_co2, a2_co2, a3_co2, a4_co2, a5_co2, &
        a1_o2, a2_o2, a3_o2, a4_o2, a5_o2, &
        a_1_n2o, a_2_n2o, a_3_n2o, a_4_n2o, &
        b_1_n2o, b_2_n2o, b_3_n2o, &
        a1_n2o, a2_n2o, a3_n2o, a4_n2o, a5_n2o


    character(len=fm_string_len) :: ice_restart_file
    character(len=fm_string_len) :: ocean_restart_file
    character(len=fm_string_len) :: IC_file

    !-----------------------------------------------------------------------
    ! Arrays for surface gas fluxes
    !-----------------------------------------------------------------------
    real, dimension(:,:), allocatable :: &
        htotallo, htotalhi,  &
        sio2, &
        co2_csurf, co2_alpha, co2_sc_no, pco2_csurf, &
        o2_csurf, o2_alpha, o2_sc_no, &
        n2o_csurf, n2o_alpha, n2o_sc_no, &
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
        b_nh4, &
        b_no3, &
        b_sil, &
        b_o2, &
        b_fe, &
        pprod_gross_2d, &
        export_prod, &
        export_inorg, &
        npp2d, &
        det_btm, &
        detfe_btm, &
        bdet_btm, &
        bdetfe_btm, &
        bdetsi_btm, &
        caco3_btm, &
        det_sed_remin, &
        detfe_sed_remin, &
        detsi_sed_remin, &
        caco3_sed_remin, &
        det_sed_denit, &
        fbury, &
        fdenit, &
        dic_intmld, &
        o2_intmld, &
        no3_intmld, &
        fe_intmld, &
        phy_intmld, &
        det_intmld, &
        pprod_gross_intmld, &
        npp_intmld, &
        radbio_intmld, &
        dic_int100, &
        o2_int100, &
        no3_int100, &
        fe_int100, &
        phy_int100, &
        det_int100, &
        pprod_gross_int100, &
        npp_int100, &
        radbio_int100, &
        zeuphot, &
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
        f_det, &
        f_detfe, &
        f_bdet, &
        f_bdetfe, &
        f_bdetsi, &
        f_doc, &
        f_don, &
        f_bac1, &
        f_bac2, &
        f_aoa, &
        f_n2o, &
        f_o2, &
        f_caco3, &
        f_fe, &
        pprod_gross, &
        zprod_gross, &
        radbio, &
        radmid, &
        radmld, &
        npp3d, &
        phy_mumax, &
        phy_mu, &
        pchl_mu, &
        pchl_lpar, &
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
        dchl_lpar, &
        dia_kni, &
        dia_kfe, &
        dia_lpar, &
        dia_lnit, &
        dia_lnh4, &
        dia_lno3, &
        dia_lfer, &
        dia_dfeupt, &
        trimumax, &
        tri_lfer, &
        tri_lpar, &
        sileqc, &
        disssi, &
        bsidiss, &
        feIII, &
        felig, &
        fecol, &
        feprecip, &
        fescaven, &
        fescadet, &
        fescabdet, &
        fecoag2det, &
        fecoag2bdet, &
        fesources, &
        fesinks, &
        phy_feupreg, &
        phy_fedoreg, &
        phygrow, &
        phydoc, &
        phyresp, &
        phymort, &
        dia_feupreg, &
        dia_fedoreg, &
        diagrow, &
        diadoc, &
        diaresp, &
        diamort, &
        zooeps, &
        zooprefbac1, &
        zooprefbac2, &
        zooprefaoa, &
        zooprefphy, &
        zooprefdia, &
        zooprefdet, &
        zoograzbac1, &
        zoograzbac2, &
        zoograzaoa, &
        zoograzphy, &
        zoograzdia, &
        zoograzdet, &
        zooresp, &
        zoomort, &
        zooexcrbac1, &
        zooexcrbac2, &
        zooexcraoa, &
        zooexcrphy, &
        zooexcrdia, &
        zooexcrdet, &
        zooegesbac1, &
        zooegesbac2, &
        zooegesaoa, &
        zooegesphy, &
        zooegesdia, &
        zooegesdet, &
        meseps, &
        mesprefbac1, &
        mesprefbac2, &
        mesprefaoa, &
        mesprefphy, &
        mesprefdia, &
        mesprefdet, &
        mesprefbdet, &
        mesprefzoo, &
        mesgrazbac1, &
        mesgrazbac2, &
        mesgrazaoa, &
        mesgrazphy, &
        mesgrazdia, &
        mesgrazdet, &
        mesgrazbdet, &
        mesgrazzoo, &
        mesresp, &
        mesmort, &
        mesexcrbac1, &
        mesexcrbac2, &
        mesexcraoa, &
        mesexcrphy, &
        mesexcrdia, &
        mesexcrdet, &
        mesexcrbdet, &
        mesexcrzoo, &
        mesegesbac1, &
        mesegesbac2, &
        mesegesaoa, &
        mesegesphy, &
        mesegesdia, &
        mesegesdet, &
        mesegesbdet, &
        mesegeszoo, &
        reminr, &
        doc1remi, &
        don1remi, &
        doc2remi, &
        don2remi, &
        detremi, &
        bdetremi, &
        pic2poc, &
        dissrat, &
        caldiss, &
        aoa_loxy, &
        aoa_lnh4, &
        aoa_yn2o, &
        aoa_mumax, &
        aoa_mu, &
        aoagrow, &
        aoaresp, &
        aoamor1, &
        aoamor2, &
        bac1grow, &
        bac1resp, &
        bac1unh4, &
        bac1ufer, &
        bac1_lnit, &
        bac1_lfer, &
        bac1_mu, &
        bac1_kdoc, &
        bac1_fanaer, &
        bac1mor1, &
        bac1mor2, &
        bac1deni, &
        bac2grow, &
        bac2resp, &
        bac2unh4, &
        bac2ufer, &
        bac2_lnit, &
        bac2_lfer, &
        bac2_mu, &
        bac2_kdoc, &
        bac2_fanaer, &
        bac2mor1, &
        bac2mor2, &
        bac2deni, &
        aox_lnh4, &
        aox_mu, &
        nitrfix, &
        ammox, &
        anammox, &
        no3_prev, &
        caco3_prev, &
        dic_correct, &
        alk_correct, &
        zw, &
        zm

    real, dimension(:,:,:,:), pointer :: &
        p_o2, p_n2o

    real, dimension(:,:,:), pointer :: &
        p_det_sediment, &
        p_detfe_sediment, &
        p_detsi_sediment, &
        p_caco3_sediment, &
        p_detbury, p_caco3bury

    real, dimension(:,:,:), pointer :: &
        p_wdet, &
        p_wdetfe, &
        p_wbdet, &
        p_wbdetfe, &
        p_wbdetsi, &
        p_wcaco3

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
        id_dic_correct = -1, &
        id_alk_correct = -1, &
        id_radbio = -1, &
        id_radmid = -1, &
        id_radmld = -1, &
        id_radbio1 = -1, &
        id_pprod_gross = -1, &
        id_phy_kni = -1, &
        id_phy_kfe = -1, &
        id_phy_lpar = -1, &
        id_pchl_lpar = -1, &
        id_phy_lnit = -1, &
        id_phy_lnh4 = -1, &
        id_phy_lno3 = -1, &
        id_phy_lfer = -1, &
        id_phy_dfeupt = -1, &
        id_dia_kni = -1, &
        id_dia_kfe = -1, &
        id_dia_lpar = -1, &
        id_dchl_lpar = -1, &
        id_dia_lnit = -1, &
        id_dia_lnh4 = -1, &
        id_dia_lno3 = -1, &
        id_dia_lfer = -1, &
        id_dia_dfeupt = -1, &
        id_trimumax = -1, &
        id_tri_lfer = -1, &
        id_tri_lpar = -1, &
        id_sileqc = -1, &
        id_disssi = -1, &
        id_bsidiss = -1, &
        id_feIII = -1, &
        id_felig = -1, &
        id_fecol = -1, &
        id_feprecip = -1, &
        id_fescaven = -1, &
        id_fescadet = -1, &
        id_fescabdet = -1, &
        id_fecoag2det = -1, &
        id_fecoag2bdet = -1, &
        id_fesources = -1, &
        id_fesinks = -1, &
        id_phy_feupreg = -1, &
        id_phy_fedoreg = -1, &
        id_phygrow = -1, &
        id_phydoc = -1, &
        id_phyresp = -1, &
        id_phymort = -1, &
        id_dia_feupreg = -1, &
        id_dia_fedoreg = -1, &
        id_diagrow = -1, &
        id_diadoc = -1, &
        id_diaresp = -1, &
        id_diamort = -1, &
        id_zooeps = -1, &
        id_zooprefbac1 = -1, &
        id_zooprefbac2 = -1, &
        id_zooprefaoa = -1, &
        id_zooprefphy = -1, &
        id_zooprefdia = -1, &
        id_zooprefdet = -1, &
        id_zoograzbac1 = -1, &
        id_zoograzbac2 = -1, &
        id_zoograzaoa = -1, &
        id_zoograzphy = -1, &
        id_zoograzdia = -1, &
        id_zoograzdet = -1, &
        id_zooresp = -1, &
        id_zoomort = -1, &
        id_zooexcrbac1 = -1, &
        id_zooexcrbac2 = -1, &
        id_zooexcraoa = -1, &
        id_zooexcrphy = -1, &
        id_zooexcrdia = -1, &
        id_zooexcrdet = -1, &
        id_zooegesbac1 = -1, &
        id_zooegesbac2 = -1, &
        id_zooegesaoa = -1, &
        id_zooegesphy = -1, &
        id_zooegesdia = -1, &
        id_zooegesdet = -1, &
        id_meseps = -1, &
        id_mesprefbac1 = -1, &
        id_mesprefbac2 = -1, &
        id_mesprefaoa = -1, &
        id_mesprefphy = -1, &
        id_mesprefdia = -1, &
        id_mesprefdet = -1, &
        id_mesprefbdet = -1, &
        id_mesprefzoo = -1, &
        id_mesgrazbac1 = -1, &
        id_mesgrazbac2 = -1, &
        id_mesgrazaoa = -1, &
        id_mesgrazphy = -1, &
        id_mesgrazdia = -1, &
        id_mesgrazdet = -1, &
        id_mesgrazbdet = -1, &
        id_mesgrazzoo = -1, &
        id_mesresp = -1, &
        id_mesmort = -1, &
        id_mesexcrbac1 = -1, &
        id_mesexcrbac2 = -1, &
        id_mesexcraoa = -1, &
        id_mesexcrphy = -1, &
        id_mesexcrdia = -1, &
        id_mesexcrdet = -1, &
        id_mesexcrbdet = -1, &
        id_mesexcrzoo = -1, &
        id_mesegesbac1 = -1, &
        id_mesegesbac2 = -1, &
        id_mesegesaoa = -1, &
        id_mesegesphy = -1, &
        id_mesegesdia = -1, &
        id_mesegesdet = -1, &
        id_mesegesbdet = -1, &
        id_mesegeszoo = -1, &
        id_reminr = -1, &
        id_doc1remi = -1, &
        id_don1remi = -1, &
        id_doc2remi = -1, &
        id_don2remi = -1, &
        id_detremi = -1, &
        id_bdetremi = -1, &
        id_pic2poc = -1, &
        id_dissrat = -1, &
        id_caldiss = -1, &
        id_aoa_loxy = -1, &
        id_aoa_lnh4 = -1, &
        id_aoa_yn2o = -1, &
        id_aoa_mumax = -1, &
        id_aoa_mu = -1, &
        id_aoagrow = -1, &
        id_aoaresp = -1, &
        id_aoamor1 = -1, &
        id_aoamor2 = -1, &
        id_bac1grow = -1, &
        id_bac1resp = -1, &
        id_bac1unh4 = -1, &
        id_bac1ufer = -1, &
        id_bac1_lnit = -1, &
        id_bac1_lfer = -1, &
        id_bac1_mu = -1, &
        id_bac1_kdoc = -1, &
        id_bac1_fanaer = -1, &
        id_bac1mor1 = -1, &
        id_bac1mor2 = -1, &
        id_bac1deni = -1, &
        id_bac2grow = -1, &
        id_bac2resp = -1, &
        id_bac2unh4 = -1, &
        id_bac2ufer = -1, &
        id_bac2_lnit = -1, &
        id_bac2_lfer = -1, &
        id_bac2_mu = -1, &
        id_bac2_kdoc = -1, &
        id_bac2_fanaer = -1, &
        id_bac2mor1 = -1, &
        id_bac2mor2 = -1, &
        id_bac2deni = -1, &
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
        id_pprod_gross_2d = -1, &
        id_export_prod = -1, &
        id_export_inorg = -1, &
        id_npp3d = -1, &
        id_npp2d = -1, &
        id_npp1 = -1, &
        id_zprod_gross = -1, &
        id_dic_intmld = -1, &
        id_o2_intmld = -1, &
        id_no3_intmld = -1, &
        id_fe_intmld = -1, &
        id_phy_intmld = -1, &
        id_det_intmld = -1, &
        id_pprod_gross_intmld = -1, &
        id_npp_intmld = -1, &
        id_radbio_intmld = -1, &
        id_dic_int100 = -1, &
        id_o2_int100 = -1, &
        id_no3_int100 = -1, &
        id_fe_int100 = -1, &
        id_phy_int100 = -1, &
        id_det_int100 = -1, &
        id_pprod_gross_int100 = -1, &
        id_npp_int100 = -1, &
        id_radbio_int100 = -1, &
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
    type(g_tracer_type), pointer :: tracer_list

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
          'Doing dynamic CaCO3 precipitation, dissolution and ballasting'
    endif

    if (do_burial) then
      write (stdoutunit,*) trim(note_header), &
          'Permanently burying organics and CaCO3 in sediments'
      if (do_conserve_tracers) then
        write (stdoutunit,*) trim(note_header), &
          'Adding back the lost NO3 and Alk due to burial to surface'
      endif
    elseif (do_conserve_tracers) then
      call mpp_error(WARNING, trim(warn_header) // &
          'do_conserve_tracers = .true. is doing nothing because do_burial = .false.')
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
    type(g_tracer_type), pointer :: tracer_list
    logical, intent(in)          :: force_update_fluxes

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
    type(g_diag_type), pointer :: diag_list ! dts: this is not actually used

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

    vardesc_temp = vardesc( &
        'dicp_vstf', 'Virtual flux of preformed dissolved inorganic carbon into ocean due to '// &
        'salinity restoring/correction', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_dicp_vstf = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'alk_vstf', 'Virtual flux of alkalinity into ocean due to salinity restoring/correction', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_alk_vstf = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    !=======================================================================
    ! Tracer and source term diagnostics
    !=======================================================================
    vardesc_temp = vardesc( &
        'dic_correct', 'Correction to DIC tracer (min limit)',  'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_dic_correct = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'alk_correct', 'Correction to Alk tracer (min limit)',  'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_alk_correct = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'radbio1', 'Photosynthetically active radiation for phytoplankton growth at surface', &
        'h', '1', 's', 'W m-2', 'f')
    wombat%id_radbio1 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
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
        'export_prod', 'Organic export through 100m', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_export_prod = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'export_inorg', 'Inorganic export through 100m', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_export_inorg = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'npp3d', 'Net primary productivity', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_npp3d = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'npp2d', 'Vertically integrated net primary productivity', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_npp2d = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
  
    vardesc_temp = vardesc( &
        'npp1', 'Net primary productivity in the first ocean layer', &
        'h', '1', 's', 'mol/kg/s', 'f')
    wombat%id_npp1 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
      
    vardesc_temp = vardesc( &
        'pprod_gross', 'Gross phytoplankton production', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_pprod_gross = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'pchl_lpar', 'Limitation of phytoplankton chlorophyll production by light', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_pchl_lpar = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'dchl_lpar', 'Limitation of microphytoplankton chlorophyll production by light', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_dchl_lpar = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'dia_dfeupt', 'Uptake of dFe by microphytoplankton', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_dia_dfeupt = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'tri_lfer', 'Limitation of trichodesmium by light', 'h', 'L', 's', '[0-1]', 'f')
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
        'fecol', 'Colloidal dissolved iron', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_fecol = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'feprecip', 'Precipitation of free Fe onto nanoparticles', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_feprecip = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fescaven', 'Scavenging of free Fe onto detritus (organic + inorganic)', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fescaven = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fescadet', 'Scavenging of free Fe onto organic detritus', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fescadet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fescabdet', 'Scavenging of free Fe onto organic big detritus', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fescabdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fecoag2det', 'Coagulation of colloidal dFe onto detritus', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fecoag2det = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fecoag2bdet', 'Coagulation of colloidal dFe onto big detritus', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fecoag2bdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'phyresp', 'Respiration of phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_phyresp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phymort', 'Mortality of phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_phymort = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'diagrow', 'Growth of microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_diagrow = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'diadoc', 'Overflow exudation of DOC by microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_diadoc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'diaresp', 'Respiration of microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_diaresp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'diamort', 'Mortality of microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_diamort = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooeps', 'Zooplankton community-wide prey capture rate coefficient', 'h', 'L', 's', 'm^6/mmolC^2/s', 'f')
    wombat%id_zooeps = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooprefbac1', 'Grazing dietary fraction of zooplankton on bacteria 1', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooprefbac1 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooprefbac2', 'Grazing dietary fraction of zooplankton on bacteria 2', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooprefbac2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'zooprefdet', 'Grazing dietary fraction of zooplankton on detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooprefdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoograzbac1', 'Grazing rate of zooplankton on bacteria 1', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoograzbac1 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoograzbac2', 'Grazing rate of zooplankton on bacteria 2', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoograzbac2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'zoograzdet', 'Grazing rate of zooplankton on detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoograzdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooresp', 'Respiration of zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooresp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoomort', 'Mortality of zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoomort = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooexcrbac1', 'Excretion rate of zooplankton eating bacteria 1', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooexcrbac1 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooexcrbac2', 'Excretion rate of zooplankton eating bacteria 2', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooexcrbac2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'zooexcrdet', 'Excretion rate of zooplankton eating detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooexcrdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooegesbac1', 'Egestion rate of zooplankton on bacteria 1', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooegesbac1 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooegesbac2', 'Egestion rate of zooplankton on bacteria 2', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooegesbac2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'zooegesdet', 'Egestion rate of zooplankton on detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooegesdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'meseps', 'mesozooplankton community-wide prey capture rate coefficient', 'h', 'L', 's', 'm^6/mmolC^2/s', 'f')
    wombat%id_meseps = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesprefbac1', 'Grazing dietary fraction of mesozooplankton on bacteria 1', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesprefbac1 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesprefbac2', 'Grazing dietary fraction of mesozooplankton on bacteria 2', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesprefbac2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'mesprefdet', 'Grazing dietary fraction of mesozooplankton on detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesprefdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesprefbdet', 'Grazing dietary fraction of mesozooplankton on big detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesprefbdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesprefzoo', 'Grazing dietary fraction of mesozooplankton on zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesprefzoo = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesgrazbac1', 'Grazing rate of mesozooplankton on bacteria 1', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesgrazbac1 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesgrazbac2', 'Grazing rate of mesozooplankton on bacteria 2', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesgrazbac2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'mesgrazdet', 'Grazing rate of mesozooplankton on detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesgrazdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesgrazbdet', 'Grazing rate of mesozooplankton on big detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesgrazbdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesgrazzoo', 'Grazing rate of mesozooplankton on zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesgrazzoo = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesresp', 'Respiration of mesozooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesresp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesmort', 'Mortality of mesozooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesmort = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesexcrbac1', 'Excretion rate of mesozooplankton eating bacteria 1', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesexcrbac1 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesexcrbac2', 'Excretion rate of mesozooplankton eating bacteria 2', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesexcrbac2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'mesexcrdet', 'Excretion rate of mesozooplankton eating detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesexcrdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesexcrbdet', 'Excretion rate of mesozooplankton eating big detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesexcrbdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesexcrzoo', 'Excretion rate of mesozooplankton eating zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesexcrzoo = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesegesbac1', 'Egestion rate of mesozooplankton on bacteria 1', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesegesbac1 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesegesbac2', 'Egestion rate of mesozooplankton on bacteria 2', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesegesbac2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'mesegesdet', 'Egestion rate of mesozooplankton on detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesegesdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesegesbdet', 'Egestion rate of mesozooplankton on big detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_mesegesdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'doc1remi', 'Remineralisation of dissolved organic carbon by bacteria #1', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_doc1remi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'don1remi', 'Remineralisation of dissolved organic nitrogen by bacteria #1', 'h', 'L', 's', 'molN/kg/s', 'f')
    wombat%id_don1remi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'doc2remi', 'Remineralisation of dissolved organic carbon by bacteria #2', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_doc2remi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'don2remi', 'Remineralisation of dissolved organic nitrogen by bacteria #2', 'h', 'L', 's', 'molN/kg/s', 'f')
    wombat%id_don2remi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'detremi', 'Hydrolysation of detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_detremi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bdetremi', 'Hydrolysation of big detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_bdetremi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'pic2poc', 'Inorganic (CaCO3) to organic carbon ratio', 'h', 'L', 's', ' ', 'f')
    wombat%id_pic2poc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dissrat', 'Dissolution rate of CaCO3', 'h', 'L', 's', '/s', 'f')
    wombat%id_dissrat = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'caldiss', 'Dissolution of CaCO3', 'h', 'L', 's', 'molCaCO3/kg/s', 'f')
    wombat%id_caldiss = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'aoa_yn2o', 'Yield of N2O produced by Ammonia Oxidizing Archaea during oxidation', 'h', 'L', 's', 'mol N / mol Biomass', 'f')
    wombat%id_aoa_yn2o = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'aoamor1', 'Linear mortality of Ammonia Oxidizing Archaea', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_aoamor1 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'aoamor2', 'Quadratic mortality of Ammonia Oxidizing Archaea', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_aoamor2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac1grow', 'Growth of facultative heterotrophic bacteria 1', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_bac1grow = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac1resp', 'Oxygen consumption of facultative heterotrophic bacteria 1', 'h', 'L', 's', 'molO2/kg/s', 'f')
    wombat%id_bac1resp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac1unh4', 'Uptake of NH4 of facultative heterotrophic bacteria 1', 'h', 'L', 's', 'molNH4/kg/s', 'f')
    wombat%id_bac1unh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac1ufer', 'Uptake of dFe of facultative heterotrophic bacteria 1', 'h', 'L', 's', 'moldFe/kg/s', 'f')
    wombat%id_bac1ufer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac1_lfer', 'Limitation of facultative heterotrophic bacteria 1 by dFe', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_bac1_lfer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac1_lnit', 'Limitation of facultative heterotrophic bacteria 1 by DON and NH4', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_bac1_lnit = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac1_mu', 'Realized growth rate of facultative heterotrophic bacteria 1', 'h', 'L', 's', '/s', 'f')
    wombat%id_bac1_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac1_kdoc', 'Half saturation coefficient for DOC uptake by facultative heterotrophic bacteria 1', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_bac1_kdoc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac1_fanaer', 'Fraction of growth supported by anaerobic metabolism', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_bac1_fanaer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac1mor1', 'Linear mortality of facultative heterotrophic bacteria 1', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_bac1mor1 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac1mor2', 'Quadratic mortality of facultative heterotrophic bacteria 1', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_bac1mor2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac1deni', 'bacterial denitrification rate (NO3 consumption)', 'h', 'L', 's', '[molN/kg/s]', 'f')
    wombat%id_bac1deni = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac2grow', 'Growth of facultative heterotrophic bacteria 2', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_bac2grow = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac2resp', 'Oxygen consumption of facultative heterotrophic bacteria 2', 'h', 'L', 's', 'molO2/kg/s', 'f')
    wombat%id_bac2resp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac2unh4', 'Uptake of NH4 of facultative heterotrophic bacteria 2', 'h', 'L', 's', 'molNH4/kg/s', 'f')
    wombat%id_bac2unh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac2ufer', 'Uptake of dFe of facultative heterotrophic bacteria 2', 'h', 'L', 's', 'moldFe/kg/s', 'f')
    wombat%id_bac2ufer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac2_lfer', 'Limitation of facultative heterotrophic bacteria 2 by dFe', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_bac2_lfer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac2_lnit', 'Limitation of facultative heterotrophic bacteria 2 by DON and NH4', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_bac2_lnit = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac2_mu', 'Realized growth rate of facultative heterotrophic bacteria 2', 'h', 'L', 's', '/s', 'f')
    wombat%id_bac2_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac2_kdoc', 'Half saturation coefficient for DOC uptake by facultative heterotrophic bacteria 2', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_bac2_kdoc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac2_fanaer', 'Fraction of growth supported by anaerobic metabolism', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_bac2_fanaer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac2mor1', 'Linear mortality of facultative heterotrophic bacteria 2', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_bac2mor1 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac2mor2', 'Quadratic mortality of facultative heterotrophic bacteria 2', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_bac2mor2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'bac2deni', 'bacterial denitrification rate (N2O consumption)', 'h', 'L', 's', '[molN2/kg/s]', 'f')
    wombat%id_bac2deni = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'pprod_gross_2d', 'Vertically integrated gross phytoplankton production', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_pprod_gross_2d = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zprod_gross', 'Gross zooplankton production', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_zprod_gross = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    ! MLD-integrated diagnostics
    !-----------------------------------------------------------------------
    vardesc_temp = vardesc( &
        'dic_intmld', 'MLD-integrated dissolved inorganic carbon', &
        'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_dic_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'o2_intmld', 'MLD-integrated dissolved oxygen', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_o2_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'no3_intmld', 'MLD-integrated nitrate', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_no3_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'fe_intmld', 'MLD-integrated iron', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_fe_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'phy_intmld', 'MLD-integrated phytoplankton', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_phy_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'det_intmld', 'MLD-integrated detritus', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_det_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'pprod_gross_intmld', 'MLD-integrated gross phytoplankton production', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_pprod_gross_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'npp_intmld', 'MLD-integrated net primary productivity', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_npp_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'radbio_intmld', 'MLD-integrated photosynthetically active radiation for phytoplankton '// &
        'growth', 'h', '1', 's', 'W m-1', 'f')
    wombat%id_radbio_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    ! 100m-integrated diagnostics
    !-----------------------------------------------------------------------
    vardesc_temp = vardesc( &
        'dic_int100', '100m-integrated dissolved inorganic carbon', &
        'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_dic_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'o2_int100', '100m-integrated dissolved oxygen', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_o2_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'no3_int100', '100m-integrated nitrate', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_no3_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'fe_int100', '100m-integrated iron', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_fe_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'phy_int100', '100m-integrated phytoplankton', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_phy_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'det_int100', '100m-integrated detritus', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_det_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'pprod_gross_int100', '100m-integrated gross phytoplankton production', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_pprod_gross_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'npp_int100', '100m-integrated net primary productivity', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_npp_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'radbio_int100', '100m-integrated photosynthetically active radiation for phytoplankton '// &
        'growth', 'h', '1', 's', 'W m-1', 'f')
    wombat%id_radbio_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zeuphot', 'Depth of the euphotic zone (%1 incident light)', 'h', '1', 's', 'm', 'f')
    wombat%id_zeuphot = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
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

    ! Coefficients for N2O solubility [1] (Weiss & Price, 1980, Marine Chemistry, 8, 347-359)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('a_1_n2o', wombat%a_1_n2o, -168.2459)
    call g_tracer_add_param('a_2_n2o', wombat%a_2_n2o, 226.0894)
    call g_tracer_add_param('a_3_n2o', wombat%a_3_n2o, 93.2817)
    call g_tracer_add_param('a_4_n2o', wombat%a_4_n2o, -1.48693)
    call g_tracer_add_param('b_1_n2o', wombat%b_1_n2o, -0.060361)
    call g_tracer_add_param('b_2_n2o', wombat%b_2_n2o, 0.033765)
    call g_tracer_add_param('b_3_n2o', wombat%b_3_n2o, -0.0051862)
    
    ! Compute the Schmidt number of N2O in seawater using the
    ! formulation proposed by Wanninkof (2014) Limnology and Oceanography: Methods, 12, 351-362.
    call g_tracer_add_param('a1_n2o', wombat%a1_n2o, 2356.2)
    call g_tracer_add_param('a2_n2o', wombat%a2_n2o, -166.38)
    call g_tracer_add_param('a3_n2o', wombat%a3_n2o, 6.3952)
    call g_tracer_add_param('a4_n2o', wombat%a4_n2o, -0.13422)
    call g_tracer_add_param('a5_n2o', wombat%a5_n2o, 0.0011506)


    ! Initial H+ concentration [mol/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('htotal_in', wombat%htotal_in, 1.e-8) ! dts: default conc from COBALT

    ! Scale factor to set lower limit of htotal range [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('htotal_scale_lo', wombat%htotal_scale_lo, 0.1)

    ! Scale factor to set upper limit of htotal range [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('htotal_scale_hi', wombat%htotal_scale_hi, 100.0)

    ! Absolute minimum of dissolved inorganic carbon [mmol/m3] for co2 sys calcs
    !-----------------------------------------------------------------------
    call g_tracer_add_param('dic_min', wombat%dic_min, 1000.0)

    ! Absolute maximum of dissolved inorganic carbon [mmol/m3] for co2 sys calcs
    !-----------------------------------------------------------------------
    call g_tracer_add_param('dic_max', wombat%dic_max, 3000.0)

    ! Absolute minimum of alkalinity [mmol/m3] for co2 sys calcs
    !-----------------------------------------------------------------------
    call g_tracer_add_param('alk_min', wombat%alk_min, 1000.0)

    ! Absolute maximum of alkalinity [mmol/m3] for co2 sys calcs
    !-----------------------------------------------------------------------
    call g_tracer_add_param('alk_max', wombat%alk_max, 3000.0)

    ! Global average surface concentration of inorganic silicate [mol/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('sio2_surf', wombat%sio2_surf, 35.0e-3 / 1035.0)

    !=======================================================================
    ! NPZD parameters
    !=======================================================================
    ! dts: note the parameter units and default values are as used by WOMBAT
    ! v3 in ACCESS-OM2 and ACCESS-ESM1.5. Unit conversions are done
    ! internally to account for the different units carried in this generic
    ! version of WOMBATmid.

    ! Initial slope of P-I curve for phytoplankton [(mg Chl m-3)-1 (W m-2)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('alphabio_phy', wombat%alphabio_phy, 2.25)

    ! Autotrophy maximum growth rate parameter a for phytoplankton [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('abioa_phy', wombat%abioa_phy, 1.0/86400.0)

    ! Autotrophy maximum growth rate parameter b for phytoplankton [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bbioa_phy', wombat%bbioa_phy, 1.050)

    ! Initial slope of P-I curve for microphytoplankton [(mg Chl m-3)-1 (W m-2)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('alphabio_dia', wombat%alphabio_dia, 2.0)

    ! Initial slope of P-I curve for trichodesmium [(mg Chl m-3)-1 (W m-2)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('alphabio_tri', wombat%alphabio_tri, 1.0)

    ! Autotrodia maximum growth rate parameter a for microphytoplankton [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('abioa_dia', wombat%abioa_dia, 1.25/86400.0)

    ! Autotrodia maximum growth rate parameter b for microphytoplankton [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bbioa_dia', wombat%bbioa_dia, 1.050)

    ! Heterotrophy maximum growth rate parameter b [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bbioh', wombat%bbioh, 1.066)

    ! Phytoplankton half saturation constant for nitrogen uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phykn', wombat%phykn, 0.3)

    ! Phytoplankton half saturation constant for iron uptake [umolFe/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phykf', wombat%phykf, 0.5)

    ! Phytoplankton minimum quota of chlorophyll to carbon [mg/mg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phyminqc', wombat%phyminqc, 0.008)

    ! Phytoplankton optimal quota of chlorophyll to carbon [mg/mg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phyoptqc', wombat%phyoptqc, 0.040)

    ! Phytoplankton optimal quota of iron to carbon [mol/mol]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phyoptqf', wombat%phyoptqf, 10e-6)

    ! Phytoplankton maximum quota of iron to carbon [mol/mol]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phymaxqf', wombat%phymaxqf, 50e-6)

    ! Phytoplankton linear mortality rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phylmor', wombat%phylmor, 0.005/86400.0)

    ! Phytoplankton quadratic mortality rate constant [m3/mmolN/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phyqmor', wombat%phyqmor, 0.05/86400.0)

    ! microphytoplankton half saturation constant for nitrogen uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diakn', wombat%diakn, 3.0)

    ! microphytoplankton half saturation constant for iron uptake [umolFe/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diakf', wombat%diakf, 1.0)

    ! microphytoplankton minimum quota of chlorodiall to carbon [mg/mg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diaminqc', wombat%diaminqc, 0.004)

    ! microphytoplankton optimal quota of chlorodiall to carbon [mg/mg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diaoptqc', wombat%diaoptqc, 0.036)

    ! microphytoplankton optimal quota of iron to carbon [mol/mol]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diaoptqf', wombat%diaoptqf, 10e-6)

    ! microphytoplankton maximum quota of iron to carbon [mol/mol]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diamaxqf', wombat%diamaxqf, 50e-6)

    ! microphytoplankton linear mortality rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('dialmor', wombat%dialmor, 0.005/86400.0)

    ! microphytoplankton quadratic mortality rate constant [m3/mmolN/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diaqmor', wombat%diaqmor, 0.05/86400.0)

    ! Chlorophyll darkness growth reduction half-saturation coefficient [W/m2]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('chlkWm2', wombat%chlkWm2, 5.0)

    ! Maximum fraction of NPP that can be routed to DOC exudation by phytoplankton [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('overflow', wombat%overflow, 0.5)

    ! Trichodesmium half saturation constant for iron uptake [umolFe/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('trikf', wombat%trikf, 0.5)

    ! Trichodesmium typical chlorophyll to carbon ratio [mg/mg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('trichlc', wombat%trichlc, 0.01)

    ! Trichodesmium typical nitrogen to carbon ratio [mol/mol]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('trin2c', wombat%trin2c, 50.0/300.0)

    ! Zooplankton carbon bulk ingestion efficiency [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooCingest', wombat%zooCingest, 0.90)

    ! Zooplankton carbon assimilation efficiency [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooCassim', wombat%zooCassim, 0.30)

    ! Zooplankton iron bulk ingestion efficiency [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooFeingest', wombat%zooFeingest, 0.20)

    ! Zooplankton iron assimilation efficiency [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooFeassim', wombat%zooFeassim, 0.90)

    ! Zooplankton fraction of excretion to DOM [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooexcrdom', wombat%zooexcrdom, 0.20)

    ! Zooplankton half saturation coefficient for linear mortality
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zookz', wombat%zookz, 0.25)

    ! Zooplankton maximum grazing rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zoogmax', wombat%zoogmax, 3.0/86400.0)

    ! Zooplankton prey capture rate constant for bacteria 1 [m6/mmol2/s]
    !  - e.g., protozoans feeding on bacteria
    !  - aim for half-saturation coefficent B1/2 = 5.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepsbac1', wombat%zooepsbac1, 0.10/86400.0)

    ! Zooplankton prey capture rate constant for bacteria 2 [m6/mmol2/s]
    !  - e.g., protozoans feeding on bacteria
    !  - aim for half-saturation coefficent B1/2 = 5.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepsbac2', wombat%zooepsbac2, 0.10/86400.0)

    ! Zooplankton prey capture rate constant for ammonia oxidizing archaea [m6/mmol2/s]
    !  - e.g., ciliates feeding on ammonia oxidizing archaea (similar size or larger than pico)
    !  - aim for half-saturation coefficent B1/2 = 5.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepsaoa', wombat%zooepsaoa, 0.25/86400.0)

    ! Zooplankton prey capture rate constant for nanophytoplankton [m6/mmol2/s]
    !  - e.g., ciliates feeding on small (nano/pico)phytoplankton
    !  - aim for half-saturation coefficent B1/2 = 2.5 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepsphy', wombat%zooepsphy, 0.50/86400.0)

    ! Zooplankton prey capture rate constant for microphytoplankton [m6/mmol2/s]
    !  - e.g., larger ciliates feeding on smaller diatoms and other microphytoplankton
    !  - aim for half-saturation coefficent B1/2 = 3.5 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepsdia', wombat%zooepsdia, 0.20/86400.0)

    ! Zooplankton prey capture rate constant for small detritus [m6/mmol2/s]
    !  - e.g., protozoa grazing on slowly sinking detrital particles
    !  - aim for half-saturation coefficent B1/2 = 5.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepsdet', wombat%zooepsdet, 0.10/86400.0)

    ! Zooplankton preference for bacteria 1 [0-1]
    ! Landry (2025) J. Plankton Res. --> find that ~100 mg C m-2 day-1 of ~500 mg C m-2 d-1
    !  of microzooplankton grazing/biomass gain comes from bacterivory
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefbac1', wombat%zprefbac1, 0.1)

    ! Zooplankton preference for bacteria 2 [0-1]
    ! Landry (2025) J. Plankton Res. --> find that ~100 mg C m-2 day-1 of ~500 mg C m-2 d-1
    !  of microzooplankton grazing/biomass gain comes from bacterivory
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefbac2', wombat%zprefbac2, 0.1)

    ! Zooplankton preference for ammonia oxidizing archaea [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefaoa', wombat%zprefaoa, 0.25)

    ! Zooplankton preference for phytoplankton [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefphy', wombat%zprefphy, 0.5)

    ! Zooplankton preference for microphytoplankton [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefdia', wombat%zprefdia, 0.025)

    ! Zooplankton preference for detritus [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefdet', wombat%zprefdet, 0.025)

    ! Zooplankton respiration rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zoolmor', wombat%zoolmor, 0.001/86400.0)

    ! Zooplankton quadratic mortality rate constant [m3/mmolN/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooqmor', wombat%zooqmor, 0.05/86400.0)

    ! Mesozooplankton carbon bulk ingestion efficiency [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesCingest', wombat%mesCingest, 0.90)

    ! Mesozooplankton carbon assimilation efficiency [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesCassim', wombat%mesCassim, 0.30)

    ! Mesozooplankton iron bulk ingestion efficiency [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesFeingest', wombat%mesFeingest, 0.20)

    ! Mesozooplankton iron assimilation efficiency [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesFeassim', wombat%mesFeassim, 0.90)

    ! Mesozooplankton fraction of excretion to DOM [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesexcrdom', wombat%mesexcrdom, 0.20)

    ! Mesozooplankton half saturation coefficient for linear mortality
    !-----------------------------------------------------------------------
    call g_tracer_add_param('meskz', wombat%meskz, 0.25)

    ! Mesozooplankton maximum grazing rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesgmax', wombat%mesgmax, 1.0/86400.0)

    ! Mesozooplankton prey capture rate constant for bacteria 1 [m6/mmol2/s]
    !  - e.g., appendicularians filter feeding on bacteria
    !  - aim for half-saturation coefficent B1/2 = 3.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsbac1', wombat%mesepsbac1, 0.11/86400.0)

    ! Mesozooplankton prey capture rate constant for bacteria 2 [m6/mmol2/s]
    !  - e.g., appendicularians filter feeding on bacteria
    !  - aim for half-saturation coefficent B1/2 = 3.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsbac2', wombat%mesepsbac2, 0.11/86400.0)

    ! Mesozooplankton prey capture rate constant for ammonia oxidizing archaea [m6/mmol2/s]
    !  - e.g., appendicularians filter feeding on ammonia oxidizing archaea
    !  - aim for half-saturation coefficent B1/2 = 3.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsaoa', wombat%mesepsaoa, 0.11/86400.0)

    ! Mesozooplankton prey capture rate constant for nanophytoplankton [m6/mmol2/s]
    !  - e.g., appendicularians filter feeding on small phytoplankton
    !  - aim for half-saturation coefficent B1/2 = 3.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsphy', wombat%mesepsphy, 0.11/86400.0)

    ! Mesozooplankton prey capture rate constant for microphytoplankton [m6/mmol2/s]
    !  - e.g., copepods preying on diatoms and other microphytoplankton
    !  - aim for half-saturation coefficent B1/2 = 2.5 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsdia', wombat%mesepsdia, 0.20/86400.0)

    ! Mesozooplankton prey capture rate constant for small detritus [m6/mmol2/s]
    !  - e.g., appendicularians filter feeding on fine detritus
    !  - aim for half-saturation coefficent B1/2 = 3.5 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsdet', wombat%mesepsdet, 0.08/86400.0)

    ! Mesozooplankton prey capture rate constant for large detritus [m6/mmol2/s]
    !  - e.g., copepods consuming sinking aggregates of marine snow
    !  - aim for half-saturation coefficent B1/2 = 10.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsbdet', wombat%mesepsbdet, 0.01/86400.0)

    ! Mesozooplankton prey capture rate constant for microzooplankton [m6/mmol2/s]
    !  - e.g., chaetognaths preying on copepods; copepods consuming ciliates
    !  - aim for half-saturation coefficent B1/2 = 5.0 mmolC/m3, where B1/2 = (gmax/eps)^(0.5)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepszoo', wombat%mesepszoo, 0.04/86400.0)
    
    ! Mesozooplankton preference for bacteria 1 [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefbac1', wombat%mprefbac1, 0.05)

    ! Mesozooplankton preference for bacteria 2 [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefbac2', wombat%mprefbac2, 0.05)

    ! Mesozooplankton preference for ammonia oxidizing archaea [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefaoa', wombat%mprefaoa, 0.1)

    ! Mesozooplankton preference for phytoplankton [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefphy', wombat%mprefphy, 0.05)

    ! Mesozooplankton preference for microphytoplankton [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefdia', wombat%mprefdia, 0.3)

    ! Mesozooplankton preference for detritus [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefdet', wombat%mprefdet, 0.05)

    ! Mesozooplankton preference for large detritus (aggregates) [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefbdet', wombat%mprefbdet, 0.1)

    ! Mesozooplankton preference for zooplankton [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefzoo', wombat%mprefzoo, 0.3)

    ! Mesozooplankton respiration rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('meslmor', wombat%meslmor, 0.001/86400.0)

    ! Mesozooplankton quadratic mortality rate constant [m3/mmolN/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesqmor', wombat%mesqmor, 0.5/86400.0)

    ! Prey switching exponent for microzooplantkon
    ! when <1, more even feeding across prey items
    ! when =1, grazing proportional to prey biomasses
    ! when >1, overweighting abundant prey and downweighting scarce prey
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zoopreyswitch', wombat%zoopreyswitch, 1.4)

    ! Prey switching exponent for mesozooplantkon
    ! when <1, more even feeding across prey items
    ! when =1, grazing proportional to prey biomasses
    ! when >1, overweighting abundant prey and downweighting scarce prey
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mespreyswitch', wombat%mespreyswitch, 1.8)

    ! Detritus remineralisation rate constant [m3/mmolC/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('detlrem', wombat%detlrem, 0.5/86400.0)
    
    ! Base detritus sinking rate coefficient [m/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('wdetbio', wombat%wdetbio, 5.0/86400.0)
    
    ! Base big detritus sinking rate coefficient [m/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('wbdetbio', wombat%wbdetbio, 25.0/86400.0)
    
    ! Detritus maximum sinking rate coefficient [m/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('wdetmax', wombat%wdetmax, 50.0/86400.0)
    
    ! Bottom thickness [m]
    !-----------------------------------------------------------------------
    ! Thickness over which tracer values are integrated to define the bottom layer
    call g_tracer_add_param('bottom_thickness', wombat%bottom_thickness, 1.0)

    ! Detritus remineralisation rate constant in sediments [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('detlrem_sed', wombat%detlrem_sed, 0.005/86400.0)
    
    ! Phytoplankton biomass threshold to scale recycling [mmolC/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phybiot', wombat%phybiot, 0.5)

    ! Microphytoplankton biomass threshold to scale recycling [mmolC/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('diabiot', wombat%diabiot, 0.5)

    ! Base CaCO3 sinking rate coefficient [m/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('wcaco3', wombat%wcaco3, 5.0/86400.0)

    ! CaCO3 remineralisation rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('caco3lrem', wombat%caco3lrem, 0.01/86400.0)

    ! CaCO3 remineralization rate constant in sediments [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('caco3lrem_sed', wombat%caco3lrem_sed, 0.01/86400.0)

    ! Ceiling of omega in the sediments (controls rate of CaCO3 dissolution) [0-1]
    ! - if == 1.0, then there may be at minimum no dissolution of CaCO3
    ! - if < 1.0, then there is always some dissolution of CaCO3 when when supersaturated
    !-----------------------------------------------------------------------
    call g_tracer_add_param('omegamax_sed', wombat%omegamax_sed, 0.7)

    ! CaCO3 inorganic fraction [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('f_inorg', wombat%f_inorg, 0.04)

    ! CaCO3 dissolution factor due to calcite undersaturation
    !-----------------------------------------------------------------------
    call g_tracer_add_param('disscal', wombat%disscal, 0.250)

    ! CaCO3 dissolution factor due to aragonite undersaturation 
    !-----------------------------------------------------------------------
    call g_tracer_add_param('dissara', wombat%dissara, 0.100)

    ! CaCO3 dissolution factor due to detritus remineralisation creating anoxic microenvironment
    !-----------------------------------------------------------------------
    call g_tracer_add_param('dissdet', wombat%dissdet, 0.200)

    ! Background concentration of iron-binding ligand [umol/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('ligand', wombat%ligand, 0.5)

    ! Fraction of dissolved iron in colloidal form [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('fcolloid', wombat%fcolloid, 0.25)

    ! Precipitation of Fe` as nanoparticles (in excess of solubility) [/d]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('knano_dfe', wombat%knano_dfe, 0.1)

    ! Scavenging of Fe` onto biogenic particles [(mmolC/m3)-1 d-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('kscav_dfe', wombat%kscav_dfe, 5e-2)

    ! Coagulation of dFe onto organic particles [(mmolC/m3)-1 d-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('kcoag_dfe', wombat%kcoag_dfe, 5e-8)

    ! Factor increase in biogenic silica dissolution caused by bacterial activity [ ]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bsi_fbac', wombat%bsi_fbac, 5.0)

    ! Half-saturation coefficient modulating increase in biogenic silica dissolution due to bacterial activity [molC/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bsi_kbac', wombat%bsi_kbac, 0.5/1.035e6)

    ! Ammonia Oxidizing Archaea half saturation constant for NH4 uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoa_knh4', wombat%aoa_knh4, 0.1)

    ! Ammonia Oxidizing Archaea diffusive uptake limit for oxygen [m3/mmolC/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoa_poxy', wombat%aoa_poxy, 275.0/86400.0)

    ! Ammonia Oxidizing Archaea biomass yield per O2 [mol O2 / mol Biomass]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoa_yoxy', wombat%aoa_yoxy, 15.5)

    ! Ammonia Oxidizing Archaea biomass yield per NH4 [mol NH4 / mol Biomass]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoa_ynh4', wombat%aoa_ynh4, 11.0)

    ! Ammonia Oxidizing Archaea biomass carbon to nitrogen ratio [mol C / mol N]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoa_C2N', wombat%aoa_C2N, 5.0)

    ! Ammonia Oxidizing Archaea biomass carbon to iron ratio [mol C / mol Fe]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoa_C2Fe', wombat%aoa_C2Fe, 1.0/20e-6)

    ! Ammonia Oxidizing Archaea linear mortality rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoalmor', wombat%aoalmor, 0.005/86400.0)

    ! Ammonia Oxidizing Archaea quadratic mortality rate constant [m3/mmolN/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoaqmor', wombat%aoaqmor, 0.05/86400.0)

    ! Facultative heterotrophic bacteria #1 maximum rate of uptake of DOC [mmol/m3/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_Vmax_doc', wombat%bac1_Vmax_doc, 6.7/86400.0)

    ! Facultative heterotrophic bacteria #1 maximum rate of uptake of NO3 [mmol/m3/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_Vmax_no3', wombat%bac1_Vmax_no3, 7.2/86400.0)

    ! Facultative heterotrophic bacteria #1 diffusive uptake limit of O2 [m3/mmolC/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_poxy', wombat%bac1_poxy, 450.0/86400.0)

    ! Facultative heterotrophic bacteria #1 half saturation constant for nitrate uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_kno3', wombat%bac1_kno3, 15.0)

    ! Facultative heterotrophic bacteria #1 minimum half saturation constant for DOC uptake [mmolC/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_kdoc_min', wombat%bac1_kdoc_min, 10.0)

    ! Facultative heterotrophic bacteria #1 maximum half saturation constant for DOC uptake [mmolC/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_kdoc_max', wombat%bac1_kdoc_max, 100.0)

    ! Facultative heterotrophic bacteria #1 half saturation constant for ammonium uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_knh4', wombat%bac1_knh4, 0.1)

    ! Facultative heterotrophic bacteria #1 half saturation constant for dissolved iron uptake [umolFe/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_kfer', wombat%bac1_kfer, 0.5)

    ! Facultative heterotrophic bacteria #1 aerobic biomass yield per DOC [mol DOC / mol Biomass]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_yaerC', wombat%bac1_yaerC, 6.7)

    ! Facultative heterotrophic bacteria #1 aerobic biomass yield per O2 [mol O2/ mol Biomass]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_yoxy', wombat%bac1_yoxy, 6.3)

    ! Facultative heterotrophic bacteria #1 anaerobic biomass yield per DOC [mol DOC / mol Biomass]
    ! NO3 --> N2O
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_yanaC', wombat%bac1_yanaC, 9.9)

    ! Facultative heterotrophic bacteria #1 aerobic biomass yield per NO3 [mol NO3/ mol Biomass]
    ! NO3 --> N2O
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_yno3', wombat%bac1_yno3, 9.8)

    ! Facultative heterotrophic bacteria #1 biomass carbon to nitrogen ratio [mol C / mol N]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_C2N', wombat%bac1_C2N, 5.0)

    ! Facultative heterotrophic bacteria #1 biomass carbon to iron ratio [mol C / mol Fe]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1_C2Fe', wombat%bac1_C2Fe, 1.0/20e-6)

    ! Facultative heterotrophic bacteria #1 linear mortality rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1lmor', wombat%bac1lmor, 0.005/86400.0)

    ! Facultative heterotrophic bacteria #1 quadratic mortality rate constant [m3/mmolN/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac1qmor', wombat%bac1qmor, 0.05/86400.0)

    ! Facultative heterotrophic bacteria #2 maximum rate of uptake of DOC [mmol/m3/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2_Vmax_doc', wombat%bac2_Vmax_doc, 6.7/86400.0)

    ! Facultative heterotrophic bacteria #2 diffusive uptake limit of O2 [m3/mmolC/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2_poxy', wombat%bac2_poxy, 450.0/86400.0)

    ! Facultative heterotrophic bacteria #2 diffusive uptake limit of N2O [m3/mmolC/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2_pn2o', wombat%bac2_pn2o, 452.0/86400.0)

    ! Facultative heterotrophic bacteria #2 minimum half saturation constant for DOC uptake [mmolC/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2_kdoc_min', wombat%bac2_kdoc_min, 10.0)

    ! Facultative heterotrophic bacteria #2 maximum half saturation constant for DOC uptake [mmolC/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2_kdoc_max', wombat%bac2_kdoc_max, 100.0)

    ! Facultative heterotrophic bacteria #2 half saturation constant for ammonium uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2_knh4', wombat%bac2_knh4, 0.1)

    ! Facultative heterotrophic bacteria #2 half saturation constant for dissolved iron uptake [umolFe/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2_kfer', wombat%bac2_kfer, 0.5)

    ! Facultative heterotrophic bacteria #2 aerobic biomass yield per DOC [mol DOC / mol Biomass]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2_yaerC', wombat%bac2_yaerC, 6.7)

    ! Facultative heterotrophic bacteria #2 aerobic biomass yield per O2 [mol O2/ mol Biomass]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2_yoxy', wombat%bac2_yoxy, 6.3)

    ! Facultative heterotrophic bacteria #2 aerobic biomass yield per N2O [mol N2/ mol Biomass]
    ! N2O --> N2
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2_yn2o', wombat%bac2_yn2o, 25.3/2.0)

    ! Facultative heterotrophic bacteria #2 anaerobic biomass yield per DOC [mol DOC / mol Biomass]
    ! N2O --> N2   |   must be >= `bac_yaerC` to ensure that N2O actually accumulates in the ocean 
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2_yanaC', wombat%bac2_yanaC, 6.701)

    ! Facultative heterotrophic bacteria #2 biomass carbon to nitrogen ratio [mol C / mol N]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2_C2N', wombat%bac2_C2N, 5.0)

    ! Facultative heterotrophic bacteria #2 biomass carbon to iron ratio [mol C / mol Fe]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2_C2Fe', wombat%bac2_C2Fe, 1.0/20e-6)

    ! Facultative heterotrophic bacteria #2 linear mortality rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2lmor', wombat%bac2lmor, 0.005/86400.0)

    ! Facultative heterotrophic bacteria #2 quadratic mortality rate constant [m3/mmolN/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bac2qmor', wombat%bac2qmor, 0.05/86400.0)

    ! Anammox bacteria half saturation constant for ammonium uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoxkn', wombat%aoxkn, 0.5)

    ! Anammox bacteria maximum growth * biomass rate [/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoxmumax', wombat%aoxmumax, 0.0025/86400.0)

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
    type(g_tracer_type), pointer :: tracer_list

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

    ! Detritus sinking velocity [m/s]
    !-----------------------------------------------------------------------
    ! Default value matches Ziehn et al 2020 but differs from Hayashida et
    ! al 2020
    call g_tracer_add_param('wdetbio', wombat%wdetbio, 5.0/86400.0)
    call g_tracer_add_param('wbdetbio', wombat%wbdetbio, 25.0/86400.0)
    call g_tracer_add_param('wcaco3', wombat%wcaco3, 5.0/86400.0) ! Based on 10µm average size

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
        flux_bottom = .true., &
        flux_virtual = .true.)

    ! Silicic acid (H4SiO4)
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Ammonium
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
        flux_gas_param = (/ as_coeff_wombatmid, 9.7561e-06 /), & ! dts: param(2) converts Pa -> atm
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

    ! Detritus
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'det', &
        longname = 'Detritus', &
        units = 'mol/kg', &
        prog = .true., &
        move_vertical = .true., &
        btm_reservoir = .true.)

    ! Detrital iron content
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'detfe', &
        longname = 'Detrital iron content', &
        units = 'mol/kg', &
        prog = .true., &
        move_vertical = .true., &
        btm_reservoir = .true.)

    ! Big detritus
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'bdet', &
        longname = 'Big detritus', &
        units = 'mol/kg', &
        prog = .true., &
        move_vertical = .true., &
        btm_reservoir = .true.)

    ! Big detrital iron content
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'bdetfe', &
        longname = 'Big detrital iron content', &
        units = 'mol/kg', &
        prog = .true., &
        move_vertical = .true., &
        btm_reservoir = .true.)

    ! Big detrital silicon content
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'bdetsi', &
        longname = 'Big detrital silicon content', &
        units = 'mol/kg', &
        prog = .true., &
        move_vertical = .true., &
        btm_reservoir = .true.)

    ! Dissolved organic carbon
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'doc', &
        longname = 'Dissolved organic carbon', &
        units = 'mol/kg', &
        prog = .true.)
    
    ! Dissolved organic nitrogen
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'don', &
        longname = 'Dissolved organic nitrogen', &
        units = 'mol/kg', &
        prog = .true.)
    
    ! Facultative heterotrophic bacteria #1
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'bac1', &
        longname = 'Facultative heterotrophic bacteria #1', &
        units = 'mol/kg', &
        prog = .true.)
    
    ! Facultative heterotrophic bacteria #2
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'bac2', &
        longname = 'Facultative heterotrophic bacteria #2', &
        units = 'mol/kg', &
        prog = .true.)
    
    ! Ammonia oxidizing archaea
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'aoa', &
        longname = 'Ammonia oxidizing archaea', &
        units = 'mol/kg', &
        prog = .true.)
    
    ! Nitrous oxide (N2O)
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'n2o', &
        longname = 'Nitrous oxide', &
        units = 'mol/kg', &
        prog = .true., &
        flux_gas = .true., &
        flux_bottom = .false., &
        flux_gas_name = 'n2o_flux', &
        flux_gas_type = 'air_sea_gas_flux_generic', &
        flux_gas_molwt = WTMCO2, & !pjb: N2O molar mass is 44.01 g/mol, same as CO2
        flux_gas_param = (/ as_coeff_wombatmid, 9.7561e-06 /), & ! dts: param(2) converts Pa -> atm
        flux_gas_restart_file = 'ocean_wombatmid_airsea_flux.res.nc', &
        flux_virtual = .false.)
    
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
        flux_gas_param = (/ as_coeff_wombatmid, 9.7561e-06 /), & ! dts: param(2) converts Pa -> atm
        flux_gas_restart_file = 'ocean_wombatmid_airsea_flux.res.nc', &
        flux_virtual = .true.)

    ! DICp (preformed Dissolved inorganic carbon)
    ! dts: Note, we use flux_virtual=.true. only to ensure that an stf array is allocated for dicp.
    ! The dicp stf is set to equal the dic stf in update_from_coupler.
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'dicp', &
        longname = 'preformed Dissolved Inorganic Carbon', &
        units = 'mol/kg', &
        prog = .true., &
        flux_virtual = .true.)

    ! DICr (remineralised dissolved inorganic carbon)
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'dicr', &
        longname = 'remineralised Dissolved Inorganic Carbon', &
        units = 'mol/kg', &
        prog = .true., &
        flux_bottom = .true.)

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
        flux_param  = (/ 1.0 /), & ! dts: conversion to mol/m2/s done in data_table
        flux_bottom = .true.)

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

    ! pjb: included here so included in restart
    call g_tracer_add(tracer_list, package_name, &
        name = 'detbury', &
        longname = 'Amount of detritus lost to burial', &
        units = 'mol m-2 s-1', &
        prog = .false.)

    ! pjb: included here so included in restart
    call g_tracer_add(tracer_list, package_name, &
        name = 'caco3bury', &
        longname = 'Amount of CaCO3 lost to burial', &
        units = 'mol m-2 s-1', &
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
    type(g_tracer_type), pointer           :: tracer_list
    integer, intent(in)                    :: ilb, jlb
    real, dimension(ilb:,jlb:), intent(in) :: salt_flux_added

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

    ! Set dicp stf equal to dic stf
    call g_tracer_set_values(tracer_list, 'dicp', 'stf', wombat%p_dic_stf, isd, jsd)

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
    type(g_tracer_type), pointer :: tracer_list
    real, intent(in)             :: dt
    integer, intent(in)          :: tau
    type(time_type), intent(in)  :: model_time

    integer                         :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, i, j
    real, dimension(:,:,:), pointer :: grid_tmask
    real                            :: orgflux
    logical                         :: used

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, &
        grid_tmask=grid_tmask)

    ! Move bottom reservoirs to sediment tracers
    !-----------------------------------------------------------------------
    call g_tracer_get_values(tracer_list, 'det', 'btm_reservoir', wombat%det_btm, isd, jsd)
    call g_tracer_get_values(tracer_list, 'detfe', 'btm_reservoir', wombat%detfe_btm, isd, jsd)
    call g_tracer_get_values(tracer_list, 'bdet', 'btm_reservoir', wombat%bdet_btm, isd, jsd)
    call g_tracer_get_values(tracer_list, 'bdetfe', 'btm_reservoir', wombat%bdetfe_btm, isd, jsd)
    call g_tracer_get_values(tracer_list, 'bdetsi', 'btm_reservoir', wombat%bdetsi_btm, isd, jsd)
    call g_tracer_get_values(tracer_list, 'caco3', 'btm_reservoir', wombat%caco3_btm, isd, jsd)

    ! Calculate burial of deposited detritus (Dunne et al., 2007)
    wombat%fbury(:,:) = 0.0
    if (do_burial) then
      do i = isc, iec
        do j = jsc, jec
          orgflux = (wombat%det_btm(i,j) + wombat%bdet_btm(i,j)) / dt * 86400 * 1e3 ! mmol C m-2 day-1
          wombat%fbury(i,j) = 0.013 + 0.53 * orgflux**2.0 / (7.0 + orgflux)**2.0  ! Eq. 3 Dunne et al. 2007
        enddo
      enddo
    endif

    call g_tracer_get_pointer(tracer_list, 'detbury', 'field', wombat%p_detbury)
    call g_tracer_get_pointer(tracer_list, 'caco3bury', 'field', wombat%p_caco3bury)
    if (do_conserve_tracers) then
      wombat%p_detbury(:,:,1) = (wombat%det_btm(:,:) + wombat%bdet_btm(:,:)) / dt * wombat%fbury(:,:)
      wombat%p_caco3bury(:,:,1) = wombat%caco3_btm(:,:) / dt * wombat%fbury(:,:)
    else
      wombat%p_detbury(:,:,1) = 0.0 
      wombat%p_caco3bury(:,:,1) = 0.0
    endif

    call g_tracer_get_pointer(tracer_list, 'det_sediment', 'field', wombat%p_det_sediment)
    wombat%p_det_sediment(:,:,1) = wombat%p_det_sediment(:,:,1) + (wombat%det_btm(:,:) + wombat%bdet_btm(:,:)) * (1.0-wombat%fbury(:,:)) ! [mol/m2]
    call g_tracer_set_values(tracer_list, 'det', 'btm_reservoir', 0.0)
    call g_tracer_set_values(tracer_list, 'bdet', 'btm_reservoir', 0.0)

    call g_tracer_get_pointer(tracer_list, 'detfe_sediment', 'field', wombat%p_detfe_sediment)
    wombat%p_detfe_sediment(:,:,1) = wombat%p_detfe_sediment(:,:,1) + (wombat%detfe_btm(:,:) + wombat%bdetfe_btm(:,:)) * (1.0-wombat%fbury(:,:)) ! [mol/m2]
    call g_tracer_set_values(tracer_list, 'detfe', 'btm_reservoir', 0.0)
    call g_tracer_set_values(tracer_list, 'bdetfe', 'btm_reservoir', 0.0)

    call g_tracer_get_pointer(tracer_list, 'detsi_sediment', 'field', wombat%p_detsi_sediment)
    wombat%p_detsi_sediment(:,:,1) = wombat%p_detsi_sediment(:,:,1) + wombat%bdetsi_btm(:,:) * (1.0-wombat%fbury(:,:)) ! [mol/m2]
    call g_tracer_set_values(tracer_list, 'bdetsi', 'btm_reservoir', 0.0)

    call g_tracer_get_pointer(tracer_list, 'caco3_sediment', 'field', wombat%p_caco3_sediment)
    wombat%p_caco3_sediment(:,:,1) =  wombat%p_caco3_sediment(:,:,1) + wombat%caco3_btm(:,:) * (1.0-wombat%fbury(:,:)) ! [mol/m2]
    call g_tracer_set_values(tracer_list, 'caco3', 'btm_reservoir', 0.0)

    ! Send diagnostics
    !-----------------------------------------------------------------------
    if (wombat%id_det_sed_depst .gt. 0) &
      used = g_send_data(wombat%id_det_sed_depst, (wombat%det_btm + wombat%bdet_btm) / dt, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_detfe_sed_depst .gt. 0) &
      used = g_send_data(wombat%id_detfe_sed_depst, (wombat%detfe_btm + wombat%bdetfe_btm) / dt, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_detsi_sed_depst .gt. 0) &
      used = g_send_data(wombat%id_detsi_sed_depst, wombat%bdetsi_btm / dt, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_caco3_sed_depst .gt. 0) &
      used = g_send_data(wombat%id_caco3_sed_depst, wombat%caco3_btm / dt, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_fbury .gt. 0) &
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
    type(g_tracer_type), pointer               :: tracer_list
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
    real                                    :: u_npz, g_npz, m_npz, Xzoo, Xmes
    real                                    :: biono3, bion2o, bionh4, biooxy, biofer, biodoc, biodon, biocaco3
    real                                    :: biophy, biodia, biozoo, biomes, biodet, biobdet, biobac1, biobac2, bioaoa
    real                                    :: biophyfe, biodiafe, biozoofe, biomesfe, zooprey, mesprey
    real                                    :: denom, wzbac1, wzbac2, wzaoa, wzphy, wzdia, wzdet, wzbdet, wzzoo, wzmes, wzsum
    real                                    :: fbc
    real                                    :: no3_bgc_change, caco3_bgc_change
    real                                    :: epsi = 1.0e-30
    real                                    :: pi = 3.14159265358979
    real                                    :: Rgas = 8.314462168 ! J/(K mol)
    integer                                 :: ichl, iter, max_iter
    real                                    :: par_phy_mldsum, par_z_mldsum
    real                                    :: chl, ndet, carb, zchl, zval, sqrt_zval, phy_chlc, dia_chlc, phi
    real                                    :: phy_limnh4, phy_limno3, phy_limdin
    real                                    :: dia_limnh4, dia_limno3, dia_limdin
    real                                    :: phy_pisl, phy_pisl2 
    real                                    :: pchl_pisl, pchl_mumin, pchl_muopt
    real                                    :: dia_pisl, dia_pisl2 
    real                                    :: dchl_pisl, dchl_mumin, dchl_muopt
    real                                    :: zooegesbac1fe, zooegesbac2fe, zooegesaoafe, zooegesphyfe, zooegesdiafe, zooegesdetfe
    real                                    :: zooassibac1fe, zooassibac2fe, zooassiaoafe, zooassiphyfe, zooassidiafe, zooassidetfe
    real                                    :: zooexcrbac1fe, zooexcrbac2fe, zooexcraoafe, zooexcrphyfe, zooexcrdiafe, zooexcrdetfe
    real                                    :: mesegesbac1fe, mesegesbac2fe, mesegesaoafe, mesegesphyfe, mesegesdiafe, mesegesdetfe, mesegesbdetfe, mesegeszoofe 
    real                                    :: mesassibac1fe, mesassibac2fe, mesassiaoafe, mesassiphyfe, mesassidiafe, mesassidetfe, mesassibdetfe, mesassizoofe
    real                                    :: mesexcrbac1fe, mesexcrbac2fe, mesexcraoafe, mesexcrphyfe, mesexcrdiafe, mesexcrdetfe, mesexcrbdetfe, mesexcrzoofe
    real                                    :: zooexcrbac1n, zooexcrbac2n, mesexcrbac1n, mesexcrbac2n,zooexcraoan, mesexcraoan
    real, dimension(:,:), allocatable       :: ek_bgr, par_bgr_mid, par_bgr_top
    real, dimension(:), allocatable         :: wsink1, wsink2, wsinkcal
    real                                    :: max_wsink
    real, dimension(4,61)                   :: zbgr
    real, dimension(3)                      :: dbgr, cbgr
    real                                    :: ztemk, I_ztemk, fe_keq, fe_par, fe_sfe, fe_tfe, partic
    real                                    :: fesol1, fesol2, fesol3, fesol4, fesol5, hp, fe3sol
    real                                    :: biof, zno3, zfermin
    real                                    :: phy_Fe2C, dia_Fe2C, zoo_Fe2C, mes_Fe2C, det_Fe2C, bdet_Fe2C, dom_N2C, dia_Si2C, bdet_Si2C
    real                                    :: phy_minqfe, phy_maxqfe
    real                                    :: dia_minqfe, dia_maxqfe
    real                                    :: zoo_slmor, mes_slmor, epsmin
    real                                    :: hco3, diss_cal, diss_ara, diss_det
    real                                    :: avedetbury, avecaco3bury
    real                                    :: dzt_bot, dzt_bot_os
    real                                    :: bac_Vdoc, bac_Voxy, bac_Vno3, bac_Vn2o, bac_muana, bac_muaer, bac_limnh4
    real                                    :: aoa_Vnh4, aoa_Voxy
    real                                    :: K_am_silica, gamma0, alphaH2O, deltaV0, spmvcorrect, disssi_temp, disssi_usat, disssi_bact
    real, dimension(:,:,:,:), allocatable   :: n_pools, c_pools, si_pools
    logical                                 :: used, converged

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
    
    call g_tracer_get_values(tracer_list, 'dic', 'field', wombat%f_dic, isd, jsd, ntau=tau, &
        positive=.true.)
    call g_tracer_get_values(tracer_list, 'no3', 'field', wombat%f_no3, isd, jsd, ntau=tau, &
        positive=.true.)
    call g_tracer_get_values(tracer_list, 'alk', 'field', wombat%f_alk, isd, jsd, ntau=tau, &
        positive=.true.)
 
    do k = 1,nk !{ 
     do j = jsc,jec; do i = isc,iec
       wombat%htotallo(i,j) = wombat%htotal_scale_lo * wombat%htotal(i,j,k)
       wombat%htotalhi(i,j) = wombat%htotal_scale_hi * wombat%htotal(i,j,k)
     enddo; enddo 

     if (k.eq.1) then !{
       call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,k), &
           Temp(:,:,k), Salt(:,:,k), &
           min(wombat%dic_max*mmol_m3_to_mol_kg, max(wombat%f_dic(:,:,k), wombat%dic_min*mmol_m3_to_mol_kg)), &
           max(wombat%f_no3(:,:,k) / 16., 1e-9), &
           wombat%sio2(:,:), &
           min(wombat%alk_max*mmol_m3_to_mol_kg, max(wombat%f_alk(:,:,k), wombat%alk_min*mmol_m3_to_mol_kg)), &
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
           Temp(:,:,k), Salt(:,:,k), &
           min(wombat%dic_max*mmol_m3_to_mol_kg, max(wombat%f_dic(:,:,k), wombat%dic_min*mmol_m3_to_mol_kg)), &
           max(wombat%f_no3(:,:,k) / 16., 1e-9), &
           wombat%sio2(:,:), &
           min(wombat%alk_max*mmol_m3_to_mol_kg, max(wombat%f_alk(:,:,k), wombat%alk_min*mmol_m3_to_mol_kg)), &
           wombat%htotallo(:,:), wombat%htotalhi(:,:), &
           wombat%htotal(:,:,k), &
           co2_calc=trim(co2_calc), zt=wombat%zw(:,:,k), &
           co2star=wombat%co2_star(:,:,k), co3_ion=wombat%co3(:,:,k), &
           omega_arag=wombat%omega_ara(:,:,k), omega_calc=wombat%omega_cal(:,:,k))

     endif !} if k.eq.1
    enddo !} do k = 1,nk

    !=======================================================================
    ! Calculate the source terms
    !=======================================================================

    wombat%dic_correct(:,:,:) = 0.0
    wombat%alk_correct(:,:,:) = 0.0
    wombat%pprod_gross(:,:,:) = 0.0
    wombat%zprod_gross(:,:,:) = 0.0
    wombat%radbio(:,:,:) = 0.0
    wombat%radmid(:,:,:) = 0.0
    wombat%radmld(:,:,:) = 0.0
    wombat%npp3d(:,:,:) = 0.0
    wombat%phy_mumax(:,:,:) = 0.0
    wombat%phy_mu(:,:,:) = 0.0
    wombat%pchl_mu(:,:,:) = 0.0
    wombat%phy_kni(:,:,:) = 0.0
    wombat%phy_kfe(:,:,:) = 0.0
    wombat%phy_lpar(:,:,:) = 0.0
    wombat%pchl_lpar(:,:,:) = 0.0
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
    wombat%dia_lpar(:,:,:) = 0.0
    wombat%dchl_lpar(:,:,:) = 0.0
    wombat%dia_lnit(:,:,:) = 0.0
    wombat%dia_lnh4(:,:,:) = 0.0
    wombat%dia_lno3(:,:,:) = 0.0
    wombat%dia_lfer(:,:,:) = 0.0
    wombat%dia_dfeupt(:,:,:) = 0.0
    wombat%tri_lfer(:,:,:) = 0.0
    wombat%tri_lpar(:,:,:) = 0.0
    wombat%trimumax(:,:,:) = 0.0
    wombat%sileqc(:,:,:) = 0.0
    wombat%disssi(:,:,:) = 0.0
    wombat%bsidiss(:,:,:) = 0.0
    wombat%feIII(:,:,:) = 0.0
    wombat%felig(:,:,:) = 0.0
    wombat%fecol(:,:,:) = 0.0
    wombat%feprecip(:,:,:) = 0.0
    wombat%fescaven(:,:,:) = 0.0
    wombat%fescadet(:,:,:) = 0.0
    wombat%fescabdet(:,:,:) = 0.0
    wombat%fesources(:,:,:) = 0.0
    wombat%fesinks(:,:,:) = 0.0
    wombat%fecoag2det(:,:,:) = 0.0
    wombat%fecoag2bdet(:,:,:) = 0.0
    wombat%phy_feupreg(:,:,:) = 0.0
    wombat%phy_fedoreg(:,:,:) = 0.0
    wombat%phygrow(:,:,:) = 0.0
    wombat%phydoc(:,:,:) = 0.0
    wombat%phyresp(:,:,:) = 0.0
    wombat%phymort(:,:,:) = 0.0
    wombat%dia_feupreg(:,:,:) = 0.0
    wombat%dia_fedoreg(:,:,:) = 0.0
    wombat%diagrow(:,:,:) = 0.0
    wombat%diadoc(:,:,:) = 0.0
    wombat%diaresp(:,:,:) = 0.0
    wombat%diamort(:,:,:) = 0.0
    wombat%zooeps(:,:,:) = 0.0
    wombat%zooprefbac1(:,:,:) = 0.0
    wombat%zooprefbac2(:,:,:) = 0.0
    wombat%zooprefaoa(:,:,:) = 0.0
    wombat%zooprefphy(:,:,:) = 0.0
    wombat%zooprefdia(:,:,:) = 0.0
    wombat%zooprefdet(:,:,:) = 0.0
    wombat%zoograzbac1(:,:,:) = 0.0
    wombat%zoograzbac2(:,:,:) = 0.0
    wombat%zoograzaoa(:,:,:) = 0.0
    wombat%zoograzphy(:,:,:) = 0.0
    wombat%zoograzdia(:,:,:) = 0.0
    wombat%zoograzdet(:,:,:) = 0.0
    wombat%zooresp(:,:,:) = 0.0
    wombat%zoomort(:,:,:) = 0.0
    wombat%zooexcrbac1(:,:,:) = 0.0
    wombat%zooexcrbac2(:,:,:) = 0.0
    wombat%zooexcraoa(:,:,:) = 0.0
    wombat%zooexcrphy(:,:,:) = 0.0
    wombat%zooexcrdia(:,:,:) = 0.0
    wombat%zooexcrdet(:,:,:) = 0.0
    wombat%zooegesbac1(:,:,:) = 0.0
    wombat%zooegesbac2(:,:,:) = 0.0
    wombat%zooegesaoa(:,:,:) = 0.0
    wombat%zooegesphy(:,:,:) = 0.0
    wombat%zooegesdia(:,:,:) = 0.0
    wombat%zooegesdet(:,:,:) = 0.0
    wombat%meseps(:,:,:) = 0.0
    wombat%mesprefbac1(:,:,:) = 0.0
    wombat%mesprefbac2(:,:,:) = 0.0
    wombat%mesprefaoa(:,:,:) = 0.0
    wombat%mesprefphy(:,:,:) = 0.0
    wombat%mesprefdia(:,:,:) = 0.0
    wombat%mesprefdet(:,:,:) = 0.0
    wombat%mesprefbdet(:,:,:) = 0.0
    wombat%mesprefzoo(:,:,:) = 0.0
    wombat%mesgrazbac1(:,:,:) = 0.0
    wombat%mesgrazbac2(:,:,:) = 0.0
    wombat%mesgrazaoa(:,:,:) = 0.0
    wombat%mesgrazphy(:,:,:) = 0.0
    wombat%mesgrazdia(:,:,:) = 0.0
    wombat%mesgrazdet(:,:,:) = 0.0
    wombat%mesgrazbdet(:,:,:) = 0.0
    wombat%mesgrazzoo(:,:,:) = 0.0
    wombat%mesresp(:,:,:) = 0.0
    wombat%mesmort(:,:,:) = 0.0
    wombat%mesexcrbac1(:,:,:) = 0.0
    wombat%mesexcrbac2(:,:,:) = 0.0
    wombat%mesexcraoa(:,:,:) = 0.0
    wombat%mesexcrphy(:,:,:) = 0.0
    wombat%mesexcrdia(:,:,:) = 0.0
    wombat%mesexcrdet(:,:,:) = 0.0
    wombat%mesexcrbdet(:,:,:) = 0.0
    wombat%mesexcrzoo(:,:,:) = 0.0
    wombat%mesegesbac1(:,:,:) = 0.0
    wombat%mesegesbac2(:,:,:) = 0.0
    wombat%mesegesaoa(:,:,:) = 0.0
    wombat%mesegesphy(:,:,:) = 0.0
    wombat%mesegesdia(:,:,:) = 0.0
    wombat%mesegesdet(:,:,:) = 0.0
    wombat%mesegesbdet(:,:,:) = 0.0
    wombat%mesegeszoo(:,:,:) = 0.0
    wombat%reminr(:,:,:) = 0.0
    wombat%doc1remi(:,:,:) = 0.0
    wombat%don1remi(:,:,:) = 0.0
    wombat%doc2remi(:,:,:) = 0.0
    wombat%don2remi(:,:,:) = 0.0
    wombat%detremi(:,:,:) = 0.0
    wombat%bdetremi(:,:,:) = 0.0
    wombat%pic2poc(:,:,:) = 0.0
    wombat%dissrat(:,:,:) = 0.0
    wombat%caldiss(:,:,:) = 0.0
    wombat%aoa_loxy(:,:,:) = 0.0
    wombat%aoa_lnh4(:,:,:) = 0.0
    wombat%aoa_yn2o(:,:,:) = 0.0
    wombat%aoa_mumax(:,:,:) = 0.0
    wombat%aoa_mu(:,:,:) = 0.0
    wombat%aoagrow(:,:,:) = 0.0
    wombat%aoaresp(:,:,:) = 0.0
    wombat%aoamor1(:,:,:) = 0.0
    wombat%aoamor2(:,:,:) = 0.0
    wombat%bac1grow(:,:,:) = 0.0
    wombat%bac1resp(:,:,:) = 0.0
    wombat%bac1unh4(:,:,:) = 0.0
    wombat%bac1ufer(:,:,:) = 0.0
    wombat%bac1_lnit(:,:,:) = 1.0
    wombat%bac1_lfer(:,:,:) = 1.0
    wombat%bac1_mu(:,:,:) = 0.0
    wombat%bac1_kdoc(:,:,:) = 10.0
    wombat%bac1_fanaer(:,:,:) = 0.0
    wombat%bac1mor1(:,:,:) = 0.0
    wombat%bac1mor2(:,:,:) = 0.0
    wombat%bac1deni(:,:,:) = 0.0
    wombat%bac2grow(:,:,:) = 0.0
    wombat%bac2resp(:,:,:) = 0.0
    wombat%bac2unh4(:,:,:) = 0.0
    wombat%bac2ufer(:,:,:) = 0.0
    wombat%bac2_lnit(:,:,:) = 1.0
    wombat%bac2_lfer(:,:,:) = 1.0
    wombat%bac2_mu(:,:,:) = 0.0
    wombat%bac2_kdoc(:,:,:) = 10.0
    wombat%bac2_fanaer(:,:,:) = 0.0
    wombat%bac2mor1(:,:,:) = 0.0
    wombat%bac2mor2(:,:,:) = 0.0
    wombat%bac2deni(:,:,:) = 0.0
    wombat%aox_lnh4(:,:,:) = 0.0
    wombat%aox_mu(:,:,:) = 0.0
    wombat%nitrfix(:,:,:) = 0.0
    wombat%ammox(:,:,:) = 0.0
    wombat%anammox(:,:,:) = 0.0
    wombat%export_prod(:,:) = 0.0
    wombat%export_inorg(:,:) = 0.0
    wombat%dic_intmld(:,:) = 0.0
    wombat%o2_intmld(:,:) = 0.0
    wombat%no3_intmld(:,:) = 0.0
    wombat%fe_intmld(:,:) = 0.0
    wombat%phy_intmld(:,:) = 0.0
    wombat%det_intmld(:,:) = 0.0
    wombat%pprod_gross_intmld(:,:) = 0.0
    wombat%npp_intmld(:,:) = 0.0
    wombat%radbio_intmld(:,:) = 0.0
    wombat%dic_int100(:,:) = 0.0
    wombat%o2_int100(:,:) = 0.0
    wombat%no3_int100(:,:) = 0.0
    wombat%fe_int100(:,:) = 0.0
    wombat%phy_int100(:,:) = 0.0
    wombat%det_int100(:,:) = 0.0
    wombat%pprod_gross_int100(:,:) = 0.0
    wombat%npp_int100(:,:) = 0.0
    wombat%radbio_int100(:,:) = 0.0
    wombat%zeuphot(:,:) = 0.0
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

    ! Allocate and initialise some multi-dimensional variables
    allocate(wsink1(nk)); wsink1(:)=0.0
    allocate(wsink2(nk)); wsink2(:)=0.0
    allocate(wsinkcal(nk)); wsinkcal(:)=0.0
    allocate(ek_bgr(nk,3)); ek_bgr(:,:)=0.0
    allocate(par_bgr_mid(nk,3)); par_bgr_mid(:,:)=0.0
    allocate(par_bgr_top(nk,3)); par_bgr_top(:,:)=0.0
    allocate(n_pools(isc:iec,jsc:jec,nk,2)); n_pools(:,:,:,:)=0.0
    allocate(c_pools(isc:iec,jsc:jec,nk,2)); c_pools(:,:,:,:)=0.0
    allocate(si_pools(isc:iec,jsc:jec,nk,2)); si_pools(:,:,:,:)=0.0

    ! Set the maximum index for euphotic depth
    ! dts: in WOMBAT v3, kmeuph and k100 are integers but here they are arrays since zw
    ! may vary spatially
    allocate(kmeuph(isc:iec, jsc:jec)); kmeuph(:,:)=1
    allocate(k100(isc:iec, jsc:jec)); k100(:,:)=1
    do j = jsc,jec; do i = isc,iec;
      nz = grid_kmt(i,j)
      do k = 1,nz
        if (wombat%zw(i,j,k) .le. 400) kmeuph(i,j)=k
        if (wombat%zw(i,j,k) .le. 100) k100(i,j)=k
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
    call g_tracer_get_values(tracer_list, 'no3', 'field', wombat%f_no3, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'nh4', 'field', wombat%f_nh4, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'sil', 'field', wombat%f_sil, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'phy', 'field', wombat%f_phy, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'pchl', 'field', wombat%f_pchl, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'phyfe', 'field', wombat%f_phyfe, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'dia', 'field', wombat%f_dia, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'dchl', 'field', wombat%f_dchl, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'diafe', 'field', wombat%f_diafe, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'diasi', 'field', wombat%f_diasi, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'zoo', 'field', wombat%f_zoo, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'zoofe', 'field', wombat%f_zoofe, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'mes', 'field', wombat%f_mes, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'mesfe', 'field', wombat%f_mesfe, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'det', 'field', wombat%f_det, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'detfe', 'field', wombat%f_detfe, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'bdet', 'field', wombat%f_bdet, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'bdetfe', 'field', wombat%f_bdetfe, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'bdetsi', 'field', wombat%f_bdetsi, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'doc', 'field', wombat%f_doc, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'don', 'field', wombat%f_don, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'bac1', 'field', wombat%f_bac1, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'bac2', 'field', wombat%f_bac2, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'aoa', 'field', wombat%f_aoa, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'n2o', 'field', wombat%f_n2o, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'o2', 'field', wombat%f_o2, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'caco3', 'field', wombat%f_caco3, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'fe', 'field', wombat%f_fe, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'dic', 'field', wombat%f_dic, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'dicr', 'field', wombat%f_dicr, isd, jsd, ntau=tau) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'alk', 'field', wombat%f_alk, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
 
    
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
    !    3.  Temperature-dependence of heterotrophy                         !
    !    4.  Light limitation of phytoplankton                              !
    !    5.  Growth of chlorophyll                                          !
    !    6.  Phytoplankton uptake of iron                                   !
    !    7.  Iron chemistry                                                 !
    !    8.  Silicic acid chemistry                                         !
    !    9.  Mortality scalings and grazing                                 !
    !    10. CaCO3 calculations                                             !
    !    11. Implicit nitrogen fixation                                     !
    !    12. Facultative heterotrophy calculations                          !
    !    13. Chemoautotroph calculations                                    !
    !    14. Sources and sinks                                              !
    !    15. Tracer tendencies                                              !
    !    16. Check for conservation by ecosystem component                  !
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
        if (max_wavelength_band(n) .lt. 710) then
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
        chl = (wombat%f_pchl(i,j,k) + wombat%f_dchl(i,j,k)) * 12.0 / mmol_m3_to_mol_kg 
        ! detritus concentration conversion from mol/kg --> mgN/m3 for look-up table
        ndet = (wombat%f_det(i,j,k) + wombat%f_bdet(i,j,k)) * 16.0/122.0 * 14.0 / mmol_m3_to_mol_kg
        ! CaCO3 concentration conversion from mol/kg --> kg/m3 for look-up table
        carb = wombat%f_caco3(i,j,k) / mmol_m3_to_mol_kg * 100.09 * 1e-3 * 1e-3 ! convert to kg/m3

        ! Attenuation coefficients given chlorophyll concentration
        zchl = max(0.05, min(10.0, chl) )
        ichl = nint( 41 + 20.0*log10(zchl) + epsi )
        ek_bgr(k,1) = (zbgr(2,ichl) + ndet * dbgr(1) + carb * cbgr(1)) * dzt(i,j,k) ! [/m * m]
        ek_bgr(k,2) = (zbgr(3,ichl) + ndet * dbgr(2) + carb * cbgr(2)) * dzt(i,j,k) ! [/m * m]
        ek_bgr(k,3) = (zbgr(4,ichl) + ndet * dbgr(3) + carb * cbgr(3)) * dzt(i,j,k) ! [/m * m]

        ! BGR light available in the water column
        if (swpar.gt.0.0) then
          if (k.eq.1) then
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

        if (swpar.gt.0.0) then
          ! Euphotic zone
          if (wombat%radmid(i,j,k).gt.(swpar*0.01) .and. wombat%radmid(i,j,k).gt.0.01) then
            wombat%zeuphot(i,j) = wombat%zw(i,j,k)
          endif
          ! Light attenuation mean over the grid cells (Eq. 19 of Baird et al., 2020 GMD)
          if (k.lt.grid_kmt(i,j)) then
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
        if (wombat%zw(i,j,k).le.hblt_depth(i,j)) then
          par_phy_mldsum = par_phy_mldsum + wombat%radbio(i,j,k) * dzt(i,j,k)
          par_z_mldsum = par_z_mldsum + dzt(i,j,k)
        endif

      enddo !}

      !--- Aggregate light in mixed layer and calculate maximum growth rates ---!
      do k = 1,grid_kmt(i,j)  !{

        ! Calculate average light level in the mixed layer
        if (wombat%zw(i,j,k).le.hblt_depth(i,j)) then
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

    wombat%no3_prev(:,:,:) = wombat%f_no3(:,:,:)
    wombat%caco3_prev(:,:,:) = wombat%f_caco3(:,:,:)

    ! Arrays for assessing conservation of mass within ecosystem component
    n_pools(:,:,:,:) = 0.0
    c_pools(:,:,:,:) = 0.0
    si_pools(:,:,:,:) = 0.0

    do tn = 1,ts_npzd  !{

      n_pools(:,:,:,1) = n_pools(:,:,:,2)
      c_pools(:,:,:,1) = c_pools(:,:,:,2)
      si_pools(:,:,:,1) = si_pools(:,:,:,2)

      do k = 1,nk; do j = jsc,jec; do i = isc,iec;

      ! Initialise some values and ratios (put into nicer units than mol/kg)
      biophy   = max(epsi, wombat%f_phy(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biodia   = max(epsi, wombat%f_dia(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biophyfe = max(epsi, wombat%f_phyfe(i,j,k))/ mmol_m3_to_mol_kg  ![mmol/m3]
      biodiafe = max(epsi, wombat%f_diafe(i,j,k))/ mmol_m3_to_mol_kg  ![mmol/m3]
      biozoo   = max(epsi, wombat%f_zoo(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biozoofe = max(epsi, wombat%f_zoofe(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biomes   = max(epsi, wombat%f_mes(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biomesfe = max(epsi, wombat%f_mesfe(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biodet   = max(epsi, wombat%f_det(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biobdet  = max(epsi, wombat%f_bdet(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biodoc   = max(epsi, wombat%f_doc(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biodon   = max(epsi, wombat%f_don(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biobac1  = max(epsi, wombat%f_bac1(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biobac2  = max(epsi, wombat%f_bac2(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      bioaoa   = max(epsi, wombat%f_aoa(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biono3   = max(epsi, wombat%f_no3(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      bion2o   = max(epsi, wombat%f_n2o(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      bionh4   = max(epsi, wombat%f_nh4(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biooxy   = max(epsi, wombat%f_o2(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biofer   = max(epsi, wombat%f_fe(i,j,k)  ) / umol_m3_to_mol_kg  ![umol/m3]
      biocaco3 = max(epsi, wombat%f_caco3(i,j,k))/ mmol_m3_to_mol_kg  ![mmol/m3]
      phy_chlc = max(epsi, wombat%f_pchl(i,j,k)) / max(epsi, wombat%f_phy(i,j,k))
      dia_chlc = max(epsi, wombat%f_dchl(i,j,k)) / max(epsi, wombat%f_dia(i,j,k))
      phy_Fe2C = max(epsi, wombat%f_phyfe(i,j,k))/ max(epsi, wombat%f_phy(i,j,k))
      dia_Fe2C = max(epsi, wombat%f_diafe(i,j,k))/ max(epsi, wombat%f_dia(i,j,k))
      zoo_Fe2C = max(epsi, wombat%f_zoofe(i,j,k))/ max(epsi, wombat%f_zoo(i,j,k))
      mes_Fe2C = max(epsi, wombat%f_mesfe(i,j,k))/ max(epsi, wombat%f_mes(i,j,k))
      det_Fe2C = max(epsi, wombat%f_detfe(i,j,k))/ max(epsi, wombat%f_det(i,j,k))
      bdet_Fe2C= max(epsi, wombat%f_bdetfe(i,j,k))/max(epsi, wombat%f_bdet(i,j,k))
      dom_N2C  = max(epsi, wombat%f_don(i,j,k))  / max(epsi, wombat%f_doc(i,j,k))
      dia_Si2C = max(epsi, wombat%f_diasi(i,j,k))/ max(epsi, wombat%f_dia(i,j,k))
      bdet_Si2C= max(epsi, wombat%f_bdetsi(i,j,k))/max(epsi, wombat%f_bdet(i,j,k))
      
    

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
      wombat%phy_kni(i,j,k) = wombat%phykn * max(0.1, max(0.0, (biophy-wombat%phybiot))**0.37)
      wombat%phy_kfe(i,j,k) = wombat%phykf * max(0.1, max(0.0, (biophy-wombat%phybiot))**0.37)
      ! Nitrogen limitation (split into NH4 and NO3 uptake limitation terms, with NH4 uptake being 5x more preferred)
      !   See Buchanan et al., (2025) Biogeosciences
      phy_limnh4 = bionh4 / (bionh4 + wombat%phy_kni(i,j,k) + epsi)
      phy_limno3 = biono3 / (biono3 + wombat%phy_kni(i,j,k) + epsi)
      phy_limdin = (biono3 + bionh4) / (biono3 + bionh4 + wombat%phy_kni(i,j,k) + epsi)
      wombat%phy_lnh4(i,j,k) = 5.0 * phy_limdin * phy_limnh4 / (phy_limno3 + 5.0 * phy_limnh4 + epsi)
      wombat%phy_lno3(i,j,k) = phy_limdin * phy_limno3 / (phy_limno3 + 5.0 * phy_limnh4 + epsi)
      wombat%phy_lnit(i,j,k) = wombat%phy_lno3(i,j,k) + wombat%phy_lnh4(i,j,k)
      ! Iron limitation (Quota model, constants from Flynn & Hipkin 1999)
      phy_minqfe = 0.00167 / 55.85 * max(wombat%phyminqc, phy_chlc)*12 + &
                   1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * wombat%phy_lnit(i,j,k) + &
                   1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * wombat%phy_lno3(i,j,k)
      wombat%phy_lfer(i,j,k) = min(1.0, max(0.0, (phy_Fe2C - phy_minqfe) / wombat%phyoptqf ))

      !!!~~~ Microphytoplankton ~~~!!!
      wombat%dia_kni(i,j,k) = wombat%diakn * max(0.1, max(0.0, (biodia-wombat%diabiot))**0.37)
      wombat%dia_kfe(i,j,k) = wombat%diakf * max(0.1, max(0.0, (biodia-wombat%diabiot))**0.37)
      ! Nitrogen limitation (split into NH4 and NO3 uptake limitation terms, with NH4 uptake being 5x more preferred)
      !   See Buchanan et al., (2025) Biogeosciences
      dia_limnh4 = bionh4 / (bionh4 + wombat%dia_kni(i,j,k) + epsi)
      dia_limno3 = biono3 / (biono3 + wombat%dia_kni(i,j,k) + epsi)
      dia_limdin = (biono3 + bionh4) / (biono3 + bionh4 + wombat%dia_kni(i,j,k) + epsi)
      wombat%dia_lnh4(i,j,k) = 5.0 * dia_limdin * phy_limnh4 / (dia_limno3 + 5.0 * dia_limnh4 + epsi)
      wombat%dia_lno3(i,j,k) = dia_limdin * dia_limno3 / (dia_limno3 + 5.0 * dia_limnh4 + epsi)
      wombat%dia_lnit(i,j,k) = wombat%dia_lno3(i,j,k) + wombat%dia_lnh4(i,j,k)
      ! Iron limitation (Quota model, constants from Flynn & Hipkin 1999)
      dia_minqfe = 0.00167 / 55.85 * max(wombat%diaminqc, dia_chlc)*12 + &
                   1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * wombat%dia_lnit(i,j,k) + &
                   1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * wombat%dia_lno3(i,j,k)
      wombat%dia_lfer(i,j,k) = min(1.0, max(0.0, (dia_Fe2C - dia_minqfe) / wombat%diaoptqf ))


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 3] Heterotrophy and remineralisation                           !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Temperature dependance of heterotrophy (applies to bact and zoo)
      fbc = wombat%bbioh ** (Temp(i,j,k))

      ! Variable rates of remineralisation
      wombat%reminr(i,j,k) = wombat%detlrem * fbc


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 4] Light limitation of phytoplankton                           !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      
      ! 1. initial slope of Photosynthesis-Irradiance curve
      ! 2. Alter the slope to account for respiration and daylength limitation
      ! 3. Light limitation 
      ! 4. Apply light and nutrient limitations to maximum growth rate
      ! 5. Determine difference in light-limited and realized growth to get DOC exudation

      !!!~~~ Phytoplankton ~~~!!!
      phy_pisl  = max(wombat%alphabio_phy * phy_chlc, wombat%alphabio_phy * wombat%phyminqc)
      phy_pisl2 = phy_pisl / ( (1. + wombat%phylmor*86400.0 * fbc) ) ! add daylength estimate here
      wombat%phy_lpar(i,j,k) = (1. - exp(-phy_pisl2 * wombat%radbio(i,j,k)))
      wombat%phy_mu(i,j,k) = wombat%phy_mumax(i,j,k) * wombat%phy_lpar(i,j,k) * & 
                             min(wombat%phy_lnit(i,j,k), wombat%phy_lfer(i,j,k))
      !!!~~~ Microphytoplankton ~~~!!!
      dia_pisl  = max(wombat%alphabio_dia * dia_chlc, wombat%alphabio_dia * wombat%diaminqc)
      dia_pisl2 = dia_pisl / ( (1. + wombat%dialmor*86400.0 * fbc) ) ! add daylength estimate here
      wombat%dia_lpar(i,j,k) = (1. - exp(-dia_pisl2 * wombat%radbio(i,j,k)))
      wombat%dia_mu(i,j,k) = wombat%dia_mumax(i,j,k) * wombat%dia_lpar(i,j,k) * & 
                             min(wombat%dia_lnit(i,j,k), wombat%dia_lfer(i,j,k))

      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 5] Growth of chlorophyll                                       !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! 1. Light limitation of chlorophyll production
      ! 2. minimum and optimal rates of chlorophyll growth
      ! 3. Calculate mg Chl m-3 s-1

      ! Reduced chlorophyll growth during extended periods of darkness
      phi = wombat%radmld(i,j,1) / (wombat%radmld(i,j,1) + wombat%chlkWm2)

      !!!~~~ Phytoplankton ~~~!!!
      pchl_pisl = phy_pisl / ( wombat%phy_mumax(i,j,k) * 86400.0 * & 
                  (1. - min(wombat%phy_lnit(i,j,k), wombat%phy_lfer(i,j,k))) + epsi )
      wombat%pchl_lpar(i,j,k) = (1. - exp(-pchl_pisl * wombat%radmld(i,j,k)))
      pchl_mumin = wombat%phyminqc * wombat%phy_mu(i,j,k) * biophy * 12.0   ![mg/m3/s] 
      pchl_muopt = wombat%phyoptqc * wombat%phy_mu(i,j,k) * biophy * 12.0   ![mg/m3/s]
      wombat%pchl_mu(i,j,k) = (pchl_muopt - pchl_mumin) * wombat%pchl_lpar(i,j,k) * &
                               min(wombat%phy_lnit(i,j,k), wombat%phy_lfer(i,j,k))
      if ( (phy_pisl * wombat%radmld(i,j,k)) .gt. 0.0 ) then
        wombat%pchl_mu(i,j,k) = phi * ( pchl_mumin + wombat%pchl_mu(i,j,k) / &
                                (phy_pisl * wombat%radmld(i,j,k)) )
      endif
      wombat%pchl_mu(i,j,k) = wombat%pchl_mu(i,j,k) / 12.0 * mmol_m3_to_mol_kg  ![mol/kg/s]

      !!!~~~ Microphytoplankton ~~~!!!
      dchl_pisl = dia_pisl / ( wombat%dia_mumax(i,j,k) * 86400.0 * & 
                  (1. - min(wombat%dia_lnit(i,j,k), wombat%dia_lfer(i,j,k))) + epsi )
      wombat%dchl_lpar(i,j,k) = (1. - exp(-dchl_pisl * wombat%radmld(i,j,k)))
      dchl_mumin = wombat%diaminqc * wombat%dia_mu(i,j,k) * biodia * 12.0   ![mg/m3/s] 
      dchl_muopt = wombat%diaoptqc * wombat%dia_mu(i,j,k) * biodia * 12.0   ![mg/m3/s]
      wombat%dchl_mu(i,j,k) = (dchl_muopt - dchl_mumin) * wombat%dchl_lpar(i,j,k) * &
                               min(wombat%dia_lnit(i,j,k), wombat%dia_lfer(i,j,k))
      if ( (dia_pisl * wombat%radmld(i,j,k)) .gt. 0.0 ) then
        wombat%dchl_mu(i,j,k) = phi * ( dchl_mumin + wombat%dchl_mu(i,j,k) / &
                                (dia_pisl * wombat%radmld(i,j,k)) )
      endif
      wombat%dchl_mu(i,j,k) = wombat%dchl_mu(i,j,k) / 12.0 * mmol_m3_to_mol_kg  ![mol/kg/s]


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 6] Phytoplankton uptake of iron                                !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! 1. Maximum iron content of phytoplankton cell
      ! 2. Ensure that dFe uptake increases or decreases in response to cell quota
      ! 3. Iron uptake of phytoplankton (reduced to 20% at night and when N is limiting)

      !!!~~~ Phytoplankton ~~~!!!
      phy_maxqfe = biophy * wombat%phymaxqf  !mmol Fe / m3
      wombat%phy_feupreg(i,j,k) = (4.0 - 4.5 * wombat%phy_lfer(i,j,k) / &
                                  (wombat%phy_lfer(i,j,k) + 0.5) )
      wombat%phy_fedoreg(i,j,k) = max(0.0, (1.0 - biophyfe/phy_maxqfe) / &
                                  abs(1.05 - biophyfe/phy_maxqfe) )
      wombat%phy_dfeupt(i,j,k) = (wombat%phy_mumax(i,j,k) * wombat%phymaxqf * &
                                  max(0.2, wombat%phy_lpar(i,j,k) * wombat%phy_lnit(i,j,k)) * &
                                  biofer / (biofer + wombat%phy_kfe(i,j,k)) * &
                                  wombat%phy_feupreg(i,j,k) * &
                                  wombat%phy_fedoreg(i,j,k) * biophy) * mmol_m3_to_mol_kg

      !!!~~~ Microphytoplankton ~~~!!!
      dia_maxqfe = biodia * wombat%diamaxqf  !mmol Fe / m3
      wombat%dia_feupreg(i,j,k) = (4.0 - 4.5 * wombat%dia_lfer(i,j,k) / &
                                  (wombat%dia_lfer(i,j,k) + 0.5) )
      wombat%dia_fedoreg(i,j,k) = max(0.0, (1.0 - biodiafe/dia_maxqfe) / &
                                  abs(1.05 - biodiafe/dia_maxqfe) )
      wombat%dia_dfeupt(i,j,k) = (wombat%dia_mumax(i,j,k) * wombat%diamaxqf * &
                                  max(0.2, wombat%dia_lpar(i,j,k) * wombat%dia_lnit(i,j,k)) * &
                                  biofer / (biofer + wombat%dia_kfe(i,j,k)) * &
                                  wombat%dia_feupreg(i,j,k) * &
                                  wombat%dia_fedoreg(i,j,k) * biodia) * mmol_m3_to_mol_kg


      !-------------------------------------------------------------------------------------!
      !-------------------------------------------------------------------------------------!
      !-------------------------------------------------------------------------------------!
      !  [Step 7] Iron chemistry (Aumont et al., 2015 GMD & Tagliabue et al., 2023 Nature)  !
      !-------------------------------------------------------------------------------------!
      !-------------------------------------------------------------------------------------!
      !-------------------------------------------------------------------------------------!

      ! Estimate solubility of Fe3+ (free Fe) in solution using temperature, pH and salinity
      ztemk = max(5.0, Temp(i,j,k)) + 273.15    ! temperature in kelvin
      I_ztemk = 1.0 / ztemk
      zval = 19.924 * Salt(i,j,k) / ( 1000. - 1.005 * Salt(i,j,k))
      sqrt_zval = sqrt(zval)
      fesol1 = 10.0**(-13.486 - 0.1856*sqrt_zval + 0.3073*zval + 5254.0*I_ztemk)
      fesol2 = 10.0**(2.517 - 0.8885*sqrt_zval + 0.2139*zval - 1320.0*I_ztemk)
      fesol3 = 10.0**(0.4511 - 0.3305*sqrt_zval - 1996.0*I_ztemk)
      fesol4 = 10.0**(-0.2965 - 0.7881*sqrt_zval - 4086.0*I_ztemk)
      fesol5 = 10.0**(4.4466 - 0.8505*sqrt_zval - 7980.0*I_ztemk)
      if (wombat%htotal(i,j,k).gt.0.0) then
        hp = wombat%htotal(i,j,k)
      else
        hp = 1.25893e-08 ! dts: =10.0**(-7.9)
      endif
      fe3sol = fesol1 * ( hp*hp*hp + fesol2*hp*hp + fesol3*hp + fesol4 + fesol5/hp ) *1e9

      ! Estimate total colloidal iron
      ! ... for now, we assume that 50% of all dFe is colloidal, and we separate this from the 
      !     equilibrium fractionation between Fe' and Fe-L below
      wombat%fecol(i,j,k) = wombat%fcolloid * biofer 

      ! Determine equilibriuim fractionation of the remain dFe (non-colloidal fraction) into Fe' and L-Fe
      fe_keq = 10.0**( 17.27 - 1565.7 * I_ztemk ) * 1e-9 ! Temperature reduces solubility
      fe_par = 0.47587 * wombat%radbio(i,j,k) ! Light increases solubility. dts: 0.47587=4.77e-7*0.5/10.0**(-6.3)
      fe_sfe = max(0.0, biofer - wombat%fecol(i,j,k))
      fe_tfe = (1.0 + fe_par) * fe_sfe
      wombat%feIII(i,j,k) = ( -( 1. + wombat%ligand * fe_keq + fe_par - fe_sfe * fe_keq ) &
                             + SQRT( ( 1. + wombat%ligand * fe_keq + fe_par - fe_sfe * fe_keq )**2 &
                                      + 4. * fe_tfe * fe_keq) ) / ( 2. * fe_keq + epsi )
      wombat%felig(i,j,k) = max(0.0, fe_sfe - wombat%feIII(i,j,k))

      ! Precipitation of Fe' (creation of nanoparticles)
      wombat%feprecip(i,j,k) = max(0.0, ( wombat%feIII(i,j,k) - fe3sol ) ) * wombat%knano_dfe/86400.0

      ! Scavenging of Fe` onto biogenic particles 
      partic = (biodet + biobdet + biocaco3)
      wombat%fescaven(i,j,k) = wombat%feIII(i,j,k) * (1e-7 + wombat%kscav_dfe * partic) / 86400.0
      wombat%fescadet(i,j,k) = wombat%fescaven(i,j,k) * biodet / (partic+epsi) 
      wombat%fescabdet(i,j,k) = wombat%fescaven(i,j,k) * biobdet / (partic+epsi) 

      ! Coagulation of colloidal Fe (umol/m3) to form sinking particles (mmol/m3)
      ! Following Tagliabue et al. (2023), make coagulation rate dependent on DOC and Phytoplankton biomass
      biof = max(1/3., biophy / (biophy + 0.03))
      ! Colloidal shunt associated with small particles and DOC (Tagliabue et al., 2023)
      if (wombat%zw(i,j,k).le.hblt_depth(i,j)) then
        zval = (      (12.*biof*biodoc + 9.*biodet) + 2.5*biodet + 128.*biof*biodoc + 725.*biodet )*wombat%kcoag_dfe
      else
        zval = ( 0.01*(12.*biof*biodoc + 9.*biodet) + 2.5*biodet + 128.*biof*biodoc + 725.*biodet )*wombat%kcoag_dfe
      endif
      wombat%fecoag2det(i,j,k) = wombat%fecol(i,j,k) * zval / 86400.0
      ! Colloidal shunt associated with big particles (Tagliabue et al., 2023)
      if (wombat%zw(i,j,k).le.hblt_depth(i,j)) then
        zval = ((      2.0 + 1.37) * biobdet )*wombat%kcoag_dfe
      else
        zval = (( 0.01*2.0 + 1.37) * biobdet )*wombat%kcoag_dfe
      endif
      wombat%fecoag2bdet(i,j,k) = wombat%fecol(i,j,k) * zval / 86400.0
      
      ! Convert the terms back to mol/kg
      wombat%feprecip(i,j,k) = wombat%feprecip(i,j,k) * umol_m3_to_mol_kg
      wombat%fescaven(i,j,k) = wombat%fescaven(i,j,k) * umol_m3_to_mol_kg
      wombat%fescadet(i,j,k) = wombat%fescadet(i,j,k) * umol_m3_to_mol_kg
      wombat%fescabdet(i,j,k) = wombat%fescabdet(i,j,k) * umol_m3_to_mol_kg
      wombat%fecoag2det(i,j,k) = wombat%fecoag2det(i,j,k) * umol_m3_to_mol_kg
      wombat%fecoag2bdet(i,j,k) = wombat%fecoag2bdet(i,j,k) * umol_m3_to_mol_kg
      wombat%feIII(i,j,k) = wombat%feIII(i,j,k) * umol_m3_to_mol_kg
      wombat%felig(i,j,k) = wombat%felig(i,j,k) * umol_m3_to_mol_kg
      wombat%fecol(i,j,k) = wombat%fecol(i,j,k) * umol_m3_to_mol_kg


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 8] Silicic acid chemistry                                      !
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
      gamma0 = 1.0 + 0.0053 * Salt(i,j,k) - 0.000034 * Salt(i,j,k)**2 ! Svenko 2014 Mar. Chem.
      alphaH2O = 0.999  ! activity of water in seawater (from TEOS-10)
      deltaV0 = -9.9 * 1e-6 ! m3/mol, (Willey 1982 Geochim. et Cosmochim. Acta) (also see Loucaides et al., 2012 Mar. Chem.)
      spmvcorrect = exp( -deltaV0/(Rgas * zval) * wombat%zm(i,j,k) * 1.0e4 ) ! 1.0e4 converts dbar to Pa
      wombat%sileqc(i,j,k) = (K_am_silica * spmvcorrect) * alphaH2O**2.0 / gamma0 ! mol/kg

      ! Dissolution of biogenic silica into silicic acid
        ! 1. Temperature effect (Arrhenius) 
        !    - activation energy of 84 kJ/mol (Greenwood et al., 2005 Aqu. Geochem.)
        !    - Greenwood et al. (2005) finds dissolution rates [1/s] at a given temperature
        !      equal to exp(20 - 10050/T), with T in Kelvin (Figure 4). We use 40ºC, which is
        !      their lowest temperature, thus 313.15 in Kelvin
        !    - NOTE: Greenwood also estimate a pH dependence, but make their rate measurements 
        !            at pH 10 - 14. We ignore this for now, but it may be worth revisiting
      disssi_temp = exp(20.0 - 10050.0/313.15) * exp( -84.0*1e3/Rgas * (1./zval - 1./313.15) ) ! [1/s]
        ! 2. Undersaturation term with an exponent of 2.0
        !    - see Eq. 2.13 and fits of this equation to ocean data in Figures 3.20 and 3.21 in
        !      Rickert, D., Dissolution kinetics of biogenic silica in marine environments, Ber. Polarforsch., 351, 2000.
        !    - From Van Cappellen et al., (2002) Global Biogeochemical Cycles:
        !      "Detailed kinetic studies of biogenic silica dissolution conducted in flow-through reactors 
        !       demonstrate that at very high degrees of undersaturation the dissolution kinetics switch 
        !       from a linear dependence on the degree of undersaturation to an exponential one [Van Cappellen 
        !       and Qiu, 1997b; Rickert, 2000]."
        !    - We therefore assume substantial undersaturation, which is the case in most of the ocean
      disssi_usat = max(0.0, (1 - wombat%f_sil(i,j,k) / wombat%sileqc(i,j,k))**2.0 )
        ! 3. Bio-interference term?
        !    - "The removal of organic or inorganic coatings enhance the reactivity by at least an order of magnitude."
        !      Ricket et al., 2002 Geochim. et Cosmochim. Acta
        !    - Diatom frustule dissolution increased by order of magnitude with bacteria (Bidle & Azam 1999 Nature) 
        !    - During a bloom off Monterey Bay, anti-biotics decreased dissolution by ~50% (Bidle et al., 2003 Limnol. Oceanogr.)
      disssi_bact = 1.0 + wombat%bsi_fbac * (wombat%f_bac1(i,j,k) + wombat%f_bac2(i,j,k)) &
                    / ( wombat%f_bac1(i,j,k) + wombat%f_bac2(i,j,k) + wombat%bsi_kbac )
        ! 4. Dissolution rate of biogenic silica (/s) composed of the above terms
      wombat%disssi(i,j,k) = disssi_temp * disssi_usat * disssi_bact


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 9] Mortality scalings and grazing                              !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! reduce linear mortality (respiration losses) of zooplankton when there is low biomass
      zoo_slmor = biozoo / (biozoo + wombat%zookz)
      mes_slmor = biomes / (biomes + wombat%meskz)
     
      !!!~~~ Zooplankton ~~~!!!
      ! Grazing function ! [1/s]
      ! normalize the prey preference kernal to reflect dietary fractions (Gentleman et al., (2003) DSRII)
      denom = (wombat%zprefbac1 + wombat%zprefbac2 + wombat%zprefaoa + wombat%zprefphy + wombat%zprefdia + wombat%zprefdet)
      wombat%zooprefbac1(i,j,k) = wombat%zprefbac1 / denom
      wombat%zooprefbac2(i,j,k) = wombat%zprefbac2 / denom
      wombat%zooprefaoa(i,j,k) = wombat%zprefaoa / denom
      wombat%zooprefphy(i,j,k) = wombat%zprefphy / denom
      wombat%zooprefdia(i,j,k) = wombat%zprefdia / denom
      wombat%zooprefdet(i,j,k) = wombat%zprefdet / denom
      
      ! Gentleman et al. (2003) DSRII 
      !   - add a switching component designed to weight the diet towards abundant prey
      !   - see their Eq. 19
      ! Emulates empirical basis of selective feeding on more abundant prey (Kiorboe et al., 2017; L&O)
      !   ... if denominator is zero, then set all preferences to 1/3 (this is a failsafe, but it should not happen)
      denom = wombat%zooprefbac1(i,j,k) + wombat%zooprefbac2(i,j,k) + wombat%zooprefaoa(i,j,k) &
              + wombat%zooprefphy(i,j,k) + wombat%zooprefdia(i,j,k) + wombat%zooprefdet(i,j,k)
      if (denom .lt. epsi) then
        wombat%zooprefbac1(i,j,k) = 1.0/6.0; wombat%zooprefbac2(i,j,k) = 1.0/6.0; wombat%zooprefaoa(i,j,k) = 1.0/6.0
        wombat%zooprefphy(i,j,k) = 1.0/6.0; wombat%zooprefdia(i,j,k) = 1.0/6.0; wombat%zooprefdet(i,j,k) = 1.0/6.0
      else
        wzbac1 = (wombat%zooprefbac1(i,j,k) * biobac1)**wombat%zoopreyswitch
        wzbac2 = (wombat%zooprefbac2(i,j,k) * biobac2)**wombat%zoopreyswitch
        wzaoa = (wombat%zooprefaoa(i,j,k) * bioaoa)**wombat%zoopreyswitch
        wzphy = (wombat%zooprefphy(i,j,k) * biophy)**wombat%zoopreyswitch
        wzdia = (wombat%zooprefdia(i,j,k) * biodia)**wombat%zoopreyswitch
        wzdet = (wombat%zooprefdet(i,j,k) * biodet)**wombat%zoopreyswitch
        wzsum = wzbac1 + wzbac2 + wzaoa + wzphy + wzdia + wzdet + epsi
        wombat%zooprefbac1(i,j,k) = wzbac1 / wzsum
        wombat%zooprefbac2(i,j,k) = wzbac2 / wzsum
        wombat%zooprefaoa(i,j,k) = wzaoa / wzsum
        wombat%zooprefphy(i,j,k) = wzphy / wzsum
        wombat%zooprefdia(i,j,k) = wzdia / wzsum
        wombat%zooprefdet(i,j,k) = wzdet / wzsum
      endif
      ! Compute sum of prey-specific Type-III terms to obtain grazing rate [1/s]
      !  - this avoids "perfect substitution" of prey types and aligns with reccommendations of Gentleman et al. (2003)
      Xzoo = (  wombat%zooepsbac1 * (wombat%zooprefbac1(i,j,k) * biobac1)**2 &
              + wombat%zooepsbac2 * (wombat%zooprefbac2(i,j,k) * biobac2)**2 &
              + wombat%zooepsaoa * (wombat%zooprefaoa(i,j,k) * bioaoa)**2 &
              + wombat%zooepsphy * (wombat%zooprefphy(i,j,k) * biophy)**2 &
              + wombat%zooepsdia * (wombat%zooprefdia(i,j,k) * biodia)**2 &
              + wombat%zooepsdet * (wombat%zooprefdet(i,j,k) * biodet)**2)
      g_npz = wombat%zoogmax * fbc * Xzoo / (wombat%zoogmax * fbc + Xzoo)
      ! find "apparent" community epsilon (prey capture rate coefficient)
      wombat%zooeps(i,j,k) = Xzoo / ( (wombat%zooprefbac1(i,j,k) * biobac1)**2 &
                                    + (wombat%zooprefbac2(i,j,k) * biobac2)**2 &
                                    + (wombat%zooprefaoa(i,j,k) * bioaoa)**2 &
                                    + (wombat%zooprefphy(i,j,k) * biophy)**2 &
                                    + (wombat%zooprefdia(i,j,k) * biodia)**2 &
                                    + (wombat%zooprefdet(i,j,k) * biodet)**2 )
      

      !!!~~~ Mesozooplankton ~~~!!!
      ! Grazing function ! [1/s]
      ! normalize the prey preference kernal to reflect dietary fractions (Gentleman et al., (2003) DSRII)
      denom = ( wombat%mprefbac1 + wombat%mprefbac2 + wombat%mprefaoa + wombat%mprefphy &
               + wombat%mprefdia + wombat%mprefdet + wombat%mprefzoo )
      wombat%mesprefbac1(i,j,k) = wombat%mprefbac1 / denom
      wombat%mesprefbac2(i,j,k) = wombat%mprefbac2 / denom
      wombat%mesprefaoa(i,j,k) = wombat%mprefaoa / denom
      wombat%mesprefphy(i,j,k) = wombat%mprefphy / denom
      wombat%mesprefdia(i,j,k) = wombat%mprefdia / denom
      wombat%mesprefdet(i,j,k) = wombat%mprefdet / denom
      wombat%mesprefbdet(i,j,k) = wombat%mprefbdet / denom
      wombat%mesprefzoo(i,j,k) = wombat%mprefzoo / denom
      denom = wombat%mesprefbac1(i,j,k) + wombat%mesprefbac2(i,j,k) + wombat%mesprefaoa(i,j,k) + wombat%mesprefphy(i,j,k) &
              + wombat%mesprefdia(i,j,k) + wombat%mesprefdet(i,j,k) + wombat%mesprefbdet(i,j,k) + wombat%mesprefzoo(i,j,k)
      ! Gentleman et al. (2003) DSRII 
      !   - add a switching component designed to weight the diet towards abundant prey
      !   - see their Eq. 19
      ! Emulates empirical basis of selective feeding on more abundant prey (Kiorboe et al., 2017; L&O)
      !   ... if denominator is zero, then set all preferences to 1/3 (this is a failsafe, but it should not happen)
      if (denom .lt. 1e-20) then
        wombat%mesprefbac1(i,j,k) = 1.0/8.0; wombat%mesprefbac2(i,j,k) = 1.0/8.0; wombat%mesprefaoa(i,j,k) = 1.0/8.0
        wombat%mesprefphy(i,j,k) = 1.0/8.0; wombat%mesprefdia(i,j,k) = 1.0/8.0; wombat%mesprefdet(i,j,k) = 1.0/8.0
        wombat%mesprefbdet(i,j,k) = 1.0/8.0; wombat%mesprefzoo(i,j,k) = 1.0/8.0
      else
        wzbac1 = (wombat%mesprefbac1(i,j,k) * biobac1)**wombat%mespreyswitch
        wzbac2 = (wombat%mesprefbac2(i,j,k) * biobac2)**wombat%mespreyswitch
        wzaoa = (wombat%mesprefaoa(i,j,k) * bioaoa)**wombat%mespreyswitch
        wzphy = (wombat%mesprefphy(i,j,k) * biophy)**wombat%mespreyswitch
        wzdia = (wombat%mesprefdia(i,j,k) * biodia)**wombat%mespreyswitch
        wzdet = (wombat%mesprefdet(i,j,k) * biodet)**wombat%mespreyswitch
        wzbdet= (wombat%mesprefbdet(i,j,k) * biobdet)**wombat%mespreyswitch
        wzzoo = (wombat%mesprefzoo(i,j,k) * biozoo)**wombat%mespreyswitch
        wzsum = wzbac1 + wzbac2 + wzaoa + wzphy + wzdia + wzdet + wzbdet + wzzoo + epsi
        wombat%mesprefbac1(i,j,k) = wzbac1 / wzsum
        wombat%mesprefbac2(i,j,k) = wzbac2 / wzsum
        wombat%mesprefaoa(i,j,k) = wzaoa / wzsum
        wombat%mesprefphy(i,j,k) = wzphy / wzsum
        wombat%mesprefdia(i,j,k) = wzdia / wzsum
        wombat%mesprefdet(i,j,k) = wzdet / wzsum
        wombat%mesprefbdet(i,j,k) = wzbdet / wzsum
        wombat%mesprefzoo(i,j,k) = wzzoo / wzsum
      endif
      ! Compute sum of prey-specific Type-III terms to obtain grazing rate [1/s]
      !  - this avoids "perfect substitution" of prey types and aligns with reccommendations of Gentleman et al. (2003)
      Xmes = (  wombat%mesepsbac1 * (wombat%mesprefbac1(i,j,k) * biobac1)**2 &
              + wombat%mesepsbac2 * (wombat%mesprefbac2(i,j,k) * biobac2)**2 &
              + wombat%mesepsaoa * (wombat%mesprefaoa(i,j,k) * bioaoa)**2 &
              + wombat%mesepsphy * (wombat%mesprefphy(i,j,k) * biophy)**2 &
              + wombat%mesepsdia * (wombat%mesprefdia(i,j,k) * biodia)**2 &
              + wombat%mesepsdet * (wombat%mesprefdet(i,j,k) * biodet)**2 &
              + wombat%mesepsdet * (wombat%mesprefbdet(i,j,k) * biobdet)**2 &
              + wombat%mesepsdet * (wombat%mesprefzoo(i,j,k) * biozoo)**2 )
      m_npz = wombat%mesgmax * fbc * Xmes / (wombat%mesgmax * fbc + Xmes)
      ! find "apparent" community epsilon (prey capture rate coefficient)
      wombat%meseps(i,j,k) = Xmes / ( (wombat%mesprefbac1(i,j,k) * biobac1)**2 &
                                    + (wombat%mesprefbac2(i,j,k) * biobac2)**2 &
                                    + (wombat%mesprefaoa(i,j,k) * bioaoa)**2 &
                                    + (wombat%mesprefphy(i,j,k) * biophy)**2 &
                                    + (wombat%mesprefdia(i,j,k) * biodia)**2 &
                                    + (wombat%mesprefdet(i,j,k) * biodet)**2 &
                                    + (wombat%mesprefbdet(i,j,k) * biobdet)**2 &
                                    + (wombat%mesprefzoo(i,j,k) * biozoo)**2 )

    
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 10] CaCO3 calculations                                         !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      if (do_caco3_dynamics) then

        ! PIC:POC ratio is a function of the substrate:inhibitor ratio, which is the 
        !  HCO3- to free H+ ions ratio (mol/umol), following Lehmann & Bach (2024). 
        !  We also add a T-dependent function to scale down CaCO3 production in waters colder 
        !  than 3 degrees C based off the observation of no E hux growth beneath this (Fielding 2013; L&O)
        hco3 = wombat%f_dic(i,j,k) - wombat%co3(i,j,k) - wombat%co2_star(i,j,k)
        wombat%pic2poc(i,j,k) = min(0.3, (wombat%f_inorg + 10.0**(-3.0 + 4.31e-6 * &
                                          hco3 / wombat%htotal(i,j,k))) * &
                                         (0.55 + 0.45 * tanh(Temp(i,j,k) - 4.0)) )

        ! The dissolution rate is a function of omegas for calcite and aragonite, as well the
        !  concentration of POC, following Kwon et al., 2024, Science Advances; Table S1
        diss_cal = (wombat%disscal * max(0.0, 1.0 - wombat%omega_cal(i,j,k))**2.2) / 86400.0
        diss_ara = (wombat%dissara * max(0.0, 1.0 - wombat%omega_ara(i,j,k))**1.5) / 86400.0
        diss_det = (wombat%dissdet * wombat%reminr(i,j,k) * (biodet**2.0 + biobdet**2.0))
        wombat%dissrat(i,j,k) = diss_cal + diss_ara + diss_det

      else
      
        wombat%pic2poc(i,j,k) = wombat%f_inorg + 0.025
        wombat%dissrat(i,j,k) = wombat%caco3lrem
      
      endif


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 11] Implicit nitrogen fixation                                 !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      
      if (do_nitrogen_fixation) then
        ! Temperature dependent maximum growth rate of Trichodesmium (Jiang et al., 2018)
        if (Temp(i,j,k).gt.17.2) then
          wombat%trimumax(i,j,k) = ( ( -3.99e-4 * Temp(i,j,k)**3.0 ) + &
                                     (  0.02685 * Temp(i,j,k)**2.0 ) + & 
                                     ( -0.555 * Temp(i,j,k) ) + 3.633 ) / 86400.0
        endif
        ! Nutrient and light limitation terms
        wombat%tri_lfer(i,j,k) = biofer / (biofer + wombat%trikf)
        wombat%tri_lpar(i,j,k) = (1. - exp(-wombat%alphabio_tri * wombat%trichlc * wombat%radbio(i,j,k)))
        ! Nitrogen fixation rate of Trichodesmium
        wombat%nitrfix(i,j,k) = wombat%trimumax(i,j,k) * (1.0 - wombat%phy_lnit(i,j,k)) &
                                * min(wombat%tri_lfer(i,j,k), wombat%tri_lpar(i,j,k)) &
                                * wombat%trin2c * 1e-6  ! 1e-6 scaler to account for biomass in mol/kg
      endif
      

      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 12] Facultative bacterial heterotrophy                         !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Assumptions:
        ! 1. We treat DOC and DON as one molecule of DOM, so they must be treated together
        ! 2. Bacteria are initially limited by DOC and electron acceptors, reflecting cellular physiology where energy aquisition is prioritized
        ! 3. Bacteria take up DON in its ratio with DOC, but may supplement their growth with available NH4 if the DOC:DON ratio is > 5:1
        ! 4. In the case that not enough NH4 is available to support DOC-limited growth, growth becomes N-limited
        ! 5. The affinity of bacteria for DOC decreases as the DOC:DON ratio increases
      
      !!!~~~ Bacterial type #1 ~~~!!!
      ! Uptake of DOC (i.e., DOC-limited growth)
      wombat%bac1_kdoc(i,j,k) = wombat%bac1_kdoc_min + 2.0 * (wombat%bac1_kdoc_max - wombat%bac1_kdoc_min) &
                               * max(0.0, min(1.0, dom_N2C / (dom_N2C + (1.0/wombat%bac1_C2N)) )) 
      bac_Vdoc = wombat%bac1_Vmax_doc * biodoc / (biodoc + wombat%bac1_kdoc(i,j,k))
      ! Aerobic growth
      bac_Voxy = biooxy * wombat%bac1_poxy
      bac_muaer = max(0.0, min( (bac_Voxy/wombat%bac1_yoxy), (bac_Vdoc/wombat%bac1_yaerC) ) ) * fbc
      ! Anaerobic growth (will always be lower than aerobic growth when DOC is limiting)
      bac_Vno3 = wombat%bac1_Vmax_no3 * biono3 / (biono3 + wombat%bac1_kno3)
      bac_muana = max(0.0, min( (bac_Vno3/wombat%bac1_yno3), (bac_Vdoc/wombat%bac1_yanaC) ) ) * fbc
      if (.not.do_wc_denitrification) bac_muana = 0.0 ! If no denitrification, anaerobic growth is zero
      ! Save occurance of anaerobic growth to array
      if (bac_muana.gt.bac_muaer) wombat%bac1_fanaer(i,j,k) = 1.0
      ! Take the maximum growth rate as the realised growth rate
      wombat%bac1_mu(i,j,k) = max(bac_muaer, bac_muana)


      !!!~~~ Bacterial type #2 ~~~!!!
      ! Uptake of DOC (i.e., DOC-limited growth)
      wombat%bac2_kdoc(i,j,k) = wombat%bac2_kdoc_min + 2.0 * (wombat%bac2_kdoc_max - wombat%bac2_kdoc_min) &
                               * max(0.0, min(1.0, dom_N2C / (dom_N2C + (1.0/wombat%bac2_C2N)) )) 
      bac_Vdoc = wombat%bac2_Vmax_doc * biodoc / (biodoc + wombat%bac2_kdoc(i,j,k))
      ! Aerobic growth
      bac_Voxy = biooxy * wombat%bac2_poxy
      bac_muaer = max(0.0, min( (bac_Voxy/wombat%bac2_yoxy), (bac_Vdoc/wombat%bac2_yaerC) ) ) * fbc
      ! Anaerobic growth (will always be lower than aerobic growth when DOC is limiting)
      bac_Vn2o = bion2o * wombat%bac2_pn2o
      bac_muana = max(0.0, min( (bac_Vn2o/wombat%bac2_yn2o), (bac_Vdoc/wombat%bac2_yanaC) ) ) * fbc
      if (.not.do_wc_denitrification) bac_muana = 0.0 ! If no denitrification, anaerobic growth is zero
      ! Save occurance of anaerobic growth to array
      if (bac_muana.gt.bac_muaer) wombat%bac2_fanaer(i,j,k) = 1.0
      ! Take the maximum growth rate as the realised growth rate
      wombat%bac2_mu(i,j,k) = max(bac_muaer, bac_muana)
      
      
      ! Determine if bacteria are limited by N or Fe
      if (wombat%bac1_mu(i,j,k)*wombat%f_bac1(i,j,k).gt.0.0) then
        ! Initial estimate of the C biomass growth, DOC and DON assimilation rate by bacteria
        wombat%bac1grow(i,j,k) = wombat%bac1_mu(i,j,k) * wombat%f_bac1(i,j,k) ! [molC/kg/s]
        wombat%doc1remi(i,j,k) = wombat%bac1_mu(i,j,k) * wombat%f_bac1(i,j,k) * wombat%bac1_yaerC * (1. - wombat%bac1_fanaer(i,j,k)) & 
                                + wombat%bac1_mu(i,j,k) * wombat%f_bac1(i,j,k) * wombat%bac1_yanaC * wombat%bac1_fanaer(i,j,k) ! [molC/kg/s]
        wombat%don1remi(i,j,k) = wombat%doc1remi(i,j,k) * dom_N2C ! [molN/kg/s]

        ! Determine degree of N limitation of bacteria and adjust growth rate
        if (bionh4.gt.1e-3) then ! First, make sure that NH4 is available for uptake
          ! Find out what they need to supplement growth (this assumes that bacteria can achieve a N biomass yield of 1.0)
          wombat%bac1unh4(i,j,k) = wombat%bac1grow(i,j,k) / wombat%bac1_C2N - wombat%don1remi(i,j,k) ! [molN/kg/s]
          ! Find limitation of NH4 uptake
          if (wombat%bac1unh4(i,j,k).gt.0.0) then ! NH4 is needed to support growth
            ! Apply uptake kinetics constraint on bacterial NH4 uptake
            bac_limnh4 = bionh4 / (bionh4 + wombat%bac1_knh4 + epsi)
            wombat%bac1unh4(i,j,k) = wombat%bac1unh4(i,j,k) * bac_limnh4 ! [molN/kg/s]
            ! Make sure bacteria don't remove more NH4 than is available (can take up as much as half at a time)
            if (wombat%bac1unh4(i,j,k).gt.wombat%f_nh4(i,j,k)*0.25/dt) then
              wombat%bac1unh4(i,j,k) = wombat%f_nh4(i,j,k)*0.25 / dt ! [molN/kg/s]
            endif
            ! Recompute N limitation of bacteria (DON + NH4 uptake)
            wombat%bac1_lnit(i,j,k) = max(0.0, min(1.0, (wombat%don1remi(i,j,k) + wombat%bac1unh4(i,j,k)) &
                                                         / (wombat%bac1grow(i,j,k) / wombat%bac1_C2N) ))
          else
            wombat%bac1unh4(i,j,k) = 0.0
            wombat%bac1_lnit(i,j,k) = 1.0
          endif
        else ! No NH4 available to supplement growth
          wombat%bac1unh4(i,j,k) = 0.0
          wombat%bac1_lnit(i,j,k) = max(0.0, min(1.0, (wombat%don1remi(i,j,k) + wombat%bac1unh4(i,j,k)) &
                                                       / (wombat%bac1grow(i,j,k) / wombat%bac1_C2N) ))
        endif

        ! Determine degree of Fe limitation of bacteria and adjust growth rate
        if (biofer.gt.1e-3) then ! First, make sure that dFe is available for uptake
          ! Apply uptake kinetics constraint on bacterial dFe uptake
          wombat%bac1_lfer(i,j,k) = biofer / (biofer + wombat%bac1_kfer + epsi)
          wombat%bac1ufer(i,j,k) = wombat%bac1grow(i,j,k) / wombat%bac1_C2Fe * wombat%bac1_lfer(i,j,k) ! [molFe/kg/s]
          ! Check that enough dFe is available to support growth (at any one timestep, bacteria can only take up half of the available dFe)
          if (wombat%bac1ufer(i,j,k).gt.wombat%f_fe(i,j,k)*0.25/dt) then
            wombat%bac1ufer(i,j,k) = wombat%f_fe(i,j,k)*0.25 / dt ! [molFe/kg/s]
            wombat%bac1_lfer(i,j,k) = wombat%bac1ufer(i,j,k) / (wombat%bac1grow(i,j,k) / wombat%bac1_C2Fe)
          endif
        else ! No dFe available
          wombat%bac1_lfer(i,j,k) = 0.0
        endif
        
        ! Adjust the growth rate to be the minumum of N-limited or dFe-limited
        wombat%bac1_mu(i,j,k) = wombat%bac1_mu(i,j,k) * min(wombat%bac1_lfer(i,j,k), wombat%bac1_lnit(i,j,k))
        
        ! Final calculation after growth rate adjustment due to possible N or Fe limitation
        wombat%bac1grow(i,j,k) = wombat%bac1_mu(i,j,k) * wombat%f_bac1(i,j,k) ! [molC/kg/s]
        wombat%doc1remi(i,j,k) = wombat%bac1_mu(i,j,k) * wombat%f_bac1(i,j,k) * wombat%bac1_yaerC * (1. - wombat%bac1_fanaer(i,j,k)) & 
                                + wombat%bac1_mu(i,j,k) * wombat%f_bac1(i,j,k) * wombat%bac1_yanaC * wombat%bac1_fanaer(i,j,k) ! [molC/kg/s]
        wombat%don1remi(i,j,k) = wombat%doc1remi(i,j,k) * dom_N2C ! [molN/kg/s]
        wombat%bac1ufer(i,j,k) = wombat%bac1grow(i,j,k) / wombat%bac1_C2Fe ! [molFe/kg/s]
        wombat%bac1resp(i,j,k) = wombat%bac1_mu(i,j,k) * wombat%f_bac1(i,j,k) * wombat%bac1_yoxy * (1. - wombat%bac1_fanaer(i,j,k)) ! [molO2/kg/s]
        wombat%bac1deni(i,j,k) = wombat%bac1_mu(i,j,k) * wombat%f_bac1(i,j,k) * wombat%bac1_yno3 * wombat%bac1_fanaer(i,j,k) ! [molNO3/kg/s]
      endif

      ! Determine if bacteria are limited by N or Fe
      if (wombat%bac2_mu(i,j,k)*wombat%f_bac2(i,j,k).gt.0.0) then
        ! Initial estimate of the C biomass growth, DOC and DON assimilation rate by bacteria
        wombat%bac2grow(i,j,k) = wombat%bac2_mu(i,j,k) * wombat%f_bac2(i,j,k) ! [molC/kg/s]
        wombat%doc2remi(i,j,k) = wombat%bac2_mu(i,j,k) * wombat%f_bac2(i,j,k) * wombat%bac2_yaerC * (1. - wombat%bac2_fanaer(i,j,k)) & 
                                + wombat%bac2_mu(i,j,k) * wombat%f_bac2(i,j,k) * wombat%bac2_yanaC * wombat%bac2_fanaer(i,j,k) ! [molC/kg/s]
        wombat%don2remi(i,j,k) = wombat%doc2remi(i,j,k) * dom_N2C ! [molN/kg/s]

        ! Determine degree of N limitation of bacteria and adjust growth rate
        if (bionh4.gt.1e-3) then ! First, make sure that NH4 is available for uptake
          ! Find out what they need to supplement growth (this assumes that bacteria can achieve a N biomass yield of 1.0)
          wombat%bac2unh4(i,j,k) = wombat%bac2grow(i,j,k) / wombat%bac2_C2N - wombat%don2remi(i,j,k) ! [molN/kg/s]
          ! Find limitation of NH4 uptake
          if (wombat%bac2unh4(i,j,k).gt.0.0) then ! NH4 is needed to support growth
            ! Apply uptake kinetics constraint on bacterial NH4 uptake
            bac_limnh4 = bionh4 / (bionh4 + wombat%bac2_knh4 + epsi)
            wombat%bac2unh4(i,j,k) = wombat%bac2unh4(i,j,k) * bac_limnh4 ! [molN/kg/s]
            ! Make sure bacteria don't remove more NH4 than is available (can take up as much as half at a time)
            if (wombat%bac2unh4(i,j,k).gt.wombat%f_nh4(i,j,k)*0.25/dt) then
              wombat%bac2unh4(i,j,k) = wombat%f_nh4(i,j,k)*0.25 / dt ! [molN/kg/s]
            endif
            ! Recompute N limitation of bacteria (DON + NH4 uptake)
            wombat%bac2_lnit(i,j,k) = max(0.0, min(1.0, (wombat%don2remi(i,j,k) + wombat%bac2unh4(i,j,k)) &
                                                         / (wombat%bac2grow(i,j,k) / wombat%bac2_C2N) ))
          else
            wombat%bac2unh4(i,j,k) = 0.0
            wombat%bac2_lnit(i,j,k) = 1.0
          endif
        else ! No NH4 available to supplement growth
          wombat%bac2unh4(i,j,k) = 0.0
          wombat%bac2_lnit(i,j,k) = max(0.0, min(1.0, (wombat%don2remi(i,j,k) + wombat%bac2unh4(i,j,k)) &
                                                       / (wombat%bac2grow(i,j,k) / wombat%bac2_C2N) ))
        endif

        ! Determine degree of Fe limitation of bacteria and adjust growth rate
        if (biofer.gt.1e-3) then ! First, make sure that dFe is available for uptake
          ! Apply uptake kinetics constraint on bacterial dFe uptake
          wombat%bac2_lfer(i,j,k) = biofer / (biofer + wombat%bac2_kfer + epsi)
          wombat%bac2ufer(i,j,k) = wombat%bac2grow(i,j,k) / wombat%bac2_C2Fe * wombat%bac2_lfer(i,j,k) ! [molFe/kg/s]
          ! Check that enough dFe is available to support growth (at any one timestep, bacteria can only take up half of the available dFe)
          if (wombat%bac2ufer(i,j,k).gt.wombat%f_fe(i,j,k)*0.25/dt) then
            wombat%bac2ufer(i,j,k) = wombat%f_fe(i,j,k)*0.25 / dt ! [molFe/kg/s]
            wombat%bac2_lfer(i,j,k) = wombat%bac2ufer(i,j,k) / (wombat%bac2grow(i,j,k) / wombat%bac2_C2Fe)
          endif
        else ! No dFe available
          wombat%bac2_lfer(i,j,k) = 0.0
        endif
        
        ! Adjust the growth rate to be the minumum of N-limited or dFe-limited
        wombat%bac2_mu(i,j,k) = wombat%bac2_mu(i,j,k) * min(wombat%bac2_lfer(i,j,k), wombat%bac2_lnit(i,j,k))
        
        ! Final calculation after growth rate adjustment due to possible N or Fe limitation
        wombat%bac2grow(i,j,k) = wombat%bac2_mu(i,j,k) * wombat%f_bac2(i,j,k) ! [molC/kg/s]
        wombat%doc2remi(i,j,k) = wombat%bac2_mu(i,j,k) * wombat%f_bac2(i,j,k) * wombat%bac2_yaerC * (1. - wombat%bac2_fanaer(i,j,k)) & 
                                + wombat%bac2_mu(i,j,k) * wombat%f_bac2(i,j,k) * wombat%bac2_yanaC * wombat%bac2_fanaer(i,j,k) ! [molC/kg/s]
        wombat%don2remi(i,j,k) = wombat%doc2remi(i,j,k) * dom_N2C ! [molN/kg/s]
        wombat%bac2ufer(i,j,k) = wombat%bac2grow(i,j,k) / wombat%bac2_C2Fe ! [molFe/kg/s]
        wombat%bac2resp(i,j,k) = wombat%bac2_mu(i,j,k) * wombat%f_bac2(i,j,k) * wombat%bac2_yoxy * (1. - wombat%bac2_fanaer(i,j,k)) ! [molO2/kg/s]
        wombat%bac2deni(i,j,k) = wombat%bac2_mu(i,j,k) * wombat%f_bac2(i,j,k) * wombat%bac2_yn2o * wombat%bac2_fanaer(i,j,k) ! [molNO3/kg/s]
      endif

      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 13] Chemoautotroph calculations                                !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      
      ! 1. find max growth rate of AOA dependent on temperature (Qin et al., 2015 PNAS)
      wombat%aoa_mumax(i,j,k) = max(0.2, 0.029 * Temp(i,j,k) - 0.147) / 86400.0
      ! 2. Limitation terms of oxygen and ammonium substrate affecting uptake
      wombat%aoa_loxy(i,j,k) = min(1.0, biooxy * wombat%aoa_poxy)
      wombat%aoa_lnh4(i,j,k) = bionh4 / (bionh4 + wombat%aoa_knh4)
      aoa_Voxy = biooxy * wombat%aoa_poxy
      aoa_Vnh4 = wombat%aoa_ynh4 * wombat%aoa_mumax(i,j,k) * wombat%aoa_lnh4(i,j,k)  ! Note: yield * max growth rate = Vmax
      ! 3. Redefine growth rate based on these limitations
      wombat%aoa_mu(i,j,k) = min( (aoa_Voxy/wombat%aoa_yoxy), (aoa_Vnh4/wombat%aoa_ynh4) )

      ! 4. Determine N2O yield from ammonia oxidation 
      !    We use the empirical relationship with O2 from Frey et al. (2023) L&O; page 433 and 434
      !    They find a maximum yield of 3% per mol NO2 produced and a baseline yield of ~0.5% in
      !    oxic conditions (i.e., when O2 is not limiting), which we note here is in excess of the 
      !    baseline yields of other studies (Ji et al., 2018; Santoro et al., 2011; Qin et al., 2017)
      wombat%aoa_yn2o(i,j,k) = min(3.0, (0.2 / (biooxy + epsi) + 0.5)) * 0.01
      ! Because Frey give yield of N2O in % per mol NO2 produced, we must solve for mol N2O per mol biomass
      !  - aNH4 + bO2 --> cBiomass + dN2O + eNO3    |    and Y = N2O produced in % of NO3 produced 
      !  - d = (a - c) * Y / (2*Y + 1)
      wombat%aoa_yn2o(i,j,k) = (wombat%aoa_ynh4 - 1.0/wombat%aoa_C2N) * wombat%aoa_yn2o(i,j,k) &
                               / (2.0 * wombat%aoa_yn2o(i,j,k) + 1.0)
      
      if (do_anammox) then
        ! Anaerobic ammonium oxidation (anammox)
        wombat%aox_lnh4(i,j,k) = bionh4 / (bionh4 + wombat%aoxkn)
        wombat%aox_mu(i,j,k) = wombat%aoxmumax * wombat%bbioh**(Temp(i,j,k)) &
                               * wombat%bac1_fanaer(i,j,k) * wombat%aox_lnh4(i,j,k)
      endif
      
    
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 14] Sources and sinks                                          !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Phytoplankton growth
      if ((wombat%f_no3(i,j,k) + wombat%f_nh4(i,j,k)) .gt. epsi) then
        wombat%phygrow(i,j,k) = wombat%phy_mu(i,j,k) * wombat%f_phy(i,j,k) ! [molC/kg/s]
        wombat%diagrow(i,j,k) = wombat%dia_mu(i,j,k) * wombat%f_dia(i,j,k) ! [molC/kg/s]
      else
        wombat%phygrow(i,j,k) = 0.0
        wombat%diagrow(i,j,k) = 0.0
      endif

      ! Excess DOC exudation (active exudation via overflow hypothesis; Fogg 1966, 1983; Williams 1990; Carlson & Hansell 2014)
        ! Up to 50% (set by `overflow`) of assimilated carbon can be exuded by phytoplankton as DOC in high light, low nutrient conditions (Thornton 2014)
        ! Some small amount of DOC is exuded via passive diffusion even in the healthiest phytoplankton (Bjornsen 1988)
        ! If too much DOC is exuded, bacterial competition for nutrients can limit phytoplankton growth (Bratbak & Thingstad, 1985; Ratnarajah et al. 2021)
        ! However, active release of DOM by mixotrophic phytoplankton can "farm" heterotrophic bacteria (Mitra et al. 2013) (NOT YET IMPLEMENTED)
      if (wombat%f_phy(i,j,k).gt.epsi) then
        wombat%phydoc(i,j,k) = wombat%phygrow(i,j,k) * max(0.02, wombat%overflow & 
                               * ( wombat%phy_lpar(i,j,k) - min(wombat%phy_lfer(i,j,k), wombat%phy_lnit(i,j,k)) ) ) ! [molC/kg/s] 
      else
        wombat%phydoc(i,j,k) = 0.0
      endif
      if (wombat%f_dia(i,j,k).gt.epsi) then
        wombat%diadoc(i,j,k) = wombat%diagrow(i,j,k) * max(0.02, wombat%overflow &
                               * ( wombat%dia_lpar(i,j,k) - min(wombat%dia_lfer(i,j,k), wombat%dia_lnit(i,j,k)) ) ) ! [molC/kg/s] 
      else
        wombat%diadoc(i,j,k) = 0.0
      endif

      ! Chemoautotrophy
      if (wombat%f_aoa(i,j,k).gt.epsi) then
        wombat%aoagrow(i,j,k) = wombat%aoa_mu(i,j,k) * wombat%f_aoa(i,j,k) ! [molC/kg/s]
        wombat%ammox(i,j,k) = wombat%aoagrow(i,j,k) * wombat%aoa_ynh4 ! [molNH4/kg/s]
        wombat%aoaresp(i,j,k) = wombat%aoagrow(i,j,k) * wombat%aoa_yoxy ! [molO2/kg/s]
      else
        wombat%aoagrow(i,j,k) = 0.0
        wombat%ammox(i,j,k) = 0.0
        wombat%aoaresp(i,j,k) = 0.0
      endif
       
      if (wombat%f_nh4(i,j,k) .gt. epsi) then
        wombat%anammox(i,j,k) = wombat%aox_mu(i,j,k) * wombat%f_nh4(i,j,k) ! [molN/kg/s]
      else
        wombat%anammox(i,j,k) = 0.0
      endif

      ! remineralisation
      if (wombat%f_det(i,j,k) .gt. epsi) then
        wombat%detremi(i,j,k) = wombat%reminr(i,j,k) / mmol_m3_to_mol_kg * wombat%f_det(i,j,k)**2.0 ! [molC/kg/s]
      else
        wombat%detremi(i,j,k) = 0.0
      endif
      if (wombat%f_bdet(i,j,k) .gt. epsi) then
        wombat%bdetremi(i,j,k) = wombat%reminr(i,j,k) / mmol_m3_to_mol_kg * wombat%f_bdet(i,j,k)**2.0 ! [molC/kg/s]
      else
        wombat%bdetremi(i,j,k) = 0.0
      endif
      ! Grazing by microzooplankton
      if (Xzoo.gt.epsi) then
        wombat%zoograzbac1(i,j,k) = g_npz * wombat%f_zoo(i,j,k) * wombat%zooepsbac1*(wombat%zooprefbac1(i,j,k)*biobac1)**2 / Xzoo ! [molC/kg/s]
        wombat%zoograzbac2(i,j,k) = g_npz * wombat%f_zoo(i,j,k) * wombat%zooepsbac2*(wombat%zooprefbac2(i,j,k)*biobac2)**2 / Xzoo ! [molC/kg/s]
        wombat%zoograzaoa(i,j,k) = g_npz * wombat%f_zoo(i,j,k) * wombat%zooepsaoa*(wombat%zooprefaoa(i,j,k)*bioaoa)**2 / Xzoo ! [molC/kg/s]
        wombat%zoograzphy(i,j,k) = g_npz * wombat%f_zoo(i,j,k) * wombat%zooepsphy*(wombat%zooprefphy(i,j,k)*biophy)**2 / Xzoo ! [molC/kg/s]
        wombat%zoograzdia(i,j,k) = g_npz * wombat%f_zoo(i,j,k) * wombat%zooepsdia*(wombat%zooprefdia(i,j,k)*biodia)**2 / Xzoo ! [molC/kg/s]
        wombat%zoograzdet(i,j,k) = g_npz * wombat%f_zoo(i,j,k) * wombat%zooepsdet*(wombat%zooprefdet(i,j,k)*biodet)**2 / Xzoo ! [molC/kg/s]
      else
        wombat%zoograzbac1(i,j,k) = 0.0
        wombat%zoograzbac2(i,j,k) = 0.0
        wombat%zoograzaoa(i,j,k) = 0.0
        wombat%zoograzphy(i,j,k) = 0.0
        wombat%zoograzdia(i,j,k) = 0.0
        wombat%zoograzdet(i,j,k) = 0.0
      endif
      ! We follow Le Mezo & Galbraith (2021) L&O - The fecal iron pump: Global impact of animals on the iron stoichiometry...
      !  - ingestion, assimilation and excretion of carbon and iron by zooplankton are calculated separately
      !  - the idea is to enrich fecal pellets in iron compared to carbon
      wombat%zooexcrbac1(i,j,k) = wombat%zoograzbac1(i,j,k) * wombat%zooCingest*(1.0 - wombat%zooCassim)
      wombat%zooexcrbac2(i,j,k) = wombat%zoograzbac2(i,j,k) * wombat%zooCingest*(1.0 - wombat%zooCassim)
      wombat%zooexcraoa(i,j,k) = wombat%zoograzaoa(i,j,k) * wombat%zooCingest*(1.0 - wombat%zooCassim)
      wombat%zooexcrphy(i,j,k) = wombat%zoograzphy(i,j,k) * wombat%zooCingest*(1.0 - wombat%zooCassim)
      wombat%zooexcrdia(i,j,k) = wombat%zoograzdia(i,j,k) * wombat%zooCingest*(1.0 - wombat%zooCassim)
      wombat%zooexcrdet(i,j,k) = wombat%zoograzdet(i,j,k) * wombat%zooCingest*(1.0 - wombat%zooCassim)
      wombat%zooegesbac1(i,j,k) = wombat%zoograzbac1(i,j,k) * (1.0-wombat%zooCingest) 
      wombat%zooegesbac2(i,j,k) = wombat%zoograzbac2(i,j,k) * (1.0-wombat%zooCingest) 
      wombat%zooegesaoa(i,j,k) = wombat%zoograzaoa(i,j,k) * (1.0-wombat%zooCingest) 
      wombat%zooegesphy(i,j,k) = wombat%zoograzphy(i,j,k) * (1.0-wombat%zooCingest) 
      wombat%zooegesdia(i,j,k) = wombat%zoograzdia(i,j,k) * (1.0-wombat%zooCingest)
      wombat%zooegesdet(i,j,k) = wombat%zoograzdet(i,j,k) * (1.0-wombat%zooCingest)
      zooegesbac1fe = wombat%zoograzbac1(i,j,k) / wombat%bac1_C2Fe * (1.0-wombat%zooFeingest)
      zooegesbac2fe = wombat%zoograzbac2(i,j,k) / wombat%bac2_C2Fe * (1.0-wombat%zooFeingest)
      zooegesaoafe = wombat%zoograzaoa(i,j,k) / wombat%aoa_C2Fe * (1.0-wombat%zooFeingest)
      zooegesphyfe = wombat%zoograzphy(i,j,k) * phy_Fe2C * (1.0-wombat%zooFeingest)
      zooegesdiafe = wombat%zoograzdia(i,j,k) * dia_Fe2C * (1.0-wombat%zooFeingest)
      zooegesdetfe = wombat%zoograzdet(i,j,k) * det_Fe2C * (1.0-wombat%zooFeingest)
      zooassibac1fe = wombat%zoograzbac1(i,j,k) / wombat%bac1_C2Fe * wombat%zooFeingest*wombat%zooFeassim
      zooassibac2fe = wombat%zoograzbac2(i,j,k) / wombat%bac2_C2Fe * wombat%zooFeingest*wombat%zooFeassim
      zooassiaoafe = wombat%zoograzaoa(i,j,k) / wombat%aoa_C2Fe * wombat%zooFeingest*wombat%zooFeassim
      zooassiphyfe = wombat%zoograzphy(i,j,k) * phy_Fe2C * wombat%zooFeingest*wombat%zooFeassim
      zooassidiafe = wombat%zoograzdia(i,j,k) * dia_Fe2C * wombat%zooFeingest*wombat%zooFeassim
      zooassidetfe = wombat%zoograzdet(i,j,k) * det_Fe2C * wombat%zooFeingest*wombat%zooFeassim
      zooexcrbac1fe = wombat%zoograzbac1(i,j,k) / wombat%bac1_C2Fe * wombat%zooFeingest*(1.0 - wombat%zooFeassim)
      zooexcrbac2fe = wombat%zoograzbac2(i,j,k) / wombat%bac2_C2Fe * wombat%zooFeingest*(1.0 - wombat%zooFeassim)
      zooexcraoafe = wombat%zoograzaoa(i,j,k) / wombat%aoa_C2Fe * wombat%zooFeingest*(1.0 - wombat%zooFeassim)
      zooexcrphyfe = wombat%zoograzphy(i,j,k) * phy_Fe2C * wombat%zooFeingest*(1.0 - wombat%zooFeassim)
      zooexcrdiafe = wombat%zoograzdia(i,j,k) * dia_Fe2C * wombat%zooFeingest*(1.0 - wombat%zooFeassim)
      zooexcrdetfe = wombat%zoograzdet(i,j,k) * det_Fe2C * wombat%zooFeingest*(1.0 - wombat%zooFeassim)
      zooexcrbac1n  = wombat%zoograzbac1(i,j,k) / wombat%bac1_C2N &
                     - (wombat%zoograzbac1(i,j,k) * wombat%zooCingest * wombat%zooCassim / (122.0/16.0)) &
                     - (wombat%zooegesbac1(i,j,k) / (122.0/16.0)) ! [molN/kg/s]
      zooexcrbac2n  = wombat%zoograzbac2(i,j,k) / wombat%bac2_C2N &
                     - (wombat%zoograzbac2(i,j,k) * wombat%zooCingest * wombat%zooCassim / (122.0/16.0)) &
                     - (wombat%zooegesbac2(i,j,k) / (122.0/16.0)) ! [molN/kg/s]
      zooexcraoan  = wombat%zoograzaoa(i,j,k) / wombat%aoa_C2N &
                     - (wombat%zoograzaoa(i,j,k) * wombat%zooCingest * wombat%zooCassim / (122.0/16.0)) &
                     - (wombat%zooegesaoa(i,j,k) / (122.0/16.0)) ! [molN/kg/s]
      
      ! Grazing by mesozooplankton
      if (Xmes.gt.epsi) then
        wombat%mesgrazbac1(i,j,k) = m_npz * wombat%f_mes(i,j,k) * wombat%mesepsbac1*(wombat%mesprefbac1(i,j,k)*biobac1)**2 / Xmes ! [molC/kg/s]
        wombat%mesgrazbac2(i,j,k) = m_npz * wombat%f_mes(i,j,k) * wombat%mesepsbac2*(wombat%mesprefbac2(i,j,k)*biobac2)**2 / Xmes ! [molC/kg/s]
        wombat%mesgrazaoa(i,j,k) = m_npz * wombat%f_mes(i,j,k) * wombat%mesepsaoa*(wombat%mesprefaoa(i,j,k)*bioaoa)**2 / Xmes ! [molC/kg/s]
        wombat%mesgrazphy(i,j,k) = m_npz * wombat%f_mes(i,j,k) * wombat%mesepsphy*(wombat%mesprefphy(i,j,k)*biophy)**2 / Xmes ! [molC/kg/s]
        wombat%mesgrazdia(i,j,k) = m_npz * wombat%f_mes(i,j,k) * wombat%mesepsdia*(wombat%mesprefdia(i,j,k)*biodia)**2 / Xmes ! [molC/kg/s]
        wombat%mesgrazdet(i,j,k) = m_npz * wombat%f_mes(i,j,k) * wombat%mesepsdet*(wombat%mesprefdet(i,j,k)*biodet)**2 / Xmes ! [molC/kg/s]
        wombat%mesgrazbdet(i,j,k) = m_npz * wombat%f_mes(i,j,k) * wombat%mesepsbdet*(wombat%mesprefbdet(i,j,k)*biobdet)**2 / Xmes ! [molC/kg/s]
        wombat%mesgrazzoo(i,j,k) = m_npz * wombat%f_mes(i,j,k) * wombat%mesepszoo*(wombat%mesprefzoo(i,j,k)*biozoo)**2 / Xmes ! [molC/kg/s]
      else
        wombat%mesgrazbac1(i,j,k) = 0.0
        wombat%mesgrazbac2(i,j,k) = 0.0
        wombat%mesgrazaoa(i,j,k) = 0.0
        wombat%mesgrazphy(i,j,k) = 0.0
        wombat%mesgrazdia(i,j,k) = 0.0
        wombat%mesgrazdet(i,j,k) = 0.0
        wombat%mesgrazbdet(i,j,k) = 0.0
        wombat%mesgrazzoo(i,j,k) = 0.0
      endif
      ! We follow Le Mezo & Galbraith (2021) L&O - The fecal iron pump: Global impact of animals on the iron stoichiometry...
      !  - ingestion, assimilation and excretion of carbon and iron by zooplankton are calculated separately
      !  - the idea is to enrich fecal pellets in iron compared to carbon
      wombat%mesexcrbac1(i,j,k) = wombat%mesgrazbac1(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcrbac2(i,j,k) = wombat%mesgrazbac2(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcraoa(i,j,k) = wombat%mesgrazaoa(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcrphy(i,j,k) = wombat%mesgrazphy(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcrdia(i,j,k) = wombat%mesgrazdia(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcrdet(i,j,k) = wombat%mesgrazdet(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcrbdet(i,j,k) = wombat%mesgrazbdet(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesexcrzoo(i,j,k) = wombat%mesgrazzoo(i,j,k) * wombat%mesCingest*(1.0 - wombat%mesCassim)
      wombat%mesegesbac1(i,j,k) = wombat%mesgrazbac1(i,j,k) * (1.0 - wombat%mesCingest) 
      wombat%mesegesbac2(i,j,k) = wombat%mesgrazbac2(i,j,k) * (1.0 - wombat%mesCingest) 
      wombat%mesegesaoa(i,j,k) = wombat%mesgrazaoa(i,j,k) * (1.0 - wombat%mesCingest) 
      wombat%mesegesphy(i,j,k) = wombat%mesgrazphy(i,j,k) * (1.0 - wombat%mesCingest) 
      wombat%mesegesdia(i,j,k) = wombat%mesgrazdia(i,j,k) * (1.0 - wombat%mesCingest)
      wombat%mesegesdet(i,j,k) = wombat%mesgrazdet(i,j,k) * (1.0 - wombat%mesCingest)
      wombat%mesegesbdet(i,j,k) = wombat%mesgrazbdet(i,j,k) * (1.0 - wombat%mesCingest)
      wombat%mesegeszoo(i,j,k) = wombat%mesgrazzoo(i,j,k) * (1.0 - wombat%mesCingest)
      mesegesbac1fe = wombat%mesegesbac1(i,j,k) / wombat%bac1_C2Fe * (1.0-wombat%mesFeingest)
      mesegesbac2fe = wombat%mesegesbac2(i,j,k) / wombat%bac2_C2Fe * (1.0-wombat%mesFeingest)
      mesegesaoafe = wombat%mesegesaoa(i,j,k) / wombat%aoa_C2Fe * (1.0-wombat%mesFeingest)
      mesegesphyfe = wombat%mesegesphy(i,j,k) * phy_Fe2C * (1.0-wombat%mesFeingest)
      mesegesdiafe = wombat%mesegesdia(i,j,k) * dia_Fe2C * (1.0-wombat%mesFeingest)
      mesegesdetfe = wombat%mesegesdet(i,j,k) * det_Fe2C * (1.0-wombat%mesFeingest)
      mesegesbdetfe = wombat%mesegesbdet(i,j,k) * bdet_Fe2C * (1.0-wombat%mesFeingest)
      mesegeszoofe = wombat%mesegeszoo(i,j,k) * zoo_Fe2C * (1.0-wombat%mesFeingest)
      mesassibac1fe = wombat%mesgrazbac1(i,j,k) / wombat%bac1_C2Fe * wombat%mesFeingest*wombat%mesFeassim
      mesassibac2fe = wombat%mesgrazbac2(i,j,k) / wombat%bac2_C2Fe * wombat%mesFeingest*wombat%mesFeassim
      mesassiaoafe = wombat%mesgrazaoa(i,j,k) / wombat%aoa_C2Fe * wombat%mesFeingest*wombat%mesFeassim
      mesassiphyfe = wombat%mesgrazphy(i,j,k) * phy_Fe2C * wombat%mesFeingest*wombat%mesFeassim
      mesassidiafe = wombat%mesgrazdia(i,j,k) * dia_Fe2C * wombat%mesFeingest*wombat%mesFeassim
      mesassidetfe = wombat%mesgrazdet(i,j,k) * det_Fe2C * wombat%mesFeingest*wombat%mesFeassim
      mesassibdetfe = wombat%mesgrazbdet(i,j,k) * bdet_Fe2C * wombat%mesFeingest*wombat%mesFeassim
      mesassizoofe = wombat%mesgrazzoo(i,j,k) * zoo_Fe2C * wombat%mesFeingest*wombat%mesFeassim
      mesexcrbac1fe = wombat%mesgrazbac1(i,j,k) / wombat%bac1_C2Fe * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrbac2fe = wombat%mesgrazbac2(i,j,k) / wombat%bac2_C2Fe * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcraoafe = wombat%mesgrazaoa(i,j,k) / wombat%aoa_C2Fe * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrphyfe = wombat%mesgrazphy(i,j,k) * phy_Fe2C * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrdiafe = wombat%mesgrazdia(i,j,k) * dia_Fe2C * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrdetfe = wombat%mesgrazdet(i,j,k) * det_Fe2C * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrbdetfe = wombat%mesgrazbdet(i,j,k) * bdet_Fe2C * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrzoofe = wombat%mesgrazzoo(i,j,k) * zoo_Fe2C * wombat%mesFeingest*(1.0 - wombat%mesFeassim)
      mesexcrbac1n  = wombat%mesgrazbac1(i,j,k) / wombat%bac1_C2N &
                     - (wombat%mesgrazbac1(i,j,k) * wombat%mesCingest * wombat%mesCassim / (122.0/16.0)) &
                     - (wombat%mesegesbac1(i,j,k) / (122.0/16.0)) ! [molN/kg/s]
      mesexcrbac2n  = wombat%mesgrazbac2(i,j,k) / wombat%bac2_C2N &
                     - (wombat%mesgrazbac2(i,j,k) * wombat%mesCingest * wombat%mesCassim / (122.0/16.0)) &
                     - (wombat%mesegesbac2(i,j,k) / (122.0/16.0)) ! [molN/kg/s]
      mesexcraoan  = wombat%mesgrazaoa(i,j,k) / wombat%aoa_C2N &
                     - (wombat%mesgrazaoa(i,j,k) * wombat%mesCingest * wombat%mesCassim / (122.0/16.0)) &
                     - (wombat%mesegesaoa(i,j,k) / (122.0/16.0)) ! [molN/kg/s]

      ! Mortality terms
      if (biophy.gt.1e-3) then
        wombat%phyresp(i,j,k) = wombat%phylmor * fbc * wombat%f_phy(i,j,k) ! [molC/kg/s]
        wombat%phymort(i,j,k) = wombat%phyqmor / mmol_m3_to_mol_kg * wombat%f_phy(i,j,k) * wombat%f_phy(i,j,k) ! [molC/kg/s]
      else
        wombat%phyresp(i,j,k) = 0.0
        wombat%phymort(i,j,k) = 0.0
      endif
      if (biodia.gt.1e-3) then
        wombat%diaresp(i,j,k) = wombat%dialmor * fbc * wombat%f_dia(i,j,k) ! [molC/kg/s]
        wombat%diamort(i,j,k) = wombat%diaqmor / mmol_m3_to_mol_kg * wombat%f_dia(i,j,k) * wombat%f_dia(i,j,k) ! [molC/kg/s]
      else
        wombat%diaresp(i,j,k) = 0.0
        wombat%diamort(i,j,k) = 0.0
      endif
      if (biozoo.gt.1e-3) then
        wombat%zooresp(i,j,k) = wombat%zoolmor * fbc * wombat%f_zoo(i,j,k) * zoo_slmor ! [molC/kg/s]
        wombat%zoomort(i,j,k) = wombat%zooqmor / mmol_m3_to_mol_kg * wombat%f_zoo(i,j,k) * wombat%f_zoo(i,j,k) ! [molC/kg/s]
      else
        wombat%zooresp(i,j,k) = 0.0
        wombat%zoomort(i,j,k) = 0.0
      endif
      if (biomes.gt.1e-3) then
        wombat%mesresp(i,j,k) = wombat%meslmor * fbc * wombat%f_mes(i,j,k) * mes_slmor ! [molC/kg/s]
        wombat%mesmort(i,j,k) = wombat%mesqmor / mmol_m3_to_mol_kg * wombat%f_mes(i,j,k) * wombat%f_mes(i,j,k) ! [molC/kg/s]
      else
        wombat%mesresp(i,j,k) = 0.0
        wombat%mesmort(i,j,k) = 0.0
      endif
      if (biobac1.gt.1e-3) then
        wombat%bac1mor1(i,j,k) = wombat%bac1lmor * fbc * wombat%f_bac1(i,j,k) ! [molC/kg/s]
        wombat%bac1mor2(i,j,k) = wombat%bac1qmor / mmol_m3_to_mol_kg * wombat%f_bac1(i,j,k) * wombat%f_bac1(i,j,k) ! [molC/kg/s]
      else
        wombat%bac1mor1(i,j,k) = 0.0
        wombat%bac1mor2(i,j,k) = 0.0
      endif
      if (biobac2.gt.1e-3) then
        wombat%bac2mor1(i,j,k) = wombat%bac2lmor * fbc * wombat%f_bac2(i,j,k) ! [molC/kg/s]
        wombat%bac2mor2(i,j,k) = wombat%bac2qmor / mmol_m3_to_mol_kg * wombat%f_bac2(i,j,k) * wombat%f_bac2(i,j,k) ! [molC/kg/s]
      else
        wombat%bac2mor1(i,j,k) = 0.0
        wombat%bac2mor2(i,j,k) = 0.0
      endif
      if (bioaoa.gt.1e-3) then
        wombat%aoamor1(i,j,k) = wombat%aoalmor * fbc * wombat%f_aoa(i,j,k) ! [molC/kg/s]
        wombat%aoamor2(i,j,k) = wombat%aoaqmor / mmol_m3_to_mol_kg * wombat%f_aoa(i,j,k) * wombat%f_aoa(i,j,k) ! [molC/kg/s]
      else
        wombat%aoamor1(i,j,k) = 0.0
        wombat%aoamor2(i,j,k) = 0.0
      endif

      ! dissolution
      if (wombat%f_caco3(i,j,k) .gt. epsi) then
        wombat%caldiss(i,j,k) = wombat%dissrat(i,j,k) * wombat%f_caco3(i,j,k) ! [mol/kg/s]
      else
        wombat%caldiss(i,j,k) = 0.0
      endif
      ! dissolution
      if (wombat%f_bdetsi(i,j,k) .gt. epsi) then
        wombat%bsidiss(i,j,k) = wombat%disssi(i,j,k) * wombat%f_bdetsi(i,j,k) ! [mol/kg/s]
      else
        wombat%bsidiss(i,j,k) = 0.0
      endif


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 15] Tracer tendencies                                          !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Nitrate equation ! [molN/kg]
      !----------------------------------------------------------------------
      wombat%f_no3(i,j,k) = wombat%f_no3(i,j,k) + dtsb * ( 0.0 & 
                            + (wombat%ammox(i,j,k) - wombat%aoagrow(i,j,k) * (1.0/wombat%aoa_C2N + 2*wombat%aoa_yn2o(i,j,k))) &
                            - wombat%bac1deni(i,j,k) ) &
                            + dtsb * 16./122. * ( 0.0 &
                            - wombat%phygrow(i,j,k) * wombat%phy_lno3(i,j,k) / ( wombat%phy_lnit(i,j,k) + epsi ) &
                            - wombat%diagrow(i,j,k) * wombat%dia_lno3(i,j,k) / ( wombat%dia_lnit(i,j,k) + epsi ) )
    
      ! Ammonium equation ! [molN/kg]
      !----------------------------------------------------------------------
      wombat%f_nh4(i,j,k) = wombat%f_nh4(i,j,k) + dtsb * ( 0.0 &
                            + zooexcrbac1n*(1.0-wombat%zooexcrdom) &
                            + mesexcrbac1n*(1.0-wombat%mesexcrdom) &
                            + zooexcrbac2n*(1.0-wombat%zooexcrdom) &
                            + mesexcrbac2n*(1.0-wombat%mesexcrdom) &
                            + zooexcraoan*(1.0-wombat%zooexcrdom) &
                            + mesexcraoan*(1.0-wombat%mesexcrdom) &
                            + wombat%nitrfix(i,j,k) &
                            + (wombat%don1remi(i,j,k) - wombat%bac1grow(i,j,k)/wombat%bac1_C2N) &
                            + (wombat%don2remi(i,j,k) - wombat%bac2grow(i,j,k)/wombat%bac2_C2N) &
                            - wombat%ammox(i,j,k) &
                            - wombat%anammox(i,j,k) ) &
                            + dtsb * 16./122. * ( 0.0 &
                            + wombat%zooresp(i,j,k) &
                            + wombat%mesresp(i,j,k) &
                            + wombat%zooexcrphy(i,j,k)*(1.0-wombat%zooexcrdom) &
                            + wombat%zooexcrdia(i,j,k)*(1.0-wombat%zooexcrdom) &
                            + wombat%zooexcrdet(i,j,k)*(1.0-wombat%zooexcrdom) &
                            + wombat%mesexcrphy(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrdia(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrdet(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrbdet(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrzoo(i,j,k)*(1.0-wombat%mesexcrdom) &
                            - wombat%phygrow(i,j,k) * wombat%phy_lnh4(i,j,k) / ( wombat%phy_lnit(i,j,k) + epsi ) &
                            - wombat%diagrow(i,j,k) * wombat%dia_lnh4(i,j,k) / ( wombat%dia_lnit(i,j,k) + epsi ) )
    
      ! Silicic acid equation ! [molSi/kg] 
      !   Microzooplankton grazing on diatoms produces clean, small, largely suspended 
      !   pieces of biogenic silica prone to rapid dissolution [Krause et al., 2010 L&O]
      !----------------------------------------------------------------------
      wombat%f_sil(i,j,k) = wombat%f_sil(i,j,k) + dtsb * ( 0.0 &
                            - wombat%diagrow(i,j,k) * 16.0/122.0 &
                            + wombat%diaresp(i,j,k) * dia_Si2C &
                            + wombat%zoograzdia(i,j,k) * dia_Si2C &
                            + wombat%bsidiss(i,j,k) )
                               
                            
      ! Nitrous oxide equation ! [molN2/kg] 
      !  pjb: note that we track N2O in units of mol N2/kg, accounting for the two N atoms
      !----------------------------------------------------------------------
      wombat%f_n2o(i,j,k) = wombat%f_n2o(i,j,k) + dtsb * ( 0.0 &
                            + wombat%aoagrow(i,j,k) * wombat%aoa_yn2o(i,j,k) &
                            + wombat%bac1deni(i,j,k)/2.0 &
                            - wombat%bac2deni(i,j,k) ) 
                              
      ! Phytoplankton equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_phy(i,j,k)  = wombat%f_phy(i,j,k) + dtsb * ( 0.0 &
                             + wombat%phygrow(i,j,k) &
                             - wombat%phyresp(i,j,k) &
                             - wombat%phymort(i,j,k) &
                             - wombat%zoograzphy(i,j,k) &
                             - wombat%mesgrazphy(i,j,k) )
      
      ! Phytoplankton chlorophyll equation ! [molChl/kg]
      !-----------------------------------------------------------------------
      wombat%f_pchl(i,j,k)  = wombat%f_pchl(i,j,k) + dtsb * ( 0.0 &
                              + wombat%pchl_mu(i,j,k) &
                              - wombat%phyresp(i,j,k) * phy_chlc &
                              - wombat%phymort(i,j,k) * phy_chlc &
                              - wombat%zoograzphy(i,j,k) * phy_chlc &
                              - wombat%mesgrazphy(i,j,k) * phy_chlc )

      ! Phytoplankton iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_phyfe(i,j,k)  = wombat%f_phyfe(i,j,k) + dtsb * ( 0.0 & 
                               + wombat%phy_dfeupt(i,j,k) &
                               - wombat%phyresp(i,j,k) * phy_Fe2C &
                               - wombat%phymort(i,j,k) * phy_Fe2C &
                               - wombat%zoograzphy(i,j,k) * phy_Fe2C &
                               - wombat%mesgrazphy(i,j,k) * phy_Fe2C )

      ! Microphytoplankton equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_dia(i,j,k)  = wombat%f_dia(i,j,k) + dtsb * ( 0.0 &
                             + wombat%diagrow(i,j,k) &
                             - wombat%diaresp(i,j,k) &
                             - wombat%diamort(i,j,k) &
                             - wombat%mesgrazdia(i,j,k) &
                             - wombat%zoograzdia(i,j,k) )
      
      ! Microphytoplankton chlorodiall equation ! [molChl/kg]
      !-----------------------------------------------------------------------
      wombat%f_dchl(i,j,k)  = wombat%f_dchl(i,j,k) + dtsb * ( 0.0 &
                              + wombat%dchl_mu(i,j,k) &
                              - wombat%diaresp(i,j,k) * dia_chlc &
                              - wombat%diamort(i,j,k) * dia_chlc &
                              - wombat%zoograzdia(i,j,k) * dia_chlc &
                              - wombat%mesgrazdia(i,j,k) * dia_chlc )

      ! Microphytoplankton iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_diafe(i,j,k)  = wombat%f_diafe(i,j,k) + dtsb * ( 0.0 & 
                               + wombat%dia_dfeupt(i,j,k) &
                               - wombat%diaresp(i,j,k) * dia_Fe2C &
                               - wombat%diamort(i,j,k) * dia_Fe2C &
                               - wombat%zoograzdia(i,j,k) * dia_Fe2C &
                               - wombat%mesgrazdia(i,j,k) * dia_Fe2C )

      ! Microphytoplankton silicon equation ! [molSi/kg] 
      !----------------------------------------------------------------------
      wombat%f_diasi(i,j,k) = wombat%f_diasi(i,j,k) + dtsb * ( 0.0 &
                              + wombat%diagrow(i,j,k) * 16.0/122.0 &
                              - wombat%diaresp(i,j,k) * dia_Si2C &
                              - wombat%diamort(i,j,k) * dia_Si2C &
                              - wombat%mesgrazdia(i,j,k) * dia_Si2C &
                              - wombat%zoograzdia(i,j,k) * dia_Si2C )
      
      ! Estimate primary productivity from phytoplankton growth ! [molC/kg/s]
      wombat%pprod_gross(i,j,k) = wombat%pprod_gross(i,j,k) + dtsb * ( 0.0 &
                                  + wombat%phygrow(i,j,k) + wombat%diagrow(i,j,k) )

      ! Net primary productivity (gross PP minus linear mortality) ! [molC/kg/s]
      wombat%npp3d(i,j,k) = wombat%npp3d(i,j,k) + dtsb * ( 0.0 &
                            + wombat%phygrow(i,j,k) - wombat%phyresp(i,j,k) &
                            + wombat%diagrow(i,j,k) - wombat%diaresp(i,j,k) )

      ! Zooplankton equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_zoo(i,j,k)  = wombat%f_zoo(i,j,k) + dtsb * ( 0.0 &
                             + wombat%zooCingest*wombat%zooCassim * wombat%zoograzbac1(i,j,k) &
                             + wombat%zooCingest*wombat%zooCassim * wombat%zoograzbac2(i,j,k) &
                             + wombat%zooCingest*wombat%zooCassim * wombat%zoograzaoa(i,j,k) &
                             + wombat%zooCingest*wombat%zooCassim * wombat%zoograzphy(i,j,k) &
                             + wombat%zooCingest*wombat%zooCassim * wombat%zoograzdia(i,j,k) &
                             + wombat%zooCingest*wombat%zooCassim * wombat%zoograzdet(i,j,k) &
                             - wombat%mesgrazzoo(i,j,k) &
                             - wombat%zooresp(i,j,k) &
                             - wombat%zoomort(i,j,k) )

      ! Zooplankton iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_zoofe(i,j,k)  = wombat%f_zoofe(i,j,k) + dtsb * ( 0.0 &
                               + zooassibac1fe &
                               + zooassibac2fe &
                               + zooassiaoafe &
                               + zooassiphyfe &
                               + zooassidiafe &
                               + zooassidetfe &
                               - wombat%mesgrazzoo(i,j,k) * zoo_Fe2C &
                               - wombat%zooresp(i,j,k) * zoo_Fe2C &
                               - wombat%zoomort(i,j,k) * zoo_Fe2C )

      ! Mesozooplankton equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_mes(i,j,k)  = wombat%f_mes(i,j,k) + dtsb * ( 0.0 &
                             + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazbac1(i,j,k) &
                             + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazbac2(i,j,k) &
                             + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazaoa(i,j,k) &
                             + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazphy(i,j,k) &
                             + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazdia(i,j,k) &
                             + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazdet(i,j,k) &
                             + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazbdet(i,j,k) &
                             + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazzoo(i,j,k) &
                             - wombat%mesresp(i,j,k) &
                             - wombat%mesmort(i,j,k) )

      ! Mesozooplankton iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_mesfe(i,j,k)  = wombat%f_mesfe(i,j,k) + dtsb * ( 0.0 &
                               + mesassibac1fe &
                               + mesassibac2fe &
                               + mesassiaoafe &
                               + mesassiphyfe &
                               + mesassidiafe &
                               + mesassidetfe &
                               + mesassibdetfe &
                               + mesassizoofe &
                               - wombat%mesresp(i,j,k) * mes_Fe2C &
                               - wombat%mesmort(i,j,k) * mes_Fe2C )

      ! Estimate secondary productivity from zooplankton growth ! [molC/kg/s]
      wombat%zprod_gross(i,j,k) = wombat%zprod_gross(i,j,k) + dtsb * ( 0.0 &
                                  + wombat%zooCingest*wombat%zooCassim * wombat%zoograzbac1(i,j,k) &
                                  + wombat%zooCingest*wombat%zooCassim * wombat%zoograzbac2(i,j,k) &
                                  + wombat%zooCingest*wombat%zooCassim * wombat%zoograzaoa(i,j,k) &
                                  + wombat%zooCingest*wombat%zooCassim * wombat%zoograzphy(i,j,k) &
                                  + wombat%zooCingest*wombat%zooCassim * wombat%zoograzdia(i,j,k) &
                                  + wombat%zooCingest*wombat%zooCassim * wombat%zoograzdet(i,j,k) &
                                  + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazbac1(i,j,k) &
                                  + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazbac2(i,j,k) &
                                  + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazaoa(i,j,k) &
                                  + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazphy(i,j,k) &
                                  + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazdia(i,j,k) &
                                  + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazdet(i,j,k) &
                                  + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazbdet(i,j,k) &
                                  + wombat%mesCingest*wombat%mesCassim * wombat%mesgrazzoo(i,j,k) )

      ! Detritus equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_det(i,j,k) = wombat%f_det(i,j,k) + dtsb * ( 0.0 & 
                            + wombat%zooegesbac1(i,j,k) &
                            + wombat%zooegesbac2(i,j,k) &
                            + wombat%zooegesaoa(i,j,k) &
                            + wombat%zooegesphy(i,j,k) &
                            + wombat%zooegesdia(i,j,k) &
                            + wombat%zooegesdet(i,j,k) &
                            + wombat%phymort(i,j,k) &
                            + wombat%zoomort(i,j,k) &
                            - wombat%zoograzdet(i,j,k) &
                            - wombat%mesgrazdet(i,j,k) &
                            - wombat%detremi(i,j,k) )

      ! Detrital iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_detfe(i,j,k) = wombat%f_detfe(i,j,k) + dtsb * ( 0.0 & 
                              + zooegesbac1fe &
                              + zooegesbac2fe &
                              + zooegesaoafe &
                              + zooegesphyfe &
                              + zooegesdiafe &
                              + zooegesdetfe &
                              + wombat%phymort(i,j,k) * phy_Fe2C &
                              + wombat%zoomort(i,j,k) * zoo_Fe2C &
                              - wombat%zoograzdet(i,j,k) * det_Fe2C &
                              - wombat%mesgrazdet(i,j,k) * det_Fe2C &
                              - wombat%detremi(i,j,k) * det_Fe2C &
                              + wombat%fescadet(i,j,k) &
                              + wombat%fecoag2det(i,j,k) )

      ! Big detritus equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_bdet(i,j,k) = wombat%f_bdet(i,j,k) + dtsb * ( 0.0 &
                             + wombat%diamort(i,j,k) &
                             + wombat%mesmort(i,j,k) &
                             + wombat%mesegesbac1(i,j,k) &
                             + wombat%mesegesbac2(i,j,k) &
                             + wombat%mesegesaoa(i,j,k) &
                             + wombat%mesegesphy(i,j,k) &
                             + wombat%mesegesdia(i,j,k) &
                             + wombat%mesegesdet(i,j,k) &
                             + wombat%mesegesbdet(i,j,k) &
                             + wombat%mesegeszoo(i,j,k) &
                             - wombat%mesgrazbdet(i,j,k) &
                             - wombat%bdetremi(i,j,k) )

      ! Big detrital iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_bdetfe(i,j,k) = wombat%f_bdetfe(i,j,k) + dtsb * ( 0.0 &
                               + wombat%diamort(i,j,k) * dia_Fe2C &
                               + wombat%mesmort(i,j,k) * mes_Fe2C &
                               + mesegesbac1fe &
                               + mesegesbac2fe &
                               + mesegesaoafe &
                               + mesegesphyfe &
                               + mesegesdiafe &
                               + mesegesdetfe &
                               + mesegesbdetfe &
                               + mesegeszoofe &
                               - wombat%mesgrazbdet(i,j,k) * bdet_Fe2C &
                               - wombat%bdetremi(i,j,k) * bdet_Fe2C &
                               + wombat%fescabdet(i,j,k) &
                               + wombat%fecoag2bdet(i,j,k) )
      
      ! Microphytoplankton silicon equation ! [molSi/kg]
      !   Copepod egestion (fecal pellets) represented 42-107% of biogenic silica export at
      !   100 metres in the spring bloom at the Antarctic Polar Front [Dagg et al., 2003 DSRII]
      !   So all mesozooplankton grazing on diatoms goes to egestion, no dissolution
      !----------------------------------------------------------------------
      wombat%f_bdetsi(i,j,k) = wombat%f_bdetsi(i,j,k) + dtsb * ( 0.0 &
                               + wombat%diamort(i,j,k) * dia_Si2C &
                               + wombat%mesgrazdia(i,j,k) * dia_Si2C &
                               - wombat%bsidiss(i,j,k)  )
      
      ! Dissolved organic carbon equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_doc(i,j,k) = wombat%f_doc(i,j,k) + dtsb * ( 0.0 &
                            + wombat%phydoc(i,j,k) &
                            + wombat%diadoc(i,j,k) &
                            + wombat%detremi(i,j,k) &
                            + wombat%bdetremi(i,j,k) &
                            + wombat%phyresp(i,j,k) &
                            + wombat%diaresp(i,j,k) &
                            + wombat%zooexcrbac1(i,j,k)*wombat%zooexcrdom &
                            + wombat%zooexcrbac2(i,j,k)*wombat%zooexcrdom &
                            + wombat%zooexcraoa(i,j,k)*wombat%zooexcrdom &
                            + wombat%zooexcrphy(i,j,k)*wombat%zooexcrdom &
                            + wombat%zooexcrdia(i,j,k)*wombat%zooexcrdom &
                            + wombat%zooexcrdet(i,j,k)*wombat%zooexcrdom &
                            + wombat%mesexcrbac1(i,j,k)*wombat%mesexcrdom &
                            + wombat%mesexcrbac2(i,j,k)*wombat%mesexcrdom &
                            + wombat%mesexcraoa(i,j,k)*wombat%mesexcrdom &
                            + wombat%mesexcrphy(i,j,k)*wombat%mesexcrdom &
                            + wombat%mesexcrdia(i,j,k)*wombat%mesexcrdom &
                            + wombat%mesexcrdet(i,j,k)*wombat%mesexcrdom &
                            + wombat%mesexcrbdet(i,j,k)*wombat%mesexcrdom &
                            + wombat%mesexcrzoo(i,j,k)*wombat%mesexcrdom &
                            - wombat%doc1remi(i,j,k) &
                            - wombat%doc2remi(i,j,k) &
                            + wombat%bac1mor1(i,j,k) &
                            + wombat%bac1mor2(i,j,k) &
                            + wombat%bac2mor1(i,j,k) &
                            + wombat%bac2mor2(i,j,k) &
                            + wombat%aoamor1(i,j,k) &
                            + wombat%aoamor2(i,j,k) )

      ! Dissolved organic nitrogen equation ! [molN/kg]
      !-----------------------------------------------------------------------
      wombat%f_don(i,j,k) = wombat%f_don(i,j,k) + dtsb * ( 0.0 &
                            - wombat%don1remi(i,j,k) &
                            - wombat%don2remi(i,j,k) &
                            + zooexcrbac1n*wombat%zooexcrdom &
                            + mesexcrbac1n*wombat%mesexcrdom &
                            + zooexcrbac2n*wombat%zooexcrdom &
                            + mesexcrbac2n*wombat%mesexcrdom &
                            + zooexcraoan*wombat%zooexcrdom &
                            + mesexcraoan*wombat%mesexcrdom &
                            + wombat%bac1mor1(i,j,k) / wombat%bac1_C2N &
                            + wombat%bac1mor2(i,j,k) / wombat%bac1_C2N &
                            + wombat%bac2mor1(i,j,k) / wombat%bac2_C2N &
                            + wombat%bac2mor2(i,j,k) / wombat%bac2_C2N &
                            + wombat%aoamor1(i,j,k) / wombat%aoa_C2N &
                            + wombat%aoamor2(i,j,k) / wombat%aoa_C2N ) &
                            + dtsb * 16./122.0 * ( 0.0 & 
                            + wombat%detremi(i,j,k) &
                            + wombat%bdetremi(i,j,k) &
                            + wombat%phyresp(i,j,k) &
                            + wombat%diaresp(i,j,k) &
                            + wombat%zooexcrphy(i,j,k)*wombat%zooexcrdom &
                            + wombat%zooexcrdia(i,j,k)*wombat%zooexcrdom &
                            + wombat%zooexcrdet(i,j,k)*wombat%zooexcrdom &
                            + wombat%mesexcrphy(i,j,k)*wombat%mesexcrdom &
                            + wombat%mesexcrdia(i,j,k)*wombat%mesexcrdom &
                            + wombat%mesexcrdet(i,j,k)*wombat%mesexcrdom &
                            + wombat%mesexcrbdet(i,j,k)*wombat%mesexcrdom &
                            + wombat%mesexcrzoo(i,j,k)*wombat%mesexcrdom )
    
      ! Heterotrophic bacteria #1 ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_bac1(i,j,k) = wombat%f_bac1(i,j,k) + dtsb * ( 0.0 &
                            + wombat%bac1grow(i,j,k) &
                            - wombat%zoograzbac1(i,j,k) &
                            - wombat%mesgrazbac1(i,j,k) &
                            - wombat%bac1mor1(i,j,k) &
                            - wombat%bac1mor2(i,j,k) )
    
      ! Heterotrophic bacteria #2 ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_bac2(i,j,k) = wombat%f_bac2(i,j,k) + dtsb * ( 0.0 &
                            + wombat%bac2grow(i,j,k) &
                            - wombat%zoograzbac2(i,j,k) &
                            - wombat%mesgrazbac2(i,j,k) &
                            - wombat%bac2mor1(i,j,k) &
                            - wombat%bac2mor2(i,j,k) )
    
      ! AOA ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_aoa(i,j,k) = wombat%f_aoa(i,j,k) + dtsb * ( 0.0 &
                            + wombat%aoagrow(i,j,k) &
                            - wombat%zoograzaoa(i,j,k) &
                            - wombat%mesgrazaoa(i,j,k) &
                            - wombat%aoamor1(i,j,k) &
                            - wombat%aoamor2(i,j,k) )
      
      ! Oxygen equation ! [molO2/kg]
      !-----------------------------------------------------------------------
      if (wombat%f_o2(i,j,k) .gt. epsi) &
        wombat%f_o2(i,j,k) = wombat%f_o2(i,j,k) - 132./122. * dtsb * ( 0.0 &
                             + wombat%zooresp(i,j,k) &
                             + wombat%mesresp(i,j,k) &
                             + wombat%zooexcrbac1(i,j,k)*(1.0-wombat%zooexcrdom) &
                             + wombat%zooexcrbac2(i,j,k)*(1.0-wombat%zooexcrdom) &
                             + wombat%zooexcraoa(i,j,k)*(1.0-wombat%zooexcrdom) &
                             + wombat%zooexcrphy(i,j,k)*(1.0-wombat%zooexcrdom) &
                             + wombat%zooexcrdia(i,j,k)*(1.0-wombat%zooexcrdom) &
                             + wombat%zooexcrdet(i,j,k)*(1.0-wombat%zooexcrdom) &
                             + wombat%mesexcrbac1(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcrbac2(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcraoa(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcrphy(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcrdia(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcrdet(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcrbdet(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcrzoo(i,j,k)*(1.0-wombat%mesexcrdom) &
                             - wombat%phygrow(i,j,k) &
                             - wombat%diagrow(i,j,k) ) &
                             - dtsb * ( 0.0 &
                             + wombat%bac1resp(i,j,k) &
                             + wombat%bac2resp(i,j,k) &
                             + wombat%aoaresp(i,j,k) )


      ! Equation for CaCO3 ! [molCaCO3/kg]
      !-----------------------------------------------------------------------
      wombat%f_caco3(i,j,k) = wombat%f_caco3(i,j,k) + dtsb * ( 0.0 &
                              + wombat%zooegesphy(i,j,k) * wombat%pic2poc(i,j,k) &
                              + wombat%mesegesphy(i,j,k) * wombat%pic2poc(i,j,k) &
                              + wombat%phymort(i,j,k) * wombat%pic2poc(i,j,k) &
                              + wombat%zoomort(i,j,k) * wombat%pic2poc(i,j,k) &
                              + wombat%mesmort(i,j,k) * wombat%pic2poc(i,j,k) &
                              - wombat%caldiss(i,j,k) )

      ! Equation for DIC ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_dic(i,j,k) = wombat%f_dic(i,j,k) + dtsb * ( 0.0 &
                            + (wombat%doc1remi(i,j,k) - wombat%bac1grow(i,j,k)) &
                            + (wombat%doc2remi(i,j,k) - wombat%bac2grow(i,j,k)) &
                            + wombat%zooresp(i,j,k) &
                            + wombat%mesresp(i,j,k) &
                            + wombat%zooexcrbac1(i,j,k)*(1.0-wombat%zooexcrdom) &
                            + wombat%zooexcrbac2(i,j,k)*(1.0-wombat%zooexcrdom) &
                            + wombat%zooexcraoa(i,j,k)*(1.0-wombat%zooexcrdom) &
                            + wombat%zooexcrphy(i,j,k)*(1.0-wombat%zooexcrdom) &
                            + wombat%zooexcrdia(i,j,k)*(1.0-wombat%zooexcrdom) &
                            + wombat%zooexcrdet(i,j,k)*(1.0-wombat%zooexcrdom) &
                            + wombat%mesexcrbac1(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrbac2(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcraoa(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrphy(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrdia(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrdet(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrbdet(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrzoo(i,j,k)*(1.0-wombat%mesexcrdom) &
                            - wombat%phygrow(i,j,k) &
                            - wombat%diagrow(i,j,k) &
                            - wombat%aoagrow(i,j,k) &
                            - wombat%phydoc(i,j,k) &
                            - wombat%diadoc(i,j,k) &
                            - wombat%zooegesphy(i,j,k) * wombat%pic2poc(i,j,k) &
                            - wombat%mesegesphy(i,j,k) * wombat%pic2poc(i,j,k) &
                            - wombat%phymort(i,j,k) * wombat%pic2poc(i,j,k) &
                            - wombat%zoomort(i,j,k) * wombat%pic2poc(i,j,k) &  
                            - wombat%mesmort(i,j,k) * wombat%pic2poc(i,j,k) &  
                            + wombat%caldiss(i,j,k) )

      ! Equation for DICr ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_dicr(i,j,k) = wombat%f_dicr(i,j,k) + dtsb * ( 0.0 &
                             + (wombat%doc1remi(i,j,k) - wombat%bac1grow(i,j,k)) &
                             + (wombat%doc2remi(i,j,k) - wombat%bac2grow(i,j,k)) &
                             + wombat%zooresp(i,j,k) &
                             + wombat%mesresp(i,j,k) &
                             + wombat%zooexcrbac1(i,j,k)*(1.0-wombat%zooexcrdom) &
                             + wombat%zooexcrbac2(i,j,k)*(1.0-wombat%zooexcrdom) &
                             + wombat%zooexcraoa(i,j,k)*(1.0-wombat%zooexcrdom) &
                             + wombat%zooexcrphy(i,j,k)*(1.0-wombat%zooexcrdom) &
                             + wombat%zooexcrdia(i,j,k)*(1.0-wombat%zooexcrdom) &
                             + wombat%zooexcrdet(i,j,k)*(1.0-wombat%zooexcrdom) &
                             + wombat%mesexcrbac1(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcrbac2(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcraoa(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcrphy(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcrdia(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcrdet(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcrbdet(i,j,k)*(1.0-wombat%mesexcrdom) &
                             + wombat%mesexcrzoo(i,j,k)*(1.0-wombat%mesexcrdom) &
                             - wombat%phygrow(i,j,k) &
                             - wombat%diagrow(i,j,k) &
                             - wombat%aoagrow(i,j,k) &
                             - wombat%phydoc(i,j,k) &
                             - wombat%diadoc(i,j,k) &
                             - wombat%zooegesphy(i,j,k) * wombat%pic2poc(i,j,k) &
                             - wombat%mesegesphy(i,j,k) * wombat%pic2poc(i,j,k) &
                             - wombat%phymort(i,j,k) * wombat%pic2poc(i,j,k) &
                             - wombat%zoomort(i,j,k) * wombat%pic2poc(i,j,k) &  
                             - wombat%mesmort(i,j,k) * wombat%pic2poc(i,j,k) &  
                             + wombat%caldiss(i,j,k) )

      ! Equation for ALK ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_alk(i,j,k) = wombat%f_alk(i,j,k) + dtsb * 16.0/122.0 * ( 0.0 &
                            + wombat%phygrow(i,j,k) * wombat%phy_lno3(i,j,k) / ( wombat%phy_lnit(i,j,k) + epsi ) &
                            + wombat%diagrow(i,j,k) * wombat%dia_lno3(i,j,k) / ( wombat%dia_lnit(i,j,k) + epsi ) &
                            + wombat%zooresp(i,j,k) &
                            + wombat%mesresp(i,j,k) &
                            + wombat%zooexcrphy(i,j,k)*(1.0-wombat%zooexcrdom) &
                            + wombat%zooexcrdia(i,j,k)*(1.0-wombat%zooexcrdom) &
                            + wombat%zooexcrdet(i,j,k)*(1.0-wombat%zooexcrdom) &
                            + wombat%mesexcrphy(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrdia(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrdet(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrbdet(i,j,k)*(1.0-wombat%mesexcrdom) &
                            + wombat%mesexcrzoo(i,j,k)*(1.0-wombat%mesexcrdom) &
                            - wombat%phygrow(i,j,k) * wombat%phy_lnh4(i,j,k) / ( wombat%phy_lnit(i,j,k) + epsi ) &
                            - wombat%diagrow(i,j,k) * wombat%dia_lnh4(i,j,k) / ( wombat%dia_lnit(i,j,k) + epsi ) ) &
                            + dtsb * ( 0.0 &
                            + (wombat%don1remi(i,j,k) - wombat%bac1grow(i,j,k)/wombat%bac1_C2N) &
                            + (wombat%don2remi(i,j,k) - wombat%bac2grow(i,j,k)/wombat%bac2_C2N) &
                            + zooexcrbac1n*(1.0-wombat%zooexcrdom) &
                            + mesexcrbac1n*(1.0-wombat%mesexcrdom) &
                            + zooexcrbac2n*(1.0-wombat%zooexcrdom) &
                            + mesexcrbac2n*(1.0-wombat%mesexcrdom) &
                            + zooexcraoan*(1.0-wombat%zooexcrdom) &
                            + mesexcraoan*(1.0-wombat%mesexcrdom) &
                            + wombat%bac1deni(i,j,k) &
                            - wombat%anammox(i,j,k) &
                            - wombat%ammox(i,j,k) &
                            - (wombat%ammox(i,j,k) - wombat%aoagrow(i,j,k)/wombat%aoa_C2N) ) &
                            + dtsb * 2.0 * ( 0.0 &
                            + wombat%caldiss(i,j,k) &
                            - wombat%zooegesphy(i,j,k) * wombat%pic2poc(i,j,k) &
                            - wombat%mesegesphy(i,j,k) * wombat%pic2poc(i,j,k) &
                            - wombat%phymort(i,j,k) * wombat%pic2poc(i,j,k) &
                            - wombat%zoomort(i,j,k) * wombat%pic2poc(i,j,k) &
                            - wombat%mesmort(i,j,k) * wombat%pic2poc(i,j,k) )

      ! Extra equation for iron ! [molFe/kg]
      !----------------------------------------------------------------------
      wombat%f_fe(i,j,k) = wombat%f_fe(i,j,k) + dtsb * ( 0.0 &
                           + wombat%detremi(i,j,k) * det_Fe2C &
                           + wombat%bdetremi(i,j,k) * bdet_Fe2C &
                           + wombat%zooresp(i,j,k) * zoo_Fe2C &
                           + wombat%mesresp(i,j,k) * mes_Fe2C &
                           + zooexcrbac1fe &
                           + zooexcrbac2fe &
                           + zooexcraoafe &
                           + zooexcrphyfe &
                           + zooexcrdiafe &
                           + zooexcrdetfe &
                           + mesexcrbac1fe &
                           + mesexcrbac2fe &
                           + mesexcraoafe &
                           + mesexcrphyfe &
                           + mesexcrdiafe &
                           + mesexcrdetfe &
                           + mesexcrbdetfe &
                           + mesexcrzoofe &
                           + wombat%phyresp(i,j,k) * phy_Fe2C &
                           + wombat%diaresp(i,j,k) * dia_Fe2C &
                           - wombat%phy_dfeupt(i,j,k) &
                           - wombat%dia_dfeupt(i,j,k) &
                           - wombat%bac1ufer(i,j,k) &
                           - wombat%bac2ufer(i,j,k) &
                           - wombat%aoagrow(i,j,k) / wombat%aoa_C2Fe &
                           + wombat%bac1mor1(i,j,k) / wombat%bac1_C2Fe &
                           + wombat%bac1mor2(i,j,k) / wombat%bac1_C2Fe &
                           + wombat%bac2mor1(i,j,k) / wombat%bac2_C2Fe &
                           + wombat%bac2mor2(i,j,k) / wombat%bac2_C2Fe &
                           + wombat%aoamor1(i,j,k) / wombat%aoa_C2Fe &
                           + wombat%aoamor2(i,j,k) / wombat%aoa_C2Fe &
                           - wombat%feprecip(i,j,k) &
                           - wombat%fescaven(i,j,k) &
                           - wombat%fecoag2det(i,j,k) &
                           - wombat%fecoag2bdet(i,j,k) )

      ! Collect dFe sources and sinks for diagnostic output
      wombat%fesources(i,j,k) = wombat%fesources(i,j,k) + dtsb * ( 0.0 &
                                + wombat%detremi(i,j,k) * det_Fe2C &
                                + wombat%bdetremi(i,j,k) * bdet_Fe2C &
                                + wombat%zooresp(i,j,k) * zoo_Fe2C &
                                + wombat%mesresp(i,j,k) * mes_Fe2C &
                                + wombat%bac1mor1(i,j,k) / wombat%bac1_C2Fe &
                                + wombat%bac1mor2(i,j,k) / wombat%bac1_C2Fe &
                                + wombat%bac2mor1(i,j,k) / wombat%bac2_C2Fe &
                                + wombat%bac2mor2(i,j,k) / wombat%bac2_C2Fe &
                                + wombat%aoamor1(i,j,k) / wombat%aoa_C2Fe &
                                + wombat%aoamor2(i,j,k) / wombat%aoa_C2Fe &
                                + zooexcrbac1fe &
                                + zooexcrbac2fe &
                                + zooexcraoafe &
                                + zooexcrphyfe &
                                + zooexcrdiafe &
                                + zooexcrdetfe &
                                + mesexcrbac1fe &
                                + mesexcrbac2fe &
                                + mesexcraoafe &
                                + mesexcrphyfe &
                                + mesexcrdiafe &
                                + mesexcrdetfe &
                                + mesexcrbdetfe &
                                + mesexcrzoofe &
                                + wombat%phyresp(i,j,k) * phy_Fe2C &
                                + wombat%diaresp(i,j,k) * dia_Fe2C)
      wombat%fesinks(i,j,k) = wombat%fesinks(i,j,k) + dtsb * ( 0.0 & 
                              + wombat%phy_dfeupt(i,j,k) &
                              + wombat%dia_dfeupt(i,j,k) &
                              + wombat%bac1ufer(i,j,k) &
                              + wombat%bac2ufer(i,j,k) &
                              + wombat%aoagrow(i,j,k) / wombat%aoa_C2Fe &
                              + wombat%feprecip(i,j,k) &
                              + wombat%fescaven(i,j,k) &
                              + wombat%fecoag2det(i,j,k) &
                              + wombat%fecoag2bdet(i,j,k)) 


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 16] Check for conservation of mass by ecosystem component      !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      n_pools(i,j,k,2) = wombat%f_no3(i,j,k) + wombat%f_nh4(i,j,k) + wombat%f_don(i,j,k) + 2*wombat%f_n2o(i,j,k) &
                          + ( wombat%f_phy(i,j,k) + wombat%f_det(i,j,k) + wombat%f_bdet(i,j,k) &
                          +   wombat%f_zoo(i,j,k) + wombat%f_mes(i,j,k) + wombat%f_dia(i,j,k) ) * 16/122.0 &
                          + ( wombat%f_bac1(i,j,k) / wombat%bac1_C2N + wombat%f_bac2(i,j,k) / wombat%bac2_C2N &
                              + wombat%f_aoa(i,j,k) / wombat%aoa_C2N )
      c_pools(i,j,k,2) = wombat%f_dic(i,j,k) + wombat%f_phy(i,j,k) + wombat%f_det(i,j,k) + wombat%f_bdet(i,j,k) + &
                         wombat%f_zoo(i,j,k) + wombat%f_mes(i,j,k) + wombat%f_caco3(i,j,k) + wombat%f_dia(i,j,k) + &
                         wombat%f_doc(i,j,k) + wombat%f_bac1(i,j,k) + wombat%f_bac2(i,j,k) + wombat%f_aoa(i,j,k)
      si_pools(i,j,k,2) = wombat%f_sil(i,j,k) + wombat%f_diasi(i,j,k) + wombat%f_bdetsi(i,j,k)
                         

      if (tn.gt.1) then
        if (do_check_n_conserve) then
          if (abs(n_pools(i,j,k,2) - n_pools(i,j,k,1)).gt.1e-16) then
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
            print *, "       N2O (molN2/kg) =", wombat%f_n2o(i,j,k)
            print *, "       PHY (molN/kg) =", wombat%f_phy(i,j,k) * 16.0 / 122.0
            print *, "       DIA (molN/kg) =", wombat%f_dia(i,j,k) * 16.0 / 122.0
            print *, "       ZOO (molN/kg) =", wombat%f_zoo(i,j,k) * 16.0 / 122.0
            print *, "       MES (molN/kg) =", wombat%f_mes(i,j,k) * 16.0 / 122.0
            print *, "       DET (molN/kg) =", wombat%f_det(i,j,k) * 16.0 / 122.0
            print *, "       BDET (molN/kg) =", wombat%f_bdet(i,j,k) * 16.0 / 122.0
            print *, "       DON (molN/kg) =", wombat%f_don(i,j,k)
            print *, "       BAC1 (molN/kg) =", wombat%f_bac1(i,j,k) / wombat%bac1_C2N
            print *, "       BAC2 (molN/kg) =", wombat%f_bac2(i,j,k) / wombat%bac2_C2N
            print *, "       AOA (molN/kg) =", wombat%f_aoa(i,j,k) / wombat%aoa_C2N
            print *, " "
            print *, "       ammox (molN/kg/s) =", wombat%ammox(i,j,k)
            print *, "       phygrow (molC/kg/s) =", wombat%phygrow(i,j,k)
            print *, "       bac1grow (molC/kg/s) =", wombat%bac1grow(i,j,k)
            print *, "       bac2grow (molC/kg/s) =", wombat%bac2grow(i,j,k)
            print *, "       diagrow (molC/kg/s) =", wombat%diagrow(i,j,k)
            print *, "       detremi (molC/kg/s) =", wombat%detremi(i,j,k)
            print *, "       bdetremi (molC/kg/s) =", wombat%bdetremi(i,j,k)
            print *, "       bac1unh4 (molN/kg/s) =", wombat%bac1unh4(i,j,k)
            print *, "       bac2unh4 (molN/kg/s) =", wombat%bac2unh4(i,j,k)
            print *, "       don1remi (molN/kg/s) =", wombat%don1remi(i,j,k)
            print *, "       don2remi (molN/kg/s) =", wombat%don2remi(i,j,k)
            print *, "       bac1deni (molN/kg/s) =", wombat%bac1deni(i,j,k)
            print *, "       bac2deni (molN/kg/s) =", wombat%bac2deni(i,j,k)
            print *, "       anammox (molN/kg/s) =", wombat%anammox(i,j,k)
            print *, "       nitrfix (molN/kg/s) =", wombat%nitrfix(i,j,k)
            print *, "       zooresp (molC/kg/s) =", wombat%zooresp(i,j,k)
            print *, "       zooexcrphy (molC/kg/s) =", wombat%zooexcrphy(i,j,k)
            print *, "       zooexcrdia (molC/kg/s) =", wombat%zooexcrdia(i,j,k)
            print *, "       zooexcrdet (molC/kg/s) =", wombat%zooexcrdet(i,j,k)
            print *, "       mesresp (molC/kg/s) =", wombat%mesresp(i,j,k)
            print *, "       mesexcrphy (molC/kg/s) =", wombat%mesexcrphy(i,j,k)
            print *, "       mesexcrdia (molC/kg/s) =", wombat%mesexcrdia(i,j,k)
            print *, "       mesexcrzoo (molC/kg/s) =", wombat%mesexcrzoo(i,j,k)
            print *, "       mesexcrdet (molC/kg/s) =", wombat%mesexcrdet(i,j,k)
            print *, "       mesexcrbdet (molC/kg/s) =", wombat%mesexcrbdet(i,j,k)
            print *, "       phyresp (molC/kg/s) =", wombat%phyresp(i,j,k)
            print *, "       diaresp (molC/kg/s) =", wombat%diaresp(i,j,k)
            print *, "--------------------------------------------"
            call mpp_error(FATAL, trim(error_header) // " Terminating run due to non-conservation of tracer")
          endif
        endif
        if (do_check_c_conserve) then
          if (abs(c_pools(i,j,k,2) - c_pools(i,j,k,1)).gt.1e-16) then
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
            print *, "       DET (molC/kg) =", wombat%f_det(i,j,k)
            print *, "       BDET (molC/kg) =", wombat%f_bdet(i,j,k)
            print *, "       BAC1 (molC/kg) =", wombat%f_bac1(i,j,k)
            print *, "       BAC2 (molC/kg) =", wombat%f_bac2(i,j,k)
            print *, "       AOA (molC/kg) =", wombat%f_aoa(i,j,k)
            print *, "       DOC (molC/kg) =", wombat%f_doc(i,j,k)
            print *, "       CaCO3 (molC/kg) =", wombat%f_caco3(i,j,k)
            print *, "       Temp =", Temp(i,j,k)
            print *, "       Salt =", Salt(i,j,k)
            print *, "       surface pCO2 =", wombat%pco2_csurf(i,j)
            print *, "       htotal =", wombat%htotal(i,j,k)
            print *, " "
            print *, "       phygrow (molC/kg/s) =", wombat%phygrow(i,j,k)
            print *, "       diagrow (molC/kg/s) =", wombat%diagrow(i,j,k)
            print *, "       bac1grow (molC/kg/s) =", wombat%bac1grow(i,j,k)
            print *, "       bac2grow (molC/kg/s) =", wombat%bac2grow(i,j,k)
            print *, "       aoagrow (molC/kg/s) =", wombat%aoagrow(i,j,k)
            print *, "       detremi (molC/kg/s) =", wombat%detremi(i,j,k)
            print *, "       bdetremi (molC/kg/s) =", wombat%bdetremi(i,j,k)
            print *, "       zooresp (molC/kg/s) =", wombat%zooresp(i,j,k)
            print *, "       zooexcrphy (molC/kg/s) =", wombat%zooexcrphy(i,j,k)
            print *, "       zooexcrdia (molC/kg/s) =", wombat%zooexcrdia(i,j,k)
            print *, "       zooexcrdet (molC/kg/s) =", wombat%zooexcrdet(i,j,k)
            print *, "       zooresp (molC/kg/s) =", wombat%zooresp(i,j,k)
            print *, "       mesexcrphy (molC/kg/s) =", wombat%mesexcrphy(i,j,k)
            print *, "       mesexcrdia (molC/kg/s) =", wombat%mesexcrdia(i,j,k)
            print *, "       mesexcrdet (molC/kg/s) =", wombat%mesexcrdet(i,j,k)
            print *, "       mesexcrbdet (molC/kg/s) =", wombat%mesexcrbdet(i,j,k)
            print *, "       mesexcrzoo (molC/kg/s) =", wombat%mesexcrzoo(i,j,k)
            print *, "       phyresp (molC/kg/s) =", wombat%phyresp(i,j,k)
            print *, "       diaresp (molC/kg/s) =", wombat%diaresp(i,j,k)
            print *, "       zooegesphy * pic2poc(i,j,k) (molC/kg/s) =", wombat%zooegesphy(i,j,k) * wombat%pic2poc(i,j,k)
            print *, "       mesegesphy * pic2poc(i,j,k) (molC/kg/s) =", wombat%mesegesphy(i,j,k) * wombat%pic2poc(i,j,k)
            print *, "       phymort * pic2poc(i,j,k) (molC/kg/s) =", wombat%phymort(i,j,k) * wombat%pic2poc(i,j,k)
            print *, "       diamort * pic2poc(i,j,k) (molC/kg/s) =", wombat%diamort(i,j,k) * wombat%pic2poc(i,j,k)
            print *, "       zoomort * pic2poc(i,j,k) (molC/kg/s) =", wombat%zoomort(i,j,k) * wombat%pic2poc(i,j,k)
            print *, "       caldiss (molC/kg/s) =", wombat%caldiss(i,j,k)
            print *, "--------------------------------------------"
            call mpp_error(FATAL, trim(error_header) // " Terminating run due to non-conservation of tracer")
          endif
        endif
        if (do_check_si_conserve) then
          if (abs(si_pools(i,j,k,2) - si_pools(i,j,k,1)).gt.1e-16) then
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
            print *, "       BDETSI (molSi/kg) =", wombat%f_bdetsi(i,j,k)
            print *, "       Temp =", Temp(i,j,k)
            print *, "       Salt =", Salt(i,j,k)
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
      no3_bgc_change = grid_tmask(i,j,k) * (wombat%f_no3(i,j,k) - wombat%no3_prev(i,j,k)) ! [mol/kg]
      caco3_bgc_change = grid_tmask(i,j,k) * (wombat%f_caco3(i,j,k) - wombat%caco3_prev(i,j,k)) ! [mol/kg]

      wombat%pprod_gross(i,j,k) = rdtts * wombat%pprod_gross(i,j,k) * grid_tmask(i,j,k) ! [mol/kg/s]
      wombat%zprod_gross(i,j,k) = rdtts * wombat%zprod_gross(i,j,k) * grid_tmask(i,j,k) ! [mol/kg/s]
      wombat%npp3d(i,j,k) = rdtts * wombat%npp3d(i,j,k) * grid_tmask(i,j,k) ! [mol/kg/s]
      wombat%fesources(i,j,k) = rdtts * wombat%fesources(i,j,k) * grid_tmask(i,j,k) ! [mol/kg/s]
      wombat%fesinks(i,j,k) = rdtts * wombat%fesinks(i,j,k) * grid_tmask(i,j,k) ! [mol/kg/s]

      if (wombat%zw(i,j,k) .le. hblt_depth(i,j)) then
        wombat%dic_intmld(i,j) = wombat%dic_intmld(i,j) + wombat%f_dic(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%o2_intmld(i,j)  = wombat%o2_intmld(i,j)  + wombat%f_o2(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%no3_intmld(i,j) = wombat%no3_intmld(i,j) + wombat%f_no3(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%fe_intmld(i,j)  = wombat%fe_intmld(i,j)  + wombat%f_fe(i,j,k)  * rho_dzt(i,j,k) ! [mol/m2]
        wombat%phy_intmld(i,j) = wombat%phy_intmld(i,j) + wombat%f_phy(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%det_intmld(i,j) = wombat%det_intmld(i,j) + wombat%f_det(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%pprod_gross_intmld(i,j) = wombat%pprod_gross_intmld(i,j) + &
            wombat%pprod_gross(i,j,k) * rho_dzt(i,j,k) ! [mol/m2/s]
        wombat%npp_intmld(i,j) = wombat%npp_intmld(i,j) + wombat%npp3d(i,j,k) * rho_dzt(i,j,k) ! [mol/m2/s]
        wombat%radbio_intmld(i,j) = wombat%radbio_intmld(i,j) + wombat%radbio(i,j,k) * dzt(i,j,k) ! [W/m]
      endif

      if (wombat%zw(i,j,k) .le. 100) then
        wombat%dic_int100(i,j) = wombat%dic_int100(i,j) + wombat%f_dic(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%o2_int100(i,j)  = wombat%o2_int100(i,j)  + wombat%f_o2(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%no3_int100(i,j) = wombat%no3_int100(i,j) + wombat%f_no3(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%fe_int100(i,j)  = wombat%fe_int100(i,j)  + wombat%f_fe(i,j,k)  * rho_dzt(i,j,k) ! [mol/m2]
        wombat%phy_int100(i,j) = wombat%phy_int100(i,j) + wombat%f_phy(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%det_int100(i,j) = wombat%det_int100(i,j) + wombat%f_det(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%pprod_gross_int100(i,j) = wombat%pprod_gross_int100(i,j) + &
            wombat%pprod_gross(i,j,k) * rho_dzt(i,j,k) ! [mol/m2/s]
        wombat%npp_int100(i,j) = wombat%npp_int100(i,j) + wombat%npp3d(i,j,k) * rho_dzt(i,j,k) ! [mol/m2/s]
        wombat%radbio_int100(i,j) = wombat%radbio_int100(i,j) + wombat%radbio(i,j,k) * dzt(i,j,k) ! [W/m]
      endif
    enddo; enddo; enddo

    ! Additional operations on tracers
    !-----------------------------------------------------------------------
    do j = jsc,jec; do i = isc,iec;
      ! mac: bottom dFe fix to 1 nM when the water is <= 200 m deep.  
      if (grid_kmt(i,j) .gt. 0) then
        k = grid_kmt(i,j)
        if (wombat%zw(i,j,k) .le. 200) wombat%f_fe(i,j,k)= umol_m3_to_mol_kg * 0.999 ! [mol/kg]
      endif
      do k = 1,nk
        ! pjb: tune minimum dissolved iron concentration to detection limit...
        !       this is essential for ensuring dFe is replenished in upper ocean and actually
        !       looks to be the secret of PISCES ability to replicate dFe limitation in the right places
        zno3 = wombat%f_no3(i,j,k) / mmol_m3_to_mol_kg
        zfermin = min( max( 3e-2 * zno3 * zno3, 5e-2), 7e-2) * umol_m3_to_mol_kg
        wombat%f_fe(i,j,k) = max(zfermin, wombat%f_fe(i,j,k)) * grid_tmask(i,j,k)
        ! pjb: set limits of some tracers and save corrections to output for budgets
        wombat%alk_correct(i,j,k) = (max(wombat%alk_min * mmol_m3_to_mol_kg, wombat%f_alk(i,j,k) ) &
                                    - wombat%f_alk(i,j,k) ) * grid_tmask(i,j,k)
        wombat%dic_correct(i,j,k) = (max(wombat%dic_min * mmol_m3_to_mol_kg, wombat%f_dic(i,j,k) ) &
                                    - wombat%f_dic(i,j,k) ) * grid_tmask(i,j,k)
        wombat%f_alk(i,j,k) = wombat%f_alk(i,j,k) + wombat%alk_correct(i,j,k)
        wombat%f_dic(i,j,k) = wombat%f_dic(i,j,k) + wombat%dic_correct(i,j,k)
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
    call g_tracer_set_values(tracer_list, 'det', 'field', wombat%f_det, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'detfe', 'field', wombat%f_detfe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'bdet', 'field', wombat%f_bdet, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'bdetfe', 'field', wombat%f_bdetfe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'bdetsi', 'field', wombat%f_bdetsi, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'doc', 'field', wombat%f_doc, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'don', 'field', wombat%f_don, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'bac1', 'field', wombat%f_bac1, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'bac2', 'field', wombat%f_bac2, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'aoa', 'field', wombat%f_aoa, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'n2o', 'field', wombat%f_n2o, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'o2', 'field', wombat%f_o2, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'caco3', 'field', wombat%f_caco3, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'fe', 'field', wombat%f_fe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'dic', 'field', wombat%f_dic, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'dicr', 'field', wombat%f_dicr, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'alk', 'field', wombat%f_alk, isd, jsd, ntau=tau)


    !-----------------------------------------------------------------------
    ! Assign spatially variable vertical movement of tracers
    !-----------------------------------------------------------------------
    call g_tracer_get_pointer(tracer_list, 'det', 'vmove', wombat%p_wdet) ! [m/s]
    call g_tracer_get_pointer(tracer_list, 'detfe', 'vmove', wombat%p_wdetfe) ! [m/s]
    call g_tracer_get_pointer(tracer_list, 'bdet', 'vmove', wombat%p_wbdet) ! [m/s]
    call g_tracer_get_pointer(tracer_list, 'bdetfe', 'vmove', wombat%p_wbdetfe) ! [m/s]
    call g_tracer_get_pointer(tracer_list, 'bdetsi', 'vmove', wombat%p_wbdetsi) ! [m/s]
    call g_tracer_get_pointer(tracer_list, 'caco3', 'vmove', wombat%p_wcaco3) ! [m/s]

    ! Variable sinking rates of organic detritus (positive for sinking when GOLDtridiag == .true.)
    !                                            (negative for sinking when IOWtridiag ==.true.)
    ! Note: sinking distances are limited in the vertdiff solver to prevent characteristics
    ! crossing within a timestep
    do j = jsc,jec; do i = isc,iec;
      if (grid_kmt(i,j).gt.0) then
        biophy = max(epsi, wombat%f_phy(i,j,1) ) / mmol_m3_to_mol_kg  ![mmol/m3]
        biodia = max(epsi, wombat%f_dia(i,j,1) ) / mmol_m3_to_mol_kg  ![mmol/m3]
        wsink1(:) = wombat%wdetbio * max(0.01, biophy - wombat%phybiot)**(0.21) 
        wsink2(:) = wombat%wbdetbio * max(0.01, biodia - wombat%diabiot)**(0.21) 
        do k=1,nk
          wsink1(k) = wsink1(k) + 10.0/86400.0 * min(1.0, & 
                      (wombat%f_caco3(i,j,k) / (wombat%f_det(i,j,k) + wombat%f_caco3(i,j,k) + epsi)))
          ! Increase sinking rate with depth to achieve power law behaviour  
          wsink1(k) = wsink1(k) + max(0.0, wombat%zw(i,j,k)/5000.0 * (wombat%wdetmax - wsink1(k)))
          wsink2(k) = wsink2(k) + max(0.0, wombat%zw(i,j,k)/5000.0 * (wombat%wdetmax - wsink2(k)))
          ! CaCO3 sinks slower than general detritus because it tends to be smaller
          wsinkcal(k) = wsink1(k) * wombat%wcaco3/wombat%wdetbio
        enddo
        wombat%p_wdet(i,j,:) = wsink1(:) 
        wombat%p_wdetfe(i,j,:) = wsink1(:)
        wombat%p_wbdet(i,j,:) = wsink2(:) 
        wombat%p_wbdetfe(i,j,:) = wsink2(:)
        wombat%p_wbdetsi(i,j,:) = wsink2(:)
        wombat%p_wcaco3(i,j,:) = wsinkcal(:)
      else
        wombat%p_wdet(i,j,:) = 0.0
        wombat%p_wdetfe(i,j,:) = 0.0
        wombat%p_wbdet(i,j,:) = 0.0
        wombat%p_wbdetfe(i,j,:) = 0.0
        wombat%p_wbdetsi(i,j,:) = 0.0
        wombat%p_wcaco3(i,j,:) = 0.0
      endif
      ! PJB: export production through 100 metres
      k = k100(i,j)
      wombat%export_prod(i,j) = (wombat%Rho_0 * wombat%p_wdet(i,j,k)) * wombat%f_det(i,j,k) + &
                                (wombat%Rho_0 * wombat%p_wbdet(i,j,k)) * wombat%f_bdet(i,j,k) ! [mol/m2/s]
      wombat%export_inorg(i,j) = (wombat%Rho_0 * wombat%p_wcaco3(i,j,k)) * wombat%f_caco3(i,j,k) ! [mol/m2/s]
    enddo; enddo

    !-----------------------------------------------------------------------
    ! Remineralisation of sediment tracers
    !-----------------------------------------------------------------------
    call g_tracer_get_pointer(tracer_list, 'det_sediment', 'field', wombat%p_det_sediment) ! [mol/m2]
    call g_tracer_get_pointer(tracer_list, 'detfe_sediment', 'field', wombat%p_detfe_sediment) ! [mol/m2]
    call g_tracer_get_pointer(tracer_list, 'detsi_sediment', 'field', wombat%p_detsi_sediment) ! [mol/m2]
    call g_tracer_get_pointer(tracer_list, 'caco3_sediment', 'field', wombat%p_caco3_sediment) ! [mol/m2]

    ! Get bottom conditions, including those that influence bottom fluxes. Bottom conditions are
    ! calculated over a layer defined by wombat%bottom_thickness (default 1 m). This is done because
    ! the bottom layers in MOM6 are usually "vanished" layers. This approach is based on what is done
    ! in COBALT v3.
    do j = jsc,jec; do i = isc,iec;
      if (grid_kmt(i,j).gt.0) then
        k_bot = 0
        dzt_bot = 0.0
        do k = grid_kmt(i,j),1,-1
            if (dzt_bot .lt. wombat%bottom_thickness) then
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
        ! relevant over one meter, and therefore can be automatically converted to mol/m3 and then
        ! subsequently converted through the mol/kg using Rho_0. With this assumption these arrays
        ! can be added together.
        ! We add these arrays together to simulate the reducing conditions of organic-rich sediments,
        ! and to calculate a lower omega for calcite, which ensures greater rates of dissolution of
        ! CaCO3 within the sediment as organic matter accumulates.
        wombat%seddic(i,j) = wombat%seddic(i,j) + wombat%p_det_sediment(i,j,1) / wombat%Rho_0
      endif
    enddo; enddo

    call FMS_ocmip2_co2calc(CO2_dope_vec, wombat%sedmask(:,:), &
        wombat%sedtemp(:,:), wombat%sedsalt(:,:), &
        min(wombat%dic_max*mmol_m3_to_mol_kg, max(wombat%seddic(:,:), wombat%dic_min*mmol_m3_to_mol_kg)), &
        max(wombat%sedno3(:,:) / 16., 1e-9), &
        wombat%sio2(:,:), & ! dts: This is currently constant, equal to wombat%sio2_surf
        min(wombat%alk_max*mmol_m3_to_mol_kg, max(wombat%sedalk(:,:), wombat%alk_min*mmol_m3_to_mol_kg)), &
        wombat%sedhtotal(:,:)*wombat%htotal_scale_lo, &
        wombat%sedhtotal(:,:)*wombat%htotal_scale_hi, &
        wombat%sedhtotal(:,:), &
        co2_calc=trim(co2_calc), zt=wombat%seddep(:,:), &
        co3_ion=wombat%sedco3(:,:), &
        omega_calc=wombat%sedomega_cal(:,:))

    do j = jsc,jec; do i = isc,iec;
      fbc = wombat%bbioh ** (wombat%sedtemp(i,j))
      wombat%det_sed_remin(i,j) = wombat%detlrem_sed * fbc * wombat%p_det_sediment(i,j,1) ! [mol/m2/s]
      wombat%detfe_sed_remin(i,j) = wombat%detlrem_sed * fbc * wombat%p_detfe_sediment(i,j,1) ! [mol/m2/s]
      wombat%detsi_sed_remin(i,j) = wombat%detlrem_sed * fbc * wombat%p_detsi_sediment(i,j,1) ! [mol/m2/s]
      if (do_caco3_dynamics) then
        wombat%caco3_sed_remin(i,j) = wombat%caco3lrem_sed * fbc * wombat%p_caco3_sediment(i,j,1) &
                                      * max((1.0 - wombat%omegamax_sed), (1.0 - wombat%sedomega_cal(i,j)))**(4.5)
      else
        wombat%caco3_sed_remin(i,j) = wombat%caco3lrem_sed * fbc * wombat%p_caco3_sediment(i,j,1) &
                                      * (1.0 - 0.2081)**(4.5)
      endif
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
      wombat%b_nh4(i,j) = -16./122. * wombat%det_sed_remin(i,j) ! [mol/m2/s]
      wombat%b_no3(i,j) = wombat%det_sed_denit(i,j) ! [molN/m2/s]
      wombat%b_o2(i,j) = -132./16. * wombat%b_nh4(i,j) * (1.0 - wombat%fdenit(i,j))! [mol/m2/s]
      wombat%b_dic(i,j) = 122./16. * wombat%b_nh4(i,j) - wombat%caco3_sed_remin(i,j) ! [mol/m2/s]
      wombat%b_dicr(i,j) = wombat%b_dic(i,j) ! [mol/m2/s]
      wombat%b_fe(i,j) = -1.0 * wombat%detfe_sed_remin(i,j) ! [mol/m2/s]
      wombat%b_sil(i,j) = -1.0 * wombat%detsi_sed_remin(i,j) ! [mol/m2/s]
      wombat%b_alk(i,j) = -2.0 * wombat%caco3_sed_remin(i,j) + wombat%b_nh4(i,j) - wombat%b_no3(i,j) ! [mol/m2/s]
    enddo; enddo


    ! Apply remineralisation rates to sediment tracers
    !-----------------------------------------------------------------------
    do j = jsc,jec; do i = isc,iec;

      if (grid_kmt(i,j) .gt. 0) then
        wombat%p_det_sediment(i,j,1) = wombat%p_det_sediment(i,j,1) - dt * wombat%det_sed_remin(i,j) ! [mol/m2]
        wombat%p_detfe_sediment(i,j,1) = wombat%p_detfe_sediment(i,j,1) - dt * wombat%detfe_sed_remin(i,j) ! [mol/m2]
        wombat%p_detsi_sediment(i,j,1) = wombat%p_detsi_sediment(i,j,1) - dt * wombat%detsi_sed_remin(i,j) ! [mol/m2]
        wombat%p_caco3_sediment(i,j,1) = wombat%p_caco3_sediment(i,j,1) - dt * wombat%caco3_sed_remin(i,j) ! [mol/m2]
      endif
    enddo; enddo

    call g_tracer_set_values(tracer_list, 'nh4', 'btf', wombat%b_nh4, isd, jsd)
    call g_tracer_set_values(tracer_list, 'no3', 'btf', wombat%b_no3, isd, jsd)
    call g_tracer_set_values(tracer_list, 'o2', 'btf', wombat%b_o2, isd, jsd)
    call g_tracer_set_values(tracer_list, 'dic', 'btf', wombat%b_dic, isd, jsd)
    call g_tracer_set_values(tracer_list, 'dicr', 'btf', wombat%b_dicr, isd, jsd)
    call g_tracer_set_values(tracer_list, 'fe', 'btf', wombat%b_fe, isd, jsd)
    call g_tracer_set_values(tracer_list, 'sil', 'btf', wombat%b_sil, isd, jsd)
    call g_tracer_set_values(tracer_list, 'alk', 'btf', wombat%b_alk, isd, jsd)


    ! Apply back burial loss of nitrogen and alkalinity to surface
    !-----------------------------------------------------------------------
    if (do_conserve_tracers) then
      call g_tracer_get_pointer(tracer_list, 'detbury', 'field', wombat%p_detbury) ! [mol/m2/s]
      call g_tracer_get_pointer(tracer_list, 'caco3bury', 'field', wombat%p_caco3bury) ! [mol/m2/s]
      call g_tracer_get_pointer(tracer_list, 'nh4', 'stf', wombat%p_nh4_stf)
      call g_tracer_get_pointer(tracer_list, 'alk', 'stf', wombat%p_alk_stf)
      ! Spread the addition around to all grid cells
      if (sum(grid_tmask(:,:,1)).gt.0.0) then
        avedetbury = sum(wombat%p_detbury(:,:,1) * grid_dat(:,:)) / sum(grid_dat(:,:) * grid_tmask(:,:,1))  ! [mol/m2/s]
        avecaco3bury = sum(wombat%p_caco3bury(:,:,1) * grid_dat(:,:)) / sum(grid_dat(:,:) * grid_tmask(:,:,1))  ! [mol/m2/s]
        wombat%p_nh4_stf(:,:) = wombat%p_nh4_stf(:,:) + avedetbury*16.0/122.0   ! [mol/m2/s]
        wombat%p_alk_stf(:,:) = wombat%p_alk_stf(:,:) - avedetbury*16.0/122.0 + avecaco3bury*2.0  ! [mol/m2/s]
      endif
    endif

    !=======================================================================
    ! Send diagnostics
    !=======================================================================

    if (wombat%id_pco2 .gt. 0) &
      used = g_send_data(wombat%id_pco2, wombat%pco2_csurf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_htotal .gt. 0) &
      used = g_send_data(wombat%id_htotal, wombat%htotal, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_omega_ara .gt. 0) &
      used = g_send_data(wombat%id_omega_ara, wombat%omega_ara, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_omega_cal .gt. 0) &
      used = g_send_data(wombat%id_omega_cal, wombat%omega_cal, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_co3 .gt. 0) &
      used = g_send_data(wombat%id_co3, wombat%co3, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_co2_star .gt. 0) &
      used = g_send_data(wombat%id_co2_star, wombat%co2_star, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_no3_vstf .gt. 0) &
      used = g_send_data(wombat%id_no3_vstf, wombat%no3_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    
    if (wombat%id_nh4_vstf .gt. 0) &
      used = g_send_data(wombat%id_nh4_vstf, wombat%nh4_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    
    if (wombat%id_dic_vstf .gt. 0) &
      used = g_send_data(wombat%id_dic_vstf, wombat%dic_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_dicp_vstf .gt. 0) &
      used = g_send_data(wombat%id_dicp_vstf, wombat%dic_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_alk_vstf .gt. 0) &
      used = g_send_data(wombat%id_alk_vstf, wombat%alk_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_dic_correct .gt. 0) &
      used = g_send_data(wombat%id_dic_correct, wombat%dic_correct, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_alk_correct .gt. 0) &
      used = g_send_data(wombat%id_alk_correct, wombat%alk_correct, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_radbio .gt. 0) &
      used = g_send_data(wombat%id_radbio, wombat%radbio, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_radmid .gt. 0) &
      used = g_send_data(wombat%id_radmid, wombat%radmid, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_radmld .gt. 0) &
      used = g_send_data(wombat%id_radmld, wombat%radmld, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_radbio1 .gt. 0) &
      used = g_send_data(wombat%id_radbio1, wombat%radbio(:,:,1), model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_pprod_gross .gt. 0) &
      used = g_send_data(wombat%id_pprod_gross, wombat%pprod_gross, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_mumax .gt. 0) &
      used = g_send_data(wombat%id_phy_mumax, wombat%phy_mumax, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_mu .gt. 0) &
      used = g_send_data(wombat%id_phy_mu, wombat%phy_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_pchl_mu .gt. 0) &
      used = g_send_data(wombat%id_pchl_mu, wombat%pchl_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_pchl_lpar .gt. 0) &
      used = g_send_data(wombat%id_pchl_lpar, wombat%pchl_lpar, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_kni .gt. 0) &
      used = g_send_data(wombat%id_phy_kni, wombat%phy_kni, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_kfe .gt. 0) &
      used = g_send_data(wombat%id_phy_kfe, wombat%phy_kfe, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_lpar .gt. 0) &
      used = g_send_data(wombat%id_phy_lpar, wombat%phy_lpar, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_lnit .gt. 0) &
      used = g_send_data(wombat%id_phy_lnit, wombat%phy_lnit, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_lnh4 .gt. 0) &
      used = g_send_data(wombat%id_phy_lnh4, wombat%phy_lnh4, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_lno3 .gt. 0) &
      used = g_send_data(wombat%id_phy_lno3, wombat%phy_lno3, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_lfer .gt. 0) &
      used = g_send_data(wombat%id_phy_lfer, wombat%phy_lfer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_dfeupt .gt. 0) &
      used = g_send_data(wombat%id_phy_dfeupt, wombat%phy_dfeupt, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_mumax .gt. 0) &
      used = g_send_data(wombat%id_dia_mumax, wombat%dia_mumax, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_mu .gt. 0) &
      used = g_send_data(wombat%id_dia_mu, wombat%dia_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dchl_mu .gt. 0) &
      used = g_send_data(wombat%id_dchl_mu, wombat%dchl_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dchl_lpar .gt. 0) &
      used = g_send_data(wombat%id_dchl_lpar, wombat%dchl_lpar, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_kni .gt. 0) &
      used = g_send_data(wombat%id_dia_kni, wombat%dia_kni, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_kfe .gt. 0) &
      used = g_send_data(wombat%id_dia_kfe, wombat%dia_kfe, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_lpar .gt. 0) &
      used = g_send_data(wombat%id_dia_lpar, wombat%dia_lpar, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_lnit .gt. 0) &
      used = g_send_data(wombat%id_dia_lnit, wombat%dia_lnit, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_lnh4 .gt. 0) &
      used = g_send_data(wombat%id_dia_lnh4, wombat%dia_lnh4, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_lno3 .gt. 0) &
      used = g_send_data(wombat%id_dia_lno3, wombat%dia_lno3, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_lfer .gt. 0) &
      used = g_send_data(wombat%id_dia_lfer, wombat%dia_lfer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_dfeupt .gt. 0) &
      used = g_send_data(wombat%id_dia_dfeupt, wombat%dia_dfeupt, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_tri_lfer .gt. 0) &
      used = g_send_data(wombat%id_tri_lfer, wombat%tri_lfer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_tri_lpar .gt. 0) &
      used = g_send_data(wombat%id_tri_lpar, wombat%tri_lpar, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_trimumax .gt. 0) &
      used = g_send_data(wombat%id_trimumax, wombat%trimumax, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_feIII .gt. 0) &
      used = g_send_data(wombat%id_feIII, wombat%feIII, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_sileqc .gt. 0) &
      used = g_send_data(wombat%id_sileqc, wombat%sileqc, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_disssi .gt. 0) &
      used = g_send_data(wombat%id_disssi, wombat%disssi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bsidiss .gt. 0) &
      used = g_send_data(wombat%id_bsidiss, wombat%bsidiss, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_felig .gt. 0) &
      used = g_send_data(wombat%id_felig, wombat%felig, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fecol .gt. 0) &
      used = g_send_data(wombat%id_fecol, wombat%fecol, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_feprecip .gt. 0) &
      used = g_send_data(wombat%id_feprecip, wombat%feprecip, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fescaven .gt. 0) &
      used = g_send_data(wombat%id_fescaven, wombat%fescaven, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fescadet .gt. 0) &
      used = g_send_data(wombat%id_fescadet, wombat%fescadet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fescabdet .gt. 0) &
      used = g_send_data(wombat%id_fescabdet, wombat%fescabdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fecoag2det .gt. 0) &
      used = g_send_data(wombat%id_fecoag2det, wombat%fecoag2det, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fecoag2bdet .gt. 0) &
      used = g_send_data(wombat%id_fecoag2bdet, wombat%fecoag2bdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fesources .gt. 0) &
      used = g_send_data(wombat%id_fesources, wombat%fesources, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fesinks .gt. 0) &
      used = g_send_data(wombat%id_fesinks, wombat%fesinks, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_feupreg .gt. 0) &
      used = g_send_data(wombat%id_phy_feupreg, wombat%phy_feupreg, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_fedoreg .gt. 0) &
      used = g_send_data(wombat%id_phy_fedoreg, wombat%phy_fedoreg, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phygrow .gt. 0) &
      used = g_send_data(wombat%id_phygrow, wombat%phygrow, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phydoc .gt. 0) &
      used = g_send_data(wombat%id_phydoc, wombat%phydoc, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phyresp .gt. 0) &
      used = g_send_data(wombat%id_phyresp, wombat%phyresp, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phymort .gt. 0) &
      used = g_send_data(wombat%id_phymort, wombat%phymort, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_feupreg .gt. 0) &
      used = g_send_data(wombat%id_dia_feupreg, wombat%dia_feupreg, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dia_fedoreg .gt. 0) &
      used = g_send_data(wombat%id_dia_fedoreg, wombat%dia_fedoreg, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_diagrow .gt. 0) &
      used = g_send_data(wombat%id_diagrow, wombat%diagrow, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_diadoc .gt. 0) &
      used = g_send_data(wombat%id_diadoc, wombat%diadoc, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_diaresp .gt. 0) &
      used = g_send_data(wombat%id_diaresp, wombat%diaresp, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_diamort .gt. 0) &
      used = g_send_data(wombat%id_diamort, wombat%diamort, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooeps .gt. 0) &
      used = g_send_data(wombat%id_zooeps, wombat%zooeps, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooprefbac1 .gt. 0) &
      used = g_send_data(wombat%id_zooprefbac1, wombat%zooprefbac1, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooprefbac2 .gt. 0) &
      used = g_send_data(wombat%id_zooprefbac2, wombat%zooprefbac2, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooprefaoa .gt. 0) &
      used = g_send_data(wombat%id_zooprefaoa, wombat%zooprefaoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooprefphy .gt. 0) &
      used = g_send_data(wombat%id_zooprefphy, wombat%zooprefphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooprefdia .gt. 0) &
      used = g_send_data(wombat%id_zooprefdia, wombat%zooprefdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooprefdet .gt. 0) &
      used = g_send_data(wombat%id_zooprefdet, wombat%zooprefdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoograzbac1 .gt. 0) &
      used = g_send_data(wombat%id_zoograzbac1, wombat%zoograzbac1, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoograzbac2 .gt. 0) &
      used = g_send_data(wombat%id_zoograzbac2, wombat%zoograzbac2, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoograzaoa .gt. 0) &
      used = g_send_data(wombat%id_zoograzaoa, wombat%zoograzaoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoograzphy .gt. 0) &
      used = g_send_data(wombat%id_zoograzphy, wombat%zoograzphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoograzdia .gt. 0) &
      used = g_send_data(wombat%id_zoograzdia, wombat%zoograzdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoograzdet .gt. 0) &
      used = g_send_data(wombat%id_zoograzdet, wombat%zoograzdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooresp .gt. 0) &
      used = g_send_data(wombat%id_zooresp, wombat%zooresp, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoomort .gt. 0) &
      used = g_send_data(wombat%id_zoomort, wombat%zoomort, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcrbac1 .gt. 0) &
      used = g_send_data(wombat%id_zooexcrbac1, wombat%zooexcrbac1, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcrbac2 .gt. 0) &
      used = g_send_data(wombat%id_zooexcrbac2, wombat%zooexcrbac2, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcraoa .gt. 0) &
      used = g_send_data(wombat%id_zooexcraoa, wombat%zooexcraoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcrphy .gt. 0) &
      used = g_send_data(wombat%id_zooexcrphy, wombat%zooexcrphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcrdia .gt. 0) &
      used = g_send_data(wombat%id_zooexcrdia, wombat%zooexcrdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcrdet .gt. 0) &
      used = g_send_data(wombat%id_zooexcrdet, wombat%zooexcrdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooegesbac1 .gt. 0) &
      used = g_send_data(wombat%id_zooegesbac1, wombat%zooegesbac1, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooegesbac2 .gt. 0) &
      used = g_send_data(wombat%id_zooegesbac2, wombat%zooegesbac2, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooegesaoa .gt. 0) &
      used = g_send_data(wombat%id_zooegesaoa, wombat%zooegesaoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooegesphy .gt. 0) &
      used = g_send_data(wombat%id_zooegesphy, wombat%zooegesphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooegesdia .gt. 0) &
      used = g_send_data(wombat%id_zooegesdia, wombat%zooegesdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooegesdet .gt. 0) &
      used = g_send_data(wombat%id_zooegesdet, wombat%zooegesdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_meseps .gt. 0) &
      used = g_send_data(wombat%id_meseps, wombat%meseps, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefbac1 .gt. 0) &
      used = g_send_data(wombat%id_mesprefbac1, wombat%mesprefbac1, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefbac2 .gt. 0) &
      used = g_send_data(wombat%id_mesprefbac2, wombat%mesprefbac2, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefaoa .gt. 0) &
      used = g_send_data(wombat%id_mesprefaoa, wombat%mesprefaoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefphy .gt. 0) &
      used = g_send_data(wombat%id_mesprefphy, wombat%mesprefphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefdia .gt. 0) &
      used = g_send_data(wombat%id_mesprefdia, wombat%mesprefdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefdet .gt. 0) &
      used = g_send_data(wombat%id_mesprefdet, wombat%mesprefdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefbdet .gt. 0) &
      used = g_send_data(wombat%id_mesprefbdet, wombat%mesprefbdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesprefzoo .gt. 0) &
      used = g_send_data(wombat%id_mesprefzoo, wombat%mesprefzoo, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesgrazbac1 .gt. 0) &
      used = g_send_data(wombat%id_mesgrazbac1, wombat%mesgrazbac1, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    
    if (wombat%id_mesgrazbac2 .gt. 0) &
      used = g_send_data(wombat%id_mesgrazbac2, wombat%mesgrazbac2, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    
    if (wombat%id_mesgrazaoa .gt. 0) &
      used = g_send_data(wombat%id_mesgrazaoa, wombat%mesgrazaoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    
    if (wombat%id_mesgrazphy .gt. 0) &
      used = g_send_data(wombat%id_mesgrazphy, wombat%mesgrazphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesgrazdia .gt. 0) &
      used = g_send_data(wombat%id_mesgrazdia, wombat%mesgrazdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesgrazdet .gt. 0) &
      used = g_send_data(wombat%id_mesgrazdet, wombat%mesgrazdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesgrazbdet .gt. 0) &
      used = g_send_data(wombat%id_mesgrazbdet, wombat%mesgrazbdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesgrazzoo .gt. 0) &
      used = g_send_data(wombat%id_mesgrazzoo, wombat%mesgrazzoo, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesresp .gt. 0) &
      used = g_send_data(wombat%id_mesresp, wombat%mesresp, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesmort .gt. 0) &
      used = g_send_data(wombat%id_mesmort, wombat%mesmort, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrbac1 .gt. 0) &
      used = g_send_data(wombat%id_mesexcrbac1, wombat%mesexcrbac1, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrbac2 .gt. 0) &
      used = g_send_data(wombat%id_mesexcrbac2, wombat%mesexcrbac2, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcraoa .gt. 0) &
      used = g_send_data(wombat%id_mesexcraoa, wombat%mesexcraoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrphy .gt. 0) &
      used = g_send_data(wombat%id_mesexcrphy, wombat%mesexcrphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrdia .gt. 0) &
      used = g_send_data(wombat%id_mesexcrdia, wombat%mesexcrdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrdet .gt. 0) &
      used = g_send_data(wombat%id_mesexcrdet, wombat%mesexcrdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrbdet .gt. 0) &
      used = g_send_data(wombat%id_mesexcrbdet, wombat%mesexcrbdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesexcrzoo .gt. 0) &
      used = g_send_data(wombat%id_mesexcrzoo, wombat%mesexcrzoo, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegesbac1 .gt. 0) &
      used = g_send_data(wombat%id_mesegesbac1, wombat%mesegesbac1, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegesbac2 .gt. 0) &
      used = g_send_data(wombat%id_mesegesbac2, wombat%mesegesbac2, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegesaoa .gt. 0) &
      used = g_send_data(wombat%id_mesegesaoa, wombat%mesegesaoa, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegesphy .gt. 0) &
      used = g_send_data(wombat%id_mesegesphy, wombat%mesegesphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegesdia .gt. 0) &
      used = g_send_data(wombat%id_mesegesdia, wombat%mesegesdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegesdet .gt. 0) &
      used = g_send_data(wombat%id_mesegesdet, wombat%mesegesdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegesbdet .gt. 0) &
      used = g_send_data(wombat%id_mesegesbdet, wombat%mesegesbdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesegeszoo .gt. 0) &
      used = g_send_data(wombat%id_mesegeszoo, wombat%mesegeszoo, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_reminr .gt. 0) &
      used = g_send_data(wombat%id_reminr, wombat%reminr, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_doc1remi .gt. 0) &
      used = g_send_data(wombat%id_doc1remi, wombat%doc1remi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_doc2remi .gt. 0) &
      used = g_send_data(wombat%id_doc2remi, wombat%doc2remi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_don1remi .gt. 0) &
      used = g_send_data(wombat%id_don1remi, wombat%don1remi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_don2remi .gt. 0) &
      used = g_send_data(wombat%id_don2remi, wombat%don2remi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_detremi .gt. 0) &
      used = g_send_data(wombat%id_detremi, wombat%detremi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bdetremi .gt. 0) &
      used = g_send_data(wombat%id_bdetremi, wombat%bdetremi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_pic2poc .gt. 0) &
      used = g_send_data(wombat%id_pic2poc, wombat%pic2poc, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dissrat .gt. 0) &
      used = g_send_data(wombat%id_dissrat, wombat%dissrat, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_caldiss .gt. 0) &
      used = g_send_data(wombat%id_caldiss, wombat%caldiss, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoa_loxy .gt. 0) &
      used = g_send_data(wombat%id_aoa_loxy, wombat%aoa_loxy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoa_lnh4 .gt. 0) &
      used = g_send_data(wombat%id_aoa_lnh4, wombat%aoa_lnh4, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoa_yn2o .gt. 0) &
      used = g_send_data(wombat%id_aoa_yn2o, wombat%aoa_yn2o, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoa_mumax .gt. 0) &
      used = g_send_data(wombat%id_aoa_mumax, wombat%aoa_mumax, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoa_mu .gt. 0) &
      used = g_send_data(wombat%id_aoa_mu, wombat%aoa_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoagrow .gt. 0) &
      used = g_send_data(wombat%id_aoagrow, wombat%aoagrow, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    
    if (wombat%id_aoaresp .gt. 0) &
      used = g_send_data(wombat%id_aoaresp, wombat%aoaresp, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoamor1 .gt. 0) &
      used = g_send_data(wombat%id_aoamor1, wombat%aoamor1, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aoamor2 .gt. 0) &
      used = g_send_data(wombat%id_aoamor2, wombat%aoamor2, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    
    if (wombat%id_bac1grow .gt. 0) &
      used = g_send_data(wombat%id_bac1grow, wombat%bac1grow, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac1resp .gt. 0) &
      used = g_send_data(wombat%id_bac1resp, wombat%bac1resp, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac1unh4 .gt. 0) &
      used = g_send_data(wombat%id_bac1unh4, wombat%bac1unh4, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac1ufer .gt. 0) &
      used = g_send_data(wombat%id_bac1ufer, wombat%bac1ufer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac1_lnit .gt. 0) &
      used = g_send_data(wombat%id_bac1_lnit, wombat%bac1_lnit, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac1_lfer .gt. 0) &
      used = g_send_data(wombat%id_bac1_lfer, wombat%bac1_lfer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac1_mu .gt. 0) &
      used = g_send_data(wombat%id_bac1_mu, wombat%bac1_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac1_fanaer .gt. 0) &
      used = g_send_data(wombat%id_bac1_fanaer, wombat%bac1_fanaer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    
    if (wombat%id_bac1mor1 .gt. 0) &
      used = g_send_data(wombat%id_bac1mor1, wombat%bac1mor1, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac1mor2 .gt. 0) &
      used = g_send_data(wombat%id_bac1mor2, wombat%bac1mor2, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    
    if (wombat%id_bac1deni .gt. 0) &
      used = g_send_data(wombat%id_bac1deni, wombat%bac1deni, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac2grow .gt. 0) &
      used = g_send_data(wombat%id_bac2grow, wombat%bac2grow, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac2resp .gt. 0) &
      used = g_send_data(wombat%id_bac2resp, wombat%bac2resp, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac2unh4 .gt. 0) &
      used = g_send_data(wombat%id_bac2unh4, wombat%bac2unh4, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac2ufer .gt. 0) &
      used = g_send_data(wombat%id_bac2ufer, wombat%bac2ufer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac2_lnit .gt. 0) &
      used = g_send_data(wombat%id_bac2_lnit, wombat%bac2_lnit, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac2_lfer .gt. 0) &
      used = g_send_data(wombat%id_bac2_lfer, wombat%bac2_lfer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac2_mu .gt. 0) &
      used = g_send_data(wombat%id_bac2_mu, wombat%bac2_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac2_fanaer .gt. 0) &
      used = g_send_data(wombat%id_bac2_fanaer, wombat%bac2_fanaer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    
    if (wombat%id_bac2mor1 .gt. 0) &
      used = g_send_data(wombat%id_bac2mor1, wombat%bac2mor1, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_bac2mor2 .gt. 0) &
      used = g_send_data(wombat%id_bac2mor2, wombat%bac2mor2, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    
    if (wombat%id_bac2deni .gt. 0) &
      used = g_send_data(wombat%id_bac2deni, wombat%bac2deni, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aox_lnh4 .gt. 0) &
      used = g_send_data(wombat%id_aox_lnh4, wombat%aox_lnh4, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_aox_mu .gt. 0) &
      used = g_send_data(wombat%id_aox_mu, wombat%aox_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    
    if (wombat%id_nitrfix .gt. 0) &
      used = g_send_data(wombat%id_nitrfix, wombat%nitrfix, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_ammox .gt. 0) &
      used = g_send_data(wombat%id_ammox, wombat%ammox, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_anammox .gt. 0) &
      used = g_send_data(wombat%id_anammox, wombat%anammox, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_pprod_gross_2d .gt. 0) then
      wombat%pprod_gross_2d = 0.0
      do k = 1,nk; do j = jsc,jec; do i = isc,iec;
        wombat%pprod_gross_2d(i,j) = wombat%pprod_gross_2d(i,j) + &
            wombat%pprod_gross(i,j,k) * rho_dzt(i,j,k) ! [mol/m2/s]
      enddo; enddo; enddo
      used = g_send_data(wombat%id_pprod_gross_2d, wombat%pprod_gross_2d, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif

    if (wombat%id_export_prod .gt. 0) &
      used = g_send_data(wombat%id_export_prod, wombat%export_prod, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_export_inorg .gt. 0) &
      used = g_send_data(wombat%id_export_inorg, wombat%export_inorg, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_npp3d .gt. 0) &
      used = g_send_data(wombat%id_npp3d, wombat%npp3d, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_npp2d .gt. 0) then
      wombat%npp2d = 0.0
      do k = 1,nk
        wombat%npp2d(isc:iec,jsc:jec) = wombat%npp2d(isc:iec,jsc:jec) + &
            wombat%npp3d(isc:iec,jsc:jec,k) * rho_dzt(isc:iec,jsc:jec,k) ! [mol/m2/s]
      enddo
      used = g_send_data(wombat%id_npp2d, wombat%npp2d, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif

    if (wombat%id_npp1 .gt. 0) &
      used = g_send_data(wombat%id_npp1, wombat%npp3d(:,:,1), model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_zprod_gross .gt. 0) &
      used = g_send_data(wombat%id_zprod_gross, wombat%zprod_gross, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_dic_intmld .gt. 0) &
      used = g_send_data(wombat%id_dic_intmld, wombat%dic_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_o2_intmld .gt. 0) &
      used = g_send_data(wombat%id_o2_intmld, wombat%o2_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_no3_intmld .gt. 0) &
      used = g_send_data(wombat%id_no3_intmld, wombat%no3_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_fe_intmld .gt. 0) &
      used = g_send_data(wombat%id_fe_intmld, wombat%fe_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_phy_intmld .gt. 0) &
      used = g_send_data(wombat%id_phy_intmld, wombat%phy_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_det_intmld .gt. 0) &
      used = g_send_data(wombat%id_det_intmld, wombat%det_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_pprod_gross_intmld .gt. 0) &
      used = g_send_data(wombat%id_pprod_gross_intmld, wombat%pprod_gross_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_npp_intmld .gt. 0) &
      used = g_send_data(wombat%id_npp_intmld, wombat%npp_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_radbio_intmld .gt. 0) &
      used = g_send_data(wombat%id_radbio_intmld, wombat%radbio_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_dic_int100 .gt. 0) &
      used = g_send_data(wombat%id_dic_int100, wombat%dic_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_o2_int100 .gt. 0) &
      used = g_send_data(wombat%id_o2_int100, wombat%o2_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_no3_int100 .gt. 0) &
      used = g_send_data(wombat%id_no3_int100, wombat%no3_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_fe_int100 .gt. 0) &
      used = g_send_data(wombat%id_fe_int100, wombat%fe_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_phy_int100 .gt. 0) &
      used = g_send_data(wombat%id_phy_int100, wombat%phy_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_det_int100 .gt. 0) &
      used = g_send_data(wombat%id_det_int100, wombat%det_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_pprod_gross_int100 .gt. 0) &
      used = g_send_data(wombat%id_pprod_gross_int100, wombat%pprod_gross_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_npp_int100 .gt. 0) &
      used = g_send_data(wombat%id_npp_int100, wombat%npp_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_radbio_int100 .gt. 0) &
      used = g_send_data(wombat%id_radbio_int100, wombat%radbio_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_det_sed_remin .gt. 0) &
      used = g_send_data(wombat%id_det_sed_remin, wombat%det_sed_remin, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_detfe_sed_remin .gt. 0) &
      used = g_send_data(wombat%id_detfe_sed_remin, wombat%detfe_sed_remin, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_detsi_sed_remin .gt. 0) &
      used = g_send_data(wombat%id_detsi_sed_remin, wombat%detsi_sed_remin, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_det_sed_denit .gt. 0) &
      used = g_send_data(wombat%id_det_sed_denit, wombat%det_sed_denit, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_fdenit .gt. 0) &
      used = g_send_data(wombat%id_fdenit, wombat%fdenit, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_caco3_sed_remin .gt. 0) &
      used = g_send_data(wombat%id_caco3_sed_remin, wombat%caco3_sed_remin, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_zeuphot .gt. 0) &
      used = g_send_data(wombat%id_zeuphot, wombat%zeuphot, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_seddep .gt. 0) &
      used = g_send_data(wombat%id_seddep, wombat%seddep, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedmask .gt. 0) &
      used = g_send_data(wombat%id_sedmask, wombat%sedmask, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedtemp .gt. 0) &
      used = g_send_data(wombat%id_sedtemp, wombat%sedtemp, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedsalt .gt. 0) &
      used = g_send_data(wombat%id_sedsalt, wombat%sedsalt, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedno3 .gt. 0) &
      used = g_send_data(wombat%id_sedno3, wombat%sedno3, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sednh4 .gt. 0) &
      used = g_send_data(wombat%id_sednh4, wombat%sednh4, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedsil .gt. 0) &
      used = g_send_data(wombat%id_sedsil, wombat%sedsil, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedo2 .gt. 0) &
      used = g_send_data(wombat%id_sedo2, wombat%sedo2, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_seddic .gt. 0) &
      used = g_send_data(wombat%id_seddic, wombat%seddic, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedalk .gt. 0) &
      used = g_send_data(wombat%id_sedalk, wombat%sedalk, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedhtotal .gt. 0) &
      used = g_send_data(wombat%id_sedhtotal, wombat%sedhtotal, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedco3 .gt. 0) &
      used = g_send_data(wombat%id_sedco3, wombat%sedco3, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_sedomega_cal .gt. 0) &
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
    type(g_tracer_type), pointer                       :: tracer_list
    real, dimension(ilb:,jlb:), intent(in)             :: SST, SSS
    real, dimension(ilb:,jlb:,:,:), intent(in)         :: rho
    integer, intent(in)                                :: ilb, jlb, tau
    real, dimension(ilb:,jlb:,:), optional, intent(in) :: dzt

    integer                                 :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, i, j
    real                                    :: sal, ST, o2_solubility, n2o_solubility
    real                                    :: tt, tk, tk100, ts, ts2, ts3, ts4, ts5
    real                                    :: mmol_m3_to_mol_kg
    real, dimension(:,:,:), pointer         :: grid_tmask
    character(len=fm_string_len), parameter :: sub_name = 'generic_WOMBATmid_set_boundary_values'

    ! Get the necessary properties
    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, &
        grid_tmask=grid_tmask)

    call g_tracer_get_pointer(tracer_list, 'o2', 'field', wombat%p_o2)
    call g_tracer_get_pointer(tracer_list, 'n2o', 'field', wombat%p_n2o)
   
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
      call g_tracer_get_values(tracer_list, 'dic', 'field', wombat%f_dic, isd, jsd, ntau=1, &
          positive=.true.)
      call g_tracer_get_values(tracer_list, 'no3', 'field', wombat%f_no3, isd, jsd, ntau=1, &
          positive=.true.)
      call g_tracer_get_values(tracer_list, 'alk', 'field', wombat%f_alk, isd, jsd, ntau=1, &
          positive=.true.)

      do j = jsc, jec; do i = isc, iec
          wombat%htotallo(i,j) = wombat%htotal_scale_lo * wombat%htotal(i,j,1)
          wombat%htotalhi(i,j) = wombat%htotal_scale_hi * wombat%htotal(i,j,1)
      enddo; enddo 

      if ((trim(co2_calc) == 'mocsy') .and. (.not. present(dzt))) then
        call mpp_error(FATAL,"mocsy method of co2_calc needs dzt to be passed to the "// &
            "FMS_ocmip2_co2calc subroutine.")
      endif

      call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,1), &
          SST(:,:), SSS(:,:), &
          min(wombat%dic_max*mmol_m3_to_mol_kg, max(wombat%f_dic(:,:,1), wombat%dic_min*mmol_m3_to_mol_kg)), &
          max(wombat%f_no3(:,:,1) / 16., 1e-9), &
          wombat%sio2(:,:), &
          min(wombat%alk_max*mmol_m3_to_mol_kg, max(wombat%f_alk(:,:,1), wombat%alk_min*mmol_m3_to_mol_kg)), &
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

    call g_tracer_get_values(tracer_list, 'n2o', 'alpha', wombat%n2o_alpha, isd, jsd)
    call g_tracer_get_values(tracer_list, 'n2o', 'csurf', wombat%n2o_csurf ,isd, jsd)

    do j=jsc,jec ; do i=isc,iec
      !-----------------------------------------------------------------------
      ! Compute the nitrous oxide saturation concentration at 1 atm total
      ! pressure in mol/kg given the temperature (t, in deg C) and 
      ! salinity (s, in permil)
      !
      ! From Weiss and Price (1980), Marine Chemistry.
      ! The formula used is from page 353, eq (13). We use coefficients from
      ! their Table 2, column 4, for "F" in units of mol/kg/atm.
      !
      ! n2o_solubility is defined between T(freezing) <= T <= 40 deg C 
      ! We impose these bounds here.
      !
      ! check value: T = 20 deg C, S = 35 permil,
      !     n2o_solubility = 0.026883995 mol kg-1 atm-1
      !-----------------------------------------------------------------------
      sal = SSS(i,j) ; ST = SST(i,j)

      sal = min(40.0, max(0.0, sal))
      tk = 273.15 + min(40.0, max(0.0, ST))
      tk100 = tk/100.0
      
      ! mol kg-1 atm-1
      n2o_solubility = exp( wombat%a_1_n2o + wombat%a_2_n2o*(100.0/tk) &
                            + wombat%a_3_n2o*log(tk100) + wombat%a_4_n2o*tk100**2 &
                            + sal * ( wombat%b_1_n2o + wombat%b_2_n2o*tk100 &
                                      + wombat%b_3_n2o*tk100**2 ) ) * grid_tmask(i,j,1)

      !-----------------------------------------------------------------------
      !  Compute the Schmidt number of N2o in seawater using the
      !  formulation proposed by Wanninkhof (2014) Limnology and Oceanography: Methods, 12, 351-362.
      !-----------------------------------------------------------------------
      wombat%n2o_sc_no(i,j) = wombat%a1_n2o + ST * (wombat%a2_n2o + ST * (wombat%a3_n2o + ST &
          * (wombat%a4_n2o + ST * wombat%a5_n2o))) * grid_tmask(i,j,1)

      wombat%n2o_alpha(i,j) = n2o_solubility * wombat%Rho_0 ! Converts from mol/kg to mol/m3
      wombat%n2o_csurf(i,j) = wombat%p_n2o(i,j,1,tau) * wombat%Rho_0 !nnz: MOM has rho(i,j,1,tau)
    enddo; enddo

    ! Set %csurf, %alpha and %sc_no for these tracers. This will mark them
    ! for sending fluxes to coupler
    call g_tracer_set_values(tracer_list, 'n2o', 'alpha', wombat%n2o_alpha, isd, jsd)
    call g_tracer_set_values(tracer_list, 'n2o', 'csurf', wombat%n2o_csurf, isd, jsd)
    call g_tracer_set_values(tracer_list, 'n2o', 'sc_no', wombat%n2o_sc_no, isd, jsd)
    
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
    allocate(wombat%sio2(isd:ied, jsd:jed)); wombat%sio2(:,:)=wombat%sio2_surf
    allocate(wombat%co2_csurf(isd:ied, jsd:jed)); wombat%co2_csurf(:,:)=0.0
    allocate(wombat%co2_alpha(isd:ied, jsd:jed)); wombat%co2_alpha(:,:)=0.0
    allocate(wombat%co2_sc_no(isd:ied, jsd:jed)); wombat%co2_sc_no(:,:)=0.0
    allocate(wombat%pco2_csurf(isd:ied, jsd:jed)); wombat%pco2_csurf(:,:)=0.0
    allocate(wombat%o2_csurf(isd:ied, jsd:jed)); wombat%o2_csurf(:,:)=0.0
    allocate(wombat%o2_alpha(isd:ied, jsd:jed)); wombat%o2_alpha(:,:)=0.0
    allocate(wombat%o2_sc_no(isd:ied, jsd:jed)); wombat%o2_sc_no(:,:)=0.0
    allocate(wombat%n2o_csurf(isd:ied, jsd:jed)); wombat%n2o_csurf(:,:)=0.0
    allocate(wombat%n2o_alpha(isd:ied, jsd:jed)); wombat%n2o_alpha(:,:)=0.0
    allocate(wombat%n2o_sc_no(isd:ied, jsd:jed)); wombat%n2o_sc_no(:,:)=0.0
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
    allocate(wombat%f_det(isd:ied, jsd:jed, 1:nk)); wombat%f_det(:,:,:)=0.0
    allocate(wombat%f_detfe(isd:ied, jsd:jed, 1:nk)); wombat%f_detfe(:,:,:)=0.0
    allocate(wombat%f_bdet(isd:ied, jsd:jed, 1:nk)); wombat%f_bdet(:,:,:)=0.0
    allocate(wombat%f_bdetfe(isd:ied, jsd:jed, 1:nk)); wombat%f_bdetfe(:,:,:)=0.0
    allocate(wombat%f_bdetsi(isd:ied, jsd:jed, 1:nk)); wombat%f_bdetsi(:,:,:)=0.0
    allocate(wombat%f_doc(isd:ied, jsd:jed, 1:nk)); wombat%f_doc(:,:,:)=0.0
    allocate(wombat%f_don(isd:ied, jsd:jed, 1:nk)); wombat%f_don(:,:,:)=0.0
    allocate(wombat%f_bac1(isd:ied, jsd:jed, 1:nk)); wombat%f_bac1(:,:,:)=0.0
    allocate(wombat%f_bac2(isd:ied, jsd:jed, 1:nk)); wombat%f_bac2(:,:,:)=0.0
    allocate(wombat%f_aoa(isd:ied, jsd:jed, 1:nk)); wombat%f_aoa(:,:,:)=0.0
    allocate(wombat%f_n2o(isd:ied, jsd:jed, 1:nk)); wombat%f_n2o(:,:,:)=0.0
    allocate(wombat%f_o2(isd:ied, jsd:jed, 1:nk)); wombat%f_o2(:,:,:)=0.0
    allocate(wombat%f_caco3(isd:ied, jsd:jed, 1:nk)); wombat%f_caco3(:,:,:)=0.0
    allocate(wombat%f_fe(isd:ied, jsd:jed, 1:nk)); wombat%f_fe(:,:,:)=0.0

    allocate(wombat%b_nh4(isd:ied, jsd:jed)); wombat%b_nh4(:,:)=0.0
    allocate(wombat%b_no3(isd:ied, jsd:jed)); wombat%b_no3(:,:)=0.0
    allocate(wombat%b_o2(isd:ied, jsd:jed)); wombat%b_o2(:,:)=0.0
    allocate(wombat%b_dic(isd:ied, jsd:jed)); wombat%b_dic(:,:)=0.0
    allocate(wombat%b_dicr(isd:ied, jsd:jed)); wombat%b_dicr(:,:)=0.0
    allocate(wombat%b_fe(isd:ied, jsd:jed)); wombat%b_fe(:,:)=0.0
    allocate(wombat%b_sil(isd:ied, jsd:jed)); wombat%b_sil(:,:)=0.0
    allocate(wombat%b_alk(isd:ied, jsd:jed)); wombat%b_alk(:,:)=0.0

    allocate(wombat%dic_correct(isd:ied, jsd:jed, 1:nk)); wombat%dic_correct(:,:,:)=0.0
    allocate(wombat%alk_correct(isd:ied, jsd:jed, 1:nk)); wombat%alk_correct(:,:,:)=0.0
    allocate(wombat%radbio(isd:ied, jsd:jed, 1:nk)); wombat%radbio(:,:,:)=0.0
    allocate(wombat%radmid(isd:ied, jsd:jed, 1:nk)); wombat%radmid(:,:,:)=0.0
    allocate(wombat%radmld(isd:ied, jsd:jed, 1:nk)); wombat%radmld(:,:,:)=0.0
    allocate(wombat%pprod_gross(isd:ied, jsd:jed, 1:nk)); wombat%pprod_gross(:,:,:)=0.0
    allocate(wombat%pprod_gross_2d(isd:ied, jsd:jed)); wombat%pprod_gross_2d(:,:)=0.0
    allocate(wombat%zprod_gross(isd:ied, jsd:jed, 1:nk)); wombat%zprod_gross(:,:,:)=0.0
    allocate(wombat%export_prod(isd:ied, jsd:jed)); wombat%export_prod(:,:)=0.0
    allocate(wombat%export_inorg(isd:ied, jsd:jed)); wombat%export_inorg(:,:)=0.0
    allocate(wombat%npp2d(isd:ied, jsd:jed)); wombat%npp2d(:,:)=0.0
    allocate(wombat%npp3d(isd:ied, jsd:jed, 1:nk)); wombat%npp3d(:,:,:)=0.0
    allocate(wombat%phy_mumax(isd:ied, jsd:jed, 1:nk)); wombat%phy_mumax(:,:,:)=0.0
    allocate(wombat%phy_mu(isd:ied, jsd:jed, 1:nk)); wombat%phy_mu(:,:,:)=0.0
    allocate(wombat%pchl_mu(isd:ied, jsd:jed, 1:nk)); wombat%pchl_mu(:,:,:)=0.0
    allocate(wombat%pchl_lpar(isd:ied, jsd:jed, 1:nk)); wombat%pchl_lpar(:,:,:)=0.0
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
    allocate(wombat%dchl_lpar(isd:ied, jsd:jed, 1:nk)); wombat%dchl_lpar(:,:,:)=0.0
    allocate(wombat%dia_kni(isd:ied, jsd:jed, 1:nk)); wombat%dia_kni(:,:,:)=0.0
    allocate(wombat%dia_kfe(isd:ied, jsd:jed, 1:nk)); wombat%dia_kfe(:,:,:)=0.0
    allocate(wombat%dia_lpar(isd:ied, jsd:jed, 1:nk)); wombat%dia_lpar(:,:,:)=0.0
    allocate(wombat%dia_lnit(isd:ied, jsd:jed, 1:nk)); wombat%dia_lnit(:,:,:)=0.0
    allocate(wombat%dia_lnh4(isd:ied, jsd:jed, 1:nk)); wombat%dia_lnh4(:,:,:)=0.0
    allocate(wombat%dia_lno3(isd:ied, jsd:jed, 1:nk)); wombat%dia_lno3(:,:,:)=0.0
    allocate(wombat%dia_lfer(isd:ied, jsd:jed, 1:nk)); wombat%dia_lfer(:,:,:)=0.0
    allocate(wombat%dia_dfeupt(isd:ied, jsd:jed, 1:nk)); wombat%dia_dfeupt(:,:,:)=0.0
    allocate(wombat%tri_lfer(isd:ied, jsd:jed, 1:nk)); wombat%tri_lfer(:,:,:)=0.0
    allocate(wombat%tri_lpar(isd:ied, jsd:jed, 1:nk)); wombat%tri_lpar(:,:,:)=0.0
    allocate(wombat%trimumax(isd:ied, jsd:jed, 1:nk)); wombat%trimumax(:,:,:)=0.0
    allocate(wombat%feIII(isd:ied, jsd:jed, 1:nk)); wombat%feIII(:,:,:)=0.0
    allocate(wombat%sileqc(isd:ied, jsd:jed, 1:nk)); wombat%sileqc(:,:,:)=0.0
    allocate(wombat%disssi(isd:ied, jsd:jed, 1:nk)); wombat%disssi(:,:,:)=0.0
    allocate(wombat%bsidiss(isd:ied, jsd:jed, 1:nk)); wombat%bsidiss(:,:,:)=0.0
    allocate(wombat%felig(isd:ied, jsd:jed, 1:nk)); wombat%felig(:,:,:)=0.0
    allocate(wombat%fecol(isd:ied, jsd:jed, 1:nk)); wombat%fecol(:,:,:)=0.0
    allocate(wombat%feprecip(isd:ied, jsd:jed, 1:nk)); wombat%feprecip(:,:,:)=0.0
    allocate(wombat%fescaven(isd:ied, jsd:jed, 1:nk)); wombat%fescaven(:,:,:)=0.0
    allocate(wombat%fescadet(isd:ied, jsd:jed, 1:nk)); wombat%fescadet(:,:,:)=0.0
    allocate(wombat%fescabdet(isd:ied, jsd:jed, 1:nk)); wombat%fescabdet(:,:,:)=0.0
    allocate(wombat%fecoag2det(isd:ied, jsd:jed, 1:nk)); wombat%fecoag2det(:,:,:)=0.0
    allocate(wombat%fecoag2bdet(isd:ied, jsd:jed, 1:nk)); wombat%fecoag2bdet(:,:,:)=0.0
    allocate(wombat%fesources(isd:ied, jsd:jed, 1:nk)); wombat%fesources(:,:,:)=0.0
    allocate(wombat%fesinks(isd:ied, jsd:jed, 1:nk)); wombat%fesinks(:,:,:)=0.0
    allocate(wombat%phy_feupreg(isd:ied, jsd:jed, 1:nk)); wombat%phy_feupreg(:,:,:)=0.0
    allocate(wombat%phy_fedoreg(isd:ied, jsd:jed, 1:nk)); wombat%phy_fedoreg(:,:,:)=0.0
    allocate(wombat%phygrow(isd:ied, jsd:jed, 1:nk)); wombat%phygrow(:,:,:)=0.0
    allocate(wombat%phydoc(isd:ied, jsd:jed, 1:nk)); wombat%phydoc(:,:,:)=0.0
    allocate(wombat%phyresp(isd:ied, jsd:jed, 1:nk)); wombat%phyresp(:,:,:)=0.0
    allocate(wombat%phymort(isd:ied, jsd:jed, 1:nk)); wombat%phymort(:,:,:)=0.0
    allocate(wombat%dia_feupreg(isd:ied, jsd:jed, 1:nk)); wombat%dia_feupreg(:,:,:)=0.0
    allocate(wombat%dia_fedoreg(isd:ied, jsd:jed, 1:nk)); wombat%dia_fedoreg(:,:,:)=0.0
    allocate(wombat%diagrow(isd:ied, jsd:jed, 1:nk)); wombat%diagrow(:,:,:)=0.0
    allocate(wombat%diadoc(isd:ied, jsd:jed, 1:nk)); wombat%diadoc(:,:,:)=0.0
    allocate(wombat%diaresp(isd:ied, jsd:jed, 1:nk)); wombat%diaresp(:,:,:)=0.0
    allocate(wombat%diamort(isd:ied, jsd:jed, 1:nk)); wombat%diamort(:,:,:)=0.0
    allocate(wombat%zooeps(isd:ied, jsd:jed, 1:nk)); wombat%zooeps(:,:,:)=0.0
    allocate(wombat%zooprefbac1(isd:ied, jsd:jed, 1:nk)); wombat%zooprefbac1(:,:,:)=0.0
    allocate(wombat%zooprefbac2(isd:ied, jsd:jed, 1:nk)); wombat%zooprefbac2(:,:,:)=0.0
    allocate(wombat%zooprefaoa(isd:ied, jsd:jed, 1:nk)); wombat%zooprefaoa(:,:,:)=0.0
    allocate(wombat%zooprefphy(isd:ied, jsd:jed, 1:nk)); wombat%zooprefphy(:,:,:)=0.0
    allocate(wombat%zooprefdia(isd:ied, jsd:jed, 1:nk)); wombat%zooprefdia(:,:,:)=0.0
    allocate(wombat%zooprefdet(isd:ied, jsd:jed, 1:nk)); wombat%zooprefdet(:,:,:)=0.0
    allocate(wombat%zoograzbac1(isd:ied, jsd:jed, 1:nk)); wombat%zoograzbac1(:,:,:)=0.0
    allocate(wombat%zoograzbac2(isd:ied, jsd:jed, 1:nk)); wombat%zoograzbac2(:,:,:)=0.0
    allocate(wombat%zoograzaoa(isd:ied, jsd:jed, 1:nk)); wombat%zoograzaoa(:,:,:)=0.0
    allocate(wombat%zoograzphy(isd:ied, jsd:jed, 1:nk)); wombat%zoograzphy(:,:,:)=0.0
    allocate(wombat%zoograzdia(isd:ied, jsd:jed, 1:nk)); wombat%zoograzdia(:,:,:)=0.0
    allocate(wombat%zoograzdet(isd:ied, jsd:jed, 1:nk)); wombat%zoograzdet(:,:,:)=0.0
    allocate(wombat%zooresp(isd:ied, jsd:jed, 1:nk)); wombat%zooresp(:,:,:)=0.0
    allocate(wombat%zoomort(isd:ied, jsd:jed, 1:nk)); wombat%zoomort(:,:,:)=0.0
    allocate(wombat%zooexcrbac1(isd:ied, jsd:jed, 1:nk)); wombat%zooexcrbac1(:,:,:)=0.0
    allocate(wombat%zooexcrbac2(isd:ied, jsd:jed, 1:nk)); wombat%zooexcrbac2(:,:,:)=0.0
    allocate(wombat%zooexcraoa(isd:ied, jsd:jed, 1:nk)); wombat%zooexcraoa(:,:,:)=0.0
    allocate(wombat%zooexcrphy(isd:ied, jsd:jed, 1:nk)); wombat%zooexcrphy(:,:,:)=0.0
    allocate(wombat%zooexcrdia(isd:ied, jsd:jed, 1:nk)); wombat%zooexcrdia(:,:,:)=0.0
    allocate(wombat%zooexcrdet(isd:ied, jsd:jed, 1:nk)); wombat%zooexcrdet(:,:,:)=0.0
    allocate(wombat%zooegesbac1(isd:ied, jsd:jed, 1:nk)); wombat%zooegesbac1(:,:,:)=0.0
    allocate(wombat%zooegesbac2(isd:ied, jsd:jed, 1:nk)); wombat%zooegesbac2(:,:,:)=0.0
    allocate(wombat%zooegesaoa(isd:ied, jsd:jed, 1:nk)); wombat%zooegesaoa(:,:,:)=0.0
    allocate(wombat%zooegesphy(isd:ied, jsd:jed, 1:nk)); wombat%zooegesphy(:,:,:)=0.0
    allocate(wombat%zooegesdia(isd:ied, jsd:jed, 1:nk)); wombat%zooegesdia(:,:,:)=0.0
    allocate(wombat%zooegesdet(isd:ied, jsd:jed, 1:nk)); wombat%zooegesdet(:,:,:)=0.0
    allocate(wombat%meseps(isd:ied, jsd:jed, 1:nk)); wombat%meseps(:,:,:)=0.0
    allocate(wombat%mesprefbac1(isd:ied, jsd:jed, 1:nk)); wombat%mesprefbac1(:,:,:)=0.0
    allocate(wombat%mesprefbac2(isd:ied, jsd:jed, 1:nk)); wombat%mesprefbac2(:,:,:)=0.0
    allocate(wombat%mesprefaoa(isd:ied, jsd:jed, 1:nk)); wombat%mesprefaoa(:,:,:)=0.0
    allocate(wombat%mesprefphy(isd:ied, jsd:jed, 1:nk)); wombat%mesprefphy(:,:,:)=0.0
    allocate(wombat%mesprefdia(isd:ied, jsd:jed, 1:nk)); wombat%mesprefdia(:,:,:)=0.0
    allocate(wombat%mesprefdet(isd:ied, jsd:jed, 1:nk)); wombat%mesprefdet(:,:,:)=0.0
    allocate(wombat%mesprefbdet(isd:ied, jsd:jed, 1:nk)); wombat%mesprefbdet(:,:,:)=0.0
    allocate(wombat%mesprefzoo(isd:ied, jsd:jed, 1:nk)); wombat%mesprefzoo(:,:,:)=0.0
    allocate(wombat%mesgrazbac1(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazbac1(:,:,:)=0.0
    allocate(wombat%mesgrazbac2(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazbac2(:,:,:)=0.0
    allocate(wombat%mesgrazaoa(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazaoa(:,:,:)=0.0
    allocate(wombat%mesgrazphy(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazphy(:,:,:)=0.0
    allocate(wombat%mesgrazdia(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazdia(:,:,:)=0.0
    allocate(wombat%mesgrazdet(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazdet(:,:,:)=0.0
    allocate(wombat%mesgrazbdet(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazbdet(:,:,:)=0.0
    allocate(wombat%mesgrazzoo(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazzoo(:,:,:)=0.0
    allocate(wombat%mesresp(isd:ied, jsd:jed, 1:nk)); wombat%mesresp(:,:,:)=0.0
    allocate(wombat%mesmort(isd:ied, jsd:jed, 1:nk)); wombat%mesmort(:,:,:)=0.0
    allocate(wombat%mesexcrbac1(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrbac1(:,:,:)=0.0
    allocate(wombat%mesexcrbac2(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrbac2(:,:,:)=0.0
    allocate(wombat%mesexcraoa(isd:ied, jsd:jed, 1:nk)); wombat%mesexcraoa(:,:,:)=0.0
    allocate(wombat%mesexcrphy(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrphy(:,:,:)=0.0
    allocate(wombat%mesexcrdia(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrdia(:,:,:)=0.0
    allocate(wombat%mesexcrdet(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrdet(:,:,:)=0.0
    allocate(wombat%mesexcrbdet(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrbdet(:,:,:)=0.0
    allocate(wombat%mesexcrzoo(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrzoo(:,:,:)=0.0
    allocate(wombat%mesegesbac1(isd:ied, jsd:jed, 1:nk)); wombat%mesegesbac1(:,:,:)=0.0
    allocate(wombat%mesegesbac2(isd:ied, jsd:jed, 1:nk)); wombat%mesegesbac2(:,:,:)=0.0
    allocate(wombat%mesegesaoa(isd:ied, jsd:jed, 1:nk)); wombat%mesegesaoa(:,:,:)=0.0
    allocate(wombat%mesegesphy(isd:ied, jsd:jed, 1:nk)); wombat%mesegesphy(:,:,:)=0.0
    allocate(wombat%mesegesdia(isd:ied, jsd:jed, 1:nk)); wombat%mesegesdia(:,:,:)=0.0
    allocate(wombat%mesegesdet(isd:ied, jsd:jed, 1:nk)); wombat%mesegesdet(:,:,:)=0.0
    allocate(wombat%mesegesbdet(isd:ied, jsd:jed, 1:nk)); wombat%mesegesbdet(:,:,:)=0.0
    allocate(wombat%mesegeszoo(isd:ied, jsd:jed, 1:nk)); wombat%mesegeszoo(:,:,:)=0.0
    allocate(wombat%reminr(isd:ied, jsd:jed, 1:nk)); wombat%reminr(:,:,:)=0.0
    allocate(wombat%doc1remi(isd:ied, jsd:jed, 1:nk)); wombat%doc1remi(:,:,:)=0.0
    allocate(wombat%don1remi(isd:ied, jsd:jed, 1:nk)); wombat%don1remi(:,:,:)=0.0
    allocate(wombat%doc2remi(isd:ied, jsd:jed, 1:nk)); wombat%doc2remi(:,:,:)=0.0
    allocate(wombat%don2remi(isd:ied, jsd:jed, 1:nk)); wombat%don2remi(:,:,:)=0.0
    allocate(wombat%detremi(isd:ied, jsd:jed, 1:nk)); wombat%detremi(:,:,:)=0.0
    allocate(wombat%bdetremi(isd:ied, jsd:jed, 1:nk)); wombat%bdetremi(:,:,:)=0.0
    allocate(wombat%pic2poc(isd:ied, jsd:jed, 1:nk)); wombat%pic2poc(:,:,:)=0.0
    allocate(wombat%dissrat(isd:ied, jsd:jed, 1:nk)); wombat%dissrat(:,:,:)=0.0
    allocate(wombat%caldiss(isd:ied, jsd:jed, 1:nk)); wombat%caldiss(:,:,:)=0.0
    allocate(wombat%aoa_loxy(isd:ied, jsd:jed, 1:nk)); wombat%aoa_loxy(:,:,:)=0.0
    allocate(wombat%aoa_lnh4(isd:ied, jsd:jed, 1:nk)); wombat%aoa_lnh4(:,:,:)=0.0
    allocate(wombat%aoa_yn2o(isd:ied, jsd:jed, 1:nk)); wombat%aoa_yn2o(:,:,:)=0.0
    allocate(wombat%aoa_mumax(isd:ied, jsd:jed, 1:nk)); wombat%aoa_mumax(:,:,:)=0.0
    allocate(wombat%aoa_mu(isd:ied, jsd:jed, 1:nk)); wombat%aoa_mu(:,:,:)=0.0
    allocate(wombat%aoagrow(isd:ied, jsd:jed, 1:nk)); wombat%aoagrow(:,:,:)=0.0
    allocate(wombat%aoaresp(isd:ied, jsd:jed, 1:nk)); wombat%aoaresp(:,:,:)=0.0
    allocate(wombat%aoamor1(isd:ied, jsd:jed, 1:nk)); wombat%aoamor1(:,:,:)=0.0
    allocate(wombat%aoamor2(isd:ied, jsd:jed, 1:nk)); wombat%aoamor2(:,:,:)=0.0
    allocate(wombat%bac1grow(isd:ied, jsd:jed, 1:nk)); wombat%bac1grow(:,:,:)=0.0
    allocate(wombat%bac1resp(isd:ied, jsd:jed, 1:nk)); wombat%bac1resp(:,:,:)=0.0
    allocate(wombat%bac1unh4(isd:ied, jsd:jed, 1:nk)); wombat%bac1unh4(:,:,:)=0.0
    allocate(wombat%bac1ufer(isd:ied, jsd:jed, 1:nk)); wombat%bac1ufer(:,:,:)=0.0
    allocate(wombat%bac1_lnit(isd:ied, jsd:jed, 1:nk)); wombat%bac1_lnit(:,:,:)=0.0
    allocate(wombat%bac1_lfer(isd:ied, jsd:jed, 1:nk)); wombat%bac1_lfer(:,:,:)=0.0
    allocate(wombat%bac1_mu(isd:ied, jsd:jed, 1:nk)); wombat%bac1_mu(:,:,:)=0.0
    allocate(wombat%bac1_kdoc(isd:ied, jsd:jed, 1:nk)); wombat%bac1_kdoc(:,:,:)=0.0
    allocate(wombat%bac1_fanaer(isd:ied, jsd:jed, 1:nk)); wombat%bac1_fanaer(:,:,:)=0.0
    allocate(wombat%bac1mor1(isd:ied, jsd:jed, 1:nk)); wombat%bac1mor1(:,:,:)=0.0
    allocate(wombat%bac1mor2(isd:ied, jsd:jed, 1:nk)); wombat%bac1mor2(:,:,:)=0.0
    allocate(wombat%bac1deni(isd:ied, jsd:jed, 1:nk)); wombat%bac1deni(:,:,:)=0.0
    allocate(wombat%bac2grow(isd:ied, jsd:jed, 1:nk)); wombat%bac2grow(:,:,:)=0.0
    allocate(wombat%bac2resp(isd:ied, jsd:jed, 1:nk)); wombat%bac2resp(:,:,:)=0.0
    allocate(wombat%bac2unh4(isd:ied, jsd:jed, 1:nk)); wombat%bac2unh4(:,:,:)=0.0
    allocate(wombat%bac2ufer(isd:ied, jsd:jed, 1:nk)); wombat%bac2ufer(:,:,:)=0.0
    allocate(wombat%bac2_lnit(isd:ied, jsd:jed, 1:nk)); wombat%bac2_lnit(:,:,:)=0.0
    allocate(wombat%bac2_lfer(isd:ied, jsd:jed, 1:nk)); wombat%bac2_lfer(:,:,:)=0.0
    allocate(wombat%bac2_mu(isd:ied, jsd:jed, 1:nk)); wombat%bac2_mu(:,:,:)=0.0
    allocate(wombat%bac2_kdoc(isd:ied, jsd:jed, 1:nk)); wombat%bac2_kdoc(:,:,:)=0.0
    allocate(wombat%bac2_fanaer(isd:ied, jsd:jed, 1:nk)); wombat%bac2_fanaer(:,:,:)=0.0
    allocate(wombat%bac2mor1(isd:ied, jsd:jed, 1:nk)); wombat%bac2mor1(:,:,:)=0.0
    allocate(wombat%bac2mor2(isd:ied, jsd:jed, 1:nk)); wombat%bac2mor2(:,:,:)=0.0
    allocate(wombat%bac2deni(isd:ied, jsd:jed, 1:nk)); wombat%bac2deni(:,:,:)=0.0
    allocate(wombat%aox_lnh4(isd:ied, jsd:jed, 1:nk)); wombat%aox_lnh4(:,:,:)=0.0
    allocate(wombat%aox_mu(isd:ied, jsd:jed, 1:nk)); wombat%aox_mu(:,:,:)=0.0
    allocate(wombat%nitrfix(isd:ied, jsd:jed, 1:nk)); wombat%nitrfix(:,:,:)=0.0
    allocate(wombat%ammox(isd:ied, jsd:jed, 1:nk)); wombat%ammox(:,:,:)=0.0
    allocate(wombat%anammox(isd:ied, jsd:jed, 1:nk)); wombat%anammox(:,:,:)=0.0
    allocate(wombat%no3_prev(isd:ied, jsd:jed, 1:nk)); wombat%no3_prev(:,:,:)=0.0
    allocate(wombat%caco3_prev(isd:ied, jsd:jed, 1:nk)); wombat%caco3_prev(:,:,:)=0.0
    allocate(wombat%det_sed_remin(isd:ied, jsd:jed)); wombat%det_sed_remin(:,:)=0.0
    allocate(wombat%det_sed_denit(isd:ied, jsd:jed)); wombat%det_sed_denit(:,:)=0.0
    allocate(wombat%det_btm(isd:ied, jsd:jed)); wombat%det_btm(:,:)=0.0
    allocate(wombat%bdet_btm(isd:ied, jsd:jed)); wombat%bdet_btm(:,:)=0.0
    allocate(wombat%fbury(isd:ied, jsd:jed)); wombat%fbury(:,:)=0.0
    allocate(wombat%fdenit(isd:ied, jsd:jed)); wombat%fdenit(:,:)=0.0
    allocate(wombat%detfe_sed_remin(isd:ied, jsd:jed)); wombat%detfe_sed_remin(:,:)=0.0
    allocate(wombat%detsi_sed_remin(isd:ied, jsd:jed)); wombat%detsi_sed_remin(:,:)=0.0
    allocate(wombat%detfe_btm(isd:ied, jsd:jed)); wombat%detfe_btm(:,:)=0.0
    allocate(wombat%bdetfe_btm(isd:ied, jsd:jed)); wombat%bdetfe_btm(:,:)=0.0
    allocate(wombat%bdetsi_btm(isd:ied, jsd:jed)); wombat%bdetsi_btm(:,:)=0.0
    allocate(wombat%caco3_sed_remin(isd:ied, jsd:jed)); wombat%caco3_sed_remin(:,:)=0.0
    allocate(wombat%caco3_btm(isd:ied, jsd:jed)); wombat%caco3_btm(:,:)=0.0
    allocate(wombat%zw(isd:ied, jsd:jed, 1:nk)); wombat%zw(:,:,:)=0.0
    allocate(wombat%zm(isd:ied, jsd:jed, 1:nk)); wombat%zm(:,:,:)=0.0

    allocate(wombat%dic_intmld(isd:ied, jsd:jed)); wombat%dic_intmld(:,:)=0.0
    allocate(wombat%o2_intmld(isd:ied, jsd:jed)); wombat%o2_intmld(:,:)=0.0
    allocate(wombat%no3_intmld(isd:ied, jsd:jed)); wombat%no3_intmld(:,:)=0.0
    allocate(wombat%fe_intmld(isd:ied, jsd:jed)); wombat%fe_intmld(:,:)=0.0
    allocate(wombat%phy_intmld(isd:ied, jsd:jed)); wombat%phy_intmld(:,:)=0.0
    allocate(wombat%det_intmld(isd:ied, jsd:jed)); wombat%det_intmld(:,:)=0.0
    allocate(wombat%pprod_gross_intmld(isd:ied, jsd:jed)); wombat%pprod_gross_intmld(:,:)=0.0
    allocate(wombat%npp_intmld(isd:ied, jsd:jed)); wombat%npp_intmld(:,:)=0.0
    allocate(wombat%radbio_intmld(isd:ied, jsd:jed)); wombat%radbio_intmld(:,:)=0.0

    allocate(wombat%dic_int100(isd:ied, jsd:jed)); wombat%dic_int100(:,:)=0.0
    allocate(wombat%o2_int100(isd:ied, jsd:jed)); wombat%o2_int100(:,:)=0.0
    allocate(wombat%no3_int100(isd:ied, jsd:jed)); wombat%no3_int100(:,:)=0.0
    allocate(wombat%fe_int100(isd:ied, jsd:jed)); wombat%fe_int100(:,:)=0.0
    allocate(wombat%phy_int100(isd:ied, jsd:jed)); wombat%phy_int100(:,:)=0.0
    allocate(wombat%det_int100(isd:ied, jsd:jed)); wombat%det_int100(:,:)=0.0
    allocate(wombat%pprod_gross_int100(isd:ied, jsd:jed)); wombat%pprod_gross_int100(:,:)=0.0
    allocate(wombat%npp_int100(isd:ied, jsd:jed)); wombat%npp_int100(:,:)=0.0
    allocate(wombat%radbio_int100(isd:ied, jsd:jed)); wombat%radbio_int100(:,:)=0.0
    allocate(wombat%zeuphot(isd:ied, jsd:jed)); wombat%zeuphot(:,:)=0.0
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
        wombat%n2o_csurf, &
        wombat%n2o_alpha, &
        wombat%n2o_sc_no, &
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
        wombat%f_det, &
        wombat%f_detfe, &
        wombat%f_bdet, &
        wombat%f_bdetfe, &
        wombat%f_bdetsi, &
        wombat%f_doc, &
        wombat%f_don, &
        wombat%f_bac1, &
        wombat%f_bac2, &
        wombat%f_aoa, &
        wombat%f_n2o, &
        wombat%f_o2, &
        wombat%f_caco3, &
        wombat%f_fe)

    deallocate( &
        wombat%b_nh4, &
        wombat%b_no3, &
        wombat%b_o2, &
        wombat%b_dic, &
        wombat%b_dicr, &
        wombat%b_fe, &
        wombat%b_sil, &
        wombat%b_alk)

    deallocate( &
        wombat%radbio, &
        wombat%radmid, &
        wombat%radmld, &
        wombat%pprod_gross, &
        wombat%pprod_gross_2d, &
        wombat%zprod_gross, &
        wombat%export_prod, &
        wombat%export_inorg, &
        wombat%npp2d, &
        wombat%npp3d, &
        wombat%phy_mumax, &
        wombat%phy_mu, &
        wombat%pchl_mu, &
        wombat%pchl_lpar, &
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
        wombat%dchl_lpar, &
        wombat%dia_kni, &
        wombat%dia_kfe, &
        wombat%dia_lpar, &
        wombat%dia_lnit, &
        wombat%dia_lnh4, &
        wombat%dia_lno3, &
        wombat%dia_lfer, &
        wombat%dia_dfeupt, &
        wombat%tri_lfer, &
        wombat%tri_lpar, &
        wombat%trimumax, &
        wombat%feIII, &
        wombat%sileqc, &
        wombat%disssi, &
        wombat%bsidiss, &
        wombat%felig, &
        wombat%fecol, &
        wombat%feprecip, &
        wombat%fescaven, &
        wombat%fescadet, &
        wombat%fescabdet, &
        wombat%fecoag2det, &
        wombat%fecoag2bdet, &
        wombat%fesources, &
        wombat%fesinks, &
        wombat%phy_feupreg, &
        wombat%phy_fedoreg, &
        wombat%phygrow, &
        wombat%phydoc, &
        wombat%phyresp, &
        wombat%phymort, &
        wombat%dia_feupreg, &
        wombat%dia_fedoreg, &
        wombat%diagrow, &
        wombat%diadoc, &
        wombat%diaresp, &
        wombat%diamort, &
        wombat%zooprefbac1, &
        wombat%zooprefbac2, &
        wombat%zooprefaoa, &
        wombat%zooprefphy, &
        wombat%zooprefdia, &
        wombat%zooprefdet, &
        wombat%zoograzbac1, &
        wombat%zoograzbac2, &
        wombat%zoograzaoa, &
        wombat%zoograzphy, &
        wombat%zoograzdia, &
        wombat%zoograzdet, &
        wombat%zooresp, &
        wombat%zoomort, &
        wombat%zooexcrbac1, &
        wombat%zooexcrbac2, &
        wombat%zooexcraoa, &
        wombat%zooexcrphy, &
        wombat%zooexcrdia, &
        wombat%zooexcrdet, &
        wombat%zooegesbac1, &
        wombat%zooegesbac2, &
        wombat%zooegesaoa, &
        wombat%zooegesphy, &
        wombat%zooegesdia, &
        wombat%zooegesdet, &
        wombat%mesprefbac1, &
        wombat%mesprefbac2, &
        wombat%mesprefaoa, &
        wombat%mesprefphy, &
        wombat%mesprefdia, &
        wombat%mesprefdet, &
        wombat%mesprefbdet, &
        wombat%mesprefzoo, &
        wombat%mesgrazbac1, &
        wombat%mesgrazbac2, &
        wombat%mesgrazaoa, &
        wombat%mesgrazphy, &
        wombat%mesgrazdia, &
        wombat%mesgrazdet, &
        wombat%mesgrazbdet, &
        wombat%mesgrazzoo, &
        wombat%mesresp, &
        wombat%mesmort, &
        wombat%mesexcrbac1, &
        wombat%mesexcrbac2, &
        wombat%mesexcraoa, &
        wombat%mesexcrphy, &
        wombat%mesexcrdia, &
        wombat%mesexcrdet, &
        wombat%mesexcrbdet, &
        wombat%mesexcrzoo, &
        wombat%mesegesbac1, &
        wombat%mesegesbac2, &
        wombat%mesegesaoa, &
        wombat%mesegesphy, &
        wombat%mesegesdia, &
        wombat%mesegesdet, &
        wombat%mesegesbdet, &
        wombat%mesegeszoo, &
        wombat%reminr, &
        wombat%doc1remi, &
        wombat%don1remi, &
        wombat%doc2remi, &
        wombat%don2remi, &
        wombat%detremi, &
        wombat%bdetremi, &
        wombat%pic2poc, &
        wombat%dissrat, &
        wombat%caldiss, &
        wombat%aoa_loxy, &
        wombat%aoa_lnh4, &
        wombat%aoa_yn2o, &
        wombat%aoa_mumax, &
        wombat%aoa_mu, &
        wombat%aoagrow, &
        wombat%aoaresp, &
        wombat%aoamor1, &
        wombat%aoamor2, &
        wombat%bac1grow, &
        wombat%bac1resp, &
        wombat%bac1unh4, &
        wombat%bac1ufer, &
        wombat%bac1_lnit, &
        wombat%bac1_lfer, &
        wombat%bac1_mu, &
        wombat%bac1_kdoc, &
        wombat%bac1_fanaer, &
        wombat%bac1mor1, &
        wombat%bac1mor2, &
        wombat%bac1deni, &
        wombat%bac2grow, &
        wombat%bac2resp, &
        wombat%bac2unh4, &
        wombat%bac2ufer, &
        wombat%bac2_lnit, &
        wombat%bac2_lfer, &
        wombat%bac2_mu, &
        wombat%bac2_kdoc, &
        wombat%bac2_fanaer, &
        wombat%bac2mor1, &
        wombat%bac2mor2, &
        wombat%bac2deni, &
        wombat%aox_lnh4, &
        wombat%aox_mu, &
        wombat%nitrfix, &
        wombat%ammox, &
        wombat%anammox, &
        wombat%no3_prev, &
        wombat%caco3_prev, &
        wombat%det_sed_remin, &
        wombat%det_sed_denit, &
        wombat%det_btm, &
        wombat%bdet_btm, &
        wombat%fbury, &
        wombat%fdenit, &
        wombat%detfe_sed_remin, &
        wombat%detsi_sed_remin, &
        wombat%detfe_btm, &
        wombat%bdetfe_btm, &
        wombat%bdetsi_btm, &
        wombat%caco3_sed_remin, &
        wombat%caco3_btm, &
        wombat%zw, &
        wombat%zm)

    deallocate( &
        wombat%dic_intmld, &
        wombat%o2_intmld, &
        wombat%no3_intmld, &
        wombat%fe_intmld, &
        wombat%phy_intmld, &
        wombat%det_intmld, &
        wombat%pprod_gross_intmld, &
        wombat%npp_intmld, &
        wombat%radbio_intmld, &
        wombat%dic_int100, &
        wombat%o2_int100, &
        wombat%no3_int100, &
        wombat%fe_int100, &
        wombat%phy_int100, &
        wombat%det_int100, &
        wombat%pprod_gross_int100, &
        wombat%npp_int100, &
        wombat%radbio_int100)

    deallocate( &
        wombat%zeuphot, &
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
