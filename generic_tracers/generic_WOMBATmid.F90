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
!  (NO3), ammonium (NH4), bio-available iron (Fe), dissolved organic carbon
!  (DOC), dissolved inorganic carbon (DIC), calcium carbonate (CaCO3), 
!  alkalinity (ALK), and oxygen (O2). Fe is carried through the all 
!  exosystem biomass pools (phytoplankton, zooplankton and detritus).
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
!  <DATA NAME="do_open_n_cycle" TYPE="logical">
!   If true, do denitrification, anammox and nitrogen fixation
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
  logical :: do_burial           = .false.  ! permanently bury organics and CaCO3 in sediments?
  logical :: do_conserve_tracers = .false. ! add back the lost NO3 and Alk due to burial to surface?
  logical :: do_open_n_cycle     = .true.  ! N cycle has denitrification, anammox and nitrogen fixation?
  logical :: do_check_n_conserve = .false. ! check that the N fluxes balance in the ecosystem
  logical :: do_check_c_conserve = .true.  ! check that the C fluxes balance in the ecosystem

  namelist /generic_wombatmid_nml/ co2_calc, do_caco3_dynamics, do_burial, do_conserve_tracers, &
                                   do_open_n_cycle, do_check_n_conserve, do_check_c_conserve

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
        trikf, &
        trichlc, &
        trin2c, &
        zooassi, &
        zooexcr, &
        zookz, &
        zoogmax, &
        zooepsmin, &
        zooepsmax, &
        zooepsrat, &
        zoooptqf, &
        zprefphy, &
        zprefdia, &
        zprefdet, &
        zoolmor, &
        zooqmor, &
        mesassi, &
        mesexcr, &
        meskz, &
        mesgmax, &
        mesepsmin, &
        mesepsmax, &
        mesepsrat, &
        mesoptqf, &
        mprefphy, &
        mprefdia, &
        mprefdet, &
        mprefbdet, &
        mprefzoo, &
        meslmor, &
        mesqmor, &
        detlrem, &
        detlrem_sed, &
        wdetbio, &
        wbdetbio, &
        wdetmax, &
        phybiot, &
        diabiot, &
        wcaco3, &
        caco3lrem, &
        caco3lrem_sed, &
        f_inorg, &
        disscal, &
        dissara, &
        dissdet, &
        ligand, &
        fcolloid, &
        knano_dfe, &
        kscav_dfe, &
        kcoag_dfe, &
        aoakn, &
        aoako, &
        aoamumax, &
        hetko, &
        hetkn, &
        hetkd, &
        hetmumax, &
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
        a1_co2, a2_co2, a3_co2, a4_co2, &
        a1_o2, a2_o2, a3_o2, a4_o2

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
        caco3_btm, &
        det_sed_remin, &
        detfe_sed_remin, &
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
        f_phy, &
        f_dia, &
        f_pchl, &
        f_dchl, &
        f_phyfe, &
        f_diafe, &
        f_zoo, &
        f_zoofe, &
        f_mes, &
        f_mesfe, &
        f_det, &
        f_detfe, &
        f_bdet, &
        f_bdetfe, &
        f_doc, &
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
        phyresp, &
        phymort, &
        dia_feupreg, &
        dia_fedoreg, &
        diagrow, &
        diaresp, &
        diamort, &
        zooeps, &
        zoograzphy, &
        zoograzdia, &
        zoograzdet, &
        zooresp, &
        zoomort, &
        zooexcrphy, &
        zooexcrdia, &
        zooexcrdet, &
        zooslopphy, &
        zooslopdia, &
        zooslopdet, &
        zooassife, &
        meseps, &
        mesgrazphy, &
        mesgrazdia, &
        mesgrazdet, &
        mesgrazbdet, &
        mesgrazzoo, &
        mesresp, &
        mesmort, &
        mesexcrphy, &
        mesexcrdia, &
        mesexcrdet, &
        mesexcrbdet, &
        mesexcrzoo, &
        messlopphy, &
        messlopdia, &
        messlopdet, &
        messlopbdet, &
        messlopzoo, &
        mesassife, &
        reminr, &
        docremi, &
        detremi, &
        bdetremi, &
        pic2poc, &
        dissrat, &
        caldiss, &
        aoa_loxy, &
        aoa_lnh4, &
        aoa_mu, &
        het_loxy, &
        het_lno3, &
        het_ldet, &
        het_mu, &
        aox_lnh4, &
        aox_mu, &
        nitrfix, &
        ammox, &
        anammox, &
        denitrif, &
        fdenitrif, &
        no3_prev, &
        caco3_prev, &
        dic_correct, &
        alk_correct, &
        zw, &
        zm

    real, dimension(:,:,:,:), pointer :: &
        p_o2

    real, dimension(:,:,:), pointer :: &
        p_det_sediment, &
        p_detfe_sediment, &
        p_caco3_sediment, &
        p_detbury, p_caco3bury

    real, dimension(:,:,:), pointer :: &
        p_wdet, &
        p_wdetfe, &
        p_wbdet, &
        p_wbdetfe, &
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
        id_phyresp = -1, &
        id_phymort = -1, &
        id_dia_feupreg = -1, &
        id_dia_fedoreg = -1, &
        id_diagrow = -1, &
        id_diaresp = -1, &
        id_diamort = -1, &
        id_zooeps = -1, &
        id_zoograzphy = -1, &
        id_zoograzdia = -1, &
        id_zoograzdet = -1, &
        id_zooresp = -1, &
        id_zoomort = -1, &
        id_zooexcrphy = -1, &
        id_zooexcrdia = -1, &
        id_zooexcrdet = -1, &
        id_zooslopphy = -1, &
        id_zooslopdia = -1, &
        id_zooslopdet = -1, &
        id_zooassife = -1, &
        id_meseps = -1, &
        id_mesgrazphy = -1, &
        id_mesgrazdia = -1, &
        id_mesgrazdet = -1, &
        id_mesgrazbdet = -1, &
        id_mesgrazzoo = -1, &
        id_mesresp = -1, &
        id_mesmort = -1, &
        id_mesexcrphy = -1, &
        id_mesexcrdia = -1, &
        id_mesexcrdet = -1, &
        id_mesexcrbdet = -1, &
        id_mesexcrzoo = -1, &
        id_messlopphy = -1, &
        id_messlopdia = -1, &
        id_messlopdet = -1, &
        id_messlopbdet = -1, &
        id_messlopzoo = -1, &
        id_mesassife = -1, &
        id_reminr = -1, &
        id_docremi = -1, &
        id_detremi = -1, &
        id_bdetremi = -1, &
        id_pic2poc = -1, &
        id_dissrat = -1, &
        id_caldiss = -1, &
        id_aoa_loxy = -1, &
        id_aoa_lnh4 = -1, &
        id_aoa_mu = -1, &
        id_het_loxy = -1, &
        id_het_lno3 = -1, &
        id_het_ldet = -1, &
        id_het_mu = -1, &
        id_aox_lnh4 = -1, &
        id_aox_mu = -1, &
        id_nitrfix = -1, &
        id_ammox = -1, &
        id_anammox = -1, &
        id_denitrif = -1, &
        id_fdenitrif = -1, &
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
        id_caco3_sed_remin = -1, &
        id_caco3_sed_depst = -1, &
        id_zeuphot = -1, &
        id_seddep = -1, &
        id_sedmask = -1, &
        id_sedtemp = -1, &
        id_sedsalt = -1, &
        id_sedno3 = -1, &
        id_sednh4 = -1, &
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

    if (do_open_n_cycle) then
      write (stdoutunit,*) trim(note_header), &
          'Doing denitrification, anammox and nitrogen fixation'
      if (do_check_n_conserve) then
        call mpp_error(FATAL, "do_check_n_conserve = .true. is going to fail because do_open_n_cycle = .false.")
      endif
    endif

    if (do_check_n_conserve) then
      write (stdoutunit,*) trim(note_header), &
          'Checking that the ecosystem model conserves nitrogen'
    endif

    if (do_check_c_conserve) then
      write (stdoutunit,*) trim(note_header), &
          'Checking that the ecosystem model conserves carbon'
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
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_detfe_sed_remin = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'detfe_sed_depst', 'Rate of deposition of detrital iron to sediment at base of water column', &
        'h', '1', 's', 'molFe/m^2/s', 'f')
    wombat%id_detfe_sed_depst = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
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
        'diaresp', 'Respiration of microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_diaresp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'diamort', 'Mortality of microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_diamort = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooeps', 'Zooplankton prey capture rate coefficient', 'h', 'L', 's', 'm^6/mmolC^2/s', 'f')
    wombat%id_zooeps = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'zooslopphy', 'Sloppy feeding of zooplankton on phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooslopphy = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooslopdia', 'Sloppy feeding of zooplankton on microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooslopdia = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooslopdet', 'Sloppy feeding of zooplankton on detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooslopdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooassife', 'Assimilation efficiency of zooplankton for Fe', 'h', 'L', 's', '', 'f')
    wombat%id_zooassife = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'meseps', 'mesozooplankton prey capture rate coefficient', 'h', 'L', 's', 'm^6/mmolC^2/s', 'f')
    wombat%id_meseps = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'messlopphy', 'Sloppy feeding of mesozooplankton on phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_messlopphy = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'messlopdia', 'Sloppy feeding of mesozooplankton on microphytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_messlopdia = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'messlopdet', 'Sloppy feeding of mesozooplankton on detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_messlopdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'messlopbdet', 'Sloppy feeding of mesozooplankton on big detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_messlopdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'messlopzoo', 'Sloppy feeding of mesozooplankton on zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_messlopzoo = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'mesassife', 'Assimilation efficiency of mesozooplankton for Fe', 'h', 'L', 's', '', 'f')
    wombat%id_mesassife = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'reminr', 'Rate of remineralisation', 'h', 'L', 's', '/s', 'f')
    wombat%id_reminr = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'docremi', 'Remineralisation of dissolved organic carbon', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_docremi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'aoa_mu', 'Realized growth rate of Ammonia Oxidizing Archaea', 'h', 'L', 's', '/s', 'f')
    wombat%id_aoa_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'het_loxy', 'Limitation of anaerobic heterotrophic bacteria by oxygen', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_het_loxy = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'het_lno3', 'Limitation of anaerobic heterotrophic bacteria by nitrate', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_het_lno3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'het_ldet', 'Limitation of anaerobic heterotrophic bacteria by detritus', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_het_ldet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'het_mu', 'Realized growth rate of anaerobic heterotrophic bacteria', 'h', 'L', 's', '/s', 'f')
    wombat%id_het_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'denitrif', 'Denitrification rate (NO3 consumption)', 'h', 'L', 's', '[mol/kg/s]', 'f')
    wombat%id_denitrif = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fdenitrif', 'Fraction of organic matter remineralised via denitrification', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_fdenitrif = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'seddep', 'Depth of the sediment', 'h', '1', 's', 'm', 'f')
    wombat%id_seddep = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedmask', 'Mask of active sediment points', 'h', '1', 's', ' ', 'f')
    wombat%id_sedmask = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedtemp', 'Temperature at the deepest grid cell', 'h', '1', 's', 'deg C', 'f')
    wombat%id_sedtemp = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedsalt', 'Salinity at the deepest grid cell', 'h', '1', 's', 'psu', 'f')
    wombat%id_sedsalt = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedno3', 'Nitrate at the deepest grid cell', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_sedno3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sednh4', 'Ammonium at the deepest grid cell', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_sednh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedo2', 'Oxygen at the deepest grid cell', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_sedo2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'seddic', 'Dissolved inorganic carbon at the deepest grid cell', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_seddic = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedalk', 'Alkalinity at the deepest grid cell', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_sedalk = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedhtotal', 'H+ ions at the deepest grid cell', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_sedhtotal = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedco3', 'CO3 ions at the deepest grid cell', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_sedco3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'sedomega_cal', 'Calcite saturation state at the deepest grid cell', 'h', '1', 's', ' ', 'f')
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
    ! formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
    ! 7373-7382).
    call g_tracer_add_param('a1_co2', wombat%a1_co2,  2073.1)
    call g_tracer_add_param('a2_co2', wombat%a2_co2, -125.62)
    call g_tracer_add_param('a3_co2', wombat%a3_co2,  3.6276)
    call g_tracer_add_param('a4_co2', wombat%a4_co2, -0.043219)

    ! Compute the Schmidt number of O2 in seawater using the
    ! formulation proposed by Keeling et al. (1998, Global Biogeochem.
    ! Cycles, 12, 141-163).
    call g_tracer_add_param('a1_o2', wombat%a1_o2, 1638.0)
    call g_tracer_add_param('a2_o2', wombat%a2_o2, -81.83)
    call g_tracer_add_param('a3_o2', wombat%a3_o2, 1.483)
    call g_tracer_add_param('a4_o2', wombat%a4_o2, -0.008004)

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

    ! Trichodesmium half saturation constant for iron uptake [umolFe/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('trikf', wombat%trikf, 0.5)

    ! Trichodesmium typical chlorophyll to carbon ratio [mg/mg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('trichlc', wombat%trichlc, 0.01)

    ! Trichodesmium typical nitrogen to carbon ratio [mol/mol]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('trin2c', wombat%trin2c, 50.0/300.0)

    ! Zooplankton assimilation efficiency [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooassi', wombat%zooassi, 0.35)

    ! Zooplankton excretion of unassimilated prey [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooexcr', wombat%zooexcr, 0.75)

    ! Zooplankton half saturation coefficient for linear mortality
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zookz', wombat%zookz, 0.25)

    ! Zooplankton maximum grazing rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zoogmax', wombat%zoogmax, 3.0/86400.0)

    ! Zooplankton minimum prey capture rate constant [m6/mmol2/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepsmin', wombat%zooepsmin, 0.10/86400.0)

    ! Zooplankton maximum prey capture rate constant [m6/mmol2/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepsmax', wombat%zooepsmax, 0.50/86400.0)

    ! Rate of transition of epsilon from nano to microzoo [per mmolC/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepsrat', wombat%zooepsrat, 0.2)

    ! Zooplankton optimal quota of iron to carbon [mol/mol]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zoooptqf', wombat%zoooptqf, 10e-6)

    ! Zooplankton preference for phytoplankton [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefphy', wombat%zprefphy, 1.0)

    ! Zooplankton preference for microphytoplankton [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefdia', wombat%zprefdia, 0.5)

    ! Zooplankton preference for detritus [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zprefdet', wombat%zprefdet, 0.25)

    ! Zooplankton respiration rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zoolmor', wombat%zoolmor, 0.01/86400.0)

    ! Zooplankton quadratic mortality rate constant [m3/mmolN/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooqmor', wombat%zooqmor, 0.1/86400.0)

    ! Mesozooplankton assimilation efficiency [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesassi', wombat%mesassi, 0.60)

    ! Mesozooplankton excretion of unassimilated prey [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesexcr', wombat%mesexcr, 0.50)

    ! Mesozooplankton half saturation coefficient for linear mortality
    !-----------------------------------------------------------------------
    call g_tracer_add_param('meskz', wombat%meskz, 0.25)

    ! Mesozooplankton maximum grazing rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesgmax', wombat%mesgmax, 1.0/86400.0)

    ! Mesozooplankton minimum prey capture rate constant [m6/mmol2/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsmin', wombat%mesepsmin, 0.01/86400.0)

    ! Mesozooplankton maximum prey capture rate constant [m6/mmol2/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsmax', wombat%mesepsmax, 0.10/86400.0)

    ! Rate of transition of epsilon from micro to meso [per mmolC/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesepsrat', wombat%mesepsrat, 0.2)

    ! Mesozooplankton optimal quota of iron to carbon [mol/mol]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesoptqf', wombat%mesoptqf, 10e-6)

    ! Mesozooplankton preference for phytoplankton [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefphy', wombat%mprefphy, 0.25)

    ! Mesozooplankton preference for microphytoplankton [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefdia', wombat%mprefdia, 1.0)

    ! Mesozooplankton preference for detritus [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefdet', wombat%mprefdet, 0.25)

    ! Mesozooplankton preference for detritus [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefbdet', wombat%mprefbdet, 0.50)

    ! Mesozooplankton preference for zooplankton [0-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mprefzoo', wombat%mprefzoo, 1.0)

    ! Mesozooplankton respiration rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('meslmor', wombat%meslmor, 0.01/86400.0)

    ! Mesozooplankton quadratic mortality rate constant [m3/mmolN/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('mesqmor', wombat%mesqmor, 0.9/86400.0)

    ! Detritus remineralisation rate constant for <180 m; value for >=180m is half of this [1/s]
    !-----------------------------------------------------------------------
    ! Detritus remineralisation rate for (>=180 m) is hard-coded to be
    ! detlrem/2
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
    
    ! Detritus remineralisation rate constant in sediments [1/s]
    !-----------------------------------------------------------------------
    ! This would normally equal detlrem, but we set the default value to be
    ! consistent with what is in ACCESS-ESM1.5 (undocumented in Ziehn et al
    ! 2020)
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
    ! Default value matches 0.001714 day-1 in Ziehn et al 2020; differs from
    ! Hayashida et al 2020
    call g_tracer_add_param('caco3lrem', wombat%caco3lrem, 0.01/86400.0)

    ! CaCO3 remineralization rate constant in sediments [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('caco3lrem_sed', wombat%caco3lrem_sed, 0.01/86400.0)

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

    ! Ammonia Oxidizing Archaea half saturation constant for NH4 uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoakn', wombat%aoakn, 0.1)

    ! Ammonia Oxidizing Archaea half saturation constant for oxygen uptake [mmolO2/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoako', wombat%aoako, 0.1)

    ! Ammonia Oxidizing Archaea maximum growth * biomass rate [/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('aoamumax', wombat%aoamumax, 0.01/86400.0)

    ! Anaerobic heterotrophic bacteria half saturation constant for oxygen uptake [mmolO2/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('hetko', wombat%hetko, 0.1)

    ! Anaerobic heterotrophic bacteria half saturation constant for nitrate uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('hetkn', wombat%hetkn, 15.0)

    ! Anaerobic heterotrophic bacteria half saturation constant for detritus uptake [mmolC/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('hetkd', wombat%hetkd, 1.0)

    ! Facultative heterotrophic bacteria maximum growth * biomass rate [/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('hetmumax', wombat%hetmumax, 0.01/86400.0)

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

    ! Dissolved organic matter
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'doc', &
        longname = 'Dissolved organic carbon', &
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

    ! Detrital irons sitting at base of column as sediment
    !-----------------------------------------------------------------------
    ! dts: included here so included in restart
    call g_tracer_add(tracer_list, package_name, &
        name = 'detfe_sediment', &
        longname = 'Detrital iron at base of column as sediment', &
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
    integer                                 :: i, j, k, n, nz
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
    real                                    :: u_npz, g_npz, m_npz, g_peffect
    real                                    :: biono3, bionh4, biooxy, biofer
    real                                    :: biophy, biodia, biozoo, biomes, biodet, biobdet, biodoc, biocaco3
    real                                    :: biophyfe, biodiafe, biozoofe, biomesfe, zooprey, mesprey
    real                                    :: fbc
    real                                    :: no3_bgc_change, caco3_bgc_change
    real                                    :: epsi = 1.0e-30
    real                                    :: pi = 3.14159265358979
    integer                                 :: ichl
    real                                    :: par_phy_mldsum, par_z_mldsum
    real                                    :: chl, zchl, zval, phy_chlc, dia_chlc
    real                                    :: phy_limnh4, phy_limno3, phy_limdin
    real                                    :: dia_limnh4, dia_limno3, dia_limdin
    real                                    :: phy_pisl, phy_pisl2 
    real                                    :: pchl_pisl, pchl_mumin, pchl_muopt
    real                                    :: dia_pisl, dia_pisl2 
    real                                    :: dchl_pisl, dchl_mumin, dchl_muopt
    real                                    :: zooslopphyfe, zooslopdiafe, zooslopdetfe, fe_deficiency
    real                                    :: zooassiphyfe, zooassidiafe, zooassidetfe
    real                                    :: zooexcrphyfe, zooexcrdiafe, zooexcrdetfe
    real                                    :: messlopphyfe, messlopdiafe, messlopdetfe, messlopbdetfe, messlopzoofe 
    real                                    :: mesassiphyfe, mesassidiafe, mesassidetfe, mesassibdetfe, mesassizoofe
    real                                    :: mesexcrphyfe, mesexcrdiafe, mesexcrdetfe, mesexcrbdetfe, mesexcrzoofe
    real, dimension(:,:), allocatable       :: ek_bgr, par_bgr_mid, par_bgr_top
    real, dimension(:), allocatable         :: wsink1, wsink2, wsinkcal
    real                                    :: max_wsink
    real, dimension(4,61)                   :: zbgr
    real                                    :: ztemk, fe_keq, fe_par, fe_sfe, fe_tfe, partic
    real                                    :: fesol1, fesol2, fesol3, fesol4, fesol5, hp, fe3sol
    real                                    :: biof, zno3, zfermin
    real                                    :: phy_Fe2C, dia_Fe2C, zoo_Fe2C, mes_Fe2C, det_Fe2C, bdet_Fe2C
    real                                    :: phy_minqfe, phy_maxqfe
    real                                    :: dia_minqfe, dia_maxqfe
    real                                    :: zoo_slmor, mes_slmor, epsmin
    real                                    :: hco3, diss_cal, diss_ara, diss_det
    real                                    :: avedetbury, avecaco3bury
    real, dimension(:,:,:,:), allocatable   :: n_pools, c_pools
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
   

    !=======================================================================
    ! Attenuation coefficients for blue, green and red light
    !=======================================================================
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
    wombat%phyresp(:,:,:) = 0.0
    wombat%phymort(:,:,:) = 0.0
    wombat%dia_feupreg(:,:,:) = 0.0
    wombat%dia_fedoreg(:,:,:) = 0.0
    wombat%diagrow(:,:,:) = 0.0
    wombat%diaresp(:,:,:) = 0.0
    wombat%diamort(:,:,:) = 0.0
    wombat%zooeps(:,:,:) = 0.0
    wombat%zoograzphy(:,:,:) = 0.0
    wombat%zoograzdia(:,:,:) = 0.0
    wombat%zoograzdet(:,:,:) = 0.0
    wombat%zooresp(:,:,:) = 0.0
    wombat%zoomort(:,:,:) = 0.0
    wombat%zooexcrphy(:,:,:) = 0.0
    wombat%zooexcrdia(:,:,:) = 0.0
    wombat%zooexcrdet(:,:,:) = 0.0
    wombat%zooslopphy(:,:,:) = 0.0
    wombat%zooslopdia(:,:,:) = 0.0
    wombat%zooslopdet(:,:,:) = 0.0
    wombat%zooassife(:,:,:) = 0.0
    wombat%meseps(:,:,:) = 0.0
    wombat%mesgrazphy(:,:,:) = 0.0
    wombat%mesgrazdia(:,:,:) = 0.0
    wombat%mesgrazdet(:,:,:) = 0.0
    wombat%mesgrazbdet(:,:,:) = 0.0
    wombat%mesgrazzoo(:,:,:) = 0.0
    wombat%mesresp(:,:,:) = 0.0
    wombat%mesmort(:,:,:) = 0.0
    wombat%mesexcrphy(:,:,:) = 0.0
    wombat%mesexcrdia(:,:,:) = 0.0
    wombat%mesexcrdet(:,:,:) = 0.0
    wombat%mesexcrbdet(:,:,:) = 0.0
    wombat%mesexcrzoo(:,:,:) = 0.0
    wombat%messlopphy(:,:,:) = 0.0
    wombat%messlopdia(:,:,:) = 0.0
    wombat%messlopdet(:,:,:) = 0.0
    wombat%messlopbdet(:,:,:) = 0.0
    wombat%messlopzoo(:,:,:) = 0.0
    wombat%mesassife(:,:,:) = 0.0
    wombat%reminr(:,:,:) = 0.0
    wombat%docremi(:,:,:) = 0.0
    wombat%detremi(:,:,:) = 0.0
    wombat%bdetremi(:,:,:) = 0.0
    wombat%pic2poc(:,:,:) = 0.0
    wombat%dissrat(:,:,:) = 0.0
    wombat%caldiss(:,:,:) = 0.0
    wombat%aoa_loxy(:,:,:) = 0.0
    wombat%aoa_lnh4(:,:,:) = 0.0
    wombat%aoa_mu(:,:,:) = 0.0
    wombat%het_loxy(:,:,:) = 0.0
    wombat%het_lno3(:,:,:) = 0.0
    wombat%het_ldet(:,:,:) = 0.0
    wombat%het_mu(:,:,:) = 0.0
    wombat%aox_lnh4(:,:,:) = 0.0
    wombat%aox_mu(:,:,:) = 0.0
    wombat%nitrfix(:,:,:) = 0.0
    wombat%ammox(:,:,:) = 0.0
    wombat%anammox(:,:,:) = 0.0
    wombat%denitrif(:,:,:) = 0.0
    wombat%fdenitrif(:,:,:) = 0.0
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
    call g_tracer_get_values(tracer_list, 'doc', 'field', wombat%f_doc, isd, jsd, ntau=tau, &
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
    !    8.  Mortality scalings and grazing                                 !
    !    9.  CaCO3 calculations                                             !
    !    10. Implicit nitrogen fixation                                     !
    !    11. Facultative heterotrophy calculations                          !
    !    12. Chemoautotroph calculations                                    !
    !    13. Sources and sinks                                              !
    !    14. Tracer tendencies                                              !
    !    15. Check for conservation by ecosystem component                  !
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

        ! Attenuation coefficients given chlorophyll concentration
        zchl = max(0.05, min(10.0, chl) )
        ichl = nint( 41 + 20.0*log10(zchl) + epsi )
        ek_bgr(k,1) = zbgr(2,ichl) * dzt(i,j,k)
        ek_bgr(k,2) = zbgr(3,ichl) * dzt(i,j,k)
        ek_bgr(k,3) = zbgr(4,ichl) * dzt(i,j,k)

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

    do tn = 1,ts_npzd  !{

      n_pools(:,:,:,1) = n_pools(:,:,:,2)
      c_pools(:,:,:,1) = c_pools(:,:,:,2)

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
      biono3   = max(epsi, wombat%f_no3(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
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
    
      !!!~~~ Phytoplankton ~~~!!!
      pchl_pisl = phy_pisl / ( wombat%phy_mumax(i,j,k) * 86400.0 * & 
                  (1. - min(wombat%phy_lnit(i,j,k), wombat%phy_lfer(i,j,k))) + epsi )
      wombat%pchl_lpar(i,j,k) = (1. - exp(-pchl_pisl * wombat%radmld(i,j,k)))
      pchl_mumin = wombat%phyminqc * wombat%phy_mu(i,j,k) * biophy * 12.0   ![mg/m3/s] 
      pchl_muopt = wombat%phyoptqc * wombat%phy_mu(i,j,k) * biophy * 12.0   ![mg/m3/s]
      wombat%pchl_mu(i,j,k) = (pchl_muopt - pchl_mumin) * wombat%pchl_lpar(i,j,k) * &
                               min(wombat%phy_lnit(i,j,k), wombat%phy_lfer(i,j,k))
      if ( (phy_pisl * wombat%radmld(i,j,k)) .gt. 0.0 ) then
        wombat%pchl_mu(i,j,k) = pchl_mumin + wombat%pchl_mu(i,j,k) / &
                                (phy_pisl * wombat%radmld(i,j,k))
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
        wombat%dchl_mu(i,j,k) = dchl_mumin + wombat%dchl_mu(i,j,k) / &
                                (dia_pisl * wombat%radmld(i,j,k))
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


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 7] Iron chemistry (Aumont et al., 2015; GMD)                   !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Estimate solubility of Fe3+ (free Fe) in solution using temperature, pH and salinity
      ztemk = max(5.0, Temp(i,j,k)) + 273.15    ! temperature in kelvin
      zval = 19.924 * Salt(i,j,k) / ( 1000. - 1.005 * Salt(i,j,k))
      fesol1 = 10**(-13.486 - 0.1856*zval**0.5 + 0.3073*zval + 5254.0/ztemk)
      fesol2 = 10**(2.517 - 0.8885*zval**0.5 + 0.2139*zval - 1320.0/ztemk)
      fesol3 = 10**(0.4511 - 0.3305*zval**0.5 - 1996.0/ztemk)
      fesol4 = 10**(-0.2965 - 0.7881*zval**0.5 - 4086.0/ztemk)
      fesol5 = 10**(4.4466 - 0.8505*zval**0.5 - 7980.0/ztemk)
      hp = 10**(-7.9)
      if (wombat%htotal(i,j,k).gt.0.0) hp = wombat%htotal(i,j,k)
      fe3sol = fesol1 * ( hp**3 + fesol2 * hp**2 + fesol3 * hp + fesol4 + fesol5 / hp ) *1e9

      ! Estimate total colloidal iron
      ! ... for now, we assume that 50% of all dFe is colloidal, and we separate this from the 
      !     equilibrium fractionation between Fe' and Fe-L below
      wombat%fecol(i,j,k) = wombat%fcolloid * biofer 

      ! Determine equilibriuim fractionation of the remain dFe (non-colloidal fraction) into Fe' and L-Fe
      fe_keq = 10**( 17.27 - 1565.7 / ztemk ) * 1e-9 ! Temperature reduces solubility
      fe_par = 4.77e-7 * wombat%radbio(i,j,k) * 0.5 / 10**(-6.3) ! Light increases solubility
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
      biodoc = 40.0 + (1.0 - min(wombat%phy_lnit(i,j,k), wombat%phy_lfer(i,j,k))) * 40.0 ! proxy of DOC (mmol/m3)
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
      !  [Step 8] Mortality scalings and grazing                              !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! reduce linear mortality (respiration losses) of zooplankton when there is low biomass
      zoo_slmor = biozoo / (biozoo + wombat%zookz)
      mes_slmor = biomes / (biomes + wombat%meskz)
     
      !!!~~~ Zooplankton ~~~!!!
      ! Grazing function ! [1/s]
      zooprey = wombat%zprefphy * biophy + wombat%zprefdia * biodia + wombat%zprefdet * biodet
      ! Epsilon (prey capture rate coefficient) is made a function of prey density (Fig 2 of Rohr et al., 2024; GRL)
      !  - scales towards lower values (microzooplankton) as prey biomass increases
      g_peffect = exp(-zooprey * wombat%zooepsrat)
      wombat%zooeps(i,j,k) = wombat%zooepsmin + (wombat%zooepsmax - wombat%zooepsmin) * g_peffect 
      g_npz = wombat%zoogmax * fbc * (wombat%zooeps(i,j,k) * zooprey*zooprey) / &
              (wombat%zoogmax * fbc + (wombat%zooeps(i,j,k) * zooprey*zooprey))

      !!!~~~ Mesozooplankton ~~~!!!
      ! Grazing function ! [1/s]
      mesprey = wombat%mprefphy * biophy + wombat%mprefdia * biodia + &
                wombat%mprefdet * biodet + wombat%mprefbdet * biobdet + wombat%mprefzoo * biozoo
      ! Epsilon (prey capture rate coefficient) is made a function of prey density (Fig 2 of Rohr et al., 2024; GRL)
      !  - scales towards lower values (mesozooplankton) as prey biomass increases
      g_peffect = exp(-mesprey * wombat%mesepsrat)
      wombat%meseps(i,j,k) = wombat%mesepsmin + (wombat%mesepsmax - wombat%mesepsmin) * g_peffect 
      m_npz = wombat%mesgmax * fbc * (wombat%meseps(i,j,k) * mesprey*mesprey) / &
              (wombat%mesgmax * fbc + (wombat%meseps(i,j,k) * mesprey*mesprey))


      ! Zooplankton Fe uptake
      fe_deficiency = wombat%zoooptqf - zoo_Fe2C
      wombat%zooassife(i,j,k) = max(0.5, min(1.5, (1.0 + 2.0 * fe_deficiency / wombat%zoooptqf ) ))
      wombat%zooassife(i,j,k) = min(0.9, max(0.1, wombat%zooassi * wombat%zooassife(i,j,k) ))
      ! Mesozooplankton Fe uptake
      fe_deficiency = wombat%mesoptqf - mes_Fe2C
      wombat%mesassife(i,j,k) = max(0.5, min(1.5, (1.0 + 2.0 * fe_deficiency / wombat%mesoptqf ) ))
      wombat%mesassife(i,j,k) = min(0.9, max(0.1, wombat%mesassi * wombat%mesassife(i,j,k) ))


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 9] CaCO3 calculations                                          !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      if (do_caco3_dynamics) then

        ! PIC:POC ratio is a function of the substrate:inhibitor ratio, which is the 
        !  HCO3- to free H+ ions ratio (mol/umol), following Lehmann & Bach (2024). 
        !  We also add a T-dependent function to scale down CaCO3 production in waters colder 
        !  than 3 degrees C based off the observation of no E hux growth beneath this (Fielding 2013; L&O)
        hco3 = wombat%f_dic(i,j,k) - wombat%co3(i,j,k) - wombat%co2_star(i,j,k)
        wombat%pic2poc(i,j,k) = min(0.3, (wombat%f_inorg + 10**(-3.0 + 4.31 * &
                                          hco3 / (wombat%htotal(i,j,k)*1e6))) * &
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
      !  [Step 10] Implicit nitrogen fixation                                 !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      
      if (do_open_n_cycle) then
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
      !  [Step 10] Anaerobic heterotroph calculations                         !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      
      if (do_open_n_cycle) then
        ! Anaerobic heterotrophic growth rate (currently only does denitrification)
        wombat%het_loxy(i,j,k) = (1.0 - biooxy**5.0 / (biooxy**5.0 + wombat%hetko))
        wombat%het_lno3(i,j,k) = biono3 / (biono3 + wombat%hetkn)
        wombat%het_ldet(i,j,k) = (biodet + biobdet) / ((biodet + biobdet) + wombat%hetkd)
        wombat%het_mu(i,j,k) = wombat%hetmumax * wombat%bbioh**(Temp(i,j,k)) &
                               * wombat%het_loxy(i,j,k) &
                               * min(wombat%het_lno3(i,j,k), wombat%het_ldet(i,j,k))
      endif
      

      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 11] Chemoautotroph calculations                                !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      
      ! Complete ammonia oxidation to nitrate with limitation from O2 and NH4 availability
      wombat%aoa_loxy(i,j,k) = biooxy / (biooxy + wombat%aoako)
      wombat%aoa_lnh4(i,j,k) = bionh4 / (bionh4 + wombat%aoakn)
      wombat%aoa_mu(i,j,k) = wombat%aoamumax * wombat%bbioh**(Temp(i,j,k)) &
                             * min(wombat%aoa_loxy(i,j,k), wombat%aoa_lnh4(i,j,k)) 
      
      if (do_open_n_cycle) then
        ! Anaerobic ammonium oxidation (anammox)
        wombat%aox_lnh4(i,j,k) = bionh4 / (bionh4 + wombat%aoxkn)
        wombat%aox_mu(i,j,k) = wombat%aoxmumax * wombat%bbioh**(Temp(i,j,k)) &
                               * wombat%het_loxy(i,j,k) * wombat%aox_lnh4(i,j,k)
      endif
      
    
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 12] Sources and sinks                                          !
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

      ! Chemoautotrophy
      if (wombat%f_nh4(i,j,k) .gt. epsi) then
        wombat%ammox(i,j,k) = wombat%aoa_mu(i,j,k) * wombat%f_nh4(i,j,k) ! [molN/kg/s]
        wombat%anammox(i,j,k) = wombat%aox_mu(i,j,k) * wombat%f_nh4(i,j,k) ! [molN/kg/s]
      else
        wombat%ammox(i,j,k) = 0.0
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
      if (wombat%f_doc(i,j,k) .gt. epsi) then
        wombat%docremi(i,j,k) = wombat%reminr(i,j,k) / mmol_m3_to_mol_kg * wombat%f_doc(i,j,k)**2.0 ! [molC/kg/s]
      else
        wombat%docremi(i,j,k) = 0.0
      endif
      
      ! Denitrification
      if (wombat%f_det(i,j,k) .gt. epsi) then
        wombat%denitrif(i,j,k) = min(wombat%het_mu(i,j,k) * wombat%f_doc(i,j,k) * 94.0 / 122.0, &
                                     0.9 * wombat%docremi(i,j,k) * 94.0 / 122.0 ) ! [molN/kg/s]
        wombat%fdenitrif(i,j,k) = wombat%denitrif(i,j,k) * 122.0/94.0 / (wombat%docremi(i,j,k) + epsi)
      else
        wombat%denitrif(i,j,k) = 0.0
        wombat%fdenitrif(i,j,k) = 0.0
      endif

      ! Grazing by microzooplankton
      if (zooprey.gt.1e-3) then
        wombat%zoograzphy(i,j,k) = g_npz * wombat%f_zoo(i,j,k) * (wombat%zprefphy*biophy)/zooprey ! [molC/kg/s]
        wombat%zoograzdia(i,j,k) = g_npz * wombat%f_zoo(i,j,k) * (wombat%zprefdia*biodia)/zooprey ! [molC/kg/s]
        wombat%zoograzdet(i,j,k) = g_npz * wombat%f_zoo(i,j,k) * (wombat%zprefdet*biodet)/zooprey ! [molC/kg/s]
      else
        wombat%zoograzphy(i,j,k) = 0.0
        wombat%zoograzdia(i,j,k) = 0.0
        wombat%zoograzdet(i,j,k) = 0.0
      endif
      wombat%zooexcrphy(i,j,k) = wombat%zoograzphy(i,j,k) * (1.0 - wombat%zooassi)*wombat%zooexcr
      wombat%zooexcrdia(i,j,k) = wombat%zoograzdia(i,j,k) * (1.0 - wombat%zooassi)*wombat%zooexcr
      wombat%zooexcrdet(i,j,k) = wombat%zoograzdet(i,j,k) * (1.0 - wombat%zooassi)*wombat%zooexcr
      wombat%zooslopphy(i,j,k) = wombat%zoograzphy(i,j,k) * (1.0 - wombat%zooassi)*(1.0-wombat%zooexcr)
      wombat%zooslopdia(i,j,k) = wombat%zoograzdia(i,j,k) * (1.0 - wombat%zooassi)*(1.0-wombat%zooexcr)
      wombat%zooslopdet(i,j,k) = wombat%zoograzdet(i,j,k) * (1.0 - wombat%zooassi)*(1.0-wombat%zooexcr)
      zooslopphyfe = wombat%zooslopphy(i,j,k) * phy_Fe2C
      zooslopdiafe = wombat%zooslopdia(i,j,k) * dia_Fe2C
      zooslopdetfe = wombat%zooslopdet(i,j,k) * det_Fe2C
      zooassiphyfe = wombat%zoograzphy(i,j,k) * wombat%zooassife(i,j,k) * phy_Fe2C 
      zooassidiafe = wombat%zoograzdia(i,j,k) * wombat%zooassife(i,j,k) * dia_Fe2C
      zooassidetfe = wombat%zoograzdet(i,j,k) * wombat%zooassife(i,j,k) * det_Fe2C
      zooexcrphyfe = wombat%zoograzphy(i,j,k)*phy_Fe2C - zooassiphyfe - zooslopphyfe
      zooexcrdiafe = wombat%zoograzdia(i,j,k)*dia_Fe2C - zooassidiafe - zooslopdiafe
      zooexcrdetfe = wombat%zoograzdet(i,j,k)*det_Fe2C - zooassidetfe - zooslopdetfe

      ! Grazing by mesozooplankton
      if (mesprey.gt.1e-3) then
        wombat%mesgrazphy(i,j,k) = m_npz * wombat%f_mes(i,j,k) * (wombat%mprefphy*biophy)/mesprey ! [molC/kg/s]
        wombat%mesgrazdia(i,j,k) = m_npz * wombat%f_mes(i,j,k) * (wombat%mprefdia*biodia)/mesprey ! [molC/kg/s]
        wombat%mesgrazdet(i,j,k) = m_npz * wombat%f_mes(i,j,k) * (wombat%mprefdet*biodet)/mesprey ! [molC/kg/s]
        wombat%mesgrazbdet(i,j,k) = m_npz * wombat%f_mes(i,j,k) * (wombat%mprefbdet*biobdet)/mesprey ! [molC/kg/s]
        wombat%mesgrazzoo(i,j,k) = m_npz * wombat%f_mes(i,j,k) * (wombat%mprefzoo*biozoo)/mesprey ! [molC/kg/s]
      else
        wombat%mesgrazphy(i,j,k) = 0.0
        wombat%mesgrazdia(i,j,k) = 0.0
        wombat%mesgrazdet(i,j,k) = 0.0
        wombat%mesgrazbdet(i,j,k) = 0.0
        wombat%mesgrazzoo(i,j,k) = 0.0
      endif
      wombat%mesexcrphy(i,j,k) = wombat%mesgrazphy(i,j,k) * (1.0 - wombat%mesassi)*wombat%mesexcr
      wombat%mesexcrdia(i,j,k) = wombat%mesgrazdia(i,j,k) * (1.0 - wombat%mesassi)*wombat%mesexcr
      wombat%mesexcrdet(i,j,k) = wombat%mesgrazdet(i,j,k) * (1.0 - wombat%mesassi)*wombat%mesexcr
      wombat%mesexcrbdet(i,j,k) = wombat%mesgrazbdet(i,j,k) * (1.0 - wombat%mesassi)*wombat%mesexcr
      wombat%mesexcrzoo(i,j,k) = wombat%mesgrazzoo(i,j,k) * (1.0 - wombat%mesassi)*wombat%mesexcr
      wombat%messlopphy(i,j,k) = wombat%mesgrazphy(i,j,k) * (1.0 - wombat%mesassi)*(1.0-wombat%mesexcr)
      wombat%messlopdia(i,j,k) = wombat%mesgrazdia(i,j,k) * (1.0 - wombat%mesassi)*(1.0-wombat%mesexcr)
      wombat%messlopdet(i,j,k) = wombat%mesgrazdet(i,j,k) * (1.0 - wombat%mesassi)*(1.0-wombat%mesexcr)
      wombat%messlopbdet(i,j,k) = wombat%mesgrazbdet(i,j,k) * (1.0 - wombat%mesassi)*(1.0-wombat%mesexcr)
      wombat%messlopzoo(i,j,k) = wombat%mesgrazzoo(i,j,k) * (1.0 - wombat%mesassi)*(1.0-wombat%mesexcr)
      messlopphyfe = wombat%messlopphy(i,j,k) * phy_Fe2C
      messlopdiafe = wombat%messlopdia(i,j,k) * dia_Fe2C
      messlopdetfe = wombat%messlopdet(i,j,k) * det_Fe2C
      messlopbdetfe = wombat%messlopbdet(i,j,k) * bdet_Fe2C
      messlopzoofe = wombat%messlopzoo(i,j,k) * zoo_Fe2C
      mesassiphyfe = wombat%mesgrazphy(i,j,k) * wombat%mesassife(i,j,k) * phy_Fe2C 
      mesassidiafe = wombat%mesgrazdia(i,j,k) * wombat%mesassife(i,j,k) * dia_Fe2C
      mesassidetfe = wombat%mesgrazdet(i,j,k) * wombat%mesassife(i,j,k) * det_Fe2C
      mesassibdetfe = wombat%mesgrazbdet(i,j,k) * wombat%mesassife(i,j,k) * bdet_Fe2C
      mesassizoofe = wombat%mesgrazzoo(i,j,k) * wombat%mesassife(i,j,k) * zoo_Fe2C
      mesexcrphyfe = wombat%mesgrazphy(i,j,k)*phy_Fe2C - mesassiphyfe - messlopphyfe
      mesexcrdiafe = wombat%mesgrazdia(i,j,k)*dia_Fe2C - mesassidiafe - messlopdiafe
      mesexcrdetfe = wombat%mesgrazdet(i,j,k)*det_Fe2C - mesassidetfe - messlopdetfe
      mesexcrbdetfe = wombat%mesgrazbdet(i,j,k)*bdet_Fe2C - mesassibdetfe - messlopbdetfe
      mesexcrzoofe = wombat%mesgrazzoo(i,j,k)*zoo_Fe2C - mesassizoofe - messlopzoofe

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

      ! dissolution
      if (wombat%f_caco3(i,j,k) .gt. epsi) then
        wombat%caldiss(i,j,k) = wombat%dissrat(i,j,k) * wombat%f_caco3(i,j,k) ! [mol/kg/s]
      else
        wombat%caldiss(i,j,k) = 0.0
      endif


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 13] Tracer tendencies                                          !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Nitrate equation ! [molN/kg]
      !----------------------------------------------------------------------
      wombat%f_no3(i,j,k) = wombat%f_no3(i,j,k) + dtsb * ( 0.0 & 
                            + wombat%ammox(i,j,k) &
                            - wombat%denitrif(i,j,k) ) &
                            + dtsb * 16./122. * ( 0.0 &
                            - wombat%phygrow(i,j,k) * wombat%phy_lno3(i,j,k) / ( wombat%phy_lnit(i,j,k) + epsi ) &
                            - wombat%diagrow(i,j,k) * wombat%dia_lno3(i,j,k) / ( wombat%dia_lnit(i,j,k) + epsi ) )
    
      ! Ammonium equation ! [molN/kg]
      !----------------------------------------------------------------------
      wombat%f_nh4(i,j,k) = wombat%f_nh4(i,j,k) + dtsb * ( 0.0 &
                            + wombat%nitrfix(i,j,k) &
                            - wombat%ammox(i,j,k) &
                            - wombat%anammox(i,j,k) ) &
                            + dtsb * 16./122. * ( 0.0 &
                            + wombat%docremi(i,j,k) &
                            + wombat%zooresp(i,j,k) &
                            + wombat%mesresp(i,j,k) &
                            - wombat%phygrow(i,j,k) * wombat%phy_lnh4(i,j,k) / ( wombat%phy_lnit(i,j,k) + epsi ) &
                            - wombat%diagrow(i,j,k) * wombat%dia_lnh4(i,j,k) / ( wombat%dia_lnit(i,j,k) + epsi ) )
                              
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
                             + wombat%zooassi * wombat%zoograzphy(i,j,k) &
                             + wombat%zooassi * wombat%zoograzdia(i,j,k) &
                             + wombat%zooassi * wombat%zoograzdet(i,j,k) &
                             - wombat%mesgrazzoo(i,j,k) &
                             - wombat%zooresp(i,j,k) &
                             - wombat%zoomort(i,j,k) )

      ! Zooplankton iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_zoofe(i,j,k)  = wombat%f_zoofe(i,j,k) + dtsb * ( 0.0 &
                               + zooassiphyfe &
                               + zooassidiafe &
                               + zooassidetfe &
                               - wombat%mesgrazzoo(i,j,k) * zoo_Fe2C &
                               - wombat%zooresp(i,j,k) * zoo_Fe2C &
                               - wombat%zoomort(i,j,k) * zoo_Fe2C )

      ! Mesozooplankton equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_mes(i,j,k)  = wombat%f_mes(i,j,k) + dtsb * ( 0.0 &
                             + wombat%mesassi * wombat%mesgrazphy(i,j,k) &
                             + wombat%mesassi * wombat%mesgrazdia(i,j,k) &
                             + wombat%mesassi * wombat%mesgrazdet(i,j,k) &
                             + wombat%mesassi * wombat%mesgrazbdet(i,j,k) &
                             + wombat%mesassi * wombat%mesgrazzoo(i,j,k) &
                             - wombat%mesresp(i,j,k) &
                             - wombat%mesmort(i,j,k) )

      ! Mesomesplankton iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_mesfe(i,j,k)  = wombat%f_mesfe(i,j,k) + dtsb * ( 0.0 &
                               + mesassiphyfe &
                               + mesassidiafe &
                               + mesassidetfe &
                               + mesassibdetfe &
                               + mesassizoofe &
                               - wombat%mesresp(i,j,k) * mes_Fe2C &
                               - wombat%mesmort(i,j,k) * mes_Fe2C )

      ! Estimate secondary productivity from zooplankton growth ! [molC/kg/s]
      wombat%zprod_gross(i,j,k) = wombat%zprod_gross(i,j,k) + dtsb * ( 0.0 &
                                  + wombat%zooassi * ( wombat%zoograzphy(i,j,k) &
                                                       + wombat%zoograzdia(i,j,k) &
                                                       + wombat%zoograzdet(i,j,k) ) &
                                  + wombat%mesassi * ( wombat%mesgrazphy(i,j,k) &
                                                       + wombat%mesgrazdia(i,j,k) &
                                                       + wombat%mesgrazdet(i,j,k) &
                                                       + wombat%mesgrazbdet(i,j,k) &
                                                       + wombat%mesgrazzoo(i,j,k) ) )

      ! Detritus equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_det(i,j,k) = wombat%f_det(i,j,k) + dtsb * ( 0.0 & 
                            + wombat%zooslopphy(i,j,k) &
                            + wombat%zooslopdia(i,j,k) &
                            + wombat%zooslopdet(i,j,k) &
                            + wombat%messlopphy(i,j,k) &
                            + wombat%messlopdia(i,j,k) &
                            + wombat%messlopdet(i,j,k) &
                            + wombat%messlopbdet(i,j,k) &
                            + wombat%messlopzoo(i,j,k) &
                            + wombat%phymort(i,j,k) &
                            + wombat%zoomort(i,j,k) &
                            - wombat%zoograzdet(i,j,k) &
                            - wombat%mesgrazdet(i,j,k) &
                            - wombat%detremi(i,j,k) )

      ! Detrital iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_detfe(i,j,k) = wombat%f_detfe(i,j,k) + dtsb * ( 0.0 & 
                              + zooslopphyfe &
                              + zooslopdiafe &
                              + zooslopdetfe &
                              + messlopphyfe &
                              + messlopdiafe &
                              + messlopdetfe &
                              + messlopbdetfe &
                              + messlopzoofe &
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
                             - wombat%mesgrazbdet(i,j,k) &
                             - wombat%bdetremi(i,j,k) )

      ! Big detrital iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_bdetfe(i,j,k) = wombat%f_bdetfe(i,j,k) + dtsb * ( 0.0 &
                               + wombat%diamort(i,j,k) * dia_Fe2C &
                               + wombat%mesmort(i,j,k) * mes_Fe2C &
                               - wombat%mesgrazbdet(i,j,k) * bdet_Fe2C &
                               - wombat%bdetremi(i,j,k) * bdet_Fe2C &
                               + wombat%fescabdet(i,j,k) &
                               + wombat%fecoag2bdet(i,j,k) )
      
      ! Dissolved organic carbon equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_doc(i,j,k) = wombat%f_doc(i,j,k) + dtsb * ( 0.0 &
                            + wombat%detremi(i,j,k) &
                            + wombat%bdetremi(i,j,k) &
                            + wombat%phyresp(i,j,k) &
                            + wombat%diaresp(i,j,k) &
                            + wombat%zooexcrphy(i,j,k) &
                            + wombat%zooexcrdia(i,j,k) &
                            + wombat%zooexcrdet(i,j,k) &
                            + wombat%mesexcrphy(i,j,k) &
                            + wombat%mesexcrdia(i,j,k) &
                            + wombat%mesexcrdet(i,j,k) &
                            + wombat%mesexcrbdet(i,j,k) &
                            + wombat%mesexcrzoo(i,j,k) &
                            - wombat%docremi(i,j,k) )

      ! Oxygen equation ! [molO2/kg]
      !-----------------------------------------------------------------------
      if (wombat%f_o2(i,j,k) .gt. epsi) &
        wombat%f_o2(i,j,k) = wombat%f_o2(i,j,k) - 132./122. * dtsb * ( 0.0 &
                             + wombat%docremi(i,j,k) * (1.0 - wombat%fdenitrif(i,j,k)) &
                             + wombat%zooresp(i,j,k) &
                             + wombat%mesresp(i,j,k) &
                             - wombat%phygrow(i,j,k) &
                             - wombat%diagrow(i,j,k) ) &
                             - 40./16. * dtsb * ( wombat%ammox(i,j,k) )


      ! Equation for CaCO3 ! [molCaCO3/kg]
      !-----------------------------------------------------------------------
      wombat%f_caco3(i,j,k) = wombat%f_caco3(i,j,k) + dtsb * ( 0.0 &
                              + wombat%zooslopphy(i,j,k) * wombat%pic2poc(i,j,k) &
                              + wombat%messlopphy(i,j,k) * wombat%pic2poc(i,j,k) &
                              + wombat%phymort(i,j,k) * wombat%pic2poc(i,j,k) &
                              + wombat%zoomort(i,j,k) * wombat%pic2poc(i,j,k) &
                              - wombat%caldiss(i,j,k) )

      ! Equation for DIC ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_dic(i,j,k) = wombat%f_dic(i,j,k) + dtsb * ( 0.0 &
                            + wombat%docremi(i,j,k) &
                            + wombat%zooresp(i,j,k) &
                            + wombat%mesresp(i,j,k) &
                            - wombat%phygrow(i,j,k) &
                            - wombat%diagrow(i,j,k) &
                            - wombat%zooslopphy(i,j,k) * wombat%pic2poc(i,j,k) &
                            - wombat%messlopphy(i,j,k) * wombat%pic2poc(i,j,k) &
                            - wombat%phymort(i,j,k) * wombat%pic2poc(i,j,k) &
                            - wombat%zoomort(i,j,k) * wombat%pic2poc(i,j,k) &  
                            + wombat%caldiss(i,j,k) )

      ! Equation for DICr ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_dicr(i,j,k) = wombat%f_dicr(i,j,k) + dtsb * ( 0.0 &
                             + wombat%docremi(i,j,k) &
                             + wombat%zooresp(i,j,k) &
                             + wombat%mesresp(i,j,k) &
                             - wombat%phygrow(i,j,k) &
                             - wombat%diagrow(i,j,k) &
                             - wombat%zooslopphy(i,j,k) * wombat%pic2poc(i,j,k) &
                             - wombat%messlopphy(i,j,k) * wombat%pic2poc(i,j,k) &
                             - wombat%phymort(i,j,k) * wombat%pic2poc(i,j,k) &
                             - wombat%zoomort(i,j,k) * wombat%pic2poc(i,j,k) &  
                             + wombat%caldiss(i,j,k) )

      ! Equation for ALK ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_alk(i,j,k) = wombat%f_alk(i,j,k) + dtsb * 16.0/122.0 * ( 0.0 &
                            + wombat%phygrow(i,j,k) * wombat%phy_lno3(i,j,k) / ( wombat%phy_lnit(i,j,k) + epsi ) &
                            + wombat%diagrow(i,j,k) * wombat%dia_lno3(i,j,k) / ( wombat%dia_lnit(i,j,k) + epsi ) &
                            + wombat%docremi(i,j,k) &
                            + wombat%zooresp(i,j,k) &
                            + wombat%mesresp(i,j,k) &
                            - wombat%phygrow(i,j,k) * wombat%phy_lnh4(i,j,k) / ( wombat%phy_lnit(i,j,k) + epsi ) &
                            - wombat%diagrow(i,j,k) * wombat%dia_lnh4(i,j,k) / ( wombat%dia_lnit(i,j,k) + epsi ) ) &
                            + dtsb * ( 0.0 &
                            + wombat%denitrif(i,j,k) &
                            - wombat%anammox(i,j,k) &
                            - 2.0 * wombat%ammox(i,j,k) ) &
                            + dtsb * 2.0 * ( 0.0 &
                            + wombat%caldiss(i,j,k) &
                            - wombat%zooslopphy(i,j,k) * wombat%pic2poc(i,j,k) &
                            - wombat%messlopphy(i,j,k) * wombat%pic2poc(i,j,k) &
                            - wombat%phymort(i,j,k) * wombat%pic2poc(i,j,k) &
                            - wombat%zoomort(i,j,k) * wombat%pic2poc(i,j,k) )

      ! Extra equation for iron ! [molFe/kg]
      !----------------------------------------------------------------------
      wombat%f_fe(i,j,k) = wombat%f_fe(i,j,k) + dtsb * ( 0.0 &
                           + wombat%detremi(i,j,k) * det_Fe2C &
                           + wombat%bdetremi(i,j,k) * bdet_Fe2C &
                           + wombat%zooresp(i,j,k) * zoo_Fe2C &
                           + wombat%mesresp(i,j,k) * mes_Fe2C &
                           + zooexcrphyfe &
                           + zooexcrdiafe &
                           + zooexcrdetfe &
                           + mesexcrphyfe &
                           + mesexcrdiafe &
                           + mesexcrdetfe &
                           + mesexcrbdetfe &
                           + mesexcrzoofe &
                           + wombat%phyresp(i,j,k) * phy_Fe2C &
                           + wombat%diaresp(i,j,k) * dia_Fe2C &
                           - wombat%phy_dfeupt(i,j,k) &
                           - wombat%dia_dfeupt(i,j,k) &
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
                                + zooexcrphyfe &
                                + zooexcrdiafe &
                                + zooexcrdetfe &
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
                              + wombat%feprecip(i,j,k) &
                              + wombat%fescaven(i,j,k) &
                              + wombat%fecoag2det(i,j,k) &
                              + wombat%fecoag2bdet(i,j,k)) 


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 14] Check for conservation of mass by ecosystem component      !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      n_pools(i,j,k,2) = wombat%f_no3(i,j,k) + wombat%f_nh4(i,j,k) + (wombat%f_phy(i,j,k) + wombat%f_det(i,j,k) + &
                          wombat%f_bdet(i,j,k) + wombat%f_zoo(i,j,k) + wombat%f_mes(i,j,k) + wombat%f_dia(i,j,k) + &
                          wombat%f_doc(i,j,k)) * 16/122.0
      c_pools(i,j,k,2) = wombat%f_dic(i,j,k) + wombat%f_phy(i,j,k) + wombat%f_det(i,j,k) + wombat%f_bdet(i,j,k) + &
                         wombat%f_zoo(i,j,k) + wombat%f_mes(i,j,k) + wombat%f_caco3(i,j,k) + wombat%f_dia(i,j,k) + &
                         wombat%f_doc(i,j,k)

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
          print *, " "
          print *, "       NO3 (molNO3/kg) =", wombat%f_no3(i,j,k)
          print *, "       NH4 (molNH4/kg) =", wombat%f_nh4(i,j,k)
          print *, "       PHY (molN/kg) =", wombat%f_phy(i,j,k) * 16.0 / 122.0
          print *, "       DIA (molN/kg) =", wombat%f_dia(i,j,k) * 16.0 / 122.0
          print *, "       ZOO (molN/kg) =", wombat%f_zoo(i,j,k) * 16.0 / 122.0
          print *, "       MES (molN/kg) =", wombat%f_mes(i,j,k) * 16.0 / 122.0
          print *, "       DET (molN/kg) =", wombat%f_det(i,j,k) * 16.0 / 122.0
          print *, "       BDET (molN/kg) =", wombat%f_bdet(i,j,k) * 16.0 / 122.0
          print *, "       DOC (molN/kg) =", wombat%f_doc(i,j,k) * 16.0 / 122.0
          print *, " "
          print *, "       ammox (molN/kg/s) =", wombat%ammox(i,j,k)
          print *, "       phygrow (molC/kg/s) =", wombat%phygrow(i,j,k)
          print *, "       diagrow (molC/kg/s) =", wombat%diagrow(i,j,k)
          print *, "       detremi (molC/kg/s) =", wombat%detremi(i,j,k)
          print *, "       bdetremi (molC/kg/s) =", wombat%bdetremi(i,j,k)
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
          print *, " "
          print *, "       DIC (molC/kg) =", wombat%f_dic(i,j,k)
          print *, "       ALK (molC/kg) =", wombat%f_alk(i,j,k)
          print *, "       PHY (molC/kg) =", wombat%f_phy(i,j,k)
          print *, "       DIA (molC/kg) =", wombat%f_dia(i,j,k)
          print *, "       ZOO (molN/kg) =", wombat%f_zoo(i,j,k)
          print *, "       MES (molN/kg) =", wombat%f_mes(i,j,k)
          print *, "       DET (molN/kg) =", wombat%f_det(i,j,k)
          print *, "       BDET (molN/kg) =", wombat%f_bdet(i,j,k)
          print *, "       DOC (molN/kg) =", wombat%f_doc(i,j,k)
          print *, "       CaCO3 (molC/kg) =", wombat%f_caco3(i,j,k)
          print *, "       Temp =", Temp(i,j,k)
          print *, "       Salt =", Salt(i,j,k)
          print *, "       surface pCO2 =", wombat%pco2_csurf(i,j)
          print *, "       htotal =", wombat%htotal(i,j,k)
          print *, " "
          print *, "       phygrow (molC/kg/s) =", wombat%phygrow(i,j,k)
          print *, "       diagrow (molC/kg/s) =", wombat%diagrow(i,j,k)
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
          print *, "       zooslopphy * pic2poc(i,j,k) (molC/kg/s) =", wombat%zooslopphy(i,j,k) * wombat%pic2poc(i,j,k)
          print *, "       messlopphy * pic2poc(i,j,k) (molC/kg/s) =", wombat%messlopphy(i,j,k) * wombat%pic2poc(i,j,k)
          print *, "       phymort * pic2poc(i,j,k) (molC/kg/s) =", wombat%phymort(i,j,k) * wombat%pic2poc(i,j,k)
          print *, "       diamort * pic2poc(i,j,k) (molC/kg/s) =", wombat%diamort(i,j,k) * wombat%pic2poc(i,j,k)
          print *, "       zoomort * pic2poc(i,j,k) (molC/kg/s) =", wombat%zoomort(i,j,k) * wombat%pic2poc(i,j,k)
          print *, "       caldiss (molC/kg/s) =", wombat%caldiss(i,j,k)
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
    call g_tracer_set_values(tracer_list, 'phy', 'field', wombat%f_phy, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'pchl', 'field', wombat%f_pchl, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'phyfe', 'field', wombat%f_phyfe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'dia', 'field', wombat%f_dia, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'dchl', 'field', wombat%f_dchl, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'diafe', 'field', wombat%f_diafe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'zoo', 'field', wombat%f_zoo, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'zoofe', 'field', wombat%f_zoofe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'mes', 'field', wombat%f_mes, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'mesfe', 'field', wombat%f_mesfe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'det', 'field', wombat%f_det, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'detfe', 'field', wombat%f_detfe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'bdet', 'field', wombat%f_bdet, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'bdetfe', 'field', wombat%f_bdetfe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'doc', 'field', wombat%f_doc, isd, jsd, ntau=tau)
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
        wombat%p_wcaco3(i,j,:) = wsinkcal(:)
      else
        wombat%p_wdet(i,j,:) = 0.0
        wombat%p_wdetfe(i,j,:) = 0.0
        wombat%p_wbdet(i,j,:) = 0.0
        wombat%p_wbdetfe(i,j,:) = 0.0
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
    call g_tracer_get_pointer(tracer_list, 'caco3_sediment', 'field', wombat%p_caco3_sediment) ! [mol/m2]

    do j = jsc,jec; do i = isc,iec;
      k = grid_kmt(i,j)
      if (k .gt. 0) then
      ! Get bottom conditions: mask, PO4, temp, salt, DIC, Alk, H+ ion
        wombat%seddep(i,j) = wombat%zw(i,j,k)
        wombat%sedmask(i,j) = grid_tmask(i,j,k)
        wombat%sedtemp(i,j) = Temp(i,j,k)
        wombat%sedsalt(i,j) = Salt(i,j,k)
        wombat%sedno3(i,j) = wombat%f_no3(i,j,k)
        wombat%sednh4(i,j) = wombat%f_nh4(i,j,k)
        wombat%sedo2(i,j) = wombat%f_o2(i,j,k)
        wombat%seddic(i,j) = wombat%f_dic(i,j,k) + wombat%p_det_sediment(i,j,1) / wombat%Rho_0  ![mol/kg] 
        wombat%sedalk(i,j) = wombat%f_alk(i,j,k)
        wombat%sedhtotal(i,j) = wombat%htotal(i,j,k)
      endif
    enddo; enddo

    call FMS_ocmip2_co2calc(CO2_dope_vec, wombat%sedmask(:,:), &
        wombat%sedtemp(:,:), wombat%sedsalt(:,:), &
        min(wombat%dic_max*mmol_m3_to_mol_kg, max(wombat%seddic(:,:), wombat%dic_min*mmol_m3_to_mol_kg)), &
        max(wombat%sedno3(:,:) / 16., 1e-9), &
        wombat%sio2(:,:), &
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
      if (do_caco3_dynamics) then
        wombat%caco3_sed_remin(i,j) = wombat%caco3lrem_sed * fbc * wombat%p_caco3_sediment(i,j,1) * &
                                            max(0.1, (1.0 - wombat%sedomega_cal(i,j)))**(4.5)
      else
        wombat%caco3_sed_remin(i,j) = wombat%caco3lrem_sed * fbc * wombat%p_caco3_sediment(i,j,1) * &
                                            max(0.1, (1.0 - 0.2081))**(4.5)
      endif
      if (do_open_n_cycle) then
        ! sedimentary denitrification (Bohlen et al., 2012 Global Biogeochemical Cycles)
        !   Stoichiometry of 94 mol NO3 used per 122 mol organic carbon oxidised (Paulmier et al, 2009 BG)
        !   Hard limit where denitrification is at maximum 90% responsible for organic matter remin
        wombat%det_sed_denit(i,j) = wombat%det_sed_remin(i,j) * min(0.9 * 94.0/122.0, &
                                    (0.083 + 0.21 * 0.98**((wombat%sedo2(i,j) - wombat%sedno3(i,j))/mmol_m3_to_mol_kg)))
        wombat%fdenit(i,j) = wombat%det_sed_denit(i,j) * 122.0/94.0 / (wombat%det_sed_remin(i,j) + epsi)
      endif

      ! Remineralisation of sediments to supply nutrient fields.
      ! btf values are positive from the water column into the sediment.
      wombat%b_nh4(i,j) = -16./122. * wombat%det_sed_remin(i,j) ! [mol/m2/s]
      wombat%b_no3(i,j) = wombat%det_sed_denit(i,j) ! [molN/m2/s]
      wombat%b_o2(i,j) = -132./16. * wombat%b_nh4(i,j) * (1.0 - wombat%fdenit(i,j))! [mol/m2/s]
      wombat%b_dic(i,j) = 122./16. * wombat%b_nh4(i,j) - wombat%caco3_sed_remin(i,j) ! [mol/m2/s]
      wombat%b_dicr(i,j) = wombat%b_dic(i,j) ! [mol/m2/s]
      wombat%b_fe(i,j) = -1.0 * wombat%detfe_sed_remin(i,j) ! [mol/m2/s]
      wombat%b_alk(i,j) = -2.0 * wombat%caco3_sed_remin(i,j) + wombat%b_nh4(i,j) - wombat%b_no3(i,j) ! [mol/m2/s]
    enddo; enddo


    ! Apply remineralisation rates to sediment tracers
    !-----------------------------------------------------------------------
    do j = jsc,jec; do i = isc,iec;

      if (grid_kmt(i,j) .gt. 0) then
        wombat%p_det_sediment(i,j,1) = wombat%p_det_sediment(i,j,1) - dt * wombat%det_sed_remin(i,j) ! [mol/m2]
        wombat%p_detfe_sediment(i,j,1) = wombat%p_detfe_sediment(i,j,1) - dt * wombat%detfe_sed_remin(i,j) ! [mol/m2]
        wombat%p_caco3_sediment(i,j,1) = wombat%p_caco3_sediment(i,j,1) - dt * wombat%caco3_sed_remin(i,j) ! [mol/m2]
      endif
    enddo; enddo

    call g_tracer_set_values(tracer_list, 'nh4', 'btf', wombat%b_nh4, isd, jsd)
    call g_tracer_set_values(tracer_list, 'no3', 'btf', wombat%b_no3, isd, jsd)
    call g_tracer_set_values(tracer_list, 'o2', 'btf', wombat%b_o2, isd, jsd)
    call g_tracer_set_values(tracer_list, 'dic', 'btf', wombat%b_dic, isd, jsd)
    call g_tracer_set_values(tracer_list, 'dicr', 'btf', wombat%b_dicr, isd, jsd)
    call g_tracer_set_values(tracer_list, 'fe', 'btf', wombat%b_fe, isd, jsd)
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

    if (wombat%id_diaresp .gt. 0) &
      used = g_send_data(wombat%id_diaresp, wombat%diaresp, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_diamort .gt. 0) &
      used = g_send_data(wombat%id_diamort, wombat%diamort, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooeps .gt. 0) &
      used = g_send_data(wombat%id_zooeps, wombat%zooeps, model_time, &
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

    if (wombat%id_zooexcrphy .gt. 0) &
      used = g_send_data(wombat%id_zooexcrphy, wombat%zooexcrphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcrdia .gt. 0) &
      used = g_send_data(wombat%id_zooexcrdia, wombat%zooexcrdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcrdet .gt. 0) &
      used = g_send_data(wombat%id_zooexcrdet, wombat%zooexcrdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooslopphy .gt. 0) &
      used = g_send_data(wombat%id_zooslopphy, wombat%zooslopphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooslopdia .gt. 0) &
      used = g_send_data(wombat%id_zooslopdia, wombat%zooslopdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooslopdet .gt. 0) &
      used = g_send_data(wombat%id_zooslopdet, wombat%zooslopdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooassife .gt. 0) &
      used = g_send_data(wombat%id_zooassife, wombat%zooassife, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_meseps .gt. 0) &
      used = g_send_data(wombat%id_meseps, wombat%meseps, model_time, &
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

    if (wombat%id_messlopphy .gt. 0) &
      used = g_send_data(wombat%id_messlopphy, wombat%messlopphy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_messlopdia .gt. 0) &
      used = g_send_data(wombat%id_messlopdia, wombat%messlopdia, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_messlopdet .gt. 0) &
      used = g_send_data(wombat%id_messlopdet, wombat%messlopdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_messlopbdet .gt. 0) &
      used = g_send_data(wombat%id_messlopbdet, wombat%messlopbdet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_messlopzoo .gt. 0) &
      used = g_send_data(wombat%id_messlopzoo, wombat%messlopzoo, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_mesassife .gt. 0) &
      used = g_send_data(wombat%id_mesassife, wombat%mesassife, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_reminr .gt. 0) &
      used = g_send_data(wombat%id_reminr, wombat%reminr, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_docremi .gt. 0) &
      used = g_send_data(wombat%id_docremi, wombat%docremi, model_time, &
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

    if (wombat%id_aoa_mu .gt. 0) &
      used = g_send_data(wombat%id_aoa_mu, wombat%aoa_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    
    if (wombat%id_het_loxy .gt. 0) &
      used = g_send_data(wombat%id_het_loxy, wombat%het_loxy, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_het_lno3 .gt. 0) &
      used = g_send_data(wombat%id_het_lno3, wombat%het_lno3, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_het_ldet .gt. 0) &
      used = g_send_data(wombat%id_het_ldet, wombat%het_ldet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_het_mu .gt. 0) &
      used = g_send_data(wombat%id_het_mu, wombat%het_mu, model_time, &
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

    if (wombat%id_denitrif .gt. 0) &
      used = g_send_data(wombat%id_denitrif, wombat%denitrif, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fdenitrif .gt. 0) &
      used = g_send_data(wombat%id_fdenitrif, wombat%fdenitrif, model_time, &
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
    real                                    :: sal, ST, o2_saturation
    real                                    :: tt, tk, ts, ts2, ts3, ts4, ts5
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
      ! Compute the Schmidt number of CO2 in seawater using the formulation
      ! presented by Wanninkhof (1992, J. Geophys. Res., 97, 7373-7382).
      !-----------------------------------------------------------------------
      ST = SST(i,j)
      wombat%co2_sc_no(i,j) = wombat%a1_co2 + ST*(wombat%a2_co2 + ST*(wombat%a3_co2 + &
          ST*wombat%a4_co2)) * grid_tmask(i,j,1)
      
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
      !
      ! o2_saturation is defined between T(freezing) <= T <= 40 deg C and
      !                                   0 permil <= S <= 42 permil
      ! We impose these bounds here.
      !
      ! check value: T = 10 deg C, S = 35 permil,
      !              o2_saturation = 0.282015 mol m-3
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

      o2_saturation = (1000.0/22391.6) * grid_tmask(i,j,1) *  & !convert from ml/l to mol m-3
          exp( wombat%a_0 + wombat%a_1*ts + wombat%a_2*ts2 + wombat%a_3*ts3 + wombat%a_4*ts4 + &
              wombat%a_5*ts5 + (wombat%b_0 + wombat%b_1*ts + wombat%b_2*ts2 + wombat%b_3*ts3 + &
              wombat%c_0*sal)*sal)

      !-----------------------------------------------------------------------
      !  Compute the Schmidt number of O2 in seawater using the
      !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
      !  Cycles, 12, 141-163).
      !-----------------------------------------------------------------------
      wombat%o2_sc_no(i,j)  = wombat%a1_o2 + ST * (wombat%a2_o2 + ST * (wombat%a3_o2 + ST * &
          wombat%a4_o2 )) * grid_tmask(i,j,1)

      ! renormalize the alpha value for atm o2
      ! data table override for o2_flux_pcair_atm is now set to 0.21
      wombat%o2_alpha(i,j) = (o2_saturation / 0.21)
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
    allocate(wombat%sio2(isd:ied, jsd:jed)); wombat%sio2(:,:)=wombat%sio2_surf
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
    allocate(wombat%f_phy(isd:ied, jsd:jed, 1:nk)); wombat%f_phy(:,:,:)=0.0
    allocate(wombat%f_pchl(isd:ied, jsd:jed, 1:nk)); wombat%f_pchl(:,:,:)=0.0
    allocate(wombat%f_phyfe(isd:ied, jsd:jed, 1:nk)); wombat%f_phyfe(:,:,:)=0.0
    allocate(wombat%f_dia(isd:ied, jsd:jed, 1:nk)); wombat%f_dia(:,:,:)=0.0
    allocate(wombat%f_dchl(isd:ied, jsd:jed, 1:nk)); wombat%f_dchl(:,:,:)=0.0
    allocate(wombat%f_diafe(isd:ied, jsd:jed, 1:nk)); wombat%f_diafe(:,:,:)=0.0
    allocate(wombat%f_zoo(isd:ied, jsd:jed, 1:nk)); wombat%f_zoo(:,:,:)=0.0
    allocate(wombat%f_zoofe(isd:ied, jsd:jed, 1:nk)); wombat%f_zoofe(:,:,:)=0.0
    allocate(wombat%f_mes(isd:ied, jsd:jed, 1:nk)); wombat%f_mes(:,:,:)=0.0
    allocate(wombat%f_mesfe(isd:ied, jsd:jed, 1:nk)); wombat%f_mesfe(:,:,:)=0.0
    allocate(wombat%f_det(isd:ied, jsd:jed, 1:nk)); wombat%f_det(:,:,:)=0.0
    allocate(wombat%f_detfe(isd:ied, jsd:jed, 1:nk)); wombat%f_detfe(:,:,:)=0.0
    allocate(wombat%f_bdet(isd:ied, jsd:jed, 1:nk)); wombat%f_bdet(:,:,:)=0.0
    allocate(wombat%f_bdetfe(isd:ied, jsd:jed, 1:nk)); wombat%f_bdetfe(:,:,:)=0.0
    allocate(wombat%f_doc(isd:ied, jsd:jed, 1:nk)); wombat%f_doc(:,:,:)=0.0
    allocate(wombat%f_o2(isd:ied, jsd:jed, 1:nk)); wombat%f_o2(:,:,:)=0.0
    allocate(wombat%f_caco3(isd:ied, jsd:jed, 1:nk)); wombat%f_caco3(:,:,:)=0.0
    allocate(wombat%f_fe(isd:ied, jsd:jed, 1:nk)); wombat%f_fe(:,:,:)=0.0

    allocate(wombat%b_nh4(isd:ied, jsd:jed)); wombat%b_nh4(:,:)=0.0
    allocate(wombat%b_no3(isd:ied, jsd:jed)); wombat%b_no3(:,:)=0.0
    allocate(wombat%b_o2(isd:ied, jsd:jed)); wombat%b_o2(:,:)=0.0
    allocate(wombat%b_dic(isd:ied, jsd:jed)); wombat%b_dic(:,:)=0.0
    allocate(wombat%b_dicr(isd:ied, jsd:jed)); wombat%b_dicr(:,:)=0.0
    allocate(wombat%b_fe(isd:ied, jsd:jed)); wombat%b_fe(:,:)=0.0
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
    allocate(wombat%phyresp(isd:ied, jsd:jed, 1:nk)); wombat%phyresp(:,:,:)=0.0
    allocate(wombat%phymort(isd:ied, jsd:jed, 1:nk)); wombat%phymort(:,:,:)=0.0
    allocate(wombat%dia_feupreg(isd:ied, jsd:jed, 1:nk)); wombat%dia_feupreg(:,:,:)=0.0
    allocate(wombat%dia_fedoreg(isd:ied, jsd:jed, 1:nk)); wombat%dia_fedoreg(:,:,:)=0.0
    allocate(wombat%diagrow(isd:ied, jsd:jed, 1:nk)); wombat%diagrow(:,:,:)=0.0
    allocate(wombat%diaresp(isd:ied, jsd:jed, 1:nk)); wombat%diaresp(:,:,:)=0.0
    allocate(wombat%diamort(isd:ied, jsd:jed, 1:nk)); wombat%diamort(:,:,:)=0.0
    allocate(wombat%zooeps(isd:ied, jsd:jed, 1:nk)); wombat%zooeps(:,:,:)=0.0
    allocate(wombat%zoograzphy(isd:ied, jsd:jed, 1:nk)); wombat%zoograzphy(:,:,:)=0.0
    allocate(wombat%zoograzdia(isd:ied, jsd:jed, 1:nk)); wombat%zoograzdia(:,:,:)=0.0
    allocate(wombat%zoograzdet(isd:ied, jsd:jed, 1:nk)); wombat%zoograzdet(:,:,:)=0.0
    allocate(wombat%zooresp(isd:ied, jsd:jed, 1:nk)); wombat%zooresp(:,:,:)=0.0
    allocate(wombat%zoomort(isd:ied, jsd:jed, 1:nk)); wombat%zoomort(:,:,:)=0.0
    allocate(wombat%zooexcrphy(isd:ied, jsd:jed, 1:nk)); wombat%zooexcrphy(:,:,:)=0.0
    allocate(wombat%zooexcrdia(isd:ied, jsd:jed, 1:nk)); wombat%zooexcrdia(:,:,:)=0.0
    allocate(wombat%zooexcrdet(isd:ied, jsd:jed, 1:nk)); wombat%zooexcrdet(:,:,:)=0.0
    allocate(wombat%zooslopphy(isd:ied, jsd:jed, 1:nk)); wombat%zooslopphy(:,:,:)=0.0
    allocate(wombat%zooslopdia(isd:ied, jsd:jed, 1:nk)); wombat%zooslopdia(:,:,:)=0.0
    allocate(wombat%zooslopdet(isd:ied, jsd:jed, 1:nk)); wombat%zooslopdet(:,:,:)=0.0
    allocate(wombat%zooassife(isd:ied, jsd:jed, 1:nk)); wombat%zooassife(:,:,:)=0.0
    allocate(wombat%meseps(isd:ied, jsd:jed, 1:nk)); wombat%meseps(:,:,:)=0.0
    allocate(wombat%mesgrazphy(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazphy(:,:,:)=0.0
    allocate(wombat%mesgrazdia(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazdia(:,:,:)=0.0
    allocate(wombat%mesgrazdet(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazdet(:,:,:)=0.0
    allocate(wombat%mesgrazbdet(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazbdet(:,:,:)=0.0
    allocate(wombat%mesgrazzoo(isd:ied, jsd:jed, 1:nk)); wombat%mesgrazzoo(:,:,:)=0.0
    allocate(wombat%mesresp(isd:ied, jsd:jed, 1:nk)); wombat%mesresp(:,:,:)=0.0
    allocate(wombat%mesmort(isd:ied, jsd:jed, 1:nk)); wombat%mesmort(:,:,:)=0.0
    allocate(wombat%mesexcrphy(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrphy(:,:,:)=0.0
    allocate(wombat%mesexcrdia(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrdia(:,:,:)=0.0
    allocate(wombat%mesexcrdet(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrdet(:,:,:)=0.0
    allocate(wombat%mesexcrbdet(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrbdet(:,:,:)=0.0
    allocate(wombat%mesexcrzoo(isd:ied, jsd:jed, 1:nk)); wombat%mesexcrzoo(:,:,:)=0.0
    allocate(wombat%messlopphy(isd:ied, jsd:jed, 1:nk)); wombat%messlopphy(:,:,:)=0.0
    allocate(wombat%messlopdia(isd:ied, jsd:jed, 1:nk)); wombat%messlopdia(:,:,:)=0.0
    allocate(wombat%messlopdet(isd:ied, jsd:jed, 1:nk)); wombat%messlopdet(:,:,:)=0.0
    allocate(wombat%messlopbdet(isd:ied, jsd:jed, 1:nk)); wombat%messlopbdet(:,:,:)=0.0
    allocate(wombat%messlopzoo(isd:ied, jsd:jed, 1:nk)); wombat%messlopzoo(:,:,:)=0.0
    allocate(wombat%mesassife(isd:ied, jsd:jed, 1:nk)); wombat%mesassife(:,:,:)=0.0
    allocate(wombat%reminr(isd:ied, jsd:jed, 1:nk)); wombat%reminr(:,:,:)=0.0
    allocate(wombat%docremi(isd:ied, jsd:jed, 1:nk)); wombat%docremi(:,:,:)=0.0
    allocate(wombat%detremi(isd:ied, jsd:jed, 1:nk)); wombat%detremi(:,:,:)=0.0
    allocate(wombat%bdetremi(isd:ied, jsd:jed, 1:nk)); wombat%bdetremi(:,:,:)=0.0
    allocate(wombat%pic2poc(isd:ied, jsd:jed, 1:nk)); wombat%pic2poc(:,:,:)=0.0
    allocate(wombat%dissrat(isd:ied, jsd:jed, 1:nk)); wombat%dissrat(:,:,:)=0.0
    allocate(wombat%caldiss(isd:ied, jsd:jed, 1:nk)); wombat%caldiss(:,:,:)=0.0
    allocate(wombat%aoa_loxy(isd:ied, jsd:jed, 1:nk)); wombat%aoa_loxy(:,:,:)=0.0
    allocate(wombat%aoa_lnh4(isd:ied, jsd:jed, 1:nk)); wombat%aoa_lnh4(:,:,:)=0.0
    allocate(wombat%aoa_mu(isd:ied, jsd:jed, 1:nk)); wombat%aoa_mu(:,:,:)=0.0
    allocate(wombat%het_loxy(isd:ied, jsd:jed, 1:nk)); wombat%het_loxy(:,:,:)=0.0
    allocate(wombat%het_lno3(isd:ied, jsd:jed, 1:nk)); wombat%het_lno3(:,:,:)=0.0
    allocate(wombat%het_ldet(isd:ied, jsd:jed, 1:nk)); wombat%het_ldet(:,:,:)=0.0
    allocate(wombat%het_mu(isd:ied, jsd:jed, 1:nk)); wombat%het_mu(:,:,:)=0.0
    allocate(wombat%aox_lnh4(isd:ied, jsd:jed, 1:nk)); wombat%aox_lnh4(:,:,:)=0.0
    allocate(wombat%aox_mu(isd:ied, jsd:jed, 1:nk)); wombat%aox_mu(:,:,:)=0.0
    allocate(wombat%nitrfix(isd:ied, jsd:jed, 1:nk)); wombat%nitrfix(:,:,:)=0.0
    allocate(wombat%ammox(isd:ied, jsd:jed, 1:nk)); wombat%ammox(:,:,:)=0.0
    allocate(wombat%anammox(isd:ied, jsd:jed, 1:nk)); wombat%anammox(:,:,:)=0.0
    allocate(wombat%denitrif(isd:ied, jsd:jed, 1:nk)); wombat%denitrif(:,:,:)=0.0
    allocate(wombat%fdenitrif(isd:ied, jsd:jed, 1:nk)); wombat%fdenitrif(:,:,:)=0.0
    allocate(wombat%no3_prev(isd:ied, jsd:jed, 1:nk)); wombat%no3_prev(:,:,:)=0.0
    allocate(wombat%caco3_prev(isd:ied, jsd:jed, 1:nk)); wombat%caco3_prev(:,:,:)=0.0
    allocate(wombat%det_sed_remin(isd:ied, jsd:jed)); wombat%det_sed_remin(:,:)=0.0
    allocate(wombat%det_sed_denit(isd:ied, jsd:jed)); wombat%det_sed_denit(:,:)=0.0
    allocate(wombat%det_btm(isd:ied, jsd:jed)); wombat%det_btm(:,:)=0.0
    allocate(wombat%bdet_btm(isd:ied, jsd:jed)); wombat%bdet_btm(:,:)=0.0
    allocate(wombat%fbury(isd:ied, jsd:jed)); wombat%fbury(:,:)=0.0
    allocate(wombat%fdenit(isd:ied, jsd:jed)); wombat%fdenit(:,:)=0.0
    allocate(wombat%detfe_sed_remin(isd:ied, jsd:jed)); wombat%detfe_sed_remin(:,:)=0.0
    allocate(wombat%detfe_btm(isd:ied, jsd:jed)); wombat%detfe_btm(:,:)=0.0
    allocate(wombat%bdetfe_btm(isd:ied, jsd:jed)); wombat%bdetfe_btm(:,:)=0.0
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
        wombat%f_phy, &
        wombat%f_pchl, &
        wombat%f_phyfe, &
        wombat%f_dia, &
        wombat%f_dchl, &
        wombat%f_diafe, &
        wombat%f_zoo, &
        wombat%f_zoofe, &
        wombat%f_mes, &
        wombat%f_mesfe, &
        wombat%f_det, &
        wombat%f_detfe, &
        wombat%f_bdet, &
        wombat%f_bdetfe, &
        wombat%f_doc, &
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
        wombat%phyresp, &
        wombat%phymort, &
        wombat%dia_feupreg, &
        wombat%dia_fedoreg, &
        wombat%diagrow, &
        wombat%diaresp, &
        wombat%diamort, &
        wombat%zoograzphy, &
        wombat%zoograzdia, &
        wombat%zoograzdet, &
        wombat%zooresp, &
        wombat%zoomort, &
        wombat%zooexcrphy, &
        wombat%zooexcrdia, &
        wombat%zooexcrdet, &
        wombat%zooslopphy, &
        wombat%zooslopdia, &
        wombat%zooslopdet, &
        wombat%zooassife, &
        wombat%mesgrazphy, &
        wombat%mesgrazdia, &
        wombat%mesgrazdet, &
        wombat%mesgrazbdet, &
        wombat%mesgrazzoo, &
        wombat%mesresp, &
        wombat%mesmort, &
        wombat%mesexcrphy, &
        wombat%mesexcrdia, &
        wombat%mesexcrdet, &
        wombat%mesexcrbdet, &
        wombat%mesexcrzoo, &
        wombat%messlopphy, &
        wombat%messlopdia, &
        wombat%messlopdet, &
        wombat%messlopbdet, &
        wombat%messlopzoo, &
        wombat%mesassife, &
        wombat%reminr, &
        wombat%docremi, &
        wombat%detremi, &
        wombat%bdetremi, &
        wombat%pic2poc, &
        wombat%dissrat, &
        wombat%caldiss, &
        wombat%aoa_loxy, &
        wombat%aoa_lnh4, &
        wombat%aoa_mu, &
        wombat%het_loxy, &
        wombat%het_lno3, &
        wombat%het_ldet, &
        wombat%het_mu, &
        wombat%aox_lnh4, &
        wombat%aox_mu, &
        wombat%nitrfix, &
        wombat%ammox, &
        wombat%anammox, &
        wombat%denitrif, &
        wombat%fdenitrif, &
        wombat%no3_prev, &
        wombat%caco3_prev, &
        wombat%det_sed_remin, &
        wombat%det_sed_denit, &
        wombat%det_btm, &
        wombat%bdet_btm, &
        wombat%fbury, &
        wombat%fdenit, &
        wombat%detfe_sed_remin, &
        wombat%detfe_btm, &
        wombat%bdetfe_btm, &
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
        wombat%sedo2, &
        wombat%seddic, &
        wombat%sedalk, &
        wombat%sedhtotal, &
        wombat%sedco3, &
        wombat%sedomega_cal)

  end subroutine user_deallocate_arrays

end module generic_WOMBATmid
