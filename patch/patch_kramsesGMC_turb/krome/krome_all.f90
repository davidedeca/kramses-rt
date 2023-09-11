
!############### MODULE ##############
module krome_commons

  implicit none

  ! *************************************************************
  !  This file has been generated with:
  !  krome 14.08.dev "Beastie Boyle" on 2016-08-08 10:48:19.
  !
  ! KROME is a nice and friendly chemistry package for a wide range of
  !  astrophysical simulations. Given a chemical network (in CSV format)
  !  it automatically generates all the routines needed to solve the kinetic
  !  of the system, modelled as system of coupled Ordinary Differential
  !  Equations.
  !  It provides different options which make it unique and very flexible.
  !  Any suggestions and comments are welcomed. KROME is an open-source
  !  package, GNU-licensed, and any improvements provided by
  !  the users is well accepted. See disclaimer below and GNU License
  !  in gpl-3.0.txt.
  !
  !  more details in http://kromepackage.org/
  !  also see https://bitbucket.org/krome/krome_stable
  !
  ! Written and developed by Tommaso Grassi
  !  tommasograssi@gmail.com,
  !  Starplan Center, Copenhagen.
  !  Niels Bohr Institute, Copenhagen.
  !
  ! Co-developer Stefano Bovino
  !  sbovino@astro.physik.uni-goettingen.de
  !  Institut fuer Astrophysik, Goettingen.
  !
  ! Others (alphabetically): D. Galli, F.A. Gianturco, T. Haugboelle,
  !  J.Prieto, J.Ramsey, D.R.G. Schleicher, D. Seifried,
  !  E. Simoncini, E. Tognelli
  !
  !
  ! KROME is provided "as it is", without any warranty.
  !  The Authors assume no liability for any damages of any kind
  !  (direct or indirect damages, contractual or non-contractual
  !  damages, pecuniary or non-pecuniary damages), directly or
  !  indirectly derived or arising from the correct or incorrect
  !  usage of KROME, in any possible environment, or arising from
  !  the impossibility to use, fully or partially, the software,
  !  or any bug or malefunction.
  !  Such exclusion of liability expressly includes any damages
  !  including the loss of data of any kind (including personal data)
  ! *************************************************************
  character(len=256),save::krome_datafolder="./data_folder/"

  integer,parameter::idx_H=1
  integer,parameter::idx_E=2
  integer,parameter::idx_Hj=3
  integer,parameter::idx_HE=4
  integer,parameter::idx_HEj=5
  integer,parameter::idx_HEjj=6
  integer,parameter::idx_Hk=7
  integer,parameter::idx_H2=8
  integer,parameter::idx_H2j=9
  integer,parameter::idx_CR=10
  integer,parameter::idx_g=11
  integer,parameter::idx_Tgas=12
  integer,parameter::idx_dummy=13
  integer,parameter::nrea=52
  integer,parameter::nmols=9
  integer,parameter::nspec=13
  integer,parameter::natoms=3
  integer,parameter::ndust=0
  integer,parameter::ndustTypes=0
  integer,parameter::nPhotoBins=10
  integer,parameter::nPhotoRea=8

  !cooling index
  integer,parameter::idx_cool_H2 = 1
  integer,parameter::idx_cool_H2GP = 2
  integer,parameter::idx_cool_atomic = 3
  integer,parameter::idx_cool_CEN = 3
  integer,parameter::idx_cool_HD = 4
  integer,parameter::idx_cool_Z = 5
  integer,parameter::idx_cool_metal = 5
  integer,parameter::idx_cool_enthalpic = 6
  integer,parameter::idx_cool_dH = 6
  integer,parameter::idx_cool_dust = 7
  integer,parameter::idx_cool_compton = 8
  integer,parameter::idx_cool_CIE = 9
  integer,parameter::idx_cool_cont = 10
  integer,parameter::idx_cool_continuum = 10
  integer,parameter::idx_cool_expansion = 11
  integer,parameter::idx_cool_exp = 11
  integer,parameter::idx_cool_ff = 12
  integer,parameter::idx_cool_bss = 12
  integer,parameter::idx_cool_custom = 13
  integer,parameter::idx_cool_CO = 14
  integer,parameter::idx_cool_ZCIE = 15
  integer,parameter::idx_cool_ZCIENOUV = 16
  integer,parameter::idx_cool_ZExtend = 17
  integer,parameter::idx_cool_gh = 18
  integer,parameter::ncools = 18

  !heating index
  integer,parameter::idx_heat_chem = 1
  integer,parameter::idx_heat_compress = 2
  integer,parameter::idx_heat_compr = 2
  integer,parameter::idx_heat_photo = 3
  integer,parameter::idx_heat_enthalpic = 4
  integer,parameter::idx_heat_dH = 4
  integer,parameter::idx_heat_photoAv = 5
  integer,parameter::idx_heat_Av = 5
  integer,parameter::idx_heat_CR = 6
  integer,parameter::idx_heat_dust = 7
  integer,parameter::idx_heat_xray = 8
  integer,parameter::idx_heat_visc = 9
  integer,parameter::idx_heat_viscous = 9
  integer,parameter::idx_heat_custom = 10
  integer,parameter::idx_heat_ZCIE = 11
  integer,parameter::idx_heat_gh = 12
  integer,parameter::nheats = 12

  real*8::arr_k(nrea)

  !commons for rate tables
  !modify ktab_n according to the requested precision
  integer,parameter::ktab_n=int(1e3)
  real*8::ktab(nrea,ktab_n),ktab_logTlow, ktab_logTup, ktab_T(ktab_n)
  real*8::inv_ktab_T(ktab_n-1), inv_ktab_idx

  !thermo toggle (when >0 do cooling/heating)
  integer::krome_thermo_toggle
  !$omp threadprivate(krome_thermo_toggle)

  !commons for implicit RHS
  integer::arr_r1(nrea)
  integer::arr_r2(nrea)
  integer::arr_p1(nrea)
  integer::arr_p2(nrea)
  integer::arr_p3(nrea)

  !commons for reduction
  integer::arr_u(nrea)
  real*8::arr_flux(nrea)

  !commons for frequency bins
  real*8::photoBinJ(nPhotoBins) !intensity per bin, eV/sr/cm2
  real*8::photoBinJ_org(nPhotoBins) !intensity per bin stored, eV/sr/cm2
  real*8::photoBinEleft(nPhotoBins) !left limit of the freq bin, eV
  real*8::photoBinEright(nPhotoBins) !right limit of the freq bin, eV
  real*8::photoBinEmid(nPhotoBins) !middle point of the freq bin, eV
  real*8::photoBinEdelta(nPhotoBins) !size of the freq bin, eV
  real*8::photoBinEidelta(nPhotoBins) !inverse of the size of the freq bin, 1/eV
  real*8::photoBinJTab(nPhotoRea,nPhotoBins) !xsecs table, cm2
  real*8::photoBinRates(nPhotoRea) !photo rates, 1/s
  real*8::photoBinHeats(nPhotoRea) !photo heating, erg/s
  real*8::photoBinEth(nPhotoRea) !energy treshold, eV
  real*8::photoPartners(nPhotoRea) !index of the photoreactants
  real*8::opacityDust(nPhotoBins) !interpolated opacity from tables
  real*8::cooling_GH,heating_GH
  !$omp threadprivate(photoBinJ,photoBinEleft,photoBinEright,photoBinEmid,photoBinEdelta, &
      !$omp    photoBinEidelta,photoBinJTab,photoBinRates,photoBinHeats,photoBinEth,photoPartners)

  !commons for dust optical properties

  !mpi rank of process. If 0, ignored
  integer::krome_mpi_rank=0, krome_omp_thread
  ! open file in the following location
  integer::krome_nfile,krome_nfile2
  !$omp threadprivate(krome_omp_thread)

  !user-defined commons variables from the reaction file
  real*8::user_crate,user_myfluxLW, user_myH2_dissociation, user_cell_size
  !$omp threadprivate(user_crate,user_myfluxLW,user_cell_size)

  !commons for anytab

  !physical commons
  real*8::phys_Tcmb
  real*8::phys_zredshift
  real*8::phys_orthoParaRatio
  real*8::phys_metallicity
  real*8::phys_Tfloor
  !$omp threadprivate(phys_Tcmb)
  !$omp threadprivate(phys_zredshift)
  !$omp threadprivate(phys_orthoParaRatio)
  !$omp threadprivate(phys_metallicity)
  !$omp threadprivate(phys_Tfloor)

  !machine precision
  real*8::krome_epsilon

  !xrayJ21 for tabulated heating and rate
  real*8::J21xray

  !total metallicity relative to solar Z/Z_solar
  real*8::total_Z
  real*8::dust2gas_ratio, dust2gas_ratio_ref

  integer,parameter::CoolZNOUVn=131,CoolZNOUVm=162
  real*8::CoolZNOUV_x(CoolZNOUVn),CoolZNOUV_y(CoolZNOUVm)
  real*8::CoolZNOUV_z(CoolZNOUVn,CoolZNOUVm),CoolZNOUV_xmul,CoolZNOUV_ymul

  !data for metal cooling from table in the presence of UV
  integer,parameter::coolZCIEn1=81
  integer,parameter::coolZCIEn2=81
  integer,parameter::coolZCIEn3=81
  real*8::coolZCIEx1(coolZCIEn1),coolZCIEx2(coolZCIEn2),coolZCIEx3(coolZCIEn3)
  real*8::coolZCIEixd1(coolZCIEn1-1),coolZCIEixd2(coolZCIEn2-1)
  real*8::coolZCIEixd3(coolZCIEn3-1)
  real*8::coolZCIEy(coolZCIEn1,coolZCIEn2,coolZCIEn3)
  real*8::heatZCIEy(coolZCIEn1,coolZCIEn2,coolZCIEn3)
  real*8::coolZCIEx1min,coolZCIEx1max
  real*8::coolZCIEx2min,coolZCIEx2max
  real*8::coolZCIEx3min,coolZCIEx3max
  real*8::coolZCIEdvn1,coolZCIEdvn2,coolZCIEdvn3

  !commons for dust tabs (cool,H2,Tdust)
  integer,parameter::dust_tab_imax=50,dust_tab_jmax=50
  real*8::dust_tab_ngas(dust_tab_imax)
  real*8::dust_tab_Tgas(dust_tab_jmax)
  real*8::dust_tab_cool(dust_tab_imax,dust_tab_jmax)
  real*8::dust_tab_heat(dust_tab_imax,dust_tab_jmax)
  real*8::dust_tab_Tdust(dust_tab_imax,dust_tab_jmax)
  real*8::dust_tab_H2(dust_tab_imax,dust_tab_jmax)
  real*8::dust_mult_ngas,dust_mult_Tgas

  !commons for exp(-a) table
  integer,parameter::exp_table_na=int(1d5)
  real*8,parameter::exp_table_aMax=1d4,exp_table_aMin=0d0
  real*8,parameter::exp_table_multa=(exp_table_na-1) &
      / (exp_table_aMax-exp_table_aMin)
  real*8,parameter::exp_table_da=1d0/exp_table_multa
  real*8::exp_table(exp_table_na)

  !stores the last evaluation of the rates in the fex
  real*8::last_coe(nrea)
  !$omp threadprivate(last_coe)

  !xsecs from file variables
  !xsec for H- -> H + E
  real*8,allocatable::xsec39_val(:)
  real*8::xsec39_Emin
  real*8::xsec39_idE
  integer::xsec39_n

  !xsec for H2 -> H2+ + E
  real*8,allocatable::xsec40_val(:)
  real*8::xsec40_Emin
  real*8::xsec40_idE
  integer::xsec40_n

  !xsec for H2+ -> H+ + H
  real*8,allocatable::xsec41_val(:)
  real*8::xsec41_Emin
  real*8::xsec41_idE
  integer::xsec41_n

  !partition function from file
  integer,parameter::zpart_nCO=641
  integer,parameter::zpart_nH2even=2000
  integer,parameter::zpart_nH2odd=2000
  real*8::zpart_CO(zpart_nCO),minpart_CO,partdT_CO
  real*8::zpart_H2even(zpart_nH2even),minpart_H2even,partdT_H2even
  real*8::zpart_H2odd(zpart_nH2odd),minpart_H2odd,partdT_H2odd

  !Habing flux for the photoelectric heating by dust
  ! and clumping factor for H2 formation
  ! on dust by Jura/Gnedin
  real*8::GHabing,Ghabing_thin,clump_factor

  ! Photo reaction rates relevant for Gnedin-Hollon cooling/heating function
  real*8::PLW,PHI,PHeI,PCVI
  !$omp threadprivate(PLW,PHI,PHeI,PCVI)
  !partition functions common vars

  ! bin averaged cross sections for Gnedin Hoolon (cm2*Hz)
  real*8,dimension(nPhotoBins)::HI_sigma = (/ 0.0,0.0,0.0,0.0,0.0008126713005254592, &
                      0.0015260648187972584,0.00539377494778203,0.0012270823637565093, &
                      0.0018177982203954513,0.0007848359841450194  /)
  real*8,dimension(nPhotoBins)::HeI_sigma = (/ 0.0,0.0,0.0,0.0,0.0,0.0,0.0, &
                      0.008278006431718244,0.01809007673951052,0.014403819587745273 /)
  real*8,dimension(nPhotoBins)::CVI_sigma = (/ 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, &
                      0.008532139279067026 /)

end module krome_commons

!############### MODULE ##############
module krome_constants
  implicit none

  !constants
  real*8,parameter::boltzmann_eV = 8.617332478d-5 !eV / K
  real*8,parameter::boltzmann_J = 1.380648d-23 !J / K
  real*8,parameter::boltzmann_erg = 1.380648d-16 !erg / K
  real*8,parameter::planck_eV = 4.135667516d-15 !eV s
  real*8,parameter::planck_J = 6.62606957d-34 !J s
  real*8,parameter::planck_erg = 6.62606957d-27 !erg s
  real*8,parameter::iplanck_eV = 1d0/planck_eV !1 / eV / s
  real*8,parameter::iplanck_J = 1d0/planck_J !1 / J / s
  real*8,parameter::iplanck_erg = 1d0/planck_erg !1 / erg / s
  real*8,parameter::gravity = 6.674d-8 !cm3 / g / s2
  real*8,parameter::e_mass = 9.10938188d-28 !g
  real*8,parameter::p_mass = 1.67262158d-24 !g
  real*8,parameter::n_mass = 1.674920d-24 !g
  real*8,parameter::ip_mass = 1d0/p_mass !1/g
  real*8,parameter::clight = 2.99792458e10 !cm/s
  real*8,parameter::pi = 3.14159265359d0 !#
  real*8,parameter::eV_to_erg = 1.60217646d-12 !eV -> erg
  real*8,parameter::ry_to_eV = 13.60569d0 !rydberg -> eV
  real*8,parameter::ry_to_erg = 2.179872d-11 !rydberg -> erg
  real*8,parameter::seconds_per_year = 365d0*24d0*3600d0 !yr -> s
  real*8,parameter::km_to_cm = 1d5 !km -> cm
  real*8,parameter::cm_to_Mpc = 1.d0/3.08d24 !cm -> Mpc
  real*8,parameter::kvgas_erg = 8.d0*boltzmann_erg/pi/p_mass !
  real*8,parameter::pre_kvgas_sqrt = sqrt(8.d0*boltzmann_erg/pi) !
  real*8,parameter::pre_planck = 2.d0*planck_erg/clight**2 !erg/cm2*s3
  real*8,parameter::exp_planck = planck_erg / boltzmann_erg !s*K
  real*8,parameter::stefboltz_erg = 5.670373d-5 !erg/s/cm2/K4
  real*8,parameter::N_avogadro = 6.0221d23 !#
  real*8,parameter::Rgas_J = 8.3144621d0 !J/K/mol
  real*8,parameter::Rgas_kJ = 8.3144621d-3 !kJ/K/mol
  real*8,parameter::hubble = 0.704d0 !dimensionless
  real*8,parameter::Omega0 = 1.0d0 !dimensionless
  real*8,parameter::Omegab = 0.0456d0 !dimensionless
  real*8,parameter::Hubble0 = 1.d2*hubble*km_to_cm*cm_to_Mpc !1/s

end module krome_constants

!############### MODULE ##############
module krome_subs
contains

  ! *************************************************************
  !  This file has been generated with:
  !  krome 14.08.dev "Beastie Boyle" on 2016-08-08 10:48:19.
  !
  ! KROME is a nice and friendly chemistry package for a wide range of
  !  astrophysical simulations. Given a chemical network (in CSV format)
  !  it automatically generates all the routines needed to solve the kinetic
  !  of the system, modelled as system of coupled Ordinary Differential
  !  Equations.
  !  It provides different options which make it unique and very flexible.
  !  Any suggestions and comments are welcomed. KROME is an open-source
  !  package, GNU-licensed, and any improvements provided by
  !  the users is well accepted. See disclaimer below and GNU License
  !  in gpl-3.0.txt.
  !
  !  more details in http://kromepackage.org/
  !  also see https://bitbucket.org/krome/krome_stable
  !
  ! Written and developed by Tommaso Grassi
  !  tommasograssi@gmail.com,
  !  Starplan Center, Copenhagen.
  !  Niels Bohr Institute, Copenhagen.
  !
  ! Co-developer Stefano Bovino
  !  sbovino@astro.physik.uni-goettingen.de
  !  Institut fuer Astrophysik, Goettingen.
  !
  ! Others (alphabetically): D. Galli, F.A. Gianturco, T. Haugboelle,
  !  J.Prieto, J.Ramsey, D.R.G. Schleicher, D. Seifried,
  !  E. Simoncini, E. Tognelli
  !
  !
  ! KROME is provided "as it is", without any warranty.
  !  The Authors assume no liability for any damages of any kind
  !  (direct or indirect damages, contractual or non-contractual
  !  damages, pecuniary or non-pecuniary damages), directly or
  !  indirectly derived or arising from the correct or incorrect
  !  usage of KROME, in any possible environment, or arising from
  !  the impossibility to use, fully or partially, the software,
  !  or any bug or malefunction.
  !  Such exclusion of liability expressly includes any damages
  !  including the loss of data of any kind (including personal data)
  ! *************************************************************

  !************************
  !compute reaction rates cm^3(n-1)/s
  function coe(n)
    use krome_commons
    use krome_constants
    use krome_user_commons
    implicit none
    real*8::coe(nrea),k(nrea),Tgas,n(nspec),kmax
    real*8::Te
    real*8::lnTe
    real*8::invT
    real*8::invTe
    real*8::small,nmax
    integer::i
    real*8::krome_fshieldHM  !preproc from coevar
    real*8::krome_fshieldH2  !preproc from coevar
    real*8::log10Tgas !preproc from coevar
    real*8::asav2 !preproc from coevar
    real*8::asav3 !preproc from coevar
    real*8::asav0 !preproc from coevar
    real*8::asav1 !preproc from coevar
    real*8::asav6 !preproc from coevar
    real*8::asav7 !preproc from coevar
    real*8::asav4 !preproc from coevar
    real*8::asav5 !preproc from coevar
    real*8::a11 !preproc from coevar
    real*8::a10 !preproc from coevar
    real*8::a12 !preproc from coevar
    real*8::a1 !preproc from coevar
    real*8::a3 !preproc from coevar
    real*8::a2 !preproc from coevar
    real*8::a5 !preproc from coevar
    real*8::T !preproc from coevar
    real*8::a7 !preproc from coevar
    real*8::a6 !preproc from coevar
    real*8::a9 !preproc from coevar
    real*8::a8 !preproc from coevar
    real*8::krome_fshieldH2j !preproc from coevar
    real*8::bsav1 !preproc from coevar
    real*8::bsav0 !preproc from coevar
    real*8::bsav3 !preproc from coevar
    real*8::bsav2 !preproc from coevar
    real*8::bsav5 !preproc from coevar
    real*8::bsav4 !preproc from coevar
    real*8::bsav7 !preproc from coevar
    real*8::bsav6 !preproc from coevar
    real*8::a4 !preproc from coevar
    !Tgas is in K
    Tgas = max(n(idx_Tgas), phys_Tcmb)
    Tgas = min(Tgas,1d9)

    !maxn initialization can be removed and small can be
    ! replaced with a proper value according to the environment
    nmax = max(maxval(n(1:nmols)),1d0)
    small = 0d0

    Te = Tgas*8.617343d-5 !Tgas in eV (eV)
    lnTe = log(Te) !ln of Te (#)
    invT = 1.d0/Tgas !inverse of T (1/K)
    invTe = 1.d0/Te !inverse of T (1/eV)

    T = Tgas
    a1 = 1.3500e-09
    a2 = 9.8493e-02
    a3 = 3.2852e-01
    a4 = 5.5610e-01
    a5 = 2.7710e-07
    a6 = 2.1826e+00
    a7 = 6.1910e-03
    a8 = 1.0461e+00
    a9 = 8.9712e-11
    a10 = 3.0424e+00
    a11 = 3.2576e-14
    a12 = 3.7741e+00
    log10Tgas = log10(Tgas)
    asav0 = -1.9153214d2
    asav1 = 4.0129114d2
    asav2 = -3.7446991d2
    asav3 = 1.9078410d2
    asav4 = -5.7263467d1
    asav5 = 1.0133210d1
    asav6 = -9.8012853d-1
    asav7 = 4.0023414d-2
    bsav0 = -8.8755774d3
    bsav1 = 1.0081246d4
    bsav2 = -4.8606622d3
    bsav3 = 1.2889659d3
    bsav4 = -2.0319575d2
    bsav5 = 1.9057493d1
    bsav6 = -9.8530668d-1
    bsav7 = 2.1675387d-2
    krome_fshieldH2 = shield_dust(n,Tgas,3.74d0)*krome_fshield(n,Tgas)
    krome_fshieldHM = shield_dust(n,Tgas,0.5d0)
    krome_fshieldH2j = shield_dust(n,Tgas,1.9d0)

    k(:) = small !inizialize coefficients

    !H + E -> H+ + E + E
    k(1) = small + (exp(-32.71396786d0+13.5365560d0&
        *lnTe-5.73932875d0*(lnTe**2)+1.56315498d0&
        *(lnTe**3)-0.28770560d0*(lnTe**4)+3.48255977d-2&
        *(lnTe**5)-2.63197617d-3*(lnTe**6)+1.11954395d-4&
        *(lnTe**7)-2.03914985d-6*(lnTe**8)))

    !H+ + E -> H
    if(Tgas.LE.5.5d3) then
      k(2) = small + (3.92d-13&
          *invTe**0.6353d0)
    end if

    !H+ + E -> H
    if(Tgas.GT.5.5d3) then
      k(3) = small + (exp(-28.61303380689232d0-0.7241125657826851d0&
          *lnTe-0.02026044731984691d0*lnTe**2-0.002380861877349834d0&
          *lnTe**3-0.0003212605213188796d0&
          *lnTe**4-0.00001421502914054107d0&
          *lnTe**5+4.989108920299513d-6*lnTe**6+5.755614137575758d-7&
          *lnTe**7-1.856767039775261d-8*lnTe**8-3.071135243196595d-9&
          *lnTe**9))
    end if

    !HE + E -> HE+ + E + E
    k(4) = small + (exp(-44.09864886d0+23.91596563d0&
        *lnTe-10.7532302d0*(lnTe**2)+3.05803875d0&
        *(lnTe**3)-0.56851189d0*(lnTe**4)+6.79539123d-2&
        *(lnTe**5)-5.00905610d-3*(lnTe**6)+2.06723616d-4&
        *(lnTe**7)-3.64916141d-6*(lnTe**8)))

    !HE+ + E -> HE
    if(Tgas.LE.9.284d3) then
      k(5) = small + (3.92d-13&
          *invTe**0.6353d0)
    end if

    !HE+ + E -> HE
    if(Tgas.GT.9.284d3) then
      k(6) = small + (1.54d-9&
          *(1.d0+0.3d0&
          /exp(8.099328789667d0&
          *invTe))&
          /(exp(40.49664394833662d0*invTe)&
          *Te**1.5d0)+3.92d-13&
          /Te**0.6353d0)
    end if

    !HE+ + E -> HE++ + E + E
    k(7) = small + (exp(-68.71040990212001d0+43.93347632635d0&
        *lnTe-18.48066993568d0*lnTe**2+4.701626486759002d0&
        *lnTe**3-0.7692466334492d0*lnTe**4+0.08113042097303d0&
        *lnTe**5-0.005324020628287001d0*lnTe**6+0.0001975705312221d0&
        *lnTe**7-3.165581065665d-6*lnTe**8))

    !HE++ + E -> HE+
    k(8) = small + (1.891d-10/(sqrt(Tgas&
        /9.37)&
        *(1.+sqrt(Tgas/9.37))**0.2476&
        *(1.+sqrt(Tgas&
        /2.774d6))**1.7524))

    !H + E -> H-
    k(9) = small + (1.4d-18*Tgas**0.928&
        *exp(-Tgas&
        /16200.))

    !H- + H -> H2 + E
    k(10) = small + (a1*(Tgas**a2+a3&
        *Tgas**a4+a5*Tgas**a6)&
        /(1.+a7*Tgas**a8+a9*Tgas**a10+a11&
        *Tgas**a12))

    !H + H+ -> H2+
    if(Tgas.LT.3d1) then
      k(11) = small + (2.10e-20&
          *(Tgas&
          /30.)**(-0.15))
    end if

    !H + H+ -> H2+
    if(Tgas.GE.3d1) then
      k(12) = small + (10**(-18.20-3.194&
          *log10Tgas+1.786*(log10Tgas)**2-0.2072&
          *(log10Tgas)**3))
    end if

    !H2+ + H -> H2 + H+
    k(13) = small + (6.0d-10)

    !H2 + H+ -> H2+ + H
    if(Tgas.GE.0.d0 .and. Tgas.LT.1.d5) then
      k(14) = small + (10**(asav0+asav1&
          *log10Tgas+asav2*(log10Tgas)**2+asav3*(log10Tgas)**3+asav4&
          *(log10Tgas)**4+asav5*(log10Tgas)**5+asav6&
          *(log10Tgas)**6+asav7*(log10Tgas)**7))
    end if

    !H2 + H+ -> H2+ + H
    if(Tgas.GE.1.d5 .and. Tgas.LE.1.d8) then
      k(15) = small + (10**(bsav0+bsav1&
          *log10Tgas+bsav2*(log10Tgas)**2+bsav3*(log10Tgas)**3+bsav4&
          *(log10Tgas)**4+bsav5*(log10Tgas)**5+bsav6&
          *(log10Tgas)**6+bsav7*(log10Tgas)**7))
    end if

    !H2 + E -> H + H-
    k(16) = small + (3.55d1*Tgas**(-2.28)&
        *exp(-46707.&
        /Tgas))

    !H2 + E -> H + H + E
    k(17) = small + (4.38d-10&
        *exp(-102000.&
        /Tgas)*Tgas**(0.35))

    !H2 + H -> H + H + H
    k(18) = small + (6.67e-12*sqrt(Tgas)&
        *exp(-(1+63593&
        /Tgas)))

    !H- + E -> H + E + E
    k(19) = small + (exp(-18.01849334273d0+2.360852208681d0&
        *lnTe-0.2827443061704d0*lnTe**2+0.01623316639567d0&
        *lnTe**3-0.03365012031362999d0*lnTe**4+0.01178329782711d0&
        *lnTe**5-0.001656194699504d0*lnTe**6+0.0001068275202678d0&
        *lnTe**7-2.631285809207d-6*lnTe**8))

    !H- + H -> H + H + E
    if(Tgas.LE.1.16d3) then
      k(20) = small + (2.56d-9&
          *Te**1.78186d0)
    end if

    !H- + H -> H + H + E
    if(Tgas.GT.1.16d3) then
      k(21) = small + (exp(-20.37260896533324d0+1.139449335841631d0&
          *lnTe-0.1421013521554148d0*lnTe**2+0.00846445538663d0&
          *lnTe**3-0.0014327641212992d0*lnTe**4+0.0002012250284791d0&
          *lnTe**5+0.0000866396324309d0*lnTe**6-0.00002585009680264d0&
          *lnTe**7+2.4555011970392d-6*lnTe**8-8.06838246118d-8&
          *lnTe**9))
    end if

    !H- + H+ -> H + H
    if(Tgas.GE.1d1 .and. Tgas.LE.1d7) then
      k(22) = small + ((2.96d-6&
          /sqrt(Tgas)-1.73d-9+2.50d-10&
          *sqrt(Tgas)-7.77d-13))
    end if

    !H- + H+ -> H2+ + E
    k(23) = small + (1.d-8*Tgas**(-0.4d0))

    !H2+ + E -> H + H
    if(Tgas.LE.6.17d2) then
      k(24) = small + (1.d-8)
    end if

    !H2+ + E -> H + H
    if(Tgas.GT.6.17d2) then
      k(25) = small + (1.32d-6&
          *Tgas**(-0.76d0))
    end if

    !H2+ + H- -> H + H2
    k(26) = small + (5.d-7*sqrt(1.d2*invT))

    !H2 + H2 -> H + H + H2
    k(27) = small + (5.996d-30&
        *Tgas**(4.1881)*(1d0+6.761d-6*Tgas)**(-5.6881)*exp(-54657.4&
        *invT))

    !HE+ + H -> HE + H+
    k(28) = small + (1.20d-15&
        *(Tgas&
        /3d2)**0.25d0)

    !HE + H+ -> HE+ + H
    if(Tgas.LE.1d4) then
      k(29) = small + (1.26d-9&
          *Tgas**(-0.75d0)*exp(-1.275d5*invT))
    end if

    !HE + H+ -> HE+ + H
    if(Tgas.GT.1d4) then
      k(30) = small + (4.d-37&
          *Tgas**(4.74d0))
    end if

    !H2 + HE+ -> HE + H + H+
    k(31) = small + (3.7d-14&
        *exp(-35.d0&
        /Tgas))

    !H2 + HE -> H + H + HE
    k(32) = small + (10**(-27.029 + 3.801&
        *log10Tgas-29487d0&
        /Tgas))

    !H2 + HE+ -> H2+ + HE
    k(33) = small + (7.2d-15)

    !H + H -> H + H+ + E
    k(34) = small + (1.2d-17*Tgas**(1.2d0)&
        *exp(-157800*invT))

    !H + HE -> HE + H+ + E
    k(35) = small + (1.75d-17&
        *Tgas**(1.3d0)*exp(-157800*invT))

    !H -> H+ + E
    k(36) = small + (photoBinRates(1))

    !HE -> HE+ + E
    k(37) = small + (photoBinRates(2))

    !HE+ -> HE++ + E
    k(38) = small + (photoBinRates(3))

    !H- -> H + E
    k(39) = small + (photoBinRates(4))

    !H2 -> H2+ + E
    k(40) = small + (photoBinRates(5))

    !H2+ -> H+ + H
    k(41) = small + (photoBinRates(6))

    !H2+ -> H+ + H+ + E
    k(42) = small + (photoBinRates(7))

    !H2 -> H + H
    k(43) = small + (photoBinRates(8))

    !H2 -> H + H
    k(44) = H2_solomonLW_ramses(user_myH2_dissociation, user_myfluxLW)&
        *krome_fshieldH2

    !H+ + E -> H
    k(45) = small + (H_recombination_on_dust(n,Tgas))

    !HE+ + E -> HE
    k(46) = small + (He_recombination_on_dust(n,Tgas))

    !H -> H+ + E
    k(47) = small + (4.6d-1*user_crate)

    !HE -> HE+ + E
    k(48) = small + (5.d-1*user_crate)

    !H2 -> H + H
    k(49) = small + (1d-1*user_crate)

    !H2 -> H+ + H-
    k(50) = small + (3d-4*user_crate)

    !H2 -> H2+ + E
    k(51) = small + (9.3d-1*user_crate)

    !H2 -> H + H+ + E
    k(52) = small + (9.3d-1*user_crate)

    coe(:) = k(:) !set coefficients to return variable

    !!uncomment below to check coefficient values
    !kmax = 1d0
    !if(maxval(k)>kmax.or.minval(k)<0d0) then
    !   print *,"***************"
    !   do i=1,size(k)
    !      if(k(i)<0d0.or.k(i)>kmax) print *,i,k(i)
    !   end do
    !end if

  end function coe

  !****************************
  !dust shielding factor
  function shield_dust(n,Tgas,gam)
    use krome_commons
    implicit none
    real*8::shield_dust,n(:),Tgas,gam,eff_d2g
    real*8::sigma_d,NHtot

    shield_dust =1.d0

#if !defined(RT)
    eff_d2g = dust2gas_ratio
    sigma_d = 2d-21*eff_d2g*gam !Richings et al. 2014
    !sigma_d = 2d-21 !Glover+2007
    !sigma_d = 4d-22 !Richings+ 2014
    !sigma_d = 4d-21 !Gnedin 2009

    NHtot = 0d0
    NHtot  = NHtot + num2col(n(idx_H),n(:))
    NHtot  = NHtot + num2col(n(idx_Hj),n(:))
    NHtot  = NHtot + 2d0 * num2col(n(idx_H2),n(:))

    shield_dust = exp(-sigma_d*NHtot)
#endif

  end function shield_dust

  !*******************
  !apply a shielding to Habing flux
  subroutine calcHabingThick(n,Tgas)
    use krome_commons
    implicit none
    real*8::getHabingThick,n(:),Tgas

    GHabing = GHabing_thin * shield_dust(n(:),Tgas,0.665d0)

  end subroutine calcHabingThick

  !*********************
  !return the ratio between the current flux an Draine's
  function get_ratioFluxDraine()
    implicit none
    real*8::get_ratioFluxDraine

    !7.95d-8 eV/cm2/sr is the integrated Draine flux
    get_ratioFluxDraine = get_integratedFlux()/7.95d-8

  end function get_ratioFluxDraine

  !**********************
  !return the curred integrated flux (eV/cm2/sr)
  ! as I(E)/E*dE
  function get_integratedFlux()
    use krome_commons
    implicit none
    integer::j
    real*8::get_integratedFlux,dE

    get_integratedFlux = 0d0
    do j=1,nPhotoBins
      dE = photoBinEdelta(j)
      get_integratedFlux = get_integratedFlux &
          + photoBinJ(j)*dE/photoBinEmid(j)
    end do

  end function get_integratedFlux

  !*******************
  !The following functions compute the recombination rate
  ! on dust for H+, He+, C+, Si+, and O+. See Weingartner&Draine 2001
  ! dust2gas_ratio, D/D_sol, default is assumed equal to Z/Z_sol
  function H_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::H_recombination_on_dust

    H_recombination_on_dust = 0d0
    GHabing = MAX(Ghabing, 1d-40)
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    H_recombination_on_dust =  1.225d-13*dust2gas_ratio &
        /(1.d0+8.074d-6*psi**(1.378)*(1.d0+5.087d2 &
        *Tgas**(0.01586)*psi**(-0.4723-1.102d-5*log(Tgas))))

  end function H_recombination_on_dust

  !******************
  function He_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::He_recombination_on_dust

    He_recombination_on_dust = 0d0
    GHabing = MAX(Ghabing, 1d-40)

    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    He_recombination_on_dust = 5.572d-14*dust2gas_ratio&
        /(1.d0+3.185d-7*psi**(1.512)*(1.d0+5.115d3&
        *Tgas**(3.903d-7)*psi**(-0.4956-5.494d-7*log(Tgas))))

  end function He_recombination_on_dust

  !*******************
  function C_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::C_recombination_on_dust

    C_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    C_recombination_on_dust = 4.558d-13*dust2gas_ratio&
        /(1.d0+6.089d-3*psi**(1.128)*(1.d0+4.331d2&
        *Tgas**(0.04845)*psi**(-0.8120-1.333d-4*log(Tgas))))

  end function C_recombination_on_dust

  !******************
  function Si_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::Si_recombination_on_dust

    Si_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    Si_recombination_on_dust = 2.166d-14*dust2gas_ratio&
        /(1.d0+5.678d-8*psi**(1.874)*(1.d0+4.375d4&
        *Tgas**(1.635d-6)*psi**(-0.8964-7.538d-5*log(Tgas))))

  end function Si_recombination_on_dust

  !********************
  function O_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,k_H
    real*8::O_recombination_on_dust

    k_H = H_recombination_on_dust(n(:),Tgas)
    O_recombination_on_dust = 0.25d0*k_H

  end function O_recombination_on_dust

  !*********************
  !This function returns the
  ! photorate of H2 occurring in the
  ! Lyman-Werner bands following the approximation
  ! provided by Glover&Jappsen 2007. Rate in 1/s.
  !Approximation valid at low-density, it assumes H2(nu = 0).
  !It also stores the rate as a common, needed for the photoheating
  function H2_solomonLW(myflux)
    use krome_commons
    use krome_constants
    implicit none
    real*8::H2_solomonLW,myflux

    !myflux is the radiation background at E = 12.87 eV
    !should be converted to erg

    H2_solomonLW = 1.38d9*myflux*eV_to_erg

  end function H2_solomonLW

  function H2_solomonLW_ramses(myH2_dissociation, myflux)
    use krome_commons
    use krome_constants
    implicit none
    real*8::H2_solomonLW_ramses,myH2_dissociation,myflux
    !real*8::JLW_1draine

    !JLW_1draine = 7.55476d-20   ! previous factor, don't know what it is
    !JLW_1draine = 5.6054d-20  ! J_LW with a draine spectrum of 1 Draine
    !JLW_1draine = 1.187e-19  ! J_LW with a flat FUV spectrum of 1 Draine
    H2_solomonLW_ramses = eV_to_erg
    H2_solomonLW_ramses = H2_solomonLW_ramses *myH2_dissociation
    H2_solomonLW_ramses = H2_solomonLW_ramses *myflux 

  end function H2_solomonLW_ramses

  !**********************
  !planck function in eV/s/cm2/Hz/sr
  ! x is the energy in eV, Tbb the black body
  ! temperature in K
  function planckBB(x,Tbb)
    use krome_constants
    implicit none
    real*8::Tbb,x,xexp,planckBB

    !exponent
    xexp = x/boltzmann_eV/Tbb

    !default value
    planckBB = 0d0

    !limit exp overflow
    if(xexp<3d2.and.x>1d-10) then
      planckBB = 2d0*x**3/planck_eV**2/clight**2 &
          / (exp(xexp)-1d0)
    end if

  end function planckBB

  !********************
  !planck function dTdust differential
  ! in eV/s/cm2/Hz/sr/K, where
  ! x is the energy in eV, Tbb the black body
  ! temperature in K
  function planckBB_dT(x,Tbb)
    use krome_constants
    real*8::a,b,x,Tbb,xexp,planckBB_dT

    b = 1d0/boltzmann_eV
    xexp = b*x/Tbb

    planckBB_dT = 0d0

    if(xexp<3d2) then
      a = 2d0/planck_eV**2/clight**2
      planckBB_dT = a*b*x**4/Tbb/Tbb * exp(xexp)/(exp(xexp)-1d0)**2
    end if

  end function planckBB_dT

  !****************************
  !tanh smoothing function that
  ! increses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_increase(xarg,xpos,slope)
    implicit none
    real*8::smooth_increase,xarg,xpos,slope

    smooth_increase = .5d0 * (tanh(slope * (xarg - xpos)) &
        + 1d0)

  end function smooth_increase

  !****************************
  !tanh smoothing function that
  ! decreses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_decrease(xarg,xpos,slope)
    implicit none
    real*8::smooth_decrease,xarg,xpos,slope

    smooth_decrease = .5d0 * (tanh(-slope * (xarg - xpos)) &
        + 1d0)

  end function smooth_decrease

  !*********************
  !sign: return 1d0 if x>=0d0,
  ! else return -1d0
  function get_sgn(x)
    implicit none
    real*8::x,get_sgn

    get_sgn = 1d0
    if(x==0d0) return
    get_sgn = x/abs(x)

  end function get_sgn

  !***********************
  !shielding function selected with -shield option
  function krome_fshield(n,Tgas)
    implicit none
    real*8::krome_fshield,n(:),Tgas

    krome_fshield = 1d0 !default shielding value

#if !defined(RT)
    !compute shielding from Richings+ 2014
    krome_fshield =  calc_H2shieldR14(n(:), Tgas)
#endif

  end function krome_fshield

  !*********************
  function conserve(n,ni)
    use krome_commons
    implicit none
    real*8::conserve(nspec),n(nspec),ni(nspec),no(nspec)
    real*8::ntot,nitot,factor

    no(:) = n(:)

    !********** H **********
    ntot = n(idx_H) &
        + n(idx_Hj) &
        + n(idx_Hk) &
        + 2d0*n(idx_H2) &
        + 2d0*n(idx_H2j)
    nitot = ni(idx_H) &
        + ni(idx_Hj) &
        + ni(idx_Hk) &
        + 2d0*ni(idx_H2) &
        + 2d0*ni(idx_H2j)
    factor = nitot/ntot
    no(idx_H) = n(idx_H) * factor
    no(idx_Hj) = n(idx_Hj) * factor
    no(idx_Hk) = n(idx_Hk) * factor
    no(idx_H2) = n(idx_H2) * factor
    no(idx_H2j) = n(idx_H2j) * factor

    !********** He **********
    ntot = n(idx_HE) &
        + n(idx_HEj) &
        + n(idx_HEjj)
    nitot = ni(idx_HE) &
        + ni(idx_HEj) &
        + ni(idx_HEjj)
    factor = nitot/ntot
    no(idx_HE) = n(idx_HE) * factor
    no(idx_HEj) = n(idx_HEj) * factor
    no(idx_HEjj) = n(idx_HEjj) * factor

    !********** E **********
    no(idx_E) = max( &
        +n(idx_Hj) &
        +n(idx_HEj) &
        +2d0*n(idx_HEjj) &
        -n(idx_Hk) &
        +n(idx_H2j), 1d-40)

    conserve(:) = 0d0
    conserve(:) = no(:)

  end function conserve

  !*************************
  !this subroutine changes the x(:) mass fractions of the species
  ! to force conservation according to the reference ref(:)
  subroutine conserveLin_x(x,ref)
    use krome_commons
    implicit none
    real*8::x(nmols),ref(natoms)
    real*8::A(natoms,natoms),B(natoms),m(nspec)

    m(:) = get_mass()
    A(:,:) = 0d0
    B(:) = ref(:)

    !charge conservation
    x(idx_E) = m(idx_E)*(+ x(idx_Hj) / m(idx_Hj) &
        + x(idx_HEj) / m(idx_HEj) &
        + 2d0*x(idx_HEjj) / m(idx_HEjj) &
        - x(idx_Hk) / m(idx_Hk) &
        + x(idx_H2j) / m(idx_H2j))
    !check if charge conservation goes wrong
    if(x(idx_E)<0d0) then
      print *,"ERROR in conserveLin, electrons < 0"
      stop
    end if

  end subroutine conserveLin_x

  !***************************
  !compute the total reference mass atom type by atom type
  function conserveLinGetRef_x(x)
    use krome_commons
    implicit none
    real*8::conserveLinGetRef_x(natoms),x(nmols)
    real*8::m(nspec)

    m(:) = get_mass()
    conserveLinGetRef_x(:) = 0d0

  end function conserveLinGetRef_x

  !***************************
  !Ref: Sasaki & Takahara (1993)
  !This function evaluate the recombination rate
  ! for H+ + e --> H + gamma and the same
  ! for D+ + e --> D + gamma
  function elec_recomb_ST93(nabund,nelec,ntot,nucleiH,Trad)
    use krome_commons
    use krome_constants
    implicit none
    real*8::nabund,nelec,Trad
    real*8::nucleiH,elec_recomb_ST93
    real*8::al,ak,rc2,r2c
    real*8::a0,b0,c0,d0,e0
    real*8::a1,b1,c1,d1,e1,f1,g1,h1
    real*8::ntot,ratio

    al = 8.227d0
    ak = 22.06d0 / (hubble  *(1d0 + phys_zredshift) &
        * sqrt(1d0 + Omega0 * phys_zredshift))
    !Rc2 evaluation
    rc2 = 8.76d-11 * (1d0 + phys_zredshift)**(-0.58)
    !R2c evaluation
    r2c = (1.80d10 * Trad)**(1.5) &
        * exp(-3.9472d4 / Trad) * rc2

    !coefficients
    a0 = nucleiH * rc2
    b0 = ak * al * nucleiH
    c0 = ak * rc2 * nucleiH * nucleiH
    d0 = r2c * exp(-1.18416d5/Trad)
    e0 = ak * r2c * nucleiH

    !polynomial terms
    a1 = -d0 * (1d0 + b0)
    b1 = d0 * (1d0 + 2d0 * b0)
    c1 = a0 + b0 * (a0 - d0)
    d1 = -a0 * b0
    e1 = a0 * c0
    f1 = 1d0 + b0 + e0
    g1 = -(b0 + e0)
    h1 = c0

    ratio = nabund / ntot

    elec_recomb_ST93 = ntot*(a1 + b1*ratio + c1*ratio**2 + d1*ratio**3 &
        + e1*ratio**4) / (f1 + g1*ratio + h1*ratio**2)

    elec_recomb_ST93 = elec_recomb_ST93 / (nabund * nelec)

  end function elec_recomb_ST93

  !***************************
  !number density to column density conversion
  function num2col(ncalc,n)
    use krome_commons
    implicit none
    real*8::num2col,ncalc,n(:),Tgas
    Tgas = n(idx_Tgas)

!    num2col = 0.5d0 * ncalc * get_jeans_length(n(:),Tgas)
    num2col = 0.5d0 * ncalc * user_cell_size

  end function num2col

  !***********************
  !column density to number density conversion
  function col2num(ncalc,n)
    use krome_commons
    implicit none
    real*8::col2num,ncalc,n(:),Tgas
    Tgas = n(idx_Tgas)

    !col2num = 2d0 * ncalc / get_jeans_length(n(:),Tgas)
    col2num = 2d0 * ncalc / user_cell_size

  end function col2num

  !**************************
  !shielding function for H2O+ and H3O+
  ! following Glover+2010 MNRAS sect 2.2 eqn.4
  function fHnOj(Av)
    implicit none
    real*8::fHnOj,Av
    if(Av.le.15d0) then
      fHnOj = exp(-2.55*Av+0.0165*Av**2)
    else
      fHnOj = exp(-2.8*Av)
    end if
  end function fHnOj

  !******************************
  !self-shielding for H2
  ! following Glover+2010 MNRAS sect2.2 eqn.6
  ! N: column density (cm-2)
  ! b: doppler broadening (cm/s)
  function fselfH2(N, b)
    implicit none
    real*8::fselfH2,N,b,x,b5

    x = N*2d-15 !normalized column density (#)
    b5 = b*1d-5 !normalized doppler broadening (#)

    fselfH2 = 0.965d0/(1+x/b5)**2 + &
        0.035d0/sqrt(1d0+x) * &
        exp(max(-8.5d-4*sqrt(1+x),-250.))

  end function fselfH2

  !********************
  subroutine load_parts()
    use krome_commons
    implicit none

  end subroutine load_parts

  !*************************
  subroutine load_part(fname,array_part,min_part,dT_part)
    character(len=*)::fname
    integer::ios,icount,i,cv
    real*8,allocatable::array_part(:),emed(:)
    real*8::min_part,dT_part,Told,array_tmp(int(1e5)),rout(2)

    open(krome_nfile,file=trim(fname),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: partition function not found"
      print *," in file "//fname
      stop
    end if

    print *,"loading partition function from "//fname
    icount = 0
    min_part = 1d99
    Told = 0d0
    do
      read(krome_nfile,*,iostat=ios) rout(:)
      if(ios<0) exit
      if(ios.ne.0) cycle
      icount = icount + 1
      min_part = min(min_part,rout(1))
      array_tmp(icount) = rout(2)
      dT_part = rout(1) - Told
      Told = rout(1)
    end do
    close(krome_nfile)

    allocate(array_part(icount),emed(icount))
    array_part(:) = array_tmp(1:icount)

  end subroutine load_part

  !**************************
  !compute 1/(gamma-1) at Tgasin using the partition function
  ! provided in the array_part with a temperature step dT_part
  ! and a minimum Tgas value min_part
  function gamma_pop(array_part,dT_part,min_part,Tgasin)
    implicit none
    real*8::array_part(:),dT_part
    real*8::min_part,Tgas,gamma_pop,Tgas2,Tgasin
    real*8::logz,logz1,logz2,emed1,emed2,Cv,inTgas,T2,T1,Cv1,Cv2
    integer::idx

    !temperature above minimum data point
    inTgas = max(Tgasin,min_part)

    !data index
    idx = (inTgas-min_part)/dT_part+1
    !corresponding Tgas
    Tgas = (idx-1)*dT_part+min_part
    !store Tgas
    T1 = Tgas

    !ln of partition functions (3 points forward)
    logz = log(array_part(idx))
    logz1 = log(array_part(idx+1))
    logz2 = log(array_part(idx+2))

    !derivative for mean energy (2 points forward)
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !derivative for 1/(gamma-1)
    Cv1 = (emed2-emed1)/dT_part

    !next point temperature
    Tgas = (idx)*dT_part+min_part
    !store Tgas
    T2 = Tgas
    !ln of partition functions
    logz = logz1
    logz1 = logz2
    logz2 = log(array_part(idx+3))

    !derivative for mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !derivative for 1/(gamma-1)
    Cv2 = (emed2-emed1)/dT_part

    !interpolation for 1/(gamma-1)
    Cv = (Cv2-Cv1)*(inTgas-T1)/(T2-T1)+Cv1

    !returns result
    gamma_pop = Cv

  end function gamma_pop

  !*****************************
  !compute 1/(gamma-1) at Tgasin using the partition function
  ! provided in the array_part with a temperature step dT_part
  ! and a minimum Tgas value min_part, for H2 with a ortho/para
  ! ratio of opratio. Needs even and odd partition functions.
  function gamma_pop_H2(array_part_even,array_part_odd,dT_part,&
        min_part,Tgasin,opratio)
    implicit none
    real*8::array_part_even(:),array_part_odd(:),dT_part,zcut(4)
    real*8::min_part,Tgas,opratio,gamma_pop_H2,Tgas2,a,b,Tgasin
    real*8::logz,logz1,logz2,emed1,emed2,Cv,inTgas,T2,T1,Cv1,Cv2
    integer::idx

    !Tgas above the data limit
    inTgas = max(Tgasin,min_part)

    !exponents for ortho/para ratio
    a = opratio/(opratio+1d0) !exponent zo
    b = 1d0-a !exponent zp

    !index in the data for the given Tgas
    idx = (inTgas-min_part)/dT_part+1
    !get the corresponding Tgas
    Tgas = (idx-1)*dT_part+min_part
    !store Tgas
    T1 = Tgas

    !needed for ortho partition function (see Boley+2007)
    zcut(1) = exp(2d0*85.4/Tgas)
    zcut(2) = exp(2d0*85.4/(Tgas+dT_part))
    zcut(3) = exp(2d0*85.4/(Tgas+2d0*dT_part))
    zcut(4) = exp(2d0*85.4/(Tgas+3d0*dT_part))

    !ln of the composite partition function
    logz = log(array_part_even(idx)**b*(3d0*array_part_odd(idx)*zcut(1))**a)
    logz1 = log(array_part_even(idx+1)**b*(3d0*array_part_odd(idx+1)*zcut(2))**a)
    logz2 = log(array_part_even(idx+2)**b*(3d0*array_part_odd(idx+2)*zcut(3))**a)
    !derivative for mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !get 1/(gamma-1) for the left point
    Cv1 = (emed2-emed1)/dT_part

    !Tgas of the right point
    Tgas = (idx)*dT_part+min_part
    !store Tgas
    T2 = Tgas
    !ln of the composite function
    logz = logz1
    logz1 = logz2
    logz2 = log(array_part_even(idx+3)**b*(3d0*array_part_odd(idx+3)*zcut(4))**a)
    !derivative for the mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !get 1/(gamma-1) for the right point
    Cv2 = (emed2-emed1)/dT_part

    !interpolation of 1/(gamma-1)
    Cv = (Cv2-Cv1)*(inTgas-T1)/(T2-T1)+Cv1

    !returns the result
    gamma_pop_H2 = Cv
  end function gamma_pop_H2

  !**************************
  !function to get the partition function
  ! of H2 at Tgas with a orto-para ratio
  ! equal to opratio
  function zfop(Tgas,opratio)
    implicit none
    real*8::Tgas,zfop,brot,ibTgas
    real*8::a,b,zo,zp,opratio
    integer::j,jmax,j1
    brot = 85.4d0 !H2 rotational constant in K
    zo = 0d0 !sum for ortho partition function
    zp = 0d0 !sum for para partition function
    jmax = 10 !number of terms in sum

    ibTgas = brot/Tgas !pre-calc

    !loop over levels
    do j=0,jmax,2 !step 2
      j1 = j + 1
      zp = zp + (2d0*j+1d0) * exp(-j*(j+1d0)*ibTgas)
      zo = zo + 3d0 * (2d0*j1+1d0) * exp(-j1*(j1+1d0)*ibTgas)
    end do

    a = opratio/(opratio+1d0) !exponent zo
    b = 1d0-a !exponent zp

    zfop = (zp**b * zo**a*exp(-2d0*ibTgas)) !final partition f

  end function zfop

  !*********************
  !get the partition function at Tgas
  ! of a diatom with rotational constant
  ! brot in K
  function zf(Tgas,brot)
    real*8::Tgas,zf,brot,z,ibTgas
    integer::j,jmax
    jmax = 10 !number of levels

    ibTgas = brot/Tgas !store
    z = 0d0
    !loop on levels
    do j=0,jmax
      z = z + (2d0*j+1d0)*exp(-j*(j+1d0)*ibTgas)
    end do

    zf = z

  end function zf

  !***********************
  !get the degrees of freedom at Tgas for
  ! the rotational component of H2 with
  ! an ortho-para ratio of opratio
  function gamma_rotop(Tgas_in,opratio)
    implicit none
    real*8::gamma_rotop,Tgas,dT,Tgas_in
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,opratio

    Tgas = max(Tgas_in,1d1)

    dT = Tgas*1d-5 !dT for derivative
    idT =  1d0/dT !stored for numeric derivative
    logp1 = log(zfop(Tgas+dT,opratio)) !store since used twice

    !derivative dlog(T)/dT = f(T)
    dlog1 = (logp1-log(zfop(Tgas,opratio)))*idT
    prot1 = dlog1*Tgas**2

    !derivative dlog(T+dT)/dT = f(T+dT)
    dlog2 = (log(zfop(Tgas+dT+dT,opratio))-logp1)*idT
    prot2 = dlog2*(Tgas+dT)**2

    !derivative df(T)/dT
    gamma_rotop = (prot2-prot1)*idT

  end function gamma_rotop

  !***********************
  !get the degrees of freedom at Tgas for
  ! the rotational component of a diatom
  ! with rotational constant brot in K
  function gamma_rot(Tgas_in,brot)
    implicit none
    real*8::gamma_rot,Tgas,dT,Tgas_in
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,brot

    Tgas = max(Tgas_in,1d1)

    dT = Tgas*1d-5 !dT for derivative
    idT =  1d0/dT !stored for numeric derivative
    logp1 = log(zf(Tgas+dT,brot)) !store since used twice

    !derivative dlog(T)/dT = f(T)
    dlog1 = (logp1-log(zf(Tgas,brot)))*idT
    prot1 = dlog1*Tgas**2

    !derivative dlog(T+dT)/dT = f(T+dT)
    dlog2 = (log(zf(Tgas+dT+dT,brot))-logp1)*idT
    prot2 = dlog2*(Tgas+dT)**2

    !derivative df(T)/dT
    gamma_rot = (prot2-prot1)*idT

  end function gamma_rot

  !*********************
  !get gamma
  function gamma_index(n)
    use krome_commons
    implicit none
    real*8::n(:),gamma_index,krome_gamma

    krome_gamma = 1.66667

    gamma_index = krome_gamma
  end function gamma_index

  !*****************************
  !get the mean molecular weight in grams
  function get_mu(n)
    use krome_commons
    use krome_constants
    implicit none
    real*8::n(:),get_mu,m(nspec)
    m(:) = get_mass()

    !ip_mass is 1/proton_mass_in_g
    get_mu = sum(n(1:nmols)*m(1:nmols)) &
        / sum(n(1:nmols)) * ip_mass

  end function get_mu

  !***************************
  !get mean molecular weight in grams
  function get_mu_rho(n,rhogas)
    use krome_commons
    use krome_constants
    implicit none
    real*8::get_mu_rho,rhogas,n(:)

    !ip_mass is 1/proton_mass_in_g
    get_mu_rho = rhogas / sum(n(1:nmols)) * ip_mass

  end function get_mu_rho

  !************************
  !get species masses (g)
  function get_mass()
    use krome_commons
    implicit none
    real*8::get_mass(nspec)

    get_mass(1) = 1.67353251819d-24	!H
    get_mass(2) = 9.10938188d-28	!E
    get_mass(3) = 1.67262158d-24	!H+
    get_mass(4) = 6.69206503638d-24	!HE
    get_mass(5) = 6.69115409819d-24	!HE+
    get_mass(6) = 6.69024316d-24	!HE++
    get_mass(7) = 1.67444345638d-24	!H-
    get_mass(8) = 3.34706503638d-24	!H2
    get_mass(9) = 3.34615409819d-24	!H2+
    get_mass(10) = 0.d0	!CR
    get_mass(11) = 0.d0	!g
    get_mass(12) = 0.d0	!Tgas
    get_mass(13) = 0.d0	!dummy

  end function get_mass

  !************************
  !get sqrt of the inverse of the masses (1/sqrt(g))
  function get_imass_sqrt()
    use krome_commons
    implicit none
    real*8::get_imass_sqrt(nspec)

    get_imass_sqrt(1) = 7.73006102111d+11	!H
    get_imass_sqrt(2) = 3.31326021505d+13	!E
    get_imass_sqrt(3) = 7.732165696d+11	!H+
    get_imass_sqrt(4) = 3.86562679981d+11	!HE
    get_imass_sqrt(5) = 3.86588992536d+11	!HE+
    get_imass_sqrt(6) = 3.86615310465d+11	!HE++
    get_imass_sqrt(7) = 7.72795806394d+11	!H-
    get_imass_sqrt(8) = 5.46597856701d+11	!H2
    get_imass_sqrt(9) = 5.46672253003d+11	!H2+
    get_imass_sqrt(10) = 0.d0	!CR
    get_imass_sqrt(11) = 0.d0	!g
    get_imass_sqrt(12) = 0.d0	!Tgas
    get_imass_sqrt(13) = 0.d0	!dummy

  end function get_imass_sqrt

  !************************
  !get inverse of the species masses (1/g)
  function get_imass()
    use krome_commons
    implicit none
    real*8::get_imass(nspec)

    get_imass(1) = 5.97538433901d+23	!H
    get_imass(2) = 1.09776932527d+27	!E
    get_imass(3) = 5.97863863505d+23	!H+
    get_imass(4) = 1.49430705554d+23	!HE
    get_imass(5) = 1.4945104915d+23	!HE+
    get_imass(6) = 1.49471398286d+23	!HE++
    get_imass(7) = 5.9721335838d+23	!H-
    get_imass(8) = 2.9876921695d+23	!H2
    get_imass(9) = 2.98850552203d+23	!H2+
    get_imass(10) = 0.d0	!CR
    get_imass(11) = 0.d0	!g
    get_imass(12) = 0.d0	!Tgas
    get_imass(13) = 0.d0	!dummy

  end function get_imass

  !************************
  !get species names
  function get_names()
    use krome_commons
    implicit none
    character*16::get_names(nspec)

    get_names(1) = "H"
    get_names(2) = "E"
    get_names(3) = "H+"
    get_names(4) = "HE"
    get_names(5) = "HE+"
    get_names(6) = "HE++"
    get_names(7) = "H-"
    get_names(8) = "H2"
    get_names(9) = "H2+"
    get_names(10) = "CR"
    get_names(11) = "g"
    get_names(12) = "Tgas"
    get_names(13) = "dummy"

  end function get_names

  !******************************
  !get the total number of H nuclei
  function get_Hnuclei(n)
    use krome_commons
    real*8::n(:),get_Hnuclei,nH

    nH = n(idx_H) + &
        n(idx_Hj) + &
        n(idx_Hk) + &
        n(idx_H2)*2d0 + &
        n(idx_H2j)*2d0
    get_Hnuclei = nH

  end function get_Hnuclei

  !***************************
  function get_zatoms()
    use krome_commons
    implicit none
    integer::get_zatoms(nspec)

    get_zatoms(1) = 1	!H
    get_zatoms(2) = 0	!E
    get_zatoms(3) = 1	!H+
    get_zatoms(4) = 2	!HE
    get_zatoms(5) = 2	!HE+
    get_zatoms(6) = 2	!HE++
    get_zatoms(7) = 1	!H-
    get_zatoms(8) = 2	!H2
    get_zatoms(9) = 2	!H2+
    get_zatoms(10) = 0	!CR
    get_zatoms(11) = 0	!g
    get_zatoms(12) = 0	!Tgas
    get_zatoms(13) = 0	!dummy

  end function get_zatoms

  !*****************************
  function fevap(T)
    implicit none
    real*8::T,fevap

    fevap = 1d0

    !if (T > 1.9d3) then
    !    fevap = 0d0
    !elseif (T < 7d2) then
    !    fevap = 1d0
    !else
    !    fevap = 0d0
    !    fevap = fevap + 0.425 * TANH(0.01973359 * (925.  - T))
    !    fevap = fevap + 0.072 * TANH(0.06906755 * (1250. - T))
    !    fevap = fevap + 0.003 * TANH(0.06906755 * (1650. - T))
    !    fevap = fevap + 0.5
    !endif

  end function fevap

  !******************************
  function get_qeff()
    use krome_commons
    implicit none
    real*8::get_qeff(nrea)

    get_qeff(:) = 0e0

  end function get_qeff

  !********************************
  function get_jeans_length(n,Tgas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::m(nspec),get_jeans_length
    m(:) = get_mass()
    rhogas = max(sum(n(1:nmols)*m(1:nmols)),1d-40)
    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length = sqrt(pi*boltzmann_erg*Tgas/rhogas&
        /p_mass/gravity/mu)

  end function get_jeans_length

  !********************************
  function get_jeans_length_rho(n,Tgas,rhogas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::get_jeans_length_rho

    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length_rho = sqrt(pi*boltzmann_erg*Tgas/rhogas&
        /p_mass/gravity/mu)

  end function get_jeans_length_rho

  !Temperature-dependent self-shielding as reported in Richings+2014.
  function calc_H2shieldR14(n,Tgas)
    use krome_commons
    use krome_constants
    implicit none
    real*8::n(nspec),Tgas,calc_H2shieldR14,N_H2,nH2
    real*8::xN_H2,b5,H_mass,bturb,btherm2
    real*8::alpha,omegaH2,Ncrit

    !check on H2 abundances to avoid weird numerical artifacts
    nH2 = max(1d-40, n(idx_H2))
    N_H2=  2d0 * num2col(nH2,n(:))
    !    N_H2 = nH2*get_jeans_length(n(:) ,Tgas)*0.5d0  !column density (cm-2)
    H_mass = p_mass+e_mass !H mass in g
    bturb = 7.1d0*km_to_cm !turbulent Doppler broadening parameter in cm/s
    btherm2 = boltzmann_erg*Tgas/H_mass !thermal Doppler broadening parameter cm/s

    !doppler broadening parameter b divided by 1d5 cm/s (#)
    b5 = ((btherm2 + bturb**2d0)**0.5)*1.d-5
    omegaH2 = 0.013d0*(1d0+(Tgas/2.7d3)**1.3)**(1.0/1.3)*exp(-(Tgas/3.9d3)**14.6)

    if(Tgas<3d3)then
      alpha = 1.4
      Ncrit = 1.3d0*(1d0+(Tgas/6d2)**0.8)
    elseif(Tgas<4d3)then
      alpha = (Tgas/4.5d3)**(-0.8)
      Ncrit = (Tgas/4.76d3)**(-3.8)
    else
      alpha = 1.1
      Ncrit = 2.d0
    endif

    xN_H2 = N_H2*1d-14/Ncrit !normalized column density (#)

    calc_H2shieldR14 = (1d0-omegaH2)/(1d0+xN_H2/b5)**alpha*exp(-5d-7*(1d0+xN_H2)) &
        +(omegaH2/sqrt(1d0+xN_H2))*exp(-8.5d-4*sqrt(1d0+xN_H2))

  end function calc_H2shieldR14

  function calc_H2shieldR14_for_ramses(nH2_in,Tgas,length)
    use krome_commons
    use krome_constants
    implicit none
    real*8::n(nspec),Tgas,calc_H2shieldR14_for_ramses,N_H2,nH2, nH2_in,length
    real*8::xN_H2,b5,H_mass,bturb,btherm2
    real*8::alpha,omegaH2,Ncrit

    !check on H2 abundances to avoid weird numerical artifacts
    nH2 = max(1d-40, nH2_in)
    N_H2=  nH2 * length
    !    N_H2 = nH2*get_jeans_length(n(:) ,Tgas)*0.5d0  !column density (cm-2)
    H_mass = p_mass+e_mass !H mass in g
    bturb = 10.d0*km_to_cm !turbulent Doppler broadening parameter in cm/s
    btherm2 = boltzmann_erg*Tgas/H_mass !thermal Doppler broadening parameter cm/s

    !doppler broadening parameter b divided by 1d5 cm/s (#)
    b5 = ((btherm2 + bturb**2d0)**0.5)*1.d-5


!    if (Tgas <10.)then
!      print*, Tgas
!      print*, Tgas/2.7d3
!      print*, (Tgas/2.7d3)**1.3
!      print*, (1d0+(Tgas/2.7d3)**1.3)**(1.0/1.3)
!      print*, (Tgas/3.9d3)**14.6
!      print*, exp(-(Tgas/3.9d3)**14.6)
!      print*, '--------'
!    endif 

    omegaH2 = 0.013d0*(1d0+(Tgas/2.7d3)**1.3)**(1.0/1.3)*exp(-(Tgas/3.9d3)**14.6)

    if(Tgas<3d3)then
      alpha = 1.4
      Ncrit = 1.3d0*(1d0+(Tgas/6d2)**0.8)
    elseif(Tgas<4d3)then
      alpha = (Tgas/4.5d3)**(-0.8)
      Ncrit = (Tgas/4.76d3)**(-3.8)
    else
      alpha = 1.1
      Ncrit = 2.d0
    endif

    xN_H2 = N_H2*1d-14/Ncrit !normalized column density (#)

    calc_H2shieldR14_for_ramses = (1d0-omegaH2)/(1d0+xN_H2/b5)**alpha*exp(-5d-7*(1d0+xN_H2)) &
        +(omegaH2/sqrt(1d0+xN_H2))*exp(-8.5d-4*sqrt(1d0+xN_H2))

  end function calc_H2shieldR14_for_ramses

  !**********************
  function troe_falloff(k0,kinf,Fc,m)
    implicit none
    real*8::troe_falloff,k0,kinf,Fc,m,rm,xexp
    rm = k0*m/kinf
    xexp = 1d0/(1d0+log10(rm)**2)
    troe_falloff = k0*m/(1d0+rm)*Fc**xexp
  end function troe_falloff

  !*************************
  function k3body(k0,kinf,Fc,nM)
    implicit none
    real*8::k3body,k0,kinf,Fc,nM
    real*8::c,n,d,Pr,xexp,F

    c = -0.4d0-0.67d0*log10(Fc)
    n = 0.75d0-1.27d0*log10(Fc)
    d = 0.14d0
    Pr = k0*nM/kinf
    xexp = (log10(Pr)+c)/(n-d*(log10(Pr)+c))
    F = 1d1**(log10(Fc)/(1d0+xexp**2))
    k3body = kinf*(Pr/(1d0+Pr)) * F

  end function k3body

  !***********************
  !see http://kida.obs.u-bordeaux1.fr/help
  function KIDA3body(ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,&
        kcFc,kdFc,npart,Tgas,pmin,pmax)
    implicit none
    real*8::ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,kcFc,kdFc
    real*8::KIDA3body,kinf,p,f,npart,Tgas,fc,fexp,invT
    real*8::k0,cc,dd,nn,pmin,pmax

    KIDA3body = 0d0

    invT = 1d0/Tgas
    k0 = ka0*(Tgas/3d2)**kb0*exp(-kc0*invT)
    kinf = kainf*(Tgas/3d2)**kbinf*exp(-kcinf*invT)

    p = k0*npart/kinf
    if(p<pmin) return
    if(p>pmax) return

    fc = (1d0-kaFc)*exp(-Tgas/kbFc) + kaFc*exp(-Tgas/kbFc) &
        + exp(-kdFc*invT)

    cc = -0.4d0 - 0.67d0 *log10(fc)
    dd = 0.14d0
    nn = 0.75d0 - 1.27d0*log10(fc)
    fexp = 1d0 + ((log10(p)+cc)/(nn-dd*(log10(p)+cc)))**2

    f = fc**(1d0/fexp)

    KIDA3body = kinf*(p/(1d0+p))*f

  end function KIDA3body

  !******************************
  !collisional ionization rate from Verner+96
  ! unit: cm3/s
  function colion_v96(Tgas,dE,P,A,X,K)
    implicit none
    real*8::colion_v96,Tgas,dE,A,X,K,U,Te,P

    Te = Tgas * 8.621738d-5 !K to eV
    U = dE / Te
    colion_v96 = A * (1d0 + P*sqrt(U)) * U**K * exp(-U) / (X+U)

  end function colion_v96

  !****************************
  !radiative recombination rates from
  ! Verner routine, standard fit, cm3/s
  function recV96(Tgas,a,b)
    implicit none
    real*8::recV96,Tgas,a,b

    recV96 = a*(1d4/Tgas)**b

  end function recV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, new fit, cm3/s
  function recNewV96(Tgas,r1,r2,r3,r4)
    implicit none
    real*8::recNewV96,Tgas,r1,r2,r3,r4,tt

    tt = sqrt(Tgas/r3)
    recNewV96 = r1/(tt*(tt + 1d0)**(1.-r2) &
        * (1d0 + sqrt(Tgas/r4))**(1.+r2))

  end function recNewV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, iron only, cm3/s
  function recFeV96(Tgas,r1,r2,r3)
    implicit none
    real*8::recFeV96,Tgas,r1,r2,r3,tt

    tt = sqrt(Tgas*1d-4)
    recFeV96 = r1/tt**(r2 + r3 + log10(tt))

  end function recFeV96

  !******************************
  !radiative recombination rates from Verner+96
  ! unit: cm3/s
  function radrec_v96(Tgas,a,b,T0,T1)
    implicit none
    real*8::Tgas,a,b,T0,T1,radrec_v96,iT0

    iT0 = 1d0/T0
    radrec_v96 = a/(sqrt(Tgas*iT0) + (1d0*sqrt(Tgas*iT0))**(1.-b) &
        * (1d0+sqrt(Tgas/T1))**(1+b))

  end function radrec_v96

  !*******************************
  !radiative recombination rates low-temp fit, Verner+96
  ! unit: cm3/s
  function radrec_low_v96(Tgas,a,b,c,d,f)
    implicit none
    real*8::Tgas,a,b,c,d,f,radrec_low_v96,t,invt

    t = Tgas*1d-4
    invt = 1d0/t

    radrec_low_v96 = 1d-12 * (a*invt + b + c*t + d*t**2) &
        * t**(-1.5) * exp(-f*invt)

    radrec_low_v96 = max(0d0,radrec_low_v96)

  end function radrec_low_v96

  !***************************
  !Collisional dissociation rate (cm-3/s) by Martin et al. 1996
  ! H2+H->H+H+H
  !NOTE: the use of this rate is suggested
  ! for high-density regime and in the presence of UV backgrounds.
  ! if necessary it must be included in the reaction file as
  ! H2,H,,H,H,H,,NONE,NONE,dissH2_Martin96(n,Tgas)
  function dissH2_Martin96(n,Tgas)
    use krome_commons
    integer::i
    real*8::n(nspec),Tgas,dissH2_Martin96
    real*8::CDrates,logTv(4),k_CIDm(21,2),k_CID,invT,logT,n_c1,n_c2,n_H
    real*8::logk_h1,logk_h2,logk_l1,logk_l2,logn_c1,logn_c2,p,logk_CID
    real*8::logT2,logT3

    !k_CID = collision-induced dissociation + dissociative tunneling

    !Collisional dissociation of H2
    k_CIDm(:,1) = (/-178.4239d0, -68.42243d0, 43.20243d0, -4.633167d0, &
        69.70086d0, 40870.38d0, -23705.70d0, 128.8953d0, -53.91334d0, &
        5.315517d0, -19.73427d0, 16780.95d0, -25786.11d0, 14.82123d0, &
        -4.890915d0, 0.4749030d0, -133.8283d0, -1.164408d0, 0.8227443d0,&
        0.5864073d0, -2.056313d0/)

    !Dissociative tunneling of H2
    k_CIDm(:,2) = (/-142.7664d0, 42.70741d0, -2.027365d0, -0.2582097d0, &
        21.36094d0, 27535.31d0, -21467.79d0, 60.34928d0, -27.43096d0, &
        2.676150d0, -11.28215d0, 14254.55d0, -23125.20d0, 9.305564d0, &
        -2.464009d0, 0.1985955d0, 743.0600d0, -1.174242d0, 0.7502286d0, &
        0.2358848d0, 2.937507d0/)

    n_H  = get_Hnuclei(n(:))
    logT = log10(Tgas)
    invT = 1.0d0/Tgas
    logT2 = logT*logT
    logT3 = logT2*logT
    logTv = (/1.d0, logT, logT2, logT3/)
    k_CID = 0.d0
    do i=1,2
      logk_h1 = k_CIDm(1,i)*logTv(1) + k_CIDm(2,i)*logTv(2) + &
          k_CIDm(3,i)*logTv(3) + k_CIDm(4,i)*logTv(4) + &
          k_CIDm(5,i)*log10(1.d0+k_CIDm(6,i)*invT)
      logk_h2 = k_CIDm(7,i)*invT
      logk_l1 = k_CIDm(8,i)*logTv(1) + k_CIDm(9,i)*logTv(2) + &
          k_CIDm(10,i)*logTv(3) + k_CIDm(11,i)*log10(1.d0+k_CIDm(12,i)*invT)
      logk_l2 = k_CIDm(13,i)*invT
      logn_c1 = k_CIDm(14,i)*logTv(1) + k_CIDm(15,i)*logTv(2) &
          + k_CIDm(16,i)*logTv(3) + k_CIDm(17,i)*invT
      logn_c2 = k_CIDm(18,i) + logn_c1
      p = k_CIDm(19,i) + k_CIDm(20,i)*exp(-Tgas/1.850d3) &
          + k_CIDm(21,i)*exp(-Tgas/4.40d2)
      n_c1 = 1d1**(logn_c1)
      n_c2 = 1d1**(logn_c2)
      logk_CID = logk_h1 - (logk_h1 - logk_l1) / (1.d0 + (n_H/n_c1)**p) &
          + logk_h2 - (logk_h2 - logk_l2) / (1.d0 + (n_H/n_c2)**p)
      k_CID = k_CID + 1.d1**logk_CID
    enddo

    dissH2_Martin96 = k_CID

  end function dissH2_Martin96

  !**********************
  !adsorpion rate Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014
  function dust_adsorption_rate(nndust,ims,stick,adust2,sqrTgas)
    use krome_constants
    implicit none
    real*8::dust_adsorption_rate,nndust,ims,stick,adust2,sqrTgas

    dust_adsorption_rate = nndust * pi * adust2 &
        * pre_kvgas_sqrt * ims * sqrTgas &
        * stick

  end function dust_adsorption_rate

  !*****************************
  !desorption rate Cazaux+2010, Hocuk+2014
  function dust_desorption_rate(fice,expEice,expEbare)
    implicit none
    real*8::dust_desorption_rate
    real*8::fice,expEice,expEbare,nu0,fbare

    nu0 = 1d12 !1/s
    fbare = 1d0 - fice
    dust_desorption_rate = nu0 * (fbare * expEbare &
        + fice * expEice)

  end function dust_desorption_rate

  !**************************
  function dust_2body_rate(p,invphi,fice,expEice1,expEice2,&
        expEbare1,expEbare2,pesc_ice,pesc_bare)
    use krome_constants
    implicit none
    real*8::fice,expEice1,expEice2,expEbare1,expEbare2,invphi
    real*8::nu0,p,dust_2body_rate,fbare,pesc_ice,pesc_bare

    !no need to calculate this if the dust is not present
    dust_2body_rate = 0d0

    fbare = 1d0-fice
    nu0 = 1d12 ! 1/s
    dust_2body_rate = fbare * (expEbare1 + expEbare2) * pesc_bare &
        + fice * (expEice1 + expEice2) * pesc_ice
    dust_2body_rate = dust_2body_rate * p * nu0 * invphi

  end function dust_2body_rate

  !*************************
  function dust_get_inv_phi(asize2,nndust)
    use krome_commons
    use krome_constants
    implicit none
    real*8::iapp2,dust_get_inv_phi(ndust),asize2(ndust)
    real*8::nndust(ndust),dephi
    integer::i

    iapp2 = (3d-8)**2 !1/cm2
    do i=1,ndust
      dust_get_inv_phi(i) = 0d0
      dephi = (4d0 * nndust(i) * pi * asize2(i))
      if(dephi.le.0d0) cycle
      dust_get_inv_phi(i) = iapp2 / dephi
    end do

  end function dust_get_inv_phi

  !****************************
  !returns an array with the sticking coefficient for each bin
  ! following Hollenbach+McKee 1979
  function dust_stick_array(Tgas,Tdust)
    use krome_commons
    implicit none
    real*8::dust_stick_array(ndust),Tgas,Tdust(ndust)
    real*8::Tg100,Td100
    integer::i

    Tg100 = Tgas * 1d-2
    do i=1,ndust
      Td100 = Tdust(i) * 1d-2
      dust_stick_array(i) = 1d0/(1d0+.4d0*sqrt(Tg100+Td100) &
          + .2d0*Tg100 + 0.08d0*Tg100**2)
    end do

  end function dust_stick_array

  !***************************
  function dust_ice_fraction_array(invphi,nH2O)
    use krome_constants
    use krome_commons
    implicit none
    integer::i
    real*8::dust_ice_fraction_array(ndust)
    real*8::invphi(ndust),nH2O(ndust)

    do i=1,ndust
      dust_ice_fraction_array(i) = min(nH2O(i) * invphi(i), 1d0)
    end do

  end function dust_ice_fraction_array

  !***********************************
  subroutine init_exp_table()
    use krome_commons
    implicit none
    integer::i
    real*8::a

    do i=1,exp_table_na
      a = (i-1)*(exp_table_aMax-exp_table_aMin)/(exp_table_na-1) + exp_table_aMin
      exp_table(i) = exp(-a)
    end do

  end subroutine init_exp_table

  !*****************************
  function get_exp_table(ain,invT)
    use krome_commons
    implicit none
    integer::ia
    real*8::get_exp_table,a,invT,ain
    real*8::x1a,f1,f2

    a = ain*invT
    a = min(a, exp_table_aMax - exp_table_da)

    ia = (a-exp_table_aMin) * exp_table_multa + 1
    ia = max(ia,1)

    x1a = (ia-1)*exp_table_da

    f1 = exp_table(ia)
    f2 = exp_table(ia+1)

    get_exp_table = (a-x1a) * exp_table_multa * (f2-f1) + f1

  end function get_exp_table

  !*****************************
  function get_Ebareice_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice_exp_array(:) = 0d0

  end function get_Ebareice_exp_array

  !*****************************
  function get_Ebareice23_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice23_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice23_exp_array(:) = 0d0

  end function get_Ebareice23_exp_array

  !************************
  !returns the binding energy for ice coated grain (K)
  function get_Ebind_ice()
    use krome_commons
    implicit none
    real*8::get_Ebind_ice(nspec)

    get_Ebind_ice(:) = 0d0

  end function get_Ebind_ice

  !************************
  !returns the binding energy for bare grain (K)
  function get_Ebind_bare()
    use krome_commons
    implicit none
    real*8::get_Ebind_bare(nspec)

    get_Ebind_bare(:) = 0d0

  end function get_Ebind_bare

  !************************
  !returns the index of the parent dust bin (0 if none)
  function get_parent_dust_bin()
    use krome_commons
    implicit none
    integer::get_parent_dust_bin(nspec)

    get_parent_dust_bin(:) = 0

  end function get_parent_dust_bin

  !***************************
  !get the index of the specie name
  function get_index(name)
    use krome_commons
    integer::get_index,i
    character*16::names(nspec)
    character*(*)::name
    names(:) = get_names()
    get_index = -1 !default index
    !loop on species to found the specie named name
    do i=1,nspec
      !when found store and break loop
      if(trim(names(i))== trim(name)) then
        get_index = i !store index
        exit
      end if
    end do

    !error if species not found
    if(get_index<0) then
      print *,"ERROR: can't find the index of ",name
      stop
    end if

  end function get_index

  !************************
  !get electrons by balancing charges
  function get_electrons(n)
    use krome_commons
    implicit none
    real*8::get_electrons,n(nspec)

    get_electrons =  + n(idx_Hj) &
        + n(idx_HEj) &
        + 2d0*n(idx_HEjj) &
        - n(idx_Hk) &
        + n(idx_H2j)
    get_electrons = max(get_electrons,0d0)

  end function get_electrons

  !************************
  !get species charges
  function get_charges()
    use krome_commons
    implicit none
    integer::get_charges(nspec)

    get_charges(1) = 0.d0 	!H
    get_charges(2) = -1.d0 	!E
    get_charges(3) = 1.d0 	!H+
    get_charges(4) = 0.d0 	!HE
    get_charges(5) = 1.d0 	!HE+
    get_charges(6) = 2.d0 	!HE++
    get_charges(7) = -1.d0 	!H-
    get_charges(8) = 0.d0 	!H2
    get_charges(9) = 1.d0 	!H2+
    get_charges(10) = 0.d0 	!CR
    get_charges(11) = 0.d0 	!g
    get_charges(12) = 0.d0 	!Tgas
    get_charges(13) = 0.d0 	!dummy

  end function get_charges

  !************************
  !get species charges
  function get_rnames()
    use krome_commons
    implicit none
    character*50::get_rnames(nrea)

    get_rnames(1) = "H + E -> H+ + E + E"
    get_rnames(2) = "H+ + E -> H"
    get_rnames(3) = "H+ + E -> H"
    get_rnames(4) = "HE + E -> HE+ + E + E"
    get_rnames(5) = "HE+ + E -> HE"
    get_rnames(6) = "HE+ + E -> HE"
    get_rnames(7) = "HE+ + E -> HE++ + E + E"
    get_rnames(8) = "HE++ + E -> HE+"
    get_rnames(9) = "H + E -> H-"
    get_rnames(10) = "H- + H -> H2 + E"
    get_rnames(11) = "H + H+ -> H2+"
    get_rnames(12) = "H + H+ -> H2+"
    get_rnames(13) = "H2+ + H -> H2 + H+"
    get_rnames(14) = "H2 + H+ -> H2+ + H"
    get_rnames(15) = "H2 + H+ -> H2+ + H"
    get_rnames(16) = "H2 + E -> H + H-"
    get_rnames(17) = "H2 + E -> H + H + E"
    get_rnames(18) = "H2 + H -> H + H + H"
    get_rnames(19) = "H- + E -> H + E + E"
    get_rnames(20) = "H- + H -> H + H + E"
    get_rnames(21) = "H- + H -> H + H + E"
    get_rnames(22) = "H- + H+ -> H + H"
    get_rnames(23) = "H- + H+ -> H2+ + E"
    get_rnames(24) = "H2+ + E -> H + H"
    get_rnames(25) = "H2+ + E -> H + H"
    get_rnames(26) = "H2+ + H- -> H + H2"
    get_rnames(27) = "H2 + H2 -> H + H + H2"
    get_rnames(28) = "HE+ + H -> HE + H+"
    get_rnames(29) = "HE + H+ -> HE+ + H"
    get_rnames(30) = "HE + H+ -> HE+ + H"
    get_rnames(31) = "H2 + HE+ -> HE + H + H+"
    get_rnames(32) = "H2 + HE -> H + H + HE"
    get_rnames(33) = "H2 + HE+ -> H2+ + HE"
    get_rnames(34) = "H + H -> H + H+ + E"
    get_rnames(35) = "H + HE -> HE + H+ + E"
    get_rnames(36) = "H -> H+ + E"
    get_rnames(37) = "HE -> HE+ + E"
    get_rnames(38) = "HE+ -> HE++ + E"
    get_rnames(39) = "H- -> H + E"
    get_rnames(40) = "H2 -> H2+ + E"
    get_rnames(41) = "H2+ -> H+ + H"
    get_rnames(42) = "H2+ -> H+ + H+ + E"
    get_rnames(43) = "H2 -> H + H"
    get_rnames(44) = "H2 -> H + H"
    get_rnames(45) = "H+ + E -> H"
    get_rnames(46) = "HE+ + E -> HE"
    get_rnames(47) = "H -> H+ + E"
    get_rnames(48) = "HE -> HE+ + E"
    get_rnames(49) = "H2 -> H + H"
    get_rnames(50) = "H2 -> H+ + H-"
    get_rnames(51) = "H2 -> H2+ + E"
    get_rnames(52) = "H2 -> H + H+ + E"

  end function get_rnames

  !*****************************
  !computes revers kinetics from reaction and
  ! product indexes
  function revKc(Tgas,ridx,pidx)
    implicit none
    real*8::revKc,Tgas
    integer::ridx(:),pidx(:),i

    revKc = 0.d0

    do i=1,size(pidx)
      revKc = revKc + revHS(Tgas,pidx(i))
    end do

    do i=1,size(ridx)
      revKc = revKc - revHS(Tgas,ridx(i))
    end do

  end function revKc

  !*****************************
  !compute H-S for species with index idx
  ! when temperature is Tgas
  function revHS(Tgas,idx)
    use krome_commons
    real*8::revHS,Tgas,Tgas2,Tgas3,Tgas4,invT,lnT,H,S
    real*8::p1(13,7), p2(13,7), Tlim(13,3), p(7)
    integer::idx

    p1(:,:) = 0.d0
    p2(:,:) = 0.d0
    Tlim(:,:) = 0.d0
    Tgas2 = Tgas * Tgas
    Tgas3 = Tgas2 * Tgas
    Tgas4 = Tgas3 * Tgas
    invT = 1d0/Tgas
    lnT = log(Tgas)

    p1(idx_H,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25473.66d0,&
        -0.44668285d0/)
    p1(idx_Hj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184021.488d0,&
        -1.14064664d0/)
    p1(idx_HE,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        0.928723974d0/)
    p1(idx_HEj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        285323.374d0,&
        1.62166556d0/)
    p1(idx_Hk,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        15976.167d0,&
        -1.1390139d0/)
    p1(idx_H2,:)  = (/2.34433112d0,&
        0.00798052075d0,&
        -1.9478151d-05,&
        2.01572094d-08,&
        -7.37611761d-12,&
        -917.935173d0,&
        0.683010238d0/)
    p1(idx_H2j,:)  = (/3.77256072d0,&
        -0.0019574659d0,&
        4.54812047d-06,&
        -2.82152141d-09,&
        5.33969209d-13,&
        178694.654d0,&
        -3.96609192d0/)
    p2(idx_H,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25473.66d0,&
        -0.44668285d0/)
    p2(idx_Hj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184021.488d0,&
        -1.14064664d0/)
    p2(idx_HE,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        0.928723974d0/)
    p2(idx_HEj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        285323.374d0,&
        1.62166556d0/)
    p2(idx_Hk,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        15976.167d0,&
        -1.1390139d0/)
    p2(idx_H2,:)  = (/2.93286575d0,&
        0.000826608026d0,&
        -1.46402364d-07,&
        1.54100414d-11,&
        -6.888048d-16,&
        -813.065581d0,&
        -1.02432865d0/)
    p2(idx_H2j,:)  = (/3.44204765d0,&
        0.000599083239d0,&
        6.69133685d-08,&
        -3.43574373d-11,&
        1.97626599d-15,&
        178650.236d0,&
        -2.79499055d0/)
    Tlim(idx_H,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_Hj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_HE,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_HEj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_Hk,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_H2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_H2j,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)

    if(Tlim(idx,2)==0.d0) then
      revHS = 0.d0
      return
    end if

    !select set of NASA polynomials using temperature
    if(Tlim(idx,1).le.Tgas .and. Tgas.le.Tlim(idx,2)) p(:) = p1(idx,:)
    if(Tlim(idx,2)<Tgas .and. Tgas.le.Tlim(idx,3)) p(:) = p2(idx,:)

    !compute NASA polynomials for enthalpy and enthropy
    H = p(1) + p(2)*0.5d0*Tgas + p(3)*Tgas2/3.d0 + p(4)*Tgas3*0.25d0 + &
        p(5)*Tgas4*0.2d0 + p(6)*invT
    S = p(1)*lnT + p(2)*Tgas + p(3)*Tgas2*0.5d0 + p(4)*Tgas3/3.d0 + &
        p(5)*Tgas4*0.25d0 + p(7)

    revHS = H - S

  end function revHS

  !******************************
  subroutine print_best_flux(n,Tgas,nbestin)
    !print the first nbestin fluxes
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,flux(nrea)
    integer::nbest,idx(nrea),i,nbestin
    character*50::name(nrea)

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
      print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux

  !******************************
  subroutine print_best_flux_frac(n,Tgas,frac)
    !print the first nbestin fluxes
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),frac
    integer::idx(nrea),i
    character*50::name(nrea)

    if(frac>1d0) then
      print *,"ERROR: fraction in krome_print_best_flux should be <=1!"
      stop
    end if

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nrea
      if(flux(idx(i))<flux(idx(1))*frac) exit
      print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux_frac

  !******************************
  subroutine print_best_flux_spec(n,Tgas,nbestin,idx_found)
    !print the first nbestin fluxes for the reactions
    ! that contains the species with index idx_found
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),maxflux
    integer::nbest,idx(nrea),i,nbestin,idx_found
    character*50::name(nrea)
    logical::found

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions
    maxflux = 0d0
    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names
    do i=1,nrea
      found = .false.
      if(arr_r1(i) == idx_found) found = .true.
      if(arr_r2(i) == idx_found) found = .true.
      if(arr_p1(i) == idx_found) found = .true.
      if(arr_p2(i) == idx_found) found = .true.
      if(arr_p3(i) == idx_found) found = .true.
      maxflux = max(maxflux,flux(i))
      if(.not.found) flux(i) = 0d0
    end do

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
      print '(I8,a1,a50,2E17.8)',idx(i)," ",name(idx(i)),flux(idx(i)),&
          flux(idx(i))/maxflux
    end do

  end subroutine print_best_flux_spec

  !*****************************
  function idx_sort(fin)
    !sorting algorithm: requires an array of real values fin
    ! and returns the sorted index list. descending.
    ! bubblesort: not very efficient, replace with what you prefer
    implicit none
    real*8::fin(:),f(size(fin)),ftmp
    integer::idx_sort(size(fin)),n,itmp,i
    logical::found

    f(:) = fin(:) !copy to local

    n = size(f)
    !init indexes
    do i=1,n
      idx_sort(i) = i
    end do

    !loop to sort
    do
      found = .false. !swapped something flag
      do i=2,n
        !> for descending, < for ascending
        if(f(i)>f(i-1)) then
          found = .true.
          !swap real value
          ftmp = f(i)
          f(i) = f(i-1)
          f(i-1) = ftmp
          !swap index
          itmp = idx_sort(i)
          idx_sort(i) = idx_sort(i-1)
          idx_sort(i-1) = itmp
        end if
      end do
      !if nothing swapped exit
      if(.not.found) exit
    end do

  end function idx_sort

  !******************************
  function get_flux(n,Tgas)
    !get the flux k*n*n*... of the rates
    use krome_commons
    implicit none
    integer::i
    integer::r1,r2
    real*8::get_flux(nrea),n(nspec),k(nrea),rrmax,Tgas

    k(:) = coe(n(:))
    rrmax = 0.d0
    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
      r1 = arr_r1(i)
      r2 = arr_r2(i)
      arr_flux(i) = k(i)*n(r1)*n(r2)
    end do
    get_flux(:) = arr_flux(:)

  end function get_flux

  !*****************************
  subroutine load_arrays()
    !load the array containing reactants
    ! and product index
    use krome_commons

    arr_r1(:) = (/1,3,3,4,5,5,5,6,1,7,1,1,9,8,8,8,8,8,7,7,7,7,7&
        ,9,9,9,8,5,4,4,8,8,8,1,1,1,4,5,7,8,9,9,8,8,3,5,1,4,8,8,8&
        ,8/)
    arr_r2(:) = (/2,2,2,2,2,2,2,2,2,1,3,3,1,3,3,2,2,1,2,1&
        ,1,3,3,2,2,7,8,1,3,3,5,4,5,1,4,13,13,13,13,13,13,13,13,13,2,2&
        ,13,13,13,13,13,13/)
    arr_p1(:) = (/3,1,1,5,4,4,6,5,7,8,9,9&
        ,8,9,9,1,1,1,1,1,1,1,9,1,1,1,1,4,5,5,4,1,9,1,4,3,5,6,1,9,3,3&
        ,1,1,1,4,3,5,1,3,9,1/)
    arr_p2(:) = (/2,13,13,2,13,13,2,13&
        ,13,2,13,13,3,1,1,7,1,1,2,1,1,1,2,1,1,8,1,3,1,1,1,1,4,3,3,2,2&
        ,2,2,2,1,3,1,1,13,13,2,2,1,7,2,3/)
    arr_p3(:) = (/2,13,13,2&
        ,13,13,2,13,13,13,13,13,13,13,13,13,2,1,2,2,2,13,13,13,13,13&
        ,8,13,13,13,3,4,13,2,2,13,13,13,13,13,13,2,13,13,13,13,13,13&
        ,13,13,13,2/)

  end subroutine load_arrays

  !********************************************
  subroutine init_anytab2D(filename,x,y,z,xmul,ymul)
    use krome_commons
    implicit none
    character(len=*)::filename
    character(len=60)::row_string
    real*8::x(:),y(:),z(:,:),rout(3),xmul,ymul
    integer::i,j,ios

    !check the size of the X input array
    if(size(x).ne.size(z,1)) then
      print *,"ERROR: in init_anytab2D x size differs from z"
      stop
    end if

    !check the size of the Y input array
    if(size(y).ne.size(z,2)) then
      print *,"ERROR: in init_anytab2D y size differs from z"
      stop
    end if

    if (krome_mpi_rank<=1) print *,"Reading tables from "//trim(filename)

    !open file and check if it exists
    open(krome_nfile,file=trim(filename),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: in init_anytab2D file ",trim(filename)," not found!"
      stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
      read(krome_nfile,'(a)') row_string
      if(row_string(1:1)/="#") exit
    end do

    !check if first line is OK
    if(scan(row_string,",")==0) then
      print *,"ERROR: file "//filename//" should"
      print *," contain the number of rows and "
      print *," columns in the format"
      print *,"  RR, CC"
      print *,row_string
      stop
    end if

    !loop to read file
    do i=1,size(x)
      do j=1,size(y)
        read(krome_nfile,*,iostat=ios) rout(:)
        y(j) = rout(2)
        z(i,j) = rout(3)
      end do
      x(i) = rout(1)
      read(krome_nfile,*,iostat=ios) !skip blanks
      if(ios.ne.0) exit
      if(ios .eq. 5008) exit
    end do
    close(krome_nfile)
 
    xmul = 1d0/(x(2)-x(1))
    ymul = 1d0/(y(2)-y(1))

  end subroutine init_anytab2D

  !******************************
  subroutine test_anytab2D(fname,x,y,z,xmul,ymul)
    use krome_commons
    implicit none
    integer::i,j
    real*8::x(:),y(:),z(:,:),xmul,ymul,xx,yy,zz
    character(len=*)::fname

    open(krome_nfile ,file=fname//".fit",status="replace")
    open(krome_nfile2,file=fname//".org",status="replace")
    do i=1,size(x)
      do j=1,size(y)
        xx = x(i)
        yy = y(j)
        zz = fit_anytab2D(x(:),y(:),z(:,:),xmul,ymul,xx,yy)
        write(krome_nfile,*) xx,yy,zz
        write(krome_nfile2,*) x(i),y(j),z(i,j)
      end do
      write(krome_nfile,*)
      write(krome_nfile2,*)
    end do
    print *,"original file wrote in ",fname//".org"
    print *,"fit test file wrote in ",fname//".fit"

    close(krome_nfile)
    close(krome_nfile2)

  end subroutine test_anytab2D

  !******************************
  function fit_anytab2D(x,y,z,xmul,ymul,xx,yy)
    real*8::fit_anytab2D,x(:),y(:),z(:,:),xmul,ymul,xx,yy
    real*8::zleft(size(x)),zright(size(x)),zl,zr
    integer::ipos,i1,i2

    ipos = (yy-y(1)) * ymul + 1
    i1 = min(max(ipos,1),size(y)-1)
    i2 = i1 + 1
    zleft(:) = z(:,i1)
    zright(:) = z(:,i2)

    zl = fit_anytab1D(x(:),zleft(:),xmul,xx)
    zr = fit_anytab1D(x(:),zright(:),xmul,xx)

    fit_anytab2D = (yy-y(i1))*ymul*(zr-zl)+zl

  end function fit_anytab2D

  !*********************
  function fit_anytab1D(x,z,xmul,xx)
    real*8::fit_anytab1D,x(:),z(:),xmul,xx,p
    integer::ipos,i1,i2

    ipos = (xx-x(1)) * xmul + 1
    i1 = min(max(ipos,1),size(x)-1)
    i2 = i1 + 1

    p = (xx-x(i1)) * xmul

    fit_anytab1D = p * (z(i2) - z(i1)) + z(i1)

  end function fit_anytab1D

  function fit_anytab2D_linlog(x,y,z,xmul,ymul,xx,yy)
    real*8::fit_anytab2D_linlog,x(:),y(:),z(:,:),xmul,ymul,xx,yy
    real*8::zleft(size(x)),zright(size(x)),zl,zr
    integer::ipos,i1,i2

    ipos = (yy-y(1)) * ymul + 1
    i1 = min(max(ipos,1),size(y)-1)
    i2 = i1 + 1
    zleft(:) = z(:,i1)
    zright(:) = z(:,i2)

    zl = fit_anytab1D_linlog(x(:),zleft(:),xmul,xx)
    zr = fit_anytab1D_linlog(x(:),zright(:),xmul,xx)

    fit_anytab2D_linlog = (yy-y(i1))*ymul*(zr-zl)+zl

  end function fit_anytab2D_linlog

  !*********************
  function fit_anytab1D_linlog(x,z,xmul,xx)
    real*8::fit_anytab1D_linlog,x(:),z(:),xmul,xx,p,z2,z1
    integer::ipos,i1,i2

    ipos = (xx-x(1)) * xmul + 1
    i1 = min(max(ipos,1),size(x)-1)
    i2 = i1 + 1

    p = (xx-x(i1)) * xmul

    z2 = z(i2)
    z1 = z(i1)
    if(z1<0d0 .and. z2<0d0) then
      z1 = log10(-z1)
      z2 = log10(-z2)
      fit_anytab1D_linlog = -1d1**(p * (z2 - z1) + z1)
      return
    end if

    if(z1>0d0 .and. z2>0d0) then
      z1 = log10(z1)
      z2 = log10(z2)
      fit_anytab1D_linlog = 1d1**(p * (z2 - z1) + z1)
      return
    end if

    fit_anytab1D_linlog = (p * (z2 - z1) + z1)

  end function fit_anytab1D_linlog

  !*****************************
  !spline interpolation at t using array  x,y (size n) as data
  function fspline(x,y,t)
    implicit none
    real*8::fspline,x(:),y(:),b(size(x)),c(size(x)),d(size(x)),t
    integer::n

    n = size(x)
    call spline(x(:),y(:),b(:),c(:),d(:),n)
    fspline = ispline(t,x(:),y(:),b(:),c(:),d(:),n)

  end function fspline

  !*******************************+
  subroutine spline(x, y, b, c, d, n)
    !======================================================================
    !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
    !  for cubic spline interpolation
    !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
    !  for  x(i) <= x <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    !  input..
    !  x = the arrays of data abscissas (in strictly increasing order)
    !  y = the arrays of data ordinates
    !  n = size of the arrays xi() and yi() (n>=2)
    !  output..
    !  b, c, d  = arrays of spline coefficients
    !  comments ...
    !  spline.f90 program is based on fortran version of program spline.f
    !  the accompanying function fspline can be used for interpolation
    !======================================================================
    implicit none
    integer::n
    real*8::x(n), y(n), b(n), c(n), d(n)
    integer::i, j, gap
    real*8::h

    gap = n-1

    !check input
    if(n<2) return
    if(n<3) then
      b(1) = (y(2)-y(1))/(x(2)-x(1)) !linear interpolation
      c(1) = 0d0
      d(1) = 0d0
      b(2) = b(1)
      c(2) = 0d0
      d(2) = 0d0
      return
    end if

    !step 1: preparation
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
      d(i) = x(i+1) - x(i)
      b(i) = 2d0*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i))/d(i)
      c(i) = c(i+1) - c(i)
    end do

    ! step 2: end conditions
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0d0
    c(n) = 0d0
    if(n.ne.3) then
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if

    ! step 3: forward elimination
    do i = 2, n
      h = d(i-1)/b(i-1)
      b(i) = b(i) - h*d(i-1)
      c(i) = c(i) - h*c(i-1)
    end do

    ! step 4: back substitution
    c(n) = c(n)/b(n)
    do j = 1, gap
      i = n-j
      c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do

    ! step 5: compute spline coefficients
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2d0*c(n))
    do i = 1, gap
      b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2d0*c(i))
      d(i) = (c(i+1) - c(i))/d(i)
      c(i) = 3d0*c(i)
    end do
    c(n) = 3d0*c(n)
    d(n) = d(n-1)
  end subroutine spline

  !*******************************
  function ispline(u, x, y, b, c, d, n)
    !======================================================================
    ! function ispline evaluates the cubic spline interpolation at point z
    ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
    ! where  x(i) <= u <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    ! input..
    ! u       = the abscissa at which the spline is to be evaluated
    ! x, y    = the arrays of given data points
    ! b, c, d = arrays of spline coefficients computed by spline
    ! n       = the number of data points
    ! output:
    ! ispline = interpolated value at point u
    !=======================================================================
    implicit none
    real*8::ispline
    integer::n
    real*8::u, x(n), y(n), b(n), c(n), d(n)
    integer::i, j, k
    real*8::dx

    ! if u is ouside the x() interval take a boundary value (left or right)
    if(u<=x(1)) then
      ispline = y(1)
      return
    end if

    if(u>=x(n)) then
      ispline = y(n)
      return
    end if

    ! binary search for for i, such that x(i) <= u <= x(i+1)
    i = 1
    j = n+1
    do while (j>i+1)
      k = (i+j)/2
      if(u<x(k)) then
        j=k
      else
        i=k
      end if
    end do

    ! evaluate spline interpolation
    dx = u - x(i)
    ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

  end function ispline

  !********************************
  !H2 formation on dust using Jura rate
  !dust2gas_ratio in terms of D_solar
  !Usually D/D_sol = Z/Z_sol
  function H2_dustJura(n,Tgas)
    use krome_commons
    use krome_user_commons
    implicit none
    real*8::n(nspec),H2_dustJura
    real*8::ntot,Tgas, fa, S, Tdust

    ntot = get_Hnuclei(n(:))

    Tdust = MIN(Tgas, 3d1)
    fa = 1d0 / (1d0 + 1d4 * EXP(-6d2/Tgas))
    S  = 1d0 / (1d0 + 4d-2*(Tgas+Tdust)**0.5 + 2d-3*Tgas + 8d-6*Tgas**2)

    H2_dustJura = n(idx_H)*ntot*3.d-17 * fa * S * (Tgas/1d2)**0.5  &
                * dust2gas_ratio*clump_factor

  end function H2_dustJura

end module krome_subs

!############### MODULE ##############
module krome_photo
contains

  !*************************
  !get the intensity of the photon flux at
  ! a given energy in eV.
  ! returned value is in eV/cm2/s/Hz
  function get_photoIntensity(energy)
    use krome_commons
    implicit none
    real*8::get_photoIntensity,energy
    integer::i

    !check if requested energy is lower than the lowest limit
    if(energy<photoBinEleft(1)) then
      get_photoIntensity = 0d0 !photoBinJ(1)
      return
    end if

    !check if requested energy is greater that the the largest limit
    if(energy>photoBinEright(nPhotoBins)) then
      get_photoIntensity = 0d0 !photoBinJ(nPhotoBins)
      return
    end if

    !look for the interval
    do i=1,nPhotoBins
      if(photoBinEleft(i).le.energy .and. photoBinEright(i).ge.energy) then
        get_photoIntensity = photoBinJ(i)
        return
      end if
    end do

    !error if nothing found
    print *,"ERROR: no interval found in get_photoIntensity"
    print *,"energy:",energy,"eV"
    stop !halt program

  end function get_photoIntensity

  !*********************
  !initialize/tabulate the bin-based xsecs
  subroutine init_photoBins()
    use krome_commons
    implicit none
    integer::i,j
    real*8::energy_eV,kk,energyL,energyR

    if(photoBinEmid(nPhotoBins)==0d0) then
      print *,"ERROR: when using photo bins you must define"
      print *," the energy interval in bins!"
      stop
    end if

    call load_xsec(trim(krome_datafolder)//"swri_H-__H_E.dat", xsec39_val, xsec39_Emin, xsec39_n, xsec39_idE)
    call load_xsec(trim(krome_datafolder)//"swri_H2__H2+_E.dat", xsec40_val, xsec40_Emin, xsec40_n, xsec40_idE)
    call load_xsec(trim(krome_datafolder)//"leiden_H2+__H+_H.dat", xsec41_val, xsec41_Emin, xsec41_n, xsec41_idE)

    !tabulate the xsecs into a bin-based array
    do j=1,nPhotoBins
      energyL = photoBinEleft(j)
      energyR = photoBinEright(j)
      energy_eV = photoBinEmid(j) !energy of the bin in eV

      !H -> H+ + E
      kk = 0d0
      if(energy_eV>1.360d+01.and.energy_eV<5.000d+04) then
        kk =  sigma_v96(energy_ev, 4.298d-01, 5.475d+04, 3.288d+01, 2.963d+00, 0.000d+00, 0.000d+00, 0.000d+00)
      endif
      photoBinJTab(1,j) = kk

      !HE -> HE+ + E
      kk = 0d0
      if(energy_eV>2.459d+01.and.energy_eV<5.000d+04) then
        kk =  sigma_v96(energy_ev, 1.361d+01, 9.492d+02, 1.469d+00, 3.188d+00, 2.039d+00, 4.434d-01, 2.136d+00)
      endif
      photoBinJTab(2,j) = kk

      !HE+ -> HE++ + E
      kk = 0d0
      if(energy_eV>5.442d+01.and.energy_eV<5.000d+04) then
        kk =  sigma_v96(energy_ev, 1.720d+00, 1.369d+04, 3.288d+01, 2.963d+00, 0.000d+00, 0.000d+00, 0.000d+00)
      endif
      photoBinJTab(3,j) = kk

      !H- -> H + E
      kk = 0d0
      if(energy_eV>0.755d0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec39_val(:), xsec39_Emin,xsec39_idE)
      photoBinJTab(4,j) = kk

      !H2 -> H2+ + E
      kk = 0d0
      if(energy_eV>15.4d0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec40_val(:), xsec40_Emin,xsec40_idE)
      photoBinJTab(5,j) = kk

      !H2+ -> H+ + H
      kk = 0d0
      if(energy_eV>2.65d0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec41_val(:), xsec41_Emin,xsec41_idE)
      photoBinJTab(6,j) = kk

      !H2+ -> H+ + H+ + E
      kk = 0d0
      if(energy_eV>3.d1.and.energy_eV<9.d1) kk = 10**(-16.926-4.528d-2*energy_eV+2.238d-4*energy_eV**2+4.245d-7*energy_eV**3)
      photoBinJTab(7,j) = kk

      !H2 -> H + H
      kk = 0d0
      if(energy_eV>14.159d0.and.energy_eV<17.7d0) kk = H2_sigmaLW(energy_eV)
      photoBinJTab(8,j) = kk

    end do

    !save interpolated xsecs to file
    call save_xsec(trim(krome_datafolder)//"swri_H-__H_E.interp",4)
    call save_xsec(trim(krome_datafolder)//"swri_H2__H2+_E.interp",5)
    call save_xsec(trim(krome_datafolder)//"leiden_H2+__H+_H.interp",6)

    !energy tresholds (eV)
    photoBinEth(1) = 1.360d+01 !H -> H+ + E
    photoBinEth(2) = 2.459d+01 !HE -> HE+ + E
    photoBinEth(3) = 5.442d+01 !HE+ -> HE++ + E
    photoBinEth(4) = 0.755d0 !H- -> H + E
    photoBinEth(5) = 15.4d0 !H2 -> H2+ + E
    photoBinEth(6) = 2.65d0 !H2+ -> H+ + H
    photoBinEth(7) = 3.d1 !H2+ -> H+ + H+ + E
    photoBinEth(8) = 14.159d0 !H2 -> H + H

  end subroutine init_photoBins

  !**********************
  !save xsecs with index idx to file
  subroutine save_xsec(fname,idx)
    use krome_commons
    implicit none
    character(len=*)::fname
    integer::idx,j
    real*8::energy_eV

    open(krome_nfile,file=trim(fname),status="replace")
    do j=1,nPhotoBins
      energy_eV = photoBinEmid(j) !energy of the bin in eV
      write(krome_nfile,*) energy_eV, photoBinJTab(idx,j)
    end do
    close(krome_nfile)

  end subroutine save_xsec

  !**********************
  !compute integrals to derive phtorates (thin)
  subroutine calc_photoBins()
    use krome_commons
    implicit none
    real*8::n(nspec)

    n(:) = 0d0
    call calc_photoBins_thick(n)

  end subroutine calc_photoBins

  !**********************
  !compute integrals to derive phtorates (thick)
  subroutine calc_photoBins_thick(n)
    use krome_commons
    use krome_constants
    use krome_subs
    implicit none
    integer::i,j
    real*8::dE,kk,Jval,E,Eth,n(:),ncol(nmols),tau
    !init rates and heating
    photoBinRates(:) = 0d0 !1/s/Hz
    photoBinHeats(:) = 0d0 !eV/s/Hz
    GHabing_thin = 0d0 !habing flux
    !loop on energy bins
    do j=1,nPhotoBins
      dE = photoBinEdelta(j) !energy interval, eV
      E = photoBinEmid(j) !energy of the bin in eV
      Jval = photoBinJ(j) !radiation intensity eV/s/cm2/sr/Hz
      if(E>=6d0.and.E<=13.6)then
        GHabing_thin = GHabing_thin + Jval * dE
      endif
      tau = 0d0
      !loop on reactions
      do i=1,nPhotoRea
        Eth = photoBinEth(i) !reaction energy treshold, eV
        if(E>Eth) then
          !approx bin integral
          kk = 4d0*pi*photoBinJTab(i,j)*Jval/E*dE * exp(-tau)
          photoBinRates(i) = photoBinRates(i) + kk
          photoBinHeats(i) = photoBinHeats(i) + kk*(E-Eth)
        else
         !if (i .eq. 6) print*, 'term', j, ': ', '0.0000000000000000'
        end if
      end do
    end do
!    print*, 'photobinrate(6)', photobinrates(6)
    !Final Habing flux
    GHabing_thin = GHabing_thin * 4d0 * pi / (1.6d-3) / planck_eV * eV_to_erg

    !converts to 1/s
    photoBinRates(:) = photoBinRates(:) * iplanck_eV

    !converts to erg/s
    photoBinHeats(:) = photoBinHeats(:) * iplanck_eV * eV_to_erg

!    if (.true.) then
!      print*, "here computing rate ..", photobinj(6), photobinrates(6)
!    endif

  end subroutine calc_photoBins_thick

  !********************
  function sigma_v96(energy_eV,E0,sigma_0,ya,P,yw,y0,y1)
    !Verner+96 cross section fit (cm2)
    real*8::sigma_v96,energy_eV,sigma_0,Fy,yw,x,y,E0
    real*8::y0,y1,ya,P
    x = energy_eV/E0 - y0
    y = sqrt(x**2 + y1**2)
    Fy = ((x - 1.d0)**2 + yw**2) *  y**(0.5*P-5.5) &
        * (1.d0+sqrt(y/ya))**(-P)
    sigma_v96 = 1d-18 * sigma_0 * Fy !cm2
  end function sigma_v96

  !********************
  function heat_v96(energy_eV,Eth,E0,sigma_0,ya,P,yw,y0,y1)
    !Heating with Verner+96 cross section fit (cm2*eV)
    use krome_constants
    real*8::heat_v96,energy_eV,sigma_0,Fy,yw,x,y,E0,Eth
    real*8::y0,y1,ya,P
    x = energy_eV/E0 - y0
    y = sqrt(x**2 + y1**2)
    Fy = ((x - 1.d0)**2 + yw**2) *  y**(0.5*P-5.5) &
        * (1.d0+sqrt(y/ya))**(-P)
    heat_v96 = 1d-18 * sigma_0 * Fy * (energy_eV - Eth) !cm2*eV
  end function heat_v96

  !************************
  !load the xsecs from file
  subroutine load_xsec(fname,xsec_val,xsec_Emin,xsec_n,xsec_idE)
    use krome_commons
    implicit none
    real*8,allocatable::xsec_val(:)
    real*8::xsec_Emin,xsec_dE,xsec_val_tmp(int(1e6)),rout(2)
    real*8::xsec_E_tmp(size(xsec_val_tmp)),xsec_idE,diff
    integer::xsec_n,ios
    character(*)::fname

    !if file already loaded skip subroutine
    if(allocated(xsec_val)) return

    xsec_n = 0 !number of lines found
    !open file
    open(krome_nfile,file=fname,status="old",iostat=ios)
    !check if file exists
    if(ios.ne.0) then
      print *,"ERROR: problems loading "//fname
      stop
    end if

    !read file line-by-line
    do
      read(krome_nfile,*,iostat=ios) rout(:) !read line
      if(ios<0) exit !eof
      if(ios/=0) cycle !skip blanks
      xsec_n = xsec_n + 1 !increase line number
      xsec_val_tmp(xsec_n) = rout(2) !read xsec value cm2
      xsec_E_tmp(xsec_n) = rout(1) !read energy value eV
      !compute the dE for the first interval
      if(xsec_n==2) xsec_dE = xsec_E_tmp(2)-xsec_E_tmp(1)
      !check if all the intervals have the same spacing
      if(xsec_n>2) then
        diff = xsec_E_tmp(xsec_n)-xsec_E_tmp(xsec_n-1)
        if(abs(diff/xsec_dE-1d0)>1d-6) then
          print *,"ERROR: spacing problem in file "//fname
          print *," energy points should be equally spaced!"
          print *,"Point number: ",xsec_n
          print *,"Found ",diff
          print *,"Should be",xsec_dE
          stop
        end if
      end if
    end do
    close(krome_nfile)

    !store the minimum energy
    xsec_Emin = xsec_E_tmp(1)
    !allocate the array with the values
    allocate(xsec_val(xsec_n))
    !copy the values from the temp array to the allocated one
    xsec_val(:) = xsec_val_tmp(1:xsec_n)
    !store the inverse of the delta energy
    xsec_idE = 1d0 / xsec_dE

  end subroutine load_xsec

  !**********************
  !return averaged xsec in the energy range [energyL,energyR]
  ! units: eV, cm2
  function xsec_interp(energyL,energyR,xsec_val,xsec_Emin,xsec_idE)
    implicit none
    real*8::xsec_interp,E0,xsecA,dE
    real*8::energy,xsec_val(:),xsec_Emin,xsec_idE,energyL,energyR
    integer::xsec_n,idx

    !xsec energy step (regular grid)
    dE = 1d0/xsec_idE
    xsecA = 0d0 !init integrated xsec
    !loop on xsec vals
    do idx=1,size(xsec_val)
      energy = (idx-1)*dE+xsec_Emin
      !if xsec energy in the interval compute area
      if(energy>=energyL .and. energy<=energyR) xsecA = xsecA &
          + xsec_val(idx)*dE
    end do

    !compute averaged xsec for the flux energy range
    xsec_interp = xsecA/(energyR-energyL)

  end function xsec_interp

  !**********************
  !linear interpolation for the photo xsec
  function xsec_interp_mid(energy,xsec_val,xsec_Emin,xsec_n,xsec_idE)
    implicit none
    real*8::xsec_interp_mid,E0
    real*8::energy,xsec_val(:),xsec_Emin,xsec_idE
    integer::xsec_n,idx

    xsec_interp_mid = 0d0
    !retrive index
    idx = (energy-xsec_Emin) * xsec_idE + 1

    !lower bound
    E0 = xsec_Emin + (idx-1)/xsec_idE

    !out of the limits is zero
    if(idx<1.or.idx>xsec_n-1) return

    !linear interpolation
    xsec_interp_mid = (energy-E0) * xsec_idE &
        * (xsec_val(idx+1)-xsec_val(idx)) + xsec_val(idx)

    !avoid negative xsec values when outside the limits
    xsec_interp_mid = max(xsec_interp_mid,0d0)

  end function xsec_interp_mid

  !************************
  function H2_sigmaLW(energy_eV)
    !H2 direct photodissociation in the Lyman-Werner bands
    ! cross-section in cm^2 fit by Abel et al. 1997 of
    ! data by Allison&Dalgarno 1969
    use krome_commons
    real*8::H2_sigmaLW,energy_eV
    real*8::sL0,sW0,sL1,sW1,fact

    !initialization
    sL0 = 0d0
    sL1 = 0d0
    sW0 = 0d0
    sW1 = 0d0

    if(energy_eV>14.675.and.energy_eV<16.820)then
      sL0 = 1d-18*1d1**(15.1289-1.05139*energy_eV)
    elseif(energy_eV>16.820.and.energy_eV<17.6d0)then
      sL0 = 1d-18*1d1**(-31.41d0+1.8042d-2*energy_eV**3-4.2339d-5*energy_eV**5)
    endif

    if(energy_eV>14.675d0.and.energy_eV<17.7d0)then
      sW0 = 1d-18*1d1**(13.5311d0-0.9182618*energy_eV)
    endif

    if(energy_eV>14.159d0.and.energy_eV<15.302d0)then
      sL1 = 1d-18*1d1**(12.0218406d0-0.819429*energy_eV)
    elseif(energy_eV>15.302d0.and.energy_eV<17.2d0)then
      sL1 = 1d-18*1d1**(16.04644d0-1.082438*energy_eV)
    endif

    if(energy_eV>14.159d0.and.energy_eV<17.2d0)then
      sW1 = 1d-18*1d1**(12.87367-0.85088597*energy_eV)
    endif

    fact = 1d0/(phys_orthoParaRatio+1d0)

    H2_sigmaLW = fact*(sL0+sW0)+(1d0-fact)*(sL1+sW1)

  end function H2_sigmaLW

end module krome_photo

!############### MODULE ##############
module krome_tabs
contains

  subroutine make_ktab()
    !build the tabs from coefficients
    use krome_commons
    use krome_subs
    use krome_photo
    integer::i,j,ierror,kwarnup(nrea),kwarndown(nrea),pblock
    real*8::kk(nrea),valmax,n(nspec)
    logical::is_rank_zero
    !temperature limits
    ktab_logTlow = log10(phys_Tfloor)
    ktab_logTup = log10(1d9)

    is_rank_zero = (krome_mpi_rank<=1)

    !loop to create tabs (it may take a while)
    valmax = 1d0
    ierror = 0 !error count
    pblock = ktab_n/10 !ouput frequency
    if(is_rank_zero) print *,"KROME: creating tabs..."
    kwarnup(:) = 0 !store warnings
    kwarndown(:) = 0 !store warnings
    !loop on temperatures
    do i=1,ktab_n
      if(mod(i,pblock)==0 .and. is_rank_zero) print *,i/pblock*10,"%"
      ktab_T(i) = 1d1**((i-1)*(ktab_logTup-ktab_logTlow)/(ktab_n-1)&
          +ktab_logTlow)
      n(:) = 0.d0
      n(idx_Tgas) = ktab_T(i)
      kk(:) = coe(n(:))
      !check for errors or discrepancies
      if((maxval(kk)>valmax.or.minval(kk)<0.d0)) then
        ierror = ierror + 1
        if(ierror==1.and. is_rank_zero) print '(a16,a5,2a11)',"",&
            "idx","Tgas","rate"
        do j=1,nrea
          if(kk(j)>valmax .and. kwarnup(j)==0) then
            kwarnup(j) = 1
            if(is_rank_zero) print '(a16,I5,2E11.3)', "WARNING: k>1.",&
                j,ktab_T(i),kk(j)
          end if
          if(kk(j)<0.d0 .and. kwarndown(j)==0) then
            kwarndown(j) = 1
            if(is_rank_zero) print *,"WARNING: k<0.d0",j,ktab_T(i),kk(j)
          end if
        end do
      end if
      ktab(:,i) = kk(:)
    end do

    !store inverse values of deltaT to speed-up at runtime
    do i=1,ktab_n-1
      inv_ktab_T(i) = 1.d0 / (ktab_T(i+1)-ktab_T(i))
    end do

    !store inverse to go fast at runtime
    inv_ktab_idx = 1.d0 / (ktab_logTup - ktab_logTlow) * (ktab_n - 1)

  end subroutine make_ktab

  !************************
  subroutine check_tabs()
    use krome_commons
    use krome_subs
    integer::i,j,pblock,ii
    real*8::kk(nrea),kktab(nrea),Tgas,kmax,n(nspec),kold(nrea),dk
    logical::is_rank_zero

    is_rank_zero = (krome_mpi_rank<=1)

    pblock = ktab_n/10 !write % every 10
    if(is_rank_zero) print *,"KROME: checking tabs..."
    !loop on tabs
    do i=1,ktab_n
      if(mod(i,pblock)==0.and.is_rank_zero) print *,i/pblock*10,"%" !output
      Tgas = 1d1**((i-1)*(ktab_logTup-ktab_logTlow)/(ktab_n-1)+ktab_logTlow)
      n(:) = 0.d0 !rates do not depends on densities
      n(idx_Tgas) = Tgas !rates depend on temperature
      kk(:) = coe(n(:)) !get rates
      kktab(:) = coe_tab(n(:)) !get rates from tabs
      kold(:) = 0.d0 !old rate value to skip discontinuities
      !loop on reactions
      do j=1,nrea
        kmax = kk(j)
        if(kmax>0.d0.and.kk(j)>0.d0) then
          dk = abs(kk(j)-kold(j))/(kold(j)+1d-40)
          if(abs(kk(j)-kktab(j))/kmax>1d-1.and.kmax>1d-12.and.dk<1d-1) then
            if(is_rank_zero) then
              print *,"ERROR: wrong rate tables"
              print *,"Rate index:",j
              print *,"Temperature:",Tgas
              print *,"Rate values:",kk(j),kktab(j)
              print *,"Error:",abs(kk(j)-kktab(j))/kmax,&
                  "(must be close to zero)"

              !dump graph
              open(krome_nfile,file=trim(krome_datafolder)//"KROME_TAB_DUMP.dat",status="replace")
              do ii=1,ktab_n
                Tgas = 1d1**((ii-1)*(ktab_logTup-ktab_logTlow)/(ktab_n-1)&
                    +ktab_logTlow)
                n(idx_Tgas) = Tgas !rates depend on temperature
                kk(:) = coe(n(:))
                kktab(:) = coe_tab(n(:))
                write(krome_nfile,'(99E12.3e3)') Tgas,kk(j),kktab(j)
              end do
              close(krome_nfile)
              print *,"Graph dump to KROME_TAB_DUMP.dat"
              print *,"gnuplot command:"
              print *," plot 'KROME_TAB_DUMP.dat' w l, '' u 1:3"
              stop
            end if
          end if
        end if
      end do
      kold(:) = kk(:)
    end do
    if(is_rank_zero) print *,"KROME: tabs are ok!"

  end subroutine check_tabs

  !***********************+
  function coe_tab(n)
    !interface to tabs
    use krome_subs
    use krome_constants
    use krome_commons
    use krome_user_commons
    integer::idx,j
    real*8::Tgas, coe_tab(nrea),n(nspec),small
    real*8::krome_fshieldHM,krome_fshieldH2,log10Tgas,asav2,asav3,asav0,asav1,asav6,asav7,asav4,asav5,a11,a10,a12,a1,a3,a2,a5,T,a7,a6,a9,a8,krome_fshieldH2j,bsav1,bsav0,bsav3,bsav2,bsav5,bsav4,bsav7,bsav6,a4

    Tgas = max(n(idx_Tgas),phys_Tcmb)
    small = 0d0

    T = Tgas
    a1 = 1.3500e-09
    a2 = 9.8493e-02
    a3 = 3.2852e-01
    a4 = 5.5610e-01
    a5 = 2.7710e-07
    a6 = 2.1826e+00
    a7 = 6.1910e-03
    a8 = 1.0461e+00
    a9 = 8.9712e-11
    a10 = 3.0424e+00
    a11 = 3.2576e-14
    a12 = 3.7741e+00
    log10Tgas = log10(Tgas)
    asav0 = -1.9153214d2
    asav1 =  4.0129114d2
    asav2 = -3.7446991d2
    asav3 =  1.9078410d2
    asav4 = -5.7263467d1
    asav5 =  1.0133210d1
    asav6 = -9.8012853d-1
    asav7 =  4.0023414d-2
    bsav0 = -8.8755774d3
    bsav1 =  1.0081246d4
    bsav2 = -4.8606622d3
    bsav3 =  1.2889659d3
    bsav4 = -2.0319575d2
    bsav5 =  1.9057493d1
    bsav6 = -9.8530668d-1
    bsav7 =  2.1675387d-2
    krome_fshieldH2  =  shield_dust(n,Tgas,3.74d0)*krome_fshield(n,Tgas)
    krome_fshieldHM  =  shield_dust(n,Tgas,0.5d0)
    krome_fshieldH2j =  shield_dust(n,Tgas,1.9d0)

    idx = (log10(Tgas)-ktab_logTlow) * inv_ktab_idx + 1
    idx = max(idx,1)
    idx = min(idx,ktab_n-1)
    coe_tab(:) = 0d0
    do j=1,nrea
      coe_tab(j) = (Tgas-ktab_T(idx)) * inv_ktab_T(idx) *&
          (ktab(j,idx+1)-ktab(j,idx)) + ktab(j,idx)
    end do

    coe_tab(36) = small + (photoBinRates(1))
    coe_tab(37) = small + (photoBinRates(2))
    coe_tab(38) = small + (photoBinRates(3))
    coe_tab(39) = small + (photoBinRates(4))
    coe_tab(40) = small + (photoBinRates(5))
    coe_tab(41) = small + (photoBinRates(6))
    coe_tab(42) = small + (photoBinRates(7))
    coe_tab(43) = small + (photoBinRates(8))
    coe_tab(44) = H2_solomonLW_ramses(user_myH2_dissociation,user_myfluxLW)*krome_fshieldH2
    coe_tab(45) = small + (H_recombination_on_dust(n,Tgas))
    coe_tab(46) = small + (He_recombination_on_dust(n,Tgas))
    coe_tab(47) = small + (4.6d-1*user_crate)
    coe_tab(48) = small + (5.d-1*user_crate)
    coe_tab(49) = small + (1d-1*user_crate)
    coe_tab(50) = small + (3d-4*user_crate)
    coe_tab(51) = small + (9.3d-1*user_crate)
    coe_tab(52) = small + (9.3d-1*user_crate)

    coe_tab(39)=coe_tab(39)*krome_fshieldHM
    coe_tab(40)=coe_tab(40)*krome_fshieldH2
    coe_tab(41)=coe_tab(41)*krome_fshieldH2j
    coe_tab(42)=coe_tab(42)*krome_fshieldH2j
    coe_tab(43)=coe_tab(43)*krome_fshieldH2

  end function coe_tab

end module krome_tabs


!############### MODULE ##############
module KROME_coolingGH
! ------------------------------------------------------------
!
! This module: Cooling and Heating Functions, table
! reader for Gnedin and Hollon cooling tables
! Language:  Fortran 77
!
!  UPDATED by Troels Haugboelle to Fortran 90, explicit kinds,
!  and encapsulated in a module
!
!  Copyright (c) 2012 Nick Gnedin
!  All rights reserved.
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions
!  are met:
!
!  Redistributions of source code must retain the above copyright
!  notice, this list of conditions and the following disclaimer.
!
!  Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer in the
!  documentation and/or other materials provided with the distribution.
!
!  Neither the name of Nick Gnedin nor the names of any contributors
!  may be used to endorse or promote products derived from this software
!  without specific prior written permission.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
!  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
!  OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ------------------------------------------------------------
private
!  Table dimensions
integer,parameter::NT=81, NX=13, NP1=24, NP2=21, NP3=16, ND=3789
!  Number of components per T-D bin
integer,parameter:: NC=6, NICH=12, NRCH=13
!  Mode of table lookup
integer::mode
!  Boundaries
integer::np(3)
!  Data and index blocks
integer:: indx(NP1,NP2,NP3)
real*4::data(NC,NT,NX,ND)
! Indices
real*8::altval(NT), altmin, altstp, &
    xval(NX), xmin, xmax, xstp, &
    qmin(3), qmax(3), qstp(3)
public :: frtInitCF, frtCFCache, frtCFGetLn, frtGetCF
public :: NICH, NRCH
contains

subroutine frtInitCF(m,fname)
implicit none
integer                      :: m
character(len=*), intent(in) :: fname
!  Internally used unit number
integer, parameter :: IOCF=97
integer :: i, j, k, id, ix, it, ic, lt, ld, lp1, lp2, lp3, lp4, lx
real*4 :: q1, q2
real*4 :: altval4(NT), xmin4, xmax4, qmin4(3), qmax4(3)

mode = m

open(unit=IOCF, file=fname, status='old', form='unformatted', err=100)
read(IOCF,err=100) lt, ld, lp1, lp2, lp3, lp4, &
    (qmin4(j),j=1,3), q1, (qmax4(j),j=1,3), q2, lx, xmin4, xmax4

if(lt.ne.NT .or. ld.ne.ND .or. lx.ne.NX .or. lp1.ne.NP1 .or. &
    lp2.ne.NP2 .or. lp3.ne.NP3 .or. lp4.ne.1 .or. ld.eq.0) then
write(0,*) 'RT::InitCF: fatal error, corrupted table:'
write(0,*) '> NT= in file: ', lt, ' in code: ', NT
write(0,*) '> NX= in file: ', lx, ' in code: ', NX
write(0,*) '> ND= in file: ', ld, ' in code: ', ND
write(0,*) '> NP1= in file: ', lp1, ' in code: ', NP1
write(0,*) '> NP2= in file: ', lp2, ' in code: ', NP2
write(0,*) '> NP3= in file: ', lp3, ' in code: ', NP3
write(0,*) '> NP4= in file: ', lp4, ' in code: ', 1
close(IOCF)
m = -1
stop
end if

qmin=qmin4; qmax=qmax4; xmin=xmin4; xmax=xmax4

np(1) = lp1
np(2) = lp2
np(3) = lp3

do i=1,3
if(np(i) .gt. 1) then
  qstp(i) = (qmax(i)-qmin(i))/(np(i)-1)
else
  qstp(i) = 1.0
end if
end do

xstp = (xmax-xmin)/(NX-1)
do i=1,NX
xval(i) = xmin + xstp*(i-1)
end do

read(IOCF,err=100) (altval4(i),i=1,NT)
altval = altval4
!  Internally use natural log
do i=1,NT
altval(i) = altval(i)*log(10.0)
end do
altmin = altval(1)
altstp = altval(2) - altval(1)

read(IOCF,err=100) (((indx(i,j,k),i=1,lp1),j=1,lp2),k=1,lp3)

do id=1,ld
read(IOCF,err=100) (((data(ic,it,ix,id),ic=1,NC),it=1,NT),ix=1,NX)
end do

do id=1,ND
do ix=1,NX
  do it=1,NT
    do ic=1,NC
      data(ic,it,ix,id) = log(1d-37+abs(data(ic,it,ix,id)))
    end do
  end do
end do
end do

close(IOCF)

if(.false.) then
write(6,*) 'RT::InitCF: Table size = ', &
    NC*(NT*NX*ND/256/1024) + (NP1*NP2*NP3/256/1024), ' MB'
endif

m = 0

return

100   m = -1

end subroutine frtInitCF
!
!  Decode the interpolated function
!
!#define IXL      ich(3)
!#define IXU      ich(4)
!#define IPP(j)   ich(4+j)
!#define WXL      rch(6)
!#define WXU      rch(7)
!#define WPL(j)   rch(7+j)
!#define WPS(j)   rch(10+j)
!
subroutine frtCFPick(it,ich,rch,cfun,hfun)
implicit none
integer :: it
integer :: ich(:)
real*8 :: rch(:)
real*8 :: cfun, hfun
!
integer :: ic, j
real*8 :: v(NC), q(8), a0, a1, a2, Z
!
do ic=1,NC
do j=1,8
  q(j) = rch(6)*data(ic,it,ich(3),ich(4+j)) + &
      rch(7)*data(ic,it,ich(4),ich(4+j))
end do
v(ic) = exp( &
    rch(10)*(rch(9)*(rch(8)*q( 1)+rch(11)*q( 2))+   &
    rch(12)*(rch(8)*q( 3)+rch(11)*q( 4))) + &
    rch(13)*(rch(9)*(rch(8)*q( 5)+rch(11)*q( 6))+   &
    rch(12)*(rch(8)*q( 7)+rch(11)*q( 8))))

end do

a0 = v(1)
a1 = v(2)
a2 = v(3)
v(2) = 2*a1 - 0.5*a2 - 1.5*a0
v(3) = 0.5*(a0+a2) - a1

a0 = v(4)
a1 = v(5)
a2 = v(6)
v(5) = 2*a1 - 0.5*a2 - 1.5*a0
v(6) = 0.5*(a0+a2) - a1

Z = rch(5)

if(mode .eq. 1) then
cfun = (Z*v(3)+v(2))*Z
hfun = (Z*v(6)+v(5))*Z
else
cfun = (Z*v(3)+v(2))*Z + v(1)
hfun = (Z*v(6)+v(5))*Z + v(4)
end if

end subroutine frtCFPick

!***********************
!  Cache some table information into arrays iCache and rCache
subroutine frtCFCache(den,Z,Plw,Ph1,Pg1,Pc6,ich,rch,ierr)
implicit none
real*8, intent(in)   :: den, Z, Plw, Ph1, Pg1, Pc6
real*8, dimension(:) :: rch
integer, dimension(:) :: ich
integer :: ierr
!
real*8  :: q(3), qh1, qg1, qc6, dl, w
integer :: j, il(3), is(3)

!  Convert from nb to nH from Cloudy models
dl = max(1.0e-10,den*(1.0-0.02*Z)/1.4)
ierr = 0

if(Plw .gt. 0.0) then
qh1 = log10(1.0e-37+Ph1/Plw)
qg1 = log10(1.0e-37+Pg1/Plw)
qc6 = log10(1.0e-37+Pc6/Plw)

q(1) = log10(1.0e-37+Plw/dl)
q(2) = 0.263*qc6 + 0.353*qh1 + 0.923*qg1
q(3) = 0.976*qc6 - 0.103*qh1 - 0.375*qg1

!  qmin, qstp, etc are boundaries of cells, not their centers
do j=1,3
  w = 0.5 + (q(j)-qmin(j))/qstp(j)
  il(j) = int(w) + 1
  if(w .gt. il(j)-0.5) then
    is(j) = il(j) + 1
  else
    is(j) = il(j) - 1
  endif
  rch(10+j) = abs(il(j)-0.5-w)
  rch(7+j) = 1 - rch(10+j)

  if(np(j) .gt. 1) then
    if(max(il(j),is(j)) .gt. np(j)) ierr =  j
    if(min(il(j),is(j)) .lt.     1) ierr = -j
  endif

  if(il(j) .lt. 1) il(j) = 1
  if(is(j) .lt. 1) is(j) = 1
  if(il(j) .gt. np(j)) il(j) = np(j)
  if(is(j) .gt. np(j)) is(j) = np(j)
enddo

else

ierr = -1
do j=1,3
  il(j) = 1
  is(j) = 1
  rch(7+j) = 1
  rch(10+j) = 0
enddo

endif

!  Density interpolation is still CIC
w = (log10(dl)-xmin)/xstp
ich(3) = int(w) + 1
if(ich(3) .lt.  1) ich(3) = 1
if(ich(3) .ge. NX) ich(3) = NX-1
ich(4) = ich(3) + 1
rch(6) = max(0.0,min(1.0,ich(3)-w))
rch(7) = 1.0 - rch(6)

!  Do not forget C-to-F77 index conversion
ich(5) = 1 + indx(il(1),il(2),il(3))
ich(6) = 1 + indx(is(1),il(2),il(3))
ich(7) = 1 + indx(il(1),is(2),il(3))
ich(8) = 1 + indx(is(1),is(2),il(3))
ich(9) = 1 + indx(il(1),il(2),is(3))
ich(10) = 1 + indx(is(1),il(2),is(3))
ich(11) = 1 + indx(il(1),is(2),is(3))
ich(12) = 1 + indx(is(1),is(2),is(3))

rch(5) = Z

!  Clear temperature cache
ich(1) = 0
ich(2) = 0
end subroutine frtCFCache

! Get the cooling and heating functions for given
!  ln(T) from the cached dat
subroutine frtCFGetLn(alt,ich,rch,cfun,hfun)
real*8 :: alt
real*8,  dimension(:) :: rch
integer, dimension(:) :: ich
real*8 :: cfun, hfun
real*8, dimension(NC) :: v
integer      :: il, iu
real*8 :: ql, qu
!
il = int((alt-altmin)/altstp*0.99999) + 1
if(il .lt.  1) il = 1
if(il .ge. NT) il = NT-1
iu = il + 1
ql = max(0.0,min(1.0,(altval(iu)-alt)/altstp))
qu = 1.0 - ql
!
!  Shift cache lines as needed
!
if(ich(1) .eq. iu) then
ich(1) = 0
ich(2) = iu
rch(3) = rch(1)
rch(4) = rch(2)
endif
if(ich(2) .eq. il) then
ich(1) = il
ich(2) = 0
rch(1) = rch(3)
rch(2) = rch(4)
endif

!  Update the cache
if(ich(1) .ne. il) then
ich(1) = il
call frtCFPick(il,ich,rch,cfun,hfun)
rch(1) = cfun
rch(2) = hfun
endif
if(ich(2) .ne. iu) then
ich(2) = iu
call frtCFPick(iu,ich,rch,cfun,hfun)
rch(3) = cfun
rch(4) = hfun
endif

cfun = ql*rch(1) + qu*rch(3)
hfun = ql*rch(2) + qu*rch(4)
end subroutine frtCFGetLn

!  Get the cooling and heating functions for T
subroutine frtGetCF(tem,den,Z,Plw,Ph1,Pg1,Pc6,cfun,hfun,ierr)
implicit none
real*8 :: tem, den, Z, Plw, Ph1, Pg1, Pc6, cfun, hfun
integer :: ierr
! Cache arrays
real*8, dimension(NRCH) :: rch
integer,      dimension(NICH) :: ich

call frtCFCache(den,Z,Plw,Ph1,Pg1,Pc6,ich,rch,ierr)
call frtCFGetLn(log(max(1.0,tem)),ich,rch,cfun,hfun)

if (ierr > 0) then
!print *,'frtGetCF: Problems with caching Gnedin and Hollon cooling table.'
!print *,'frtGetCF:  WARNING. Error code :', ierr
!print *,'frtGetCF: Input variables ntot, Tgas, PLW, PHI, PHeI, PCVI: ', den, tem, Plw, Ph1, Pg1, Pc6
!stop
end if

end subroutine frtGetCF

end module KROME_coolingGH

!############### MODULE ##############
module krome_dust
contains

  !***********************
  subroutine init_dust_tabs()
    use krome_commons
    use krome_subs
    implicit none

    call init_anytab2D("dust_table_cool.dat",dust_tab_ngas(:),&
        dust_tab_Tgas(:), dust_tab_cool(:,:), dust_mult_ngas, &
        dust_mult_Tgas)
    call init_anytab2D("dust_table_Tdust.dat",dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_Tdust(:,:), dust_mult_ngas, &
        dust_mult_Tgas)
    call init_anytab2D("dust_table_H2.dat",dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_H2(:,:), dust_mult_ngas, &
        dust_mult_Tgas)

  end subroutine init_dust_tabs

end module krome_dust

!############### MODULE ##############
module KROME_cooling
  ! *************************************************************
  !  This file has been generated with:
  !  krome 14.08.dev "Beastie Boyle" on 2016-08-08 10:48:19.
  !
  ! KROME is a nice and friendly chemistry package for a wide range of
  !  astrophysical simulations. Given a chemical network (in CSV format)
  !  it automatically generates all the routines needed to solve the kinetic
  !  of the system, modelled as system of coupled Ordinary Differential
  !  Equations.
  !  It provides different options which make it unique and very flexible.
  !  Any suggestions and comments are welcomed. KROME is an open-source
  !  package, GNU-licensed, and any improvements provided by
  !  the users is well accepted. See disclaimer below and GNU License
  !  in gpl-3.0.txt.
  !
  !  more details in http://kromepackage.org/
  !  also see https://bitbucket.org/krome/krome_stable
  !
  ! Written and developed by Tommaso Grassi
  !  tommasograssi@gmail.com,
  !  Starplan Center, Copenhagen.
  !  Niels Bohr Institute, Copenhagen.
  !
  ! Co-developer Stefano Bovino
  !  sbovino@astro.physik.uni-goettingen.de
  !  Institut fuer Astrophysik, Goettingen.
  !
  ! Others (alphabetically): D. Galli, F.A. Gianturco, T. Haugboelle,
  !  J.Prieto, J.Ramsey, D.R.G. Schleicher, D. Seifried,
  !  E. Simoncini, E. Tognelli
  !
  !
  ! KROME is provided "as it is", without any warranty.
  !  The Authors assume no liability for any damages of any kind
  !  (direct or indirect damages, contractual or non-contractual
  !  damages, pecuniary or non-pecuniary damages), directly or
  !  indirectly derived or arising from the correct or incorrect
  !  usage of KROME, in any possible environment, or arising from
  !  the impossibility to use, fully or partially, the software,
  !  or any bug or malefunction.
  !  Such exclusion of liability expressly includes any damages
  !  including the loss of data of any kind (including personal data)
  ! *************************************************************
  integer,parameter::coolTab_n=int(1e2)
  integer,parameter::nZrate=0
  real*8::coolTab(nZrate,coolTab_n),coolTab_logTlow, coolTab_logTup
  real*8::coolTab_T(coolTab_n),inv_coolTab_T(coolTab_n-1),inv_coolTab_idx
contains

  !*******************
  function cooling(n,inTgas)
    use krome_commons
    implicit none
    real*8::n(:),inTgas,cooling,Tgas

    Tgas = inTgas
    cooling = sum(get_cooling_array(n(:),Tgas))

  end function cooling

  !*******************************
  function get_cooling_array(n, Tgas)
    use krome_commons
    implicit none
    real*8::n(:), Tgas
    real*8::get_cooling_array(ncools),cools(ncools)
    real*8::f1,f2,smooth

    !returns cooling in erg/cm3/s
    f1 = 1d0
    f2 = 1d0

    cools(:) = 0d0

    cools(idx_cool_H2) = cooling_H2(n(:), Tgas) - cooling_H2(n(:), phys_Tfloor)

    !this parameter controls the smoothness of the
    ! merge between the two cooling functions
    smooth = 2.d-3

    !smoothing functions | f1+f2=1
!    f1 = (tanh(smooth*(Tgas-1d4))+1.d0)*0.5d0
!    f2 = (tanh(smooth*(-Tgas+1d4))+1.d0)*0.5d0
    
    ! Only compute cooling terms if it really matters
!    if (f2 > 1d-20) then
      cools(idx_cool_atomic)  = f2 * cooling_Atomic(n(:), Tgas)
      cools(idx_cool_compton) = f2 * cooling_compton(n(:), Tgas)
      cools(idx_cool_ff)      = f2 * cooling_ff(n(:), Tgas)
      cools(idx_cool_ZCIENOUV)    = f2 * (cooling_Z_CIENOUV(n(:), Tgas) &
                                        - cooling_Z_CIENOUV(n(:), phys_Tfloor))
!    endif

!    if (f1 > 1d-20)then
!      cools(idx_cool_GH) = f1 * cooling_GH
!    endif

    !cools(idx_cool_GH) = cooling_GH

    cools(idx_cool_H2) = cooling_H2(n(:), Tgas)
    cools(idx_cool_custom) = cooling_custom(n(:),Tgas)

    get_cooling_array(:) = cools(:)

    end function get_cooling_array

    !***************************
    !Metal line cooling CIE
    ! method: CLOUDY 10, NOUV
    ! tables kindly provided by Sijing Shen.
    function cooling_Z_CIENOUV(n,inTgas)
      use krome_commons
      use krome_subs
      use krome_user_commons
      implicit none
      real*8::cooling_Z_CIENOUV,n(:),inTgas
      real*8::cH,Tgas,xLd,logcH

      cooling_Z_CIENOUV = 0d0
      cH = get_Hnuclei(n(:))

      !check if the abundance is close to zero to
      !avoid weird log evaluation
      if(cH.lt.1d-20)return

      Tgas = log10(inTgas)
      logcH = log10(cH)

      xLd = fit_anytab2D(CoolZNOUV_x(:), CoolZNOUV_y(:), CoolZNOUV_z(:,:), &
           CoolZNOUV_xmul, CoolZNOUV_ymul,logcH,Tgas)

      cooling_Z_CIENOUV = 10**xLd * cH * cH * total_Z

    end function cooling_Z_CIENOUV


  !***************************
  !Metal line cooling CIE
  ! method: CLOUDY 10, including the
  ! extragalactic (quasars + galaxies) UV flux
  ! by Haardt&Madau 2012.
  !Tables kindly provided by Sijing Shen.
  function cooling_Z_CIE(n,inTgas)
    use krome_commons
    use krome_subs
    implicit none
    integer,parameter::imax=coolZCIEn1
    integer,parameter::jmax=coolZCIEn2
    integer,parameter::kmax=coolZCIEn3
    integer::i,j,k
    real*8::cooling_Z_CIE,n(:),inTgas,Tgas
    real*8::v1,v2,v3,prev1,prev2,cH
    real*8::vv1,vv2,vv3,vv4,vv12,vv34,xLd
    real*8::x1(imax),x2(jmax),x3(kmax)
    real*8::ixd1(imax-1),ixd2(jmax-1),ixd3(kmax-1)
    real*8::v1min,v1max,v2min,v2max,v3min,v3max
    real*8,parameter::eps=1d-5

    Tgas = inTgas
    cooling_Z_CIE = 0d0

    !local copy of limits
    v1min = coolZCIEx1min
    v1max = coolZCIEx1max
    v2min = coolZCIEx2min
    v2max = coolZCIEx2max
    v3min = coolZCIEx3min
    v3max = coolZCIEx3max

    !local copy of variables arrays
    x1(:) = coolZCIEx1(:)
    x2(:) = coolZCIEx2(:)
    x3(:) = coolZCIEx3(:)

    ixd1(:) = coolZCIEixd1(:)
    ixd2(:) = coolZCIEixd2(:)
    ixd3(:) = coolZCIEixd3(:)

    !local variables
    cH = get_Hnuclei(n(:))

    !check if the abundance is close to zero to
    !avoid weird log evaluation
    if(cH.lt.1d-20)return

    v1 = Tgas           !Tgas
    v2 = cH             !total H number density
    v3 = phys_zredshift !redshift is linear

    !logs of variables
    v1 = log10(v1)
    v2 = log10(v2)

    !check limits
    if(v1>=v1max) v1 = v1max*(1d0-eps)
    if(v2>=v2max) v2 = v2max*(1d0-eps)
    if(v3>=v3max) v3 = v3max*(1d0-eps)

    if(v1<v1min) return
    if(v2<v2min) return
    if(v3<v3min) return

    !gets position of variable in the array
    i = (v1-v1min)*coolZCIEdvn1+1
    j = (v2-v2min)*coolZCIEdvn2+1
    k = (v3-v3min)*coolZCIEdvn3+1

    !precompute shared variables
    prev1 = (v1-x1(i))*ixd1(i)
    prev2 = (v2-x2(j))*ixd2(j)

    !linear interpolation on x1 for x2,x3
    vv1 = prev1 * (coolZCIEy(k,j,i+1) - &
        coolZCIEy(k,j,i)) + coolZCIEy(k,j,i)
    !linear interpolation on x1 for x2+dx2,x3
    vv2 = prev1 * (coolZCIEy(k,j+1,i+1) - &
        coolZCIEy(k,j+1,i)) + coolZCIEy(k,j+1,i)
    !linear interpolation on x2 for x3
    vv12 = prev2 * (vv2 - vv1) + vv1

    !linear interpolation on x1 for x2,x3+dx3
    vv3 = prev1 * (coolZCIEy(k+1,j,i+1) - &
        coolZCIEy(k+1,j,i)) + coolZCIEy(k+1,j,i)
    !linear interpolation on x1 for x2+dx2,x3+dx3
    vv4 = prev1 * (coolZCIEy(k+1,j+1,i+1) - &
        coolZCIEy(k+1,j+1,i)) + coolZCIEy(k+1,j+1,i)
    !linear interpolation on x2 for x3+dx3
    vv34 = prev2 * (vv4 - vv3) + vv3

    !linear interpolation on x3
    xLd = (v3-x3(k))*ixd3(k)*(vv34 - &
        vv12) + vv12

    !Z cooling in erg/s/cm3
    cooling_Z_CIE = 1d1**xLd * cH * cH * total_Z

  end function cooling_Z_CIE

  !************************
  subroutine init_coolingZCIE()
    use krome_commons
    implicit none
    integer::ios,iout(3),i
    real*8::rout(5)

    if(krome_mpi_rank<=1) then
      print *,"load Z_CIE2012 cooling..."
    endif

!    per test
!    print *,'process' , krome_mpi_rank

    open(krome_nfile,file=trim(krome_datafolder)//"coolZ_CIE2012.dat",status="old",iostat=ios)
    !check if file exists
    if(ios.ne.0) then
      print *,"ERROR: problems loading "//trim(krome_datafolder)//"coolZ_CIE2012.dat !"
      stop
    end if

    do
      read(krome_nfile,*,iostat=ios) iout(:),rout(:) !read line
      if(ios<0) exit !eof
      if(ios .eq. 5008) exit !gfortran eof
      if(ios/=0) cycle !skip blanks
      coolZCIEx1(iout(1)) = rout(1)
      coolZCIEx2(iout(2)) = rout(2)
      coolZCIEx3(iout(3)) = rout(3)
      coolZCIEy(iout(3),iout(2),iout(1)) = rout(4)
      heatZCIEy(iout(3),iout(2),iout(1)) = rout(5)
    end do

    close(krome_nfile)

    !store inverse of the differences
    ! to speed up interpolation
    do i=1,coolZCIEn1-1
      coolZCIEixd1(i) = 1d0/(coolZCIEx1(i+1)-coolZCIEx1(i))
    end do
    do i=1,coolZCIEn2-1
      coolZCIEixd2(i) = 1d0/(coolZCIEx2(i+1)-coolZCIEx2(i))
    end do
    do i=1,coolZCIEn3-1
      coolZCIEixd3(i) = 1d0/(coolZCIEx3(i+1)-coolZCIEx3(i))
    end do

    coolZCIEx1min = minval(coolZCIEx1)
    coolZCIEx1max = maxval(coolZCIEx1)
    coolZCIEx2min = minval(coolZCIEx2)
    coolZCIEx2max = maxval(coolZCIEx2)
    coolZCIEx3min = minval(coolZCIEx3)
    coolZCIEx3max = maxval(coolZCIEx3)

    coolZCIEdvn1 = (coolZCIEn1-1)/(coolZCIEx1max-coolZCIEx1min)
    coolZCIEdvn2 = (coolZCIEn2-1)/(coolZCIEx2max-coolZCIEx2min)
    coolZCIEdvn3 = (coolZCIEn3-1)/(coolZCIEx3max-coolZCIEx3min)

  end subroutine init_coolingZCIE

  !*****************************
  function cooling_custom(n,Tgas)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::n(:),Tgas,cooling_custom

    cooling_custom = 0d0

  end function cooling_custom

  !**********************************
  function kpla(n,Tgas)
    !Planck opacity mean fit (Lenzuni+1996)
    !only denisity dependent (note that the
    ! fit provided by Lenzuni is wrong)
    ! valid for T<3e3 K
    use krome_subs
    use krome_commons
    implicit none
    real*8::kpla,rhogas,Tgas,n(:),y
    real*8::a0,a1,m(nspec)

    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols)) !g/cm3

    kpla = 0.d0
    !opacity is zero under 1e-12 g/cm3
    if(rhogas<1d-12) return

    !fit coefficients
    a0 = 1.000042d0
    a1 = 2.14989d0

    !log density cannot exceed 0.5 g/cm3
    y = log10(min(rhogas,0.5d0))

    kpla = 1d1**(a0*y + a1) !fit density only

  end function kpla

  !*****************************
  function coolingChem(n,Tgas)
    implicit none
    real*8::coolingChem,n(:),Tgas

    !note that this function is a dummy.
    ! For chemical cooling you should see
    ! heatingChem function in krome_heating.f90

    coolingChem = 0.d0

  end function coolingChem

  !*******************************
  function cooling_compton(n, Tgas)
    !compton cooling erg/cm3/s from Cen1992
    use krome_user_commons
    use krome_commons
    real*8::cooling_compton,n(:),Tgas

    !note that redhsift is a common variable and
    ! should be provided by the user, otherwise the default is zero
    cooling_compton = 5.65d-36 * (1.d0 + phys_zredshift)**4 &
        * (Tgas - 2.73d0 * (1.d0 + phys_zredshift)) * n(idx_e) !erg/s/cm3

  end function cooling_compton

  !*****************
  !sigmoid function with x0 shift and s steepness
  function sigmoid(x,x0,s)
    implicit none
    real*8::sigmoid,x,x0,s

    sigmoid = 1d1/(1d1+exp(-s*(x-x0)))

  end function sigmoid

  !*******************
  !window function for H2 cooling to smooth limits
  function wCool(logTgas,logTmin,logTmax)
    implicit none
    real*8::wCool,logTgas,logTmin,logTmax,x

    x = (logTgas-logTmin)/(logTmax-logTmin)
    wCool = 1d1**(2d2*(sigmoid(x,-2d-1,5d1)*sigmoid(-x,-1.2d0,5d1)-1d0))
    if(wCool<1d-199) wCool = 0d0
    if(wCool>1d0) then
      print *,"ERROR: wCool>1"
      stop
    end if

  end function wCool

  !ALL THE COOLING FUNCTIONS ARE FROM GLOVER & ABEL, MNRAS 388, 1627, 2008
  !FOR LOW DENSITY REGIME: CONSIDER AN ORTHO-PARA RATIO OF 3:1
  !UPDATED TO THE DATA REPORTED BY GLOVER 2015, MNRAS
  !EACH SINGLE FUNCTION IS IN erg/s
  !FINAL UNITS = erg/cm3/s
  !*******************************
  function cooling_H2(n, Tgas)
    use krome_commons
    use krome_subs
    real*8::n(:),Tgas
    real*8::temp,logt3,logt,cool,cooling_H2,T3
    real*8::LDL,HDLR,HDLV,HDL
    real*8::logt32,logt33,logt34,logt35,logt36,logt37,logt38
    real*8::dump14,fH2H,fH2e,fH2H2,fH2Hp,fH2He,w14,w24
    integer::i
    character*16::names(nspec)

    temp = Tgas
    cooling_H2 = 0d0
    !if(temp<2d0) return

    T3 = temp * 1.d-3
    logt3 = log10(T3)
    logt = log10(temp)
    cool = 0d0

    logt32 = logt3 * logt3
    logt33 = logt32 * logt3
    logt34 = logt33 * logt3
    logt35 = logt34 * logt3
    logt36 = logt35 * logt3
    logt37 = logt36 * logt3
    logt38 = logt37 * logt3

    w14 = wCool(logt, 1d0, 4d0)
    w24 = wCool(logt, 2d0, 4d0)

    !//H2-H
    if(temp<=1d2) then
      fH2H = 1.d1**(-16.818342D0 +3.7383713D1*logt3 &
          +5.8145166D1*logt32 +4.8656103D1*logt33 &
          +2.0159831D1*logt34 +3.8479610D0*logt35 )*n(idx_H)
    elseif(temp>1d2 .and. temp<=1d3) then
      fH2H = 1.d1**(-2.4311209D1 +3.5692468D0*logt3 &
          -1.1332860D1*logt32 -2.7850082D1*logt33 &
          -2.1328264D1*logt34 -4.2519023D0*logt35)*n(idx_H)
    elseif(temp>1.d3.and.temp<=6d3) then
      fH2H = 1d1**(-2.4311209D1 +4.6450521D0*logt3 &
          -3.7209846D0*logt32 +5.9369081D0*logt33 &
          -5.5108049D0*logt34 +1.5538288D0*logt35)*n(idx_H)
    else
      fH2H = 1.862314467912518E-022*wCool(logt,1d0,log10(6d3))*n(idx_H)
    end if
    cool = cool + fH2H

    !//H2-Hp
    if(temp>1d1.and.temp<=1d4) then
      fH2Hp = 1d1**(-2.2089523d1 +1.5714711d0*logt3 &
          +0.015391166d0*logt32 -0.23619985d0*logt33 &
          -0.51002221d0*logt34 +0.32168730d0*logt35)*n(idx_Hj)
    else
      fH2Hp = 1.182509139382060E-021*n(idx_Hj)*w14
    endif
    cool = cool + fH2Hp

    !//H2-H2
    fH2H2 = w24*1d1**(-2.3962112D1 +2.09433740D0*logt3 &
        -.77151436D0*logt32 +.43693353D0*logt33 &
        -.14913216D0*logt34 -.033638326D0*logt35)*n(idx_H2) !&
        cool = cool + fH2H2

    !//H2-e
    fH2e = 0d0
    if(temp<=5d2) then
      fH2e = 1d1**(min(-2.1928796d1 + 1.6815730d1*logt3 &
          +9.6743155d1*logt32 +3.4319180d2*logt33 &
          +7.3471651d2*logt34 +9.8367576d2*logt35 &
          +8.0181247d2*logt36 +3.6414446d2*logt37 &
          +7.0609154d1*logt38,3d1))*n(idx_e)
    elseif(temp>5d2)  then
      fH2e = 1d1**(-2.2921189D1 +1.6802758D0*logt3 &
          +.93310622D0*logt32 +4.0406627d0*logt33 &
          -4.7274036d0*logt34 -8.8077017d0*logt35 &
          +8.9167183*logt36 + 6.4380698*logt37 &
          -6.3701156*logt38)*n(idx_e)
    end if
    cool = cool + fH2e*w24

    !//H2-He
    if(temp>1d1.and.temp<=1d4)then
      fH2He = 1d1**(-2.3689237d1 +2.1892372d0*logt3&
          -.81520438d0*logt32 +.29036281d0*logt33 -.16596184d0*logt34 &
          +.19191375d0*logt35)*n(idx_He)
    else
      fH2He = 1.002560385050777E-022*n(idx_He)*w14
    endif
    cool = cool + fH2He

    !check error
    if(cool>1.d30) then
      print *,"ERROR: cooling >1.d30 erg/s/cm3"
      print *,"cool (erg/s/cm3): ",cool
      names(:) = get_names()
      do i=1,size(n)
        print '(I3,a18,E11.3)',i,names(i),n(i)
      end do
      stop
    end if

    !this to avoid negative, overflow and useless calculations below
    if(cool<=0d0) then
      cooling_H2 = 0d0
      return
    end if

    !high density limit from HM79, GP98 below Tgas = 2d3
    !UPDATED TO THE DATA BY GLOVER 2015 for high temperature corrections, MNRAS
    !IN THE HIGH DENSITY REGIME LAMBDA_H2 = LAMBDA_H2(LTE) = HDL
    !the following mix of functions ensures the right behaviour
    ! at low temperatures (Tgas < 10 K) and at high temperatures (T > 2000 K) by
    ! using both the original Hollenbach and the new Glover data, merging them
    ! in a smooth way.
    if(temp.lt.2d3)then
      HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*exp(-(0.13/t3)**3)+&
          3.e-24*exp(-0.51/t3)) !erg/s
      HDLV = (6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3)) !erg/s
      HDL  = HDLR + HDLV !erg/s
    elseif(temp>=2d3 .and. temp<=1d4)then
      HDL = 1d1**(-2.0584225d1 + 5.0194035*logt3 &
          -1.5738805*logt32 -4.7155769*logt33 &
          +2.4714161*logt34 +5.4710750*logt35 &
          -3.9467356*logt36 -2.2148338*logt37 &
          +1.8161874*logt38)
    else
      dump14 = 1d0 / (1d0 + exp(min((temp-3d4)*2d-4,3d2)))
      HDL = 5.531333679406485E-019*dump14
    endif

    LDL = cool !erg/s
    if (HDL==0.) then
      cooling_H2 = 0.d0
    else
      cooling_H2 = n(idx_H2)/(1.d0/HDL+1.d0/LDL)  !erg/cm3/s
    endif

  end function cooling_H2

  !Atomic COOLING  Cen ApJS, 78, 341, 1992
  !UNITS = erg/s/cm3
  !*******************************
  function cooling_Atomic(n, Tgas)
    use krome_commons
    use krome_subs
    real*8::Tgas,cooling_atomic,n(:)
    real*8::temp,T5,cool

    temp = max(Tgas,10.d0) !K
    T5 = temp/1.d5 !K
    cool = 0.d0 !erg/cm3/s

    !COLLISIONAL IONIZATION: H, He, He+, He(2S)
    cool = cool+ 1.27d-21*sqrt(temp)/(1.d0+sqrt(T5))&
        *exp(-1.578091d5/temp)*n(idx_e)*n(idx_H)
    cool = cool+ 9.38d-22*sqrt(temp)/(1.d0+sqrt(T5))&
        *exp(-2.853354d5/temp)*n(idx_e)*n(idx_He)
    cool = cool+ 4.95d-22*sqrt(temp)/(1.d0+sqrt(T5))&
        *exp(-6.31515d5/temp)*n(idx_e)*n(idx_Hej)
    cool = cool+ 5.01d-27*temp**(-0.1687)/(1.d0+sqrt(T5))&
        *exp(-5.5338d4/temp)*n(idx_e)**2*n(idx_Hej)

    !RECOMBINATION: H+, He+,He2+
    cool = cool+ 8.7d-27*sqrt(temp)*(temp/1.d3)**(-0.2)&
        /(1.d0+(temp/1.d6)**0.7)*n(idx_e)*n(idx_Hj)
    cool = cool+ 1.55d-26*temp**(0.3647)*n(idx_e)*n(idx_Hej)
    cool = cool+ 3.48d-26*sqrt(temp)*(temp/1.d3)**(-0.2)&
        /(1.d0+(temp/1.d6)**0.7)*n(idx_e)*n(idx_Hejj)

    !DIELECTRONIC RECOMBINATION: He
    cool = cool+ 1.24d-13*temp**(-1.5)*exp(-4.7d5/temp)&
        *(1.d0+0.3d0*exp(-9.4d4/temp))*n(idx_e)*n(idx_Hej)

    !COLLISIONAL EXCITATION:
    !H(all n), He(n=2,3,4 triplets), He+(n=2)
    cool = cool+ 7.5d-19/(1.d0+sqrt(T5))*exp(-1.18348d5/temp)*n(idx_e)*n(idx_H)
    cool = cool+ 9.1d-27*temp**(-.1687)/(1.d0+sqrt(T5))&
        *exp(-1.3179d4/temp)*n(idx_e)**2*n(idx_Hej)
    cool = cool+ 5.54d-17*temp**(-.397)/(1.d0+sqrt(T5))&
        *exp(-4.73638d5/temp)*n(idx_e)*n(idx_Hej)

    cooling_atomic = max(cool, 0.d0)  !erg/cm3/s

  end function cooling_Atomic

  !**************************
  !free-free cooling (bremsstrahlung for all ions)
  ! using mean Gaunt factor value (Cen+1992)
  function cooling_ff(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(:),Tgas,cool,cooling_ff,gaunt_factor,bms_ions

    gaunt_factor = 1.5d0 !mean value

    !BREMSSTRAHLUNG: all ions
    bms_ions = +n(idx_Hj) +n(idx_HEj) +4.d0*n(idx_HEjj)
    cool = 1.42d-27*gaunt_factor*sqrt(Tgas)&
        *bms_ions*n(idx_e)

    cooling_ff = max(cool, 0.d0)  !erg/cm3/s

  end function cooling_ff

!************** GNEDIN HOLLON COOLING FUNCTION
  subroutine gnedin_hollon(n, Tgas)
  use krome_commons
  use krome_constants
  use krome_subs
  use krome_coolingGH
  implicit none
  real*8::n(:),Tgas
  real*8, dimension(nPhotoBins)::Nphot
  real*8:: PLW_last=0d0,PHI_last=0d0,PHeI_last=0d0, &
      PCVI_last=0d0,ntot_last=0d0
  real*8,  save :: rch(NRCH)
  integer, save :: ich(NICH)
  !$omp threadprivate(PLW_last,PHI_last,PHeI_last,PCVI_last,ntot_last,rch,ich)
  real*8::log10Tgas,ntot,cfun,hfun
  logical, save :: first_call=.true.
  integer       :: ierr

  Nphot(:) = photoBinJ(:) * 4. * pi / photoBinEmid(:)

  PLW  = H2_solomonLW_ramses(user_myH2_dissociation,user_myfluxLW)
  PHI  = DOT_PRODUCT(Nphot, HI_sigma)
  PHEI = DOT_PRODUCT(Nphot, HeI_sigma)
  PCVI = DOT_PRODUCT(Nphot, CVI_sigma)

  ! Check if we have to load table -- only done by one thread
  ! Do the double if-block so that we only enter the critical region
  ! if first_call has not happened or is in progress
  !----------------------------------
  if (first_call) then
  !$omp critical
  if (first_call) then
    ierr = 0
    !call frtInitCF(ierr,'cf_table.I2.dat')
    if(ierr .ne. 0) then
      print *,'Error in reading Gnedin and Hollon'
      print *,' cooling table data file: ', ierr
      print *,'Maybe data file cf_table.I2.dat'
      print *,' does not exist in current directory ?'
      stop
    endif
  end if
  first_call = .false.
  !$omp end critical
  end if

  !ntot = sum(n(1:nmols))
  ntot = 1.d0 * ( n(idx_HK) + n(idx_H) + n(idx_Hj) )       &
       + 4.d0 * ( n(idx_HE) + n(idx_HEj) + n(idx_HEjj) )   &
       + 2.d0 * ( n(idx_H2) + n(idx_H2j) )

  if (Tgas > 0d0) then       ! LIMIT WHERE GH IS RELIABLE: BEFORE IT WAS 1D3, I CHANGED TO 0
  !check input values to see if we have to regenerate the cache
  if (abs(ntot-ntot_last) > 1d-6*ntot .or. &
      PLW  .ne. PLW_last  .or. &
      PHI  .ne. PHI_last  .or. &
      PHeI .ne. PHeI_last .or. &
      PCVI .ne. PCVI_last) then
  call frtCFCache(ntot, 1d0, PLW, PHI, PHeI, PCVI, ich, rch, ierr)
  ntot_last=ntot
  PLW_last=PLW
  PHI_last=PHI
  PHeI_last=PHeI
  PCVI_last=PCVI
  if (ierr > 0) then
    !print *,'Problems with caching Gnedin and Hollon cooling table.'
    !print *,' WARNING. Error code :', ierr
    !print *,'Input variables ntot, Tgas, PLW, PHI, PHeI, PCVI: ', ntot, Tgas, PLW, PHI, PHeI, PCVI
    !stop
  end if
  end if

  call frtGetCF(Tgas, ntot, 1d0, Plw, Phi, Phei, Pcvi, cfun, hfun, ierr)

  else
  cfun = 0d0
  hfun = 0d0
  endif

  cooling_GH = SUM(n(1:nmols)) * ntot * cfun
  heating_GH = SUM(n(1:nmols)) * ntot * hfun

  end subroutine

  !***********************
  subroutine mylin2(a,b)
    !solve Ax=B analytically for a 2-levels system
    implicit none
    integer,parameter::n=2
    real*8::a(n,n),b(n),c(n),iab

    !uncomment this: safer but slower function
    !if(a(2,2)==a(2,1)) then
    !   print *,"ERROR: a22=a21 in mylin2"
    !   stop
    !end if
    iab = b(1)/(a(2,2)-a(2,1))
    c(1) = a(2,2) * iab
    c(2) = -a(2,1) * iab
    b(:) = c(:)

  end subroutine mylin2

  !************************
  subroutine mylin3(a,b)
    !solve Ax=B analytically for a 3-levels system
    implicit none
    integer,parameter::n=3
    real*8::iab,a(n,n),b(n),c(n)

    !uncomment this: safer but slower function
    !if(a(2,2)==a(2,3)) then
    !   print *,"ERROR: a22=a23 in mylin3"
    !   stop
    !end if

    !uncomment this: safer but slower
    !if(a(2,1)*a(3,2)+a(2,2)*a(3,3)+a(2,3)*a(3,1) == &
        !     a(2,1)*a(3,3)+a(2,2)*a(3,1)+a(2,3)*a(3,2)) then
    !   print *,"ERROR: division by zero in mylin3"
    !   stop
    !end if

    iab = b(1) / (a(2,1)*(a(3,3)-a(3,2)) + a(2,2)*(a(3,1)-a(3,3)) &
        + a(2,3)*(a(3,2)-a(3,1)))
    c(1) = (a(2,3)*a(3,2)-a(2,2)*a(3,3)) * iab
    c(2) = -(a(2,3)*a(3,1)-a(2,1)*a(3,3)) * iab
    c(3) = (a(3,1)*a(2,2)-a(2,1)*a(3,2)) * iab
    b(:) = c(:)

  end subroutine mylin3

  !************************************
  subroutine plot_cool(n)
    !routine to plot cooling at runtime
    use krome_commons
    implicit none
    real*8::n(:),Tgas,Tmin,Tmax
    real*8::cool_atomic,cool_H2,cool_HD,cool_tot, cool_totGP,cool_H2GP
    real*8::cool_dH,cool_Z
    integer::i,imax
    imax = 1000
    Tmin = log10(1d1)
    Tmax = log10(1d8)
    print *,"plotting cooling..."
    open(krome_nfile,file="KROME_cooling_plot.dat",status="replace")
    do i=1,imax
      Tgas = 1d1**(i*(Tmax-Tmin)/imax+Tmin)
      cool_H2 = 0.d0
      cool_H2GP = 0.d0
      cool_HD = 0.d0
      cool_atomic = 0.d0
      cool_Z = 0.d0
      cool_dH = 0.d0
      cool_H2 = cooling_H2(n(:),Tgas)
      cool_atomic = cooling_atomic(n(:),Tgas)
      cool_tot = cool_H2 + cool_atomic + cool_HD + cool_Z + cool_dH
      cool_totGP = cool_H2GP + cool_atomic + cool_HD + cool_Z + cool_dH
      write(krome_nfile,'(99E12.3e3)') Tgas, cool_tot, cool_totGP, cool_H2, &
          cool_atomic, cool_HD, cool_H2GP, cool_Z, cool_dH
    end do
    close(krome_nfile)
    print *,"done!"

  end subroutine plot_cool

!***********************************
  !routine to dump cooling in unit krome_nfile
  subroutine dump_cool(n,Tgas)
    use krome_commons
    implicit none
    real*8::Tgas,n(:),cools(ncools)

    cools(:) = get_cooling_array(n(:),Tgas)
    open(krome_nfile,file=trim(krome_datafolder)//"krome_dump_cool.dat",status="replace")
    write(krome_nfile,'(99E14.5e3)') Tgas, sum(cools), cools(:)
    close(krome_nfile)
  end subroutine dump_cool

end module KROME_cooling


!############### MODULE ##############
module KROME_heating
contains

  ! *************************************************************
  !  This file has been generated with:
  !  krome 14.08.dev "Beastie Boyle" on 2016-08-08 10:48:19.
  !
  ! KROME is a nice and friendly chemistry package for a wide range of
  !  astrophysical simulations. Given a chemical network (in CSV format)
  !  it automatically generates all the routines needed to solve the kinetic
  !  of the system, modelled as system of coupled Ordinary Differential
  !  Equations.
  !  It provides different options which make it unique and very flexible.
  !  Any suggestions and comments are welcomed. KROME is an open-source
  !  package, GNU-licensed, and any improvements provided by
  !  the users is well accepted. See disclaimer below and GNU License
  !  in gpl-3.0.txt.
  !
  !  more details in http://kromepackage.org/
  !  also see https://bitbucket.org/krome/krome_stable
  !
  ! Written and developed by Tommaso Grassi
  !  tommasograssi@gmail.com,
  !  Starplan Center, Copenhagen.
  !  Niels Bohr Institute, Copenhagen.
  !
  ! Co-developer Stefano Bovino
  !  sbovino@astro.physik.uni-goettingen.de
  !  Institut fuer Astrophysik, Goettingen.
  !
  ! Others (alphabetically): D. Galli, F.A. Gianturco, T. Haugboelle,
  !  J.Prieto, J.Ramsey, D.R.G. Schleicher, D. Seifried,
  !  E. Simoncini, E. Tognelli
  !
  !
  ! KROME is provided "as it is", without any warranty.
  !  The Authors assume no liability for any damages of any kind
  !  (direct or indirect damages, contractual or non-contractual
  !  damages, pecuniary or non-pecuniary damages), directly or
  !  indirectly derived or arising from the correct or incorrect
  !  usage of KROME, in any possible environment, or arising from
  !  the impossibility to use, fully or partially, the software,
  !  or any bug or malefunction.
  !  Such exclusion of liability expressly includes any damages
  !  including the loss of data of any kind (including personal data)
  ! *************************************************************

  !************************
  function heating(n,inTgas,k,nH2dust)
    implicit none
    real*8::n(:), Tgas, inTgas, k(:), nH2dust
    real*8::heating

    Tgas = inTgas
    heating = sum(get_heating_array(n(:),Tgas,k(:), nH2dust))

  end function heating

  !*******************************
  function get_heating_array(n, Tgas, k, nH2dust)
    use krome_commons
    implicit none
    real*8::n(:), Tgas, k(:), nH2dust
    real*8::get_heating_array(nheats),heats(nheats)
    real*8::smooth,f1,f2
    !returns heating in erg/cm3/s

    f1 = 1.d0
    f2 = 1.d0

!    f2 = (tanh(smooth*(-Tgas+1d4))+1.d0)*0.5d0

    smooth = 2.d3

    heats(:) = 0.d0

!    if(f2 > 1.d-20)then  ! T < 1e4 K
      heats(idx_heat_chem)  = f2 * heatingChem(n(:), Tgas, k(:), nH2dust)
      heats(idx_heat_photo) = f2 * photo_heating(n(:))
      !heats(idx_heat_ZCIE)  = f2 * heat_ZCIE(n(:),Tgas)
!    endif

!    if(f1 > 1.d-20)then  ! T > 1e4 K
!       heats(idx_heat_GH) = f1 * heating_GH
!    endif

    !heats(idx_heat_GH) = heating_GH

    heats(idx_heat_dust)  = heat_netPhotoDust(n(:),Tgas)
    heats(idx_heat_CR)    = heat_CR(n(:),Tgas,k(:))
    heats(idx_heat_custom) = heat_custom(n(:),Tgas)

    get_heating_array(:) = heats(:)

  end function get_heating_array

  !*************************
  function heat_custom(n,Tgas)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::n(:),Tgas,heat_custom

    heat_custom = 0d0

  end function heat_custom

  function heat_ZCIE(n,inTgas)
    use krome_commons
    use krome_subs
    implicit none
    integer,parameter::imax=coolZCIEn1
    integer,parameter::jmax=coolZCIEn2
    integer,parameter::kmax=coolZCIEn3
    integer::i,j,k
    real*8::heat_ZCIE,n(:),inTgas,Tgas
    real*8::v1,v2,v3,prev1,prev2,cH
    real*8::vv1_h,vv2_h,vv3_h,vv4_h,vv12_h,vv34_h,xGd
    real*8::x1(imax),x2(jmax),x3(kmax)
    real*8::ixd1(imax-1),ixd2(jmax-1),ixd3(kmax-1)
    real*8::v1min,v1max,v2min,v2max,v3min,v3max
    real*8,parameter::eps=1d-5

    Tgas = inTgas
    heat_ZCIE = 0d0

    !local copy of limits
    v1min = coolZCIEx1min
    v1max = coolZCIEx1max
    v2min = coolZCIEx2min
    v2max = coolZCIEx2max
    v3min = coolZCIEx3min
    v3max = coolZCIEx3max

    !local copy of variables arrays
    x1(:) = coolZCIEx1(:)
    x2(:) = coolZCIEx2(:)
    x3(:) = coolZCIEx3(:)

    ixd1(:) = coolZCIEixd1(:)
    ixd2(:) = coolZCIEixd2(:)
    ixd3(:) = coolZCIEixd3(:)

    !local variables
    cH = get_Hnuclei(n(:))

    !check if the abundance is close to zero to
    !avoid weird log evaluation
    if(cH.lt.1d-20)return

    v1 = Tgas           !Tgas
    v2 = cH             !total H number density
    v3 = phys_zredshift !redshift is linear

    !logs of variables
    v1 = log10(v1)
    v2 = log10(v2)

    !check limits
    if(v1>=v1max) v1 = v1max*(1d0-eps)
    if(v2>=v2max) v2 = v2max*(1d0-eps)
    if(v3>=v3max) v3 = v3max*(1d0-eps)

    if(v1<v1min) return
    if(v2<v2min) return
    if(v3<v3min) return

    !gets position of variable in the array
    i = (v1-v1min)*coolZCIEdvn1+1
    j = (v2-v2min)*coolZCIEdvn2+1
    k = (v3-v3min)*coolZCIEdvn3+1

    !precompute shared variables
    prev1 = (v1-x1(i))*ixd1(i)
    prev2 = (v2-x2(j))*ixd2(j)

    !linear interpolation on x1 for x2,x3
    vv1_h = prev1 * (heatZCIEy(k,j,i+1) - &
        heatZCIEy(k,j,i)) + heatZCIEy(k,j,i)
    !linear interpolation on x1 for x2+dx2,x3
    vv2_h = prev1 * (heatZCIEy(k,j+1,i+1) - &
        heatZCIEy(k,j+1,i)) + heatZCIEy(k,j+1,i)
    !linear interpolation on x2 for x3
    vv12_h = prev2 * (vv2_h - vv1_h) + vv1_h

    !linear interpolation on x1 for x2,x3+dx3
    vv3_h = prev1 * (heatZCIEy(k+1,j,i+1) - &
        heatZCIEy(k+1,j,i)) + heatZCIEy(k+1,j,i)
    !linear interpolation on x1 for x2+dx2,x3+dx3
    vv4_h = prev1 * (heatZCIEy(k+1,j+1,i+1) - &
        heatZCIEy(k+1,j+1,i)) + heatZCIEy(k+1,j+1,i)
    !linear interpolation on x2 for x3+dx3
    vv34_h = prev2 * (vv4_h - vv3_h) + vv3_h

    !linear interpolation on x3
    xGd = (v3-x3(k))*ixd3(k)*(vv34_h - &
        vv12_h) + vv12_h

    !Z cooling in erg/s/cm3
    heat_ZCIE = 1d1**xGd * cH * cH * total_Z

  end function heat_ZCIE

  !***************************
  function heat_netPhotoDust(n,Tgas)
    !photoelectric effect from dust in erg/s/cm3
    !including the recombination cooling and a generic radiation flux
    !eq. 42 and 44 in Bakes&Tielens, 1994
    ! dust2gas_ratio is D/D_sol, default assumes D/D_sol = Z/Z_sol
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    integer::i
    real*8::heat_netPhotoDust,n(:),Tgas,ntot,eps
    real*8::psi,recomb_cool,bet

    ntot = get_Hnuclei(n(:))
    bet = 0.735d0*(Tgas)**(-0.068)

    if(n(idx_e)>0d0) then
      psi = GHabing * sqrt(Tgas) / n(idx_e)
    else
      psi = 0d0
    end if

    !grains recombination cooling
    recomb_cool = 4.65d-30*Tgas**0.94*psi**bet &
        * n(idx_e)*n(idx_H)

    eps = 4.9d-2 / (1d0 + 4d-3 * psi**.73) + &
        3.7d-2 * (Tgas * 1d-4)**.7 / (1d0 + 2d-4 * psi)

    !net photoelectric heating
    heat_netPhotoDust = (1.3d-24*eps*GHabing*ntot-recomb_cool)*dust2gas_ratio

  end function heat_netPhotoDust

  !***************************
  function heat_CR(n,Tgas,k)
    !heating from cosmic rays erg/s/cm3
    use krome_commons
    implicit none
    real*8::heat_CR,n(:),Tgas,Hfact,k(:)

    Hfact = 3.20435313d-11 !erg

    heat_CR = 0d0

    !H -> H+ + E
    heat_CR = heat_CR + k(47) * n(idx_H) * Hfact

    !HE -> HE+ + E
    heat_CR = heat_CR + k(48) * n(idx_HE) * Hfact

    !H2 -> H + H
    heat_CR = heat_CR + k(49) * n(idx_H2) * Hfact

    !H2 -> H+ + H-
    heat_CR = heat_CR + k(50) * n(idx_H2) * Hfact

    !H2 -> H2+ + E
    heat_CR = heat_CR + k(51) * n(idx_H2) * Hfact

    !H2 -> H + H+ + E
    heat_CR = heat_CR + k(52) * n(idx_H2) * Hfact

  end function heat_CR

  !**************************
  function photo_heating(n)
    !photo heating in erg/cm3/s using bin-based
    ! approach. Terms are computed in the
    ! krome_photo module
    use krome_commons
    use krome_constants
    implicit none
    real*8::photo_heating,n(:)

    photo_heating = 0.d0
    photo_heating = photoBinHeats(1) * n(idx_H) &
        + photoBinHeats(2) * n(idx_HE) &
        + photoBinHeats(3) * n(idx_HEj) &
        + photoBinHeats(4) * n(idx_Hk) &
        + photoBinHeats(5) * n(idx_H2) &
        + photoBinHeats(6) * n(idx_H2j) &
        + photoBinHeats(7) * n(idx_H2j) &
        + photoBinHeats(8) * n(idx_H2)

  end function photo_heating

  !H2 FORMATION HEATING and other exo/endothermic
  ! processes (including H2 on dust) in erg/cm3/s
  !krome builds the heating/cooling term according
  ! to the chemical network employed
  !*******************************
  function heatingChem(n, Tgas, k, nH2dust)
    use krome_constants
    use krome_commons
    use krome_dust
    use krome_subs
    implicit none
    real*8::heatingChem, n(:), Tgas,k(:),nH2dust
    real*8::h2heatfac,HChem,yH,yH2, denom
    real*8::ncr,ncrn,ncrd1,ncrd2,dd,n2H,small,nmax
    dd = get_Hnuclei(n(:))

    !replace small according to the desired enviroment
    ! and remove nmax if needed
    nmax = maxval(n(1:nmols))
    small = 0d0

    heatingChem = 0.d0

    ncrn  = 1.0d6*(Tgas**(-0.5d0))
    ncrd1 = 1.6d0*exp(-(4.0d2/Tgas)**2)
    ncrd2 = 1.4d0*exp(-1.2d4/(Tgas+1.2d3))

    yH = n(idx_H)/dd   !dimensionless
    yH2= n(idx_H2)/dd  !dimensionless

    denom = max((ncrd1*yH+ncrd2*yH2), 1e-10)
    ncr = ncrn/denom      !1/cm3
    h2heatfac = 1.0d0/(1.0d0+ncr/dd)     !dimensionless

    HChem = 0.d0 !inits chemical heating
    n2H = n(idx_H) * n(idx_H)

    !H- + H -> H2 + E (heating)
    HChem = HChem + k(10) * (3.53d0*h2heatfac*n(idx_Hk)*n(idx_H))
    !H2+ + H -> H2 + H+ (heating)
    HChem = HChem + k(13) * (1.83d0*h2heatfac*n(idx_H2j)*n(idx_H))
    !H2 + E -> H + H + E (cooling)
    HChem = HChem + k(17) * (-4.48d0*n(idx_H2)*n(idx_E))
    !H2 + H -> H + H + H (cooling)
    HChem = HChem + k(18) * (-4.48d0*n(idx_H2)*n(idx_H))
    !H2 + H2 -> H + H + H2 (cooling)
    HChem = HChem + k(27) * (-4.48d0*n(idx_H2)*n(idx_H2))
    HChem = HChem + k(44) * (18.7d0*h2heatfac + 0.4d0)*n(idx_H2)

    HChem = HChem + nH2dust * (4.2d0*h2heatfac + 0.2d0)

    heatingChem = HChem * eV_to_erg  !erg/cm3/s

  end function heatingChem

end module KROME_heating

!############### MODULE ##############
module krome_stars

end module krome_stars

!############### MODULE ##############
module krome_ode
contains

  ! *************************************************************
  !  This file has been generated with:
  !  krome 14.08.dev "Beastie Boyle" on 2016-08-08 10:48:19.
  !
  ! KROME is a nice and friendly chemistry package for a wide range of
  !  astrophysical simulations. Given a chemical network (in CSV format)
  !  it automatically generates all the routines needed to solve the kinetic
  !  of the system, modelled as system of coupled Ordinary Differential
  !  Equations.
  !  It provides different options which make it unique and very flexible.
  !  Any suggestions and comments are welcomed. KROME is an open-source
  !  package, GNU-licensed, and any improvements provided by
  !  the users is well accepted. See disclaimer below and GNU License
  !  in gpl-3.0.txt.
  !
  !  more details in http://kromepackage.org/
  !  also see https://bitbucket.org/krome/krome_stable
  !
  ! Written and developed by Tommaso Grassi
  !  tommasograssi@gmail.com,
  !  Starplan Center, Copenhagen.
  !  Niels Bohr Institute, Copenhagen.
  !
  ! Co-developer Stefano Bovino
  !  sbovino@astro.physik.uni-goettingen.de
  !  Institut fuer Astrophysik, Goettingen.
  !
  ! Others (alphabetically): D. Galli, F.A. Gianturco, T. Haugboelle,
  !  J.Prieto, J.Ramsey, D.R.G. Schleicher, D. Seifried,
  !  E. Simoncini, E. Tognelli
  !
  !
  ! KROME is provided "as it is", without any warranty.
  !  The Authors assume no liability for any damages of any kind
  !  (direct or indirect damages, contractual or non-contractual
  !  damages, pecuniary or non-pecuniary damages), directly or
  !  indirectly derived or arising from the correct or incorrect
  !  usage of KROME, in any possible environment, or arising from
  !  the impossibility to use, fully or partially, the software,
  !  or any bug or malefunction.
  !  Such exclusion of liability expressly includes any damages
  !  including the loss of data of any kind (including personal data)
  ! *************************************************************

  subroutine fex(neq,tt,nin,dn)
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_cooling
    use krome_heating
    use krome_tabs
    use krome_photo
    implicit none
    integer::neq,idust
    real*8::tt,dn(neq),n(neq),k(nrea),krome_gamma
    real*8::gamma,Tgas,vgas,ntot,nH2dust,nd,nin(neq)
    real*8::rr
    integer::i,r1,r2,p1,p2,p3

    n(:) = nin(:)

    nH2dust = 0.d0
    n(idx_CR) = 1.d0
    n(idx_g)  = 1.d0
    n(idx_dummy) = 1.d0

    dn(:) = 0.d0 !initialize differentials
    n(idx_Tgas) = max(n(idx_tgas),phys_Tfloor)
    n(idx_Tgas) = min(n(idx_tgas),1d9)
    Tgas = n(idx_Tgas) !get temperature

    call calcHabingThick(n(:),Tgas)

    dust2gas_ratio = dust2gas_ratio_ref * fevap(Tgas)

    k(:) = coe_tab(n(:)) !compute coefficients

    nH2dust = nH2dust + H2_dustJura(n(:),Tgas)

    !H
    dn(idx_H) = &
        -k(1)*n(idx_H)*n(idx_E) &
        +k(2)*n(idx_Hj)*n(idx_E) &
        +k(3)*n(idx_Hj)*n(idx_E) &
        -k(9)*n(idx_H)*n(idx_E) &
        -k(10)*n(idx_Hk)*n(idx_H) &
        -k(11)*n(idx_H)*n(idx_Hj) &
        -k(12)*n(idx_H)*n(idx_Hj) &
        -k(13)*n(idx_H2j)*n(idx_H) &
        +k(14)*n(idx_H2)*n(idx_Hj) &
        +k(15)*n(idx_H2)*n(idx_Hj) &
        +k(16)*n(idx_H2)*n(idx_E) &
        +2.d0*k(17)*n(idx_H2)*n(idx_E) &
        -k(18)*n(idx_H2)*n(idx_H) &
        +3.d0*k(18)*n(idx_H2)*n(idx_H) &
        +k(19)*n(idx_Hk)*n(idx_E) &
        -k(20)*n(idx_Hk)*n(idx_H) &
        +2.d0*k(20)*n(idx_Hk)*n(idx_H) &
        -k(21)*n(idx_Hk)*n(idx_H) &
        +2.d0*k(21)*n(idx_Hk)*n(idx_H) &
        +2.d0*k(22)*n(idx_Hk)*n(idx_Hj) &
        +2.d0*k(24)*n(idx_H2j)*n(idx_E) &
        +2.d0*k(25)*n(idx_H2j)*n(idx_E) &
        +k(26)*n(idx_H2j)*n(idx_Hk) &
        +2.d0*k(27)*n(idx_H2)*n(idx_H2) &
        -k(28)*n(idx_HEj)*n(idx_H) &
        +k(29)*n(idx_HE)*n(idx_Hj) &
        +k(30)*n(idx_HE)*n(idx_Hj) &
        +k(31)*n(idx_H2)*n(idx_HEj) &
        +2.d0*k(32)*n(idx_H2)*n(idx_HE) &
        -2.d0*k(34)*n(idx_H)*n(idx_H) &
        +k(34)*n(idx_H)*n(idx_H) &
        -k(35)*n(idx_H)*n(idx_HE) &
        -k(36)*n(idx_H) &
        +k(39)*n(idx_Hk) &
        +k(41)*n(idx_H2j) &
        +2.d0*k(43)*n(idx_H2) &
        +2.d0*k(44)*n(idx_H2) &
        +k(45)*n(idx_Hj)*n(idx_E) &
        -k(47)*n(idx_H) &
        +2.d0*k(49)*n(idx_H2) &
        +k(52)*n(idx_H2) - 2d0*nH2dust

   ! if (user_myfluxLW .le. 6.d-59 .and. k(39) .ge.1.d-20) then
      !print*, -k(1), +k(2), +k(3), -k(9), -k(10), -k(11), -k(12), -k(13), +k(14), +k(15), +k(16), k(17), +k(18), &
      !       k(19), k(20), k(21), k(22), k(24), k(25), k(26), k(27), -k(28), k(29), k(30), k(31), k(32), -k(34), &
      !       -k(35), -k(36), k(39), k(41), k(43), k(44), k(45), -k(47), k(49), k(52)
 !     print*, "user_myfluxLW",  user_myfluxLW
    !  print*, "rates", k(39), k(41), photobinj(4), photobinj(6)
    !endif

    !E
    dn(idx_E) = &
        -k(1)*n(idx_H)*n(idx_E) &
        +2.d0*k(1)*n(idx_H)*n(idx_E) &
        -k(2)*n(idx_Hj)*n(idx_E) &
        -k(3)*n(idx_Hj)*n(idx_E) &
        -k(4)*n(idx_HE)*n(idx_E) &
        +2.d0*k(4)*n(idx_HE)*n(idx_E) &
        -k(5)*n(idx_HEj)*n(idx_E) &
        -k(6)*n(idx_HEj)*n(idx_E) &
        -k(7)*n(idx_HEj)*n(idx_E) &
        +2.d0*k(7)*n(idx_HEj)*n(idx_E) &
        -k(8)*n(idx_HEjj)*n(idx_E) &
        -k(9)*n(idx_H)*n(idx_E) &
        +k(10)*n(idx_Hk)*n(idx_H) &
        -k(16)*n(idx_H2)*n(idx_E) &
        -k(17)*n(idx_H2)*n(idx_E) &
        +k(17)*n(idx_H2)*n(idx_E) &
        -k(19)*n(idx_Hk)*n(idx_E) &
        +2.d0*k(19)*n(idx_Hk)*n(idx_E) &
        +k(20)*n(idx_Hk)*n(idx_H) &
        +k(21)*n(idx_Hk)*n(idx_H) &
        +k(23)*n(idx_Hk)*n(idx_Hj) &
        -k(24)*n(idx_H2j)*n(idx_E) &
        -k(25)*n(idx_H2j)*n(idx_E) &
        +k(34)*n(idx_H)*n(idx_H) &
        +k(35)*n(idx_H)*n(idx_HE) &
        +k(36)*n(idx_H) &
        +k(37)*n(idx_HE) &
        +k(38)*n(idx_HEj) &
        +k(39)*n(idx_Hk) &
        +k(40)*n(idx_H2) &
        +k(42)*n(idx_H2j) &
        -k(45)*n(idx_Hj)*n(idx_E) &
        -k(46)*n(idx_HEj)*n(idx_E) &
        +k(47)*n(idx_H) &
        +k(48)*n(idx_HE) &
        +k(51)*n(idx_H2) &
        +k(52)*n(idx_H2)

    !H+
    dn(idx_Hj) = &
        +k(1)*n(idx_H)*n(idx_E) &
        -k(2)*n(idx_Hj)*n(idx_E) &
        -k(3)*n(idx_Hj)*n(idx_E) &
        -k(11)*n(idx_H)*n(idx_Hj) &
        -k(12)*n(idx_H)*n(idx_Hj) &
        +k(13)*n(idx_H2j)*n(idx_H) &
        -k(14)*n(idx_H2)*n(idx_Hj) &
        -k(15)*n(idx_H2)*n(idx_Hj) &
        -k(22)*n(idx_Hk)*n(idx_Hj) &
        -k(23)*n(idx_Hk)*n(idx_Hj) &
        +k(28)*n(idx_HEj)*n(idx_H) &
        -k(29)*n(idx_HE)*n(idx_Hj) &
        -k(30)*n(idx_HE)*n(idx_Hj) &
        +k(31)*n(idx_H2)*n(idx_HEj) &
        +k(34)*n(idx_H)*n(idx_H) &
        +k(35)*n(idx_H)*n(idx_HE) &
        +k(36)*n(idx_H) &
        +k(41)*n(idx_H2j) &
        +2.d0*k(42)*n(idx_H2j) &
        -k(45)*n(idx_Hj)*n(idx_E) &
        +k(47)*n(idx_H) &
        +k(50)*n(idx_H2) &
        +k(52)*n(idx_H2)

    !HE
    dn(idx_HE) = &
        -k(4)*n(idx_HE)*n(idx_E) &
        +k(5)*n(idx_HEj)*n(idx_E) &
        +k(6)*n(idx_HEj)*n(idx_E) &
        +k(28)*n(idx_HEj)*n(idx_H) &
        -k(29)*n(idx_HE)*n(idx_Hj) &
        -k(30)*n(idx_HE)*n(idx_Hj) &
        +k(31)*n(idx_H2)*n(idx_HEj) &
        -k(32)*n(idx_H2)*n(idx_HE) &
        +k(32)*n(idx_H2)*n(idx_HE) &
        +k(33)*n(idx_H2)*n(idx_HEj) &
        -k(35)*n(idx_H)*n(idx_HE) &
        +k(35)*n(idx_H)*n(idx_HE) &
        -k(37)*n(idx_HE) &
        +k(46)*n(idx_HEj)*n(idx_E) &
        -k(48)*n(idx_HE)

    !HE+
    dn(idx_HEj) = &
        +k(4)*n(idx_HE)*n(idx_E) &
        -k(5)*n(idx_HEj)*n(idx_E) &
        -k(6)*n(idx_HEj)*n(idx_E) &
        -k(7)*n(idx_HEj)*n(idx_E) &
        +k(8)*n(idx_HEjj)*n(idx_E) &
        -k(28)*n(idx_HEj)*n(idx_H) &
        +k(29)*n(idx_HE)*n(idx_Hj) &
        +k(30)*n(idx_HE)*n(idx_Hj) &
        -k(31)*n(idx_H2)*n(idx_HEj) &
        -k(33)*n(idx_H2)*n(idx_HEj) &
        +k(37)*n(idx_HE) &
        -k(38)*n(idx_HEj) &
        -k(46)*n(idx_HEj)*n(idx_E) &
        +k(48)*n(idx_HE)

    !HE++
    dn(idx_HEjj) = &
        +k(7)*n(idx_HEj)*n(idx_E) &
        -k(8)*n(idx_HEjj)*n(idx_E) &
        +k(38)*n(idx_HEj)

    !H-
    dn(idx_Hk) = &
        +k(9)*n(idx_H)*n(idx_E) &
        -k(10)*n(idx_Hk)*n(idx_H) &
        +k(16)*n(idx_H2)*n(idx_E) &
        -k(19)*n(idx_Hk)*n(idx_E) &
        -k(20)*n(idx_Hk)*n(idx_H) &
        -k(21)*n(idx_Hk)*n(idx_H) &
        -k(22)*n(idx_Hk)*n(idx_Hj) &
        -k(23)*n(idx_Hk)*n(idx_Hj) &
        -k(26)*n(idx_H2j)*n(idx_Hk) &
        -k(39)*n(idx_Hk) &
        +k(50)*n(idx_H2)

    !H2
    dn(idx_H2) = &
        +k(10)*n(idx_Hk)*n(idx_H) &
        +k(13)*n(idx_H2j)*n(idx_H) &
        -k(14)*n(idx_H2)*n(idx_Hj) &
        -k(15)*n(idx_H2)*n(idx_Hj) &
        -k(16)*n(idx_H2)*n(idx_E) &
        -k(17)*n(idx_H2)*n(idx_E) &
        -k(18)*n(idx_H2)*n(idx_H) &
        +k(26)*n(idx_H2j)*n(idx_Hk) &
        -2.d0*k(27)*n(idx_H2)*n(idx_H2) &
        +k(27)*n(idx_H2)*n(idx_H2) &
        -k(31)*n(idx_H2)*n(idx_HEj) &
        -k(32)*n(idx_H2)*n(idx_HE) &
        -k(33)*n(idx_H2)*n(idx_HEj) &
        -k(40)*n(idx_H2) &
        -k(43)*n(idx_H2) &
        -k(44)*n(idx_H2) &
        -k(49)*n(idx_H2) &
        -k(50)*n(idx_H2) &
        -k(51)*n(idx_H2) &
        -k(52)*n(idx_H2) + nH2dust

!    if (user_myfluxLW .le. 6.d-59) then
!     print*, 'H2', user_myfluxLW, dn(idx_H2), &
!        +k(10), &
!        +k(13), &
!        -k(14), &
!        -k(15), &
!        -k(16), &
!        -k(17), &
!        -k(18), &
!        +k(26), &
!        -2.d0*k(27), &
!        +k(27), &
!        -k(31), &
!        -k(32), &
!        -k(33), &
!        -k(40), &
!        -k(43), &
!        -k(44), &
!        -k(49), &
!        -k(50), &
!        -k(51), &
!        -k(52),
!    endif

    !H2+
    dn(idx_H2j) = &
        +k(11)*n(idx_H)*n(idx_Hj) &
        +k(12)*n(idx_H)*n(idx_Hj) &
        -k(13)*n(idx_H2j)*n(idx_H) &
        +k(14)*n(idx_H2)*n(idx_Hj) &
        +k(15)*n(idx_H2)*n(idx_Hj) &
        +k(23)*n(idx_Hk)*n(idx_Hj) &
        -k(24)*n(idx_H2j)*n(idx_E) &
        -k(25)*n(idx_H2j)*n(idx_E) &
        -k(26)*n(idx_H2j)*n(idx_Hk) &
        +k(33)*n(idx_H2)*n(idx_HEj) &
        +k(40)*n(idx_H2) &
        -k(41)*n(idx_H2j) &
        -k(42)*n(idx_H2j) &
        +k(51)*n(idx_H2)

    !CR
    dn(idx_CR) = &
        0.d0

    !g
    dn(idx_g) = 0.d0

    !Tgas
    dn(idx_Tgas) = 0.d0

    !dummy
    dn(idx_dummy) = 0.d0

    krome_gamma = gamma_index(n(:))

   ! call gnedin_hollon(n(:), Tgas)

    dn(idx_Tgas) = (heating(n(:), Tgas, k(:), nH2dust) &
        - cooling(n(:), Tgas)  ) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))

    last_coe(:) = k(:)

  end subroutine fex

  !***************************
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use krome_commons
    use krome_subs
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    implicit none
    integer::neq, j, ian, jan, r1, r2, p1, p2, p3, i
    real*8::tt, n(neq), pdj(neq), dr1, dr2, kk,k(nrea),Tgas
    real*8::nn(neq),dn0,dn1,dnn,nH2dust,dn(neq),krome_gamma

    nH2dust = 0.d0
    Tgas = n(idx_Tgas)

    krome_gamma = gamma_index(n(:))

    k(:) = last_coe(:) !get rate coefficients

    if(j==1) then
      pdj(1) =  &
          +3.d0*k(18)*n(idx_H2)  &
          -k(9)*n(idx_E)  &
          -k(13)*n(idx_H2j)  &
          -k(10)*n(idx_Hk)  &
          -k(12)*n(idx_Hj)  &
          -k(21)*n(idx_Hk)  &
          -k(1)*n(idx_E)  &
          -k(35)*n(idx_HE)  &
          -4.d0*k(34)*n(idx_H)  &
          -k(28)*n(idx_HEj)  &
          -k(11)*n(idx_Hj)  &
          -k(18)*n(idx_H2)  &
          -k(20)*n(idx_Hk)  &
          +2.d0*k(20)*n(idx_Hk)  &
          -k(36)  &
          +2.d0*k(34)*n(idx_H)  &
          -k(47)  &
          +2.d0*k(21)*n(idx_Hk)
      pdj(2) =  &
          +k(36)  &
          +k(47)  &
          +k(10)*n(idx_Hk)  &
          +2.d0*k(34)*n(idx_H)  &
          -k(1)*n(idx_E)  &
          -k(9)*n(idx_E)  &
          +2.d0*k(1)*n(idx_E)  &
          +k(20)*n(idx_Hk)  &
          +k(21)*n(idx_Hk)  &
          +k(35)*n(idx_HE)
      pdj(3) =  &
          +k(36)  &
          -k(11)*n(idx_Hj)  &
          -k(12)*n(idx_Hj)  &
          +2.d0*k(34)*n(idx_H)  &
          +k(28)*n(idx_HEj)  &
          +k(13)*n(idx_H2j)  &
          +k(1)*n(idx_E)  &
          +k(47)  &
          +k(35)*n(idx_HE)
      pdj(4) =  &
          +k(28)*n(idx_HEj)  &
          +k(35)*n(idx_HE)  &
          -k(35)*n(idx_HE)
      pdj(5) =  &
          -k(28)*n(idx_HEj)
      pdj(7) =  &
          +k(9)*n(idx_E)  &
          -k(21)*n(idx_Hk)  &
          -k(20)*n(idx_Hk)  &
          -k(10)*n(idx_Hk)
      pdj(8) =  &
          +k(10)*n(idx_Hk)  &
          -k(18)*n(idx_H2)  &
          +k(13)*n(idx_H2j)
      pdj(9) =  &
          -k(13)*n(idx_H2j)  &
          +k(11)*n(idx_Hj)  &
          +k(12)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(1)*1d-3
      if(dnn>0.d0) then
        nn(1) = n(1) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==2) then
      pdj(1) =  &
          -k(1)*n(idx_H)  &
          +k(45)*n(idx_Hj)  &
          +k(3)*n(idx_Hj)  &
          +2.d0*k(17)*n(idx_H2)  &
          -k(9)*n(idx_H)  &
          +2.d0*k(25)*n(idx_H2j)  &
          +k(16)*n(idx_H2)  &
          +k(19)*n(idx_Hk)  &
          +k(2)*n(idx_Hj)  &
          +2.d0*k(24)*n(idx_H2j)
      pdj(2) =  &
          -k(1)*n(idx_H)  &
          -k(25)*n(idx_H2j)  &
          -k(3)*n(idx_Hj)  &
          -k(7)*n(idx_HEj)  &
          -k(6)*n(idx_HEj)  &
          -k(4)*n(idx_HE)  &
          -k(9)*n(idx_H)  &
          -k(46)*n(idx_HEj)  &
          +2.d0*k(1)*n(idx_H)  &
          -k(5)*n(idx_HEj)  &
          -k(24)*n(idx_H2j)  &
          -k(45)*n(idx_Hj)  &
          +2.d0*k(7)*n(idx_HEj)  &
          +k(17)*n(idx_H2)  &
          -k(16)*n(idx_H2)  &
          -k(17)*n(idx_H2)  &
          -k(2)*n(idx_Hj)  &
          -k(19)*n(idx_Hk)  &
          +2.d0*k(4)*n(idx_HE)  &
          -k(8)*n(idx_HEjj)  &
          +2.d0*k(19)*n(idx_Hk)
      pdj(3) =  &
          -k(2)*n(idx_Hj)  &
          -k(3)*n(idx_Hj)  &
          -k(45)*n(idx_Hj)  &
          +k(1)*n(idx_H)
      pdj(4) =  &
          -k(4)*n(idx_HE)  &
          +k(6)*n(idx_HEj)  &
          +k(5)*n(idx_HEj)  &
          +k(46)*n(idx_HEj)
      pdj(5) =  &
          -k(7)*n(idx_HEj)  &
          -k(6)*n(idx_HEj)  &
          -k(5)*n(idx_HEj)  &
          -k(46)*n(idx_HEj)  &
          +k(8)*n(idx_HEjj)  &
          +k(4)*n(idx_HE)
      pdj(6) =  &
          +k(7)*n(idx_HEj)  &
          -k(8)*n(idx_HEjj)
      pdj(7) =  &
          -k(19)*n(idx_Hk)  &
          +k(16)*n(idx_H2)  &
          +k(9)*n(idx_H)
      pdj(8) =  &
          -k(16)*n(idx_H2)  &
          -k(17)*n(idx_H2)
      pdj(9) =  &
          -k(25)*n(idx_H2j)  &
          -k(24)*n(idx_H2j)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(2)*1d-3
      if(dnn>0.d0) then
        nn(2) = n(2) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==3) then
      pdj(1) =  &
          +2.d0*k(22)*n(idx_Hk)  &
          +k(29)*n(idx_HE)  &
          +k(45)*n(idx_E)  &
          +k(3)*n(idx_E)  &
          -k(11)*n(idx_H)  &
          +k(30)*n(idx_HE)  &
          -k(12)*n(idx_H)  &
          +k(15)*n(idx_H2)  &
          +k(2)*n(idx_E)  &
          +k(14)*n(idx_H2)
      pdj(2) =  &
          +k(23)*n(idx_Hk)  &
          -k(45)*n(idx_E)  &
          -k(3)*n(idx_E)  &
          -k(2)*n(idx_E)
      pdj(3) =  &
          -k(2)*n(idx_E)  &
          -k(3)*n(idx_E)  &
          -k(29)*n(idx_HE)  &
          -k(11)*n(idx_H)  &
          -k(30)*n(idx_HE)  &
          -k(15)*n(idx_H2)  &
          -k(14)*n(idx_H2)  &
          -k(12)*n(idx_H)  &
          -k(23)*n(idx_Hk)  &
          -k(45)*n(idx_E)  &
          -k(22)*n(idx_Hk)
      pdj(4) =  &
          -k(30)*n(idx_HE)  &
          -k(29)*n(idx_HE)
      pdj(5) =  &
          +k(29)*n(idx_HE)  &
          +k(30)*n(idx_HE)
      pdj(7) =  &
          -k(23)*n(idx_Hk)  &
          -k(22)*n(idx_Hk)
      pdj(8) =  &
          -k(15)*n(idx_H2)  &
          -k(14)*n(idx_H2)
      pdj(9) =  &
          +k(23)*n(idx_Hk)  &
          +k(11)*n(idx_H)  &
          +k(15)*n(idx_H2)  &
          +k(14)*n(idx_H2)  &
          +k(12)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(3)*1d-3
      if(dnn>0.d0) then
        nn(3) = n(3) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==4) then
      pdj(1) =  &
          +k(29)*n(idx_Hj)  &
          +2.d0*k(32)*n(idx_H2)  &
          -k(35)*n(idx_H)  &
          +k(30)*n(idx_Hj)
      pdj(2) =  &
          +k(37)  &
          +k(48)  &
          +k(35)*n(idx_H)  &
          +2.d0*k(4)*n(idx_E)  &
          -k(4)*n(idx_E)
      pdj(3) =  &
          +k(35)*n(idx_H)  &
          -k(30)*n(idx_Hj)  &
          -k(29)*n(idx_Hj)
      pdj(4) =  &
          +k(32)*n(idx_H2)  &
          +k(35)*n(idx_H)  &
          -k(48)  &
          -k(29)*n(idx_Hj)  &
          -k(32)*n(idx_H2)  &
          -k(37)  &
          -k(35)*n(idx_H)  &
          -k(4)*n(idx_E)  &
          -k(30)*n(idx_Hj)
      pdj(5) =  &
          +k(29)*n(idx_Hj)  &
          +k(37)  &
          +k(48)  &
          +k(4)*n(idx_E)  &
          +k(30)*n(idx_Hj)
      pdj(8) =  &
          -k(32)*n(idx_H2)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(4)*1d-3
      if(dnn>0.d0) then
        nn(4) = n(4) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==5) then
      pdj(1) =  &
          -k(28)*n(idx_H)  &
          +k(31)*n(idx_H2)
      pdj(2) =  &
          -k(6)*n(idx_E)  &
          +2.d0*k(7)*n(idx_E)  &
          +k(38)  &
          -k(7)*n(idx_E)  &
          -k(46)*n(idx_E)  &
          -k(5)*n(idx_E)
      pdj(3) =  &
          +k(28)*n(idx_H)  &
          +k(31)*n(idx_H2)
      pdj(4) =  &
          +k(31)*n(idx_H2)  &
          +k(46)*n(idx_E)  &
          +k(6)*n(idx_E)  &
          +k(5)*n(idx_E)  &
          +k(33)*n(idx_H2)  &
          +k(28)*n(idx_H)
      pdj(5) =  &
          -k(33)*n(idx_H2)  &
          -k(31)*n(idx_H2)  &
          -k(6)*n(idx_E)  &
          -k(28)*n(idx_H)  &
          -k(5)*n(idx_E)  &
          -k(7)*n(idx_E)  &
          -k(46)*n(idx_E)  &
          -k(38)
      pdj(6) =  &
          +k(7)*n(idx_E)  &
          +k(38)
      pdj(8) =  &
          -k(33)*n(idx_H2)  &
          -k(31)*n(idx_H2)
      pdj(9) =  &
          +k(33)*n(idx_H2)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(5)*1d-3
      if(dnn>0.d0) then
        nn(5) = n(5) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==6) then
      pdj(2) =  &
          -k(8)*n(idx_E)
      pdj(5) =  &
          +k(8)*n(idx_E)
      pdj(6) =  &
          -k(8)*n(idx_E)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(6)*1d-3
      if(dnn>0.d0) then
        nn(6) = n(6) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==7) then
      pdj(1) =  &
          +k(26)*n(idx_H2j)  &
          -k(20)*n(idx_H)  &
          +2.d0*k(22)*n(idx_Hj)  &
          +2.d0*k(20)*n(idx_H)  &
          +k(39)  &
          -k(10)*n(idx_H)  &
          -k(21)*n(idx_H)  &
          +2.d0*k(21)*n(idx_H)  &
          +k(19)*n(idx_E)
      pdj(2) =  &
          +k(10)*n(idx_H)  &
          +k(23)*n(idx_Hj)  &
          +k(20)*n(idx_H)  &
          +k(39)  &
          -k(19)*n(idx_E)  &
          +k(21)*n(idx_H)  &
          +2.d0*k(19)*n(idx_E)
      pdj(3) =  &
          -k(23)*n(idx_Hj)  &
          -k(22)*n(idx_Hj)
      pdj(7) =  &
          -k(20)*n(idx_H)  &
          -k(22)*n(idx_Hj)  &
          -k(26)*n(idx_H2j)  &
          -k(10)*n(idx_H)  &
          -k(19)*n(idx_E)  &
          -k(21)*n(idx_H)  &
          -k(23)*n(idx_Hj)  &
          -k(39)
      pdj(8) =  &
          +k(26)*n(idx_H2j)  &
          +k(10)*n(idx_H)
      pdj(9) =  &
          -k(26)*n(idx_H2j)  &
          +k(23)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(7)*1d-3
      if(dnn>0.d0) then
        nn(7) = n(7) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==8) then
      pdj(1) =  &
          +4.d0*k(27)*n(idx_H2)  &
          +2.d0*k(32)*n(idx_HE)  &
          +k(14)*n(idx_Hj)  &
          +2.d0*k(43)  &
          -k(18)*n(idx_H)  &
          +2.d0*k(17)*n(idx_E)  &
          +2.d0*k(44)  &
          +k(15)*n(idx_Hj)  &
          +k(16)*n(idx_E)  &
          +k(31)*n(idx_HEj)  &
          +3.d0*k(18)*n(idx_H)  &
          +k(52)  &
          +2.d0*k(49)
      pdj(2) =  &
          +k(52)  &
          +k(17)*n(idx_E)  &
          +k(51)  &
          -k(17)*n(idx_E)  &
          -k(16)*n(idx_E)  &
          +k(40)
      pdj(3) =  &
          +k(50)  &
          +k(52)  &
          -k(15)*n(idx_Hj)  &
          +k(31)*n(idx_HEj)  &
          -k(14)*n(idx_Hj)
      pdj(4) =  &
          +k(33)*n(idx_HEj)  &
          -k(32)*n(idx_HE)  &
          +k(32)*n(idx_HE)  &
          +k(31)*n(idx_HEj)
      pdj(5) =  &
          -k(31)*n(idx_HEj)  &
          -k(33)*n(idx_HEj)
      pdj(7) =  &
          +k(50)  &
          +k(16)*n(idx_E)
      pdj(8) =  &
          +2.d0*k(27)*n(idx_H2)  &
          -k(50)  &
          -4.d0*k(27)*n(idx_H2)  &
          -k(43)  &
          -k(31)*n(idx_HEj)  &
          -k(51)  &
          -k(18)*n(idx_H)  &
          -k(52)  &
          -k(40)  &
          -k(32)*n(idx_HE)  &
          -k(44)  &
          -k(49)  &
          -k(17)*n(idx_E)  &
          -k(16)*n(idx_E)  &
          -k(33)*n(idx_HEj)  &
          -k(15)*n(idx_Hj)  &
          -k(14)*n(idx_Hj)
      pdj(9) =  &
          +k(40)  &
          +k(15)*n(idx_Hj)  &
          +k(51)  &
          +k(33)*n(idx_HEj)  &
          +k(14)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(8)*1d-3
      if(dnn>0.d0) then
        nn(8) = n(8) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==9) then
      pdj(1) =  &
          +2.d0*k(24)*n(idx_E)  &
          +2.d0*k(25)*n(idx_E)  &
          -k(13)*n(idx_H)  &
          +k(41)  &
          +k(26)*n(idx_Hk)
      pdj(2) =  &
          -k(24)*n(idx_E)  &
          -k(25)*n(idx_E)  &
          +k(42)
      pdj(3) =  &
          +k(13)*n(idx_H)  &
          +k(41)  &
          +2.d0*k(42)
      pdj(7) =  &
          -k(26)*n(idx_Hk)
      pdj(8) =  &
          +k(13)*n(idx_H)  &
          +k(26)*n(idx_Hk)
      pdj(9) =  &
          -k(24)*n(idx_E)  &
          -k(26)*n(idx_Hk)  &
          -k(41)  &
          -k(25)*n(idx_E)  &
          -k(13)*n(idx_H)  &
          -k(42)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(9)*1d-3
      if(dnn>0.d0) then
        nn(9) = n(9) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==10) then
      pdj(12) = 0.d0
    elseif(j==11) then
      pdj(12) = 0.d0
    elseif(j==12) then
      !use fex to compute temperature-dependent Jacobian
      dnn = n(idx_Tgas)*1d-3
      nn(:) = n(:)
      nn(idx_Tgas) = n(idx_Tgas) + dnn
      call fex(neq,tt,nn(:),dn(:))
      do i=1,neq-1
        pdj(i) = dn(i) / dnn
      end do
    elseif(j==13) then
      pdj(12) = 0.d0
    end if

    return
  end subroutine jes

  !*************************
  subroutine jex(neq,t,n,ml,mu,pd,npd)
    use krome_commons
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_subs
    implicit none
    real*8::n(neq),pd(neq,neq),t,k(nrea),dn0,dn1,dnn,Tgas
    real*8::krome_gamma,nn(neq),nH2dust
    integer::neq,ml,mu,npd

    Tgas = n(idx_Tgas)
    npd = neq
    k(:) = coe_tab(n(:))
    pd(:,:) = 0d0
    krome_gamma = gamma_index(n(:))

    !d[H_dot]/d[H]
    pd(1,1) =  &
        +3.d0*k(18)*n(idx_H2)  &
        -k(9)*n(idx_E)  &
        -k(13)*n(idx_H2j)  &
        -k(10)*n(idx_Hk)  &
        -k(12)*n(idx_Hj)  &
        -k(21)*n(idx_Hk)  &
        -k(1)*n(idx_E)  &
        -k(35)*n(idx_HE)  &
        -4.d0*k(34)*n(idx_H)  &
        -k(28)*n(idx_HEj)  &
        -k(11)*n(idx_Hj)  &
        -k(18)*n(idx_H2)  &
        -k(20)*n(idx_Hk)  &
        +2.d0*k(20)*n(idx_Hk)  &
        -k(36)  &
        +2.d0*k(34)*n(idx_H)  &
        -k(47)  &
        +2.d0*k(21)*n(idx_Hk)

    !d[E_dot]/d[H]
    pd(2,1) =  &
        +k(36)  &
        +k(47)  &
        +k(10)*n(idx_Hk)  &
        +2.d0*k(34)*n(idx_H)  &
        -k(1)*n(idx_E)  &
        -k(9)*n(idx_E)  &
        +2.d0*k(1)*n(idx_E)  &
        +k(20)*n(idx_Hk)  &
        +k(21)*n(idx_Hk)  &
        +k(35)*n(idx_HE)

    !d[H+_dot]/d[H]
    pd(3,1) =  &
        +k(36)  &
        -k(11)*n(idx_Hj)  &
        -k(12)*n(idx_Hj)  &
        +2.d0*k(34)*n(idx_H)  &
        +k(28)*n(idx_HEj)  &
        +k(13)*n(idx_H2j)  &
        +k(1)*n(idx_E)  &
        +k(47)  &
        +k(35)*n(idx_HE)

    !d[HE_dot]/d[H]
    pd(4,1) =  &
        +k(28)*n(idx_HEj)  &
        +k(35)*n(idx_HE)  &
        -k(35)*n(idx_HE)

    !d[HE+_dot]/d[H]
    pd(5,1) =  &
        -k(28)*n(idx_HEj)

    !d[H-_dot]/d[H]
    pd(7,1) =  &
        +k(9)*n(idx_E)  &
        -k(21)*n(idx_Hk)  &
        -k(20)*n(idx_Hk)  &
        -k(10)*n(idx_Hk)

    !d[H2_dot]/d[H]
    pd(8,1) =  &
        +k(10)*n(idx_Hk)  &
        -k(18)*n(idx_H2)  &
        +k(13)*n(idx_H2j)

    !d[H2+_dot]/d[H]
    pd(9,1) =  &
        -k(13)*n(idx_H2j)  &
        +k(11)*n(idx_Hj)  &
        +k(12)*n(idx_Hj)

    !d[Tgas_dot]/d[H]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(1)*1d-3
    if(dnn>0.d0) then
      nn(1) = n(1) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,1) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[E]
    pd(1,2) =  &
        -k(1)*n(idx_H)  &
        +k(45)*n(idx_Hj)  &
        +k(3)*n(idx_Hj)  &
        +2.d0*k(17)*n(idx_H2)  &
        -k(9)*n(idx_H)  &
        +2.d0*k(25)*n(idx_H2j)  &
        +k(16)*n(idx_H2)  &
        +k(19)*n(idx_Hk)  &
        +k(2)*n(idx_Hj)  &
        +2.d0*k(24)*n(idx_H2j)

    !d[E_dot]/d[E]
    pd(2,2) =  &
        -k(1)*n(idx_H)  &
        -k(25)*n(idx_H2j)  &
        -k(3)*n(idx_Hj)  &
        -k(7)*n(idx_HEj)  &
        -k(6)*n(idx_HEj)  &
        -k(4)*n(idx_HE)  &
        -k(9)*n(idx_H)  &
        -k(46)*n(idx_HEj)  &
        +2.d0*k(1)*n(idx_H)  &
        -k(5)*n(idx_HEj)  &
        -k(24)*n(idx_H2j)  &
        -k(45)*n(idx_Hj)  &
        +2.d0*k(7)*n(idx_HEj)  &
        +k(17)*n(idx_H2)  &
        -k(16)*n(idx_H2)  &
        -k(17)*n(idx_H2)  &
        -k(2)*n(idx_Hj)  &
        -k(19)*n(idx_Hk)  &
        +2.d0*k(4)*n(idx_HE)  &
        -k(8)*n(idx_HEjj)  &
        +2.d0*k(19)*n(idx_Hk)

    !d[H+_dot]/d[E]
    pd(3,2) =  &
        -k(2)*n(idx_Hj)  &
        -k(3)*n(idx_Hj)  &
        -k(45)*n(idx_Hj)  &
        +k(1)*n(idx_H)

    !d[HE_dot]/d[E]
    pd(4,2) =  &
        -k(4)*n(idx_HE)  &
        +k(6)*n(idx_HEj)  &
        +k(5)*n(idx_HEj)  &
        +k(46)*n(idx_HEj)

    !d[HE+_dot]/d[E]
    pd(5,2) =  &
        -k(7)*n(idx_HEj)  &
        -k(6)*n(idx_HEj)  &
        -k(5)*n(idx_HEj)  &
        -k(46)*n(idx_HEj)  &
        +k(8)*n(idx_HEjj)  &
        +k(4)*n(idx_HE)

    !d[HE++_dot]/d[E]
    pd(6,2) =  &
        +k(7)*n(idx_HEj)  &
        -k(8)*n(idx_HEjj)

    !d[H-_dot]/d[E]
    pd(7,2) =  &
        -k(19)*n(idx_Hk)  &
        +k(16)*n(idx_H2)  &
        +k(9)*n(idx_H)

    !d[H2_dot]/d[E]
    pd(8,2) =  &
        -k(16)*n(idx_H2)  &
        -k(17)*n(idx_H2)

    !d[H2+_dot]/d[E]
    pd(9,2) =  &
        -k(25)*n(idx_H2j)  &
        -k(24)*n(idx_H2j)

    !d[Tgas_dot]/d[E]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(2)*1d-3
    if(dnn>0.d0) then
      nn(2) = n(2) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,2) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[H+]
    pd(1,3) =  &
        +2.d0*k(22)*n(idx_Hk)  &
        +k(29)*n(idx_HE)  &
        +k(45)*n(idx_E)  &
        +k(3)*n(idx_E)  &
        -k(11)*n(idx_H)  &
        +k(30)*n(idx_HE)  &
        -k(12)*n(idx_H)  &
        +k(15)*n(idx_H2)  &
        +k(2)*n(idx_E)  &
        +k(14)*n(idx_H2)

    !d[E_dot]/d[H+]
    pd(2,3) =  &
        +k(23)*n(idx_Hk)  &
        -k(45)*n(idx_E)  &
        -k(3)*n(idx_E)  &
        -k(2)*n(idx_E)

    !d[H+_dot]/d[H+]
    pd(3,3) =  &
        -k(2)*n(idx_E)  &
        -k(3)*n(idx_E)  &
        -k(29)*n(idx_HE)  &
        -k(11)*n(idx_H)  &
        -k(30)*n(idx_HE)  &
        -k(15)*n(idx_H2)  &
        -k(14)*n(idx_H2)  &
        -k(12)*n(idx_H)  &
        -k(23)*n(idx_Hk)  &
        -k(45)*n(idx_E)  &
        -k(22)*n(idx_Hk)

    !d[HE_dot]/d[H+]
    pd(4,3) =  &
        -k(30)*n(idx_HE)  &
        -k(29)*n(idx_HE)

    !d[HE+_dot]/d[H+]
    pd(5,3) =  &
        +k(29)*n(idx_HE)  &
        +k(30)*n(idx_HE)

    !d[H-_dot]/d[H+]
    pd(7,3) =  &
        -k(23)*n(idx_Hk)  &
        -k(22)*n(idx_Hk)

    !d[H2_dot]/d[H+]
    pd(8,3) =  &
        -k(15)*n(idx_H2)  &
        -k(14)*n(idx_H2)

    !d[H2+_dot]/d[H+]
    pd(9,3) =  &
        +k(23)*n(idx_Hk)  &
        +k(11)*n(idx_H)  &
        +k(15)*n(idx_H2)  &
        +k(14)*n(idx_H2)  &
        +k(12)*n(idx_H)

    !d[Tgas_dot]/d[H+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(3)*1d-3
    if(dnn>0.d0) then
      nn(3) = n(3) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,3) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[HE]
    pd(1,4) =  &
        +k(29)*n(idx_Hj)  &
        +2.d0*k(32)*n(idx_H2)  &
        -k(35)*n(idx_H)  &
        +k(30)*n(idx_Hj)

    !d[E_dot]/d[HE]
    pd(2,4) =  &
        +k(37)  &
        +k(48)  &
        +k(35)*n(idx_H)  &
        +2.d0*k(4)*n(idx_E)  &
        -k(4)*n(idx_E)

    !d[H+_dot]/d[HE]
    pd(3,4) =  &
        +k(35)*n(idx_H)  &
        -k(30)*n(idx_Hj)  &
        -k(29)*n(idx_Hj)

    !d[HE_dot]/d[HE]
    pd(4,4) =  &
        +k(32)*n(idx_H2)  &
        +k(35)*n(idx_H)  &
        -k(48)  &
        -k(29)*n(idx_Hj)  &
        -k(32)*n(idx_H2)  &
        -k(37)  &
        -k(35)*n(idx_H)  &
        -k(4)*n(idx_E)  &
        -k(30)*n(idx_Hj)

    !d[HE+_dot]/d[HE]
    pd(5,4) =  &
        +k(29)*n(idx_Hj)  &
        +k(37)  &
        +k(48)  &
        +k(4)*n(idx_E)  &
        +k(30)*n(idx_Hj)

    !d[H2_dot]/d[HE]
    pd(8,4) =  &
        -k(32)*n(idx_H2)

    !d[Tgas_dot]/d[HE]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(4)*1d-3
    if(dnn>0.d0) then
      nn(4) = n(4) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,4) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[HE+]
    pd(1,5) =  &
        -k(28)*n(idx_H)  &
        +k(31)*n(idx_H2)

    !d[E_dot]/d[HE+]
    pd(2,5) =  &
        -k(6)*n(idx_E)  &
        +2.d0*k(7)*n(idx_E)  &
        +k(38)  &
        -k(7)*n(idx_E)  &
        -k(46)*n(idx_E)  &
        -k(5)*n(idx_E)

    !d[H+_dot]/d[HE+]
    pd(3,5) =  &
        +k(28)*n(idx_H)  &
        +k(31)*n(idx_H2)

    !d[HE_dot]/d[HE+]
    pd(4,5) =  &
        +k(31)*n(idx_H2)  &
        +k(46)*n(idx_E)  &
        +k(6)*n(idx_E)  &
        +k(5)*n(idx_E)  &
        +k(33)*n(idx_H2)  &
        +k(28)*n(idx_H)

    !d[HE+_dot]/d[HE+]
    pd(5,5) =  &
        -k(33)*n(idx_H2)  &
        -k(31)*n(idx_H2)  &
        -k(6)*n(idx_E)  &
        -k(28)*n(idx_H)  &
        -k(5)*n(idx_E)  &
        -k(7)*n(idx_E)  &
        -k(46)*n(idx_E)  &
        -k(38)

    !d[HE++_dot]/d[HE+]
    pd(6,5) =  &
        +k(7)*n(idx_E)  &
        +k(38)

    !d[H2_dot]/d[HE+]
    pd(8,5) =  &
        -k(33)*n(idx_H2)  &
        -k(31)*n(idx_H2)

    !d[H2+_dot]/d[HE+]
    pd(9,5) =  &
        +k(33)*n(idx_H2)

    !d[Tgas_dot]/d[HE+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(5)*1d-3
    if(dnn>0.d0) then
      nn(5) = n(5) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,5) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[HE++]
    pd(2,6) =  &
        -k(8)*n(idx_E)

    !d[HE+_dot]/d[HE++]
    pd(5,6) =  &
        +k(8)*n(idx_E)

    !d[HE++_dot]/d[HE++]
    pd(6,6) =  &
        -k(8)*n(idx_E)

    !d[Tgas_dot]/d[HE++]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(6)*1d-3
    if(dnn>0.d0) then
      nn(6) = n(6) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,6) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[H-]
    pd(1,7) =  &
        +k(26)*n(idx_H2j)  &
        -k(20)*n(idx_H)  &
        +2.d0*k(22)*n(idx_Hj)  &
        +2.d0*k(20)*n(idx_H)  &
        +k(39)  &
        -k(10)*n(idx_H)  &
        -k(21)*n(idx_H)  &
        +2.d0*k(21)*n(idx_H)  &
        +k(19)*n(idx_E)

    !d[E_dot]/d[H-]
    pd(2,7) =  &
        +k(10)*n(idx_H)  &
        +k(23)*n(idx_Hj)  &
        +k(20)*n(idx_H)  &
        +k(39)  &
        -k(19)*n(idx_E)  &
        +k(21)*n(idx_H)  &
        +2.d0*k(19)*n(idx_E)

    !d[H+_dot]/d[H-]
    pd(3,7) =  &
        -k(23)*n(idx_Hj)  &
        -k(22)*n(idx_Hj)

    !d[H-_dot]/d[H-]
    pd(7,7) =  &
        -k(20)*n(idx_H)  &
        -k(22)*n(idx_Hj)  &
        -k(26)*n(idx_H2j)  &
        -k(10)*n(idx_H)  &
        -k(19)*n(idx_E)  &
        -k(21)*n(idx_H)  &
        -k(23)*n(idx_Hj)  &
        -k(39)

    !d[H2_dot]/d[H-]
    pd(8,7) =  &
        +k(26)*n(idx_H2j)  &
        +k(10)*n(idx_H)

    !d[H2+_dot]/d[H-]
    pd(9,7) =  &
        -k(26)*n(idx_H2j)  &
        +k(23)*n(idx_Hj)

    !d[Tgas_dot]/d[H-]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(7)*1d-3
    if(dnn>0.d0) then
      nn(7) = n(7) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,7) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[H2]
    pd(1,8) =  &
        +4.d0*k(27)*n(idx_H2)  &
        +2.d0*k(32)*n(idx_HE)  &
        +k(14)*n(idx_Hj)  &
        +2.d0*k(43)  &
        -k(18)*n(idx_H)  &
        +2.d0*k(17)*n(idx_E)  &
        +2.d0*k(44)  &
        +k(15)*n(idx_Hj)  &
        +k(16)*n(idx_E)  &
        +k(31)*n(idx_HEj)  &
        +3.d0*k(18)*n(idx_H)  &
        +k(52)  &
        +2.d0*k(49)

    !d[E_dot]/d[H2]
    pd(2,8) =  &
        +k(52)  &
        +k(17)*n(idx_E)  &
        +k(51)  &
        -k(17)*n(idx_E)  &
        -k(16)*n(idx_E)  &
        +k(40)

    !d[H+_dot]/d[H2]
    pd(3,8) =  &
        +k(50)  &
        +k(52)  &
        -k(15)*n(idx_Hj)  &
        +k(31)*n(idx_HEj)  &
        -k(14)*n(idx_Hj)

    !d[HE_dot]/d[H2]
    pd(4,8) =  &
        +k(33)*n(idx_HEj)  &
        -k(32)*n(idx_HE)  &
        +k(32)*n(idx_HE)  &
        +k(31)*n(idx_HEj)

    !d[HE+_dot]/d[H2]
    pd(5,8) =  &
        -k(31)*n(idx_HEj)  &
        -k(33)*n(idx_HEj)

    !d[H-_dot]/d[H2]
    pd(7,8) =  &
        +k(50)  &
        +k(16)*n(idx_E)

    !d[H2_dot]/d[H2]
    pd(8,8) =  &
        +2.d0*k(27)*n(idx_H2)  &
        -k(50)  &
        -4.d0*k(27)*n(idx_H2)  &
        -k(43)  &
        -k(31)*n(idx_HEj)  &
        -k(51)  &
        -k(18)*n(idx_H)  &
        -k(52)  &
        -k(40)  &
        -k(32)*n(idx_HE)  &
        -k(44)  &
        -k(49)  &
        -k(17)*n(idx_E)  &
        -k(16)*n(idx_E)  &
        -k(33)*n(idx_HEj)  &
        -k(15)*n(idx_Hj)  &
        -k(14)*n(idx_Hj)

    !d[H2+_dot]/d[H2]
    pd(9,8) =  &
        +k(40)  &
        +k(15)*n(idx_Hj)  &
        +k(51)  &
        +k(33)*n(idx_HEj)  &
        +k(14)*n(idx_Hj)

    !d[Tgas_dot]/d[H2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(8)*1d-3
    if(dnn>0.d0) then
      nn(8) = n(8) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,8) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[H2+]
    pd(1,9) =  &
        +2.d0*k(24)*n(idx_E)  &
        +2.d0*k(25)*n(idx_E)  &
        -k(13)*n(idx_H)  &
        +k(41)  &
        +k(26)*n(idx_Hk)

    !d[E_dot]/d[H2+]
    pd(2,9) =  &
        -k(24)*n(idx_E)  &
        -k(25)*n(idx_E)  &
        +k(42)

    !d[H+_dot]/d[H2+]
    pd(3,9) =  &
        +k(13)*n(idx_H)  &
        +k(41)  &
        +2.d0*k(42)

    !d[H-_dot]/d[H2+]
    pd(7,9) =  &
        -k(26)*n(idx_Hk)

    !d[H2_dot]/d[H2+]
    pd(8,9) =  &
        +k(13)*n(idx_H)  &
        +k(26)*n(idx_Hk)

    !d[H2+_dot]/d[H2+]
    pd(9,9) =  &
        -k(24)*n(idx_E)  &
        -k(26)*n(idx_Hk)  &
        -k(41)  &
        -k(25)*n(idx_E)  &
        -k(13)*n(idx_H)  &
        -k(42)

    !d[Tgas_dot]/d[H2+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(9)*1d-3
    if(dnn>0.d0) then
      nn(9) = n(9) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,9) = (dn1-dn0)/dnn
    end if

    !d[Tgas_dot]/d[CR]
    pd(12,10) = 0.d0

    !d[Tgas_dot]/d[g]
    pd(12,11) = 0.d0

    !d[H_dot]/d[Tgas]
    pd(1,12) = 0.d0

    !d[E_dot]/d[Tgas]
    pd(2,12) = 0.d0

    !d[H+_dot]/d[Tgas]
    pd(3,12) = 0.d0

    !d[HE_dot]/d[Tgas]
    pd(4,12) = 0.d0

    !d[HE+_dot]/d[Tgas]
    pd(5,12) = 0.d0

    !d[HE++_dot]/d[Tgas]
    pd(6,12) = 0.d0

    !d[H-_dot]/d[Tgas]
    pd(7,12) = 0.d0

    !d[H2_dot]/d[Tgas]
    pd(8,12) = 0.d0

    !d[H2+_dot]/d[Tgas]
    pd(9,12) = 0.d0

    !d[CR_dot]/d[Tgas]
    pd(10,12) = 0.d0

    !d[g_dot]/d[Tgas]
    pd(11,12) = 0.d0

    !d[Tgas_dot]/d[Tgas]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(12)*1d-3
    if(dnn>0.d0) then
      nn(12) = n(12) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,12) = (dn1-dn0)/dnn
    end if

    !d[dummy_dot]/d[Tgas]
    pd(13,12) = 0.d0

    !d[Tgas_dot]/d[dummy]
    pd(12,13) = 0.d0

  end subroutine jex

end module krome_ode

!############### MODULE ##############
module krome_user
  implicit none

  ! *************************************************************
  !  This file has been generated with:
  !  krome 14.08.dev "Beastie Boyle" on 2016-08-08 10:48:19.
  !
  ! KROME is a nice and friendly chemistry package for a wide range of
  !  astrophysical simulations. Given a chemical network (in CSV format)
  !  it automatically generates all the routines needed to solve the kinetic
  !  of the system, modelled as system of coupled Ordinary Differential
  !  Equations.
  !  It provides different options which make it unique and very flexible.
  !  Any suggestions and comments are welcomed. KROME is an open-source
  !  package, GNU-licensed, and any improvements provided by
  !  the users is well accepted. See disclaimer below and GNU License
  !  in gpl-3.0.txt.
  !
  !  more details in http://kromepackage.org/
  !  also see https://bitbucket.org/krome/krome_stable
  !
  ! Written and developed by Tommaso Grassi
  !  tommasograssi@gmail.com,
  !  Starplan Center, Copenhagen.
  !  Niels Bohr Institute, Copenhagen.
  !
  ! Co-developer Stefano Bovino
  !  sbovino@astro.physik.uni-goettingen.de
  !  Institut fuer Astrophysik, Goettingen.
  !
  ! Others (alphabetically): D. Galli, F.A. Gianturco, T. Haugboelle,
  !  J.Prieto, J.Ramsey, D.R.G. Schleicher, D. Seifried,
  !  E. Simoncini, E. Tognelli
  !
  !
  ! KROME is provided "as it is", without any warranty.
  !  The Authors assume no liability for any damages of any kind
  !  (direct or indirect damages, contractual or non-contractual
  !  damages, pecuniary or non-pecuniary damages), directly or
  !  indirectly derived or arising from the correct or incorrect
  !  usage of KROME, in any possible environment, or arising from
  !  the impossibility to use, fully or partially, the software,
  !  or any bug or malefunction.
  !  Such exclusion of liability expressly includes any damages
  !  including the loss of data of any kind (including personal data)
  ! *************************************************************

  integer,parameter::KROME_idx_H = 1	!H
  integer,parameter::KROME_idx_E = 2	!E
  integer,parameter::KROME_idx_Hj = 3	!H+
  integer,parameter::KROME_idx_HE = 4	!HE
  integer,parameter::KROME_idx_HEj = 5	!HE+
  integer,parameter::KROME_idx_HEjj = 6	!HE++
  integer,parameter::KROME_idx_Hk = 7	!H-
  integer,parameter::KROME_idx_H2 = 8	!H2
  integer,parameter::KROME_idx_H2j = 9	!H2+
  integer,parameter::KROME_idx_CR = 10	!CR
  integer,parameter::KROME_idx_g = 11	!g
  integer,parameter::KROME_idx_Tgas = 12	!Tgas
  integer,parameter::KROME_idx_dummy = 13	!dummy

  integer,parameter::krome_idx_cool_H2 = 1
  integer,parameter::krome_idx_cool_H2GP = 2
  integer,parameter::krome_idx_cool_atomic = 3
  integer,parameter::krome_idx_cool_CEN = 3
  integer,parameter::krome_idx_cool_HD = 4
  integer,parameter::krome_idx_cool_Z = 5
  integer,parameter::krome_idx_cool_metal = 5
  integer,parameter::krome_idx_cool_enthalpic = 6
  integer,parameter::krome_idx_cool_dH = 6
  integer,parameter::krome_idx_cool_dust = 7
  integer,parameter::krome_idx_cool_compton = 8
  integer,parameter::krome_idx_cool_CIE = 9
  integer,parameter::krome_idx_cool_cont = 10
  integer,parameter::krome_idx_cool_continuum = 10
  integer,parameter::krome_idx_cool_expansion = 11
  integer,parameter::krome_idx_cool_exp = 11
  integer,parameter::krome_idx_cool_ff = 12
  integer,parameter::krome_idx_cool_bss = 12
  integer,parameter::krome_idx_cool_custom = 13
  integer,parameter::krome_idx_cool_CO = 14
  integer,parameter::krome_idx_cool_ZCIE = 15
  integer,parameter::krome_idx_cool_ZCIENOUV = 16
  integer,parameter::krome_idx_cool_ZExtend = 17
  integer,parameter::krome_idx_cool_gh = 18
  integer,parameter::krome_ncools = 18

  integer,parameter::krome_idx_heat_chem = 1
  integer,parameter::krome_idx_heat_compress = 2
  integer,parameter::krome_idx_heat_compr = 2
  integer,parameter::krome_idx_heat_photo = 3
  integer,parameter::krome_idx_heat_enthalpic = 4
  integer,parameter::krome_idx_heat_dH = 4
  integer,parameter::krome_idx_heat_photoAv = 5
  integer,parameter::krome_idx_heat_Av = 5
  integer,parameter::krome_idx_heat_CR = 6
  integer,parameter::krome_idx_heat_dust = 7
  integer,parameter::krome_idx_heat_xray = 8
  integer,parameter::krome_idx_heat_visc = 9
  integer,parameter::krome_idx_heat_viscous = 9
  integer,parameter::krome_idx_heat_custom = 10
  integer,parameter::krome_idx_heat_ZCIE = 11
  integer,parameter::krome_nheats = 11

  integer,parameter::krome_nrea=52
  integer,parameter::krome_nmols=9
  integer,parameter::krome_nspec=13
  integer,parameter::krome_natoms=3
  integer,parameter::krome_ndust=0
  integer,parameter::krome_ndustTypes=0
  integer,parameter::krome_nPhotoBins=10
  integer,parameter::krome_nPhotoRates=8

  real*8,parameter::krome_boltzmann_eV = 8.617332478d-5 !eV / K
  real*8,parameter::krome_boltzmann_J = 1.380648d-23 !J / K
  real*8,parameter::krome_boltzmann_erg = 1.380648d-16 !erg / K
  real*8,parameter::krome_planck_eV = 4.135667516d-15 !eV s
  real*8,parameter::krome_planck_J = 6.62606957d-34 !J s
  real*8,parameter::krome_planck_erg = 6.62606957d-27 !erg s
  real*8,parameter::krome_iplanck_eV = 1d0/krome_planck_eV !1 / eV / s
  real*8,parameter::krome_iplanck_J = 1d0/krome_planck_J !1 / J / s
  real*8,parameter::krome_iplanck_erg = 1d0/krome_planck_erg !1 / erg / s
  real*8,parameter::krome_gravity = 6.674d-8 !cm3 / g / s2
  real*8,parameter::krome_e_mass = 9.10938188d-28 !g
  real*8,parameter::krome_p_mass = 1.67262158d-24 !g
  real*8,parameter::krome_n_mass = 1.674920d-24 !g
  real*8,parameter::krome_ip_mass = 1d0/krome_p_mass !1/g
  real*8,parameter::krome_clight = 2.99792458e10 !cm/s
  real*8,parameter::krome_pi = 3.14159265359d0 !#
  real*8,parameter::krome_eV_to_erg = 1.60217646d-12 !eV -> erg
  real*8,parameter::krome_ry_to_eV = 13.60569d0 !rydberg -> eV
  real*8,parameter::krome_ry_to_erg = 2.179872d-11 !rydberg -> erg
  real*8,parameter::krome_seconds_per_year = 365d0*24d0*3600d0 !yr -> s
  real*8,parameter::krome_km_to_cm = 1d5 !km -> cm
  real*8,parameter::krome_cm_to_Mpc = 1.d0/3.08d24 !cm -> Mpc
  real*8,parameter::krome_kvgas_erg = 8.d0*krome_boltzmann_erg/krome_pi/krome_p_mass !
  real*8,parameter::krome_pre_kvgas_sqrt = sqrt(8.d0*krome_boltzmann_erg/krome_pi) !
  real*8,parameter::krome_pre_planck = 2.d0*krome_planck_erg/krome_clight**2 !erg/cm2*s3
  real*8,parameter::krome_exp_planck = krome_planck_erg / krome_boltzmann_erg !s*K
  real*8,parameter::krome_stefboltz_erg = 5.670373d-5 !erg/s/cm2/K4
  real*8,parameter::krome_N_avogadro = 6.0221d23 !#
  real*8,parameter::krome_Rgas_J = 8.3144621d0 !J/K/mol
  real*8,parameter::krome_Rgas_kJ = 8.3144621d-3 !kJ/K/mol
  real*8,parameter::krome_hubble = 0.704d0 !dimensionless
  real*8,parameter::krome_Omega0 = 1.0d0 !dimensionless
  real*8,parameter::krome_Omegab = 0.0456d0 !dimensionless
  real*8,parameter::krome_Hubble0 = 1.d2*krome_hubble*krome_km_to_cm*krome_cm_to_Mpc !1/s

contains

  subroutine krome_set_user_cell_size(argset)
    use krome_commons
    implicit none
    real*8 :: argset
    user_cell_size = argset
  end subroutine krome_set_user_cell_size

  !*******************
  function krome_get_user_cell_size()
    use krome_commons
    implicit none
    real*8 :: krome_get_user_cell_size
    krome_get_user_cell_size = user_cell_size
  end function krome_get_user_cell_size

  !*******************
  subroutine krome_set_user_crate(argset)
    use krome_commons
    implicit none
    real*8 :: argset
    user_crate = argset
  end subroutine krome_set_user_crate

  !*******************
  function krome_get_user_crate()
    use krome_commons
    implicit none
    real*8 :: krome_get_user_crate
    krome_get_user_crate = user_crate
  end function krome_get_user_crate

  !*******************
  subroutine krome_set_user_myfluxLW(argset)
    use krome_commons
    implicit none
    real*8 :: argset
    user_myfluxLW = argset
  end subroutine krome_set_user_myfluxLW

  !*******************
  function krome_get_user_myfluxLW()
    use krome_commons
    implicit none
    real*8 :: krome_get_user_myfluxLW
    krome_get_user_myfluxLW = user_myfluxLW
  end function krome_get_user_myfluxLW

  !*******************
  subroutine krome_set_user_myH2_dissociation(argset)
    use krome_commons
    implicit none
    real*8 :: argset
    user_myH2_dissociation = argset
  end subroutine krome_set_user_myH2_dissociation

  !*******************
  function krome_get_user_myH2_dissociation()
    use krome_commons
    implicit none
    real*8 :: krome_get_user_myH2_dissociation
    krome_get_user_myH2_dissociation = user_myH2_dissociation
  end function krome_get_user_myH2_dissociation

  !************************
  !returns the Tdust averaged over the number density
  ! as computed in the tables
  function krome_get_table_Tdust(x,Tgas)
    use krome_commons
    use krome_subs
    implicit none
    real*8 :: Tgas
    real*8 :: x(nmols), krome_get_table_Tdust
    real*8::ntot

    ntot = sum(x)
    krome_get_table_Tdust = 1d1**fit_anytab2D(dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_Tdust(:,:), dust_mult_ngas, &
        dust_mult_Tgas, &
        log10(ntot), log10(Tgas))

  end function krome_get_table_Tdust

  !**********************
  !convert from MOCASSIN abundances to KROME
  ! xmoc: MOCASSIN matrix (note: cm-3, real*4),
  ! imap: matrix position index map, integer
  ! returns KROME abundances (cm-3, real*8)
  function krome_convert_xmoc(xmoc,imap)
    use krome_commons
    use krome_subs
    implicit none
    real*4 :: xmoc(:,:)
    real*8 :: krome_convert_xmoc(nmols),x(nmols)
    integer :: imap(:)
    real*8::n(nspec)

    x(:) = 0d0

    x(idx_H) = xmoc(imap(1), 1)
    x(idx_Hj) = xmoc(imap(1), 2)
    x(idx_HE) = xmoc(imap(2), 1)
    x(idx_HEj) = xmoc(imap(2), 2)
    x(idx_HEjj) = xmoc(imap(2), 3)

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0
    x(idx_e) = get_electrons(n(:))
    krome_convert_xmoc(:) = x(:)

  end function krome_convert_xmoc

  !*************************
  !convert from KROME abundances to MOCASSIN
  ! xmoc: KROME matrix (cm-3, real*4),
  ! imap: matrix position index map, integer
  ! xmoc (out), matrix MOCASSIN abundances (cm-3, real*4)
  subroutine krome_return_xmoc(x,imap,xmoc)
    use krome_commons
    implicit none
    real*8 :: x(nmols)
    real*4 :: xmoc(:,:)
    integer :: imap(:)

    xmoc(:,:) = 0d0

    xmoc(imap(1), 1) = x(idx_H)
    xmoc(imap(1), 2) = x(idx_Hj)
    xmoc(imap(2), 1) = x(idx_HE)
    xmoc(imap(2), 2) = x(idx_HEj)
    xmoc(imap(2), 3) = x(idx_HEjj)

  end subroutine krome_return_xmoc

  !**********************
  !convert number density (cm-3) into column
  ! density (cm-2) using the specific density
  ! column method (see help for option
  ! -columnDensityMethod)
  ! num is the number density, x(:) is the species
  ! array, Tgas is the gas temperature
  ! If the method is not JEANS, x(:) and Tgas
  ! are dummy variables
  function krome_num2col(num,x,Tgas)
    use krome_subs
    use krome_commons
    implicit none
    real*8 :: x(nmols),krome_num2col
    real*8 :: Tgas,num
    real*8::n(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    krome_num2col = num2col(num,n(:))

  end function krome_num2col

  !***********************
  !print on screen the current values of the phys variables
  subroutine krome_print_phys_variables()
    use krome_commons
    implicit none

    print *, "Tcmb:", phys_Tcmb
    print *, "zredshift:", phys_zredshift
    print *, "orthoParaRatio:", phys_orthoParaRatio
    print *, "metallicity:", phys_metallicity
    print *, "Tfloor:", phys_Tfloor

  end subroutine krome_print_phys_variables

  !*******************
  subroutine krome_set_Tcmb(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_Tcmb = arg
  end subroutine krome_set_Tcmb

  !*******************
  function krome_get_Tcmb()
    use krome_commons
    implicit none
    real*8 :: krome_get_Tcmb
    krome_get_Tcmb = phys_Tcmb
  end function krome_get_Tcmb

  !*******************
  subroutine krome_set_zredshift(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_zredshift = arg
  end subroutine krome_set_zredshift

  !*******************
  function krome_get_zredshift()
    use krome_commons
    implicit none
    real*8 :: krome_get_zredshift
    krome_get_zredshift = phys_zredshift
  end function krome_get_zredshift

  !*******************
  subroutine krome_set_orthoParaRatio(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_orthoParaRatio = arg
  end subroutine krome_set_orthoParaRatio

  !*******************
  function krome_get_orthoParaRatio()
    use krome_commons
    implicit none
    real*8 :: krome_get_orthoParaRatio
    krome_get_orthoParaRatio = phys_orthoParaRatio
  end function krome_get_orthoParaRatio

  !*******************
  subroutine krome_set_metallicity(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_metallicity = arg
  end subroutine krome_set_metallicity

  !*******************
  function krome_get_metallicity()
    use krome_commons
    implicit none
    real*8 :: krome_get_metallicity
    krome_get_metallicity = phys_metallicity
  end function krome_get_metallicity

  !*******************
  subroutine krome_set_Tfloor(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_Tfloor = arg
  end subroutine krome_set_Tfloor

  !*******************
  function krome_get_Tfloor()
    use krome_commons
    implicit none
    real*8 :: krome_get_Tfloor
    krome_get_Tfloor = phys_Tfloor
  end function krome_get_Tfloor

  !*****************************
  !dump the data for restart (UNDER DEVELOPEMENT!)
  !arguments: the species array and the gas temperature
  subroutine krome_store(x,Tgas,dt)
    use krome_commons
    implicit none
    integer::i
    real*8 :: x(nmols)
    real*8 :: Tgas,dt

    open(krome_nfile,file=trim(krome_datafolder)//"krome_dump.dat",status="replace")
    !dump temperature
    write(krome_nfile,*) Tgas
    write(krome_nfile,*) dt
    !dump species
    do i=1,nmols
      write(krome_nfile,*) x(i)
    end do
    close(krome_nfile)

  end subroutine krome_store

  !*****************************
  !restore the data from a dump (UNDER DEVELOPEMENT!)
  !arguments: the species array and the gas temperature
  subroutine krome_restore(x,Tgas,dt)
    use krome_commons
    implicit none
    integer::i
    real*8 :: x(nmols)
    real*8 :: Tgas,dt

    open(krome_nfile,file=trim(krome_datafolder)//"krome_dump.dat",status="old")
    !restore temperature
    read(krome_nfile,*) Tgas
    read(krome_nfile,*) dt
    !restore species
    do i=1,nmols
      read(krome_nfile,*) x(i)
    end do
    close(krome_nfile)

  end subroutine krome_restore

  !****************************
  !switch on the thermal calculation
  subroutine krome_thermo_on()
    use krome_commons
    krome_thermo_toggle = 1
  end subroutine krome_thermo_on

  !****************************
  !switch off the thermal calculation
  subroutine krome_thermo_off()
    use krome_commons
    krome_thermo_toggle = 0
  end subroutine krome_thermo_off

  !************************
  ! prepares tables for cross sections and
  ! photorates
  subroutine krome_calc_photobins()
    use krome_photo
    call calc_photobins()
  end subroutine krome_calc_photobins

  !****************************
  ! set the energy per photo bin
  ! eV/cm2/sr
  subroutine krome_set_photoBinJ(phbin)
    use krome_commons
    use krome_photo
    implicit none
    real*8 :: phbin(nPhotoBins)
    photoBinJ(:) = phbin(:)
    photoBinJ_org(:) = phbin(:) !for restore

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBinJ

  !*************************
  ! set the energy (frequency) of the photobin
  ! as left-right limits in eV
  subroutine krome_set_photobinE_lr(phbinleft,phbinright)
    use krome_commons
    use krome_photo
    implicit none
    real*8 :: phbinleft(nPhotoBins),phbinright(nPhotoBins)
    photoBinEleft(:) = phbinleft(:)
    photoBinEright(:) = phbinright(:)
    photoBinEmid(:) = 0.5d0*(phbinleft(:)+phbinright(:))
    photoBinEdelta(:) = phbinright(:)-phbinleft(:)
    photoBinEidelta(:) = 1d0/photoBinEdelta(:)

    !initialize xsecs table
    call init_photoBins()

  end subroutine krome_set_photobinE_lr

  !*******************************
  !set the energy (eV) of the photobin according
  ! to MOCASSIN way (position and width array)
  subroutine krome_set_photobinE_moc(binPos,binWidth)
    use krome_commons
    use krome_photo
    implicit none
    real*8 :: binPos(nPhotoBins),binWidth(nPhotoBins)

    photoBinEleft(:) = binPos(:)-binWidth(:)/2d0
    photoBinEright(:) = binPos(:)+binWidth(:)/2d0
    photoBinEmid(:) = binPos(:)
    photoBinEdelta(:) = photoBinEright(:)-photoBinEleft(:)
    photoBinEidelta(:) = 1d0/photoBinEdelta(:)

    !initialize xsecs table
    call init_photoBins()

  end subroutine krome_set_photobinE_moc

  !********************************
  ! set the energy (eV) of the photobin
  ! linearly from lowest to highest energy value
  ! in eV
  subroutine krome_set_photobinE_lin(lower,upper)
    use krome_commons
    use krome_photo
    implicit none
    real*8 :: lower,upper
    real*8::dE
    integer::i
    dE = abs(upper-lower)/nPhotoBins
    do i=1,nPhotoBins
      photoBinEleft(i) = dE*(i-1) + lower
      photoBinEright(i) = dE*i + lower
      photoBinEmid(i) = 0.5d0*(photoBinEleft(i)+photoBinEright(i))
    end do
    photoBinEdelta(:) = photoBinEright(:)-photoBinEleft(:)
    photoBinEidelta(:) = 1d0/photoBinEdelta(:)

    !initialize xsecs table
    call init_photoBins()

  end subroutine krome_set_photobinE_lin

  !********************************
  ! set the energy (eV) of the photobin
  ! logarithmically from lowest to highest energy value
  ! in eV
  subroutine krome_set_photobinE_log(lower,upper)
    use krome_commons
    use krome_photo
    implicit none
    real*8 :: lower,upper
    real*8::dE,logup,loglow
    integer::i
    if(lower.ge.upper) then
      print *,"ERROR: in  krome_set_photobinE_log lower >= upper limit!"
      stop
    end if
    loglow = log10(lower)
    logup = log10(upper)
    dE = 1d1**(abs(logup-loglow)/nPhotoBins)
    do i=1,nPhotoBins
      photoBinEleft(i) = 1d1**((i-1)*(logup-loglow)/nPhotoBins + loglow)
      photoBinEright(i) = 1d1**(i*(logup-loglow)/nPhotoBins + loglow)
      photoBinEmid(i) = 0.5d0*(photoBinEleft(i)+photoBinEright(i))
    end do
    photoBinEdelta(:) = photoBinEright(:)-photoBinEleft(:)
    photoBinEidelta(:) = 1d0/photoBinEdelta(:)

    !initialize xsecs table
    call init_photoBins()

  end subroutine krome_set_photobinE_log

  !*********************************
  !returns an array containing the flux for each photo bin
  ! in eV/cm2/sr
  function krome_get_photoBinJ()
    use krome_commons
    real*8 :: krome_get_photoBinJ(nPhotoBins)
    krome_get_photoBinJ(:) = photoBinJ(:)
  end function krome_get_photoBinJ

  !*********************************
  function krome_get_photoBinE_left()
    !returns an array of size krome_nPhotoBins with the
    ! left energy limits (eV)
    use krome_commons
    real*8 :: krome_get_photoBinE_left(nPhotoBins)
    krome_get_photoBinE_left(:) = photoBinEleft(:)
  end function krome_get_photoBinE_left

  !*********************************
  !returns an array of size krome_nPhotoBins with the
  ! right energy limits (eV)
  function krome_get_photoBinE_right()
    use krome_commons
    real*8 :: krome_get_photoBinE_right(nPhotoBins)
    krome_get_photoBinE_right(:) = photoBinEright(:)
  end function krome_get_photoBinE_right

  !*********************************
  !returns an array of size krome_nPhotoBins with the
  ! middle energy values (eV)
  function krome_get_photoBinE_mid()
    use krome_commons
    real*8 :: krome_get_photoBinE_mid(nPhotoBins)
    krome_get_photoBinE_mid(:) = photoBinEmid(:)
  end function krome_get_photoBinE_mid

  !*********************************
  !returns an array of size krome_nPhotoBins with the
  ! bin span (eV)
  function krome_get_photoBinE_delta()
    use krome_commons
    real*8 :: krome_get_photoBinE_delta(nPhotoBins)
    krome_get_photoBinE_delta(:) = photoBinEdelta(:)
  end function krome_get_photoBinE_delta

  !*********************************
  !returns an array of size krome_nPhotoBins with the
  ! inverse of the bin span (1/eV)
  function krome_get_photoBinE_idelta()
    use krome_commons
    real*8 :: krome_get_photoBinE_idelta(nPhotoBins)
    krome_get_photoBinE_idelta(:) = photoBinEidelta(:)
  end function krome_get_photoBinE_idelta

  !*********************************
  !returns an array of size krome_nPhotoBins with the
  ! integrated photo rates (1/s)
  function krome_get_photoBin_rates()
    use krome_commons
    real*8 :: krome_get_photoBin_rates(nPhotoRea)
    krome_get_photoBin_rates(:) = photoBinRates(:)
  end function krome_get_photoBin_rates

  !*********************************
  !returns an array of size krome_nPhotoBins containing
  ! the cross section (cm2) of the idx-th photoreaction
  function krome_get_xsec(idx)
    use krome_commons
    implicit none
    real*8 :: krome_get_xsec(nPhotoBins)
    integer :: idx

    krome_get_xsec(:) = photoBinJTab(idx,:)

  end function krome_get_xsec

  !*********************************
  !returns an array of size krome_nPhotoBins with the
  ! integrated photo heatings (erg/s)
  function krome_get_photoBin_heats()
    use krome_commons
    implicit none
    real*8 :: krome_get_photoBin_heats(nPhotoRea)
    krome_get_photoBin_heats(:) = photoBinHeats(:)

  end function krome_get_photoBin_heats

  !****************************
  !multiply all the bins by a factor real*8 xscale
  subroutine krome_photoBin_scale(xscale)
    use krome_commons
    use krome_photo
    implicit none
    real*8 :: xscale

    photoBinJ(:) = photoBinJ(:) * xscale

    !compute rates
    call calc_photobins()

  end subroutine krome_photoBin_scale

  !****************************
  !multiply all the bins by a real*8 array xscale(:)
  ! of size krome_nPhotoBins
  subroutine krome_photoBin_scale_array(xscale)
    use krome_commons
    use krome_photo
    implicit none
    real*8 :: xscale(nPhotoBins)

    photoBinJ(:) = photoBinJ(:) * xscale(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_photoBin_scale_array

  !********************************
  !restore the original flux (i.e. undo any rescale).
  ! the flux is automatically stored by the functions
  ! that set the flux, or by the function
  ! krome_photoBin_store()
  subroutine krome_photoBin_restore()
    use krome_commons
    implicit none

    photoBinJ(:) = photoBinJ_org(:)

  end subroutine krome_photoBin_restore

  !**********************
  !store flux to be restored with the subroutine
  ! krome_photoBin_restore later
  subroutine krome_photoBin_store()
    use krome_commons
    implicit none

    photoBinJ_org(:) = photoBinJ(:)

  end subroutine krome_photoBin_store

  !********************************
  !load the radiation bins from the file fname
  ! data should be a 3-column file with
  ! energy Left (eV), energy Right (eV)
  ! intensity (eV/cm2/sr).
  ! This subroutine sets also the bin-size
  subroutine krome_load_photoBin_file(fname) !! !! not yet callable from C
    use krome_commons
    implicit none
    integer::ios,icount
    character(len=*) :: fname
    real*8::tmp_El(nPhotoBins),tmp_Er(nPhotoBins)
    real*8::rout(3),tmp_J(nPhotoBins)

    !open file and check for errors
    open(krome_nfile,file=fname,status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: problem opening "//fname//"!"
      print *," (e.g. file not found)"
      stop
    end if

    icount = 0 !count valid line
    !loop on file
    do
      read(krome_nfile,*,iostat=ios) rout(:)
      if(ios==-1) exit !EOF
      if(ios.ne.0) cycle !skip comments
      icount = icount + 1
      if(icount>nPhotoBins) exit !can't load more than nPhotoBins
      tmp_El(icount) = rout(1) !energy L eV
      tmp_Er(icount) = rout(2) !energy R eV
      !check if left interval is before right
      if(tmp_El(icount)>tmp_Er(icount)) then
        print *,"ERROR: in file "//fname//" left"
        print *, " interval larger than right one!"
        print *,tmp_El(icount),tmp_Er(icount)
        stop
      end if
      tmp_J(icount) = rout(3) !intensity eV/cm2/sr
    end do
    close(krome_nfile)

    !file data lines should be the same number of the photobins
    if(icount/=nPhotoBins) then
      print *,"ERROR: the number of data lines in the file"
      print *," "//fname//" should be equal to the number of"
      print *," photobins ",nPhotoBins
      print *,"Found",icount
      stop
    end if

    !initialize inteval and indensity according to data
    call krome_set_photobinE_lr(tmp_El(:),tmp_Er(:))
    call krome_set_photoBinJ(tmp_J(:))

  end subroutine krome_load_photoBin_file

  !**********************************
  !this subroutine set a flux HM in the energy limits
  ! as argument
  subroutine krome_set_photoBin_HMlog(lower_in,upper_in)
    use krome_commons
    use krome_photo
    use krome_subs
    implicit none
    real*8::z(59),energy(500),HM(59,500)
    real*8::z_mul,energy_mul,x,lower,upper
    real*8,parameter::limit_lower = 0.1237d0
    real*8,parameter::limit_upper = 4.997d7
    real*8,parameter::limit_redshift = 15.660d0
    real*8, optional :: lower_in,upper_in
    integer::i

    lower = limit_lower
    upper = limit_upper
    if(present(lower_in)) lower = lower_in
    if(present(upper_in)) upper = upper_in

    if(phys_zredshift>limit_redshift) then
      print *,"ERROR: redshift out of range in HM"
      print *,"redshift:",phys_zredshift
      print *,"limit:",limit_redshift
      stop
    end if

    if(lower<limit_lower .or. upper>limit_upper) then
      print *,"ERROR: upper or lower limit out of range in HM."
      print *,"lower limit (eV):",limit_lower
      print *,"upper limit (eV):",limit_upper
      stop
    end if

    call krome_set_photoBinE_log(lower,upper)

    call init_anytab2D(trim(krome_datafolder)//"krome_HMflux.dat", z(:), energy(:), &
        HM(:,:), z_mul, energy_mul)

    do i=1,nPhotoBins
      x = log10(photoBinEmid(i)) !log(eV)
      photoBinJ(i) = 1d1**fit_anytab2D(z(:), energy(:), HM(:,:), &
          z_mul, energy_mul, phys_zredshift, x)
    end do

    photoBinJ_org(:) = photoBinJ(:)

    call calc_photobins()

  end subroutine krome_set_photoBin_HMlog

  !**********************************
  !this subroutine set a flux HM in the energy limits
  ! as argument
  subroutine krome_set_photoBin_HMCustom(lower_in,upper_in,additive)
    use krome_commons
    use krome_photo
    use krome_subs
    implicit none
    real*8::z(59),energy(500),HM(59,500)
    real*8::z_mul,energy_mul,x,lower,upper
    real*8::photoTmpJ(nPhotoBins)
    real*8,parameter::limit_lower = 0.1237d0
    real*8,parameter::limit_upper = 4.997d7
    real*8,parameter::limit_redshift = 15.660d0
    logical,optional::additive
    logical::add
    real*8, optional :: lower_in,upper_in
    integer::i

    lower = limit_lower
    upper = limit_upper
    if(present(lower_in)) lower = lower_in
    if(present(upper_in)) upper = upper_in

    add = .false.
    if(present(additive)) add = additive

    if(phys_zredshift>limit_redshift) then
      print *,"ERROR: redshift out of range in HM"
      print *,"redshift:",phys_zredshift
      print *,"limit:",limit_redshift
      stop
    end if

    if(lower<limit_lower .or. upper>limit_upper) then
      print *,"ERROR: upper or lower limit out of range in HM."
      print *,"lower limit (eV):",limit_lower
      print *,"upper limit (eV):",limit_upper
      stop
    end if

    call init_anytab2D(trim(krome_datafolder)//"krome_HMflux.dat", z(:), energy(:), &
        HM(:,:), z_mul, energy_mul)

    do i=1,nPhotoBins
      x = log10(photoBinEmid(i)) !log(eV)
      photoTmpJ(i) = 1d1**fit_anytab2D(z(:), energy(:), HM(:,:), &
          z_mul, energy_mul, phys_zredshift, x)
    end do

    !add flux to already-present flux if optional argument
    if(add) then
      photoBinJ(:) = photoBinJ(:) + photoTmpJ(:)
    else
      photoBinJ(:) = photoTmpJ(:)
    end if

    photoBinJ_org(:) = photoBinJ(:)

    call calc_photobins()

  end subroutine krome_set_photoBin_HMCustom

  !**********************************
  !set the flux as a black body with temperature Tbb (K)
  ! in the range lower to upper (eV). the spacing is linear
  subroutine krome_set_photoBin_BBlin(lower,upper,Tbb)
    use krome_commons
    use krome_constants
    use krome_photo
    use krome_subs
    implicit none
    real*8 :: lower,upper,Tbb
    real*8::x
    integer::i

    call krome_set_photoBinE_lin(lower,upper)

    !eV/cm2/sr
    do i=1,nPhotoBins
      x = photoBinEmid(i) !eV
      photoBinJ(i) = planckBB(x,Tbb)
    end do
    photoBinJ_org(:) = photoBinJ(:)

    call calc_photobins()

  end subroutine krome_set_photoBin_BBlin

  !**********************************
  !set the flux as a black body with temperature Tbb (K)
  ! in the range lower to upper (eV). the spacing is logarithmic
  subroutine krome_set_photoBin_BBlog(lower,upper,Tbb)
    use krome_commons
    use krome_constants
    use krome_photo
    use krome_subs
    implicit none
    real*8 :: lower,upper,Tbb
    real*8::x,xmax,xexp,Jlim
    integer::i

    !limit for the black body intensity to check limits
    Jlim = 1d-3

    call krome_set_photoBinE_log(lower,upper)

    !eV/cm2/sr
    do i=1,nPhotoBins
      x = photoBinEmid(i) !eV
      photoBinJ(i) = planckBB(x,Tbb)
    end do
    photoBinJ_org(:) = photoBinJ(:)

    !uncomment this below for additional control
    !!$    !find the maximum using Wien's displacement law
    !!$    xmax = Tbb/2.8977721d-1 * clight * planck_eV !eV
    !!$
    !!$    if(xmax<lower) then
    !!$       print *,"WARNING: maximum of the Planck function"
    !!$       print *," is below the lowest energy bin!"
    !!$       print *,"max (eV)",xmax
    !!$       print *,"lowest (eV)",lower
    !!$       print *,"Tbb (K)",Tbb
    !!$    end if
    !!$
    !!$    if(xmax>upper) then
    !!$       print *,"WARNING: maximum of the Planck function"
    !!$       print *," is above the highest energy bin!"
    !!$       print *,"max (eV)",xmax
    !!$       print *,"highest (eV)",upper
    !!$       print *,"Tbb (K)",Tbb
    !!$    end if
    !!$
    !!$    if(photoBinJ(1)>Jlim) then
    !!$       print *,"WARNING: lower bound of the Planck function"
    !!$       print *," has a flux of (ev/cm2/s/Hz/sr)",photoBinJ(1)
    !!$       print *," which is larger than the limit Jlim",Jlim
    !!$       print *,"Tbb (K)",Tbb
    !!$    end if
    !!$
    !!$    if(photoBinJ(nPhotoBins)>Jlim) then
    !!$       print *,"WARNING: upper bound of the Planck function"
    !!$       print *," has a flux of (ev/cm2/s/Hz/sr)",photoBinJ(nPhotoBins)
    !!$       print *," which is larger than the limit Jlim",Jlim
    !!$       print *,"Tbb (K)",Tbb
    !!$    end if

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_BBlog

  !*************************************
  !set the BB spectrum and the limits using bisection
  subroutine krome_set_photoBin_BBlog_auto(Tbb)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8 :: Tbb
    real*8::xlow,xup,eps,xmax,J0,J1,x0,x1,xm,Jm
    eps = 1d-6

    !Rayleigh–Jeans approximation for the minimum energy
    xlow = planck_eV*clight*sqrt(.5d0/Tbb/boltzmann_eV*eps)

    !find energy of the Wien maximum (eV)
    xmax = Tbb / 2.8977721d-1 * clight * planck_eV

    !bisection to find the maximum
    x0 = xmax
    x1 = 2.9d2*Tbb*boltzmann_eV
    J0 = planckBB(x0,Tbb) - eps
    J1 = planckBB(x1,Tbb) - eps
    if(J0<0d0.or.J1>0d0) then
      print *,"ERROR: problems with auto planck bisection!"
      stop
    end if

    do
      xm = 0.5d0*(x0+x1)
      Jm = planckBB(xm,Tbb) - eps
      if(Jm>0d0) x0 = xm
      if(Jm<0d0) x1 = xm
      if(abs(Jm)<eps*1d-3) exit
    end do
    xup = xm

    !initialize BB radiation using the values found
    call krome_set_photoBin_BBlog(xlow,xup,Tbb)

  end subroutine krome_set_photoBin_BBlog_auto

  !*********************************
  !return the ratio between the current flux an Draine's
  function krome_get_ratioFluxDraine()
    use krome_subs
    implicit none
    real*8::krome_get_ratioFluxDraine

    krome_get_ratioFluxDraine = get_ratioFluxDraine()

  end function krome_get_ratioFluxDraine

  !**********************************
  !set the flux as Draine's function
  ! in the range lower to upper (eV). the spacing is linear
  subroutine krome_set_photoBin_draineLin(lower,upper)
    use krome_commons
    use krome_photo
    use krome_constants
    real*8 :: upper,lower
    real*8::x
    integer::i

    call krome_set_photoBinE_lin(lower,upper)

    do i=1,nPhotoBins
      x = photoBinEmid(i) !eV
      !eV/cm2/sr
      if(x<13.6d0.and.x>5d0) then
        photoBinJ(i) = (1.658d6*x - 2.152d5*x**2 + 6.919d3*x**3) &
            * x *planck_eV
      else
        photoBinJ(i) = 0d0
      end if
    end do

    photoBinJ_org(:) = photoBinJ(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_draineLin

  !**************************
  !set the flux as Draine's function
  ! in the range lower to upper (eV). the spacing is logarithmic
  subroutine krome_set_photoBin_draineLog(lower,upper)
    use krome_commons
    use krome_photo
    use krome_constants
    real*8 :: upper,lower
    real*8::x
    integer::i

    call krome_set_photoBinE_log(lower,upper)

    do i=1,nPhotoBins
      x = photoBinEmid(i) !eV
      !eV/cm2/sr/s/Hz
      if(x<13.6d0.and.x>5d0) then
        photoBinJ(i) = (1.658d6*x - 2.152d5*x**2 + 6.919d3*x**3) &
            * x *planck_eV
      else
        photoBinJ(i) = 0d0
      end if
    end do

    photoBinJ_org(:) = photoBinJ(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_draineLog

  !**************************
  !set the flux as Draine's function
  ! in the range lower to upper (eV). the spacing is custom
  ! Note: you have to set the binning first
  subroutine krome_set_photoBin_draineCustom(lower,upper)
    use krome_commons
    use krome_photo
    use krome_constants
    real*8 :: upper,lower
    real*8::x
    integer::i

    if(maxval(photoBinEmid)==0d0) then
      print *,"ERROR: not initialized binning in draineCustom!"
      stop
    end if

    do i=1,nPhotoBins
      x = photoBinEmid(i) !eV
      !eV/cm2/sr/s/Hz
      if(x<13.6d0.and.x>5d0) then
        photoBinJ(i) = (1.658d6*x - 2.152d5*x**2 + 6.919d3*x**3) &
            * x *planck_eV
      else
        photoBinJ(i) = 0d0
      end if
    end do

    photoBinJ_org(:) = photoBinJ(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_draineCustom

  !**************************
  !set the flux as power-law (J21-style)
  ! in the range lower to upper (eV). the spacing is linear
  subroutine krome_set_photoBin_J21lin(lower,upper)
    use krome_commons
    use krome_photo
    real*8 :: upper,lower

    call krome_set_photoBinE_lin(lower,upper)
    photoBinJ(:) = 6.2415d-10 * (13.6d0/photoBinEmid(:)) !eV/cm2/s/Hz/sr
    photoBinJ_org(:) = photoBinJ(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_J21lin

  !**************************
  !set the flux as power-law (J21-style)
  ! in the range lower to upper (eV). the spacing is logarithmic
  subroutine krome_set_photoBin_J21log(lower,upper)
    use krome_commons
    use krome_photo
    real*8 :: upper,lower

    call krome_set_photoBinE_log(lower,upper)
    photoBinJ(:) = 6.2415d-10 * (13.6d0/photoBinEmid(:)) !eV/cm2/s/Hz/sr
    photoBinJ_org(:) = photoBinJ(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_J21log

  !*****************************
  !get the opacity exp(-tau) correpsonding the to x(:)
  ! chemical composition. The column density
  ! is computed using the expression in the
  ! num2col(x) function.
  ! An array of size krome_nPhotoBins is returned.
  function krome_get_opacity(x,Tgas)
    use krome_commons
    use krome_constants
    use krome_photo
    use krome_subs
    implicit none
    real*8 :: x(nmols),krome_get_opacity(nPhotoBins)
    real*8 :: Tgas
    real*8::tau,n(nspec)
    integer::i,j,idx

    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    !loop on frequency bins
    do j=1,nPhotoBins
      tau = 0d0
      !loop on species
      do i=1,nPhotoRea
        !calc opacity as column_density * cross_section
        idx = photoPartners(i)
        tau = tau + num2col(x(idx),n(:)) * photoBinJTab(i,j)
      end do
      krome_get_opacity(j) = tau !store
    end do

  end function krome_get_opacity

  !*****************************
  !get the opacity exp(-tau) correpsonding to the x(:)
  ! chemical composition. The column density
  ! is computed using the size of the cell (csize)
  ! An array of size krome_nPhotoBins is returned.
  function krome_get_opacity_size(x,Tgas,csize)
    use krome_commons
    use krome_constants
    use krome_photo
    use krome_subs
    use krome_dust
    implicit none
    real*8 :: x(nmols),krome_get_opacity_size(nPhotoBins)
    real*8 :: Tgas,csize
    real*8::n(nspec),energy,tau
    integer::i,j,idx

    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    !loop on frequency bins
    do j=1,nPhotoBins
      tau = 0d0
      !loop on species
      do i=1,nPhotoRea
        !calc opacity as column_density * cross_section
        !where column_density is density*cell_size
        idx = photoPartners(i)
        tau = tau + x(idx) * photoBinJTab(i,j)
      end do

      krome_get_opacity_size(j) = tau * csize !store
    end do

  end function krome_get_opacity_size

  !*****************************
  !get the opacity exp(-tau) correpsonding to the x(:)
  ! chemical composition. The column density
  ! is computed using the size of the cell (csize).
  ! Dust is included using dust-to-gas mass ratio (d2g).
  ! You should load the dust tables with the subroutine
  ! krome_load_dust_opacity(fileName).
  ! An array of size krome_nPhotoBins is returned.
  function krome_get_opacity_size_d2g(x,Tgas,csize,d2g)
    use krome_commons
    use krome_constants
    use krome_photo
    use krome_subs
    use krome_dust
    implicit none
    real*8 :: x(nmols),krome_get_opacity_size_d2g(nPhotoBins)
    real*8 :: Tgas,csize,d2g
    real*8::n(nspec),energy,tau,m(nspec),mgas
    integer::i,j,idx

    m(:) = get_mass()
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    mgas = sum(n(1:nmols)*m(1:nmols))

    !loop on frequency bins
    do j=1,nPhotoBins
      tau = 0d0
      !loop on species
      do i=1,nPhotoRea
        !calc opacity as column_density * cross_section
        !where column_density is density*cell_size
        idx = photoPartners(i)
        tau = tau + x(idx) * photoBinJTab(i,j)
      end do

      !sum dust opacity from interpolated table
      tau = tau + d2g*mgas*opacityDust(j)

      krome_get_opacity_size_d2g(j) = tau * csize !store
    end do

  end function krome_get_opacity_size_d2g

  !*******************************
  !load a frequency-dependent opacity table stored in fname file,
  ! column 1 is energy or wavelenght in un units of unitEnergy
  ! (default eV), column 2 is opacity in cm2/g.
  ! opacity is interpolated over the current photo-binning.
  subroutine krome_load_opacity_table(fname,unitEnergy)
    use krome_commons
    use krome_constants
    implicit none
    integer,parameter::ntmp=int(1e5)
    character(len=*)::fname
    character(len=*),optional::unitEnergy
    character*10::eunit
    integer::ios,icount,iR,iL,i,j
    real*8::wl,opac,fL,fR,kk,dE
    real*8::wls(ntmp),opacs(ntmp)
    real*8,allocatable::energy(:),kappa(:)

    !read energy unit optional argument
    eunit = "eV" !default is eV
    if(present(unitEnergy)) then
      eunit = trim(unitEnergy)
    end if

    !read form file
    open(krome_nfile,file=trim(fname),status="old",iostat=ios)
    !error if problems reading file
    if(ios/=0) then
      print *,"ERROR: problem while loading "//trim(fname)
      stop
    end if
    icount = 0
    !loop on file lines
    do
      !read wavelength and opacity
      read(krome_nfile,*,iostat=ios) wl,opac
      if(ios/=0) exit
      icount = icount + 1
      wls(icount) = wl
      opacs(icount) = opac
    end do
    close(krome_nfile)

    !allocate arrays
    allocate(energy(icount), kappa(icount))
    !copy temp arrays into allocated arrays, converting units
    if(trim(eunit)=="eV") then
      !eV->eV (default)
      kappa(:) = opacs(1:icount)
      energy(:) = wls(1:icount)
    elseif(trim(eunit)=="micron") then
      !micron->eV
      kappa(:) = opacs(1:icount)
      energy(:) = planck_eV*clight/(wls(1:icount)*1d-4)
    else
      print *,"ERROR: in load opacity table energy unit unknow",trim(eunit)
      stop
    end if

    !reverse array if necessary
    if(energy(2)<energy(1)) then
      energy(:) = energy(size(energy):1:-1)
      kappa(:) = kappa(size(kappa):1:-1)
    end if

    !check lower limit
    if(photoBinEleft(1)<energy(1)) then
      print *,"ERROR: dust table "//trim(fname)//" lower bound"
      print *,photoBinEleft(1), "<", energy(1)
      stop
    end if

    !check upper limit
    if(photoBinEright(nPhotoBins)>energy(size(energy))) then
      print *,"ERROR: dust table "//trim(fname)//" upper bound"
      print *,photoBinEright(nPhotoBins), ">", energy(size(energy))
      stop
    end if

    !interpolate on current energy distribution
    do j=1,nPhotoBins
      do i=2,size(energy)
        !find left bound position
        if(photoBinEleft(j)>energy(i-1) &
            .and. photoBinEleft(j)<energy(i)) then
        dE = energy(i)-energy(i-1)
        fL = (photoBinEleft(j)-energy(i-1))/dE &
            * (kappa(i)-kappa(i-1)) + kappa(i-1)
        iL = i
      end if

      !find right bound position
      if(photoBinEright(j)>energy(i-1) &
          .and. photoBinEright(j)<energy(i)) then
      dE = energy(i)-energy(i-1)
      fR = (photoBinEright(j)-energy(i-1))/dE &
          * (kappa(i)-kappa(i-1)) + kappa(i-1)
      iR = i
    end if
  end do

  !sum opacity for the given photo bin
  kk = 0d0
  !if there are other opacity points in between left and right limits
  if(iR-iL>0) then
    kk = kk + (energy(iL)-photoBinEleft(j))*(fL+kappa(iL))/2d0
    kk = kk + (photoBinEright(j)-energy(iR-1))*(fR+kappa(iR-1))/2d0
    !sum points in between
    do i=iL,iR-2
      kk = kk + (energy(i+1)-energy(i))*(kappa(i+1)+kappa(i))/2d0
    end do
  elseif(iR==iL) then
    !no opacity points in between
    kk = kk + (fL+fR)*(photoBinEright(j)-photoBinEleft(j))/2d0
  else
    print *,"ERROR: dust opacity interpolation error, iR-iL<0!"
    print *,"iR,iL:",iR,iL
    stop
  end if

  !copy to common and scale to bin size
  dE = photoBinEright(j)-photoBinEleft(j)
  opacityDust(j) = kk/dE

end do

!dump interpolated opacity
open(krome_nfile,file=trim(krome_datafolder)//"opacityDust.interp",status="replace")
do j=1,nPhotoBins
  write(krome_nfile,*) photoBinEmid(j),opacityDust(j)
end do
close(23)

!dump original opacity file (as loaded by krome)
open(krome_nfile,file=trim(krome_datafolder)//"opacityDust.org",status="replace")
do i=1,size(energy)
  write(krome_nfile,*) energy(i),kappa(i)
end do
close(krome_nfile)

end subroutine krome_load_opacity_table

!***************************
!alias for coe in krome_subs
! returns the coefficient array of size krome_nrea
! for a given Tgas
function krome_get_coef(Tgas)
use krome_commons
use krome_subs
real*8 :: krome_get_coef(nrea),Tgas
real*8::n(nspec)
n(:) = 0d0
n(idx_Tgas) = Tgas

krome_get_coef(:) = coe(n(:))

end function krome_get_coef

!****************************
!get the mean molecular weight from
! mass fractions
function krome_get_mu_x(xin)
use krome_commons
implicit none
real*8 :: xin(nmols), krome_get_mu_x
real*8::n(nmols)
n(:) = krome_x2n(xin(:),1d0)
krome_get_mu_x = krome_get_mu(n(:))
end function krome_get_mu_x

!****************************
!return the adiabatic index from mass fractions
! and temperature in K
function krome_get_gamma_x(xin,inTgas)
use krome_commons
implicit none
real*8 :: inTgas
real*8 :: xin(nmols), krome_get_gamma_x
real*8::x(nmols),Tgas,rhogas

Tgas = inTgas
x(:) = krome_x2n(xin(:),1d0)
krome_get_gamma_x = krome_get_gamma(x(:),Tgas)

end function krome_get_gamma_x

!***************************
!normalize mass fractions and
! set charge to zero
subroutine krome_consistent_x(x)
use krome_commons
use krome_constants
implicit none
real*8 :: x(nmols)
real*8::isumx,sumx,xerr,imass(nmols),ee

!1. charge consistency
imass(:) = krome_get_imass()

x(idx_e) = 0.d0

ee = sum(krome_get_charges()*x(:)*imass(:))
ee = max(ee*e_mass,0d0)
x(idx_e) = ee

!2. mass fraction consistency
sumx = sum(x)

!NOTE: uncomment here if you want some additional control
!conservation error threshold: rise an error if above xerr
!xerr = 1d-2
!if(abs(sum-1d0)>xerr) then
!   print *,"ERROR: some problem with conservation!"
!   print *,"|sum(x)-1|=",abs(sum-1d0)
!   stop
!end if

isumx = 1d0/sumx
x(:) = x(:) * isumx

end subroutine krome_consistent_x

!*********************
!return an array sized krome_nmols containing
! the mass fractions (#), computed from the number
! densities (1/cm3) and the total density in g/cm3
function krome_n2x(n,rhogas)
use krome_commons
implicit none
real*8 :: n(nmols)
real*8 :: rhogas
real*8 :: krome_n2x(nmols)

krome_n2x(:) = n(:) * krome_get_mass() / rhogas

end function krome_n2x

!********************
!return an array sized krome_nmols containing
! the number densities (1/cm3), computed from the mass
! fractions and the total density in g/cm3
function krome_x2n(x,rhogas)
use krome_commons
implicit none
real*8 :: x(nmols),rhogas,krome_x2n(nmols)

!compute densities from fractions
krome_x2n(:) = rhogas * x(:) * krome_get_imass()

end function krome_x2n

!*******************
!do only cooling and heating
subroutine krome_thermo(x,Tgas,dt)
use krome_commons
use krome_cooling
use krome_heating
use krome_subs
use krome_tabs
use krome_constants
implicit none
real*8 :: x(nmols)
real*8 :: Tgas,dt
real*8::n(nspec),nH2dust,dTgas,k(nrea),krome_gamma

nH2dust = 0d0
n(:) = 0d0
n(idx_Tgas) = Tgas
n(1:nmols) = x(:)
k(:) = coe_tab(n(:)) !compute coefficients
krome_gamma = gamma_index(n(:))

!dTgas = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
!    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))

!Tgas = Tgas + dTgas*dt !update gas

end subroutine krome_thermo

!*************************
!get heating (erg/cm3/s) for a given species
! array x(:) and Tgas
function krome_get_heating(x,inTgas)
use krome_heating
use krome_subs
use krome_commons
implicit none
real*8 :: inTgas
real*8 :: x(nmols), krome_get_heating
real*8::Tgas,k(nrea),nH2dust,n(nspec)
n(1:nmols) = x(:)
Tgas = inTgas
n(idx_Tgas) = Tgas
k(:) = coe(n(:))
nH2dust = 0d0
krome_get_heating = heating(n(:),Tgas,k(:),nH2dust)
end function krome_get_heating

!*****************************
! get an array containing individual heatings (erg/cm3/s)
! the array has size krome_nheats. see heatcool.gps
! for index list
function krome_get_heating_array(x,inTgas)
use krome_heating
use krome_subs
use krome_commons
implicit none
real*8::n(nspec),Tgas,k(nrea),nH2dust
real*8 :: x(nmols),krome_get_heating_array(nheats),inTgas

n(:) = 0d0
n(1:nmols) = x(:)
n(idx_Tgas) = inTgas
!#KROME_Tdust_copy
k(:) = coe(n(:))
Tgas = inTgas
nH2dust = 0d0
krome_get_heating_array(:) = get_heating_array(n(:),Tgas,k(:),nH2dust)

end function krome_get_heating_array

!*************************
!get cooling (erg/cm3/s) for x(:) species array
! and Tgas
function krome_get_cooling(x,inTgas)
use krome_cooling
use krome_commons
implicit none
real*8 :: inTgas
real*8 :: x(nmols), krome_get_cooling
real*8::Tgas,n(nspec)
n(1:nmols) = x(:)
Tgas = inTgas
n(idx_Tgas) = Tgas
krome_get_cooling = cooling(n,Tgas)
end function krome_get_cooling

!*****************************
! get an array containing individual coolings (erg/cm3/s)
! the array has size krome_ncools. see heatcool.gps
! for index list
function krome_get_cooling_array(x,inTgas)
use krome_cooling
use krome_commons
implicit none
real*8::n(nspec),Tgas
real*8 :: x(nmols),krome_get_cooling_array(ncools),inTgas

n(:) = 0d0
n(1:nmols) = x(:)
n(idx_Tgas) = inTgas
!#KROME_Tdust_copy
Tgas = inTgas
krome_get_cooling_array(:) = get_cooling_array(n(:),Tgas)

end function krome_get_cooling_array

!******************
!alias of plot_cool
subroutine krome_plot_cooling(n)
use krome_cooling
implicit none
real*8 :: n(krome_nmols)

call plot_cool(n(:))

end subroutine krome_plot_cooling

!************************
!conserve the total amount of nucleii,
! alias for conserveLin_x in subs
subroutine krome_conserveLin_x(x,ref)
use krome_commons
use krome_subs
implicit none
real*8 :: x(nmols),ref(natoms)

call conserveLin_x(x(:),ref(:))

end subroutine krome_conserveLin_x

!************************
!conserve the total amount of nucleii,
! alias for conserveLin_x in subs
function krome_conserveLinGetRef_x(x)
use krome_commons
use krome_subs
implicit none
real*8 :: x(nmols),krome_conserveLinGetRef_x(natoms)

krome_conserveLinGetRef_x(:) = &
    conserveLinGetRef_x(x(:))

end function krome_conserveLinGetRef_x

!*************************
!force conservation to array x(:)
!using xi(:) as initial abundances.
!alias for conserve in krome_subs
function krome_conserve(x,xi)
use krome_subs
implicit none
real*8 :: x(krome_nmols),xi(krome_nmols),krome_conserve(krome_nmols)
real*8::n(krome_nspec),ni(krome_nspec)

n(:) = 0d0
ni(:) = 0d0
n(1:krome_nmols) = x(1:krome_nmols)
ni(1:krome_nmols) = xi(1:krome_nmols)
n(:) = conserve(n(:), ni(:))
krome_conserve(:) = n(1:krome_nmols)

end function krome_conserve

!***************************
!get the adiabatic index for x(:) species abundances
! and Tgas.
! alias for gamma_index in krome_subs
function krome_get_gamma(x,Tgas)
use krome_subs
use krome_commons
real*8 :: Tgas
real*8 :: x(nmols), krome_get_gamma
real*8::n(nspec)
n(:) = 0.d0
n(1:nmols) = x(:)
n(idx_Tgas) = Tgas
krome_get_gamma = gamma_index(n(:))
end function krome_get_gamma

!***************************
!get an integer array containing the atomic numbers Z
! of the spcecies.
! alias for get_zatoms
function krome_get_zatoms()
use krome_subs
use krome_commons
implicit none
integer :: krome_get_zatoms(nmols)
integer::zatoms(nspec)

zatoms(:) = get_zatoms()
krome_get_zatoms(:) = zatoms(1:nmols)

end function krome_get_zatoms

!****************************
!get the mean molecular weight from
! number density and mass density.
! alias for get_mu in krome_subs module
function krome_get_mu(x)
use krome_commons
use krome_subs
implicit none
real*8 :: x(nmols), krome_get_mu
real*8::n(1:nspec)
n(:) = 0d0
n(1:nmols) = x(:)
krome_get_mu = get_mu(n(:))
end function krome_get_mu

!***************************
!get the names of the reactions as a
! character*50 array of krome_nrea
! elements
function krome_get_rnames() !! !! cannot yet be called from C
use krome_commons
use krome_subs
implicit none
character*50 :: krome_get_rnames(nrea)

krome_get_rnames(:) = get_rnames()

end function krome_get_rnames

!*****************
!get an array of double containing the masses in g
! of the species.
! alias for get_mass in krome_subs
function krome_get_mass()
use krome_subs
use krome_commons
implicit none
real*8::tmp(nspec)
real*8 :: krome_get_mass(nmols)
tmp(:) = get_mass()
krome_get_mass = tmp(1:nmols)
end function krome_get_mass

!*****************
!get an array of double containing the inverse
! of the mass (1/g) of the species
!alias for get_imass in krome_subs
function krome_get_imass()
use krome_subs
use krome_commons
implicit none
real*8::tmp(nspec)
real*8 :: krome_get_imass(nmols)
tmp(:) = get_imass()
krome_get_imass = tmp(1:nmols)
end function krome_get_imass

!***********************
!get the total number of H nuclei
function krome_get_Hnuclei(x)
use krome_commons
use krome_subs
real*8::n(nspec)
real*8 :: krome_get_Hnuclei, x(nmols)
n(:) = 0d0
n(1:nmols) = x(:)

krome_get_Hnuclei = get_Hnuclei(n(:))

end function krome_get_Hnuclei

!*****************
!get an array of size krome_nmols containing the
! charges of the species.
! alias for get_charges
function krome_get_charges()
use krome_subs
use krome_commons
implicit none
real*8::tmp(nspec)
real*8 :: krome_get_charges(nmols)
tmp(:) = get_charges()
krome_get_charges = tmp(1:nmols)
end function krome_get_charges

!*****************
!get an array of character*16 and size krome_nmols
! containing the names of all the species.
! alias for get_names
function krome_get_names() !!  !! cannot yet be called from C
use krome_subs
use krome_commons
implicit none
character*16 :: krome_get_names(nmols)
character*16::tmp(nspec)
tmp(:) = get_names()
krome_get_names = tmp(1:nmols)
end function krome_get_names

!*****************
!get the index of the species with name name.
! alias for get_index
function krome_get_index(name) !! !! cannot yet be called from C
use krome_subs
implicit none
integer :: krome_get_index
character*(*) :: name
krome_get_index = get_index(name)
end function krome_get_index

!*******************
!get the total density of the gas in g/cm3
! giving all the number densities n(:)
function krome_get_rho(n)
use krome_commons
real*8 :: krome_get_rho, n(nmols)
real*8::m(nmols)
m(:) = krome_get_mass()
krome_get_rho = sum(m(:)*n(:))
end function krome_get_rho

!*************************
!scale the abundances of the metals contained in n(:)
! to Z according to Asplund+2009.
! note that this applies only to neutral atoms.
subroutine krome_scale_Z(n,Z)
use krome_commons
real*8 :: n(nmols)
real*8 :: Z
real*8::Htot

Htot =  &
    + n(idx_H) &
    + n(idx_Hj) &
    + n(idx_Hk) &
    + 2d0 * n(idx_H2) &
    + 2d0 * n(idx_H2j)

end subroutine krome_scale_Z

!*************************
!set the total metallicity
! in terms of Z/Z_solar
subroutine krome_set_Z(xarg)
use krome_commons
real*8 :: xarg

total_Z = xarg

end subroutine krome_set_Z

!*************************
!set D is in terms of D_solar (D/D_sol).
subroutine krome_set_dust_to_gas(xarg)
use krome_commons
real*8 :: xarg

dust2gas_ratio_ref = xarg

end subroutine

subroutine krome_get_dust_to_gas(xarg)
use krome_commons
real*8 :: xarg

xarg = dust2gas_ratio_ref

end subroutine krome_get_dust_to_gas

!*************************
!fraction of dust not evaporated
function krome_get_fevap(Tgas)
use krome_commons
use krome_subs, ONLY:fevap

implicit none
real*8 :: Tgas, krome_get_fevap

krome_get_fevap = fevap(Tgas)

end function krome_get_fevap

!*************************
!set the clumping factor
subroutine krome_set_clump(xarg)
use krome_commons
real*8 :: xarg

clump_factor = xarg

end subroutine krome_set_clump

! AP85
! Asano 2014, D(Z/Zsun)/Dsun at 5Gyr (fig 13)
function asano_dust_to_gas(metal_Z)
! Z in  is in term of Z_solar
! D out is in term of D_solar.
use krome_commons
implicit none
real*8 ::metal_Z,asano_dust_to_gas
real*8::x
integer ::i
! spostare in commons, come parametri
real*8::asano_x_low = -1.4
real*8::asano_x_up  = -0.55
real*8,dimension(1:2)::asano_c_low = (/1.21563124, -1.73748871/)
real*8,dimension(1:2)::asano_c_up  = (/1.13122257,  0.02788859/)
real*8,dimension(1:8)::asano_c_mid = (/234.54077964, 1636.69098226, 4771.67185654,&
    & 7504.18772305, 6841.94288515, 3597.86313533, 1009.51484418, 116.35450858/)

asano_dust_to_gas = 0d0

if(metal_Z.gt.0) then
  x = log10(metal_Z)
  if(x.ge.asano_x_up) then
    asano_dust_to_gas = asano_c_up(1)*x + asano_c_up(2)
  elseif(x.le.asano_x_low) then
    asano_dust_to_gas = asano_c_low(1)*x + asano_c_low(2)
  else
    do i=1,8
      asano_dust_to_gas = asano_dust_to_gas*x + asano_c_mid(i)
    enddo
  endif
  asano_dust_to_gas = 10d0**asano_dust_to_gas
endif

end function

!***********************
!get the number of electrons assuming
! total neutral charge (cations-anions)
function krome_get_electrons(x)
use krome_commons
use krome_subs
real*8 :: x(nmols), krome_get_electrons
real*8::n(nspec)
n(1:nmols) = x(:)
n(nmols+1:nspec) = 0d0
krome_get_electrons = get_electrons(n(:))
end function krome_get_electrons

!**********************
!print on screen the first nbest highest reaction fluxes
subroutine krome_print_best_flux(xin,Tgas,nbest)
use krome_subs
use krome_commons
implicit none
real*8 :: xin(nmols)
real*8 :: Tgas
real*8::x(nmols),n(nspec)
integer :: nbest
n(1:nmols) = xin(:)
n(idx_Tgas) = Tgas
call print_best_flux(n,Tgas,nbest)

end subroutine krome_print_best_flux

!*********************
!print only the highest fluxes greater than a fraction frac
! of the maximum flux
subroutine krome_print_best_flux_frac(xin,Tgas,frac)
use krome_subs
use krome_commons
implicit none
real*8 :: xin(nmols)
real*8 :: Tgas,frac
real*8::n(nspec)
n(1:nmols) = xin(:)
n(idx_Tgas) = Tgas
call print_best_flux_frac(n,Tgas,frac)

end subroutine krome_print_best_flux_frac

!**********************
!print the highest nbest fluxes for reactions involving
!a given species using the index idx_find (e.g. krome_idx_H2)
subroutine krome_print_best_flux_spec(xin,Tgas,nbest,idx_find)
use krome_subs
use krome_commons
implicit none
real*8 :: xin(nmols)
real*8 :: Tgas
real*8::n(nspec)
integer :: nbest,idx_find
n(1:nmols) = xin(:)
n(idx_Tgas) = Tgas
call print_best_flux_spec(n,Tgas,nbest,idx_find)
end subroutine krome_print_best_flux_spec

!*******************************
!get an array of size krome_nrea with
! the fluxes of all the reactions in cm-3/s
function krome_get_flux(n,Tgas)
use krome_commons
use krome_subs
real*8 :: krome_get_flux(nrea),n(nmols),Tgas
real*8::x(nspec)
x(:) = 0.d0
x(1:nmols) = n(:)
x(idx_Tgas) = Tgas
krome_get_flux(:) = get_flux(x(:), Tgas)
end function krome_get_flux

!*****************************
!store the fluxes to the file unit ifile
! using the chemical composition x(:), and the
! gas temperature Tgas. xvar is th value of an
! user-defined independent variable that
! can be employed for plots.
! the file columns are as follow
! rate number, xvar, absolute flux,
!  flux/maxflux, flux fraction wrt total,
!  reaction name (*50 string)
subroutine krome_explore_flux(x,Tgas,ifile,xvar)
use krome_commons
use krome_subs
implicit none
real*8 :: x(nmols)
real*8 :: Tgas,xvar
real*8::flux(nrea),fluxmax,sumflux,n(nspec)
integer :: ifile
integer::i
character*50::rname(nrea)

!get reaction names
rname(:) = get_rnames()
n(:) = 0d0
n(1:nmols) = x(:)
n(idx_Tgas) = Tgas
!get fluxes
flux(:) = get_flux(n(:), Tgas)
fluxmax = maxval(flux) !maximum flux
sumflux = sum(flux) !sum of all the fluxes
!loop on reactions
do i=1,nrea
  write(ifile,'(I8,4E17.8e3,a3,a50)') i,xvar,flux(i),&
      flux(i)/fluxmax, flux(i)/sumflux," ",rname(i)
end do
write(ifile,*)

end subroutine krome_explore_flux

!*********************
!get nulcear qeff for the reactions
function krome_get_qeff()
use krome_commons
use krome_subs
implicit none
real*8 :: krome_get_qeff(nrea)

krome_get_qeff(:) = get_qeff()

end function krome_get_qeff

!************************
!dump all the evaluation of the coefficient rates in
! the file funit, in the range inTmin, inTmax, using
! imax points
subroutine krome_dump_rates(inTmin,inTmax,imax,funit)
use krome_commons
use krome_subs
implicit none
integer::i,j
integer :: funit,imax
real*8 :: inTmin,inTmax
real*8::Tmin,Tmax,Tgas,k(nrea),n(nspec)

Tmin = log10(inTmin)
Tmax = log10(inTmax)

n(:) = 1d-40
do i=1,imax
  Tgas = 1d1**((i-1)*(Tmax-Tmin)/(imax-1)+Tmin)
  n(idx_Tgas) = Tgas
  k(:) = coe(n(:))
  do j=1,nrea
    write(funit,'(E17.8e3,I8,E17.8e3)') Tgas,j,k(j)
  end do
  write(funit,*)
end do

end subroutine krome_dump_rates

!************************
!print species informations on screen
subroutine krome_get_info(x, Tgas)
use krome_commons
use krome_subs
integer::i,charges(nspec)
real*8 :: x(nmols)
real*8 :: Tgas
real*8::masses(nspec)
character*16::names(nspec)

names(:) = get_names()
charges(:) = get_charges()
masses(:) = get_mass()

print '(a4,a10,a11,a5,a11)',"#","Name","m (g)","Chrg","x"
do i=1,size(x)
  print '(I4,a10,E11.3,I5,E11.3)',i," "//names(i),masses(i),charges(i),x(i)
end do
print '(a30,E11.3)'," sum",sum(x)

print '(a14,E11.3)',"Tgas",Tgas
end subroutine krome_get_info

end module krome_user

!############### MODULE ##############
module krome_reduction
contains

!**************************
function fex_check(n,Tgas)
use krome_commons
use krome_tabs
implicit none
integer::i
integer::r1,r2
real*8::fex_check,n(nspec),k(nrea),rrmax,Tgas

k(:) = coe_tab(n(:))
rrmax = 0.d0
n(idx_dummy) = 1.d0
n(idx_g) = 1.d0
n(idx_CR) = 1.d0
do i=1,nrea
  r1 = arr_r1(i)
  r2 = arr_r2(i)
  arr_flux(i) = k(i)*n(r1)*n(r2)
  rrmax = max(rrmax, arr_flux(i))
end do
fex_check = rrmax

end function fex_check

end module krome_reduction

!############### MODULE ##############
module krome_main

integer::krome_call_to_fex
!$omp threadprivate(krome_call_to_fex)

contains

! *************************************************************
!  This file has been generated with:
!  krome 14.08.dev "Beastie Boyle" on 2016-08-08 10:48:19.
!
! KROME is a nice and friendly chemistry package for a wide range of
!  astrophysical simulations. Given a chemical network (in CSV format)
!  it automatically generates all the routines needed to solve the kinetic
!  of the system, modelled as system of coupled Ordinary Differential
!  Equations.
!  It provides different options which make it unique and very flexible.
!  Any suggestions and comments are welcomed. KROME is an open-source
!  package, GNU-licensed, and any improvements provided by
!  the users is well accepted. See disclaimer below and GNU License
!  in gpl-3.0.txt.
!
!  more details in http://kromepackage.org/
!  also see https://bitbucket.org/krome/krome_stable
!
! Written and developed by Tommaso Grassi
!  tommasograssi@gmail.com,
!  Starplan Center, Copenhagen.
!  Niels Bohr Institute, Copenhagen.
!
! Co-developer Stefano Bovino
!  sbovino@astro.physik.uni-goettingen.de
!  Institut fuer Astrophysik, Goettingen.
!
! Others (alphabetically): D. Galli, F.A. Gianturco, T. Haugboelle,
!  J.Prieto, J.Ramsey, D.R.G. Schleicher, D. Seifried,
!  E. Simoncini, E. Tognelli
!
!
! KROME is provided "as it is", without any warranty.
!  The Authors assume no liability for any damages of any kind
!  (direct or indirect damages, contractual or non-contractual
!  damages, pecuniary or non-pecuniary damages), directly or
!  indirectly derived or arising from the correct or incorrect
!  usage of KROME, in any possible environment, or arising from
!  the impossibility to use, fully or partially, the software,
!  or any bug or malefunction.
!  Such exclusion of liability expressly includes any damages
!  including the loss of data of any kind (including personal data)
! *************************************************************

!********************************
!KROME main (interface to the solver library)
subroutine krome(x,Tgas,dt )
use krome_commons
use krome_subs
use krome_ode
use krome_reduction
use krome_dust
real*8 :: Tgas,dt
real*8 :: x(nmols)
real*8 :: rhogas
real*8::mass(nspec),n(nspec),tloc,xin
real*8::rrmax,totmass,n_old(nspec),ni(nspec),invTdust(ndust)
integer::icount,i,ierr,icount_max

!DLSODES variables
integer,parameter::meth=2 !1=adam, 2=BDF
integer::neq(1),itol,itask,istate,iopt,lrw,liw,mf
integer::iwork(131)
real*8::atol(nspec),rtol(nspec)
real*8::rwork(520)
logical::got_error,equil

!****************************
!init DLSODES (see DLSODES manual)
call XSETF(0)!toggle solver verbosity
got_error = .false.
neq = nspec !number of eqns
liw = size(iwork)
lrw = size(rwork)
iwork(:) = 0
rwork(:) = 0.d0
itol = 4 !both tolerances are scalar
rtol(:) = 1.000000d-04 !relative tolerance
atol(:) = 1.000000d-10 !absolute tolerance
icount_max = 100 !maximum number of iterations

itask = 1
iopt = 0

!MF=
!  = 222 internal-generated JAC and sparsity
!  = 121 user-provided JAC and internal generated sparsity
!  =  22 internal-generated JAC but sparsity user-provided
!  =  21 user-provided JAC and sparsity
MF = 222
!end init DLSODES
!****************************

ierr = 0 !error flag, zero==OK!
n(:) = 0.d0 !initialize densities

n(1:nmols) = x(:)

n(idx_Tgas) = Tgas !put temperature in the input array

icount = 0 !count solver iterations
istate = 1 !init solver state
tloc = 0.d0 !set starting time

!store initial values
ni(:) = n(:)

n_old(:) = -1d99
krome_call_to_fex = 0
do
  icount = icount + 1
  !solve ODE
  CALL DLSODES(fex, NEQ(:), n(:), tloc, dt, ITOL, RTOL, ATOL,&
      ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JES, MF)
  krome_call_to_fex = krome_call_to_fex + IWORK(12)
  !check DLSODES exit status
  if(istate==2) then
    exit !sucsessful integration
  elseif(istate==-1) then
    istate = 1 !exceeded internal max iterations
  elseif(istate==-5 .or. istate==-4) then
    istate = 3 !wrong sparsity recompute
  elseif(istate==-3) then
    n(:) = ni(:)
    istate = 1
  else
    got_error = .true.
  end if

  if(got_error.or.icount>icount_max) then
    if (krome_mpi_rank>0) then
      print *,krome_mpi_rank,"ERROR: wrong solver exit status!"
      print *,krome_mpi_rank,"istate:",istate
      print *,krome_mpi_rank,"iter count:",icount
      print *,krome_mpi_rank,"max iter count:",icount_max
      print *,krome_mpi_rank,"SEE KROME_ERROR_REPORT file"
    else
      print *,"ERROR: wrong solver exit status!"
      print *,"istate:",istate
      print *,"iter count:",icount
      print *,"max iter count:",icount_max
      print *,"SEE KROME_ERROR_REPORT file"
    end if
    call krome_dump(n(:), rwork(:), iwork(:), ni(:))
    stop
  end if

end do

!avoid negative species
do i=1,nspec
  n(i) = max(n(i),1d-40)
end do

n(:) = conserve(n(:),ni(:))

!returns to user array
x(:) = n(1:nmols)

Tgas = n(idx_Tgas) !get new temperature

end subroutine krome

!*********************************
!integrates to equilibrium using constant temperature
subroutine krome_equilibrium(x,Tgas)
use krome_ode
use krome_subs
use krome_commons
use krome_constants
implicit none
integer::mf,liw,lrw,itol,meth,iopt,itask,istate,neq(1)
integer::i,imax
real*8 :: Tgas
real*8 :: x(nmols)
real*8 :: rhogas
real*8::tloc,n(nspec),mass(nspec),ni(nspec)
real*8::dt,xin
integer::iwork(131)
real*8::atol(nspec),rtol(nspec)
real*8::rwork(520)
real*8::ertol,eatol,max_time,t_tot
logical::converged

integer, save :: ncall=0
integer, parameter :: ncall_print_frequency=20000
integer :: ncallp
integer::charges(nspec)
real*8::masses(nspec)
character*16::names(nspec)

call XSETF(0)!toggle solver verbosity
meth = 2
neq = nspec !number of eqns
liw = size(iwork)
lrw = size(rwork)
iwork(:) = 0
rwork(:) = 0.d0
itol = 4 !both tolerances are scalar
rtol(:) = 1d-6 !relative tolerance
atol(:) = 1d-20 !absolute tolerance

! Switches to decide when equilibrium has been reached
ertol = 1d-5  ! relative min change in a species
eatol = 1d-12 ! absolute min change in a species
max_time=seconds_per_year*5d8 ! max time we will be integrating for

!for DLSODES options see its manual
iopt = 0
itask = 1
istate = 1

mf = 222 !internally evaluated sparsity and jacobian
tloc = 0d0 !initial time

n(:) = 0.d0 !initialize densities
!copy into array
n(nmols+1:) = 0d0
n(1:nmols) = x(:)

n(idx_Tgas) = Tgas

!store previous values
ni(:) = n(:)

imax = 1000

dt = seconds_per_year * 100.
t_tot = dt
converged = .false.
do while (.not. converged)
  do i=1,imax
    !solve ODE
    CALL DLSODES(fcn_tconst, NEQ(:), n(:), tloc, dt, ITOL, RTOL, ATOL,&
        ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, jcn_dummy, MF)
    if(istate==2) then
      exit
    else
      istate=1
    end if
  end do
  !check errors
  if(istate.ne.2) then
    print *,"ERROR: no equilibrium found!"
    stop
  end if

  !avoid negative species
  do i=1,nspec
    n(i) = max(n(i),0.d0)
  end do

  n(:) = conserve(n(:),ni(:))
  ! check if we have converged by comparing the error in any species with an relative abundance above eatol
  converged = maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) .lt. ertol &
      .or. t_tot .gt. max_time

  ! Increase integration time by a reasonable factor
  if (.not. converged) then
    dt = dt * 3.
    t_tot = t_tot + dt
    ni = n
  endif
enddo
!returns to user array
x(:) = n(1:nmols)

if (t_tot > max_time .and. &
    maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) > 0.2) then
print *, 'krome_equilibrium: Did not converge in ', max_time / seconds_per_year, ' years.'
print *, 'Tgas :', Tgas
names(:) = get_names()
charges(:) = get_charges()
masses(:) = get_mass()

print '(a4,a10,a11,a5,a16)',"#","Name","m (g)","Chrg","  Current / Last"
do i=1,nmols
  print '(I4,a10,E11.3,I5,2E14.6,E11.3)',i," "//names(i),masses(i),charges(i),n(i),ni(i),abs(n(i) - ni(i)) / max(n(i),eatol*sum(n(1:nmols)))
end do
print '(a30,2E14.6)'," sum",sum(n(1:nmols)),sum(ni(1:nmols))
print *, 'Fractional error :', maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols))))
print *, 'Absolute and relative floors:', eatol, ertol
endif

! Print info ever so often
!$omp critical
ncall=ncall+1
ncallp = ncall
!$omp end critical

if (modulo(ncallp,ncall_print_frequency)==0) then
print *, 'Found equilibrium for ', ncallp, ' cells.'
endif

end subroutine krome_equilibrium

!********************
!dummy jacobian
subroutine jcn_dummy()
implicit none
end subroutine jcn_dummy

!*******************
!dn/dt where dT/dt=0
subroutine fcn_tconst(n,tt,x,f)
use krome_commons
use krome_ode
implicit none
integer::n,ierr
real*8::x(n),f(n),tt
call fex(n,tt,x(:),f(:))
f(idx_Tgas) = 0d0
end subroutine fcn_tconst

!*******************************
subroutine krome_dump(n,rwork,iwork,ni)
use krome_commons
use krome_subs
use krome_tabs
use krome_reduction
use krome_ode
integer::fnum,i,iwork(:),idx(nrea),j
real*8::n(:),rwork(:),rrmax,k(nrea),kmax,rperc,kperc,dn(nspec),tt,ni(:)
character*16::names(nspec),FMTi,FMTr
character*50::rnames(nrea),fname,prex
integer,save::mx_dump=1000 ! max nr of reports before terminating

fnum = 99+krome_mpi_rank
if (krome_mpi_rank>0) then
write(fname,'(a,i5.5)') "KROME_ERROR_REPORT_",krome_mpi_rank
else
fname = "KROME_ERROR_REPORT"
endif
open(fnum,FILE=trim(fname),status="replace")
tt = 0d0
names(:) = get_names()
rnames(:) = get_rnames()
call fex(nspec,tt,n(:),dn(:))

write(fnum,*) "KROME ERROR REPORT"
write(fnum,*)
!SPECIES
write(fnum,*) "Species abundances"
write(fnum,*) "**********************"
write(fnum,'(a5,a20,3a12)') "#","name","qty","dn/dt","ninit"
write(fnum,*) "**********************"
do i=1,nspec
write(fnum,'(I5,a20,3E12.3e3)') i,names(i),n(i),dn(i),ni(i)
end do
write(fnum,*) "**********************"

!F90 FRIENDLY RESTART
write(fnum,*)
write(fnum,*) "**********************"
write(fnum,*) "F90-friendly species"
write(fnum,*) "**********************"
do i=1,nspec
write(prex,'(a,i3,a)') "x(",i,") = "
write(fnum,*) trim(prex),ni(i),"!"//names(i)
end do

write(fnum,*) "**********************"

!RATE COEFFIECIENTS
k(:) = coe_tab(n(:))
idx(:) = idx_sort(k(:))
kmax = maxval(k)
write(fnum,*)
write(fnum,*) "Rate coefficients (sorted) at Tgas",n(idx_Tgas)
write(fnum,*) "**********************"
write(fnum,'(a5,2a12,a10)') "#","k","k %","  name"
write(fnum,*) "**********************"
do j=1,nrea
i = idx(j)
kperc = 0.d0
if(kmax>0.d0) kperc = k(i)*1d2/kmax
write(fnum,'(I5,2E12.3e3,a2,a50)') i,k(i),kperc,"  ", rnames(i)
end do
write(fnum,*) "**********************"
write(fnum,*)

!FLUXES
call load_arrays
rrmax = fex_check(n(:), n(idx_Tgas))
idx(:) = idx_sort(arr_flux(:))
write(fnum,*)
write(fnum,*) "Reaction magnitude (sorted) [k*n1*n2*n3*...]"
write(fnum,*) "**********************"
write(fnum,'(a5,2a12,a10)') "#","flux","flux %","  name"
write(fnum,*) "**********************"
do j=1,nrea
i = idx(j)
rperc = 0.d0
if(rrmax>0.d0) rperc = arr_flux(i)*1d2/rrmax
write(fnum,'(I5,2E12.3e3,a2,a50)') i,arr_flux(i),rperc,"  ",rnames(i)
end do
write(fnum,*) "**********************"
write(fnum,*)

!SOLVER
FMTr = "(a30,E16.7e3)"
FMTi = "(a30,I10)"
write(fnum,*) "Solver-related information:"
write(fnum,FMTr) "step size last",rwork(11)
write(fnum,FMTr) "step size attempt",rwork(12)
write(fnum,FMTr) "time current",rwork(13)
write(fnum,FMTr) "tol scale factor",rwork(14)
write(fnum,FMTi) "numeber of steps",iwork(11)
write(fnum,FMTi) "call to fex",iwork(12)
write(fnum,FMTi) "call to jex",iwork(13)
write(fnum,FMTi) "last order used",iwork(14)
write(fnum,FMTi) "order attempt",iwork(15)
write(fnum,FMTi) "idx largest error",iwork(16)
write(fnum,FMTi) "RWORK size required",iwork(17)
write(fnum,FMTi) "IWORK size required",iwork(18)
write(fnum,FMTi) "NNZ in Jac",iwork(19)
write(fnum,FMTi) "extra fex to compute jac",iwork(20)
write(fnum,FMTi) "number of LU decomp",iwork(21)
write(fnum,FMTi) "base address in RWORK",iwork(22)
write(fnum,FMTi) "base address of IAN",iwork(23)
write(fnum,FMTi) "base address of JAN",iwork(24)
write(fnum,FMTi) "NNZ in lower LU",iwork(25)
write(fnum,FMTi) "NNZ in upper LU",iwork(21)
write(fnum,*) "See DLSODES manual for further details on Optional Outputs"
write(fnum,*)
write(fnum,*) "END KROME ERROR REPORT"
write(fnum,*)
close(fnum)

mx_dump = mx_dump - 1
if (mx_dump==0) stop

end subroutine krome_dump

!********************************
subroutine krome_init()
use krome_commons
use krome_tabs
use krome_subs
use krome_reduction
use krome_dust
use krome_cooling
use krome_user

real*8,dimension(10)::flux_for_krome_init

!init phys common variables
!$omp parallel
phys_Tcmb = 2.73d0
phys_zredshift = 0d0
phys_orthoParaRatio = 3d0
phys_metallicity = 0d0
phys_Tfloor = 2.73d0
!$omp end parallel

!krome_datafolder=trim("./data_folder/")

!init metallicity default
!assuming solar
total_Z = 1d0

!default D/D_sol = Z/Z_sol
!assuming linear scaling
dust2gas_ratio_ref = total_Z

!default clumping factor for
! H2 formation on dust by Jura/Gnedin
clump_factor = 1d0

!default for thermo toggle is ON
!$omp parallel
krome_thermo_toggle = 1
!$omp end parallel

call load_arrays

flux_for_krome_init = (/0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 /)
call krome_set_photoBinJ(flux_for_krome_init)
call krome_set_user_myfluxLW(0d0)
call krome_set_user_myH2_dissociation(2.07d9)
call krome_set_user_crate(1.0d-17)

call make_ktab()
call check_tabs()

!initialize metal CIE cooling
!call init_coolingZCIE()

!initialize metal CIE cooling no UV case
call init_anytab2D(trim(krome_datafolder)//"coolZ_CIE2012NOUV.dat",CoolZNOUV_x(:), &
      CoolZNOUV_y(:), CoolZNOUV_z(:,:), CoolZNOUV_xmul, &
      CoolZNOUV_ymul)
!call test_anytab2D(trim(krome_datafolder)//"coolZ_CIE2012NOUV.dat",CoolZNOUV_x(:), &
!      CoolZNOUV_y(:), CoolZNOUV_z(:,:), CoolZNOUV_xmul, &
!      CoolZNOUV_ymul)

!initialize the table for exp(-a/T) function
call init_exp_table()

call load_parts()

!init photo reactants indexes
photoPartners(1) = idx_H
photoPartners(2) = idx_HE
photoPartners(3) = idx_HEj
photoPartners(4) = idx_Hk
photoPartners(5) = idx_H2
photoPartners(6) = idx_H2j
photoPartners(7) = idx_H2j
photoPartners(8) = idx_H2

!get machine precision
krome_epsilon = epsilon(0d0)

end subroutine krome_init

!****************************
function krome_get_coe(x,Tgas)
!krome_get_coe: public interface to obtain rate coefficients
use krome_commons
use krome_subs
use krome_tabs
implicit none
real*8 :: krome_get_coe(nrea), x(nmols), Tgas
real*8::n(nspec)

n(:) = 0d0
n(1:nmols) = x(:)
n(idx_Tgas) = Tgas
krome_get_coe(:) = coe_tab(n(:))

end function krome_get_coe

!****************************
function krome_get_coeT(Tgas)
!krome_get_coeT: public interface to obtain rate coefficients
! with argument Tgas only
use krome_commons
use krome_subs
use krome_tabs
implicit none
real*8 :: krome_get_coeT(nrea),Tgas
real*8::n(nspec)
n(idx_Tgas) = Tgas
krome_get_coeT(:) = coe_tab(n(:))
end function krome_get_coeT

end module krome_main
