!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  !-------------------------------------------
  ! This routine builds the initial AMR grid
  !-------------------------------------------
  integer::ilevel

  if(myid==1)write(*,*)'Building initial AMR grid'
  init=.true.

  ! Base refinement
  do ilevel=1,levelmin
     call flag
     call refine
  end do

  ! Further refinements if necessary
  do ilevel=levelmin+1,nlevelmax
     if(initfile(levelmin).ne.' '.and.initfile(ilevel).eq.' ')exit
     if(hydro)call init_flow
     first_time = .false.
#ifdef RT
     if(rt)call rt_init_flow
#endif
     if(ivar_refine==0)call init_refmap
     call flag
     call refine
     if(nremap>0)call load_balance
     if(numbtot(1,ilevel)==0)exit
  end do

  ! Final pass to initialize the flow
  init=.false.
  if(hydro)call init_flow
#ifdef RT
  if(rt)call rt_init_flow
#endif

end subroutine init_refine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_2
  !--------------------------------------------------------------
  ! This routine builds additional refinements to the
  ! the initial AMR grid for filetype ne 'grafic'
  !--------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use dice_commons
#ifdef RT
  use rt_hydro_commons
#endif
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel,i,j, ivar
  real(dp)::eps_star2
  integer:: i_conv
  
  if(filetype.eq.'grafic')return

  do i=levelmin,nlevelmax+1
     ! DICE------
     do ilevel=levelmin-1,1,-1
        if(pic)call merge_tree_fine(ilevel)
     enddo
     ! ----------
     call refine_coarse

     do ilevel=1,nlevelmax
        call build_comm(ilevel)
        call make_virtual_fine_int(cpu_map(1),ilevel)
        call refine_fine(ilevel)
        ! DICE------
        if(pic)call make_tree_fine(ilevel)
        ! ----------
        if(hydro)call init_flow_fine(ilevel)
        ! DICE------
        if(pic)then
           call kill_tree_fine(ilevel)
           call virtual_tree_fine(ilevel)
        endif
        ! ----------
#ifdef RT
        if(rt)call rt_init_flow_fine(ilevel)
#endif
     end do

     ! DICE------
     do ilevel=nlevelmax-1,levelmin,-1
        if(pic)call merge_tree_fine(ilevel)
     enddo
    ! ----------

    if(nremap>0)call load_balance

    do ilevel=levelmin,nlevelmax
        if(pic)call make_tree_fine(ilevel)
        if(poisson)call rho_fine(ilevel,2)
        if(pic)then
           call kill_tree_fine(ilevel)
           call virtual_tree_fine(ilevel)
        endif
     end do

     do ilevel=nlevelmax,levelmin,-1
        if(pic)call merge_tree_fine(ilevel)
        if(hydro)then
           call upload_fine(ilevel)
#ifdef SOLVERmhd
           do ivar=1,nvar+3
#else
           do ivar=1,nvar
#endif
              call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
           end do
#else
           end do
#endif
           if(simple_boundary)call make_boundary_hydro(ilevel)
        endif
#ifdef RT
        if(rt)then
           call rt_upload_fine(ilevel)
           do ivar=1,nrtvar
              call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
           end do
           if(simple_boundary)call rt_make_boundary_hydro(ilevel)
        end if
#endif
     end do

     do ilevel=nlevelmax,1,-1
        call flag_fine(ilevel,2)
     end do
     call flag_coarse

  end do

  ! DICE------
  do ilevel=levelmin-1,1,-1
    if(pic)call merge_tree_fine(ilevel)
  enddo
 call kill_gas_part(1)

  do ilevel=1,nlevelmax
     if(pic)then
        call make_tree_fine(ilevel)
        call kill_tree_fine(ilevel)
        call virtual_tree_fine(ilevel)
     endif
  end do
  do ilevel=nlevelmax,levelmin,-1
     call merge_tree_fine(ilevel)
  end do
!  deallocate(uthp)
 if(sf_virial)then
     eps_star2=eps_star
     eps_star=0d0
     do ilevel=nlevelmax,levelmin,-1
        call star_formation(ilevel)
     enddo
     eps_star=eps_star2
  endif
  dice_init=.false.
  ! ----------


!### DICE: KROME EQUILIBRIUM STATES 


do i_conv = 1, n_krome_init
  if (krome_init_equilibrium)then
   do ilevel=nlevelmax,1,-1
     if(myid==1) print*, 'running krome on IC (',i_conv,'/', n_krome_init, ') at lev', ilevel
     call dice_init_chem(ilevel, dt_krome_init)
   end do
  endif
enddo

do ilevel=nlevelmax,1,-1
     call upload_fine(ilevel)
end do



!#if defined(RT)
!  if(rt_is_init_xion) then
!     if(myid==1) write(*,*) 'Initializing ionization states from T profile'
!     do ilevel=nlevelmax,1,-1
!        call rt_init_xion(ilevel)
!        call upload_fine(ilevel)
!     end do
!  endif
!#endif

end subroutine init_refine_2
!################################################################
!################################################################
!################################################################
!################################################################

!################################################################
subroutine kill_gas_part(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  include 'mpif.h'
  integer::ilevel
  !--------------------------------------------------------
  ! This subroutine removes the gas particles
  ! initially present in the gadget1 DICE output
  !--------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu,info
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  logical,dimension(1:nvector)::ok=.true.
  integer::npart_all
  integer,dimension(1:ncpu)::npart_cpu,npart_cpu_all

  npart_cpu = 0
  npart_all = 0

  if(numbtot(1,ilevel)==0)return
  ! Gather gas particles.
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        ! Count gas particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).eq.1)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
           npart_cpu(myid)=npart_cpu(myid)+npart2
        endif
        ! Gather gas particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only gas particles
              if(idp(ipart).eq.1)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 call remove_list(ind_part,ind_grid_part,ok,ip)
                 call add_free_cond(ind_part,ok,ip)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)then
        call remove_list(ind_part,ind_grid_part,ok,ip)
        call add_free_cond(ind_part,ok,ip)
     end if
  end do

#ifndef WITHOUTMPI
  ! Give an array of number of gas on each cpu available to all cpus
  call MPI_ALLREDUCE(npart_cpu,npart_cpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  npart_all=sum(npart_cpu_all(1:ncpu))
  if(npart_all>0) then
     if(myid==1) then
        write(*,'(A50)')"__________________________________________________"
        write(*,'(A,I15)')' Gas particles deleted ->',npart_all
        write(*,'(A50)')"__________________________________________________"
     endif
  endif
npart_cpu(myid)=npart
#ifndef WITHOUTMPI
  ! Give an array of number of gas on each cpu available to all cpus
  call MPI_ALLREDUCE(npart_cpu,npart_cpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  npart_all=sum(npart_cpu_all(1:ncpu))
  if(npart_all>0) then
     if(myid==1) then
        write(*,'(A50)')"__________________________________________________"
        write(*,'(A,I15)')' Remaining particles ->',npart_all
        write(*,'(A50)')"__________________________________________________"
     endif
  endif

  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        ! Count gas particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).gt.0) idp(ipart)=idp(ipart)-1
              ipart=next_part  ! Go to next particle
           end do
           npart_cpu(myid)=npart_cpu(myid)+npart2
        endif
     end do
  end do


111 format('   Entering kill_gas_part for level ',I2)
!---------------------------------------------
end subroutine



SUBROUTINE dice_init_chem(ilevel, dt_krome_init)

!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer:: ilevel
  integer:: ncache,i,igrid,ngrid
  real(dp):: dt_krome_init
  integer,dimension(1:nvector),save:: ind_grid
!-------------------------------------------------------------------------
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! Do the initialization by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call dice_init_chem_vsweep(ind_grid, ngrid, ilevel, dt_krome_init)
  end do

111 format('   Entering dice_init_chem for level',i2)

END SUBROUTINE dice_init_chem


!#####################

SUBROUTINE dice_init_chem_vsweep(ind_grid, ngrid, ilevel, dt_krome_init)

  use amr_commons
  use hydro_commons
  use rt_parameters
  use krome_user
  use krome_photo
  use krome_main !mandatory
  use krome_commons, ONLY:photoPartners,photoBinJTab
  use krome_subs, ONLY: calc_H2shieldR14_for_ramses

  implicit none
  include 'mpif.h'
  
  integer::ngrid, ilevel
  integer,dimension(1:nvector)::ind_grid
  integer::i, ind, iskip, idim, nleaf,nx_loc
  integer::i_photo_rea,idx, i_bin
  real(dp)::H2shield
  real(dp)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  integer,dimension(1:nvector),save::ind_cell, ind_leaf
  real(dp)::density, T2, ekk, err, emag, x, mu
  real(dp),dimension(krome_nmols)::unoneq   
  real(kind=8)::Z_for_krome,dust_for_krome,size_loc_cm,energy_norm
  real(dp),dimension(1:krome_nmols)::krome_scale,inv_krome_scale
  real(dp)::Zsolar,mu_noneq_old
  integer:: i_krome
  integer,dimension(1:krome_nmols) ::krome_id_vec
  real(dp),dimension(krome_nPhotoBins)::flux_for_krome
  real(dp)::scale_np,scale_fp
  real(dp)::dx,dx_loc,scale
  real(dp),dimension(nGroups),save::phAbs
  real(dp)::dt_krome_init
#if NENER>0
  integer::irad
#endif


!-------------------------------------------------------------------------

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  !skip_loc=(/0.0d0,0.0d0,0.0d0/)
  !if(ndim>0)skip_loc(1)=dble(icoarse_min)
  !if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  !if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call krome_units(scale_d,krome_scale,inv_krome_scale)
  call rt_units(scale_np, scale_fp)
  size_loc_cm   = dx_loc*scale_l

  ! id shortcuts
  do i_krome =1,krome_nmols
    krome_id_vec(i_krome)=krome_off+i_krome
  enddo

  call krome_set_zredshift(1.d0/aexp -1.d0)
  call krome_set_Tcmb(10d0)
  call krome_set_Tfloor(10d0)
  call krome_set_user_cell_size(size_loc_cm)
  call krome_set_user_crate(CR_ionisation_rate)
  call krome_set_user_myH2_dissociation(H2_dissociation_rate)

  call krome_set_zredshift(0d0)

! Loop over cells in each oct
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do
     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do
     if(nleaf .eq. 0) cycle

     do i=1,nleaf

        ! Compute rho
        density = MAX(uold(ind_leaf(i),1),smallr)   !       Mass density of gas
        ! Compute pressure from energy density
        T2 = uold(ind_leaf(i),ndim+2)          ! Energy density (kin+heat)
        ekk = 0.0d0                            !            Kinetic energy

        do idim=1,ndim
           ekk = ekk+0.5*uold(ind_leaf(i),idim+1)**2/density
        end do
        err = 0.0d0
#if NENER>0
        do irad=0,nener-1
           err = err+uold(ind_leaf(i),inener+irad)
        end do
#endif
        emag = 0.0d0
#ifdef SOLVERmhd
        do idim=1,ndim
           emag=emag+0.125d0*(uold(ind_leaf(i),idim+ndim+2)+uold(ind_leaf(i),idim+nvar))**2
        end do
#endif
        T2 = (gamma-1.0)*(T2-ekk-err-emag)     !     Gamma is ad. exponent

        ! now T2 is pressure (in user units)   !    (relates p and energy)
        ! Compute T2=T/mu in Kelvin from pressure
        T2 = scale_T2*T2/density

        ! set metallicity
        if(metal) then
          Zsolar = (uold(ind_leaf(i),imetal)/density)/0.02
        else
          Zsolar = z_ave
        endif

        call krome_set_Z(Zsolar)
        call krome_set_dust_to_gas(Zsolar)

        ! set cosmological initial abundances
        ! \rho_X [units dens]
        do i_krome =1,krome_nmols
          unoneq(i_krome)     =  uold(ind_leaf(i),krome_id_vec(i_krome))
        end do 
        ! n_X from code units of RAMSES to KROME [cm-3]
        unoneq(1:krome_nmols) = unoneq(1:krome_nmols)*krome_scale(1:krome_nmols)
        ! T/mu [K] -> T [K]
        mu_noneq_old          = krome_get_mu(unoneq(1:krome_nmols))
        T2                    = T2 * mu_noneq_old

  phAbs(1:nGroups)=0.

     do i_bin=1,krome_nPhotoBins
       phAbs(i_bin) = 0
       do i_photo_rea=1,krome_nPhotoRates
         idx      = photoPartners(i_photo_rea)
         ! density * cross_section
         phAbs(i_bin) = phAbs(i_bin) + unoneq(idx) * photoBinJTab(i_photo_rea,i_bin)
       enddo
       ! rate [s-1]   !should be no H2 in initial conditions
       if(groupL0(i_bin) .ge. 11.2d0 .and.  groupL1(i_bin) .le. 13.6d0) then
         ! f_shield = exp(-sigma n dx)
         H2shield = max(calc_H2shieldR14_for_ramses(unoneq(KROME_idx_H2), T2, krome_get_user_cell_size()), tiny(1.d0))
         phAbs(i_bin) = phAbs(i_bin) - log(H2shield) / krome_get_user_cell_size()
       endif
       phAbs(i_bin) = phAbs(i_bin) * rt_c_cgs
     enddo

  flux_for_krome(1:krome_nPhotoBins) = 0d0

  !!if UV background, add it here to krome calculation
  if (rt_uv_background)then
      flux_for_krome(1:krome_nPhotoBins) = flux_for_krome(1:krome_nPhotoBins) &
                     + rt_n_background(1:krome_nPhotoBins)/rt_c_cgs/scale_Np  &
                     * EXP(-phAbs(1:krome_nPhotoBins) / rt_c_cgs * krome_get_user_cell_size())
  endif

   


  ! eV/cm3
  flux_for_krome(:) = flux_for_krome(:)*krome_get_photoBinE_mid()
  ! ev/cm2/s
  flux_for_krome(:) = flux_for_krome(:)*rt_c_cgs
  ! ev/cm2/s/Hz
  flux_for_krome(:) = flux_for_krome(:)/(krome_iplanck_eV*krome_get_photoBinE_delta())
  ! ev/cm2/s/Hz/sr
  flux_for_krome(:) = flux_for_krome(:)/(4*krome_pi)

  ! input is in eV/s/cm2/sr/Hz
  call krome_set_photoBinJ(flux_for_krome(:))
  call krome_set_user_myfluxLW(get_photoIntensity(12.87d0))

        ! compute chemical equlibrium keeping the temperature fixed
        ! even tough -- as Bovinsky uses say -- chemical equilibrium does not
        ! exists
        ! note that convergence by krome is defined if the relative error in ANY
        ! species is below a set accuracy
        ! this is typically an overkill requirement for IC, particular when
        ! idealized and considering that T does not change

        !call krome_equilibrium(unoneq,T2)
        call krome_Tconst(unoneq, T2, dt_krome_init)

  if (T2 < 10d0) then
    T2 = 10d0
  endif
  if (T2 > 1e7) then
    T2 = 1e7
  endif

        ! [cm-3] -> [unit dens]
        unoneq(:) = unoneq(:)*inv_krome_scale(:)

        do i_krome =1,krome_nmols
           uold(ind_leaf(i),krome_id_vec(i_krome)) = unoneq(i_krome)
        end do
      end do

  end do
  ! End loop over cells
END SUBROUTINE dice_init_chem_vsweep



