!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flow
  use amr_commons
  use hydro_commons, ONLY: nvar, uold
  implicit none

  integer::ilevel,ivar

  if(verbose)write(*,*)'Entering init_flow'
  do ilevel=nlevelmax,1,-1
     if(ilevel>=levelmin)call init_flow_fine(ilevel)
     call upload_fine(ilevel)
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     end do
     if(simple_boundary)call make_boundary_hydro(ilevel)
  end do
  if(verbose)write(*,*)'Complete init_flow'

end subroutine init_flow

subroutine init_krome_and85

  use amr_commons
  !KROME: add the main module to call krome_init
  use krome_main
#if defined(avg_field)
  use krome_user, ONLY: krome_set_clump,krome_set_photoBin_draineLog, krome_photoBin_store
#elif defined(RT)
  use krome_user, ONLY: krome_set_clump,krome_set_photobinE_lr,krome_nPhotoBins
  use rt_parameters, ONLY: groupL0,groupL1,group_egy,nGroups
#else
  use krome_user, ONLY: krome_set_clump
#endif
  use krome_commons, ONLY: krome_mpi_rank,krome_nfile,krome_nfile2

  implicit none

#if defined(RT)
  real*8,dimension(nGroups):: photo_right_for_krome
  integer :: i
#endif

  krome_mpi_rank = myid
  krome_nfile    = ncpu   + 22 + myid
  krome_nfile2   = ncpu*2 + 42 + myid

  !KROME: initialize KROME. This call is mandatory.
  call krome_init()

  if(myid.eq.1) then
    write(*,*) "krome init done, performing additional settings"
  endif

  !avg_sfr    = 0.d0
  krome_chem = .true.

#if defined(RT) && defined(avg_field)
  print* 'Incompatible option in the makefile'
  print* 'aborting'
  call clean_stop
#endif


#if defined(RT)
  if(krome_nPhotoBins .ne. nGroups) then
    write(*,*) 'n of photobins from the namelist    ',nGroups
    write(*,*) 'n of photobins for krome compilation',krome_nPhotoBins
    write(*,*) 'Aborting'
    call clean_stop
  endif

  photo_right_for_krome(1:krome_nPhotoBins)=groupL1(1:krome_nPhotoBins)
  ! ramses_rt takes 0 as infty ...
  if(photo_right_for_krome(krome_nPhotoBins) .eq.0) then
    photo_right_for_krome(krome_nPhotoBins) = 2*group_egy(krome_nPhotoBins)-groupL0(krome_nPhotoBins)
  endif

  if(myid.eq.1) then
    write(*,*) 'Setting krome radiation in ',krome_nPhotoBins,' bins:'
    write(*,*) '  nbin | E_left [eV] | E_right [eV]'
    do i=1,krome_nPhotoBins
      write(*,*)'  ',i,groupL0(i),photo_right_for_krome(i)
    enddo
  endif

  call krome_set_photobinE_lr(groupL0(1:krome_nPhotoBins),photo_right_for_krome(1:krome_nPhotoBins))
#endif

end subroutine init_krome_and85
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flow_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
#if defined(krome_on)
  use amr_parameters, ONLY: krome_off
  use krome_user, ONLY: krome_nmols,KROME_idx_H,KROME_idx_E,KROME_idx_Hj,&
    &                   KROME_idx_HE,KROME_idx_HEj,KROME_idx_HEjj,KROME_idx_Hk,&
    &                   KROME_idx_H2,KROME_idx_H2j
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info,info2,dummy_io
#endif
  integer::ilevel

  integer::i,icell,igrid,ncache,iskip,ngrid,ilun
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count
  integer ,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,rr,vx,vy,vz,ek,ei,pp,xx1,xx2,xx3,dx_loc,scale,xval
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector)       ,save::vv
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar),save::uu

  real(dp),allocatable,dimension(:,:,:)::init_array
  real(kind=4),allocatable,dimension(:,:)  ::init_plane

  logical::error,ok_file1,ok_file2,ok_file3,ok_file
  character(LEN=80)::filename, foldername
  character(LEN=80),dimension(NVAR)::varnames
  character(LEN=5)::nchar,ncharvar, fmt_lev
  character(LEN=2)::str_lev
  real(dp)::xmin,xmax,ymin,ymax,zmin,zmax

#if defined(krome_on)
  integer::i_krome
  character*4,dimension(1:krome_nmols)::  krome_names_arr
  real(kind=8),dimension(1:krome_nmols):: krome_init_var,krome_check_init
  real(kind=8)::krome_init_temp,krome_init_temp_check
#endif
  integer,parameter::tag=1107

#if NDIM==1
    real(dp),dimension(2**file_ncell(ilevel))::file_var
#endif
#if NDIM==2
    real(dp),dimension(2**file_ncell(ilevel),2**file_ncell(ilevel))::file_var
#endif
#if NDIM==3
    !real,dimension(file_ncell(ilevel),file_ncell(ilevel),file_ncell(ilevel))::file_var
    real, allocatable :: file_var(:,:,:)
#endif
  

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

#if defined(krome_on)

  if(krome_nmols .ne.9) then
    write(*,*)'Mismatch in init file'
    write(*,*)'aborting'
    call clean_stop
  endif

  ! galli & palla 1998 @ z = 100
  krome_init_var(KROME_idx_H   ) = 0.752471d+00 !H
  krome_init_var(KROME_idx_E   ) = 0.608748d-07 !E
  krome_init_var(KROME_idx_Hj  ) = 0.111772d-03 !H+
  krome_init_var(KROME_idx_HE  ) = 0.247527d+00 !HE
  krome_init_var(KROME_idx_HEj ) = 0.287201d-26 !HE+
  krome_init_var(KROME_idx_HEjj) = 0.334857d-52 !HE++
  krome_init_var(KROME_idx_Hk  ) = 0.626462d-11 !H-
  krome_init_var(KROME_idx_H2  ) = 0.559307d-06 !H2
  krome_init_var(KROME_idx_H2j ) = 0.728777d-13 !H2+

  krome_init_temp   = 0.120995d+03

  krome_names_arr(KROME_idx_H   ) = 'H   '
  krome_names_arr(KROME_idx_E   ) = 'E   '
  krome_names_arr(KROME_idx_Hj  ) = 'H+  '
  krome_names_arr(KROME_idx_HE  ) = 'HE  '
  krome_names_arr(KROME_idx_HEj ) = 'HE+ '
  krome_names_arr(KROME_idx_HEjj) = 'HE++'
  krome_names_arr(KROME_idx_Hk  ) = 'H-  '
  krome_names_arr(KROME_idx_H2  ) = 'H2  '
  krome_names_arr(KROME_idx_H2j ) = 'H2+ '
#endif

  ! Conversion factor from user units to cgs units

  allocate(file_var(file_ncell(ilevel),file_ncell(ilevel),file_ncell(ilevel)))

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  ncache=active(ilevel)%ngrid

  !--------------------------------------
  ! Compute initial conditions from files
  !--------------------------------------

  !!! COSMO
  if(cosmo)then

   filename=TRIM(initfile(ilevel))//'/ic_d'
   INQUIRE(file=filename,exist=ok_file1)
   if(multiple)then
     filename=TRIM(initfile(ilevel))//'/dir_deltab/ic_deltab.00001'
     INQUIRE(file=filename,exist=ok_file2)
   else
     filename=TRIM(initfile(ilevel))//'/ic_deltab'
     INQUIRE(file=filename,exist=ok_file2)
   endif
   ok_file = ok_file1 .or. ok_file2

  end if

  if(cosmo .and. ok_file)then

     !-------------------------------------------------------------------------
     ! First step: compute level boundaries in terms of initial condition array
     !-------------------------------------------------------------------------
     if(ncache>0)then
     i1_min=n1(ilevel)+1; i1_max=0
     i2_min=n2(ilevel)+1; i2_max=0
     i3_min=n3(ilevel)+1; i3_max=0
     do ind=1,twotondim
        do i=1,ncache
           igrid=active(ilevel)%igrid(i)
           xx1=xg(igrid,1)+xc(ind,1)-skip_loc(1)
           xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
           xx2=xg(igrid,2)+xc(ind,2)-skip_loc(2)
           xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
           xx3=xg(igrid,3)+xc(ind,3)-skip_loc(3)
           xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
           i1_min=MIN(i1_min,int(xx1)+1)
           i1_max=MAX(i1_max,int(xx1)+1)
           i2_min=MIN(i2_min,int(xx2)+1)
           i2_max=MAX(i2_max,int(xx2)+1)
           i3_min=MIN(i3_min,int(xx3)+1)
           i3_max=MAX(i3_max,int(xx3)+1)
        end do
     end do
     error=.false.
     if(i1_min<1.or.i1_max>n1(ilevel))error=.true.
     if(i2_min<1.or.i2_max>n2(ilevel))error=.true.
     if(i3_min<1.or.i3_max>n3(ilevel))error=.true.
     if(error) then
        write(*,*)'Some grid are outside initial conditions sub-volume'
        write(*,*)'for ilevel=',ilevel
        write(*,*)i1_min,i1_max
        write(*,*)i2_min,i2_max
        write(*,*)i3_min,i3_max
        write(*,*)n1(ilevel),n2(ilevel),n3(ilevel)
        call clean_stop
     end if
     endif

     !-----------------------------------------
     ! Second step: read initial condition file
     !-----------------------------------------
     ! Allocate initial conditions array
     if(ncache>0)allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
     allocate(init_plane(1:n1(ilevel),1:n2(ilevel)))
     ! Loop over input variables
     do ivar=1,nvar
        if(cosmo)then
           ! Read baryons initial overdensity and displacement at a=aexp
           if(multiple)then
              call title(myid,nchar)
              if(ivar==1)filename=TRIM(initfile(ilevel))//'/dir_deltab/ic_deltab.'//TRIM(nchar)
              if(ivar==2)filename=TRIM(initfile(ilevel))//'/dir_velcx/ic_velcx.'//TRIM(nchar)
              if(ivar==3)filename=TRIM(initfile(ilevel))//'/dir_velcy/ic_velcy.'//TRIM(nchar)
              if(ivar==4)filename=TRIM(initfile(ilevel))//'/dir_velcz/ic_velcz.'//TRIM(nchar)
              if(ivar==5)filename=TRIM(initfile(ilevel))//'/dir_tempb/ic_tempb.'//TRIM(nchar)
           else
              if(ivar==1)filename=TRIM(initfile(ilevel))//'/ic_deltab'
              if(ivar==2)filename=TRIM(initfile(ilevel))//'/ic_velcx'
              if(ivar==3)filename=TRIM(initfile(ilevel))//'/ic_velcy'
              if(ivar==4)filename=TRIM(initfile(ilevel))//'/ic_velcz'
              if(ivar==5)filename=TRIM(initfile(ilevel))//'/ic_tempb'
           endif
        else
           ! Read primitive variables
           if(ivar==1)filename=TRIM(initfile(ilevel))//'/ic_d'
           if(ivar==2)filename=TRIM(initfile(ilevel))//'/ic_u'
           if(ivar==3)filename=TRIM(initfile(ilevel))//'/ic_v'
           if(ivar==4)filename=TRIM(initfile(ilevel))//'/ic_w'
           if(ivar==5)filename=TRIM(initfile(ilevel))//'/ic_p'
        endif
        call title(ivar,ncharvar)
        if(ivar>5)then
           call title(ivar-5,ncharvar)
           filename=TRIM(initfile(ilevel))//'/ic_pvar_'//TRIM(ncharvar)
        endif

        INQUIRE(file=filename,exist=ok_file3)
        if(ok_file3)then
           ! Reading the existing file
           if(myid==1)write(*,*)'Reading file '//TRIM(filename)
           if(multiple)then
              ilun=ncpu+myid+10

              ! Wait for the token
#ifndef WITHOUTMPI
              if(IOGROUPSIZE>0) then
                 if (mod(myid-1,IOGROUPSIZE)/=0) then
                    call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                         & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
                 end if
              endif
#endif

              open(ilun,file=filename,form='unformatted')
              rewind ilun
              read(ilun) ! skip first line
              do i3=1,n3(ilevel)
                 read(ilun) ((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                 if(ncache>0)then
                    if(i3.ge.i3_min.and.i3.le.i3_max)then
                       init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                            & init_plane(i1_min:i1_max,i2_min:i2_max)
                    end if
                 endif
              end do
              close(ilun)
              ! Send the token

#ifndef WITHOUTMPI
              if(IOGROUPSIZE>0) then
                 if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
                    dummy_io=1
                    call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
                         & MPI_COMM_WORLD,info2)
                 end if
              endif
#endif
           else
              if(myid==1)then
                 open(10,file=filename,form='unformatted')
                 rewind 10
                 read(10) ! skip first line
              endif
              do i3=1,n3(ilevel)
                 if(myid==1)then
                    read(10) ((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                 else
                    init_plane=0.0
                 endif
                 buf_count=n1(ilevel)*n2(ilevel)
#ifndef WITHOUTMPI
                 call MPI_BCAST(init_plane,buf_count,MPI_REAL,0,MPI_COMM_WORLD,info)
#endif
                 if(ncache>0)then
                    if(i3.ge.i3_min.and.i3.le.i3_max)then
                       init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                            & init_plane(i1_min:i1_max,i2_min:i2_max)
                    end if
                 endif
              end do
              if(myid==1)close(10)
           endif
        else
           ! If file doesn't exist, initialize variable to default value
           ! In most cases, this is zero (you can change that if necessary)
           if(myid==1)write(*,*)'File '//TRIM(filename)//' not found'
           if(myid==1)write(*,*)'Initialize corresponding variable to default value'
           if(ncache>0)then
              init_array=0d0
              ! Default value for metals
              if(cosmo.and.ivar==imetal.and.metal)init_array=z_ave*0.02 ! from solar units
              ! Default value for ionization fraction
              if(cosmo)xval=sqrt(omega_m)/(h0/100.*omega_b) ! From the book of Peebles p. 173
              if(cosmo.and.ivar==ixion.and.aton)init_array=1.2d-5*xval
#if defined(krome_on)
            ! sarebbero da leggere da file in funzione di z (anche T2_start)
            if((ivar.ge. krome_off+1) .and.(ivar.le.krome_off+9)) then
              init_array = krome_init_var(ivar-krome_off)
              if(myid==1) then
                if(ivar.eq. krome_off+1) then
                   print *,'Initializing krome:'
                endif
                print *,'                   ',krome_names_arr(ivar-krome_off),' to',krome_init_var(ivar-krome_off)
              endif
            endif
#endif
           endif
        endif

        if(ncache>0)then

        ! For cosmo runs, rescale initial conditions to code units
        if(cosmo)then
           ! Compute approximate average temperature in K
           if(.not. cooling)T2_start=1.356d-2/aexp**2
#if defined(krome_on)
          T2_start = krome_init_temp
          if((ivar==ndim+2) .and. (myid==1)) then
            print *,'  Setting init T to ',krome_init_temp,' (via precomputed table from Galli&Palla+98)'
          endif
#endif
           if(ivar==1)init_array=(1.0+dfact(ilevel)*init_array)*omega_b/omega_m
           if(ivar==2)init_array=dfact(ilevel)*vfact(1)*dx_loc/dxini(ilevel)*init_array/vfact(ilevel)
           if(ivar==3)init_array=dfact(ilevel)*vfact(1)*dx_loc/dxini(ilevel)*init_array/vfact(ilevel)
           if(ivar==4)init_array=dfact(ilevel)*vfact(1)*dx_loc/dxini(ilevel)*init_array/vfact(ilevel)
           if(ivar==ndim+2)init_array=(1.0+init_array)*T2_start/scale_T2
        endif

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ncache
              igrid=active(ilevel)%igrid(i)
              icell=igrid+iskip
              xx1=xg(igrid,1)+xc(ind,1)-skip_loc(1)
              xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
              xx2=xg(igrid,2)+xc(ind,2)-skip_loc(2)
              xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
              xx3=xg(igrid,3)+xc(ind,3)-skip_loc(3)
              xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
              i1=int(xx1)+1
              i1=int(xx1)+1
              i2=int(xx2)+1
              i2=int(xx2)+1
              i3=int(xx3)+1
              i3=int(xx3)+1
              ! Scatter to corresponding primitive variable
              uold(icell,ivar)=init_array(i1,i2,i3)
           end do
        end do
        ! End loop over cells
        endif
     end do
     ! End loop over input variables

     ! Deallocate initial conditions array
     if(ncache>0)deallocate(init_array)
     deallocate(init_plane)

     !----------------------------------------------------------------
     ! For cosmology runs: compute pressure, prevent negative density
     !----------------------------------------------------------------
     if(cosmo)then
        ! Loop over grids by vector sweeps
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           end do
           ! Loop over cells
           do ind=1,twotondim
              ! Gather cell indices
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              ! Prevent negative density
              do i=1,ngrid
                 rr=max(uold(ind_cell(i),1),0.1*omega_b/omega_m)
                 uold(ind_cell(i),1)=rr
              end do
              ! Compute pressure from temperature and density
              do i=1,ngrid
                 uold(ind_cell(i),ndim+2)=uold(ind_cell(i),1)*uold(ind_cell(i),ndim+2)
              end do
           end do
           ! End loop over cells
        end do
        ! End loop over grids
     end if

     !---------------------------------------------------
     ! Third step: compute initial conservative variables
     !---------------------------------------------------
     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        vy=0.0
        vz=0.0
        ! Loop over cells
        do ind=1,twotondim
           ! Gather cell indices
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           ! Compute total energy density
           do i=1,ngrid
              rr=uold(ind_cell(i),1)
              vx=uold(ind_cell(i),2)
#if NDIM>1
              vy=uold(ind_cell(i),3)
#endif
#if NDIM>2
              vz=uold(ind_cell(i),4)
#endif
              pp=uold(ind_cell(i),ndim+2)
              ek=0.5d0*(vx**2+vy**2+vz**2)
              ei=pp/(gamma-1.0)
              vv(i)=ei+rr*ek
           end do
           ! Scatter to corresponding conservative variable
           do i=1,ngrid
              uold(ind_cell(i),ndim+2)=vv(i)
           end do
           ! Compute momentum density
           do ivar=1,ndim
              do i=1,ngrid
                 rr=uold(ind_cell(i),1)
                 vx=uold(ind_cell(i),ivar+1)
                 vv(i)=rr*vx
              end do
              ! Scatter to corresponding conservative variable
              do i=1,ngrid
                 uold(ind_cell(i),ivar+1)=vv(i)
              end do
           end do
#if NVAR > NDIM + 2
           ! Compute passive variable density
           do ivar=ndim+3,nvar
              do i=1,ngrid
                 rr=uold(ind_cell(i),1)
                 uold(ind_cell(i),ivar)=rr*uold(ind_cell(i),ivar)
              end do
           enddo
#endif
        end do
        ! End loop over cells

     end do
     ! End loop over grids

  !-------------------------------------------------------
  ! Compute initial conditions from subroutine condinit
  !-------------------------------------------------------

  !!! IF NOT COSMO
  else

    foldername=TRIM(initfile(levelmin))

    varnames(1) = TRIM('d')
    varnames(2) = TRIM('vx')
#if NDIM>1
    varnames(3) = TRIM('vy')
#endif
#if NDIM>2
    varnames(4) = TRIM('vz')
#endif
    varnames(ndim+2)  = TRIM('P')
if(metal) then
    varnames(ndim+3)  = TRIM('Z')
endif
    varnames(krome_off+1)  = TRIM('x1')
    varnames(krome_off+2)  = TRIM('x2')
    varnames(krome_off+3)  = TRIM('x3')
    varnames(krome_off+4)  = TRIM('x4')
    varnames(krome_off+5)  = TRIM('x5')
    varnames(krome_off+6)  = TRIM('x6')
    varnames(krome_off+7)  = TRIM('x7')
    varnames(krome_off+8)  = TRIM('x8')
    varnames(krome_off+9)  = TRIM('x9')


!     !find max and minimum of refined regions
!     xmin = boxlen
!     ymin = boxlen
!     zmin = boxlen
!     xmax = 0.
!     ymax = 0.
!     zmax = 0.
!
!     do igrid=1,ncache,nvector
!        ngrid=MIN(nvector,ncache-igrid+1)
!        do i=1,ngrid
!           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
!        end do
!        ! Loop over cells
!        do ind=1,twotondim
!           ! Gather cell indices
!           iskip=ncoarse+(ind-1)*ngridmax
!           do i=1,ngrid
!              ind_cell(i)=iskip+ind_grid(i)
!           end do
!           ! Gather cell centre positions
!           do idim=1,ndim
!              do i=1,ngrid
!                 xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
!              end do
!           end do
!           ! Rescale position from code units to user units
!           do idim=1,ndim
!              do i=1,ngrid
!                 xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
!              end do
!           end do
!           if (xx(i, 1) > xmax) xmax=xx(i,1)
!           if (xx(i, 2) > ymax) ymax=xx(i,2)
!           if (xx(i, 3) > zmax) zmax=xx(i,3)
!           if (xx(i, 1) < xmin) xmin=xx(i,1)
!           if (xx(i, 2) < ymin) ymin=xx(i,2)
!           if (xx(i, 3) < zmin) zmin=xx(i,3)
!        enddo    
!     enddo
!     
!
!    if(myid==1) then 
!       print*, '.. reading condinit, level', ilevel
!       print*, '   xmin, xmax : ', xmin, xmax
!       print*, '   ymin, ymax : ', ymin, zmax
!       print*, '   zmin, zmax : ', ymin, zmax
!    endif

    !!! LOOP OVER VARIABLES
    do ivar=1,nvar

      ok_file = .false.
      file_var = 0.d0
      if (ilevel < 10)then
        fmt_lev = '(I1)'
      else
        fmt_lev= '(I2)'
      endif

      write (str_lev,fmt_lev) ilevel
      if (ilevel < 10) str_lev = TRIM('0') // TRIM(str_lev)
      filename = TRIM(foldername) // TRIM('/ic_') // TRIM(varnames(ivar))//TRIM('_') // TRIM(str_lev) // TRIM('.dat')
      INQUIRE(file=filename, exist=ok_file)

      if(ok_file)then
        open(10,file=filename,form='unformatted')
        read(10) file_var
        close(10)
      endif

     ! Loop over grids by vector sweeps

     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        ! Loop over cells
        do ind=1,twotondim
           ! Gather cell indices
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           ! Gather cell centre positions
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              end do
           end do
           ! Rescale position from code units to user units
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
              end do
           end do
           ! Call initial condition routine
           call condinit(xx,uu,dx_loc,ngrid,ilevel, ivar, ok_file, file_var)

           ! Scatter variables
           do i=1,ngrid
                 if(ivar==1)then
                   uold(ind_cell(i),1)=uu(i,1)
                 endif
                 ! velocity -> momentum
                 if(ivar==2)then
                   uold(ind_cell(i),2)=uold(ind_cell(i),1)*uu(i,2) 
                 endif
#if NDIM>1
                 if(ivar==3)then
                   uold(ind_cell(i),3)=uold(ind_cell(i),1)*uu(i,3) 
                 endif
#endif
#if NDIM>2
                 if(ivar==4)then
                   uold(ind_cell(i),4)=uold(ind_cell(i),1)*uu(i,4) 
                 endif
#endif
                 ! kinetic energy
                 if(ivar==ndim+2)then
                   uold(ind_cell(i),ndim+2)=0.0d0
                   uold(ind_cell(i),ndim+2)=uold(ind_cell(i),ndim+2)+0.5/uold(ind_cell(i),1)*uold(ind_cell(i),2)**2
#if NDIM>1
                   uold(ind_cell(i),ndim+2)=uold(ind_cell(i),ndim+2)+0.5/uold(ind_cell(i),1)*uold(ind_cell(i),3)**2
#endif
#if NDIM>2
                   uold(ind_cell(i),ndim+2)=uold(ind_cell(i),ndim+2)+0.5/uold(ind_cell(i),1)*uold(ind_cell(i),4)**2
#endif
                 ! thermal pressure -> total fluid energy
                   uold(ind_cell(i),ndim+2)=uold(ind_cell(i),ndim+2)+uu(i,ndim+2)/(gamma-1.0d0)
                 endif
#if NENER>0
                 ! radiative pressure -> radiative energy
                 ! radiative energy -> total fluid energy
                 if ((ivar>=ndim+2+1) .and. (ivar<=ndim+2+nener))then
                    uold(ind_cell(i),ivar)=uu(i,ivar)/(gamma_rad(ivar-2-ndim)-1.0d0)
                    uold(ind_cell(i),ndim+2)=uold(ind_cell(i),ndim+2)+uold(ind_cell(i),ivar)
                 endif
#endif
#if NVAR>NDIM+2+NENER
                 ! passive scalars
                 if ((ivar>=ndim+3+nener) .and. (ivar<=nvar))then
                    uold(ind_cell(i),ivar)=uold(ind_cell(i),1)*uu(i,ivar)
                 endif
#endif
           end do

        end do
        ! End loop over cells
     end do
     ! End loop over grids

   end do ! End loop over variables
  end if

  deallocate(file_var)

111 format('   Entering init_flow_fine for level ',I2)

end subroutine init_flow_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine region_condinit(x,q,dx,nn,ilevel, ivar, ok_file, file_var)
  use amr_commons, ONLY: myid
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn,ilevel,nline
  real(dp)::dx
  real(dp),dimension(1:nvector,1:nvar)::q
  real(dp),dimension(1:nvector,1:ndim)::x

#if NDIM==1
  real(dp),dimension(file_ncell(ilevel))::file_var
#endif
#if NDIM==2
  real(dp),dimension(file_ncell(ilevel),file_ncell(ilevel))::file_var
#endif
#if NDIM==3
  real,dimension(file_ncell(ilevel),file_ncell(ilevel),file_ncell(ilevel))::file_var
#endif

  logical::error,ok_file
  integer::i,k,ii,kk
  real(dp)::vol,r,xn,yn,zn,en
  integer::ivar

  ! Set some (tiny) default values in case n_region=0
  if (ivar==1) q(1:nn,1)=smallr
  if (ivar==2) q(1:nn,2)=0.0d0
#if NDIM>1
  if (ivar==3) q(1:nn,3)=0.0d0
#endif
#if NDIM>2
  if (ivar==4) q(1:nn,4)=0.0d0
#endif
  if (ivar==ndim+2) q(1:nn,ndim+2)=smallr*smallc**2/gamma
#if NVAR > NDIM + 2
  if ((ivar>=ndim+3) .and. (ivar<=nvar)) q(1:nn,ivar)=0.0d0
#endif

  if (ok_file)then
    do i=1,nn
       xn=0.0d0; yn=0.0d0; zn=0.0d0

       xn = (x(i,1) - file_x_offset(ilevel)) / dx + 0.5 
       yn = (x(i,2) - file_y_offset(ilevel)) / dx + 0.5
       zn = (x(i,3) - file_z_offset(ilevel)) / dx + 0.5

       if (xn < 1 .or. xn > file_ncell(ilevel))then
           print*, 'WARNING: refined cell not in ic', ilevel, xn, x(i,1), file_x_offset(ilevel), dx, file_ncell(ilevel)
       endif
       if (yn < 1 .or. yn > file_ncell(ilevel))then
           print*, 'WARNING: refined cell not in ic', ilevel, yn, x(i,2), file_y_offset(ilevel), dx, file_ncell(ilevel)
       endif
       if (zn < 1 .or. zn > file_ncell(ilevel))then
           print*, 'WARNING: refined cell not in ic', ilevel, zn, x(i,3), file_z_offset(ilevel), dx, file_ncell(ilevel)
       endif
 
#if NDIM==1
       q(i, ivar) = file_var(INT(xn))
#endif
#if NDIM==2
       q(i, ivar) = file_var(INT(xn), INT(yn))
#endif
#if NDIM==3
       q(i, ivar) = file_var(INT(xn), INT(yn), INT(zn))
#endif

    end do

  else ! no file, loading ic from namelist

  ! Loop over initial conditions regions
   do k=1,nregion
     ! For "square" regions only:
     if(region_type(k) .eq. 'square')then
        ! Exponent of choosen norm
        en=exp_region(k)
        do i=1,nn
           ! Compute position in normalized coordinates
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-x_center(k))/length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-y_center(k))/length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-z_center(k))/length_z(k)
#endif
           ! Compute cell "radius" relative to region center
           if(exp_region(k)<10)then
              r=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              r=max(xn,yn,zn)
           end if
           ! If cell lies within region,
           ! REPLACE primitive variables by region values
           if(r<1.0)then

             if (ivar==1) q(i,1)=d_region(k)

             if (ivar==2) q(i,2)=u_region(k)
#if NDIM>1
             if (ivar==3) q(i,3)=v_region(k)
#endif
#if NDIM>2
             if (ivar==4) q(i,4)=w_region(k)
#endif
             if (ivar==ndim+2) q(i,ndim+2)=p_region(k)
#if NENER>0
             if (ivar >= ndim+2+1 .and. ivar <= ndim+2+NENER) then
                 q(i,ivar)=prad_region(k,ivar-ndim-2)
             endif
#endif
#if NVAR>NDIM+2+NENER
             if (ivar >= ndim+3+nener .and. ivar <= nvar) then
                 q(i,ivar)=var_region(k,ivar-NDIM-2-NENER)
             endif
#endif
           end if

        end do
     end if

     ! For "point" regions only:
     if(region_type(k) .eq. 'point')then
        ! Volume elements
        vol=dx**ndim
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=1.0; yn=1.0; zn=1.0
           xn=max(1.0-abs(x(i,1)-x_center(k))/dx,0.0_dp)
#if NDIM>1
           yn=max(1.0-abs(x(i,2)-y_center(k))/dx,0.0_dp)
#endif
#if NDIM>2
           zn=max(1.0-abs(x(i,3)-z_center(k))/dx,0.0_dp)
#endif
           r=xn*yn*zn
           ! If cell lies within CIC cloud,
           ! ADD to primitive variables the region values
           if (ivar==1) q(i,1)=q(i,1)+d_region(k)*r/vol
           if (ivar==2) q(i,2)=q(i,2)+u_region(k)*r
#if NDIM>1
           if (ivar==3) q(i,3)=q(i,3)+v_region(k)*r
#endif
#if NDIM>2
           if (ivar==4) q(i,4)=q(i,4)+w_region(k)*r
#endif
           if (ivar==ndim+2) q(i,ndim+2)=q(i,ndim+2)+p_region(k)*r/vol
#if NENER>0
           if (ivar >= ndim+2+1 .and. ivar <= ndim+2+NENER) then
               q(i,ivar)=prad_region(k,ivar-ndim-2)*r/vol
           endif
#endif
#if NVAR>NDIM+2+NENER
           if (ivar >= ndim+3+nener .and. ivar <= nvar) then
               q(i,ivar)=var_region(k,ivar-NDIM-2-NENER)
           endif
#endif
        end do
     end if

   end do ! end loop for regions

 end if ! end if for file or regions

  return
end subroutine region_condinit
