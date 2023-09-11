subroutine write_screen
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use rt_hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  integer::igrid,jgrid,ind,icpu
  integer::i,icell,ncell,ilevel,ncache
  integer::icellmin,nx_loc
  real(dp)::dx,scale,smallp,ddd,ppp
  real(dp)::sss1, sss2, sss3, sss4,sss5,sss6,sss7,sss8,sss9
  character(LEN=5)::nchar
  character(LEN=80)::filename
  integer,dimension(:),allocatable::ind_grid,ind_cell,ind_sort,ll,ll_all
  real(qdp),dimension(:),allocatable::rr,et,ei,dd,uu,mm,gg,dtot
  real(qdp),dimension(:),allocatable::ss1, ss2, ss3,ss4,ss5,ss6,ss7,ss8,ss9
  real(qdp),dimension(:),allocatable::rr_all,et_all,ei_all
  real(qdp),dimension(:),allocatable::dd_all,uu_all,mm_all,gg_all,dtot_all
  real(qdp),dimension(:),allocatable::ss1_all, ss2_all,ss3_all,ss4_all,ss5_all,ss6_all,ss7_all,ss8_all,ss9_all
#if NENER>0
  integer::irad
  real(qdp),dimension(:,:),allocatable::prad_all,prad
#endif
#if defined(RT)
  integer::irt
  real(qdp),dimension(:,:),allocatable:: ff, ff_all   !flux
#endif

  integer,dimension(1:ncpu)::iskip,ncell_loc,ncell_all

  if(ndim>1)return

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

  ncell=0
  do ilevel=1,nlevelmax
     ncache=numbl(myid,ilevel)
     if(ncache > 0)then
        allocate(ind_grid(1:ncache),ind_cell(1:ncache))
        ! Gather all grids
        igrid=headl(myid,ilevel)
        do jgrid=1,ncache
           ind_grid(jgrid)=igrid
           igrid=next(igrid)
        end do
        ! Count leaf cells
        do ind=1,twotondim
           do i=1,ncache
              ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
           end do
           do i=1,ncache
              if(son(ind_cell(i))== 0)then
                 ncell=ncell+1
              end if
           end do
        end do
        deallocate(ind_grid, ind_cell)
     end if
  end do

  ncell_loc=0
  ncell_all=0
  ncell_loc(myid)=ncell
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ncell_loc,ncell_all,ncpu,MPI_INTEGER,MPI_SUM,&
       & MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ncell_all=ncell_loc
#endif

  ncell=0
  iskip=0
  do icpu=1,ncpu
     iskip(icpu)=ncell
     ncell=ncell+ncell_all(icpu)
  end do

  if(myid==1)write(*,114)ncell

  if(ncell>0)then

  allocate(rr(1:ncell),mm(1:ncell),dd(1:ncell),dtot(1:ncell),ss1(1:ncell),ss2(1:ncell),ss3(1:ncell),ss4(1:ncell),ss5(1:ncell),ss6(1:ncell),ss7(1:ncell),ss8(1:ncell),ss9(1:ncell))
  allocate(et(1:ncell),ei(1:ncell))
  allocate(uu(1:ncell),ll(1:ncell),gg(1:ncell))
  allocate(rr_all(1:ncell),mm_all(1:ncell),dd_all(1:ncell),dtot_all(1:ncell),ss1_all(1:ncell),ss2_all(1:ncell),ss3_all(1:ncell),ss4_all(1:ncell),ss5_all(1:ncell),ss6_all(1:ncell),ss7_all(1:ncell),ss8_all(1:ncell),ss9_all(1:ncell))
  allocate(et_all(1:ncell),ei_all(1:ncell))
  allocate(uu_all(1:ncell),ll_all(1:ncell),gg_all(1:ncell))

  rr=0.0D0; mm=0.0D0; dd=0.0D0; dtot=0.0D0; et=0.0D0
  ei=0.0D0; uu=0.0D0; gg=0.0D0; ll=0; ss1=0.0D0; ss2=0.0D0;ss3=0.0D0;ss4=0.0D0;ss5=0.0D0;ss6=0.0D0;ss7=0.0D0;ss8=0.0D0;ss9=0.0D0
  rr_all=0.0D0; mm_all=0.0D0; dd_all=0.0D0; dtot_all=0.0D0; et_all=0.0D0
  ei_all=0.0D0; uu_all=0.0D0; gg_all=0.0D0; ll_all=0; ss1_all=0.0D0; ss2_all=0.0D0; ss3_all=0.0D0; ss4_all=0.0D0; ss5_all=0.0D0; ss6_all=0.0D0; ss7_all=0.0D0; ss8_all=0.0D0; ss9_all=0.0D0
#if NENER>0
  allocate(prad(1:ncell,1:nener),prad_all(1:ncell,1:nener))
  prad=0.0D0; prad_all=0.0D0
#endif
#if defined(RT)
  allocate(ff(1:ncell,1:ngroups),ff_all(1:ncell,1:ngroups))
  ff=0.0D0; ff_all=0.0D0
#endif

  icell=iskip(myid)
  do ilevel=1,nlevelmax
     icellmin=icell
     ncache=numbl(myid,ilevel)
     if(ncache > 0)then
        dx=0.5D0**ilevel
        allocate(ind_grid(1:ncache),ind_cell(1:ncache))
        ! Gather all grids
        igrid=headl(myid,ilevel)
        do jgrid=1,ncache
           ind_grid(jgrid)=igrid
           igrid=next(igrid)
        end do
        ! Gather variables
        icell=icellmin
        do ind=1,twotondim
           do i=1,ncache
              ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
           end do
           do i=1,ncache
              if(son(ind_cell(i))==0)then
                 icell=icell+1
                 rr(icell)=xg(ind_grid(i),1)+(dble(ind)-1.5D0)*dx
                 ll(icell)=ilevel
              end if
           end do
        end do
        if(hydro)then
           icell=icellmin
           do ind=1,twotondim
              do i=1,ncache
                 ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
              end do
              do i=1,ncache
                 if(son(ind_cell(i))==0)then
                    icell=icell+1
                    dd(icell)=uold(ind_cell(i),1)
                    ss1(icell)=uold(ind_cell(i),4+NENER)/dd(icell)
                    ss2(icell)=uold(ind_cell(i),5+NENER)/dd(icell)
                    ss3(icell)=uold(ind_cell(i),6+NENER)/dd(icell)
                    ss4(icell)=uold(ind_cell(i),7+NENER)/dd(icell)
                    ss5(icell)=uold(ind_cell(i),8+NENER)/dd(icell)
                    ss6(icell)=uold(ind_cell(i),9+NENER)/dd(icell)
                    ss7(icell)=uold(ind_cell(i),10+NENER)/dd(icell)
                    ss8(icell)=uold(ind_cell(i),11+NENER)/dd(icell)
                    ss9(icell)=uold(ind_cell(i),12+NENER)/dd(icell)
                    
                    mm(icell)=dd(icell)
                    uu(icell)=uold(ind_cell(i),2)/dd(icell)
                    et(icell)=uold(ind_cell(i),3)
                    ei(icell)=et(icell)/dd(icell)-0.5d0*uu(icell)**2
                    ei(icell)=ei(icell)*dd(icell)
#if NENER>0
                      do irad=1,nener
                        ei(icell)=ei(icell)-uold(ind_cell(i),3+irad)
                        prad(icell,irad)=(gamma_rad(irad)-1.0d0)*uold(ind_cell(i),3+irad)
                      end do
#endif
#if defined(RT)
                      do irt=1,ngroups
                        ff(icell, irt)=rtuold(ind_cell(i),iGroups(irt)) * rt_c
                      end do
#endif
                 end if
              end do
           end do
        end if
        deallocate(ind_grid, ind_cell)
     end if
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rr,rr_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mm,mm_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dd,dd_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ss1,ss1_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ss2,ss2_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ss3,ss3_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ss4,ss4_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ss5,ss5_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ss6,ss6_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ss7,ss7_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ss8,ss8_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ss9,ss9_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dtot,dtot_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(et,et_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ei,ei_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(uu,uu_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(gg,gg_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ll,ll_all,ncell,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(rr,rr_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  rr=rr_all; mm=mm_all; dd=dd_all; dtot=dtot_all; et=et_all
  ei=ei_all; uu=uu_all; gg=gg_all; ll=ll_all; ss1=ss1_all; ss2=ss2_all; ss3=ss3_all
  ss4=ss4_all; ss5=ss5_all; ss6=ss6_all; ss7=ss7_all; ss8=ss8_all; ss9=ss9_all
#if NENER>0
  call MPI_ALLREDUCE(prad,prad_all,ncell*nener,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  prad=prad_all
#endif
#if defined(RT)
  call MPI_ALLREDUCE(ff,ff_all,ncell*ngroups,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  ff = ff_all
#endif

#endif

  if(myid==1)then
     write(*,*)'================================================'
     write(*,*)'t=', t

#if NENER>0

#if !defined(RT)
     if(poisson)then
        write(*,*)'lev, x, d, u, Pnt, P, H, H2, HII'
     else
        write(*,*)'lev, x ,d ,u, Pnt, P, H, H2, HII'
     endif
#endif
#if defined(RT)
     if(poisson)then
        write(*,*)'lev, x, d, u, Pnt, P, H, H2, HII, flux'
     else
        write(*,*)'lev, x ,d ,u, Pnt, P, H, H2, HII, flux'
    endif
#endif

#else

#if !defined(RT)
     if(poisson)then
        write(*,*)'lev, x, d, u, P, H, H2, HII'
     else
        write(*,*)'lev, x ,d ,u, P, H, H2, HII'
     endif
#endif
#if defined(RT)
     if(poisson)then
        write(*,*)'lev, x, d, u, P, HI, e, HII, HeI, HeII, HeIII, H-, H2, H2+, u1, u2, u3, u4, u5, u6, u7, u8, u9, u10'
     else
        write(*,*)'lev, x, d, u, P, HI, e, HII, HeI, HeII, HeIII, H-, H2, H2+, u1, u2, u3, u4, u5, u6, u7, u8, u9, u10'
     endif
#endif

#endif
     ! Sort radius
     allocate(ind_sort(1:ncell))
     call quick_sort(rr,ind_sort,ncell)
     ! Write results to screen
     smallp=smallc**2/gamma
     nx_loc=icoarse_max-icoarse_min+1
     scale=boxlen/dble(nx_loc)
     ! Prevent underflow for velocity
     do i=1,ncell
        if(ABS(uu(i))<smallc)uu(i)=0.0D0
     end do
     call title(ifout,nchar)
     filename = 'output_'//TRIM(nchar)
     open(unit=33, file=filename, form='formatted')
     do i=1,ncell
        ddd=MAX(dd(ind_sort(i)),smallr)
        ppp=MAX((gamma-1.0)*ei(ind_sort(i)),ddd*smallp)
        sss1=ss1(ind_sort(i))
        sss2=ss2(ind_sort(i))
        sss3=ss3(ind_sort(i))
        sss4=ss4(ind_sort(i))
        sss5=ss5(ind_sort(i))
        sss6=ss6(ind_sort(i))
        sss7=ss7(ind_sort(i))
        sss8=ss8(ind_sort(i))
        sss9=ss9(ind_sort(i))
        write(33,113) &
            & ll(ind_sort(i)),  &
             & (rr(i)-dble(icoarse_min))*scale, &
             & ddd , &
             & uu(ind_sort(i)), &
#if NENER>0
             & prad(ind_sort(i),1), &

#endif
             & ppp, &
             & sss1, &
             & sss2, &
#if !defined(RT)
             & sss3
#endif
#if defined(RT)
             & sss3, sss4, sss5, sss6, sss7, sss8, sss9,&
             & ff(ind_sort(i), 1), ff(ind_sort(i), 2), ff(ind_sort(i), 3), &
             & ff(ind_sort(i), 4), ff(ind_sort(i), 5), ff(ind_sort(i), 6), &
             & ff(ind_sort(i), 7), ff(ind_sort(i), 8), ff(ind_sort(i), 9), &
             & ff(ind_sort(i), 10)
#endif

     end do
     deallocate(ind_sort)
     close(33)
     write(*,*)'================================================'
  end if

  ! Deallocate local arrays
  deallocate(mm,rr,dd,dtot,et,ei,uu,ll,gg,ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,ss9)
  deallocate(mm_all,rr_all,dd_all,dtot_all,et_all,ei_all,uu_all,ll_all,gg_all,ss1_all, ss2_all, ss3_all,ss4_all, ss5_all, ss6_all, ss7_all, ss8_all, ss9_all)
#if NENER>0
  deallocate(prad,prad_all)
#endif
#if defined(RT)
  deallocate(ff, ff_all)
#endif
  end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

113 format(i3,1x,1pe12.5,1x,23(1pe12.5,1x))
114 format(' Output ',i5,' cells')

end subroutine write_screen
