module interpolate_amr
 use kernels, only:radkernel2, radkernel,cnormk3D,wallint
 use amr_parameters
 implicit none
 integer, parameter :: doub_prec = kind(0.d0)
 integer, parameter :: sing_prec = kind(1.0)
 public :: interpolate !, interpolate_vec

contains

subroutine interpolate(x,y,z,hh,weight,dat,npart,datsmooth,xg,yg,zg,myid) 

 integer, intent(in) :: npart
 real(doub_prec), intent(in), dimension(npart) :: x,y,z,hh,weight,dat
 real(doub_prec):: datnorm

 integer :: nn,j, i, ipix, jpix, kpix,myid
 real(doub_prec) :: xi,yi,zi,hi,hi1,hi21,radkern,wab,q2,const,dyz2,dz2,boxlen
 integer :: nx,ny,nz,xmin, xmax, ymin, ymax, zmin, zmax
 real(doub_prec) :: cellsize
 real(doub_prec) :: term,termnorm,dy,dz,ypix,zpix,xpixi,pixwidthmax
 real(doub_prec), dimension(1:nvector) :: xg, yg, zg
 real(doub_prec), dimension(1:nvector), intent(out) :: datsmooth

 real :: pixint, wint
 logical, parameter :: exact_rendering = .true.

 datsmooth = 0.
 datnorm   = 0.

 const = cnormk3D

 do i=1,npart

  hi = hh(i)
  xi = x(i)  
  yi = y(i)
  zi = z(i)

  hi1 = 1./hi
  hi21 = hi1*hi1
  radkern = radkernel*hi  !2*smoothing_length

  do j=1,nvector

  ! rough selection.. a particle could still not be in the kernel
  if( xg(j) < xi - radkern ) cycle
  if( yg(j) < yi - radkern ) cycle
  if( zg(j) < zi - radkern ) cycle
  if( xg(j) > xi + radkern ) cycle
  if( yg(j) > yi + radkern ) cycle
  if( zg(j) > zi + radkern ) cycle

  !termnorm = const*weight(i)
  term = termnorm*dat(i)

  q2 = 0.
  q2 = q2 + (xg(j) - xi)**2
  q2 = q2 + (yg(j) - yi)**2
  q2 = q2 + (zg(j) - zi)**2
  q2 = q2 * hi21

  !! only case .not. exact_rendering available
  if (q2 < radkernel2) then

    wab = wkernel(real(q2))
    datsmooth(j) = datsmooth(j) + term*wab    
    !datnorm   = datnorm   + termnorm*wab 

  endif
  
 enddo
 enddo 

! where(datnorm > tiny(datnorm))
!     datsmooth = datsmooth/datnorm
! end where


end subroutine interpolate


!------------------------------------------------------------
! interface to kernel routine to avoid problems with openMP
!-----------------------------------------------------------
real function wkernel(q2)
 use kernels, only:wfunc
 real, intent(in) :: q2

 wkernel = wfunc(q2)

end function wkernel

end module interpolate_amr
