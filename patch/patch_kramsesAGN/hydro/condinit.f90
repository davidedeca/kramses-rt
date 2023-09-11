!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn,ilevel, ivar, ok_file, file_var)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer::ilevel
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  logical::ok_file
  integer::myid, i

#if NDIM==1
    real(dp),dimension(2**lev_file)::file_var
#endif
#if NDIM==2
    real(dp),dimension(2**lev_file,2**lev_file)::file_var
#endif
#if NDIM==3
    real,dimension(2**lev_file,2**lev_file,2**lev_file)::file_var
#endif

  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
#if NENER>0 || NVAR>NDIM+2+NENER
  integer::ivar
#endif
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables

  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn,ilevel, ivar, ok_file, file_var)

  u(1:nn,ivar)=q(1:nn,ivar)
  
  ! Add here, if you wish, some user-defined initial conditions
  ! ........

  ! Convert primitive to conservative variables
  ! density -> density


end subroutine condinit
