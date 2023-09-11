subroutine backup_poisson(filename, filename_desc)
  use amr_commons
  use hydro_commons, only: gamma
  use poisson_commons
  use dump_utils, only : dump_header_info, generic_dump, dim_keys
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  character(len=80), intent(in) :: filename, filename_desc

  integer::i,ivar,idim,ncache,ind,ilevel,igrid,iskip,istart,ibound
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar, ncharcpu
  character(LEN=80)::fileloc
  integer :: unit_out, unit_info

  logical :: dump_info_flag
  character(len=100) :: field_name
  integer :: info_var_count

#ifndef WITHOUTMPI
  integer,parameter::tag=1123
  integer::dummy_io,info2
#endif

  if(verbose)write(*,*)'Entering backup_poisson'
  if(verbose)write(*,*)filename_desc
  !ilun=ncpu+myid+10

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

 ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif


  open(newunit=unit_out,file=fileloc,form='unformatted')

  if (myid == 1) then
     open(newunit=unit_info, file=filename_desc, form='formatted')
     call dump_header_info(unit_info)
     info_var_count = 1
     dump_info_flag = .true.
  else
     dump_info_flag = .false.
  end if

  write(unit_out)ncpu
  write(unit_out)ndim+1
  write(unit_out)ndim
  write(unit_out)nlevelmax
  write(unit_out)nboundary
  write(unit_out)gamma

  do ilevel= 1, nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        write(unit_out)ilevel
        write(unit_out)ncache
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              ! Write potential
              do i=1,ncache
                 xdp(i)=phi(ind_grid(i)+iskip)
              end do
              field_name='potential'
              call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              ! Write force
              do ivar=1,ndim
                 do i=1,ncache
                    xdp(i)=f(ind_grid(i)+iskip,ivar)
                 end do
                 field_name = 'force_' // dim_keys(ivar)
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do
              dump_info_flag=.false.
           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do
  close(unit_out)
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

end subroutine backup_poisson





