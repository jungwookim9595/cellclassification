!----------------------------------------------------------------------------------
!
!     This program creates segmented solids in a STL file
!     input : STL file contains many solids.
!     output: STL files contain a one solid.
!
!     Jungwoo Kim.(2023.02)
!     Updated for real(8) consistency (2026-03-23)
!
!----------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
!
!     count_solids : count number of solids in original STL file.
!
!---------------------------------------------------------------------------------------
subroutine count_solids(orgstlfile, nlsfile)
  implicit none

  character(len=80), intent(in)  :: orgstlfile
  integer, intent(out)           :: nlsfile

  character(len=80) :: cname
  integer :: stat

  nlsfile = 0

  open(11, file=trim(orgstlfile))

  do
     read(11, *, iostat=stat) cname
     if (stat /= 0) exit
     if (trim(cname) == 'endsolid') then
        nlsfile = nlsfile + 1
     end if
  end do
  close(11)

  return
end subroutine count_solids

!---------------------------------------------------------------------------------------
!
!     devide_file : read original STL file and write STL files for each solids.
!                   + save STL file names in array filestl
!
!---------------------------------------------------------------------------------------
subroutine devide_file(orgstlfile, nlsfile, filestl)
  implicit none

  character(len=80), intent(in)                      :: orgstlfile
  integer, intent(in)                                :: nlsfile
  character(len=80), intent(out), dimension(nlsfile) :: filestl

  ! local variables
  integer, dimension(nlsfile) :: nelem
  integer :: solidnum, nvert
  character(len=80) :: solidnum2str
  real(8), allocatable, dimension(:,:) :: dpts, nvts

  call count_element(orgstlfile, nlsfile, nelem)

  open(123, file=trim(orgstlfile))

  do solidnum = 1, nlsfile

     nvert = nelem(solidnum) * 3
     allocate(dpts(nvert, 3), nvts(nelem(solidnum), 3))

     call read_stl_solid(nvert, nelem(solidnum), dpts, nvts)

     write(solidnum2str, *) solidnum
     solidnum2str = adjustl(solidnum2str)
     solidnum2str = trim(solidnum2str) // '.stl'
     filestl(solidnum) = 'devidedsolid/solid' // trim(solidnum2str)

     call write_stl_solid(filestl(solidnum), nvert, nelem(solidnum), dpts, nvts)

     deallocate(dpts, nvts)

  end do

  close(123)

end subroutine devide_file

!---------------------------------------------------------------------------------------
!
!     count_Element : count number of elements at each solids.
!
!---------------------------------------------------------------------------------------
subroutine count_element(orgstlfile, nlsfile, nelem)
  implicit none

  character(len=80), intent(in)  :: orgstlfile
  integer, intent(in)           :: nlsfile
  integer, intent(out)          :: nelem(nlsfile)

  ! local variables
  integer :: line, solidnum, stat
  character(len=80) :: cname

  line = 0
  solidnum = 0
  open(11, file=trim(orgstlfile))

  do
     read(11, *, iostat=stat) cname
     if (stat /= 0) exit
     line = line + 1
     if (trim(cname) == 'endsolid') then
        solidnum = solidnum + 1
        nelem(solidnum) = (line - 2) / 7
        line = 0
     end if
  end do

  close(11)

  if (solidnum /= nlsfile) then
     print *, 'count element error in soliddevider'
     stop
  end if

  return
end subroutine count_element

!---------------------------------------------------------------------------------------
!
!     read_stl_solid : read stl file from original stl file
!
!---------------------------------------------------------------------------------------
subroutine read_stl_solid(nvert, nelem, dpts, nvts)
  implicit none

  integer, intent(in)                            :: nvert, nelem
  real(8), intent(out), dimension(nvert, 3)      :: dpts
  real(8), intent(out), dimension(nelem, 3)      :: nvts

  ! local variables
  character(len=80) :: cname
  integer :: i, nn

  read(123, *)
  nn = 0
  do i = 1, nelem
     read(123, *) cname, cname, nvts(i, 1), nvts(i, 2), nvts(i, 3)
     read(123, *) cname
     nn = nn + 1
     read(123, *) cname, dpts(nn, 1), dpts(nn, 2), dpts(nn, 3)
     nn = nn + 1
     read(123, *) cname, dpts(nn, 1), dpts(nn, 2), dpts(nn, 3)
     nn = nn + 1
     read(123, *) cname, dpts(nn, 1), dpts(nn, 2), dpts(nn, 3)
     read(123, *) cname
     read(123, *) cname
  end do
  read(123, *) cname

  if (cname /= 'endsolid') then
     print *, 'error read solid in soliddevider'
     stop
  end if

  return
end subroutine read_stl_solid

!---------------------------------------------------------------------------------------
!
!     write_stl_solid : write each solid's stl file
!
!---------------------------------------------------------------------------------------
subroutine write_stl_solid(filename, nvert, nelem, dpts, nvts)
  implicit none

  character(len=80), intent(in) :: filename
  integer, intent(in) :: nvert, nelem
  real(8), dimension(nvert, 3), intent(in) :: dpts
  real(8), dimension(nelem, 3), intent(in) :: nvts

  ! local variables
  integer :: i, nn

  open(321, file=trim(filename))

  write(321, '(A13)') 'solid devided'
  nn = 0
  do i = 1, nelem
     write(321, '(2X,A5,1X,A6,3E)') 'facet', 'normal', nvts(i, 1), nvts(i, 2), nvts(i, 3)
     write(321, '(4X,A10)') 'outer loop'
     nn = nn + 1
     write(321, '(6X,A,3E)') 'vertex', dpts(nn, 1), dpts(nn, 2), dpts(nn, 3)
     nn = nn + 1
     write(321, '(6X,A,3E)') 'vertex', dpts(nn, 1), dpts(nn, 2), dpts(nn, 3)
     nn = nn + 1
     write(321, '(6X,A,3E)') 'vertex', dpts(nn, 1), dpts(nn, 2), dpts(nn, 3)
     write(321, '(4X,A8)') 'end loop'
     write(321, '(2X,A9)') 'end facet'
  end do

  write(321, '(A8)') 'endsolid'

  close(321)

end subroutine write_stl_solid
