!----------------------------------------------------------------------------------
!
!     This program creates segmented solids in a STL file
!     input : STL file contains many solids.
!     output: STL files contain a one solid.
! 
!     Jungwoo Kim.(2023.02)
!
!----------------------------------------------------------------------------------

      ! program debug
      ! implicit none

      ! character*80      :: filename
      ! integer           :: nlsfile
      ! character*80, allocatable, dimension(:)     ::  filestl

      ! filename = 'sinchon.stl'

      ! call count_solids(filename, nlsfile)

      ! print *, 'nlsfile=',nlsfile
      ! allocate(filestl(nlsfile))
      ! call devide_file(filename, nlsfile, filestl)
      ! deallocate(filestl)
      ! end program

!---------------------------------------------------------------------------------------
! 
!     count_solids : count number of solids in original STL file.
!
!---------------------------------------------------------------------------------------
      subroutine count_solids(orgstlfile, nlsfile)

      implicit none

      character*80, intent( in)                       ::  orgstlfile
      character*80                                    ::  cname
      integer     , intent(out)                       ::  nlsfile
      integer                                         ::  stat

      nlsfile = 0

      open(11,file=trim(orgstlfile))

      do
        read(11,*,iostat=stat) cname
        if(stat.NE.0) exit
        if(trim(cname).EQ.'endsolid') then
          nlsfile = nlsfile + 1
        endif
      enddo
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

      character*80, intent( in)                         ::  orgstlfile
      integer     , intent( in)                         ::  nlsfile

      character*80, intent(out), dimension(nlsfile)     ::  filestl

 !    local variables
      integer                  , dimension(nlsfile)     ::  nelem
      integer                                           ::  solidnum, nvert
      character*80                                      ::  solidnum2str
      double precision, allocatable, dimension(:,:)     ::  dpts, nvts


      call count_element(orgstlfile, nlsfile, nelem)

      open(123, FILE=trim(orgstlfile))

      do solidnum = 1, nlsfile

        nvert = nelem(solidnum)*3
        allocate(dpts(nvert,3), nvts(nelem(solidnum),3))

        call read_stl_solid(nvert, nelem(solidnum), dpts, nvts)

        write(solidnum2str,*) solidnum
        solidnum2str = adjustl(solidnum2str)
        solidnum2str = trim(solidnum2str)//'.stl'
        filestl(solidnum)= 'devidedsolid/solid'//trim(solidnum2str)

        call write_stl_solid(filestl(solidnum),nvert,nelem(solidnum),dpts,nvts)

        deallocate(dpts, nvts)

      enddo

      close(123)

      end subroutine devide_file

!---------------------------------------------------------------------------------------
!
!     count_Element : count number of elements at each solids.
!
!---------------------------------------------------------------------------------------
      subroutine count_element(orgstlfile, nlsfile, nelem)

      implicit none

      character*80, intent( in)                       ::  orgstlfile
      integer     , intent( in)                       ::  nlsfile

      integer     , intent(out)                       ::  nelem(nlsfile)
 !    local variables
      integer       ::  line, solidnum
      integer       ::  stat
      character*80  ::  cname

      line=0
      solidnum = 0
      open(11,file=trim(orgstlfile))

      do
        read(11,*,iostat=stat) cname
        if(stat.NE.0) exit
        line=line+1
        if(trim(cname).EQ.'endsolid') then
            solidnum = solidnum + 1
            nelem(solidnum) = (line-2)/7
            line = 0
        endif
      enddo

      close(11)

      if(solidnum.NE.nlsfile) then
        print*, 'count element error in soliddevider'
        stop
      endif

      return
      end subroutine count_element

!---------------------------------------------------------------------------------------
! 
!     read_stl_solid : read stl file from original stl file
!
!---------------------------------------------------------------------------------------
      subroutine read_stl_solid(nvert, nelem, dpts, nvts)
      implicit none

      integer         , intent( in)                             ::  nvert, nelem

      double precision, intent(out), dimension(nvert,3)         ::  dpts
      double precision, intent(out), dimension(nelem,3)         ::  nvts
 !    local variables
      character*80                                              ::  cname
      integer                                                   ::  i, nn

      read(123,*)
      nn = 0
      do i=1,nelem
      read(123,*) cname,cname,nvts(i,1),nvts(i,2),nvts(i,3)
      read(123,*) cname
      nn = nn + 1
      read(123,*) cname,dpts(nn,1),dpts(nn,2),dpts(nn,3)
      nn = nn + 1
      read(123,*) cname,dpts(nn,1),dpts(nn,2),dpts(nn,3)
      nn = nn + 1
      read(123,*) cname,dpts(nn,1),dpts(nn,2),dpts(nn,3)
      read(123,*) cname
      read(123,*) cname
      enddo
      read(123,*) cname
      
      if(cname.NE.'endsolid') then
        print*, 'error read solid in soliddevider'
        stop
      endif

      return
      end subroutine read_stl_solid

!---------------------------------------------------------------------------------------
! 
!     write_stl_solid : write each solid's stl file
!
!---------------------------------------------------------------------------------------
      subroutine write_stl_solid(filename,nvert,nelem,dpts,nvts)
      implicit none

      character*80                                         ::  filename
      integer                                              ::  nvert, nelem
      double precision, dimension(nvert,3)                 ::  dpts
      double precision, dimension(nelem,3)                 ::  nvts
      
 !    local variables
      integer                                                   ::  i, nn

      open(321, FILE=trim(filename))

      write(321,'(A13)') 'solid devided'
      nn = 0
      do i = 1,nelem
      write(321,'(2X,A5,1X,A6,3E)') 'facet','normal',nvts(i,1),nvts(i,2),nvts(i,3)
      write(321,'(4X,A10)') 'outer loop'
      nn = nn +1
      write(321,'(6X,A,3E)') 'vertex',dpts(nn,1),dpts(nn,2),dpts(nn,3)
      nn = nn +1
      write(321,'(6X,A,3E)') 'vertex',dpts(nn,1),dpts(nn,2),dpts(nn,3)
      nn = nn +1
      write(321,'(6X,A,3E)') 'vertex',dpts(nn,1),dpts(nn,2),dpts(nn,3)
      write(321,'(4X,A8)') 'end loop'
      write(321,'(2X,A9)') 'end facet'
      enddo

      write(321,'(A8)') 'endsolid'

      close(321)
      
      end subroutine write_stl_solid
