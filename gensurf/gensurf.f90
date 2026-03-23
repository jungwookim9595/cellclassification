!-----------------------------------------------------------------------------------
!
!     Program GenSurf.f90
!
!     This program is designed to generate the node surface points from TecPlot FE data
!     Input data: vertices, triangle information (index for triangle elements' vertices
!                                             and index for connected triangle elements)
!     Output data: position, normal vectors, velocity (optional)
!                  at vertices, edges, and centers of the elements (triangles)
!
!     07/01/08
!
!     Jung-Il Choi, Ph.D.,
!     Research Assistant Professor
!     Dept. Mechanical & Aerospace Engineering
!     North Carolina State University
!
!     Modernized to Fortran 90 free-format (2026-03-23)
!
!-----------------------------------------------------------------------------------

!     Define Data Structures MODULE

module GenSurf
  implicit none

  integer, parameter :: nshare = 500

  type vertex
     real(8), dimension(3) :: xyz, xyz0, norm, vel, vel0, acc
     real(8) :: area, power
     integer, dimension(nshare) :: share
     integer :: open
  end type

  type element
     real(8), dimension(3) :: xyz, norm
     real(8) :: area, power
     integer, dimension(3) :: vert, conn, edge
  end type

  type edges
     real(8), dimension(3) :: xyz, norm
     real(8) :: area, power
     integer, dimension(2) :: vert, elem
     integer :: open
  end type

  type surfpar
     integer :: nvert, nelem, nedge, nedgea, iedge
  end type

  type vertexorg
     real(8), dimension(3) :: xyz
     integer, dimension(nshare) :: share
  end type

  type elementorg
     integer, dimension(3) :: vert, conn
  end type

contains

  !---------------------------------------------------------------------------------------
  !
  !     STR2NUM : String to number
  !
  !---------------------------------------------------------------------------------------
  function str2num(string) result(ival)
    implicit none
    character(len=100), intent(in) :: string
    integer :: ival

    integer :: lens, istr, iend, i, inum

    lens = len(trim(string))
    istr = lens
    iend = 0
    do i = 1, lens
       if (ichar(string(i:i)) >= 48 .and. ichar(string(i:i)) < 58) then
          istr = min(i, istr)
          iend = max(i, iend)
       end if
    end do

    ival = 0
    do i = istr, iend
       inum = ichar(string(i:i)) - 48
       ival = ival + inum * 10**(iend - i)
    end do
  end function str2num

end module GenSurf


!---------------------------------------------------------------------------------------
!
!     Normal_Correction : Correction for Normal Vectors considering Connectivity
!                         Correct the order of the element storage!
!
!---------------------------------------------------------------------------------------
subroutine normal_correction(paramx, vert, elem, para, indx1, indx2, nc_corr)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertex), intent(inout) :: vert(paramx%nvert)
  type(element), intent(inout) :: elem(paramx%nelem)
  integer, intent(in) :: indx1, indx2
  integer, intent(inout) :: nc_corr

  integer :: indxv1(3), indxv2(3)
  integer :: m1, m2, n1, n2, n3, j, iflg
  integer :: mt1, mt2, nn3, nt1, nt2, nt3
  logical :: found

  indxv1(:) = elem(indx1)%vert(:)
  indxv2(:) = elem(indx2)%vert(:)

  ! find sharing edge
  found = .false.
  do m1 = 1, 3
     m2 = mod(m1, 3) + 1
     n1 = indxv1(m1)
     n2 = indxv1(m2)

     iflg = 0
     do j = 1, 3
        if (indxv2(j) == n1) iflg = iflg + 1
        if (indxv2(j) == n2) iflg = iflg + 1
        if (indxv2(j) /= n1 .and. indxv2(j) /= n2) n3 = indxv2(j)
     end do
     if (iflg == 2) then
        mt1 = m1
        mt2 = m2
        nn3 = n3
        found = .true.
        exit  ! replaced goto 10
     end if
  end do

  ! continue after finding sharing edge
  nt1 = 4 - mt1
  nt2 = 4 - mt2
  nt3 = 6 - nt1 - nt2

  if (indxv2(nt1) /= indxv1(mt1) .or. indxv2(nt2) /= indxv1(mt2)) nc_corr = nc_corr + 1
  elem(indx2)%vert(nt1) = indxv1(mt1)
  elem(indx2)%vert(nt2) = indxv1(mt2)
  elem(indx2)%vert(nt3) = nn3

  return
end subroutine normal_correction


!---------------------------------------------------------------------------------------
!
!     Normal_Correction_Org : Correction for Normal Vectors (original types)
!
!---------------------------------------------------------------------------------------
subroutine normal_correction_org(paramx, vert, elem, para, indx1, indx2, nc_corr)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertexorg), intent(inout) :: vert(paramx%nvert)
  type(elementorg), intent(inout) :: elem(paramx%nelem)
  integer, intent(in) :: indx1, indx2
  integer, intent(inout) :: nc_corr

  integer :: indxv1(3), indxv2(3)
  integer :: m1, m2, n1, n2, n3, j, iflg
  integer :: mt1, mt2, nn3, nt1, nt2, nt3
  logical :: found

  indxv1(:) = elem(indx1)%vert(:)
  indxv2(:) = elem(indx2)%vert(:)

  ! find sharing edge
  found = .false.
  do m1 = 1, 3
     m2 = mod(m1, 3) + 1
     n1 = indxv1(m1)
     n2 = indxv1(m2)

     iflg = 0
     do j = 1, 3
        if (indxv2(j) == n1) iflg = iflg + 1
        if (indxv2(j) == n2) iflg = iflg + 1
        if (indxv2(j) /= n1 .and. indxv2(j) /= n2) n3 = indxv2(j)
     end do
     if (iflg == 2) then
        mt1 = m1
        mt2 = m2
        nn3 = n3
        found = .true.
        exit  ! replaced goto 10
     end if
  end do

  nt1 = 4 - mt1
  nt2 = 4 - mt2
  nt3 = 6 - nt1 - nt2

  if (indxv2(nt1) /= indxv1(mt1) .or. indxv2(nt2) /= indxv1(mt2)) nc_corr = nc_corr + 1
  elem(indx2)%vert(nt1) = indxv1(mt1)
  elem(indx2)%vert(nt2) = indxv1(mt2)
  elem(indx2)%vert(nt3) = nn3

  return
end subroutine normal_correction_org


!---------------------------------------------------------------------------------------
!
!     Conn_Correction : Correction for wrong connectivity
!
!---------------------------------------------------------------------------------------
subroutine conn_correction(paramx, vert, elem, para)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertex), intent(inout) :: vert(paramx%nvert)
  type(element), intent(inout) :: elem(paramx%nelem)

  integer, allocatable :: iconn(:,:)
  integer :: nn, m1, m2, m3, mv1, mv2, n, nt, n1, n2, n3, nv1, nv2

  allocate(iconn(para%nelem, 3))
  iconn = 0
  do nn = 1, para%nelem
     do m1 = 1, 3
        m2 = mod(m1, 3) + 1
        m3 = 6 - m1 - m2
        mv1 = elem(nn)%vert(m1)
        mv2 = elem(nn)%vert(m2)

        do n = 1, 3
           nt = elem(nn)%conn(n)

           if (nt /= 0) then
              do n1 = 1, 3
                 n2 = mod(n1, 3) + 1
                 n3 = 6 - n1 - n2
                 nv1 = elem(nt)%vert(n1)
                 nv2 = elem(nt)%vert(n2)
                 if (mv1 == nv2 .and. mv2 == nv1) then
                    iconn(nn, m3) = nt
                    iconn(nt, n3) = nn
                 end if
              end do
           end if
        end do
     end do
  end do

  do nn = 1, para%nelem
     elem(nn)%conn(:) = iconn(nn, :)
  end do

  deallocate(iconn)

  return
end subroutine conn_correction


!---------------------------------------------------------------------------------------
!
!     Shared_Vertices : Find shared vertices
!
!---------------------------------------------------------------------------------------
subroutine shared_vertices(paramx, vert, elem, para)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertex), intent(inout) :: vert(paramx%nvert)
  type(element), intent(in) :: elem(paramx%nelem)

  ! local variable
  integer, allocatable :: nc(:)
  integer :: nn, mm, np

  allocate(nc(para%nvert))
  nc = 0
  do nn = 1, para%nelem
     do mm = 1, 3
        np = elem(nn)%vert(mm)
        nc(np) = nc(np) + 1
        vert(np)%share(nc(np) + 1) = nn
     end do
  end do

  do nn = 1, para%nvert
     vert(nn)%share(1) = nc(nn)  ! Total number of sharing points
     if (nc(nn) == 0) write(6, *) 'No sharing vertex', nn
  end do

  deallocate(nc)

  return
end subroutine shared_vertices


!---------------------------------------------------------------------------------------
!
!     Shared_Vertices_Org : Find shared vertices (original types)
!
!---------------------------------------------------------------------------------------
subroutine shared_vertices_org(paramx, vert, elem, para)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertexorg), intent(inout) :: vert(paramx%nvert)
  type(elementorg), intent(in) :: elem(paramx%nelem)

  ! local variable
  integer, allocatable, dimension(:) :: nc
  integer :: nn, mm, np

  allocate(nc(para%nvert))

  nc = 0
  do nn = 1, para%nelem
     do mm = 1, 3
        np = elem(nn)%vert(mm)
        nc(np) = nc(np) + 1
        vert(np)%share(nc(np) + 1) = nn
     end do
  end do

  do nn = 1, para%nvert
     vert(nn)%share(1) = nc(nn)  ! Total number of sharing points
     if (nc(nn) == 0) write(6, *) 'No sharing vertex', nn
  end do

  deallocate(nc)

  return
end subroutine shared_vertices_org


!---------------------------------------------------------------------------------------
!
!     Face_Normal : Calculate outward-pointing normal vectors at the face of elements
!
!---------------------------------------------------------------------------------------
subroutine face_normal(paramx, vert, elem, para)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertex), intent(inout) :: vert(paramx%nvert)
  type(element), intent(inout) :: elem(paramx%nelem)

  ! local variable
  integer :: np(3)
  real(8) :: xp1(3), xp2(3), xp3(3)
  real(8) :: onethird, xpnorm
  integer :: nn

  onethird = 1.0d0 / 3.0d0

  do nn = 1, para%nelem
     np(:) = elem(nn)%vert(:)
     xp2(:) = vert(np(3))%xyz(:) - vert(np(2))%xyz(:)
     xp1(:) = vert(np(2))%xyz(:) - vert(np(1))%xyz(:)

     xp3(1) = xp1(2) * xp2(3) - xp1(3) * xp2(2)
     xp3(2) = xp1(3) * xp2(1) - xp1(1) * xp2(3)
     xp3(3) = xp1(1) * xp2(2) - xp1(2) * xp2(1)

     xpnorm = sqrt(xp3(1)**2 + xp3(2)**2 + xp3(3)**2)
     elem(nn)%xyz(:) = onethird * (vert(np(1))%xyz(:) + vert(np(2))%xyz(:) + vert(np(3))%xyz(:))
     elem(nn)%norm(:) = xp3(:) / xpnorm
     elem(nn)%area = 0.5d0 * xpnorm
  end do

  return
end subroutine face_normal


!---------------------------------------------------------------------------------------
!
!     Find_Connectivity : Find connectivities of each elements
!
!---------------------------------------------------------------------------------------
subroutine find_connectivity(paramx, vert, elem, para)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertex), intent(in) :: vert(paramx%nvert)
  type(element), intent(inout) :: elem(paramx%nelem)

  integer, allocatable :: indx(:,:)
  integer :: np(3), mp(3)
  integer :: nn, mm, ll, k, nc, ip1, ip2

  allocate(indx(1000, 3))

  do nn = 1, para%nelem
     np(:) = elem(nn)%vert(:)
     mp(:) = vert(np(:))%share(1)
     do k = 1, 3
        do mm = 1, mp(k)
           indx(mm, k) = vert(np(k))%share(mm + 1)
        end do
     end do

     elem(nn)%conn(:) = 0
     nc = 0
     do ip1 = 1, 3
        ip2 = mod(ip1, 3) + 1
        do mm = 1, mp(ip1)
           if (indx(mm, ip1) /= nn) then
              do ll = 1, mp(ip2)
                 if (indx(mm, ip1) == indx(ll, ip2)) then
                    nc = nc + 1
                    elem(nn)%conn(nc) = indx(mm, ip1)
                 end if
              end do
           end if
        end do
     end do
     if (nc == 0) write(6, *) 'element w/o connection is found', nn
     if (nc > 3) write(6, *) 'wrong connection is found', nn
  end do

  deallocate(indx)

  return
end subroutine find_connectivity


!---------------------------------------------------------------------------------------
!
!     Find_Connectivity_Org : Find connectivities of each elements (original types)
!
!---------------------------------------------------------------------------------------
subroutine find_connectivity_org(paramx, vert, elem, para)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertexorg), intent(in) :: vert(paramx%nvert)
  type(elementorg), intent(inout) :: elem(paramx%nelem)

  integer, allocatable :: indx(:,:)
  integer :: np(3), mp(3)
  integer :: nn, mm, ll, k, nc, ip1, ip2

  allocate(indx(nshare, 3))

  do nn = 1, para%nelem
     np(:) = elem(nn)%vert(:)
     mp(:) = vert(np(:))%share(1)
     do k = 1, 3
        do mm = 1, mp(k)
           indx(mm, k) = vert(np(k))%share(mm + 1)
        end do
     end do

     elem(nn)%conn(:) = 0
     nc = 0
     do ip1 = 1, 3
        ip2 = mod(ip1, 3) + 1
        do mm = 1, mp(ip1)
           if (indx(mm, ip1) /= nn) then
              do ll = 1, mp(ip2)
                 if (indx(mm, ip1) == indx(ll, ip2)) then
                    nc = nc + 1
                    elem(nn)%conn(nc) = indx(mm, ip1)
                 end if
              end do
           end if
        end do
     end do
     if (nc == 0) write(6, *) 'element w/o connection is found', nn
     if (nc > 3) write(6, *) 'wrong connection is found', nn
  end do

  deallocate(indx)

  return
end subroutine find_connectivity_org


!---------------------------------------------------------------------------------------
!
!     Find_Segment : Find internal segments and correct order of vertices in elements
!
!---------------------------------------------------------------------------------------
subroutine find_segment(paramx, vert, elem, ielemseg, para)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertexorg), intent(in) :: vert(paramx%nvert)
  type(elementorg), intent(inout) :: elem(paramx%nelem)
  integer, intent(out) :: ielemseg(paramx%nelem)

  ! local variables
  integer, allocatable, dimension(:) :: ivisit, indx_list1, indx_list2
  integer :: nn, mm, nc, nc_corr, nseg, nlist1, nlist2
  integer :: indx, nt

  allocate(ivisit(para%nelem))
  allocate(indx_list1(para%nelem))
  allocate(indx_list2(para%nelem))

  ivisit = 0
  ielemseg = 0
  nseg = 0
  nc_corr = 0

  do nn = 1, para%nelem
     if (ivisit(nn) == 0) then
        nseg = nseg + 1
        ivisit(nn) = 1
        ielemseg(nn) = nseg

        ! BFS to find all connected elements
        ! Normal correction is done as each element is discovered (BFS order)
        nlist1 = 1
        indx_list1(1) = nn

        do while (nlist1 > 0)
           nlist2 = 0
           do mm = 1, nlist1
              indx = indx_list1(mm)
              do nc = 1, 3
                 nt = elem(indx)%conn(nc)
                 if (nt /= 0) then
                    if (ivisit(nt) == 0) then
                       ! Correct vertex ordering of newly discovered element
                       ! relative to its already-processed BFS parent
                       call normal_correction_org(paramx, vert, elem, para, indx, nt, nc_corr)
                       ivisit(nt) = 1
                       ielemseg(nt) = nseg
                       nlist2 = nlist2 + 1
                       indx_list2(nlist2) = nt
                    end if
                 end if
              end do
           end do
           nlist1 = nlist2
           indx_list1(1:nlist1) = indx_list2(1:nlist1)
        end do
     end if
  end do

  deallocate(ivisit, indx_list1, indx_list2)

  return
end subroutine find_segment


!---------------------------------------------------------------------------------------
!
!     Shared_Edge : Find the elements sharing edge
!
!---------------------------------------------------------------------------------------
subroutine shared_edge(paramx, vert, elem, edge, para)
  use GenSurf
  implicit none
  type(surfpar), intent(inout) :: paramx, para
  type(vertex), intent(in) :: vert(paramx%nvert)
  type(element), intent(inout) :: elem(paramx%nelem)
  type(edges), intent(inout) :: edge(paramx%nedge)

  integer, allocatable, dimension(:,:) :: ivisit
  integer :: nn, m1, m2, m3, n1, n2, nt, ne, mm, n3
  integer :: iopt_check

  allocate(ivisit(para%nelem, 3))

  ivisit = 0
  ne = 0
  do nn = 1, para%nelem
     do m1 = 1, 3
        m2 = mod(m1, 3) + 1
        m3 = 6 - m1 - m2
        n1 = elem(nn)%vert(m1)
        n2 = elem(nn)%vert(m2)
        nt = elem(nn)%conn(m3)

        if (nt /= 0) then
           do mm = 1, 3
              if (elem(nt)%conn(mm) == nn) n3 = mm
           end do

           if (ivisit(nn, m3) == 0) then
              ne = ne + 1
              edge(ne)%vert(1) = n1
              edge(ne)%vert(2) = n2
              edge(ne)%elem(1) = nn
              edge(ne)%elem(2) = nt

              ivisit(nn, m3) = 1
              ivisit(nt, n3) = 1
              elem(nn)%edge(m3) = ne
              elem(nt)%edge(n3) = ne
           end if

        else if (nt == 0) then

           if (ivisit(nn, m3) == 0) then
              ne = ne + 1
              edge(ne)%vert(1) = n1
              edge(ne)%vert(2) = n2
              edge(ne)%elem(1) = nn
              edge(ne)%elem(2) = nt

              ivisit(nn, m3) = 1
              elem(nn)%edge(m3) = ne
           end if

        end if

     end do
  end do

  para%nedgea = ne

  deallocate(ivisit)

  ! Check
  iopt_check = 0
  if (iopt_check == 1) then
     write(6, *) 'Total edge number', ne
     write(6, *) 'edge'
     do nn = 1, para%nedgea
        write(6, '(30I6)') nn, edge(nn)%elem(1), edge(nn)%elem(2), edge(nn)%vert(1), edge(nn)%vert(2)
        write(6, '(I6,30E15.7)') nn, vert(edge(nn)%vert(1))%xyz(1), vert(edge(nn)%vert(1))%xyz(2), &
             vert(edge(nn)%vert(1))%xyz(3), vert(edge(nn)%vert(2))%xyz(1), &
             vert(edge(nn)%vert(2))%xyz(2), vert(edge(nn)%vert(2))%xyz(3)
     end do
  end if

  return
end subroutine shared_edge


!---------------------------------------------------------------------------------------
!
!     Edge_Normal : Calculate outward-pointing pseudo-normal vectors at edges
!                   Treatment of thin body surface is included
!
!---------------------------------------------------------------------------------------
subroutine edge_normal(paramx, vert, elem, edge, para)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertex), intent(inout) :: vert(paramx%nvert)
  type(element), intent(in) :: elem(paramx%nelem)
  type(edges), intent(inout) :: edge(paramx%nedge)

  ! local variable
  real(8) :: xp1(3), xp2(3), xp3(3), xpm(3), xn1(3), xn2(3), xnp(3), xnm3(3)
  real(8) :: xnpmag, area, vecdot, xnsign
  integer :: ne, n1, n2, nn, nt, n0, mm, m3, ioptedge

  vert(:)%open = 0  ! connected vertex
  ioptedge = para%iedge

  do ne = 1, para%nedgea

     edge(ne)%open = 0  ! connected edge

     n1 = edge(ne)%vert(1)
     n2 = edge(ne)%vert(2)
     nn = edge(ne)%elem(1)
     nt = edge(ne)%elem(2)

     if (nt /= 0) then
        xp1(:) = vert(n1)%xyz(:)
        xp2(:) = vert(n2)%xyz(:)
        xpm(:) = 0.5d0 * (xp1(:) + xp2(:))

        xn1(:) = elem(nn)%norm(:)
        xn2(:) = elem(nt)%norm(:)
        xnp(:) = 0.5d0 * (xn1(:) + xn2(:))
        xnpmag = sqrt(xnp(1)**2 + xnp(2)**2 + xnp(3)**2)

        xnp(:) = xnp(:) / xnpmag
        area = 0.5d0 * (elem(nn)%area + elem(nt)%area)

        edge(ne)%xyz(:) = xpm(:)
        edge(ne)%norm(:) = xnp(:)
        edge(ne)%area = area

     else if (nt == 0) then

        xp1(:) = vert(n1)%xyz(:)
        xp2(:) = vert(n2)%xyz(:)
        xpm(:) = 0.5d0 * (xp1(:) + xp2(:))

        xn1(:) = elem(nn)%norm(:)  ! plane normal vector
        xn2(:) = xp1(:) - xp2(:)   ! edge vector

        call calc_normal_thin_edge(xn1, xn2, xnp)

        do mm = 1, 3
           m3 = elem(nn)%vert(mm)
           if (m3 /= n1 .and. m3 /= n2) n0 = m3  ! not shared points at the edge
        end do

        xp3(:) = vert(n0)%xyz(:)

        xnm3(:) = xpm(:) - xp3(:)
        vecdot = xnm3(1) * xnp(1) + xnm3(2) * xnp(2) + xnm3(3) * xnp(3)
        xnsign = sign(1.0d0, vecdot)

        xnp(:) = xnsign * xnp(:)

        ! edge normal vector back to face normal vector
        ! if(ioptedge==0) nothing - use tangential vector
        if (ioptedge == 2) xnp(:) = elem(nn)%norm(:) + xnp(:)  ! use plane normal + tangential vector
        if (ioptedge == 1) xnp(:) = elem(nn)%norm(:)            ! use plane normal vector

        xnpmag = sqrt(xnp(1)**2 + xnp(2)**2 + xnp(3)**2)

        xnp(:) = xnp(:) / xnpmag

        area = elem(nn)%area
        edge(ne)%xyz(:) = xpm(:)
        edge(ne)%norm(:) = xnp(:)
        edge(ne)%area = area
        edge(ne)%open = 1  ! disconnected edge

        vert(n1)%open = 1  ! disconnected vertex
        vert(n2)%open = 1  ! disconnected vertex

     end if

  end do

  return
end subroutine edge_normal


!---------------------------------------------------------------------------------------
!
!     Calc_Normal_Thin_Edge : Calculate outward-pointing pseudo-normal vectors at edges
!                             for thin body
!                             a*c =0, b*c=0, |c|=1, find c vector
!
!---------------------------------------------------------------------------------------
subroutine calc_normal_thin_edge(xn1, xn2, xn3)
  implicit none
  real(8), intent(in) :: xn1(3), xn2(3)
  real(8), intent(out) :: xn3(3)
  ! local variable
  real(8) :: xna(3)
  real(8) :: xnamx, c2, c3, c2a, c3a, camax, a1, a2, a3
  integer :: i, indx1, indx2, indx3

  xna(:) = abs(xn1(:))
  xnamx = max(max(xna(1), xna(2)), xna(3))

  do i = 1, 3
     if (xnamx == xna(i)) then
        indx1 = i
        if (indx1 == 1) indx2 = 2
        if (indx1 == 2) indx2 = 3
        if (indx1 == 3) indx2 = 1
     end if
  end do
  indx3 = 6 - indx1 - indx2

  c2 = -xn1(indx2) * xn2(indx1) / xn1(indx1) + xn2(indx2)
  c3 = -xn1(indx3) * xn2(indx1) / xn1(indx1) + xn2(indx3)

  c2a = abs(c2)
  c3a = abs(c3)
  camax = max(c2a, c3a)
  if (camax == c2a) then
     a2 = -c3 / c2
     a1 = -xn1(indx2) / xn1(indx1) * a2 - xn1(indx3) / xn1(indx1)
     xn3(indx3) = sqrt(1.0d0 / (a1**2 + a2**2 + 1.0d0))
     xn3(indx2) = a2 * xn3(indx3)
     xn3(indx1) = a1 * xn3(indx3)
  else
     a3 = -c2 / c3
     a1 = -xn1(indx2) / xn1(indx1) - xn1(indx3) / xn1(indx1) * a3
     xn3(indx2) = sqrt(1.0d0 / (a1**2 + a3**2 + 1.0d0))
     xn3(indx3) = a3 * xn3(indx2)
     xn3(indx1) = a1 * xn3(indx2)
  end if

  return
end subroutine calc_normal_thin_edge


!---------------------------------------------------------------------------------------
!
!     Pseudo_Normal : Calculate outward-pointing pseudo-normal vectors at vertices
!                     Angle-weighted pseudo-normal vectors
!
!---------------------------------------------------------------------------------------
subroutine pseudo_normal(paramx, vert, elem, edge, para)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertex), intent(inout) :: vert(paramx%nvert)
  type(element), intent(in) :: elem(paramx%nelem)
  type(edges), intent(in) :: edge(paramx%nedge)

  ! local variable
  integer :: np(3), nconn(3)
  real(8) :: xp1(3), xp2(3), pnorm(3)
  real(8) :: area, suma, x12dot, vecmag1, vecmag2, angle, pnormmag, xnormmag
  integer :: nn, mm, indx, ipp, ip, iconn, ne, ndetect, nconnp

  do nn = 1, para%nvert

     pnorm(:) = 0.0d0
     area = 0.0d0
     suma = 0.0d0

     do mm = 1, vert(nn)%share(1)

        ! angle computation
        indx = vert(nn)%share(mm + 1)
        np(:) = elem(indx)%vert(:)

        if (nn == np(1)) then
           xp1(:) = vert(np(2))%xyz(:) - vert(np(1))%xyz(:)
           xp2(:) = vert(np(3))%xyz(:) - vert(np(1))%xyz(:)
        else if (nn == np(2)) then
           xp1(:) = vert(np(1))%xyz(:) - vert(np(2))%xyz(:)
           xp2(:) = vert(np(3))%xyz(:) - vert(np(2))%xyz(:)
        else if (nn == np(3)) then
           xp1(:) = vert(np(1))%xyz(:) - vert(np(3))%xyz(:)
           xp2(:) = vert(np(2))%xyz(:) - vert(np(3))%xyz(:)
        end if

        x12dot = xp1(1) * xp2(1) + xp1(2) * xp2(2) + xp1(3) * xp2(3)
        vecmag1 = sqrt(xp1(1)**2 + xp1(2)**2 + xp1(3)**2)
        vecmag2 = sqrt(xp2(1)**2 + xp2(2)**2 + xp2(3)**2)
        angle = acos(x12dot / vecmag1 / vecmag2)

        pnorm(:) = pnorm(:) + angle * elem(indx)%norm(:)
        suma = suma + angle
        area = area + angle * elem(indx)%area
     end do

     pnorm(:) = pnorm(:) / suma
     pnormmag = sqrt(pnorm(1)**2 + pnorm(2)**2 + pnorm(3)**2)

     vert(nn)%norm(:) = pnorm(:) / pnormmag
     vert(nn)%area = area / suma

  end do

  ! re-define pseudo-normal vectors for disconnected edges
  if (para%nedgea /= para%nelem * 3 / 2) then  ! This must be checked!!!

     do nn = 1, para%nvert

        ! detect vertex involving disconnected edges
        ndetect = 1  ! indicator for disconnected edge, if ndetect=0
        do mm = 1, vert(nn)%share(1)
           indx = vert(nn)%share(mm + 1)
           np(:) = elem(indx)%vert(:)
           nconn(:) = elem(indx)%conn(:)
           ndetect = ndetect * nconn(1) * nconn(2) * nconn(3)
        end do

        if (ndetect == 0) then  ! do it again
           pnorm(:) = 0.0d0
           do mm = 1, vert(nn)%share(1)
              indx = vert(nn)%share(mm + 1)
              np(:) = elem(indx)%vert(:)
              nconn(:) = elem(indx)%conn(:)
              nconnp = nconn(1) * nconn(2) * nconn(3)

              if (nconnp == 0) then
                 do ip = 1, 3
                    if (nn == np(ip)) ipp = ip
                 end do
                 do iconn = 1, 3  ! find edge information
                    if (nconn(iconn) == 0 .and. ipp /= iconn) then
                       ne = elem(indx)%edge(iconn)
                       pnorm(:) = pnorm(:) + edge(ne)%norm(:)
                    end if
                 end do
              end if
           end do

           xnormmag = sqrt(pnorm(1)**2 + pnorm(2)**2 + pnorm(3)**2)
           vert(nn)%norm(:) = pnorm(:) / xnormmag
        end if
     end do  ! end of vertex

  end if

  return
end subroutine pseudo_normal


!---------------------------------------------------------------------------------------
!
!     Define_Element : Read ASCII STL Files and define # of vertices and elements
!
!---------------------------------------------------------------------------------------
subroutine define_element(filestl, nelem)
  implicit none
  character(len=80), intent(in) :: filestl
  integer, intent(out) :: nelem

  character(len=80) :: cname
  integer :: line

  ! Check number of triangle elements
  line = 0
  open(11, file=trim(filestl))
  do
     read(11, *) cname
     line = line + 1
     if (trim(cname) == 'endsolid') exit  ! replaced goto 1
  end do
  close(11)

  ! Read STL data file and define domain size
  nelem = (line - 2) / 7

  return
end subroutine define_element


!---------------------------------------------------------------------------------------
!
!     Read_STL : Read ASCII STL Files and store positions of vertices in datapts
!
!---------------------------------------------------------------------------------------
subroutine read_stl(filestl, datapts, nelemorg, nvertorg)
  implicit none
  integer, intent(in) :: nelemorg, nvertorg
  real(8), intent(out) :: datapts(nvertorg, 3)
  character(len=80), intent(in) :: filestl

  character(len=80) :: cname
  real(8) :: dummy1, dummy2, dummy3
  integer :: nn, i

  open(11, file=trim(filestl))
  read(11, *)
  nn = 0
  do i = 1, nelemorg
     read(11, *) cname, cname, dummy1, dummy2, dummy3
     read(11, *) cname
     nn = nn + 1
     read(11, *) cname, datapts(nn, 1), datapts(nn, 2), datapts(nn, 3)
     nn = nn + 1
     read(11, *) cname, datapts(nn, 1), datapts(nn, 2), datapts(nn, 3)
     nn = nn + 1
     read(11, *) cname, datapts(nn, 1), datapts(nn, 2), datapts(nn, 3)
     read(11, *) cname
     read(11, *) cname
  end do
  close(11)

  return
end subroutine read_stl


!---------------------------------------------------------------------------------------
!
!     Find_Duplications : Find duplicated vertices and elements
!
!---------------------------------------------------------------------------------------
subroutine find_duplications(datapts, indxmapping, ielemmapping, paramx, para, nudvert, nudelem)
  use GenSurf
  implicit none
  type(surfpar), intent(inout) :: paramx, para
  integer, intent(out) :: indxmapping(paramx%nvert, 2)
  integer, intent(out) :: ielemmapping(paramx%nelem)
  real(8), intent(in) :: datapts(para%nvert, 3)
  integer, intent(out) :: nudvert, nudelem

  ! local variable
  integer, allocatable, dimension(:,:) :: elem2vert, vert2share, indxshare
  integer, allocatable, dimension(:) :: ncshare
  integer, allocatable, dimension(:) :: indxarray, indxarray2
  real(8), allocatable, dimension(:) :: distarray, distarray2
  integer, allocatable, dimension(:) :: indxvisit, idupelem
  integer, allocatable, dimension(:,:) :: idup

  integer :: np(3), mp(3)
  real(8) :: qp(3)
  integer :: nsmx, nsmx2
  integer :: nn, mm, ll, ncount, npl, indxtag, indxdup
  integer :: k, m1, m2, m3

  allocate(elem2vert(para%nelem, 3))
  allocate(vert2share(para%nvert, nshare))
  allocate(ncshare(para%nvert))
  allocate(indxshare(para%nelem, 3))
  allocate(indxvisit(para%nvert))
  allocate(idupelem(para%nelem))
  allocate(idup(para%nelem, 3))

  ! set up search space
  nsmx = min(50, para%nvert)
  nsmx2 = para%nvert  ! entire search

  allocate(indxarray(nsmx))
  allocate(distarray(nsmx))
  allocate(indxarray2(nsmx2))
  allocate(distarray2(nsmx2))

  call wldst_create_tree(datapts(1, 1), datapts(1, 2), datapts(1, 3), para%nvert)

  indxvisit = 0
  ncount = 0

  do nn = 1, para%nvert
     if (indxvisit(nn) == 0) then
        qp(:) = datapts(nn, :)
        ncount = ncount + 1
        indxvisit(nn) = ncount
        call wldst_getdists_xyz(qp(1), qp(2), qp(3), indxarray, distarray, nsmx)
        if (distarray(nsmx) /= 0.0d0) then
           do mm = 1, nsmx
              if (distarray(mm) == 0.0d0) indxvisit(indxarray(mm) + 1) = ncount
           end do
        else  ! extended search
           if (nsmx /= nsmx2) then
              call wldst_getdists_xyz(qp(1), qp(2), qp(3), indxarray2, distarray2, nsmx2)
              do mm = 1, nsmx2
                 if (distarray2(mm) == 0.0d0) indxvisit(indxarray2(mm) + 1) = ncount
              end do
           end if
        end if
     end if
  end do

  call wldst_delete_tree()

  deallocate(indxarray, distarray, indxarray2, distarray2)

  ! Find data mapping between un-duplicate vertices and original data
  nudvert = ncount

  ! indxmapping(:,1) : indx from new set to original set
  ! indxmapping(:,2) : indx from original set to new set
  indxmapping = para%nvert
  do nn = 1, para%nvert
     indxtag = indxvisit(nn)
     indxmapping(indxtag, 1) = min(indxmapping(indxtag, 1), nn)
     indxmapping(nn, 2) = indxtag
  end do

  ! Check duplicate elements

  ! Define shared vertices
  ncshare = 0
  do nn = 1, para%nelem
     mm = 3 * (nn - 1)
     elem2vert(nn, 1) = indxmapping(mm + 1, 2)
     elem2vert(nn, 2) = indxmapping(mm + 2, 2)
     elem2vert(nn, 3) = indxmapping(mm + 3, 2)
     do ll = 1, 3
        npl = elem2vert(nn, ll)
        ncshare(npl) = ncshare(npl) + 1
        vert2share(npl, ncshare(npl) + 1) = nn
     end do
  end do

  do nn = 1, nudvert
     vert2share(nn, 1) = ncshare(nn)
  end do

  idupelem = 0
  do nn = 1, para%nelem
     np(:) = elem2vert(nn, :)
     mp(:) = vert2share(np(:), 1)
     do k = 1, 3
        do mm = 1, mp(k)
           indxshare(mm, k) = vert2share(np(k), mm + 1)
        end do
     end do
     do m1 = 1, mp(1)
        do m2 = 1, mp(2)
           indxdup = indxshare(m1, 1)
           if (indxdup /= nn .and. indxshare(m2, 2) == indxdup) then
              do m3 = 1, mp(3)
                 if (indxshare(m3, 3) == indxdup) idupelem(indxdup) = 1
              end do
           end if
        end do
     end do
  end do

  ncount = 0
  ielemmapping = 0
  do nn = 1, para%nelem
     if (idupelem(nn) == 0) then
        ncount = ncount + 1
        ielemmapping(nn) = ncount
     end if
  end do
  nudelem = ncount

  deallocate(elem2vert, vert2share, ncshare, indxshare, indxvisit, idupelem, idup)

  return
end subroutine find_duplications


!---------------------------------------------------------------------------------------
!
!     Convert2IBSurface : Read ASCII STL Files
!                         Find Duplicated Vertices & Elements
!                         Find Shared Vertices and Edges
!                         Find Connectivity
!                         Find Segments
!                         Write IBSurface in Tecplot, Binary, new STL, Nodesurf formats
!
!---------------------------------------------------------------------------------------
subroutine convert2ibsurface(filestl, paramx, vert, elem, para, indxmapping, ielemmapping, ielemseg, nlsf)
  use GenSurf
  implicit none
  type(surfpar), intent(inout) :: paramx, para
  type(vertexorg), intent(inout) :: vert(paramx%nvert)
  type(elementorg), intent(inout) :: elem(paramx%nelem)
  integer, intent(inout) :: indxmapping(paramx%nvert, 2)
  integer, intent(inout) :: ielemmapping(paramx%nelem)
  integer, intent(inout) :: ielemseg(paramx%nelem)
  integer, intent(out) :: nlsf
  character(len=80), intent(in) :: filestl

  real(8), allocatable, dimension(:,:) :: datapts
  integer :: nelemorg, nvertorg, nudvert, nudelem
  integer :: nn, mm, maxshare, ielemmap, indx

  nelemorg = para%nelem
  nvertorg = para%nvert

  allocate(datapts(para%nvert, 3))

  ! Read STL data
  call read_stl(filestl, datapts, nelemorg, nvertorg)
  write(6, *) 'Original parameters, #vert, #elem', nvertorg, nelemorg

  ! Find duplication of vertices and index mapping
  call find_duplications(datapts, indxmapping, ielemmapping, paramx, para, nudvert, nudelem)

  ! Make data structure of triangle elements
  para%nvert = nudvert  ! re-define para%nvert using un-duplicate vertices
  para%nelem = nudelem  ! re-define para%nelem using un-duplicate elements
  para%nedge = nudelem * 3 / 2
  write(6, *) 'Undupliate #vert, #elem', para%nvert, para%nelem

  do nn = 1, para%nvert
     indx = indxmapping(nn, 1)
     vert(nn)%xyz(:) = datapts(indx, :)
  end do

  do nn = 1, nelemorg
     ielemmap = ielemmapping(nn)
     if (ielemmap /= 0) then
        mm = 3 * (nn - 1)
        elem(ielemmap)%vert(1) = indxmapping(mm + 1, 2)
        elem(ielemmap)%vert(2) = indxmapping(mm + 2, 2)
        elem(ielemmap)%vert(3) = indxmapping(mm + 3, 2)
     end if
  end do

  ! Find share vertices
  call shared_vertices_org(paramx, vert, elem, para)
  maxshare = 1
  do nn = 1, para%nvert
     maxshare = max(maxshare, vert(nn)%share(1))
  end do
  write(6, *) 'Maximum shared vertices', maxshare

  ! Find connectivity
  call find_connectivity_org(paramx, vert, elem, para)

  ! Find segment
  call find_segment(paramx, vert, elem, ielemseg, para)
  do nn = 1, para%nelem
     nlsf = max(1, ielemseg(nn))
     if (ielemseg(nn) == 0) then
        write(6, *) 'Un-assigned elements are detected!!', nn
     end if
  end do
  write(6, *) 'Number of levelsets', nlsf

  deallocate(datapts)

  return
end subroutine convert2ibsurface


!---------------------------------------------------------------------------------------
!
!     IBSurface_Seg_Param : Define segmented ibsurface parameters
!
!---------------------------------------------------------------------------------------
subroutine ibsurface_seg_param(paramx, vertorg, elemorg, paraorg, &
     indxmapping, ielemmapping, ielemseg, ncelem, ncvert, &
     ielem_seg2org, ivert_seg2org, ivert_org2seg, nls, nlsmap, nlsf, nlsfile)
  use GenSurf
  implicit none
  integer, intent(in) :: nls, nlsfile
  type(surfpar), intent(in) :: paramx, paraorg(nlsfile)
  type(vertexorg), intent(in) :: vertorg(paramx%nvert, nlsfile)
  type(elementorg), intent(in) :: elemorg(paramx%nelem, nlsfile)
  integer, intent(in) :: indxmapping(paramx%nvert, 2, nlsfile)
  integer, intent(in) :: ielemmapping(paramx%nelem, nlsfile)
  integer, intent(in) :: ielemseg(paramx%nelem, nlsfile)
  integer, intent(in) :: nlsmap(nls), nlsf(nlsfile)
  integer, intent(out) :: ncelem(nls), ncvert(nls)
  integer, intent(out) :: ielem_seg2org(paramx%nelem, nls)
  integer, intent(out) :: ivert_seg2org(paramx%nvert, nls)
  integer, intent(out) :: ivert_org2seg(paramx%nvert, nlsfile)

  integer, allocatable, dimension(:,:) :: ivisit
  integer :: nlsindex, nf, nn, mm, iseg, nv

  allocate(ivisit(paramx%nvert, nls))

  ! re-define elements and vertices in each segment
  nlsindex = 0
  ncelem = 0
  ncvert = 0
  ivisit = 0

  do nf = 1, nlsfile

     if (nf /= 1) nlsindex = nlsindex + nlsf(nf - 1)

     do nn = 1, paraorg(nf)%nelem
        iseg = ielemseg(nn, nf) + nlsindex
        ncelem(iseg) = ncelem(iseg) + 1
        ielem_seg2org(ncelem(iseg), iseg) = nn
        do mm = 1, 3
           nv = elemorg(nn, nf)%vert(mm)
           if (ivisit(nv, iseg) == 0) then
              ivisit(nv, iseg) = 1
              ncvert(iseg) = ncvert(iseg) + 1
              ivert_seg2org(ncvert(iseg), iseg) = nv
              ivert_org2seg(nv, nf) = ncvert(iseg)
           end if
        end do
     end do

  end do

  deallocate(ivisit)

  return
end subroutine ibsurface_seg_param


!---------------------------------------------------------------------------------------
!
!     IBSurface_Seg : Define segmented ibsurface parameters (vert,elem,edge)
!
!---------------------------------------------------------------------------------------
subroutine ibsurface_seg(paraorgmx, vertorg, elemorg, paraorg, &
     paramx, vert, elem, edge, para, &
     ielem_seg2org, ivert_seg2org, ivert_org2seg, nls, nlsmap, nlsf, nlsfile)
  use GenSurf
  implicit none
  integer, intent(in) :: nls, nlsfile
  type(surfpar), intent(in) :: paraorgmx, paraorg(nlsfile)
  type(vertexorg), intent(in) :: vertorg(paraorgmx%nvert, nlsfile)
  type(elementorg), intent(in) :: elemorg(paraorgmx%nelem, nlsfile)
  integer, intent(in) :: ielem_seg2org(paraorgmx%nelem, nls)
  integer, intent(in) :: ivert_seg2org(paraorgmx%nvert, nls)
  integer, intent(in) :: ivert_org2seg(paraorgmx%nvert, nlsfile)

  ! Segmented objects
  type(surfpar), intent(inout) :: paramx, para(nls)
  type(vertex), intent(inout) :: vert(paramx%nvert, nls)
  type(element), intent(inout) :: elem(paramx%nelem, nls)
  type(edges), intent(inout) :: edge(paramx%nedge, nls)
  integer, intent(in) :: nlsmap(nls), nlsf(nlsfile)

  integer :: np(3), mp(3)
  integer :: ll, nn, mm, nf

  do ll = 1, nls
     nf = nlsmap(ll)
     do nn = 1, para(ll)%nvert
        mm = ivert_seg2org(nn, ll)
        vert(nn, ll)%xyz(:) = vertorg(mm, nf)%xyz(:)
     end do
     do nn = 1, para(ll)%nelem
        mm = ielem_seg2org(nn, ll)
        np(:) = elemorg(mm, nf)%vert(:)
        mp(:) = ivert_org2seg(np(:), nf)
        elem(nn, ll)%vert(:) = mp(:)
     end do
  end do

  ! All calculations for the elements (shared vertices/edges, face-normal,edge-normal,pseudo-normal)
  ! Check the normal vectors

  do ll = 1, nls
     para(ll)%iedge = 2  ! initialize before edge_normal (default value from output_bin)
     call shared_vertices(paramx, vert(1, ll), elem(1, ll), para(ll))
     call find_connectivity(paramx, vert(1, ll), elem(1, ll), para(ll))
     call conn_correction(paramx, vert(1, ll), elem(1, ll), para(ll))
     call shared_edge(paramx, vert(1, ll), elem(1, ll), edge(1, ll), para(ll))
     call face_normal(paramx, vert(1, ll), elem(1, ll), para(ll))
     call edge_normal(paramx, vert(1, ll), elem(1, ll), edge(1, ll), para(ll))
     call pseudo_normal(paramx, vert(1, ll), elem(1, ll), edge(1, ll), para(ll))
  end do

  return
end subroutine ibsurface_seg


!---------------------------------------------------------------------------------------
!
!     Scale_IB_Data : Scaling Data or Shift IB Position
!
!---------------------------------------------------------------------------------------
subroutine scale_ibsurf(paramx, vert, para, scale)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertex), intent(inout) :: vert(paramx%nvert)
  real(8), intent(in) :: scale

  integer :: nn
  real(8) :: xmax, xmin, ymax, ymin, zmax, zmin, xctr, yctr, zctr

  do nn = 1, para%nvert
     vert(nn)%xyz(1) = scale * vert(nn)%xyz(1)
     vert(nn)%xyz(2) = scale * vert(nn)%xyz(2)
     vert(nn)%xyz(3) = scale * vert(nn)%xyz(3)
  end do

  xmax = -1.0d30;  xmin = 1.0d30
  ymax = -1.0d30;  ymin = 1.0d30
  zmax = -1.0d30;  zmin = 1.0d30

  do nn = 1, para%nvert
     xmax = max(xmax, vert(nn)%xyz(1))
     xmin = min(xmin, vert(nn)%xyz(1))
     ymax = max(ymax, vert(nn)%xyz(2))
     ymin = min(ymin, vert(nn)%xyz(2))
     zmax = max(zmax, vert(nn)%xyz(3))
     zmin = min(zmin, vert(nn)%xyz(3))
  end do

  xctr = 0.5d0 * (xmax + xmin)
  yctr = 0.5d0 * (ymax + ymin)
  zctr = 0.5d0 * (zmax + zmin)

  return
end subroutine scale_ibsurf


!---------------------------------------------------------------------------------------
!
!     Read_ibsurf_bin : Read IBSurface in binary Format w or w/o elem_vert
!
!---------------------------------------------------------------------------------------
subroutine read_ibsurf_bin(filein, paramx, vert, elem, para, ibseq, nls, nseq, ibscale, ibfactor, xyzoff)
  use GenSurf
  implicit none
  integer, intent(in) :: nls, nseq
  type(surfpar), intent(in) :: paramx
  type(surfpar), intent(inout) :: para(nls)
  type(vertex), intent(inout) :: vert(paramx%nvert, nls)
  type(element), intent(inout) :: elem(paramx%nelem, nls)
  integer, intent(in) :: ibseq(nls), ibscale(nls)
  real(8), intent(in) :: ibfactor(nls), xyzoff(nls, 3)
  character(len=80), intent(in) :: filein

  character(len=80) :: fout, foutl
  integer :: nn, ll, nls_read
  integer :: ibcheck
  real(8) :: xmax, xmin, ymax, ymin, zmax, zmin, xctr, yctr, zctr

  write(foutl, '(I10)') nseq
  fout = trim(filein) // trim(adjustl(foutl)) // '.bin'

  open(10, file=trim(fout), form='unformatted')

  if (nseq == 0) then
     read(10) nls_read
     do ll = 1, nls
        read(10) para(ll)
     end do

     do ll = 1, nls
        read(10) (vert(nn, ll)%xyz, nn = 1, para(ll)%nvert)
     end do
     do ll = 1, nls
        read(10) (elem(nn, ll)%vert, nn = 1, para(ll)%nelem)
     end do

  else

     do ll = 1, nls
        if (ibseq(ll) == 1) then
           read(10) (vert(nn, ll)%xyz, nn = 1, para(ll)%nvert)
        end if
     end do

  end if
  close(10)

  ! scaling and translation
  do ll = 1, nls
     if (ibscale(ll) == 1 .and. (ibseq(ll) == 1 .or. nseq == 0)) then
        do nn = 1, para(ll)%nvert
           vert(nn, ll)%xyz(:) = ibfactor(ll) * vert(nn, ll)%xyz(:) + xyzoff(ll, :)
        end do
     end if
  end do

  ibcheck = 1

  if (ibcheck == 1) then
     do ll = 1, nls
        if (ibseq(ll) == 1 .or. nseq == 0) then

           xmax = -1.0d30;  xmin = 1.0d30
           ymax = -1.0d30;  ymin = 1.0d30
           zmax = -1.0d30;  zmin = 1.0d30

           do nn = 1, para(ll)%nvert
              xmax = max(xmax, vert(nn, ll)%xyz(1))
              xmin = min(xmin, vert(nn, ll)%xyz(1))
              ymax = max(ymax, vert(nn, ll)%xyz(2))
              ymin = min(ymin, vert(nn, ll)%xyz(2))
              zmax = max(zmax, vert(nn, ll)%xyz(3))
              zmin = min(zmin, vert(nn, ll)%xyz(3))
           end do

           xctr = 0.5d0 * (xmax + xmin)
           yctr = 0.5d0 * (ymax + ymin)
           zctr = 0.5d0 * (zmax + zmin)

        end if
     end do
  end if

  return
end subroutine read_ibsurf_bin


!---------------------------------------------------------------------------------------
!
!     ibsurf_conn : Calculate connectivity, shared_vertices, etc.
!
!---------------------------------------------------------------------------------------
subroutine ibsurf_conn(paramx, vert, elem, edge, para, ibseq, nls, nseq)
  use GenSurf
  implicit none
  integer, intent(in) :: nls, nseq
  type(surfpar), intent(inout) :: paramx
  type(surfpar), intent(inout) :: para(nls)
  type(vertex), intent(inout) :: vert(paramx%nvert, nls)
  type(element), intent(inout) :: elem(paramx%nelem, nls)
  type(edges), intent(inout) :: edge(paramx%nedge, nls)
  integer, intent(in) :: ibseq(nls)

  real(8) :: powervalue
  integer :: ll

  powervalue = 1.0d0 / 7.0d0

  if (nseq == 0) then
     do ll = 1, nls
        call shared_vertices(paramx, vert(1, ll), elem(1, ll), para(ll))
        call find_connectivity(paramx, vert(1, ll), elem(1, ll), para(ll))
        call conn_correction(paramx, vert(1, ll), elem(1, ll), para(ll))
        call shared_edge(paramx, vert(1, ll), elem(1, ll), edge(1, ll), para(ll))
        call face_normal(paramx, vert(1, ll), elem(1, ll), para(ll))
        call edge_normal(paramx, vert(1, ll), elem(1, ll), edge(1, ll), para(ll))
        call pseudo_normal(paramx, vert(1, ll), elem(1, ll), edge(1, ll), para(ll))

        vert(:, :)%power = powervalue
        elem(:, :)%power = powervalue
        edge(:, :)%power = powervalue
     end do

  else

     do ll = 1, nls
        if (ibseq(ll) == 1) then
           call face_normal(paramx, vert(1, ll), elem(1, ll), para(ll))
           call edge_normal(paramx, vert(1, ll), elem(1, ll), edge(1, ll), para(ll))
           call pseudo_normal(paramx, vert(1, ll), elem(1, ll), edge(1, ll), para(ll))
        end if
     end do

  end if

  return
end subroutine ibsurf_conn


!---------------------------------------------------------------------------------------
!
!     calc_ibvel : calc_ibvel and ibacc
!
!---------------------------------------------------------------------------------------
subroutine calc_ibvel(paramx, vert, elem, edge, para, ibseq, nls, nseq, dtadvec, velmax, accmax, posmax, posmin)
  use GenSurf
  implicit none
  integer, intent(in) :: nls, nseq
  type(surfpar), intent(in) :: paramx, para(nls)
  type(vertex), intent(inout) :: vert(paramx%nvert, nls)
  type(element), intent(in) :: elem(paramx%nelem, nls)
  type(edges), intent(in) :: edge(paramx%nedge, nls)
  integer, intent(in) :: ibseq(nls)
  real(8), intent(in) :: dtadvec
  real(8), intent(out) :: velmax(nls), accmax(nls), posmax(nls, 3), posmin(nls, 3)

  integer :: ll, nn
  real(8) :: velmag, accmag

  posmin = 1.0d30
  posmax = -1.0d30
  velmax = 0.0d0
  accmax = 0.0d0

  if (nseq == 0) then

     do ll = 1, nls
        do nn = 1, para(ll)%nvert
           vert(nn, ll)%xyz0 = vert(nn, ll)%xyz
           vert(nn, ll)%vel = 0.0d0
           vert(nn, ll)%vel0 = vert(nn, ll)%vel
           vert(nn, ll)%acc = 0.0d0

           velmag = sqrt(vert(nn, ll)%vel(1)**2 + vert(nn, ll)%vel(2)**2 + vert(nn, ll)%vel(3)**2)
           accmag = sqrt(vert(nn, ll)%acc(1)**2 + vert(nn, ll)%acc(2)**2 + vert(nn, ll)%acc(3)**2)
           velmax(ll) = max(velmax(ll), velmag)
           accmax(ll) = max(accmax(ll), accmag)
           posmax(ll, 1:3) = max(posmax(ll, 1:3), vert(nn, ll)%xyz(1:3))
           posmin(ll, 1:3) = min(posmin(ll, 1:3), vert(nn, ll)%xyz(1:3))
        end do
     end do

  else

     do ll = 1, nls
        do nn = 1, para(ll)%nvert
           vert(nn, ll)%vel = (vert(nn, ll)%xyz - vert(nn, ll)%xyz0) / dtadvec
           vert(nn, ll)%acc = (vert(nn, ll)%vel - vert(nn, ll)%vel0) / dtadvec
           vert(nn, ll)%vel0 = vert(nn, ll)%vel
           vert(nn, ll)%xyz0 = vert(nn, ll)%xyz

           velmag = sqrt(vert(nn, ll)%vel(1)**2 + vert(nn, ll)%vel(2)**2 + vert(nn, ll)%vel(3)**2)
           accmag = sqrt(vert(nn, ll)%acc(1)**2 + vert(nn, ll)%acc(2)**2 + vert(nn, ll)%acc(3)**2)

           velmax(ll) = max(velmax(ll), velmag)
           accmax(ll) = max(accmax(ll), accmag)
           posmax(ll, 1:3) = max(posmax(ll, 1:3), vert(nn, ll)%xyz(1:3))
           posmin(ll, 1:3) = min(posmin(ll, 1:3), vert(nn, ll)%xyz(1:3))
        end do
     end do

     write(6, *) 'Maximum velocity', velmax
     write(6, *) 'Maximum acceleration', accmax

  end if

  return
end subroutine calc_ibvel


!---------------------------------------------------------------------------------------
!
!     IBSurface_Sequence : Read ASCII STL Files
!                          Find Duplicated Vertices & Elements
!
!---------------------------------------------------------------------------------------
subroutine ibsurface_sequence(fileseq, paraorgmx, paraorg, indxmapping, ivert_seg2org, nlsf, nlsfseq, nlsfile, &
     paramx, vert, elem, edge, para, nls, ibseq, nseq, nlsfscale, fscalenls, icheck)
  use GenSurf
  implicit none
  integer, intent(in) :: nlsfile, nls, nseq
  type(surfpar), intent(in) :: paraorgmx, paraorg(nlsfile)
  type(surfpar), intent(in) :: paramx, para(nls)
  type(vertex), intent(inout) :: vert(paramx%nvert, nls)
  type(element), intent(inout) :: elem(paramx%nelem, nls)
  type(edges), intent(inout) :: edge(paramx%nedge, nls)
  integer, intent(in) :: indxmapping(paraorgmx%nvert, 2, nlsfile)
  integer, intent(in) :: ivert_seg2org(paraorgmx%nvert, nls)
  integer, intent(in) :: nlsf(nlsfile), nlsfseq(nlsfile), nlsfscale(nlsfile), ibseq(nls)
  real(8), intent(in) :: fscalenls(nlsfile)
  integer, intent(in) :: icheck
  character(len=80), intent(in) :: fileseq(nlsfile)

  integer :: nlsmn(nlsfile), nlsmx(nlsfile)
  real(8), allocatable, dimension(:,:) :: datapts
  character(len=80) :: filestl, foutl
  integer :: nf, ll, nn, mm, indx

  allocate(datapts(paraorgmx%nvert, 3))

  nlsmn(1) = 1
  nlsmx(1) = nlsf(1)

  do nf = 2, nlsfile
     nlsmn(nf) = nlsmn(nf - 1) + 1
     nlsmx(nf) = nlsmx(nf - 1) + nlsf(nf)
  end do

  do nf = 1, nlsfile

     if (nlsfseq(nf) == 1) then

        ! Read STL data
        write(foutl, '(I10)') nseq
        filestl = trim(fileseq(nf)) // trim(adjustl(foutl)) // '.stl'

        call read_stl(filestl, datapts, paraorg(nf)%nelem, paraorgmx%nvert)
        write(6, *) 'Original parameters, #seq #vert,#elem', nseq, paraorg(nf)%nelem, paraorg(nf)%nvert

        do ll = nlsmn(nf), nlsmx(nf)
           do nn = 1, para(ll)%nvert
              mm = ivert_seg2org(nn, ll)
              indx = indxmapping(mm, 1, nf)
              vert(nn, ll)%xyz(:) = datapts(indx, :)
           end do
        end do

        if (nlsfscale(nf) == 1) then
           do ll = nlsmn(nf), nlsmx(nf)
              call scale_ibsurf(paramx, vert(1, ll), para(ll), fscalenls(nf))
           end do
        end if

     end if  ! end of nlsfseq

  end do  ! end of nf

  if (icheck == 1) call ibsurf_conn(paramx, vert, elem, edge, para, ibseq, nls, nseq)

  deallocate(datapts)

  return
end subroutine ibsurface_sequence


!---------------------------------------------------------------------------------------
!
!     Output_STL : Write IBSurface in STL Format
!
!---------------------------------------------------------------------------------------
subroutine output_stl(paramx, vert, elem, para, l)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertex), intent(in) :: vert(paramx%nvert)
  type(element), intent(in) :: elem(paramx%nelem)
  integer, intent(in) :: l

  ! local variable
  integer :: np(3)
  character(len=80) :: fout, foutl
  integer :: nn

  write(foutl, '(I10)') l
  fout = 'output/Output_STL' // trim(adjustl(foutl)) // '.stl'

  open(10, file=trim(fout))
  write(10, '(A)') 'solid object'

  do nn = 1, para%nelem
     np(:) = elem(nn)%vert(:)
     write(10, '(A,3E15.7)') '  facet normal ', elem(nn)%norm(1), elem(nn)%norm(2), elem(nn)%norm(3)
     write(10, '(A)') '    outer loop'
     write(10, '(A,3E15.7)') '      vertex ', vert(np(1))%xyz(1), vert(np(1))%xyz(2), vert(np(1))%xyz(3)
     write(10, '(A,3E15.7)') '      vertex ', vert(np(2))%xyz(1), vert(np(2))%xyz(2), vert(np(2))%xyz(3)
     write(10, '(A,3E15.7)') '      vertex ', vert(np(3))%xyz(1), vert(np(3))%xyz(2), vert(np(3))%xyz(3)
     write(10, '(A)') '    endloop'
     write(10, '(A)') '  endfacet'
  end do
  write(10, '(A)') 'endsolid object'
  close(10)

  return
end subroutine output_stl


!---------------------------------------------------------------------------------------
!
!     Output_Node : Write IBSurface in Scatter Format
!
!---------------------------------------------------------------------------------------
subroutine output_node(paramx, vert, elem, edge, para, l)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertex), intent(in) :: vert(paramx%nvert)
  type(element), intent(in) :: elem(paramx%nelem)
  type(edges), intent(in) :: edge(paramx%nedge)
  integer, intent(in) :: l

  ! local variables
  integer :: np(3)
  real(8) :: xg(3), xn(3)
  character(len=80) :: fout, foutl
  real(8) :: onethird
  integer :: nn

  onethird = 1.0d0 / 3.0d0

  write(foutl, '(I10)') l
  fout = 'output/Output_Node' // trim(adjustl(foutl)) // '.dat'

  open(20, file=trim(fout))

  write(20, *) 'zone t="T', para%nelem, '"'
  do nn = 1, para%nelem
     np(:) = elem(nn)%vert(:)
     xg(:) = onethird * (vert(np(1))%xyz(:) + vert(np(2))%xyz(:) + vert(np(3))%xyz(:))
     xn(:) = elem(nn)%norm(:)
     write(20, '(30E15.7)') xg, xn
  end do
  write(20, *) 'zone t="V', para%nvert, '"'
  do nn = 1, para%nvert
     xg(:) = vert(nn)%xyz(:)
     xn(:) = vert(nn)%norm(:)
     write(20, '(30E15.7)') xg, xn
  end do
  write(20, *) 'zone t="E', para%nedgea, '"'
  do nn = 1, para%nedgea
     xg(:) = edge(nn)%xyz(:)
     xn(:) = edge(nn)%norm(:)
     write(20, '(30E15.7)') xg, xn
  end do
  close(20)

  return
end subroutine output_node


!---------------------------------------------------------------------------------------
!
!     Output_Custom : Write IBSurface in Custom format
!
!---------------------------------------------------------------------------------------
subroutine output_custom(paramx, vert, elem, edge, para, l)
  use GenSurf
  implicit none
  type(surfpar), intent(in) :: paramx, para
  type(vertex), intent(in) :: vert(paramx%nvert)
  type(element), intent(in) :: elem(paramx%nelem)
  type(edges), intent(in) :: edge(paramx%nedge)
  integer, intent(in) :: l

  ! local variables
  integer :: np(3)
  real(8) :: xg(3), xn(3)
  character(len=80) :: fout, foutl
  real(8) :: onethird
  integer :: nn, mm

  onethird = 1.0d0 / 3.0d0

  write(foutl, '(I10)') l
  fout = 'output/Output' // trim(adjustl(foutl)) // '.dat'

  open(10, file=trim(fout))

  write(10, '(A,I10)') 'nvert= ', para%nvert
  do nn = 1, para%nvert
     write(10, '(A)') '  vertloop'
     write(10, '(A,3E15.7)') '    vertex ', vert(nn)%xyz(1), vert(nn)%xyz(2), vert(nn)%xyz(3)
     write(10, '(A,3E15.7)') '    pseudonormal ', vert(nn)%norm(1), vert(nn)%norm(2), vert(nn)%norm(3)
     write(10, '(A,I3)') '      nshare= ', vert(nn)%share(1)
     do mm = 1, vert(nn)%share(1)
        write(10, '(A,I10)') '        ', vert(nn)%share(mm + 1)
     end do
     write(10, '(A)') '  endvertloop'
  end do
  write(10, '(A)') 'endvertex'

  write(10, '(A,I10)') 'nedge= ', para%nedgea
  do nn = 1, para%nedgea
     write(10, '(A)') '  edgeloop'
     write(10, '(A,3E15.7)') '    edgecenter ', edge(nn)%xyz(1), edge(nn)%xyz(2), edge(nn)%xyz(3)
     write(10, '(A,3E15.7)') '    edgenormal ', edge(nn)%norm(1), edge(nn)%norm(2), edge(nn)%norm(3)
     write(10, '(A,2I10)') '      vertindx ', edge(nn)%vert(1), edge(nn)%vert(2)
     write(10, '(A)') '  endedgeloop'
  end do
  write(10, '(A)') 'endedge'

  write(10, '(A,I10)') 'nelem= ', para%nelem
  do nn = 1, para%nelem
     write(10, '(A)') '  elemloop'
     write(10, '(A,3E15.7)') '    facetcenter ', elem(nn)%xyz(1), elem(nn)%xyz(2), elem(nn)%xyz(3)
     write(10, '(A,3E15.7)') '    facetnormal ', elem(nn)%norm(1), elem(nn)%norm(2), elem(nn)%norm(3)
     write(10, '(A,3I10)') '      vertindx ', elem(nn)%vert(1), elem(nn)%vert(2), elem(nn)%vert(3)
     write(10, '(A,3I10)') '      edgeindx ', elem(nn)%edge(1), elem(nn)%edge(2), elem(nn)%edge(3)
     write(10, '(A)') '  endelemloop'
  end do
  write(10, '(A)') 'endelem'

  close(10)

  return
end subroutine output_custom


!---------------------------------------------------------------------------------------
!
!     Output_TecPlot : Write IBSurface in TecPlot Format
!
!---------------------------------------------------------------------------------------
subroutine output_tecplot(paramx, vert, elem, para, nls, nseq)
  use GenSurf
  implicit none
  integer, intent(in) :: nls, nseq
  type(surfpar), intent(in) :: paramx, para(nls)
  type(vertex), intent(in) :: vert(paramx%nvert, nls)
  type(element), intent(in) :: elem(paramx%nelem, nls)

  ! local variable
  integer :: np(3)
  character(len=80) :: fout, foutl
  character(len=80) :: fvert, felem, fedge, fzone, fnls
  integer :: nn, mm, ll, mc, icount

  write(foutl, '(I10)') nseq
  fout = 'Output-Tec/Output_Tec' // trim(adjustl(foutl)) // '.dat'

  open(10, file=trim(fout))
  write(10, '(A)') ' TITLE= "Finite-Element Data"'
  write(10, '(A)') ' VRIABLES = "X" "Y" "Z"'

  do ll = 1, nls

     write(fnls, '(I10)') ll
     fzone = 'Object' // trim(adjustl(fnls)) // ''

     icount = 0
     do nn = 1, para(ll)%nelem
        do mm = 1, 3
           if (elem(nn, ll)%conn(mm) /= 0) icount = icount + 1
        end do
     end do

     write(fvert, '(I10)') para(ll)%nvert
     write(felem, '(I10)') para(ll)%nelem
     write(fedge, '(I10)') icount

     write(10, '(A,A,A)') ' ZONE T="', trim(adjustl(fzone)), '"'
     write(10, '(A,A,A)') ' STRANDID=', trim(adjustl(fnls)), ', SOLUTIONTIME=0'
     write(10, '(A,A,A,A,A)') ' N=', trim(adjustl(fvert)), ', E=', trim(adjustl(felem)), ', ZONETYPE=FETriangle'
     write(10, '(A)') ' DATAPACKING=POINT'
     write(10, '(A,A)') ' FACENEIGHBORCONNECTIONS=', trim(adjustl(fedge))
     write(10, '(A)') ' FACENEIGHBORMODE=LOCALONETOONE'
     write(10, '(A)') ' DT=(DOUBLE DOUBLE DOUBLE )'

     do nn = 1, para(ll)%nvert
        write(10, '(3(1PE22.9E3))') vert(nn, ll)%xyz(1), vert(nn, ll)%xyz(2), vert(nn, ll)%xyz(3)
     end do
     do nn = 1, para(ll)%nelem
        write(10, *) elem(nn, ll)%vert(1), elem(nn, ll)%vert(2), elem(nn, ll)%vert(3)
     end do

     do nn = 1, para(ll)%nelem
        mc = 0
        do mm = 1, 3
           if (elem(nn, ll)%conn(mm) /= 0) then
              mc = mc + 1
              write(10, *) nn, mc, elem(nn, ll)%conn(mm)
           end if
        end do
     end do

  end do

  close(10)

  return
end subroutine output_tecplot


!---------------------------------------------------------------------------------------
!
!     Output_bin : Write IBSurface in binary Format w or w/o elem_vert
!
!---------------------------------------------------------------------------------------
subroutine output_bin(fileout, paramx, vert, elem, para, ibseq, nls, nseq)
  use GenSurf
  implicit none
  integer, intent(in) :: nls, nseq
  type(surfpar), intent(in) :: paramx
  type(surfpar), intent(inout) :: para(nls)
  type(vertex), intent(in) :: vert(paramx%nvert, nls)
  type(element), intent(in) :: elem(paramx%nelem, nls)
  integer, intent(in) :: ibseq(nls)
  character(len=80), intent(in) :: fileout

  character(len=80) :: fout, foutl
  integer :: nn, ll

  para(:)%iedge = 2  ! default

  write(foutl, '(I10)') nseq
  fout = trim(fileout) // trim(adjustl(foutl)) // '.bin'

  open(10, file=trim(fout), form='unformatted')

  if (nseq == 0) then

     write(10) nls
     do ll = 1, nls
        write(10) para(ll)
     end do
     do ll = 1, nls
        write(10) (vert(nn, ll)%xyz, nn = 1, para(ll)%nvert)
     end do
     do ll = 1, nls
        write(10) (elem(nn, ll)%vert, nn = 1, para(ll)%nelem)
     end do

  else

     do ll = 1, nls
        if (ibseq(ll) == 1) then
           write(10) (vert(nn, ll)%xyz, nn = 1, para(ll)%nvert)
        end if
     end do

  end if

  close(10)

  return
end subroutine output_bin