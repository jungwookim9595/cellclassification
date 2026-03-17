C-----------------------------------------------------------------------------------
C
C     Program GenSurf.f 
C
C     This program is designed to generate the node surface points from TecPlot FE data
C     Input data: vertices, triangle information (index for triange elements' vertices
C                                             and index for connected triangle elements)
C     Output data: position, normal vectors, velocity (optional)
C                  at vertices, edges, and centers of the elements (triangles)
C
C     07/01/08
C
C     Jung-Il Choi, Ph.D.,
C     Research Assistant Professor
C     Dept. Mechanical & Aerospace Engineering
C     North Carolina State University
C
C-----------------------------------------------------------------------------------

C     Define Data Strctures MODULE

      Module GenSurf
      integer, parameter :: nshare = 500
      type vertex
           real,dimension(3) :: xyz,xyz0,norm,vel,vel0,acc
           real :: area,power
           integer,dimension(nshare) :: share
           integer :: open
      end type
      type element
           real,dimension(3) :: xyz,norm
           real :: area,power
           integer,dimension(3) :: vert,conn,edge
      end type
      type edges
           real,dimension(3) :: xyz,norm
           real :: area,power
           integer, dimension(2) :: vert,elem
           integer :: open
      end type
      type surfpar
           integer :: nvert,nelem,nedge,nedgea,iedge
      end type
      type vertexorg
           real,dimension(3) :: xyz
           integer,dimension(nshare) :: share
      end type
      type elementorg
           integer,dimension(3) :: vert,conn
      end type
      end Module

C---------------------------------------------------------------------------------------
C 
C     STR2NUM : String to number
C
C---------------------------------------------------------------------------------------

      function str2num(string)
      character*100 string
      integer str2num

      lens = len(trim(string))
      istr = lens
      iend = 0
      do i=1,lens
         if(ichar(string(i:i)).ge.48.and.ichar(string(i:i)).lt.58) then
            istr = min(i,istr)
            iend = max(i,iend)
         endif
      enddo

      str2num = 0
      do i=istr,iend
         inum = ichar(string(i:i))-48
         str2num  = str2num + inum*10**(iend-i)
      enddo
      end function



C---------------------------------------------------------------------------------------
C 
C     Normal_Correction : Correction for Normal Vectors considering Connectivity
C                         Correct the order of the element storage!
C
C---------------------------------------------------------------------------------------
      subroutine normal_correction(paramx,vert,elem,para,indx1,indx2,nc_corr)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertex) :: vert(paramx%nvert)
      type(element):: elem(paramx%nelem)

      integer indxv1(3),indxv2(3)
      
      indxv1(:) = elem(indx1)%vert(:)
      indxv2(:) = elem(indx2)%vert(:)

C---- find sharing edge
      do m1=1,3
         m2 = mod(m1,3)+1
         n1 = indxv1(m1)
         n2 = indxv1(m2)

         iflg = 0
         do j=1,3
            if(indxv2(j).eq.n1) iflg = iflg+1
            if(indxv2(j).eq.n2) iflg = iflg+1
            if(indxv2(j).ne.n1.and.indxv2(j).ne.n2) n3 = indxv2(j)
         enddo
         if(iflg.eq.2) then 
            mt1 = m1
            mt2 = m2
            nn3 = n3
            goto 10
         endif
      enddo

10    continue
      nt1 = 4 - mt1
      nt2 = 4 - mt2
      nt3 = 6 - nt1 - nt2

      if(indxv2(nt1).ne.indxv1(mt1).or.indxv2(nt2).ne.indxv1(mt2)) nc_corr = nc_corr + 1
      elem(indx2)%vert(nt1) = indxv1(mt1)
      elem(indx2)%vert(nt2) = indxv1(mt2)
      elem(indx2)%vert(nt3) = nn3

      return
      end

C---------------------------------------------------------------------------------------
C 
C     Conn_Correction : Correction for wrong connectivity
C
C---------------------------------------------------------------------------------------
      subroutine conn_correction(paramx,vert,elem,para)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertex) :: vert(paramx%nvert)
      type(element):: elem(paramx%nelem)

      integer iconn(para%nelem,3)
      
      iconn = 0
      do nn=1,para%nelem
      do m1=1,3
         m2 = mod(m1,3)+1
         m3 = 6 - m1 - m2
         mv1 = elem(nn)%vert(m1)
         mv2 = elem(nn)%vert(m2)

         do n=1,3
         nt = elem(nn)%conn(n)

         if(nt.ne.0) then     
            do n1=1,3
               n2 = mod(n1,3)+1
               n3 = 6 - n1 - n2
               nv1 = elem(nt)%vert(n1)
               nv2 = elem(nt)%vert(n2)
               if(mv1.eq.nv2.and.mv2.eq.nv1) then
                  iconn(nn,m3) = nt
                  iconn(nt,n3) = nn
               endif
            enddo
         endif
         enddo
      enddo
      enddo

      do nn=1,para%nelem
         elem(nn)%conn(:) = iconn(nn,:)
      enddo

      return
      end



C---------------------------------------------------------------------------------------
C 
C     Shared_Vertices : Find shared vertices
C
C---------------------------------------------------------------------------------------
      subroutine shared_vertices(paramx,vert,elem,para)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertex)  :: vert(paramx%nvert)
      type(element) :: elem(paramx%nelem)

c     local variable
      integer nc(para%nvert)

      nc = 0
      do nn=1,para%nelem
         do mm=1,3
           np = elem(nn)%vert(mm)
           nc(np) = nc(np) + 1
           vert(np)%share(nc(np)+1) = nn
         enddo
      enddo

      do nn=1,para%nvert
         vert(nn)%share(1) = nc(nn) ! Total number of sharing points
         if(nc(nn).eq.0) write(6,*) 'No sharing vertex',nn
      enddo

      return
      end

C---------------------------------------------------------------------------------------
C 
C     Face_Normal : Calculate outward-pointing normal vectors at the face of elements
C
C---------------------------------------------------------------------------------------
      subroutine face_normal(paramx,vert,elem,para)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertex)  :: vert(paramx%nvert)
      type(element) :: elem(paramx%nelem)

c     local variable
      integer np(3)
      real xp1(3),xp2(3),xp3(3)

      onethird = 1.0/3.0

      do nn=1,para%nelem
         np(:) = elem(nn)%vert(:)
         xp2(:) = vert(np(3))%xyz(:)-vert(np(2))%xyz(:)
         xp1(:) = vert(np(2))%xyz(:)-vert(np(1))%xyz(:)

         xp3(1) = xp1(2)*xp2(3) - xp1(3)*xp2(2)
         xp3(2) = xp1(3)*xp2(1) - xp1(1)*xp2(3)
         xp3(3) = xp1(1)*xp2(2) - xp1(2)*xp2(1)

         xpnorm = sqrt(xp3(1)**2+xp3(2)**2+xp3(3)**2)
         elem(nn)%xyz(:)  = onethird*(vert(np(1))%xyz(:)+vert(np(2))%xyz(:)+vert(np(3))%xyz(:))
         elem(nn)%norm(:) = xp3(:)/xpnorm
         elem(nn)%area = 0.5*xpnorm
      enddo

      return
      end

C---------------------------------------------------------------------------------------
C 
C     Pseudo_Normal : Calculate outward-pointing pseudo-normal vectors at vertices
C                     Angle-weighted pseudo-normal vectors
C
C---------------------------------------------------------------------------------------
      subroutine pseudo_normal(paramx,vert,elem,edge,para)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertex)  :: vert(paramx%nvert)
      type(element) :: elem(paramx%nelem)
      type(edges)   :: edge(paramx%nedge)

c     local variable
      integer np(3),nconn(3)
      real xp1(3),xp2(3),pnorm(3)

      do nn=1,para%nvert

         pnorm(:) = 0.0
         area = 0.0
         suma = 0.0

         do mm=1,vert(nn)%share(1)

        !angle computation
         indx = vert(nn)%share(mm+1)
         np(:) = elem(indx)%vert(:)

         if(nn.eq.np(1)) then
            xp1(:) = vert(np(2))%xyz(:) - vert(np(1))%xyz(:)
            xp2(:) = vert(np(3))%xyz(:) - vert(np(1))%xyz(:)
         elseif(nn.eq.np(2)) then
            xp1(:) = vert(np(1))%xyz(:) - vert(np(2))%xyz(:)
            xp2(:) = vert(np(3))%xyz(:) - vert(np(2))%xyz(:)
         elseif(nn.eq.np(3)) then
            xp1(:) = vert(np(1))%xyz(:) - vert(np(3))%xyz(:)
            xp2(:) = vert(np(2))%xyz(:) - vert(np(3))%xyz(:)
         endif

         x12dot= xp1(1)*xp2(1)+xp1(2)*xp2(2)+xp1(3)*xp2(3)
         vecmag1=sqrt(xp1(1)**2+xp1(2)**2+xp1(3)**2)
         vecmag2=sqrt(xp2(1)**2+xp2(2)**2+xp2(3)**2)
         angle=acos(x12dot/vecmag1/vecmag2)
 
         pnorm(:)=pnorm(:)+angle*elem(indx)%norm(:)
         suma = suma + angle
         area = area + angle*elem(indx)%area
         enddo

         pnorm(:) = pnorm(:)/suma
         pnormmag = sqrt(pnorm(1)**2+pnorm(2)**2+pnorm(3)**2)

         vert(nn)%norm(:) = pnorm(:)/pnormmag
         vert(nn)%area = area/suma

      enddo

c**** re-define pseudo-normal vectors for disconnected edges
      if(para%nedgea.ne.para%nelem*3/2) then ! This must be checked!!!

      do nn=1,para%nvert

c *** detect vertex involving disconnected edges
      ndetect = 1 ! indicator for disconnected edge, if ndetect=0
      do mm=1,vert(nn)%share(1)
         indx = vert(nn)%share(mm+1)
         np(:) = elem(indx)%vert(:)
         nconn(:) = elem(indx)%conn(:)
         ndetect= ndetect*nconn(1)*nconn(2)*nconn(3) 
      enddo

      if(ndetect.eq.0) then  ! do it again
         pnorm(:) = 0.0
         do mm=1,vert(nn)%share(1)
            indx = vert(nn)%share(mm+1)
            np(:) = elem(indx)%vert(:)
            nconn(:) = elem(indx)%conn(:)
            nconnp= nconn(1)*nconn(2)*nconn(3) 

            if(nconnp.eq.0) then
               do ip=1,3
                  if(nn.eq.np(ip)) ipp = ip
               enddo
               do iconn=1,3 !find edge information
                  if(nconn(iconn).eq.0.and.ipp.ne.iconn) then 
                     ne = elem(indx)%edge(iconn)
                     pnorm(:) = pnorm(:) + edge(ne)%norm(:)
                  endif
               enddo
            endif
         enddo

         xnormmag = sqrt(pnorm(1)**2+pnorm(2)**2+pnorm(3)**2)
         vert(nn)%norm(:) = pnorm(:)/xnormmag
      endif
      enddo ! end of vertex

      endif

      return
      end

C---------------------------------------------------------------------------------------
C 
C     Shared_Edge : Find the elements sharing edge
C
C---------------------------------------------------------------------------------------
      subroutine shared_edge(paramx,vert,elem,edge,para)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertex)  :: vert(paramx%nvert)
      type(element) :: elem(paramx%nelem)
      type(edges)   :: edge(paramx%nedge)

      integer,allocatable,dimension(:,:)::ivisit

      allocate(ivisit(para%nelem,3))

      ivisit = 0
      ne = 0
      do nn=1,para%nelem
      do m1=1,3
         m2 = mod(m1,3)+1
         m3 = 6-m1-m2
         n1 = elem(nn)%vert(m1)
         n2 = elem(nn)%vert(m2)
         nt = elem(nn)%conn(m3)

         if(nt.ne.0) then
         do mm=1,3
            if(elem(nt)%conn(mm).eq.nn) n3 = mm
         enddo

         if(ivisit(nn,m3).eq.0) then
            ne = ne + 1
            edge(ne)%vert(1) = n1
            edge(ne)%vert(2) = n2
            edge(ne)%elem(1) = nn
            edge(ne)%elem(2) = nt

            ivisit(nn,m3) = 1
            ivisit(nt,n3) = 1
            elem(nn)%edge(m3) = ne
            elem(nt)%edge(n3) = ne
         endif

         elseif(nt.eq.0)  then

         if(ivisit(nn,m3).eq.0) then
            ne = ne + 1
            edge(ne)%vert(1) = n1
            edge(ne)%vert(2) = n2
            edge(ne)%elem(1) = nn
            edge(ne)%elem(2) = nt

            ivisit(nn,m3) = 1
            elem(nn)%edge(m3) = ne
         endif

         endif

      enddo
      enddo

      para%nedgea = ne

      deallocate(ivisit)

C---- Check ------------------------------------
      iopt_check = 0
      if(iopt_check.eq.1) then
         write(6,*) 'Total edge number',ne
         write(6,*) 'edge'
         do nn=1,para%nedgea
c            write(6,999) nn,edge(nn)%vert(1:2)
            write(6,999) nn,edge(nn)%elem(1:2),edge(nn)%vert(1:2)
            write(6,998) nn,vert(edge(nn)%vert(1))%xyz,vert(edge(nn)%vert(2))%xyz
         enddo
c         write(6,*) 'elem'
c         do nn=1,para%nelem
c            write(6,999) nn,elem(nn)%edge
c         enddo
      endif
999   format(30i6)
998   format(i6,30e15.7)

      return
      end

C---------------------------------------------------------------------------------------
C 
C     Edge_Normal : Calculate outward-pointing pseudo-normal vectors at edges
C                   Treatment of thin body surface is included
C
C---------------------------------------------------------------------------------------
      subroutine edge_normal(paramx,vert,elem,edge,para)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertex)  :: vert(paramx%nvert)
      type(element) :: elem(paramx%nelem)
      type(edges)   :: edge(paramx%nedge)

c     local variable
      real xp1(3),xp2(3),xp3(3),xpm(3),xn1(3),xn2(3),xnp(3),xnm3(3)

      vert(:)%open = 0 ! connected vortex
      ioptedge = para%iedge


      do ne=1,para%nedgea

         edge(ne)%open = 0 ! connected edge

         n1 = edge(ne)%vert(1)
         n2 = edge(ne)%vert(2)
         nn = edge(ne)%elem(1)
         nt = edge(ne)%elem(2)

         if(nt.ne.0) then
         xp1(:) = vert(n1)%xyz(:)
         xp2(:) = vert(n2)%xyz(:)
         xpm(:) = 0.5*(xp1(:)+xp2(:))

         xn1(:) = elem(nn)%norm(:)
         xn2(:) = elem(nt)%norm(:)
         xnp(:) = 0.5*(xn1(:)+xn2(:))
         xnpmag = sqrt(xnp(1)**2+xnp(2)**2+xnp(3)**2)

         xnp(:) = xnp(:)/xnpmag
         area = 0.5*(elem(nn)%area + elem(nt)%area)

         edge(ne)%xyz(:)  = xpm(:)
         edge(ne)%norm(:) = xnp(:)
         edge(ne)%area    = area

         elseif(nt.eq.0) then

         xp1(:) = vert(n1)%xyz(:)
         xp2(:) = vert(n2)%xyz(:)
         xpm(:) = 0.5*(xp1(:)+xp2(:))
 
         xn1(:) = elem(nn)%norm(:) ! plane normal vector
         xn2(:) = xp1(:)-xp2(:)    ! edge vector

         call calc_normal_thin_edge(xn1,xn2,xnp)

         do mm=1,3
            m3 = elem(nn)%vert(mm)
            if(m3.ne.n1.and.m3.ne.n2) n0 = m3 ! not shared points at the edge but in elements
         enddo

         xp3(:) = vert(n0)%xyz(:)

         xnm3(:) = xpm(:)-xp3(:)
         vecdot = xnm3(1)*xnp(1)+xnm3(2)*xnp(2)+xnm3(3)*xnp(3)
         xnsign = sign(1.0,vecdot)

         xnp(:) = xnsign*xnp(:)

C ***  edge normal vector back to face normal vector
c        !modified edge_normal vector considering face_normal vector
        !if(ioptedge.eq.0) nothing - use tangential vecotr
         if(ioptedge.eq.2) xnp(:) = elem(nn)%norm(:)+xnp(:) ! use plane normal + tangential vector
         if(ioptedge.eq.1) xnp(:) = elem(nn)%norm(:) ! use plane normal vector

         xnpmag = sqrt(xnp(1)**2+xnp(2)**2+xnp(3)**2)

         xnp(:) = xnp(:)/xnpmag

         area = elem(nn)%area
         edge(ne)%xyz(:)  = xpm(:)
         edge(ne)%norm(:) = xnp(:)
         edge(ne)%area    = area
         edge(ne)%open    = 1 ! disconnected edge

         vert(n1)%open = 1 ! disconnected vertex
         vert(n2)%open = 1 ! disconnected vertex

         endif

      enddo

      return
      end

C---------------------------------------------------------------------------------------
C 
C     Calc_Normal_Thin_Edge : Calculate outward-pointing pseudo-normal vectors at edges
C                             for thin body
C                             a*c =0, b*c=0, |c|=1, find c vector
C
C---------------------------------------------------------------------------------------
      subroutine calc_normal_thin_edge(xn1,xn2,xn3)

      real xn1(3),xn2(3),xn3(3)
c     local variable
      real xna(3)

      xna(:) = abs(xn1(:))
      xnamx = max(max(xna(1),xna(2)),xna(3))

      do i=1,3
         if(xnamx.eq.xna(i)) then
            indx1 = i
            if(indx1.eq.1) indx2 = 2
            if(indx1.eq.2) indx2 = 3
            if(indx1.eq.3) indx2 = 1
         endif
      enddo
      indx3 = 6 - indx1 - indx2

      c2 = -xn1(indx2)*xn2(indx1)/xn1(indx1)+xn2(indx2)
      c3 = -xn1(indx3)*xn2(indx1)/xn1(indx1)+xn2(indx3)

      c2a = abs(c2)
      c3a = abs(c3)
      camax = max(c2a,c3a)
      if(camax.eq.c2a) then
         a2 = -c3/c2
         a1 = -xn1(indx2)/xn1(indx1)*a2-xn1(indx3)/xn1(indx1)
         xn3(indx3) = sqrt(1.0/(a1**2+a2**2+1.0))
         xn3(indx2) = a2*xn3(indx3)
         xn3(indx1) = a1*xn3(indx3)
      else
         a3 = -c2/c3
         a1 = -xn1(indx2)/xn1(indx1)-xn1(indx3)/xn1(indx1)*a3
         xn3(indx2) = sqrt(1.0/(a1**2+a3**2+1.0))
         xn3(indx3) = a3*xn3(indx2)
         xn3(indx1) = a1*xn3(indx2)
      endif

      return
      end


C---------------------------------------------------------------------------------------
C 
C     Convert2IBSurface : Read ASCII STL Files
C                         Find Duplicated Vertices & Elements 
C                         Find Shared Vertices and Edges
C                         Find Connectivity
C                         Find Segments
C                         Write IBSurface in Tecplot, Binary, new STL, Nodesurf formats
C
C---------------------------------------------------------------------------------------
      subroutine convert2ibsurface(filestl,paramx,vert,elem,para,indxmapping,ielemmapping,ielemseg,nlsf)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertexorg) :: vert(paramx%nvert)
      type(elementorg):: elem(paramx%nelem)

      integer indxmapping(paramx%nvert,2)
      integer ielemmapping(paramx%nelem)
      integer ielemseg(paramx%nelem)

      real,allocatable,dimension(:,:) ::  datapts

      character*80 filestl

      nelemorg = para%nelem
      nvertorg = para%nvert

      allocate(datapts(para%nvert,3))

C---- Read STL data 
      call read_stl(filestl,datapts,nelemorg,nvertorg)
      write(6,*) 'Original parameters, #vert, #elem', nvertorg, nelemorg

C---- Find duplication of vertices and index mapping
      call find_duplications(datapts,indxmapping,ielemmapping,paramx,para,nudvert,nudelem)

C---- Make data structure of triangle elements
      para%nvert = nudvert ! re-define para%nvert using un-duplicate vertices
      para%nelem = nudelem ! re-define para%nvert using un-duplicate elements
      para%nedge = nudelem*3/2
      write(6,*) 'Undupliate #vert, #elem',para%nvert, para%nelem

      do nn=1,para%nvert
         indx = indxmapping(nn,1)
         vert(nn)%xyz(:) = datapts(indx,:)
      enddo

      do nn=1,nelemorg
         ielemmap = ielemmapping(nn)
         if(ielemmap.ne.0) then
         mm = 3*(nn-1)
         elem(ielemmap)%vert(1) = indxmapping(mm+1,2)
         elem(ielemmap)%vert(2) = indxmapping(mm+2,2)
         elem(ielemmap)%vert(3) = indxmapping(mm+3,2)
         endif
      enddo

C---- Find share vertices
      call shared_vertices_org(paramx,vert,elem,para)
      maxshare = 1
      do nn=1,para%nvert
         maxshare = max(maxshare,vert(nn)%share(1))
      enddo
      write(6,*) 'Maximum shared vertices', maxshare
999   format(100i6)

C---- Find connectivity
      call find_connectivity_org(paramx,vert,elem,para)

C---- Find segment
      call find_segment(paramx,vert,elem,ielemseg,para)
      do nn=1,para%nelem
         nlsf = max(1,ielemseg(nn))
         if(ielemseg(nn).eq.0) then
         write(6,*) 'Un-assigned elements are detected!!',nn
         end if
      enddo
      write(6,*) 'Number of levelsets', nlsf

      deallocate(datapts)

      return
      end


C---------------------------------------------------------------------------------------
C 
C     Define_Element : Read ASCII STL Files and define # of vertices and elements
C
C---------------------------------------------------------------------------------------
      subroutine define_element(filestl,nelem)
      character*80 filestl,cname

C     Check number of triangle elements
      line=0
      open(11,file=trim(filestl))
1     read(11,*) cname
      line=line+1
      if(trim(cname).ne.'endsolid') goto 1
      close(11)

C     Read STL data file and define domain size
      nelem=(line-2)/7

      return
      end

C---------------------------------------------------------------------------------------
C 
C     Read_STL : Read ASCII STL Files and store positions of vertices in datapts
C
C---------------------------------------------------------------------------------------
      subroutine read_stl(filestl,datapts,nelemorg,nvertorg)
      real datapts(nvertorg,3)
      character*80 filestl,cname

      open(11,file=trim(filestl))
      read(11,*)
      nn = 0
      do i=1,nelemorg
      read(11,*) cname,cname,dummy1,dummy2,dummy3
      read(11,*) cname
      nn = nn + 1
      read(11,*) cname,datapts(nn,1),datapts(nn,2),datapts(nn,3)
      nn = nn + 1
      read(11,*) cname,datapts(nn,1),datapts(nn,2),datapts(nn,3)
      nn = nn + 1
      read(11,*) cname,datapts(nn,1),datapts(nn,2),datapts(nn,3)
      read(11,*) cname
      read(11,*) cname
      enddo
      close(11)

      return
      end

C---------------------------------------------------------------------------------------
C 
C     Find_Duplications : Find duplicated vertices and elements
C
C---------------------------------------------------------------------------------------
      subroutine find_duplications(datapts,indxmapping,ielemmapping,paramx,para,nudvert,nudelem)
      USE GenSurf
      type(surfpar):: paramx,para

      integer indxmapping(paramx%nvert,2)
      integer ielemmapping(paramx%nelem)

      real datapts(para%nvert,3)

c     local variable
      integer, allocatable, dimension(:,:) :: elem2vert,vert2share,indxshare
      integer, allocatable, dimension(:) :: ncshare

      integer, allocatable, dimension(:) :: indxarray,indxarray2
      real, allocatable, dimension(:) :: distarray,distarray2

      integer, allocatable, dimension(:) :: indxvisit,idupelem
      integer, allocatable, dimension(:,:) :: idup

      integer np(3),mp(3)
      real qp(3)

      allocate(elem2vert(para%nelem,3))
      allocate(vert2share(para%nvert,nshare))
      allocate(ncshare(para%nvert))
      allocate(indxshare(para%nelem,3))
      allocate(indxvisit(para%nvert))
      allocate(idupelem(para%nelem))
      allocate(idup(para%nelem,3))

C *** set up search space
      nsmx = min(50,para%nvert)
      nsmx2= para%nvert ! entire search

      allocate(indxarray(nsmx))
      allocate(distarray(nsmx))
      allocate(indxarray2(nsmx2))
      allocate(distarray2(nsmx2))

      call wldst_create_tree(datapts(1,1),datapts(1,2),datapts(1,3),para%nvert)

      indxvisit = 0
      ncount = 0

      do nn=1,para%nvert
         if(indxvisit(nn).eq.0) then
            qp(:) = datapts(nn,:)
            ncount = ncount + 1
            indxvisit(nn) = ncount
            call wldst_getdists_xyz(qp(1),qp(2),qp(3),indxarray,distarray,nsmx)
            if(distarray(nsmx).ne.0.0) then
               do mm=1,nsmx
                  if(distarray(mm).eq.0.0) indxvisit(indxarray(mm)+1) = ncount
               enddo
            else  ! extended search
               if(nsmx.ne.nsmx2) then
               call wldst_getdists_xyz(qp(1),qp(2),qp(3),indxarray2,distarray2,nsmx2)
               do mm=1,nsmx2
                  if(distarray2(mm).eq.0.0) indxvisit(indxarray2(mm)+1) = ncount
               enddo
               endif
            endif
         endif
      enddo

      call wldst_delete_tree

      deallocate(indxarray,distarray,indxarray2,distarray2)

C---- Find data mapping between un-duplicate vertices and original data
      nudvert = ncount

c     indxmapping(:,1) : indx from new set to original set
c     indxmapping(:,2) : indx from original set to new set

      indxmapping = para%nvert
      do nn=1,para%nvert
         indxtag = indxvisit(nn)
         indxmapping(indxtag,1) = min(indxmapping(indxtag,1),nn)
         indxmapping(nn,2) = indxtag
      enddo

c *** Check duplicate elements

C --- Define shared vertices
      ncshare = 0
      do nn=1,para%nelem
         mm = 3*(nn-1)
         elem2vert(nn,1) = indxmapping(mm+1,2)
         elem2vert(nn,2) = indxmapping(mm+2,2)
         elem2vert(nn,3) = indxmapping(mm+3,2)
         do ll=1,3
            npl = elem2vert(nn,ll)
            ncshare(npl) = ncshare(npl) + 1
            vert2share(npl,ncshare(npl)+1) = nn
         enddo
      enddo

      do nn=1,nudvert
         vert2share(nn,1) = ncshare(nn)
      enddo

      idupelem = 0
      do nn=1,para%nelem
         np(:) = elem2vert(nn,:)
         mp(:) = vert2share(np(:),1)
         do k=1,3
         do mm=1,mp(k)
            indxshare(mm,k) = vert2share(np(k),mm+1)
         enddo
         enddo
         do m1=1,mp(1)
         do m2=1,mp(2)
            indxdup = indxshare(m1,1)
            if(indxdup.ne.nn.and.indxshare(m2,2).eq.indxdup) then
               do m3=1,mp(3)
                  if(indxshare(m3,3).eq.indxdup) idupelem(indxdup) = 1
               enddo
            endif
         enddo
         enddo
      enddo

      ncount = 0
      ielemmapping = 0
      do nn=1,para%nelem
         if(idupelem(nn).eq.0) then
         ncount = ncount + 1
         ielemmapping(nn) = ncount
         endif
      enddo
      nudelem = ncount

      deallocate(elem2vert,vert2share,ncshare,indxshare,indxvisit,idupelem,idup)

      return
      end


C---------------------------------------------------------------------------------------
C
C     Shared_Vertices : Find shared vertices
C
C---------------------------------------------------------------------------------------
      subroutine shared_vertices_org(paramx,vert,elem,para)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertexorg)  :: vert(paramx%nvert)
      type(elementorg) :: elem(paramx%nelem)

c     local variable
      integer,allocatable,dimension(:) :: nc

      allocate(nc(para%nvert))

      nc = 0
      do nn=1,para%nelem
         do mm=1,3
           np = elem(nn)%vert(mm)
           nc(np) = nc(np) + 1
           vert(np)%share(nc(np)+1) = nn
         enddo
      enddo

      do nn=1,para%nvert
         vert(nn)%share(1) = nc(nn) ! Total number of sharing points
         if(nc(nn).eq.0) write(6,*) 'No sharing vertex',nn
      enddo

      deallocate(nc)

      return
      end

C---------------------------------------------------------------------------------------
C
C     Find_Connectivity : Find connectivities of each elements
C
C---------------------------------------------------------------------------------------
      subroutine find_connectivity_org(paramx,vert,elem,para)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertexorg) :: vert(paramx%nvert)
      type(elementorg):: elem(paramx%nelem)

      integer indx(nshare,3)
      integer np(3),mp(3)

      do nn=1,para%nelem
         np(:) = elem(nn)%vert(:)
         mp(:) = vert(np(:))%share(1)
         do k=1,3
         do mm=1,mp(k)
            indx(mm,k) = vert(np(k))%share(mm+1)
         enddo
         enddo

         elem(nn)%conn(:) = 0
         nc = 0
         do ip1=1,3
            ip2 = mod(ip1,3)+1
            do mm=1,mp(ip1)
               if(indx(mm,ip1).ne.nn) then
               do ll=1,mp(ip2)
                  if(indx(mm,ip1).eq.indx(ll,ip2)) then
                     nc = nc + 1
                     elem(nn)%conn(nc) = indx(mm,ip1)
                  endif
               enddo
               endif
             enddo
         enddo
         if(nc.eq.0) write(6,*) 'element w/o connection is found',nn
         if(nc.gt.3) write(6,*) 'wrong connection is found',nn
      enddo


      return
      end



C---------------------------------------------------------------------------------------
C 
C     Find_Connectivity : Find connectivities of each elements
C
C---------------------------------------------------------------------------------------
      subroutine find_connectivity(paramx,vert,elem,para)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertex) :: vert(paramx%nvert)
      type(element):: elem(paramx%nelem)

      integer indx(1000,3)
      integer np(3),mp(3)

      do nn=1,para%nelem
         np(:) = elem(nn)%vert(:)
         mp(:) = vert(np(:))%share(1)
         do k=1,3
         do mm=1,mp(k)
            indx(mm,k) = vert(np(k))%share(mm+1)
         enddo
         enddo

         elem(nn)%conn(:) = 0
         nc = 0
         do ip1=1,3
            ip2 = mod(ip1,3)+1
            do mm=1,mp(ip1)
               if(indx(mm,ip1).ne.nn) then
               do ll=1,mp(ip2)
                  if(indx(mm,ip1).eq.indx(ll,ip2)) then
                     nc = nc + 1
                     elem(nn)%conn(nc) = indx(mm,ip1)
                  endif
               enddo
               endif
             enddo
         enddo
         if(nc.eq.0) write(6,*) 'element w/o connection is found',nn 
         if(nc.gt.3) write(6,*) 'wrong connection is found',nn 
      enddo
 
      return
      end

C---------------------------------------------------------------------------------------
C
C     Find_Segment : Find internal segments and correct order of the vertices in the elements
C
C---------------------------------------------------------------------------------------
      subroutine find_segment(paramx,vert,elem,ielemseg,para)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertexorg) :: vert(paramx%nvert)
      type(elementorg):: elem(paramx%nelem)

      integer ielemseg(paramx%nelem)

c     local variables
      integer, allocatable, dimension(:) :: ivisit,indx_list1,indx_list2

      allocate(ivisit(para%nelem))
      allocate(indx_list1(para%nelem))
      allocate(indx_list2(para%nelem))

      inseg = 1
      ivisit = 0
      nc_corr = 0

      nloopmx = 1000000
      nlist1 = 1
      indx_list1(1) = 1
      ivisit(1) = 1

30    continue
      nloop = 0
10    nlist2 = 0
      do n=1,nlist1
         indx1 = indx_list1(n)
         ivisit(indx1) = 1
         ielemseg(indx1) = inseg
         do m=1,3
            indx2 = elem(indx1)%conn(m)
            if(indx2.ne.0) then
               if(ivisit(indx2).eq.0) then
                  ivisit(indx2) = 1
                  nlist2 = nlist2 + 1
                  indx_list2(nlist2) = indx2
                  call normal_correction_org(paramx,vert,elem,para,indx1,indx2,nc_corr)
               endif
            endif
         enddo
      enddo

      if(nlist2.eq.0) goto 20
      indx_list1(1:nlist2) = indx_list2(1:nlist2)
      nlist1 = nlist2
      nloop = nloop + 1

      if(nloop.gt.nloopmx) goto 25
      goto 10

25    write(6,*) 'failed to find segment in Find_Segment'

20    continue

      iflg = 0
      do n=1,para%nelem
         if(ivisit(n).eq.0) then
            iflg = 1
            nlist1 = 1
            indx_list1(1) = n
            ivisit(n) = 1
            inseg = inseg + 1
            goto 30
         endif
      enddo

      if(iflg.eq.1) write(6,*) 'Wrong in check connection'
      if(nc_corr.ne.0) write(6,*) 'Number of normal-corrected elements',nc_corr

      deallocate(ivisit,indx_list1,indx_list2)

      return
      end


C---------------------------------------------------------------------------------------
C
C     Normal_Correction : Correction for Normal Vectors considering Connectivity
C                         Correct the order of the element storage!
C
C---------------------------------------------------------------------------------------
      subroutine normal_correction_org(paramx,vert,elem,para,indx1,indx2,nc_corr)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertexorg) :: vert(paramx%nvert)
      type(elementorg):: elem(paramx%nelem)

      integer indxv1(3),indxv2(3)

      indxv1(:) = elem(indx1)%vert(:)
      indxv2(:) = elem(indx2)%vert(:)

C---- find sharing edge
      do m1=1,3
         m2 = mod(m1,3)+1
         n1 = indxv1(m1)
         n2 = indxv1(m2)

         iflg = 0
         do j=1,3
            if(indxv2(j).eq.n1) iflg = iflg+1
            if(indxv2(j).eq.n2) iflg = iflg+1
            if(indxv2(j).ne.n1.and.indxv2(j).ne.n2) n3 = indxv2(j)
         enddo
         if(iflg.eq.2) then
            mt1 = m1
            mt2 = m2
            nn3 = n3
            goto 10
         endif
      enddo

10    continue
      nt1 = 4 - mt1
      nt2 = 4 - mt2
      nt3 = 6 - nt1 - nt2

      if(indxv2(nt1).ne.indxv1(mt1).or.indxv2(nt2).ne.indxv1(mt2)) nc_corr = nc_corr + 1
      elem(indx2)%vert(nt1) = indxv1(mt1)
      elem(indx2)%vert(nt2) = indxv1(mt2)
      elem(indx2)%vert(nt3) = nn3

      return
      end

C---------------------------------------------------------------------------------------
C
C     IBSurface_Seg_Param : Define segmented ibsurface parameters 
C
C---------------------------------------------------------------------------------------
      subroutine ibsurface_seg_param(paramx,vertorg,elemorg,paraorg,
     c                               indxmapping,ielemmapping,ielemseg,ncelem,ncvert,
     c                               ielem_seg2org,ivert_seg2org,ivert_org2seg,nls,nlsmap,nlsf,nlsfile)
      USE GenSurf
      type(surfpar):: paramx,paraorg(nlsfile)
      type(vertexorg) :: vertorg(paramx%nvert,nlsfile)
      type(elementorg):: elemorg(paramx%nelem,nlsfile)

      integer indxmapping(paramx%nvert,2,nlsfile)
      integer ielemmapping(paramx%nelem,nlsfile)
      integer ielemseg(paramx%nelem,nlsfile)
      integer nlsmap(nls),nlsf(nlsfile)

      integer ielem_seg2org(paramx%nelem,nls)
      integer ivert_seg2org(paramx%nvert,nls)
      integer ivert_org2seg(paramx%nvert,nlsfile)

      integer ncelem(nls),ncvert(nls)

      integer,allocatable,dimension(:,:) :: ivisit

      allocate(ivisit(paramx%nvert,nls))

c      if(nls.eq.1) then
c
c         ncelem(nls) = paraorg(1)%nelem
c         ncvert(nls) = paraorg(1)%nvert
c
c      else

C ***    re-define elements and vertices in each segment
         nlsindex = 0
         ncelem = 0
         ncvert = 0
         ivisit = 0

         do nf=1,nlsfile

         if(nf.ne.1) nlsindex = nlsindex + nlsf(nf-1)

         do nn=1,paraorg(nf)%nelem
            iseg = ielemseg(nn,nf) + nlsindex
            ncelem(iseg) = ncelem(iseg) + 1
            ielem_seg2org(ncelem(iseg),iseg) = nn
            do mm=1,3
               nv = elemorg(nn,nf)%vert(mm)
               if(ivisit(nv,iseg).eq.0) then
                  ivisit(nv,iseg) = 1
                  ncvert(iseg) = ncvert(iseg) + 1
                  ivert_seg2org(ncvert(iseg),iseg) = nv
                  ivert_org2seg(nv,nf) = ncvert(iseg)
               endif
             enddo
         enddo

         enddo

c      endif

      deallocate(ivisit)

      return
      end

C---------------------------------------------------------------------------------------
C
C     IBSurface_Seg : Define segmented ibsurface parameters (vert,elem,edge)
C
C---------------------------------------------------------------------------------------
      subroutine ibsurface_seg(paraorgmx,vertorg,elemorg,paraorg,
     c                         paramx,vert,elem,edge,para,
     c                         ielem_seg2org,ivert_seg2org,ivert_org2seg,nls,nlsmap,nlsf,nlsfile)
      USE GenSurf
      type(surfpar):: paraorgmx,paraorg(nlsfile)
      type(vertexorg) :: vertorg(paraorgmx%nvert,nlsfile)
      type(elementorg):: elemorg(paraorgmx%nelem,nlsfile)

      integer ielem_seg2org(paraorgmx%nelem,nls)
      integer ivert_seg2org(paraorgmx%nvert,nls)
      integer ivert_org2seg(paraorgmx%nvert,nlsfile)

c     Segmented objects
      type(surfpar):: paramx,para(nls)
      type(vertex) :: vert(paramx%nvert,nls)
      type(element):: elem(paramx%nelem,nls)
      type(edges)  :: edge(paramx%nedge,nls)

      integer nlsmap(nls),nlsf(nlsfile)

      integer np(3),mp(3)

      do ll=1,nls
         nf = nlsmap(ll)
         do nn=1,para(ll)%nvert
            mm = ivert_seg2org(nn,ll)
            vert(nn,ll)%xyz(:) = vertorg(mm,nf)%xyz(:)
         enddo
         do nn=1,para(ll)%nelem
            mm = ielem_seg2org(nn,ll)
            np(:) = elemorg(mm,nf)%vert(:)
            mp(:) = ivert_org2seg(np(:),nf)
            elem(nn,ll)%vert(:) = mp(:)
         enddo
      enddo

C *** All calculations for the elements (shared vertices/edges, face-normal,edge-normal,pseudo-normal)
C *** Check the normal vectors

      do ll=1,nls

         call shared_vertices(paramx,vert(1,ll),elem(1,ll),para(ll))
         call find_connectivity(paramx,vert(1,ll),elem(1,ll),para(ll))
         call conn_correction(paramx,vert(1,ll),elem(1,ll),para(ll))
         call shared_edge(paramx,vert(1,ll),elem(1,ll),edge(1,ll),para(ll))
         call face_normal(paramx,vert(1,ll),elem(1,ll),para(ll))
         call edge_normal(paramx,vert(1,ll),elem(1,ll),edge(1,ll),para(ll))
         call pseudo_normal(paramx,vert(1,ll),elem(1,ll),edge(1,ll),para(ll))

      enddo


      return
      end

C---------------------------------------------------------------------------------------
C
C     Ouput_STL : Write IBSurface in STL Format
C
C---------------------------------------------------------------------------------------
      subroutine output_stl(paramx,vert,elem,para,l)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertex) :: vert(paramx%nvert)
      type(element):: elem(paramx%nelem)

c     local variable
      integer np(3)
      character*80 fout,foutl

      write(foutl,990) l
990   format(i10)

      fout='output/Output_STL'//trim(adjustl(foutl))//'.stl'

      open(10,file=trim(fout))
      write(10,993)

      do nn=1,para%nelem
         np(:) = elem(nn)%vert(:)
         write(10,996) elem(nn)%norm(:)
         write(10,997)
         write(10,995) vert(np(1))%xyz(:)
         write(10,995) vert(np(2))%xyz(:)
         write(10,995) vert(np(3))%xyz(:)
         write(10,998)
         write(10,999)
      enddo
      write(10,994)
      close(10)

999   format('  endfacet')
998   format('    endloop')
997   format('    outer loop')
996   format('  facet normal ',3e15.7)
995   format('      vertex ',3e15.7)
994   format('endsolid object')
993   format('solid object')

      return
      end

C---------------------------------------------------------------------------------------
C
C     Ouput_Node : Write IBSurface in Scatter Format
C
C---------------------------------------------------------------------------------------
      subroutine output_node(paramx,vert,elem,edge,para,l)
      USE GenSurf

      type(surfpar):: paramx,para
      type(vertex)  :: vert(paramx%nvert)
      type(element) :: elem(paramx%nelem)
      type(edges)   :: edge(paramx%nedge)

c     local variables
      integer np(3)
      real xg(3),xn(3)
      character*80 fout,foutl

      onethird = 1.0/3.0

      write(foutl,990) l
990   format(i10)

      fout='output/Output_Node'//trim(adjustl(foutl))//'.dat'

      open(20,file=trim(fout))

      write(20,*) 'zone t="T',para%nelem,'"'
      do nn=1,para%nelem
         np(:) = elem(nn)%vert(:)
         xg(:) = onethird*(vert(np(1))%xyz(:)+vert(np(2))%xyz(:)+vert(np(3))%xyz(:))
         xn(:) = elem(nn)%norm(:)
         write(20,999) xg,xn
      enddo
      write(20,*) 'zone t="V',para%nvert,'"'
      do nn=1,para%nvert
         xg(:) = vert(nn)%xyz(:)
         xn(:) = vert(nn)%norm(:)
         write(20,999) xg,xn
      enddo
      write(20,*) 'zone t="E',para%nedgea,'"'
      do nn=1,para%nedgea
         xg(:) = edge(nn)%xyz(:)
         xn(:) = edge(nn)%norm(:)
         write(20,999) xg,xn
      enddo
      close(20)

999   format(30e15.7)

      return
      end

C---------------------------------------------------------------------------------------
C
C     Ouput_Custom : Write IBSurface in Custom format
C
C---------------------------------------------------------------------------------------
      subroutine output_custom(paramx,vert,elem,edge,para,l)
      USE GenSurf

      type(surfpar) :: paramx,para
      type(vertex)  :: vert(paramx%nvert)
      type(element) :: elem(paramx%nelem)
      type(edges)   :: edge(paramx%nedge)

c     local variables
      integer np(3)
      real xg(3),xn(3)
      character*80 fout,foutl

      onethird = 1.0/3.0

      write(foutl,990) l
990   format(i10)

      fout='output/Output'//trim(adjustl(foutl))//'.dat'

      open(10,file=trim(fout))

      write(10,770) para%nvert
      do nn=1,para%nvert
         write(10,771)
         write(10,772) vert(nn)%xyz(:)
         WRITE(10,773) vert(nn)%norm(:)
         write(10,774) vert(nn)%share(1)
         do mm = 1,vert(nn)%share(1)
            write(10,775) vert(nn)%share(mm+1)
         enddo
         write(10,776)
      enddo
      write(10,777)

      write(10,880) para%nedgea
      do nn=1,para%nedgea
         write(10,881)
         write(10,882) edge(nn)%xyz(:)
         write(10,883) edge(nn)%norm(:)
         write(10,884) edge(nn)%vert(:)
         write(10,885)
      enddo
      write(10,886)

      write(10,991) para%nelem
      do nn=1,para%nelem
         write(10,992)
         write(10,993) elem(nn)%xyz(:)
         write(10,994) elem(nn)%norm(:)
         write(10,995) elem(nn)%vert(:)
         write(10,996) elem(nn)%edge(:)
         write(10,997)
      enddo
      write(10,998)

      close(10)

770   format('nvert= ',I10)
771   format('  vertloop')
772   format('    vertex ',3e15.7)
773   format('    pseudonormal ',3e15.7)
774   format('      nshare= ',I3)
775   format('        ',I10)
776   format('  endvertloop')
777   format('endvertex')

880   format('nedge= ',I10)
881   format('  edgeloop')
882   format('    edgecenter ',3e15.7)
883   format('    edgenormal ',3e15.7)
884   format('      vertindx ',2I10)
885   format('  endedgeloop')
886   format('endedge')

991   format('nelem= ',I10)
992   format('  elemloop')
993   format('    facetcenter ', 3e15.7)
994   format('    facetnormal ', 3e15.7)
995   format('      vertindx ', 3I10)
996   format('      edgeindx ', 3I10)
997   format('  endelemloop')
998   format('endelem')


      return
      end

C---------------------------------------------------------------------------------------
C
C     Ouput_TecPlot : Write IBSurface in TecPlot Format
C
C---------------------------------------------------------------------------------------
      subroutine output_tecplot(paramx,vert,elem,para,nls,nseq)
      USE GenSurf
      type(surfpar):: paramx,para(nls)
      type(vertex) :: vert(paramx%nvert,nls)
      type(element):: elem(paramx%nelem,nls)

c     local variable
      integer np(3)
      character*80 fout,foutl
      character*80 fvert,felem,fedge,fzone,fnls

      write(foutl,990) nseq
990   format(i10)
      fout='Output-Tec/Output_Tec'//trim(adjustl(foutl))//'.dat'

      open(10,file=trim(fout))
      write(10,991)
      write(10,992)

      do ll=1,nls

      write(fnls,990) ll
      fzone = 'Object'//trim(adjustl(fnls))//''

      icount = 0
      do nn=1,para(ll)%nelem
         do mm=1,3
            if(elem(nn,ll)%conn(mm).ne.0) icount = icount + 1
         enddo
      enddo

      write(fvert,990) para(ll)%nvert
      write(felem,990) para(ll)%nelem
      write(fedge,990) icount

      write(10,993) trim(adjustl(fzone))
      write(10,994) trim(adjustl(fnls))
      write(10,995) trim(adjustl(fvert)),trim(adjustl(felem))
      write(10,996)
      write(10,997) trim(adjustl(fedge))
      write(10,998)
      write(10,999)

      do nn=1,para(ll)%nvert
         write(10,980) vert(nn,ll)%xyz
      enddo
      do nn=1,para(ll)%nelem
         write(10,*) elem(nn,ll)%vert
      enddo

      do nn=1,para(ll)%nelem
         mc = 0
         do mm=1,3
            if(elem(nn,ll)%conn(mm).ne.0) then
               mc = mc + 1
               write(10,*) nn,mc,elem(nn,ll)%conn(mm)
            endif
         enddo
      enddo

      enddo

      close(10)

991   format(' TITLE= "Finite-Element Data"')
992   format(' VRIABLES = "X" "Y" "Z"')
993   format(' ZONE T="',A,'"')
994   format(' STRANDID=',A,', SOLUTIONTIME=0')
995   format(' N=',A,', E=',A,', ZONETYPE=FETriangle')
996   format(' DATAPACKING=POINT')
997   format(' FACENEIGHBORCONNECTIONS=',A)
998   format(' FACENEIGHBORMODE=LOCALONETOONE')
999   format(' DT=(DOUBLE DOUBLE DOUBLE )')
980   format(3(1pe22.9e3))


      return
      end

C---------------------------------------------------------------------------------------
C
C     Ouput_bin : Write IBSurface in binary Format w or w/o elem_vert
C
C---------------------------------------------------------------------------------------
      subroutine output_bin(fileout,paramx,vert,elem,para,ibseq,nls,nseq)
      USE GenSurf
      type(surfpar):: paramx,para(nls)
      type(vertex) :: vert(paramx%nvert,nls)
      type(element):: elem(paramx%nelem,nls)
      integer ibseq(nls)

      character*80 fileout,fout,foutl

      para(:)%iedge = 2 ! default

      write(foutl,990) nseq
990   format(i10)

      fout=trim(fileout)//trim(adjustl(foutl))//'.bin'

      open(10,file=trim(fout),form='unformatted')

      if(nseq.eq.0) then

         write(10) nls
         do ll=1,nls
            write(10) para(ll)
         enddo
         do ll=1,nls
            write(10) (vert(nn,ll)%xyz,nn=1,para(ll)%nvert)
         enddo
         do ll=1,nls
            write(10) (elem(nn,ll)%vert,nn=1,para(ll)%nelem)
         enddo
 
      else

         do ll=1,nls
            if(ibseq(ll).eq.1) then
            write(10) (vert(nn,ll)%xyz,nn=1,para(ll)%nvert)
            endif
         enddo

      endif
      

      close(10)


      return
      end



C---------------------------------------------------------------------------------------
C 
C     IBSurface_Sequence : Read ASCII STL Files
C                          Find Duplicated Vertices & Elements 
C
C---------------------------------------------------------------------------------------
      subroutine ibsurface_sequence(fileseq,paraorgmx,paraorg,indxmapping,ivert_seg2org,nlsf,nlsfseq,nlsfile,
     c                              paramx,vert,elem,edge,para,nls,ibseq,nseq,nlsfscale,fscalenls,icheck)
 
      USE GenSurf
      type(surfpar):: paraorgmx,paraorg(nlsfile)
      type(surfpar):: paramx,para(nls)
      type(vertex) :: vert(paramx%nvert,nls)
      type(element):: elem(paramx%nelem,nls)
      type(edges)  :: edge(paramx%nedge,nls)

      integer indxmapping(paraorgmx%nvert,2,nlsfile)
      integer ivert_seg2org(paraorgmx%nvert,nls)

      integer nlsf(nlsfile),nlsfseq(nlsfile),nlsfscale(nlsfile),ibseq(nls)
      integer nlsmn(nlsfile),nlsmx(nlsfile)
      real fscalenls(nlsfile)

      real,allocatable,dimension(:,:) :: datapts

      character*80 fileseq(nlsfile)
      character*80 filestl,foutl

      allocate(datapts(paraorgmx%nvert,3))

      nlsmn(1) = 1
      nlsmx(1) = nlsf(1)

      do nf=2,nlsfile
         nlsmn(nf) = nlsmn(nf-1) + 1
         nlsmx(nf) = nlsmx(nf-1) + nlsf(nf)
      enddo

      do nf=1,nlsfile

      if(nlsfseq(nf).eq.1) then

C---- Read STL data 
      write(foutl,990) nseq
990   format(i10)
      filestl=trim(fileseq(nf))//trim(adjustl(foutl))//'.stl'

      call read_stl(filestl,datapts,paraorg(nf)%nelem,paraorgmx%nvert)
      write(6,*) 'Original parameters, #seq #vert,#elem',nseq,paraorg(nf)%nelem,paraorg(nf)%nvert

      do ll=nlsmn(nf),nlsmx(nf)
         do nn=1,para(ll)%nvert
            mm = ivert_seg2org(nn,ll)
            indx = indxmapping(mm,1,nf)
            vert(nn,ll)%xyz(:) = datapts(indx,:)
         enddo
      enddo


      if(nlsfscale(nf).eq.1) then
      do ll=nlsmn(nf),nlsmx(nf)
         call scale_ibsurf(paramx,vert(1,ll),para(ll),fscalenls(nf))
      enddo
      endif

      endif ! end of nlsfseq

      enddo ! end of nf

      if(icheck.eq.1) call ibsurf_conn(paramx,vert,elem,edge,para,ibseq,nls,nseq)
 
      deallocate(datapts)

      return
      end


C---------------------------------------------------------------------------------------
C 
C     Scale_IB_Data : Scaling Data or Shift IB Position
C
C---------------------------------------------------------------------------------------
      subroutine scale_ibsurf(paramx,vert,para,scale)
      USE GenSurf
      type(surfpar):: paramx,para
      type(vertex) :: vert(paramx%nvert)

      do nn=1,para%nvert
         vert(nn)%xyz(1) =  scale*vert(nn)%xyz(1)
         vert(nn)%xyz(2) =  scale*vert(nn)%xyz(2)
         vert(nn)%xyz(3) =  scale*vert(nn)%xyz(3)
      enddo

      xmax = -1e30
      xmin =  1e30
      ymax = -1e30
      ymin =  1e30
      zmax = -1e30
      zmin =  1e30

      do nn=1,para%nvert
         xmax = max(xmax,vert(nn)%xyz(1))
         xmin = min(xmin,vert(nn)%xyz(1))
         ymax = max(ymax,vert(nn)%xyz(2))
         ymin = min(ymin,vert(nn)%xyz(2))
         zmax = max(zmax,vert(nn)%xyz(3))
         zmin = min(zmin,vert(nn)%xyz(3))
      enddo

      xctr = 0.5*(xmax+xmin)
      yctr = 0.5*(ymax+ymin)
      zctr = 0.5*(zmax+zmin)

c      write(6,*) 'Scaled IBsurf xctr,xmin,xmax,yctr,ymin,ymax,zctr,zmax,zmin'
c      write(6,999) xctr,xmin,xmax,yctr,ymin,ymax,zctr,zmax,zmin
999   format(30e15.7)

      return
      end


C---------------------------------------------------------------------------------------
C
C     Read_ibsurf_bin : Read IBSurface in binary Format w or w/o elem_vert
C
C---------------------------------------------------------------------------------------
      subroutine read_ibsurf_bin(filein,paramx,vert,elem,para,ibseq,nls,nseq,ibscale,ibfactor,xyzoff)
      USE GenSurf
      type(surfpar):: paramx,para(nls)
      type(vertex) :: vert(paramx%nvert,nls)
      type(element):: elem(paramx%nelem,nls)
      integer ibseq(nls),ibscale(nls)
      real ibfactor(nls),xyzoff(nls,3)

      character*80 filein,fout,foutl

      write(foutl,990) nseq
990   format(i10)

      fout=trim(filein)//trim(adjustl(foutl))//'.bin'

      open(10,file=trim(fout),form='unformatted')

      if(nseq.eq.0) then
         read(10) nls
         do ll=1,nls
            read(10) para(ll)
         enddo

         do ll=1,nls
            read(10) (vert(nn,ll)%xyz,nn=1,para(ll)%nvert)
         enddo
         do ll=1,nls
            read(10) (elem(nn,ll)%vert,nn=1,para(ll)%nelem)
         enddo

      else

         do ll=1,nls
            if(ibseq(ll).eq.1) then
            read(10) (vert(nn,ll)%xyz,nn=1,para(ll)%nvert)
            endif
         enddo

      endif
      close(10)

C *** scaling and translation
      do ll=1,nls
         if(ibscale(ll).eq.1.and.(ibseq(ll).eq.1.or.nseq.eq.0)) then
         do nn=1,para(ll)%nvert
            vert(nn,ll)%xyz(:) = ibfactor(ll)*vert(nn,ll)%xyz(:) + xyzoff(ll,:)
         enddo
         endif
      enddo
C ***

      ibcheck = 1

      if(ibcheck.eq.1) then
      do ll=1,nls
      if(ibseq(ll).eq.1.or.nseq.eq.0) then

      xmax = -1e30
      xmin =  1e30
      ymax = -1e30
      ymin =  1e30
      zmax = -1e30
      zmin =  1e30

      do nn=1,para(ll)%nvert
         xmax = max(xmax,vert(nn,ll)%xyz(1))
         xmin = min(xmin,vert(nn,ll)%xyz(1))
         ymax = max(ymax,vert(nn,ll)%xyz(2))
         ymin = min(ymin,vert(nn,ll)%xyz(2))
         zmax = max(zmax,vert(nn,ll)%xyz(3))
         zmin = min(zmin,vert(nn,ll)%xyz(3))
      enddo

      xctr = 0.5*(xmax+xmin)
      yctr = 0.5*(ymax+ymin)
      zctr = 0.5*(zmax+zmin)

c      write(6,*) 'Scaled IBsurf xctr,xmin,xmax,yctr,ymin,ymax,zctr,zmax,zmin'
c      write(6,999) xctr,xmin,xmax,yctr,ymin,ymax,zctr,zmax,zmin
999   format(30e15.7)
      endif
      enddo

      endif

      return
      end


C--------------------------------------------------------------------------------------
C
C     ibsurf_conn : Calcuate connectivity, shared_vertices, etc.
C
C---------------------------------------------------------------------------------------
      subroutine ibsurf_conn(paramx,vert,elem,edge,para,ibseq,nls,nseq)
      USE GenSurf
      type(surfpar):: paramx,para(nls)
      type(vertex) :: vert(paramx%nvert,nls)
      type(element):: elem(paramx%nelem,nls)
      type(edges)  :: edge(paramx%nedge,nls)
      integer ibseq(nls)

      powervalue = 1.0/7.0
c      powervalue = 1.0
      if(nseq.eq.0) then
         do ll=1,nls
            call shared_vertices(paramx,vert(1,ll),elem(1,ll),para(ll))
            call find_connectivity(paramx,vert(1,ll),elem(1,ll),para(ll))
            call conn_correction(paramx,vert(1,ll),elem(1,ll),para(ll))
            call shared_edge(paramx,vert(1,ll),elem(1,ll),edge(1,ll),para(ll))
            call face_normal(paramx,vert(1,ll),elem(1,ll),para(ll))
            call edge_normal(paramx,vert(1,ll),elem(1,ll),edge(1,ll),para(ll))
            call pseudo_normal(paramx,vert(1,ll),elem(1,ll),edge(1,ll),para(ll))
 
            vert(:,:)%power = powervalue
            elem(:,:)%power = powervalue
            edge(:,:)%power = powervalue
         enddo

      else

      do ll=1,nls
         if(ibseq(ll).eq.1) then
         call face_normal(paramx,vert(1,ll),elem(1,ll),para(ll))
         call edge_normal(paramx,vert(1,ll),elem(1,ll),edge(1,ll),para(ll))
         call pseudo_normal(paramx,vert(1,ll),elem(1,ll),edge(1,ll),para(ll))
         endif
      enddo
    
      endif

      return
      end


C---------------------------------------------------------------------------------------
C
C     calc_ibvel : calc_ibvel and ibacc
C
C---------------------------------------------------------------------------------------
      subroutine calc_ibvel(paramx,vert,elem,edge,para,ibseq,nls,nseq,dtadvec,velmax,accmax,posmax,posmin)
      USE GenSurf
      type(surfpar):: paramx,para(nls)
      type(vertex) :: vert(paramx%nvert,nls)
      type(element):: elem(paramx%nelem,nls)
      type(edges)  :: edge(paramx%nedge,nls)
      integer ibseq(nls)

      real velmax(nls),accmax(nls),posmax(nls,3),posmin(nls,3)

      posmin = 1e30
      posmax =-1e30
      velmax = 0.0
      accmax = 0.0

      if(nseq.eq.0) then

         do ll=1,nls
         do nn=1,para(ll)%nvert
            vert(nn,ll)%xyz0 = vert(nn,ll)%xyz
            vert(nn,ll)%vel = 0.0
            vert(nn,ll)%vel0 = vert(nn,ll)%vel
            vert(nn,ll)%acc = 0.0

            velmag = sqrt(vert(nn,ll)%vel(1)**2+vert(nn,ll)%vel(2)**2+vert(nn,ll)%vel(3)**2)
            accmag = sqrt(vert(nn,ll)%acc(1)**2+vert(nn,ll)%acc(2)**2+vert(nn,ll)%acc(3)**2)
            velmax(ll) = max(velmax(ll),velmag)
            accmax(ll) = max(accmax(ll),accmag)
            posmax(ll,1:3) = max(posmax(ll,1:3),vert(nn,ll)%xyz(1:3))
            posmin(ll,1:3) = min(posmin(ll,1:3),vert(nn,ll)%xyz(1:3))

         enddo
         enddo

      else

         do ll=1,nls
         do nn=1,para(ll)%nvert

            vert(nn,ll)%vel = (vert(nn,ll)%xyz - vert(nn,ll)%xyz0)/dtadvec
            vert(nn,ll)%acc = (vert(nn,ll)%vel - vert(nn,ll)%vel0)/dtadvec
            vert(nn,ll)%vel0 = vert(nn,ll)%vel
            vert(nn,ll)%xyz0 = vert(nn,ll)%xyz

            velmag = sqrt(vert(nn,ll)%vel(1)**2+vert(nn,ll)%vel(2)**2+vert(nn,ll)%vel(3)**2)
            accmag = sqrt(vert(nn,ll)%acc(1)**2+vert(nn,ll)%acc(2)**2+vert(nn,ll)%acc(3)**2)

            velmax(ll) = max(velmax(ll),velmag)
            accmax(ll) = max(accmax(ll),accmag)
            posmax(ll,1:3) = max(posmax(ll,1:3),vert(nn,ll)%xyz(1:3))
            posmin(ll,1:3) = min(posmin(ll,1:3),vert(nn,ll)%xyz(1:3))

         enddo
         enddo

         write(6,*) 'Maximum velocity',velmax
         write(6,*) 'Maximum acceleration',accmax

      endif

999   format(30e15.7)
      return
      end

