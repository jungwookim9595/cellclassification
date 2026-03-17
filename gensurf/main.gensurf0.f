C----------------------------------------------------------------------------------
C
C     This program creates segmented immersed surfaces in STL, TecPlot, Node, Bin formats
C     from STL file sequences
C
C     Jung-Il Choi, Ph.D.,
C     Research Assistant Professor
C     Dept. Mechanical & Aerospace Engineering
C     North Carolina State University
C
C----------------------------------------------------------------------------------
      Program main
      USE GenSurf

      type(surfpar) :: ibparaorgmx
      type(surfpar),allocatable,dimension(:):: ibparaorg
      type(vertexorg),allocatable,dimension(:,:) :: ibvertorg
      type(elementorg),allocatable,dimension(:,:):: ibelemorg

      integer, allocatable,dimension(:,:,:) :: indxmapping
      integer, allocatable,dimension(:,:) :: ielemmapping,ielemseg

      integer, allocatable,dimension(:,:) :: ielem_seg2org,ivert_seg2org,ivert_org2seg

      type(surfpar) :: ibparamx
      type(surfpar),allocatable,dimension(:)  :: ibpara
      type(vertex),allocatable,dimension(:,:) :: ibvert
      type(element),allocatable,dimension(:,:):: ibelem
      type(edges),allocatable,dimension(:,:)  :: ibedge

      character*80                           :: orgstlfile
      integer                                :: nlsfile, nsolid, solidnum
      character*80, allocatable,dimension(:) :: filestl,fileseq, filesolidname
      integer, allocatable,dimension(:) :: nlsf,nlsfseq,nlsfscale,ibseq
      integer, allocatable,dimension(:) :: nlsmap,ncelem,ncvert
      real, allocatable,dimension(:) :: fscalenls

      character*80 fileout

      write(6,*) 'input STL filename'
      read (6,*) orgstlfile


      call count_solids(orgstlfile, nsolid)
      print *, 'total number of solids :', nsolid

      allocate(filesolidname(nsolid))
      call devide_file(orgstlfile, nsolid, filesolidname)

      nlsfile = 1
      allocate(filestl(nlsfile))
      allocate(nlsf(nlsfile))
      allocate(ibparaorg(nlsfile))

      allocate(fileseq(nlsfile))
      allocate(nlsfseq(nlsfile))
      allocate(nlsfscale(nlsfile))
      allocate(fscalenls(nlsfile))

      do solidnum = 1,nsolid
      filestl(1) = filesolidname(solidnum)
      
      nlsfseq = 0
   
      nlsfscale = 0
      fscalenls = 1.0

      nseqmx = 0
      fileout = 'output/ibsurf'

C---- Define the maximum number of the element
      ibparaorgmx%nvert = 1
      ibparaorgmx%nelem = 1
      ibparaorgmx%nedge = 1

      do nn=1,nlsfile
         call Define_Element(filestl(nn),nelemorg)
         nvertorg = 3*nelemorg
         nedgeorg = 3*nelemorg

         ibparaorg(nn)%nvert = nvertorg
         ibparaorg(nn)%nelem = nelemorg
         ibparaorg(nn)%nedge = nedgeorg

         ibparaorgmx%nvert = max(ibparaorgmx%nvert,nvertorg)
         ibparaorgmx%nelem = max(ibparaorgmx%nelem,nelemorg)
         ibparaorgmx%nedge = max(ibparaorgmx%nedge,nedgeorg)
      enddo

      allocate(indxmapping(ibparaorgmx%nvert,2,nlsfile))
      allocate(ielemmapping(ibparaorgmx%nelem,nlsfile))
      allocate(ielemseg(ibparaorgmx%nelem,nlsfile))

      allocate(ibvertorg(ibparaorgmx%nvert,nlsfile))
      allocate(ibelemorg(ibparaorgmx%nelem,nlsfile))
      nls = 0
      do nn=1,nlsfile
         call convert2ibsurface(filestl(nn),ibparaorgmx,ibvertorg(1,nn),ibelemorg(1,nn),ibparaorg(nn),
     c                          indxmapping(1,1,nn),ielemmapping(1,nn),ielemseg(1,nn),nlsf(nn))
         nls = nls + nlsf(nn)
      enddo
      write(6,*) 'Total # of levelsets',nls

      allocate(nlsmap(nls))
      nc = 0
      do nn=1,nlsfile
         do mm=1,nlsf(nn)
            nc = nc + 1
            nlsmap(nc) = nn
         enddo
      enddo

      allocate(ielem_seg2org(ibparaorgmx%nelem,nls))
      allocate(ivert_seg2org(ibparaorgmx%nvert,nls))
      allocate(ivert_org2seg(ibparaorgmx%nvert,nlsfile))

      allocate(ncelem(nls))
      allocate(ncvert(nls))

      call ibsurface_seg_param(ibparaorgmx,ibvertorg,ibelemorg,ibparaorg,
     c                         indxmapping,ielemmapping,ielemseg,ncelem,ncvert,
     c                         ielem_seg2org,ivert_seg2org,ivert_org2seg,nls,nlsmap,nlsf,nlsfile)

      allocate(ibpara(nls))

      ibparamx%nvert = 1
      ibparamx%nelem = 1
      ibparamx%nedge = 1

      do ll=1,nls
         ibpara(ll)%nelem = ncelem(ll)
         ibpara(ll)%nvert = ncvert(ll)
         ibpara(ll)%nedge = ncelem(ll)*3
         ibparamx%nvert = max(ibparamx%nvert,ibpara(ll)%nvert)
         ibparamx%nelem = max(ibparamx%nelem,ibpara(ll)%nelem)
         ibparamx%nedge = max(ibparamx%nedge,ibpara(ll)%nedge)
      enddo

      allocate(ibvert(ibparamx%nvert,nls))
      allocate(ibelem(ibparamx%nelem,nls))
      allocate(ibedge(ibparamx%nedge,nls))

      call ibsurface_seg(ibparaorgmx,ibvertorg,ibelemorg,ibparaorg,
     c                   ibparamx,ibvert,ibelem,ibedge,ibpara,
     c                   ielem_seg2org,ivert_seg2org,ivert_org2seg,nls,nlsmap,nlsf,nlsfile)

      deallocate(ibvertorg,ibelemorg,ielemmapping,ielemseg)
      deallocate(ielem_seg2org,ivert_org2seg)
      deallocate(ncelem,ncvert)
      deallocate(nlsmap)

      write(6,*) 'Segmented levelset info.'
      do ll=1,nls
         write(6,*) ll,ibpara(ll)%nvert,ibpara(ll)%nelem,ibpara(ll)%nedgea
      enddo

      allocate(ibseq(nls))
      ibseq = 0
c      ibseq(1) = 1

C *** Segmentation and dignostics are completed!!!

      write(6,*)
      write(6,*) 'Now reading sequences'
      write(6,*)

      icheck = 1
      icheckfile = 1
      do nseq=0,nseqmx
         write(6,*) 'Reading sequence',nseq

         call ibsurface_sequence(fileseq,ibparaorgmx,ibparaorg,indxmapping,ivert_seg2org,nlsf,nlsfseq,nlsfile,
     c                           ibparamx,ibvert,ibelem,ibedge,ibpara,nls,ibseq,nseq,nlsfscale,fscalenls,icheck)

         call output_bin(fileout,ibparamx,ibvert,ibelem,ibpara,ibseq,nls,nseq)

         if(nseq.eq.0.and.icheckfile.eq.1) then
            call output_stl(ibparamx,ibvert(1,1),ibelem(1,1),ibpara(1),solidnum)
            call output_node(ibparamx,ibvert(1,1),ibelem(1,1),ibedge(1,1),ibpara(1),solidnum)
            call output_custom(ibparamx,ibvert(1,1),ibelem(1,1),ibedge(1,1),ibpara(1),solidnum)
         !do ll=1,nls
         !   call output_stl(ibparamx,ibvert(1,ll),ibelem(1,ll),ibpara(ll),ll)
         !   call output_node(ibparamx,ibvert(1,ll),ibelem(1,ll),ibedge(1,ll),ibpara(ll),ll)
         !   !call output_particle_method(ibparamx,ibvert(1,ll),ibelem(1,ll),ibedge(1,ll),ibpara(ll),ll)
         !enddo
         !   !call output_particle_method_all(ibparamx,ibvert,ibelem,ibedge,ibpara,nls)
         !   call output_tecplot(ibparamx,ibvert,ibelem,ibpara,nls,nseq)
         endif

      enddo

      deallocate(indxmapping,ivert_seg2org)
      deallocate(ibpara,ibvert,ibelem,ibedge)
      deallocate(ibseq)

      enddo

      deallocate(filesolidname)
      deallocate(filestl, nlsf, ibparaorg, fileseq, nlsfseq, nlsfscale, fscalenls)

      stop
      end

