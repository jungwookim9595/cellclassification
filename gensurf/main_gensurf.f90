!----------------------------------------------------------------------------------
!
!     This program creates segmented immersed surfaces in STL, TecPlot, Node, Bin formats
!     from STL file sequences
!
!     Jung-Il Choi, Ph.D.,
!     Research Assistant Professor
!     Dept. Mechanical & Aerospace Engineering
!     North Carolina State University
!
!     Modernized to Fortran 90 free-format (2026-03-23)
!
!----------------------------------------------------------------------------------
program main
  use GenSurf
  implicit none

  type(surfpar) :: ibparaorgmx
  type(surfpar), allocatable, dimension(:) :: ibparaorg
  type(vertexorg), allocatable, dimension(:,:) :: ibvertorg
  type(elementorg), allocatable, dimension(:,:) :: ibelemorg

  integer, allocatable, dimension(:,:,:) :: indxmapping
  integer, allocatable, dimension(:,:) :: ielemmapping, ielemseg

  integer, allocatable, dimension(:,:) :: ielem_seg2org, ivert_seg2org, ivert_org2seg

  type(surfpar) :: ibparamx
  type(surfpar), allocatable, dimension(:) :: ibpara
  type(vertex), allocatable, dimension(:,:) :: ibvert
  type(element), allocatable, dimension(:,:) :: ibelem
  type(edges), allocatable, dimension(:,:) :: ibedge

  character(len=80) :: orgstlfile
  integer :: nlsfile, nsolid, solidnum
  character(len=80), allocatable, dimension(:) :: filestl, fileseq, filesolidname
  integer, allocatable, dimension(:) :: nlsf, nlsfseq, nlsfscale, ibseq
  integer, allocatable, dimension(:) :: nlsmap, ncelem, ncvert
  real(8), allocatable, dimension(:) :: fscalenls

  character(len=80) :: fileout
  integer :: nls, nseq, nseqmx, nc, nn, mm, ll, nf
  integer :: icheck, icheckfile

  write(6, *) 'input STL filename'
  read(5, *) orgstlfile

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

  do solidnum = 1, nsolid

     filestl(1) = filesolidname(solidnum)
     nlsfseq = 0
     nlsfscale = 0
     fscalenls = 1.0d0
     nseqmx = 0

     fileout = 'output/ibsurf'

     ! Define element count
     do nf = 1, nlsfile
        call define_element(filestl(nf), ibparaorg(nf)%nelem)
        ibparaorg(nf)%nvert = ibparaorg(nf)%nelem * 3
        ibparaorg(nf)%nedge = ibparaorg(nf)%nelem * 3
     end do

     ibparaorgmx%nvert = 1
     ibparaorgmx%nelem = 1
     ibparaorgmx%nedge = 1

     do nf = 1, nlsfile
        ibparaorgmx%nvert = max(ibparaorgmx%nvert, ibparaorg(nf)%nvert)
        ibparaorgmx%nelem = max(ibparaorgmx%nelem, ibparaorg(nf)%nelem)
        ibparaorgmx%nedge = max(ibparaorgmx%nedge, ibparaorg(nf)%nedge)
     end do

     allocate(ibvertorg(ibparaorgmx%nvert, nlsfile))
     allocate(ibelemorg(ibparaorgmx%nelem, nlsfile))
     allocate(indxmapping(ibparaorgmx%nvert, 2, nlsfile))
     allocate(ielemmapping(ibparaorgmx%nelem, nlsfile))
     allocate(ielemseg(ibparaorgmx%nelem, nlsfile))

     ! Convert STL to IBSurface
     do nf = 1, nlsfile
        call convert2ibsurface(filestl(nf), ibparaorgmx, ibvertorg(1, nf), ibelemorg(1, nf), &
             ibparaorg(nf), indxmapping(1, 1, nf), ielemmapping(1, nf), ielemseg(1, nf), nlsf(nf))
     end do
     
     ! Compute total number of level sets
     nls = 0
     do nf = 1, nlsfile
        nls = nls + nlsf(nf)
     end do

     allocate(nlsmap(nls))
     nc = 0
     do nn = 1, nlsfile
        do mm = 1, nlsf(nn)
           nc = nc + 1
           nlsmap(nc) = nn
        end do
     end do

     allocate(ielem_seg2org(ibparaorgmx%nelem, nls))
     allocate(ivert_seg2org(ibparaorgmx%nvert, nls))
     allocate(ivert_org2seg(ibparaorgmx%nvert, nlsfile))

     allocate(ncelem(nls))
     allocate(ncvert(nls))

     write(6,*) 'DEBUG: calling ibsurface_seg_param'; flush(6)
     call ibsurface_seg_param(ibparaorgmx, ibvertorg, ibelemorg, ibparaorg, &
          indxmapping, ielemmapping, ielemseg, ncelem, ncvert, &
          ielem_seg2org, ivert_seg2org, ivert_org2seg, nls, nlsmap, nlsf, nlsfile)
     write(6,*) 'DEBUG: ibsurface_seg_param done, ncelem(1)=',ncelem(1),' ncvert(1)=',ncvert(1); flush(6)

     allocate(ibpara(nls))

     ibparamx%nvert = 1
     ibparamx%nelem = 1
     ibparamx%nedge = 1

     do ll = 1, nls
        ibpara(ll)%nelem = ncelem(ll)
        ibpara(ll)%nvert = ncvert(ll)
        ibpara(ll)%nedge = ncelem(ll) * 3
        ibparamx%nvert = max(ibparamx%nvert, ibpara(ll)%nvert)
        ibparamx%nelem = max(ibparamx%nelem, ibpara(ll)%nelem)
        ibparamx%nedge = max(ibparamx%nedge, ibpara(ll)%nedge)
     end do

     write(6,*) 'DEBUG: ibparamx nvert/nelem/nedge=',ibparamx%nvert,ibparamx%nelem,ibparamx%nedge; flush(6)
     allocate(ibvert(ibparamx%nvert, nls))
     write(6,*) 'DEBUG: ibvert allocated'; flush(6)
     allocate(ibelem(ibparamx%nelem, nls))
     write(6,*) 'DEBUG: ibelem allocated'; flush(6)
     allocate(ibedge(ibparamx%nedge, nls))
     write(6,*) 'DEBUG: ibedge allocated, calling ibsurface_seg'; flush(6)

     call ibsurface_seg(ibparaorgmx, ibvertorg, ibelemorg, ibparaorg, &
          ibparamx, ibvert, ibelem, ibedge, ibpara, &
          ielem_seg2org, ivert_seg2org, ivert_org2seg, nls, nlsmap, nlsf, nlsfile)
     write(6,*) 'DEBUG: ibsurface_seg done'; flush(6)

     write(6,*) 'DEBUG: deallocating org arrays'; flush(6)
     deallocate(ibvertorg, ibelemorg, ielemmapping, ielemseg)
     deallocate(ielem_seg2org, ivert_org2seg)
     deallocate(ncelem, ncvert)
     deallocate(nlsmap)
     write(6,*) 'DEBUG: deallocations done'; flush(6)

     write(6, *) 'Segmented levelset info.'
     do ll = 1, nls
        write(6, *) ll, ibpara(ll)%nvert, ibpara(ll)%nelem, ibpara(ll)%nedgea
     end do

     allocate(ibseq(nls))
     ibseq = 0

     ! Segmentation and diagnostics are completed!!!

     write(6, *)
     write(6, *) 'Now reading sequences'
     write(6, *)

     icheck = 1
     icheckfile = 1
     do nseq = 0, nseqmx
        write(6, *) 'Reading sequence', nseq

        call ibsurface_sequence(fileseq, ibparaorgmx, ibparaorg, indxmapping, ivert_seg2org, nlsf, nlsfseq, nlsfile, &
             ibparamx, ibvert, ibelem, ibedge, ibpara, nls, ibseq, nseq, nlsfscale, fscalenls, icheck)

        call output_bin(fileout, ibparamx, ibvert, ibelem, ibpara, ibseq, nls, nseq)

        if (nseq == 0 .and. icheckfile == 1) then
           call output_stl(ibparamx, ibvert(1, 1), ibelem(1, 1), ibpara(1), solidnum)
           call output_node(ibparamx, ibvert(1, 1), ibelem(1, 1), ibedge(1, 1), ibpara(1), solidnum)
           call output_custom(ibparamx, ibvert(1, 1), ibelem(1, 1), ibedge(1, 1), ibpara(1), solidnum)
        end if

     end do

     deallocate(indxmapping, ivert_seg2org)
     deallocate(ibpara, ibvert, ibelem, ibedge)
     deallocate(ibseq)

  end do

  deallocate(filesolidname)
  deallocate(filestl, nlsf, ibparaorg, fileseq, nlsfseq, nlsfscale, fscalenls)

  stop
end program main
