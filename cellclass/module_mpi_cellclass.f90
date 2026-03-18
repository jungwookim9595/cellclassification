module surfvar
  implicit none

  integer                                             ::  nshare = 200
  type ibvert
    double precision, dimension(3)                    ::  xyz, norm
    integer         , dimension(200)                  ::  share ! dimension = nshare
  end type ibvert
  
  type ibedge
    double precision, dimension(3)                    ::  xyz, norm
    integer         , dimension(2)                    ::  vert
  end type ibedge
  
  type ibelem
    double precision, dimension(3)                    ::  xyz, norm
    integer         , dimension(3)                    ::  vert, edge
  end type ibelem
  
  type ibpara
    integer                                           ::  nvert, nedge, nelem
  end type ibpara
  
end module surfvar

module mpi_cellclass
    use mpi
    use global
    use mpi_subdomain
    use mpi_topology
    use surfvar

    implicit none

    integer                                           ::  nstlf    = 1
    integer                                           ::  nsearch  = 1
    double precision                                  ::  dist_CAD = 0.01

    type(ibpara)                                      ::  para
    type(ibvert), allocatable, dimension(:)           ::  vert
    type(ibelem), allocatable, dimension(:)           ::  elem
    type(ibedge), allocatable, dimension(:)           ::  edge

    double precision , allocatable, dimension(:,:,:)  :: sgned_g, heavi_g

    contains

    subroutine var_surf_make
        implicit none

        allocate( vert(para%nvert), elem(para%nelem), edge(para%nedge) )

    end subroutine var_surf_make

    subroutine var_surf_clean
        deallocate( vert, elem, edge )
    end subroutine var_surf_clean

    subroutine var_cc_make
        implicit none

        allocate( sgned_g(0:n1sub,0:n2sub,0:n3sub), heavi_g(0:n1sub,0:n2sub,0:n3sub) )

        sgned_g = 999.0
    end subroutine var_cc_make

    subroutine var_cc_clean
        deallocate( sgned_g, heavi_g )
    end subroutine var_cc_clean

    subroutine check_npt_total(filestl)
        implicit none

        character(len=80)   ::  filestl

        integer             ::  i, j
        character(len=10)   ::  cname
        integer   ::  dummy

        open(1,file = trim(filestl))

          read(1,*) cname, para%nvert
          do i = 1,para%nvert
              read(1,*) 
              read(1,*) 
              read(1,*) 
              read(1,*) cname, dummy
              do j = 1, dummy
                read(1,*)
              enddo
              read(1,*)
          enddo

          read(1,*) cname
          if(trim(cname).NE.'endvertex') then
            write(*,*) 'read vertex error'
            stop
          endif

          read(1,*) cname, para%nedge
          do i = 1,para%nedge
              read(1,*) 
              read(1,*) 
              read(1,*) 
              read(1,*) 
              read(1,*) 
          enddo

          read(1,*) cname
          if(trim(cname).NE.'endedge') then
            write(*,*) 'read edge error'
            stop
          endif

          read(1,*) cname, para%nelem
          do i = 1,para%nelem
              read(1,*)
              read(1,*)
              read(1,*)
              read(1,*)
              read(1,*)
              read(1,*)
          enddo

          read(1,*) cname
          if(trim(cname).NE.'endelem') then
            write(*,*) 'read element error'
            stop
          endif

        close(1)

    end subroutine check_npt_total

    subroutine read_stl(filestl)
        implicit none

        character(len=80)   :: filestl

        integer             ::  i, j, nn, mm
        character(len=10)   ::  cname

        open(2,file = trim(filestl))

          read(2,*)
          do i = 1,para%nvert
              read(2,*) 
              read(2,*) cname, vert(i)%xyz(1), vert(i)%xyz(2), vert(i)%xyz(3)
              read(2,*) cname, vert(i)%norm(1), vert(i)%norm(2), vert(i)%norm(3)
              read(2,*) cname, vert(i)%share(1)
              do j = 1,vert(i)%share(1)
                read(2,*) vert(i)%share(j+1)
              enddo
              read(2,*)
          enddo

          read(2,*) cname
          if(trim(cname).NE.'endvertex') then
            write(*,*) 'read vertex error'
            stop
          endif

          read(2,*)
          do i = 1,para%nedge
              read(2,*) 
              read(2,*) cname, edge(i)%xyz(1), edge(i)%xyz(2), edge(i)%xyz(3)
              read(2,*) cname, edge(i)%norm(1), edge(i)%norm(2), edge(i)%norm(3)
              read(2,*) cname, edge(i)%vert(1), edge(i)%vert(2)
              read(2,*) 
          enddo
          
          read(2,*) cname
          if(trim(cname).NE.'endedge') then
            write(*,*) 'read edge error'
            stop
          endif

          read(2,*) 
          do i = 1,para%nelem
              read(2,*)
              read(2,*) cname, elem(i)%xyz(1), elem(i)%xyz(2), elem(i)%xyz(3)
              read(2,*) cname, elem(i)%norm(1), elem(i)%norm(2), elem(i)%norm(3)
              read(2,*) cname, elem(i)%vert(1), elem(i)%vert(2), elem(i)%vert(3)
              read(2,*) cname, elem(i)%edge(1), elem(i)%edge(2), elem(i)%edge(3)
              read(2,*)
          enddo
          
          read(2,*) cname
          if(trim(cname).NE.'endelem') then
            write(*,*) 'read element error'
            stop
          endif

        close(2)

    end subroutine read_stl
    
     
    !=========================================================================
    ! Phase 1: Find narrowband cells using triangle AABB intersection
    !=========================================================================
    subroutine find_narrowband_cells(band_mask, bandwidth)
      implicit none
 
      logical, dimension(0:n1sub,0:n2sub,0:n3sub), intent(out) :: band_mask
      double precision, intent(in)                              :: bandwidth
 
      integer          :: ielem, i_lo, i_hi, j_lo, j_hi, k_lo, k_hi
      double precision :: xmin, xmax, ymin, ymax, zmin, zmax
      double precision :: p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z
 
      band_mask = .false.
 
      do ielem = 1, para%nelem
        ! Get triangle vertex coordinates
        p1x = vert(elem(ielem)%vert(1))%xyz(1)
        p1y = vert(elem(ielem)%vert(1))%xyz(2)
        p1z = vert(elem(ielem)%vert(1))%xyz(3)
        p2x = vert(elem(ielem)%vert(2))%xyz(1)
        p2y = vert(elem(ielem)%vert(2))%xyz(2)
        p2z = vert(elem(ielem)%vert(2))%xyz(3)
        p3x = vert(elem(ielem)%vert(3))%xyz(1)
        p3y = vert(elem(ielem)%vert(3))%xyz(2)
        p3z = vert(elem(ielem)%vert(3))%xyz(3)
 
        ! Compute AABB expanded by bandwidth
        xmin = min(p1x, p2x, p3x) - bandwidth
        xmax = max(p1x, p2x, p3x) + bandwidth
        ymin = min(p1y, p2y, p3y) - bandwidth
        ymax = max(p1y, p2y, p3y) + bandwidth
        zmin = min(p1z, p2z, p3z) - bandwidth
        zmax = max(p1z, p2z, p3z) + bandwidth
 
        ! Find grid index ranges using binary search
        call find_grid_index(x1_sub, n1sub, xmin, i_lo)
        call find_grid_index(x1_sub, n1sub, xmax, i_hi)
        call find_grid_index(x2_sub, n2sub, ymin, j_lo)
        call find_grid_index(x2_sub, n2sub, ymax, j_hi)
        call find_grid_index(x3_sub, n3sub, zmin, k_lo)
        call find_grid_index(x3_sub, n3sub, zmax, k_hi)
 
        ! Clamp to subdomain bounds
        i_lo = max(i_lo, 0);  i_hi = min(i_hi, n1sub)
        j_lo = max(j_lo, 0);  j_hi = min(j_hi, n2sub)
        k_lo = max(k_lo, 0);  k_hi = min(k_hi, n3sub)
 
        ! Mark cells in the narrowband
        if (i_lo <= i_hi .and. j_lo <= j_hi .and. k_lo <= k_hi) then
          band_mask(i_lo:i_hi, j_lo:j_hi, k_lo:k_hi) = .true.
        endif
      enddo
 
    end subroutine find_narrowband_cells
 
    !=========================================================================
    ! Utility: Binary search for grid index (find leftmost index where x_arr(idx) <= val)
    !=========================================================================
    subroutine find_grid_index(x_arr, nmax, val, idx)
      implicit none
 
      integer, intent(in)                        :: nmax
      double precision, dimension(0:nmax), intent(in) :: x_arr
      double precision, intent(in)               :: val
      integer, intent(out)                       :: idx
 
      integer :: lo, hi, mid
 
      ! Handle out-of-bounds
      if (val <= x_arr(0)) then
        idx = 0
        return
      endif
      if (val >= x_arr(nmax)) then
        idx = nmax
        return
      endif
 
      ! Binary search
      lo = 0
      hi = nmax
      do while (hi - lo > 1)
        mid = (lo + hi) / 2
        if (x_arr(mid) <= val) then
          lo = mid
        else
          hi = mid
        endif
      enddo
      idx = lo
 
    end subroutine find_grid_index
 
    !=========================================================================
    ! Phase 2: Compute SDF only for narrowband cells (reuses ANN logic)
    !=========================================================================
    subroutine compute_sdf_narrowband(sgned, band_mask)
      implicit none
 
      double precision, dimension(0:n1sub,0:n2sub,0:n3sub), intent(out) :: sgned
      logical, dimension(0:n1sub,0:n2sub,0:n3sub), intent(in)          :: band_mask
 
      integer                                               ::  i, j, k, LL, mm, id
      integer         , allocatable, dimension(:)           ::  indexarray, indexarray2, idmapping, idmapping2
      double precision, allocatable, dimension(:)           ::  distarray, distarray2
      double precision                                      ::  dist, distance, sign_SDF
      integer                                               ::  restriction
      double precision, dimension(3)                        ::  real_clst_pt, real_clst_norm, clst_pt, clst_norm
      double precision, dimension(3)                        ::  qp
      double precision, dimension(para%nvert)               ::  ptx, pty, ptz
      double precision                                      ::  check_dist
      integer                                               ::  band_count, total_count, loop_count
 
      allocate(indexarray(nsearch), idmapping(nsearch), distarray(nsearch))
      allocate(indexarray2(para%nvert), idmapping2(para%nvert), distarray2(para%nvert))
 
      ! Create ANN-tree
      ptx(:) = vert(:)%xyz(1)
      pty(:) = vert(:)%xyz(2)
      ptz(:) = vert(:)%xyz(3)
      call wldst_create_tree(ptx(1), pty(1), ptz(1), para%nvert)
 
      ! Initialize: band exterior gets large positive value (undetermined)
      sgned = 999.0d0
 
      band_count = 0
      total_count = (n1sub+1) * (n2sub+1) * (n3sub+1)
      loop_count = 0
      do k = 0, n3sub
      do j = 0, n2sub
      do i = 0, n1sub
 
        if (.not. band_mask(i,j,k)) cycle   ! Skip non-narrowband cells
 
        band_count = band_count + 1
        distance = 999.0d0
 
        ! Get distarray & index array from ANN-tree
        qp(1) = x1_sub(i)
        qp(2) = x2_sub(j)
        qp(3) = x3_sub(k)
        call wldst_getdists_xyz(qp(1), qp(2), qp(3), indexarray, distarray, nsearch)
 
        restriction = 0
        do LL = 1, nsearch
          if (distarray(LL) .LT. distarray(1)) then
            write(*,*) 'ANN search error'
            stop
          else
            check_dist = (ptx(indexarray(1)+1)-ptx(indexarray(LL)+1))**2 &
            &           +(pty(indexarray(1)+1)-pty(indexarray(LL)+1))**2 &
            &           +(ptz(indexarray(1)+1)-ptz(indexarray(LL)+1))**2
            check_dist = sqrt(check_dist)
          endif
          if (check_dist .LT. 2*dist_CAD) then
            restriction = restriction + 1
            idmapping(restriction) = LL
          endif
        enddo
 
        if (restriction .EQ. nsearch) then ! extend search
          restriction = 0
          call wldst_getdists_xyz(qp(1), qp(2), qp(3), indexarray2, distarray2, para%nvert)
          do mm = 1, para%nvert
            if (distarray2(mm) .LT. distarray2(1)) then
              write(*,*) 'ANN extend search error'
              stop
            else
              check_dist = (ptx(indexarray2(1)+1)-ptx(indexarray2(mm)+1))**2 &
              &           +(pty(indexarray2(1)+1)-pty(indexarray2(mm)+1))**2 &
              &           +(ptz(indexarray2(1)+1)-ptz(indexarray2(mm)+1))**2
              check_dist = sqrt(check_dist)
            endif
            if (check_dist .LT. 2*dist_CAD) then
              restriction = restriction + 1
              idmapping2(restriction) = mm
            endif
          enddo
        endif
 
        ! Find real closest point & distance
        if (restriction .LT. nsearch) then
          do LL = 1, restriction
            id = idmapping(LL)
            do mm = 1, vert(indexarray(id)+1)%share(1)
              call triangle_distance(qp, vert(indexarray(id)+1)%share(mm+1), dist, clst_pt, clst_norm)
              if (dist .LT. distance) then
                distance = dist
                real_clst_pt(:)   = clst_pt(:)
                real_clst_norm(:) = clst_norm(:)
              endif
            enddo
          enddo
        else
          do LL = 1, restriction
            id = idmapping2(LL)
            do mm = 1, vert(indexarray2(id)+1)%share(1)
              call triangle_distance(qp, vert(indexarray2(id)+1)%share(mm+1), dist, clst_pt, clst_norm)
              if (dist .LT. distance) then
                distance = dist
                real_clst_pt(:)   = clst_pt(:)
                real_clst_norm(:) = clst_norm(:)
              endif
            enddo
          enddo
        endif
 
        ! Get signed distance function
        sign_SDF     = (x1_sub(i)-real_clst_pt(1))*(real_clst_norm(1)) &
        &             +(x2_sub(j)-real_clst_pt(2))*(real_clst_norm(2)) &
        &             +(x3_sub(k)-real_clst_pt(3))*(real_clst_norm(3))
        if (abs(sign_SDF) > 0.0d0) then
          sign_SDF = sign_SDF / abs(sign_SDF)
        else
          sign_SDF = 1.0d0
        endif
 
        sgned(i,j,k) = sign_SDF * distance
        
        loop_count = loop_count + 1

        if (mod(loop_count,100)==0) then
          if(myrank==0) write(*,*) '  [Narrowband] Loop count: ', loop_count, ' / Total: ', total_count
        endif
      enddo
      enddo
      enddo
 
      if (myrank == 0) then
        write(*,'(A,I10,A,I10,A,F6.2,A)') '  [Narrowband] Band cells: ', band_count, &
          ' / Total: ', total_count, ' (', 100.0d0*dble(band_count)/dble(total_count), '%)'
      endif
 
      call wldst_delete_tree
      deallocate(indexarray, indexarray2, distarray, distarray2)
 
    end subroutine compute_sdf_narrowband
 
    !=========================================================================
    ! Phase 3: Fast Sweeping sign propagation for band exterior
    ! Propagates SDF sign from narrowband boundary outward using 8-direction sweeps
    !=========================================================================
    subroutine propagate_sign_sweep(sgned, band_mask)
      implicit none
 
      double precision, dimension(0:n1sub,0:n2sub,0:n3sub), intent(inout) :: sgned
      logical, dimension(0:n1sub,0:n2sub,0:n3sub), intent(in)             :: band_mask
 
      integer :: i, j, k, sweep, isweep
      integer :: i_start, i_end, i_step
      integer :: j_start, j_end, j_step
      integer :: k_start, k_end, k_step
      double precision :: neighbor_val, large_val
      integer :: num_sweeps
      logical :: changed
 
      large_val = 998.0d0  ! Threshold to detect unset values
 
      ! Perform 8-direction sweeps to propagate sign from band boundary
      ! Each sweep goes in one of 8 diagonal directions
      num_sweeps = 2  ! Usually 2 passes of 8-sweep are sufficient
 
      do isweep = 1, num_sweeps
      do sweep = 1, 8
 
        select case(sweep)
        case(1); i_start=0;     i_end=n1sub; i_step=1;  j_start=0;     j_end=n2sub; j_step=1;  k_start=0;     k_end=n3sub; k_step=1
        case(2); i_start=n1sub; i_end=0;     i_step=-1; j_start=0;     j_end=n2sub; j_step=1;  k_start=0;     k_end=n3sub; k_step=1
        case(3); i_start=0;     i_end=n1sub; i_step=1;  j_start=n2sub; j_end=0;     j_step=-1; k_start=0;     k_end=n3sub; k_step=1
        case(4); i_start=n1sub; i_end=0;     i_step=-1; j_start=n2sub; j_end=0;     j_step=-1; k_start=0;     k_end=n3sub; k_step=1
        case(5); i_start=0;     i_end=n1sub; i_step=1;  j_start=0;     j_end=n2sub; j_step=1;  k_start=n3sub; k_end=0;     k_step=-1
        case(6); i_start=n1sub; i_end=0;     i_step=-1; j_start=0;     j_end=n2sub; j_step=1;  k_start=n3sub; k_end=0;     k_step=-1
        case(7); i_start=0;     i_end=n1sub; i_step=1;  j_start=n2sub; j_end=0;     j_step=-1; k_start=n3sub; k_end=0;     k_step=-1
        case(8); i_start=n1sub; i_end=0;     i_step=-1; j_start=n2sub; j_end=0;     j_step=-1; k_start=n3sub; k_end=0;     k_step=-1
        end select
 
        do k = k_start, k_end, k_step
        do j = j_start, j_end, j_step
        do i = i_start, i_end, i_step
 
          ! Only process unset cells (band exterior)
          if (band_mask(i,j,k)) cycle
          if (abs(sgned(i,j,k)) < large_val) cycle
 
          ! Check 6-connected neighbors for a known sign
          ! X-neighbors
          if (i+1 <= n1sub) then
            neighbor_val = sgned(i+1,j,k)
            if (abs(neighbor_val) < large_val) then
              if (neighbor_val > 0.0d0) then
                sgned(i,j,k) = large_val
              else
                sgned(i,j,k) = -large_val
              endif
              cycle
            endif
          endif
          if (i-1 >= 0) then
            neighbor_val = sgned(i-1,j,k)
            if (abs(neighbor_val) < large_val) then
              if (neighbor_val > 0.0d0) then
                sgned(i,j,k) = large_val
              else
                sgned(i,j,k) = -large_val
              endif
              cycle
            endif
          endif
          ! Y-neighbors
          if (j+1 <= n2sub) then
            neighbor_val = sgned(i,j+1,k)
            if (abs(neighbor_val) < large_val) then
              if (neighbor_val > 0.0d0) then
                sgned(i,j,k) = large_val
              else
                sgned(i,j,k) = -large_val
              endif
              cycle
            endif
          endif
          if (j-1 >= 0) then
            neighbor_val = sgned(i,j-1,k)
            if (abs(neighbor_val) < large_val) then
              if (neighbor_val > 0.0d0) then
                sgned(i,j,k) = large_val
              else
                sgned(i,j,k) = -large_val
              endif
              cycle
            endif
          endif
          ! Z-neighbors
          if (k+1 <= n3sub) then
            neighbor_val = sgned(i,j,k+1)
            if (abs(neighbor_val) < large_val) then
              if (neighbor_val > 0.0d0) then
                sgned(i,j,k) = large_val
              else
                sgned(i,j,k) = -large_val
              endif
              cycle
            endif
          endif
          if (k-1 >= 0) then
            neighbor_val = sgned(i,j,k-1)
            if (abs(neighbor_val) < large_val) then
              if (neighbor_val > 0.0d0) then
                sgned(i,j,k) = large_val
              else
                sgned(i,j,k) = -large_val
              endif
              cycle
            endif
          endif
 
        enddo
        enddo
        enddo
 
      enddo ! sweep
      enddo ! isweep
 
    end subroutine propagate_sign_sweep
 
    !=========================================================================
    ! Phase 4: Ghost cell communication for SDF (3-direction MPI exchange)
    !=========================================================================
    subroutine ghostcell_communication_sdf()
      implicit none
 
      integer                                       ::  i, j, k
      integer                                       ::  buffsize
      integer                                       ::  row
      double precision, allocatable, dimension(:)   ::  sendeast, sendwest
      double precision, allocatable, dimension(:)   ::  recveast, recvwest
      integer                                       ::  request(4)
      integer                                       ::  ierr
 
      ! X-direction
      buffsize = (n2sub+1)*(n3sub+1)
      allocate(sendeast(0:buffsize-1), sendwest(0:buffsize-1))
      allocate(recveast(0:buffsize-1), recvwest(0:buffsize-1))
 
      do k = 0, n3sub
      do j = 0, n2sub
          row           = (n2sub+1)*k + j
          sendeast(row) = sgned_g(n1msub, j, k)
          sendwest(row) = sgned_g(1,      j, k)
      enddo
      enddo
 
      call MPI_Isend(sendeast, buffsize, MPI_DOUBLE_PRECISION, comm_1d_x1%east_rank, 111, comm_1d_x1%mpi_comm, request(1), ierr)
      call MPI_Irecv(recvwest, buffsize, MPI_DOUBLE_PRECISION, comm_1d_x1%west_rank, 111, comm_1d_x1%mpi_comm, request(2), ierr)
      call MPI_Isend(sendwest, buffsize, MPI_DOUBLE_PRECISION, comm_1d_x1%west_rank, 222, comm_1d_x1%mpi_comm, request(3), ierr)
      call MPI_Irecv(recveast, buffsize, MPI_DOUBLE_PRECISION, comm_1d_x1%east_rank, 222, comm_1d_x1%mpi_comm, request(4), ierr)
      call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
 
      do k = 0, n3sub
      do j = 0, n2sub
          row              = (n2sub+1)*k + j
          sgned_g(0,     j, k) = recvwest(row)
          sgned_g(n1sub, j, k) = recveast(row)
      enddo
      enddo
 
      deallocate(sendeast, sendwest, recveast, recvwest)
 
      ! Y-direction
      buffsize = (n1sub+1)*(n3sub+1)
      allocate(sendeast(0:buffsize-1), sendwest(0:buffsize-1))
      allocate(recveast(0:buffsize-1), recvwest(0:buffsize-1))
 
      do k = 0, n3sub
      do i = 0, n1sub
          row           = (n1sub+1)*k + i
          sendeast(row) = sgned_g(i, n2msub, k)
          sendwest(row) = sgned_g(i, 1,      k)
      enddo
      enddo
 
      call MPI_Isend(sendeast, buffsize, MPI_DOUBLE_PRECISION, comm_1d_x2%east_rank, 333, comm_1d_x2%mpi_comm, request(1), ierr)
      call MPI_Irecv(recvwest, buffsize, MPI_DOUBLE_PRECISION, comm_1d_x2%west_rank, 333, comm_1d_x2%mpi_comm, request(2), ierr)
      call MPI_Isend(sendwest, buffsize, MPI_DOUBLE_PRECISION, comm_1d_x2%west_rank, 444, comm_1d_x2%mpi_comm, request(3), ierr)
      call MPI_Irecv(recveast, buffsize, MPI_DOUBLE_PRECISION, comm_1d_x2%east_rank, 444, comm_1d_x2%mpi_comm, request(4), ierr)
      call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
 
      do k = 0, n3sub
      do i = 0, n1sub
          row              = (n1sub+1)*k + i
          sgned_g(i, 0,     k) = recvwest(row)
          sgned_g(i, n2sub, k) = recveast(row)
      enddo
      enddo
 
      deallocate(sendeast, sendwest, recveast, recvwest)
 
      ! Z-direction
      buffsize = (n1sub+1)*(n2sub+1)
      allocate(sendeast(0:buffsize-1), sendwest(0:buffsize-1))
      allocate(recveast(0:buffsize-1), recvwest(0:buffsize-1))
 
      do j = 0, n2sub
      do i = 0, n1sub
          row           = (n1sub+1)*j + i
          sendeast(row) = sgned_g(i, j, n3msub)
          sendwest(row) = sgned_g(i, j, 1     )
      enddo
      enddo
 
      call MPI_Isend(sendeast, buffsize, MPI_DOUBLE_PRECISION, comm_1d_x3%east_rank, 555, comm_1d_x3%mpi_comm, request(1), ierr)
      call MPI_Irecv(recvwest, buffsize, MPI_DOUBLE_PRECISION, comm_1d_x3%west_rank, 555, comm_1d_x3%mpi_comm, request(2), ierr)
      call MPI_Isend(sendwest, buffsize, MPI_DOUBLE_PRECISION, comm_1d_x3%west_rank, 666, comm_1d_x3%mpi_comm, request(3), ierr)
      call MPI_Irecv(recveast, buffsize, MPI_DOUBLE_PRECISION, comm_1d_x3%east_rank, 666, comm_1d_x3%mpi_comm, request(4), ierr)
      call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
 
      do j = 0, n2sub
      do i = 0, n1sub
          row              = (n1sub+1)*j + i
          sgned_g(i, j, 0     ) = recvwest(row)
          sgned_g(i, j, n3sub ) = recveast(row)
      enddo
      enddo
 
      deallocate(sendeast, sendwest, recveast, recvwest)
 
    end subroutine ghostcell_communication_sdf

    !=========================================================================
    ! Modified: get_signed_distance_function_global with narrowband optimization
    !=========================================================================
    subroutine get_signed_distance_function_global
      implicit none
 
      integer                                              ::  solidnum, i, k, j, ierr
      character(len=80)                                    ::  solidnum2str, filestl
      double precision, dimension(0:n1sub,0:n2sub,0:n3sub) ::  sgned
      logical, dimension(0:n1sub,0:n2sub,0:n3sub)          ::  band_mask
      double precision                                     ::  bandwidth, dx_max
      double precision                                     ::  t_band, t_sdf, t_sweep, t_tmp

      do solidnum = 1, nstlf

        write(solidnum2str,*) solidnum
        solidnum2str = adjustl(solidnum2str)
        solidnum2str = trim(solidnum2str)//'.dat'
        filestl= '../gensurf/output/Output'//trim(solidnum2str)

        call check_npt_total(filestl)

        if(myrank == 0) then
          write(*,*) '=============================================='
          write(*,*) 'ib parameters, #vert, #elem, #edge', para%nvert, para%nelem, para%nedge, 'in ', solidnum
        endif
        
        call var_surf_make
        call read_stl(filestl)
        
        ! Compute bandwidth from grid spacing
        dx_max = 0.0d0
        do i = 1, n1msub
          dx_max = max(dx_max, dx1_sub(i))
        enddo
        do j = 1, n2msub
          dx_max = max(dx_max, dx2_sub(j))
        enddo
        do k = 1, n3msub
          dx_max = max(dx_max, dx3_sub(k))
        enddo
        bandwidth = dble(N_band) * dx_max
 
        if(myrank == 0) then
          write(*,'(A,F12.6,A,F12.6)') '  [Narrowband] dx_max = ', dx_max, ', bandwidth = ', bandwidth
        endif
 
        ! Phase 1: Find narrowband cells
        t_tmp = MPI_Wtime()
        call find_narrowband_cells(band_mask, bandwidth)
        t_band = MPI_Wtime() - t_tmp

         if(myrank == 0) then
          write(*,'(A,F10.3,A)') '  [Timing] Narrowband detection: ', t_band, ' sec'
        endif

        ! Phase 2: Compute SDF in narrowband only
        t_tmp = MPI_Wtime()
        call compute_sdf_narrowband(sgned, band_mask)
        t_sdf = MPI_Wtime() - t_tmp

         if(myrank == 0) then
          write(*,'(A,F10.3,A)') '  [Timing] SDF computation:      ', t_sdf, ' sec'
        endif

        ! Phase 3: Propagate sign to exterior cells
        t_tmp = MPI_Wtime()
        call propagate_sign_sweep(sgned, band_mask)
        t_sweep = MPI_Wtime() - t_tmp
 
        if(myrank == 0) then 
          write(*,'(A,F10.3,A)') '  [Timing] Sign propagation:     ', t_sweep, ' sec'
        endif
 
        call var_surf_clean

        do k=0,n3sub
        do j=0,n2sub
        do i=0,n1sub
          sgned_g(i,j,k) = min(sgned(i,j,k),sgned_g(i,j,k))
        enddo
        enddo
        enddo

      enddo
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

    end subroutine get_signed_distance_function_global

    !=========================================================================
    ! Phase 5: Improved Heaviside function (SDF sign-based + smoothed option)
    !=========================================================================

    subroutine get_heaviside_function
      implicit none
      integer          :: i, j, k
      double precision :: epsilon_h, dx_local

      if(myrank==0) then
          write(*,*) ' Computing Heaviside Function...'
      endif
      heavi_g = 0.0d0
 
      if (use_smoothed_heaviside) then
        ! Smoothed Heaviside function for better interface representation
        do k = 1, n3msub
        do j = 1, n2msub
        do i = 1, n1msub
          ! Local grid spacing for smoothing width
          dx_local = max(dx1_sub(i), dx2_sub(j), dx3_sub(k))
          epsilon_h = 1.5d0 * dx_local
 
          if (sgned_g(i,j,k) > epsilon_h) then
            heavi_g(i,j,k) = 1.0d0
          elseif (sgned_g(i,j,k) < -epsilon_h) then
            heavi_g(i,j,k) = 0.0d0
          else
            ! Smooth transition using sine function
            heavi_g(i,j,k) = 0.5d0 * (1.0d0 + sgned_g(i,j,k)/epsilon_h &
                             + sin(PI * sgned_g(i,j,k)/epsilon_h) / PI)
          endif
        enddo
        enddo
        enddo
      else
        ! Sharp Heaviside based on SDF sign
        do k = 1, n3msub
        do j = 1, n2msub
        do i = 1, n1msub
          if (sgned_g(i,j,k) > 0.0d0) then
            heavi_g(i,j,k) = 1.0d0   ! Fluid
          else
            heavi_g(i,j,k) = 0.0d0   ! Solid
          endif
        enddo
        enddo
        enddo
      endif

    end subroutine get_heaviside_function

    !=========================================================================
    ! Main cell classification driver (modified call sequence)
    !=========================================================================

    subroutine cell_classification
      implicit none

      double precision :: t_total, t_sdf, t_ghost, t_heavi, t_tmp
      integer :: ierr
      t_total = MPI_Wtime()
 
      call var_cc_make
 
      ! SDF computation with narrowband optimization
      t_tmp = MPI_Wtime()
      call get_signed_distance_function_global
      t_sdf = MPI_Wtime() - t_tmp
 
      ! Phase 4: Ghost cell communication
      t_tmp = MPI_Wtime()
      call ghostcell_communication_sdf
      t_ghost = MPI_Wtime() - t_tmp
 
      ! Phase 5: Heaviside function
      t_tmp = MPI_Wtime()
      call get_heaviside_function
      t_heavi = MPI_Wtime() - t_tmp
 
      call save_cc
 
      if(myrank == 0) then
        write(*,*) '=============================================='
        write(*,'(A,F10.3,A)') ' [Total] SDF computation:      ', t_sdf, ' sec'
        write(*,'(A,F10.3,A)') ' [Total] Ghost cell comm:      ', t_ghost, ' sec'
        write(*,'(A,F10.3,A)') ' [Total] Heaviside function:   ', t_heavi, ' sec'
        write(*,'(A,F10.3,A)') ' [Total] Cell classification:  ', MPI_Wtime() - t_total, ' sec'
        write(*,*) '=============================================='
      endif
 
      call var_cc_clean
 
    end subroutine cell_classification
 
    !=========================================================================
    ! Original SDF computation (kept for reference/validation)
    !=========================================================================

    subroutine get_signed_distance_function_ANN(sgned)
      implicit none

      integer                                               ::  i, j, K, LL ,mm, id
      integer         , allocatable, dimension(:)           ::  indexarray, indexarray2, idmapping, idmapping2
      double precision, allocatable, dimension(:)           ::  distarray, distarray2
      double precision                                      ::  dist, distance, sign_SDF
      integer                                               ::  restriction
      double precision, dimension(3)                        ::  real_clst_pt, real_clst_norm, clst_pt, clst_norm
      double precision, dimension(3)                        ::  qp
      double precision, dimension(para%nvert)               ::  ptx, pty, ptz
      double precision                                      ::  check_dist
      
      double precision, dimension(0:n1sub,0:n2sub,0:n3sub)  ::  sgned

      allocate(indexarray(nsearch), idmapping(nsearch), distarray(nsearch))
      allocate(indexarray2(para%nvert), idmapping2(para%nvert), distarray2(para%nvert))

      ! create ANN-tree
      ptx(:) = vert(:)%xyz(1)
      pty(:) = vert(:)%xyz(2)
      ptz(:) = vert(:)%xyz(3)
      call wldst_create_tree(ptx(1), pty(1), ptz(1), para%nvert)

      do k = 0,n3sub
      do j = 0,n2sub
      do i = 0,n1sub
        distance = 999
        ! get distarray & index array from ANN-tree
        qp(1) = x1_sub(i)
        qp(2) = x2_sub(j)
        qp(3) = x3_sub(k)
        call wldst_getdists_xyz(qp(1), qp(2), qp(3), indexarray, distarray, nsearch)
        
        restriction = 0
        do LL = 1,nsearch
          if(distarray(LL).LT.distarray(1)) then
            write(*,*) 'ANN search error'
            stop
          else
            check_dist = (ptx(indexarray(1)+1)-ptx(indexarray(LL)+1))**2 &
            &           +(pty(indexarray(1)+1)-pty(indexarray(LL)+1))**2 &
            &           +(ptz(indexarray(1)+1)-ptz(indexarray(LL)+1))**2
            check_dist = sqrt(check_dist)
          endif
          if(check_dist.LT.2*dist_CAD) then
            restriction = restriction + 1
            idmapping(restriction) = LL
          endif
        enddo

        if(restriction.EQ.nsearch) then !extend search
          restriction = 0
          call wldst_getdists_xyz(qp(1), qp(2), qp(3), indexarray2, distarray2, para%nvert)
          do mm = 1,para%nvert
            if(distarray2(mm).LT.distarray2(1)) then
              write(*,*) 'ANN extend search error'
              stop
            else
              check_dist = (ptx(indexarray2(1)+1)-ptx(indexarray2(mm)+1))**2 &
              &           +(pty(indexarray2(1)+1)-pty(indexarray2(mm)+1))**2 &
              &           +(ptz(indexarray2(1)+1)-ptz(indexarray2(mm)+1))**2
              check_dist = sqrt(check_dist)
            endif
            if(check_dist.LT.2*dist_CAD) then
              restriction = restriction + 1
              idmapping2(restriction) = mm
            endif
          enddo
        endif

        ! find real closest point & distance
        if(restriction.LT.nsearch) then
          do LL = 1, restriction
            id = idmapping(LL)
            do mm = 1, vert(indexarray(id)+1)%share(1)
              call triangle_distance(qp, vert(indexarray(id)+1)%share(mm+1), dist, clst_pt, clst_norm)
              if(dist.LT.distance) then
                distance = dist
                real_clst_pt(:)   = clst_pt(:)
                real_clst_norm(:) = clst_norm(:)
              endif
            enddo
          enddo
        else
          do LL = 1, restriction
            id = idmapping2(LL)
            do mm = 1, vert(indexarray2(id)+1)%share(1)
              call triangle_distance(qp, vert(indexarray2(id)+1)%share(mm+1), dist, clst_pt, clst_norm)
              if(dist.LT.distance) then
                distance = dist
                real_clst_pt(:)   = clst_pt(:)
                real_clst_norm(:) = clst_norm(:)
              endif
            enddo
          enddo
        endif

        ! get_signed_distance_function
        sign_SDF     = (x1_sub(i)-real_clst_pt(1))*(real_clst_norm(1)) &
        &             +(x2_sub(j)-real_clst_pt(2))*(real_clst_norm(2)) &
        &             +(x3_sub(k)-real_clst_pt(3))*(real_clst_norm(3))
        sign_SDF = sign_SDF/abs(sign_SDF)

        sgned(i,j,k) = sign_SDF*distance
      enddo
      enddo
        ! if(myrank==0) then
        !     if(k==int(n3sub*0.1 )) write(*,*) '   10%..   '
        !     if(k==int(n3sub*0.25)) write(*,*) '   25%..   '
        !     if(k==int(n3sub*0.5 )) write(*,*) '   50%..   '
        !     if(k==int(n3sub*0.75)) write(*,*) '   75%..   '
        !     if(k==int(n3sub     )) write(*,*) '   100%..  '
        ! endif
      enddo

      call wldst_delete_tree
      deallocate(indexarray, indexarray2, distarray, distarray2)

    end subroutine get_signed_distance_function_ANN

    subroutine triangle_distance(qp, indx, distance, clst_pt, clst_norm)
      implicit none
      
      double precision, dimension(3)                  ::  qp
      integer                                         ::  indx

      double precision                                ::  distance
      double precision, dimension(3)                  ::  clst_pt, clst_norm

      ! Local variables
      double precision, dimension(3)                  ::  p1, p2, p3
      double precision                                ::  alpha, beta, gamma, delta, dist1, dist2, dist3
      double precision, parameter                     ::  dbl_epsilon = EPSILON(1.0d0)

      p1(:) = vert(elem(indx)%vert(1))%xyz(:)
      p2(:) = vert(elem(indx)%vert(2))%xyz(:)
      p3(:) = vert(elem(indx)%vert(3))%xyz(:)
      
      call get_ratio_of_area(p1, p2, p3, qp, elem(indx)%norm, alpha, beta, gamma)

      if(abs(alpha+beta+gamma-1.0d0).LT.dbl_epsilon) then
        distance = (qp(1)-elem(indx)%xyz(1))*elem(indx)%norm(1) &
        &         +(qp(2)-elem(indx)%xyz(2))*elem(indx)%norm(2) &
        &         +(qp(3)-elem(indx)%xyz(3))*elem(indx)%norm(3)
        distance = sqrt(distance)

        clst_pt(:)   = elem(indx)%xyz(:)
        clst_norm(:) = elem(indx)%norm(:)
      else
        call get_delta(p1, p2, qp, delta)
        if(delta.GT.0.0d0.AND.delta.LT.1.0d0) then
          call edge_dist(p1, p2, qp, distance)
          call edge_mdpt(elem(indx)%vert(1), elem(indx)%vert(2), elem(indx)%edge, clst_pt, clst_norm)
        else
          call get_delta(p2, p3, qp, delta)
          if(delta.GT.0.0d0.AND.delta.LT.1.0d0) then
            call edge_dist(p2, p3, qp, distance)
            call edge_mdpt(elem(indx)%vert(2), elem(indx)%vert(3), elem(indx)%edge, clst_pt, clst_norm)
          else
            call get_delta(p3, p1, qp, delta) 
            if(delta.GT.0.0d0.AND.delta.LT.1.0d0) then
              call edge_dist(p3, p1, qp, distance)
              call edge_mdpt(elem(indx)%vert(3), elem(indx)%vert(2), elem(indx)%edge, clst_pt, clst_norm)
            else
              dist1 = sqrt( (qp(1)-p1(1))**2 + (qp(2)-p1(2))**2 + (qp(3)-p1(3))**2 )
              dist2 = sqrt( (qp(1)-p2(1))**2 + (qp(2)-p2(2))**2 + (qp(3)-p2(3))**2 )
              dist3 = sqrt( (qp(1)-p3(1))**2 + (qp(2)-p3(2))**2 + (qp(3)-p3(3))**2 )
              distance = min(dist1, dist2, dist3)
              if(distance.EQ.dist1) then
                clst_pt(:)  = p1(:)
                clst_norm(:)= vert(elem(indx)%vert(1))%norm(:)
              elseif(distance.EQ.dist2) then
                clst_pt(:)  = p2(:)
                clst_norm(:)= vert(elem(indx)%vert(2))%norm(:)
              elseif(distance.EQ.dist3) then
                clst_pt(:)  = p3(:)
                clst_norm(:)= vert(elem(indx)%vert(3))%norm(:)
              endif
            endif
          endif
        endif
      endif
    end subroutine triangle_distance

    subroutine get_ratio_of_area(p1, p2, p3, qp, norm, alpha, beta, gamma)
      implicit none

      double precision, dimension(3)                  ::  p1, p2, p3, qp, norm
      double precision                                ::  alpha, beta, gamma

      ! Local variables
      double precision, dimension(3)                  ::  vec1, vec2, q
      double precision                                ::  cpn1, cpn2

      ! cal q
      q(1) = qp(1) - ( (norm(1)*(qp(1)-p1(1))) + (norm(2)*(qp(2)-p1(2))) + (norm(3)*(qp(3)-p1(3))) )*norm(1)
      q(2) = qp(2) - ( (norm(1)*(qp(1)-p1(1))) + (norm(2)*(qp(2)-p1(2))) + (norm(3)*(qp(3)-p1(3))) )*norm(2)
      q(3) = qp(3) - ( (norm(1)*(qp(1)-p1(1))) + (norm(2)*(qp(2)-p1(2))) + (norm(3)*(qp(3)-p1(3))) )*norm(3)

      ! cal alpha
      vec1(:) = p2(:) - p1(:)
      vec2(:) =  q(:) - p1(:)
      cpn1    = (vec1(2)*vec2(3)-vec1(3)*vec2(2))**2 &
      &        +(vec1(3)*vec2(1)-vec1(1)*vec2(3))**2 &
      &        +(vec1(1)*vec2(2)-vec1(2)*vec2(1))**2
      cpn1    = sqrt(cpn1)
      
      vec1(:) = p2(:) - p1(:)
      vec2(:) = p3(:) - p1(:)
      cpn2    = (vec1(2)*vec2(3)-vec1(3)*vec2(2))**2 &
      &        +(vec1(3)*vec2(1)-vec1(1)*vec2(3))**2 &
      &        +(vec1(1)*vec2(2)-vec1(2)*vec2(1))**2
      cpn2    = sqrt(cpn2)
      
      alpha = cpn1/cpn2

      ! cal beta
      vec1(:) = p2(:) - p3(:)
      vec2(:) =  q(:) - p2(:)
      cpn1    = (vec1(2)*vec2(3)-vec1(3)*vec2(2))**2 &
      &        +(vec1(3)*vec2(1)-vec1(1)*vec2(3))**2 &
      &        +(vec1(1)*vec2(2)-vec1(2)*vec2(1))**2
      cpn1    = sqrt(cpn1)
      
      vec1(:) = p2(:) - p1(:)
      vec2(:) = p3(:) - p1(:)
      cpn2    = (vec1(2)*vec2(3)-vec1(3)*vec2(2))**2 &
      &        +(vec1(3)*vec2(1)-vec1(1)*vec2(3))**2 &
      &        +(vec1(1)*vec2(2)-vec1(2)*vec2(1))**2
      cpn2    = sqrt(cpn2)

      beta = cpn1/cpn2

      ! cal gamma
      vec1(:) = p3(:) - p1(:)
      vec2(:) =  q(:) - p3(:)
      cpn1    = (vec1(2)*vec2(3)-vec1(3)*vec2(2))**2 &
      &        +(vec1(3)*vec2(1)-vec1(1)*vec2(3))**2 &
      &        +(vec1(1)*vec2(2)-vec1(2)*vec2(1))**2
      cpn1    = sqrt(cpn1)
      
      vec1(:) = p2(:) - p1(:)
      vec2(:) = p3(:) - p1(:)
      cpn2    = (vec1(2)*vec2(3)-vec1(3)*vec2(2))**2 &
      &        +(vec1(3)*vec2(1)-vec1(1)*vec2(3))**2 &
      &        +(vec1(1)*vec2(2)-vec1(2)*vec2(1))**2
      cpn2    = sqrt(cpn2)

      gamma = cpn1/cpn2

    end subroutine get_ratio_of_area

    subroutine get_delta(p1, p2, qp, delta)
      implicit none

      double precision, dimension(3)                      ::  p1, p2, qp
      double precision                                    ::  delta

      delta = (qp(1) - p1(1))*(p2(1) - p1(1)) + (qp(2) - p1(2))*(p2(2) - p1(2)) + (qp(3) - p1(3))*(p2(3) - p1(3))
      delta = delta/sqrt( (p2(1)-p1(1))**2 + (p2(2)-p1(2))**2 + (p2(3)-p1(3))**2 )

    end subroutine get_delta

    subroutine edge_dist(p1, p2, qp, distance)
      implicit none

      double precision, dimension(3)                        ::  p1, p2, qp
      double precision                                      ::  distance

      ! Local variables
      double precision, dimension(3)                        ::  vec1, vec2
      double precision                                      ::  cpn

      vec1(:) = p2(:) - p1(:)
      vec2(:) = p1(:) - qp(:)
      cpn     = (vec1(2)*vec2(3)-vec1(3)*vec2(2))**2 &
      &        +(vec1(3)*vec2(1)-vec1(1)*vec2(3))**2 &
      &        +(vec1(1)*vec2(2)-vec1(2)*vec2(1))**2
      cpn     = sqrt(cpn)

      distance = cpn/sqrt( (p2(1)-p1(1))**2 + (p2(2)-p1(2))**2 + (p2(3)-p1(3))**2)

    end subroutine edge_dist

    subroutine edge_mdpt(p1, p2, indx, clst_pt, clst_norm)
      implicit none

      integer                                     ::  p1, p2
      integer         , dimension(3)              ::  indx
      double precision, dimension(3)              ::  clst_pt, clst_norm

      ! Local variables
      integer                                     ::  ev1, ev2
      integer                                     ::  check

      ev1 = edge(indx(1))%vert(1)
      ev2 = edge(indx(1))%vert(2)

      check = 0
      if(p1.EQ.ev1.AND.p2.EQ.ev2) then
        clst_pt(:)    = edge(indx(1))%xyz(:)
        clst_norm(:)  = edge(indx(1))%norm(:)
        check = check +1
      elseif(p1.EQ.ev2.AND.p2.EQ.ev1) then
        clst_pt(:)    = edge(indx(1))%xyz(:)
        clst_norm(:)  = edge(indx(1))%norm(:)
        check = check +1
      endif

      ev1 = edge(indx(2))%vert(1)
      ev2 = edge(indx(2))%vert(2)

      if(p1.EQ.ev1.AND.p2.EQ.ev2) then
        clst_pt(:)    = edge(indx(2))%xyz(:)
        clst_norm(:)  = edge(indx(2))%norm(:)
        check = check +1
      elseif(p1.EQ.ev2.AND.p2.EQ.ev1) then
        clst_pt(:)    = edge(indx(2))%xyz(:)
        clst_norm(:)  = edge(indx(2))%norm(:)
        check = check +1
      endif

      ev1 = edge(indx(3))%vert(1)
      ev2 = edge(indx(3))%vert(2)

      if(p1.EQ.ev1.AND.p2.EQ.ev2) then
        clst_pt(:)    = edge(indx(3))%xyz(:)
        clst_norm(:)  = edge(indx(3))%norm(:)
        check = check +1
      elseif(p1.EQ.ev2.AND.p2.EQ.ev1) then
        clst_pt(:)    = edge(indx(3))%xyz(:)
        clst_norm(:)  = edge(indx(3))%norm(:)
        check = check +1
      endif

      if(check.NE.1) then
        write(*,*) 'edge_mdpt error, check =',check
      endif

    end subroutine edge_mdpt
    subroutine save_cc
      implicit none
      integer             ::  i, j, K
      character(len=100)  ::  fileout
 
      if(myrank==0) then
          write(*,*) ' Writing outputs..'
      endif
 
      write(fileout, '(1A33,3(A1,I1),1A4)') 'output/result/SDF/out_signed_dist','_',&
                                            comm_1d_x1%myrank,'_',comm_1d_x2%myrank,'_',comm_1d_x3%myrank,'.dat'
 
      open(4,file=trim(fileout),action='write')
          write(4,*) 'variables="X","Y","Z","Phi"'
          write(4,*) 'zone i=',n1msub,', j=',n2msub,', k=',n3msub
          do k=1,n3msub
          do j=1,n2msub
          do i=1,n1msub
              write(4,*) x1_sub(i),x2_sub(j),x3_sub(k),sgned_g(i,j,k)
          enddo
          enddo
          enddo
      close(4)
 
      write(fileout, '(1A25,3(A1,I1),1A4)') 'output/result/G/out_heavi','_',&
                                            comm_1d_x1%myrank,'_',comm_1d_x2%myrank,'_',comm_1d_x3%myrank,'.dat'
 
      open(5,file=trim(fileout),action='write')
          write(5,*) 'variables="X","Y","Z","G"'
          write(5,*) 'zone i=',n1msub,', j=',n2msub,', k=',n3msub
          do k=1,n3msub
          do j=1,n2msub
          do i=1,n1msub
              write(5,*) x1_sub(i),x2_sub(j),x3_sub(k),heavi_g(i,j,k)
          enddo
          enddo
          enddo
      close(5)
 
      if(myrank==0) then
          write(*,*) '=============================================='
          write(*,*) ' Cell Classification is Done.'
          write(*,*) '=============================================='
      endif
 
    end subroutine save_cc
 
    subroutine consensus_algorithm(sgned)
      implicit none

      double precision, dimension(0:n1sub,0:n2sub,0:n3sub)  ::  sgned

      ! Local variables
      integer                                               ::  No, Ni
      integer                                               ::  Lx, Ly, Lz, i, j ,k
      integer                                               ::  corr_sign, nc, check_corr
      integer                                               ::  count = 0

      do 
        nc = 0

        do k = 1,n3msub
        do j = 1,n2msub
        do i = 1,n1msub

          do Lz = -1,1
          do Ly = -1,1
          do Lx = -1,1
            
            No = 0
            Ni = 0

            if(sgned(i+Lx,j+Ly,k+Lz).LT.0.0d0) then
              Ni = Ni + 1
            else
              No = No + 1
            endif

          enddo
          enddo
          enddo

          corr_sign = (No - Ni)/abs(No - Ni)
          check_corr = corr_sign*sgned(i,j,k)/abs(sgned(i,j,k))
          sgned(i,j,k) = corr_sign*abs(sgned(i,j,k))

          if(check_corr.LT.0) nc = nc + 1

        enddo
        enddo
        enddo
      
        if(nc.EQ.0) then
          if(myrank==0) then
            print*, 'iter=', count, 'nc=',nc
          endif
          exit
        else
          count = count + 1
          if(myrank==0) then
            print*, 'iter=', count, 'nc=',nc
          endif
        endif
          
      enddo

      if(myrank==0) then
        write(*,*) 'end correct sign'
      endif

    end subroutine consensus_algorithm

    subroutine debug_read(sgned)
      implicit none
      integer             ::  i, j, K
      character(len=100)  ::  fileout, dummy
      double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub)::sgned

      write(fileout, '(1A33,3(A1,I1),1A4)') 'output/debug/out_signed_dist','_',&
                                            comm_1d_x1%myrank,'_',comm_1d_x2%myrank,'_',comm_1d_x3%myrank,'.dat'

      open(4,file=trim(fileout),action='read')
          read(4,*)
          read(4,*) 
          do k=0,n3sub
          do j=0,n2sub
          do i=0,n1sub
              read(4,*) x1_sub(i),x2_sub(j),x3_sub(k),sgned(i,j,k)
          enddo
          enddo
          enddo
      close(4)      
    end subroutine debug_read

    subroutine debug_save_cc(sgned)
      implicit none
      integer             ::  i, j, K
      character(len=100)  ::  fileout
      double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub)::sgned

      write(fileout, '(1A33,3(A1,I1),1A4)') 'output/debug/out_signed_dist','_',&
                                            comm_1d_x1%myrank,'_',comm_1d_x2%myrank,'_',comm_1d_x3%myrank,'.dat'

      open(4,file=trim(fileout),action='write')
          write(4,*) 'variables="X","Y","Z","Phi"'
          write(4,*) 'zone i=',n1sub,', j=',n2sub,', k=',n3sub
          do k=0,n3sub
          do j=0,n2sub
          do i=0,n1sub
              write(4,*) x1_sub(i),x2_sub(j),x3_sub(k),sgned(i,j,k)
          enddo
          enddo
          enddo
      close(4)

    end subroutine debug_save_cc
end module mpi_cellclass