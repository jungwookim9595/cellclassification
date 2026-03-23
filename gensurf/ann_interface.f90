!---------------------------------------------------------------------------------------
!
!     Module ann_interface
!
!     iso_c_binding interface for ANN kd-tree C++ functions in walldist.cpp
!
!     walldist.cpp functions are wrapped with extern "C" so no name mangling occurs.
!     The function names in walldist.cpp already include trailing underscores
!     (e.g. wldst_create_tree_), so bind(C, name="...") matches them directly.
!
!     All arguments in walldist.cpp are passed by pointer (C convention for
!     Fortran interop), which matches Fortran's default pass-by-reference.
!     Therefore NO 'value' attribute is used.
!
!     Functions available in walldist.cpp:
!       wldst_init_          : Initialize tree from file
!       wldst_create_tree_   : Create tree from x,y,z arrays
!       wldst_delete_tree_   : Delete tree
!       wldst_getdist_       : Get min distance (query_pt array, dist)
!       wldst_getdist_xyz_   : Get min distance + nearest index (x,y,z,indx,dist)
!       wldst_getdist_xyzs_  : Get signed distance (x,y,z,normal,nsize,indx,dist,Sign)
!       wldst_getdists_      : Batch distance for 3D array (x,y,z,d,ii,jj,kk)
!       wldst_getdists_xyz_  : k-nearest search (x,y,z,indx,dist,knear)
!
!     Created: 2026-03-23
!     Updated: 2026-03-23 — removed nonexistent wldst_ann_search,
!              added all actual walldist.cpp functions
!
!---------------------------------------------------------------------------------------
module ann_interface
  use iso_c_binding
  implicit none

  interface

     !------------------------------------------------------------------------
     ! wldst_init_ : Initialize ANN tree from surface points file
     !------------------------------------------------------------------------
     subroutine wldst_init(filename, length) bind(C, name="wldst_init_")
       use iso_c_binding
       implicit none
       character(kind=c_char), intent(in) :: filename(*)
       integer(c_int), intent(in) :: length
     end subroutine wldst_init

     !------------------------------------------------------------------------
     ! wldst_create_tree_ : Create kd-tree from separate x, y, z arrays
     !   px, py, pz : coordinate arrays (length = *n)
     !   n          : number of points (passed by pointer)
     !------------------------------------------------------------------------
     subroutine wldst_create_tree(px, py, pz, n) bind(C, name="wldst_create_tree_")
       use iso_c_binding
       implicit none
       real(c_double), intent(in) :: px(*)
       real(c_double), intent(in) :: py(*)
       real(c_double), intent(in) :: pz(*)
       integer(c_int), intent(in) :: n
     end subroutine wldst_create_tree

     !------------------------------------------------------------------------
     ! wldst_delete_tree_ : Delete kd-tree and free memory
     !------------------------------------------------------------------------
     subroutine wldst_delete_tree() bind(C, name="wldst_delete_tree_")
       use iso_c_binding
       implicit none
     end subroutine wldst_delete_tree

     !------------------------------------------------------------------------
     ! wldst_getdist_ : Get minimum distance to nearest surface point
     !   query_pt : 3-element double array (query coordinates)
     !   dist     : output minimum distance (scalar)
     !------------------------------------------------------------------------
     subroutine wldst_getdist(query_pt, dist) bind(C, name="wldst_getdist_")
       use iso_c_binding
       implicit none
       real(c_double), intent(in) :: query_pt(*)
       real(c_double), intent(out) :: dist
     end subroutine wldst_getdist

     !------------------------------------------------------------------------
     ! wldst_getdist_xyz_ : Get nearest distance and index from x,y,z scalars
     !   x, y, z  : query point coordinates (passed by pointer)
     !   indx     : output index of nearest point (0-based C index)
     !   dist     : output minimum distance
     !------------------------------------------------------------------------
     subroutine wldst_getdist_xyz(x, y, z, indx, dist) bind(C, name="wldst_getdist_xyz_")
       use iso_c_binding
       implicit none
       real(c_double), intent(in) :: x, y, z
       integer(c_int), intent(out) :: indx
       real(c_double), intent(out) :: dist
     end subroutine wldst_getdist_xyz

     !------------------------------------------------------------------------
     ! wldst_getdist_xyzs_ : Get signed distance using normal vectors
     !   x, y, z  : query point coordinates
     !   normal   : normal vector array (flattened)
     !   nsize    : number of surface points
     !   indx     : output nearest index (0-based)
     !   dist     : output minimum distance
     !   sgn      : output sign (+1.0 or -1.0, inside/outside)
     !------------------------------------------------------------------------
     subroutine wldst_getdist_xyzs(x, y, z, normal, nsize, indx, dist, sgn) &
          bind(C, name="wldst_getdist_xyzs_")
       use iso_c_binding
       implicit none
       real(c_double), intent(in) :: x, y, z
       real(c_double), intent(in) :: normal(*)
       integer(c_int), intent(in) :: nsize
       integer(c_int), intent(out) :: indx
       real(c_double), intent(out) :: dist
       real(c_double), intent(out) :: sgn
     end subroutine wldst_getdist_xyzs

     !------------------------------------------------------------------------
     ! wldst_getdists_ : Batch distance computation for 3D Fortran array
     !   x, y, z    : 3D coordinate arrays (with ghost cells, Fortran layout)
     !   d          : output 3D distance array
     !   ii, jj, kk : array dimensions
     !------------------------------------------------------------------------
     subroutine wldst_getdists(x, y, z, d, ii, jj, kk) bind(C, name="wldst_getdists_")
       use iso_c_binding
       implicit none
       real(c_double), intent(in) :: x(*), y(*), z(*)
       real(c_double), intent(out) :: d(*)
       integer(c_int), intent(in) :: ii, jj, kk
     end subroutine wldst_getdists

     !------------------------------------------------------------------------
     ! wldst_getdists_xyz_ : k-nearest neighbor search
     !   x, y, z  : query point coordinates (passed by pointer)
     !   indx     : output array of k nearest neighbor indices (0-based)
     !   dist     : output array of k nearest distances (Euclidean)
     !   knear    : number of nearest neighbors to find (passed by pointer)
     !
     !   This is the function used by find_duplications for vertex
     !   deduplication via k-nearest neighbor search.
     !------------------------------------------------------------------------
     subroutine wldst_getdists_xyz(x, y, z, indx, dist, knear) &
          bind(C, name="wldst_getdists_xyz_")
       use iso_c_binding
       implicit none
       real(c_double), intent(in) :: x, y, z
       integer(c_int), intent(out) :: indx(*)
       real(c_double), intent(out) :: dist(*)
       integer(c_int), intent(in) :: knear
     end subroutine wldst_getdists_xyz

  end interface

end module ann_interface
