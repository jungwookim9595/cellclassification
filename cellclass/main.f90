!!======================================================================================================================
!> @file        main.f90
!> @brief       This file contains the main subroutines for an example problem of PaScaL_TCS.
!> @details     The target example problem is a three-dimensional(3D) Rayleigh-Bernard natural convection problem
!>              in a channel domain with the boundary conditions of vertically different wall temperature and 
!>              horizontally periodic boundaries. The non-Oberbeck–Boussinesq effect is considered. 
!> @author      
!>              - Kiha Kim (k-kiha@yonsei.ac.kr), Department of Computational Science & Engineering, Yonsei University
!>              - Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
!>              - Jung-Il Choi (jic@yonsei.ac.kr), Department of Computational Science & Engineering, Yonsei University
!>
!> @date        October 2022
!> @version     1.0
!> @par         Copyright
!>              Copyright (c) 2022 Kiha Kim and Jung-Il choi, Yonsei University and 
!>              Ji-Hoon Kang, Korea Institute of Science and Technology Information, All rights reserved.
!> @par         License     
!>              This project is release under the terms of the MIT License (see LICENSE in )
!======================================================================================================================

!>
!> @brief       Main execution program for the example problem of PaScaL_TCS.
!>

program main
  
    use mpi
    use mpi_topology
    use global
    use mpi_subdomain
    use mpi_cellclass

    implicit none

    integer :: ierr
    
    call MPI_Init(ierr)
    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)

    if(myrank==0) write(*,*) '[Main] cell classification start'
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
    call global_inputpara()
    if(myrank==0) write(*,*) '[Main] Read input parameters!'

    call mpi_topology_make()
    call mpi_subdomain_make()
    call mpi_subdomain_mesh()
    call mpi_subdomain_indices()

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call cell_classification    

    call mpi_subdomain_indices_clean()
    call mpi_subdomain_clean()

    call mpi_topology_clean()
    call MPI_FINALIZE(ierr)

    if(myrank==0) write(*,*) '[Main] cell classification complete! '
end
