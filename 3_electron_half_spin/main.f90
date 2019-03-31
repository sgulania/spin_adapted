!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Spin Adapted using Young Frame 
! Authors - Sahil Gulania and James Daniel Whitfield
! University of Southern California and Dartmouth College
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Work flow 
! 1. Read the system information - 
!    No. electron (n_el), Spin of the system (sp_deg), Nuclear
!    energy (nuc_repul), number of basis (n_bas)
! 2. Read 1e and 2e atomic integrals 
! 3. Generate young frame dependent on n_el, sp_deg  
! 4. Using young frame generate irreducible representation
! 5. Use irreps to generate wavefunction
! 6. Evaluate Hamiltonian matrix and perform diagonalization
! 7. Output - Spectrum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Present code use irrpes directly and generate wavefunction
! Using the wavefunction evaluate Hamiltonian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
USE variables

implicit none 

call read_write

! call generate_irrep

! call generate_wavefunction

! call hamiltonian

call young_frame_S3

end program main

