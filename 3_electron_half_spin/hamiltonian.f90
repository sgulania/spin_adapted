!--------------------------------------------------------------------------------
! Subroutine for evaluating Hamiltonian matrix
! Hardcoded for N=3, S=1/2 but can be easily extened
! fof general case
! Input - Wavefunction 
! Output - Spectrum
!--------------------------------------------------------------------------------

subroutine hamiltonian(cof_phi_1,cof_phi_2,phi_1,phi_2,list_phi_1,list_phi_2,n_phi_1,n_phi_2)
USE variables
implicit none


integer i,j
!--------------------------------------------------------------------------------
! Wavefunction information

integer n_phi_1,n_phi_2
integer list_phi_1(n_weyl+1),list_phi_2(n_weyl+1)
integer phi_1(n_phi_1,n_el), phi_2(n_phi_2,n_el)
real*8 cof_phi_1(n_phi_1),cof_phi_2(n_phi_2)
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
! Variable for orthogonalization and digonalization
real*8 H_f(n_weyl,n_weyl), S_f(n_weyl,n_weyl), W_f(n_weyl)
!--------------------------------------------------------------------------------

integer i11,i11n,j11,j11n,i21,i21n,j21,j21n  !dummy variables
real*8 val,over                              !storing hamiltonian value and overlap

!--------------------------------------------------------------------------------
! Using the list of phi_1 and phi_2 to evaluate 
! Hamiltonian element 
!--------------------------------------------------------------------------------

do i=1,n_weyl
    do j=1,n_weyl

       i11  = list_phi_1(i)+1
       i11n = list_phi_1(i+1)

       j11  = list_phi_1(j)+1
       j11n = list_phi_1(j+1)
      
       i21  = list_phi_2(i)+1
       i21n = list_phi_2(i+1)
       
       j21  = list_phi_2(j)+1
       j21n = list_phi_2(j+1)
       
       
       call one_element(phi_1(i11:i11n,:),phi_2(i21:i21n,:),&
                        phi_1(j11:j11n,:),phi_2(j21:j21n,:),&
                    cof_phi_1(i11:i11n),cof_phi_2(i21:i21n),&
                    cof_phi_1(j11:j11n),cof_phi_2(j21:j21n),&
                    i11n-i11+1,i21n-i21+1,j11n-j11+1,j21n-j21+1,val,over)

      S_f(i,j)=over
      if(abs(over).lt.1.D-8) then
        S_f(i,j)=0.d0
      end if

      H_f(i,j)= val

      if(abs(H_f(i,j)).lt.1.D-8) then
       H_f(i,j)=0.d0
      end if
    
    end do
end do

!open (unit = 23, file = "Hamiltonian_Matrix.out") 
! Printing Hamiltonian matrix for 3-electron wavefunction obtained using                      
! Young frame formalism 
!do i=1,n_weyl
! write(23,"(9f13.8)")H_f(i,:)
!enddo


open (unit = 24, file = "Overlap.out")

! Printing overlap matrix for 3-electron wavefunction obtained using
! Young frame formalism
do i=1,n_weyl
 write(24,"(9f13.8)")S_f(i,:)
enddo

! Orthogonalizing the basis and doing Hamiltonian transformation
! accordingly 

call ortho_transformation(S_f,H_f,n_weyl)

open (unit = 25, file = "Transformed_Hamiltonian_Matrix.out") 

!Digonalizing the transformed Hamiltonian matrix
do i=1,n_weyl                                                                             
 write(25,"(9f13.8)")H_f(i,:)                                                             
enddo  


call diagonalization(H_f,n_weyl,W_f)

 write(6,*) "Spectrum for Sz=1/2"
do i=1,n_weyl
 write(6,"(I5,9f13.8)")i,W_f(i)+nuc_repul
!10 format (4I5,9f13.8)
enddo


end subroutine hamiltonian
