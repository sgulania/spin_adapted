!--------------------------------------------------------------------------------
! Hard-coded subroutine for 3 electron and S=1/2
! Generate the wavefunction and call subroutine
! for computing the Hamiltonian matrix. In this case
! the wave function is given by
!   Psi = phi_1 X theta_1 + phi_2 X theta_2 , 
! where theta_1 and theta_2 are orthonormal. 
!--------------------------------------------------------------------------------

subroutine young_frame_S3
  USE variables
  implicit none

  integer i,j
  integer, dimension(:,:), allocatable  :: weyl,phi_1,phi_2
  real*8,  dimension(:), allocatable  :: cof_phi_1,cof_phi_2
  integer, dimension(:), allocatable  ::  list_phi_1,list_phi_2 
  integer n_phi_1,n_phi_2

!--------------------------------------------------------------------------------
! weyl - contains ordered elements of each weyl tabuleau
! cof_phi_1 - set of coefficients for each permuatation in phi_1
! cof_phi_2 - set of coefficients for each permuatation in phi_2
! list_phi_1 - array list for phi_1
! list_phi_2 - array list for phi_2
! n_phi_1 - number of elements in phi_1
! n_phi_2 - number of elements in phi_2
!--------------------------------------------------------------------------------

  open(unit = 104, file = 'weyl.dat', status = 'old', action = 'read')
  read(104,*) n_weyl  

   allocate ( weyl(n_weyl,n_el),list_phi_1(n_weyl+1),list_phi_2(n_weyl+1), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

   n_phi_1 = n_weyl*6
   n_phi_2 = n_weyl*4

   allocate (phi_1(n_phi_1,n_el),phi_2(n_phi_2,n_el), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
 
   allocate (cof_phi_1(n_phi_1),cof_phi_2(n_phi_2), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"


   do i = 1, n_weyl
    read(104,*) weyl(i,1),weyl(i,2),weyl(i,3)   
   end do
   close(104)

  list_phi_1(1) =0
  list_phi_2(1) =0

!--------------------------------------------------------------------------------
! Generating phi_1
do i=1,n_weyl
j=((i-1)*6)
j=j+1
cof_phi_1(j) = 1.d0; phi_1(j,1)= weyl(i,1);phi_1(j,2)= weyl(i,2);phi_1(j,3)=weyl(i,3) ! e
j=j+1
cof_phi_1(j) =-1.d0; phi_1(j,1)= weyl(i,2);phi_1(j,2)= weyl(i,1);phi_1(j,3)=weyl(i,3) ! 12
j=j+1
cof_phi_1(j) = 0.5d0; phi_1(j,1)= weyl(i,3);phi_1(j,2)= weyl(i,2);phi_1(j,3)=weyl(i,1) ! 13
j=j+1
cof_phi_1(j) = 0.5d0; phi_1(j,1)= weyl(i,1);phi_1(j,2)= weyl(i,3);phi_1(j,3)=weyl(i,2) ! 23
j=j+1
cof_phi_1(j) =-0.5d0; phi_1(j,1)= weyl(i,3);phi_1(j,2)= weyl(i,1);phi_1(j,3)=weyl(i,2) ! 123
j=j+1
cof_phi_1(j) =-0.5d0; phi_1(j,1)= weyl(i,2);phi_1(j,2)= weyl(i,3);phi_1(j,3)=weyl(i,1) ! 132
list_phi_1(i+1) =j
end do

!--------------------------------------------------------------------------------
! Generating phi_2
do i=1,n_weyl
j=((i-1)*4)
j=j+1
cof_phi_2(j) = sqrt(3.d0)/2.d0; phi_2(j,1)= weyl(i,3);phi_2(j,2)= weyl(i,2);phi_2(j,3)=weyl(i,1) ! 13
j=j+1
cof_phi_2(j) =-sqrt(3.d0)/2.d0; phi_2(j,1)= weyl(i,1);phi_2(j,2)= weyl(i,3);phi_2(j,3)=weyl(i,2) ! 23
j=j+1
cof_phi_2(j) =-sqrt(3.d0)/2.d0; phi_2(j,1)= weyl(i,3);phi_2(j,2)= weyl(i,1);phi_2(j,3)=weyl(i,2) ! 123
j=j+1
cof_phi_2(j) = sqrt(3.d0)/2.d0; phi_2(j,1)= weyl(i,2);phi_2(j,2)= weyl(i,3);phi_2(j,3)=weyl(i,1) ! 132
list_phi_2(i+1) =j
end do

!write(6,*) "Pass Young"
call hamiltonian(cof_phi_1,cof_phi_2,phi_1,phi_2,list_phi_1,list_phi_2,n_phi_1,n_phi_2)

end subroutine young_frame_S3
