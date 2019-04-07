!--------------------------------------------------------------------------------
! Subroutine for evaluating one element of Hamiltonian Matrix
!--------------------------------------------------------------------------------
subroutine one_element(p1,p2,s1,s2,cp1,cp2,cs1,cs2,np1,np2,ns1,ns2,val,over)
  USE variables
  implicit none
  
  integer i,j,k
  integer np1,np2,ns1,ns2
  real*8 cp1(np1),cp2(np2),cs1(ns1),cs2(ns2),val
  integer p1(np1,n_el),p2(np2,n_el),s1(ns1,n_el),s2(ns2,n_el)
  real*8 v,val1,val2,nm,norm_aa1,norm_aa2,norm_bb1,norm_bb2,over

!--------------------------------------------------------------------------------
!Hamiltonian value for spatial wavefunction assuming 
!p and s are two spatial parts for orthogonal spin eignefunction

!p1 \times theta1 + s1 \times theta2 -> bra  
!p2 \times theta1 + s2 \times theta2  -> ket
! <theta1|theta2> = 0 -> orthogonal
!np1 - number of element in p1
!np2 - number of element in p2
!ns1 - number of element in s1                                                            
!ns2 - number of element in s2
!cp1 - coefficieng of p1                                                            
!cp2 - coefficieng of p2                                                            
!cs1 - coefficieng of s1                                                           
!cs2 - coefficieng of s2 
!val = hamiltonian = <bra|H|ket>
!over = overlap = <brak|ket>
!--------------------------------------------------------------------------------


!write(6,*)s1(ns1,:)
!write(6,*)np1,ns1,np2,ns2

! do i=1,np1  
!  write(6,*)cp1(i),p1(i,:)
! end do

 norm_aa1=0.d0;
 norm_aa2=0.d0;

 norm_bb1=0.d0;
 norm_bb2=0.d0;

! Computing norm square of phi_11(i)

 do i=1,np1
  do k=1,np1
!  write(6,*)p1(i,:),p1(k,:)
!  write(6,*)cp1(i),cp1(k)
  call norm1(p1(i,:),p1(k,:),nm)
  norm_aa1=norm_aa1+cp1(i)*cp1(k)*nm
  end do
 end do

!write(6,*) norm_aa1

! Computing norm square of phi_11(j)



 do i=1,ns1
  do k=1,ns1
  call norm1(s1(i,:),s1(k,:),nm)
  norm_aa2=norm_aa2+cs1(i)*cs1(k)*nm
  end do
 end do

! do i=1,ns1
! write(6,*)cs1(i),s1(i,:)
! end do
! write(6,*) norm_aa2

!stop
! Computing the <phi_11(i)|H|phi_11(j)>

val1=0
 do i=1,np1
   do j=1,ns1
       call ham(p1(i,:),s1(j,:),v)
       val1=val1+cp1(i)*cs1(j)*v
   end do
 end do


! Computing norm square of phi_12(i)

 do i=1,np2
  do k=1,np2
   call norm1(p2(i,:),p2(k,:),nm)
   norm_bb1=norm_bb1+cp2(i)*cp2(k)*nm
  end do
 end do

! Computing norm square of phi_12(j)

 do i=1,ns2
  do k=1,ns2
  call norm1(s2(i,:),s2(k,:),nm)
  norm_bb2=norm_bb2+cs2(i)*cs2(k)*nm
  end do
 end do

! Computing <phi_12(i)|H|phi_12(j)>

val2=0;
  do i=1,np2
   do j=1,ns2
       call ham(p2(i,:),s2(j,:),v)
       val2=val2+cp2(i)*cs2(j)*v
   end do
 end do

 over=0.d0
 ! Computing overlap between psi_i and psi_j
    do i=1,np1
      do k=1,ns1
         call norm1(p1(i,:),s1(k,:),nm)
         over=over+cp1(i)*cs1(k)*nm
!         write(6,*) over         
      end do
    end do
!     over = 0.d0

    do i=1,np2
      do k=1,ns2
         call norm1(p2(i,:),s2(k,:),nm)
         over=over+cp2(i)*cs2(k)*nm
!         write(6,*) over
      end do
    end do
! write(6,*) i,j,over

! This is more crucial step 

!
!1. This assuming the theta1 and theta2 are already orthonormal and normalizing it with norms
!   of individual part
!
   val= (val1+val2)
!   (sqrt((norm_aa1+norm_bb1))*sqrt((norm_aa2+norm_bb2)))
!   write (6,*) val,val1,val2!,norm_aa1,norm_bb1,norm_aa2,norm_bb2,"val"

!2. This assuming the alpha1 and alpha2 are already orthonormal and using the direct expression
!   the book by Ruben Paunz
!
! val=(val1+val2)/(6.d0)

!3. This assuming the alpha1 and alpha2 are not normalized and get constructed from conjugate young
!   tabulae
!   theta1= sqrt(3)[beta*alpha*alpha - alpha*alpha*beta]
!   theta2= [2*alpha*alpha*beta - beta*lpha*alpha - alpha*beta*alpha ]
!   <alpha1|alpha1> = 6 ; <alpha2|alpha2>=6
!
!  val= (val1*6.d0+val2*3.d0)/(sqrt((norm_aa1*6.d0+norm_bb1*3.d0))*sqrt((norm_aa2*6.d0+norm_bb2*3.d0)))
 
end subroutine one_element
