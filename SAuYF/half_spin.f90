!Authors - James Daniel Whitfield - Dartmouth College - http://physics.dartmouth.edu/people/james-daniel-whitfield
!        - Sahil Gulania - University of Souther California 
program half_spin
!implicit double precision (A-Z)
integer i,j,k,l,a,b,c,d,ii
integer aa(9,3),a1(81),b1(81),c1(81),d1(81)
integer p11(54,4)
integer p12(36,4)
real*8  Ov(3,3), HOv(3,3),e1(81),S1(3,3)
real*8  EEOv(3,3,3,3)
!common  S(3,3),H(3,3),EE(3,3,3,3)
real*8  S(3,3),H(3,3),EE(3,3,3,3)
real*8  W(3),WORK(8),W_f(9),WORK_f(26),X1(3,3),W1(3),WORK1(8)
integer INFO, LWORK,LWORK_f,LWORK1
real*8  val,H_f(9,9)
real*8  Xsymm(3,3),H1(3,3)
open(unit = 100, file = 'weyl.dat', status = 'old', action = 'read')
 do i = 1,  9
 read(100,*) aa(i,1),aa(i,2),aa(i,3)
end do

! do i = 1,8
! write(6,*) a(i,1),a(i,2),a(i,3)
! end do
!

!------------------------------------------------------------------
! Generating all the phi_11
do i=1,9
j=((i-1)*6)
j=j+1
p11(j,1) = 2; p11(j,2)= aa(i,1);p11(j,3)= aa(i,2);p11(j,4)=aa(i,3)
j=j+1
p11(j,1) =-2; p11(j,2)= aa(i,2);p11(j,3)= aa(i,1);p11(j,4)=aa(i,3)
j=j+1
p11(j,1) = 1; p11(j,2)= aa(i,3);p11(j,3)= aa(i,2);p11(j,4)=aa(i,1)
j=j+1
p11(j,1) = 1; p11(j,2)= aa(i,1);p11(j,3)= aa(i,3);p11(j,4)=aa(i,2)
j=j+1
p11(j,1) = -1; p11(j,2)= aa(i,3);p11(j,3)= aa(i,1);p11(j,4)=aa(i,2)
j=j+1
p11(j,1) = -1;p11(j,2)= aa(i,2);p11(j,3)= aa(i,3);p11(j,4)=aa(i,1)
end do

!do i=1,6
!write(6,*) p11(i,:)
!end do

!-----------------------------------------------------------------
! Generating phi_12
do i=1,9
j=((i-1)*4)
j=j+1
p12(j,1) =1; p12(j,2)= aa(i,3);p12(j,3)= aa(i,2);p12(j,4)=aa(i,1)
j=j+1
p12(j,1) =1; p12(j,2)= aa(i,1);p12(j,3)= aa(i,3);p12(j,4)=aa(i,2)
j=j+1
p12(j,1) =-1; p12(j,2)= aa(i,3);p12(j,3)= aa(i,1);p12(j,4)=aa(i,2)
j=j+1
p12(j,1) =-1; p12(j,2)= aa(i,2);p12(j,3)= aa(i,3);p12(j,4)=aa(i,1)

end do

!do i=1,4
!write(6,*) p12(i,:)
!end do

!----------------------------------------------------------------------
! loading the integrals - Overlap , One electrong, two electron 

open(unit = 101, file = 'ham_ov.dat', status = 'old', action = 'read')
 do i = 1,3
  read(101,*) Ov(i,1),Ov(i,2),Ov(i,3)
 end do

 do i = 1,3 
  read(101,*) HOv(i,1),HOv(i,2),HOv(i,3)
 end do

! do i = 1,3
!  write(6,*) Ov(i,1),Ov(i,2),Ov(i,3)
! end do

! do i = 1,3
!  write(6,*) HOv(i,1),HOv(i,2),HOv(i,3)
! end do



open(unit = 102, file = '2_ele.dat', status = 'old', action = 'read')

! do i = 1,81
!  read(102,*) a1(i),b1(i),c1(i),d1(i),e1(i)
!  EEOv(a1(i),c1(i),b1(i),d1(i))=e1(i)
!  EEOv(a1(i),c1(i),b1(i),d1(i))=0.d0
! end do

do i=1,81
 read(102,*) e1(i)
end do

ii=0
do i=1,3;
 do j=1,3;
  do k=1,3;
   do l=1,3;
       ii=ii+1;
       EEOv(i,j,k,l)=e1(ii);
!       EEOv(i,j,k,l)=0.d0
   end do
  end do
 end do
end do


! do i = 1,81
!  write(6,*) EEOv(a1(i),b1(i),c1(i),d1(i))
! end do

! Orthoganizlation of basis using canonical orthogonalization
! [U1,D1]=eig(S1);

!do i=1,3
!   write(6,*)Ov(i,:)
!end do

S1=Ov

LWORK=8
call DSYEV( 'V', 'U', 3, Ov, 3, W, WORK, LWORK, INFO )

! do i=1,3
!   write(6,*)S1(i,:)
! end do

!-------------------------------------------------------------
! Working in orthogonal basis 

 do i=1,3
    do  j=1,3
      X1(i,j)=Ov(i,j)/sqrt(W(j))
    end do 
 end do 

! Symmetric Orthogonalization
! do i=1,3
!    do  j=1,3
!      X1(i,j)=Ov(i,j)/sqrt(W(j))
!    end do 
! end do 
!Xsymm=matmul(X1,transpose(Ov))
!
!X1=Xsymm
!X1=0.d0

!do i=1,3
!X1(i,i)=1.d0
!enddo 


S=matmul(transpose(X1),matmul(S1,X1))
H=matmul(transpose(X1),matmul(HOv,X1))


do i=1,3
!write(6,*)S(i,:)
enddo 
do i=1,3
  do j=1,3
    do k=1,3
      do l=1,3

          EE(i,j,k,l)=0.d0

         do a=1,3
          do b=1,3
           do c=1,3
            do d=1,3
              EE(i,j,k,l)=EE(i,j,k,l)+X1(a,i)*X1(b,j)*X1(c,k)*X1(d,l)*EEOv(a,b,c,d)
            end do 
           end do 
          end do 
         end do 

       end do 
     end do 
   end do 
 end do 


!--------------------------------------------------------------------------------
! computing the Hessian matrix

do i=1,9
    do j=1,9
       call one_element(p11(6*(i-1)+1:6*(i-1)+6,:),p12(4*(i-1)+1:4*(i-1)+4,:), &
                        p11(6*(j-1)+1:6*(j-1)+6,:),p12(4*(j-1)+1:4*(j-1)+4,:),val,S,H,EE)

      H_f(i,j)=val
      
      if(abs(H_f(i,j)).lt.1.D-8) then
       H_f(i,j)=0.d0
      end if

    end do
end do

do i=1,9
! write(23,10)H_f(i,:)
enddo 

!--------------------------------------------------------------------------------
!Digonalizing the matrix

LWORK_f=26
call DSYEV( 'V', 'U', 9, H_f, 9, W_f, WORK_f, LWORK_f, INFO )

 write(6,*) "Spectrum for Sz=1/2"
do i=1,9
 write(6,10)W_f(i)
10 format(8f13.8)
enddo

!H1=H
!call DSYEV( 'V', 'U', 3, H1, 3, W, WORK, LWORK, INFO )
!write(6,*)W


end program half_spin


subroutine one_element(a1,b1,a2,b2,val,S,H,EE)
integer i,j,k,l
integer a1(6,4),a2(6,4),b1(4,4),b2(4,4)
integer aa1(6,3),aa2(6,3),bb1(4,3),bb2(4,3)
real*8 v,val,val1,val2,nm,norm_aa1,norm_aa2,norm_bb1,norm_bb2
real*8  S(3,3),H(3,3),EE(3,3,3,3)
!--------------------------------------------------------------------------------
! This subroutines compute 
!  part <PHI(i)|H|PHI(j)>  where 
! PHI(i)=phi_11(i)*alpha_1 + phi_12(i)*alpha_2
! PHI(i)=phi_11(j)*alpha_1 + phi_12(j)*alpha_2
! and <alpha_1|alpha_2> = 0 i.e orthogonal
!--------------------------------------------------------------------------------

aa1=a1(:,2:4)   ! storing the components of phi_11(i) 
aa2=a2(:,2:4)   ! storing the components of phi_12(j)


 norm_aa1=0.d0;
 norm_aa2=0.d0;

! Computing norm square of phi_11(i)

 do i=1,6
  do k=1,6
  call norm1(aa1(i,:),aa1(k,:),nm,S,H,EE)
  norm_aa1=norm_aa1+a1(i,1)*a1(k,1)*nm/4.d0
  end do
 end do

! Computing norm square of phi_11(j)

 do i=1,6
  do k=1,6
   call norm1(aa2(i,:),aa2(k,:),nm,S,H,EE)
   norm_aa2=norm_aa2+a2(i,1)*a2(k,1)*nm/4.d0
  end do
 end do

! Computing the <phi_11(i)|H|phi_11(j)>

val1=0
 do i=1,6
   do j=1,6
       call ham(aa1(i,:),aa2(j,:),v,S,H,EE)
       val1=val1+a1(i,1)*a2(j,1)*v/4.d0
   end do
 end do



bb1=b1(:,2:4);    ! storing the components of phi_12(i)
bb2=b2(:,2:4);    ! storing the components of phi_12(j)


 norm_bb1=0.d0;
 norm_bb2=0.d0;

! Computing norm square of phi_12(i)

 do i=1,4
  do k=1,4
   call norm1(bb1(i,:),bb1(k,:),nm,S,H,EE)
   norm_bb1=norm_bb1+b1(i,1)*b1(k,1)*nm*3.d0/4.d0
  end do 
 end do 

! Computing norm square of phi_12(j)

 do i=1,4
  do k=1,4
  call norm1(bb2(i,:),bb2(k,:),nm,S,H,EE)
  norm_bb2=norm_bb2+b2(i,1)*b2(k,1)*nm*3.d0/4.d0
  end do
 end do

! Computing <phi_12(i)|H|phi_12(j)>

val2=0;
  do i=1,4
   do j=1,4
       call ham(bb1(i,:),bb2(j,:),v,S,H,EE)
       val2=val2+b1(i,1)*b2(j,1)*v*3.d0/4.d0
   end do
 end do


! This is more crucial step 
!
!1. This assuming the alpha1 and alpha2 are already orthonormal and normalizing it with norms
!   of individual part
!
   val= (val1+val2)/(sqrt((norm_aa1+norm_bb1))*sqrt((norm_aa2+norm_bb2)))

!2. This assuming the alpha1 and alpha2 are already orthonormal and using the direct expression
!   the book by Ruben Paunz
!
! val=(val1+val2)/(6.d0)

!3. This assuming the alpha1 and alpha2 are not normalized and get constructed from conjugate young
!   tabulae
!   alpha1= sqrt(3)[beta*alpha*alpha - alpha*alpha*beta]
!   alpha2= [2*alpha*alpha*beta - beta*lpha*alpha - alpha*beta*alpha ]
!   <alpha1|alpha1> = 6 ; <alpha2|alpha2>=6
!
!  val= (val1*6.d0+val2*3.d0)/(sqrt((norm_aa1*6.d0+norm_bb1*3.d0))*sqrt((norm_aa2*6.d0+norm_bb2*3.d0)))
end subroutine one_element 


subroutine ham(a,p,v,S,H,EE)
!common S(3,3),H(3,3),EE(3,3,3,3)
real*8  S(3,3),H(3,3),EE(3,3,3,3)
integer a(3),p(3)
real*8  v
! this subroutine computes <123|H|456>

v = H(a(1),p(1))*S(a(2),p(2))*S(a(3),p(3))+ &
    H(a(2),p(2))*S(a(1),p(1))*S(a(3),p(3))+ &
    H(a(3),p(3))*S(a(1),p(1))*S(a(2),p(2))+ &
    EE(a(1),a(2),p(1),p(2))*S(a(3),p(3))+   &
    EE(a(1),a(3),p(1),p(3))*S(a(2),p(2))+   &
    EE(a(2),a(3),p(2),p(3))*S(a(1),p(1));

end subroutine ham

subroutine norm1(a,p,nm,S,H,EE)

integer a(3),p(3)
!common S(3,3),H(3,3),EE(3,3,3,3)
real*8  S(3,3),H(3,3),EE(3,3,3,3)
real*8 nm
! this subroutine computes the overalp <123|456>
nm=S(a(1),p(1))*S(a(2),p(2))*S(a(3),p(3));

end subroutine norm1


