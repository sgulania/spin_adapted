subroutine elec_full_ci
USE variables                                                                             
implicit none                                                                             
integer a(3),p(3)                                                                         
real*8  v          
real*8 H2(3,3)
real*8 W_f(3),WORK_f(3*3-1)
integer LWORK_F, INFO,i  
real*8 Energy


H2(1,1) = 2.d0*H(1,1)+EE(1,1,1,1);
H2(1,2) = EE(1,1,2,2);
H2(2,1) = H2(1,2);
H2(1,3) = EE(1,1,3,3);
H2(3,1) = H2(1,3);
H2(2,3) = EE(2,2,3,3);
H2(3,2) = H2(2,3);
H2(2,2) = 2.d0*H(2,2) + EE(2,2,2,2);
H2(3,3) = 2.d0*H(3,3) + EE(3,3,3,3);

print*,EE(1,1,1,1),H(1,1),EE(2,2,2,2),H(2,2),EE(3,3,3,3),H(3,3)
print*,H2(1,1),H2(2,2),H2(3,3)

LWORK_f=3*3-1
call DSYEV( 'V', 'U', 3, H2, 3, W_f, WORK_f, LWORK_f, INFO)                    

Energy = 0.d0
do i=1,3                                                                            
! Energy = Energy+ 3.d0*W_f(i)-3.d0*H(i,i)                                               
print*, W_f(i)
enddo      


call DSYEV( 'V', 'U', 3, H, 3, W_f, WORK_f, LWORK_f, INFO)                               
                                                                                          
Energy = 0.d0                                                                             
do i=1,3                                                                                  
! Energy = Energy+ 3.d0*W_f(i)-3.d0*H(i,i)                                                
print*, W_f(i)                                                                            
enddo

print*, 'Energy', Energy

end subroutine
