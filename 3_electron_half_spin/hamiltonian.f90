!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine for evaluating Hamiltonian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hamiltonian(cof_p11,cof_p21,p11,p21,list_p11,list_p21,n_p11,n_p21)
USE variables
implicit none

integer i,j,k,l,m,n,o
integer n_p11,n_p21
integer list_p11(n_weyl+1),list_p21(n_weyl+1),LWORK_f,LWORK_S,INFO
real*8 cof_p11(n_p11),cof_p21(n_p21)
integer p11(n_p11,n_el), p21(n_p21,n_el)
real*8 H_f(n_weyl,n_weyl), S_f(n_weyl,n_weyl), W_f(n_weyl),WORK_f(3*n_weyl-1)
real*8 Ov(n_weyl,n_weyl), W_S(n_weyl),WORK_S(3*n_weyl-1) 
real*8 X_S(n_weyl,n_weyl), H_S(n_weyl,n_weyl)
real*8 X_symm(n_weyl,n_weyl)
integer i11,i11n,j11,j11n,i21,i21n,j21,j21n
real*8 val,over



!write(6,*) n_weyl
!stop
!do i=1,34
!write(6,*) cof_p21(i)
!end do

!do i=1,list_p11(n_weyl+1)
!write(6,*) i, cof_p11(i),p11(i,:)
!end do

do i=1,n_weyl
    do j=1,n_weyl

 !      write(6,*) i,j
       i11  = list_p11(i)+1
       i11n = list_p11(i+1)

       j11  = list_p11(j)+1
       j11n = list_p11(j+1)
      
       i21  = list_p21(i)+1
       i21n = list_p21(i+1)
       
       j21  = list_p21(j)+1
       j21n = list_p21(j+1)
       
  !     write(6,*) i11,i11n,j11,j11n,i21,i21n,j21,j21n

  !     write(6,*) p11(i11:i11n,:)
       
       call one_element(p11(i11:i11n,:),p21(i21:i21n,:),&
                        p11(j11:j11n,:),p21(j21:j21n,:),&
                    cof_p11(i11:i11n),cof_p21(i21:i21n),&
                    cof_p11(j11:j11n),cof_p21(j21:j21n),&
                    i11n-i11+1,i21n-i21+1,j11n-j11+1,j21n-j21+1,val,over)
      S_f(i,j)=over
      if(abs(over).lt.1.D-8) then
        S_f(i,j)=0.d0
      end if
      !write(6,*) S_f(i,j),i,j     
      !if(i.eq.1 .and. j.eq.1) then
      ! write(6,*) i,j,list_p11(i)+1
      ! write(6,*) "out"
      ! write(6,*) val,"val out" 
      !end if

      H_f(i,j)= val
      !H_f(j,i)= H_f(i,j)

      if(abs(H_f(i,j)).lt.1.D-8) then
       H_f(i,j)=0.d0
      end if
    
!write(6,*) "Ham"
    end do
end do


do i=1,n_weyl
 write(24,"(9f13.8)")S_f(i,:)
enddo

do i=1,n_weyl
 write(23,"(9f13.8)")H_f(i,:)
enddo



 Ov=S_f                                                                                    
                                                                                          
 LWORK_S=3*n_weyl-1                                                                          
 call DSYEV( 'V', 'U', n_weyl, Ov, n_weyl, W_S, WORK_S, LWORK_S, INFO )                           
                                                                                          
! do i=1,n_weyl                                                                             
!    do  j=1,n_weyl                                                                         
!      X_S(i,j)=Ov(i,j)/sqrt(W_S(j))                                                          
!    end do                                                                                
! end do                                                                                   


! Symmetric Orthogonalization                                                             
 do i=1,n_weyl                                                                                
    do  j=1,n_weyl                                                                           
      X_S(i,j)=Ov(i,j)/sqrt(W_S(j))                                                         
    end do                                                                               
 end do                                                                                  
X_symm=matmul(X_S,transpose(Ov))                                                           
!                                                                                         
X_S=X_symm                                     


do i=1,n_weyl                                                                             
 write(27,*)W_S(i)                                                             
enddo
                            
do i=1,n_weyl                                                                             
 write(26,"(9f13.8)")X_S(i,:)                                                            
enddo                                                              
                                                                                          
H_S=matmul(transpose(X_S),matmul(H_f,X_S))                                                     

!Digonalizing the matrix

do i=1,n_weyl                                                                             
 write(25,"(9f13.8)")H_S(i,:)                                                             
enddo  


LWORK_f=3*n_weyl-1
call DSYEV( 'V', 'U', n_weyl, H_S, n_weyl, W_f, WORK_f, LWORK_f, INFO)

 write(6,*) "Spectrum for Sz=1/2"
do i=1,n_weyl
 write(6,"(I5,9f13.8)")i,W_f(i)+nuc_repul
!10 format (4I5,9f13.8)
enddo


end subroutine hamiltonian
