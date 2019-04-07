!------------------------------------------------------------- 
! Subroutine for orthogonalizing the basis and 
! transforming the matrix accordingly
!------------------------------------------------------------- 

subroutine ortho_transformation(overlap,matrix,dimen)
 implicit none
 
 integer i,j
 integer INFO, LWORK, dimen
 real*8  overlap(dimen,dimen),matrix(dimen,dimen) 
 real*8  X1(dimen,dimen),W(dimen),WORK(3*dimen-1),S1(dimen,dimen)


!-------------------------------------------------------------
! Working in orthogonal basis 

 S1=overlap

 LWORK=3*dimen-1
 call DSYEV( 'V', 'U', dimen, overlap, dimen, W, WORK, LWORK, INFO )

! do i=1,dimen
! print*, 'overlap', W(i)
! end do 
 do i=1,dimen
    do  j=1,dimen
      X1(i,j)=overlap(i,j)/sqrt(W(j))
    end do
 end do


! Symmetric Orthogonalization
! do i=1,3
!    do  j=1,3
!      X1(i,j)=overlap(i,j)/sqrt(W(j))
!    end do 
! end do 
!Xsymm=matmul(X1,transpose(overlap))
!
!X1=Xsymm
!X1=0.d0

!do i=1,3
!X1(i,i)=1.d0
!enddo 

overlap=matmul(transpose(X1),matmul(S1,X1))
matrix =matmul(transpose(X1),matmul(matrix,X1))

end subroutine 
