subroutine irrep_generate(young,u,v,f_lam)
integer*8 i,j
integer*8 u,v,f_lam
integer*8 young(f_lam,u+v)
integer*8 order 
real*8 mat(u+v-1,f_lam*(f_lam+1)/2) ! storing upper triangular of iirep
real*8 d
integer*8 m(u+v-1,f_lam)
integer*8 k
integer*8 flag
order = u+v
mat=0


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Computing diagonal terms of Irrep for elementary
!transposition using young algorithm
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

do i=1,order-1
!print*,"Diagonal for U(",i,i+1,")" 
  flag = 1
  do k=1,f_lam
    call distance(i,i+1,order,young(k,:),u,v,d)
     mat(i,flag) = -1.d0/d
     flag = flag+1 + f_lam -k
  end do
end do



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Finding action of elementary transposition on standard
! young frame
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

m=0  !initialize

do i=2,order-1
      call action_perm_young(i,i+1,young,f_lam,u,v,m(i,:))
      !print*,m(i,:)
end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Computing upper off diagonal  terms of Irrep for elementary
!transposition using young algorithm
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
do i=2,order-1

!      print*,"off diagonal terms ofU(",i,i+1,")"
     
  flag = 0
  do j=1,f_lam
     flag = flag + 1
    do k=j+1,f_lam
      flag = flag +1 
       
      if (m(i,j).eq.k) then
          call distance(i,i+1,order,young(j,:),u,v,d)
          mat(i,flag) = sqrt(1-1/(d*d))
      endif

    enddo 
  end do

end do

open (unit = 23, file = "Irrep.dat")

do i=1,u+v-1
  write(23,*) i,i+1
  do k=1,f_lam*(f_lam+1)/2 
     write(23,*) mat(i,k)
  end do
end do
close (23) 

end subroutine irrep_generate


