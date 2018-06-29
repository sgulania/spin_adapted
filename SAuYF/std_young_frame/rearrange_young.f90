!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Subroutine for arrangment of young frames using Yamanouchi symbols.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine rearrange_young(young,u,v,f_lam)
integer*8 u,v,f_lam
integer*8 young(f_lam,u+v)
integer*8 num(f_lam,u+v)
integer*8 i,k,j,ord
integer*8 tmp(u+v)

!cccccccccccccccccccccccccccccccc
!assigning Yamanouchi symbols.
!cccccccccccccccccccccccccccccccc
do i=1,f_lam

 do k=1,u
   num(i,young(i,k))=1
 enddo

 do k=u+1,u+v
   num(i,young(i,k))=2
 end do

!print*,num(i,:)
enddo

!cccccccccccccccccccccccccccccccccccccccccccccccc
!arrange the Yamanouchi numbers along with young
!cccccccccccccccccccccccccccccccccccccccccccccccc
202 continue

do i=1,f_lam
 do j=i+1,f_lam

    ord = u+v      
201 continue
    if (num(i,ord).lt.num(j,ord)) then
        tmp(:)  = num(j,:)
        num(j,:) = num (i,:)
        num(i,:) = tmp(:)
        tmp(:) = young (j,:)
        young (j,:) = young (i,:)
        young (i,:) = tmp(:)
        go to 202

    else if (num(i,ord).eq.num(j,ord)) then
        ord = ord-1
        go to 201
    else  
       go to 203
    end if

203 continue

 enddo
enddo



end subroutine rearrange_young
