!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Subrouting for computing action of elementary transposition
!on standard young frame
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine action_perm_young(p,q,young,f_lam,u,v,m1) 
integer*8 p,q,f_lam,u,v
integer*8 young(f_lam,u+v)
integer*8 i,j,l,loc_p,loc_q,tmp(u+v),tm,m1(f_lam)

do i = 1, f_lam

!    determining the position of (p,q) in young frame T_{i}
    do j = 1, u+v
        if (young(i,j) .eq. p) then
          loc_p = j
        exit
       endif
    end do

    do j = 1, u+v
        if (young(i,j) .eq. q) then
          loc_q = j
          exit
        endif
    end do
   
!   action of permutation on young frame T_{i} 
    tmp(:) = young(i,:)

    tm = tmp(loc_p)
    tmp(loc_p) = tmp(loc_q)
    tmp(loc_q) = tm

!   determining if (p,q)T_{i} = T_{j} and storing it in m1
    do j=1,f_lam
  
       flag = 0

         do k=1,u+v
            if (young(j,k).eq.tmp(k)) then
               flag = flag + 1
            end if
         end do
   
          if (flag.eq.u+v) then
!            print*, "Permutation action successful on = ", p,q,j
            m1(i)=j
          end if 

     end do
end do

end subroutine action_perm_young
