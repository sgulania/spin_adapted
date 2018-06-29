!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Subroutine for determining distance (d) between two elements in 
!pu and pv in standard young fram using axial distance rule 
!Step moved in left and down are counted as positive
!and negative when moved right and up
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine distance (pu,pv,i_max,a,u,v,d)
integer*8 i,j
integer*8 i_max
integer*8 pu,pv,u,v
integer*8 loc_u, loc_v
integer*8 a(i_max)
real*8 d

do i = 1, i_max
    if (a(i) .eq. pu) then
        loc_u = i
        exit
    endif
end do

do i = 1, i_max                                                                    
    if (a(i) .eq. pv) then 
        loc_v = i         
        exit                                                                 
    endif                                                             
end do

!d = loc_u -loc_v

if (loc_u .le. u .and. loc_v .le. u) then
    d = loc_u-loc_v
    
else if (  loc_u .gt. u .and. loc_v .gt. u ) then
    d = loc_u-loc_v

else if (loc_u .le. u   .and. loc_v .gt.u ) then
    d = loc_u - (loc_v -u) + 1
else 
    d = - (loc_v - (loc_u -u) + 1)

endif 

end subroutine distance 
