program generate
implicit none 
integer*8 i,j,k,l
real*8 a,b,c,d
integer*8 el, ba, sp, u_el
integer*8 u, v

write(*,*) "Number of electron"
read(*,*) el

write(*,*) "Spin multiplicity of the system"
read(*,*) sp

write(*,*) "No. of basis"
read(*,*) ba

! writing in output file fort.22
write(22,*) "Number of electron = ", el

write(22,*) "Spin multiplicity of the system = ", sp

write(22,*) "No. of basis = ", ba 


!print*, el,sp,ba

! No. of unpaired electron = spin multiplicity - 1

u_el = sp - 1

! Because electrons are fermions with only two posible spin alpha
! and beta corresponding to spin 1/2 and -1/2, the young shape is always 
! [u,v], where u = (el + u_el)/2  and v = (el - u_el)/2 


u = (el+u_el)/2
v = (el-u_el)/2
 
if (u+v.ne.el) then 

print*, "!!!!!!! Read about spin multiplicity before running this program"

stop 
end if 


!print*,u_el,u,v

call young_frame(u,v)

end program 
