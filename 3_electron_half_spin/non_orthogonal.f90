subroutine non_orthogonal(Ov,HOv,EEOv)
 USE variables
 implicit none
 
 real*8  Ov(n_bas,n_bas),HOv(n_bas,n_bas),EEOv(n_bas,n_bas,n_bas,n_bas)  


 allocate ( S(n_bas,n_bas),H(n_bas,n_bas),EE(n_bas,n_bas,n_bas,n_bas), &
          STAT = AllocateStatus)

!-------------------------------------------------------------
! Working in one electron non-orthogonal basis 

S = Ov
H = HOv
EE = EEOv

end subroutine 
