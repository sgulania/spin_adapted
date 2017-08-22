subroutine young_frame(u,v)
implicit none

integer*8 u,v
integer*8 i,j,k,l
integer*8 ii,jj,m,kk
integer*8 a(u),b(v)
integer*8 hka(u),hkb(v),f_lam
integer*16 h1,h2
integer*8 a1(u),b1(v),a11(u),b11(v)
integer*8 flag,flag2 
integer*8 pp
integer*8 alpha
integer*8, dimension(:), allocatable  :: comb
integer*8, dimension(:,:), allocatable  :: rem
integer*8, dimension(:,:), allocatable  :: remb
integer*8, dimension(:,:), allocatable :: young
integer AllocateStatus
integer*8 v1,num
integer*16  tabl(15)

! a = first row of young frame 
! b = second row of young frame
! u = size (a)
! v = size (b)
! hka = hooklength of first row
! hkb = hooklength of second row
! f_lam = no. of standard young tabulae

print*,"u and v", u,v
do i = 1,u
  a(i) = i
end do

do j=1,v
  b(j) = u+j
end do 

! Computing number of standard yound tabuleaux using
! hook length formula f(lam) = N! /(h1*h2* ... hn)


! hook lengths are stored in hka and hkb 

if (v .ne. 0) then
 
  do i=1,v
    hka(i) = 2 + (u-i) 
  end do 

  do i=v+1,u
   hka(i) = 1 + (u-i)
  end do 

 else 
 
  do i=1,u
    hka(i) = 1 + (u-i) 
  end do

end if

do j=1,v
 hkb(j) = 1 + (v-j)
end do 

!write(22,*) "hook lengths"

!write(22,*) hka
!write(22,*) hkb

! computing number of standard young tabulaex

h1=1 
do i=1,u
  h1=h1*hka(i)
end do  

h2=1
do j=1,v
  h2=h2*hkb(j)
end do

if (u+v.le.20) then
f_lam = Factorial(u+v)/(h1*h2)

print*,"Fatcorial of ", u+v,"=", Factorial (u+v)
print*, "No. of standard yound tabulae = ", f_lam
write(22,*), "No. of standard yound tabulae = ", f_lam 

end if 

do i=1,13
read(51,*) num,tabl(i)
end do

if (u+v.gt.20) then
f_lam = tabl(u+v-20)/(h1*h2)
print*,"Fatcorial of ", u+v,"=", tabl(u+v-20)
print*, "No. of standard yound tabulae = ", f_lam
write(22,*), "No. of standard yound tabulae = ", f_lam
end if

!allocate (young(f_lam,u+v))

allocate ( young(f_lam,u+v), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

! generating all the standard young tabulae

write(22,*) "Standard Young Frame = 1"

write(22,*) a
write(22,*) b

write(22,*) " "

young(1,1:size(a)) = a(:)
young(1,size(a)+1:size(a)+size(b))=b(:)

m=2

!print*, young(1,:),"above"

a1 = a
b1 = b

! no of ways to pull down numbers from above row to below row

v1 = v
if (u.eq.v) then 
 v1 =v-1
end if

pp=0
do i=1,v1

  k = Factorial (u-1)/(Factorial(i)*Factorial (u-1-i))
  l = Factorial (v)/(Factorial(i)*Factorial (v-i))
 
!  print*, k,l,k*l,"k and l"
! pp= pp+k*l

  alpha = 1
!  print*, u-1,i,"comb"
  allocate ( comb(i))
  Call gen (alpha,u-1,i,comb,1)
  deallocate (comb)
  
  CLOSE ( 23, STATUS='KEEP') 
  allocate (rem(k,i))

  do j=1,k
  read(23,*) rem(j,:)  
  end do

  do j=1,k
!  write(*,*) rem(j,:)
  end do

!----------------------------------------

  alpha = 1
  allocate ( comb(i))
  Call gen (alpha,v,i,comb,2)
  deallocate (comb)

  CLOSE ( 24, STATUS='KEEP')
  allocate (remb(l,i))

  do j=1,l
  read(24,*) remb(j,:)
  end do

  do j=1,l
!  write(*,*) remb(j,:)
  end do



!----------------------------------------



 do kk = 1,k 
  do jj = 1,l
  do j=1,i
   
!   print*, rem(kk,j)+1, remb(jj,j),jj,"inside"
!   print*, a1(rem(kk,j)+1), b1(remb(jj,j)),"inside"
   Call Swap(a1(rem(kk,j)+1),b1(remb(jj,j)))
  end do
!  stop
  Call Sort (a1,u)
  Call Sort (b1,v)

!-------------------------------
  flag = 1
 do j=1,v
  if (a1(j).gt.b1(j)) then
   flag=flag*0
  end if
 end do

 if (flag.eq.1) then

 young(m,1:size(a1)) = a1(:)
 young(m,size(a1)+1:size(a1)+size(b1))=b1(:)
!print*,young(m,:)
 write(22,*) "Standard Young Frame = ",m

 write(22,*) a1
 write(22,*) b1

 write(22,*) " "
 m=m+1
 end if

 a1 = a
 b1 = b
  
!-----------------------------

 end do 
 end do
  deallocate (rem)
  deallocate (remb)
  CLOSE ( 23, STATUS='DELETE' )
  CLOSE ( 24, STATUS='DELETE' )
!print*, young(1,:) ,"loop",i
end do

stop
print*, "Satndard Young Tabulae constructed, checking degenracy"
!print*, young(1,:)

!stop
flag2 = 0
do j=1,f_lam
 do k=j+1,f_lam
  
flag = 0

  do i=1,u+v
   if (young(j,i).eq.young(k,i)) then
   flag = flag + 1
   end if
  end do
   
  if (flag.eq.u+v) then
   flag2 =flag2+1
   print*, "degeneracy number = ", j,k
   print*, young(j,1:u)
   print*, young(j,u+1:u+v)
  end if 

 end do
end do

print*,"No. of degeneracy = ", flag2
deallocate (young)


!-------------------------------------------------------------------

CONTAINS

! http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90

! --------------------------------------------------------------------
! INTEGER*8 FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

   INTEGER*8 FUNCTION  FindMinimum(x, Start, End)
      IMPLICIT  NONE
      INTEGER*8, DIMENSION(1:), INTENT(IN) :: x
      INTEGER*8, INTENT(IN)                :: Start, End
      INTEGER*8                            :: Minimum
      INTEGER*8                            :: Location
      INTEGER*8                            :: i

      Minimum  = x(Start)               ! assume the first is the min
      Location = Start                  ! record its position
      DO i = Start+1, End               ! start with next elements
         IF (x(i) < Minimum) THEN       !   if x(i) less than the min?
            Minimum  = x(i)             !      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location            ! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      INTEGER*8, INTENT(INOUT) :: a, b
      INTEGER*8                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, Size)
      IMPLICIT  NONE
      INTEGER*8, DIMENSION(1:), INTENT(INOUT) :: x
      INTEGER*8, INTENT(IN)                   :: Size
      INTEGER*8                               :: i
      INTEGER*8                               :: Location

      DO i = 1, Size-1                      ! except for the last
         Location = FindMinimum(x, i, Size) ! find min from this to last
         CALL  Swap(x(i), x(Location)) ! swap this and the minimum
      END DO
   END SUBROUTINE  Sort




! Stanford - https://web.stanford.edu/class/me200c/tutorial_90/08_subprograms.html
RECURSIVE FUNCTION Factorial(n)  RESULT(Fact)

IMPLICIT NONE
INTEGER*8 :: Fact
INTEGER*8, INTENT(IN) :: n

IF (n == 0) THEN
   Fact = 1
ELSE
   Fact = n * Factorial(n-1)
END IF

END FUNCTION Factorial


end subroutine 

recursive subroutine gen (m,n_max,m_max,comb,l)

    implicit none
    integer*8 m
    integer*8 n, n_max,m_max
    integer*8 comb(m_max)
    integer l
    
    if (m > m_max) then
     if (l.eq.1) then
      write (23,*) comb
      write(223,*) comb
     end if 
      
     if (l.eq.2) then
      write(24,*) comb
      write(224,*) comb
     end if 
    else
      do n = 1, n_max
        if ((m == 1) .or. (n > comb (m - 1))) then
          comb (m) = n
          call gen (m + 1,n_max,m_max,comb,l)
        end if
      end do
    end if

end subroutine gen



