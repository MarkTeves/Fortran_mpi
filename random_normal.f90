subroutine random_stdnormal(x)
    ! Box-Muller transformation
   implicit none
   real,intent(out) :: x
   real,parameter :: pi=3.14159265
   real :: u1,u2
   call random_stduniform(u1)
   call random_stduniform(u2)
   x = sqrt(-2*log(u1))*cos(2*pi*u2)
end subroutine random_stdnormal
    
subroutine random_stduniform(u)
   implicit none
   real,intent(out) :: u
   real :: r
   call random_number(r)
   u = 1 - r
end subroutine random_stduniform