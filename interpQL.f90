module interpQL
    contains 

subroutine interp1D_QL(sigma,nx,x,f,xi,outval)
implicit none
integer, intent(in)		:: nx
real, intent(in)	:: x(nx),f(nx),xi
real ::  xin
real :: interp_QL
real :: sigma
real, intent(out) :: outval
!integer :: llc(1), n
integer :: klo, khi, k
integer :: llc(1)

  ! Define the box lower left front corner 
  ! llc=(/find_grid_point(x, nx, xi), find_grid_point(y, ny, yi), find_grid_point(z, nz, zi)/)
        
! X term

       klo = 1
       khi = nx
       do while (khi-klo>1)
          k=(khi+klo)/2
          if (x(k)>xi) then
            khi=k
          else
            klo=k
          end if
        end do
        llc(1) = klo
    

  ! Normalize points to make the length exactly 1 in each direction
   xin = (xi - x(llc(1))) / ( x(llc(1)+1) - x(llc(1)))

   ! Formula  
!Vxyz = &
! V000 * (1 - x) * (1 - y)  +
! V100 *    x    * (1 - y)  + 



   ! CRRA
   interp_QL = &
     f(llc(1)    )   * (1-xin)  + &
     f(llc(1)+1  )      *   xin     

outval = ((interp_QL**(1-sigma))/(1-sigma))    
end subroutine    
    

subroutine interp2D_QL(sigma,nx,x,ny,y,f,xi,yi,outval)
implicit none
integer, intent(in)		:: nx,ny
real, intent(in)	:: x(nx),y(ny),f(nx,ny),xi,yi
real ::  xin, yin
real :: interp_QL
real :: sigma
integer :: klo, khi, k
real, intent(out) :: outval
integer :: llc(2), n

          
  ! Define the box lower left front corner 
  ! llc=(/find_grid_point(x, nx, xi), find_grid_point(y, ny, yi), find_grid_point(z, nz, zi)/)
        
! X term
       klo = 1
       khi = nx
       do while (khi-klo>1)
          k=(khi+klo)/2
          if (x(k)>xi) then
            khi=k
          else
            klo=k
          end if
        end do
        llc(1) = klo        
        
! Y term
       klo = 1
       khi = ny
       do while (khi-klo>1)
          k=(khi+klo)/2
          if (y(k)>yi) then
            khi=k
          else
            klo=k
          end if
        end do
        llc(2) = klo
       

  ! Normalize points to make the length exactly 1 in each direction
   xin = (xi - x(llc(1))) / ( x(llc(1)+1) - x(llc(1)))
   yin = (yi - y(llc(2))) / ( y(llc(2)+1) - y(llc(2)))

   ! Formula  
!Vxyz = &
! V000 * (1 - x) * (1 - y)  +
! V100 *    x    * (1 - y)  + 
! V010 * (1 - x) *    y     + 
! V110 *    x    *    y      


   ! CRRA
   interp_QL = &
     f(llc(1)  , llc(2)    )   * (1-xin) * (1-yin) + &
     f(llc(1)+1, llc(2)    )   *   xin   * (1-yin) + &
     f(llc(1)  , llc(2)+1  )   * (1-xin) *  yin    + &
     f(llc(1)+1, llc(2)+1  )   *   xin   *  yin  

outval  = ((interp_QL**(1-sigma))/(1-sigma))


end subroutine
    
    
subroutine interp3D_QL(sigma,nx,x,ny,y,nz,z,f,xi,yi,zi,outval)
implicit none
integer, intent(in)		:: nx,ny,nz
real, intent(in)	:: x(nx),y(ny),z(nz),f(nx,ny,nz),xi,yi,zi 
real ::  xin, yin, zin
real :: interp_QL
real :: sigma
integer :: klo, khi, k
real, intent(out) :: outval
integer :: llc(3), n

          
  ! Define the box lower left front corner 
  ! llc=(/find_grid_point(x, nx, xi), find_grid_point(y, ny, yi), find_grid_point(z, nz, zi)/)
        
! X term

       klo = 1
       khi = nx
       do while (khi-klo>1)
          k=(khi+klo)/2
          if (x(k)>xi) then
            khi=k
          else
            klo=k
          end if
        end do
        llc(1) = klo  
! Y term

       klo = 1
       khi = ny
       do while (khi-klo>1)
          k=(khi+klo)/2
          if (y(k)>yi) then
            khi=k
          else
            klo=k
          end if
        end do
        llc(2) = klo    

! Z term

       klo = 1
       khi = nz
       do while (khi-klo>1)
          k=(khi+klo)/2
          if (z(k)>zi) then
            khi=k
          else
            klo=k
          end if
        end do
        llc(3) = klo    
        
  
  ! Normalize points to make the length exactly 1 in each direction
   xin = (xi - x(llc(1))) / ( x(llc(1)+1) - x(llc(1)))
   yin = (yi - y(llc(2))) / ( y(llc(2)+1) - y(llc(2)))
   zin = (zi - z(llc(3))) / ( z(llc(3)+1) - z(llc(3)))

   ! Formula  
!Vxyz = &
! V000 * (1 - x) * (1 - y) * (1 - z) +
! V100 *    x    * (1 - y) * (1 - z) + 
! V010 * (1 - x) *    y    * (1 - z) + 
! V001 * (1 - x) * (1 - y) *    z    +
! V011 * (1 - x) *    y    *    z    + 
! V101 *    x    * (1 - y) *    z    + 
! V110 *    x    *    y    * (1 - z) + 
! V111 *    x    *    y    *    z

   ! CRRA
   interp_QL = &
     f(llc(1)  , llc(2)  , llc(3)  )   * (1-xin) * (1-yin)*(1-zin) + &
     f(llc(1)+1, llc(2)  , llc(3)  )   *   xin   * (1-yin)*(1-zin) + &
     f(llc(1)  , llc(2)+1, llc(3)  )   * (1-xin) *  yin  *(1-zin)  + &
     f(llc(1)  , llc(2)  , llc(3)+1)   * (1-xin) * (1-yin)*  zin   + &
     f(llc(1)  , llc(2)+1, llc(3)+1)   * (1-xin) *  yin  *  zin    + &
     f(llc(1)+1, llc(2)  , llc(3)+1)   *   xin   * (1-yin)*  zin   + &
     f(llc(1)+1, llc(2)+1, llc(3)  )   *   xin   *  yin  *(1-zin)  + &
     f(llc(1)+1, llc(2)+1, llc(3)+1)   *   xin   *  yin  *  zin   

outval= ((interp_QL**(1-sigma))/(1-sigma))


end subroutine
!---------------------------------------------------------------------------------------------------------!


!--------------------------------------------------------------------------------------------------------------------------------
! This function uses 4-d linear interpolation to estimate the quasi-linear version of the value
! of a function f at point x , y, z, q 
!--------------------------------------------------------------------------------------------------------------------------------  
subroutine interp4D_QL(sigma,nx,x,ny,y,nz,z,nq,q,f,xi,yi,zi,qi,outval)
implicit none
integer, intent(in)		:: nx,ny,nz,nq
real, intent(in)	:: x(nx),y(ny),z(nz),q(nq),f(nx,ny,nz,nq),xi,yi,zi,qi 
real ::  xin, yin, zin, qin
real :: interp_QL
real :: sigma
integer :: klo, khi, k
real, intent(out) :: outval
integer :: llc(4), n
      
        
! X term

       klo = 1
       khi = nx
       do while (khi-klo>1)
          k=(khi+klo)/2
          if (x(k)>xi) then
            khi=k
          else
            klo=k
          end if
        end do
        llc(1) = klo  
! Y term

       klo = 1
       khi = ny
       do while (khi-klo>1)
          k=(khi+klo)/2
          if (y(k)>yi) then
            khi=k
          else
            klo=k
          end if
        end do
        llc(2) = klo    

! Z term

       klo = 1
       khi = nz
       do while (khi-klo>1)
          k=(khi+klo)/2
          if (z(k)>zi) then
            khi=k
          else
            klo=k
          end if
        end do
        llc(3) = klo
 
! Q term

       klo = 1
       khi = nq
       do while (khi-klo>1)
          k=(khi+klo)/2
          if (q(k)>qi) then
            khi=k
          else
            klo=k
          end if
        end do
        llc(4) = klo
        

  
  ! Normalize points to make the length exactly 1 in each direction
   xin = (xi - x(llc(1))) / (x(llc(1)+1) - x(llc(1)))
   yin = (yi - y(llc(2))) / (y(llc(2)+1) - y(llc(2)))
   zin = (zi - z(llc(3))) / (z(llc(3)+1) - z(llc(3)))
   qin = (qi - q(llc(4))) / (q(llc(4)+1) - q(llc(4)))


   ! Formula  
! Vxyzq = &
! V0000 * (1 - x) * (1 - y) * (1 - z) * (1 - q)+
! V1000 *    x    * (1 - y) * (1 - z) * (1 - q)+ 
! V0100 * (1 - x) *    y    * (1 - z) * (1 - q)+ 
! V0010 * (1 - x) * (1 - y) *    z    * (1 - q) +
! V0110 * (1 - x) *    y    *    z    * (1 - q)+ 
! V1010 *    x    * (1 - y) *    z    * (1 - q) + 
! V1100 *    x    *    y    * (1 - z) * (1 - q) + 
! V1110 *    x    *    y    *    z    * (1 - q) +
! V0001 * (1 - x) * (1 - y) * (1 - z) * (1 - q)+
! V1001 *    x    * (1 - y) * (1 - z) * (1 - q)+ 
! V0101 * (1 - x) *    y    * (1 - z) * (1 - q)+ 
! V0011 * (1 - x) * (1 - y) *    z    * (1 - q) +
! V0111 * (1 - x) *    y    *    z    * (1 - q)+ 
! V1011 *    x    * (1 - y) *    z    * (1 - q) + 
! V1101 *    x    *    y    * (1 - z) * (1 - q) + 
! V1111 *    x    *    y    *    z    * (1 - q)   

    interp_QL = &
     f(llc(1)  , llc(2)  , llc(3)   , llc(4)   )   * (1-xin) * (1-yin) * (1-zin) * (1-qin) + &
     f(llc(1)+1, llc(2)  , llc(3)   , llc(4)   )   *   xin   * (1-yin) * (1-zin) * (1-qin) + &
     f(llc(1)  , llc(2)+1, llc(3)   , llc(4)   )   * (1-xin) *  yin    * (1-zin) * (1-qin) + &
     f(llc(1)  , llc(2)  , llc(3)+1 , llc(4)   )   * (1-xin) * (1-yin) *   zin   * (1-qin) + &
     f(llc(1)  , llc(2)+1, llc(3)+1 , llc(4)   )   * (1-xin) *  yin    *   zin   * (1-qin) + &
     f(llc(1)+1, llc(2)  , llc(3)+1 , llc(4)   )   *   xin   * (1-yin) *   zin   * (1-qin) + &
     f(llc(1)+1, llc(2)+1, llc(3)   , llc(4)   )   *   xin   *  yin    * (1-zin) * (1-qin) + &
     f(llc(1)+1, llc(2)+1, llc(3)+1 , llc(4)   )   *   xin   *  yin    *   zin   * (1-qin) + & 
     f(llc(1)  , llc(2)  , llc(3)   , llc(4)+1 )   * (1-xin) * (1-yin) * (1-zin) *   qin  + &
     f(llc(1)+1, llc(2)  , llc(3)   , llc(4)+1 )   *   xin   * (1-yin) * (1-zin) *   qin  + &
     f(llc(1)  , llc(2)+1, llc(3)   , llc(4)+1 )   * (1-xin) *  yin    * (1-zin) *   qin  + &
     f(llc(1)  , llc(2)  , llc(3)+1 , llc(4)+1 )   * (1-xin) * (1-yin) *   zin   *   qin  + &
     f(llc(1)  , llc(2)+1, llc(3)+1 , llc(4)+1 )   * (1-xin) *  yin    *   zin   *   qin  + &
     f(llc(1)+1, llc(2)  , llc(3)+1 , llc(4)+1 )   *   xin   * (1-yin) *   zin   *   qin  + &
     f(llc(1)+1, llc(2)+1, llc(3)   , llc(4)+1 )   *   xin   *  yin    * (1-zin) *   qin  + &
     f(llc(1)+1, llc(2)+1, llc(3)+1 , llc(4)+1 )   *   xin   *  yin    *   zin   *   qin 

     outval = ((interp_QL**(1-sigma))/(1-sigma))

end subroutine

end module 