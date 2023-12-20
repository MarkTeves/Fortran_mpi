module bspline_sub_module
    contains 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine db1ink(x,nx,fcn,kx,iknot,tx,bcoef,iflag)

    implicit none

    integer,intent(in)                  :: nx     !! Number of \(x\) abcissae
    integer,intent(in)                  :: kx     !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    real,dimension(:),intent(in)        :: x      !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real,dimension(:),intent(in)        :: fcn    !! `(nx)` array of function values to interpolate. `fcn(i)` should
                                                      !! contain the function value at the point `x(i)`
    integer,intent(in)                  :: iknot  !! knot sequence flag:
                                                      !!
                                                      !! * 0 = knot sequence chosen by [[db1ink]].
                                                      !! * 1 = knot sequence chosen by user.
    real,dimension(:),intent(inout)     :: tx     !! The `(nx+kx)` knots in the \(x\) direction
                                                      !! for the spline interpolant:
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db1ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real,dimension(:),intent(out)       :: bcoef  !! `(nx)` array of coefficients of the b-spline interpolant.
    integer,intent(out)                 :: iflag  !! status flag:
                                                      !!
                                                      !! * 0 = successful execution.
                                                      !! * 2 = `iknot` out of range.
                                                      !! * 3 = `nx` out of range.
                                                      !! * 4 = `kx` out of range.
                                                      !! * 5 = `x` not strictly increasing.
                                                      !! * 6 = `tx` not non-decreasing.
                                                      !! * 700 = `size(x)` \( \ne \) `size(fcn,1)`.
                                                      !! * 706 = `size(x)` \( \ne \) `nx`.
                                                      !! * 712 = `size(tx)` \( \ne \) `nx+kx`.
                                                      !! * 800 = `size(x)` \( \ne \) `size(bcoef,1)`.

    logical :: status_ok
    real,dimension(:),allocatable :: work   !! work array of dimension `2*kx*(nx+1)`

    !check validity of inputs

    !call check_inputs(  iknot,&
    !                    iflag,&
    !                    nx=nx,&
    !                    kx=kx,&
    !                    x=x,&
    !                    f1=fcn,&
    !                    bcoef1=bcoef,&
    !                    tx=tx,&
    !                    status_ok=status_ok)

    !if (status_ok) then

        !choose knots
        if (iknot == 0) then
            call dbknot(x,nx,kx,tx)
        end if

        allocate(work(2*kx*(nx+1)))

        !construct b-spline coefficients
        call dbtpcf(x,nx,fcn,nx,1,tx,kx,bcoef,work,iflag)

        deallocate(work)

    !end if

    end subroutine db1ink
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine db1val(xval,idx,tx,nx,kx,bcoef,f,iflag,inbvx,w0,extrap)

    implicit none

    integer,intent(in)               :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)               :: nx       !! the number of interpolation points in \(x\).
                                                     !! (same as in last call to [[db1ink]])
    integer,intent(in)               :: kx       !! order of polynomial pieces in \(x\).
                                                     !! (same as in last call to [[db1ink]])
    real,intent(in)                  :: xval     !! \(x\) coordinate of evaluation point.
    real,dimension(nx+kx),intent(in) :: tx       !! sequence of knots defining the piecewise polynomial
                                                     !! in the \(x\) direction. (same as in last call to [[db1ink]])
    real,dimension(nx),intent(in)    :: bcoef    !! the b-spline coefficients computed by [[db1ink]].
    real,intent(out)                 :: f        !! interpolated value
    integer,intent(out)              :: iflag    !! status flag:
                                                     !!
                                                     !! * \( = 0 \)   : no errors
                                                     !! * \( \ne 0 \) : error
    integer,intent(inout)            :: inbvx    !! initialization parameter which must be set
                                                     !! to 1 the first time this routine is called,
                                                     !! and must not be changed by the user.
    real,dimension(3*kx),intent(inout) :: w0  !! work array
    logical,intent(in),optional          :: extrap   !! if extrapolation is allowed
                                                     !! (if not present, default is False)

    f = 0.0

    iflag = check_value(xval,tx,1,extrap); if (iflag/=0) return

    call dbvalu(tx,bcoef,nx,kx,idx,xval,inbvx,w0,iflag,f,extrap)

    end subroutine db1val
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dbknot(x,n,k,t)

    implicit none

    integer,intent(in)             :: n  !! dimension of `x`
    integer,intent(in)             :: k
    real,dimension(:),intent(in)   :: x
    real,dimension(:),intent(out)  :: t

    integer :: i, j, ipj, npj, ip1, jstrt
    real :: rnot

    !put k knots at each endpoint
    !(shift right endpoints slightly -- see pg 350 of reference)
    rnot = x(n) + 0.1*( x(n)-x(n-1) )
    do j=1,k
        t(j)   = x(1)
        npj    = n + j
        t(npj) = rnot
    end do

    !distribute remaining knots

    if (mod(k,2) == 1)  then

        !case of odd k --  knots between data points

        i = (k-1)/2 - k
        ip1 = i + 1
        jstrt = k + 1
        do j=jstrt,n
            ipj = i + j
            t(j) = 0.5*( x(ipj) + x(ipj+1) )
        end do

    else

        !case of even k --  knots at data points

        i = (k/2) - k
        jstrt = k+1
        do j=jstrt,n
            ipj = i + j
            t(j) = x(ipj)
        end do

    end if

    end subroutine dbknot
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dbvalu(t,a,n,k,ideriv,x,inbv,work,iflag,val,extrap)

    implicit none

    real,intent(out)             :: val     !! the interpolated value
    integer,intent(in)           :: n       !! number of b-spline coefficients.
                                                !! (sum of knot multiplicities-`k`)
    real,dimension(:),intent(in) :: t       !! knot vector of length `n+k`
    real,dimension(n),intent(in) :: a       !! b-spline coefficient vector of length `n`
    integer,intent(in)           :: k       !! order of the b-spline, `k >= 1`
    integer,intent(in)           :: ideriv  !! order of the derivative, `0 <= ideriv <= k-1`.
                                                !! `ideriv = 0` returns the b-spline value
    real,intent(in)              :: x       !! argument, `t(k) <= x <= t(n+1)`
    integer,intent(inout)        :: inbv    !! an initialization parameter which must be set
                                                !! to 1 the first time [[dbvalu]] is called.
                                                !! `inbv` contains information for efficient processing
                                                !! after the initial call and `inbv` must not
                                                !! be changed by the user.  distinct splines require
                                                !! distinct `inbv` parameters.
    real,dimension(:),intent(inout) :: work !! work vector of length at least `3*k`
    integer,intent(out)          :: iflag   !! status flag:
                                                !!
                                                !! * 0: no errors
                                                !! * 401: `k` does not satisfy `k` \( \ge \) 1
                                                !! * 402: `n` does not satisfy `n` \( \ge \) `k`
                                                !! * 403: `ideriv` does not satisfy 0 \( \le \) `ideriv` \(<\) `k`
                                                !! * 404: `x` is not greater than or equal to `t(k)`
                                                !! * 405: `x` is not less than or equal to `t(n+1)`
                                                !! * 406: a left limiting value cannot be obtained at `t(k)`
    logical,intent(in),optional :: extrap   !! if extrapolation is allowed
                                            !! (if not present, default is False)

    integer :: i,iderp1,ihi,ihmkmj,ilo,imk,imkpj,ipj,&
               ip1,ip1mj,j,jj,j1,j2,kmider,kmj,km1,kpk,mflag
    real :: fkmj
    real :: xt
    logical :: extrapolation_allowed  !! if extrapolation is allowed

    val = 0.0

    if (k<1) then
        iflag = 401  ! dbvalu - k does not satisfy k>=1
        return
    end if

    if (n<k) then
        iflag = 402  ! dbvalu - n does not satisfy n>=k
        return
    end if

    if (ideriv<0 .or. ideriv>=k) then
        iflag = 403  ! dbvalu - ideriv does not satisfy 0<=ideriv<k
        return
    end if

    if (present(extrap)) then
        extrapolation_allowed = extrap
    else
        extrapolation_allowed = .false.
    end if

    ! make a temp copy of x (for computing the
    ! interval) in case extrapolation is allowed
    if (extrapolation_allowed) then
        if (x<t(k)) then
            xt = t(k)
        else if (x>t(n+1)) then
            xt = t(n+1)
        else
            xt = x
        end if
    else
        xt = x
    end if

    kmider = k - ideriv

    ! find *i* in (k,n) such that t(i) <= x < t(i+1)
    ! (or, <= t(i+1) if t(i) < t(i+1) = t(n+1)).

    km1 = k - 1
    call dintrv(t, n+1, xt, inbv, i, mflag)
    if (xt<t(k)) then
        iflag = 404  ! dbvalu - x is not greater than or equal to t(k)
        return
    end if

    if (mflag/=0) then

        if (xt>t(i)) then
            iflag = 405  ! dbvalu - x is not less than or equal to t(n+1)
            return
        end if

        do
            if (i==k) then
                iflag = 406  ! dbvalu - a left limiting value cannot be obtained at t(k)
                return
            end if
            i = i - 1
            if (xt/=t(i)) exit
        end do

    end if

    ! difference the coefficients *ideriv* times
    ! work(i) = aj(i), work(k+i) = dp(i), work(k+k+i) = dm(i), i=1.k

    imk = i - k
    do j=1,k
        imkpj = imk + j
        work(j) = a(imkpj)
    end do

    if (ideriv/=0) then
        do j=1,ideriv
            kmj = k - j
            fkmj = real(kmj)
            do jj=1,kmj
                ihi = i + jj
                ihmkmj = ihi - kmj
                work(jj) = (work(jj+1)-work(jj))/(t(ihi)-t(ihmkmj))*fkmj
            end do
        end do
    end if

    ! compute value at *x* in (t(i),(t(i+1)) of ideriv-th derivative,
    ! given its relevant b-spline coeff. in aj(1),...,aj(k-ideriv).

    if (ideriv/=km1) then
        ip1 = i + 1
        kpk = k + k
        j1 = k + 1
        j2 = kpk + 1
        do j=1,kmider
            ipj = i + j
            work(j1) = t(ipj) - x
            ip1mj = ip1 - j
            work(j2) = x - t(ip1mj)
            j1 = j1 + 1
            j2 = j2 + 1
        end do
        iderp1 = ideriv + 1
        do j=iderp1,km1
            kmj = k - j
            ilo = kmj
            do jj=1,kmj
                work(jj) = (work(jj+1)*work(kpk+ilo)+work(jj)*&
                            work(k+jj))/(work(kpk+ilo)+work(k+jj))
                ilo = ilo - 1
            end do
        end do
    end if

    iflag = 0
    val = work(1)

    end subroutine dbvalu
    
    
    function check_value(x,t,i,extrap) result(iflag)

    implicit none

    integer                      :: iflag   !! returns 0 if value is OK, otherwise returns `600+i`
    real,intent(in)              :: x       !! the value to check
    integer,intent(in)           :: i       !! 1=x, 2=y, 3=z, 4=q, 5=r, 6=s
    real,dimension(:),intent(in) :: t       !! the knot vector
    logical,intent(in),optional      :: extrap  !! if extrapolation is allowed
                                                !! (if not present, default is False)

    logical :: allow_extrapolation  !! if extrapolation is allowed

    if (present(extrap)) then
        allow_extrapolation = extrap
    else
        allow_extrapolation = .false.
    end if

    if (allow_extrapolation) then
        ! in this case all values are OK
        iflag = 0
    else
        if (x<t(1) .or. x>t(size(t))) then
            iflag = 600 + i  ! value out of bounds (601, 602, etc.)
        else
            iflag = 0
        end if
    end if

    end function check_value
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dbtpcf(x,n,fcn,ldf,nf,t,k,bcoef,work,iflag)

    integer,intent(in)                :: n  !! dimension of `x`
    integer,intent(in)                :: nf
    integer,intent(in)                :: ldf
    integer,intent(in)                :: k
    real,dimension(:),intent(in)      :: x
    real,dimension(ldf,nf),intent(in) :: fcn
    real,dimension(:),intent(in)      :: t
    real,dimension(nf,n),intent(out)  :: bcoef
    real,dimension(*),intent(out)     :: work   !! work array of size >= `2*k*(n+1)`
    integer,intent(out)               :: iflag  !!   0: no errors
                                                    !! 301: n should be >0

    integer :: i, j, m1, m2, iq, iw

    ! check for null input

    if (nf > 0)  then

        ! partition work array
        m1 = k - 1
        m2 = m1 + k
        iq = 1 + n
        iw = iq + m2*n+1

        ! compute b-spline coefficients

        ! first data set

        call dbintk(x,fcn,t,n,k,work,work(iq),work(iw),iflag)
        if (iflag == 0) then
            do i=1,n
                bcoef(1,i) = work(i)
            end do

            !  all remaining data sets by back-substitution

            if (nf == 1)  return
            do j=2,nf
                do i=1,n
                    work(i) = fcn(i,j)
                end do
                call dbnslv(work(iq),m2,n,m1,m1,work)
                do i=1,n
                    bcoef(j,i) = work(i)
                end do
            end do
        end if

    else
        !write(error_unit,'(A)') 'dbtpcf - n should be >0'
        iflag = 301
    end if

    end subroutine dbtpcf
    
    subroutine dintrv(xt,lxt,xx,ilo,ileft,mflag,extrap)

    implicit none

    integer,intent(in)             :: lxt    !! length of the `xt` vector
    real,dimension(:),intent(in)   :: xt     !! a knot or break point vector of length `lxt`
    real,intent(in)                :: xx     !! argument
    integer,intent(inout)          :: ilo    !! an initialization parameter which must be set
                                                 !! to 1 the first time the spline array `xt` is
                                                 !! processed by dintrv. `ilo` contains information for
                                                 !! efficient processing after the initial call and `ilo`
                                                 !! must not be changed by the user.  distinct splines
                                                 !! require distinct `ilo` parameters.
    integer,intent(out)            :: ileft  !! largest integer satisfying `xt(ileft)` \( \le \) `x`
    integer,intent(out)            :: mflag  !! signals when `x` lies out of bounds
    logical,intent(in),optional        :: extrap !! if extrapolation is allowed
                                                 !! (if not present, default is False)

    integer :: ihi, istep, middle
    real :: x

    x = get_temp_x_for_extrap(xx,xt(1),xt(lxt),extrap)

    ihi = ilo + 1
    if ( ihi>=lxt ) then
        if ( x>=xt(lxt) ) then
            mflag = -2
            ileft = lxt
            return
        end if
        if ( lxt<=1 ) then
            mflag = -1
            ileft = 1
            return
        end if
        ilo = lxt - 1
        ihi = lxt
    end if

    if ( x>=xt(ihi) ) then

        ! now x >= xt(ilo). find upper bound
        istep = 1
        do
            ilo = ihi
            ihi = ilo + istep
            if ( ihi>=lxt ) then
                if ( x>=xt(lxt) ) then
                    mflag = -2
                    ileft = lxt
                    return
                end if
                ihi = lxt
            else if ( x>=xt(ihi) ) then
                istep = istep*2
                cycle
            end if
            exit
        end do

    else

        if ( x>=xt(ilo) ) then
            mflag = 0
            ileft = ilo
            return
        end if
        ! now x <= xt(ihi). find lower bound
        istep = 1
        do
            ihi = ilo
            ilo = ihi - istep
            if ( ilo<=1 ) then
                ilo = 1
                if ( x<xt(1) ) then
                    mflag = -1
                    ileft = 1
                    return
                end if
            else if ( x<xt(ilo) ) then
                istep = istep*2
                cycle
            end if
            exit
        end do

    end if

    ! now xt(ilo) <= x < xt(ihi). narrow the interval
    do
        middle = (ilo+ihi)/2
        if ( middle==ilo ) then
            mflag = 0
            ileft = ilo
            return
        end if
        ! note. it is assumed that middle = ilo in case ihi = ilo+1
        if ( x<xt(middle) ) then
            ihi = middle
        else
            ilo = middle
        end if
    end do

    end subroutine dintrv
    
    function get_temp_x_for_extrap(x,tmin,tmax,extrap) result(xt)

    implicit none

    real,intent(in) :: x    !! variable value
    real,intent(in) :: tmin !! first knot vector element for b-splines
    real,intent(in) :: tmax !! last knot vector element for b-splines
    real            :: xt   !! The value returned (it will either
                                !! be `tmin`, `x`, or `tmax`)
    logical,intent(in),optional :: extrap  !! if extrapolation is allowed
                                           !! (if not present, default is False)

    logical :: extrapolation_allowed  !! if extrapolation is allowed

    if (present(extrap)) then
        extrapolation_allowed = extrap
    else
        extrapolation_allowed = .false.
    end if

    if (extrapolation_allowed) then
        if (x<tmin) then
            xt = tmin
        else if (x>tmax) then
            ! Put it just inside the upper bound.
            ! This is sort of a hack to get
            ! extrapolation to work.
            xt = tmax - 2.0*spacing(tmax)
        else
            xt = x
        end if
    else
        xt = x
    end if

    end function get_temp_x_for_extrap
    
    subroutine dbintk(x,y,t,n,k,bcoef,q,work,iflag)

    implicit none

    integer,intent(in)            :: n      !! number of data points, n >= k
    real,dimension(n),intent(in)  :: x      !! vector of length n containing data point abscissa
                                                !! in strictly increasing order.
    real,dimension(n),intent(in)  :: y      !! corresponding vector of length n containing data
                                                !! point ordinates.
    real,dimension(*),intent(in)  :: t      !! knot vector of length n+k
                                                !! since t(1),..,t(k) <= x(1) and t(n+1),..,t(n+k)
                                                !! >= x(n), this leaves only n-k knots (not
                                                !! necessarily x(i) values) interior to (x(1),x(n))
    integer,intent(in)            :: k      !! order of the spline, k >= 1
    real,dimension(n),intent(out) :: bcoef  !! a vector of length n containing the b-spline coefficients
    real,dimension(*),intent(out) :: q      !! a work vector of length (2*k-1)*n, containing
                                                !! the triangular factorization of the coefficient
                                                !! matrix of the linear system being solved.  the
                                                !! coefficients for the interpolant of an
                                                !! additional data set (x(i),yy(i)), i=1,...,n
                                                !! with the same abscissa can be obtained by loading
                                                !! yy into bcoef and then executing
                                                !! call dbnslv(q,2k-1,n,k-1,k-1,bcoef)
    real,dimension(*),intent(out) :: work   !! work vector of length 2*k
    integer,intent(out)           :: iflag  !! *   0: no errors.
                                                !! * 100: k does not satisfy k>=1.
                                                !! * 101: n does not satisfy n>=k.
                                                !! * 102: x(i) does not satisfy x(i)<x(i+1) for some i.
                                                !! * 103: some abscissa was not in the support of the.
                                                !! corresponding basis function and the system is singular.
                                                !! * 104: the system of solver detects a singular system.
                                                !! although the theoretical conditions for a solution were satisfied.

    integer :: iwork, i, ilp1mx, j, jj, km1, kpkm2, left,lenq, np1
    real :: xi
    logical :: found

    if (k<1) then
        !write(error_unit,'(A)') 'dbintk - k does not satisfy k>=1'
        iflag = 100
        return
    end if

    if (n<k) then
        !write(error_unit,'(A)') 'dbintk - n does not satisfy n>=k'
        iflag = 101
        return
    end if

    jj = n - 1
    if (jj/=0) then
        do i=1,jj
            if (x(i)>=x(i+1)) then
                !write(error_unit,'(A)') 'dbintk - x(i) does not satisfy x(i)<x(i+1) for some i'
                iflag = 102
                return
            end if
        end do
    end if

    np1 = n + 1
    km1 = k - 1
    kpkm2 = 2*km1
    left = k
    ! zero out all entries of q
    lenq = n*(k+km1)
    do i=1,lenq
        q(i) = 0.0
    end do

    ! loop over i to construct the n interpolation equations
    do i=1,n

        xi = x(i)
        ilp1mx = min(i+k,np1)
        ! find left in the closed interval (i,i+k-1) such that
        !         t(left) <= x(i) < t(left+1)
        ! matrix is singular if this is not possible
        left = max(left,i)
        if (xi<t(left)) then
            !write(error_unit,'(A)') 'dbintk - some abscissa was not in the support of the'//&
            !             ' corresponding basis function and the system is singular'
            iflag = 103
            return
        end if
        found = .false.
        do
            found = (xi<t(left+1))
            if (found) exit
            left = left + 1
            if (left>=ilp1mx) exit
        end do
        if (.not. found) then
            left = left - 1
            if (xi>t(left+1)) then
                !write(error_unit,'(A)') 'dbintk - some abscissa was not in the support of the'//&
                !             ' corresponding basis function and the system is singular'
                iflag = 103
                return
            end if
        end if
        ! the i-th equation enforces interpolation at xi, hence
        ! a(i,j) = b(j,k,t)(xi), all j. only the  k  entries with  j =
        ! left-k+1,...,left actually might be nonzero. these  k  numbers
        ! are returned, in  bcoef (used for temp.storage here), by the
        ! following
        call dbspvn(t, k, k, 1, xi, left, bcoef, work, iwork, iflag)
        if (iflag/=0) return

        ! we therefore want  bcoef(j) = b(left-k+j)(xi) to go into
        ! a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since
        ! a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q
        ! as a two-dim. array , with  2*k-1  rows (see comments in
        ! dbnfac). in the present program, we treat  q  as an equivalent
        ! one-dimensional array (because of fortran restrictions on
        ! dimension statements) . we therefore want  bcoef(j) to go into
        ! entry
        !     i -(left+j) + 2*k + ((left+j) - k-1)*(2*k-1)
        !            = i-left+1 + (left -k)*(2*k-1) + (2*k-2)*j
        ! of q.
        jj = i - left + 1 + (left-k)*(k+km1)
        do j=1,k
            jj = jj + kpkm2
            q(jj) = bcoef(j)
        end do

    end do

    ! obtain factorization of a, stored again in q.
    call dbnfac(q, k+km1, n, km1, km1, iflag)

    if (iflag==1) then !success
        ! solve  a*bcoef = y  by backsubstitution
        do i=1,n
            bcoef(i) = y(i)
        end do
        call dbnslv(q, k+km1, n, km1, km1, bcoef)
        iflag = 0
    else  !failure
        !write(error_unit,'(A)') 'dbintk - the system of solver detects a singular system'//&
        !             ' although the theoretical conditions for a solution were satisfied'
        iflag = 104
    end if

    end subroutine
    
    
    subroutine dbnslv(w,nroww,nrow,nbandl,nbandu,b)

    integer,intent(in) :: nroww   !! describes the lu-factorization of a banded matrix a of order `nrow`
                                      !! as constructed in [[dbnfac]].
    integer,intent(in) :: nrow    !! describes the lu-factorization of a banded matrix a of order `nrow`
                                      !! as constructed in [[dbnfac]].
    integer,intent(in) :: nbandl  !! describes the lu-factorization of a banded matrix a of order `nrow`
                                      !! as constructed in [[dbnfac]].
    integer,intent(in) :: nbandu  !! describes the lu-factorization of a banded matrix a of order `nrow`
                                      !! as constructed in [[dbnfac]].
    real,dimension(nroww,nrow),intent(in) :: w !! describes the lu-factorization of a banded matrix a of
                                                   !! order `nrow` as constructed in [[dbnfac]].
    real,dimension(nrow),intent(inout) :: b  !! * **in**: right side of the system to be solved
                                                 !! * **out**: the solution x, of order nrow

    integer :: i, j, jmax, middle, nrowm1

    middle = nbandu + 1
    if (nrow/=1) then

        nrowm1 = nrow - 1
        if (nbandl/=0) then

            ! forward pass
            ! for i=1,2,...,nrow-1, subtract right side(i)*(i-th column of l)
            !                       from right side (below i-th row).
            do i=1,nrowm1
                jmax = min(nbandl,nrow-i)
                do j=1,jmax
                    b(i+j) = b(i+j) - b(i)*w(middle+j,i)
                end do
            end do

        end if

        ! backward pass
        ! for i=nrow,nrow-1,...,1, divide right side(i) by i-th diagonal
        !                          entry of u, then subtract right side(i)*(i-th column
        !                          of u) from right side (above i-th row).
        if (nbandu<=0) then
            ! a is lower triangular.
            do i=1,nrow
                b(i) = b(i)/w(1,i)
            end do
            return
        end if

        i = nrow
        do
            b(i) = b(i)/w(middle,i)
            jmax = min(nbandu,i-1)
            do j=1,jmax
                b(i-j) = b(i-j) - b(i)*w(middle-j,i)
            end do
            i = i - 1
            if (i<=1) exit
        end do

    end if

    b(1) = b(1)/w(middle,1)

    end subroutine dbnslv
    
    subroutine dbspvn(t,jhigh,k,index,x,ileft,vnikx,work,iwork,iflag)

    implicit none

    real,dimension(*),intent(in)  :: t        !! knot vector of length `n+k`, where
                                                  !! `n` = number of b-spline basis functions
                                                  !! `n` = sum of knot multiplicities-`k`
                                                  !! dimension `t(ileft+jhigh)`
    integer,intent(in)            :: jhigh    !! order of b-spline, `1 <= jhigh <= k`
    integer,intent(in)            :: k        !! highest possible order
    integer,intent(in)            :: index    !! index = 1 gives basis functions of order `jhigh`
                                                  !!       = 2 denotes previous entry with `work`, `iwork`
                                                  !!         values saved for subsequent calls to
                                                  !!         dbspvn.
    real,intent(in)               :: x        !! argument of basis functions, `t(k) <= x <= t(n+1)`
    integer,intent(in)            :: ileft    !! largest integer such that `t(ileft) <= x < t(ileft+1)`
    real,dimension(k),intent(out) :: vnikx    !! vector of length `k` for spline values.
    real,dimension(*),intent(out) :: work     !! a work vector of length `2*k`
    integer,intent(out)           :: iwork    !! a work parameter.  both `work` and `iwork` contain
                                                  !! information necessary to continue for `index = 2`.
                                                  !! when `index = 1` exclusively, these are scratch
                                                  !! variables and can be used for other purposes.
    integer,intent(out)           :: iflag    !! *   0: no errors
                                                  !! * 201: `k` does not satisfy `k>=1`
                                                  !! * 202: `jhigh` does not satisfy `1<=jhigh<=k`
                                                  !! * 203: `index` is not 1 or 2
                                                  !! * 204: `x` does not satisfy `t(ileft)<=x<=t(ileft+1)`

    integer :: imjp1, ipj, jp1, jp1ml, l
    real :: vm, vmprev

    ! content of j, deltam, deltap is expected unchanged between calls.
    ! work(i) = deltap(i),
    ! work(k+i) = deltam(i), i = 1,k

    if (k<1) then
        !write(error_unit,'(A)') 'dbspvn - k does not satisfy k>=1'
        iflag = 201
        return
    end if
    if (jhigh>k .or. jhigh<1) then
        !write(error_unit,'(A)') 'dbspvn - jhigh does not satisfy 1<=jhigh<=k'
        iflag = 202
        return
    end if
    if (index<1 .or. index>2) then
        !write(error_unit,'(A)') 'dbspvn - index is not 1 or 2'
        iflag = 203
        return
    end if
    if (x<t(ileft) .or. x>t(ileft+1)) then
        !write(error_unit,'(A)') 'dbspvn - x does not satisfy t(ileft)<=x<=t(ileft+1)'
        iflag = 204
        return
    end if

    iflag = 0

    if (index==1) then
        iwork = 1
        vnikx(1) = 1.0
        if (iwork>=jhigh) return
    end if

    do
        ipj = ileft + iwork
        work(iwork) = t(ipj) - x
        imjp1 = ileft - iwork + 1
        work(k+iwork) = x - t(imjp1)
        vmprev = 0.0
        jp1 = iwork + 1
        do l=1,iwork
            jp1ml = jp1 - l
            vm = vnikx(l)/(work(l)+work(k+jp1ml))
            vnikx(l) = vm*work(l) + vmprev
            vmprev = vm*work(k+jp1ml)
        end do
        vnikx(jp1) = vmprev
        iwork = jp1
        if (iwork>=jhigh) exit
    end do

    end subroutine dbspvn
    
    
    subroutine dbnfac(w,nroww,nrow,nbandl,nbandu,iflag)

    integer,intent(in) :: nroww   !! row dimension of the work array w. must be >= nbandl + 1 + nbandu.
    integer,intent(in) :: nrow    !! matrix order
    integer,intent(in) :: nbandl  !! number of bands of a below the main diagonal
    integer,intent(in) :: nbandu  !! number of bands of a above the main diagonal
    integer,intent(out) :: iflag  !! indicating success(=1) or failure (=2)
    real,dimension(nroww,nrow),intent(inout) :: w  !! work array. See header for details.

    integer :: i, ipk, j, jmax, k, kmax, middle, midmk, nrowm1
    real :: factor, pivot

    iflag = 1
    middle = nbandu + 1   ! w(middle,.) contains the main diagonal of a.
    nrowm1 = nrow - 1

    if (nrowm1 < 0) then
        iflag = 2
        return
    else if (nrowm1 == 0) then
        if (w(middle,nrow)==0.0) iflag = 2
        return
    end if

    if (nbandl<=0) then
        ! a is upper triangular. check that diagonal is nonzero .
        do i=1,nrowm1
            if (w(middle,i)==0.0) then
                iflag = 2
                return
            end if
        end do
        if (w(middle,nrow)==0.0) iflag = 2
        return
    end if

    if (nbandu<=0) then
        ! a is lower triangular. check that diagonal is nonzero and
        ! divide each column by its diagonal.
        do i=1,nrowm1
            pivot = w(middle,i)
            if (pivot==0.0) then
                iflag = 2
                return
            end if
            jmax = min(nbandl,nrow-i)
            do j=1,jmax
                w(middle+j,i) = w(middle+j,i)/pivot
            end do
        end do
        return
    end if

    ! a is not just a triangular matrix. construct lu factorization
    do i=1,nrowm1
        ! w(middle,i)  is pivot for i-th step .
        pivot = w(middle,i)
        if (pivot==0.0) then
            iflag = 2
            return
        end if
        ! jmax is the number of (nonzero) entries in column i
        ! below the diagonal.
        jmax = min(nbandl,nrow-i)
        ! divide each entry in column i below diagonal by pivot.
        do j=1,jmax
            w(middle+j,i) = w(middle+j,i)/pivot
        end do
        ! kmax is the number of (nonzero) entries in row i to
        ! the right of the diagonal.
        kmax = min(nbandu,nrow-i)
        ! subtract a(i,i+k)*(i-th column) from (i+k)-th column
        ! (below row i).
        do k=1,kmax
            ipk = i + k
            midmk = middle - k
            factor = w(midmk,ipk)
            do j=1,jmax
                w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor
            end do
        end do
    end do

    ! check the last diagonal entry.
    if (w(middle,nrow)==0.0) iflag = 2

    end subroutine dbnfac
end module
    