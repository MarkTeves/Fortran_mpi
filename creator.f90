program creator

    use mpi
    use mpi_base
    ! use mpi_constants
    use mpi_sizeofs
    use bspline_sub_module
    use interpQL

    implicit none
    !Auxiliar stuff:
    character(8000)  :: filename = ' '
    character(8000) :: path
    real :: infinity=-1e+10
    real :: rnd_num

    !Clock stuff
    real :: startTime, midTime1, midTime2, stopTime, timeElapsed

    !Switches
    integer :: solve_nodebt = 1
    integer :: cluster = 1 

    integer :: ierror, rank, size, len
    character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: version
    !Lifetime
    integer :: t_death = 45 
    integer :: t_born = 25
    integer :: total_periods
    integer :: t
    integer :: age
    
    ! Preference parameters
    real :: gamma = 2
    real :: bqt = 10.0
    real :: delta = 0.98

    !Income 
    !real :: c = -2.170042+2.700381
    !real :: b1 =  0.16818 !Age coefficients
    !real :: b2 = -0.0323371/10
    !real :: b3 = 0.0019704/100
    real :: mu_income = 0.02   !!!!!! Still need to calibrate this one !!!!!!!!
    real :: start_inc = 21.360 !Start income of deterministic part on CGM

    real ::  sig_yp = 0.158 !permanent income
    real ::  rho_yp = 0.946 !persistence perm inc
    real ::  sig_yt = 0.145

    real, allocatable :: g(:,:) !Deterministic income component
    real, allocatable :: Tmyp(:,:) !Permanent income transition matrix 
    real, allocatable :: gyp(:,:) !Permanent income grid
    real, allocatable :: gyt(:,:) !Temporary shocks !(allow this to be correlated with bcycle)

    real :: income_tax = 0.2

    real, allocatable  :: gyt_levels(:,:)
    real, allocatable :: jensen_term(:,:)

    ! Grids for state variables
    real, allocatable :: gcash(:,:)
    real, allocatable :: g_debt_long(:,:,:,:)
    real, allocatable :: g_maturity(:,:) 
    real :: grid_curvature 

    !Debt limits and cost
    real :: DTI_lim = 4.5, DTI
    real :: DTI_init = 3.5
    real :: LTV_lim = 0.9, LTV
    !real :: PTI_lim = 0.28, PTI 
    real, allocatable :: loan_premia(:,:)
    real :: default_penalty = 0.00

    real :: prepayment_penalty = 0.010, prepayment_penalty2
    !real :: prepayment_penalty = 0.000, prepayment_penalty2
    real :: origination_fee = 0.005, origination_fee2

    integer :: underwater_refi = 0 !set to 1 to allow underwater refi (with no eqt extraction). zero otherwise

    !Auxiliar debt variables
    real, allocatable :: annuity_factor(:,:)
    real :: interest_payment, principal_payment, payment_long
    real, allocatable :: tot_payment(:,:,:)
    real, allocatable :: long_debt_repay(:,:,:,:)
    integer :: debt_maturity, debt_maturity2

    !Grid limits
    real :: cash_max=200.0, cash_min=0.5
    real :: debt_long_min, debt_long_max


    !Grid sizes
    integer :: n_cash = 70.0
    integer :: n_yp = 5 !permanent income 
    integer :: n_rf = 3  
    integer :: n_bcycle = 2 !Normal times vs recession indicator (for asset value)
    integer :: n_long_debt = 10 !Number of points for debt
    integer :: n_maturity = 3

    !Rf calibration, transition matrix & grid
    real :: mu_rf = 0.012
    real :: sig_rf  = 0.018
    real :: rho_rf  = 0.825
    real :: term_premium
    real, allocatable :: grf(:,:)
    real, allocatable :: grf_long(:,:)
    real :: rf_long_unc  !Unconditional long risk free rate in logs
    real, allocatable :: Tmrf(:,:), Tmrf_steady(:,:), Tmrf_ss(:) !First matrix is transition matrix. Last two matrices are just auxiliars. 
    integer :: nsim = 200000 !Simulate to get expected long-term rate
    real, allocatable :: r_sim(:,:), r_sim2(:,:)
    real :: shock
    real :: rf_unc !Unconditional rf in logs for rental cost

    !Asset value calibration
    real :: mu_h = 0.016 !mu_h = 0.017
    real :: delta_h = 0.2
    real :: A0 
    real :: rc !rental cost
    real, allocatable :: Tmbcycle(:,:)
    real :: rental_premium = 0.02
    real :: prop_tax = 0.015, prop_maintenance = 0.025
    real :: A_value, delta_A !Auxiliar values
    real :: asset_sale_fee = 0.06
    real :: he, he_net !home equity asset equity
    real :: asset_purchase !Value of asset when purchased at t=1

    !Consumption choice
    real, allocatable :: gc(:,:)
    real ::  c_min = 0.01, c_max = 200
    integer :: n_c = 201

    !Quadrature
    integer :: nqp = 3 
    real, allocatable :: weig(:,:)
    real, allocatable :: grid(:,:)
    real, allocatable :: weig_yt(:,:) !Weight matrix for temporary shock

    !Invest 
    real, allocatable :: invest2(:,:), nw2(:,:)
    real, allocatable :: invest3(:,:,:), nw3(:,:,:)

    !Choice refinancing
    real, allocatable :: refinancing_choice(:,:,:)

    !Utility and auxiliars
    real, allocatable :: u(:,:)
    real, allocatable :: v2(:,:), v3(:,:,:)
    real :: auxv, auxv1
    integer, dimension(1) :: pt1
    integer, dimension(2) :: pt2
    integer, dimension(3) :: pt3
    real :: v_cashout, v_keep, v_default



    !Policy Functions
    ! real, allocatable :: c_nodebt(:,:,:,:,:), v_nodebt(:,:,:,:,:)
    ! real, allocatable :: c_debt(:,:,:,:,:,:,:,:,:), v_debt(:,:,:,:,:,:,:,:,:), default_pf(:,:,:,:,:,:,:,:,:), refinancing_pf(:,:,:,:,:,:,:,:,:)
    ! real, allocatable :: debt_long_pf(:,:,:,:,:,:,:,:,:), maturity_pf(:,:,:,:,:,:,:,:,:)


    !Counters
    integer :: indt, indcash, indrf, indyp, indbcycle !outer loop
    integer :: indyp2, indrf2, indbcycle2, indyt !inner loop
    integer :: indc
    integer :: indsim
    integer :: indlongdebt, indmaturity, indloan_rf, indloanper
    integer :: indlongdebt2, indmaturity2, indloan_rf2, indloanper2


        ! Calculate the number of iterations per process
    integer :: iter_per_process
    integer :: start_index, end_index

    !For tensor product splines
    integer :: Kxord = 3, Kyord = 3, Kzord = 3
    real, allocatable :: Xknot(:), Yknot(:), Zknot(:)
    real, allocatable :: bscoef4(:,:,:,:), bscoef8(:,:,:,:,:,:,:,:)
    real, allocatable :: w0(:)
    integer :: iflag
    integer :: inbvx = 1

    !Auxiliar variables
    real :: inc, cash
    real :: aux, sub
    real, allocatable :: aux_grid(:,:)
    integer :: scenario


    allocate(loan_premia(n_rf,1))


        
        n_maturity = 1
        allocate(g_maturity(n_maturity,1))
        g_maturity(1,1) =	30
        term_premium = 0.01/20
        dti_init 	=	3.5
            
        loan_premia(1,1) 	= -0.02
        loan_premia(2,1) 	= 0.0025	
        loan_premia(3,1) 	= 0.0400	
            
        prepayment_penalty 	=	0.01
        start_inc = 21.360

        path	=	'/output_path_here/'  

        


    !external BSINT, BS2IN, BSNAK
    !real BSVAL
    !real

    ! Initialize MPI    
    call MPI_INIT(ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
    call MPI_GET_LIBRARY_VERSION(version, len, ierror)

    startTime = MPI_WTIME()
    total_periods = t_death - t_born + 1
    iter_per_process = total_periods / size



  ! Calculate the starting and ending indices for each process
  
  start_index = rank * iter_per_process + 1
  end_index = (rank + 1) * iter_per_process

    !Quadrature
    allocate(weig(nqp,1))
    allocate(grid(nqp,1))

    if (nqp .eq. 2) then
        weig(1,1) = 0.5
        weig(2,1) = 0.5
        grid(1,1) = -1.0
        grid(2,1) = 1.0
    elseif (nqp .eq. 3) then
        weig(1,1)= 0.1666666666666
        weig(2,1)= 0.6666666666666
        weig(3,1)= 0.1666666666666
        grid(1,1)= -1.73205080756887
        grid(2,1)=  0.0
        grid(3,1)=  1.73205080756887
    else
        write(*,*) "error"

    end if


    !Create grids
    allocate(gcash(n_cash,1))
    allocate(gyp(n_yp, 1))
    allocate(Tmyp(n_yp,n_yp))
    allocate(grf(n_rf, 1))
    allocate(grf_long(n_rf, n_maturity))
    allocate(Tmrf(n_rf,n_rf))
    allocate(gc(n_c,1))
    allocate(Tmbcycle(n_bcycle, n_bcycle))

    grid_curvature = 0.7
    call linspace(cash_min**grid_curvature, cash_max**grid_curvature, n_cash, gcash)
    gcash = gcash**(1/grid_curvature)

    grid_curvature = 0.7
    call linspace(c_min**grid_curvature, c_max**grid_curvature, n_c, gc)
    gc = gc**(1/grid_curvature)

    call rouwenhorst(rho_rf,sig_rf,n_rf,grf,Tmrf) !Double checked with Kurt Mitman's Matlab rouwenhorst code that this is alright
    grf = grf + mu_rf
    call rouwenhorst(rho_yp,sig_yp,n_yp,gyp,Tmyp) 

    allocate(Tmrf_steady(n_rf, n_rf))
    allocate(Tmrf_ss(n_rf))
     Tmrf_steady = Tmrf
     do indsim = 1,10000
        Tmrf_steady = matmul(Tmrf,Tmrf_steady)
     end do
     Tmrf_ss = Tmrf_steady(1,:)
     rf_unc = sum(Tmrf_ss(:)*grf(:,1)) !Unconditional long-term interest rate to define mortgage annuity

    deallocate(Tmrf_steady, Tmrf_ss)

    !Need to calibrate this
    Tmbcycle(1,1) = 0.80! Rec | Rec
    Tmbcycle(1,2) = 0.20! Exp | Rec
    Tmbcycle(2,1) = 0.09! Rec | Exp
    Tmbcycle(2,2) = 0.91! Exp | Exp

    !Setup income deterministic part and idiosyncratic shocks
    allocate(g(total_periods,1)) 
    do indt = 1, total_periods
        age = indt + t_born - 1
        !g(indt,1) = exp(c + b1 * age + b2*age**2 + b3*age**3)
        g(indt,1) = start_inc*(1+mu_income)**(indt-1)
    end do

    ! Temporary shock 
    allocate(gyt(nqp,n_bcycle))
    allocate(weig_yt(nqp, n_bcycle))

    !Temporary shock in recession and expansion
    gyt(1,1) = sig_yt*1
    gyt(2,1) = sig_yt*(-1)
    gyt(3,1) = log(0.5) !50% drop on the worst realization of temporary income

    weig_yt(3,1) = 0.063 !Probability of large income drop in worst realization (recession)
    weig_yt(2,1) = (1-weig_yt(3,1))/2 !Probability of the other 2 states is identical
    weig_yt(1,1) = (1-weig_yt(3,1))/2

    gyt(1,2) = sig_yt*1
    gyt(2,2) = sig_yt*(-1)
    gyt(3,2) = log(0.5) !50% drop on the worst realization of temporary income

    weig_yt(3,2) = 0.055 !Probability of large income drop in worst realization (expansion)
    weig_yt(2,2) = (1-weig_yt(3,2))/2 !Probability of the other 2 states is identical
    weig_yt(1,2) = (1-weig_yt(3,2))/2

    !Make jensen adjustment
    allocate(jensen_term(1,n_bcycle))
    allocate(gyt_levels(nqp,n_bcycle))

    jensen_term(1,1) = sum(weig_yt(:,1)*exp(gyt(:,1)))
    jensen_term(1,2) = sum(weig_yt(:,2)*exp(gyt(:,2)))

    gyt_levels(:,1) = exp(gyt(:,1))/jensen_term(1,1)
    gyt_levels(:,2) = exp(gyt(:,2))/jensen_term(1,2)




    !A0 = DTI_lim*(g(1,1))-0.05 ! -1 just to make sure DTI constraint of 3.5 on the middle point of permanent income is not binding
    A0 = DTI_init*(g(1,1))-0.05 !

    !Allocate splines
    allocate(Xknot(n_cash+Kxord))
    !allocate(Yknot(n_debt+Kyord))
    !allocate(Zknot(+Kzord))
    allocate(bscoef4(n_cash,n_yp, n_rf, n_bcycle))

    Xknot = 0.0
    !Yknot = 0.0
    !Zknot = 0.0
    !call BSNAK (n_cash, gcash(:,1) , Kxord , Xknot)
    !call BSNAK (, ) , Kyord , Yknot)
    ! call BSNAK ( ,  , Kzord , Zknot)

    call dbknot(gcash(:,1) ,n_cash, Kxord ,  Xknot)

        !This code is prepared to use 3 types of splines;

    !(1) Type 1: Standard spline/splint cubic splines from CGM. In practice these work pretty poorly specially on the "no debt" problem when the agent has no cash to pay the rent
    !Further it cannot cannot be extended to multiple dimensions easily             
    !(2) Type 2: Tensor product splines. These splines work pretty well, and can accomodate multiple dimensions, but are not very fast. 
    !(3) Type 3: Quasi linear interpolation. This is a trick that I learned from Taha Choukmane at MIT. First transform the value function into a linear function; 
        !        Do a linear interpolation and then transform the problem back. This works pretty well for a CRRA utility function.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solving no debt problem !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Allocate value/policy functions
    allocate(c_nodebt(n_cash, n_yp, n_rf, n_bcycle, total_periods))
    allocate(v_nodebt(n_cash, n_yp, n_rf, n_bcycle, total_periods))
    allocate(invest2(n_c,1))
    allocate(nw2(n_c,1))
    allocate(u(n_c,1))
    allocate(v2(n_c,1))
    allocate(w0(3*Kxord))

    if (solve_nodebt .eq. 1) then
        !Terminal period
        t = total_periods
        do indyp = 1, n_yp
            do indrf = 1, n_rf
                do indbcycle = 1, n_bcycle    
                    v_nodebt(:,indyp,indrf,indbcycle,t) = bqt*(gcash(:,1)**(1-gamma))/(1-gamma)
                    c_nodebt(:,indyp,indrf,indbcycle,t) = gcash(:,1)
                end do
            end do
        end do

        u = (gc**(1-gamma))/(1-gamma)
        !Other periods
        do indt = 1, total_periods - 1
            t = total_periods - indt
            write(*,*) t
            do indyp2 = 1, n_yp
                do indrf2 = 1, n_rf
                    do indbcycle2 = 1, n_bcycle


                        !(1)
                        !call  spline(gcash, v_nodebt(:,indyp2,indrf2,indbcycle2,t+1),n_cash,gamma,bscoef4(:,indyp2,indrf2,indbcycle2))
                        !(2)            
                        !call db1ink(gcash(:,1),n_cash,v_nodebt(:,indyp2,indrf2,indbcycle2,t+1),Kxord,0,Xknot,bscoef4(:,indyp2,indrf2,indbcycle2),iflag)
                        !(3)
                        bscoef4(:,indyp2,indrf2,indbcycle2) = (((1-gamma)*v_nodebt(:,indyp2,indrf2,indbcycle2,t+1))**(1/(1-gamma)))
                
                    end do
                end do
            end do
            do indyp = 1, n_yp
                     
                !(1)
                !call splint(gcash,v_nodebt(:,indyp2,indrf2,indbcycle2,t+1),bscoef4(:,indyp2,indrf2,indbcycle2),n_cash,inc,auxv)
                
                !(2)
                !inbvx = 1
                !call db1val(cash,0,Xknot,n_cash,Kxord,bscoef4(:,indyp2,indrf2,indbcycle2),auxv,iflag,inbvx,w0)
                
                !(3)
                call interp1D_QL(gamma, n_cash, gcash(:,1), bscoef4(:,indyp2,indrf2,indrf1,indrf0))    !write(*,'(a,i0,a,i0)') "Time period:", t, " indyp:", indyp
                do indrf = 1, n_rf
                    do indbcycle = 1, n_bcycle
                            !$omp parallel do default(private) shared(indt,t,indyp, indrf, indbcycle, Tmbcycle, mu_h, delta_h, A0, grf, rental_premium,  &
                            !$omp& gc, gcash, n_c, n_yp, n_rf, n_bcycle, n_cash, total_periods, g, gyp, gyt, gyt_levels, nqp, Kxord, Xknot, bscoef4, Tmyp, weig, Tmrf, delta, infinity, &
                            !$omp& u, c_nodebt, v_nodebt, cash_min, cash_max,gamma, prop_maintenance, prop_tax, weig_yt, income_tax, rf_unc)

                        ! Perform the parallel computation using MPI
                        integer :: i    
                        do i = start_index, end_index
                            indcash = i
                            
                            !The line below is wrong, replaced with if statement
                            !delta_A = Tmbcycle(indbcycle,1)*(1+mu_h)*(1-delta_h) + Tmbcycle(indbcycle,2)*(1+mu_h)
                            if (indbcycle .eq. 1) then
                                delta_A = Tmbcycle(indbcycle,1)*(1+mu_h) + Tmbcycle(indbcycle,2)*(1+mu_h)*(1/(1-delta_h))
                            else 
                                delta_A = Tmbcycle(indbcycle,1)*(1+mu_h)*(1-delta_h) + Tmbcycle(indbcycle,2)*(1+mu_h)
                            end if        
                            
                            if (indbcycle .eq. 1) then
                                A_value = A0*((1+mu_h)**(t-1))*(1-delta_h)
                            else
                                A_value = A0*((1+mu_h)**(t-1))
                            end if
                            
                            
                            rc = ( (exp(grf(indrf,1))-1) - (delta_A - 1) + prop_tax + prop_maintenance + rental_premium)*A_value
                            !rc = ( (exp(rf_unc)-1) - (delta_A - 1) + prop_tax + prop_maintenance + rental_premium)*A_value
                            invest2(:,1) = gcash(indcash,1) - gc(:,1) - rc
                            
                            nw2(:,1) = invest2(:,1) * (exp(grf(indrf,1)))
                            
                            aux = maxval(nw2) !If agents have no money for rent, give the minimum possible subsidy so that they can aford it.
                            if (aux .le. gcash(1,1)) then
                                sub = -aux + gcash(1,1)
                                nw2(:,1) = nw2(:,1) + sub
                            end if
                            
                            
                            v2(:,:) = 0.0
                            
                            do indc = 1, n_c
                                
                                auxv1 = 0.0    
                                if (nw2(indc,1) .ge. 0.0) then 
                                    do indyp2 = 1, n_yp
                                        do indyt = 1, nqp    
                                            do indrf2 = 1,n_rf
                                                do indbcycle2 = 1, n_bcycle             
                                                    inc = g(t,1)*exp(gyp(indyp2,1))*gyt_levels(indyt,indbcycle2)
                                                    cash = nw2(indc,1) + inc*(1-income_tax)
                                                    cash = max(cash, cash_min) 
                                                    cash = min(cash, cash_max)
                                        indbcycle2), cash, auxv)
                                                    
                                                    auxv1 = auxv1+auxv*Tmyp(indyp,indyp2)*weig_yt(indyt,indbcycle2)*Tmrf(indrf,indrf2)*Tmbcycle(indbcycle, indbcycle2)
                                                end do 
                                            end do
                                        end do 
                                    end do
                                    v2(indc,1) = u(indc,1) + delta*auxv1
                                else
                                    v2(indc,1) = infinity
                                end if
                            end do            
                            
                            pt2 = maxloc(v2)
                            c_nodebt(indcash,indyp,indrf,indbcycle,t) = gc(pt2(1),1)
                            v_nodebt(indcash,indyp,indrf,indbcycle,t) = v2(pt2(1),1)
                            
                        end do
                        !$omp end parallel do
                        call MPI_BARRIER( MPI_COMM_WORLD, ierror)
                    end do
                end do
            end do
            
            midTime2 = MPI_WTIME()
            write(*,"(A,F15.1,A)") "Total time elapsed",midTime2-startTime," seconds." 
        end do
        
        midTime2 = MPI_WTIME()
        write(*,"(A,F15.1,A)") "Finished in",midTime2-startTime," seconds." 

        write(filename,'(a)')trim(path)//'c_nodebt.txt'
        open(unit=25,file=filename,status='replace',action='write',form='formatted')
        write(25,*) c_nodebt(:,:,:,:,:)
        close(unit=25)
        
        write(filename,'(a)')trim(path)//'v_nodebt.txt'
        open(unit=25,file=filename,status='replace',action='write',form='formatted')
        write(25,*) v_nodebt(:,:,:,:,:)
        close(unit=25)
    else
        write(filename,'(a)')trim(path)//'c_nodebt.txt'
        open(unit=25,file=filename, status='old', action='read')
        read(25,*) c_nodebt(:,:,:,:,:)
        close(unit=25)
        
        write(filename,'(a)')trim(path)//'v_nodebt.txt'
        open(unit=25,file=filename, status='old', action='read')
        read(25,*) v_nodebt(:,:,:,:,:)
        close(unit=25)
    end if     
    deallocate(w0)  



    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Solving debt problem !
    !!!!!!!!!!!!!!!!!!!!!!!!
    
    !Simulate to get expectation of long-term interest rate

       do indmaturity = 1,n_maturity

        debt_maturity = g_maturity(indmaturity,1)    
            
        allocate(r_sim(debt_maturity,nsim))
        allocate(r_sim2(nsim,1))


        do indrf = 1,n_rf
            do t = 1,debt_maturity
                do indsim = 1, nsim
                        
                        if (t .eq. 1) then
                            r_sim(t,indsim) = grf(indrf,1)
                        else 
                            call random_stdnormal(shock)
                            r_sim(t,indsim) = mu_rf * (1-rho_rf)+ rho_rf * r_sim(t-1,indsim) + sig_rf*shock
                        end if
                end do 
            end do
            r_sim2(:,1) = sum(r_sim, dim=1)/debt_maturity+term_premium*debt_maturity
            grf_long(indrf,indmaturity) = sum(r_sim2)/nsim
        end do

        deallocate(r_sim, r_sim2)
    end do

    debt_maturity = maxval(g_maturity)

    !Create grid for long term debt
    allocate(g_debt_long(n_maturity,debt_maturity,n_long_debt, n_rf)) 
    allocate(long_debt_repay(n_maturity,debt_maturity,n_long_debt, n_rf))
    allocate(tot_payment(n_maturity,n_long_debt, n_rf))
    allocate(annuity_factor(n_maturity,n_rf))

    annuity_factor = 0.0
    g_debt_long = 0.0
    long_debt_repay = 0.0
    tot_payment = 0.0

    debt_long_max = LTV_lim*A0*((1+mu_h)**(total_periods))
    debt_long_min = 0.0


    do indmaturity = 1,n_maturity

        debt_maturity = g_maturity(indmaturity,1)  

        do indrf = 1, n_rf
            call linspace(debt_long_min,debt_long_max,n_long_debt, g_debt_long(indmaturity,1,:,indrf))
        end do

        do indrf = 1,n_rf
            do indt = 1, debt_maturity 
                annuity_factor(indmaturity,indrf) = annuity_factor(indmaturity,indrf) + 1/((exp(grf_long(indrf,indmaturity))+loan_premia(indrf,1))**indt)
            end do
        end do

        do indrf = 1,n_rf
            do indlongdebt = 1,n_long_debt
                tot_payment(indmaturity,indlongdebt,indrf) = g_debt_long(indmaturity,1,indlongdebt,indrf)/annuity_factor(indmaturity,indrf)
            end do
        end do

        do indrf = 1,n_rf
            do indlongdebt = 1,n_long_debt
                do t = 2,debt_maturity
                    interest_payment = g_debt_long(indmaturity,t-1,indlongdebt,indrf)*(exp(grf_long(indrf,indmaturity))+loan_premia(indrf,1)-1)
                    principal_payment = (tot_payment(indmaturity,indlongdebt,indrf) - interest_payment)
                    g_debt_long(indmaturity,t,indlongdebt,indrf) = g_debt_long(indmaturity,t-1,indlongdebt,indrf) - principal_payment
                    long_debt_repay(indmaturity,t-1,indlongdebt,indrf) = principal_payment
                end do
                interest_payment = g_debt_long(indmaturity,debt_maturity,indlongdebt,indrf)*(exp(grf_long(indrf,indmaturity))+loan_premia(indrf,1)-1)
                principal_payment = (tot_payment(indmaturity,indlongdebt,indrf) - interest_payment)
                long_debt_repay(indmaturity,debt_maturity,indlongdebt,indrf) = principal_payment
            end do
        end do

    end do


    !Save some variables to open in matlab
    write(filename,'(a)')trim(path)//'g_debt_long'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) g_debt_long 
    close(unit=25)
    write(filename,'(a)')trim(path)//'g_maturity'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) g_maturity 
    close(unit=25)
    write(filename,'(a)')trim(path)//'grf_long'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) grf_long 
    close(unit=25)
    write(filename,'(a)')trim(path)//'long_debt_repay'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) long_debt_repay 
    close(unit=25)



    !Allocate value/policy functions
    debt_maturity = maxval(g_maturity)
    allocate(c_debt(n_cash, n_yp, n_rf, n_bcycle, n_long_debt, debt_maturity, n_maturity, n_rf, total_periods)) !First point on debt_maturity dimension is when loan is taken. 
    allocate(v_debt(n_cash, n_yp, n_rf, n_bcycle, n_long_debt, debt_maturity, n_maturity, n_rf, total_periods))
    allocate(default_pf(n_cash, n_yp, n_rf, n_bcycle, n_long_debt, debt_maturity, n_maturity, n_rf, total_periods))
    allocate(refinancing_pf(n_cash, n_yp, n_rf, n_bcycle, n_long_debt, debt_maturity, n_maturity, n_rf, total_periods))
    allocate(debt_long_pf(n_cash, n_yp, n_rf, n_bcycle, n_long_debt, debt_maturity, n_maturity, n_rf, total_periods))
    allocate(maturity_pf(n_cash, n_yp, n_rf, n_bcycle, n_long_debt, debt_maturity, n_maturity, n_rf, total_periods))

    allocate(invest3(n_c,n_long_debt,n_maturity))
    allocate(refinancing_choice(n_c,n_long_debt,n_maturity))
    allocate(nw3(n_c,n_long_debt,n_maturity))
    allocate(v3(n_c,n_long_debt,n_maturity))
    allocate(w0(3*Kxord))


    ! Read in policy functions (useful to debug)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !write(filename,'(a)')trim(path)//'c_debt.txt'
    !open(unit=25,file=filename,status='old',action='read')
    !read(25,*) c_debt(:,:,:,:,:,:,:,:,:)
    !close(unit=25)
    !
    !write(filename,'(a)')trim(path)//'v_debt.txt'
    !open(unit=25,file=filename,status='old',action='read')
    !read(25,*) v_debt(:,:,:,:,:,:,:,:,:)
    !close(unit=25)
    !
    !write(filename,'(a)')trim(path)//'default_pf.txt'
    !open(unit=25,file=filename,status='old',action='read')
    !read(25,*) default_pf(:,:,:,:,:,:,:,:,:)
    !close(unit=25)
    !
    !write(filename,'(a)')trim(path)//'debt_long_pf.txt'
    !open(unit=25,file=filename,status='old',action='read')
    !read(25,*) debt_long_pf(:,:,:,:,:,:,:,:,:)
    !close(unit=25)
    !
    !write(filename,'(a)')trim(path)//'maturity_pf.txt'
    !open(unit=25,file=filename,status='old',action='read')
    !read(25,*) maturity_pf(:,:,:,:,:,:,:,:,:)
    !close(unit=25)
    !
    !write(filename,'(a)')trim(path)//'refinancing_pf.txt'
    !open(unit=25,file=filename,status='old',action='read')
    !read(25,*) refinancing_pf(:,:,:,:,:,:,:,:,:)
    !close(unit=25)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !Terminal period
    t = total_periods

    do indyp = 1, n_yp
        do indrf = 1, n_rf
            do indbcycle = 1, n_bcycle  
                do indlongdebt = 1,n_long_debt
                    do indmaturity = 1,n_maturity
                        debt_maturity = g_maturity(indmaturity,1)
                        debt_maturity2 = min(debt_maturity,t)    
                        do indloanper = 1,debt_maturity2 
                            do indloan_rf = 1,n_rf
                                ! integer :: j    
                                do i = start_index, end_index
                                    indcash = i
                                    
                                    if (indbcycle .eq. 1) then
                                        A_value = A0*((1+mu_h)**(t-1))*(1-delta_h)
                                    else
                                        A_value = A0*((1+mu_h)**(t-1))
                                    end if
                                    
                                    interest_payment = (g_debt_long(indmaturity,indloanper,indlongdebt,indloan_rf))*(exp(grf_long(indloan_rf,indmaturity))+loan_premia(indloan_rf,1)-1)
                                    principal_payment = g_debt_long(indmaturity,indloanper,indlongdebt,indloan_rf)
                                    he =  A_value - principal_payment !Need to repay all debt (bullet part plus remaining principal of non-bullet part)
                                    if (indloanper .eq. debt_maturity) then 
                                        prepayment_penalty2 = 0.0
                                    else
                                        prepayment_penalty2 = prepayment_penalty*g_debt_long(indmaturity,indloanper,indlongdebt,indloan_rf)
                                    end if
                                    he_net = A_value*(1-asset_sale_fee) - principal_payment - prepayment_penalty2
                                    
                                    !Check if agent:
                                    !(1) Makes last payment on the loan (keep as is)
                                    cash = gcash(indcash,1) + he - interest_payment
                                    if (cash .le. 0.0) then
                                        v_keep = infinity
                                    else
                                        v_keep = bqt*(cash**(1-gamma))/(1-gamma)
                                    end if 
                                    
                                    !(2) Sells the asset (comission for sale) and moves to rental market (on the last period he will never do this since cashout on the last period is dominated by bequesting the home equity)
                                    if (he_net - interest_payment .ge. 0.0) then
                                        cash = gcash(indcash,1) + he_net - interest_payment
                                        pt1 = minloc(abs(gcash(:,1)-cash))
                                        v_cashout = v_nodebt(pt1(1),indyp,indrf,indbcycle,t)
                                    else
                                        v_cashout = infinity
                                    end if

                                    !(iii) Defaults
                                    v_default = v_nodebt(indcash,indyp,indrf,indbcycle,t) - default_penalty
                                    
                                    if ((v_keep .ge. v_default) .and. (v_keep .ge. v_cashout)) then 
                                        v_debt(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper, indmaturity, indloan_rf, total_periods) = v_keep
                                        default_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper, indmaturity, indloan_rf, total_periods) = 1.0
                                        c_debt(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper, indmaturity, indloan_rf, total_periods) = gcash(indcash,1)
                                    elseif ((v_cashout .ge. v_keep) .and. (v_cashout .ge. v_default)) then
                                        v_debt(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper, indmaturity, indloan_rf, total_periods) = v_cashout
                                        default_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper, indmaturity, indloan_rf, total_periods) = 2.0
                                    else 
                                        v_debt(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper, indmaturity, indloan_rf, total_periods) = v_default
                                        default_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper, indmaturity, indloan_rf, total_periods) = 3.0
                                        
                                    end if
                                    
                                end do
                                call MPI_BARRIER( MPI_COMM_WORLD, ierror)    
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do



    u = (gc**(1-gamma))/(1-gamma)

    debt_maturity = maxval(g_maturity)
    allocate(bscoef8(n_cash, n_yp, n_rf, n_bcycle, n_long_debt, debt_maturity, n_maturity, n_rf))

    !Other periods
    do indt = 1, total_periods - 1
        t = total_periods - indt
        write(*,*) t
        do indyp2 = 1, n_yp
            do indrf2 = 1, n_rf
                do indbcycle2 = 1, n_bcycle
                    do indlongdebt2 = 1,n_long_debt
                        do indmaturity2 = 1,n_maturity            
                            debt_maturity = g_maturity(indmaturity2,1)
                            debt_maturity2 = min(debt_maturity,t)
                            do indloanper2 = 1,debt_maturity2         
                                do indloan_rf2 = 1,n_rf
                                    
                                    !(1)
                                    !call spline(gcash(:,1),v_debt(:,indyp2,indrf2,indbcycle2, indlongdebt2, indloanper2, indmaturity2, indloan_rf2, t+1), n_cash, gamma,bscoef8(:,indyp2,indrf2,indbcycle2, indlongdebt2, indloanper2, indmaturity2, indloan_rf2))
                                    !(2)
                                    !call db1ink(gcash(:,1),n_cash,v_debt(:,indyp2,indrf2,indbcycle2, indlongdebt2, indloanper2, indmaturity2, indloan_rf2, t+1),Kxord,0,Xknot,bscoef8(:,indyp2,indrf2,indbcycle2, indlongdebt2, indloanper2, indmaturity2, indloan_rf2),iflag)
                                    !(3)
                                    bscoef8(:,indyp2,indrf2,indbcycle2, indlongdebt2, indloanper2, indmaturity2, indloan_rf2) = (((1-gamma)*v_debt(:,indyp2,indrf2,indbcycle2, indlongdebt2, indloanper2, indmaturity2, indloan_rf2, t+1))**(1/(1-gamma)))
                                    
                                    
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
        
        do indyp = 1, n_yp
        do indrf = 1, n_rf
        do indbcycle = 1, n_bcycle
        
        do indmaturity = 1,n_maturity
        debt_maturity = g_maturity(indmaturity,1)
        debt_maturity2 = min(debt_maturity,t)
        do indloanper = 1,debt_maturity2        
            midTime2 = MPI_WTIME()
            write(*,'(a,i0,a,i0,a,i0,a,i0)') "Time period:", t, " indyp:", indyp, " indrf:", indrf, " indbcycle:", indbcycle
            write(*,"(A,F15.1,A)") "Total time elapsed",midTime2-startTime," seconds." 
            !PTI_lim
            !$omp parallel default(private) shared(indt, total_periods, t, n_yp, n_rf, n_bcycle, n_long_debt, debt_maturity, debt_maturity2, n_maturity, gcash, n_cash, Kxord, Xknot, bscoef8, iflag,  &
            !$omp& indyp, indrf, indbcycle, indmaturity, indloanper,  &
            !$omp& A0, mu_h, delta_h, prepayment_penalty, origination_fee, g_debt_long, g_maturity, grf_long, loan_premia, long_debt_repay, &
            !$omp& gc, g, gyp, gyt, gyt_levels, LTV_lim, DTI_lim, grf, n_c, nqp, cash_min, cash_max, u, delta, &
            !$omp& v_debt, default_pf, debt_long_pf, maturity_pf, refinancing_pf, c_debt, &
            !$omp& asset_sale_fee, default_penalty, underwater_refi, &
            !$omp& Tmyp, weig, Tmrf, Tmbcycle, infinity, v_nodebt, tot_payment, gamma, prop_maintenance, prop_tax, weig_yt, income_tax)
            ! $omp do collapse(3)
        do indlongdebt = 1,n_long_debt
            do indloan_rf = 1,n_rf
                do i = start_index, end_index
                    indcash = i

                    
                    
                    if (indbcycle .eq. 1) then
                        A_value = A0*((1+mu_h)**(t-1))*(1-delta_h)
                    else
                        A_value = A0*((1+mu_h)**(t-1))
                    end if
                    
                    !if (t .eq. 1) then !Make agent purchase asset on first period
                    !    asset_purchase = A_value
                    !    !asset_purchase = 0.0
                    !else
                    !    asset_purchase = 0.0
                    !end if
                    !
                    
                    if (t > 1) then
                        
                        
                        ! no refinancing nor prepayment occuring (same initial loan amount, same maturity and loan is not at maturity)
                        if (indloanper .ne. debt_maturity) then                 
                            interest_payment = (g_debt_long(indmaturity,indloanper,indlongdebt,indloan_rf))*(exp(grf_long(indloan_rf,indmaturity))+loan_premia(indloan_rf,1)-1)
                            principal_payment = long_debt_repay(indmaturity,indloanper,indlongdebt,indloan_rf) 
                            invest2(:,1) = gcash(indcash,1) - gc(:,1) - principal_payment - interest_payment -  A_value * (prop_maintenance + prop_tax)
                        else !If loan is at maturity it has to be refinanced - even if to a "zero loan" below
                            invest2(:,1) = infinity
                        end if 
                        
                        do indlongdebt2 = 1,n_long_debt !This allows to refinance to the same loan amount
                            do indmaturity2 = 1,n_maturity           
                                
                                if ((indloanper .eq. debt_maturity) .or. (indlongdebt .eq. 1)) then ! (i) no leverage outstanding or previous loan reached maturity (no prepayment cost but origination fee);
                                    prepayment_penalty2 = 0.0
                                    origination_fee2 = origination_fee*g_debt_long(indmaturity2,1,indlongdebt2,indrf)
                                else !(ii) prepayment cost if changing loan amount/maturity and loan outstanding>0 or not at maturity AND origination fee.
                                    prepayment_penalty2 = prepayment_penalty*g_debt_long(indmaturity,indloanper,indlongdebt,indloan_rf)
                                    origination_fee2 = origination_fee*g_debt_long(indmaturity2,1,indlongdebt2,indrf)
                                end if

                                !Check restrictions first 
                                interest_payment = (g_debt_long(indmaturity2,1,indlongdebt2,indrf))*(exp(grf_long(indrf,indmaturity2))+loan_premia(indloan_rf,1)-1) !(This is identical to using total_payment on long_loan + interest on bullet part)
                                principal_payment = long_debt_repay(indmaturity2,1,indlongdebt2,indrf)
                                !PTI = (principal_payment+interest_payment)/(g(t,1)*exp(gyp(indyp,1)))  
                                LTV = (g_debt_long(indmaturity2,1,indlongdebt2,indrf))/A_value
                                DTI = (g_debt_long(indmaturity2,1,indlongdebt2,indrf))/(g(t,1)*exp(gyp(indyp,1)))
                        
                                !if ((PTI < PTI_lim) .and. (LTV < LTV_lim) .and. (DTI < DTI_lim)) then           
                                if ((LTV .le. LTV_lim) .and. (DTI .le. DTI_lim)) then               
                            
                                    !Pay all the old debts
                                    interest_payment = (g_debt_long(indmaturity,indloanper,indlongdebt,indloan_rf))*(exp(grf_long(indloan_rf,indmaturity))+loan_premia(indloan_rf,1)-1)
                                    principal_payment = g_debt_long(indmaturity,indloanper,indlongdebt,indloan_rf)
                            
                        
                                    invest3(:,indlongdebt2,indmaturity2) =  gcash(indcash,1) - gc(:,1) !consume
                                    invest3(:,indlongdebt2,indmaturity2) =  invest3(:,indlongdebt2,indmaturity2) - principal_payment - interest_payment !Pay all old debts
                                    invest3(:,indlongdebt2,indmaturity2) =  invest3(:,indlongdebt2,indmaturity2) + g_debt_long(indmaturity2,1,indlongdebt2,indrf) - prepayment_penalty2  - origination_fee2 - A_value * (prop_maintenance + prop_tax)!Get new loans (this allows for prepayment)
                                elseif ((underwater_refi .eq. 1) .and. (indlongdebt2 .le. indlongdebt)) then !Agents allowed to refinance to original loan amount
                                    interest_payment = (g_debt_long(indmaturity,indloanper,indlongdebt,indloan_rf))*(exp(grf_long(indloan_rf,indmaturity))+loan_premia(indloan_rf,1)-1)
                                    principal_payment = g_debt_long(indmaturity,indloanper,indlongdebt,indloan_rf)
                                
                                                                    
                                    invest3(:,indlongdebt2,indmaturity2) =  gcash(indcash,1) - gc(:,1) !consume
                                    invest3(:,indlongdebt2,indmaturity2) =  invest3(:,indlongdebt,indmaturity2) - principal_payment - interest_payment !Pay all old debts
                                    invest3(:,indlongdebt2,indmaturity2) =  invest3(:,indlongdebt,indmaturity2) + g_debt_long(indmaturity2,1,indlongdebt2,indrf) - prepayment_penalty2  - origination_fee2 - A_value * (prop_maintenance + prop_tax)!Get new loan with identical size to old one, but reset interest rate.
                                else
                                    invest3(:,indlongdebt2,indmaturity2) = infinity
                                end if
                                                            
                            end do
                        end do
                    else !First period choice 
                        do indmaturity2 = 1,n_maturity
                            do indlongdebt2 = 1,n_long_debt       
                                !No refinancing cost on first period
                                prepayment_penalty2 = 0.0
                                origination_fee2 = origination_fee*g_debt_long(indmaturity2,1,indlongdebt2,indrf)
                                
                                !Check restrictions
                                interest_payment = (g_debt_long(indmaturity2,1,indlongdebt2,indrf))*(exp(grf_long(indrf,indmaturity2))+loan_premia(indrf,1)-1) 
                                principal_payment = long_debt_repay(indmaturity2,1,indlongdebt2,indrf)
                                !PTI = (principal_payment+interest_payment)/(g(t,1)*exp(gyp(indyp,1)))  
                                LTV = (g_debt_long(indmaturity2,1,indlongdebt2,indrf))/A_value
                                DTI = (g_debt_long(indmaturity2,1,indlongdebt2,indrf))/(g(t,1)*exp(gyp(indyp,1)))
                                
                                !Assume agents have no mortgage debt prior to period 1. So always start agents on the simulation at first point of indlongdebt
                                !if ((PTI < PTI_lim) .and. (LTV < LTV_lim) .and. (DTI < DTI_lim)) then           
                                if ((LTV .le. LTV_lim) .and. (DTI .le. DTI_lim)) then               
                                    invest3(:,indlongdebt2,indmaturity2) =  gcash(indcash,1) - gc(:,1) !consume
                                    invest3(:,indlongdebt2,indmaturity2) =  invest3(:,indlongdebt2,indmaturity2) - A_value + g_debt_long(indmaturity2,1,indlongdebt2,indrf) - prepayment_penalty2 - origination_fee2 !Get new loan and pay for the asset and origination fee
                                    invest3(:,indlongdebt2,indmaturity2) =  invest3(:,indlongdebt2,indmaturity2) - A_value * (prop_maintenance + prop_tax)
                                else 
                                    invest3(:,indlongdebt2,indmaturity2) = infinity
                                end if 
                            end do
                        end do
                        invest2(:,1) = infinity !First period there's no loan outstanding. So making normal payment not possible
                    end if
                        
            
                    nw3(:,:,:) = invest3(:,:,:) * (exp(grf(indrf,1)))
                    nw2(:,1) = invest2(:,1) * (exp(grf(indrf,1)))
                    
                    v3(:,:,:) = 0.0
                    v2(:,:) = 0.0
                    
                    !This is refinancing decision
                    do indlongdebt2 = 1,n_long_debt
                        do indmaturity2 = 1,n_maturity
                            do indc = 1, n_c

                                indloan_rf2 = indrf
                                indloanper2 = 1
                                
                                auxv1 = 0.0    
                                if (nw3(indc,indlongdebt2,indmaturity2) .ge. 0.0) then 
                                    do indyp2 = 1, n_yp
                                        do indyt = 1, nqp    
                                            do indrf2 = 1,n_rf
                                                do indbcycle2 = 1, n_bcycle             
                                                    inc = g(t,1)*exp(gyp(indyp2,1))*gyt_levels(indyt,indbcycle2)
                                                    cash = nw3(indc,indlongdebt2,indmaturity2) + inc*(1-income_tax)
                                                    cash = max(cash, cash_min) 
                                                    cash = min(cash, cash_max)
                                                    
                                                    !(1)
                                                    !call splint(gcash(:,1), v_debt(:,indyp2,indrf2,indbcycle2, indlongdebt2, indloanper2, indmaturity2, indloan_rf2, t+1), &
                                                    !            bscoef8(:,indyp2,indrf2,indbcycle2, indlongdebt2, indloanper2, indmaturity2, indloan_rf2), n_cash, cash,1,1,auxv)
                                                    
                                                    !(2)
                                                    !inbvx = 1
                                                    !call db1val(cash,0,Xknot,n_cash,Kxord,bscoef8(:,indyp2,indrf2,indbcycle2, indlongdebt2, indloanper2, indmaturity2, indloan_rf2),auxv,iflag,inbvx,w0)
                                                    
                                                    !(3)
                                                    call interp1D_QL(gamma, n_cash, gcash(:,1), bscoef8(:,indyp2,indrf2,indbcycle2, indlongdebt2, indloanper2, indmaturity2, indloan_rf2), cash, auxv) 

                                                    auxv1 = auxv1+auxv*Tmyp(indyp,indyp2)*weig_yt(indyt,indbcycle2)*Tmrf(indrf,indrf2)*Tmbcycle(indbcycle, indbcycle2)
                                                end do 
                                            end do
                                        end do 
                                    end do
                                    v3(indc,indlongdebt2,indmaturity2) = u(indc,1) + delta*auxv1
                                else
                                    v3(indc,indlongdebt2,indmaturity2) = infinity
                                end if
                            end do
                        end do
                    end do
                    
                    !This is non-refinancing decision
                    do indc = 1, n_c

                        indloan_rf2 = indloan_rf
                        indloanper2 = indloanper+1
                                                
                        auxv1 = 0.0    
                        if (nw2(indc,1) .ge. 0.0) then 
                            do indyp2 = 1, n_yp
                                do indyt = 1, nqp    
                                    do indrf2 = 1,n_rf
                                        do indbcycle2 = 1, n_bcycle             
                                            inc = g(t,1)*exp(gyp(indyp2,1))*gyt_levels(indyt,indbcycle2)
                                            cash = nw2(indc,1) + inc*(1-income_tax)
                                            cash = max(cash, cash_min) 
                                            cash = min(cash, cash_max)
                                            
                                            !(1)
                                            !call splint(gcash(:,1), v_debt(:,indyp2,indrf2,indbcycle2, indlongdebt, indloanper2, indmaturity, indloan_rf2, t+1), &
                                            !            bscoef8(:,indyp2,indrf2,indbcycle2, indlongdebt, indloanper2, indmaturity, indloan_rf2), n_cash, cash,1,1,auxv)
                                            
                                            !(2)
                                            !inbvx = 1
                                            !call db1val(cash,0,Xknot,n_cash,Kxord,bscoef8(:,indyp2,indrf2,indbcycle2, indlongdebt, indloanper2, indmaturity, indloan_rf2),auxv,iflag,inbvx,w0)
                                            
                                            !(3)
                                            call interp1D_QL(gamma, n_cash, gcash(:,1), bscoef8(:,indyp2,indrf2,indbcycle2, indlongdebt, indloanper2, indmaturity, indloan_rf2), cash, auxv) 

                                            auxv1 = auxv1+auxv*Tmyp(indyp,indyp2)*weig_yt(indyt,indbcycle2)*Tmrf(indrf,indrf2)*Tmbcycle(indbcycle, indbcycle2)
                                        end do 
                                    end do
                                end do 
                            end do
                            v2(indc,1) = u(indc,1) + delta*auxv1
                        else
                            v2(indc,1) = infinity
                        end if
                    end do

                        

                    
                    pt3 = maxloc(v3)
                    pt2 = maxloc(v2)
                    
                    if (v3(pt3(1),pt3(2),pt3(3)) .ge.  v2(pt2(1),1)) then
                    
                        c_debt(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = gc(pt3(1),1)
                        refinancing_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = 1
                        debt_long_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = g_debt_long(pt3(3),1,pt3(2),1) !The last index here does not matter since this pf only keeps track of the initial loan amount outstanding which is independent of interest rate
                        maturity_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = g_maturity(pt3(3),1)
                        v_debt(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = v3(pt3(1),pt3(2),pt3(3))
                        
                    else
                        c_debt(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = gc(pt2(1),1)
                        refinancing_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = 0
                        debt_long_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = g_debt_long(indmaturity,1,indlongdebt,1) !No change in loan amount
                        maturity_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = g_maturity(indmaturity,1) !no change in maturity
                        v_debt(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) =  v2(pt2(1),1)
                        
                    end if 
                    
                    !Check default decision
                    !Check if agent:
                    !(1) Does not default
                    v_keep =  v_debt(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t)
                    
                    !Cash out 
                    interest_payment = (g_debt_long(indmaturity,indloanper,indlongdebt,indloan_rf))*(exp(grf_long(indloan_rf,indmaturity))+loan_premia(indloan_rf,1)-1)
                    principal_payment = g_debt_long(indmaturity,indloanper,indlongdebt,indloan_rf)
                    
                    he =  A_value - principal_payment 
                                                    
                    if (indloanper .eq. debt_maturity) then 
                        prepayment_penalty2 = 0.0
                    else
                        prepayment_penalty2 = prepayment_penalty*g_debt_long(indmaturity,indloanper,indlongdebt,indloan_rf)
                    end if
                    he_net = A_value*(1-asset_sale_fee) - principal_payment - prepayment_penalty2
                    
                    if (he_net - interest_payment .ge. 0.0) then
                        cash = gcash(indcash,1) + he_net - interest_payment
                        pt1 = minloc(abs(gcash(:,1)-cash))
                        v_cashout = v_nodebt(pt1(1),indyp,indrf,indbcycle,t)
                    else
                        v_cashout = infinity
                    end if
                    if (t .eq. 1) then !agent has no asset so no cashout possible
                        v_cashout = infinity
                    end if
                    
                        
                        
                    
                    
                    !(iii) Defaults
                    
                    v_default = v_nodebt(indcash,indyp,indrf,indbcycle,t) - default_penalty
                    if (t .eq. 1) then 
                        !on the first period default indicator tells me if agent prefers to move to rental market 
                        !(he always starts with zero debt)
                        !hopefully rental premium is enough to give agent incentive to lever up and buy the asset
                        !v_default = v_nodebt(indcash,indyp,indrf,indbcycle,t) 
                        
                        !No default possible on first period
                        v_default = infinity
                    end if
                
                    if ((v_keep .ge. v_default) .and. (v_keep .ge. v_cashout)) then !keep loan
                        v_debt(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = v_keep
                        default_pf(indcash,indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t)  = 1.0
                    elseif ((v_cashout .ge. v_keep) .and. (v_cashout .ge. v_default)) then !cash out loan 
                        v_debt(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = v_cashout
                        default_pf(indcash,indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = 2.0
                        debt_long_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = 0.0
                        maturity_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = 0.0
                        refinancing_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = 0.0
                    else  !default loan
                        v_debt(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = v_default
                        default_pf(indcash,indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = 3.0
                        debt_long_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t)  = 0.0
                        maturity_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = 0.0
                        refinancing_pf(indcash, indyp, indrf, indbcycle, indlongdebt, indloanper,  indmaturity, indloan_rf, t) = 0.0
                    
                    end if

                end do
                call MPI_BARRIER( MPI_COMM_WORLD, ierror) 
            end do
        end do
        !$omp end do
        !$omp end parallel
        end do
        end do
        end do
        end do
        end do


        
        midTime2 = MPI_WTIME()
        write(*,"(A,F15.1,A)") "Total time elapsed",midTime2-startTime," seconds." 
        
        ! write(filename,'(a,i0,a)') trim(path)//'c_debt', t, '.txt'
        ! open(unit=25,file=filename,status='replace',action='write',form='formatted')
        ! write(25,*) c_debt(:,:,:,:,:,:,:,:,t)
        ! close(unit=25)
        ! 
        ! write(filename,'(a,i0,a)')trim(path)//'v_debt', t, '.txt'
        !open(unit=25,file=filename,status='replace',action='write',form='formatted')
        ! write(25,*) v_debt(:,:,:,:,:,:,:,:,t)
        ! close(unit=25)
        ! 
        ! write(filename,'(a,i0,a)')trim(path)//'debt_long_pf', t, '.txt'
        ! open(unit=25,file=filename,status='replace',action='write',form='formatted')
        ! write(25,*) debt_long_pf(:,:,:,:,:,:,:,:,t)
        ! close(unit=25)
        ! 
        ! write(filename,'(a,i0,a)')trim(path)//'maturity_pf', t, '.txt'
        ! open(unit=25,file=filename,status='replace',action='write',form='formatted')
        ! write(25,*) maturity_pf(:,:,:,:,:,:,:,:,t)
        ! close(unit=25)
        ! 
        ! write(filename,'(a,i0,a)')trim(path)//'default_pf', t, '.txt'
        ! open(unit=25,file=filename,status='replace',action='write',form='formatted')
        ! write(25,*) default_pf(:,:,:,:,:,:,:,:,t)
        ! close(unit=25)
        ! 
        ! write(filename,'(a,i0,a)')trim(path)//'refinancing_pf', t, '.txt'
        ! open(unit=25,file=filename,status='replace',action='write',form='formatted')
        ! write(25,*) refinancing_pf(:,:,:,:,:,:,:,:,t)
        ! close(unit=25)
        
    end do

!Save time elapsed    
    write(filename,'(a)')trim(path)//'time_elapsed'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) midTime2-startTime 
    close(unit=25)


    !Save policy functions for matlab    
    write(filename,'(a)')trim(path)//'c_debt.txt'
    open(unit=25,file=filename,status='replace',action='write',form='formatted')
    write(25,*) c_debt(:,:,:,:,:,:,:,:,:)
    close(unit=25)

    write(filename,'(a)')trim(path)//'v_debt.txt'
    open(unit=25,file=filename,status='replace',action='write',form='formatted')
    write(25,*) v_debt(:,:,:,:,:,:,:,:,:)
    close(unit=25)

    write(filename,'(a)')trim(path)//'default_pf.txt'
    open(unit=25,file=filename,status='replace',action='write',form='formatted')
    write(25,*) default_pf(:,:,:,:,:,:,:,:,:)
    close(unit=25)

    write(filename,'(a)')trim(path)//'debt_long_pf.txt'
    open(unit=25,file=filename,status='replace',action='write',form='formatted')
    write(25,*) debt_long_pf(:,:,:,:,:,:,:,:,:)
    close(unit=25)

    write(filename,'(a)')trim(path)//'maturity_pf.txt'
    open(unit=25,file=filename,status='replace',action='write',form='formatted')
    write(25,*) maturity_pf(:,:,:,:,:,:,:,:,:)
    close(unit=25)

    write(filename,'(a)')trim(path)//'refinancing_pf.txt'
    open(unit=25,file=filename,status='replace',action='write',form='formatted')
    write(25,*) refinancing_pf(:,:,:,:,:,:,:,:,:)
    close(unit=25)

    midTime2 = MPI_WTIME()
    write(*,"(A,F15.1,A)") "Total time elapsed",midTime2-startTime," seconds." 

    write(filename,'(a)')trim(path)//'total_time_elapsed.txt'
    open(unit=25,file=filename,status='replace',action='write',form='formatted')
    write(25,*) midTime2-startTime
    close(unit=25)
    

    !Save variables and grids to open in matlab
    write(filename,'(a)')trim(path)//'t_death'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) t_death 
    close(unit=25)
    write(filename,'(a)')trim(path)//'t_born'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) t_born 
    close(unit=25)
    write(filename,'(a)')trim(path)//'gcash'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) gcash 
    close(unit=25)
    write(filename,'(a)')trim(path)//'n_cash'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) n_cash 
    close(unit=25)
    write(filename,'(a)')trim(path)//'g_maturity'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) g_maturity 
    close(unit=25)
    write(filename,'(a)')trim(path)//'n_yp'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) n_yp 
    close(unit=25)
    write(filename,'(a)')trim(path)//'n_rf'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) n_rf 
    close(unit=25)
    write(filename,'(a)')trim(path)//'n_bcycle'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) n_bcycle 
    close(unit=25)
    write(filename,'(a)')trim(path)//'n_maturity'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) n_maturity
    close(unit=25)
    write(filename,'(a)')trim(path)//'n_long_debt'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) n_long_debt
    close(unit=25)
    !write(filename,'(a)')trim(path)//'c'
    !open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
    !    write(25,*) c
    !close(unit=25)
    !write(filename,'(a)')trim(path)//'b1'
    !open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
    !    write(25,*) b1
    !close(unit=25)
    !write(filename,'(a)')trim(path)//'b2'
    !open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
    !    write(25,*) b2
    !close(unit=25)
    !write(filename,'(a)')trim(path)//'b3'
    !open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
    !    write(25,*) b3
    !close(unit=25)
    write(filename,'(a)')trim(path)//'sig_yp'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) sig_yp
    close(unit=25)
    write(filename,'(a)')trim(path)//'rho_yp'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) rho_yp
    close(unit=25)
    write(filename,'(a)')trim(path)//'sig_yt'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) sig_yt
    close(unit=25)
    write(filename,'(a)')trim(path)//'g'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) g
    close(unit=25)
    write(filename,'(a)')trim(path)//'Tmyp'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) Tmyp
    close(unit=25)
    write(filename,'(a)')trim(path)//'gyp'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) gyp
    close(unit=25)
    write(filename,'(a)')trim(path)//'gyt'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) gyt
    close(unit=25)
    write(filename,'(a)')trim(path)//'DTI_lim'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) DTI_lim
    close(unit=25)
    write(filename,'(a)')trim(path)//'LTV_lim'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) LTV_lim
    close(unit=25)
    !write(filename,'(a)')trim(path)//'PTI_lim'
    !open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
    !    write(25,*) PTI_lim
    !close(unit=25)
    write(filename,'(a)')trim(path)//'loan_premia'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) loan_premia
    close(unit=25)
    write(filename,'(a)')trim(path)//'prepayment_penalty'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) prepayment_penalty
    close(unit=25)
    write(filename,'(a)')trim(path)//'origination_fee'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) origination_fee
    close(unit=25)
    write(filename,'(a)')trim(path)//'cash_min'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) cash_min
    close(unit=25)
    write(filename,'(a)')trim(path)//'cash_max'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) cash_max
    close(unit=25)
    write(filename,'(a)')trim(path)//'mu_rf'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) mu_rf
    close(unit=25)
    write(filename,'(a)')trim(path)//'sig_rf'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) sig_rf
    close(unit=25)
    write(filename,'(a)')trim(path)//'rho_rf'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) rho_rf
    close(unit=25)
    write(filename,'(a)')trim(path)//'term_premium'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) term_premium
    close(unit=25)
    write(filename,'(a)')trim(path)//'mu_h'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) mu_h
    close(unit=25)
    write(filename,'(a)')trim(path)//'delta_h'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) delta_h
    close(unit=25)
    write(filename,'(a)')trim(path)//'A0'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) A0
    close(unit=25)
    write(filename,'(a)')trim(path)//'rental_premium'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) rental_premium
    close(unit=25)
    write(filename,'(a)')trim(path)//'asset_sale_fee'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) asset_sale_fee
    close(unit=25)
    write(filename,'(a)')trim(path)//'c_min'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) c_min
    close(unit=25)
    write(filename,'(a)')trim(path)//'c_max'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) c_max
    close(unit=25)
    write(filename,'(a)')trim(path)//'gc'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) gc
    close(unit=25)
    write(filename,'(a)')trim(path)//'n_c'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) n_c
    close(unit=25)
    write(filename,'(a)')trim(path)//'nqp'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) nqp
    close(unit=25)
    write(filename,'(a)')trim(path)//'weig'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) weig
    close(unit=25)
    write(filename,'(a)')trim(path)//'grid'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) grid
    close(unit=25)
    write(filename,'(a)')trim(path)//'A0'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) A0
    close(unit=25)
    write(filename,'(a)')trim(path)//'Tmbcycle'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) Tmbcycle
    close(unit=25)
    write(filename,'(a)')trim(path)//'Tmrf'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) Tmrf
    close(unit=25)
    write(filename,'(a)')trim(path)//'g'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) g
    close(unit=25)
    write(filename,'(a)')trim(path)//'grf'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) grf
    close(unit=25)
    write(filename,'(a)')trim(path)//'prop_tax'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) prop_tax
    close(unit=25)
    write(filename,'(a)')trim(path)//'prop_maintenance'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) prop_maintenance
    close(unit=25)
    write(filename,'(a)')trim(path)//'weig_yt'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) weig_yt
    close(unit=25)
    write(filename,'(a)')trim(path)//'income_tax'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) income_tax
    close(unit=25)   
    write(filename,'(a)')trim(path)//'gyt_levels'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) gyt_levels
    close(unit=25)  
    write(filename,'(a)')trim(path)//'gamma'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) gamma
    close(unit=25)  
    write(filename,'(a)')trim(path)//'bqt'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) bqt
    close(unit=25)  
    write(filename,'(a)')trim(path)//'delta'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) delta
    close(unit=25)
    write(filename,'(a)')trim(path)//'start_inc'
    open(UNIT=25,FILE=filename,STATUS='replace',ACTION='write',form='formatted')
        write(25,*) start_inc
    close(unit=25)

    deallocate(gcash,gyp,Tmyp,grf,grf_long,Tmrf,gc,Tmbcycle,g,gyt,weig_yt,jensen_term,gyt_levels,Xknot,bscoef4,c_nodebt,v_nodebt,invest2,nw2,u,v2,w0,g_debt_long,long_debt_repay,tot_payment,annuity_factor,g_maturity,c_debt,v_debt,default_pf,refinancing_pf,debt_long_pf,maturity_pf,invest3,refinancing_choice,nw3,v3,bscoef8)
    deallocate(grid, weig)           

    call MPI_FINALIZE(ierror)


end program