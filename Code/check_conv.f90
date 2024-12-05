
      subroutine check_conv(av,g,d_avg,d_max)

!     This subroutine checks the residuals in all primary flow variables and
!     prints their values, you should not need to change this subroutine

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(in) :: g(:)
      real, intent(out) :: d_avg, d_max
      real :: d_avg_array(av%nn), d_max_array(av%nn)
      real, allocatable :: dro(:,:), droe(:,:), drovx(:,:), drovy(:,:)
      integer :: ijx_max(2), ijy_max(2), ij_max(2), ncells
      real :: dro_max(av%nn), drovx_max(av%nn), drovy_max(av%nn), droe_max(av%nn), dro_avg(av%nn),&
      drovx_avg(av%nn), drovy_avg(av%nn), droe_avg(av%nn), flow_ratio
      character(len=100) :: fmt_step
            integer :: n
      do n=1,av%nn
      
      allocate(dro(g(n)%ni,g(n)%nj),droe(g(n)%ni,g(n)%nj),drovx(g(n)%ni,g(n)%nj),drovy(g(n)%ni,g(n)%nj))
      

!     Get the number of cells from the size of the residual arrays
      ncells = size(g(n)%dro)

!     Use "abs" to make all residual values positive and store locally
      dro = abs(g(n)%dro); droe = abs(g(n)%droe);
      drovx = abs(g(n)%drovx); drovy = abs(g(n)%drovy);

!     Calculate the mean changes for each variable
      dro_avg(n) = sum(abs(dro)) / (ncells * av%ro_ref)
      droe_avg(n) = sum(abs(droe)) / (ncells * av%roe_ref)
      drovx_avg(n) = sum(abs(drovx)) / (ncells * av%rov_ref)
      drovy_avg(n) = sum(abs(drovy)) / (ncells * av%rov_ref)

!     Find the maximum value of change for the momenta and the positions
      dro_max(n) = maxval(dro) / av%ro_ref; droe_max = maxval(droe) / av%roe_ref;
      ijx_max = maxloc(drovx); ijy_max = maxloc(drovy);
      drovx_max(n) = drovx(ijx_max(1),ijx_max(2)) / av%rov_ref
      drovy_max(n) = drovy(ijy_max(1),ijy_max(2)) / av%rov_ref

!     Store single values as the maximum of either the x or y-momentum
      if(drovx_avg(n) > drovy_avg(n)) then
          d_max_array(n) = drovx_max(n); d_avg_array(n) = drovx_avg(n); ij_max = ijx_max;
      else
          d_max_array(n) = drovy_max(n); d_avg_array(n) = drovy_avg(n); ij_max = ijy_max;
      end if


     

          
          
          deallocate(dro,droe,drovx,drovy)
      end do
      !     Write a short human readable output summary to the screen.
      write(6,*) 'Time step number ', av%nstep
      fmt_step = '(a,e10.3,a,i4,a,i4,a,e10.3)'
      write(*,fmt_step) '   d_max =', maxval(d_max_array), ' at i =', ij_max(1), ', j =', &
          ij_max(2), ', d_avg =', maxval(d_avg_array)
          
          
      !     Write the average and maximum changes in the primary variables to unit 3
!     for convergenge plotting
      write(3,'(i13,8e15.6)') av%nstep, maxval(dro_avg), maxval(droe_avg), maxval(drovx_avg), &
          maxval(drovy_avg), maxval(dro_max), maxval(droe_max), maxval(drovx_max), maxval(drovy_max)

      end subroutine check_conv


