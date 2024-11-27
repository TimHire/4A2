      
      subroutine set_secondary(av,g,bcs,constant_enthalpy)

!     This subroutine calculates the secondary flow variables from the primary
!     ones at every node

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g(:)
      type(t_bconds), intent(in) :: bcs
      logical, intent(in) :: constant_enthalpy
      integer :: n
      do n=1,av%nn

!     Define any further variables you may need
!     INSERT
      
      

!     The primary flow variables are "ro", "roe", "rovx" and "rovy", these are 
!     the conserved quantities within the Euler equations. Write the code to
!     calculate the secondary flow variables which are the velocity components
!     "vx" and "vy", the static pressure "p" and the stagnation enthalpy
!     "hstag". These are needed at every timestep, there is no need for any 
!     loops as the operations can be performed elementwise, although you may
!     wish to define some intermediate variables to improve readability.
!     INSERT
      g(n)%vx = g(n)%rovx / g(n)%ro
      g(n)%vy = g(n)%rovy / g(n)%ro
      if (constant_enthalpy) then
      g(n)%p = g(n)%ro * av%rgas * (bcs%tstag - 0.5 * hypot(g(n)%vx, g(n)%vy)**2 / av%cp)
      g(n)%hstag = av%cp * bcs%tstag
      g(n)%roe = g(n)%hstag * g(n)%ro - g(n)%p
      else
      g(n)%p = (av%rgas/av%cv) * (g(n)%roe - g(n)%ro * 0.5 * hypot(g(n)%vx, g(n)%vy)**2)
      g(n)%hstag = (g(n)%roe + g(n)%p) / g(n)%ro
      end if
      
      end do

      end subroutine set_secondary


