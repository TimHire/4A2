      
      subroutine apply_bconds(av,g,bcs,constant_enthalpy)

!     This subroutine applies both the inlet and outlet boundary conditions, as
!     it modifies both the primary and secondary flow variables they must be
!     calculated first

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(inout) :: bcs
      logical, intent(in) :: constant_enthalpy

!     Declare the other variables you need here
!     INSERT
!      real :: t(g%nj), v(g%nj), p(g%nj)
      real :: t(1:av%ni,1:av%nj), v(1:av%ni,1:av%nj)


!     At the inlet boundary the change in density is driven towards "rostag",
!     which is then used to obtain the other flow properties to match the
!     specified stagnation pressure, temperature and flow angle. 

!     To help prevent instabilities forming at the inlet boundary condition the 
!     changes in inlet density are relaxed by a factor "rfin" normally set to 
!     0.25 but it can be reduced further.

!     It is also worth checking if "ro" is greater than "rostag" and limiting 
!     the values to be slightly less than "rostag". This can prevent the solver 
!     crashing during severe transients.
      if(av%nstep == 1) then
          bcs%ro = g%ro(1,:)
      else
          bcs%ro = bcs%rfin * g%ro(1,:) + (1 - bcs%rfin) * bcs%ro
      endif
      bcs%ro = min(bcs%ro,0.9999 * bcs%rostag)

!     Calculate "p(1,:)", "rovx(1,:)", "rovy(1,:)" and "roe(1,:)" from the inlet 
!     "ro(:)", "pstag", "tstag" and "alpha". Also set "vx(1,:)", "vy(1,:)" and 
!     "hstag(1,:)"
!     INSERT
!     Need to calculate the velocities for each cell at the inlet
      t(1,:) = bcs%tstag * (bcs%ro / bcs%rostag) ** (av%gam - 1)
      v(1,:) = sqrt(2 * av%cp * (bcs%tstag - t(1,:)))

      g%p(1,:) = bcs%pstag * (bcs%ro / bcs%rostag)**av%gam
      g%vx(1,:) = v(1,:) * cos(bcs%alpha)
      g%vy(1,:) = v(1,:) * sin(bcs%alpha)
      g%rovx_start(1,:) = g%vx(1,:) * bcs%ro
      g%rovy_start(1,:) = g%vy(1,:) * bcs%ro

      if (.not. constant_enthalpy) then
           g%roe_start(1,:) = bcs%ro * (av%cv * t(1,:) + 0.5 * v(1,:)**2)   ! constant stagnation enthalpy
           g%hstag(1,:) = (g%roe(1,:) + g%p(1,:) ) / bcs%ro 
!          g%hstag(1,:) = av%cp * bcs%tstag
      end if

!     For the outlet boundary condition set the value of "p(ni,:)" to the
!     specified value of static pressure "p_out" in "bcs"
!     INSERT
      g%p(av%ni,:) = bcs%p_out

      end subroutine apply_bconds


