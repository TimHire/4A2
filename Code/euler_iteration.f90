
      subroutine euler_iteration(av,g,p,constant_enthalpy)

!     This subroutine calculates the fluxes into each cell and then sums them to
!     update the primary flow properties

!     Explicitly declare the required variables
      use types
      use flux_stencil
      use smooth_stencil
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g(:)
      type(t_match), intent(inout) :: p(:)
      logical, intent(in) :: constant_enthalpy
      real, allocatable :: mass_i(:,:), flux_i(:,:)
      real, allocatable :: mass_j(:,:), flux_j(:,:)
      integer :: i, j, ni, nj
      integer :: n
      do n=1,av%nn
      

      allocate(mass_i(g(n)%ni,g(n)%nj-1));  allocate(flux_i(g(n)%ni,g(n)%nj-1));
      allocate(mass_j(g(n)%ni - 1,g(n)%nj));  allocate(flux_j(g(n)%ni - 1,g(n)%nj));

!     Get the block size and store locally for convenience
      ni = g(n)%ni; nj = g(n)%nj

!     Setup the continuity equation by calculating the mass flow through
!     the facets in both the i and j-directions. Store these values in
!     "mass_i" and "mass_j"
!     INSERT
      mass_i = ((g(n)%rovx(:,1:nj-1) + g(n)%rovx(:,2:nj)) * g(n)%lx_i(:,1:nj-1) +&
               (g(n)%rovy(:,1:nj-1) + g(n)%rovy(:,2:nj)) * g(n)%ly_i(:,1:nj-1)) / 2
               
      mass_j = ((g(n)%rovx(1:ni-1,:) + g(n)%rovx(2:ni,:)) * g(n)%lx_j(1:ni-1,:) +&
               (g(n)%rovy(1:ni-1,:) + g(n)%rovy(2:ni,:)) * g(n)%ly_j(1:ni-1,:)) / 2
     
!     Apply the wall boundary condition by checking that two nodes at the
!     end of a facet are both on a wall, if so then the appropriate mass
!     flow array is set to have zero flow through that facet
      where(g(n)%wall(1:ni-1,:) .and. g(n)%wall(2:ni,:)) mass_j = 0 
      where(g(n)%wall(:,1:nj-1) .and. g(n)%wall(:,2:nj)) mass_i = 0 

!     Update the density with mass fluxes by calling "sum_fluxes"
!     INSERT
!     Parsing in out mass functions so also need to parse in related cell property ro and its derivative dro
      call sum_fluxes(av,mass_i,mass_j,g(n)%area,g(n)%ro_start,g(n)%dro)


      if (.not. constant_enthalpy) then
!     Setup the conservation of energy equation by calculated the enthalpy flux
!     and storing the values in "flux_i" and "flux_j", you will need "mass_i"
!     and "mass_j" from before
!     INSERT
      flux_i = mass_i * (g(n)%hstag(:,1:nj-1) + g(n)%hstag(:,2:nj)) / 2
      flux_j = mass_j * (g(n)%hstag(1:ni-1,:) + g(n)%hstag(2:ni,:)) / 2

!     Update the internal energy with enthalpy fluxes
!     INSERT
      call sum_fluxes(av,flux_i,flux_j,g(n)%area,g(n)%roe_start,g(n)%droe)
      end if



!     After discussing with Felix decided needed to do pressure times lenth for each section before averaging

!     Setup the x-momentum equation including momentum flux and pressure forces
!     INSERT
      flux_i = (mass_i * (g(n)%vx(:,1:nj-1) + g(n)%vx(:,2:nj)) + (g(n)%p(:,1:nj-1) + g(n)%p(:,2:nj) ) * g(n)%lx_i ) / 2
      flux_j = (mass_j * (g(n)%vx(1:ni-1,:) + g(n)%vx(2:ni,:)) + (g(n)%p(1:ni-1,:) + g(n)%p(2:ni,:) ) * g(n)%lx_j ) / 2

!     Update the x-momentum with momentum flux
!     INSERT
      call sum_fluxes(av,flux_i,flux_j,g(n)%area,g(n)%rovx_start,g(n)%drovx)

!     Setup the y-momentum equation including momentum flux and pressure forces
!     INSERT
      flux_i = (mass_i * (g(n)%vy(:,1:nj-1) + g(n)%vy(:,2:nj)) + (g(n)%p(:,1:nj-1) + g(n)%p(:,2:nj) ) * g(n)%ly_i ) / 2
      flux_j = (mass_j * (g(n)%vy(1:ni-1,:) + g(n)%vy(2:ni,:)) + (g(n)%p(1:ni-1,:) + g(n)%p(2:ni,:) ) * g(n)%ly_j ) / 2

!     Update the y-momentum with momentum flux
!     INSERT
      call sum_fluxes(av,flux_i,flux_j,g(n)%area,g(n)%rovy_start,g(n)%drovy)
     
     
!     Reset the variables to be used again in the solver loop
      g(n)%ro = g(n)%ro_start
      if (.not. constant_enthalpy) then
      g(n)%roe = g(n)%roe_start
      end if 
      g(n)%rovx = g(n)%rovx_start
      g(n)%rovy = g(n)%rovy_start
      
      

      
      
      
     
      !     Add artificial viscosity by smoothing all of the primary flow variables
      call smooth_array(av,g(n)%ro)
      if (.not. constant_enthalpy) then
      call smooth_array(av,g(n)%roe)
      end if
      call smooth_array(av,g(n)%rovx)
      call smooth_array(av,g(n)%rovy)
      
      deallocate(mass_i, flux_i, mass_j, flux_j)
      end do

      

      end subroutine euler_iteration


