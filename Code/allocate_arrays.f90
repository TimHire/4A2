      
      subroutine allocate_arrays(av,g,bcs)

!     Allocate memory for all arrays in the grid and bcs datatype, this has been
!     completed for the basic solver. If you add further arrays to the code in
!     the extensions you will need to allocate them here.

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), allocatable, intent(inout) :: g(:)
      type(t_bconds), intent(inout) :: bcs
      integer :: ni, nj
      
      allocate(g(1))

!     Get the size of the mesh and store locally for convenience
      ni = av%ni; nj = av%nj;

!     Copy the mesh size to the grid datatype
      g(1)%ni = ni; g(1)%nj = nj;

!     Wall flag array
      allocate(g(1)%wall(ni,nj))

!     Arrays to store static conditions at the inlet plane
      allocate(bcs%ro(nj),bcs%p(nj))

!     Node coordinates in the mesh
      allocate(g(1)%x(ni,nj),g(1)%y(ni,nj))

!     Cell areas and projected side lengths at the centre of each respectively
      allocate(g(1)%area(ni-1,nj-1),g(1)%lx_i(ni,nj-1),g(1)%ly_i(ni,nj-1), &
          g(1)%lx_j(ni-1,nj),g(1)%ly_j(ni-1,nj))

!     Primary flow variables in the mesh
      allocate(g(1)%ro(ni,nj),g(1)%rovx(ni,nj),g(1)%rovy(ni,nj),g(1)%roe(ni,nj))

!     Primary flow variables in the mesh_start
      allocate(g(1)%ro_start(ni,nj),g(1)%rovx_start(ni,nj),g(1)%rovy_start(ni,nj),g(1)%roe_start(ni,nj))

!     Cell centred primary increments
      allocate(g(1)%dro(ni-1,nj-1),g(1)%drovx(ni-1,nj-1), &
          g(1)%drovy(ni-1,nj-1),g(1)%droe(ni-1,nj-1))

!     Secondary variables stored at the nodes for convenience
      allocate(g(1)%p(ni,nj),g(1)%hstag(ni,nj),g(1)%vx(ni,nj),g(1)%vy(ni,nj))


      end subroutine allocate_arrays


