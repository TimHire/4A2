      
      subroutine generate_mesh(geom,g,geometric_mesh)

!     Create cells of the mesh to cover the domain defined by geometry curves,
!     the values of the node coordinates, x(i,j) and y(i,j) are stored in "g"

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_geometry), intent(in) :: geom
      type(t_grid), intent(inout) :: g(:)
      logical, intent(in) :: geometric_mesh
      real :: si_a(geom%ni_a), si_b(geom%ni_b)
      real, allocatable :: si(:), sj(:), si1(:), sj1(:), si2(:), sj2(:)
      integer :: ni, nj

!     Declare integers or any extra variables you need here
!     INSERT
      integer :: i, x
      real :: ri, rj, percent, middle

!     Get the size of the mesh and store locally for convenience
      ni = g(1)%ni; nj = g(1)%nj;
      allocate(si(ni)); allocate(sj(nj));
      allocate(si1(ni)); allocate(sj1(nj));
      allocate(si2(ni)); allocate(sj2(nj));

!     Calculate the non-dimensional curve lengths of the geometry input and
!     generate linearly spaced points in i-direction at desired mesh resolution
      call dist(geom%x_a,geom%y_a,1,si_a)
      call dist(geom%x_b,geom%y_b,1,si_b)
!      
      
      if (geometric_mesh) then
!     Splitting up the x direction into 2 sections and then appending the meshes together
      percent = 0.3
      ri = 0.07
!     Need to add some offset so less small cells on the right side
      x = nint(percent * ni) + nint(0.1 * ni)
      call geometric(0.0, percent, x, 1.0-ri, si1)
!     Get the smallest interval of the first section
      middle = si1(x) - si1(x-1)
!     Need to ensure there is a space between the two appended arrays which is arbitrary
      call geometric(percent+middle, 1.0, ni-x, 1.0+ri, si2)
      si(:x) = si1(:x)
      si(x+1:) = si2(:ni-x)
      else
      call linspace(0.0,1.0,si)
      end if

!     Interpolate the geometry curves to the required resolution in the 
!     i-direction, this allows the mesh to be refined without altering the 
!     geometry definition file, the data is stored at "j = 1" and "j = nj"
      call interp(si_a,geom%x_a,si,g(1)%x(:,1))
      call interp(si_a,geom%y_a,si,g(1)%y(:,1))
      call interp(si_b,geom%x_b,si,g(1)%x(:,nj))
      call interp(si_b,geom%y_b,si,g(1)%y(:,nj))

!     Calculate the coordinates of all the intermediate points within the mesh.
!     Create a new vector of non-dimensional spacings in the j-direction using 
!     "linspace", loop over the mesh in the i-direction and calculate the
!     intermediate coordinates from a weighted sum of the two boundaries
!      call linspace(0.0,1.0,sj)
      
      if (geometric_mesh) then
!     Define the common ratio for the j direction geometric sequence
      rj = 1.06
      call geometric (0.0, 1.0, nj, rj, sj)
      else
      call linspace(0.0,1.0,sj)
      end if 
         
      do i=1, ni
!         Only considering the interior points as the outside points are defined by the input geometry
!	  Use a weighted sum of the values at each of the walls along with the current linspace value for proportions
      	  g(1)%x(i,2:nj-1) = (1.0 - sj(2:nj-1))*g(1)%x(i,1) + sj(2:nj-1)*g(1)%x(i,nj)
      	  g(1)%y(i,2:nj-1) = (1.0 - sj(2:nj-1))*g(1)%y(i,1) + sj(2:nj-1)*g(1)%y(i,nj)
      end do

!     In all of the test cases for the basic solver the the "j = 1" and "j = nj"
!     boundaries are walls, for the extensions you may need to return to this
!     and communicate the position of the walls to the solver in a more 
!     flexible way. The "wall" variable is an "ni * nj" logical array where 
!     "true" indicates the node is on a wall.
      g(1)%wall = .false.
      g(1)%wall(:,[1,g(1)%nj]) = .true.

!     Print that the mesh has been created
      write(6,*) 'Interpolated mesh from the bounding geometry curves'
      write(6,*)

      end subroutine generate_mesh


