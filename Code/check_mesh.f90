      
      subroutine check_mesh(g,av)

!     Check the cell area and facet length calculations before attempting to
!     solve the flow, make sure you do this for both the "bump" and "bend" test
!     cases

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_grid), intent(inout) :: g(:)
      type(t_appvars), intent(in) :: av
      real :: area_min, dx_error, dy_error, tol, min_area, x_error, y_error
      integer :: ni, nj, n
      
      do n=1, av%nn

!     Get the size of the mesh and store locally for convenience
      ni = g(n)%ni; nj = g(n)%nj;

!     Exact checking of floating point numbers never goes well, define a
!     small tolerance to use for a comparative operation instead
      tol = 1e-4 * g(n)%l_min

!     Check that all of the cell areas are positive, either with the intrinsic
!     "minval" function or with nested do loops. Print the output to the screen
!     and flag negative numbers as an error with an if statement to "stop" the
!     program
!     INSERT
      min_area = minval(g(n)%area)
      if (min_area < 0) then
          write(6,*) 'Negative cell area found, stopping program'
          stop
      end if 
      	  

!     Next check that the sum of the edge vectors around every quadrilateral is 
!     very nearly zero in both the x and y-coordinate directions. You can
!     complete this with some elementwise addition of the arrays and use of the
!     "maxval" and "abs" intrinsic functions.
!     INSERT
!     Consider the absolute difference in the lengths of the x and y sides separately
      x_error = maxval(abs(g(n)%lx_i(1:ni-1,1:nj-1) + g(n)%lx_j(1:ni-1,1:nj-1) - g(n)%lx_i(2:ni,1:nj-1) - g(n)%lx_j(1:ni-1,2:nj)))
      y_error = maxval(abs(g(n)%ly_i(1:ni-1,1:nj-1) + g(n)%ly_j(1:ni-1,1:nj-1) - g(n)%ly_i(2:ni,1:nj-1) - g(n)%ly_j(1:ni-1,2:nj)))
      
      if (x_error > tol .or. y_error > tol) then
          write(6,*) 'Error of', max(x_error, y_error), 'has been found which is greater than the tolerance of', tol
          stop
      end if
      	

!     It may be worthwhile to complete some other checks, the prevous call to
!     the "write_output" subroutine has written a file that you can read and
!     postprocess using the Python script plot_mesh.py. This program also has
!     access to all of the mesh parameters used within the solver that you could
!     inspect graphically.


      end do
!     Print a blank line
      write(6,*)
      

      end subroutine check_mesh
