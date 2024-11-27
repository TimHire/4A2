
      program solver

!     The main body of the CFD solver, calls all subroutines to march towards a
!     flow solution

!     Change to the directory you want to run your case within and execute with
!     "path_to_solver/solver.x < input_casename.txt" to run with screen output
!     or "path_to_solver/solver.x < input_casename.txt > log_casename.txt &" to
!     run in the background

!     Use derived data types defined in a separate module
      use types

!     Don't use historical implicit variable naming
      implicit none

!     Explicitly declare the required variables
      type(t_appvars) :: av
      type(t_bconds) :: bcs
      type(t_match) :: p
      type(t_geometry) :: geom
      type(t_grid) :: g
!	GETTING AN ERROR IF G IS NOT DEFINED AS A RAY BUT THIS SHOULD NOT BE RIGHT!!!!!!!!!!!!!!!!
      real :: d_max = 1, d_avg = 1
      integer :: nstep, nconv = 5, ncheck = 5
      integer :: nrkut, nrkuts = 4             ! set nkruts = 1 to turn of the Runge/Kutta improvement
      integer :: n
      logical :: geometric_mesh = .true.       ! to turn on/off geometric mesh improvement
      logical :: constant_enthalpy = .true.    ! to turn on/off the constant enthalpy improvement


!   THIS IS THE LAST CHANGE SO KNOW HAS BEEN SAVED


!     Read in the data on the run settings
      call read_settings(av,bcs)
      
!     Determine whether to generate the mesh within this Fortran program or read
!     it directly from a binary file written in Python
      if(av%ni /= -1) then
      !    av%nn = 1

!         Now the size of the grid is known, the space in memory can be 
!         allocated within the grid type
          call allocate_arrays(av,g,bcs)

!         Read in the case geometry
          call read_geom(av,geom)

!         Set up the mesh coordinates, interpolated between the geometry curves
          call generate_mesh(geom,g,geometric_mesh)

      else 
	 
!         Read the mesh coordinates directly from file - used for extension

          call read_mesh(av,g,bcs,p)

      end if
      
      !open(unit=15,file='g.txt')
      !write(15) g(1)
      
     ! close(15)
          
    !  do n=1,av%nn
      

!     Calculate cell areas and facet lengths
      call calc_areas(g)

!     Optional output call to inspect the mesh you have generated
      call write_output(av,g,n,1)

!     Check that the areas and projected lengths are correct
      call check_mesh(g)

!     Calculate the initial guess of the flowfield in the domain. There are two
!     options that can be chosen with the input argument "guesstype":
!         1. Uniform flow properties when "guesstype = 1", this is completed
!            for you already, it will allow you to get the solver started but
!            convergence will take more iterations.
!         2. A 1D varying flowfield when "guesstype = 2", assuming isentropic
!            flow in the i-direction allows a calculation of a better
!            approximation to the converged flowfield and so the time to
!            solution will be reduced. You will need to complete this option.
      call flow_guess(av,g,bcs,2)

!     Optional output call to inspect the initial guess of the flowfield
      call write_output(av,g,n,2)
      
      call set_timestep(av,g,bcs)
      
    !end do


!     Open file to store the convergence history. This is human readable during
!     a long run by using "tail -f conv_example.csv" in a terminal window
      open(unit=3,file='conv_' // av%casename // '.csv')

!     Initialise the "stopit" file, during long runs you can request an output
!     is written by setting the value to 1, or to terminate the calculation if
!     set to 2
      open(unit=11,file='stopit')
      write(11,*) 0; close(11);

!     Start the time stepping do loop for "nsteps". This is now the heart of the
!     program, you should aim to program anything inside this loop to operate as
!     efficiently as you can.
       
      do nstep = 1, av%nsteps
    !  do n=1,av%nn
!     Iterate through all the meshes and run an iteration one at a time before passing the information between

!         Smooth between the meshes



!         Update record of nstep to use in subroutines
          av%nstep = nstep

          ! For the Runge-Kutta improvement
          g%ro_start = g%ro
          g%roe_start = g%roe
          g%rovx_start = g%rovx
          g%rovy_start = g%rovy
          
          do nrkut = 1, nrkuts                  
          	  av%dt = av%dt_total / (1 + nrkuts - nrkut)
          	  
!          	  write(6,*) 'Current timestep of ', av%dt_total,'iterations'

	!         Calculate secondary flow variables used in conservation equations
		  call set_secondary(av,g,bcs,constant_enthalpy)

	!         Apply inlet and outlet values at the boundaries of the domain
		  call apply_bconds(av,g,bcs,constant_enthalpy)

	!         Perform the timestep to update the primary flow variables
		  call euler_iteration(av,g,constant_enthalpy)
		  
        
        		  
 	  end do

!         Write out summary every "nconv" steps and update "davg" and "dmax" 
          if(mod(av%nstep,nconv) == 0) then
              call check_conv(av,g,d_avg,d_max)
          end if

!         Check the solution hasn't diverged or a stop has been requested 
!         every "ncheck" steps
          if(mod(av%nstep,ncheck) == 0) then
              call check_stop(av,g)
          end if

!         Stop marching if converged to the desired tolerance "conlim"
          if(d_max < av%d_max .and. d_avg < av%d_avg) then
              write(6,*) 'Calculation converged in', nstep,'iterations'
              exit
          end if

     ! end do
      end do
      
      
    !  do n=1,av%nn

!     Calculation finished, call "write_output" to write the final, not 
!     necessarily converged flowfield
      write(6,*) 'Calculation completed after', av%nstep,'iterations'
      call write_output(av,g,n,3)
     ! end do
!
!     Close open convergence history file
      close(3)

      end program solver


