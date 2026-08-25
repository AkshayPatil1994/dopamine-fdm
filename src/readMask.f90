!> IBM setup: reads precomputed cell-centre SDF and builds ghost-cell stencil lists for U,V,W (active when ibm_input_mode==1, unused when 0)
Module ibmSetup

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use ibm
  Use monitor, Only : open_force_csv

  ! prevent implicit typing
  Implicit None

Contains

	! Load the precomputed SDF file and build ghost-cell lists.
	! Must be called after xg/yg/zg are initialised (grid setup complete).
	Subroutine readSDF

        	Call setup_ibm
        	If ( nsampling > 0 ) Call open_force_csv

	End Subroutine readSDF

End Module ibmSetup
