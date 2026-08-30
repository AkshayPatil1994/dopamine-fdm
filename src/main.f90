!			Main Dopamine Program
Program dopamine

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use input_output
  Use initialization
  Use time_integration
  Use monitor
  Use finalization
  Use reynolds_stress_budget
  Use probe_output
  Use synthetic_eddy_method, Only : accumulate_ti_rescale, apply_ti_rescale
  
  ! prevent implicit typing
  Implicit None

  ! initialize everything and read input file and input flow field
  Call initialize

  ! initialise Reynolds stress budget module (allocates accumulators)
  Call init_rsb

  ! initialise slice/line probe output (resolve grid indices, open files)
  Call init_probes

  ! small summary of input parameters
  Call summary
     
  ! temporal loop: nsteps>0 runs a fixed step count, nsteps<0 runs until t >= sim_end_time
  istep = 0
  Do

     istep = istep + 1

     ! time step (pressure forcing update handled inside )
     Call compute_time_step_RK3

     ! accumulate Reynolds stress statistics (no-op when rsb_active==0)
     Call accumulate_rsb

     ! write RSB window average when a full window has been accumulated
     If ( rsb_active == 1 .And. istep >= rsb_nstart .And. &
          n_accum > 0 .And. Mod(istep - rsb_nstart + rsb_freq, rsb_freq) == 0 ) Then
        Call output_rsb
     End If

     ! accumulate SEM inflow TI-rescale statistics (no-op when ti_rescale_active==0)
     Call accumulate_ti_rescale

     ! nudge the injected SEM target profile when a full window has been accumulated
     If ( ti_rescale_active == 1 .And. istep > ti_rescale_nstart .And. &
          Mod(istep - ti_rescale_nstart, ti_rescale_freq) == 0 ) Then
        Call apply_ti_rescale
     End If

     ! output some key values
     Call output_monitor

     ! write snapshot if needed
     Call output_data

     ! write 2-D slice and 1-D line probes (gated on slice_freq/line_freq)
     Call write_probes

     ! stopping criteria: step-count based (nsteps>0) or time based (nsteps<0)
     If ( nsteps > 0 ) Then
        If ( istep >= nsteps ) Exit
     Else
        ! tolerance guards against t landing one ULP short of sim_end_time after dt-snapping in compute_time_step_RK3
        If ( t >= sim_end_time - 1d-10 ) Exit
     End If

  End Do

  ! finalize stuff
  Call finalize_probes
  Call finalize_rsb
  Call finalize

End program dopamine
