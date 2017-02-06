Program congugate_gradient

  use numbers_module
  Use base_data
  Use io
  Use init
  Use cg

  Implicit None

  Real( wp ) :: tolerance

  Integer :: n
  Integer :: iterations

  Call read_input( n, tolerance )
  Call allocate_arrays( n )
  Call initialize_arrays( a, b )
  Call generate_symm_pos_def( a )
  Call cg_solve( a, b, tolerance, iterations, x )
  Call write_results( a, b, iterations, x )
  Call deallocate_arrays

End Program congugate_gradient
