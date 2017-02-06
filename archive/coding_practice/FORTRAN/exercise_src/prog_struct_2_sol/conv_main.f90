program dr

  Use numbers_module
  Use convert

  Implicit None

  Real( wp ) :: hour, min, sec, deg

  Call convert_time( 3661.0_wp, hour, min, sec )
  Write( *, * ) hour, min, sec
  deg = 90.0_wp
  Call deg_to_rad( deg )
  Write( *, * ) deg
  Call rad_to_deg( deg )
  Write( *, * ) deg


End program dr
