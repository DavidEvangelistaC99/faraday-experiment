      parameter(nspline=30) ! limited to 40 at present (see below) (was 29)
      parameter(norder=4)
      parameter(delr=39.0) ! this times nspline ~ top altitude
      parameter(t0=800.0)
      parameter(t1=3000.0) ! was 3750.0
      real ta(nspline+norder),bcoef(nspline,5)
      common/spline/ta,bcoef
