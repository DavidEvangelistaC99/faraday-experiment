      parameter(nl=16,nrange=75) ! was 95
      parameter(nstate=5,npen=5)
      real plag(nl,nrange+nl),plag_errors(nl,nrange+nl),sconst
      real densp(nrange+nl),tep(nrange+nl),trp(nrange+nl),
     & tip(nrange+nl),hfp(nrange+nl),hefp(nrange+nl),altp(nrange+nl),
     & r0,dr
      real edensp(nrange+nl),etep(nrange+nl),etip(nrange+nl),
     & ehfp(nrange+nl),ehefp(nrange+nl)
      common /fpa/densp,tep,trp,tip,hfp,hefp,altp,r0,dr,wl
      common /data/plag,plag_errors
      common /sys/sconst
      common /errs/edensp,etep,etip,ehfp,ehefp
