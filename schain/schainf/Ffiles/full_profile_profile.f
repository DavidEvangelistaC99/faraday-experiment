c
c     $Id: full_profile.f,v 1.3 2013/09/12 19:11:48 daveh Exp daveh $ DLH 
c      

      subroutine profile(acf_sum,acf_err,power,en,alag,thb2,bfm2,
     & ote,ete,oti,eti,oph,eph,ophe,ephe,range2,ut,
     &     NHTS,NACF,IBITS,acf_avg_real,status)
c
      complex acf_sum(4,NHTS,IBITS)
      real acf_avg_real(NHTS,IBITS)
      real acf_err(NHTS,IBITS),range2(NHTS)
      real power(NHTS),en(NHTS),alag(IBITS),thb2(NHTS),bfm2(NHTS)
      real ote(NHTS),ete(NHTS),oti(NHTS),eti(NHTS),oph(NHTS),eph(NHTS)
      real ophe(NHTS),ephe(NHTS)
      integer NHTS,NACF,IBITS
      real status(1)
      common /chisq/chi2
c
c     arguments above, internal variables below: have to make these match
c
      include 'fpa.h'
      include 'bfield.h'
c
      common /mode/imode
      common /utime/uttime
c
      include 'spline.h'
      include 'fitter.h'
c
c     functions to minimize are all lagged profiles plus penalties
c     variables to minimize them with are all spline values
c

      parameter(nfun=nrange*nl+nstate*(nspline-2)+npen)
      parameter(nvar=nspline*nstate)
      parameter(lwa=nfun*nvar+nstate*nvar+nfun)
      real fvec(nfun),wa(lwa),xe(nvar)
      integer iwa(nvar)
c      integer lwa
c
      external fcn_lpreg
c
      pi=4.0*atan(1.0)
      uttime=ut
c      write(*,*) "ut: ",ut
c      write(*,*) "NHTS: ",NHTS
c      write(*,*) "NACF: ",NACF
c      write(*,*) "IBITS: ",IBITS
c      call exit
c      write(*,*) "acf_sum: ", acf_sum(1,111,11)      
c      write(*,*) "Starting Profile Process"
c     set up basic physical parameters
c      dens=639747.44
c      dens=509664.781
c      te=3971.98926
c      write(*,*) "acf_sum: ", acf_sum(1,31,2)
c      call    exit

      wl=3.0
      ak=2.0*pi/wl

c     analysis starts at gate 15

      ir0=15

      r0=range2(ir0)
      dr=range2(ir0+1)-range2(ir0)

      do i=1,nfield ! nrange+10   memory hole
         bfld_prof(i)=bfm2(i+ir0-1)
         alpha_prof(i)=thb2(i+ir0-1)*pi/180.0
      end do

      do i=1,nrange+nl ! careful ... check nfield compared to nrange+nl
         densp(i)=en(i+ir0-1)
         do j=1,nl
c            write(*,*) acf_sum(1,i+ir0-1,j)
            plag(j,i)=acf_sum(1,i+ir0-1,j)
            
            plag_errors(j,i)=acf_err(i+ir0-1,j)
c            write(*,*) j,i,plag_errors(j,i)
         end do
      end do      
c      call exit
c     set some fitting conventions
c      write(*,*) "plag: ",plag(2,2)
      imode=2
      nion=3

      ven=0.0
      do i=1,nion
         vin(i)=0.0
      end do
      wi(1)=16
      wi(2)=1
      wi(3)=4

c     system constant

      sconst=1.0

c     set up splines: DLH need to modify for different nu. splines

      do i=1,nspline+norder
         ta(i)=r0+float(i-1-2)*delr ! 2-knot offset for 4th order spline
      end do

c     initial guess - set density, He+ profiles here
         
      do i=1,nspline ! DLH modify for different no. splines 
         bcoef(i,1)=log10(abs(densp(1+int(float(i-1)*delr/dr))))
                                ! need good guess for density profile 
         tmp=0.03
         bcoef(i,5)=atanh(2.0*tmp-1.0) ! not so critical
      end do        

c      write(*,*) bcoef
c      call exit

c     grid search for other params to help set initial guess
c     could include Ne, He+ in grid search too
c      write(*,fmt='("before GRID")')
      call grid()
c      write(*,fmt='("After GRID")')
      write(*,*) "chi2 after grid:",chi2
c      write(*,fmt='("Just After chi")')
c      call exit
      if(chi2.le.3.5) then

c      write(*,fmt='("After chi")')
c      call exit
      tol=1.0e-9

c      call lmdif1(fcn_lpreg,nfun-nl*nl,nvar,bcoef,fvec,tol,
c     &     info,iwa,wa,lwa)

      info=0
      m=nfun-nl*nl
      n=nvar
      zero=0.0
      factor=1.0e-2   ! rather small initial step size

      maxfev = 200*(n + 1)
      ftol = tol
      xtol = tol
      gtol = zero
      epsfcn = zero
      mode = 1
      nprint = 0
      mp5n = m + 5*n
c      write(*,*) wa(5*n+1)
c      call exit
      write(*,*) "before densp:", densp(1),bcoef(1,2)
      call lmdif(fcn_lpreg,m,n,bcoef,fvec,ftol,xtol,gtol,
     &     maxfev,epsfcn,wa(1),
     *           mode,factor,nprint,info,nfev,wa(mp5n+1),m,iwa,
     *           wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))

c      write(*,fmt='("After lmdif")')
c      write(*,*) info
      write(*,*) "densp_after:", densp(1),bcoef(1,2)
c      call exit

c
c     error handling; evalate solution point first, then Jacobian
c
      iflag=0
c      write(*,fmt='("before fcn_lpreg")')
      call fcn_lpreg(m,n,bcoef,fvec,iflag) ! common sys constant throughout
c      write(*,fmt='("After fcn_lpreg")')
      nm1=n*m+1
      call fdjac2(fcn_lpreg,m,n,bcoef,fvec,wa,m,iflag,epsfcn,wa(nm1))

      do i=1,n
         err=0.0
         xe(i)=0.0
         do j=1,m-npen
c         do j=2+(nspline-2)*nstate,m-npen ! skip equations unrelated to data
            err=err+MIN(wa(j+(i-1)*m)**2,10.0)
         end do
         if(err.gt.0.0) xe(i)=sqrt(err**(-1))
      end do
    
      call propagate(xe)

      do i=1,nrange
         ete(i+ir0)=etep(i)
         eti(i+ir0)=etip(i)
         eph(i+ir0)=ehfp(i)
         ephe(i+ir0)=ehefp(i)
      end do
      
c
c     output parameters
c

      do i=1,nrange
         altp(i)=r0+dr*float(i-1)
         call get_spline(altp(i),densp(i),tep(i),tip(i),
     &        hfp(i),hefp(i))
         ote(i+ir0)=MAX(tep(i),tip(i))
         oti(i+ir0)=tip(i)
         oph(i+ir0)=hfp(i)
         ophe(i+ir0)=hefp(i)
      end do
c
c     theoretical lag products and other quantities for plotting
c

      call lagp(plag,wl,r0,dr,nl,nrange)
      write(*,*) "plag_plag",REAL(plag(12,58))
      write(*,*) "plag_plag",REAL(plag(12,72))
c      write(*,*) "plag_plag",AIMAG(plag(12,72))
      do i=1,nrange                        !
         en(i+ir0-1)=densp(i)
         do j=1,nl
            acf_avg_real(i+ir0-1,j)=REAL(plag(j,i))
c            acf_avg_imag(i+ir0-1,j)=AIMAG(acf_sum(1,i+ir0-1,j))
         end do
      end do      
      write(*,*) "acf",acf_avg_real(72,12)
c      write(*,*) "acf",acf_avg_imag(72,12)
      end if
      status(1)=chi2
      
c      write(*,*) "status: ",status

 5    return
      end


      subroutine grid()
c
c     minimize chi-squared statistic through grid search
c     just optmize temperature and H+ profiles
c
      include 'spline.h'
      include 'fpa.h'
      include 'fitter.h'
c
      real plag2(nl,nrange)
      real fvec(nrange*nl)
      real chmin,uttime
      integer lmin,mmin,nmin
      integer iloop,iskip,l1,l2,m1,m2,n1,n2
c
      common /chisq/chi2
      common /utime/uttime
c
      chmin=1.0e10

      do i=1,nrange
         altp(i)=r0+dr*float(i-1)
      end do

      goto 3
      call newplot(.true.)
      call graphinit()
      call erase

      call setorigin(6.5,3.5,90.0)
      call linearaxis("L",altp(1),altp(nrange),1.5,5.0,2,0.0,50.0)
      call linearaxis("B",0.0,5000.0,1.5,5.0,2,0.0,500.0)
      call drawaxis("T",2)
      call drawaxis("R",2)
      call labelaxis("B","(f5.0)",.1,2,0.0,2000.0)
      call labelaxis("L","(f5.0)",.1,2,0.0,100.0)
      call ctitle("L","Altitude (km)",.125,2)
      call ctitle("B","T\\se\\u (K)",.125,2)

      call setorigin(6.5,7.5,90.0)
      call linearaxis("L",altp(1),altp(nrange),1.5,5.0,2,0.0,50.0)
      call linearaxis("B",0.0,1.0,1.5,5.0,2,0.0,0.25)
      call drawaxis("T",2)
      call drawaxis("R",2)
      call labelaxis("B","(f5.1)",.1,2,0.0,0.5)
      call labelaxis("L","(f5.0)",.1,2,0.0,100.0)
      call ctitle("B","H\\S+\\U frac.",.125,2)

 3    continue

      do iloop=1,2
         if(iloop.eq.1) then
            l1=2
            l2=26
            m1=2
            m2=10
            n1=2
            n2=16
            iskip=2
         else
            l1=lmin-2
            l2=lmin+2
            m1=mmin-2
            m2=mmin+2
            n1=nmin-2
            n2=nmin+2
            iskip=1
         end if

      do l=l1,l2,iskip ! temperature search
         do k=1,nspline

            tmp=-2.5+float(k-1)*float(l-1)*1.3e-2 
c            if(uttime.ge.2.0.and.uttime.lt.10.0) then
c               tmp=-1.5+float(k-1)*float(l-1)*8.0e-3 
c            else
c               tmp = -0.75 -0.2*float(l) + 3.5e-1*float(k-1) ! was 3.5
c               tmp=-1.5 + 3.0e-2*float(l-1) 
c     &              *(float(k-1)**3.0/float(nspline-1)**2.0)
c            end if

            bcoef(k,2)=tmp
            if(uttime.ge.10.0.and.uttime.le.14.0) then ! was 10.5
               bcoef(k,3)=tmp-0.1 
            else
               bcoef(k,3)=tmp-0.02
            end if
         end do

         do m=m1,m2,iskip ! H+ scale height search DLH 2015
            
            do n=n1,n2,iskip ! H+ crossing  height search
	      
               do k=1,nspline
                  tmp=float(k-nspline/4-n)*(float(m-1)*5.0e-2+0.1)
                  bcoef(k,4)=tmp
               end do
               
               do i=1,nrange
                  call get_spline(altp(i),densp(i),tep(i),tip(i),
     &                 hfp(i),hefp(i))

               end do
	       
               goto 4
               call setorigin(6.5,3.5,90.0)
               call setscale("L",altp(1),altp(nrange),1)
               call setscale("B",0.0,5000.0,1)
               call setclip("L",250.0,1500.0)
               call plotline("B,L",tep,altp,nrange,2,0)
              
               if(l.eq.1) then
                  call setorigin(6.5,7.5,90.0)
                  call setscale("L",altp(1),altp(nrange),1)
                  call setscale("B",0.0,1.0,1)
                  call setclip("L",250.0,1500.0)
                  call plotline("B,L",hfp,altp,nrange,2,0)
               end if
 4             continue
	       
c     calculate chi-2
c               write(*,*) "wl: ",wl, "r0: ",r0
c               write(*,*) "dr: ",dr,"nl: ",n,"nrange", nrange
c	       write(*,*) "before lagp inside grid"
               call lagp(plag2,wl,r0,dr,nl,nrange)
c               write(*,*) "plag2: ",plag2(12,72)
c               write(*,*) "plag2: ",plag2
c               call exit	       
               call get_scale(plag2)
c               write(*,*) "plag: ",plag(13,70)
c               write(*,*) "plag2: ",plag2(13,70)
c	       write(*,*) "plag_errors: ",plag_errors(13,70)
c	       write(*,*) "sconst: ",sconst              
c               call exit
c               write(*,*) "cal chi2.1: ", chi2
c               call exit
               chi2=0.0
               ipt=0
c               write(*,*) "cal chi2.1: ", chi2
               do j=nl+1,nrange
                  do i=1,nl
                     ipt=ipt+1
                     fvec(ipt)=(plag(i,j)-plag2(i,j)*sconst)/
     &                    plag_errors(i,j)
c                     write(*,*) j,i,plag(i,j),plag2(i,j),
c     &                    plag_errors(i,j),sconst
c		     call exit
                     chi2=chi2+fvec(ipt)**2
c                     write(*,*) "cal chi2.1: ", chi2
                  end do
c                  call exit
               end do
c               write(*,*) "cal chi2.2: ", chi2
c               call exit
               chi2=sqrt(chi2/float(nrange*nl))!+float(l)*0.005 ! DLH??
               if(chi2.lt.chmin) then
                  lmin=l
                  mmin=m
                  nmin=n
                  chmin=chi2
               end if

               write(*,*) "A",l,m,n,chi2,chmin,"B" ! best
c               call exit
c	        write(*,fmt='("intermediate")')
            end do
         end do
      end do

      end do
c      write(*,fmt='("befire flush")')
c      call where()
      write(22,*) lmin,mmin,nmin,chmin
      call flush(22)
c      write(*,fmt='("after flush")')
c     best result here

      do k=1,nspline
         tmp=-2.5+float(k-1)*float(lmin-1)*1.3e-2 
c         if(uttime.ge.2.0.and.uttime.lt.10.0) then
c         tmp=-1.5+float(k-1)*float(lmin-1)*8.0e-3 
c         else
c            tmp = -0.75 -0.2*float(l) + 3.5e-1*float(k-1) ! was 3.5
c         tmp=-1.5 + 3.0e-2*float(lmin-1) ! new family of temps 
c     &     *(float(k-1)**3.0/float(nspline-1)**2.0)
c         end if
         bcoef(k,2)=tmp
         if(uttime.ge.10.0.and.uttime.le.14.0) then  ! was 10.5
            bcoef(k,3)=tmp-0.1
         else
            bcoef(k,3)=tmp-0.02
         end if
         tmp=float(k-nspline/4-nmin)*(float(mmin-1)*5.0e-2+0.1)
         bcoef(k,4)=tmp
      end do

c     reset system constant

      do i=1,nrange
         call get_spline(altp(i),densp(i),tep(i),tip(i),
     &        hfp(i),hefp(i))
      end do
c      write(*,fmt='("------------before lagp------------")')
c      write(*,*) chi2
c      call exit
      call lagp(plag2,wl,r0,dr,nl,nrange)   
c      write(*,fmt='("after lagp")')      
      call get_scale(plag2)
c      write(*,*) chi2
      chi2=chmin
      return
      end

      subroutine propagate(xe)
c
c     propagate errors through splines
c
      include 'fpa.h'
      include 'spline.h'
      parameter(nvar=nspline*nstate)
      real xe(nvar)


      do i=1,nrange

         write(*,*) "densp:", densp(i),bcoef(1,2)
c         call exit
c                  
         edensp(i)=bvalue(ta,xe(1+0*nspline),nspline,norder,altp(i),0)
     &        *densp(i)*2.30
         x=bvalue(ta,bcoef(1,2),nspline,norder,altp(i),0)
         etep(i)=bvalue(ta,xe(1+1*nspline),nspline,norder,altp(i),0)
     &        *(t1/2.0)/cosh(x)**2 ! DLH 10/14 was 3500

         write(*,*) i,altp(i),x,
     &        bvalue(ta,xe(1+1*nspline),nspline,norder,altp(i),0),
     &        (3500.0/2.0)/cosh(x)**2,xe(1+i/4+1*nspline)

         x=bvalue(ta,bcoef(1,3),nspline,norder,altp(i),0)
         etip(i)=bvalue(ta,xe(1+2*nspline),nspline,norder,altp(i),0)
     &        *(t1/2.0)/cosh(x)**2  ! DLH 10/14 was 3500
         x=bvalue(ta,bcoef(1,4),nspline,norder,altp(i),0)
         ehfp(i)=bvalue(ta,xe(1+3*nspline),nspline,norder,altp(i),0)
     &        *(1.0/2.0)/cosh(x)**2
         x=bvalue(ta,bcoef(1,5),nspline,norder,altp(i),0)
         ehefp(i)=bvalue(ta,xe(1+4*nspline),nspline,norder,altp(i),0)
     &        *(1.0/2.0)/cosh(x)**2

         write(*,*) altp(i),edensp(i),etep(i),etip(i),ehfp(i),ehefp(i)
c         call exit
      end do

c      stop

      return
      end

      subroutine fcn_lpreg(m,n,x,fvec,iflag)
c
      include 'spline.h'
      include 'fpa.h'
      include 'fitter.h'
c
      integer m,n,iflag,iflaglast
      real x(n),fvec(m),pen(nstate+npen)
      real plag2(nl,nrange)
      real regu(nstate+npen)
      data regu/1.0e1,5.0e1,10.0e1,1.0e2,5.0e1,
     &     1.0e1,2.0e0,1.0e0,5.0e-3,1.0e0/ 

c     2019 
      common /iflaglast/iflaglast
      common /chisq/chi2
      common /utime/uttime
c
c      write(*,*)uttime," 3 *"

c     zero out, copy over some arrays

c      write(*,*) "Starting fcn_lpreg"
      write(*,*) "x:",x(2)

      do j=1,nstate+npen
         pen(j)=0.0
      end do

      ipt=0
      do j=1,nstate
         do i=1,nspline
            ipt=ipt+1
            bcoef(i,j)=x(ipt)
c            write(*,*)j,i,bcoef(i,j)
         end do
      end do

c     populate profiles from splines
      write(*,*) "densp1:",densp(1)
      do i=1,nrange
         altp(i)=r0+dr*float(i-1)
c         write(*,*) "densp1:",densp(i)
         call get_spline(altp(i),densp(i),tep(i),tip(i),
     &        hfp(i),hefp(i))
c         write(*,*) "densp2:",densp(i)
c         call exit
      end do
      write(*,*) "densp2:",densp(1)
c     regularize spline coefficients, establish boundary conditions
c      write(*,*) "After spline_fcn_lpreg"
      ipt=0
      do i=2,nspline-1
         do j=1,nstate

            fac=1.2  ! DLH 7/14 could be bendy (or else get spurious He+)
            if((uttime.lt.10.0.or.uttime.gt.15.0).and.j.eq.2) fac=4.0
            
            ipt=ipt+1
            tmp=(bcoef(i-1,j)-2.*bcoef(i,j)+bcoef(i+1,j))
            fvec(ipt)=tmp*regu(j)*fac
            pen(j)=pen(j)+fvec(ipt)**2
         end do
      end do

c     rescale for clearer diagnostic (just for printing)

      do i=1,nstate
         pen(i)=sqrt(pen(i)/float(nspline)) 
      end do
       
c      write(*,*) x

c     compute model prediction error to minimize in LS sense
c      write(*,*) "before lagp_fcn_lpreg"
      call lagp(plag2,wl,r0,dr,nl,nrange)
c      write(*,*) "before lagp_fcn_lpreg"
c     adjust system constant
      write(*,*) "flgas:",iflag,iflaglast
      if(iflag.eq.1.and.iflaglast.eq.2) then
         sclast=sconst
         call get_scale(plag2)
         sconst=0.1*sconst+0.9*sclast
      end if
      iflaglast=iflag
         
      chi2=0.0
      do j=nl+1,nrange
         do i=1,nl
            ipt=ipt+1
            fvec(ipt)=(plag(i,j)-plag2(i,j)*sconst)/plag_errors(i,j)
            chi2=chi2+fvec(ipt)**2
         end do
      end do
      chi2=sqrt(chi2/float(nrange*nl))

c     add penalties for Tr, composition, He+, lower boundary values, etc.
c     have to adjust penalties for different conditions (e.g. sunrise, night)

      do i=1,nrange
         pen(nstate+1)=pen(nstate+1)+hefp(i)**2
         
         if(tep(i).gt.tip(i)) then ! 3/2015 500 -> 1000 ?
            if(uttime.ge.10.0.and.uttime.le.15.0) then ! was 10.5, 14->15
               scale=500.0  ! dlh 7/14 was 1500, but this defeated He+, so 500
            else           
               scale=50.0  
            end if
         else
            scale=10.0           
         end if 

         pen(nstate+2)=pen(nstate+2)+((tip(i)-tep(i))/scale)**2
         pen(nstate+3)=pen(nstate+3)+exp((hfp(i)+hefp(i)-1.0)/2.0)**2

      end do

      pen(nstate+4)=(tip(1)-700.0)**2*float(nrange) ! or 800?
      pen(nstate+5)=(hfp(nrange)-1.0)**2*float(nrange)

      do i=1,npen
         tmp=sqrt(pen(nstate+i))*regu(nstate+i)
         ipt=ipt+1
         fvec(ipt)=tmp
         pen(nstate+i)=tmp/sqrt(float(nrange))
      end do
      
c      write(*,*) "chi2:",chi2," sys const:",sconst,
c     &     " flag:",iflag," penalties:"
c      write(*,*)(pen(i),i=1,5)
c      write(*,*)(pen(i),i=6,10)
c      write(*,*) "ut:",uttime
c      write(*,*) fvec(222)

      return
      end

      subroutine get_scale(plag2)     
c
      include 'fpa.h'
      real plag2(nl,nrange)
c      
      a=0.
      b=0.
c      do j=3*nrange/4,nrange
      do j=nl+1,nrange 
c         write(*,*) j,plag(1,j),plag2(1,j)
         a=a+plag(1,j)*plag2(1,j)/plag_errors(1,j)**2
         b=b+plag(1,j)*plag(1,j)/plag_errors(1,j)**2
      end do

      sconst=b/a
      
      write(*,*)"sconst_inside_get_scape:",sconst,a,b
c      stop
      
      return
      end
