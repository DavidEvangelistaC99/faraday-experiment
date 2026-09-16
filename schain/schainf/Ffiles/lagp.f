      subroutine lagp(plag,wl,r0,dr,nl,nrange)
c
c     predict lagged product matrix using state parameters stored
c     as cubic b-splines
c     read weights for numeric integration from a file
c
c     faster version - assumes samples taken at half-integer lags, ranges
c
      include 'fitter.h'
      include 'spline.h'
      include 'bfield.h'
      integer nl,nrange
      real plag(nl,nrange),taus(nl)
      parameter(nw=1000)
      integer ilag(nw)
      real weight(nw),dtau(nw),drange(nw)
      logical first
      real store(2*nl,2*nrange,2)

      character(1024) :: fqual_temp
      character(:), allocatable :: fqual
      character bfmodel_str*8

      data first/.true./
c
      
c      write(*,*) "dens_before: ",dens
c      call exit
c     zero out array
c      write(*,*) "Starting lagp"
      do i=1,nrange
         do j=1,nl
            plag(j,i)=0.0
         end do
      end do

c
c     compute all required acfs
c
      do i=1,2*nrange
         ii=(i-1)/2+1
         alt=r0+float(i-1)*dr*0.5  ! half integer ranges
c         write(*,*) "dens_before: ",dens
c         call exit
c         write(*,*) "alt: ",alt
c	 write(*,*) "dens: ",dens 
c	 write(*,*) "te: ",te
c	 write(*,*) "ti1: ",ti1
c	 write(*,*) "hf: ",hf
c	 write(*,*) "hef: ",hef
c         call exit
         call get_spline(alt,dens,te,ti1,hf,hef)
c         write(*,*) "dens_after: ",dens
c         call exit
         do k=1,nion
            ti(k)=ti1
         end do
         fi(2)=hf
         fi(3)=hef
         fi(1)=1.0-hf-hef         

         do j=1,2*nl
            tau=float(j-1)*(dr/1.5e5)*0.5  ! half integer lags
            
c     just consider a single azimuth angle for now

c            call acf2(wl,tau,te,ti,fi,ven,vin,wi,nion,
c     &           alpha_prof(ii)+1.74e-2,dens,bfld_prof(ii),rho)       
c            store(j,i,1)=rho*dens*(100.0/alt)**2

c            call acf2(wl,tau,te,ti,fi,ven,vin,wi,nion,
c     &           alpha_prof(ii)-1.74e-2,dens,bfld_prof(ii),rho)       
c            store(j,i,2)=rho*dens*(100.0/alt)**2

            call acf2(wl,tau,te,ti,fi,ven,vin,wi,nion,
     &           alpha_prof(ii),dens,bfld_prof(ii),rho)       
            store(j,i,1)=rho*dens*(100.0/alt)**2
c            write(*,*) "store", store(j,i,1)
c            call exit
         end do
      end do

c     construct lagged products
c
c      write(*,*) "before weights"
c      write(*,*) first 
c      write(*,*) "weight before: ",weight
c      write(*,*) "BEFORE GET_PATH"
      call get_path(fqual_temp)
c      write(*,*) "L_BEF: ", fqual_temp, "L_BEF_end"
      fqual = TRIM(fqual_temp)
      bfmodel_str = 'bfmodel/'
      write(*,*) "bfmodel_str: ", bfmodel_str
c      fqual(1:len_trim(fqual)-len_trim(bfmodel_str))
      fqual = fqual(1:index(fqual, bfmodel_str)-1)//'weights.dat'
      write(*,*) "fqual: ", fqual
c      write(*,*) "Final fqual",fqual,"FINAL"
c      call exit
      if(first) then
         open(unit=25,file=fqual,
     &    status='old')
         read (25,*) nwt
         do i=1,nwt
            read(25,*) weight(i),dtau(i),drange(i),ilag(i)
c            write(*,*) "ilag: ",ilag(i)
         end do
c         flush(25)
         close(25)
c         first=.false.
      end if

c      write(*,*) "after weights"
 
c      write(*,*) "weight(3): ",weight(3)
c      call exit
c      write(*,*) "nl: ",nl,"nrange: ",nrange,"nwt: ",nwt

      do i=1+nl,nrange  ! don't worry about redundancy here
         
         do j=1,nwt
c            write(*,*) "j: ",j
            alt=r0+(float(i-1)+drange(j))*dr
            tau=abs(dtau(j)*dr/1.5e5)
            lag=ilag(j)+1
c            write(*,*) tau

            k=1+2.0*(tau/(dr/1.5e5))
            l=1+2.0*(alt-r0)/dr

            do m=1,1 ! 2
               
c               call exit
c               write(*,*) i
c               call exit
c               if (i .eq. nrange) then
c               	  if (j .eq. 227) then

c              	   write(*,*) "lag: ",lag
c              	   write(*,*) "ilag: ",ilag(j)
c               call exit
c                   write(*,*) "plag: ",plag(lag,i)
c                   write(*,*) "weight: ",weight
c
c                   write(*,*) "weight: ",weight(j)
c                   write(*,*) "store: ",store(k,l,m)
c                   call exit
c	       	  end if
c               end if
               plag(lag,i)=plag(lag,i)+weight(j)*store(k,l,m)
c               if (i .eq. nrange) then
c	       write(*,*) "lag: ",lag				
c               write(*,*) "plag: ",plag(lag,i)
c                write(*,*) "plag: ",plag(lag,i)
c                call exit
c               end if
            end do            
c            write(*,*) "plag: ",plag(lag,i)
c             write(*,*) "plag: ",plag(12,72) 
c            call exit
         end do
c         call exit
c      write(*,*) "plag: ",plag(12,72) 
      end do
c      write(*,*) "plag: ",plag(12,72)
c      call exit
c      write(*,*) "End LAGP"
      return
      end

      subroutine lagp_old(plag,wl,r0,dr,nl,nrange)
c
c     predict lagged product matrix using state parameters stored
c     as cubic b-splines
c     read weights for numeric integration from a file
c
c     general version - samples can be anywhere
c
      include 'fitter.h'
      include 'spline.h'
      include 'bfield.h'
      integer nl,nrange
      real plag(nl,nrange),taus(nl)
      parameter(nw=1000)
      integer ilag(nw)
      real weight(nw),dtau(nw),drange(nw)
      logical first
      data first/.true./
c

c     zero out array

      do i=1,nrange
         do j=1,nl
            plag(j,i)=0.0
         end do
      end do

c
c     construct lagged products
c

      if(first) then
         open(unit=25,file='weights.dat',status='old')
         read (25,*) nwt
         do i=1,nwt
            read(25,*) weight(i),dtau(i),drange(i),ilag(i)
         end do
         close(25)
         first=.false.
      end if
      
c     can avoid recalculating splines and/or ACFs when parameters repeat

      altp=-1.0
      taup=-1.0

      do i=1+nl,nrange  
         do j=1,nwt
            alt=r0+(float(i-1)+drange(j))*dr
            ii=(alt-r0)/dr+1

            if(alt.ne.altp) then
               call get_spline(alt,dens,te,ti1,hf,hef)

c     write(*,*) alt,dens,te,ti1,hf,hef

               do k=1,nion
                  ti(k)=ti1
               end do
               fi(2)=hf
               fi(3)=hef
               fi(1)=1.0-hf-hef
            end if

            tau=abs(dtau(j)*dr/1.5e5)
            lag=ilag(j)+1

c            write(*,*) i,alpha_prof(i),bfld_prof(i)

            if(tau.ne.taup.or.alt.ne.altp) then

c     two- or three- point Gauss Hermite quadrature rule

c               call acf2(wl,tau,te,ti,fi,ven,vin,wi,nion,
c     &              alpha_prof(ii)+1.74e-2,dens,bfld_prof(ii),rho1)

c               call acf2(wl,tau,te,ti,fi,ven,vin,wi,nion,
c     &              alpha_prof(ii)-1.74e-2,dens,bfld_prof(ii),rho2)

c               call acf2(wl,tau,te,ti,fi,ven,vin,wi,nion,
c     &              alpha_prof(ii),dens,bfld_prof(ii),rho1)

c               call acf2(wl,tau,te,ti,fi,ven,vin,wi,nion,
c     &              alpha_prof(ii)+1.9e-2,dens,bfld_prof(ii),rho2)

c               call acf2(wl,tau,te,ti,fi,ven,vin,wi,nion,
c     &              alpha_prof(ii)-1.9e-2,dens,bfld_prof(ii),rho3)

               call acf2(wl,tau,te,ti,fi,ven,vin,wi,nion,
     &              alpha_prof(ii),dens,bfld_prof(ii),rho)

c               rho=(rho1+rho2)/2.0
c               rho=(rho2+rho3)*0.29541+rho1*1.18164)/1.77246              

            end if

            plag(lag,i)=plag(lag,i)+weight(j)*rho*dens*(100.0/alt)**2

            altp=alt
            taup=tau

         end do
      end do

      return
      end

      function atanh(x)
c
      real atanh,x
c
      atanh=log(sqrt((1.+x)/(1.-x)))
      return
      end

      subroutine get_spline(alt,dens,te,ti,hf,hef)
c
c     routines for handling cubic b-spline interpolation
c     gets values consistent with stored coefficients
c

      include 'spline.h'

c
c     bspline values from 200-1500 km in 15-km intervals
c     must specify five knots above ceiling
c     note offset ... start accessing splines above bottom two knots
c
c     five banks of splines, for Ne, Te, Ti, H+, and He+ versus altitude
c     need to initialize ta somewhere
c
c      write(*,*) "ta: ",ta
c      write(*,*) "bcoef(1,1): ",bcoef(1,1)
c      write(*,*) "nspline: ",nspline
c      write(*,*) "norder: ",norder
c      write(*,*) "alt: ",alt
      dens=bvalue(ta,bcoef(1,1),nspline,norder,alt,0)
c      write(*,*) "dens_inside_get_spline: ",dens
c      write(*,*) dens
c      call exit
      dens=10.0**MAX(dens,2.0)
c      write(*,*) dens
c      call exit
      te=bvalue(ta,bcoef(1,2),nspline,norder,alt,0)
      te=t0+t1*(1.0+tanh(te))/2.0  ! DLH 10/14 was 3500
c      tr=bvalue(ta,bcoef(1,3),nspline,norder,alt,0)
c      tr=exp(tr)
c      ti=te/tr
      ti=bvalue(ta,bcoef(1,3),nspline,norder,alt,0)
      ti=t0+t1*(1.0+tanh(ti))/2.0 ! DLH 10/14 was 3500
      hf=bvalue(ta,bcoef(1,4),nspline,norder,alt,0)
      hf=(1.0+tanh(hf))/2.0
      hef=bvalue(ta,bcoef(1,5),nspline,norder,alt,0)
      hef=(1.0+tanh(hef))/2.0
      return
      end


      real function bvalue ( t, bcoef, n, k, x, jderiv )
c  from  * a practical guide to splines *  by c. de boor    
calls  interv
c
calculates value at  x  of  jderiv-th derivative of spline from b-repr.
c  the spline is taken to be continuous from the right, EXCEPT at the
c  rightmost knot, where it is taken to be continuous from the left.
c
c******  i n p u t ******
c  t, bcoef, n, k......forms the b-representation of the spline  f  to
c        be evaluated. specifically,
c  t.....knot sequence, of length  n+k, assumed nondecreasing.
c  bcoef.....b-coefficient sequence, of length  n .
c  n.....length of  bcoef  and dimension of spline(k,t),
c        a s s u m e d  positive .
c  k.....order of the spline .
c
c  w a r n i n g . . .   the restriction  k .le. kmax (=20)  is imposed
c        arbitrarily by the dimension statement for  aj, dl, dr  below,
c        but is  n o w h e r e  c h e c k e d  for.
c
c  x.....the point at which to evaluate .
c  jderiv.....integer giving the order of the derivative to be evaluated
c        a s s u m e d  to be zero or positive.
c
c******  o u t p u t  ******
c  bvalue.....the value of the (jderiv)-th derivative of  f  at  x .
c
c******  m e t h o d  ******
c     The nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo-
c  cated with the aid of  interv . The  k  b-coeffs of  f  relevant for
c  this interval are then obtained from  bcoef (or taken to be zero if
c  not explicitly available) and are then differenced  jderiv  times to
c  obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
c  Precisely, with  j = jderiv, we have from x.(12) of the text that
c
c     (d**j)f  =  sum ( bcoef(.,j)*b(.,k-j,t) )
c
c  where
c                   / bcoef(.),                     ,  j .eq. 0
c                   /
c    bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1)
c                   / ----------------------------- ,  j .gt. 0
c                   /    (t(.+k-j) - t(.))/(k-j)
c
c     Then, we use repeatedly the fact that
c
c    sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
c  with
c                 (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1)
c    a(.,x)  =    ---------------------------------------
c                 (x - t(.))      + (t(.+m-1) - x)
c
c  to write  (d**j)f(x)  eventually as a linear combination of b-splines
c  of order  1 , and the coefficient for  b(i,1,t)(x)  must then be the
c  desired number  (d**j)f(x). (see x.(17)-(19) of text).
c
      parameter (kmax = 40)
      integer jderiv,k,n,   i,ilo,imk,j,jc,jcmin,jcmax,jj,kmj,km1
     *                     ,mflag,nmi,jdrvp1      
c      integer kmax
C     real bcoef(n),t(1),x,   aj(20),dl(20),dr(20),fkmj
      real bcoef(n),x,   aj(kmax),dl(kmax),dr(kmax),fkmj
      dimension t(n+k)
c former fortran standard made it impossible to specify the length of  t
c  precisely without the introduction of otherwise superfluous addition-
c  al arguments.
      bvalue = 0.
      if (jderiv .ge. k)                go to 99
c
c  *** Find  i   s.t.   1 .le. i .lt. n+k   and   t(i) .lt. t(i+1)   and
c      t(i) .le. x .lt. t(i+1) . If no such i can be found,  x  lies
c      outside the support of  the spline  f , hence  bvalue = 0.
c      (The asymmetry in this choice of  i  makes  f  rightcontinuous, except
c      at  t(n+k) where it is leftcontinuous.)
      call interv ( t, n+k, x, i, mflag )
      if (mflag .ne. 0)                 go to 99
c  *** if k = 1 (and jderiv = 0), bvalue = bcoef(i).
      km1 = k - 1
      if (km1 .gt. 0)                   go to 1
      bvalue = bcoef(i)
                                        go to 99
c
c  *** store the k b-spline coefficients relevant for the knot interval
c     (t(i),t(i+1)) in aj(1),...,aj(k) and compute dl(j) = x - t(i+1-j),
c     dr(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable
c     from input to zero. set any t.s not obtainable equal to t(1) or
c     to t(n+k) appropriately.
    1 jcmin = 1
      imk = i - k
      if (imk .ge. 0)                   go to 8
      jcmin = 1 - imk
      do 5 j=1,i
    5    dl(j) = x - t(i+1-j)
      do 6 j=i,km1
         aj(k-j) = 0.
    6    dl(j) = dl(i)
                                        go to 10
    8 do 9 j=1,km1
    9    dl(j) = x - t(i+1-j)
c
   10 jcmax = k
      nmi = n - i
      if (nmi .ge. 0)                   go to 18
      jcmax = k + nmi
      do 15 j=1,jcmax
   15    dr(j) = t(i+j) - x
      do 16 j=jcmax,km1
         aj(j+1) = 0.
   16    dr(j) = dr(jcmax)
                                        go to 20
   18 do 19 j=1,km1
   19    dr(j) = t(i+j) - x
c
   20 do 21 jc=jcmin,jcmax
   21    aj(jc) = bcoef(imk + jc)
c
c               *** difference the coefficients  jderiv  times.
      if (jderiv .eq. 0)                go to 30
      do 23 j=1,jderiv
         kmj = k-j
         fkmj = float(kmj)
         ilo = kmj
         do 23 jj=1,kmj
            aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
   23       ilo = ilo - 1
c
c  *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
c     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
   30 if (jderiv .eq. km1)              go to 39
      jdrvp1 = jderiv + 1     
      do 33 j=jdrvp1,km1
         kmj = k-j
         ilo = kmj
         do 33 jj=1,kmj
            aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
   33       ilo = ilo - 1
   39 bvalue = aj(1)
c
   99                                   return
      end
      subroutine interv ( xt, lxt, x, left, mflag )
c  from  * a practical guide to splines *  by C. de Boor    
computes  left = max( i :  xt(i) .lt. xt(lxt) .and.  xt(i) .le. x )  .
c
c******  i n p u t  ******
c  xt.....a real sequence, of length  lxt , assumed to be nondecreasing
c  lxt.....number of terms in the sequence  xt .
c  x.....the point whose location with respect to the sequence  xt  is
c        to be determined.
c
c******  o u t p u t  ******
c  left, mflag.....both integers, whose value is
c
c   1     -1      if               x .lt.  xt(1)
c   i      0      if   xt(i)  .le. x .lt. xt(i+1)
c   i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
c   i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
c
c        In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
c        indicates that  x  lies outside the CLOSED interval
c        xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
c        intervals is due to the decision to make all pp functions cont-
c        inuous from the right, but, by returning  mflag = 0  even if
C        x = xt(lxt), there is the option of having the computed pp function
c        continuous from the left at  xt(lxt) .
c
c******  m e t h o d  ******
c  The program is designed to be efficient in the common situation that
c  it is called repeatedly, with  x  taken from an increasing or decrea-
c  sing sequence. This will happen, e.g., when a pp function is to be
c  graphed. The first guess for  left  is therefore taken to be the val-
c  ue returned at the previous call and stored in the  l o c a l  varia-
c  ble  ilo . A first check ascertains that  ilo .lt. lxt (this is nec-
c  essary since the present call may have nothing to do with the previ-
c  ous call). Then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
c  ilo  and are done after just three comparisons.
c     Otherwise, we repeatedly double the difference  istep = ihi - ilo
c  while also moving  ilo  and  ihi  in the direction of  x , until
c                      xt(ilo) .le. x .lt. xt(ihi) ,
c  after which we use bisection to get, in addition, ilo+1 = ihi .
c  left = ilo  is then returned.
c
      integer left,lxt,mflag,   ihi,ilo,istep,middle
      real x,xt(lxt)
      data ilo /1/
      save ilo  
      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt
c
   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100
c
c              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
      istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)                go to 35
         if (x .ge. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      if (x .lt. xt(1))                 go to 90
                                        go to 50
c              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         if (ihi .ge. lxt)              go to 45
         if (x .lt. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 if (x .ge. xt(lxt))               go to 110
      ihi = lxt
c
c           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      if (x .lt. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50
c**** set output and return.
   90 mflag = -1
      left = 1
                                        return
  100 mflag = 0
      left = ilo
                                        return
  110 mflag = 1
	  if (x .eq. xt(lxt)) mflag = 0
      left = lxt
  111 if (left .eq. 1)                  return
	  left = left - 1
	  if (xt(left) .lt. xt(lxt))        return
										go to 111
      end
