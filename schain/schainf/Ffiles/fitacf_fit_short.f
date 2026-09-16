
      subroutine fit(wl,taup,rhop,covar,cinv,sigma2p,paramp,ebp,
     &     bfldp,alphap,densp,alt,time,ifitp,ist,nl)
c
c     subroutine to fit measured ACF 
c     wavelength wl (m),lags taup (s), normalized acf rhop,
c     experimental variances sigma2p, covariances covar,
c     fit parameters params, error bars ebp, magnetic field bfldp (gauss),
c     alphap B-field angle (radians), density densp (gcs),
c     altitude alt (km), time (LT hours), nl lags
c     ifitp determines which parameters are fit (see below)
c
      real tol,pi,wl
      parameter(nlmax=100,npmax=10,lwa=2000,tol=1.0e-5)
      external fcn

      real covar(nl,nl),cinv(nl,nl)
      real ev(nlmax*nlmax),ap(nlmax*(nlmax+1)/2)
      real fv1(nlmax),fv2(nlmax),det(2),w(nlmax)

      real taup(nl),rhop(nl),sigma2p(nl),paramp(npmax)
      real p2(npmax),ebp(npmax),alphap,bfldp,densp
      real tau(nlmax),rho(nlmax),sigma2(nlmax),params(npmax)
      real fvec(nlmax),wa(lwa),wa2(lwa),eb(npmax)

      integer iwa(npmax),ifitp(npmax),ifit(npmax)

      include 'fitter.h'
      common /mode/imode

      common/fitter/tau,rho,sigma2,params,ifit
      common/trans/ev      
C-----

c     set ifit(1-5) to unity to fit:
c     1 normalization
c     2 Te 
c     3 Ti
c     4 H+
c     5 He+

      imode=2
c     
      pi=4.0*atan(1.0)
      ak=2.0*pi/wl

      wi(1)=16
      wi(2)=1
      wi(3)=4
      nion=3
      
c
c     invert covariances and find eigenvalues and eigenvectors
c
      l = 0
      do 20 j = 1, nl
         do 10 i = 1, j
            l = l + 1
            ap(l) = covar(i,j)
 10      continue
 20   continue      
      
      call sppfa(ap,nl,info)
      call sppdi(ap,nl,det,01)

      l = 0
      do 40 j = 1, nl
         do 30 i = 1, j
            l = l + 1
            cinv(i,j)=ap(l)
            cinv(j,i)=ap(l)
 30      continue
 40   continue   
      
c
c     transformation matrix is inverse (transpose) of eigenvectors
c
      
      matz=1
      call rs(nl,nl,cinv,w,matz,ev,fv1,fv2,ierr)      
      
      do i=1,nl
         sigma2(i)=1.0/w(i)
      end do

c      
c     extract desired parameters
c
      iparm=0
      do i=1,5
         eb(i)=0.0
         if(ifitp(i).eq.1) then 
            iparm=iparm+1
            p2(iparm)=paramp(i)            
         end if
         ifit(i)=ifitp(i)
         params(i)=paramp(i)
      end do
      np=iparm
      
      alpha=alphap
      dens=densp
      do i=1,nl
        tau(i)=taup(i)
        rho(i)=rhop(i)
        sigma2p(i)=sigma2(i)
      end do

      bfld=bfldp

c
c     no. equations is no. lags - do nlls fit
c
      
      call lmdif1(fcn,nl,np,p2,fvec,tol,info,iwa,wa,lwa)
      ist=info
      
c
c     generate error bars here
c
      call fdjac2(fcn,nl,np,p2,fvec,wa,nl,iflag,0.0e0,wa2)

      do i=1,np
         err=0.0
         do j=1,nl
            err=err+(wa(j+(i-1)*nl))**2
         end do
         if(err.gt.0.0) eb(i)=sqrt(err**-1)
      end do

c     reorder results
      iparm=0
      do i=1,5
         if(ifit(i).eq.1) then 
            iparm=iparm+1
            paramp(i)=p2(iparm)
            ebp(i)=eb(iparm)
         else            
            ebp(i)=0.0
            paramp(i)=params(i)
         end if
      end do

 9    continue
c      write(*,*) dens
c      write(*,*) te
      
      return
      end


      subroutine fcn(m,n,x,fvec,iflag)
c
c     provides m functions (fvec) in n variables (x) to minimize
c     in least-squares sense
c
      
      parameter(nlmax=100,npmax=10)
      integer m,n,iflag
      real fvec(m)
      
      integer ifit(npmax)
      real x(n),fv(nlmax),anorm,params(npmax),ev(nlmax*nlmax)
      real tau(nlmax),rho(nlmax),sigma2(nlmax)
      real chisq,acf

      include 'fitter.h'

      common/fitter/tau,rho,sigma2,params,ifit
      common /errs/ chisq
      common /trans/ev

c     1 0 0 0 0    fit zero lag (otherwise set to default)
c     0 1 0 0 0    fit Te (otherwise default)
c     0 0 1 0 0    fit Ti (otherwise set to Te)
c     0 0 0 1 0    fit H+ (otherwize default)
c     0 0 0 0 1    fit He+ (otherwise default)

 8    continue

c
c     ignore collisions
c
c      write(*,*) "starting fcn"
      ven=0.0
      do i=1,3
         vin(i)=0.0
      end do

      iparm=0

c
c     zero lag normalization constant
c
      if(ifit(1).eq.1) then
         iparm=iparm+1
         c=x(iparm)
      else
         c=params(1)
      end if
c
c     Te
c
      if(ifit(2).eq.1) then
         iparm=iparm+1
         te=x(iparm)
      else
         te=params(2)
      end if
c
c     Ti - default is Te rather than initial value
c
      if(ifit(3).eq.1) then
         iparm=iparm+1
         do i=1,3
            ti(i)=x(iparm)
         end do
      else
         do i=1,3
            ti(i)=te
         end do
         params(3)=te
      end if
  
c     three-ion plasma

      nion=3
c
c     composition H+ first
c
      fi(1)=1.0
      if(ifit(4).eq.1) then
         iparm=iparm+1
         fi(2)=(x(iparm))
         fi(1)=fi(1)-(x(iparm))
      else
         fi(2)=params(4)
         fi(1)=fi(1)-params(4)
      end if
c
c     He+
c
      if(ifit(5).eq.1) then
         iparm=iparm+1
         fi(3)=(x(iparm))
         fi(1)=fi(1)-(x(iparm))
      else
         fi(3)=params(5)
         fi(1)=fi(1)-params(5)
      end if

c      write(*,*) x,bfld,alpha,dens,ak,m
     
      call gaussq(0.0,anorm)
      anorm=anorm/c
 
      if(anorm.eq.0.0.or.m.eq.0)return
      chisq=0.0
      do i=1,m
        call gaussq(tau(i),acf)
        fv(i)=(acf/anorm-rho(i))
      end do

c
c     transform to space where s2 is diagonal
c     (are i and j transposed below?)
c      
      do i=1,m
         fvec(i)=0.0
         do j=1,m
            fvec(i)=fvec(i)+fv(j)*ev(j+(i-1)*m)/sqrt(sigma2(i))
         end do 
         chisq=chisq+fvec(i)**2
      end do

      chisq=chisq/float(m)

c      write(*,*) fvec
c      write(*,*) chisq
c      stop
      
      return
      end

c
c     not really fond of above code
c


      complex function cj_ion(theta,psi)
c     
      real theta,phi,psi,alpha
      complex z
      complex*16 zz
      z=zz(dcmplx(-theta,psi))
      cj_ion=z*cmplx(0.0,-1.0)
      return
      end

      complex function cj_electron(theta,phi,psi,alpha)
c
c     Theta, phi, and psi are the normalized frequency, gyrofrequency,
c     and collision frequency. Alpha is angle between wavevector and
c     magnetic field in radians.
c
      parameter(nterms=10)
      real theta,phi,psi,alpha
      real arg
      complex cj,cy(0:nterms)
      complex z
      complex*16 zz
      integer m,nz,ierr
      integer imode
      common/mode/imode
c
      if(imode.eq.3) then

         arg=0.5*(sin(alpha)/phi)**2

c     calculate modified Bessel functions using amos library
      
         call cbesi(cmplx(arg,0.0),0.0,2,nterms+1,cy,nz,ierr)

         cj=cmplx(0.0,0.0)
     
         do m=-nterms,nterms
            z=zz(dcmplx(-(theta-float(m)*phi)/cos(alpha),
     &           psi/cos(alpha)))
            cj=cj+z*cy(iabs(m))
         end do
         
         cj_electron=cj*cmplx(0.0,-1.0/cos(alpha))

      else

         z=zz(dcmplx(-theta/cos(alpha),psi/cos(alpha)))
         cj_electron=z*cmplx(0.0,-1.0/cos(alpha))
         cj_electron=cj_electron*(1.0-(sin(alpha)/phi)**2/2.0)

      end if

      return
      end

      complex function y_ion(theta,psi)
c
      real theta,phi,psi,alpha
      complex cj, cj_ion
      cj=cj_ion(theta,psi)
      y_ion=cj*cmplx(theta,-psi)+cmplx(0.0,1.0)
      y_ion=y_ion/(1-psi*cj)
      return
      end

      complex function y_electron(theta,phi,psi,alpha)
c
      real theta,phi,psi,alpha
      complex cj, cj_electron
      cj=cj_electron(theta,phi,psi,alpha)
      y_electron=cj*cmplx(theta,-psi)+cmplx(0.0,1.0)
      y_electron=y_electron/(1.0-psi*cj)
      return
      end

      real function spect1(omega)
c
c     Function to generate IS ion-line spectrum 
c     Provisions for nimax ions
c     No provisions for drifts
c
c     imode=1   Farley B-field treatment for electrons 
c           2   use Mike Sulzer's model for ye
c           3   calculate ye using sums of Bessel functions
c
      include 'fitter.h'

      real omega,thetae,thetai(nimax),psie,psii(nimax),phi,p
      real tr(nimax),vti(nimax),vte,omegae
      real bk,em,e,dlf,dl,pi
      real alpha2,densmks,freq
      complex ye,yed,yet,yi(nimax),sum1,sum2,sumdl
      complex y_ion, y_electron, y_esum
      integer i,j,k,imode
      common/mode/ imode
      data bk,em,e,dlf/1.38e-23,9.1e-31,1.6e-19,4772.9/
      
c     fi is ion fraction, wi is ion atomic weight 
c     dens is electron density (cgs), alpha is angle
c     between k and the magnetic field (radians)

      pi=4.0*atan(1.0)

      if(omega.eq.0.0) omega=1.0

      omegae=e*bfld*1.0e-4/em
      densmks=dens*1.0e6

      sum1=0.0
      sum2=0.0
      
      do i=1,nion
        tr(i)=te/ti(i)
        vti(i)=sqrt(bk*ti(i)/(1.67e-27*wi(i)))
        thetai(i)=(omega/ak)/(sqrt(2.0)*vti(i))
        psii(i)=(vin(i)/ak)/(sqrt(2.0)*vti(i))
        yi(i)=y_ion(thetai(i),psii(i))
        sum1=sum1+fi(i)*tr(i)*yi(i)
        sum2=sum2+fi(i)*yi(i)
      end do
      dl=ak**2*dlf*te/densmks
c      write(*,fmt='("Before YE")')
c      write(*,*) imode
c      call exit

      if(imode.eq.1.or.imode.eq.3) then
      
         vte=sqrt(bk*te/em)
         thetae=(omega/ak)/(sqrt(2.0)*vte)
         phi=(omegae/ak)/(sqrt(2.0)*vte)
         psie=(ven/ak)/(sqrt(2.0)*vte)
         ye=y_electron(thetae,phi,psie,alpha)
	 write(*,fmt='("AFTER YE")')
c         call exit
  
      else if(imode.eq.2) then
c
c     use Mike Sulzer's library here: alpha2 is angle off perp in degrees
c
         
         freq=omega/(2.0*pi)
         alpha2=abs(pi/2.0-alpha)*180.0/pi
c         write(*,*) "ye: ", ye
         call collision(densmks, te, freq, alpha2, ye)
c         write(*,fmt='("AFTER COLLISION")')
c         call exit
         ye=ye*omega+cmplx(0.0,1.0)        

      end if
      yed=ye+cmplx(0.0,dl)  
	
      p=(cabs(ye))**2*real(sum2)+cabs(sum1+cmplx(0.0,dl))**2*real(ye)
      p=p/(cabs(yed+sum1))**2
      spect1=p*2.0e0/(omega*pi)
      write(*,*) "spect1:",spect1
      return
      end

      subroutine acf2(wl, tau, te1, ti1, fi1, ven1, vin1, wi1,
     &     alpha1, dens1, bfld1, acf, nion1)
c
c     computes autocorrelation function for given plasma parameters
c     by integrating real spectrum
c     tau in sec., alpha in radians, density in cgs, bfield in cgs
c     scattering wavelength (wl) in meters
c
      include 'fitter.h'

      real wl,tau,te1,ti1(nion1),fi1(nion1),ven1,vin1(nion1),alpha1,
     &     dens1,bfld1,acf
      real pi
      integer nion1
      integer wi1(nion1)
      integer i,j,k,imode
      common /mode/imode
c
      write(*,*) "INITIAL acf:",wl,tau,te1,ti1,fi1,ven1,vin1,wi1,alpha1
      write(*,*) "INITIAL acf:",dens1, bfld1, acf, nion1
c      write(*,fmt='("INIT")')
      pi=4.0*atan(1.0)
c
c     copy arguments to common block
c
      ak=2.0*pi/wl
      imode=2
      write(*,*) "imode:",imode

      nion=nion1
      alpha=alpha1
      te=te1
      ven=ven1
c      write(*,fmt='("INIT2")')
      do i=1,nion
c      write(*,fmt='("INIT2.5")')
        ti(i)=ti1(i)
        fi(i)=fi1(i)
        vin(i)=vin1(i)
        wi(i)=wi1(i)
      end do
c      write(*,fmt='("INIT3")')
      dens=dens1
      bfld=bfld1

c      write(*,*) wl,alpha1,bfld1,dens
c      call exit
c      write(*,fmt='("Before Gauss")')
      call gaussq(tau,acf)

      write(*,*) "FINAL acf:",acf

c      write(*,fmt='("After Gauss")')
      return
      end

      subroutine gaussq(tau,acf)
c
c     Computes cosine or sine transform of given real function
c     Uses quadpack, tau in sec.
c     
      real a,abserr,epsabs,acf,tau,work,chebmo
      integer ier,integr,iwork,last,leniw,lenw,limit,limlst,
     *     lst,maxp1,neval
      integer knorm,ksave,momcom,nrmom
      dimension iwork(300),work(1625),chebmo(61,25)
      external spect1
      include 'fitter.h'
c
c      write(*,*) "inside gaussq"
      a = 0.0
      b = 2.5e4*ak  ! upper integration limit open to debate (was 1.5, now 2.5)
      integr = 1
      epsabs = 1.0e-4
      limlst = 50
      limit = 100
      leniw = limit*2+limlst
      maxp1 = 61
      lenw = leniw*2+maxp1*25

c      write(*,*) leniw,lenw
      
      nrmom=0
      ksave=0
      momcom=0
      write(*,*) "Before qc25f:",acf,imode
c     much faster, more robust
c      write(*,*) "acf_in: ",acf
      call qc25f(spect1,a,b,tau,integr,nrmom,maxp1,ksave,acf,
     &     abserr,neval,resabs,resasc,momcom,chebmo)
c      write(*,*) "acf_out: ",acf
c      call qawf(spect1,a,tau,integr,epsabs,acf,abserr,neval,
c     &     ier,limlst,lst,leniw,maxp1,lenw,iwork,work)
      write(*,*) "After qc25f:",acf
c      call exit
      return
      end
