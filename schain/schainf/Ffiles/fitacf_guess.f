      subroutine guess(acf,tau,npts,zero,amin,te,tr)
c
c     find zero crossing (zero), depth of minimum (amin), height of maximum
c

      real acf(npts),tau(npts)
c      write(*,*) 'acf: ',acf
c      write(*,*) 'tau: ',tau
c      read(*,*) xx

      zero=0.0
      amin=1.0
      tmin=0.0
      jmin=0

      do i=npts,2,-1
         if(acf(i)*acf(i-1).lt.0.0) then
            zero=(tau(i-1)*acf(i)-tau(i)*acf(i-1))/(acf(i)-acf(i-1))
         end if
         if(acf(i).lt.amin) then
            amin=acf(i)
            jmin=i
         end if
      end do

      if(jmin.gt.0) then
         call parab1(tau(jmin-1),acf(jmin-1),a,b,c)
         tmin=-b/(2.0*a)
         amin=c+tmin*(b+tmin*a)
      end if
      
      tr=cdtr1(-amin)
      te=czte1(zero*1000.0,tr)
c      write(*,*) 'zero: ',zero
c      write(*,*) 'amin: ',amin
c      write(*,*) 'te: ',te
c      write(*,*) 'tr: ',tr
c      read(*,*) xx

      return
      end
      
      subroutine parab1(x,y,a,b,c)
C-----
      dimension x(3),y(3)
      delta=x(1)-x(2)
      a=(y(1)-2.*y(2)+y(3))/(2.*delta*delta)
      b=(y(1)-y(2))/delta - a*(x(1)+x(2))
      c=y(1)-a*x(1)*x(1)-b*x(1)
      return
      end

      real function cdtr1(depth)
C-----convert depth to te/ti ratio
      dimension tr(4)
C     according to the curve published in farley et al 1967
c     modified for 2004 conditions on axis
c      data tr/7.31081,3.53286,5.92271,.174/
      data tr/9.5,4.0,8.5,.3/
      data nt/4/
      cdtr1=tr(1)
      do i=2,nt
        cdtr1=cdtr1*depth + tr(i)
      end do
      return
      end

      real function czte1(zlag,tr)
C-----convert zero crossing point to te
C     according to the curve published in farley et al 1967
c     modified for 2004 conditions on axis
      dimension dt(4)
c      data dt/0.00945025,-0.0774338,.203626,.812397/,nd/4/
      data dt/0.00945025,-0.0774338,.2,0.9/,nd/4/
      data t0/1000./
      tr1=min(abs(tr),5.)
      if(zlag .eq. 0)then
        czte1=1000000.
      else
        dt0=dt(1)
        do  i=2,nd
          dt0=dt0*tr1 + dt(i)
        end do
        czte1=t0*(dt0/zlag)**2
      end if
      return
      end

