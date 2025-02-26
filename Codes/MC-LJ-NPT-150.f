*****************************************************************************
*          Molecular Dynamics code to simulate at NVE collectivity 
*          a system of N atoms of Ar inside a cubic box, interacting
*          through a truncated Lennard-Jones force field at 2.5 sigma. 
*          Leap-frog Verlet integration algorithm.
*
*****************************************************************************
      PROGRAM MC_LJ_NPT_150
      implicit double precision(a-h,o-z)
      double precision mass
      REAL*8 :: nid
    

c     1. Defining dimensions
 
      dimension r(3,1000),vinf(3,1000)
      dimension g(100), gaux(100)

c     2. Reading data and computing related quantities

      open(1,file='leap-lj.data',status='old')
      read(1,*) nconf
      read(1,*) natoms
      read(1,*) mass
      read(1,*) sigma,epsil
      read(1,*) deltat 
      close(1)
      
      pi = 3.14159265359
      nf = 3*natoms-3 !number of degrees of freedom, since we impose the center of masses not to leave the box, hence we have -3 degrees of freedom
      rc = 2.5d0 !sigma = range of the potential in reduced units
      rc1 = rc*sigma !! Not in reduced units, do not use it
      

c     3. Reading initial configuration (positions, velocities) in A and A/ps

      open(2,file='last_config_P_150.dat',status='old') ! We read the last configuration of the last rho (if we want to calculate rho=0.7, we use last conf of rho=0.8)
      do is = 1,natoms
         read(2,*) (r(l,is),l=1,3)
         read(2,*) (vinf(l,is),l=1,3) !OBS: Now we read the velocities for completenness, even though it's not necessary.
      end do         
      read(2,*) boxlength
      close(2)
      
c     Opening files to write results
c
      open(3,file='thermodynamics_P_150.dat',status='unknown')
      open(4,file='gr_P_150.dat',status='unknown')
      open(5,file='last_config_P_150.dat',status='unknown')

      do ii=1,100
        g(ii) = 0
        gaux(ii) = 0
      enddo

c     4. Change to reduced units
      
      
      call reduced(natoms,r,vinf,boxlength,deltat,epsil,sigma,
     &mass,uvel)

c     PSEUDO CODI NPT

      ncycl = 1000000
      nsamp = natoms+1
      ndumped = 0
      ncountervol=0
      delta= 0.15
      beta=1/(2.d0) !2 is the temperature of the isotherm. Since we are working on the (N,P,T) collectivity, T remains constant, hence the MC method will make the system evolve towards thr equilibibriuum temperaature.
      vmax = 0.035
      pressure = 15.d0
      ppot = 0.d0
      delg = boxlength/(2*DFLOAT(100))
      ncounttest = 0
      nconf = ncycl/dfloat(nsamp)
      boxlengthorig=boxlength
      rcg=boxlength/2
      do icycl=1,ncycl !perform ncycl MC cycles
        nran = int(rand()*(natoms+1))+1
        if (nran.le.natoms) then
          call mcmove(r,rc,boxlength,natoms,beta,ndumped,delta) !! displace a particle
        else
          call mcvol(r,boxlength,vmax,beta,boxlengthn,natoms,pressure
     &               ,nchange) !Attempts to change the volume
          if (nchange.eq.1) then
              boxlength=boxlengthn
              ncountervol=ncountervol+1
          endif
          ncounttest = ncounttest +1
        endif
        if (mod(icycl,nsamp).eq.0) then
          call sample(natoms,r,boxlength,boxlengthorig
     &                ,rc,epot,gaux,delg,ppot) !! sample averages
          rho = dfloat(natoms)/(boxlength**3.d0)
          ppot=ppot+rho*2.d0
          do i=1,100  !100 = nhis = total number of bins
            vb = ((i+1)**3-i**3)*delg**3 !volume
            nid = (4.0/3.0)*(4.d0*DATAN(1.d0))*vb*rho !number of ideal gas particles in vb
            gaux(i) = gaux(i) /(dfloat(natoms)*nid)  !normalize g(r) and divide over nconf since we are averaging over all the times that we have "added up" g(r)
            g(i) = g(i) + gaux(i)
            gaux(i)=0
          enddo
          deltappot = 16.d0/3.d0*pi*(rho**2)*
     &            (2.d0/3.d0*(1/rc)**9-(1/rc)**3)  
   
          deltaE = (8.d0*pi/3.d0)*natoms*rho*
     &         (1.d0/3.d0*(1/rc)**9-(1/rc)**3)
          ppot = ppot + deltappot
          epot = epot + deltaE
          write(3,*) icycl/nsamp, ppot, epot, rho, boxlength, g(1)
        endif
      enddo
      nacceptedvol = dfloat(ncountervol)/(dfloat(ncounttest))*100.d0 ! We compute the percentage of accepted vol changes. Important to readjust volmax?
      naccepted = dfloat(ncycl - ndumped)/dfloat(ncycl)*100.d0 ! We compute the percentage of accepted trial moves. Important to readjust delta
      write(3,*) naccepted,nacceptedvol,rho,ncountervol,0,0
      write(3,*) ncounttest,ndumped,boxlength,0,0,0
      do jj = 1, 100
        g(jj)=g(jj)/nconf
      enddo
      do ii = 1, 99 !! to save all the components of g(r)
         write(4,*) g(ii), delg*(ii-0.5) 
      enddo
      
      do is = 1,natoms
        do l = 1,3
          r(l,is)=sigma*r(l,is)
          vinf(l,is)=uvel*vinf(l,is)
        end do   
      end do
      
      do is = 1,natoms
         write(5,*) (r(l,is),l=1,3)
         write(5,*) (vinf(l,is),l=1,3)
      end do   
      boxlength = boxlength*sigma
      write(5,*) boxlength
      
      close(3)
      close(4)
      close(5)
      end
      

*********************************************************
*********************************************************
c              subroutine reduced
*********************************************************
*********************************************************

      subroutine reduced(natoms,r,vinf,boxlength,deltat,
     &epsil,sigma,mass,uvel)
      implicit double precision(a-h,o-z)
      double precision mass
      dimension r(3,1000),vinf(3,1000)

      rgas = 8.314472673d0 !J/(mol*K) 
      utime = sigma*dsqrt(mass/epsil)*dsqrt(10.d0/rgas)
      uvel = sigma/utime !unit of velocity, expressed in A/ps

      boxlength = boxlength/sigma
      deltat = deltat/utime
      do is = 1,natoms
         do l = 1,3
            r(l,is) = r(l,is)/sigma
            vinf(l,is) = vinf(l,is)/uvel 
         end do
      end do

      return
      end


*********************************************************
*********************************************************
c              ener
*********************************************************
*********************************************************

      subroutine ener(r,rc,boxlength,natoms,iio,epot)
      implicit double precision(a-h,o-z)
      dimension r(3,1000)
      dimension rij(3),g(100)
      epot = 0.d0 
      do is = 1,natoms ! since we only need the diference between old and new configurations this means we compute only the enrgy of the affected particle.
        if (is.ne.iio) then
            rr2 = 0.d0
            do l = 1,3
              rijl = r(l,iio) - r(l,is)
              rij(l) = rijl - boxlength*dnint(rijl/boxlength)
              rr2 = rr2 + rij(l)*rij(l)
            end do
            rr = dsqrt(rr2)
            if (rr.lt.rc) then
              ynvrr2 = 1.d0/rr2
              ynvrr6 = ynvrr2*ynvrr2*ynvrr2
              ynvrr12 = ynvrr6*ynvrr6
              forcedist = 24.d0*(2.d0*ynvrr12-ynvrr6)*ynvrr2
              epot = epot + 4.d0*(ynvrr12-ynvrr6)  
            end if
        end if
      end do
      return
      end
      


*********************************************************
*********************************************************
c              mcmove
*********************************************************
*********************************************************
      
      
      SUBROUTINE mcmove(r,rc,boxlength,natoms,beta,ndumped,delta) !! Attempts to displace a particle
      implicit double precision(a-h,o-z)
      dimension ro(3)
      dimension r(3,1000)
      iio=int(rand()*natoms)+1 !! select a particle at random
      
      call ener(r,rc,boxlength,natoms,iio,epoto) !! energy old configuration
      do l=1,3 
        ro(l)=r(l,iio)
        r(l,iio) = ro(l)+ (rand()-0.5)*delta !! give particle random displacement
      !!boundary conditions
         if (r(l,iio).lt.0) r(l,iio) = r(l,iio) + boxlength
         if (r(l,iio).gt.boxlength) r(l,iio) = r(l,iio) - boxlength
      enddo
      call ener(r,rc,boxlength,natoms,iio,epotn) !! energy new configuration
      if(rand().gt.exp(-beta*(epotn-epoto))) then !! inverse acceptance law
        ndumped = ndumped +1 !! Counter of how many configurations are being rejected. We need it to know whether the chosen Delta parameter is the one that makes de simulation robust or not.

        do l=1,3 
          r(l,iio) = ro(l) !! give particle random displacement !! reverse the replacement r(o) by rn
        enddo
      endif
      
      return
      end
      
*********************************************************
*********************************************************
c              subroutine sample
*********************************************************
*********************************************************

      subroutine sample(natoms,r,boxlength,boxlengthorig,
     &                  rc,epot,g,delg,pvirial)
      implicit double precision(a-h,o-z)
      dimension r(3,1000)
      dimension g(100)
      epot = 0.d0 
      pvirial = 0.d0
      !ppot=0.d0
c      atom-atom interactions
      do is = 1,natoms-1
         do js = is+1,natoms
            call ljMC(is,js,r,boxlength,boxlenthorig,rc
     &       ,pot,g,delg,ppot)
            epot = epot + pot
            pvirial= pvirial + ppot
         end do
      end do
      return
      end

*********************************************************
*********************************************************
c        subroutine Lennard-Jones-Montecarlo
*********************************************************
*********************************************************

      subroutine ljMC(is,js,r,boxlength,boxlengthorig,
     &            rc,pot,g,delg,ppot)
      implicit double precision(a-h,o-z)
      dimension r(3,1000)
      dimension rij(3),g(100)

      rr2 = 0.d0
      pot = 0.d0
      ppot = 0.d0
      do l = 1,3
         rijl = r(l,js) - r(l,is)
         rij(l) = rijl - boxlength*dnint(rijl/boxlength)
         rr2 = rr2 + rij(l)*rij(l)
      end do

      rr = dsqrt(rr2)
      if (rr.lt.boxlength/2) then ! only within half the box length
        ig = int(rr/delg)
        g(ig) = g(ig) +2 ! contribution for particle i and j. We add two because we are only looking at each pair once, namely particle 1 and 4. But 4 also has 1 at a ditance r. Since we'll divide by num of particles later on, we add the contribution of all the pairs, at once
      endif

       !! LJ potential computation
      if (rr.lt.rc) then
         ynvrr2 = 1.d0/rr2
         ynvrr6 = ynvrr2*ynvrr2*ynvrr2
         ynvrr12 = ynvrr6*ynvrr6
         forcedist = 24.d0*(2.d0*ynvrr12-ynvrr6)*ynvrr2
         pot = 4.d0*(ynvrr12-ynvrr6)  
         do l = 1,3 
            ppot = ppot + 1.0/(3.d0*boxlength**3)*(rij(l)*rij(l)
     &      *forcedist)
         end do

      end if

      return
      end
      
*********************************************************
*********************************************************
c        subroutine mcvol
*********************************************************
*********************************************************

      SUBROUTINE mcvol(r,boxlength,vmax,beta,boxlengthn,natoms
     & ,p,nchange) !Attempts to change the volume
      implicit double precision(a-h,o-z)
      dimension r(3,1000)
      nchange=1
      call toterg(r,rc,boxlength,natoms,eno) !total energy old configuration
      vo = boxlength**3.d0 !determine old volume
      vnln=log(vo)+(rand()-0.5)*vmax !perform random walk in ln V
      vn = exp(vnln)
      boxlengthn=vn**(1.d0/3.d0) !new bow length
      do i=1,natoms
        do l=1,3
          r(l,i) = r(l,i)*boxlengthn/boxlength !rescale center of mass
        enddo
      enddo
      call toterg(r,rc,boxlengthn,natoms,enn) !total energy new configuration
      arg= -beta*((enn-eno)+p*(vn-vo)
     &     -(natoms+1)*log(vn/vo)/beta) !appropriate weight function!
      if (rand().gt.exp(arg)) then !acceptance rule
        do i=1,natoms !if REJECTED restore the old positions
          do l=1,3
            r(l,i) = r(l,i)*boxlength/boxlengthn !rescale center of mass
          enddo
        enddo
        boxlengthn=boxlength !we go back to the previous boxlength
        nchange=0 !comprobar si la condicion esta bien
      endif
      return
      end
*********************************************************
*********************************************************
c              toterg
*********************************************************
*********************************************************

      subroutine toterg(r,rc,boxlength,natoms,epot)
      implicit double precision(a-h,o-z)
      dimension r(3,1000)
      dimension rij(3)
      epot = 0.d0 
      do is = 1,natoms-1
         do js = is+1,natoms
            rr2 = 0.d0
            do l = 1,3
              rijl = r(l,js) - r(l,is)
              rij(l) = rijl - boxlength*dnint(rijl/boxlength)
              rr2 = rr2 + rij(l)*rij(l)
            enddo
            rr = dsqrt(rr2)
            if (rr.lt.rc) then
              ynvrr2 = 1.d0/rr2
              ynvrr6 = ynvrr2*ynvrr2*ynvrr2
              ynvrr12 = ynvrr6*ynvrr6
              !forcedist = 24.d0*(2.d0*ynvrr12-ynvrr6)*ynvrr2
              epot = epot + 4.d0*(ynvrr12-ynvrr6)  
            endif
          enddo
      enddo
      return
      end
