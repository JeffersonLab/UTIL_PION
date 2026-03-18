! tshift_piff.for

c this program uses kinematics for (e,e'pi) L/T separation experiment
c to calculate the implied shift in t if the MM is wrong by some amount
c - keeping Q2, W, and meson CM-angle constant, calculate the change in
c   meson kinematics implied by the MMshift
c - meson momentum and meson lab angle changed slightly

c Conventions:
c MMshift = difference between correct and expt values
c   MMexpt=MMcorrect+MMshift (i.e. opposite sign of the MMshift applied to data)
c tEXP=tCORRECT+tSHIFT (i.e. apply the opposite sign of tSHIFT to data)

      implicit none

      real*8 mpi,mpig,mn,mng,mp,mpg,me,pi
      real*8 nu,nug,pgam,pgamg,pgamx,pgamy,pgamz,thgam
      real*8 tinc,tincg,einc,eincg,escat,escatg,tscat,tscatg
      real*8 pinc,pincx,pincz,pscat,thscat,pscatx,pscatz
      real*8 thpi0,thpicm,epi0,epicm0,ppi0,ppicm0,ppix0,ppiy0,ppiz0
      real*8 thpi1,epi1,epicm1,ppi1,ppicm1,ppix1,ppiy1,ppiz1
      real*8 mmshift,mmshgev,tshift,tshiftg
      real*8 ppig,phix,thpilab0,thpilab1
      real*8 betacm,gammacm,pgamtemp,ppitemp
      real*8 ep1cm,pp1cm,egamcm,pgamcm,p1plus,gamplus,piplus
      real*8 t0,t1,w,wg,q2g,ang

c stuff for linux_suppl.inc
      real*8 sind, cosd, tand
      real*8 asind, acosd, atand, atan2d
      external sind, cosd, tand
      external asind, acosd, atand, atan2d

      me  = 0.511
      mpi = 139.57
      mpig = mpi/1000.
      mp  = 938.27
      mpg = mp/1000.
      mn  = 939.57
      mng = mn/1000.
      
      pi=3.14156

 1000 write(6,100)
 100  format(' Enter Q^2, W in GeV')
      read(5,*)q2g,wg
      if (q2g.eq.0.0) goto 200
      w = wg*1000.
      
      write(6,101)
 101  format(' Enter theta_pq CM (deg)')
      read(5,*)thpicm
c      thpicm=0
      
      if (thpicm.eq.0.) then
         write(6,102)
 102     format(/' Pion will be emitted in parallel kinematics')
         phix=0.
      elseif (thpicm.ne.0.) then
         write(6,105)
 105     format(' Enter PHI angle in degrees (0->360)')
         read(5,*)phix
         if (phix.lt.90. .or. phix.ge. 270.) write(6,103)
 103     format(/' Pion will be emitted forward of q vector')
         if (phix.ge.90. .and. phix.lt. 270.) write(6,104)
 104     format(/' Pion will be emitted rearward of q vector')
      endif

      write(6,110)
 110  format(' Enter the MM shift (MeV) >0 means data is high')
      read(5,*)mmshift      
      mmshgev = mmshift/1000.
         
      nug = (wg**2 + q2g - mpg**2)/(2.*mpg)
      nu  = nug * 1000.
      
      pgamg = sqrt( nug**2 + q2g)
      pgam  = pgamg * 1000.

! find pion lab momentum
! first, find speed of virtual photon+proton c.m. frame
      betacm  = pgam/(nu+mp)
      gammacm = (nu+mp)/w
         
      epicm0 = (w**2 + mpi**2 - (mn**2)) / (2.*w)
      epicm1 = (w**2 + mpi**2 - ((mn+mmshift)**2)) / (2.*w)
      ppicm0 = sqrt( epicm0**2 - mpi**2 )
      ppicm1 = sqrt( epicm1**2 - mpi**2 )
c      ppicmg = ppicm/1.e3
      
! now transform to lab frame wrt q-vector
      epi0 = gammacm* (betacm*ppicm0*cosd(thpicm) + epicm0)
      epi1 = gammacm* (betacm*ppicm1*cosd(thpicm) + epicm1)
      ppi0 = sqrt(epi0**2 - mpi**2)
      ppi1 = sqrt(epi1**2 - mpi**2)
!      ppig=ppi/1000.
      thpilab0 = asind(sind(thpicm)*ppicm0/ppi0)
      thpilab1 = asind(sind(thpicm)*ppicm1/ppi1)

      write(6,*)' Beam energy (MeV) '
      read(5,*)tinc
      tincg=tinc/1000.
      einc=tinc+me
      eincg=einc/1.e3
      pinc=sqrt(einc**2-me**2)
      pincz=pinc
      pincx=0.
         
      escat=einc-nu
      escatg=escat/1.e3
      if (escat.lt.me) goto 200
      tscat=escat-me
      tscatg=tscat/1000.
      pscat=sqrt(escat**2-me**2)
            
      ang=(pinc**2+pscat**2-pgam**2)/2./pinc/pscat
      if (ang.lt.-1. .or. ang.gt.1.) goto 200
      thscat=acosd(ang)
      pscatz = pscat*cosd(thscat)
      pscatx = pscat*sind(thscat)
            
      pgamz = pincz-pscatz
      pgamx = pincx-pscatx
      pgamy = 0.

      pgamtemp = sqrt(pgamx**2 + pgamy**2 + pgamz**2)
      if (abs(pgamtemp-pgam).gt.0.2) write(6,150)pgam,
     *     pgamtemp
 150  format(' PGam disagreement ',2f10.3)      
      thgam = atand(pgamx/pgamz)
            
! now rotate to lab frame wrt e- beam
      ppix0 = ppi0*( cosd(thpilab0)*sind(thgam) -
     *     sind(thpilab0)*cosd(thgam)*cosd(phix) )
      ppiy0 = ppi0*sind(thpilab0)*sind(phix)
      ppiz0 = ppi0*( cosd(thpilab0)*cosd(thgam) +
     *     sind(thpilab0)*sind(thgam)*cosd(phix) )
      ppix1 = ppi1*( cosd(thpilab1)*sind(thgam) -
     *     sind(thpilab1)*cosd(thgam)*cosd(phix) )
      ppiy1 = ppi1*sind(thpilab1)*sind(phix)
      ppiz1 = ppi1*( cosd(thpilab1)*cosd(thgam) +
     *     sind(thpilab1)*sind(thgam)*cosd(phix) )
            
      ppitemp = sqrt(ppix0**2 + ppiy0**2 + ppiz0**2)
      if (abs(ppitemp-ppi0).gt.0.2) write(6,160)ppi0,ppitemp
 160  format(' Ppi disagreement ',2f10.2)

! general equation for t for non-parallel kinematics
! equation is actually for -t
      t0 = (pgamx-ppix0)**2 +(pgamy-ppiy0)**2 
     *     +(pgamz-ppiz0)**2 -(nu-epi0)**2
      t1 = (pgamx-ppix1)**2 +(pgamy-ppiy1)**2 
     *     +(pgamz-ppiz1)**2 -(nu-epi0)**2
!      tg = t/1.e6
      tshift=t1-t0
      tshiftg=tshift/1.e6
      write(6,165)mmshift,tshiftg
 165  format(' MM shift of ',f5.3,' MeV gives ',f8.5,' GeV^2 shift')
         
! equation for u: between initial proton and outgoing pion
!      u = mp**2 + mpi**2 -2.*mp*epi
!      ug= u/1.e6
            
 200  continue
      end

      include 'linux_suppl.inc'
