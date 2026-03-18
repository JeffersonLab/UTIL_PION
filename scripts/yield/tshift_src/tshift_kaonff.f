! tshift_kaon.for

c this program uses kinematics for (e,e'K) L/T separation experiment
c to calculate the implied shift in t if the MM is wrong by some amount
c - keeping Q2, W, and meson CM-angle constant, calculate the change in
c   meson kinematics implied by the MMshift
c - meson momentum and meson lab angle changed slightly

c Conventions:
c MMshift = difference between correct and expt values
c   MMexpt=MMcorrect+MMshift (i.e. opposite sign of the MMshift applied to data)
c tEXP=tCORRECT+tSHIFT (i.e. apply the opposite sign of tSHIFT to data)
      
      implicit none

      real*8 mkp,mkpg,mL,mLg,mp,mpg,me,pi
      real*8 nu,nug,pgam,pgamg,pgamx,pgamy,pgamz,thgam
      real*8 tinc,tincg,einc,eincg,escat,escatg,tscat,tscatg
      real*8 pinc,pincx,pincz,pscat,thscat,pscatx,pscatz
      real*8 thkp0,thkpcm,ekp0,ekpcm0,pkp0,pkpcm0,pkpx0,pkpy0,pkpz0
      real*8 thkp1,ekp1,ekpcm1,pkp1,pkpcm1,pkpx1,pkpy1,pkpz1
      real*8 mmshift,mmshgev,tshift,tshiftg
      real*8 pkpg,phix,thkplab0,thkplab1
      real*8 betacm,gammacm,pgamtemp,pkptemp
      real*8 ep1cm,pp1cm,egamcm,pgamcm,p1plus,gamplus,piplus
      real*8 t0,t1,w,wg,q2g,ang

c stuff for linux_suppl.inc
      real*8 sind, cosd, tand
      real*8 asind, acosd, atand, atan2d
      external sind, cosd, tand
      external asind, acosd, atand, atan2d

      me  = 0.511
      mkp = 493.677
      mkpg = mkp/1000.
      mp  = 938.27
      mpg = mp/1000.
      mL  = 1115.683
      mLg = mL/1000.
      
      pi=3.14156

 1000 write(6,100)
 100  format(' Enter Q^2, W in GeV')
      read(5,*)q2g,wg
      if (q2g.eq.0.0) goto 200
      w = wg*1000.
      
      write(6,101)
 101  format(' Enter theta_pq CM (deg)')
      read(5,*)thkpcm
c      thkpcm=0
      
      if (thkpcm.eq.0.) then
         write(6,102)
 102     format(/' Kaon will be emitted in parallel kinematics')
         phix=0.
      elseif (thkpcm.ne.0.) then
         write(6,105)
 105     format(' Enter PHI angle in degrees (0->360)')
         read(5,*)phix
         if (phix.lt.90. .or. phix.ge. 270.) write(6,103)
 103     format(/' Kaon will be emitted forward of q vector')
         if (phix.ge.90. .and. phix.lt. 270.) write(6,104)
 104     format(/' Kaon will be emitted rearward of q vector')
      endif

      write(6,110)
 110  format(' Enter the MM shift (MeV) >0 means data is high')
      read(5,*)mmshift
      mmshgev = mmshift/1000.
         
      nug = (wg**2 + q2g - mpg**2)/(2.*mpg)
      nu  = nug * 1000.
      
      pgamg = sqrt( nug**2 + q2g)
      pgam  = pgamg * 1000.

! find kaon lab momentum
! first, find speed of virtual photon+proton c.m. frame
      betacm  = pgam/(nu+mp)
      gammacm = (nu+mp)/w
         
      ekpcm0 = (w**2 + mkp**2 - (mL**2)) / (2.*w)
      ekpcm1 = (w**2 + mkp**2 - ((mL+mmshift)**2)) / (2.*w)
      pkpcm0 = sqrt( ekpcm0**2 - mkp**2 )
      pkpcm1 = sqrt( ekpcm1**2 - mkp**2 )
c      pkpcmg = pkpcm/1.e3
      
! now transform to lab frame wrt q-vector
      ekp0 = gammacm* (betacm*pkpcm0*cosd(thkpcm) + ekpcm0)
      ekp1 = gammacm* (betacm*pkpcm1*cosd(thkpcm) + ekpcm1)
      pkp0 = sqrt(ekp0**2 - mkp**2)
      pkp1 = sqrt(ekp1**2 - mkp**2)
!      pkpg=pkp/1000.
      thkplab0 = asind(sind(thkpcm)*pkpcm0/pkp0)
      thkplab1 = asind(sind(thkpcm)*pkpcm1/pkp1)

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
      pkpx0 = pkp0*( cosd(thkplab0)*sind(thgam) -
     *     sind(thkplab0)*cosd(thgam)*cosd(phix) )
      pkpy0 = pkp0*sind(thkplab0)*sind(phix)
      pkpz0 = pkp0*( cosd(thkplab0)*cosd(thgam) +
     *     sind(thkplab0)*sind(thgam)*cosd(phix) )
      pkpx1 = pkp1*( cosd(thkplab1)*sind(thgam) -
     *     sind(thkplab1)*cosd(thgam)*cosd(phix) )
      pkpy1 = pkp1*sind(thkplab1)*sind(phix)
      pkpz1 = pkp1*( cosd(thkplab1)*cosd(thgam) +
     *     sind(thkplab1)*sind(thgam)*cosd(phix) )
      
      pkptemp = sqrt(pkpx0**2 + pkpy0**2 + pkpz0**2)
      if (abs(pkptemp-pkp0).gt.0.2) write(6,160)pkp0,pkptemp
 160  format(' Pkp disagreement ',2f10.2)
           
! general equation for t for non-parallel kinematics
! equation is actually for -t
      t0 = (pgamx-pkpx0)**2 +(pgamy-pkpy0)**2 
     *     +(pgamz-pkpz0)**2 -(nu-ekp0)**2
      t1 = (pgamx-pkpx1)**2 +(pgamy-pkpy1)**2 
     *     +(pgamz-pkpz1)**2 -(nu-ekp0)**2
!      tg = t/1.e6
      tshift=t1-t0
      tshiftg=tshift/1.e6
      write(6,165)mmshift,tshiftg
 165  format(' MM shift of ',f5.3,' MeV gives ',f8.5,' GeV^2 shift')
         
! equation for u: between initial proton and outgoing kaon
!      u = mp**2 + mkp**2 -2.*mp*ekp
!      ug= u/1.e6
            
 200  continue
      end

      include 'linux_suppl.inc'
