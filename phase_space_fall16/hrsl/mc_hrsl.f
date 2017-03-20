	subroutine mc_hrsl (p_spec, th_spec, dpp, x, y, z, dxdz, dydz,
     >		x_fp, dx_fp, y_fp, dy_fp, m2,
     >		ms_flag, wcs_flag, decay_flag, resmult, fry, ok_spec, pathlen,
     >          collimator)

C+______________________________________________________________________________
!
! Monte-Carlo of HRSL spectrometer (USED TO BE 'ELECTRON ARM').
!
!
! Author: David Meekins April 2000
! based on code by David Potterveld
!
! Modification History:
!
!  units for this file are percents, cm, and mrads.
!
C-______________________________________________________________________________

	implicit 	none

	include 'struct_hrsl.inc'
	include 'apertures_hrsl.inc'

	include '../g_dump_all_events.inc'
	include '../constants.inc'
	include '../spectrometers.inc'

C Spectrometer definitions - for double arm monte carlo compatability

	integer*4 spectr
	parameter (spectr = 4)	!hrs-left is spec #4.

C Math constants

	real*8 d_r,r_d
	parameter (d_r = pi/180.)
	parameter (r_d = 180./pi)

C The arguments

	real*8	x,y,z				!(cm)
	real*8	dpp				!delta p/p (%)
	real*8 	dxdz,dydz			!X,Y slope in spectrometer
	real*8	x_fp,y_fp,dx_fp,dy_fp		!Focal plane values to return
	real*8	p_spec, th_spec			!spectrometer setting
	real*8  fry				!vertical position@tgt (+y=down)
	real*8  pathlen
	logical	ms_flag				!mult. scattering flag
	logical	wcs_flag			!wire chamber smearing flag
	logical	decay_flag			!check for particle decay
	logical	ok_spec				!true if particle makes it
	integer*4 collimator                    !position of collimator

C Option for "wide open" configuration
        
	logical use_open /.false./
        
C Option for "Large collimator" configuration
        
	logical use_coll /.false./

C Option for sieve slit. Previously, it just took the particles 
! and forced trajectory to put the event in the nearest hole. This would
! mess up the physics distributions somewhat.
! This codes checks to see if events actually go through the sieve holes.
! It checks both the front and back of the holes. Of course, a lot of events
! will get stopped, but this is the real way to do it.
	
	logical use_sieve /.false./ !deprecated, same as external sieve.	
	logical use_ext_sieve /.false./ !use a 12GeV external sieve slit
        logical use_gmp_sieve /.false./ !use a GMp-style sieve slit
	logical sieve_flag /.false./ !set if event passes through sieve

C Collimator (rectangle) dimensions and offsets.

	real*8  h_entr,v_entr	!horiz. and vert. 1/2 gap of fixed slit
	real*8  h_exit,v_exit	!horiz. and vert. 1/2 gap of fixed slit
	real*8  y_off,x_off	!horiz. and vert. position offset of slit
	real*8  z_off		!offset in distance from target to front of sli+

! Offsets...
	parameter (y_off  = 0.0)
	parameter (x_off  = 0.0)
	parameter (z_off  = 0.0)

! z-position of important apertures.
!	real*8 z_entr,z_exit
!	parameter (z_entr = 110.9e0 + z_off)	!nominally 1.109 m
!	parameter (z_exit = z_entr + 8.0e0)	!8.0 cm thick

C Local declarations.

	integer*4	chan/1/,n_classes

	logical*4	first_time_here/.true./

	real*8 dpp_recon,dth_recon,dph_recon	!reconstructed quantities
	real*8 y_recon
	real*8 p,m2				!More kinematic variables.
	real*8 xt,yt,rt,dxt,dyt                 !temporaries
	integer ii,jj		                !counters 
	real*8 resmult				!DC resolution factor (unused)
	real*8 zdrift,ztmp				

	logical dflag			!has particle decayed?
	logical	ok

	real*8 grnd

	save		!Remember it all!

C ================================ Executable Code =============================


! Initialize ok_spec to .flase., restet decay flag

	ok_spec = .false.
	dflag = .false.			!particle has not decayed yet
	lSTOP_trials = lSTOP_trials + 1
	xt = th_spec    !avoid 'unused variable' error for th_spec

! Collimator Option
        if (first_time_here) then
           if(collimator .eq. 0) then
              use_open = .true.
           endif
           
           if(collimator .eq. 1) then
              use_coll = .true.
           endif
           
           if(collimator .eq. 2) then
              use_sieve = .true.
           endif

	   if(collimator .eq. 3) then
              use_ext_sieve = .true.
           endif

           if(collimator .eq. 4) then
              use_gmp_sieve = .true.
           endif

	   ! No collimator - wide open
	   if (use_open .or. use_sieve .or. use_ext_sieve .or. use_gmp_sieve) then
              h_entr = 99.
              v_entr = 99.
              h_exit = 99.
              v_exit = 99.
           endif
           
	   ! Large collimator: (hallaweb.jlab.org/news/minutes/collimator-distance.html)
           if (use_coll) then
              h_entr = 3.145
              v_entr = 6.090
              h_exit = 3.335    !0.1mm narrower than 'hadron arm'.
              v_exit = 6.485
           endif
        endif	

! Save spectrometer coordinates.

	xs = x
	ys = y
	zs = z
	dxdzs = dxdz
	dydzs = dydz

! particle momentum

	dpps = dpp
	p = p_spec*(1.+dpps/100.)

C Read in transport coefficients.

	if (first_time_here) then
	  call transp_init(spectr,n_classes)
	  close (unit=chan)
	  if (n_classes.ne.12) stop 'MC_HRSL, wrong number of transport classes'
	  first_time_here = .false.
	endif

! Begin transporting particle.
! Do transformations, checking against apertures.

!
! Front of external Collimator/Sieve
	zdrift = 107.07
	ztmp = zdrift
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen)
	if ( (abs(ys-y_off).gt.h_entr) .or. (abs(xs-x_off).gt.v_entr) ) then
	   lSTOP_col_entr = lSTOP_col_entr + 1
	   stop_where=1.
	   x_stop=xs
	   y_stop=ys
	   goto 500
	endif

	if (use_sieve .or. use_ext_sieve) then
	   sieve_flag = .false.
	   
	   do ii = 0,20
	      do jj = 0,20
		 xt = 2.50*(jj-10)   !coordinates of hole center
		 yt = 1.25*(ii-10)
		 
		 dxt = xs - xt	     !distances for hole center
		 dyt = ys - yt
		 rt = 0.2	     !distance from center of hole (r=2.0mm)
		 if((abs(xt)<0.01 .and. abs(yt)<0.01) .or. (abs(xt+5.00)<0.01 .and. abs(yt+1.25)<0.01)) then
		    rt = 0.3
		 endif
		 
		 if(sqrt(dxt*dxt+dyt*dyt).lt.rt) then
		    sieve_flag = .true.
		 endif
	      enddo
	   enddo
	   
	   if(.not.sieve_flag) then
	      lSTOP_col_entr = lSTOP_col_entr + 1
	      stop_where=1.
	      x_stop=xs
	      y_stop=ys
	      goto 500
	   endif 
	endif

	if(use_gmp_sieve) then
	   sieve_flag = .false.
	   
           do ii = 0,20
              do jj = 0,20
                 xt = 2.50*(jj-10)     !coordinates of hole center
                 yt = 1.25*(ii-10)
		 
                 dxt = xs - xt	       !distances for hole center                   
		 dyt = ys - yt
		 rt = 0.2	       !distance from center of hole (r=2.0mm)
		 if((abs(xt)<0.01 .and. abs(yt)<0.01) .or. (abs(xt+5.00)<0.01 .and. abs(yt+1.25)<0.01)) then
                    rt = 0.3
                 endif
		 
                 if(sqrt(dxt*dxt+dyt*dyt).lt.rt) then
		    sieve_flag = .true.
                 endif
	      enddo
	   enddo

	   do ii = 0,20
              do jj = 0,20
                 xt = 2.50*(jj-10) + 1.25     !coordinates of hole center
                 yt = 1.25*(ii-10) + 0.625
                 
		 dxt = xs - xt	              !distances for hole center
                 dyt = ys - yt
                 rt = 0.2	              !distance from center of hole (r=2.0mm)
                 if(sqrt(dxt*dxt+dyt*dyt).lt.rt) then
                    sieve_flag = .true.
                 endif
              enddo
	   enddo

           if(.not.sieve_flag) then
              lSTOP_col_entr = lSTOP_col_entr + 1
              stop_where=1.
              x_stop=xs
              y_stop=ys
              goto 500
	   endif
	endif

!
! Back of external Collimator/Sieve
	zdrift = 2.54    !External Sieve is 1" thick
	ztmp = ztmp + zdrift
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen)
	if ( (abs(ys-y_off).gt.h_exit) .or. (abs(xs-x_off).gt.v_exit) ) then
	   lSTOP_col_exit = lSTOP_col_exit + 1
	   stop_where=2.
	   x_stop=xs
	   y_stop=ys
	   goto 500
	endif

	if (use_sieve .or. use_ext_sieve) then
	   sieve_flag = .false.
	   
	   do ii = 0,20
	      do jj = 0,20
		 xt = 2.50*(jj-10)   !coordinates of hole center
		 yt = 1.25*(ii-10)
		 
		 dxt = xs - xt	     !distances for hole center
		 dyt = ys - yt
		 rt = 0.2	     !distance from center of hole (r=2.0mm)
		 if((abs(xt)<0.01 .and. abs(yt)<0.01) .or. (abs(xt+5.00)<0.01 .and. abs(yt+1.25)<0.01)) then
		    rt = 0.3
		 endif
		 
		 if(sqrt(dxt*dxt+dyt*dyt).lt.rt) then
		    sieve_flag = .true.
		 endif
	      enddo
	   enddo
	   
	   if(.not.sieve_flag) then
	      lSTOP_col_exit = lSTOP_col_exit + 1
	      stop_where=2.
	      x_stop=xs
	      y_stop=ys
	      goto 500
	   endif 
	endif

	if(use_gmp_sieve) then
	   sieve_flag = .false.
	   
           do ii = 0,20
              do jj = 0,20
                 xt = 2.50*(jj-10)     !coordinates of hole center
                 yt = 1.25*(ii-10)
		 
                 dxt = xs - xt	       !distances for hole center                   
		 dyt = ys - yt
		 rt = 0.2	       !distance from center of hole (r=2.0mm)
		 if((abs(xt)<0.01 .and. abs(yt)<0.01) .or. (abs(xt+5.00)<0.01 .and. abs(yt+1.25)<0.01)) then
                    rt = 0.3
                 endif
		 
                 if(sqrt(dxt*dxt+dyt*dyt).lt.rt) then
		    sieve_flag = .true.
                 endif
	      enddo
	   enddo

	   do ii = 0,20
              do jj = 0,20
                 xt = 2.50*(jj-10) + 1.25     !coordinates of hole center
                 yt = 1.25*(ii-10) + 0.625
                 
		 dxt = xs - xt	              !distances for hole center
                 dyt = ys - yt
                 rt = 0.2	              !distance from center of hole (r=2.0mm)
                 if(sqrt(dxt*dxt+dyt*dyt).lt.rt) then
                    sieve_flag = .true.
                 endif
              enddo
	   enddo

           if(.not.sieve_flag) then
              lSTOP_col_exit = lSTOP_col_exit + 1
              stop_where=2.
              x_stop=xs
              y_stop=ys
              goto 500
	   endif
	endif


! Circular entrance to spectrometer (where vacuum begins)
	zdrift = 115.73 - ztmp
	ztmp = 115.73
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen)
	if (sqrt(xs*xs+ys*ys).gt.9.92) then
	  lSTOP_spec_entr = lSTOP_spec_entr + 1
	  stop_where=3.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Aperture before Q1 (SOS Quad). This is the long pipe going through magnet.
! (can only check this if next transformation is DRIFT).
        zdrift = 120.64 - ztmp
        ztmp = 120.64
        call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
        if (sqrt(xs*xs+ys*ys).gt.10.185) then
          lSTOP_Q1_in = lSTOP_Q1_in + 1
          stop_where=4.
          x_stop=xs
          y_stop=ys
          goto 500
        endif

! Go to Q1 (SOS Quad) IN  mag bound.

	if (.not.adrift(spectr,1)) write(6,*) 'Transformation #1 is NOT a drift'
	zdrift = driftdist(spectr,1) - ztmp
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay

C Pipe going through Q1 is a little smaller than r_Q1, but use r_Q1 for now
	if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	  lSTOP_Q1_in = lSTOP_Q1_in + 1
	  stop_where=5.	
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Check aperture at 2/3 of Q1 (SOS Quad).

	call transp(spectr,2,decay_flag,dflag,m2,p,46.66666667e0,pathlen)
	if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	  lSTOP_Q1_mid = lSTOP_Q1_mid + 1
	  stop_where=6.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Go to Q1 (SOS Quad) OUT mag boundary.

	call transp(spectr,3,decay_flag,dflag,m2,p,23.33333333e0,pathlen)
	if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	  lSTOP_Q1_out = lSTOP_Q1_out + 1
	  stop_where=7.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Apertures after Q1(SOS Quad), before Q2 (can only check this if next trans. is DRIFT).
	  
	zdrift =  260.92 - 242.05               !Q1 (SOS Quad) exit is z=242.05
        ztmp = zdrift                           !distance from Q1 (SOS Quad) exit
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (sqrt(xs*xs+ys*ys).gt.12.395) then
	  lSTOP_Q1_out = lSTOP_Q1_out + 1
	  stop_where=8.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

	zdrift =  313.77 - 260.92 
        ztmp = ztmp + zdrift                    !distance from Q1 (SOS Quad) exit
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (sqrt(xs*xs+ys*ys).gt.29.53) then
	  lSTOP_Q2_in = lSTOP_Q2_in + 1
	  stop_where=9.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Go to Q2 IN  mag bound.

	if (.not.adrift(spectr,4)) write(6,*) 'Transformation #4 is NOT a drift'
	zdrift = driftdist(spectr,4) - ztmp
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	  lSTOP_Q2_in = lSTOP_Q2_in + 1
	  stop_where=10.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Check aperture at 2/3 of Q2.

	call transp(spectr,5,decay_flag,dflag,m2,p,121.77333333e0,pathlen)
	if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	  lSTOP_Q2_mid = lSTOP_Q2_mid + 1
	  stop_where=11.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Go to Q2 OUT mag boundary.

	call transp(spectr,6,decay_flag,dflag,m2,p,60.88666667e0,pathlen)
	if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	  lSTOP_Q2_out = lSTOP_Q2_out + 1
	  stop_where=12.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Apertures after Q2, before D1 (can only check this if next trans. is DRIFT).

	zdrift = 609.664 - 553.020		!Q2 exit is z=553.02
	ztmp = zdrift				!distance from Q2 exit
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (sqrt(xs*xs+ys*ys).gt.30.0073) then
	  lSTOP_Q2_out = lSTOP_Q2_out + 1
	  stop_where=13.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

	zdrift = 641.800 - 609.664
	ztmp = ztmp + zdrift			!distance from Q2 exit.
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (sqrt(xs*xs+ys*ys).gt.30.0073) then
	  lSTOP_Q2_out = lSTOP_Q2_out + 1
	  stop_where=14.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

	zdrift = 819.489 - 641.800
	ztmp = ztmp + zdrift			!distance from Q2 exit.
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (abs(xs).gt.50.0 .or. abs(ys).gt.15.0) then
	  lSTOP_D1_in = lSTOP_D1_in + 1
	  stop_where=15.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Go to D1 IN magnetic boundary, Find intersection with rotated aperture plane.

	if (.not.adrift(spectr,7)) write(6,*) 'Transformation #7 is NOT a drift'
	zdrift = driftdist(spectr,7) - ztmp
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	xt=xs
	yt=ys
	call rotate_haxis(-30.0,xt,yt)
	if (abs(xt-2.500).gt.52.5) then		! -50 < x < +55
	  lSTOP_D1_in = lSTOP_D1_in + 1	
	  stop_where=16.
	  x_stop=xt
	  y_stop=yt
	  goto 500
	endif
	if ( (abs(yt)+0.01861*xt) .gt. 12.5 ) then !tan(1.066) ~ 0.01861
	  lSTOP_D1_in = lSTOP_D1_in +1
	  stop_where=16.
	  x_stop=xt
	  y_stop=yt
	  goto 500
	endif

! Go to D1 OUT magnetic boundary.
! Find intersection with rotated aperture plane.

	call transp(spectr,8,decay_flag,dflag,m2,p,659.73445725e0,pathlen)
	xt=xs
	yt=ys
	call rotate_haxis(30.0,xt,yt)
	if (abs(xt-2.500).gt.52.5) then		! -50 < x < +55
	  lSTOP_D1_out = lSTOP_D1_out + 1
	  stop_where=17.
	  x_stop=xt
	  y_stop=yt
	  goto 500
	endif

	if ( (abs(yt)+0.01861*xt) .gt. 12.5 ) then !tan(1.066) ~ 0.01861
	  lSTOP_D1_out = lSTOP_D1_out + 1
	  stop_where = 17.
	  x_stop=xt
	  y_stop=yt
	  goto 500
	endif 


! Apertures after D1, before Q3 (can only check this if next trans. is DRIFT).

	zdrift = 1745.33546 - 1655.83446	!D1 exit is z=1655.83446
	ztmp = zdrift				!distance from D1 exit
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (sqrt(xs*xs+ys*ys).gt.30.3276) then
	  lSTOP_D1_out = lSTOP_D1_out + 1
	  stop_where=18.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif
	if (abs(xs).gt.50.0 .or. abs(ys).gt.15.0) then
	  lSTOP_D1_out = lSTOP_D1_out + 1
	  stop_where=18.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

	zdrift = 1759.00946 - 1745.33546
	ztmp = ztmp + zdrift			!distance from D1 exit
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (sqrt(xs*xs+ys*ys).gt.30.3276) then
	  lSTOP_Q3_in = lSTOP_Q3_in + 1
	  stop_where=19.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Go to Q3 IN  mag bound.

	if (.not.adrift(spectr,9)) write(6,*) 'Transformation #9 is NOT a drift'
	zdrift = driftdist(spectr,9) - ztmp
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	  lSTOP_Q3_in = lSTOP_Q3_in + 1	
	  stop_where=20.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Check aperture at 2/3 of Q3.

	call transp(spectr,10,decay_flag,dflag,m2,p,121.7866667e0,pathlen)
	if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	  lSTOP_Q3_mid = lSTOP_Q3_mid + 1
	  stop_where=21.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Go to Q3 OUT mag boundary.

	call transp(spectr,11,decay_flag,dflag,m2,p,60.89333333e0,pathlen)
	if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	  lSTOP_Q3_out = lSTOP_Q3_out + 1
	  stop_where=22.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Apertures after Q3 (can only check this if next trans. is DRIFT).
        zdrift = 2054.03446 - 1997.76446        !Q3 exit is z=1997.76446
        ztmp = zdrift
        call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
        if ((sqrt(xs*xs+ys*ys).gt.35.3219)) then
          lSTOP_Q3_out = lSTOP_Q3_out + 1
          stop_where=23.
          x_stop=xs
          y_stop=ys
          goto 500
        endif

        zdrift = 2080.38746 - 2054.03446 
        ztmp = ztmp + zdrift                    !distance from Q3 exit
        call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
        if (abs(xs).gt.35.56 .or. abs(ys).gt.17.145) then
          lSTOP_Q3_out = lSTOP_Q3_out + 1
          stop_where=24.
          x_stop=xs
          y_stop=ys
          goto 500
        endif

! Vacuum window is 15.522cm before FP (which is at VDC1)
! The window is in the horrizontal plane.

	zdrift = 2327.47246 - 2080.38746
	ztmp = ztmp + zdrift			!distance from Q3 exit
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay

	xt=xs
        yt=ys
        call rotate_haxis(45.0e0,xt,yt)

	if (abs(xt).gt.99.76635 .or. abs(yt).gt.17.145) then
	  lSTOP_Q3_out = lSTOP_Q3_out + 1
	  stop_where=25.
	  x_stop=xs   !Keep as transport
	  y_stop=ys   !Keep as transport
	  goto 500
	endif

! If we get this far, the particle is in the hut.

	lSTOP_hut = lSTOP_hut + 1	

	if (.not.adrift(spectr,12)) write(6,*) 'Transformation #12 is NOT a drift'

	zdrift = driftdist(spectr,12) - ztmp    !distance left to go.
	call mc_hrsl_hut(m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		decay_flag,dflag,resmult,ok,-zdrift,pathlen)
	if (.not.ok) goto 500

! replace xs,ys,... with 'tracked' quantities.  
	xs=x_fp
	ys=y_fp
	dxdzs=dx_fp
	dydzs=dy_fp

! Apply offset to y_fp (detectors offset w.r.t optical axis).
! In the Hall A Analyser, the offset is taken out for recon during conversion
! to rotated (focal plane) coordinates, but NOT for y_fp (transport) in ntuple;
! so we do not apply it to ys (which goes to recon), but do shift it for y_fp.
! This offset is taken as the negative of the y000 term in the reconstruction
! database (vdc). We also apply this offset to the detectors in the hut. But in general 
! the dectectors don't limit the acceptance, so it's not a huge effect.


	y_fp = y_fp - 0.807		!VDC center is +8.07mm for HRSL.

C Reconstruct target quantities.

	call mc_hrsl_recon(dpp_recon,dth_recon,dph_recon,y_recon,fry)
	if (.not.ok) then
	  write(6,*) 'mc_hrsl_recon returned ok=',ok
	  goto 500
	endif

C Fill output to return to main code
	dpp = dpp_recon
	dxdz = dph_recon
	dydz = dth_recon
	y = y_recon
	
	ok_spec = .true.
	lSTOP_successes = lSTOP_successes + 1

C We are done with this event, whether GOOD or BAD.

 500	continue

C ALL done!

	return
	end
