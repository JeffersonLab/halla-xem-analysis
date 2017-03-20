	subroutine mc_hrsl_hut (m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		decay_flag,dflag,resmult,ok_hut,zinit,pathlen)

C----------------------------------------------------------------------
C
C Monte-Carlo of HRSL detector hut.
C
C	The particle is stepped through the detector (using project), and
C	multiple scattering is applied for each detector or air gap.
C	If particle decay in enabled, then project.f also checks for
C	decay of particle.  The particle starts at z=zinit.  This
C	needs to be before the first mult. scattering point (the exit window)
C	or the decay distance is negative, and things are BAD.
C
C----------------------------------------------------------------------

	implicit 	none

	include 'struct_hrsl.inc'
	include '../spectrometers.inc'
	include '../g_dump_all_events.inc'

C Math constants

	real*8 pi,d_r,r_d

	parameter (pi = 3.141592654)
	parameter (d_r = pi/180.)
	parameter (r_d = 180./pi)

	real*8 gauss1			!external functions

C all parameters, later to take from .parm files
C----------------------------------------------------------------------
C HRSL_MATERIALS
C CTP parameter file containing the materials of all the HRSL detectors.
C For all materials AFTER the bend only have to do multiple scattering.
C     radlen = 1 radiation length (in cm)
C     thick  = thickness in cm
C In case a "+" is added to the comment, the thickness is a guess only.
C----------------------------------------------------------------------
C spectrometer exit window, 0.1 mm of titanium (X0=3.56cm)
C See JLab-TN-00-001.
	real*8 hfoil_exit_radlen,hfoil_exit_thick
	parameter (hfoil_exit_radlen = 3.56)
	parameter (hfoil_exit_thick  = 0.01)

C spectrometer air gaps
	real*8 hair_radlen
	parameter (hair_radlen = 30420.)


C VDC Information comes from "K.G. Fissum, et al., Nucl. Instr. and Meth. 
C A 474 (2001) 108;"
C According to JLab-TN-00-001, there is ~15cm of air vacuum exit window 
C and center of VDC 1. Here it's set to 15.522cm.

C Aluminum Cage around VDCs. Put at entrance of VDC 1 and exit of
C VDC 2. Cage is 25um thick.
	real*8 hdc_alcage_radlen,hdc_alcage_thick
	parameter (hdc_alcage_radlen = 8.897)
	parameter (hdc_alcage_thick = 0.0025)

C Copper entrance/exit window for each chamber. 125um thick.
	real*8 hdc_cufoil_radlen,hdc_cufoil_thick
	parameter (hdc_cufoil_radlen = 1.436)
	parameter (hdc_cufoil_thick = 0.0125)
	
C Gas windows near entrance/exit of each chamber (put at entrance/exit).
C They are made of 6um thick mylar.
	real*8 hdc_gwin_radlen,hdc_gwin_thick
	parameter (hdc_gwin_radlen = 28.54)
	parameter (hdc_gwin_thick = 0.0006)

C Cathode planes located made of 6um thich mylar. They are located
C at 26.5mm, 52.5mm, and 78.5 mm relative to the chamber entrance.
C (Each chamber is 105 mm thick.) For convenience, we put the first
C and third planes at the chamber entrance and exit, respectively.
	real*8 hdc_cath_radlen,hdc_cath_thick
	parameter (hdc_cath_radlen = 28.54)
	parameter (hdc_cath_thick = 0.0006)

C Gold plating for cathode planes (double sided for c2).
C Thickness is 850A (0.0850um).
	real*8 hdc_gpl_radlen,hdc_gpl_thick
	parameter (hdc_gpl_radlen = 0.9415)
	parameter (hdc_gpl_thick = 8.5E-6)

C Chamber gas, 32/68 ethane/argon
C Each chamber is 105mm thick and we divide it in two for
C the particle transport. The U and V planes are located 39.5mm
C from the chamber entrance and exit, respectively. But for simplicity
C we are keeping everything in the transport coordinate system for tracking. 
C We will determine the tracked x,y position at the center of each
C chamber. Z=0 will be at the center of chamber 1. This is consistent with
C what is done in COSY. In the Hall A Analyzer, however, Z=0 is at the 
C U1 plane (so vdc 1 center is at z = 1.3*sqrt(2) cm in transport coor.). 
C Maybe someone can improve this at some point...
C Also note that the sense wires are 20nm thick tungsten, but
C we assume the particle does not pass through them.
	real*8 hdc_radlen,hdc_thick
	parameter (hdc_radlen = 15815.0)
	!parameter (hdc_thicka  = 3.95)
	!parameter (hdc_thickb  = 1.3)
	parameter (hdc_thick = 5.25)


C There is also one straw chamber (taken from former FPP detector)
C on the spectrometer. The chamber is oriented along the transport
C coordinate system. The chamber has six layers of straw tubes. I put
C the total thickness for all the materials here...

C 100um thick mylar surrounds each tube. The particle will pass through 
C the layer twice. So the total is 0.01cm x 2 x 6 = 0.12cm.
	real*8 hsc_myl_radlen,hsc_myl_thick
	parameter (hsc_myl_radlen = 28.54)
	parameter (hsc_myl_thick = 0.12)

C 10um thick aluminum surrounds each tube. The particle will pass through 
C the layer twice. So the total is 0.001cm x 2 x 6 = 0.012cm.
	real*8 hsc_al_radlen,hsc_al_thick
	parameter (hsc_al_radlen = 8.897)
	parameter (hsc_al_thick = 0.012)

C Chamber gas, 32/68 ethane/argon (same as VDCs)
C Each tube is 10.44mm thick, but the tubes are arranged such that the centers
C of the tubes are only 9.5mm apart. So, the total amount of material is
C 0.95mm x 6 = 5.7cm.
	real*8 hsc_radlen,hsc_thick
	parameter (hsc_radlen = 15815.0)
	parameter (hsc_thick = 5.7)


C hodoscopes (Bicron BC-402 equivalent, I think)
C Base material is  polyvinyltolunene 
	real*8 hscin_radlen
	parameter (hscin_radlen =  42.4)


C Cherenkov entrance foil is 125um (~5mil) thick Al.
	real*8 hcer_entr_radlen,hcer_entr_thick
	parameter (hcer_entr_radlen = 8.897)
	parameter (hcer_entr_thick  = 0.0125)

C Cherenkov, 1 atm of CO2
	real*8 hcer_radlen
	parameter (hcer_radlen = 19650.0)

C Cherenkov mirrors are composite structures, composed of
C ~1mm thick plexiglass + 13mm phenolic honeycomb + 70nm Al +
C 30nm MgF2. In any case, it is supposed to be 5.5e-3 radiation
C lengths thick.
	real*8 hcer_mir_radlen,hcer_mir_thick
	parameter (hcer_mir_radlen = 254.5)
	parameter (hcer_mir_thick  = 1.4)

C Cherenkov exit foil is 75um Tedlar
	real*8 hcer_exit_radlen,hcer_exit_thick
	parameter (hcer_exit_radlen = 25.9)
	parameter (hcer_exit_thick  = 0.0075)


C Both layers of the pion rejector are made of SF-5 type lead glass. The radiation length
C for this material is 2.55cm.
	real*8 hcal_radlen
	parameter (hcal_radlen = 2.55)

C There is 25mm of Al before the PRL1
	real*8 hcal_al_radlen,hcal_al_thick
	parameter (hcal_al_radlen = 8.897)
	parameter (hcal_al_thick  = 2.5)

C---------------------------------------------------------------------------

C Wire chamber resolutions (sigma)

	real*8 hdc_sigma(1:2)/ 0.0225,0.0225 /

C Wire plane positions, construct hdc_zpos array using these parameters

	!integer*4 hdc_nr_cham,hdc_nr_plan
	!parameter (hdc_nr_cham = 1)
	!parameter (hdc_nr_plan = 2)

	real*8 hdc_1_zpos,hdc_1_left,hdc_1_right,hdc_1_top,hdc_1_bot
	real*8 hdc_1x_offset,hdc_1y_offset
	real*8 hdc_2_zpos,hdc_2_left,hdc_2_right,hdc_2_top,hdc_2_bot
	real*8 hdc_2x_offset,hdc_2y_offset
	!real*8 hdc_del_plane

C Drift chamber 1 is the focal plane, so shift all zpos values by 25cm
	parameter (hdc_1_zpos = 0.0)
	parameter (hdc_2_zpos = 47.4) !z-position in transport coor.
	!parameter (hdc_del_plane = hdc_thick)
	parameter (hdc_1_left  =  14.4)
	parameter (hdc_1_right = -14.4)
	parameter (hdc_1y_offset = 0.807)
	parameter (hdc_1_top   = -105.9)
	parameter (hdc_1_bot   =  105.9)
	parameter (hdc_1x_offset = 0.000)
	parameter (hdc_2_left  =  14.4)
	parameter (hdc_2_right = -14.4)
	parameter (hdc_2y_offset = 0.807)
	parameter (hdc_2_top   = -105.9)
	parameter (hdc_2_bot   =  105.9)
	parameter (hdc_2x_offset = 0.000)

C Straw Chamber position
	real *8 hsc_zpos,hsc_left,hsc_right,hsc_top,hsc_bottom
	real *8 hsc_y_offset,hsc_x_offset
	
	parameter (hsc_zpos = 151.0) !entrance to chamber; remember z=0 is vdc 1 center
	parameter (hsc_left = 30.0)
	parameter (hsc_right = -30.0)
	parameter (hsc_y_offset = -2.69)
	parameter (hsc_top = -104.5)
	parameter (hsc_bottom = 104.5)
	parameter (hsc_x_offset = -9.5)

C Scintillator positions and thicknesses

	real*8 hscin1_zpos,hscin1_left,hscin1_right,hscin1_top,hscin1_bottom
	real*8 hscin1_thick,hscin1_y_offset,hscin1_x_offset

	real*8 hscin0_zpos,hscin0_left,hscin0_right,hscin0_top,hscin0_bottom
	real*8 hscin0_thick,hscin0_y_offset,hscin0_x_offset
	
	real*8 hscin2_zpos,hscin2_left,hscin2_right,hscin2_top,hscin2_bottom
	real*8 hscin2_thick,hscin2_y_offset,hscin2_x_offset

	parameter (hscin1_zpos =  135.2)
	parameter (hscin1_thick =  0.5)
	parameter (hscin1_left =  18.0)
	parameter (hscin1_right = -18.0)
	parameter (hscin1_y_offset = 0.0)
	parameter (hscin1_top = -88.0)
	parameter (hscin1_bottom = 88.0)
	parameter (hscin1_x_offset = 0.0)
	
	parameter (hscin0_zpos =  171.3)
	parameter (hscin0_thick =  1.0)
	parameter (hscin0_left =  12.5)
	parameter (hscin0_right = -12.5)
	parameter (hscin0_y_offset = 0.31)
	parameter (hscin0_top = -85.0)
	parameter (hscin0_bottom = 85.0)
	parameter (hscin0_x_offset = -14.0)

	parameter (hscin2_zpos =  316.1)
	parameter (hscin2_thick =  5.0)
	parameter (hscin2_left =  21.6)
	parameter (hscin2_right = -21.6)
	parameter (hscin2_y_offset = 0.0)
	parameter (hscin2_top = -111.8)
	parameter (hscin2_bottom = 111.8)
	parameter (hscin2_x_offset = -12.1)

C Cherenkov position

	real*8 hcer_entr_z,hcer_left,hcer_right,hcer_top,hcer_bottom
	real*8 hcer_thick,hcer_y_offset,hcer_x_offset

	parameter (hcer_entr_z =  179.2)
	parameter (hcer_thick =  128.0)
	parameter (hcer_left =  30.2)
	parameter (hcer_right = -30.2)
	parameter (hcer_y_offset = 0.31)
	parameter (hcer_top = -122.3)
	parameter (hcer_bottom = 122.3)
	parameter (hcer_x_offset = -14.0)

C Calorimeter position

	real*8 hcal_prl1_zpos,hcal_prl1_left,hcal_prl1_right,hcal_prl1_top,hcal_prl1_bottom
	real*8 hcal_prl1_thick,hcal_prl1_y_offset,hcal_prl1_x_offset
	
	real*8 hcal_prl2_left,hcal_prl2_right,hcal_prl2_top,hcal_prl2_bottom
	real*8 hcal_prl2_thick,hcal_prl2_y_offset,hcal_prl2_x_offset

	parameter (hcal_prl1_zpos = 456.2) !may not be exactly mid-point
	parameter (hcal_prl1_thick = 14.75)
	parameter (hcal_prl1_left = 32.5)
	parameter (hcal_prl1_right = -32.5)
	parameter (hcal_prl1_y_offset = 0.0)
	parameter (hcal_prl1_top = -125.4)
	parameter (hcal_prl1_bottom = 125.4)
	parameter (hcal_prl1_x_offset = -20.0)

	!parameter (hcal_prl2_zpos = 475.2) !Not used right now
	parameter (hcal_prl2_thick = 15.0)
	parameter (hcal_prl2_left = 32.5)
	parameter (hcal_prl2_right = -32.5)
	parameter (hcal_prl2_y_offset = 2.3)
	parameter (hcal_prl2_top = -127.5)
	parameter (hcal_prl2_bottom = 127.5)
	parameter (hcal_prl2_x_offset = -11.0)

C The arguments

	real*8	p,m2			!momentum and mass of particle
	real*8	xt,yt			!temporary variables.
	real*8	x_fp,y_fp,dx_fp,dy_fp	!Focal plane values to return
	real*8	xcal,ycal		!Position of track at calorimeter.
	real*8	zinit			!Initial z-position (Not at F.P.)
	real*8	pathlen
	logical ms_flag			!mult. scattering flag.
	logical wcs_flag		!wire chamber smearing flag
	logical decay_flag		!check for decay
	logical ok_hut			!true if particle makes it

C Local declarations.

	integer*4 i

	logical dflag				!has particle decayed?

	real*8	resmult
	real*8	tmpran1,tmpran2			!temporary random numbers
	real*8	radw,drift

	real*8 nsig_max
	parameter(nsig_max=99.0e0)	!max #/sigma for gaussian ran #s.

C These have to be real*4 for the CERNLIB lfit routine

	!real*4	badf				!temporaries
	real*4  xfp4,yfp4,dxfp4,dyfp4           !real*4 versions of fp track
	real*4	xdc(2),ydc(2),zdc(2)		!positions at d.c. planes

C ================================ Executable Code =============================

C Initialize ok_hut to zero

	ok_hut = .false.

C Initialize the xdc and ydc arrays to zero

	do i=1,2
	  xdc(i) = 0.
	  ydc(i) = 0.
	enddo
	resmult = 1.0

C------------------------------------------------------------------------------C
C                           Top of loop through hut                            C
C------------------------------------------------------------------------------C

C Scatter in spectrometer exit foil, which is located at zinit.
C As usual, neglect effect of nonzero dydzs and dxdzs on radw (which means
C the particle goes through the foil at 45 deg).

	radw = (hfoil_exit_thick*sqrt(2.))/hfoil_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C Go to first drift chamber set
C For simplicity, call air MS (probably negligeable) after calling
C 'project' instead of 1/2 way through.

 	drift = hdc_1_zpos - (hdc_thick*sqrt(2.)) - zinit
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

C First Chamber

	!Entrance Foils/Windows (includes Al cage)
	radw = hdc_alcage_thick/hdc_alcage_radlen
	radw = radw + (hdc_cufoil_thick/hdc_cufoil_radlen)
	radw = radw + (hdc_gwin_thick/hdc_gwin_radlen)
	radw = radw + (hdc_cath_thick/hdc_cath_radlen)
	radw = radw + (hdc_gpl_thick/hdc_gpl_radlen)
	radw = radw*sqrt(2.) !for 45 degree particle
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

        !Chamber Gas
	drift = hdc_thick*sqrt(2.) !conversion to transport coor.
	radw = drift/hdc_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	
        !C2 Plane (mutiple scattering only)
	radw = hdc_cath_thick/hdc_cath_radlen
	radw = radw + (2.*hdc_gpl_thick/hdc_gpl_radlen)
	radw = radw*sqrt(2.)	!for 45 degree particle
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	
	!Get Tracking Information
	if(wcs_flag) then
	   tmpran1 = gauss1(nsig_max)
	   tmpran2 = gauss1(nsig_max)
	else
	   tmpran1 = 0.
	   tmpran2 = 0.
	endif
	xdc(1) = xs + hdc_sigma(1)*tmpran1*resmult
	ydc(1) = ys + hdc_sigma(1)*tmpran2*resmult
	zdc(1) = hdc_1_zpos

	!Chamber Gas
	drift = hdc_thick*sqrt(2.) !conversion to transport coor.
	radw = drift/hdc_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	!Exit Foils/Windows
	radw = hdc_cufoil_thick/hdc_cufoil_radlen
	radw = radw + (hdc_gwin_thick/hdc_gwin_radlen)
	radw = radw + (hdc_cath_thick/hdc_cath_radlen)
	radw = radw + (hdc_gpl_thick/hdc_gpl_radlen)
	radw = radw*sqrt(2.) !for 45 degree particle
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	!rotate 45 degrees to compare to VDCs.  CHECK SIGN AND SIZE OF ROTATAION!!!
	xt=xs
	yt=ys
	call rotate_haxis(45.0e0,xt,yt)

	if (xt.gt.(hdc_1_bot+hdc_1x_offset) .or.
     >      xt.lt.(hdc_1_top+hdc_1x_offset) .or.
     >      yt.gt.(hdc_1_left+hdc_1y_offset) .or.
     >      yt.lt.(hdc_1_right+hdc_1y_offset) ) then
	  lSTOP_dc1 = lSTOP_dc1 + 1
	  stop_where=26.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

C First Chamber Done

C Drift to second chamber
C Note that hdc_2_zpos is already in transport coordinate system
	drift = hdc_2_zpos - hdc_1_zpos - ( 2.*hdc_thick*sqrt(2.) )
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

C Second Chamber

	!Entrance Foils/Windows
	radw = hdc_cufoil_thick/hdc_cufoil_radlen
	radw = radw + (hdc_gwin_thick/hdc_gwin_radlen)
	radw = radw + (hdc_cath_thick/hdc_cath_radlen)
	radw = radw + (hdc_gpl_thick/hdc_gpl_radlen)
	radw = radw*sqrt(2.) !for 45 degree particle
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

        !Chamber Gas
	drift = hdc_thick*sqrt(2.) !conversion to transport coor.
	radw = drift/hdc_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	
        !C2 Plane (mutiple scattering only)
	radw = hdc_cath_thick/hdc_cath_radlen
	radw = radw + (2.*hdc_gpl_thick/hdc_gpl_radlen)
	radw = radw*sqrt(2.)	!for 45 degree particle
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	!Get Tracking Information
	if(wcs_flag) then
	   tmpran1 = gauss1(nsig_max)
	   tmpran2 = gauss1(nsig_max)
	else
	   tmpran1 = 0.
	   tmpran2 = 0.
	endif
	xdc(2) = xs + hdc_sigma(2)*tmpran1*resmult
	ydc(2) = ys + hdc_sigma(2)*tmpran2*resmult
	zdc(2) = hdc_2_zpos

	!Chamber Gas
	drift = hdc_thick*sqrt(2.) !conversion to transport coor.
	radw = drift/hdc_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	!Exit Foils/Windows (includes Al cage)
	radw = hdc_alcage_thick/hdc_alcage_radlen
	radw = radw + (hdc_cufoil_thick/hdc_cufoil_radlen)
	radw = radw + (hdc_gwin_thick/hdc_gwin_radlen)
	radw = radw + (hdc_cath_thick/hdc_cath_radlen)
	radw = radw + (hdc_gpl_thick/hdc_gpl_radlen)
	radw = radw*sqrt(2.) !for 45 degree particle
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	
	!rotate 45 degrees to compare to VDCs. CHECK SIGN AND SIZE OF ROTATAION!!!
        xt=xs
        yt=ys
        call rotate_haxis(45.0e0,xt,yt)

	if (xt.gt.(hdc_2_bot+hdc_2x_offset) .or.
     >      xt.lt.(hdc_2_top+hdc_2x_offset) .or.
     >      yt.gt.(hdc_2_left+hdc_2y_offset) .or.
     >      yt.lt.(hdc_2_right+hdc_2y_offset) ) then
	  lSTOP_dc2 = lSTOP_dc2 + 1
	  stop_where=27.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

C Second Chamber Done

C 'fit' track to give new focal plane values

	xfp4 = xdc(1)
	yfp4 = ydc(1)
	dxfp4 = (xdc(2)-xdc(1))/(zdc(2)-zdc(1))
	dyfp4 = (ydc(2)-ydc(1))/(zdc(2)-zdc(1))

	x_fp = dble(xfp4)
	y_fp = dble(yfp4)
	dx_fp = dble(dxfp4)
	dy_fp = dble(dyfp4)

C finished VDCs, drift to unused hodoscope (s1)
	drift = hscin1_zpos - hdc_2_zpos - (hdc_thick*sqrt(2.))
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

C Check if particle will pass through s1 hodoscope.
C Otherwise pass through air.
C Just multiple scatter, no drift
	
	if (xs.lt.(hscin1_bottom+hscin1_x_offset) .or.
     >      xs.gt.(hscin1_top+hscin1_x_offset) .or.
     >      ys.lt.(hscin1_left+hscin1_y_offset) .or.
     >      ys.gt.(hscin1_right+hscin1_y_offset) ) then

	   radw = hscin1_thick/hscin_radlen
	   if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	   
	else
	   !Air
	   radw = hscin1_thick/hair_radlen
	   if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	endif

C finish s1 hodoscope, drift to straw chamber
	drift = hsc_zpos - hscin1_zpos
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

C Check if particle will pass through straw chamber.
C Otherwise pass through air.

	drift = hsc_thick
	
	if (xs.lt.(hsc_bottom+hsc_x_offset) .or.
     >      xs.gt.(hsc_top+hsc_x_offset) .or.
     >      ys.lt.(hsc_left+hsc_y_offset) .or.
     >      ys.gt.(hsc_right+hsc_y_offset) ) then

	   !Chamber Gas
	   radw = drift/hsc_radlen
	   call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	   if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	   !Straws
	   radw = hsc_myl_thick/hsc_myl_radlen
	   radw = radw + (hsc_al_thick/hsc_al_radlen)
	   if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	   
	else
	   !Air
	   radw = drift/hair_radlen
	   call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	   if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	endif

C finished straw chamber, drift to first hodoscope (s0)

	drift = hscin0_zpos - hsc_zpos - hsc_thick
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	if (xs.gt.(hscin0_bottom+hscin0_x_offset) .or.
     >      xs.lt.(hscin0_top+hscin0_x_offset) .or.
     >      ys.gt.(hscin0_left+hscin0_y_offset) .or.
     >      ys.lt.(hscin0_right+hscin0_y_offset) ) then
	   lSTOP_s0 = lSTOP_s0 + 1
	   stop_where=28.
	   x_stop=xs
	   y_stop=ys
	   goto 500
	endif

	radw = hscin0_thick/hscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C finished first hodoscope, drift to cherenkov

	drift = hcer_entr_z - hscin0_zpos
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

C Check if particle is inside Cherenkov. The (x,y) size of
C the cherenkov changes with Z. The minimum is near the 
C entrance, which is where we check.
*       Fix Me!!! Maybe should protect with a pid flag...
*       -------
	if (xs.gt.(hcer_bottom+hcer_x_offset) .or.
     >      xs.lt.(hcer_top+hcer_x_offset) .or.
     >      ys.gt.(hcer_left+hcer_y_offset) .or.
     >      ys.lt.(hcer_right+hcer_y_offset) ) then
	   lSTOP_cer = lSTOP_cer + 1
	   stop_where=29.
	   x_stop=xs
	   y_stop=ys
	   goto 500
	endif

	!Entrance Window
	radw = hcer_entr_thick/hcer_entr_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	!Radiator (CO2 Gas at 1 atm)
	drift = hcer_thick
	radw = drift/hcer_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	!Mirror and Exit Window
	radw = hcer_mir_thick/hcer_mir_radlen
	radw = radw + (hcer_exit_thick/hcer_exit_radlen)
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C drift to second hodoscope (s2m)
C Since s2m is 5cm thick, we'll treat it as an extended
C multiple scatterer. We'll check the location of the particle
C at the back.

	drift = hscin2_zpos - (hscin2_thick/2.) - hcer_entr_z - hcer_thick
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	
	drift = hscin2_thick
	radw = drift/hscin_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	if (xs.gt.(hscin2_bottom+hscin2_x_offset) .or.
     >      xs.lt.(hscin2_top+hscin2_x_offset) .or.
     >      ys.gt.(hscin2_left+hscin2_y_offset) .or.
     >      ys.lt.(hscin2_right+hscin2_y_offset) ) then
	   lSTOP_s2 = lSTOP_s2 + 1
	   stop_where=30.
	   x_stop=xs
	   y_stop=ys
	   goto 500
	endif

C Drift to front of Calorimeter.
	drift = hcal_prl1_zpos - (hcal_prl1_thick/2.) - hcal_al_thick -  hscin2_zpos - (hscin2_thick/2.)
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	! Aluminum before Pre-Shower
	drift = hcal_al_thick
	radw = drift/hcal_al_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	
	!prl1 layer
	drift = hcal_prl1_thick
	radw = drift/hcal_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

C Check particle at back of prl1
C Right now particle must be at least 2.5cm from any edge.
*       Fix Me!!! Maybe should protect with a pid flag...
*       -------
	if (xs.gt.(hcal_prl1_bottom+hcal_prl1_x_offset-2.5) .or.
     >      xs.lt.(hcal_prl1_top+hcal_prl1_x_offset+2.5) .or.
     >      ys.gt.(hcal_prl1_left+hcal_prl1_y_offset-2.5) .or.
     >      ys.lt.(hcal_prl1_right+hcal_prl1_y_offset+2.5) ) then
	   lSTOP_prl1 = lSTOP_prl1 + 1
	   stop_where=31.
	   x_stop=xs
	   y_stop=ys
	   goto 500
	endif

	!prl2 layer
	drift = hcal_prl2_thick
	radw = drift/hcal_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

C Check particle at back of Shower
C Right now particle must be at least 2.5cm from any edge.
*       Fix Me!!! Maybe should protect with a pid flag...
*       -------
	if (xs.gt.(hcal_prl2_bottom+hcal_prl2_x_offset-2.5) .or.
     >      xs.lt.(hcal_prl2_top+hcal_prl2_x_offset+2.5) .or.
     >      ys.gt.(hcal_prl2_left+hcal_prl2_y_offset-2.5) .or.
     >      ys.lt.(hcal_prl2_right+hcal_prl2_y_offset+2.5) ) then
	   lSTOP_prl2 = lSTOP_prl2 + 1
	   stop_where=32.
	   x_stop=xs
	   y_stop=ys
	   goto 500
	endif	
	
	ok_hut = .true.

C We are done with this event, whether GOOD or BAD.

 500   continue

C ALL done!

	return
	end
