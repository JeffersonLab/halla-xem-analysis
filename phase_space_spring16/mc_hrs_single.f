	program mc_hrs_single

C CHANGES FOR OPTICS TESTING:
C 1. Remove target multiple scattering/energy loss

C+______________________________________________________________________________
!
! Monte-Carlo of HRS spectrometer using uniform illumination.
!   This version uses TRANSPORT right-handed coordinate system.
!
! Author: David Potterveld, March-1993
!
! Modification History:
!
!  11-Aug-1993	(D. Potterveld) Modified to use new transformation scheme in
!		which each transformation begins at the pivot.
!
!  19-AUG-1993  (D. Potterveld) Modified to use COSY INFINITY transformations.
C-______________________________________________________________________________

	implicit none

	include 'hrsr/struct_hrsr.inc'
	include 'hrsl/struct_hrsl.inc'
	include 'spectrometers.inc'
	include 'g_dump_all_events.inc'
	include 'constants.inc'

C HBOOK/NTUPLE common block and parameters.
	integer*4	pawc_size
	parameter	(pawc_size = 80000)
	common		/pawc/ hbdata(pawc_size)
	integer*4	hbdata
	character*8	hut_nt_names(18)/
     >			'hsxfp', 'hsyfp', 'hsxpfp', 'hsypfp',
     >			'hsytari', 'hsdeltai', 'hsyptari', 'hsxptari',
     >			'hsytar', 'hsdelta', 'hsyptar', 'hsxptar','fry',
     >                  'frx', 'ok_spec', 'stopwhen', 'x_stop','y_stop'/
	real*4		hut(18)

C Local declarations.
	integer*4	i,
     >			chanin	/1/,
     >			chanout	/2/,
     >			n_trials,trial,
     >			tmp_int

	logical*4	iss

C Event limits, topdrawer limits, physics quantities
	real*8 gen_lim(8)			!M.C. phase space limits.

	real*8 gen_lim_up(3)
	real*8 gen_lim_down(3)

	real*8 cut_dpp,cut_dth,cut_dph,cut_z	!cuts on reconstructed quantities
	real*8 xoff,yoff,zoff                   !Beam offsets (z target offset)
	real*8 spec_xoff,spec_yoff,spec_zoff    !Spectrometer offsets
	real*8 spec_xpoff, spec_ypoff           !Spectrometer angle offsets
	real*8 cos_ts,sin_ts			!cos and sin of spectrometer angle
	real*8 th_ev,cos_ev,sin_ev		!cos and sin of event angle
	real*8 x,y,z,dpp,dxdz,dydz,t1,t2,t3,t4	!temporaries
	real*8 musc_targ_len			!target length for multiple scattering
	real*8 m2				!particle mass squared.
	real*8 rad_len_cm			!conversion r.l. to cm for target
	real*8 pathlen				!path length through spectrometer.
C DJG Variables used for target calcs
	real*8 t,atmp,btmp,ctmp,z_can
	real*8 side_path,costmp,th_can,s_Al
	real*8 forward_path,s_target

C Miscellaneous
	logical*4 ok_spec			!indicates whether event makes it in MC
c	integer*4 hit_calo                      !flag for hitting the calorimeter

C Initial and reconstructed track quantities.
	real*8 dpp_init,dth_init,dph_init,xtar_init,ytar_init,ztar_init
	real*8 dpp_recon,dth_recon,dph_recon,ytar_recon,ztar_recon
	real*8 x_fp,y_fp,dx_fp,dy_fp		!at focal plane
	real*8 frx,fry,fr1,fr2
	real*8 p_spec,th_spec			!spectrometer setting
	real*8 resmult
	real*8 good_evt

C Control flags (from input file)
	integer*4 p_flag			!particle identification
	logical*4 ms_flag
	logical*4 wcs_flag
	integer*4 col_flag
	logical*4 gen_evts_file_flag
	logical*4 use_left_arm

C Hardwired control flags.
	logical*4 hut_ntuple	/.true./
	logical*4 decay_flag	/.false./

	real*8	dpp_var(2),dth_var(2),dph_var(2),ztg_var(2)
	real*4	stime,etime,zero

	character*132	str_line

C Random Number Generation
	integer*4 rnd_seed_flag
	character*80 random_state_file
	logical restorerndstate

C Function definitions.

	integer*4	last_char
	logical*4	rd_int,rd_real
	real*8          grnd,gauss1

        character*80 rawname, filename
	real*4  secnds

        integer iquest
        common/quest/iquest(100)

	save		!Remember it all!

C ================================ Executable Code =============================

C Initialize
C xiaochao:
C using SIMC unstructured version
C
	rSTOP_trials	= 0
	rSTOP_col_entr  = 0
	rSTOP_col_exit  = 0
	rSTOP_spec_entr = 0
	rSTOP_Q1_in	= 0
	rSTOP_Q1_mid	= 0
	rSTOP_Q1_out	= 0
	rSTOP_Q2_in	= 0
	rSTOP_Q2_mid	= 0
	rSTOP_Q2_out	= 0
	rSTOP_Q3_in	= 0
	rSTOP_Q3_mid	= 0
	rSTOP_Q3_out	= 0
	rSTOP_D1_in	= 0
	rSTOP_D1_out	= 0
	rSTOP_hut	= 0
	rSTOP_dc1	= 0
	rSTOP_dc2	= 0
	rSTOP_s0	= 0
	rSTOP_cer       = 0
	rSTOP_s2	= 0
	rSTOP_ps	= 0
	rSTOP_sh        = 0
	rSTOP_successes	= 0

	lSTOP_trials	= 0
	lSTOP_col_entr  = 0
	lSTOP_col_exit  = 0
	lSTOP_spec_entr = 0
	lSTOP_box_entr  = 0
	lSTOP_Q1_in	= 0
	lSTOP_Q1_mid	= 0
	lSTOP_Q1_out	= 0
	lSTOP_Q2_in	= 0
	lSTOP_Q2_mid	= 0
	lSTOP_Q2_out	= 0
	lSTOP_Q3_in	= 0
	lSTOP_Q3_mid	= 0
	lSTOP_Q3_out	= 0
	lSTOP_D1_in	= 0
	lSTOP_D1_out	= 0
	lSTOP_hut	= 0
	lSTOP_dc1	= 0
	lSTOP_dc2	= 0
	lSTOP_s0	= 0
	lSTOP_cer       = 0
	lSTOP_s2	= 0
	lSTOP_prl1	= 0
	lSTOP_prl2      = 0
	lSTOP_successes	= 0

C Open setup file.

	write(*,*)'Enter input filename (assumed to be in infiles dir)'
	read(*,1968) rawname
 1968	format(a)
	filename = 'infiles/'//rawname(1:last_char(rawname))//'.inp'
	write(6,*) filename,'opened'
	open(unit=chanin,status='old',file=filename)

C Initialize HBOOK/NTUPLE if used.
	if (hut_ntuple) then
	  call hlimit(pawc_size)
	  filename = 'worksim/'//rawname(1:last_char(rawname))//'.rzdat'
!	  call hropen(30,'HUT',filename,'N',1024,i)
	  iquest(10) = 256000
	  iquest(10) = 510000
! see for example
!   http://wwwasd.web.cern.ch/wwwasd/cgi-bin/listpawfaqs.pl/7
! the file size is limited to ~260M no matter how I change iquest !
	  call hropen(30,'HUT',filename,'NQ',4096,i) !CERNLIB
 
	  if (i.ne.0) then
	    write(6,*),'HROPEN error: istat = ',i
	    stop
	  endif
	  call hbookn(1,'HUT NTUPLE',18,'HUT',10000,hut_nt_names)
	endif	   

C Open Output file.
	filename = 'outfiles/'//rawname(1:last_char(rawname))//'.hist'
	open (unit=chanout,status='unknown',file=filename)

C Read in real*8's from setup file

	str_line = '!'

C Strip off header

	do while (str_line(1:1).eq.'!')
	  write(6,*) str_line(1:last_char(str_line))
	  read (chanin,1001) str_line
	enddo

! Read data lines.

! N_TRIALS:
c	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_int(str_line,n_trials)
	if (.not.iss) stop 'ERROR (ntrials) in setup!'

! Spectrometer momentum:
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,p_spec)
	if (.not.iss) stop 'ERROR (Spec momentum) in setup!'

! Spectrometer angle:
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,th_spec)
	if (.not.iss) stop 'ERROR (Spec theta) in setup!'
	th_spec = abs(th_spec) / degrad
	cos_ts = cos(th_spec)
	sin_ts = sin(th_spec)

! M.C. limits (half width's for dp,th,ph, full width's for x,y,z)
	do i=1,3
	  read (chanin,1001) str_line
	  write(6,*) str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	  gen_lim_down(i) = gen_lim(i)
	  read (chanin,1001) str_line
	  write(6,*) str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	  gen_lim_up(i) = gen_lim(i)
	enddo
!	do i=1,3
!	   write(*,*)'gen_lim_up/down = ',gen_lim_up(i),' ',gen_lim_down(i)
!	enddo

	do i = 4,6
	  read (chanin,1001) str_line
	  write(6,*) str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	enddo

! Raster size

	do i=7,8
	   read (chanin,1001) str_line
	   write(6,*) str_line(1:last_char(str_line))
	   iss = rd_real(str_line,gen_lim(i))
	   if (.not.iss) stop 'ERROR (Fast Raster) in setup'
	enddo

! Cuts on reconstructed quantities
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dpp)) stop 'ERROR (CUT_DPP) in setup!'

	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dth)) stop 'ERROR (CUT_DTH) in setup!'

	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dph)) stop 'ERROR (CUT_DPH) in setup!'

	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_z)) stop 'ERROR (CUT_Z) in setup!'

! Read in radiation length of target material in cm
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,rad_len_cm)) stop 'ERROR (RAD_LEN_CM) in setup!'

! Beam and target offsets
	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,xoff)
	if(.not.iss) stop 'ERROR (xoff) in setup!'

	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,yoff)
	if(.not.iss) stop 'ERROR (yoff) in setup!'

	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,zoff)
	if(.not.iss) stop 'ERROR (zoff) in setup!'

! Spectrometer offsets
	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_xoff)
	if(.not.iss) stop 'ERROR (spect. xoff) in setup!'

	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_yoff)
	if(.not.iss) stop 'ERROR (spect. yoff) in setup!'

	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_zoff)
	if(.not.iss) stop 'ERROR (spect. zoff) in setup!'

	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_xpoff)
	if(.not.iss) stop 'ERROR (spect. xpoff) in setup!'

	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_ypoff)
	if(.not.iss) stop 'ERROR (spect. ypoff) in setup!'

! read in flag for particle type.
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,p_flag)) stop 'ERROR: p_flag in setup file!'

! Read in flag for aerogel usage in the HMS.
!	read (chanin,1001) str_line
!	write(6,*) str_line(1:last_char(str_line))
!	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: use_aer in setup file!'
!	if (tmp_int.eq.1) use_aer = .true.

! Read in flag for multiple scattering.
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: ms_flag in setup file!'
	if (tmp_int.eq.1) ms_flag = .true.

! Read in flag for wire chamber smearing.
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: wcs_flag in setup file!'
	if (tmp_int.eq.1) wcs_flag = .true.

! Read in flag for dumping ALL events into the HUT NTUPLE (note that
! the ...recon quantities will be ill defined but the FAIL_ID could be
! used to tell...
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR:dump_all_in_ntuple in setup file!'
	if (tmp_int.eq.1) dump_all_in_ntuple = .true.

! Read in flag that sets the collimator option.
	read (chanin,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,col_flag)) stop 'ERROR:col_flag in setup file!'

! Read in flag that sets random number generation option
	read (chanin,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,rnd_seed_flag)) stop 'ERROR:rnd_seed_flag in setup file!'

! Read in flag for using generated events file.
        read (chanin,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: gen_evts_file_flag in setup file!'
	if (tmp_int.eq.1) gen_evts_file_flag = .true.

! Read in flag for which spectrometer to use.
        read (chanin,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: use_left_arm in setup file!'
	if (tmp_int.eq.1) use_left_arm = .true.

C Set particle masses.
	m2 = me2			!default to electron
	if(p_flag.eq.0) then
	  m2 = me2
	else if(p_flag.eq.1) then
	  m2 = mp2
	else if(p_flag.eq.2) then
	  m2 = md2
	else if(p_flag.eq.3) then
	  m2 = mpi2
	else if(p_flag.eq.4) then
	  m2 = mk2
	endif

C Initialize the random number generator
	random_state_file = 'infiles/rand.dat'

        if(rnd_seed_flag.eq.1) then
           if(restorerndstate(random_state_file)) then
              write(6,'(1x,''Random state restored from '',a)')
     >             random_state_file(1:index(random_state_file,' ')-1)
	   endif
	endif

C Initialize (open) generated events file to ID 15
        if(gen_evts_file_flag) then
           open (unit=15,file='infiles/generated_events.dat')
	    write(6,*) 'Taking input from file infiles/generated_events.dat!!!'
        endif

C Print out which arm we will be using
        if(use_left_arm) then
           write(6,*) 'Will be Running Events through LEFT HRS!'
        else
           write(6,*) 'Will be Running Events through RIGHT HRS!'
        endif

C------------------------------------------------------------------------------C
C                           Top of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

	zero=sngl(0.0)
	stime = secnds(zero)
c	print *,'Enter total number of trials'
c	read *,n_trials
	do trial = 1,n_trials

        if(use_left_arm) then
	   if(mod(trial,5000).eq.0) write(*,*)'event #: ',trial,'       successes: ',lSTOP_successes
	else
	   if(mod(trial,5000).eq.0) write(*,*)'event #: ',trial,'       successes: ',rSTOP_successes
	endif

C Pick starting point within target. Z is picked uniformly, X and Y are
C chosen as truncated gaussians, using numbers picked above.
C Units are cm.

	  x = gauss1(3.0) * gen_lim(4) / 6.0			!beam width
	  y = gauss1(3.0) * gen_lim(5) / 6.0			!beam height
	  z = (grnd() - 0.5) * gen_lim(6)                       !along target

C DJG Assume flat raster
	  fr1 = (grnd() - 0.5) * gen_lim(7)   !raster x
	  fr2 = (grnd() - 0.5) * gen_lim(8)   !raster y

	  fry = -fr2  !+y = up, but fry needs to be positive when pointing down

c ... Actually, it seems that +x is left based on
c ... the transformations to spectrometer coordinates
c ... below. (Barak S. April, 2016)
	  frx = -fr1  !+x = right, but frx is left

	  x = x + fr1
	  y = y + fr2

	  x = x + xoff
	  y = y + yoff
	  z = z + zoff

C Pick scattering angles and DPP from independent, uniform distributions.
C dxdz and dydz in HMS TRANSPORT coordinates.

C April 23, 2016 (E. Cohen): Read generated events from a file
          if(gen_evts_file_flag) then
             read(15,*) z, dpp, dxdz, dydz
             dxdz = dxdz/1000.
             dydz = dydz/1000.   
	  else
	     dpp  = grnd()*(gen_lim_up(1)-gen_lim_down(1)) + gen_lim_down(1)
	     dydz = grnd()*(gen_lim_up(2)-gen_lim_down(2))/1000. + gen_lim_down(2)/1000.
	     dxdz = grnd()*(gen_lim_up(3)-gen_lim_down(3))/1000. + gen_lim_down(3)/1000.
	  endif

C Transform from target to HRS (TRANSPORT) coordinates.
C We do the transformation assuming +x is beam left looking
C downstream. See above comments. Barak Schmookler, Dec 2016
          if(use_left_arm) then
             xs    = -y
             ys    = x * cos_ts - z * sin_ts
             zs    = z * cos_ts + x * sin_ts
          else
             xs    = -y
             ys    = x * cos_ts + z * sin_ts
             zs    = z * cos_ts - x * sin_ts
          endif

C DJG Apply spectrometer offsets
C DJG If the spectrometer if too low (positive x offset) a particle
C DJG at "x=0" will appear in the spectrometer to be a little high
C DJG so we should subtract the offset
	  xs = xs - spec_xoff
	  ys = ys - spec_yoff
	  zs = zs - spec_zoff

C Version for spectrometer on the left-hand side:
!	  xs    = -y
!	  ys    = x * cos_ts - z * sin_ts
!	  zs    = z * cos_ts + x * sin_ts

	  dpps  = dpp
	  dxdzs = dxdz
	  dydzs = dydz

C DJG Apply spectrometer angle offsets
	  dxdzs = dxdzs - spec_xpoff/1000.0
	  dydzs = dydzs - spec_ypoff/1000.0


C Save init values for later.
	  xtar_init = xs
	  ytar_init = ys
	  ztar_init = z
	  dpp_init = dpp
	  dth_init = dydzs*1000.		!mr
	  dph_init = dxdzs*1000.		!mr

C Drift back to zs = 0, the plane through the target center
	  xs = xs - zs * dxdzs
	  ys = ys - zs * dydzs
	  zs = 0.0

	  cos_ev = (cos_ts+dydzs*sin_ts)/sqrt(1+dydzs**2+dxdzs**2)
	  th_ev = acos(cos_ev)
	  sin_ev = sin(th_ev)

C Calculate multiple scattering length of target
C-----------------------------------------------------------------------------C
C Version for Hall A (modified by Barak Schmookler 5/23/16)

c ... For GMP Cryo Target:
c ... width: 2.630", wall thickness: 0.203 mm Al
c ... the target has a tip of radius 1.315" and thickness 0.147 mm Al
c ... Although any target length can be set in input file, it's really 15 cm 

c ... So make Beer can with same dimensions:
c ... width: 2.630", wall thickness: 0.203 mm Al
c ... entrance/exit thickness 0.147 mm Al 

! ... compute distances travelled for 12GeV LHRS (RHRS)
! ... 16.0 mil of Al-foil for the target chamber exit foil
! ... 10.62" (13.95") of air between scattering chamber and HRS vacuum
! ... 12.0 mil Kapton for spectrometer entrance (Use mylar, since 
! .....                 X0=28.6cm for Kapton, X0=28.7cm for Mylar) 

!Require Liquid target to be bigger than...2cm
	  if (abs(gen_lim(6)).gt.2.) then   
	     forward_path = (gen_lim(6)/2.-z) / abs(cos_ev)
	     s_target = forward_path
	     
	     side_path = (1.315*2.54) / abs(sin_ev)
	     if (forward_path.lt.side_path) then
		s_Al = 0.0147 / abs(cos_ev)
	     else
		s_target = side_path
		s_Al = 0.0203 / abs(sin_ev)
	     endif
	     musc_targ_len = s_target/rad_len_cm + s_Al/8.89
	  else
!Solid Target
	     musc_targ_len = abs(gen_lim(6)/2. + zoff - z)/rad_len_cm/abs(cos_ev)
	  endif

!Scattering before spectrometer vacuum (assume 12 GeV RHRS)
	  musc_targ_len = musc_targ_len + .016*2.54/8.89 +
     >           13.95*2.54/30420 + .012*2.54/28.7

C Begin transporting particle.  
	  if (ms_flag) call musc(m2,p_spec*(1.+dpps/100.),musc_targ_len,dydzs,dxdzs)

	  if(use_left_arm) then     
             call mc_hrsl(p_spec, th_spec, dpps, xs, ys, zs, dxdzs, dydzs,
     >          x_fp, dx_fp, y_fp, dy_fp, m2,
     >          ms_flag, wcs_flag, decay_flag, resmult, fry, ok_spec, pathlen,
     >          col_flag)
          else       
             call mc_hrsr(p_spec, th_spec, dpps, xs, ys, zs, dxdzs, dydzs,
     >          x_fp, dx_fp, y_fp, dy_fp, m2,
     >          ms_flag, wcs_flag, decay_flag, resmult, fry, ok_spec, pathlen,
     >          col_flag)
          endif

! Option for dumping all events is implemented
	  if (ok_spec .or. dump_all_in_ntuple) then !Success, increment arrays
!	  if (ok_spec) then !Success, increment arrays
	    dpp_recon = dpps
	    dth_recon = dydzs*1000.		!mr
	    dph_recon = dxdzs*1000.		!mr
	    ytar_recon = + ys

            !rough calculation of z-target
            if(use_left_arm) then
               ztar_recon = -ytar_recon/sin_ts
            else
               ztar_recon = ytar_recon/sin_ts
            endif

	    good_evt = 0
	    if(ok_spec) then
	       good_evt = 1
	    endif

C Output NTUPLE entry.

	    if (hut_ntuple) then
	      hut(1) = x_fp
	      hut(2) = y_fp
	      hut(3) = dx_fp
	      hut(4) = dy_fp
	      hut(5) = ytar_init
	      hut(6) = dpp_init
	      hut(7) = dth_init/1000.
	      hut(8) = dph_init/1000.
	      hut(9) = ytar_recon
	      hut(10)= dpp_recon
	      hut(11)= dth_recon/1000.
	      hut(12)= dph_recon/1000.
	      hut(13) = fry
	      hut(14) = frx
	      hut(15) = good_evt
	      hut(16) = stop_where
	      hut(17) = x_stop
	      hut(18) = y_stop
!	      hut(13)= hit_calo 
	      call hfn(1,hut)
	    endif

C Cut on reconstructed quantities.
	    if(ok_spec) then
	       if ((abs(dpp_recon).gt.cut_dpp) .or.
     >             (abs(dth_recon).gt.cut_dth) .or.
     >		   (abs(dph_recon).gt.cut_dph) .or.
     >		   (abs(ztar_recon).gt.cut_z)) then
		  goto 500	!quit if failed
	       endif

C Compute sums for calculating reconstruction variances.
	       dpp_var(1) = dpp_var(1) + (dpp_recon - dpp_init)
	       dth_var(1) = dth_var(1) + (dth_recon - dth_init)
	       dph_var(1) = dph_var(1) + (dph_recon - dph_init)
	       ztg_var(1) = ztg_var(1) + (ztar_recon - ztar_init)

	       dpp_var(2) = dpp_var(2) + (dpp_recon - dpp_init)**2
	       dth_var(2) = dth_var(2) + (dth_recon - dth_init)**2
	       dph_var(2) = dph_var(2) + (dph_recon - dph_init)**2
	       ztg_var(2) = ztg_var(2) + (ztar_recon - ztar_init)**2
	    endif		!Incremented the arrays
	 
	 endif
	   
C We are done with this event, whether GOOD or BAD.
C Loop for remainder of trials.

500	  continue
	enddo				!End of M.C. loop

C------------------------------------------------------------------------------C
C                           End of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

C Save the random number state, but only if good events were generated.
	if(rnd_seed_flag.eq.1) then
	   call saverndstate(random_state_file)
           write(6,'(1x,''Random state saved to '',a)')
     >             random_state_file(1:index(random_state_file,' ')-1)
	endif

C close generated events file....
        if(gen_evts_file_flag) then
           close(15)
        endif

C Elapsed time
	etime = secnds(stime)
	write(6,*) 'Elapsed time = ',etime,' seconds'
	write(6,*) ' '

C Close NTUPLE file.

	call hrout(1,i,' ')
	call hrend('HUT')

	write (chanout,1002)
	write (chanout,1003) p_spec,th_spec*degrad
        write (chanout,1004) (gen_lim(i),i=1,6)

	write (chanout,1005) n_trials

C Indicate where particles are lost in spectrometer.
	if(use_left_arm) then
           write (chanout,1016)
     >  lSTOP_col_entr,lSTOP_col_exit,
     >  lSTOP_spec_entr,lSTOP_box_entr,
     >  lSTOP_Q1_in,lSTOP_Q1_mid,lSTOP_Q1_out,
     >  lSTOP_Q2_in,lSTOP_Q2_mid,lSTOP_Q2_out,
     >  lSTOP_Q3_in,lSTOP_Q3_mid,lSTOP_Q3_out,
     >  lSTOP_D1_in,lSTOP_D1_out

           write (chanout,1007)
     >  lSTOP_trials,lSTOP_hut,lSTOP_dc1,lSTOP_dc2,lSTOP_s0,
     >  lSTOP_cer,lSTOP_s2,lSTOP_prl1,lSTOP_prl2,
     >  lSTOP_successes,lSTOP_successes
        
        else
           write (chanout,1015)
     >  rSTOP_col_entr,rSTOP_col_exit,
     >  rSTOP_spec_entr,
     >  rSTOP_Q1_in,rSTOP_Q1_mid,rSTOP_Q1_out,
     >  rSTOP_Q2_in,rSTOP_Q2_mid,rSTOP_Q2_out,
     >  rSTOP_Q3_in,rSTOP_Q3_mid,rSTOP_Q3_out,
     >  rSTOP_D1_in,rSTOP_D1_out

           write (chanout,1006)
     >  rSTOP_trials,rSTOP_hut,rSTOP_dc1,rSTOP_dc2,rSTOP_s0,
     >  rSTOP_cer,rSTOP_s2,rSTOP_ps,rSTOP_sh,
     >  rSTOP_successes,rSTOP_successes
        
        endif

C Compute reconstruction resolutions.
	if(use_left_arm) then
           if (lSTOP_successes.eq.0) lSTOP_successes=1
	   t1 = sqrt(max(0.,dpp_var(2)/lSTOP_successes - (dpp_var(1)/lSTOP_successes)**2))
	   t2 = sqrt(max(0.,dth_var(2)/lSTOP_successes - (dth_var(1)/lSTOP_successes)**2))
	   t3 = sqrt(max(0.,dph_var(2)/lSTOP_successes - (dph_var(1)/lSTOP_successes)**2))
	   t4 = sqrt(max(0.,ztg_var(2)/lSTOP_successes - (ztg_var(1)/lSTOP_successes)**2))

	   write (chanout,1011) dpp_var(1)/lSTOP_successes,t1,dth_var(1)/lSTOP_successes,
     >          t2,dph_var(1)/lSTOP_successes,t3,ztg_var(1)/lSTOP_successes,t4

	   write(6,*) lSTOP_trials,' Trials',lSTOP_successes,' Successes'
	   write (6,1011) dpp_var(1)/lSTOP_successes,t1,dth_var(1)/lSTOP_successes,
     >          t2,dph_var(1)/lSTOP_successes,t3,ztg_var(1)/lSTOP_successes,t4

        else
           if (rSTOP_successes.eq.0) rSTOP_successes=1
	   t1 = sqrt(max(0.,dpp_var(2)/rSTOP_successes - (dpp_var(1)/rSTOP_successes)**2))
	   t2 = sqrt(max(0.,dth_var(2)/rSTOP_successes - (dth_var(1)/rSTOP_successes)**2))
	   t3 = sqrt(max(0.,dph_var(2)/rSTOP_successes - (dph_var(1)/rSTOP_successes)**2))
	   t4 = sqrt(max(0.,ztg_var(2)/rSTOP_successes - (ztg_var(1)/rSTOP_successes)**2))
           
	   write (chanout,1011) dpp_var(1)/rSTOP_successes,t1,dth_var(1)/rSTOP_successes,
     >          t2,dph_var(1)/rSTOP_successes,t3,ztg_var(1)/rSTOP_successes,t4

	   write(6,*) rSTOP_trials,' Trials',rSTOP_successes,' Successes'
	   write (6,1011) dpp_var(1)/rSTOP_successes,t1,dth_var(1)/rSTOP_successes,
     >          t2,dph_var(1)/rSTOP_successes,t3,ztg_var(1)/rSTOP_successes,t4
        endif

C ALL done!

	stop ' '

C =============================== Format Statements ============================

1001	format(a)
1002	format('!',/,'! Uniform illumination Monte-Carlo results')
1003	format('!',/'! Spectrometer setting:',/,'!',/,
     >  g18.8,' =  P  spect (MeV)',/,
     >  g18.8,' =  TH spect (deg)')

1004	format('!',/'! Monte-Carlo limits:',/,'!',/,
     >  g18.8,' =  GEN_LIM(1) - DP/P                    (half width,% )',/,
     >  g18.8,' =  GEN_LIM(2) - Theta                   (half width,mr)',/,
     >  g18.8,' =  GEN_LIM(3) - Phi                     (half width,mr)',/,
     >  g18.8,' =  GEN_LIM(4) - HORIZ (full width of 3 sigma cutoff,cm)',/,
     >  g18.8,' =  GEN_LIM(5) - VERT  (full width of 3 sigma cutoff,cm)',/,
     >  g18.8,' =  GEN_LIM(6) - Z                       (Full width,cm)')

!inp     >	,/,
!inp     >	g18.8,' =  Hor. 1/2 gap size (cm)',/,
!inp     >	g18.8,' =  Vert. 1/2 gap size (cm)')

1005	format('!',/,'! Summary:',/,'!',/,
     >  i8,' Monte-Carlo trials:')

1006	format(i8,' Initial Trials',/
     >  i8,' Trials made it to the hut',/
     >  i8,' Trial cut in dc1',/
     >  i8,' Trial cut in dc2',/
     >  i8,' Trial cut in scin0',/
     >  i8,' Trial cut in cer',/
     >  i8,' Trial cut in scin2',/
     >  i8,' Trial cut in pre-shower',/
     >  i8,' Trial cut in shower',/
     >  i8,' Trials made it thru the detectors and were reconstructed',/
     >  i8,' Trials passed all cuts and were histogrammed.',/
     >  )

1007	format(i8,' Initial Trials',/
     >  i8,' Trials made it to the hut',/
     >  i8,' Trial cut in dc1',/
     >  i8,' Trial cut in dc2',/
     >  i8,' Trial cut in scin0',/
     >  i8,' Trial cut in cer',/
     >  i8,' Trial cut in scin2',/
     >  i8,' Trial cut in prl1',/
     >  i8,' Trial cut in prl2',/
     >  i8,' Trials made it thru the detectors and were reconstructed',/
     >  i8,' Trials passed all cuts and were histogrammed.',/
     >  )


1008	format(i8)
1009	format(1x,i4,g18.8,i8)
1010	format(a,i8)
1011	format(
     >  'DPP ave error, resolution = ',2g18.8,' %',/,
     >  'DTH ave error, resolution = ',2g18.8,' mr',/,
     >  'DPH ave error, resolution = ',2g18.8,' mr',/,
     >  'ZTG ave error, resolution = ',2g18.8,' cm')

1012	format(1x,16i4)


1015	format(/,
     >  i8,' stopped in COL (Sieve) Entrance',/
     >  i8,' stopped in COL (Sieve) Exit',/
     >  i8,' stopped in Spectrometer Entrance',/
     >  i8,' stopped in Q1 ENTRANCE',/
     >  i8,' stopped in Q1 MIDPLANE',/
     >  i8,' stopped in Q1 EXIT',/
     >  i8,' stopped in Q2 ENTRANCE',/
     >  i8,' stopped in Q2 MIDPLANE',/
     >  i8,' stopped in Q2 EXIT',/
     >  i8,' stopped in Q3 ENTRANCE',/
     >  i8,' stopped in Q3 MIDPLANE',/
     >  i8,' stopped in Q3 EXIT',/
     >  i8,' stopped in D1 ENTRANCE',/
     >  i8,' stopped in D1 EXIT',/
     >  )

1016	format(/,
     >  i8,' stopped in COL (Sieve) Entrance',/
     >  i8,' stopped in COL (Sieve) Exit',/
     >  i8,' stopped in Spectrometer Entrance',/
     >  i8,' stopped in BOX Entrance',/
     >  i8,' stopped in Q1 ENTRANCE',/
     >  i8,' stopped in Q1 MIDPLANE',/
     >  i8,' stopped in Q1 EXIT',/
     >  i8,' stopped in Q2 ENTRANCE',/
     >  i8,' stopped in Q2 MIDPLANE',/
     >  i8,' stopped in Q2 EXIT',/
     >  i8,' stopped in Q3 ENTRANCE',/
     >  i8,' stopped in Q3 MIDPLANE',/
     >  i8,' stopped in Q3 EXIT',/
     >  i8,' stopped in D1 ENTRANCE',/
     >  i8,' stopped in D1 EXIT',/
     >  )

1100	format('!',79('-'),/,'! ',a,/,'!')
1200	format(/,'! ',a,' Coefficients',/,/,
     >  (5(g18.8,','))
     >  )
1300	format(/,'! ',a,' Coefficient uncertainties',/,/,
     >  (5(g18.8,','))
     >  )

	end