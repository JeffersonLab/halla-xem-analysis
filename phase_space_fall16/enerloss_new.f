	subroutine enerloss_new(len,dens,zeff,aeff,epart,mpart,typeflag,Eloss)

	implicit none

	real*8 thick,len,dens,zeff,aeff,epart,mpart,Eloss
	real*8 x,chsi,lambda,gauss1,Eloss_mp,gamma,tau
	real*8 denscorr,CO,hnup,log10bg,I,beta,Eloss_mp_new
	real*8 Eloss_save,Ke
	real*8 grnd
	real*4 uni,ranlan

	integer typeflag          !1=normal eloss (picked from distribution)
                                  !2=min eloss
	                          !3=max eloss
                                  !4=most probable eloss
        integer numerr
        data numerr /0/

	real*8 me,econ
        parameter(me=0.51099906)
	parameter(econ=0.577216)  !Euler's Constant

	thick = len*dens
	gamma=epart/mpart
        beta = sqrt(1.-1./gamma**2)
	tau = gamma - 1.

C Use Ionization potential from Leo. Barak Schmookler, Nov. 2016	
	if(zeff.lt.1.5) then	!Ionization potential in MeV
	   I = 21.8e-06
	elseif(zeff.lt.13) then
	   I = ((12.*zeff) + 7.)*1.0e-06
	else
	   I = ( (9.76*zeff) + (58.8/(zeff**0.19)) )*1.0e-06
           !I = (16.*zeff**0.9)*1.0e-06
	endif

C Use Sternheimer's parameterization. Also see Hall C engine document.
C Barak Schmookler, Nov. 2016

	hnup = 28.816e-06*sqrt(dens*zeff/aeff) !plasma frequency
	log10bg = log(beta*gamma)/log(10.)
	CO= 2.*log(hnup) - 2.*log(I) - 1.0

	if(log10bg.lt.0.) then
	  denscorr=0.
	elseif(log10bg.lt.3.) then
	  denscorr= CO + 2.*log(10.)*log10bg + abs(CO/27.)*(3.-log10bg)**3
	else
	  denscorr= CO + 2.*log(10.)*log10bg
	endif
	
	if (thick.le.0.) then
	  Eloss = 0.
	else
!	  Eloss_mp = 0.1536e-03 * zeff/aeff * thick * ( 19.26 +
!     &          log(thick/dens) )

	  Eloss_mp_new = 0.1536e-03 * zeff/aeff *thick/beta**2* (
     &	       log(me/I**2) + 1. + log(2.) - econ + 2.*log(gamma*beta) + 
     &	       log(0.1536*zeff/aeff*thick/beta**2)-beta**2-denscorr)


C Average energy loss for electrons. Used in Hall A Analyzer.
C Barak Schmookler, Nov. 2016

!	  Eloss_mp_new = 0.1536e-03 * zeff/aeff *thick/beta**2* (
!     &          log((2+tau)*tau**2/(2*(I/me)**2)) + 1. - beta**2 +
!     &          (tau**2/8. - (2*tau+1)*log(2.))/((tau**2 + 1)**2) -
!     &          denscorr )


c	  write(6,*) 'ELOSS',Eloss_mp,Eloss_mp_new 
! ........ convert to MeV, the unit of choice in THIS program
! ........ (cf. EVCOIN where GeV prevail)
	  Eloss_mp = Eloss_mp_new*1000.
	  chsi = 0.307075/2.*zeff/aeff*thick/beta**2
	  if(typeflag.eq.1)then
	     !x=abs(gauss1(10.0e0))
	     !if(x.gt.0.0) then            
		!lambda = -2.0*log(x)
	     !else
		!lambda = 100000.
	     !endif 
	     
	     uni=grnd()
	     lambda = dble(ranlan(uni))
	     !if(lambda.gt.100.0) lambda = 100.0 !cutoff
	  elseif(typeflag.eq.2)then
	     x=3
	     lambda = -2.0*log(x)
	  elseif(typeflag.eq.3)then
	     x=0.0067
	     lambda = -2.0*log(x)
	  elseif(typeflag.eq.4)then
	     x=1
	     lambda = -2.0*log(x)
          endif
	  
	  Eloss = lambda*chsi+eloss_mp
	endif

        if (eloss.gt.(epart-mpart)) then
	   Eloss_save = eloss
	   Ke = epart-mpart
	   eloss=(epart-mpart)-0.0000001
	   numerr=numerr+1
	   if (numerr.le.10) then
	      if (numerr.eq.1) write(6,*)'--------------------------------'
	      write(6,*) 'Eloss>Total KE; forcing Eloss=KE'
	      write(6,*) 'This occurs in mode',typeflag
	      write(6,*) 'Energy is',epart,'Kinetic Energy is',ke
	      write(6,*) 'Eloss is',eloss_save
	      write(6,*) '--------------------------------'
	      if (numerr.eq.10) write(6,*) '     FURTHER ELOSS ERRORS SUPPRESSED'
	   endif
        endif

	return
	end
