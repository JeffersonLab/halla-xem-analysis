! CONSTANTS.INC

! Note: All angles are in radians
!	All distances are in cm
!	All energies (momenta,masses) are in MeV
!	All deltas are in percent
!	All densities (thicknesses) are in g/cm3 (g/cm2)
!	All B fields are in kG
!	All roads lead to infinity

! ... some physical constants
	real*8 Me, Me2, Mp, Mp2, Mn, Mn2, Mpi, Mpi2, Mk, Mk2, Md, Md2
	real*8 Mlambda, Msigma0, Msigma_minus
	real*8 amu, hbarc, pi, twopi, alpha, alpi, degrad, euler
	parameter (Me=0.51099906)
	parameter (Me2=Me**2)
	parameter (Mp=938.27231)
	parameter (Mp2=Mp**2)
	parameter (Mn=939.56563)
	parameter (Mn2=Mn**2)
	parameter (Mpi=139.56995)
	parameter (Mpi2=Mpi**2)
        parameter (Mk=493.677)
        parameter (Mk2=Mk**2)
        parameter (Md=1875.613)
        parameter (Md2=Md**2)
	parameter (Mlambda=1115.68)
	parameter (Msigma0=1192.64)
	parameter (Msigma_minus=1197.45)
	parameter (amu=931.49432)
	parameter (hbarc=197.327053)
	parameter (pi=3.141592653589793)
	parameter (twopi=2.*pi)
	parameter (alpha=1./137.0359895)
	parameter (alpi=alpha/pi)
	parameter (degrad=180./pi)
	parameter (euler=0.577215665)

! ... some constants for energy loss
        real*8          rho_Al, Z_Al, A_Al, X0_Al, X0_cm_Al
        parameter       (rho_Al = 2.70)
        parameter       (Z_Al = 13.)
        parameter       (A_Al = 26.98)
        parameter       (X0_Al = 24.01)
        parameter       (X0_cm_Al = X0_Al/rho_Al)
        real*8          rho_mylar, Z_mylar, A_mylar, X0_mylar, X0_cm_mylar
        parameter       (rho_mylar = 1.39)
        parameter       (Z_mylar = 4.545)
        parameter       (A_mylar = 8.735)
        parameter       (X0_mylar = 39.95)
        parameter       (X0_cm_mylar = X0_mylar/rho_mylar)
        real*8          rho_kevlar, Z_kevlar, A_kevlar, X0_kevlar, X0_cm_kevlar
        parameter       (rho_kevlar = 0.74)
        parameter       (Z_kevlar = 2.67)
        parameter       (A_kevlar = 4.67)
        parameter       (X0_kevlar = 55.2)
        parameter       (X0_cm_kevlar = X0_kevlar/rho_kevlar)
        real*8          rho_air, Z_air, A_air, X0_air, X0_cm_air
        parameter       (rho_air = 0.00121)
        parameter       (Z_air = 7.2)   !80% N2, 20% O2
        parameter       (A_air = 14.4)
        parameter       (X0_air = 36.66)
        parameter       (X0_cm_air = X0_air/rho_air)
        real*8          inch_cm
        parameter       (inch_cm = 2.54)

