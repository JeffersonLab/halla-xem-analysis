! STRUCT_HRSL.INC
!
! Define some record structures
!
! ... where do we loose events?
!	structure /STOP_EVENTS_HRSL/
!		integer*4	trials,slit_hor,slit_vert,
!     >				Q1_in,Q1_mid,Q1_out,
!     >				Q2_in,Q2_mid,Q2_out,
!     >				Q3_in,Q3_mid,Q3_out,
!     >				D1_in,D1_out,
!     >				hut,dc1,dc2,s1,s2,cal,successes
!	end structure
!
! ... a general common block
!	common /spectrometer4/	lSTOP
!
!	record /STOP_EVENTS_HRSL/	lSTOP


! NEW, NON-STRUCTURE VERSION TO MAXIMIZE MACHINE INDEPENDENCE FOR
! STAND ALONE CODES.

	integer*4 lSTOP_trials
	integer*4 lSTOP_col_entr,lSTOP_col_exit
	integer*4 lSTOP_spec_entr
	integer*4 lSTOP_Q1_in,lSTOP_Q1_mid,lSTOP_Q1_out
	integer*4 lSTOP_Q2_in,lSTOP_Q2_mid,lSTOP_Q2_out
	integer*4 lSTOP_Q3_in,lSTOP_Q3_mid,lSTOP_Q3_out
	integer*4 lSTOP_D1_in,lSTOP_D1_out
	integer*4 lSTOP_hut,lSTOP_dc1,lSTOP_dc2
	integer*4 lSTOP_s0,lSTOP_cer,lSTOP_s2
	integer*4 lSTOP_prl1,lSTOP_prl2,lSTOP_successes
	
	common/spectrometer4/ lSTOP_trials,
     >	         lSTOP_col_entr,lSTOP_col_exit,
     >		 lSTOP_spec_entr,
     >		 lSTOP_Q1_in,lSTOP_Q1_mid,lSTOP_Q1_out,
     >		 lSTOP_Q2_in,lSTOP_Q2_mid,lSTOP_Q2_out,
     >		 lSTOP_Q3_in,lSTOP_Q3_mid,lSTOP_Q3_out,
     >		 lSTOP_D1_in,lSTOP_D1_out,
     >		 lSTOP_hut,lSTOP_dc1,lSTOP_dc2,
     >           lSTOP_s0,lSTOP_cer,lSTOP_s2,
     >           lSTOP_prl1,lSTOP_prl2,lSTOP_successes

