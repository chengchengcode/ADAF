

module constant
implicit none
real*8,parameter	::	pi=3.14159265d0
real*8,parameter	::	c=3.d10
real*8,parameter	::	h=6.62606896d-27
real*8,parameter	::	pc=3.08d18
real*8,parameter	::	M_sun=1.989d33
real*8,parameter	::	G=6.672d-8
real*8,parameter	::	m_H=1.673d-24
real*8,parameter	::	m_p=1.673d-24
real*8,parameter	::	sigma_t=6.652d-25
real*8,parameter	::	sigma_SB = 5.67d-5
real*8,parameter	::	e=4.8032d-10
real*8,parameter	::	m_e=9.11d-28
real*8,parameter	::	alpha_f=1./137.

real*8,parameter	::	yr=3.15d7
real*8,parameter	::	kpc=1000.*pc
real*8,parameter	::	km=100000.
real*8,parameter	::	k_bol=1.381d-16
end module constant


module vars_bh
	use constant
	implicit none
	real*8	::	mm_bh=1.d8
	real*8	::	mm_dot=1.d-5
	real*8	::	M_bh,M_dot
	real*8	::	r_g,r_sch
	real*8	::	M_dot_edd
end module vars_bh

module vars_ADAF
	use constant
	use vars_bh
	implicit none
	
	real*8		::		alpha=0.1
	real*8		::		beta=0.5
	real*8		::		gamma=5./3.
	
	real*8		::		mu_i=1.23
	real*8		::		mu_e=1.14

	real*8		::		T_i,T_e
	real*8		::		theta_i 
	real*8		::		theta_e 
	real*8		::		sigma,w_i,w_e,n_i,n_e
!	real*8		::		r
	real*8		::		W,cs,P_ADAF
	

	real*8		::		sonic_jump
	real*8		::		sonic_crit

	real*8		::		l_in
	real*8		::		omega
!	real*8		::		
end module vars_adaf

module vars_clouds
	use constant
	use vars_bh
	implicit none
	real*8		::		sigma_cl,v_r_cl,v_phi_cl
	real*8		::		k_c_cl=0.
	real*8		::		k_c_term=0.
	real*8		::		cloud_1=0.
	real*8		::		cloud_2=0.
	real*8		::		cloud_3=0.
!
	
	real*8		::		v_r_2
	
	real*8		::		mass_ratio = 1.
	real*8		::		mass_ratio_factor = 1.!2.
	real*8		::		M_cl = 4.d23
	real*8		::		R_cl = 0.!5.d10
	real*8		::		R_cl_out = 1.d12
	real*8		::		cloud_size_index = 0.5
	real*8		::		T_cl = 1.d4
	real*8		::		tao_cl,R_cl_max,R_cl_min, n_cl,qq_cl
	
	real*8		::		rr_tidal
	
	real*8		::		gamma_phi= 0.!-1.d-2
	real*8		::		gamma_R  = 0.!1.d-2
	
	real*8		::		det_VR2 = 0.
	real*8		::		det_VPhi2 = 0.
	real*8		::		v_phi_ratio_v_R = 0.
end module vars_clouds

module vars_boundary
    implicit none
    real*8      ::      Omega_boundary_factor = 0.8
    real*8      ::      T_out_factor = 0.1
end module vars_boundary
	

subroutine ADAF_solve()
	use constant
	use vars_bh
	use vars_adaf
	use vars_clouds
    	implicit none 
      	real*8	::		a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,b51,b52,b53,b54,&
					b61,b62,b63,b64,b65,c1,c2,c3,c4,c5,c6,dc1,dc2,dc3,dc4,dc5,dc6
      	parameter(a2=0.2d0,a3=0.3d0,a4=0.6d0,a5=1.0d0&
     		,a6=0.875d0,b21=0.2d0,b31=0.075d0,b32=0.225d0&
     		,b41=0.3d0,b42=-0.9d0, b43=1.2d0&
     		,b51=-11.0d0/54.0d0,b52=2.5d0,b53=-70.0d0/27.0d0&
     		,b54=35.0d0/27.0d0,b61=1631.0d0/55296.0d0&
     		,b62=175.0d0/512.0d0,b63=575.0d0/13824.0d0&
     		,b64=44275.0d0/110592.0d0,b65=253.0d0/4096.0d0,&
			c1=37.0d0/378.0d0,c2=0.0d0,c3=250.0d0/621.0d0,&
			c4=125.0d0/594.0d0,c5=0.0d0,c6=512.0d0/1771.0d0,&
			dc1=c1-2825.0d0/27648.0d0,dc2=0.0d0,&
			dc3=c3-18575.0d0/48384.0d0,dc4=c4-13525.0d0/55296.0d0,&
			dc5=-277.0d0/14336.0d0,dc6=c6-0.25d0)    	
	  	real*8 safety,pgrow,pshrnk,errcon,eacc_test
      	parameter(safety=0.9d0,pgrow=-0.20d0,pshrnk=-0.25d0,errcon=1.89d-4,&
			eacc_test=1.0d-5)
		integer,parameter	::		N_step=50000000
		integer		::		i_test,i_test_0,i_mu,i,j
		real*8		::		k1(3),k2(3),k3(3),k4(3),k5(3),k6(3)	
  		real*8 		::		x_now,wip,wep,sigp,sig(N_step),wi(N_step),we(N_step)
		real*8		::		y1_0,y2_0,y3_0		
		real*8		::		x_test_0,x_test,y_test_0(3),y_test(3),yerr_test(3),df_test(3),&
						k1_test(3),k2_test(3),k3_test(3),k4_test(3),k5_test(3),k6_test(3),&
							h_step_test,df0_test(3)
		real*8		::		x_in,x_out,y_in(3),y(3)
		real*8		::		errmax_test,htemp_test

		real*8		::		drstep_max=5.d-3
		real*8		::		drstep_min=1.d-9
		
		real*8		::		ww_i,ww_e,ssigma,rr
		real*8		::		W_i_gas,W_e_gas,W_gas,sigma_gas,r
		real*8		::		v_r,c_s 
		
		real*8		::		dim_W,dim_sigma,dim_l_in,dim_Q,dim_omega

		real*8		::		SIGMA_GAS_0,SSIGMA_0,W_I_GAS_0,WW_I_0,W_E_GAS_0,WW_E_0,DSIGMA_DR,F_R,F_PHI,V_PHI_2,V_R_2_0
		real*8		::		omega_0,x_last
!		integer		::		stat,error_allocate=0
	
		integer*8,parameter	::		N_radius = 1900
	
		real*8		::		flux(200,N_radius),r_array(N_radius),tau_e_array(N_radius),T_e_array(N_radius)
		real*8		::		sigma_array(N_radius),W_i_array(N_radius),W_e_array(N_radius)
	
		real*8		::		mu(200),f_r_mu(N_radius,100),r_test(N_radius-1)!,f_mu(200)
		
		real*8,external	::		height,omega_k,integ_for_q,YITA_COMPTON,flux_nu,flux_nu_syn,flux_nu_bre
		
		open(unit=21,file='r_sigma_W_i_W_e.txt')
!		
!		sigma_array=999.
!		W_i_array=999.
!		W_e_array=999.
!		
!	
!		write(*,*)r
!		stop
!	


	open(unit=35,file='./spec/nu.txt')
	do i =1,200,1
		read(35,*)mu(i)
	end do
	close(35)

		dim_W=M_dot_edd*c/r_sch
		dim_sigma=M_dot_edd/c/r_sch
		dim_omega=c/r_sch
		dim_l_in=c*r_sch
		dim_Q=M_dot_edd*c*c/r_sch/r_sch

		call INITIALIZATION_func(x_in,x_out,y_in)
! 
	
	open(unit=39,file='./spec/r_array.txt')	
	
	open(unit = 58, file= 'spec_r_mu.txt')

		write(*,*)x_in,'x_out=',exp(x_out),y_in
!		write(*,*)'x_out=',exp(x_in)
!		stop

		x_now=x_in		
		h_step_test=(x_out-x_now)/N_step
		sig(1)=y_in(1)
		wi(1)=y_in(2)
		we(1)=y_in(3)

	ssigma_0 =  sig(1)	
	ww_i_0 	 =  wi(1)
	ww_e_0   =  we(1)
	v_R_2_0  = (M_dot/2./pi/ssigma_0/dim_sigma/exp(x_in)/r_sch)**2.

	sigma_gas = ssigma_0*dim_sigma
	W_i_gas   =	ww_i_0*dim_W
	W_e_gas   =	ww_e_0*dim_W
	W_gas=(W_i_gas+W_e_gas)/beta
	omega_0 = l_in/(exp(x_in)*r_sch)**2.+2.*pi*(alpha/mass_ratio_factor)*W_gas/M_dot

	write(*,*)'omega_ratio',omega_0/omega_k(exp(x_in)*r_sch)
!	stop
	
	x_last = x_in

do 	i_test=2,N_step,1
	
!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	
!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	

!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	
!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	!!!!!!!!!!	

		i_test_0=i_test-1
	
40		sigp=sig(i_test_0)
		wip=wi(i_test_0)
		wep=we(i_test_0)
	
		x_test_0=x_now
		y_test_0(1)=sigp
		y_test_0(2)=wip
		y_test_0(3)=wep
	
	call deriva_test(df_test,x_test_0,y_test_0)
		k1_test=h_step_test*df_test
 		x_test=x_test_0+a2*h_step_test
		y_test=y_test_0+b21*k1_test	

		x_now=x_test
		sigp=y_test(1)
		wip=y_test(2)
		wep=y_test(3)
		
	call deriva_test(df_test,x_test,y_test)
		k2_test=h_step_test*df_test
		x_test=x_test_0+a3*h_step_test
		y_test=y_test_0+b31*k1+b32*k2

		x_now=x_test
		sigp=y_test(1)
		wip=y_test(2)
		wep=y_test(3)
				
		
	call deriva_test(df_test,x_test,y_test)
		k3_test=h_step_test*df_test
		x_test=x_test_0+a4*h_step_test
		y_test=y_test_0+b41*k1_test+b42*k2_test+b43*k3_test
	
		x_now=x_test
		sigp=y_test(1)
		wip=y_test(2)
		wep=y_test(3)
		
		
	call deriva_test(df_test,x_test,y_test) 
		k4_test=h_step_test*df_test
		x_test=x_test_0+a5*h_step_test
		y_test=y_test_0+b51*k1_test+b52*k2_test+b53*k3_test+b54*k4_test

		x_now=x_test
		sigp=y_test(1)
		wip=y_test(2)
		wep=y_test(3)


	call deriva_test(df_test,x_test,y_test)
		k5_test=h_step_test*df_test
		x_test=x_test_0+a6*h_step_test
		y_test=y_test_0+b61*k1_test+b62*k2_test+b63*k3_test+b64*k4_test+b65*k5_test

		x_now=x_test
		sigp=y_test(1)
		wip=y_test(2)
		wep=y_test(3)


	call deriva_test(df_test,x_test,y_test)
		k6_test=h_step_test*df_test
		y_test=c1*k1_test+c3*k3_test+c4*k4_test+c6*k6_test

		df0_test=y_test/h_step_test
		yerr_test=dc1*k1_test+dc3*k3_test+dc4*k4_test+dc5*k5_test+dc6*k6_test 

		errmax_test=0.d0
		errmax_test=max(dabs(yerr_test(1)),dabs(yerr_test(2)))
		errmax_test=max(errmax_test,dabs(yerr_test(3)))
		errmax_test=errmax_test/eacc_test

!		write(*,*)x_test,y_test,errmax_test,yerr_test
!stop

	if(errmax_test.gt.1.0d0)then
    	htemp_test=safety*h_step_test*(errmax_test**pshrnk)/10.0
    	htemp_test=sign(max(dabs(htemp_test),0.1*dabs(h_step_test)),h_step_test)
       	if(dabs(htemp_test).lt.drstep_min)then
       		write(*,*)'hstep underflow'
       		htemp_test=-drstep_min
       		stop
       	end if
    	h_step_test=htemp_test
    	write(*,*)'goto 40'
		goto 40
	else
   		if(errmax_test.gt.errcon)then
  			htemp_test=safety*h_step_test*(errmax_test**pgrow)
   		else
   				htemp_test=5.0*h_step_test
   		end if
   		if(dabs(htemp_test).gt.drstep_max)then
   			htemp_test=-drstep_max 
   		end if


		x_test=x_test_0+h_step_test
		sig(i_test)=sig(i_test_0)+y_test(1)
		wi(i_test)=wi(i_test_0)+y_test(2)
		we(i_test)=we(i_test_0)+y_test(3)

		x_now  = x_test
!		sigp = y_test(1)
!		wip  = y_test(2)
!		wep	 = y_test(3)

		
		rr=exp(x_test)
		ssigma=sig(i_test)
		ww_i=wi(i_test)
		ww_e=we(i_test)
		
		r=rr*r_sch
		sigma_gas = ssigma*dim_sigma
		W_i_gas   =	ww_i*dim_W
		W_e_gas   =	ww_e*dim_W
		W_gas=(W_i_gas+W_e_gas)/beta
	
	T_i=W_i_gas*mu_i*m_p/sigma_gas/k_bol
	T_e=W_e_gas*mu_e*m_p/sigma_gas/k_bol
	
	v_r=-M_dot/2./pi/sigma_gas/r
	c_s=sqrt((3.0*gamma-1.0+2.0*(gamma-1.0)*alpha*alpha)/(gamma+1.0)*W_gas/sigma_gas )

	
	if(-v_r/c_s>sonic_crit .and. - v_r/c_s<(2.-sonic_crit))then
!		pause
		write(*,*)'sony'!,x_test,x_test_0,i_test,h_step_test
!		write(*,*) sig(i_test),sig(i_test_0),df0_test(1),(h_step_test*10.)
		
!		pause
!		call crosssonic(i_test_0,df0_test,-1.5d-1)
!			stop
			x_test		=	x_test_0-sonic_jump               
			sig(i_test)	=	sig(i_test_0)+df0_test(1)*(-sonic_jump)
			wi(i_test)	= 	 wi(i_test_0)+df0_test(2)*(-sonic_jump)
			we(i_test)	= 	 we(i_test_0)+df0_test(3)*(-sonic_jump)

			rr=exp(x_test)
			ssigma=sig(i_test)
			ww_i=wi(i_test)
			ww_e=we(i_test)

			r=rr*r_sch
			sigma_gas = ssigma*dim_sigma
			W_i_gas   =	ww_i*dim_W
			W_e_gas   =	ww_e*dim_W
			W_gas=(W_i_gas+W_e_gas)/beta

			T_i=W_i_gas*mu_i*m_p/sigma_gas/k_bol
			T_e=W_e_gas*mu_e*m_p/sigma_gas/k_bol

			v_r=-M_dot/2./pi/sigma_gas/r
			c_s=sqrt( ((3.0*gamma-1.0+2.0*(gamma-1.0)*alpha*alpha)/(gamma+1.0))*W_gas/sigma_gas )
			x_now  = x_test
			
	end if
	write(11,*)log10(exp(x_test)),log10(sigma_gas)
	write(12,*)log10(exp(x_test)),W_i_gas
	write(13,*)log10(exp(x_test)),W_e_gas
	sigma_array(i_test-1) = sigma_gas
	W_i_array(i_test-1) = W_i_gas
	W_e_array(i_test-1) = W_e_gas
	write(21,*)exp(x_test)*r_sch,sigma_gas,W_i_gas,W_e_gas	
	write(14,*)log10(exp(x_test)),log10(T_i),log10(T_e)
!	stop
	write(15,*)log10(exp(x_test)),log10(-v_r/c),log10(c_s/c)
	write(16,*)log10(exp(x_test)),height(exp(x_test)*r_sch,ssigma*dim_sigma,ww_i*dim_W,ww_e*dim_W)/r,&
			sigma_gas*sigma_T/m_p
	
	omega=l_in/r**2.+2.*pi*(alpha/mass_ratio_factor)*W_gas/M_dot+2.*pi*cloud_1/r**2./M_dot

	write(17,*)log10(exp(x_test)),log10(omega*exp(x_test)*r_sch/c),log10(omega_k(exp(x_test)*r_sch)*exp(x_test)*r_sch/c)
	write(18,*)log10(exp(x_test)),-v_r/c_s
!	tao_cl = 500.	
	P_ADAF = W_gas/height(exp(x_test)*r_sch,ssigma*dim_sigma,ww_i*dim_W,ww_e*dim_W)/sqrt(2.*pi)	
	R_cl_max = 4.4d10*(mm_bh/1.d8)*(T_cl/1.d4)**0.5*(rr/10.)**1.5*2.
	
	R_cl_max = sqrt(4*pi**2.*(rr*r_sch)**3.*k_bol*T_cl/G/M_bh/m_p)
	
	R_cl_min = (5.4**(T_e/1.d10)**3.5/8./(P_ADAF/k_bol/T_cl/1.d14)**2.)**0.5*1.d10
	n_cl = P_ADAF/k_bol/T_cl	
	
	qq_cl = ((sqrt(8./pi))*sigma_gas*det_VR2**(3./2.)*f_R - (sqrt(8./pi))*sigma_gas*det_VPhi2**(3./2.)*f_phi) &
				/height(exp(x_test)*r_sch,ssigma*dim_sigma,ww_i*dim_W,ww_e*dim_W)

	M_cl = (P_ADAF/k_bol/T_cl)*m_p*(4.*pi/3.)*R_cl**3.		
	tao_cl = (P_ADAF/k_bol/T_cl)*sigma_T*R_cl
!	R_cl = (R_cl_max*R_cl_min)**0.5
	R_cl = (5.4d8*(T_e/1.d10)**3.5/(-1.d-26*n_cl*(1.-n_cl*(1.d7)*exp(-118400./(T_cl+1.d3))+1.4d-2*sqrt(T_cl)*exp(-92./T_cl))))**0.5*1.d10

	R_cl = R_cl_out*(rr/1.d3)**(cloud_size_index)

	rr_tidal = (G*M_bh*m_p*R_cl**2./3./k_bol/T_cl)**(1./3.)/r_sch


		write(19,'(F8.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3,ES17.3)')(exp(x_test)),&!log10(exp(x_test)),&
			(P_ADAF/k_bol/T_cl)*m_p*(4.*pi/3.)*R_cl**3.,&
			(P_ADAF/k_bol/T_cl)*sigma_T*R_cl,&
			8.d8*(P_ADAF/k_bol/T_cl/1.d14)**2.,&		! 	lambda_line cooling rate erg/s/cm3
			5.4d8*(T_e/1.d10)**3.5/(R_cl/1.d10)**2.,&	!	Heating  rate		erg/s/cm3
			R_cl_max/r_sch,&			!	Rcl max
			R_cl_min/r_sch,&	! Rcl min
			(P_ADAF/k_bol/T_cl)*sigma_T*R_cl_max*2.,&	!	tao_max
			(P_ADAF/k_bol/T_cl)*sigma_T*R_cl_min*2.,&	!	tao_min
			(P_ADAF/k_bol/T_cl)*sigma_T*(R_cl_max*R_cl_min)**0.5,&
			sigma_SB*T_cl**4.*(n_cl*sigma_T),&
			qq_cl,&
			(5.4d8*(T_e/1.d10)**3.5/(sigma_SB*T_cl**4.*n_cl*sigma_T + 8.8d8*(P_ADAF/k_bol/T_cl/1.d14)**2.*(T_cl/1.d4)**(-1.3)))**0.5*1.d10,&
			n_cl*sigma_T*(5.4d8*(T_e/1.d10)**3.5/(sigma_SB*T_cl**4.*n_cl*sigma_T + 8.8d8*(P_ADAF/k_bol/T_cl/1.d14)**2.*(T_cl/1.d4)**(-1.3)))**0.5*1.d10,&
			n_cl*sigma_T*R_cl,&
			3.3d16*((T_e/1.d8)**(7./4.)/(P_ADAF/1.d-2)/10.),&
		 	(G*M_bh*m_p*R_cl_max**2./4./pi/k_bol/T_cl)**(1./3.)/r_sch,&
			8.*(R_cl_max/1.d11)**(2./3.),&
			sqrt(4.*pi*k_bol*T_cl*(rr*r_sch)**3./G/M_bh/m_p)/r_sch,&
			R_cl/r_sch,&
			R_cl_out*(rr/1.d3)**(0.9)/r_sch,&
			n_cl*sigma_T*R_cl,&
			n_cl*sigma_T*R_cl_out*(rr/1.d3)**(0.9),&
			rr_tidal
			
			if ((G*M_bh*m_p*R_cl_max**2./4./pi/k_bol/T_cl)**(1./3.)/r_sch > exp(x_test))then
!				stop
			end if
			
	if(R_cl_max< 1.d12)then
		write(52,*)	log10(exp(x_test)),	R_cl_max, (P_ADAF/k_bol/T_cl)*sigma_T*R_cl_max
		tao_cl = (P_ADAF/k_bol/T_cl)*sigma_T*R_cl_max
	else
		write(52,*)	log10(exp(x_test)),	1.d12, (P_ADAF/k_bol/T_cl)*sigma_T*1.d12
		tao_cl = (P_ADAF/k_bol/T_cl)*sigma_T*1.d12
	end if
			
	
!	if (W_i_gas/height(exp(x_test)*r_sch,ssigma*dim_sigma,ww_i*dim_W,ww_e*dim_W)/sqrt(2.*pi)<G*M_bh*(m_p/sigma_T)*tao_cl*R_cl/4./pi/r**3.)then
!		write(*,*)exp(x_test)
!		stop
!	end if
	
	
!,(sigma_gas)/m_p*k_bol*T_i

!	p_ADAF = (sigma_gas/height(exp(x_test)*r_sch,ssigma*dim_sigma,ww_i*dim_W,ww_e*dim_W))/m_p*k_bol*T_i

	write(*,*)	-v_r/c_s,i_test,log10(exp(x_test)) 

	h_step_test=htemp_test

	end if
	
		sigma_gas_0 = ssigma_0*dim_sigma
		W_i_gas_0   =	ww_i_0*dim_W
		W_e_gas_0   =	ww_e_0*dim_W
!		W_gas_0=(W_i_gas_0+W_e_gas_0)/beta
!	
!	T_i_0=W_i_gas_0*mu_i*m_p/sigma_gas_0/k_bol
!	T_e_0=W_e_gas_0*mu_e*m_p/sigma_gas_0/k_bol
!	
!	v_r_0=M_dot/2./pi/sigma_gas_0/r
!	c_s_0=sqrt((3.0*gamma-1.0+2.0*(gamma-1.0)*alpha*alpha)/(gamma+1.0)*W_gas_0/sigma_gas_0 )
!	
	dsigma_dR = (sigma_gas-sigma_gas_0)/(x_test-x_last)/exp(x_test)/r_sch
!	
!	l

	f_R = gamma_R/r_sch
	f_phi = gamma_phi/r_sch

	v_phi_2 = (omega*r)**2.+v_r*(2.*r*omega + r**2.*(omega-omega_0)/(x_test-x_last)/exp(x_test)/r_sch)/r/f_phi
	
!	write(*,*)v_phi_2, h_step_test,dsigma_dR


	if(-v_r/c_s>sonic_crit .and. - v_r/c_s<(2.-sonic_crit))then
		v_R_2 = v_R_2_0 + r_sch*exp(x_test)*(-sonic_jump)*(&
					- v_R_2_0*dsigma_dR/sigma_gas - v_R_2_0/r+f_R*(v_R_2_0) + v_phi_2/R - G*M_bh/R/R - f_R*(M_dot/2./pi/sigma_gas/r)**2.)
			else
	v_R_2 = v_R_2_0 + r_sch*exp(x_test)*(x_test-x_last)*(&
				- v_R_2_0*dsigma_dR/sigma_gas - v_R_2_0/r+f_R*(v_R_2_0) + v_phi_2/R - G*M_bh/R/R - f_R*(M_dot/2./pi/sigma_gas/r)**2.)
	end if

	det_VR2 = v_R_2 - v_r**2.
	det_VPhi2 = v_phi_2 - (r*omega)**2.
	v_phi_ratio_v_R = -r*omega/v_r
	
	write(30,'(es17.7,es17.7,es17.7,es17.7,es17.7,es17.7,es17.7,es17.7)')&
	v_R_2,dsigma_dR,-v_R_2_0*dsigma_dR/sigma_gas, v_R_2_0/r,f_R*(v_R_2_0),v_phi_2/R, G*M_bh/R/R ,f_R*(M_dot/2./pi/sigma_gas/r)**2.
	
	write(31,*)log10(exp(x_test)), log10(v_R_2**0.5/c), log10(v_phi_2**0.5/c)
		
	write(32,*)log10(exp(x_test)), omega,((omega-omega_0)/h_step_test/exp(x_test))!dsigma_dR
	write(33,*)log10(r/r_sch),log10(v_phi_2/c/c),log10(v_R_2/c/c)
	write(34,*)log10(exp(x_test)),h_step_test,-(x_test-x_last)
	
	
	
	write(51,'(F17.7,F17.7,F17.7,F17.7,F17.7,F17.7,F17.7,F17.7,F17.7,F17.7,F17.7,F17.7)')&
	log10(exp(x_test)),log10(sigma_gas),log10(W_i_gas),log10(W_e_gas),&
	log10(T_i),log10(T_e),log10(-v_r/c),log10(c_s/c),log10(omega*r/c),log10(omega_k(r)*r/c),&
	log10(v_R_2**0.5/c), log10(v_phi_2**0.5/c)
	
!	if(exp(x_test).le.8.)then
	if(x_test.le.x_out)then
		write(*,*)'x.le.x_in',i_test
		exit
!	else
!		if(error_allocate .eq. 0)then
!		allocate(ssigma_array(i_test),stat=error_allocate)
!		ssigma_array(i_test) = sigma_gas/dim_sigma
!		write(*,*)ssigma_array(i_test)
!		end if
	end if

ssigma_0 = ssigma
ww_i_0 	 = ww_i
ww_e_0   = ww_e
v_R_2_0 = v_R_2
omega_0 = omega	

x_last = x_test

r_test(i_test-1) = r

!do i_mu = 1, 100
!	mu = 10.**(9. + i_mu*(22. - 9.)/100.)
!	f_r_mu(i_test,i_mu) = integ_for_q(mu,r,sigma_gas,W_i_gas,W_e_gas)/2./yita_compton(mu,r,sigma_gas,W_i_gas,W_e_gas)
!	write(58,*) r_test(i_test-1), mu, f_r_mu(i_test,i_mu)
!end do

r_array(i_test-1)=r
T_e_array(i_test-1)=T_e
tau_e_array(i_test-1)=sigma_gas*sigma_t/m_p
do i=1,200,1
	flux(i,i_test-1)=flux_nu(mu(i),r,sigma_gas,W_i_gas,W_e_gas)
end do


end do

	close(58)
	close(21)
	
!	write(*,*)f_r_mu
	
!	f_mu = 0.
!do i_mu = 1,100
!	do i_test = 2,1358
!	f_mu(i_mu) = f_mu(i_mu) + f_r_mu(i_test,i_mu)*2.*pi*r_test(i_test)*(-1)*(r_test(i_test)-r_test(i_test-1))
!	write(23,*)r_test(i_test),r_test(i_test)-r_test(i_test-1),2.*pi*r_test(i_test)*(-1)*(r_test(i_test)-r_test(i_test-1))
!	end do	
!	write(21,*)10.**(9. + i_mu*(22. - 9.)/100.), f_mu(i_mu)!*10.**(9. + i_mu*(22. - 9.)/100.)
!end do

do i=1,N_radius-1,1
	if(mod(i,10) .eq. 1)then
	write(39,*)r_array(i)
	end if
end do

open(unit=34,file='spec_input.txt')
do j=1,N_radius-1,1
	if(mod(j,10) .eq. 1)then
	do i=1,200,1
		write(34,*)T_e_array(j),tau_e_array(j),flux(i,j)
	end do
	end if
end do
close(34)



end subroutine ADAF_solve


function beta_gamma(gamma_factor)
	implicit none
	real*8		::		beta_gamma,gamma_factor
	beta_gamma = 1.-1./gamma_factor**2.
	return
end function beta_gamma

subroutine test_omega_integ()
	implicit none
	real*8		::		omega_incoming, gamma_factor,omega_ave,omega_ave_2,R_omega_gamma
        
	omega_incoming = 20901553236.7367
	gamma_factor = 1.35043649799882
	call omega_integ(omega_incoming, gamma_factor,omega_ave,omega_ave_2,R_omega_gamma)
	write(*,*)omega_ave,omega_ave_2,R_omega_gamma
!	stop
	return
end subroutine test_omega_integ

subroutine omega_integ(omega_incoming, gamma_factor,omega_ave,omega_ave_2,R_omega_gamma)
	use constant
	implicit none
	integer*8,parameter	::	n = 100
	integer*8	::		i_mu,i
	real*8		::		xx(n),ww(n)
	real*8		::	omega_incoming,omega_ave,gamma_factor,omega_alpha_2
	real*8		::	BETA_GAMMA,R_omega_gamma,R_alpha,omega_mean,omega_mean_2,omega_ave_2
	real*8		::	x1,x2,mu,omega_alpha,alpha,x,xp_over_x,d_sigma_d_alpha,P_mu_alpha_omega_gamma
	
	beta_gamma = sqrt(1.-1./gamma_factor**2.)

	x1=-1.
	x2=+1.
    call gauleg(x1,x2,xx,ww,n)
	omega_ave = 0.
	omega_ave_2 = 0.
	R_omega_gamma = 0.
	do i_mu = 1,n,1
		mu = xx(i_mu)
		
		omega_alpha = 0.
		omega_alpha_2 = 0.
		R_alpha = 0.
		do i = 1,n,1
		alpha = xx(i)	
		
		x = gamma_factor*omega_incoming*(1.-beta_gamma*mu)
		xp_over_x = 1./(1.+(1.-alpha)*x)
		d_sigma_d_alpha = (3./8.)*sigma_T*(xp_over_x**2.)*(1./xp_over_x+xp_over_x-1.+alpha**2.)
		P_mu_alpha_omega_gamma = c*(1.-beta_gamma*mu)*d_sigma_d_alpha
		omega_mean = gamma_factor*omega_incoming*(gamma_factor*(1.-beta_gamma*mu) + gamma_factor*beta_gamma*alpha*(mu - beta_gamma))*xp_over_x
		omega_mean_2 = gamma_factor*omega_incoming**2.&
			*((gamma_factor*(1.-beta_gamma*mu)+gamma_factor*beta_gamma*alpha*(mu-beta_gamma))**2.+ 0.5*beta_gamma**2.*(1.-alpha**2.)*(1.-beta_gamma**2.))&
			/(1.+gamma_factor*(1.-beta_gamma*mu)*(1-alpha)*omega_incoming)**2.

		omega_alpha = omega_alpha + ww(i)*omega_mean * d_sigma_d_alpha
		omega_alpha_2 = omega_alpha_2 + ww(i)*omega_mean_2 * d_sigma_d_alpha

		R_alpha = R_alpha + ww(i)*d_sigma_d_alpha

		end do
		omega_ave = omega_ave + ww(i_mu)*omega_alpha*(1-beta_gamma*mu)*c/2.
		omega_ave_2 = omega_ave_2 + ww(i_mu)*omega_alpha_2*(1-beta_gamma*mu)*c/2.
		R_omega_gamma = R_omega_gamma + ww(i_mu)*R_alpha*(1-beta_gamma*mu)*c/2.		
	end do

	omega_ave = omega_ave/R_omega_gamma
	omega_ave_2 = omega_ave_2/R_omega_gamma
	
!    call gauleg(x1,x2,xx,ww,n)

!    result=0.
!    do i=1,n,1
!    result=result+ww(i)*func(xx(i))
!    end do
	
	return
end subroutine omega_integ

subroutine n_omega_in_out()
	use constant
	implicit none
	integer*8,parameter	::	n = 100
	integer*8	::		i_omega,i_gamma
	real*8		::		xx_gamma(n),ww_gamma(n),xx_omega(n),ww_omega(n),x1,x2
	real*8		::		omega_out,f_mu_out,omega_int,f_mu_alpha,gamma_factor
	real*8		::		omega_ave,R_omega_gamma,P_omega_omega_gamma
	real*8		::		delta_omega_2,omega_ave_2,D_omega_omega,beta_gamma
	real*8		::		omega_min,omega_max,theta_e
	
	real*8		::		T_e
	
	real*8,external	::	BESSK,f_mu_input

!	T_i=W_i_in*mu_i*m_p/k_bol/sigma_in
!	T_e=W_e_in*mu_e*m_p/k_bol/sigma_in
	
!	theta_i=k_bol*T_i/m_p/c/c
	T_e = 10.**9.6
	theta_e=k_bol*T_e/m_e/c/c
	
!	n_e=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_e
!	n_i=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_i
	

	omega_min = 0.

	omega_out = 10.**10.

	write(*,*)'n_omega_in_out',omega_out,f_mu_out


!	write(*,*)omega_out

	x1=20.7233
	x2=50.6569
    call gauleg(x1,x2,xx_omega,ww_omega,n)

	omega_min = exp(x1)
	omega_max = exp(x2)

    x1=dlog10(1.001d0)
    x2=dlog10(50.0d0)	
    call gauleg(x1,x2,xx_gamma,ww_gamma,n)

	f_mu_out = 0.
	do i_omega = 1, n,1
		omega_int = exp(xx_omega(i_omega))
		f_mu_alpha = 0.
		do i_gamma = 1, n,1
			gamma_factor = exp(xx_gamma(i_gamma))			
			beta_gamma = sqrt(1.-1./gamma_factor**2.)
			
			write(*,*)omega_int	
			call omega_integ(omega_int,gamma_factor,omega_ave,omega_ave_2,R_omega_gamma)	
			write(*,*)'omega_int',omega_int,gamma_factor,omega_ave,omega_ave_2,R_omega_gamma
			call test_omega_integ(omega_int,gamma_factor)
!			stop			
			delta_omega_2=omega_ave_2 - omega_ave**2.

			if((sqrt(3.*delta_omega_2) < (omega_ave - omega_min)) .and. (sqrt(3.*delta_omega_2) < (omega_max - omega_ave)))then
				D_omega_omega = sqrt(3.*delta_omega_2)
			else if	((omega_ave - omega_min) < (sqrt(3.*delta_omega_2))&
				.and.&
				(omega_ave - omega_min) < (omega_max - omega_ave))then
				D_omega_omega = (omega_ave - omega_min)
			else
				D_omega_omega = (omega_max - omega_ave)
			end if
			
			if(D_omega_omega > abs(omega_out - omega_ave))then
				P_omega_omega_gamma = 1./2./D_omega_omega
			else
				P_omega_omega_gamma = 0.
			end if

			write(*,*)D_omega_omega,P_omega_omega_gamma

			f_mu_alpha = f_mu_alpha + ww_gamma(i_gamma)*P_omega_omega_gamma*(R_omega_gamma/c/sigma_T)&
					*(gamma_factor**2.*beta_gamma*exp(-gamma_factor/theta_e)/theta_e/bessk(2,1/theta_e))&
					*(f_mu_input(omega_int)/h/omega_int)*omega_int*gamma_factor
			write(*,*)'f_mu_input',	ww_gamma(i_gamma),P_omega_omega_gamma*(R_omega_gamma/c/sigma_T)&
							*(gamma_factor**2.*beta_gamma*exp(-gamma_factor/theta_e)/theta_e/bessk(2,1/theta_e))&
							*(f_mu_input(omega_int)/h/omega_int)*omega_int*gamma_factor
					
		end do
!		stop
		f_mu_out = f_mu_out + ww_omega(i_omega)*f_mu_alpha
	end do

		f_mu_out = f_mu_out!*(log(10.)**2.)*(10.**(-0.8)/m_p)*(0.5*1.d15)*sigma_T

	write(*,*)'n_omega_in_out',omega_out,f_mu_out

return
end subroutine n_omega_in_out



function f_mu_input(omega_int)
	implicit none
	real*8	::		f_mu_input,omega_int
	
	f_mu_input = 1.d40
	return
end function f_mu_input



program none
		implicit none
		call INITIALIZATION_bh()
		call INITIALIZATION_clouds()
!		call test_omega_integ()
!		call n_omega_in_out()
!		stop
		call ADAF_solve()
!		call n_omega_in_out()
!		call cloud_solve()
end program none
	
subroutine 	cloud_solve()
	use constant
	use vars_bh
	use vars_adaf
	use vars_clouds
	implicit none
	call interpolation()
	return
end subroutine cloud_solve
				
subroutine interpolation()
use constant
use vars_bh
use vars_adaf
use vars_clouds
implicit none

!write(*,*)sigma_array
!write(*,*)W_i_array
!write(*,*)W_e_array

return
end subroutine interpolation
				
					
subroutine deriva_test(df_test,x_test,y_test)	
use constant
use vars_bh
use vars_adaf
implicit none
	real*8		::		x_test,y_test(3),df_test(3)
	real*8		::		coefficient_a(3,3),coefficient_c(3)
	real*8		::		W_i_in,W_e_in,sigma_in,r
	
!	write(*,*)'derivs',exp(x),y
!	stop
	
	
	
	r=exp(x_test)*r_sch
	sigma_in=y_test(1)*M_dot_edd/c/r_sch
	W_i_in=y_test(2)*M_dot_edd*c/r_sch
	W_e_in=y_test(3)*M_dot_edd*c/r_sch

	call coefficient(r,sigma_in,W_i_in,W_e_in,coefficient_a,coefficient_c)
	call tri_linear(coefficient_a,df_test,coefficient_c)
	
	return
end subroutine deriva_test

subroutine derivs(x,y,dydx)
use constant
use vars_bh
use vars_adaf
implicit none
	real*8		::		x,y(3),dydx(3)
	real*8		::		coefficient_a(3,3),coefficient_c(3)
	real*8		::		W_i_in,W_e_in,sigma_in,r
	
!	write(*,*)'derivs',exp(x),y
!	stop
	
	r=exp(x)*r_sch
	sigma_in=y(1)*M_dot_edd/c/r_sch
	W_i_in=y(2)*M_dot_edd*c/r_sch
	W_e_in=y(3)*M_dot_edd*c/r_sch

	call coefficient(r,sigma_in,W_i_in,W_e_in,coefficient_a,coefficient_c)
	call tri_linear(coefficient_a,dydx,coefficient_c)
	
	return
end subroutine derivs

subroutine	coefficient(r,sigma_in,W_i_in,W_e_in,coefficient_a,coefficient_c)
	use constant
	use vars_bh
	use vars_adaf
	use vars_clouds
	implicit none
	real*8	::	coefficient_a(3,3)
	real*8	::	coefficient_c(3)
	real*8	::	W_in,W_i_in,W_e_in,sigma_in,r,rr
	real*8	::	ww_i,ww_e,ww,ssigma,x,ll_in,oomega,oomega_k,qq_rad,llambda_ie,ddlogomegakdr,f_R,f_phi
	real*8	::	dim_W,dim_sigma,dim_l_in,dim_Q,dim_omega

	real*8,external		::		q_rad,lambda_ie,omega_k,dlogomegakdr,height

	f_R = gamma_R/r_sch
	f_phi = gamma_phi/r_sch

!	initialization
!	l_in,	Omega,	Omega_k,	dlnomega_kdr
!	qrad,	lambda_ie
!
	W_in=(W_i_in+W_e_in)/beta
	
!	dimensionless

	omega=l_in/r**2.+2.*pi*(alpha/mass_ratio_factor)*W_in/M_dot!+2.*pi*cloud_1/r**2./M_dot

	dim_W=M_dot_edd*c/r_sch
	dim_sigma=M_dot_edd/c/r_sch
	dim_omega=c/r_sch
	dim_l_in=c*r_sch
	dim_Q=M_dot_edd*c*c/r_sch/r_sch
	
!	write(*,*)'coeffecient: sigma_in,W_i_in,W_e_in',sigma_in,W_i_in,W_e_in

	ww_i=W_i_in/dim_W
	ww_e=W_e_in/dim_W
	ww=W_in/dim_W
	ssigma=sigma_in/dim_sigma
	ll_in=l_in/dim_l_in
	rr=r/r_sch
	oomega=omega/dim_omega
!		write(*,*)'omega',omega
	oomega_k=omega_k(r)/dim_omega
	qq_rad=q_rad(r,sigma_in,W_i_in,W_e_in)/dim_Q
	
	T_i=W_i_in*mu_i*m_p/sigma_in/k_bol
	T_e=W_e_in*mu_e*m_p/sigma_in/k_bol

	llambda_ie=lambda_ie(r,sigma_in,W_i_in,W_e_in)/dim_Q
	ddlogomegakdr=dlogomegakdr(r)*r_sch
	
!	write(*,*)'coeffecient:q',dim_Q,qq_rad,llambda_ie
!	stop

	x=log(rr)

	
	coefficient_a(1,1)=-mm_dot**2./4./pi**2./exp(3.*x)/ssigma**2.
	coefficient_a(1,2)=1./beta/exp(x)
	coefficient_a(1,3)=1./beta/exp(x)
	
	coefficient_a(2,1)=-ww_i*mm_dot*(3.*gamma-1.)/ssigma**2./exp(x)/2./(gamma-1.)
	coefficient_a(2,2)=mm_dot*(gamma+1.)/ssigma/exp(x)/2./(gamma-1.)&
						-2.*pi*exp(x)*alpha**2.*ww*2.*pi/mm_dot/beta/mass_ratio_factor
	coefficient_a(2,3)=-2.*pi*exp(x)*alpha**2.*ww*2.*pi/mm_dot/beta/mass_ratio_factor

	coefficient_a(3,1)=-ww_e*mm_dot*(3.*gamma-1.)/ssigma**2./exp(x)/2./(gamma-1.)
	coefficient_a(3,2)=0.
	coefficient_a(3,3)=mm_dot*(gamma+1.)/ssigma/exp(x)/2./(gamma-1.)

	cloud_2 = ssigma*gamma_R*(v_R_2/c/c - (mm_dot/2./pi/ssigma/rr)**2.)

	write(41,*)r,ssigma*gamma_R*((v_R_2/c/c) - (mm_dot/2./pi/ssigma/rr)**2.)
	write(42,*)r,v_R_2/c/c, (mm_dot/2./pi/ssigma/rr)**2.
!!!!!!!!!!

	qq_cl = (-f_R*det_VR2-f_phi*det_VPhi2*v_phi_ratio_v_R)*M_dot/2./pi/r/dim_Q
	
	qq_cl = (sqrt(8./pi))*sigma_in*det_VR2**(3./2.)*f_R/dim_Q - (sqrt(8./pi))*sigma_in*det_VPhi2**(3./2.)*f_phi/dim_Q
!	write(*,*)det_VR2,det_VPhi2

	


!if(x<log(rr_tidal))then
	cloud_2 = 0.
	qq_cl = 0.
!end if

	
	
!	sigma_in,-sigma_in*k_bol*T_i/height(exp(x)*r_sch,sigma_in,W_i_in,W_e_in)/m_p/dim_Q/qq_cl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	coefficient_c(1)=	mm_dot**2./4./pi**2./exp(3.*x)/ssigma&
						+ssigma*exp(x)*(oomega**2.-oomega_k**2.)&
						-ww*ddlogomegakdr!&
!						-cloud_2
												
	coefficient_c(2)=	mm_dot*ww_i*ddlogomegakdr/ssigma&
						-2.*pi*alpha*ww*2.*ll_in/exp(x)&
						+2.*pi*exp(x)*llambda_ie!&
!						+2.*pi*exp(x)*qq_cl*0.05

	coefficient_c(3)=	mm_dot*ww_e*ddlogomegakdr/ssigma&
						+2.*pi*exp(x)*qq_rad&
						-2.*pi*exp(x)*llambda_ie!&
!						-2.*pi*exp(x)*qq_cl*0.5
						
!		write(*,*)'coeffi_c:3',mm_dot*ww_e*ddlogomegakdr/ssigma,2.*pi*exp(x)*qq_rad&
!								,2.*pi*exp(x)*llambda_ie+k_c_term
	write(99,'(es17.7,es17.7,es17.7,es17.7,es17.7,es17.7)')&
		log10(r/r_sch),&
		(dim_Q*qq_rad),&
		(dim_Q*llambda_ie),&
		(dim_Q*qq_cl*0.1),&
		-((mm_dot*ww_i*ddlogomegakdr/ssigma)/(2.*pi*exp(x)))*dim_Q,&!-(dim_Q*qq_cl)-(dim_Q*llambda_ie)
		-((mm_dot*ww_i*ddlogomegakdr/ssigma)/(2.*pi*exp(x)))*dim_Q+(dim_Q*qq_cl*0.1)-(dim_Q*llambda_ie)
		
		
		
end subroutine coefficient

subroutine	tri_linear(aa,bb,cc)
	implicit none
	real*8		::	aa(3,3),bb(3),cc(3)
	real*8		::	delta,delta1,delta2,delta3
	              
!	write(*,*)'aa',aa
!	write(*,*)'cc',cc
	           
	delta=aa(1,1)*aa(2,2)*aa(3,3)&
		+ aa(1,2)*aa(2,3)*aa(3,1)&
		+ aa(2,1)*aa(3,2)*aa(1,3)&
		- aa(1,3)*aa(2,2)*aa(3,1)&
		- aa(1,2)*aa(2,1)*aa(3,3)&
		- aa(1,1)*aa(2,3)*aa(3,2)
		
	delta1=	  cc(1)*aa(2,2)*aa(3,3)&
			+ aa(1,2)*aa(2,3)*cc(3)&
			+ cc(2)*aa(3,2)*aa(1,3)&
			- aa(1,3)*aa(2,2)*cc(3)&
			- aa(1,2)*cc(2)*aa(3,3)&
			- cc(1)*aa(2,3)*aa(3,2)
			
    delta2=   aa(1,1)*cc(2)*aa(3,3)&
			+ cc(1)*aa(2,3)*aa(3,1)&
			+ aa(2,1)*cc(3)*aa(1,3)&
			- aa(1,3)*cc(2)*aa(3,1)&
			- cc(1)*aa(2,1)*aa(3,3)&
			- aa(1,1)*aa(2,3)*cc(3)		
		
	delta3=   aa(1,1)*aa(2,2)*cc(3)&
			+ aa(1,2)*cc(2)*aa(3,1)&
			+ aa(2,1)*aa(3,2)*cc(1)&
			- cc(1)*aa(2,2)*aa(3,1)&
			- aa(1,2)*aa(2,1)*cc(3)&
			- aa(1,1)*cc(2)*aa(3,2)

!	write(*,*)"delta",delta,delta1,delta2,delta3
	

	

			bb(1)=delta1/delta
			bb(2)=delta2/delta
			bb(3)=delta3/delta
!	delta_ex=delta
return
end subroutine tri_linear

subroutine INITIALIZATION_func(x_in,x_out,y_in)
use constant
	implicit none
	real*8		::		x_in,x_out,y_in(3)
	real*8		::		rr_out

	call boundary_condition_initial(y_in(1),y_in(2),y_in(3),rr_out)
	x_in=log(rr_out)
	x_out=log(1.01)
	
	write(*,*)log10(y_in(2)*c*c*1.23*m_p/y_in(1)/k_bol)
!stop
	return
end subroutine INITIALIZATION_func

subroutine	boundary_condition_initial(ssigma_out,ww_i_out,ww_e_out,rr_out)

	use constant
	use vars_bh
	use vars_adaf
	use vars_clouds
	use vars_boundary
	implicit none
	real*8		::			dim_W,dim_sigma,dim_l_in,dim_Q,dim_omega
	real*8		::			ll_in,rr_out,ww_out,ww_i_out,ww_e_out,ssigma_out
	real*8		::			R_out,cloud_1_out,T_vir_out,W_out
	real*8		::			W_i_out,W_e_out,T_i_out,T_e_out,sigma_out,l_out,omega_out

	real*8,external		::		height,omega_k

	r_out=1.d4*r_sch
	
	dim_W=M_dot_edd*c/r_sch
	dim_sigma=M_dot_edd/c/r_sch
	dim_omega=c/r_sch
	dim_l_in=c*r_sch
	dim_Q=M_dot_edd*c*c/r_sch/r_sch
	
	ll_in=l_in/dim_l_in
	rr_out=r_out/r_sch
!	write(*,*)'rr_out',rr_out,r_out,r_sch
!	stop
	cloud_1_out=0.
	
	ww_out = Omega_boundary_factor*mm_dot/2./pi/(alpha/mass_ratio_factor)/(rr_out-1.)/sqrt(2.*rr_out)&
		-mm_dot*ll_in/2./pi/(alpha/mass_ratio_factor)/rr_out**2.
	
	ssigma_out=2.*beta*ww_out*rr_out/T_out_factor/(gamma-1.)*mu_i*mu_e/(mu_i+mu_e)
	
	ww_i_out=beta*ww_out*mu_e/(mu_i+mu_e)
	ww_e_out=beta*ww_out*mu_i/(mu_i+mu_e)

	W_i_out=ww_i_out*dim_W
	W_e_out=ww_e_out*dim_W
	sigma_out=ssigma_out*dim_sigma
	
	T_i_out=W_i_out*mu_i*m_p/sigma_out/k_bol
	T_e_out=W_e_out*mu_e*m_p/sigma_out/k_bol

!write(*,*)0.8*1./(rr_out-1.)/sqrt(2.*rr_out)*dim_omega,(ww_i_out+ww_e_out)/beta,sigma_out/dim_sigma,T_i_out,T_e_out,W_i_out+W_e_out!,0.1*(gamma-1.)*G*M_bh*m_p/k_bol/r_out

	omega_out=0.8*Omega_k(r_out)
	l_out=omega_out*r_out**2.
	W_out=M_dot*(l_out)/2./pi/r_out**2./(alpha/mass_ratio_factor)-M_dot*(l_in)/2./pi/r_out**2./(alpha/mass_ratio_factor)
	
	T_vir_out=(gamma-1.)*G*M_bh*m_p/k_bol/r_out
	T_i_out=0.1*T_vir_out
	T_e_out=0.1*T_vir_out
	
	write(*,*)log10(T_i_out)
!	stop
	
	sigma_out=m_p*beta*W_out*mu_i*mu_e/(mu_i+mu_e)/k_bol/T_i_out
	write(*,*)'omega ratio',omega_out*r_out/(omega_k(r_out)*r_out)
	
	write(*,*)W_out/dim_W,sigma_out/dim_sigma,T_i_out,T_e_out,W_out*beta
	write(*,*)r_out
!	stop

	return
end subroutine	boundary_condition_initial

subroutine	INITIALIZATION_bh()
	use constant
	use vars_bh
	use vars_adaf
	use vars_clouds
	implicit none
!	M_dot_edd=32.*pi*c*r_sch/0.34

	M_dot_edd=1.39d18*mm_bh
	M_dot=mm_dot*M_dot_edd


	M_bh=mm_bh*m_sun
!	r_g=G*M_bh/c/c
	r_sch=2.*G*M_bh/c/c

	
!	rr_tidal = (G*M_bh*m_p*R_cl**2./3./k_bol/T_cl)**(1./3.)/r_sch

	
!	M_dot=mm_dot*M_dot_edd
	
!	for Li's Boundary & Li's BH mass, beta	
!	l_in=0.8999999*r_sch*c

! 	for Manmoto's boundary & Li's mass, beta
!	l_in=0.904*r_sch*c
	
!	for M_dot_edd=32.*pi*c*r_sch/0.34
!	l_in=0.554*r_sch*c
!

!	for beta = 0.5	
!	l_in=0.5434*r_sch*c

!	l_in = 0.855*r_sch*c

!	gamma_R = 1.d-3
!	gamma_phi = 1.d-3
	l_in = 0.743*r_sch*c

!l_in = 0.535*r_sch*c

l_in = 0.9182115*r_sch*c
!l_in = 0.
!l_in = 1.1009*r_sch*c

!!!!!!!!!!!!!!!!!!!!!!!!

!	f_cl = 0.1:0.9
l_in = 0.8193*r_sch*c

!	f_cl = 0.8:0.2
l_in = 0.81857*r_sch*c

!	f_cl = 0.2:0.8
	l_in = 0.8129695*r_sch*c
	
	
	
!	f_cl = 0.5 :
	l_in = 0.81334146855*r_sch*c
	sonic_jump=0.01d-1
	sonic_crit = 0.9999!994


	l_in = 0.9317*r_sch*c
	sonic_jump=0.5d-1
	sonic_crit = 0.99!99

!0.8 : 0.3
l_in = 1.005*r_sch*c
sonic_jump=0.4d-1
sonic_crit = 0.98!994

l_in = 0.8565*r_sch*c

l_in = 0.59*r_sch*c

l_in=0.536*r_sch*c

!l_in = 2*r_sch*c

sonic_jump=1.d-1
sonic_crit = 0.98!994


end subroutine	INITIALIZATION_bh



subroutine INITIALIZATION_clouds()
	use constant
	use vars_bh
	use vars_adaf
	use vars_clouds
	implicit none
	
	k_c_cl=0.
	cloud_1=0.
	cloud_2=0.
	cloud_3=0.
	return
end subroutine INITIALIZATION_clouds


function omega_k(r)
	use constant
	use vars_bh
	implicit none
	real*8					::	omega_k
	real*8,intent(in)		::	r
	omega_k= sqrt(G*M_bh/r/(r-r_sch)**2.)
	return
end function omega_k



function dlogomegakdr(r)
	use constant
	use vars_bh
	implicit none
	real*8		::		dlogomegakdr,r
	dlogomegakdr=-(3.*r-r_sch)/(2.*r)/(r-r_sch)
	return
end function dlogomegakdr


function height(r,sigma_in,W_i_in,W_e_in)
	use constant
	use vars_bh
	use vars_adaf
	implicit none
	real*8		::		height,r
	real*8		::		sigma_in,W_i_in,W_e_in,W_in
	
	real*8,external	::		omega_k
	
	W_in=(W_i_in+W_e_in)/beta
	cs=(W_in/sigma_in)**0.5

	height=cs/omega_k(r)
!	write(*,*)'height',height,cs/c,omega_k(r)*r_sch/c,r/r_sch,sigma_in,W_i_in,W_e_in
!	stop
	return
end function height


subroutine qgaus(func,x1,x2,result)
    implicit none
    real*8    x1,x2,result,s,t
    integer    i,n
    
    parameter(n=200)		!高斯积分的节点数
    real*8    xx(n),ww(n)
    
    real*8    func
    external    func
    call gauleg(x1,x2,xx,ww,n)

    result=0.
    do i=1,n,1
    result=result+ww(i)*func(xx(i))
    end do
    return
end subroutine qgaus

SUBROUTINE gauleg(x1,x2,x,w,n)
    implicit none
    INTEGER n
    REAL*8 x1,x2,x(n),w(n)
    real*8 EPS
    PARAMETER (EPS=3.d-14) 
    INTEGER i,j,m
    real*8 p1,p2,p3,pp,xl,xm,z,z1
    m=(n+1)/2    
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
  do  i=1,m
    z=cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
110    continue
    p1=1.d0
    p2=0.d0
    do  j=1,n
    p3=p2
    p2=p1
    p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
    enddo
    pp=n*(z*p1-p2)/(z*z-1.d0)
    z1=z
    z=z1-p1/pp                      
    if(abs(z-z1).gt.EPS)then
!    write(*,*)’goto 110′
    goto 110
    end if
    	x(i)=xm-xl*z                 
     	x(n+1-i)=xm+xl*z                 
     	w(i)=2.d0*xl/((1.d0-z*z)*pp*pp) 
    	w(n+1-i)=w(i) 
    enddo
    return
END SUBROUTINE gauleg

FUNCTION gammp(a,x)
	      REAL*8 a,gammp,x
	!CU    USES gcf,gser
	      REAL*8 gammcf,gamser,gln
	      if(x.lt.0..or.a.le.0.)then
	      write(*,*)'bad arguments in gammp'
	      stop
	      end if
	      if(x.lt.a+1.)then
	        call gser(gamser,a,x,gln)
	        gammp=gamser
	      else
	        call gcf(gammcf,a,x,gln)
	        gammp=1.-gammcf
	      endif
	      return
END FUNCTION gammp
!C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.
!C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.


      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
!CU    USES gammln
      INTEGER n

      	REAL*8 ap,del,sum,gammln
!		external	gammln

      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)then
        write(*,*)'x < 0 in gser'
        end if
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      write(*,*)'a too large, ITMAX too small in gser'
      stop
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END SUBROUTINE gser
!C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.

      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
!CU    USES gammln
      INTEGER i
      	REAL*8 an,b,c,d,del,h,gammln
!		external	gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      write(*,*)'a too large, ITMAX too small in gcf',x,a
      
	  stop
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END SUBROUTINE gcf
	
!	C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.

!	c================================================
!	c Bessel functions 
!	c================================================
FUNCTION bessk(n,x)
		implicit none
	      INTEGER n
	      REAL*8 bessk,x
!		USES bessk0,bessk1
	      INTEGER j
	      	REAL*8 bk,bkm,bkp,tox
			real*8		bessk0,bessk1
			external	bessk0,bessk1	
	      if (n.lt.2)then
	      write(*,*)'bad argument n in bessk'
	      stop
	      end if
	      tox=2.0/x
	      bkm=bessk0(x)
	      bk=bessk1(x)

!do 11 j=1,n-1
!	        bkp=bkm+j*tox*bk
!	        bkm=bk
!	        bk=bkp
!11    continue

do  j=1,n-1
	        bkp=bkm+j*tox*bk
	        bkm=bk
	        bk=bkp
end do

	      bessk=bk
	      return
END FUNCTION bessk

FUNCTION bessk0(x)
		implicit none
		REAL*8 bessk0,x
!		   USES bessi0
	      REAL*8,external	::	 bessi0
	      real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
	      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
	      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,&
			0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
	      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,&
			-0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
	      if (x.le.2.0) then
	        y=x*x/4.0
	        bessk0=(-log(x/2.0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
	      else
	        y=(2.0/x)
	        bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))))
	      endif
	      return
END FUNCTION bessk0

FUNCTION bessk1(x)
	implicit none
	REAL*8 bessk1,x
!	USES bessi1
    real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
    SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
    DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,&
		-0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
    DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,&
  		0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/

	    REAL*8		bessi1
		external	bessi1
		
	if (x.le.2.0) then
	    y=x*x/4.0
		bessk1=(log(x/2.0)*bessi1(x))+(1.0/x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
	else
	    y=2.0/x
	    bessk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))))
	endif
	return
END	FUNCTION bessk1

FUNCTION bessi0(x)
	implicit none
	REAL*8 bessi0,x
	REAL*8 ax
	real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
	SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
	DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,&
		1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
	DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,&
		0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,&
				-0.1647633d-1,0.392377d-2/
	if (abs(x).lt.3.75) then
		y=(x/3.75)**2
	    bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
	else
	    ax=abs(x)
	    y=3.75/ax
	    bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
	endif
	      return
END FUNCTION bessi0

FUNCTION bessi1(x)
 REAL*8 bessi1,x
 REAL*8 ax
 real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
 SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
 DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,&
	0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
 DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,&
	-0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,-0.2895312d-1,&
		0.1787654d-1,-0.420059d-2/
if (abs(x).lt.3.75) then
   y=(x/3.75)**2
   bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
 else
   ax=abs(x)
   y=3.75/ax
   bessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
   if(x.lt.0.)bessi1=-bessi1
endif

 return
END FUNCTION bessi1



FUNCTION gammln(xx)
		implicit none
      	REAL*8 gammln,xx
      	INTEGER j
      	real*8 ser,stp,tmp,x,y,cof(6)
!      SAVE cof,stp
!      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
!	24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
!	-.5395239384953d-5,2.5066282746310005d0/
	
	cof(1)=76.18009172947146d0
	cof(2)=-86.50532032941677d0
	cof(3)=24.01409824083091d0
	cof(4)=-1.231739572450155d0
	cof(5)=0.1208650973866179d-2
	cof(6)=-0.5395239384953d-5
	
	stp=2.5066282746310005d0
	
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
	end do
      gammln=tmp+log(stp*ser/x)

      return
END FUNCTION

SUBROUTINE rk4(y,dydx,n,x,h,yout)
implicit none
INTEGER n,NMAX
REAL*8 h,x,dydx(n),y(n),yout(n)
EXTERNAL derivs
PARAMETER (NMAX=3)
INTEGER i
REAL*8 h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
hh=h*0.5
h6=h/6.
xh=x+hh

!write(*,*)x,y,dydx,h,n
!stop

do i=1,n
yt(i)=y(i)+hh*dydx(i)
end do
!write(*,*)'xhyt',xh,yt,dyt
!stop
call derivs(xh,yt,dyt)
do i=1,n
  yt(i)=y(i)+hh*dyt(i)
end do
	call derivs(xh,yt,dym)
do i=1,n
  yt(i)=y(i)+h*dym(i)
  dym(i)=dyt(i)+dym(i)
end do
call derivs(x+h,yt,dyt)
do  i=1,n
  yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
end do

!write(*,*)'rk4:y',y,x
!write(*,*)'rk4:dydx',dydx
!write(*,*)'rk4:yout',yout
!write(*,*)'yout',yout,x
!stop

	return
END


!	radiation process
 function lambda_ie(r,sigma_in,W_i_in,W_e_in)
	use constant
	use vars_bh
	use vars_adaf
 implicit none
real*8				::		r,sigma_in,W_i_in,W_e_in
 real*8 			::		radheat,lambda_ie,lambda_ie_small!ne,ni,
 real*8,external	::	 	bessk,bessk0,bessk1,height
 real*8 			::		tempi,tempe,tempie,temp
	
	T_i=W_i_in*mu_i*m_p/k_bol/sigma_in
	T_e=W_e_in*mu_e*m_p/k_bol/sigma_in

	theta_i=k_bol*T_i/m_p/c/c
	theta_e=k_bol*T_e/m_e/c/c

	n_i=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/mu_i/m_p
	n_e=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/mu_e/m_p

 if((theta_e.ge.1.0d-2).and.(theta_i.ge.1.0d-2))then
 lambda_ie_small=5.61d-32*n_e*n_i*(T_i-T_e)&
        /(bessk(2,1.0/theta_e)*bessk(2,1.0/theta_i))&
        *((2.0*(theta_e+theta_i)**2+1.0)/(theta_e+theta_i)&
        *bessk1((theta_e+theta_i)/theta_e/theta_i)&
        +2.0*bessk0((theta_e+theta_i)/theta_e/theta_i))     

lambda_ie=lambda_ie_small*sqrt(pi)*height(r,sigma_in,W_i_in,W_e_in)
!write(*,*)lambda_ie_small,'1'
 return
 end if
 
 temp=theta_e*theta_i/(theta_e+theta_i) 
      
 if((theta_e.ge.1.0d-2).and.(theta_i.lt.1.0d-2))then
 lambda_ie_small=5.61d-32*n_e*n_i*(T_i-T_e)&
	/bessk(2,1.0/theta_e)*dexp(-1.0/theta_e)&
	*dsqrt(theta_e/(theta_e+theta_i))&
   	*((2.0*(theta_e+theta_i)**2+1)/(theta_e+theta_i)*&
	(1.0+3.0/8.0*temp-15.0/128.0*temp**2+15.0*21.0/6.0/8.0**3*temp**3) &
	+2.0*(1.0-1.0/8.0*temp+9.0/128.0*temp**2&
      -9.0*25.0/6.0/8.0**3*temp**3))&
 	/(1.0+15.0/8.0*theta_i+15.0*7.0/128.0*theta_i**2&
   -15.0*7.0*9.0/6.0/8.0**3*theta_i**3)        
lambda_ie=lambda_ie_small*sqrt(pi)*height(r,sigma_in,W_i_in,W_e_in)  

!write(*,*)lambda_ie_small,'2'
return
 end if 
 
 if((theta_e.lt.1.0d-2).and.(theta_i.ge.1.0d-2))then
 lambda_ie_small=5.61d-32*n_e*n_i*(T_i-T_e)&
	/bessk(2,1.0/theta_i)*dexp(-1.0/theta_i)&
	*dsqrt(theta_i/(theta_e+theta_i))&
	*((2.0*(theta_e+theta_i)**2+1)/(theta_e+theta_i)*&
	(1.0+3.0/8.0*temp-15.0/128.0*temp**2+15.0*21.0/6.0/8.0**3*temp**3) &
	+2.0*(1.0-1.0/8.0*temp+9.0/128.0*temp**2&
    -9.0*25.0/6.0/8.0**3*temp**3)) &
 	/(1.0+15.0/8.0*theta_e+15.0*7.0/128.0*theta_e**2&
   	-15.0*7.0*9.0/6.0/8.0**3*theta_e**3)           
lambda_ie=lambda_ie_small*sqrt(pi)*height(r,sigma_in,W_i_in,W_e_in)

!write(*,*)lambda_ie_small,'3'


 return
 end if          
 
 if((theta_e.lt.1.0d-2).and.(theta_i.lt.1.0d-2))then
	 lambda_ie_small=5.61d-32*n_e*n_i*(T_i-T_e)&
	   *dsqrt(2.0/pi)/dsqrt(theta_e+theta_i)&
	   *((2.0*(theta_e+theta_i)**2+1)/(theta_e+theta_i)*&
	   (1.0+3.0/8.0*temp-15.0/128.0*temp**2+15.0*21.0/6.0/8.0**3*temp**3) &
	   +2.0*(1.0-1.0/8.0*temp+9.0/128.0*temp**2&
	         -9.0*25.0/6.0/8.0**3*temp**3))&
	    /(1.0+15.0/8.0*theta_i+15.0*7.0/128.0*theta_i**2&
	      -15.0*7.0*9.0/6.0/8.0**3*theta_i**3)  &
	    /(1.0+15.0/8.0*theta_e+15.0*7.0/128.0*theta_e**2&
	      -15.0*7.0*9.0/6.0/8.0**3*theta_e**3) 
	lambda_ie=lambda_ie_small*sqrt(pi)*height(r,sigma_in,W_i_in,W_e_in)
!	           write(*,*)lambda_ie_small,'4'
			
return
 end if    


 end function lambda_ie

function ki_mu(mu,r,sigma_in,W_i_in,W_e_in)
use constant
use vars_bh
use vars_adaf
	implicit none
	real*8				::		r,W_i_in,W_e_in,sigma_in
	real*8				::		mu,ki_mu
!	real*8				::		W_i,W_e,sigma
	real*8				::		ki,ki_brem,ki_syn
!	real*8				::		T_i,T_e,theta_i,theta_e
!	real*8				::		n_i,n_e
	real*8				::		B,vars_in_ifun,k2e
	real*8				::		f_theta,g_factor
	real*8				::		q_brem,q_ei,q_ee
	real*8,external		::		syn_i,height,bessk

!W_i_in=W_i
!W_e_in=W_e
!sigma_in=sigma

T_i=W_i_in*mu_i*m_p/k_bol/sigma_in
T_e=W_e_in*mu_e*m_p/k_bol/sigma_in

theta_i=k_bol*T_i/m_p/c/c
theta_e=k_bol*T_e/m_e/c/c

n_e=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_e
n_i=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_i

B=sqrt(8.*pi*(1.-beta)*(W_i_in+W_e_in)/beta/sqrt(2.*pi)/height(r,sigma_in,W_i_in,W_e_in))

vars_in_ifun=4.*pi*m_e*c*mu/3./e/B/theta_e**2.
k2e=bessk(2,1./theta_e)
if (k2e<1.d-150)then
	write(*,*)"syn_ki","bessk(2,1.0/thetae).lt.1.0d-150"
	k2e=1.d-150
end if

ki_syn=4.43d-30*4.*pi*n_e*mu*syn_i(vars_in_ifun)/k2e

!	bremsstraglung
!	f_factor, q_ee
	if(theta_e>1.)then
	f_theta=(9.*theta_e/2./pi)*(log(0.48+1.123*theta_e)+1.5)+2.3*theta_e*(log(1.123*theta_e)+1.28)
		else
	f_theta=4.*(2.*theta_e/pi**3.)**0.5*(1.+1.78*theta_e**1.34)+&
			1.73*theta_e**1.5*(1+1.1*theta_e+theta_e**2.-1.25*theta_e**2.5)

	end if

!	Gaunt factor
	if (k_bol*T_e/h/mu<1.)then
		g_factor=h*(3.*k_bol*T_e/pi/h/mu)**0.5/k_bol/T_e
	else
		g_factor=h*sqrt(3.)*log(4.*k_bol*T_e/0.6979/h/mu)/k_bol/T_e/pi
	end if
	q_brem=1.48d-22*n_e**2.*f_theta
	ki_brem=q_brem*g_factor*exp(-h*mu/k_bol/T_e)	
	ki_mu=ki_brem+ki_syn
	return
end function ki_mu





function syn_i(x)
implicit none
real*8		::		syn_i,x
syn_i=4.0505*(1.+0.4/x**0.25+0.5316/x**0.5)*exp(-1.8899*x**(1./3.))/x**(1./6.)
return
end function syn_i

function yita_compton(mu,r,sigma_in,W_i_in,W_e_in)
use constant
use vars_bh
use vars_adaf
implicit none
real*8		yita_compton,mu
real*8		AA,j_m,s,yita_max,tau_es
real*8		sigma_in,W_i_in,W_e_in,r
real*8,external		::		syn_i,height,gammp

T_i=W_i_in*mu_i*m_p/k_bol/sigma_in
T_e=W_e_in*mu_e*m_p/k_bol/sigma_in

theta_i=k_bol*T_i/m_p/c/c
theta_e=k_bol*T_e/m_e/c/c

n_e=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_e
n_i=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_i

tau_es=sqrt(2.*pi)*n_e*sigma_t*height(r,sigma_in,W_i_in,W_e_in)
s=tau_es+tau_es**2.

AA=1.+4.*theta_e+16.*theta_e**2.
yita_max=3.*k_bol*T_e/h/mu

if(yita_max>1.)then
	j_m=log(yita_max)/log(AA)
	yita_compton=exp(s*(AA-1.))*(1.-gammp(j_m+1.,AA*s))+yita_max*gammp(j_m+1.,s)
	yita_compton=max(1.0d0,yita_compton)
	else
		j_m=0.0
		yita_compton=1.
end if

return
end function yita_compton

function integ_for_q(mu,r,sigma_in,W_i_in,W_e_in)
	use constant
	use vars_bh
	use vars_adaf
	use vars_clouds
	implicit none
	real*8	::	integ_for_q,mu
	real*8	::	black_body
	real*8	::	kappa_mu,tau_mu_star,f_mu
	real*8	::	sigma_in,W_i_in,W_e_in,r
	real*8	::	temp_nu,temp_bnu,B_mu_Tcl,f_c
	real*8,external		::		ki_mu,height
	real*8,external		::		yita_compton

	T_i=W_i_in*mu_i*m_p/k_bol/sigma_in
	T_e=W_e_in*mu_e*m_p/k_bol/sigma_in

	theta_i=k_bol*T_i/m_p/c/c
	theta_e=k_bol*T_e/m_e/c/c


	n_e=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_e
	n_i=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_i
 
	temp_nu=h*mu/k_bol/T_e
	temp_bnu=2.0*h/c**2      

	if(temp_nu.lt.100.0d0)then
	black_body=temp_bnu*mu**3/(dexp(temp_nu)-1.0d0)
	else
	black_body=temp_bnu*mu**3/(dexp(temp_nu/3.0)**3.0d0-1.0d0)
	end if
    
    if(black_body.gt.1.d-100)then
	   kappa_mu=ki_mu(mu,r,sigma_in,W_i_in,W_e_in)/4./pi/black_body
	   tau_mu_star=(pi**0.5/2.)*kappa_mu*height(r,sigma_in,W_i_in,W_e_in)
	
    	if(tau_mu_star.gt.1.0d-6)then
			f_mu=2.*pi*black_body*(1.-exp(-2.*sqrt(3.)*tau_mu_star))/sqrt(3.)
    	else
		f_mu=2.0*pi/dsqrt(3.0d0)*black_body*(2.0*dsqrt(3.0d0)*tau_mu_star&
				-6.0*tau_mu_star**2.+4.0*dsqrt(3.0d0)*tau_mu_star**3.)
    	end if 
    else 	
	f_mu=dsqrt(pi)*height(r,sigma_in,W_i_in,W_e_in)*ki_mu(mu,r,sigma_in,W_i_in,W_e_in)
    end if

	f_c = mass_ratio*sigma_in*R_cl**2./M_cl
	
!	B_mu_Tcl = 2.*h*mu**3./c/c/(exp(h*mu/k_bol/T_cl)-1.)
	temp_nu=h*mu/k_bol/T_cl
	temp_bnu=2.0*h/c**2      
	if(temp_nu.lt.100.0d0)then
	B_mu_Tcl=temp_bnu*mu**3/(dexp(temp_nu)-1.0d0)
	else
	B_mu_Tcl=temp_bnu*mu**3/(dexp(temp_nu/3.0)**3.0d0-1.0d0)
	end if
	
	
	if(r/r_sch<8.)then
		f_c = 0.
	end if	
	
	
	if(tao_cl<1.)then
		f_c = 0.
	end if
	
	
	f_mu = f_mu !+ f_c*B_mu_Tcl

    integ_for_q=yita_compton(mu,r,sigma_in,W_i_in,W_e_in)*f_mu*2.
	return
end function integ_for_q



function q_rad(r,sigma_in,W_i_in,W_e_in)
implicit none
	real*8	::	q_rad
	real*8	::	x1,x2,result,r,mu_temp
	integer,parameter	::		N_gauss=400
	integer		::	i
	real*8	::	xx(N_gauss),ww(N_gauss)
	real*8	::	sigma_in,W_i_in,W_e_in
	
	real*8,external		::		integ_for_q_test
	real*8,external		::		integ_for_q,test_input
	x1=20.7233
	x2=50.6569
	
	call gauleg(x1,x2,xx,ww,N_gauss)
	q_rad=0.

    do i=1,n_gauss,1
		q_rad=q_rad+ww(i)*exp(xx(i))*integ_for_q(exp(xx(i)),r,sigma_in,W_i_in,W_e_in)
    end do
	return
end function q_rad


function flux_nu(mu,r,sigma_in,W_i_in,W_e_in)
	use constant
	use vars_bh
	use vars_adaf
	use vars_clouds
	implicit none
	real*8	::	flux_nu
	real*8	::	integ_for_q,mu
	real*8	::	black_body
	real*8	::	kappa_mu,tau_mu_star,f_mu
	real*8	::	sigma_in,W_i_in,W_e_in,r
	real*8	::	temp_nu,temp_bnu,f_c,B_mu_Tcl
	real*8,external		::		ki_mu,height
	real*8,external		::		yita_compton

	T_i=W_i_in*mu_i*m_p/k_bol/sigma_in
	T_e=W_e_in*mu_e*m_p/k_bol/sigma_in

	theta_i=k_bol*T_i/m_p/c/c
	theta_e=k_bol*T_e/m_e/c/c


n_e=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_e
n_i=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_i
 
temp_nu=h*mu/k_bol/T_e
temp_bnu=2.0*h/c**2      

if(temp_nu.lt.100.0d0)then
black_body=temp_bnu*mu**3/(dexp(temp_nu)-1.0d0)
else
black_body=temp_bnu*mu**3/(dexp(temp_nu/3.0)**3.0d0-1.0d0)
end if
    
    if(black_body.gt.1.d-100)then
	   kappa_mu=ki_mu(mu,r,sigma_in,W_i_in,W_e_in)/4./pi/black_body
	   tau_mu_star=(pi**0.5/2.)*kappa_mu*height(r,sigma_in,W_i_in,W_e_in)
	
    	if(tau_mu_star.gt.1.0d-6)then
			flux_nu=2.*pi*black_body*(1.-exp(-2.*sqrt(3.)*tau_mu_star))/sqrt(3.)
    	else
			flux_nu=2.0*pi/dsqrt(3.0d0)*black_body*(2.0*dsqrt(3.0d0)*tau_mu_star&
				-6.0*tau_mu_star**2.+4.0*dsqrt(3.0d0)*tau_mu_star**3.)
    	end if 
    else 	
		flux_nu=dsqrt(pi)*height(r,sigma_in,W_i_in,W_e_in)&
				*ki_mu(mu,r,sigma_in,W_i_in,W_e_in)
    end if

  
	f_c = mass_ratio*sigma_in*R_cl**2./M_cl
	
!	B_mu_Tcl = 2.*h*mu**3./c/c/(exp(h*mu/k_bol/T_cl)-1.)
	temp_nu=h*mu/k_bol/T_cl
	temp_bnu=2.0*h/c**2      
	if(temp_nu.lt.100.0d0)then
	B_mu_Tcl=temp_bnu*mu**3/(dexp(temp_nu)-1.0d0)
	else
	B_mu_Tcl=temp_bnu*mu**3/(dexp(temp_nu/3.0)**3.0d0-1.0d0)
	end if
	
	
	if(tao_cl<1.)then
		f_c = 0.
	end if
	
	if(r/r_sch<rr_tidal)then
		f_c = 0.
	end if	
	
	
	flux_nu = flux_nu !+ f_c*B_mu_Tcl

	return
end function flux_nu


function flux_nu_syn(mu,r,sigma_in,W_i_in,W_e_in)
	use constant
	use vars_bh
	use vars_adaf
	use vars_clouds
	implicit none
	real*8	::	flux_nu_syn
	real*8	::	integ_for_q,mu
	real*8	::	black_body
	real*8	::	kappa_mu,tau_mu_star,f_mu
	real*8	::	sigma_in,W_i_in,W_e_in,r
	real*8	::	temp_nu,temp_bnu,f_c,B_mu_Tcl
	real*8,external		::		ki_mu_syn,height
	real*8,external		::		yita_compton

	T_i=W_i_in*mu_i*m_p/k_bol/sigma_in
	T_e=W_e_in*mu_e*m_p/k_bol/sigma_in

	theta_i=k_bol*T_i/m_p/c/c
	theta_e=k_bol*T_e/m_e/c/c


n_e=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_e
n_i=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_i
 
temp_nu=h*mu/k_bol/T_e
temp_bnu=2.0*h/c**2      

if(temp_nu.lt.100.0d0)then
black_body=temp_bnu*mu**3/(dexp(temp_nu)-1.0d0)
else
black_body=temp_bnu*mu**3/(dexp(temp_nu/3.0)**3.0d0-1.0d0)
end if
    
    if(black_body.gt.1.d-100)then
	   kappa_mu=ki_mu_syn(mu,r,sigma_in,W_i_in,W_e_in)/4./pi/black_body
	   tau_mu_star=(pi**0.5/2.)*kappa_mu*height(r,sigma_in,W_i_in,W_e_in)
	
    	if(tau_mu_star.gt.1.0d-6)then
			flux_nu_syn=2.*pi*black_body*(1.-exp(-2.*sqrt(3.)*tau_mu_star))/sqrt(3.)
    	else
			flux_nu_syn=2.0*pi/dsqrt(3.0d0)*black_body*(2.0*dsqrt(3.0d0)*tau_mu_star&
				-6.0*tau_mu_star**2.+4.0*dsqrt(3.0d0)*tau_mu_star**3.)
    	end if 
    else 	
		flux_nu_syn=dsqrt(pi)*height(r,sigma_in,W_i_in,W_e_in)&
				*ki_mu_syn(mu,r,sigma_in,W_i_in,W_e_in)
    end if

  
	f_c = mass_ratio*sigma_in*R_cl**2./M_cl
	
!	B_mu_Tcl = 2.*h*mu**3./c/c/(exp(h*mu/k_bol/T_cl)-1.)
	temp_nu=h*mu/k_bol/T_cl
	temp_bnu=2.0*h/c**2      
	if(temp_nu.lt.100.0d0)then
	B_mu_Tcl=temp_bnu*mu**3/(dexp(temp_nu)-1.0d0)
	else
	B_mu_Tcl=temp_bnu*mu**3/(dexp(temp_nu/3.0)**3.0d0-1.0d0)
	end if
	
	flux_nu_syn = flux_nu_syn! + f_c*B_mu_Tcl

	return
end function flux_nu_syn

function ki_mu_syn(mu,r,sigma_in,W_i_in,W_e_in)
use constant
use vars_bh
use vars_adaf
	implicit none
	real*8				::		r,W_i_in,W_e_in,sigma_in
	real*8				::		mu,ki_mu_syn
!	real*8				::		W_i,W_e,sigma
	real*8				::		ki,ki_brem,ki_syn
!	real*8				::		T_i,T_e,theta_i,theta_e
!	real*8				::		n_i,n_e
	real*8				::		B,vars_in_ifun,k2e
	real*8				::		f_theta,g_factor
	real*8				::		q_brem,q_ei,q_ee
	real*8,external		::		syn_i,height,bessk

!W_i_in=W_i
!W_e_in=W_e
!sigma_in=sigma

T_i=W_i_in*mu_i*m_p/k_bol/sigma_in
T_e=W_e_in*mu_e*m_p/k_bol/sigma_in

theta_i=k_bol*T_i/m_p/c/c
theta_e=k_bol*T_e/m_e/c/c

n_e=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_e
n_i=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_i

B=sqrt(8.*pi*(1.-beta)*(W_i_in+W_e_in)/beta/sqrt(2.*pi)/height(r,sigma_in,W_i_in,W_e_in))

vars_in_ifun=4.*pi*m_e*c*mu/3./e/B/theta_e**2.
k2e=bessk(2,1./theta_e)
if (k2e<1.d-150)then
	write(*,*)"syn_ki","bessk(2,1.0/thetae).lt.1.0d-150"
	k2e=1.d-150
end if

ki_syn=4.43d-30*4.*pi*n_e*mu*syn_i(vars_in_ifun)/k2e

!	bremsstraglung
!	f_factor, q_ee
	if(theta_e>1.)then
	f_theta=(9.*theta_e/2./pi)*(log(0.48+1.123*theta_e)+1.5)+2.3*theta_e*(log(1.123*theta_e)+1.28)
		else
	f_theta=4.*(2.*theta_e/pi**3.)**0.5*(1.+1.78*theta_e**1.34)+&
			1.73*theta_e**1.5*(1+1.1*theta_e+theta_e**2.-1.25*theta_e**2.5)

	end if

!	Gaunt factor
	if (k_bol*T_e/h/mu<1.)then
		g_factor=h*(3.*k_bol*T_e/pi/h/mu)**0.5/k_bol/T_e
	else
		g_factor=h*sqrt(3.)*log(4.*k_bol*T_e/0.6979/h/mu)/k_bol/T_e/pi
	end if
	q_brem=1.48d-22*n_e**2.*f_theta
	ki_brem=q_brem*g_factor*exp(-h*mu/k_bol/T_e)	
	ki_mu_syn=ki_syn!+ki_brem
	return
end function ki_mu_syn


function flux_nu_bre(mu,r,sigma_in,W_i_in,W_e_in)
	use constant
	use vars_bh
	use vars_adaf
	use vars_clouds
	implicit none
	real*8	::	flux_nu_bre
	real*8	::	integ_for_q,mu
	real*8	::	black_body
	real*8	::	kappa_mu,tau_mu_star,f_mu
	real*8	::	sigma_in,W_i_in,W_e_in,r
	real*8	::	temp_nu,temp_bnu,f_c,B_mu_Tcl
	real*8,external		::		ki_mu_bre,height
	real*8,external		::		yita_compton

	T_i=W_i_in*mu_i*m_p/k_bol/sigma_in
	T_e=W_e_in*mu_e*m_p/k_bol/sigma_in

	theta_i=k_bol*T_i/m_p/c/c
	theta_e=k_bol*T_e/m_e/c/c


n_e=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_e
n_i=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_i
 
temp_nu=h*mu/k_bol/T_e
temp_bnu=2.0*h/c**2      

if(temp_nu.lt.100.0d0)then
black_body=temp_bnu*mu**3/(dexp(temp_nu)-1.0d0)
else
black_body=temp_bnu*mu**3/(dexp(temp_nu/3.0)**3.0d0-1.0d0)
end if
    
    if(black_body.gt.1.d-100)then
	   kappa_mu=ki_mu_bre(mu,r,sigma_in,W_i_in,W_e_in)/4./pi/black_body
	   tau_mu_star=(pi**0.5/2.)*kappa_mu*height(r,sigma_in,W_i_in,W_e_in)
	
    	if(tau_mu_star.gt.1.0d-6)then
			flux_nu_bre=2.*pi*black_body*(1.-exp(-2.*sqrt(3.)*tau_mu_star))/sqrt(3.)
    	else
			flux_nu_bre=2.0*pi/dsqrt(3.0d0)*black_body*(2.0*dsqrt(3.0d0)*tau_mu_star&
				-6.0*tau_mu_star**2.+4.0*dsqrt(3.0d0)*tau_mu_star**3.)
    	end if 
    else 	
		flux_nu_bre=dsqrt(pi)*height(r,sigma_in,W_i_in,W_e_in)&
				*ki_mu_bre(mu,r,sigma_in,W_i_in,W_e_in)
    end if

  
	f_c = mass_ratio*sigma_in*R_cl**2./M_cl
	
!	B_mu_Tcl = 2.*h*mu**3./c/c/(exp(h*mu/k_bol/T_cl)-1.)
	temp_nu=h*mu/k_bol/T_cl
	temp_bnu=2.0*h/c**2      
	if(temp_nu.lt.100.0d0)then
	B_mu_Tcl=temp_bnu*mu**3/(dexp(temp_nu)-1.0d0)
	else
	B_mu_Tcl=temp_bnu*mu**3/(dexp(temp_nu/3.0)**3.0d0-1.0d0)
	end if
	
	flux_nu_bre = flux_nu_bre! + f_c*B_mu_Tcl

	return
end function flux_nu_bre

function ki_mu_bre(mu,r,sigma_in,W_i_in,W_e_in)
use constant
use vars_bh
use vars_adaf
	implicit none
	real*8				::		r,W_i_in,W_e_in,sigma_in
	real*8				::		mu,ki_mu_bre
!	real*8				::		W_i,W_e,sigma
	real*8				::		ki,ki_brem,ki_syn
!	real*8				::		T_i,T_e,theta_i,theta_e
!	real*8				::		n_i,n_e
	real*8				::		B,vars_in_ifun,k2e
	real*8				::		f_theta,g_factor
	real*8				::		q_brem,q_ei,q_ee
	real*8,external		::		syn_i,height,bessk

!W_i_in=W_i
!W_e_in=W_e
!sigma_in=sigma

T_i=W_i_in*mu_i*m_p/k_bol/sigma_in
T_e=W_e_in*mu_e*m_p/k_bol/sigma_in

theta_i=k_bol*T_i/m_p/c/c
theta_e=k_bol*T_e/m_e/c/c

n_e=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_e
n_i=sigma_in/sqrt(2.*pi)/Height(r,sigma_in,W_i_in,W_e_in)/m_p/mu_i

B=sqrt(8.*pi*(1.-beta)*(W_i_in+W_e_in)/beta/sqrt(2.*pi)/height(r,sigma_in,W_i_in,W_e_in))

vars_in_ifun=4.*pi*m_e*c*mu/3./e/B/theta_e**2.
k2e=bessk(2,1./theta_e)
if (k2e<1.d-150)then
	write(*,*)"syn_ki","bessk(2,1.0/thetae).lt.1.0d-150"
	k2e=1.d-150
end if

ki_syn=4.43d-30*4.*pi*n_e*mu*syn_i(vars_in_ifun)/k2e

!	bremsstraglung
!	f_factor, q_ee
	if(theta_e>1.)then
	f_theta=(9.*theta_e/2./pi)*(log(0.48+1.123*theta_e)+1.5)+2.3*theta_e*(log(1.123*theta_e)+1.28)
		else
	f_theta=4.*(2.*theta_e/pi**3.)**0.5*(1.+1.78*theta_e**1.34)+&
			1.73*theta_e**1.5*(1+1.1*theta_e+theta_e**2.-1.25*theta_e**2.5)

	end if

!	Gaunt factor
	if (k_bol*T_e/h/mu<1.)then
		g_factor=h*(3.*k_bol*T_e/pi/h/mu)**0.5/k_bol/T_e
	else
		g_factor=h*sqrt(3.)*log(4.*k_bol*T_e/0.6979/h/mu)/k_bol/T_e/pi
	end if
	q_brem=1.48d-22*n_e**2.*f_theta
	ki_brem=q_brem*g_factor*exp(-h*mu/k_bol/T_e)	
	ki_mu_bre=ki_brem
	return
end function ki_mu_bre

