!> This module contains  a number of mineral and rock physics
  !!utility routines. There are a number of derived types defined
  !! also. fit contains the data for polynomial fitting parameters.
  !! It is used by functions to convert temperature into Vs, Vp,
  !! and density. Type melt contains parameters for use in EOS
  !! for melt. Composition contains information on the bulk
  !! composition of the reference mantle and the potential
  !! temperature. Type mantle contains the information from
  !! from mineral physics on the physical properties of the
  !! reference mantle for a given value of the variable composition.
  !! Finally, unit cell contains the information on the seismic
  !! signature of the reference mantle (subscript sol), the observed
  !! seismic signature (subscript obs), and the effective
  !! values (subscript eff). If the difference between reference and
  !! observed values can be explained by melting, then the
  !! melt fraction and melt fraction dependent contiguity
  !! are calculated. Temperature is the potential (not actual)
  !! temperature for the unit cell. For example, it can be 
  !! calculated from the depression/uplift of the transition zone
  !! boundaries in later modules. See the description of each
  !! function within the function body.
  !! Created by Saswata Hier-Majumder, August, 2012. !<
MODULE microgeodynamics
  USE global
  TYPE fit
     REAL(sp)::p1,p2,p3,p4,p5,p6
  END TYPE fit
  TYPE melt
     REAL(sp)::K0,Kp,rho0,dihedral,aspect
  END TYPE melt
  TYPE composition
     INTEGER::potential_temperature,basalt_fraction
  END TYPE composition
  TYPE mantle
     TYPE(composition)::comp
     REAL(sp),DIMENSION(352)::P,Z,rho,vs,vp,K,G,nu,T,X
  END TYPE mantle
  TYPE unit_cell
     REAL(sp)::vs_sol,vp_sol,K_sol,G_sol,nu,Kl,rhol,rho,pressure
     REAL(sp)::vs_eff,vp_eff,K_eff,G_eff
     REAL(sp)::vs_obs,vp_obs,K_obs,G_obs,amplitude
     REAL(sp)::melt_fraction,contiguity,temperature
  END TYPE unit_cell

CONTAINS

  FUNCTION clear_unit_cell(cell)
    IMPLICIT NONE
    TYPE(unit_cell),INTENT(in)::cell
    TYPE(unit_cell)::clear_unit_cell

    !> As the name suggests, it zeros out the elements
    !!of the unit cell.!<
    clear_unit_cell%vs_sol=0.0_sp
    clear_unit_cell%vp_sol=0.0_sp
    clear_unit_cell%K_sol=0.0_sp
    clear_unit_cell%G_sol=0.0_sp
    clear_unit_cell%nu=0.0_sp
    clear_unit_cell%Kl=0.0_sp
    clear_unit_cell%rhol=0.0_sp
    clear_unit_cell%rho=0.0_sp
    clear_unit_cell%pressure=0.0_sp

    clear_unit_cell%vs_obs=0.0_sp
    clear_unit_cell%vp_obs=0.0_sp
    clear_unit_cell%K_obs=0.0_sp
    clear_unit_cell%G_obs=0.0_sp

    clear_unit_cell%vs_eff=0.0_sp
    clear_unit_cell%vp_eff=0.0_sp
    clear_unit_cell%K_eff=0.0_sp
    clear_unit_cell%G_eff=0.0_sp

    clear_unit_cell%melt_fraction=0.0_sp
    clear_unit_cell%contiguity=0.0_sp
    clear_unit_cell%temperature=0.0_sp

  END FUNCTION clear_unit_cell

  !> This function sets the properties of the melt
  !! See Table 1 of Wimert and Hier-Majumder(2012) for details
  !!The choices are 
  !!1 = peridotite melt, (2273 K) Guillot and Sator (2007)
  !!2 = MORB,(2073 K) Guillot and Sator  (2007)
  !!3 = peridotite melt, Ohtani and Maeda (2001)
  !!4 = MORB, Ohtani and Maeda (2001)
  !!5= peridotite+5%CO2 from Ghosh etal (2007)
  !!Any other value defaults to 1
  !! The optional value of dihedral angle
  !! is the second input. !<
  FUNCTION set_melt(i,theta)
    IMPLICIT NONE
    TYPE(melt)::set_melt,m
    INTEGER,INTENT(in)::i
    REAL(sp),OPTIONAL,INTENT(in)::theta
    REAL(sp)::temp
    IF(PRESENT(theta)) THEN
       temp=theta
    ELSE
       temp=25.0_sp
    END IF

    IF(i==2) THEN
       m%K0=15.5e9_sp
       m%Kp=7.2_sp
       m%rho0=2590.0_sp
       m%dihedral=temp
    ELSEIF (i==3) THEN
       m%K0=24.9e9_sp
       m%Kp=6.4_sp
       m%rho0=2610.0_sp
       m%dihedral=temp
    ELSEIF(i==4) THEN
       m%K0=18.1e9_sp
       m%Kp=5.5_sp
       m%rho0=2590.0_sp
       m%dihedral=temp
    ELSEIF(i==5) THEN !Fig 4
       m%K0=24.9e9_sp
       m%Kp=5.1_sp
       m%rho0=2670.0_sp
       m%dihedral=temp
    ELSE
       m%K0=16.5e9_sp
       m%Kp=7.2_sp
       m%rho0=2610.0_sp
       m%dihedral=temp
    END IF
    set_melt%K0=m%K0
    set_melt%kp=m%Kp
    set_melt%rho0=m%rho0
    set_melt%dihedral= m%dihedral

  END FUNCTION set_melt

  FUNCTION vbw(dihedral,meltfrac)
    IMPLICIT NONE
    REAL(sp),INTENT(in)::dihedral,meltfrac
    REAL(sp)::vbw
    REAL(sp),DIMENSION(3)::bss,pss,bsl,psl
    REAL(sp)::b,p,Agg,Agm

    !> This function returns the contiguity
    !!of a partially molten unit cell as a function
    !!of melt volume fraction and dihedral angle,
    !!using the parametrization of von Bargen and
    !!Waff (1986). Notice that the parameter Agg
    !!shouldn't become zero at zero melt fraction,
    !!as erroneously indicated in their article.
    !! it is fixed by subtracting it from pi to match their
    !! Figure 10. Doesn't work
    !! beyond melt volume fraction fo 0.18. !<
    bss=(/8.16_sp, -7.71e-2_sp, 1.03e-3_sp/)
    pss=(/0.424_sp, 9.95e-4_sp, 8.6645e-6_sp/)
    bsl=(/12.86_sp, -7.85e-2_sp, 1.0043e-3_sp/)
    psl =(/0.43_sp, 8.63e-5_sp, 2.41e-5_sp/)

    b=bss(3)*dihedral**2+bss(2)*dihedral+bss(1) 
    p=pss(3)*dihedral**2+pss(2)*dihedral+pss(1)

    Agg=PI-b*meltfrac**p 
    IF(Agg<0.0_sp) THEN
       !print*,'VBW: error in contiguity calculation'
    END IF

    b=bsl(3)*dihedral**2+bsl(2)*dihedral+bsl(1) 
    p=psl(3)*dihedral**2+psl(2)*dihedral+psl(1) 

    Agm=b*meltfrac**p 


    vbw=2.0_sp*Agg/(2.0_sp*Agg+Agm)

  END FUNCTION vbw

  FUNCTION hmrb(dihedral,meltfrac)
    IMPLICIT NONE
    REAL(sp),INTENT(in)::dihedral,meltfrac
    REAL(sp)::hmrb,temp1,chi1,chi2,Agm,Agg,temp2
    !> This function returns the two dimensional
    !! measurement of contiguity from Hier-Majumder
    !! et al (2006). This is not recommended as it typically
    !! returns contiguity values higher than the 3D models. !<
    ! 
    temp2=dihedral*PI/90.0_sp  
    temp1=SQRT(COS(temp2)**2-SIN(temp2)*COS(temp2)-0.25_sp*PI+temp2)

    chi1=ABS(0.5_sp*PI -2.0_sp*temp2)/temp1

    chi2=ABS(COS(temp2)-SIN(temp2))/temp1
    Agm=chi1*SQRT(meltfrac)
    Agg= 1.0_sp-chi2*SQRT(meltfrac)
    hmrb=Agg/(Agg+Agm)

  END FUNCTION hmrb

  FUNCTION whm(dihedral,meltfrac)
    IMPLICIT NONE
    REAL(sp),INTENT(in)::dihedral,meltfrac
    REAL(sp)::whm
    REAL(sp),DIMENSION(6)::p
    !> This function returns the contiguity as a function of
    !!melt fraction. The dihedral angle, even taken as an input
    !!is currently not used, as the model of Wimert and
    !!Hier-Majumder is valid for a constant dihedral angle of
    !!approximately 30 degreees. Doesn't work well beyond melt
    !!volume fraction of 0.25.!<
    p(1)=-8065.0_sp
    p(2)=6149.0_sp
    p(3)=-1778.0_sp
    p(4)=249.0_sp
    p(5)=-19.77_sp
    p(6)=1.0_sp

    whm=p(1)*meltfrac**5+p(2)*meltfrac**4+p(3)*meltfrac**3&
         &+p(4)*meltfrac**2+p(5)*meltfrac+p(6)
    IF (whm<=0.0_sp) THEN
       ! PRINT*,'whm: error in calculation'
    END IF
  END FUNCTION whm

  SUBROUTINE vel2mod(vs,vp,rho,G,K,nu)
    IMPLICIT NONE
    REAL(sp),INTENT(in)::vs,vp,rho
    REAL(sp),INTENT(out)::G,K,nu
    !> This utility subroutine calculates shear modulus G,
    !!bulk modulus K, and Poisson's ratio nu, from known
    !!shear and P wave velocities. The use of SI units are not
    !!necessary, but strongly recommended, especially if
    !!using the functions related to melting. !<

    G = rho*vs**2
    K = rho*(vp**2-1.33_sp*vs**2)
    nu=(3.0_sp*K-2.0_sp*G)/(6.0_sp*K+2.0_sp*G)
  END SUBROUTINE vel2mod

  SUBROUTINE mod2vel(G,K,rho,vs,vp)
    IMPLICIT NONE
    REAL(sp),INTENT(in)::G,K,rho
    REAL(sp),INTENT(out)::vs,vp
    !> This subroutine calculates Vs and Vp from the moduli.
    !!Here, too, the use of SI units are strongly recommended.!<
    vs=SQRT(G/rho)
    vp=SQRT(K/rho+1.33_sp*G/rho)

  END SUBROUTINE mod2vel

  function elastic_film(aspect,meltfrac,K1,mu1,rho1,K2,rho2)
    use global
    REAL(sp),INTENT(in)::aspect,meltfrac,K1,mu1,K2,rho1,rho2
    REAL(sp),DIMENSION(4)::elastic_film
    real(sp)::rho_av,term1,term2,vs0,vp0,vs,vp,Keff,Geff,mu2,term3,term4
    !> This function calculates Vs/Vs0,Vp/VP0, K/K0, and G/G0
    !! Based on the model of Walsh, 1969. On input, the aspect
    !! ratio is the ratio between the minor and major axis 
    !! of the elliptical melt inclusion.
    !! Subscript 1 is the stronger phase.
    !! On return, the four members in the array
    !! elastic_melt are respectively
    !!Vs/Vs0,Vp/VP0, K/K0, and G/G0. Please only enter SI units. !<

    call mod2vel(mu1,K1,rho1,vs0,vp0)

    mu2=0.0_sp ! Shear modulus of melt
    term3=3.0_sp*aspect*PI*mu1*(3.0_sp*K1+mu1)/(3.0_sp*K1+4.0_sp*mu1)
    term4=3.0_sp*aspect*PI*mu1*(3.0_sp*K1+2.0_sp*mu1)/(3.0_sp*K1+4.0_sp*mu1)
    term1=meltfrac*(K1-K2)*(3.0_sp*K1+4.0_sp*mu2)/&
         &(3.0_sp*K2+4.0_sp*mu2+term3)/K1
    term2=meltfrac*(1+8.0_sp*mu1/(4.0_sp*mu2+term4)+&
         &2.0_sp*(3.0_sp*K2+2.0_sp*mu2+2.0_sp*mu1)/&
         &(3.0_sp*K2+4.0_sp*mu2+term3))*(mu1-mu2)/mu1/5.0_sp
    Keff=1.0_sp/(1.0_sp+term1)
    Geff=1.0_sp/(1.0_sp+term2)
    
    rho_av=(1.0_sp-meltfrac)*rho1+meltfrac*rho2
    call mod2vel(Geff*mu1,Keff*K1,rho_av,vs,vp)
    elastic_film(1)=vs/vs0
    elastic_film(2)=vp/vp0
    elastic_film(3)=Keff
    elastic_film(4)=Geff
  end function elastic_film

  FUNCTION elastic_melt(contiguity,meltfrac,K,G,nu,rho,Kl,rhol)
    REAL(sp),INTENT(in)::contiguity,meltfrac,K,G,nu,Kl,rho,rhol
    REAL(sp),DIMENSION(4)::elastic_melt
    REAL(sp)::vs1,vp1,K1,G1
    REAL(sp)::a1,a2,a3,b1,b2,b3
    REAL(sp)::m,n,h,gg,KboverK,bet,rhobaroverrho,rr0
    !> This function calculates Vs/Vs0,Vp/VP0, K/K0, and G/G0
    !! as a function of melt fraction, contiguity, and
    !! known solid bulk modulus,
    !! solid shear modulus, Poisson's ratio, solid density,
    !! liquid bulk modulus, and
    !! liquid density. 
    !! On return, the four members in the array
    !! elastic_melt are respectively
    !!Vs/Vs0,Vp/VP0, K/K0, and G/G0. Please only enter SI units. !<

    a1=1.8625_sp+0.52594_sp*nu-4.8397_sp*nu**2
    a2=4.5001_sp-6.1551_sp*nu-4.3634_sp*nu**2
    a3=-5.6512_sp+6.9159_sp*nu+29.595_sp*nu**2-58.96_sp*nu**3
    b1=1.6122_sp+0.13527_sp*nu
    b2=4.5869_sp+3.6086_sp*nu
    b3=-7.5395_sp-4.8676_sp*nu-4.3182_sp*nu**2
    !Equation A5
    m=a1*contiguity+a2*(1.0_sp-contiguity)+&
         &a3*((1.0_sp-contiguity)**1.5_sp)*contiguity 
    !Equation A6
    n=b1*contiguity+b2*(1.0_sp-contiguity)+&
         &b3*contiguity*(1.0_sp-contiguity)**2 

    h=1.0_sp-(1.0_sp-contiguity)**m    !Eq A3
    gg=1.0_sp-(1.0_sp-contiguity)**n   ! Eq A4

    !normalized bulk modulus of the skeletal framework eq 27, H-m 2008
    KboverK=(1.0_sp-meltfrac)*h        !Eq A1 divided by k
    G1=(1.0_sp-meltfrac)*gg            !Eq A2 divided by mu

    bet=K/Kl
    gam=G/K
    rr0=rhol/rho

    K1=KboverK+((1.0_sp-KboverK)**2)/(1.0_sp-meltfrac-KboverK+meltfrac*bet)
    rhobaroverrho=1.0_sp-meltfrac+rr0*meltfrac

    vs1=SQRT(G1)/SQRT(rhobaroverrho)

    vp1=SQRT(K1+4.0_sp*gam*G1/3.0_sp)&
         &/SQRT(1.0_sp+4.0_sp*gam/3.0_sp)/SQRT(rhobaroverrho)

    elastic_melt(1)=vs1
    elastic_melt(2)=vp1
    elastic_melt(3)=K1
    elastic_melt(4)=G1
  END FUNCTION elastic_melt

  SUBROUTINE vinet(K0,rho,Kp,rho0,P,K)
    IMPLICIT NONE
    REAL(sp),INTENT(in)::K0,rho,Kp,rho0
    REAL(sp),INTENT(out)::P,K
    REAL(sp)::z,term1,term2,term3,dpdz
    !> Returns the bulk modulus, K, and pressure,P,
    !!from known values of density, rho, surface
    !!bulk modulus, K0, and the pressure derivative
    !!of the bulk modulus, Kp. It uses the Vinet EOS. !<
    z=rho/rho0
    P=3*K0*z**(2.0_sp/3.0_sp)*(1.0_sp-z**(-1.0_sp/3.0_sp))&
         &*EXP(1.5_sp*(Kp-1.0_sp)*(1.0_sp-z**(-1.0_sp/3.0_sp)))
    term1=3.0_sp*K0*EXP(1.5_sp*(Kp-1.0_sp)*(1.0_sp-z**(-1.0_sp/3.0_sp)))
    term2=z**(-1.0_sp/3.0_sp)*(2.0_sp-z**(-1.0_sp/3.0_sp))/3.0_sp
    term3=0.5_sp*z**(-2.0_sp/3.0_sp)*(1.0_sp-z**(-1.0_sp/3.0_sp))*(Kp-1.0_sp)
    dpdz=term1*(term2+term3)
    K=z*dpdz
  END SUBROUTINE vinet

  SUBROUTINE bm3(K0,rho,Kp,rho0,P,K)
    IMPLICIT NONE
    REAL(sp),INTENT(in)::K0,rho,Kp,rho0
    REAL(sp),INTENT(out)::P,K
    REAL(sp)::z,term1,term2,term3,dpdz
    !> Returns the bulk modulus, K, and pressure,P,
    !!from known values of density, rho, surface
    !!bulk modulus, K0, and the pressure derivative
    !!of the bulk modulus, Kp. It uses the third order Birch
    !!Murnaghan EOS. !<
    z=rho/rho0
    P=1.5_sp*K0*((rho/rho0)**(7.0_sp/3.0_sp)&
         &-(rho/rho0)**(5.0_sp/3.0_sp))*&
         &(1.0_sp+0.75_sp*(Kp-4.0_sp)*&
         &((rho/rho0)**(2.0_sp/3.0_sp)-1.0_sp))
    term1=0.75_sp*K0*(Kp-4.0_sp)*(z**2-z**(4.0_sp/3.0_sp))
    term2=1.0_sp+0.75_sp*(Kp-4.0_sp)*(z**(2.0_sp/3.0_sp)-1.0_sp)
    term3=0.5_sp*K0*(7.0_sp*z**(4.0_sp/3.0_sp)-5.0_sp*z**(2.0_sp/3.0_sp))
    dpdz=term1+term2*term3
    K=z*dpdz
  END SUBROUTINE bm3

  FUNCTION load_composition(comp,fname)
    IMPLICIT NONE
    TYPE(composition),INTENT(in)::comp
    CHARACTER(len=*),INTENT(in)::fname
    TYPE(mantle)::load_composition
    INTEGER::ii,filesize

    REAL(sp)::P,Z,T,rho,X,vs,vp,S
    load_composition%comp=comp
    filesize=line_numbers(fname,10000)
    OPEN(10,file=fname,form='formatted')

    DO ii=1,filesize
       READ(10,'(8(f10.2,1x))'),P,Z,T,rho,X,vs,vp,S
       load_composition%P(ii)=P*1.0e9_sp! Convert to Pa
       load_composition%Z(ii)=Z*1.0e3_sp! Convert to m
       load_composition%rho(ii)=rho*1.0e3_sp! Convert to kg/m^3
       load_composition%vs(ii)=vs*1.0e3_sp! Convert to m/s
       load_composition%vp(ii)=vp*1.0e3_sp! Convert to m/s
       load_composition%T(ii)=T           !potential temperature
       load_composition%X(ii)=X           !basalt fraction
       ! Calculate the bulk and shear modulus
       CALL vel2mod(vs*1.0e3_sp,vp*1.0e3_sp,rho*1.0e3_sp,&
            &load_composition%G(ii),&
            &load_composition%K(ii),load_composition%nu(ii))
    END DO
    CLOSE(10)
  END FUNCTION load_composition

  FUNCTION bisection(x1,x2,eps,cell,magma)
    IMPLICIT NONE
    REAL(sp),INTENT(in)::x1,x2,eps
    TYPE(unit_cell),INTENT(in)::cell
    TYPE(melt),INTENT(in)::magma
    REAL(sp)::bisection,dx,f,fmid,xmid
    REAL(sp),DIMENSION(2)::temp
    INTEGER::ii,maxit,jj
    maxit=200
    temp=funcd(x2,cell,magma)
    fmid=temp(1)
    temp=funcd(x1,cell,magma)
    f=temp(1)
    IF(f*fmid>0) THEN
       PRINT*,'root not bracketted in bisection method'
       bisection=x1
       RETURN
    END IF
    IF (f<0.0_sp) THEN
       bisection=x1
       dx=x2-x1
    ELSE
       bisection=x2
       dx=x1-x2
    END IF
    DO jj=1,maxit
       dx=0.5_sp*dx
       xmid=bisection+dx
       temp=funcd(xmid,cell,magma)
       fmid=temp(1)
       IF (fmid<0.0_sp) THEN
          bisection=xmid
       END IF

       IF(ABS(dx)<eps.OR.fmid==0.0_sp) RETURN

    END DO

  END FUNCTION bisection



  FUNCTION nonlinear_solve(x1,x2,eps,cell,magma)
    IMPLICIT NONE
    REAL(sp),INTENT(in)::x1,x2,eps
    TYPE(unit_cell),INTENT(in)::cell
    TYPE(melt),INTENT(in)::magma
    REAL(sp)::sol,nonlinear_solve,xlo,xhi,dx,dxold,temp
    REAL(sp),DIMENSION(2)::flo,fhi,f
    INTEGER::ii,maxit,jj

    !> This function finds the root of a known, nonlinear function funcd
    !!between values x1 and x2, within precision 1.0e-3
    !! This function uses a combination of bisection 
    !!method with Newton Raphson method.
    !<
    !! First get the values of the function and its derivative
    flo=funcd(x1,cell,magma)
    fhi=funcd(x2,cell,magma)
!!$

    !! Does the function cross zero within the provided bounds?
    IF((flo(1)<0.0_sp.AND.fhi(1)<0.0_sp).&
         &OR.(flo(1)>0.0_sp.AND.fhi(1)>0.0_sp)) THEN

       IF((flo(1)<0.0_sp.AND.fhi(1)<0.0_sp)) THEN
          nonlinear_solve=0.0_sp
          !print*,'observed velocity faster than prediction, no melt'
       ELSEIF((flo(1)>0.0_sp.AND.fhi(1)>0.0_sp)) THEN
          nonlinear_solve=0.26
          PRINT*,'congratulations! you found a magma ocean!'
       END IF

       RETURN    
    ELSEIF(flo(1)==0.0_sp) THEN
       nonlinear_solve=x1 ! The lower limit is the solution
    ELSEIF(fhi(1)==0.0_sp) THEN
       nonlinear_solve=x2 ! The upper limit is the solution
    ELSEIF(flo(1)<0.0_sp) THEN
       xlo=x1 ! Set it such that f(xlo)<0
       xhi=x2
    ELSE
       xhi=x1
       xlo=x2
    END IF
    nonlinear_solve=0.5_sp*(x1+x2) ! Start with bisection
    dxold=ABS(x2-x1)   !stepsize before last stepsize
    dx=dxold
    ! Evaluate function at the current location
    f=funcd(nonlinear_solve,cell,magma)

    maxit=10000
    DO jj=1,maxit
       !! Bisect if the Newton iteration is
       !! out of range or not decreasing fast enough
       IF(((nonlinear_solve-xhi)*f(2)-f(1))*&
            &((nonlinear_solve-xlo)*f(2)-f(1))>0.0_sp.OR.&
            &ABS(2.0*f(1))>ABS(dxold*f(2))) THEN
          dxold=dx
          dx=0.5_sp*(xhi-xlo)
          nonlinear_solve=xlo+dx

          IF(xlo==nonlinear_solve) RETURN
       ELSE
          dxold=dx
          dx=f(1)/f(2)
          temp=nonlinear_solve
          nonlinear_solve=nonlinear_solve-dx

          IF(temp==nonlinear_solve) RETURN
       END IF
       IF(ABS(dx)<eps) THEN

          RETURN
       END IF
       f=funcd(nonlinear_solve,cell,magma)      
       IF(f(1)<0.0_sp) THEN
          xlo=nonlinear_solve
       ELSE
          xhi=nonlinear_solve
       END IF
    END DO
  END FUNCTION nonlinear_solve


  FUNCTION funcd(x,cell,magma)
    IMPLICIT NONE
    REAL(sp),INTENT(in)::x
    TYPE(unit_cell),INTENT(in)::cell
    TYPE(melt),INTENT(in)::magma
    REAL(sp),DIMENSION(2)::funcd
    REAL(sp)::dx,cont
    REAL(sp),DIMENSION(4)::effective
    !> This function inverts returns vs_calc-vs_obs
    !! For inverting v_s to melt volume fraction using 
    !! the nonlinear Newton-Raphson/bisection algorithm.
    !! This function calculates the first derivative numerically
    !! f' = (f(x+dx)-f(x))/dx !<

    ! Change the following line to use another microgeodynamic mdel

    !cont=whm(magma%dihedral,x)
    cont=vbw(magma%dihedral,x)

    effective= elastic_melt(cont,x,cell%K_sol,cell%G_sol&
         &,cell%nu,cell%rho,cell%Kl,cell%rhol)

    funcd(1)=effective(1)*cell%vs_sol-cell%vs_obs
    !funcd(1)=effective(1)-cell%vp_obs
    dx=5.0e-4_sp
    !cont=whm(magma%dihedral,x+dx)
    cont=vbw(magma%dihedral,x+dx)
    effective= elastic_melt(cont,x+dx,cell%K_sol,cell%G_sol&
         &,cell%nu,cell%rho,cell%Kl,cell%rhol)
    funcd(2)=(effective(1)*cell%vs_sol-cell%vs_obs-funcd(1))/dx
    !funcd(2)=(effective(1)-cell%vp_obs-funcd(1))/dx

  END FUNCTION funcd

  FUNCTION forward_melting(melt_fraction,magma,mymantle,compfile,outfile,&
       &geometry_option)
    IMPLICIT NONE
    TYPE(mantle),INTENT(in)::mymantle
    TYPE(melt),INTENT(in)::magma
    REAL(sp),INTENT(in)::melt_fraction
    integer,intent(in)::geometry_option
    REAL(sp)::cont,Kl,rhol,pressure
    real(sp),dimension(:,:),pointer::forward_melting
    REAL(sp),DIMENSION(4)::effective
    CHARACTER(len=*),INTENT(in)::compfile,outfile
    TYPE(mantle)::load_composition
    INTEGER::ii,filesize
    !> Geometry option of 1 corresponds to melt tubes anything else 
    !! corresponds to melt films. <!

    OPEN(1,file=outfile,form='formatted')
    filesize=line_numbers(compfile,10000)
    cont=vbw(magma%dihedral,melt_fraction)
    allocate(forward_melting(filesize,4))
   
    do ii=1,filesize
       rhol=magma%rho0
       CALL bm3(magma%K0,rhol,magma%Kp,magma%rho0,pressure,Kl)
       if(geometry_option==1) then
       effective=elastic_melt(cont,melt_fraction,mymantle%K(ii),mymantle%G(ii)&
            &,mymantle%nu(ii),mymantle%rho(ii),Kl,rhol)
       else
          effective=elastic_film(magma%aspect,melt_fraction,mymantle%K(ii),&
               &mymantle%G(ii),mymantle%rho(ii),Kl,rhol)
         
       end if
       forward_melting(ii,1:4)=effective(1:4)
       WRITE(1,'(4(1x,f12.5))'),effective(1:4)
    end do


    CLOSE(1)

  END FUNCTION forward_melting

  FUNCTION temp2vs(T,basalt_fraction)!<model)
    REAL(sp),INTENT(in)::T
    INTEGER,INTENT(in)::basalt_fraction
    !type(mantle),intent(in)::model
    TYPE(fit),DIMENSION(3)::f
    REAL(sp)::temp2vs,temptol
    INTEGER::ii,n
    !> This function returns the Vs at a given temperature and basalt fraction
    !!based on Xu et al. (2008)!<
    f=temperaturefit(basalt_fraction)
    temp2vs=f(2)%p1*T**2+f(2)%p2*T+f(2)%p3
    temp2vs=temp2vs*1.0e3_sp
  END FUNCTION temp2vs

  FUNCTION temp2vp(T,basalt_fraction)
    !> This function returns the Vp at a given temperature and basalt fraction
    !!based on Xu et al. (2008)!<
    REAL(sp),INTENT(in)::T
    INTEGER,INTENT(in)::basalt_fraction
    TYPE(fit),DIMENSION(3)::f
    REAL(sp)::temp2vp
    f=temperaturefit(basalt_fraction)
    temp2vp=f(3)%p1*T**2+f(3)%p2*T+f(3)%p3
    temp2vp=temp2vp*1.0e3_sp
  END FUNCTION temp2vp

  FUNCTION temp2rho(T,basalt_fraction)
    !> This function returns the density at a given temperature and basalt fraction
    !!based on Xu et al. (2008)!<
    REAL(sp),INTENT(in)::T
    INTEGER,INTENT(in)::basalt_fraction
    TYPE(fit),DIMENSION(3)::f
    REAL(sp)::temp2rho
    f=temperaturefit(basalt_fraction)
    temp2rho=f(1)%p1*T**2+f(1)%p2*T+f(1)%p3
    temp2rho=temp2rho*1.0e3_sp
  END FUNCTION temp2rho


  FUNCTION depth2dT(depth,mtzdepth)
    IMPLICIT NONE
    REAL(sp),INTENT(in)::depth,mtzdepth
    REAL(sp)::depth2dT,dzdT,dPdT,dPdz
    !>This function returns the rise or fall in  temperature in K
    !! corresponding to a depression or elevation of the 410 km disocntinuity.
    !! According to Katsura et al (Section 4.2,2006) dzdT=0.1km/K
    !! dPdT=4 MPa/K from Katsura et al
    !! dPdT=3.1 MPa/K from Houser and Williams,2010
    !! dPdT=2.5 MPa/K Katsura and Ito,1989
    !! dT=1900-1600 K
    !! dz =430-400 km
    !! dP=14.5-13.4 GPa
    !<
    dPdz=(14.5_sp-13.4_sp)*1e3_sp/30.0_sp !! MPa/km
    dPdT=3.1 !!MPa/km
    dzdT=dPdT/dPdz
    depth2dT=(depth-mtzdepth)/dzdT
  END FUNCTION depth2dT

  FUNCTION temperaturefit(basalt_fraction)
    IMPLICIT NONE
    INTEGER,INTENT(in)::basalt_fraction
    TYPE(fit),DIMENSION(3)::f,temperaturefit
    INTEGER::ii,n
    REAL(sp),DIMENSION(22,10)::par
    REAL(sp)::X
    REAL(sp),DIMENSION(22,10)::coeffs
    X=REAL(basalt_fraction,kind=sp)/100.0_sp

    coeffs(1,:)=(/0.0_sp,   1.62E-07_sp, -0.00073476_sp,  &
         &   4.3039_sp,   4.35E-07_sp,-0.0021304_sp,    &
         & 7.0076_sp,   6.42E-07_sp, -0.0030743_sp,     11.953_sp/)
    coeffs(2,:)=(/ 0.05_sp,   1.55E-07_sp, -0.00070334_sp,   &
         &  4.2853_sp,   4.15E-07_sp, -0.0020299_sp,     6.9023_sp, &
         &  6.17E-07_sp, -0.0029462_sp,      11.83_sp/)
    coeffs(3,:)=(/ 0.1_sp,   1.46E-07_sp, -0.00067004_sp, &
         &    4.2649_sp,   4.01E-07_sp, -0.0019522_sp,   &
         &  6.8155_sp,   5.94E-07_sp, -0.0028289_sp,     11.716_sp/)
    coeffs(4,:)=(/0.15_sp,   1.33E-07_sp, -0.00061938_sp,  &
         &   4.2297_sp,   3.69E-07_sp, -0.0018162_sp,  &
         &  6.6814_sp,   5.42E-07_sp,  -0.002618_sp,     11.526_sp/)
    coeffs(5,:)=(/0.18_sp,  1.23E-07_sp, -0.00056978_sp,  &
         &   4.1939_sp,   3.36E-07_sp, -0.0016798_sp,     &
         & 6.563_sp,   4.97E-07_sp, -0.0024148_sp,     11.346_sp/)
    coeffs(6,:)=(/0.2_sp,   1.22E-07_sp, -0.00057739_sp,  &
         &   4.2013_sp,   3.39E-07_sp, -0.0016885_sp,    &
         & 6.5548_sp,   4.97E-07_sp, -0.0024321_sp,     11.357_sp/)
    coeffs(7,:)=(/0.25_sp,   1.12E-07_sp, -0.00053771_sp,&
         &     4.1745_sp,   3.07E-07_sp, -0.0015563_sp,  &
         &   6.4254_sp,   4.52E-07_sp, -0.0022473_sp,     11.189_sp/)
    coeffs(8,:)=(/0.3_sp,  1.02E-07_sp, -0.0004991_sp,  &
         &   4.1484_sp,   2.75E-07_sp,  -0.0014213_sp,    &
         & 6.2948_sp,   4.07E-07_sp, -0.0020626_sp,     11.022_sp/)
    coeffs(9,:)=(/0.35_sp,   9.39E-08_sp, -0.00046484_sp, &
         &    4.1254_sp,   2.44E-07_sp, -0.0012891_sp,    &
         & 6.1672_sp,   3.65E-07_sp, -0.0018862_sp,     10.862_sp/)
    coeffs(10,:)=(/0.4_sp,   8.65E-08_sp, -0.0004324_sp,  &
         &   4.1026_sp,   2.12E-07_sp, -0.0011527_sp,    &
         & 6.0348_sp,   3.23E-07_sp, -0.0017067_sp,     10.695_sp/)
    coeffs(11,:)=(/0.45_sp,   7.86E-08_sp, -0.00039711_sp,&
         &     4.0761_sp,   1.77E-07_sp, -0.0010066_sp,  &
         &   5.8934_sp,   2.78E-07_sp, -0.0015126_sp,     10.514_sp/)
    coeffs(12,:)=(/ 0.5_sp,   7.11E-08_sp, -0.00036279_sp, &
         &    4.0501_sp,   1.44E-07_sp, -0.00086483_sp,    &
         & 5.7567_sp,   2.34E-07_sp, -0.0013235_sp,     10.338_sp/)
    coeffs(13,:)=(/0.55_sp,   6.36E-08_sp, -0.00032808_sp, &
         &    4.0237_sp,   1.11E-07_sp, -0.00072765_sp,   &
         &  5.6253_sp,   1.91E-07_sp, -0.0011384_sp,     10.166_sp/)
    coeffs(14,:)=(/0.6_sp,   5.61E-08_sp, -0.00029206_sp,  &
         &   3.9949_sp,   8.06E-08_sp, -0.00059562_sp,    &
         & 5.4966_sp,   1.50E-07_sp, -0.00095617_sp,     9.9932_sp/)
    coeffs(15,:)=(/0.65_sp,   4.66E-08_sp, -0.00024911_sp,&
         &     3.9601_sp,   5.03E-08_sp, -0.0004655_sp,  &
         &   5.3701_sp,   1.07E-07_sp, -0.00076796_sp,   &
         &  9.8152_sp/)
    coeffs(16,:)=(/0.7_sp,   2.87E-08_sp, -0.00018026_sp, &
         &    3.9065_sp,   1.59E-08_sp, -0.00032465_sp,   &
         &  5.2383_sp,   5.08E-08_sp, -0.00054152_sp,     9.6117_sp/)
    coeffs(17,:)=(/0.75_sp,   9.38E-09_sp, -0.00010847_sp,  &
         &    3.851_sp,  -1.04E-08_sp, -0.00021256_sp,    &
         & 5.1297_sp,   8.97E-10_sp, -0.00033815_sp,     9.4268_sp/)
    coeffs(18,:)=(/0.8_sp,  -5.48E-09_sp,  -5.38E-05_sp,  &
         &   3.8096_sp,  -3.18E-08_sp, -0.00011513_sp,    &
         & 5.0327_sp,  -3.97E-08_sp, -0.00016716_sp,      9.268_sp/)
    coeffs(19,:)=(/0.85_sp,  -1.80E-08_sp,  -1.31E-05_sp, &
         &    3.7954_sp,  -5.42E-08_sp,  -4.69E-05_sp,    &
         & 5.0235_sp,  -6.82E-08_sp,  -7.56E-05_sp,     9.2518_sp/)
    coeffs(20,:)=(/0.9_sp,  -3.04E-08_sp,   2.48E-05_sp,  &
         &   3.7895_sp,  -7.45E-08_sp,   2.26E-06_sp,     &
         &5.0545_sp,  -9.19E-08_sp,  -1.52E-05_sp,     9.2915_sp/)
    coeffs(21,:)=(/0.95_sp, -4.28E-08_sp,    6.23E-05_sp , &
         &   3.7849_sp,  -9.53E-08_sp,   5.22E-05_sp,    &
         & 5.0896_sp,  -1.16E-07_sp,   4.58E-05_sp,     9.3371_sp/)
    coeffs(22,:)=(/1.0_sp,  -5.52E-08_sp,   9.97E-05_sp,   &
         &  3.7808_sp,  -1.17E-07_sp, 0.00010339_sp,     &
         &5.1261_sp,   1.41E-07_sp,  0.00010825_sp,     9.3848_sp/)

    DO ii=1,22
       IF(X==coeffs(ii,1)) THEN
          !Density
          f(1)%p1=coeffs(ii,2)
          f(1)%p2=coeffs(ii,3)
          f(1)%p3=coeffs(ii,4)
          !Vs
          f(2)%p1=coeffs(ii,5)
          f(2)%p2=coeffs(ii,6)
          f(2)%p3=coeffs(ii,7)
          !Vp
          f(3)%p1=coeffs(ii,8)
          f(3)%p2=coeffs(ii,9)
          f(3)%p3=coeffs(ii,10)
       END IF
    END DO
    temperaturefit=f
  END FUNCTION temperaturefit

  FUNCTION PREM_vs_vp_rho(depth)
    ! returns vs vp and density in SI units. Enter depth in km
    IMPLICIT NONE
    REAL(sp),INTENT(in)::depth
    REAL(sp),DIMENSION(3)::PREM_vs_vp_rho
    REAL(sp)::vs,vp,rho,x,r

    r=6371.0_sp-depth
    x=r/6371.0_sp

    IF (r>3480.0_sp.AND.r<=3630.0_sp) THEN

       rho=7.9565_sp-6.4761_sp*x+5.5283_sp*x**2-3.0807_sp*x**3
       vp=15.3891_sp-5.3181_sp*x+5.5242_sp*x**2-2.5514_sp*x**3
       vs=6.9254_sp+1.4672_sp*x-2.0834_sp*x**2+0.9783_sp*x**3

    ELSEIF(r>3630.0_sp.AND.r<=5600.0_sp) THEN
       rho=7.9565_sp-6.4761_sp*x+5.5283_sp*x**2-3.0807_sp*x**3
       vp=24.952_sp-40.4673_sp*x+51.4832_sp*x**2-26.6419_sp*x**3
       vs=11.671_sp-13.7818_sp*x+17.4575_sp*x**2-9.2777_sp*x**3

    ELSEIF(r>5600.0_sp.AND.r<=5701.0_sp) THEN
       rho=7.9565_sp-6.4761_sp*x+5.5283_sp*x**2-3.0807_sp*x**3
       vp=29.2766_sp-23.6027_sp*x+5.5242_sp*x**2-2.5514_sp*x**3
       vs=22.3459_sp-17.2473_sp*x-2.0834_sp*x**2+0.9783_sp*x**3
    ELSEIF(r>5701.0_sp.AND.r<=5771.0_sp) THEN
       rho=5.3197_sp-1.4836_sp*x
       vp=19.0957_sp-9.8672_sp*x
       vs=9.9839_sp-4.9324_sp*x
    ELSEIF(r>5771.0_sp.AND.r<=5971.0_sp) THEN
       rho=11.2494_sp-8.0298_sp*x
       vp=39.7027_sp-32.6166_sp*x
       vs=22.3512_sp-18.5856_sp*x
    ELSEIF(r>5971.0_sp.AND.r<=6151.0_sp) THEN
       rho=7.1089_sp-3.8045_sp*x
       vp=20.3926_sp-12.2569_sp*x
       vs=8.9496_sp-4.4597_sp*x
    ELSEIF(r>6151.0_sp.AND.r<=6291.0_sp) THEN
       PRINT*,'Depth too shallow'
    END IF
    PREM_vs_vp_rho(1)=vs*1.0e3_sp
    PREM_vs_vp_rho(2)=vp*1.0e3_sp
    PREM_vs_vp_rho(3)=rho*1.0e3_sp
  END FUNCTION PREM_vs_vp_rho

END MODULE microgeodynamics
