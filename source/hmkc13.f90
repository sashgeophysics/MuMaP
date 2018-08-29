PROGRAM hmkc13
  USE global
  USE microgeodynamics
  USE regional
  !TYPE(unit_cell)::PREM
  TYPE(region)::HAWAII,CORAL_SEA
  TYPE(composition)::HAWAII_REF_COMP
  TYPE(mantle)::HAWAII_REF_MANTLE
  LOGICAL::MTZ_TOPO_NORMALIZE_HAWAII,MTZ_TOPO_NORMALIZE_CORAL_SEA
  INTEGER::ii,jj,kk
  INTEGER,DIMENSION(10)::basalt,potential
  REAL(sp),DIMENSION(10)::dihedral_angle
  INTEGER :: clock_start,clock_end,clock_rate
  REAL(sp) :: elapsed_time,theta
  real(sp),dimension(3)::prem,temp
 

  CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
  CALL SYSTEM_CLOCK(COUNT=clock_start) ! Start timing

  ! Create the parameter space for basalt fractions
  basalt=(/00,05,10,15,18,20,25,30,35,40/)
  ! Create the parameter space for potential temperature

  potential=linspace(1400,2300,10)
  ! Create the parameter space for dihedral angles
  dihedral_angle=linspace(5.0_sp,30.0_sp,10)
  MTZ_TOPO_NORMALIZE_HAWAII=.TRUE.
  MTZ_TOPO_NORMALIZE_CORAL_SEA=.FALSE.
  DO ii=1,10
     DO jj=1,10
        ! Define composition space 
        HAWAII_REF_COMP=composition(potential_temperature&
             &=potential(ii),basalt_fraction=basalt(jj))
        ! Load solid composition for 352 km from Xu et al (2008)
        HAWAII_REF_MANTLE=load_composition(HAWAII_REF_COMP&
             &,'../input/352.dat')          
        ! Load the regional data
        HAWAII=load_seismo_variable_temp(HAWAII_REF_COMP,3&
             &,'../input/&
             &khm13_input.csv',435.0_sp,MTZ_TOPO_NORMALIZE_HAWAII,velocity_contrast_410) 
        CORAL_SEA=load_seismo_variable_temp(HAWAII_REF_COMP,3&
             &,'../input/&
             &Coral_Sea.csv',420.0_sp,MTZ_TOPO_NORMALIZE_CORAL_SEA,impedance_contrast_410)
        DO kk=1,10
           ! Set the melt composition
           HAWAII%MAGMA=set_melt(5,dihedral_angle(kk)) 
           CORAL_SEA%MAGMA=set_melt(5,dihedral_angle(kk))
           ! Calculate melt volume fraction from the observed data
           HAWAII=melt_calculate(HAWAII)
           CORAL_SEA=melt_calculate(CORAL_SEA)
           ! Write the data for various formats
           CALL vtk_write(HAWAII,HAWAII_REF_MANTLE,'../data/HI/')
           CALL data_table_output(HAWAII,HAWAII_REF_MANTLE,'../data/HI/')
           CALL melt_gmt_output(HAWAII,HAWAII_REF_MANTLE,'../data/HI/')
           CALL vtk_write(CORAL_SEA,HAWAII_REF_MANTLE,'../data/CS/')
           CALL data_table_output(CORAL_SEA,HAWAII_REF_MANTLE,'../data/CS/')
           CALL melt_gmt_output(CORAL_SEA,HAWAII_REF_MANTLE,'../data/CS/')
           !CALL inversion_test(HAWAII)
           !call compare_prediction_observation(HAWAII)
        END DO
     END DO
  END DO
  CALL SYSTEM_CLOCK(COUNT=clock_end) ! Stop timing
  ! Calculate the elapsed time in seconds:
  elapsed_time=REAL((clock_end-clock_start)/clock_rate,sp)


  PRINT*,'elapsed_time',elapsed_time
  PRINT*,'# of iterations',(ii-1)*(jj-1)*(kk-1)


CONTAINS 


  SUBROUTINE inversion_test(reg)
    IMPLICIT NONE
    TYPE(region),INTENT(in)::reg
    TYPE(region)::reg2
    INTEGER::ii
    REAL(sp)::cont
    REAL(sp),DIMENSION(4)::temp


    reg2=reg

    reg2%CELL(1:reg%NGRID)%melt_fraction=linspace_r(0.001_sp,&
         &0.1_sp,reg%NGRID)

!!$    reg2%CELL(1:reg%NGRID)%K_sol=655.6e9_sp!GPa
!!$    reg2%CELL(1:reg%NGRID)%G_sol=293.8e9_sp!GPa
!!$    reg2%CELL(1:reg%NGRID)%vp_sol=13716.6_sp! m/s
!!$    reg2%CELL(1:reg%NGRID)%vs_sol=7264.66_sp
!!$    reg2%CELL(1:reg%NGRID)%rho=5566.45_sp
!!$    reg2%CELL(1:reg%NGRID)%nu=0.3051_sp
!!$
!!$    reg2%CELL(1:reg%NGRID)%Kl=452.98_sp!GPa
!!$    reg2%CELL(1:reg%NGRID)%rhol=2.4_sp*2690.0_sp
    PRINT*,'================================================='
    PRINT*,'Error testing for melt fraction inversion'
    PRINT*,'================================================='
    WRITE(*,'(5(1x,a))'),'|Original', '  |Nonlinear', '|Error',  '|  Bisection', '|   Error|'
    PRINT*,'================================================='
    DO ii=1,reg%NGRID
       cont=whm(reg2%MAGMA%dihedral,reg%CELL(ii)%melt_fraction)
       temp=elastic_melt(cont,reg2%CELL(ii)%melt_fraction&
            &,reg2%CELL(ii)%K_sol,reg2%CELL(ii)%G_sol,reg2%CELL(ii)%nu,&
            &reg2%CELL(ii)%rho,reg2%CELL(ii)%Kl,reg2%CELL(ii)%rhol)
       reg2%CELL(ii)%vs_eff=temp(1)
       reg2%CELL(ii)%vp_eff=temp(2)
       reg2%CELL(ii)%K_eff=temp(3)
       reg2%CELL(ii)%G_eff=temp(4)
       reg2%CELL(ii)%vs_obs=nonlinear_solve(0.0001_sp,0.45_sp,1.0e-4_sp,reg2%CELL(ii),reg2%MAGMA)
       reg2%CELL(ii)%vp_obs=bisection(0.0001_sp,0.45_sp,1.0e-4_sp,reg2%CELL(ii),reg2%MAGMA)
       WRITE(*,'(5(1x,es10.2))'),reg2%CELL(ii)%melt_fraction,reg2%CELL(ii)%vs_obs,ABS(reg2%CELL(ii)%vs_obs-reg2%CELL(ii)%melt_fraction),reg2%CELL(ii)%vp_obs,ABS(reg2%CELL(ii)%vp_obs-reg2%CELL(ii)%melt_fraction)
    END DO
  END SUBROUTINE inversion_test

  SUBROUTINE compare_prediction_observation(reg)
    IMPLICIT NONE
    TYPE(region),INTENT(in)::reg
    TYPE(region)::reg2
    INTEGER::ii
    REAL(sp)::cont
    REAL(sp),DIMENSION(4)::temp
    reg2=reg
    !reg2%CELL(1:reg%NGRID)%melt_fraction=0.02_sp
    PRINT*,'================================================='
    WRITE(*,'(3(1x,a))'),'|melt', '  |predicted ', '|observed|'
    PRINT*,'================================================='
    DO ii=1,reg%NGRID
       cont=vbw(reg%MAGMA%dihedral,reg%CELL(ii)%melt_fraction)
       temp=elastic_melt(cont,reg%CELL(ii)%melt_fraction&
            &,reg%CELL(ii)%K_sol,reg%CELL(ii)%G_sol,reg%CELL(ii)%nu,&
            &reg%CELL(ii)%rho,reg%CELL(ii)%Kl,reg%CELL(ii)%rhol)
       reg2%CELL(ii)%vs_eff=temp(1)
       reg2%CELL(ii)%vp_eff=temp(2)
       reg2%CELL(ii)%K_eff=temp(3)
       reg2%CELL(ii)%G_eff=temp(4)
       WRITE(*,'(3(1x,es10.2))'),reg%CELL(ii)%melt_fraction,reg%CELL(ii)%vs_obs/reg%CELL(ii)%vs_sol,temp(1)
    END DO
  end SUBROUTINE compare_prediction_observation
END PROGRAM hmkc13
