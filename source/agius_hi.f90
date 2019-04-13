PROGRAM agius_hi
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
  REAL(sp),DIMENSION(3)::prem,temp
  REAL(sp),DIMENSION(5)::potential_temperature
  real(sp),dimension(7)::dihedral

  potential_temperature=(/1400.0_sp,1500.0_sp,1600.0_sp,1700.0_sp,1800.0_sp/)
  basalt = (/10,20,30,40,50,60,70,80,90,99/)
  MTZ_TOPO_NORMALIZE_HAWAII=.TRUE.
  ! Create the parameter space for dihedral angles
  dihedral_angle=linspace(5.0_sp,30.0_sp,10)
  do kk=1,10
     do ii=1,5
        ! Define composition space 
        HAWAII_REF_COMP=composition(potential_temperature&
             &=potential_temperature(ii),basalt_fraction=basalt(kk))
        ! Load solid composition for 352 km from Xu et al (2008)
        HAWAII_REF_MANTLE=load_composition(HAWAII_REF_COMP&
             &,'../input/352.dat')          
        ! Load the regional data
        ! The first call uses temperature calculated from MTZ thickness
        HAWAII=load_seismo_hmt(HAWAII_REF_COMP,3&
        !     &,'../input/Hawaii_temperature_agius.csv',MTZ_TOPO_NORMALIZE_HAWAII,velocity_contrast_410,14)
        !This call uses temperature calculated from 410 totpography
        !HAWAII=load_seismo_hmt(HAWAII_REF_COMP,3&
             &,'../input/Hawaii_temperature_from410topo_agius.csv',MTZ_TOPO_NORMALIZE_HAWAII,velocity_contrast_410,14)
        ! loop through dihedral angle values
        dihedral=(/10.0_sp,15.0_sp,20.0_sp,25.0_sp,30.0_sp,35.0_sp,40.0_sp/)
        do jj=1,7
           ! Set the melt composition
           HAWAII%MAGMA=set_melt(5,dihedral(jj)) 

           PRINT*,  'HMT: initialization successful'
           ! Calculate melt volume fraction from the observed data
           HAWAII=melt_calculate(HAWAII)

           ! Write the data for various formats
           !CALL vtk_write(HAWAII,HAWAII_REF_MANTLE,'../data/HMT/')
           CALL data_table_output(HAWAII,HAWAII_REF_MANTLE,'../data/AGIUS/')
           !CALL melt_gmt_output(HAWAII,HAWAII_REF_MANTLE,'../data/AGIUS/')
        end do

     end do
  end do


CONTAINS 




END PROGRAM agius_hi
