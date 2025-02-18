!> This module contains data types and functions applicable
  !!to the region of interest. Each location in the region
  !!is assigned a unit cell, which contains the physical properties.
  !!The parameters for the equation of state of the melt phase are
  !!stored in the derived type MAGMA. It is assumed that only one
  !!kind of melt exists in the region. Three dimensional location
  !!of each point in the region are stored in the array LOC. Typically
  !!this data is stored as lat,lon, depth (km). The variable composition
  !!is also regional, it contains the basalt fraction
  !!of the bulk ocmposition and the average potential temperature
  !!of the region.
  !! Created by Saswata Hier-Majumder August, 2012.!<
MODULE regional
  USE global
  USE microgeodynamics

  TYPE region
     INTEGER::NGRID
     TYPE(composition)::comp
     TYPE(unit_cell),DIMENSION(:),POINTER::CELL
     TYPE(melt)::MAGMA
     REAL(sp),DIMENSION(:,:),POINTER::LOC
     REAL(sp),DIMENSION(:,:),POINTER::TOPO
     real(SP)::MTZ                      ! Depth to transition zone in km
  END TYPE region


CONTAINS

  FUNCTION load_seismo_variable_temp(c,dim,fname,mtzdepth,NORM,func,depth_in)
    IMPLICIT NONE
    TYPE(composition),INTENT(in)::c
    real(sp),intent(in)::mtzdepth
    INTEGER,INTENT(in)::dim
    logical,intent(in)::NORM
    CHARACTER(len=*)::fname
    REAL(sp),INTENT(in),OPTIONAL::depth_in
    REAL(sp)::z,impedance350,prem_impedance
    REAL(sp),DIMENSION(2)::impedance_vel
    REAL(sp),DIMENSION(:,:),POINTER::DATA
    REAL(sp),DIMENSION(:,:),POINTER::x
    REAL(sp),DIMENSION(3)::prem_value
    TYPE(region)::load_seismo_variable_temp
    INTEGER::ii
    interface
       function func(int)
         use global
         IMPLICIT NONE
         INTEGER,INTENT(in)::int
         REAL(sp)::func
       end function func
    end interface

    !> Reads in input data from seismological models.
    !!The format of the input data
    !! should be  the following:
    !! Lat Lon Depth Vs Vp
    !! If the data set is two dimensional,
    !! then an optional argument, depth_in (km) 
    !! needs to be supplied. This feature is not
    !! yet fully implemented. !<

    load_seismo_variable_temp%comp=c
    load_seismo_variable_temp%NGRID=line_numbers(fname,10000)
    load_seismo_variable_temp%MTZ=mtzdepth

    ALLOCATE(DATA(load_seismo_variable_temp%NGRID,11),x(load_seismo_variable&
         &_temp%NGRID,3))
    OPEN(10,file=fname,status='unknown',form='formatted')
    DO ii=1,load_seismo_variable_temp%NGRID
       READ(10,'(11(1x,es10.3))'),DATA(ii,:)
    END DO
    CLOSE(10)

    ALLOCATE(load_seismo_variable_temp%CELL(load_seismo_variable_temp%NGRID),&
         &load_seismo_variable_temp%LOC(load_seismo_variable_temp%NGRID,3),&
         &load_seismo_variable_temp%TOPO(load_seismo_variable_temp%NGRID,2))

    !> Here we load the observed variables into the region !<
    DO ii=1,load_seismo_variable_temp%NGRID
       ! Location in Lat,Lon,Depth(km)
       load_seismo_variable_temp%LOC(ii,1:3)=DATA(ii,1:3)

       !Depth of top of LVL and 410 in (km)
       load_seismo_variable_temp%TOPO(ii,1:2)=DATA(ii,3:4)

       !>calculate potential temperature from depression of 410
       !!Using function depth2dT in module microgeodynamics!<
       load_seismo_variable_temp%CELL(ii)%temperature=DBLE(load_seismo_variable&
            &_temp%comp&
            &%potential_temperature)+depth2dT(DATA(ii,4),mtzdepth)

       !> Convert the temperature to  Shear wave velocity
       !!using the function temp2vs in module microgeodynamics!<
       load_seismo_variable_temp%CELL(ii)%vs_sol=temp2vs(load_seismo_variable&
            &_temp%CELL(ii)&
            &%temperature,load_seismo_variable_temp%comp%basalt_fraction)

       !> Convert the temperature to  P wave velocity
       !!using the function temp2vp in module microgeodynamics!<
       load_seismo_variable_temp%CELL(ii)%vp_sol=temp2vp(load_seismo_variable&
            &_temp%CELL(ii)&
            &%temperature,&
            &load_seismo_variable_temp%comp%basalt_fraction)

       !> Convert the temperature to density using the function
       !! temp2rho in module microgeodynamics!<
       load_seismo_variable_temp%CELL(ii)%rho=temp2rho(load_seismo_variable&
            &_temp%CELL(ii)&
            &%temperature,&
            &load_seismo_variable_temp%comp%basalt_fraction)
       !Calculate the corresponding solid K,G, and nu
       CALL vel2mod(load_seismo_variable_temp%CELL(ii)%vs_sol,load_seismo&
            &_variable_temp%CELL(ii)&
            &%vp_sol,&
            &load_seismo_variable_temp%CELL(ii)%rho,load_seismo_variable&
            &_temp%CELL(ii)%G_sol,&
            &load_seismo_variable_temp%CELL(ii)%K_sol,load_seismo_variable&
            &_temp%CELL(ii)%nu )
       !> Calculate observed Shear wave velocity from impedance
       !!contrast observation!<
       impedance_vel=func(1)
       if (NORM .eqv. .TRUE.) then
       impedance350=DATA(ii,7)*impedance_vel(1)
       else
          impedance350=DATA(ii,8)
       end if
       ! Impedance contrast from the PREM model across 400
       !Calculate vs from the impedance contrast and PREM value
       prem_value=PREM_vs_vp_rho(350.0_sp)
          load_seismo_variable_temp%CELL(ii)%vs_obs=prem_value(1)&
           &*(1.0+impedance350)   
       !Impedance contrast at 350
       load_seismo_variable_temp%CELL(ii)%vp_obs=1.0+impedance350
    END DO

    DATA=>NULL()
    x=>NULL()
  END FUNCTION load_seismo_variable_temp

  FUNCTION load_seismo_hmt(c,dim,fname,NORM,func,ncols)
    IMPLICIT NONE
    TYPE(composition),INTENT(in)::c
    INTEGER,INTENT(in)::dim
    logical,intent(in)::NORM
    CHARACTER(len=*)::fname
    character(len=24)::colnum,fmt
    INTEGER,INTENT(in)::ncols ! The number ofcolumns in the input file
    REAL(sp)::z,impedance350,prem_impedance
    REAL(sp),DIMENSION(2)::impedance_vel
    REAL(sp),DIMENSION(:,:),POINTER::DATA,temp1
    REAL(sp),DIMENSION(:,:),POINTER::x
    REAL(sp),DIMENSION(3)::prem_value
    TYPE(region)::load_seismo_hmt
    INTEGER::ii
    INTERFACE
       FUNCTION func(int)
         USE global
         IMPLICIT NONE

         INTEGER,INTENT(in)::int
         REAL(sp)::func
       END FUNCTION func
    END INTERFACE

    !> Reads in input data from seismological models.
    !!The format of the input data
    !! should be  the following:
    !! Lat Lon 350 410  Rnorm dT
    !! If there are any additional column
    !! data added afterwards. Those should be
    !! properly commented. !<

    load_seismo_hmt%comp=c
    !load_seismo_hmt%NGRID=line_numbers(fname,10000)
    ! This value is only for Agius Hawaii data, change it
    load_seismo_hmt%NGRID=1490
    !> Read the number of columns into a character and 
    !! format the input !<

    WRITE(colnum,'(i3)'),ncols
    fmt="("//trim(colnum)//"(1x,e12.5))"
    ALLOCATE(DATA(load_seismo_hmt%NGRID,ncols),&
         &temp1(ncols,load_seismo_hmt%NGRID),&
         &x(load_seismo_hmt%NGRID,3))
    OPEN(10,file=fname,status='unknown',form='formatted')    
    READ(10,fmt),temp1
    
    CLOSE(10)
    DATA=transpose(temp1)
    print*,'load_seismo_hmt: max(amp),min(amp)',maxval(DATA(:,5)),minval(DATA(:,5))
    ALLOCATE(load_seismo_hmt%CELL(load_seismo_hmt%NGRID),&
         &load_seismo_hmt%LOC(load_seismo_hmt%NGRID,3),&
         &load_seismo_hmt%TOPO(load_seismo_hmt%NGRID,2))

    !> Here we load the observed variables into the region !<
    DO ii=1,load_seismo_hmt%NGRID
       ! Location in Lat,Lon,Depth(km)h
       load_seismo_hmt%LOC(ii,1:3)=DATA(ii,1:3)
       
       !Depth of top of LVL and 410 in (km)
       load_seismo_hmt%TOPO(ii,1:2)=DATA(ii,3:4)
       

       !>calculate potential temperature from input data
       !!which should contain the temperature anomaly <!
       !> The input data contains dT for claperon slopes 0.5-4.5MPa/K
       !! with increments of 0.5 starting from column 5 !>

       ! For 3 MPa/K read from column 11
       ! For 4 MPa/K read from column 13
       load_seismo_hmt%CELL(ii)%temperature=DBLE(load_seismo_hmt&
            &%comp%potential_temperature) + DATA(ii,11) + 0.3_sp*350.0_sp
       !> Convert the temperature to  Shear wave velocity
       !!using the function temp2vs in module microgeodynamics!<
       load_seismo_hmt%CELL(ii)%vs_sol=temp2vs(load_seismo_hmt%CELL(ii)&
            &%temperature,load_seismo_hmt%comp%basalt_fraction)

       !> Convert the temperature to  P wave velocity
       !!using the function temp2vp in module microgeodynamics!<
       load_seismo_hmt%CELL(ii)%vp_sol=temp2vp(load_seismo_hmt%CELL(ii)&
            &%temperature,&
            &load_seismo_hmt%comp%basalt_fraction)

       !> Convert the temperature to density using the function
       !! temp2rho in module microgeodynamics!<
       load_seismo_hmt%CELL(ii)%rho=temp2rho(load_seismo_hmt%CELL(ii)&
            &%temperature,&
            &load_seismo_hmt%comp%basalt_fraction)
       !Calculate the corresponding solid K,G, and nu
       CALL vel2mod(load_seismo_hmt%CELL(ii)%vs_sol,load_seismo&
            &_hmt%CELL(ii)&
            &%vp_sol,&
            &load_seismo_hmt%CELL(ii)%rho,load_seismo_hmt%CELL(ii)%G_sol,&
            &load_seismo_hmt%CELL(ii)%K_sol,load_seismo_hmt%CELL(ii)%nu )
       !> Calculate observed Shear wave velocity from impedance
       !!contrast observation!<
       impedance_vel=func(1)
       IF (NORM .eqv. .TRUE.) THEN
          impedance350=DATA(ii,5)*impedance_vel(1)
       ELSE
          impedance350=DATA(ii,5)
       END IF
       ! Impedance contrast from the PREM model across 400
       !Calculate vs from the impedance contrast and PREM value
       prem_value=PREM_vs_vp_rho(350.0_sp)
       load_seismo_hmt%CELL(ii)%vs_obs=prem_value(1)&
            &*(1.0+impedance350)   
       !Impedance contrast at 350
       load_seismo_hmt%CELL(ii)%vp_obs=1.0+impedance350
       ! Load the raw amplitude data into the cell
       load_seismo_hmt%CELL(ii)%amplitude=DATA(ii,5)
    END DO
    
    DATA=>NULL()
    temp1=>NULL()
    x=>NULL()
  END FUNCTION load_seismo_hmt

  FUNCTION melt_calculate(reg)
    IMPLICIT NONE
    TYPE(region),INTENT(in)::reg
    TYPE(region):: melt_calculate
    INTEGER::ii
    REAL(sp),DIMENSION(4)::temp
    REAL(sp)::cont
    !> This function returns the bulk modulus and pressure from
    !!the properties of melt set in function melt_set in module
    !!microgeodynamics. Either the Vinet or a third order Birch-
    !!Murnaghan EOS can be used for the calculation. Just replace
    !!vinet with bm3 in the following line.
    !!Since both EOS are implicit equations, the density at depth
    !!must be known. If not known ahead of time, run a few trials
    !!and read in the pressure. Once the pressure is satisfactory
    !! for the region of interest, then use the prefactor.
    !!This process of iteration can also be automated either within 
    !!the function or by calling an external function.
    !!Once the bulk modulus of the melt is known, it uses a nonlinear
    !!root finder to calculate the melt volume fraction based on the
    !!discrepancy between the predicted and observed Shear wave
    !!velocities.  Currently, two solvers are implemented, a
    !!bisection algorithm and a combined bisection/ Newton-Raphson
    !!algorithm. The combined algorithm is called nonlinear_solver.
    !!The bisection algortihm is called bisection. They both take
    !!the same arguments. See more on these routines in the module
    !!microgeodynamics.!<
    melt_calculate=reg
    DO ii=1,melt_calculate%NGRID
       !! Calculate elastic constants of the melt
       IF(melt_calculate%MAGMA%K0==18.1e9_sp) THEN !MORB
          melt_calculate%CELL(ii)%rhol=melt_calculate%MAGMA%rho0*1.3275_sp
          CALL vinet(melt_calculate%MAGMA%K0,&
               &melt_calculate%CELL(ii)%rhol,melt_calculate%MAGMA%Kp,&
               &melt_calculate%MAGMA%rho0,&
               &melt_calculate%CELL(ii)%pressure,melt_calculate&
               &%CELL(ii)%Kl)
       ELSE !carbonated peridotite
          melt_calculate%CELL(ii)%rhol=melt_calculate%MAGMA%rho0*1.282_sp
          CALL bm3(melt_calculate%MAGMA%K0,&
               &melt_calculate%CELL(ii)%rhol,melt_calculate%MAGMA%Kp,&
               &melt_calculate%MAGMA%rho0,&
               &melt_calculate%CELL(ii)%pressure,melt_calculate&
               &%CELL(ii)%Kl)
       END IF

       !! Calculate melt volume fraction from the known vs
       melt_calculate%CELL(ii)%melt_fraction=nonlinear_solve(1.0e-6_sp,&
            &0.1_sp,&
            &1.0e-9_sp,&
            &melt_calculate%CELL(ii),melt_calculate%MAGMA)

       cont=vbw(melt_calculate%MAGMA%dihedral,&
            &melt_calculate%CELL(ii)%melt_fraction)
       temp=elastic_melt(cont,melt_calculate%CELL(ii)%melt_fraction&
            &,melt_calculate%CELL(ii)%K_sol,melt_calculate%CELL(ii)%G_sol,&
            &melt_calculate%CELL(ii)%nu,&
            &melt_calculate%CELL(ii)%rho,melt_calculate%CELL(ii)%Kl,&
            &melt_calculate%CELL(ii)%rhol)
       melt_calculate%CELL(ii)%vs_eff=temp(1) ! back calculated vs

    END DO

  END FUNCTION melt_calculate


  SUBROUTINE vtk_write(reg,mod,dir)
    IMPLICIT NONE
    TYPE(region),INTENT(in)::reg
    TYPE(mantle),INTENT(in)::mod
    CHARACTER(len=11),intent(in)::dir
    CHARACTER(len=26)::fname
    character(len=100)::dirf
    INTEGER::ii,N,NNODES

    !> This subroutine writes the regional data into a vtk file.
    !!Currently, the data is written as unstructured grid. !<
    fname=dir//filename(reg%MAGMA%dihedral,reg%comp%basalt_fraction,reg%comp%&
         &potential_temperature,'.vtk')

    N=reg%NGRID
    OPEN(10,file=fname,form='formatted')
    WRITE(10,'(''# vtk DataFile Version 2.0'')')
    WRITE(10,'(''Seismology data'')')
    WRITE(10,'(''ASCII'')')
    WRITE(10,'(''           '')')
    WRITE(10,'(''DATASET UNSTRUCTURED_GRID'')')

    WRITE(10,'("POINTS",i7,1x "float")'),N
    DO ii=1,N ! write the position of the vertices
       WRITE(10,'(3(1x,f20.7))'), reg%LOC(ii,2),reg%LOC(ii,1),0.0_sp
    END DO
    WRITE(10,'(''           '')')
    WRITE(10,'("POINT_DATA",i7)'),N

    WRITE(10,'("SCALARS LVL float ")')
    WRITE(10,'("LOOKUP_TABLE default")')
    DO  ii=1,N
       WRITE(10,'(1x,f12.5)'),reg%TOPO(ii,1)
    END DO
    WRITE(10,'("SCALARS TZ float ")')
    WRITE(10,'("LOOKUP_TABLE default")')
    DO  ii=1,N
       WRITE(10,'(1x,f12.5)'),reg%TOPO(ii,2)
    END DO
    WRITE(10,'("SCALARS Temperature float ")')
    WRITE(10,'("LOOKUP_TABLE default")')
    DO  ii=1,N
       WRITE(10,'(1x,f12.5)'),reg%CELL(ii)%temperature
    END DO
    WRITE(10,'("SCALARS vs_obs float ")')
    WRITE(10,'("LOOKUP_TABLE default")')
    DO  ii=1,N
       WRITE(10,'(1x,f12.5)'),reg%CELL(ii)%vs_obs
    END DO
    WRITE(10,'("SCALARS vs_sol float ")')
    WRITE(10,'("LOOKUP_TABLE default")')
    DO  ii=1,N
       WRITE(10,'(1x,f12.5)'),reg%CELL(ii)%vs_sol
    END DO
    WRITE(10,'("SCALARS melt float ")')
    WRITE(10,'("LOOKUP_TABLE default")')
    DO  ii=1,N
       WRITE(10,'(1x,f12.5)'),reg%CELL(ii)%melt_fraction
    END DO

    CLOSE(10)
  END SUBROUTINE vtk_write

  SUBROUTINE data_table_output(reg,mod,dir)
    IMPLICIT NONE
    TYPE(region),INTENT(in)::reg
    TYPE(mantle),INTENT(in)::mod
    CHARACTER(len=11),intent(in)::dir
    CHARACTER(len=26)::fname
    INTEGER::ii,N,NNODES
    fname=dir//filename(reg%MAGMA%dihedral,reg%comp%basalt_fraction,&
         &reg%comp%potential_temperature,'.csv')
    N=reg%NGRID
    OPEN(10,file=fname,form='formatted')
!!$    WRITE(10,'(a,i2)'),'//Basalt percent: ',reg%comp%basalt_fraction
!!$    WRITE(10,'(a,i4)'),'//Potential temperature: ',reg%comp&
!!$         &%potential_temperature
!!$    WRITE(10,'(a,es10.3)'),'//Pressure (Pa):',reg%CELL(1)%pressure
!!$    WRITE(10,'(a,es10.3)'),'//Melt Density (kg/m^3):',reg%CELL(1)%rhol
!!$
!!$    WRITE(10,'(12(1x,a))'),'    //Latitude', '   Longitude',&
!!$         &'    LVL(km)','    410(km)',&
!!$         &' Potential T (K)',' Vs_obs(m/s)', ' Vs_sol(m/s)','Melt fraction',&
!!$         & 'K/K0', 'G/G0','REsidual', 'Amplitude'
    print*,'data_table_output: max(amp),min(amp)',maxval(reg%CELL%amplitude),minval(reg%CELL%amplitude)
    DO ii=1,N
       WRITE(10,'(10(1x,f12.5),2(1x,es16.7),(1x,f12.5))') reg%LOC(ii,1),&
            &reg%LOC(ii,2),reg%TOPO(ii,1),&
            &reg%TOPO(ii,2),reg%CELL(ii)%temperature&
            &,reg%CELL(ii)%vs_obs, reg%CELL(ii)&
            &%vs_sol, reg%CELL(ii)%melt_fraction,&
            &reg%CELL(ii)%K_sol,reg%CELL(ii)%G_sol,&
            &reg%CELL(ii)%vs_obs/reg%CELL(ii)%vs_sol-&
            &reg%CELL(ii)%vs_eff,reg%CELL(ii)%amplitude

    END DO
    CLOSE(10)
  END SUBROUTINE data_table_output


  SUBROUTINE melt_gmt_output(reg,mod,dir)
    IMPLICIT NONE
    TYPE(region),INTENT(in)::reg
    TYPE(mantle),INTENT(in)::mod
    CHARACTER(len=11),intent(in)::dir
    CHARACTER(len=26)::fname
    INTEGER::ii,N,NNODES
    fname=dir//filename(reg%MAGMA%dihedral,reg%comp%basalt_fraction,&
         &reg%comp%potential_temperature,'.dat')
    N=reg%NGRID
    OPEN(10,file=fname,form='formatted')

    DO ii=1,N
       WRITE(10,'(3(1x,f12.5))') reg%LOC(ii,1),&
            &reg%LOC(ii,2), reg%CELL(ii)%melt_fraction
    END DO
    CLOSE(10)
  END SUBROUTINE melt_gmt_output

  FUNCTION impedance_contrast_410(int)
    IMPLICIT NONE
    INTEGER,INTENT(in)::int
    REAL(sp)::impedance_contrast_410
    REAL(sp)::vstop,vsbot,rhotop,rhobot
    rhotop=3543.25_sp
    rhobot=3723.78_sp
    IF (int==1) THEN  
       !PREM model
       vstop=4769.89_sp
       vsbot=4932.59_sp
    ELSE          
       ! IASP
       vstop=4870.0_sp
       vsbot=5070.0_sp
    END IF
    impedance_contrast_410=(rhobot*vsbot-rhotop*vstop)/rhotop/vstop
  END FUNCTION impedance_contrast_410

  FUNCTION velocity_contrast_410(int)
    IMPLICIT NONE
    INTEGER,INTENT(in)::int
    REAL(sp)::velocity_contrast_410
    REAL(sp)::vstop,vsbot,rhotop,rhobot
    rhotop=3543.25_sp
    rhobot=3723.78_sp
    IF (int==1) THEN  
       !PREM model
       vstop=4769.89_sp
       vsbot=4932.59_sp

    ELSE          

       ! IASP
       vstop=4870.0_sp
       vsbot=5070.0_sp
    END IF

    velocity_contrast_410=(vsbot-vstop)/vstop

  END FUNCTION velocity_contrast_410



END MODULE regional
