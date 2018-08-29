PROGRAM LAB
  USE global
  USE microgeodynamics
  USE regional
  !TYPE(unit_cell)::PREM
  TYPE(region)::HAWAII
  TYPE(composition)::HAWAII_REF_COMP
  TYPE(mantle)::HAWAII_REF_MANTLE
  INTEGER::ii,jj,kk,filesize
  INTEGER :: clock_start,clock_end,clock_rate
  REAL(sp) :: elapsed_time,theta,dihedral_angle
  REAL(sp),DIMENSION(3)::prem,temp
  REAL(sp),DIMENSION(:,:),POINTER::melting_column

  CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
  CALL SYSTEM_CLOCK(COUNT=clock_start) ! Start timing

  ! Define composition space 
  HAWAII_REF_COMP=composition(potential_temperature&
       &=1600,basalt_fraction=18)
  ! Load solid composition for 352 km from Xu et al (2008)
  HAWAII_REF_MANTLE=load_composition(HAWAII_REF_COMP&
       &,'../input/LAB_1600.csv')          
  dihedral_angle=5.0_sp
  ! Set the melt composition
  HAWAII%MAGMA=set_melt(5,dihedral_angle)
  
  filesize=line_numbers('../input/LAB_1600.csv',10000)
  ALLOCATE(melting_column(filesize,4))
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1600.csv','../data/LAB_T1600_d05_phi001.dat',1)
  melting_column=>forward_melting(0.02_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1600.csv','../data/LAB_T1600_d05_phi02.dat',1)
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1500.csv','../data/LAB_T1500_d05_phi001.dat',1)
  melting_column=>forward_melting(0.02_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1500.csv','../data/LAB_T1500_d05_phi02.dat',1)
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1300.csv','../data/LAB_T1300_d05_phi001.dat',1)
  melting_column=>forward_melting(0.02_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1300.csv','../data/LAB_T1300_d05_phi02.dat',1)
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1200.csv','../data/LAB_T1200_d05_phi001.dat',1)
  melting_column=>forward_melting(0.02_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1200.csv','../data/LAB_T1200_d05_phi02.dat',1)
  dihedral_angle=25.0_sp
  HAWAII%MAGMA=set_melt(5,dihedral_angle)
  
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1600.csv','../data/LAB_T1600_d25_phi001.dat',1)
  melting_column=>forward_melting(0.02_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1600.csv','../data/LAB_T1600_d25_phi02.dat',1)
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1500.csv','../data/LAB_T1500_d25_phi001.dat',1)
  melting_column=>forward_melting(0.02_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1500.csv','../data/LAB_T1500_d25_phi02.dat',1)
   melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1300.csv','../data/LAB_T1300_d25_phi001.dat',1)
  melting_column=>forward_melting(0.02_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1300.csv','../data/LAB_T1300_d25_phi02.dat',1)
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1200.csv','../data/LAB_T1200_d25_phi001.dat',1)
  melting_column=>forward_melting(0.02_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1200.csv','../data/LAB_T1200_d25_phi02.dat',1)
  ! Now for melt films
  HAWAII%MAGMA%aspect=0.005_sp
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1600.csv','../data/LAB_T1600_005_phi001.dat',2)
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1500.csv','../data/LAB_T1500_005_phi001.dat',2)
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1300.csv','../data/LAB_T1300_005_phi001.dat',2)
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1200.csv','../data/LAB_T1200_005_phi001.dat',2)
  HAWAII%MAGMA%aspect=0.002_sp
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1600.csv','../data/LAB_T1600_002_phi001.dat',2)
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1500.csv','../data/LAB_T1500_002_phi001.dat',2)
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1300.csv','../data/LAB_T1300_002_phi001.dat',2)
  melting_column=>forward_melting(0.001_sp,HAWAII%MAGMA,HAWAII_REF_MANTLE,&
       &'../input/LAB_1200.csv','../data/LAB_T1200_002_phi001.dat',2)
  melting_column=>null()
  CALL SYSTEM_CLOCK(COUNT=clock_end) ! Stop timing
  ! Calculate the elapsed time in seconds:
  elapsed_time=REAL((clock_end-clock_start)/clock_rate,sp)


  PRINT*,'elapsed_time',elapsed_time
  PRINT*,'# of iterations',(ii-1)*(jj-1)*(kk-1)

contains

END PROGRAM LAB
