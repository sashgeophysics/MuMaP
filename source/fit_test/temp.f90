program temp
  implicit none

  INTEGER, PARAMETER :: SP = KIND(1.0D0)
  INTEGER,dimension(22)::basalt_fraction
  real(sp),DIMENSION(9)::f,temperaturefit
  INTEGER::ii,n,jj
  REAL(sp),DIMENSION(22,10)::par
  REAL(sp)::X
 
  !> This function interpolates VS, Vp, and density
  !! as a function of potential temperature for a given
  !! basalt fraction. The fut us determined by a fifthe order
  !! polynomial. Given the basalt fraction, the function
  !! returns 3 sets of six coefficients, for Vs,Vp, and rho.
  !! The fits are done on the data of Xu et al !<


  ! Enter coefficients for polynomial fit for vs   
  basalt_fraction=(/0,5,10,15,18,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100/)
  OPEN(10,file='../input/xu_fit_param.csv',form='formatted')
  
  n=22
  DO ii=1,22
     READ(10,'(10(1x,es10.3))'),par(ii,:)
  END DO
  CLOSE(10)
  do jj=1,22
     X=DBLE(basalt_fraction(jj))/100.0_sp
  DO ii=1,22
     IF(X==par(ii,1)) THEN
        !Density
        f(1)=par(ii,2)
        f(2)=par(ii,3)
        f(3)=par(ii,4)
        !Vs
        f(4)=par(ii,5)
        f(5)=par(ii,6)
        f(6)=par(ii,7)
        !Vp
        f(7)=par(ii,8)
        f(8)=par(ii,9)
        f(9)=par(ii,10)
     END IF
  END DO
  print*,'basalt_fraction',basalt_fraction(jj)
     !write(*,'(9(1x,es10.3))'),f
  end do
end program temp
