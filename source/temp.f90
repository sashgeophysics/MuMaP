program temp

  IMPLICIT NONE
  INTEGER, PARAMETER :: SP = KIND(1.0D0)
  CHARACTER(len=40)::fname
  character(len=24)::colnum,fmt
  REAL(sp)::z,impedance350,prem_impedance
  REAL(sp),DIMENSION(2)::impedance_vel
  REAL(sp),DIMENSION(:,:),POINTER::DATA,temp1
  REAL(sp),DIMENSION(:,:),POINTER::x
  REAL(sp),DIMENSION(3)::prem_value
  INTEGER::ii,ncols,NGRID
  fname='../input/&
       &Hawaii_temperature_agius.csv'
  NGRID=line_numbers(fname,10000)
  NGRID=1681
  ncols=14
  print *, 'Nrid is',NGRID
  WRITE(colnum,'(i3)'),ncols
  fmt="("//trim(colnum)//"(1x,e12.5))"
  ALLOCATE(DATA(NGRID,ncols),&
       &temp1(ncols,NGRID),&
       &x(NGRID,3))
  OPEN(10,file=fname,status='unknown',form='formatted')    
  READ(10,fmt),temp1

  CLOSE(10)
  DATA=transpose(temp1)
  x(:,1)=DATA(:,1)
  x(:,2)=DATA(:,10)
  x(:,3)=DATA(:,10)+0.3_sp*350.0_sp
  write(*,fmt),x(:,2)
  print *, '#######################'
  write(*,fmt),x(:,3)
  print*,'load_seismo_hmt: max(amp),min(amp)',maxval(DATA(:,5)),minval(DATA(:,5))

contains
  FUNCTION line_numbers(fname,lmax)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)::fname
    INTEGER,INTENT(in)::lmax
    INTEGER::line_numbers,ii,lineno,ios
    REAL::junk
    lineno=0
    OPEN(10,file=fname,form='formatted')
    DO ii=1,lmax
       READ(10,*,IOSTAT=ios), junk
       IF(ios/=0) EXIT
       lineno=lineno+1
    END DO
    CLOSE(10)
    !PRINT*,'linenumber:',lineno
    line_numbers=lineno
  END FUNCTION line_numbers
end program temp
