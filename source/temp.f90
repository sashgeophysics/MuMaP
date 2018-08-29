program temp
  implicit none
  INTEGER::dim
  CHARACTER(len=72)::fname,fname2
  character(len=24)::colnum,fmt
  INTEGER::ncols,NGRID ! The number ofcolumns in the input file

  REAL,DIMENSION(:,:),POINTER::DATA,temp1
  REAL,DIMENSION(:,:),POINTER::x
  REAL,DIMENSION(3)::prem_value

  INTEGER::ii
  
  fname='../input/mydata.csv'
  fname2='hmtdata.csv'
  NGRID=line_numbers(fname,10000)
  ncols=14
  WRITE(colnum,'(i3)'),ncols
  
  fmt="("//trim(colnum)//"(1x,f10.3))"
  ALLOCATE(DATA(NGRID,ncols),temp1(ncols,NGRID),x(3,NGRID))
  OPEN(10,file=fname,status='unknown',form='formatted')
  OPEN(20,file=fname2,status='unknown',form='formatted')
  OPEN(20,file='x.dat',status='unknown',form='formatted')
  
  READ(10,fmt),temp1
  DATA=transpose(temp1)
  DO ii=1,NGRID
     write(20,fmt),DATA(ii,:)
     write(*,fmt),DATA(ii,:)
  END DO
  CLOSE(10)
  close(20)
  deallocate(DATA,x)
contains
  !> This is a utility function to count the number of lines in a file
  !On input, filename is the name of the file and lmax is a guess
  ! about how many lines there are. If you are not sure, provide a
  ! large number, say, 10000. !<
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
