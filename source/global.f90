MODULE global
  !USE LAPACK95
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: DP = KIND(1.0)
  INTEGER, PARAMETER :: SP = KIND(1.0D0)
  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
  TYPE sparse
     REAL(sp),DIMENSION(:),POINTER::A
     INTEGER,DIMENSION(:),POINTER::IA,JA
     INTEGER::m,n
  END TYPE sparse
  INTERFACE linspace
     MODULE PROCEDURE linspace_r, linspace_i
  END INTERFACE linspace

CONTAINS
  PURE FUNCTION linspace_r(a,b,n)
    IMPLICIT NONE
    REAL(sp),INTENT(in)::a,b
    INTEGER(I4B),INTENT(in)::n
    INTEGER(I4B)::ii
    REAL(sp),DIMENSION(n)::linspace_r
    DO ii=1,n
       linspace_r(ii)=a+(REAL(ii,kind=sp)-1.0_sp)/(REAL(n,kind=sp)-1.0_sp)*(b-a)
    END DO
  END FUNCTION linspace_r
  PURE FUNCTION linspace_i(a,b,n)
    IMPLICIT NONE
    INTEGER(I4B),INTENT(in)::a,b
    INTEGER(I4B),INTENT(in)::n
    INTEGER(I4B)::ii
    INTEGER(I4B),DIMENSION(n)::linspace_i
    DO ii=1,n
       linspace_i(ii)=a+(ii-1)*(b-a)/(n-1)
    END DO
  END FUNCTION linspace_i

  PURE FUNCTION eye(n)
    INTEGER(I4B),INTENT(in)::n
    REAL(sp),DIMENSION(n,n)::eye
    INTEGER(I4B)::ii,jj
    eye=0.0_sp
    DO ii=1,n
       DO jj=1,n
          IF(ii==jj)THEN
             eye(ii,jj)=1.0_sp
          END IF
       END DO
    END DO
  END FUNCTION eye

  FUNCTION average(array)
    IMPLICIT NONE
    REAL(sp),DIMENSION(:),INTENT(in)::array
    REAL(sp)::average
    INTEGER::int
    int=SIZE(array)
    average=SUM(array)/REAL(int,kind=sp)

  END FUNCTION average
  FUNCTION lusol(KPACK,FPACK)

    IMPLICIT NONE
    REAL(sp),DIMENSION(:,:),INTENT(in)::KPACK
    REAL(sp),DIMENSION(:),INTENT(in)::FPACK
    REAL(sp),DIMENSION(SIZE(FPACK,1))::lusol
    REAL(SP),DIMENSION(SIZE(KPACK,1),SIZE(KPACK,2))::A
    !This is a wrapper for solving a general system of linear
    !equations using the intel mkl library
    !KPACK*lusol=FPACK
    INTEGER(I4B)::n,lda,info,ldb,nrhs
    INTEGER(I4B),DIMENSION(SIZE(KPACK,1))::ipiv
    CHARACTER(len=1)::trans
    n=SIZE(KPACK,1)
    lda=n
    lusol=FPACK   
    A=KPACK
    !LU factorize the matrix first
    CALL dgetrf(n,n,A,lda,ipiv,info)    
    trans='N'
    nrhs=1
    ldb=n
    !Then solve   
    CALL dgetrs(trans,n,nrhs,A,lda,ipiv,lusol,ldb,info)    
    PRINT*,' lusol:error',SUM(ABS(MATMUL(KPACK,lusol)-FPACK))

  END FUNCTION lusol


  FUNCTION tridag(a,b,c,r)
    IMPLICIT NONE
    REAL(sp),DIMENSION(:),INTENT(in)::a,b,c,r
    REAL(sp),DIMENSION(SIZE(b,1)-1)::du2
    REAL(sp),DIMENSION(SIZE(r,1))::tridag,sol
    INTEGER(I4B)::n,lda,info,ldb,nrhs
    INTEGER(I4B),DIMENSION(SIZE(b,1))::ipiv
    CHARACTER(len=1)::trans
    n=SIZE(b,1)
    sol=r
    !LU factorize the matrix first
    CALL dgttrf(n,a,b,c,du2,ipiv,info)
    trans='N'
    nrhs=1
    ldb=n
    !Then solve   

    CALL dgttrs(trans,n,nrhs,a,b,c,du2,ipiv,sol,ldb,info)
    tridag=sol
  END FUNCTION tridag



  FUNCTION dnscsr(inmat,tolin)
    REAL(sp),DIMENSION(:,:),INTENT(in)::inmat
    REAL(sp),OPTIONAL::tolin
    TYPE(sparse)::dnscsr
    INTEGER::Nz,m,n,job(6),info
    REAL(sp)::tol
    REAL(sp),DIMENSION(SIZE(inmat,1),SIZE(inmat,2))::mat
    !! set the tolerance for zero search
    IF(PRESENT(tolin)) THEN
       tol=tolin
    ELSE
       tol=1.0e-8_sp
    END IF
    mat=inmat
    !Zero out the entries that are below tolerance level
    WHERE(ABS(inmat)<tol)
       mat=0.0_sp
    END WHERE
    ! Count the number of nonzero entries in the incoming matrix
    Nz=sparscalc(mat,tol)
    !Number of rows in the incoming matrix
    m=SIZE(mat,1)
    !Number of columns in the incoming matrix
    n=SIZE(mat,2)
    ! Allocate the arrays for sparse storage
    ALLOCATE(dnscsr%A(Nz),dnscsr%IA(m+1),dnscsr%JA(Nz))

    ! Set the job preference
    job(1)=0 ! Convert regular matrix to CSR
    job(2)=1 ! Use one-based matrix notation
    job(3)=1 ! Use one-based indexing in CSR format
    job(4)=2 ! Input the whole matrix
    job(5)=Nz! Maximum number of nonzero elements
    job(6)=1 ! Return Ia, Ja, and AA

    CALL mkl_ddnscsr(job,m,n,mat,m,dnscsr%A,dnscsr%JA,dnscsr%IA,info)

  END FUNCTION dnscsr

  FUNCTION sparscalc(mat,tolin)
    ! This function calculates the number of nonzero entries in 
    ! a sparse matrix. By default, numbers smaller than 1.0e-8 
    ! are considered zero, but the default can be overriden by
    ! supplying the optional second paramter

    IMPLICIT NONE
    REAL(sp),DIMENSION(:,:),INTENT(in)::mat
    REAL(sp),OPTIONAL::tolin
    REAL(sp)::tol
    INTEGER::nn,sparscalc
    INTEGER,DIMENSION(SIZE(mat,1),SIZE(mat,2))::tmp
    tmp=0
    IF(PRESENT(tolin)) THEN
       tol=tolin
    ELSE
       tol=1.0e-8_sp
    END IF

    WHERE(ABS(mat)>tol)
       tmp=1
    END WHERE
    sparscalc=SUM(tmp)
  END FUNCTION sparscalc

  FUNCTION inverse(MAT)
    IMPLICIT NONE
    REAL(sp),DIMENSION(:,:),INTENT(in)::MAT
    REAL(SP),DIMENSION(SIZE(MAT,1),SIZE(MAT,2))::A,inverse
    !This is a wrapper for inverting a matrix 
    INTEGER(I4B)::n,lda,info,lwork
    INTEGER(I4B),DIMENSION(SIZE(MAT,1))::ipiv
    REAL(sp),DIMENSION(:),ALLOCATABLE::work
    n=SIZE(MAT,1)
    lda=n
    inverse=MAT
    A=MAT
    !LU factorize the matrix first
    CALL dgetrf(n,n,A,lda,ipiv,info)    
    lwork=n*64
    ALLOCATE(work(lwork))
    !Then invert  
    CALL dgetri( n, A, lda, ipiv, work, lwork, info )
    DEALLOCATE(work)
    inverse=A
    PRINT*,'inverse: maximum error',MAXVAL(eye(n)-MATMUL(MAT,A))

  END FUNCTION inverse

  SUBROUTINE eigen_val(A,wr,wi)
    IMPLICIT NONE
    REAL(sp),DIMENSION(:,:),INTENT(in)::A
    REAL(sp),DIMENSION(SIZE(A,1)),INTENT(out)::wr,wi
    REAL(sp),DIMENSION(SIZE(A,1))::scale,work,w
    REAL(sp),DIMENSION(SIZE(A,1)-1)::tau
    REAL(SP),DIMENSION(SIZE(A,1),SIZE(A,2))::z,H
    INTEGER::ilo,ihi,info,n,lda,ldz,ldh,lwork
    CHARACTER(len=1)::job,compz
    job='B'
    n=SIZE(A,1)
    lda=n
    lwork=n
    ldz=n
    ldh=n
    H=A
    ! Balance the incoming matrix
    CALL dgebal(job, n, H, lda, ilo, ihi, scale, info)
    ! Reduce the matrix to upper-Hysenberg form
    CALL dgehrd(n, ilo, ihi, H, lda, tau, work, lwork, info)
    ! Calculate eigenvalues
    job='E'   ! Only eigenvalues needed
    compz='N' ! No Schur vectors are needed
    z=0.0_sp  ! It is not referenced when compz='N', but just set to zero
    CALL dhseqr(job, compz, n, ilo, ihi, H, ldh, wr, wi, z, ldz, work, lwork, info)
    !hseqr(H, w, ilo,ihi,z,job,compz,info)
    ! In the output w is the array of eigenvalues calculated by the routine

  END SUBROUTINE eigen_val

  FUNCTION eigentest(k)
    IMPLICIT NONE
    INTEGER,INTENT(in)::k
    INTEGER::eigentest
    REAL(sp),DIMENSION(8,8)::mat
    REAL(sp),DIMENSION(8)::analytical,wr,wi,imag,rel


    INTEGER::ii
    mat(1,:)=(/611.0_sp,196.0_sp,-192.0_sp,407.0_sp,-8.0_sp,-52.0_sp,-49.0_sp,29.0_sp/)
    mat(2,:)=(/196.0_sp,899.0_sp,113.0_sp,-192.0_sp,-71.0_sp,-43.0_sp,-8.0_sp,-44.0_sp/)
    mat(3,:)=(/-192.0_sp,113.0_sp,899.0_sp,196.0_sp,61.0_sp,49.0_sp,8.0_sp,52.0_sp/)
    mat(4,:)=(/407.0_sp,-192.0_sp,196.0_sp,611.0_sp,8.0_sp,44.0_sp,59.0_sp,-23.0_sp/)
    mat(5,:)=(/-8.0_sp,-71.0_sp,61.0_sp,8.0_sp,411.0_sp,-599.0_sp,208.0_sp,208.0_sp/)
    mat(6,:)=(/-52.0_sp,-43.0_sp,49.0_sp,44.0_sp,-599.0_sp,411.0_sp,208.0_sp,208.0_sp/)
    mat(7,:)=(/-49.0_sp,-8.0_sp,8.0_sp,59.0_sp,208.0_sp,208.0_sp,99.0_sp,-911.0_sp/)
    mat(8,:)=(/29.0_sp,-44.0_sp,52.0_sp,-23.0_sp,208.0_sp,208.0_sp,-911.0_sp,99.0_sp/)

    PRINT*,'eigentest: input matrix'
    DO ii=1,8
       WRITE(*,'(8(1x,f11.4))'),mat(ii,:)
    END DO

    CALL eigen_val(mat,wr,wi)


    analytical=(/ 0.0_sp, 1020.0_sp, 510.0_sp-100.0_sp*SQRT(26.0_sp), &
         & 510.0_sp+100.0_sp*SQRT(26.0_sp), 10.0_sp*SQRT(10405.0_sp), &
         &-10.0_sp*SQRT(10405.0_sp),  1000.0_sp, 1000.0_sp/)


    PRINT*,'Calculated eigenvalues: analytical, real, imaginary'
    DO ii=1,8
       WRITE(*,'(3(1x,f11.4))'),analytical(ii),wr(ii),wi(ii)
    END DO

    mat(1,:)=(/10.0_sp,-65.0_sp,74.0_sp,42.0_sp,-22.0_sp,45.0_sp,51.0_sp,4.0_sp/)
    mat(2,:)=(/25.0_sp,-69.0_sp,-18.0_sp,-18.0_sp,-9.0_sp,24.0_sp,-1.0_sp,79.0_sp/)
    mat(3,:)=(/21.0_sp,-90.0_sp,45.0_sp,-64.0_sp,67.0_sp,-35.0_sp,49.0_sp,-42.0_sp/)
    mat(4,:)=(/97.0_sp,80.0_sp,0.0_sp,74.0_sp,-77.0_sp,65.0_sp,86.0_sp,-72.0_sp/)
    mat(5,:)=(/-79.0_sp,-97.0_sp,-65.0_sp,2.0_sp,0.0_sp,41.0_sp,47.0_sp,-5.0_sp/)
    mat(6,:)=(/77.0_sp,-86.0_sp,69.0_sp,32.0_sp,-50.0_sp,0.0_sp,79.0_sp,95.0_sp/)
    mat(7,:)=(/65.0_sp,-49.0_sp,-46.0_sp,-32.0_sp,-69.0_sp,88.0_sp,35.0_sp,-82.0_sp/)
    mat(8,:)=(/3.0_sp,27.0_sp,-25.0_sp,3.0_sp,-14.0_sp,42.0_sp,-38.0_sp,-3.0_sp/)

    PRINT*,'eigentest: input matrix'
    DO ii=1,8
       WRITE(*,'(8(1x,f11.4))'),mat(ii,:)
    END DO

    CALL eigen_val(mat,wr,wi)

    rel=(/140.01_sp,140.01_sp,-127.04_sp,-50.56_sp,&
         &-50.56_sp,19.66_sp,19.66_sp,0.8_sp/)

    imag=(/55.68_sp,-55.68_sp,0.0_sp,48.15_sp,&
         &-48.15_sp,69.28_sp,-69.28_sp,0.0_sp/)

    PRINT*,'Eigenvalues: real(analytical), imaginary(analytical), real, imaginary'
    DO ii=1,8
       WRITE(*,'(4(1x,f11.4))'),rel(ii),imag(ii),wr(ii),wi(ii)
    END DO


    eigentest=1
  END FUNCTION eigentest

  FUNCTION filename(dr,b,T,str)
    real(sp),intent(in)::dr
    INTEGER(I4B),INTENT(in)::b,T
    integer::d
    CHARACTER(len=4),INTENT(in)::str
    CHARACTER(len=3)::dc,bc
    CHARACTER(len=5)::tc
    CHARACTER(len=40)::frmt
    CHARACTER(len=15)::filename

    d=aint(dr)
    IF(T<10) THEN
       frmt='("T000",i1)'
       WRITE(tc,frmt),T
    ELSEIF((T<100).AND.(T>=10)) THEN
       frmt='("T00",i2)'
       WRITE(tc,frmt),T       
    ELSEIF((T<1000).AND.(T>=100)) THEN
       frmt='("T0",i3)'
       WRITE(tc,frmt),T
    ELSE
       frmt='("T",i4)'
       WRITE(tc,frmt),T
    END IF

    IF(b<10) THEN
       frmt='("b0",i1)'
       WRITE(bc,frmt),b
    ELSE
       frmt='("b",i2)'
       WRITE(bc,frmt),b      
    END IF

    IF(d<10) THEN
       frmt='("d0",i1)'
       WRITE(dc,frmt),d
    ELSE
       frmt='("d",i2)'
       WRITE(dc,frmt),d      
    END IF

    frmt='(a,a,a,a)'
    WRITE(filename,frmt),dc,bc,tc,str


  END FUNCTION filename

  

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

END MODULE global
