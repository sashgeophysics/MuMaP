MODULE global
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

  pure function eye(n)
    integer(I4B),intent(in)::n
    real(sp),dimension(n,n)::eye
    integer(I4B)::ii,jj
    eye=0.0_sp
    do ii=1,n
       do jj=1,n
          if(ii==jj)then
             eye(ii,jj)=1.0_sp
          end if
       end do
    end do
  end function eye

  FUNCTION lusol(KPACK,FPACK)

    IMPLICIT NONE
    REAL(sp),DIMENSION(:,:),INTENT(in)::KPACK
    REAL(sp),DIMENSION(:),INTENT(in)::FPACK
    REAL(sp),DIMENSION(SIZE(FPACK,1))::lusol
    real(SP),dimension(size(KPACK,1),size(KPACK,2))::A
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
    print*,' lusol:error',sum(abs(matmul(KPACK,lusol)-FPACK))

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

  FUNCTION gmres(KPACK,FPACK)

    IMPLICIT NONE
    REAL(sp),DIMENSION(:,:),INTENT(in)::KPACK
    REAL(sp),DIMENSION(:),INTENT(in)::FPACK
    REAL(sp),DIMENSION(SIZE(FPACK,1))::gmres,rhs,soln
    type(sparse)::sprs
    integer::SIZ
    SIZ=128
    rhs=FPACK
    soln=0.0
    sprs=dnscsr(KPACK,1.0e-6_sp)
    call fgmres_nopc(sprs%IA,sprs%JA,sprs%A,rhs,soln,SIZE(FPACK,1),SIZ)
    gmres=soln
  end FUNCTION gmres
  
  SUBROUTINE fgmres_nopc(IA,JA,A,RHS,COMPUTED_SOLUTION,N,SIZ)

    INTEGER,DIMENSION(:),INTENT(inout)::IA,JA
    REAL(sp),DIMENSION(:),INTENT(inout)::A,COMPUTED_SOLUTION,RHS
    INTEGER,INTENT(in)::N,SIZ
    INTEGER, PARAMETER::EXPECTED_ITERCOUNT=5

    INTEGER::IPAR(SIZ),ITERCOUNT,RCI_REQUEST, I
    REAL(sp)::DPAR(SIZ),TMP(N*(2*N+1)+(N*(N+9))/2+1)
    REAL(sp):: DVAR,EXPECTED_SOLUTION(N)

    !---------------------------------------------------------------------------
    ! Initialize variables and the right hand side through matrix-vector product
    !---------------------------------------------------------------------------
    CALL MKL_DCSRGEMV('N', N, A, IA, JA, EXPECTED_SOLUTION, RHS)

    !---------------------------------------------------------------------------
    ! Initialize the initial guess
    !---------------------------------------------------------------------------
    DO I=1,N
       COMPUTED_SOLUTION(I)=1.0_sp
    ENDDO
    !---------------------------------------------------------------------------
    ! Initialize the solver
    !---------------------------------------------------------------------------
    CALL DFGMRES_INIT(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST, IPAR, DPAR, TMP)
    !IF (RCI_REQUEST.NE.0) GOTO 999
    !---------------------------------------------------------------------------
    ! Set the desired parameters:
    ! LOGICAL parameters:
    ! do residual stopping test
    ! do not request for the user defined stopping test
    ! do the check of the norm of the next generated vector automatically
    ! DOUBLE PRECISION parameters
    ! set the relative tolerance to 1.0D-3 instead of default value 1.0D-6
    !---------------------------------------------------------------------------
    IPAR(9)=1
    IPAR(10)=0
    IPAR(12)=1
    DPAR(1)=1.0E-3_sp
    !---------------------------------------------------------------------------
    ! Check the correctness and consistency of the newly set parameters
    !---------------------------------------------------------------------------
    CALL DFGMRES_CHECK(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST, IPAR, DPAR, TMP)

    ! Reverse Communication starts here
    !---------------------------------------------------------------------------
1   CALL DFGMRES(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST, IPAR, DPAR, TMP)
    !---------------------------------------------------------------------------
    ! If RCI_REQUEST=0, then the solution was found with the required precision
    !---------------------------------------------------------------------------
    IF (RCI_REQUEST.EQ.0) GOTO 3
    !---------------------------------------------------------------------------
    ! If RCI_REQUEST=1, then compute the vector A*TMP(IPAR(22))
    ! and put the result in vector TMP(IPAR(23))
    !---------------------------------------------------------------------------
    IF (RCI_REQUEST.EQ.1) THEN
       CALL MKL_DCSRGEMV('N',N, A, IA, JA, TMP(IPAR(22)), TMP(IPAR(23)))
       GOTO 1
       !---------------------------------------------------------------------------
       ! If RCI_REQUEST=anything else, then DFGMRES subroutine failed
       ! to compute the solution vector: COMPUTED_SOLUTION(N)
       !---------------------------------------------------------------------------
    ELSE
       GOTO 999
    ENDIF
    !---------------------------------------------------------------------------
    ! Reverse Communication ends here
    ! Get the current iteration number and the FGMRES solution (DO NOT FORGET to
    ! call DFGMRES_GET routine as COMPUTED_SOLUTION is still containing
    ! the initial guess!)
    !---------------------------------------------------------------------------
3   CALL DFGMRES_GET(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST, IPAR, DPAR, TMP, ITERCOUNT)
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Print solution vector: COMPUTED_SOLUTION(N) and
    ! the number of iterations: ITERCOUNT
    !---------------------------------------------------------------------------
    

    !---------------------------------------------------------------------------
    ! Release internal MKL memory that might be used for computations
    ! NOTE: It is important to call the routine below to avoid memory leaks
    ! unless you disable MKL Memory Manager
    !---------------------------------------------------------------------------
    CALL MKL_FREEBUFFERS

    DVAR=DNRM2(N,EXPECTED_SOLUTION,1)

    !---------------------------------------------------------------------------
    ! Release internal MKL memory that might be used for computations
    ! NOTE: It is important to call the routine below to avoid memory leaks
    ! unless you disable MKL Memory Manager
    !---------------------------------------------------------------------------
999 WRITE( *,'(A,A,I5)') 'This example FAILED as the solver has',&
         & ' returned the ERROR code', RCI_REQUEST
    CALL MKL_FREEBUFFERS
    ! STOP 1
  END SUBROUTINE fgmres_nopc

FUNCTION dnscsr(inmat,tolin)
REAL(sp),DIMENSION(:,:),INTENT(in)::inmat
real(sp),optional::tolin
TYPE(sparse)::dnscsr
integer::Nz,m,n,job(6),info
REAL(sp)::tol
REAL(sp),dimension(size(inmat,1),size(inmat,2))::mat
!! set the tolerance for zero search
IF(PRESENT(tolin)) THEN
     tol=tolin
  ELSE
     tol=1.0e-8_sp
  END IF
mat=inmat
!Zero out the entries that are below tolerance level
where(abs(inmat)<tol)
mat=0.0_sp
end where
! Count the number of nonzero entries in the incoming matrix
Nz=sparscalc(mat,tol)
!Number of rows in the incoming matrix
m=size(mat,1)
!Number of columns in the incoming matrix
n=size(mat,2)
! Allocate the arrays for sparse storage
ALLOCATE(dnscsr%A(Nz),dnscsr%IA(m+1),dnscsr%JA(Nz))

! Set the job preference
job(1)=0 ! Convert regular matrix to CSR
job(2)=1 ! Use one-based matrix notation
job(3)=1 ! Use one-based indexing in CSR format
job(4)=2 ! Input the whole matrix
job(5)=Nz! Maximum number of nonzero elements
job(6)=1 ! Return Ia, Ja, and AA

call mkl_ddnscsr(job,m,n,mat,m,dnscsr%A,dnscsr%JA,dnscsr%IA,info)

END FUNCTION dnscsr

FUNCTION sparscalc(mat,tolin)
  ! This function calculates the number of nonzero entries in 
  ! a sparse matrix. By default, numbers smaller than 1.0e-8 
  ! are considered zero, but the default can be overriden by
  ! supplying the optional second paramter
  
  implicit none
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

  FUNCTION filename(f,str)
    INTEGER(I4B),INTENT(in)::f
    character(len=3),intent(in)::str
    CHARACTER(len=40)::frmt
    CHARACTER(len=18)::filename
    
    IF(f<10) THEN
       frmt='("../data/",a3,"00",i1,".vtk")'
       WRITE(filename,frmt),str,f
    ELSEIF((f<100).AND.(f>=10)) THEN
       frmt='("../data/",a3,"0",i2,".vtk")'
       WRITE(filename,frmt),str,f
    ELSEIF((f<1000).AND.(f>=100)) THEN
       frmt='("../data/",a3,i3,".vtk")'
       WRITE(filename,frmt),str,f
    END IF
  END FUNCTION filename

 
END MODULE global
