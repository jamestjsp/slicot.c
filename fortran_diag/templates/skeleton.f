C
C SPDX-License-Identifier: BSD-3-Clause
C
C {{ROUTINE}} DIAGNOSTIC PROGRAM - FORTRAN VERSION
C
C Instrumented version of T{{ROUTINE}} that prints all intermediate
C values for comparison with C implementation.
C
C DEVELOPER TODO:
C 1. Update parameters (N, M, etc.) based on routine signature
C 2. Update READ statements to match data file format
C 3. Update CALL statement with correct routine signature
C 4. Add WRITE statements for output matrices/vectors
C
C Refer to: SLICOT-Reference/examples/T{{ROUTINE}}.f
C
C     .. Parameters ..
      INTEGER           NIN, NOUT
      PARAMETER         ( NIN = 5, NOUT = 6 )
      INTEGER           NMAX
      PARAMETER         ( NMAX = 20 )
C
C     TODO: Update leading dimensions based on routine requirements
      INTEGER           LDA
      PARAMETER         ( LDA = NMAX )
C
C     TODO: Update LDWORK based on routine requirements
      INTEGER           LDWORK
      PARAMETER         ( LDWORK = MAX( 1, 4*NMAX ) )
C
C     .. Local Scalars ..
C     TODO: Add routine-specific parameters (e.g., DICO, FACT, TRANS)
      INTEGER           I, INFO, J, N
C
C     .. Local Arrays ..
C     TODO: Add routine-specific arrays (A, B, Q, etc.)
      DOUBLE PRECISION  DWORK(LDWORK)
C
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C
C     .. External Subroutines ..
      EXTERNAL          {{ROUTINE}}
C
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C
C     .. Executable Statements ..
C
      WRITE ( NOUT, FMT = '(A)' )
     $  '=== {{ROUTINE}} FORTRAN DIAGNOSTIC ==='
      WRITE ( NOUT, FMT = '(A)' ) ''
C
C     Skip the heading in the data file and read the data.
      READ ( NIN, FMT = '()' )
C
C     TODO: Read parameters from data file
C     Example: READ ( NIN, FMT = * ) N, M, DICO, FACT, TRANS
C     READ ( NIN, FMT = * ) <parameters>
C
C     TODO: Print parameters
C     Example: WRITE ( NOUT, FMT = '(A,I3)' ) 'N = ', N
C
C     TODO: Validate parameter ranges
C     IF ( N.LT.0 .OR. N.GT.NMAX ) THEN
C        WRITE ( NOUT, FMT = '(A,I5)' ) 'ERROR: N out of range: ', N
C        STOP
C     END IF
C
C     TODO: Read input matrices
C     Example: READ ( NIN, FMT = * ) ( ( A(I,J), J = 1,N ), I = 1,N )
C     Note: Pay attention to loop order (row-wise vs column-wise)
C     Note: Handle conditional reads (IF LSAME(FLAG, 'X'))
C
C     TODO: Print input matrices with high precision
C     Example:
C     WRITE ( NOUT, FMT = '(/,A)' ) 'INPUT A matrix:'
C     DO 10 I = 1, N
C        WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
C    $        ( A(I,J), J = 1,N )
C  10 CONTINUE
C
C     TODO: Call the routine
C     Example: CALL {{ROUTINE}}( N, M, A, LDA, ..., DWORK, LDWORK, INFO )
C
C     TODO: Print return status
C     WRITE ( NOUT, FMT = '(/,A,I5)' ) 'INFO = ', INFO
C
C     TODO: Print output matrices/vectors
C     Example:
C     WRITE ( NOUT, FMT = '(/,A)' ) 'OUTPUT A matrix:'
C     DO 20 I = 1, N
C        WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
C    $        ( A(I,J), J = 1,N )
C  20 CONTINUE
C
      STOP
      END
