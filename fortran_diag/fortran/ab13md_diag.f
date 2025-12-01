C
C SPDX-License-Identifier: BSD-3-Clause
C
C AB13MD DIAGNOSTIC PROGRAM - FORTRAN VERSION
C
C Instrumented version of TAB13MD that prints all intermediate
C values for comparison with C implementation.
C
C     .. Parameters ..
      INTEGER           NIN, NOUT
      PARAMETER         ( NIN = 5, NOUT = 6 )
      INTEGER           NMAX, MMAX
      PARAMETER         ( NMAX = 10, MMAX = 10 )
      INTEGER           LDZ
      PARAMETER         ( LDZ = NMAX )
      INTEGER           LIWORK
      PARAMETER         ( LIWORK = MAX( 4*MMAX-2, NMAX ) )
      INTEGER           LDWORK
      PARAMETER         ( LDWORK = 2*NMAX*NMAX*MMAX - NMAX*NMAX +
     $                            9*MMAX*MMAX + NMAX*MMAX + 11*NMAX +
     $                            33*MMAX - 11 )
      INTEGER           LZWORK
      PARAMETER         ( LZWORK = 6*NMAX*NMAX*MMAX + 12*NMAX*NMAX +
     $                            6*MMAX + 6*NMAX - 3 )
C     .. Local Scalars ..
      INTEGER           I, INFO, J, M, N
      DOUBLE PRECISION  BOUND
C     .. Local Arrays ..
      INTEGER           ITYPE(MMAX), IWORK(LIWORK), NBLOCK(MMAX)
      DOUBLE PRECISION  D(NMAX), DWORK(LDWORK), G(NMAX), X(2*MMAX-1)
      COMPLEX*16        Z(LDZ,NMAX), ZWORK(LZWORK)
C     .. External Subroutines ..
      EXTERNAL          AB13MD
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, DBLE, DIMAG
C
C     .. Executable Statements ..
C
      WRITE ( NOUT, FMT = '(A)' )
     $  '=== AB13MD FORTRAN DIAGNOSTIC ==='
      WRITE ( NOUT, FMT = '(A)' ) ''
C
C     Skip the heading in the data file and read the data.
      READ ( NIN, FMT = '()' )
      READ ( NIN, FMT = * ) N, M
C
      WRITE ( NOUT, FMT = '(A,I3)' ) 'N = ', N
      WRITE ( NOUT, FMT = '(A,I3)' ) 'M = ', M
C
      IF ( N.LT.0 .OR. N.GT.NMAX ) THEN
         WRITE ( NOUT, FMT = '(A,I5)' ) 'ERROR: N out of range: ', N
         STOP
      END IF
      IF ( M.LT.0 .OR. M.GT.MMAX ) THEN
         WRITE ( NOUT, FMT = '(A,I5)' ) 'ERROR: M out of range: ', M
         STOP
      END IF
C
      READ ( NIN, FMT = * ) ( NBLOCK(I), I = 1, M )
      READ ( NIN, FMT = * ) ( ITYPE(I), I = 1, M )
      READ ( NIN, FMT = * ) ( ( Z(I,J), J = 1,N ), I = 1,N )
C
      WRITE ( NOUT, FMT = '(/,A)' ) 'NBLOCK:'
      WRITE ( NOUT, FMT = '(10I5)' ) ( NBLOCK(I), I = 1, M )
C
      WRITE ( NOUT, FMT = '(/,A)' ) 'ITYPE:'
      WRITE ( NOUT, FMT = '(10I5)' ) ( ITYPE(I), I = 1, M )
C
      WRITE ( NOUT, FMT = '(/,A)' ) 'INPUT Z matrix (real parts):'
      DO 10 I = 1, N
         WRITE ( NOUT, FMT = '(10(1X,E16.9))' )
     $        ( DBLE(Z(I,J)), J = 1,N )
   10 CONTINUE
C
      WRITE ( NOUT, FMT = '(/,A)' ) 'INPUT Z matrix (imag parts):'
      DO 20 I = 1, N
         WRITE ( NOUT, FMT = '(10(1X,E16.9))' )
     $        ( DIMAG(Z(I,J)), J = 1,N )
   20 CONTINUE
C
C     Call the routine
      CALL AB13MD( 'N', N, Z, LDZ, M, NBLOCK, ITYPE, X, BOUND, D, G,
     $             IWORK, DWORK, LDWORK, ZWORK, LZWORK, INFO )
C
      WRITE ( NOUT, FMT = '(/,A,I5)' ) 'INFO = ', INFO
C
      IF ( INFO.EQ.0 ) THEN
         WRITE ( NOUT, FMT = '(/,A)' ) 'BOUND:'
         WRITE ( NOUT, FMT = '(E23.16)' ) BOUND
C
         WRITE ( NOUT, FMT = '(/,A)' ) 'D vector:'
         WRITE ( NOUT, FMT = '(5(1X,E16.9))' ) ( D(I), I = 1, N )
C
         WRITE ( NOUT, FMT = '(/,A)' ) 'G vector:'
         WRITE ( NOUT, FMT = '(5(1X,E16.9))' ) ( G(I), I = 1, N )
C
         WRITE ( NOUT, FMT = '(/,A)' ) 'X vector (scaling):'
         WRITE ( NOUT, FMT = '(5(1X,E16.9))' ) ( X(I), I = 1, M-1+
     $        (M-1) )
      ELSE
         WRITE ( NOUT, FMT = '(A,I3)' ) 'ERROR: INFO = ', INFO
      END IF
C
      STOP
      END
