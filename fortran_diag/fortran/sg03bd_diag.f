C
C SPDX-License-Identifier: BSD-3-Clause
C
C SG03BD DIAGNOSTIC PROGRAM - FORTRAN VERSION
C
C Instrumented version of TSG03BD that prints all intermediate
C values for comparison with C implementation.
C
C     .. Parameters ..
      INTEGER           NIN, NOUT
      PARAMETER         ( NIN = 5, NOUT = 6 )
      INTEGER           NMAX
      PARAMETER         ( NMAX = 20 )
      INTEGER           LDA, LDB, LDE, LDQ, LDZ
      PARAMETER         ( LDA = NMAX, LDB = NMAX, LDE = NMAX,
     $                    LDQ = NMAX, LDZ = NMAX )
      INTEGER           LDWORK
      PARAMETER         ( LDWORK = MAX( 1, 4*NMAX, 6*NMAX-6 ) )
C     .. Local Scalars ..
      CHARACTER*1       DICO, FACT, TRANS
      DOUBLE PRECISION  SCALE
      INTEGER           I, INFO, J, N, M
C     .. Local Arrays ..
      DOUBLE PRECISION  A(LDA,NMAX), ALPHAI(NMAX), ALPHAR(NMAX),
     $                  B(LDB,NMAX), BETA(NMAX), DWORK(LDWORK),
     $                  E(LDE,NMAX), Q(LDQ,NMAX), Z(LDZ,NMAX)
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          SG03BD
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      WRITE ( NOUT, FMT = '(A)' )
     $  '=== SG03BD FORTRAN DIAGNOSTIC ==='
      WRITE ( NOUT, FMT = '(A)' ) ''
C     Skip the heading in the data file and read the data.
      READ ( NIN, FMT = '()' )
      READ ( NIN, FMT = * ) N, M, DICO, FACT, TRANS

      WRITE ( NOUT, FMT = '(A,I3)' ) 'N = ', N
      WRITE ( NOUT, FMT = '(A,I3)' ) 'M = ', M
      WRITE ( NOUT, FMT = '(A,A1)' ) 'DICO = ', DICO
      WRITE ( NOUT, FMT = '(A,A1)' ) 'FACT = ', FACT
      WRITE ( NOUT, FMT = '(A,A1)' ) 'TRANS = ', TRANS

      IF ( N.LT.0 .OR. N.GT.NMAX ) THEN
         WRITE ( NOUT, FMT = '(A,I5)' ) 'ERROR: N out of range: ', N
         STOP
      ELSE IF ( M.LT.0 .OR. M.GT.NMAX ) THEN
         WRITE ( NOUT, FMT = '(A,I5)' ) 'ERROR: M out of range: ', M
         STOP
      ELSE
         READ ( NIN, FMT = * ) ( ( A(I,J), J = 1,N ), I = 1,N )
         READ ( NIN, FMT = * ) ( ( E(I,J), J = 1,N ), I = 1,N )
         IF ( LSAME( FACT, 'F' ) ) THEN
            READ ( NIN, FMT = * ) ( ( Q(I,J), J = 1,N ), I = 1,N )
            READ ( NIN, FMT = * ) ( ( Z(I,J), J = 1,N ), I = 1,N )
         END IF
         IF ( LSAME( TRANS, 'T' ) ) THEN
            READ ( NIN, FMT = * ) ( ( B(I,J), J = 1,M ), I = 1,N )
         ELSE
            READ ( NIN, FMT = * ) ( ( B(I,J), J = 1,N ), I = 1,M )
         END IF

C        Print input matrices with high precision
         WRITE ( NOUT, FMT = '(/,A)' ) 'INPUT A matrix:'
         DO 10 I = 1, N
            WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $           ( A(I,J), J = 1,N )
   10    CONTINUE

         WRITE ( NOUT, FMT = '(/,A)' ) 'INPUT E matrix:'
         DO 15 I = 1, N
            WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $           ( E(I,J), J = 1,N )
   15    CONTINUE

         IF ( LSAME( FACT, 'F' ) ) THEN
            WRITE ( NOUT, FMT = '(/,A)' ) 'INPUT Q matrix:'
            DO 17 I = 1, N
               WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $              ( Q(I,J), J = 1,N )
   17       CONTINUE

            WRITE ( NOUT, FMT = '(/,A)' ) 'INPUT Z matrix:'
            DO 18 I = 1, N
               WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $              ( Z(I,J), J = 1,N )
   18       CONTINUE
         END IF

         WRITE ( NOUT, FMT = '(/,A)' ) 'INPUT B matrix:'
         IF ( LSAME( TRANS, 'T' ) ) THEN
            DO 20 I = 1, N
               WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $              ( B(I,J), J = 1,M )
   20       CONTINUE
         ELSE
            DO 25 I = 1, M
               WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $              ( B(I,J), J = 1,N )
   25       CONTINUE
         END IF

C        Call SG03BD
         CALL SG03BD( DICO, FACT, TRANS, N, M, A, LDA, E, LDE, Q, LDQ,
     $                Z, LDZ, B, LDB, SCALE, ALPHAR, ALPHAI, BETA,
     $                DWORK, LDWORK, INFO )

C        Print results
         WRITE ( NOUT, FMT = '(/,A,I3)' ) 'INFO = ', INFO

         IF ( INFO.EQ.0 ) THEN
            WRITE ( NOUT, FMT = '(/,A,E23.16)' ) 'SCALE = ', SCALE

            WRITE ( NOUT, FMT = '(/,A)' )
     $        'OUTPUT U matrix (Cholesky factor):'
            DO 30 I = 1, N
               WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $              ( B(I,J), J = 1,N )
   30       CONTINUE

            WRITE ( NOUT, FMT = '(/,A)' ) 'ALPHAR eigenvalues:'
            WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $           ( ALPHAR(I), I = 1,N )

            WRITE ( NOUT, FMT = '(/,A)' ) 'ALPHAI eigenvalues:'
            WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $           ( ALPHAI(I), I = 1,N )

            WRITE ( NOUT, FMT = '(/,A)' ) 'BETA values:'
            WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $           ( BETA(I), I = 1,N )

            WRITE ( NOUT, FMT = '(/,A)' )
     $        'OUTPUT A matrix (Schur form):'
            DO 40 I = 1, N
               WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $              ( A(I,J), J = 1,N )
   40       CONTINUE

            WRITE ( NOUT, FMT = '(/,A)' )
     $        'OUTPUT E matrix (triangular):'
            DO 50 I = 1, N
               WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $              ( E(I,J), J = 1,N )
   50       CONTINUE

            WRITE ( NOUT, FMT = '(/,A)' )
     $        'OUTPUT Q matrix (orthogonal):'
            DO 60 I = 1, N
               WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $              ( Q(I,J), J = 1,N )
   60       CONTINUE

            WRITE ( NOUT, FMT = '(/,A)' )
     $        'OUTPUT Z matrix (orthogonal):'
            DO 70 I = 1, N
               WRITE ( NOUT, FMT = '(20(1X,E23.16))' )
     $              ( Z(I,J), J = 1,N )
   70       CONTINUE
         ELSE
            WRITE ( NOUT, FMT = '(A)' )
     $        'SG03BD failed - no output matrices'
         END IF
      END IF

      STOP
      END
