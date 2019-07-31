*> \brief \b SDOT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       REAL SX(*),SY(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    SDOT forms the dot product of two vectors.
*>    uses unrolled loops for increments equal to one.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] SX
*> \verbatim
*>          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of SX
*> \endverbatim
*>
*> \param[in] SY
*> \verbatim
*>          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>         storage spacing between elements of SY
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup single_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE RIADDC(C,SX,N)
*     SY=SY-SX
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER C, N
*     ..
*     .. Array Arguments ..
      INTEGER SX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               SX(I)=SX(I)+C
            END DO
            IF (N.LT.4) THEN
*               RSDOT=STEMP
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            SX(I)=SX(I)+C
            SX(I+1)=SX(I+1)+C
            SX(I+2)=SX(I+2)+C
            SX(I+3)=SX(I+3)+C
         END DO
      RETURN
      END
      SUBROUTINE RISUM(SX,N,ANS)
      INTEGER ANS,N
      INTEGER SX(*)
      INTEGER I,M,MP1
      ANS = 0
      IF (N.LE.0 ) RETURN
         M = MOD(N,6)
         IF (M.NE.0) THEN
            DO I = 1,M
               ANS = ANS + SX(I)
            END DO
            IF (N.LT.6) THEN
               RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,6
            ANS = ANS + (SX(I)) + (SX(I+1)) +
     $              (SX(I+2)) + (SX(I+3)) +
     $              (SX(I+4)) + (SX(I+5))
         END DO
      RETURN
      END
      SUBROUTINE RSADD(SX,SY,N)
      INTEGER N
      REAL SX(*),SY(*)
      INTEGER I,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               SY(I)=SY(I)+SX(I)
            END DO
            IF (N.LT.4) THEN
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            SY(I)=SY(I)+SX(I)
            SY(I+1)=SY(I+1)+SX(I+1)
            SY(I+2)=SY(I+2)+SX(I+2)
            SY(I+3)=SY(I+3)+SX(I+3)
         END DO
      RETURN
      END
      SUBROUTINE RSCOSVEC(SX,N)
      INTEGER N
      REAL SX(*)
      INTRINSIC SQRT
      INTEGER I,M,MP1
      INTRINSIC MOD, COS
      IF (N.LE.0) RETURN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               SX(I)=COS(SX(I))
            END DO
            IF (N.LT.4) THEN
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            SX(I)=COS(SX(I))
            SX(I+1)=COS(SX(I+1))
            SX(I+2)=COS(SX(I+2))
            SX(I+3)=COS(SX(I+3))
         END DO
      RETURN
      END
      SUBROUTINE RSCSUB(VAL,SX,N)
      INTEGER N
      REAL VAL, SX(*)
      INTEGER I,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               SX(I)=VAL-SX(I)
            END DO
            IF (N.LT.4) THEN
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            SX(I)=VAL-SX(I)
            SX(I+1)=VAL-SX(I+1)
            SX(I+2)=VAL-SX(I+2)
            SX(I+3)=VAL-SX(I+3)
         END DO
      RETURN
      END
      SUBROUTINE RSDIFF(SX,N)
      INTEGER N
      REAL SX(*)
      INTEGER I,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               SX(I)=SX(I+1)-SX(I)
            END DO
            IF (N.LT.4) THEN
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            SX(I)   = SX(I+1)  -SX(I)
            SX(I+1) = SX(I+2)-SX(I+1)
            SX(I+2) = SX(I+3)-SX(I+2)
            SX(I+3) = SX(I+4)-SX(I+3)
         END DO
      RETURN
      END
      REAL FUNCTION RSDOT(N,SX,INCX,SY,INCY)
      INTEGER INCX,INCY,N
      REAL SX(*),SY(*)
      REAL STEMP
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      STEMP = 0.0e0
      RSDOT = 0.0e0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               STEMP = STEMP + SX(I)*SY(I)
            END DO
            IF (N.LT.5) THEN
               RSDOT=STEMP
            RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
          STEMP = STEMP + SX(I)*SY(I) + SX(I+1)*SY(I+1) +
     $            SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            STEMP = STEMP + SX(IX)*SY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RSDOT = STEMP
      RETURN
      END
      SUBROUTINE RSLOGVEC(SX,N)
      INTEGER N
      REAL SX(*)
      INTRINSIC LOG
      INTEGER I,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               SX(I)=LOG(SX(I))
            END DO
            IF (N.LT.4) THEN
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            SX(I)=LOG(SX(I))
            SX(I+1)=LOG(SX(I+1))
            SX(I+2)=LOG(SX(I+2))
            SX(I+3)=LOG(SX(I+3))
         END DO
      RETURN
      END
      SUBROUTINE  RSMAXIDX(SX,N, MAXVAL, MAXIDX)
      INTEGER N, MAXIDX
      REAL SX(*), MAXVAL
      INTEGER I,IX
      MAXIDX = 0
      IF (N.LT.1 ) RETURN
      MAXVAL=SX(1)
      MAXIDX = 0
      IF (N.EQ.1) RETURN      
         
         DO I = 2,N
            IF ((SX(I)).GT.MAXVAL) THEN
               MAXIDX = I-1
               MAXVAL = (SX(I))
            END IF
         END DO
 
      RETURN
      END
      SUBROUTINE  RSMINIDX(SX,N, MINVAL, MINIDX)
      INTEGER N, MINIDX
      REAL SX(*), MINVAL
      INTEGER I,IX
      IF (N.LT.1 ) RETURN
      MINVAL=SX(1)
      MINIDX = 0
      IF (N.EQ.1) RETURN      
         
         DO I = 2,N
            IF ((SX(I)).LT.MINVAL) THEN
               MINIDX = I-1
               MINVAL = (SX(I))
            END IF
         END DO
 
      RETURN
      END
      SUBROUTINE RSMNSD(SX,N,MEAN, STD)
      REAL MEAN,STD
      INTEGER N
      REAL SX(*)
      INTRINSIC  FLOAT,SQRT
      INTEGER I,M,MP1
      REAL TMP1, TMP2, FN
      TMP1=0.0e0
      TMP2=0.0e0
      IF (N.LE.0 ) RETURN
         M = MOD(N,6)
         IF (M.NE.0) THEN
            DO I = 1,M
               TMP1 = TMP1 + SX(I)
               TMP2 = TMP2 + SX(I)*SX(I)
            END DO
            IF (N.LT.6) THEN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,6
            TMP1 = TMP1 + (SX(I)) + (SX(I+1)) +
     $              (SX(I+2)) + (SX(I+3)) +
     $              (SX(I+4)) + (SX(I+5))
            TMP2 = TMP2 + (SX(I)*SX(I)) + (SX(I+1)*SX(I+1)) +
     $              (SX(I+2)*SX(I+2)) + (SX(I+3)*SX(I+3)) +
     $              (SX(I+4)*SX(I+4)) + (SX(I+5)*SX(I+5))
         END DO
         
         FN=FLOAT(N)
         MEAN=TMP1/FN
         STD=SQRT((TMP2-FN*MEAN*MEAN)/(FN-1))
      RETURN
      END
      SUBROUTINE RSPOWVEC(SX,POW,N)
      INTEGER N
      REAL POW
      REAL SX(*)
      INTRINSIC SQRT
      INTEGER I,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               SX(I)= SX(I)**POW
            END DO
            IF (N.LT.4) THEN
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            SX(I)=(SX(I))**POW
            SX(I+1)=(SX(I+1))**POW
            SX(I+2)=(SX(I+2))**POW
            SX(I+3)=(SX(I+3))**POW
         END DO
      RETURN
      END
      SUBROUTINE RSSCAL(N,SA,SX,INCX)
      REAL SA
      INTEGER INCX,N
      REAL SX(*)
      INTEGER I,M,MP1,NINCX
      INTRINSIC MOD
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               SX(I) = SA*SX(I)
            END DO
            IF (N.LT.5) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            SX(I) = SA*SX(I)
            SX(I+1) = SA*SX(I+1)
            SX(I+2) = SA*SX(I+2)
            SX(I+3) = SA*SX(I+3)
            SX(I+4) = SA*SX(I+4)
         END DO
      ELSE
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            SX(I) = SA*SX(I)
         END DO
      END IF
      RETURN
      END
      SUBROUTINE RSSINVEC(SX,N)
      INTEGER N
      REAL SX(*)
      INTRINSIC SIN
      INTEGER I,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               SX(I)=SIN(SX(I))
            END DO
            IF (N.LT.4) THEN
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            SX(I)=SIN(SX(I))
            SX(I+1)=SIN(SX(I+1))
            SX(I+2)=SIN(SX(I+2))
            SX(I+3)=SIN(SX(I+3))
         END DO
      RETURN
      END
      SUBROUTINE RSSQRVEC(SX,SY,N)
      INTEGER N
      REAL SX(*),SY(*)
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               SY(I)=SX(I)*SX(I)
            END DO
            IF (N.LT.4) THEN
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            SY(I)=SX(I)*SX(I)
            SY(I+1)=SX(I+2)*SX(I+1)
            SY(I+2)=SX(I+2)*SX(I+2)
            SY(I+3)=SX(I+3)*SX(I+3)     
         END DO
      RETURN
      END
      SUBROUTINE RSSQRTVEC(SX,N)
      INTEGER N
      REAL SX(*)
      INTRINSIC SQRT
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               SX(I)=SQRT(SX(I))
            END DO
            IF (N.LT.4) THEN
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            SX(I)=SQRT(SX(I))
            SX(I+1)=SQRT(SX(I+1))
            SX(I+2)=SQRT(SX(I+2))
            SX(I+3)=SQRT(SX(I+3))
         END DO
      RETURN
      END
      SUBROUTINE RSSUB(SX,SY, SZ, N)
      INTEGER N
      REAL SX(*),SY(*),SZ(*)
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               SZ(I)=SY(I)-SX(I)
            END DO
            IF (N.LT.4) THEN
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            SZ(I)=SY(I)-SX(I)
            SZ(I+1)=SY(I+1)-SX(I+1)
            SZ(I+2)=SY(I+2)-SX(I+2)
            SZ(I+3)=SY(I+3)-SX(I+3)
         END DO
      RETURN
      END
      SUBROUTINE RSSUB_I(SX,SY,N)
      INTEGER N
      REAL SX(*),SY(*)
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               SY(I)=SY(I)-SX(I)
            END DO
            IF (N.LT.4) THEN
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            SY(I)=SY(I)-SX(I)
            SY(I+1)=SY(I+1)-SX(I+1)
            SY(I+2)=SY(I+2)-SX(I+2)
            SY(I+3)=SY(I+3)-SX(I+3)
         END DO
      RETURN
      END
      SUBROUTINE RSSUBC(VAL,SX,N)
      INTEGER N
      REAL VAL, SX(*)
      INTEGER I,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               SX(I)=SX(I)-VAL
            END DO
            IF (N.LT.4) THEN
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            SX(I)=SX(I)-VAL
            SX(I+1)=SX(I+1)-VAL
            SX(I+2)=SX(I+2)-VAL
            SX(I+3)=SX(I+3)-VAL    
         END DO
      RETURN
      END
      SUBROUTINE RSSUM(SX,N,ANS)
      REAL ANS
      INTEGER N
      REAL SX(*)
      INTEGER I,M,MP1
      ANS = 0.0e0
      IF (N.LE.0 ) RETURN
         M = MOD(N,6)
         IF (M.NE.0) THEN
            DO I = 1,M
               ANS = ANS + SX(I)
            END DO
            IF (N.LT.6) THEN
               RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,6
            ANS = ANS + (SX(I)) + (SX(I+1)) +
     $              (SX(I+2)) + (SX(I+3)) +
     $              (SX(I+4)) + (SX(I+5))
         END DO
      RETURN
      END
      SUBROUTINE RSSUMLOG(SX,N, ANS)
      INTEGER N
      REAL SX(*), ANS
      INTEGER I,M,MP1
      INTRINSIC MOD, LOG      
      ANS=0.0e0
      IF (N.LE.0) RETURN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
              ANS=ANS+LOG(SX(I))
            END DO
            IF (N.LT.4) THEN
                RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
            ANS=ANS+LOG(SX(I))
            ANS=ANS+LOG(SX(I+1))
            ANS=ANS+LOG(SX(I+2))
            ANS=ANS+LOG(SX(I+3))
         END DO
      RETURN
      END
      LOGICAL FUNCTION RLSAME(CA,CB)
      CHARACTER CA,CB
      INTRINSIC ICHAR
      INTEGER INTA,INTB,ZCODE
      RLSAME = CA .EQ. CB
      IF (RLSAME) RETURN
      ZCODE = ICHAR('Z')
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
          IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
          IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
          IF (INTA.GE.129 .AND. INTA.LE.137 .OR.
     +        INTA.GE.145 .AND. INTA.LE.153 .OR.
     +        INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
          IF (INTB.GE.129 .AND. INTB.LE.137 .OR.
     +        INTB.GE.145 .AND. INTB.LE.153 .OR.
     +        INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
          IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
          IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
      RLSAME = INTA .EQ. INTB
      END
      SUBROUTINE RSGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,
     + B,LDB,BETA,C,LDC)
      REAL ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
      REAL A(LDA,*),B(LDB,*),C(LDC,*)
      LOGICAL RLSAME
      EXTERNAL RLSAME
      EXTERNAL RXERBLA
      INTRINSIC MAX
      REAL TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL NOTA,NOTB
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      NOTA = RLSAME(TRANSA,'N')
      NOTB = RLSAME(TRANSB,'N')
      IF (NOTA) THEN
          NROWA = M
          NCOLA = K
      ELSE
          NROWA = K
          NCOLA = M
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.RLSAME(TRANSA,'C')) .AND.
     +    (.NOT.RLSAME(TRANSA,'T'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.RLSAME(TRANSB,'C')) .AND.
     +         (.NOT.RLSAME(TRANSB,'T'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
          CALL RXERBLA('SGEMM ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
      IF (NOTB) THEN
          IF (NOTA) THEN
              DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                      TEMP = ALPHA*B(L,J)
                      DO 70 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (NOTA) THEN
              DO 170 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 130 I = 1,M
                          C(I,J) = ZERO
  130                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 140 I = 1,M
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  END IF
                  DO 160 L = 1,K
                      TEMP = ALPHA*B(J,L)
                      DO 150 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  150                 CONTINUE
  160             CONTINUE
  170         CONTINUE
          ELSE
              DO 200 J = 1,N
                  DO 190 I = 1,M
                      TEMP = ZERO
                      DO 180 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
      RETURN
      END
      SUBROUTINE RSGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      REAL ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
      REAL A(LDA,*),X(*),Y(*)
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
      LOGICAL RLSAME
      EXTERNAL RLSAME
      EXTERNAL RXERBLA
      INTRINSIC MAX
      INFO = 0
      IF (.NOT.RLSAME(TRANS,'N') .AND. .NOT.RLSAME(TRANS,'T') .AND.
     +    .NOT.RLSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL RXERBLA('SGEMV ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
      IF (RLSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (RLSAME(TRANS,'N')) THEN
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  TEMP = ALPHA*X(JX)
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 100 J = 1,N
                  TEMP = ZERO
                  DO 90 I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  100         CONTINUE
          ELSE
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO 110 I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
      RETURN
      END
      SUBROUTINE RSSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      REAL ALPHA,BETA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
      REAL A(LDA,*),X(*),Y(*)
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
      LOGICAL RLSAME
      EXTERNAL RLSAME
      EXTERNAL RXERBLA
      INTRINSIC MAX
      INFO = 0
      IF (.NOT.RLSAME(UPLO,'U') .AND. .NOT.RLSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 5
      ELSE IF (INCX.EQ.0) THEN
          INFO = 7
      ELSE IF (INCY.EQ.0) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL RXERBLA('SSYMV ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (N-1)*INCY
      END IF
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,N
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,N
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,N
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,N
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (RLSAME(UPLO,'U')) THEN
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  DO 50 I = 1,J - 1
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   50             CONTINUE
                  Y(J) = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2
   60         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 80 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  IX = KX
                  IY = KY
                  DO 70 I = 1,J - 1
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
                      IX = IX + INCX
                      IY = IY + INCY
   70             CONTINUE
                  Y(JY) = Y(JY) + TEMP1*A(J,J) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
   80         CONTINUE
          END IF
      ELSE
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 100 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  Y(J) = Y(J) + TEMP1*A(J,J)
                  DO 90 I = J + 1,N
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   90             CONTINUE
                  Y(J) = Y(J) + ALPHA*TEMP2
  100         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 120 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  Y(JY) = Y(JY) + TEMP1*A(J,J)
                  IX = JX
                  IY = JY
                  DO 110 I = J + 1,N
                      IX = IX + INCX
                      IY = IY + INCY
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
      RETURN
      END
      SUBROUTINE RSSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
      REAL ALPHA,BETA
      INTEGER K,LDA,LDC,N
      CHARACTER TRANS,UPLO
      REAL A(LDA,*),C(LDC,*)
      LOGICAL RLSAME
      EXTERNAL RLSAME
      EXTERNAL RXERBLA
      INTRINSIC MAX
      REAL TEMP
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      IF (RLSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = RLSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.RLSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.RLSAME(TRANS,'N')) .AND.
     +         (.NOT.RLSAME(TRANS,'T')) .AND.
     +         (.NOT.RLSAME(TRANS,'C'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL RXERBLA('SSYRK ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
      IF (RLSAME(TRANS,'N')) THEN
          IF (UPPER) THEN
              DO 130 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 100 I = 1,J
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                  END IF
                  DO 120 L = 1,K
                      IF (A(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*A(J,L)
                          DO 110 I = 1,J
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 150 I = J,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  END IF
                  DO 170 L = 1,K
                      IF (A(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*A(J,L)
                          DO 160 I = J,N
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
          IF (UPPER) THEN
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP = ZERO
                      DO 190 L = 1,K
                          TEMP = TEMP + A(L,I)*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP = ZERO
                      DO 220 L = 1,K
                          TEMP = TEMP + A(L,I)*A(L,J)
  220                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
      RETURN
      END
      SUBROUTINE RSTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
      REAL A(LDA,*),X(*)
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
      REAL TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOUNIT
      LOGICAL RLSAME
      EXTERNAL RLSAME
      EXTERNAL RXERBLA
      INTRINSIC MAX
      INFO = 0
      IF (.NOT.RLSAME(UPLO,'U') .AND. .NOT.RLSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.RLSAME(TRANS,'N') .AND. .NOT.RLSAME(TRANS,'T') .AND.
     +         .NOT.RLSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.RLSAME(DIAG,'U') .AND. .NOT.RLSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL RXERBLA('STRMV ',INFO)
          RETURN
      END IF
      IF (N.EQ.0) RETURN
      NOUNIT = RLSAME(DIAG,'N')
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (RLSAME(TRANS,'N')) THEN
          IF (RLSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 10 I = 1,J - 1
                              X(I) = X(I) + TEMP*A(I,J)
   10                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX
                  DO 40 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 30 I = 1,J - 1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX + INCX
   30                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX + INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 50 I = N,J + 1,-1
                              X(I) = X(I) + TEMP*A(I,J)
   50                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   60             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 70 I = N,J + 1,-1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX - INCX
   70                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX - INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
          IF (RLSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 100 J = N,1,-1
                      TEMP = X(J)
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 90 I = J - 1,1,-1
                          TEMP = TEMP + A(I,J)*X(I)
   90                 CONTINUE
                      X(J) = TEMP
  100             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 120 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 110 I = J - 1,1,-1
                          IX = IX - INCX
                          TEMP = TEMP + A(I,J)*X(IX)
  110                 CONTINUE
                      X(JX) = TEMP
                      JX = JX - INCX
  120             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140 J = 1,N
                      TEMP = X(J)
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 130 I = J + 1,N
                          TEMP = TEMP + A(I,J)*X(I)
  130                 CONTINUE
                      X(J) = TEMP
  140             CONTINUE
              ELSE
                  JX = KX
                  DO 160 J = 1,N
                      TEMP = X(JX)
                      IX = JX
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 150 I = J + 1,N
                          IX = IX + INCX
                          TEMP = TEMP + A(I,J)*X(IX)
  150                 CONTINUE
                      X(JX) = TEMP
                      JX = JX + INCX
  160             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
      SUBROUTINE RSTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      REAL ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
      REAL A(LDA,*),B(LDB,*)
      LOGICAL RLSAME
      EXTERNAL RLSAME
      EXTERNAL RXERBLA
      INTRINSIC MAX
      REAL TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      LSIDE = RLSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOUNIT = RLSAME(DIAG,'N')
      UPPER = RLSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.RLSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.RLSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.RLSAME(TRANSA,'N')) .AND.
     +         (.NOT.RLSAME(TRANSA,'T')) .AND.
     +         (.NOT.RLSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.RLSAME(DIAG,'U')).AND.(.NOT.RLSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL RXERBLA('STRSM ',INFO)
          RETURN
      END IF
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
      IF (LSIDE) THEN
          IF (RLSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 60 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      END IF
                      DO 50 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      END IF
                      DO 90 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 130 J = 1,N
                      DO 120 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          DO 110 K = 1,I - 1
                              TEMP = TEMP - A(K,I)*B(K,J)
  110                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  120                 CONTINUE
  130             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          DO 140 K = I + 1,M
                              TEMP = TEMP - A(K,I)*B(K,J)
  140                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (RLSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 210 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 170 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  170                     CONTINUE
                      END IF
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 200 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  200                     CONTINUE
                      END IF
  210             CONTINUE
              ELSE
                  DO 260 J = N,1,-1
                      IF (ALPHA.NE.ONE) THEN
                          DO 220 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  220                     CONTINUE
                      END IF
                      DO 240 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              DO 230 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 250 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 310 K = N,1,-1
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
                      DO 290 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 280 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  280                         CONTINUE
                          END IF
  290                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 300 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  300                     CONTINUE
                      END IF
  310             CONTINUE
              ELSE
                  DO 360 K = 1,N
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 320 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  320                     CONTINUE
                      END IF
                      DO 340 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 330 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  330                         CONTINUE
                          END IF
  340                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 350 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  350                     CONTINUE
                      END IF
  360             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
      SUBROUTINE RXERBLA( SRNAME, INFO )
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
      INTRINSIC          LEN_TRIM
      RETURN 
 9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
      END
      INTEGER FUNCTION RILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME, TWOSTAGE
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*16
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
      INTEGER            IEEECK, IPARMQ, IPARAM2STAGE
      EXTERNAL           IEEECK, IPARMQ, IPARAM2STAGE
          
	 IF( ISPEC .LE. 0) THEN
		  RILAENV = -1
		  RETURN
	  END IF
      SELECT CASE (ISPEC)
		CASE (1:3)	
   10 CONTINUE
      RILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $             I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
      TWOSTAGE = LEN( SUBNAM ).GE.11
     $           .AND. SUBNAM( 11: 11 ).EQ.'2'
      GO TO ( 50, 60, 70 )ISPEC
   50 CONTINUE
      NB = 1
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'QR ') THEN
            IF( N3 .EQ. 1) THEN
               IF( SNAME ) THEN
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               ELSE
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               END IF
            ELSE
               IF( SNAME ) THEN
                  NB = 1
               ELSE
                  NB = 1
               END IF
            END IF
         ELSE IF( C3.EQ.'LQ ') THEN
            IF( N3 .EQ. 2) THEN
               IF( SNAME ) THEN
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               ELSE
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               END IF
            ELSE
               IF( SNAME ) THEN
                  NB = 1
               ELSE
                  NB = 1
               END IF
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( TWOSTAGE ) THEN
                  NB = 192
               ELSE
                  NB = 64
               END IF
            ELSE
               IF( TWOSTAGE ) THEN
                  NB = 192
               ELSE
                  NB = 64
               END IF
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( TWOSTAGE ) THEN
               NB = 192
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF ( C3.EQ.'EVC' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      ELSE IF( C2.EQ.'GG' ) THEN
         NB = 32
         IF( C3.EQ.'HD3' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         END IF
      END IF
      RILAENV = NB
      RETURN
   60 CONTINUE
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'GG' ) THEN
         NBMIN = 2
         IF( C3.EQ.'HD3' ) THEN
            NBMIN = 2
         END IF
      END IF
      RILAENV = NBMIN
      RETURN
   70 CONTINUE
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'GG' ) THEN
         NX = 128
         IF( C3.EQ.'HD3' ) THEN
            NX = 128
         END IF
      END IF
      RILAENV = NX
      RETURN
	  CASE (4)
   80 CONTINUE
      RILAENV = 6
      RETURN
	  CASE (5)
   90 CONTINUE
      RILAENV = 2
      RETURN
	  CASE (6)
  100 CONTINUE
      RILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
	  CASE (7)
  110 CONTINUE
      RILAENV = 1
      RETURN
	  CASE (8)
  120 CONTINUE
      RILAENV = 50
      RETURN
	  CASE (9)
  130 CONTINUE
      RILAENV = 25
      RETURN
	  CASE (10)
  140 CONTINUE
      RILAENV = 1
      IF( RILAENV.EQ.1 ) THEN
         RILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
	  CASE (11)
  150 CONTINUE
      RILAENV = 1
      IF( RILAENV.EQ.1 ) THEN
         RILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
	  CASE (12:16)
  160 CONTINUE
      RILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
	  END SELECT
      END
      LOGICAL FUNCTION RSISNAN( SIN )
      REAL, INTENT(IN) :: SIN
      LOGICAL SLAISNAN
      EXTERNAL SLAISNAN
      RSISNAN = SLAISNAN(SIN,SIN)
      RETURN
      END
      SUBROUTINE RSPOTRF( UPLO, N, A, LDA, INFO )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
      REAL               A( LDA, * )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      LOGICAL            UPPER
      INTEGER            J, JB, NB
      LOGICAL            RLSAME
      INTEGER            RILAENV
      EXTERNAL           RLSAME, RILAENV
      EXTERNAL           RSGEMM, RSPOTRF2, RSSYRK, RSTRSM, RXERBLA
      INTRINSIC          MAX, MIN
      INFO = 0
      UPPER = RLSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.RLSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL RXERBLA( 'SPOTRF', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 )
     $   RETURN
      NB = RILAENV( 1, 'SPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
         CALL RSPOTRF2( UPLO, N, A, LDA, INFO )
      ELSE
         IF( UPPER ) THEN
            DO 10 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
               CALL RSSYRK( 'Upper', 'Transpose', JB, J-1, -ONE,
     $                     A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL RSPOTRF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
                  CALL RSGEMM('Transpose', 'No transpose', JB, N-J-JB+1,
     $                        J-1, -ONE, A( 1, J ), LDA, A( 1, J+JB ),
     $                        LDA, ONE, A( J, J+JB ), LDA )
                  CALL RSTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit',
     $                        JB, N-J-JB+1, ONE, A( J, J ), LDA,
     $                        A( J, J+JB ), LDA )
               END IF
   10       CONTINUE
         ELSE
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
               CALL RSSYRK( 'Lower', 'No transpose', JB, J-1, -ONE,
     $                     A( J, 1 ), LDA, ONE, A( J, J ), LDA )
               CALL RSPOTRF2( 'Lower', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
                  CALL RSGEMM('No transpose', 'Transpose', N-J-JB+1, JB,
     $                        J-1, -ONE, A( J+JB, 1 ), LDA, A( J, 1 ),
     $                        LDA, ONE, A( J+JB, J ), LDA )
                  CALL RSTRSM('Right', 'Lower', 'Transpose', 'Non-unit',
     $                        N-J-JB+1, JB, ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
   30 CONTINUE
      INFO = INFO + J - 1
   40 CONTINUE
      RETURN
      END
      RECURSIVE SUBROUTINE RSPOTRF2( UPLO, N, A, LDA, INFO )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
      REAL               A( LDA, * )
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO=0.0E+0 )
      LOGICAL            UPPER
      INTEGER            N1, N2, IINFO
      LOGICAL            RLSAME, RSISNAN
      EXTERNAL           RLSAME, RSISNAN
      EXTERNAL           RSSYRK, RSTRSM, RXERBLA
      INTRINSIC          MAX, SQRT
      INFO = 0
      UPPER = RLSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.RLSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL RXERBLA( 'SPOTRF2', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 )
     $   RETURN
      IF( N.EQ.1 ) THEN
         IF( A( 1, 1 ).LE.ZERO.OR.RSISNAN( A( 1, 1 ) ) ) THEN
            INFO = 1
            RETURN
         END IF
         A( 1, 1 ) = SQRT( A( 1, 1 ) )
      ELSE
         N1 = N/2
         N2 = N-N1
         CALL RSPOTRF2( UPLO, N1, A( 1, 1 ), LDA, IINFO )
         IF ( IINFO.NE.0 ) THEN
            INFO = IINFO
            RETURN
         END IF
         IF( UPPER ) THEN
            CALL RSTRSM( 'L', 'U', 'T', 'N', N1, N2, ONE,
     $                  A( 1, 1 ), LDA, A( 1, N1+1 ), LDA )
            CALL RSSYRK( UPLO, 'T', N2, N1, -ONE, A( 1, N1+1 ), LDA,
     $                  ONE, A( N1+1, N1+1 ), LDA )
            CALL RSPOTRF2( UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO )
            IF ( IINFO.NE.0 ) THEN
               INFO = IINFO + N1
               RETURN
            END IF
         ELSE
            CALL RSTRSM( 'R', 'L', 'T', 'N', N2, N1, ONE,
     $                  A( 1, 1 ), LDA, A( N1+1, 1 ), LDA )
            CALL RSSYRK( UPLO, 'N', N2, N1, -ONE, A( N1+1, 1 ), LDA,
     $                  ONE, A( N1+1, N1+1 ), LDA )
            CALL RSPOTRF2( UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO )
            IF ( IINFO.NE.0 ) THEN
               INFO = IINFO + N1
               RETURN
            END IF
         END IF
      END IF
      RETURN
      END
      SUBROUTINE RSPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
      REAL               A( LDA, * ), B( LDB, * )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      LOGICAL            UPPER
      LOGICAL            RLSAME
      EXTERNAL           RLSAME
      EXTERNAL           RSTRSM, RXERBLA
      INTRINSIC          MAX
      INFO = 0
      UPPER = RLSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.RLSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL RXERBLA( 'SPOTRS', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
      IF( UPPER ) THEN
         CALL RSTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
         CALL RSTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
         CALL RSTRSM( 'Left', 'Lower', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
         CALL RSTRSM( 'Left', 'Lower', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
      END IF
      RETURN
      END
      SUBROUTINE RSTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,
     $                   INFO )
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
      REAL               A( LDA, * ), B( LDB, * )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      LOGICAL            NOUNIT
      LOGICAL            RLSAME
      EXTERNAL           RLSAME
      EXTERNAL           RSTRSM, RXERBLA
      INTRINSIC          MAX
      INFO = 0
      NOUNIT = RLSAME( DIAG, 'N' )
      IF( .NOT.RLSAME( UPLO, 'U' ) .AND. .NOT.RLSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.RLSAME( TRANS, 'N' ) .AND. .NOT.
     $       RLSAME( TRANS, 'T' ) .AND. .NOT.RLSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.RLSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL RXERBLA( 'STRTRS', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 )
     $   RETURN
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
      END IF
      INFO = 0
      CALL RSTRSM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B,
     $            LDB )
      RETURN
      END
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
      INTEGER            ISPEC
      REAL               ONE, ZERO
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
      IEEECK = 1
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
      IF( ISPEC.EQ.0 )
     $   RETURN
      NAN1 = POSINF + NEGINF
      NAN2 = POSINF / NEGINF
      NAN3 = POSINF / POSINF
      NAN4 = POSINF*ZERO
      NAN5 = NEGINF*NEGZRO
      NAN6 = NAN5*ZERO
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
      RETURN
      END
      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,
     $                   ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14,
     $                   NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
      INTEGER            NH, NS
      INTEGER            I, IC, IZ
      CHARACTER          SUBNAM*6
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR.
     $    ( ISPEC.EQ.IACC22 ) ) THEN
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 )
     $      NS = 4
         IF( NH.GE.60 )
     $      NS = 10
         IF( NH.GE.150 )
     $      NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
         IF( NH.GE.590 )
     $      NS = 64
         IF( NH.GE.3000 )
     $      NS = 128
         IF( NH.GE.6000 )
     $      NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
      IF( ISPEC.EQ.INMIN ) THEN
         IPARMQ = NMIN
      ELSE IF( ISPEC.EQ.INIBL ) THEN
         IPARMQ = NIBBLE
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
         IPARMQ = NS
      ELSE IF( ISPEC.EQ.INWIN ) THEN
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
         IPARMQ = 0
         SUBNAM = NAME
         IC = ICHAR( SUBNAM( 1: 1 ) )
         IZ = ICHAR( 'Z' )
         IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
            IF( IC.GE.97 .AND. IC.LE.122 ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( IC.GE.97 .AND. IC.LE.122 )
     $               SUBNAM( I: I ) = CHAR( IC-32 )
               END DO
            END IF
         ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
            IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $          ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $          ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC+64 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $                ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $                ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $                I ) = CHAR( IC+64 )
               END DO
            END IF
         ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
            IF( IC.GE.225 .AND. IC.LE.250 ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( IC.GE.225 .AND. IC.LE.250 )
     $               SUBNAM( I: I ) = CHAR( IC-32 )
               END DO
            END IF
         END IF
         IF( SUBNAM( 2:6 ).EQ.'GGHRD' .OR.
     $       SUBNAM( 2:6 ).EQ.'GGHD3' ) THEN
            IPARMQ = 1
            IF( NH.GE.K22MIN )
     $         IPARMQ = 2
         ELSE IF ( SUBNAM( 4:6 ).EQ.'EXC' ) THEN
            IF( NH.GE.KACMIN )
     $         IPARMQ = 1
            IF( NH.GE.K22MIN )
     $         IPARMQ = 2
         ELSE IF ( SUBNAM( 2:6 ).EQ.'HSEQR' .OR.
     $             SUBNAM( 2:5 ).EQ.'LAQR' ) THEN
            IF( NS.GE.KACMIN )
     $         IPARMQ = 1
            IF( NS.GE.K22MIN )
     $         IPARMQ = 2
         END IF
      ELSE
         IPARMQ = -1
      END IF
      END
      LOGICAL FUNCTION SLAISNAN( SIN1, SIN2 )
      REAL, INTENT(IN) :: SIN1, SIN2
      SLAISNAN = (SIN1.NE.SIN2)
      RETURN
      END
