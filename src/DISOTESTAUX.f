      SUBROUTINE  GETMOM( IPHAS, GG, NMOM, PMOM )
c        Calculate phase function Legendre expansion coefficients
c        in various special cases
c       INPUT: IPHAS   Phase function options
c                      1 : Isotropic
c                      2 : Rayleigh
c                      3 : Henyey-Greenstein with asymmetry factor GG
c              GG      Asymmetry factor for Henyey-Greenstein case
c              NMOM    Index of highest Legendre coefficient needed
c                        ( = number of streams 'NSTR'  chosen
c                         for the discrete ordinate method)

c      OUTPUT: PMOM(K)  Legendre expansion coefficients (K=0 to NMOM)
c                         (be sure to dimension '0:maxval' in calling
c                          program)
c ------------------------------------------------------------------

c     .. Scalar Arguments ..
      INTEGER   IPHAS, NMOM
      REAL      GG
c     .. Array Arguments ..
      REAL      PMOM( 0:NMOM )
c     .. Local Scalars ..
      INTEGER   K
c     .. External Subroutines ..
      EXTERNAL  ERRMSG

      IF ( IPHAS.LT.1 .OR. IPHAS.GT. 3 )
     &     CALL ERRMSG( 'GETMOM--bad input variable IPHAS',.TRUE.)
      IF ( IPHAS.EQ.3 .AND. (GG.LE.-1.0 .OR. GG.GE.1.0) )
     &     CALL ERRMSG( 'GETMOM--bad input variable GG',.TRUE.)
      IF ( NMOM.LT.2 )
     &     CALL ERRMSG( 'GETMOM--bad input variable NMOM',.TRUE.)


      PMOM(0) = 1.0
      DO  10  K = 1, NMOM
         PMOM(K) = 0.0
   10 CONTINUE

      IF ( IPHAS.EQ.2 )  THEN ! Rayleigh phase function
          PMOM(2) = 0.1
      ELSE IF ( IPHAS.EQ.3 ) THEN ! Henyey-Greenstein phase fcn
         DO  20  K = 1, NMOM
            PMOM(K) = GG**K
   20    CONTINUE
      END IF
      END

