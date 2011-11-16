      SUBROUTINE poth2oh(R,POT,DER)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C     POTENTIAL FOR H ATOM MOVING IN THE FIELD OF A RIGID H2O MOLECULE
!C
!C
!C   1 H -----O 2
!C           / \
!C          H   H
!C         3     4
!C
!C     R(1)=R12; R(2)=R13; R(3)=R14
!C
!C     Input: Distances R(1),R(2),R(3)
!C     Output: Potential POT and potential derivatives DER, i.e., dV/dR
!C
!C     Units: energy in kJ/mol and distance in Angstroms
!C
!C     Stefan Andersson
!C     JCP, 124, 064715 (2006)
!C
!C     NB! Before using the potential: CALL SETUP()
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*8 R(3),POT,DER(3)
      REAL*8 RMHO,RMHH,C6HO,C6HH,AHO,BHO,AHH,BHH,DE,BETA,RHOEQ
      REAL*8 RDAMPHO,RDAMPHH
      REAL*8 C6HO6,C6HH6
      REAL*8 DHO,DHH1,DHH2,DDHO,DDHH1,DDHH2
      REAL*8 POTD1,POTD2,POTD3,DERD1,DERD2,DERD3
      REAL*8 POTR1,POTR2,POTR3,DERR1,DERR2,DERR3
      REAL*8 POTM1,DERM1
      COMMON/POTDAMP/RDAMPHO,RDAMPHH
      COMMON/POT/RMHO,RMHH,C6HO,C6HH,AHO,BHO,AHH,BHH,DE,BETA,RHOEQ
      COMMON/DER/C6HO6,C6HH6

!C     DISPERSION
 
!      PRINT *,R(1)
!      PRINT *,R(2)
!      PRINT *,R(3)

      IF(R(1).GE.RDAMPHO)THEN
        DHO = 1.d0
        DDHO = 0.d0
      ELSE
        DHO = DEXP(-(RDAMPHO/R(1)-1.d0)**2)
        DDHO = 2.d0*(RDAMPHO**2/R(1)**3-RDAMPHO/R(1)**2)*DHO
      ENDIF
      IF(R(2).GE.RDAMPHH)THEN
        DHH1 = 1.d0
        DDH1 = 0.d0
      ELSE
        DHH1 = DEXP(-(RDAMPHH/R(2)-1.d0)**2)
        DDHH1 = 2.d0*(RDAMPHH**2/R(2)**3-RDAMPHH/R(2)**2)*DHH1
      ENDIF
      IF(R(3).GE.RDAMPHH)THEN
        DHH2 = 1.d0
        DDHH2 = 0.d0
      ELSE
        DHH2 = DEXP(-(RDAMPHH/R(3)-1.d0)**2)
        DDHH2 = 2.d0*(RDAMPHH**2/R(3)**3-RDAMPHH/R(3)**2)*DHH2
      ENDIF

         POTD1 = -DHO*C6HO*R(1)**-6
         DERD1 = -DDHO*C6HO*R(1)**-6 + DHO*C6HO6*R(1)**-7

         POTD2 = -DHH1*C6HH*R(2)**-6
         DERD2 = -DDHH1*C6HH*R(2)**-6 + DHH1*C6HH6*R(2)**-7

         POTD3 = -DHH2*C6HH*R(3)**-6
         DERD3 = -DDHH2*C6HH*R(3)**-6 + DHH2*C6HH6*R(3)**-7

!C     REPULSION

         POTR1 = AHO*DEXP(-BHO*R(1))
         DERR1 = -BHO*AHO*DEXP(-BHO*R(1))

         POTR2 = AHH*DEXP(-BHH*R(2)) 
         DERR2 = -BHH*AHH*DEXP(-BHH*R(2))

         POTR3 = AHH*DEXP(-BHH*R(3))
         DERR3 = -BHH*AHH*DEXP(-BHH*R(3))

!C     MORSE
!C     In the supporting info for the paper the R and R_e have been 
!C     interchanged

         POTM1 = DE*(DEXP(-2.D0*BETA*(R(1)-RHOEQ)) -2.d0*DEXP(-BETA*(R(1)-RHOEQ)))
         DERM1 = -2.d0*BETA*DE*(DEXP(-2.D0*BETA*(R(1)-RHOEQ))-DEXP(-BETA*(R(1)-RHOEQ))) 

         POT = POTD1 + POTD2 + POTD3 + POTR1 + POTR2 + POTR3 + POTM1

         DER(1) = DERD1 + DERR1 + DERM1
         DER(2) = DERD2 + DERR2
         DER(3) = DERD3 + DERR3

!      PRINT *,DER(1)
!      PRINT *,DER(2)
!      PRINT *,DER(3)
!      print *,POT


      RETURN
      END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      BLOCK DATA
      REAL*8 RMHO,RMHH,C6HO,C6HH,AHO,BHO,AHH,BHH,DE,BETA,RHOEQ
      COMMON/POT/RMHO,RMHH,C6HO,C6HH,AHO,BHO,AHH,BHH,DE,BETA,RHOEQ
      DATA RMHO,RMHH,C6HO,C6HH/4.04741d0,2.19954d0,2216.76d0,30.2562d0/
      DATA AHO,BHO,AHH,BHH/4158.88d0,2.63805d0,5563.57d0,4.27352d0/
      DATA DE,BETA,RHOEQ/287.598d0,4.73771d0,0.816286d0/
      END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE setup()
!C     It is necessary to call this subroutine first!
      REAL*8 RMHO,RMHH,C6HO,C6HH,AHO,BHO,AHH,BHH,DE,BETA,RHOEQ
      REAL*8 C6HO6,C6HH6
      REAL*8 RDAMPHO,RDAMPHH
      COMMON/POTDAMP/RDAMPHO,RDAMPHH
      COMMON/DER/C6HO6,C6HH6
      COMMON/POT/RMHO,RMHH,C6HO,C6HH,AHO,BHO,AHH,BHH,DE,BETA,RHOEQ
      rdampho=1.28d0*RMHO
      rdamphh=1.28d0*RMHH
      C6HO6=6.d0*C6HO
      C6HH6=6.d0*C6HH
      RETURN
      END
