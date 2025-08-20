      PROGRAM TOVSOLVER                   
C     *********************************************************
C     
C     RMF BSP parameter set
c      for p-anisotropic slow rotating isotropic NS, 24 Nov 2015
C     R. L. Bowers and E. P. T. Liang, Astrophys. J 88, 657 (1974)
C     LBL=-2 and 2             
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
C      INTEGER I, J, IM, NS, LI, N, IL, IN    
C      DIMENSION YA(10), EK(5,10), Y(10)
       INTEGER IL
C--------------------------------------------------------
C result save in file X.dat 
C--------------------
      
      OPEN (unit=8,STATUS='unknown',FILE='data1_LBL_K15.6_L.dat')
      OPEN (unit=9,STATUS='unknown',FILE='data2_LBL_K15.6_L.dat')

c Set constants: 
c HC: hbar c, 
c PI, 
c GS: Newton constant in m/(m^3 MeV fm^-3)
c MSS: solar mass in m^3 MeV fm^-3
      HC  = 197.327D0
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
  
C---------------------------------------------------------------------
C  PCC is pressure in center
C  TUCT is trial initial condition for metric component nu
c Start looping for different central pressure
c       DO 10 IL=600,2,-1
        DO 10 IL=2,600,1
C        DO 10 IL=400,400,1   ! For one profile
        
       PCC=1.0D0*IL
       TUCT=0.1D-8

c       PCC=330.D0
C------------------------------------------------------------------------
c Initial steps to calculate the correct initial condition for 
C metric component nu
C   ROS is radius in km
C   GMOS is mass in solar mass at r=R=ROS
C   MNURT is metric component nu at r=R=ROS
C   Rot is the correct initial condition for metric component nu
C
      CALL TOVI(PCC,TUCT,ROS,GMOS,MNURT)


C----------------------------------------------------------------------
C In subroutine TOVI the calculation for pressure and mass had been
C done, but initial condition for nu is NOT determined yet 
C Since Exp(nu) -> Exp(nu+constant) does not change the EoM, then
C we can obtain nu(0) -> nu(0)=nu(0)-(nu(R)-Log(1-2GM/R))

      RNS=ROS*1.D3
      DBS=1.D0-2.D0*GS*GMOS*MSS/RNS

      KC=DLOG(DBS)-MNURT
      ROT= TUCT + KC
      
      
C----------------------------------------
C  Check metric nu
C----------------------------------------- 
C   MAMA is metric component nu in R calculate using initial correct nu
c      CALL TOVI(PCC,ROT,ROS,GMOS,MAMA)
c      CMNUR=DLOG(1.D0-2.D0*GS*GMOS*MSS/RNS)
c      WRITE(*,*)MAMA,MNURT,CMNUR,PCC,ROT,ROS,GMOS,MAMA
C--------------------------------------------------
c Moment of Inertia related properties
c----------------------------------------------------
C   ROS2 is radius in km calculate using initial correct nu
C   GMOS2 is mass in solar mass calculate using initial correct nu
C   OMEGA is rotation frequency and KAPPA is nedded to calculate moment
C   of inertia
C------------------------------------------------------
C In subroutine TOVMI, boundary conditions for both 
C omega and kappa are already satisfied using 
C multiplication by a constant zeta to both of them,
C which is similar to what we have done to nu
C------------------------------------------------------
C   MOMIN is I/MR^2 dimensionless of NS
C   MI is moment of inertia in 10^45 g cm^2 of NS
       PMIN=1.0D-9
       CALL TOVMI(PCC,ROT,PMIN,ROS2,GMOS2,MNURT,OMEGA,KAPPA)
       RNS2=ROS2*1.D3
       MOMIN=KAPPA/(GS*GMOS*MSS*RNS2*RNS2)
       MI=MOMIN*1.98892D33*1.D10*GMOS*ROS2*ROS2/1.0D45
      
c-----------------------------------------------------------
c     CRUST Properties
c----------------------------------------------------------
C     PT pressure at core-crust transition
C     RNST is core radius
c     MOMINC I/MR^2 dimesionless of the core of NS
c     MIC is moment of Inertia    in 10^45 g cm^2 of NS core
c     MICR,MGCR,RPMIC,RPMGC are crust related properties
c--------------------------------------------------------------       
       PT=2.863D-1
       CALL TOVMI(PCC,ROT,PT,RT,GMC,MNUC,OMEGAC,KAPPAC)
       RNST=RT*1.D3
       MOMINC=KAPPAC/(GS*GMC*MSS*RNST*RNST)
       MIC=MOMINC*1.98892D33*1.D10*GMC*RT*RT/1.0D45
C------------------------------------------------------------
       MICR=MI-MIC
       MGCR=GMOS2-GMC
       RPMIC=MICR/MI*1D2
       RPMGC= MGCR/GMOS2*1D2
    

       !WRITE(*,*)IL,ROS,GMOS2,MOMIN
       WRITE(8,*)IL,PCC,ROS2,RT,GMOS2,MOMIN,MI
       WRITE(9,*)IL,GMOS2,MGCR,MICR,RPMIC,RPMGC

 10   CONTINUE

       
      STOP
      END

C----------------------------------------------------------
C TOV Inertia Moment 
C--------------------------------------------
           
      SUBROUTINE TOVMI(PCC,TUCT,PMIN,ROS,GMOS,MNURT,OMEGA,KAPPA)
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IN, IL
      DIMENSION YA(10), EK(5,10), Y(10)
      HC  = 197.327D0
      PI  = 3.14159265358979D0  
      GS=1.325D-12
      MSS=1.1155D15


C     IM = NUMBER OF EQUATIONS (6)
C     IN =  NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS (5)


      IM=6
      IN=IM-1

C---------------------------------------------------
C     Y(1)=Pressure, 
C     Y(2)=Star Mass,
C     Y(3)=Metric Nu
C     Y(4)=Omega, 
C     Y(5)=Kappa,
C----------------------------------------------------
C     Y(6)=Energy Density, 

         
      Y(1)=PCC
      Y(3)=TUCT

      Y(2)=0.1D-8    
      Y(4)=0.1D-8
      Y(5)=0.1D-8

     
c--------------------------------------------------------------------------
c Compute initial Energy density
c---------------------------------------------------------------------      

      
      Y(6)=FED(PCC)


     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D-1
      NS=8   ! recommended 32
      XL=30.0D3

      H=PU/NS
      XP=1.0D-3
      HH=H/(2.0D0)

C     LINE NUMBER INITIALIZATION

      LI=0
 

 28   LI=LI+1

C     XB=OLD radius, XP=NEW radius, XM=MIDPOINT radius

      DO N=1,NS
         XB=XP
         XP=XP+H
         XM=XB+HH


C    COMPUTE K1, L1, M1, N1, O1

         J=1
         DO I=1,IM
            YA(I)=Y(I)
         END DO
         XA=XB

   
         CALL FUNCX(EK,J,YA,XA,H)

C    COMPUTE K2, L2, M2, N2, O2

         J=2
         
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
            P0=YA(1)
c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(6)=ED
c---------------------------------------------------------------------
                     
         XA=XM

         CALL FUNCX(EK,J,YA,XA,H)

C    COMPUTE K3, L3, M3, N3, O3

         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
     

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(6)=ED    
c---------------------------------------------------------------------- 
                               
         XA=XM


         CALL FUNCX(EK,J,YA,XA,H)

C    COMPUTE K4, L4, M4, N4, O4

         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(6)=ED

c---------------------------------------------------------------------
                 

         XA=XP

         CALL FUNCX(EK,J,YA,XA,H)

C    4-TH ORDER RUNGE KUTTA SCHEME

         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
          P0=Y(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon density parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        Y(6)=ED
c---------------------------------------------------------------------  
                   
          
       END DO




       PS=Y(1)
         
c       PT=2.863D-1 
c       PMIN=1.0D-9 
      
       IF (PS .GT. PMIN  ) GOTO 28

C --------------------------------------------------
C The initial condition from omet and kapat are not
C satisfied but since: 
C   omet'(r) is a linear function of kapat and 
C   kapat'(r) is a linear function of omet,
C then if zeta is a constant, changing both
C omet -> zeta*omet and kapat -> zeta*kapat
C will not change each EoM. Suppose at r=R we
C obtain omet(R)=(1-2GI/R^3)/zeta and
C kapat(R)=GI/zeta, then the boundary conditions
C will be satisfied by changing both 
C omet(R) -> zeta*omet(R) and kapat -> zeta*kapat(R)
C with zeta=1/(omet(R)+2 kapat(R)/R^3)

C----------------------------------------------------
C In subroutine FUNCMI, I had implemented zeta
C
C      ROS=(XP/1.D3)
C      GMOS=Y(2)
C      MNURT=Y(3)
      OMET=Y(4)
      KAPAT=Y(5)
      ZETA=1.D0/(OMET+2.D0*KAPAT/(XP*XP*XP))
C      OMEGA=ZETA*OMET
C      KAPPA=ZETA*KAPAT


       CALL FUNCMI(XP,Y,MPHY,OMTPHY,KPTPHY)   

 
      
      ROS=(XP/1.D3)
      GMOS=MPHY
      MNURT=Y(3)
      OMEGA=OMTPHY
      KAPPA=KPTPHY
      
      IL=PCC
      MOMINEFF=ZETA*KAPAT/(GS*Y(2)*MSS*XP*XP)
      MOMINPHY=KPTPHY/(GS*GMOS*MSS*XP*XP)
      WRITE(*,*) Y(2),GMOS,MOMINEFF,MOMINPHY

   
      RETURN
      END


C--------------------------------------------------------------
C  TOV initial condition for nu
C----------------------------------------------------------

      SUBROUTINE TOVI(PCC,TUCT,ROS,GMOS,MNURT)
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IN, IL
      DIMENSION YA(10), EK(5,10), Y(10)

      HC  = 197.327D0
      PI  = 3.14159265358979D0  
      GS=1.325D-12
      MSS=1.1155D15

C     IM = NUMBER OF EQUATIONS
C------------------------------------------
C     IN =  NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS
C     Y(1)=Pressure, 
C     Y(2)=Star Mass,
C     Y(3)=Matric Nu,  
C----------------------------------------------
C     Y(4)=Energy Density, 


      IM=4
      IN=IM-1  
      
      XP=1.0D-3

      Y(1)=PCC
      Y(2)=0.1D-8
      Y(3)=TUCT

      
c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      
    
      Y(4)=FED(PCC)
     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D-1
      NS=8   ! recommended 16
      XL=30.0D3

      H=PU/NS
      XP=1.0D-3
      HH=H/(2.0D0)

C     LINE NUMBER INITIALIZATION

      LI=0
 

 28   LI=LI+1

C     XB=OLD radius, XP=NEW radius, XM=MIDPOINT radius

      DO N=1,NS
         XB=XP
         XP=XP+H
         XM=XB+HH

C    COMPUTE K1, L1, M1

         J=1
         DO I=1,IM
            YA(I)=Y(I)
         END DO
         XA=XB
   
         CALL FUNCT(EK,J,YA,XA,H)
         
C    COMPUTE K2, L2, M2

         J=2
         
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
         
        P0=YA(1)
c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(4)=ED
c---------------------------------------------------------------------

        
         XA=XM

         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K3, L3, M3

         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
     

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(4)=ED 
                 
c---------------------------------------------------------------------- 
           
         XA=XM

         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K4, L4, M4

         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(4)=ED

c---------------------------------------------------------------------
                

         XA=XP

         CALL FUNCT(EK,J,YA,XA,H)

C    4-TH ORDER RUNGE KUTTA SCHEME

         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
          P0=Y(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon density parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        Y(4)=ED
c---------------------------------------------------------------------  
                  
          
       END DO

  

       PS=Y(1)
       PMIN=1.0D-9
c       PMIN=2.0D-5
     
      IF (PS .GT. PMIN  ) GOTO 28
      
      CALL FUNCM(XP,Y,MPHY)

      ROS=(XP/1.D3)
      GMOS=MPHY
      MNURT=Y(3)
      
      IL=PCC
      WRITE(*,*) IL
      WRITE(*,*) Y(2),GMOS
           
      
      RETURN
      END


       SUBROUTINE FUNCX(EK,J,YA,XA,H)
C     *********************************************************
C     DEFINES TOV EQUATIONS AIP BL 
C     P-unisotropic used to calculate Moment of inertia
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(5,10), YA(10)
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
C ----------------------------------
c    model parameter !! Mmax=2.08 Ms
  
      !LBL=1.82D0   
      LBL=0.D0 
      K=15.6D6
      L=1.D0
C--------------------------------------       
      PRESS=YA(1)
      MASST=YA(2)
      MNU=YA(3)
      OMGT=YA(4)
      KAPT=YA(5)
      EDEN=YA(6) 
      EXPMNU=EXP(MNU/2.D0)
      EXPMMNU=EXP(-MNU/2.D0)
      
C sigma anisotropik pressure
      SGM=-(LBL*GS*EDEN*EDEN*XA*XA/3.D0)
     &        *(1.D0+3.D0*PRESS/EDEN)*(1.D0+PRESS/EDEN)
     &        /(1.D0-2.D0*GS*MASST*MSS/XA)
     
      A=DSQRT(L+8.D0*PI*GS*K*EDEN)
      B=DSQRT(L-8.D0*PI*GS*K*PRESS)
      C=DSQRT(L-8.D0*PI*GS*K*(PRESS-SGM))
      
      EDENEFF=(A*A-B*B-2.D0*C*C+2.D0*A*B*C*C)
     &       /(16.D0*PI*GS*K*A*B*C*C)
      PRESSEFF=(A*A-B*B+2.D0*C*C-2.D0*A*B*C*C)
     &        /(16.D0*PI*GS*K*A*B*C*C)
    
      DSGMDP=-(LBL*GS*EDEN*EDEN*XA*XA/3.D0)
     &        *((3.D0/EDEN)*(1.D0+PRESS/EDEN)
     &        +(1.D0+3.D0*PRESS/EDEN)*(1.D0/EDEN))
     &        /(1.D0-2.D0*GS*MASST*MSS/XA)
      
C Persamaan pressure
      SOS=1.D0/DEDP(PRESS)
      XXX=(1.D0/(B*C*C)+B/(A*A*C*C)-2.D0/(A*A*B))
     &       *1.D0/(4.D0*A*SOS)
     &      +(A/(B*B*C*C)+1.D0/(A*C*C)+2.D0/(A*B*B))
     &       *1.D0/(4.D0*B)
     &      +(A/B-B/A)*1.D0/(2.D0*C*C*C*C)
     &       *(1.D0+DSGMDP)
     
      EK(J,1)=H*(-EDENEFF*(1.D0+PRESSEFF/EDENEFF)/2.D0
     &        *(2.D0*GS*MASST*MSS/(XA*XA)
     &        *(1.D0+4.D0*PI*XA*XA*XA*PRESSEFF/(MSS*MASST))
     &        /(1.D0-2.D0*GS*MASST*MSS/XA))
     &        -(2.D0/XA)*SGM/(A*B*C*C))
     &        /XXX
     
C Persamaan massa
      EK(J,2)=H*4.D0*PI*XA*XA*EDENEFF/MSS
      
C Persamaan nu
      EK(J,3)=H*2.D0*GS*MASST*MSS/(XA*XA)
     &        *(1.D0+4.D0*PI*XA*XA*XA*PRESSEFF/(MSS*MASST))
     &        /(1.D0-2.D0*GS*MASST*MSS/XA)

      SGMEFF=SGM/(A*B*C*C)

C Persamaan omega tilde = omet
      EK(J,4)=H*6.D0*EXPMNU/DSQRT(1.D0-2.D0*GS*MASST*MSS/XA)
     &     *KAPT/(XA*XA*XA*XA)
     
     
C Persamaan kappa tilde = kapat
      EK(J,5)=H*(8.D0*PI*GS/3.D0)*(XA*XA*XA*XA)*EXPMMNU
     &        *EDENEFF*(1.D0+PRESSEFF/EDENEFF)
     &        *(1.D0-SGMEFF
     &          /(EDENEFF*(1.D0+PRESSEFF/EDENEFF)))
     &        *OMGT/DSQRT(1.D0-2.D0*GS*MASST*MSS/XA)
 
      
      RETURN
      END


      SUBROUTINE FUNCT(EK,J,YA,XA,H)
C     *********************************************************
C     DEFINES TOV EQUATIONS AIP BL 
C     P-anisotropic used to calculate initial condition for nu
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(5,10), YA(10)
      HC  = 197.327D0
      PI  = 3.14159265358979D0     
      GS=1.325D-12
      MSS=1.1155D15
C ----------------------------------
c    model parameter !! Mmax=2.08 Ms
  
      !LBL=1.82D0   
      LBL=0.D0 
      K=15.6D6
      L=1.D0
c--------------------------------------------      
      PRESS=YA(1)
      MASST=YA(2)
      MNU=YA(3)
      EDEN=YA(4) 
      
C sigma anisotropik pressure
      SGM=-(LBL*GS*EDEN*EDEN*XA*XA/3.D0)
     &        *(1.D0+3.D0*PRESS/EDEN)*(1.D0+PRESS/EDEN)
     &        /(1.D0-2.D0*GS*MASST*MSS/XA)
     
      A=DSQRT(L+8.D0*PI*GS*K*EDEN)
      B=DSQRT(L-8.D0*PI*GS*K*PRESS)
      C=DSQRT(L-8.D0*PI*GS*K*(PRESS-SGM))
      
      EDENEFF=(A*A-B*B-2.D0*C*C+2.D0*A*B*C*C)
     &       /(16.D0*PI*GS*K*A*B*C*C)
      PRESSEFF=(A*A-B*B+2.D0*C*C-2.D0*A*B*C*C)
     &        /(16.D0*PI*GS*K*A*B*C*C)
    
      DSGMDP=-(LBL*GS*EDEN*EDEN*XA*XA/3.D0)
     &        *((3.D0/EDEN)*(1.D0+PRESS/EDEN)
     &        +(1.D0+3.D0*PRESS/EDEN)*(1.D0/EDEN))
     &        /(1.D0-2.D0*GS*MASST*MSS/XA)
      
C Persamaan pressure
      SOS=1.D0/DEDP(PRESS)
      XXX=(1.D0/(B*C*C)+B/(A*A*C*C)-2.D0/(A*A*B))
     &       *1.D0/(4.D0*A*SOS)
     &      +(A/(B*B*C*C)+1.D0/(A*C*C)+2.D0/(A*B*B))
     &       *1.D0/(4.D0*B)
     &      +(A/B-B/A)*1.D0/(2.D0*C*C*C*C)
     &       *(1.D0+DSGMDP)
     
      EK(J,1)=H*(-EDENEFF*(1.D0+PRESSEFF/EDENEFF)/2.D0
     &        *(2.D0*GS*MASST*MSS/(XA*XA)
     &        *(1.D0+4.D0*PI*XA*XA*XA*PRESSEFF/(MSS*MASST))
     &        /(1.D0-2.D0*GS*MASST*MSS/XA))
     &        -(2.D0/XA)*SGM/(A*B*C*C))
     &        /XXX
     
C Persamaan massa
      EK(J,2)=H*4.D0*PI*XA*XA*EDENEFF/MSS
      
C Persamaan nu
      EK(J,3)=H*2.D0*GS*MASST*MSS/(XA*XA)
     &        *(1.D0+4.D0*PI*XA*XA*XA*PRESSEFF/(MSS*MASST))
     &        /(1.D0-2.D0*GS*MASST*MSS/XA)
      
C Semua persamaan di EK harus dikali H

   
      RETURN
      END

 
C---------------------------------------------------------------------
C  Energy density as a function of pressure
C------------------------------------------------------------------
      FUNCTION FED(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)

        IF ( P0 .GT. 50.D0 ) THEN


        FED=752.1779439721782D0 - 40.108620203769D0*P0 + 
     -  1.4752000864096837D0*P0**2 - 
     -  0.027036299898738337D0*P0**3 + 
     -  0.0003125837415657684D0*P0**4 - 
     -  2.441653734444954D-6*P0**5 + 
     -  1.3308776414303D-8*P0**6 - 
     -  5.1246996017212506D-11*P0**7 + 
     -  1.3886907993633848D-13*P0**8 - 
     -  2.590695651603915D-16*P0**9 + 
     -  3.166994006343349D-19*P0**10 - 
     -  2.283247911494129D-22*P0**11 + 
     -  7.357464623695704D-26*P0**12



  

        ELSE IF ( P0 .GT.  2.863D-1 .AND. P0 .LE. 50.D0 ) THEN


        FED=65.94376092512762D0 + 57.952630512664356D0*P0 - 
     -  15.437797576854498D0*P0**2 + 
     -  2.7511762122977346D0*P0**3 - 
     -  0.299643408613215D0*P0**4 + 
     -  0.02089083577091044D0*P0**5 - 
     -  0.0009679766327308126D0*P0**6 + 
     -  0.000030437894581131215D0*P0**7 - 
     -  6.521271499180282D-7*P0**8 + 
     -  9.372308503275227D-9*P0**9 - 
     -  8.643173476154485D-11*P0**10 + 
     -  4.621195967220976D-13*P0**11 - 
     -  1.0889239514378402D-15*P0**12
     
        ELSE IF (P0 .GT. 4.99313436D-4 .AND. P0 .LE. 2.863D-1) THEN


        FED=0.05015663787134234D0 + 836.2363942486941D0*P0 - 
     -  9315.146969977652D0*P0**2 + 79689.1930322726D0*P0**3 - 
     -  412197.6475732246D0*P0**4 + 1.116366190255507D6*P0**5 - 
     -  1.1988188657021397D6*P0**6


        ELSE
        
        FED=0.00020663104786863406D0 + 985.7550962048012D0*P0 - 
     -  6.452649548410687D6*P0**2 + 4.045493650683396D10*P0**3 - 
     -  1.2422017897554384D14*P0**4 + 1.765493095354931D17*P0**5 - 
     -  9.15294170121406D19*P0**6

        
        END IF

   
      RETURN
      END


      SUBROUTINE FUNCM(XA,YA,MPHY)
C     *********************************************************
C     Definition for physical mass MPHY
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)   
      DIMENSION YA(10)
      HC  = 197.327D0
      PI  = 3.14159265358979D0     
      GS=1.325D-12
      MSS=1.1155D15
C ----------------------------------
      !LBL=1.82D0   
      LBL=0.D0 
      K=15.6D6
      L=1.D0
      
      PRESS=YA(1)
      EDEN=YA(4)
      MASST=YA(2)
      
      SGM=-(LBL*GS*EDEN*EDEN*XA*XA/3.D0)
     &        *(1.D0+3.D0*PRESS/EDEN)*(1.D0+PRESS/EDEN)
     &        /(1.D0-2.D0*GS*MASST*MSS/XA)

      A=DSQRT(L+8.D0*PI*GS*K*EDEN)
      B=DSQRT(L-8.D0*PI*GS*K*PRESS)
      C=DSQRT(L-8.D0*PI*GS*K*(PRESS-SGM))
      
      MPHY=XA/(2.D0*GS*MSS*DSQRT(A*B))
     &     *(1.D0+A*C*C/(B*B)*(2.D0*GS*MASST*MSS/XA-1.D0))
     

   
      RETURN
      END


      SUBROUTINE FUNCMI(XA,YA,MPHY,OMTPHY,KPTPHY)
C     *********************************************************
C     Definition for physical mass MPHY and moment of inertia MIPHY
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)   
      DIMENSION YA(10)
      EXTERNAL FED,DEDP
      
      HC  = 197.327D0
      PI  = 3.14159265358979D0     
      GS=1.325D-12
      MSS=1.1155D15
C ----------------------------------
      !LBL=1.82D0   
      LBL=0.D0 
      K=15.6D6
      L=1.D0

      PRESS=YA(1)
      MASST=YA(2)
      MNU=YA(3)
      OMET=YA(4)
      KAPAT=YA(5)
      EDEN=YA(6)
      
      
      SGM=-(LBL*GS*EDEN*EDEN*XA*XA/3.D0)
     &        *(1.D0+3.D0*PRESS/EDEN)*(1.D0+PRESS/EDEN)
     &        /(1.D0-2.D0*GS*MASST*MSS/XA)

      A=DSQRT(L+8.D0*PI*GS*K*EDEN)
      B=DSQRT(L-8.D0*PI*GS*K*PRESS)
      C=DSQRT(L-8.D0*PI*GS*K*(PRESS-SGM))
      
      MPHY=XA/(2.D0*GS*MSS*DSQRT(A*B))
     &     *(1.D0+A*C*C/(B*B)*(2.D0*GS*MASST*MSS/XA-1.D0))

C ----------------------------------
C To convert omega & kappa to physical ones,
C I need to fix them by zeta first
      ZETA=1.D0/(OMET+2.D0*KAPAT/(XA*XA*XA))
      
      
      SOS=1.D0/DEDP(PRESS)
      XXX=(1.D0/(B*C*C)+B/(A*A*C*C)-2.D0/(A*A*B))
     &       *1.D0/(4.D0*A*SOS)
     &      +(A/(B*B*C*C)+1.D0/(A*C*C)+2.D0/(A*B*B))
     &       *1.D0/(4.D0*B)
     &      +(A/B-B/A)*1.D0/(2.D0*C*C*C*C)
     &       *(1.D0+DSGMDP)
      PERSP=-(-EDENEFF*(1.D0+PRESSEFF/EDENEFF)/2.D0
     &        *(2.D0*GS*MASST*MSS/(XA*XA)
     &        *(1.D0+4.D0*PI*XA*XA*XA*PRESSEFF/(MSS*MASST))
     &        /(1.D0-2.D0*GS*MASST*MSS/XA))
     &        -(2.D0/XA)*SGM/(A*B*C*C))
     &        /XXX
     
      DSGMDP=-(LBL*GS*EDEN*EDEN*XA*XA/3.D0)
     &        *((3.D0/EDEN)*(1.D0+PRESS/EDEN)
     &        +(1.D0+3.D0*PRESS/EDEN)*(1.D0/EDEN))
     &        /(1.D0-2.D0*GS*MASST*MSS/XA)
      PERSSGM=DSGMDP*PERSP
      
      MIA=KAPAT/GS*ZETA
      
      SOS=1.D0/DEDP(PRESS)
      
      MIPHY=MIA/(B*B*C*C)
     &     +4.D0*PI*K*XA*XA*XA*XA/3.D0
     &     *OMET*ZETA
     &     *(1.D0/(SOS*A*A*B*B*C*C)*PERSP
     &     +1.D0/(B*B*C*C*C*C)*(PERSP-PERSSGM))
           
      OMTPHY=A*A/(C*C)*ZETA*OMET 
      KPTPHY=GS*MIPHY
      
      
C      !Cek hasilnya

   
      RETURN
      END
      
C-----------------------------------------------------------------------
C     DE/DP AS A FUNCTION OF (P0) for NS (ANTO'S VERSION)
C-----------------------------------------------------------------------
      FUNCTION DEDP(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL FED
 60   h  = 1.D-6
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      C1=FED(x3)
      C2=FED(x4)
      MKK = (FED(x4)-8.D0*FED(x2)+8.D0*FED(x1)-FED(x3))/(12.D0*h)
      IF (ABS(MKK) .GE. 1.D-2/h) THEN  ! very rough
         xa=xa-(xa*h*1.D3)
         GOTO 60   
      ENDIF
      DEDP = MKK
c      WRITE(*,*)xa,DEDP
      RETURN
      END
 
c---------------------------------------------------------------------
C end of the code
