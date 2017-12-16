C      SMAP2000.FOR
C
      DIMENSION T(5000),BT(5000),Q(5000),E(5000),TSET(5000),ID(5000),
     /         QQ(5000),NPN(5000),NSD(5000),NQD(5000),C(5000),CR(5000),
     /    AS(5000),IOPT(5000),IOPE(5000),IOPQ(5000),CL(5000),ICT(5000)
      DIMENSION TD(100),QD(100),ED(100),LMS(100),LDS(100),JDN(100),
     /    JTD(100),JQD(100),JED(100),QE(5000),BF(5000),HTR(5000)
      DIMENSION NP(5000,40),NKR(5000,40),HR(5000,40),
     /         S1(5000),S2(5000),TT(5000),ER(5000),FLR(5000),IL(5000)
      DIMENSION BW(5000),BS(5000),SS(5000),WW(5000),EV(5000),DPH(5000),
     /   BDPH(5000),MLT(5000),SAT(5000),BST(5000),NC(5000),smf(5000) 
      DIMENSION ITMP(24),ISLR(24),IPRE(24)
      DIMENSION IX0(24),IRDV(24),IRFV(24),IUNR(24),IFKO(24),IVEL(24)
      DIMENSION CFD(100),CSCA(100),NFD(100),NSCA(100),
     /          CTIME(5,5),NTIME(5,5),DPD(100),SSD(100),TSD(100) 
      COMMON /POSITION/RD,LD,PIDO,PKED
      character fn*64
C
      SSE=0
      RT=0
      IFDX=0
      ISCAX=0
      ETSX=0
      QF=0
      qf2=0
      IDR=2
      irt=0
      ssmft=0
      DO 140 I=1,100
        CFD(I)=0
        CSCA(I)=0
  140 CONTINUE
      DO 1 I=1,5
        DO 2 J=1,5
          CTIME(I,J)=0
    2   CONTINUE
    1 CONTINUE 
      RD=3.1416/180
      PIDO=43.04
      PKED=141.33
c ****************************************************
      PCTH=0.
      pctl=-0.2
      TVM=(PCTH+PCTL)/2
      CTH=1.67
c ****************************************************
      OPEN(2,FILE='file.net')
      OPEN(7,FILE='op.csv')
      OPEN(4,FILE='smap.run')
      OPEN(9,FILE='sekisetu.csv')
C      WRITE(*,*)'DTM,LPR0,CF'
      READ(4,*)DTM,LPR0,CF,CK,LPSX
      IDX=INT(1/DTM+.4)
      LPR=1
      DLT=.0001
C      WRITE(*,*)'start LD, stop LY,LD'
      READ(4,*)LST,LYE,LEND
      READ(4,*)ND
      DO 7 I=1,ND
        READ(4,*)LMS(I),LDS(I)
    7 CONTINUE 
      LD=LST-1
      LYR=0
c
      READ(4,*)MODE,fn
      READ(4,*)QB,WF,TEMP,I1,I2,TMOFF,iroff,irset
      READ(4,*)A1,B1,C1,A2,B2,C2,A3,B3,C3,TAOFF
      READ(4,*)FLO,HIX,AR0,CVEL,DPH0,YMID
c
      IF (MODE.EQ.1) THEN
        OPEN(3,FILE=fn,FORM='FORMATTED',
     +         ACCESS='DIRECT',RECL=80) 
      END IF
      IF (MODE.EQ.2) THEN
        OPEN(3,FILE=fn,FORM='FORMATTED',
     +         ACCESS='DIRECT',RECL=80) 
      END IF
      IF (MODE.EQ.3) THEN
        OPEN(3,FILE=fn) 
      END IF
      T00=0
      IF (MODE.EQ.0) THEN
        READ(4,*)T0,RX,SUN,RNT,pre,VEL,TME,DPH0,YMID
      END IF
      PCT=0
      BS0=DPH0*YMID
      CLOSE(4)
C
      ITN=0
      IQN=0
      IEN=0
      IDN=0
   10 READ(2,*)IPX 
      DO 20 IP=1,IPX
        READ(2,*)DM,BT(IP),ID(IP),C(IP),CL(IP)
        READ(2,*)NSD(IP),AS(IP),NQD(IP),QQ(IP)
        READ(2,*)IOPT(IP),IOPQ(IP),IOPE(IP)
        IF (IOPT(IP).EQ.1) THEN
          ITN=ITN+1
          JTD(ITN)=IP
        END IF
        IF (IOPQ(IP).EQ.1) THEN
          IQN=IQN+1
          JQD(IQN)=IP
        END IF
        IF (IOPE(IP).EQ.1) THEN
          IEN=IEN+1
          JED(IEN)=IP
        END IF          
        IF (NSD(IP).EQ.1) THEN
          AT=AT+AS(IP)
          BS(IP)=BS0
          DPH(IP)=DPH0
          BDPH(IP)=DPH0
          IDN=IDN+1
          JDN(IDN)=IP
        END IF  
        READ(2,*)NPN(IP)
          DO 30 J=1,NPN(IP)
            READ(2,*)NP(IP,J),NKR(IP,J),HR(IP,J)
   30     CONTINUE
        FLR(IP)=0
   20 CONTINUE
      CLOSE(2)
      ih=ipx+1     
C
      WRITE(7,600)'MT,MD,TM,T0,X0,SUN,RNT,VEL,PRE,FD,SCA,WA,EVA,SM,SE,
     &flow',
     &  (JTD(I),I=1,ITN),(JQD(I),I=1,IQN),(JED(I),I=1,IEN) 
  600 FORMAT(A55,',',200(I4,','))
      WRITE(9,601)'MT,MD,TM',(JDN(I),I=1,IDN),
     &             (JDN(I),I=1,IDN),(JDN(I),I=1,IDN)
  601 FORMAT(A8,',',200(I4,','))
C
C     *DAYLOOP
   40 LD=LD+1
      IF (LD.GT.365) LD=1
      IF (LD.EQ.LST) LY=LY+1 
      WRITE(*,*)LY,LD
C
C     *WEATHER DATA
      IF (MODE.EQ.1) THEN
        READ(3,500,REC=4*LD-3)(ITMP(I),I=1,24),MM,MT,MD
        READ(3,500,REC=4*LD-2)(IVEL(I),I=1,24),MM,MT,MD
        READ(3,500,REC=4*LD-1)(ISLR(I),I=1,24),MM,MT,MD
        READ(3,500,REC=4*LD)(IPRE(I),I=1,24),MM,MT,MD
  500   FORMAT(24I3,I4,2I2)
      END IF
      IF (MODE.EQ.2) THEN
        READ(3,510,REC=7*LD-6)(ITMP(I),I=1,24),MM,MT,MD,MY,N
        READ(3,510,REC=7*LD-5)(IX0(I),I=1,24),MM,MT,MD,MY,N
        READ(3,510,REC=7*LD-4)(IRDV(I),I=1,24),MM,MT,MD,MY,N
        READ(3,510,REC=7*LD-3)(IRFV(I),I=1,24),MM,MT,MD,MY,N
        READ(3,510,REC=7*LD-2)(IUNR(I),I=1,24),MM,MT,MD,MY,N
        READ(3,510,REC=7*LD-1)(IFKO(I),I=1,24),MM,MT,MD,MY,N
        READ(3,510,REC=7*LD)(IVEL(I),I=1,24),MM,MT,MD,MY,N
  510   FORMAT(24I3,3I2,2I1)
      END IF
      IDS=0
      DO 3 I=1,ND
        IF ((MT.EQ.LMS(I)).AND.(MD.EQ.LDS(I))) IDS=1
    3 CONTINUE
      IF (MODE.EQ.0) THEN
        MD=LD
        IDS=1
      END IF
C
C     *TIMELOOP
      DO 95 ITM=1,24
      T00=T0
      if (mode.eq.3) then
        read(3,*)MM,MT,MD,mhr,ta,vp,vel0,sun,pre,unr
        sun=sun/4.186*1000
        RNT=45
        IF (VP.GE.0) THEN
            x0=0.622*vp/(1013.25-vp)
            RH=vp*760./1013.25
            IF (UNR.GE.0) THEN
              RNT=4.88E-08*
     &            (Ta+273.16)**4*(1-.62*UNR/10)*(.49-.076*SQRT(RH))
            END IF
          ELSE
            X0=0.7*FUNX(Ta)
        END IF
        DO 33 I=1,ND
          IF ((MT.EQ.LMS(I)).AND.(MD.EQ.LDS(I))) IDS=1
   33   CONTINUE
      end if
      DO 96 IDT=1,IDX  
      TM=ITM-1+DTM*IDT
      IF (MODE.EQ.0) THEN
        IF (TM.GT.TME) THEN
          FS=0
          FR=0
          PRE=0
        END IF
        X0=RX*FUNX(T0)
      END IF
      IF (MODE.EQ.1) THEN
        Ta=ITMP(ITM)
        Ta=Ta/10-50
        VEL=IVEL(ITM)
        VEL0=VEL/10.
        SUN=ISLR(ITM)
        SUN=SUN/4.186*10
        RX=.7
        X0=RX*FUNX(T0)
        PRE=IPRE(ITM)
        pre=pre/10
        RNT=45
      END IF
      IF (MODE.EQ.2) THEN
        Ta=ITMP(ITM)
        Ta=Ta/10-50
        RDV=IRDV(ITM)
        RFV=IRFV(ITM)
        X0=IX0(ITM)
        X0=X0/10000
        UNR=IUNR(ITM)
        FKO=IFKO(ITM)
        VEL=IVEL(ITM)
        VEL0=VEL/10.
        FKO=FKO*22.5-180
        RH=760*X0/(X0+.622)
        RNT=4.88E-08*(Ta+273.16)**4*(1-.62*UNR/10)*(.49-.076*SQRT(RH))
        CALL SUNPO(TM,SH,CH,SA,CA)
        CALL QDS(SH,CH,SA,CA,RDV,RFV,SUN)
        pre=0
      END IF  
      T0=T00+(Ta-T00)*DTM*IDT
      VEL=VEL0*CVEL
      IF (T0.LT.0) THEN
          FS=PRE
          FR=0
      END IF
      IF ((T0.GE.0).AND.(T0.LT.2)) THEN
          FS=PRE*(2-T0)/2
          FR=PRE*T0/2
      END IF  
      IF (T0.GE.2) THEN
          FS=0
          FR=PRE
      END IF
      YG=1000/(.091*T0**2-1.81*T0+9.47)
      IF (YG.LT.50) YG=50
C
C     *SET 
      SCA=0
      ssmf=0
      FD=0
      WA=0
      SE=0
      EROT=0
      EVA=0
      ILP=0
      ERX=0
      L1=0
      L2=0
      if (idr.eq.2) then
        pmx=0
        ptl=0
        ptm=0
      end if
      if (pre.gt.0) then
        ptl=ptl+pre*dtm
        ptm=ptm+dtm
      end if
      if (pre.gt.pmx) pmx=pre
      if (a1.eq.99) then
        a2=pmx*2
        if (a2.lt.2) then
          a2=2
        end if
        b1=a2
        b2=a2
      end if          
      if (a1.eq.999) then
        a2=pmx+2.
c        if (a2.lt.3.) then
c          a2=3.
c        end if
        a3=a2
        b1=a2
        b2=a2
        b3=a2
      end if 
      TMR=TMR+DTM
      IF (PRE.GT.0) TMR=0
c      IF ((MT.LE.3).OR.(MT.GE.12)) THEN
      IF ((MT.LE.4).OR.(MT.GE.11)) THEN
        IF (PRE.GT.0) THEN
                LEVEL=1
          ELSE
                IF ((BS(I2)+BW(I2)).GT.0) THEN
                    LEVEL=2
                  ELSE
                    LEVEL=3
                END IF
        END IF
        IF (LEVEL.EQ.1) THEN
          TN=A1+B1*T0
          TF=TN+C1
        END IF
        IF (LEVEL.EQ.2) THEN
          TN=A2+B2*T0
          TF=TN+C2
        END IF
        IF (LEVEL.EQ.3) THEN
          TN=A3+B3*T0
          TF=TN+C3
        END IF
        IF (BT(i1).LT.TN) IDR=1
        IF (BT(i1).GT.TF) IDR=2
        IF (TMR.GT.TMOFF+0.0001) IDR=2
        IF (T0.GT.TAOFF) IDR=2
        if (idr.eq.2) irt=-1
        imo=0
        if (idr.eq.1) then
          irt=irt+1
          ia=mod(irt,irset*iroff)
          if ((ia.ge.0).and.(ia.lt.iroff)) imo=1
          if ((ia.ge.iroff).and.(ia.lt.iroff*2)) imo=2
          if ((ia.ge.iroff*2).and.(ia.lt.iroff*3)) imo=3
        end if
        NTIME(LEVEL,IDR)=NTIME(LEVEL,IDR)+1
      END IF      
      DO 50 IP=1,IPX
        E(IP)=0
        Q(IP)=0
        NC(IP)=0
        MLT(IP)=0
        ICT(IP)=0
        BST(IP)=SAT(IP)
        IF (ID(IP).EQ.2) THEN
          TSET(IP)=T0
          ICT(IP)=1
        END IF
        IF (ID(IP).EQ.9) THEN 
          TSET(IP)=BT(IP)
          ICT(IP)=1 
        END IF
        IF (NQD(IP).EQ.1) Q(IP)=QQ(IP)
        IF (NQD(IP).EQ.11) THEN
            IF (IDR.EQ.1) THEN
              IH=IP
              Q(IP)=QB
              ICT(IP)=0
            END IF
        END IF
        IF (IDR.EQ.1) then
          if (ID(IP).EQ.10) then
            q(ip)=QB
            EROT=EROT+Q(IP)
          END IF
          if (ID(IP).EQ.11) then
            q(ip)=QB*0.5
            if (imo.EQ.1) Q(IP)=QB*2
            EROT=EROT+Q(IP)
          END IF
          if (ID(IP).EQ.12) then
            q(ip)=QB*0.5
            if (imo.EQ.2) Q(IP)=QB*2
            EROT=EROT+Q(IP)
          END IF
          if (ID(IP).EQ.13) then
            q(ip)=QB*0.5
            if (imo.EQ.3) Q(IP)=QB*2
            EROT=EROT+Q(IP)
          END IF
        END IF
        DO 88 J=1,NPN(IP)
          IF (NKR(IP,J).EQ.11) THEN
            HR(IP,J)=0
            IF (IDR.EQ.1) THEN
              HR(IP,J)=WF
            END IF
            flow=hr(ip,j)
          END IF
   88   CONTINUE
        CR(IP)=C(IP)
        IF (CL(IP).GT.0) THEN
          IF (BT(IP).GT.TVM) THEN
              CR(IP)=C(IP)+CL(IP)
            ELSE
              CR(IP)=C(IP)+CL(IP)*.5
          END IF
          Q(IP)=FLR(IP)
          IF (IL(IP).EQ.1) THEN
            TIC=BT(IP)
            PIC=(PCTH-TIC)/(PCTH-PCTL)
            CR(IP)=C(IP)+CL(IP)*(1-pic)+CL(IP)*0.5*pic
     $             +CL(IP)*80/ABS(PCTH-PCTL)
          END IF
        END IF
        IF (NSD(IP).EQ.1) THEN
          AR=.8-30*BDPH(IP)
          IF (AR.LT.AR0) THEN
              AR=AR0
          END IF
          SAT(IP)=T0+(AR*SUN-.9*RNT)/(FUNA(VEL)+4)
          EV(IP)=0
          TS=BT(IP)
          LPS=0
          DPS=BDPH(IP)
          IF ((BS(IP)+FS).GT.0) THEN
              IF ((TS.GT.0).OR.(SAT(IP).GT.0)) THEN
                MLT(IP)=1
              END IF
          END IF
          IF ((BW(IP)+FR).GT.0) THEN
C              IF ((TS.LT.0).OR.(SAT(IP).LT.0)) THEN
              IF (TS.LT.0) THEN
                MLT(IP)=1
              END IF
          END IF
          IF (((BS(IP)+FS).GT.0).AND.((BW(IP)+FR).GT.0)) THEN
              MLT(IP)=1
          END IF
          HI=0
          IF (BDPH(IP).GT.0) THEN
              HI=BW(IP)/(1000-BS(IP)/BDPH(IP))
            ELSE
              HI=BW(IP)/1000
          END IF
          IF (HI.GT.HIX) THEN
              HI=HIX
          END IF
          WAT=BW(IP)+FR*DTM
C
   49     SM=0
          EV(IP)=0
          BF(IP)=0
          IF (MLT(IP).EQ.1) THEN
            DH=DPS-HI
            IF ((DH.LE.0).OR.(SAT(IP).GT.0)) THEN
              DH=0
              TSV=0
              EV(IP)=4*FUNA(VEL)*(FUNX(TSV)-X0)
              WAT=BW(IP)+(FR-EV(IP))*DTM
              IF (WAT.LT.0) THEN
                EV(IP)=BW(IP)/DTM+FR
                WAT=0
              END IF
            END IF
            HTRM=1/(1/(FUNA(VEL)+4)+DH/0.08)
            SM0=(200*TS+HTRM*SAT(IP)-590*EV(IP))*DTM/80
            SM=SM0
            BF(IP)=1
            IF (SM0.LT.-WAT) THEN
              SM=-WAT
              BF(IP)=SM/SM0
            END IF
            IF (SM0.GT.(BS(IP)+FS*DTM)) THEN
              SM=BS(IP)+FS*DTM 
              BF(IP)=SM/SM0
            END IF
          END IF
          IF (MLT(IP).EQ.0) THEN
            TSV=BT(IP)
            EV(IP)=4*FUNA(VEL)*(FUNX(TSV)-X0)
            IF (EV(IP).GT.(BW(IP)/DTM+FR)) THEN
              EV(IP)=BW(IP)/DTM+FR
            END IF
          END IF
          HTR(IP)=1/(1/(FUNA(VEL)+4)+DPS/0.08)
          QE(IP)=-590*EV(IP)*AS(IP)*(1-BF(IP))
          SS(IP)=BS(IP)-SM+FS*DTM
          WW(IP)=BW(IP)+SM+(FR-EV(IP))*DTM
          smf(ip)=SM
          IF (SS(IP).GT.0) THEN
              IF (SM.GT.0) THEN
                 G=(BS(IP)+FS*DTM)/(BDPH(IP)+FS*DTM/YG)
                 DPH(IP)=SS(IP)/G
                ELSE
                 DPH(IP)=BDPH(IP)+FS*DTM/YG-SM/916
              END IF
            ELSE
                DPH(IP)=0 
          END IF
          DDPS=(DPH(IP)+BDPH(IP))/2-DPS
          IF ((ABS(DDPS).GT.0.001).AND.(LPS.LT.10)) THEN
            DPS=0.7*DDPS+DPS
            LPS=LPS+1
            GOTO 49
          END IF
        END IF
   50 CONTINUE
C
C     *MATLIX
      DO 70 IP=1,IPX
        T(IP)=BT(IP)
   70 CONTINUE
   71 LPS=0
   76 DO 72 IP=1,IPX
          IF (ICT(IP).EQ.1) THEN 
            T(IP)=TSET(IP)
            E(IP)=0
          END IF
          S1(IP)=(Q(IP)+QE(IP)+E(IP))*DTM+BT(IP)*CR(IP)
          S2(IP)=CR(IP) 
          DO 74 J=1,NPN(IP)
              S1(IP)=S1(IP)+HR(IP,J)*(CF*T(NP(IP,J))+(1-CF)*
     &               (BT(NP(IP,J))-BT(IP)))*DTM
              S2(IP)=S2(IP)+CF*HR(IP,J)*DTM
   74     CONTINUE  
          IF (NSD(IP).EQ.1) THEN
            F=BF(IP)
            S1(IP)=S1(IP)+((1-F)*HTR(IP)*SAT(IP)-F*200*(1-CF)*BT(IP)
     &        -(1-F)*HTR(IP)*(1-CF)*BT(IP))*DTM*AS(IP)
            S2(IP)=S2(IP)+(F*200+(1-F)*HTR(IP))*DTM*AS(IP)*CF
          END IF  
          IF (ICT(IP).EQ.1) THEN
               E(IP)=(S2(IP)*TSET(IP)-S1(IP))/DTM
            ELSE
              IF (CF.LT.0.01) THEN
                  T(IP)=S1(IP)/S2(IP)
                  ERX=0
                ELSE
                  TT(IP)=S1(IP)/S2(IP)
                  ER(IP)=TT(IP)-T(IP)
                  AER=ABS(ER(IP))
                  IF (AER.GT.ERX) THEN 
                    ERX=AER
                    IER=IP
                  END IF  
                  T(IP)=CK*ER(IP)+T(IP)
              END IF
          END IF
C        if (e(ip).lt.0)  write(*,*)e(ip),idr,ilp
   72   CONTINUE
        LPS=LPS+1
        IF (LPS.GT.LPSX) THEN
          WRITE(*,*)MT,MD,TM,ERX,IER,T(IER)
          GOTO 77
        END IF
        IF (ERX.GT.DLT) THEN
          ERX=0
          GOTO 76
        END IF
        IF (LPS.GT.LPMX) LPMX=LPS
   77   IF (ILP.lt.2) THEN
          IF ((T(IH).GT.TEMP).and.(Q(IH).GT.0)) THEN
            ILP=ilp+1
            TSET(IH)=TEMP
            ICT(IH)=1
            Q(IH)=0
            ERX=0
            GOTO 71
          END IF  
          IF (E(IH).LT.0) THEN
            ILP=ilp+1
            ICT(IH)=0
            E(IH)=0
            ERX=0
            GOTO 71
          END IF
        END IF    
C
        DO 90 IP=1,IPX
        IF (CL(IP).GT.0) THEN
          IL(IP)=0
          FLR(IP)=0
          IF ((T(IP).LT.PCTH).AND.(T(IP).GT.PCTL)) THEN
            IL(IP)=1
          END IF
          IF ((BT(IP).GT.PCTH).AND.(T(IP).LT.PCTH)) THEN
            FLR(IP)=(C(IP)+CL(IP))*(T(IP)-PCTH)
            T(IP)=PCTH
            IL(IP)=1
          END IF
          IF ((BT(IP).GT.PCTL).AND.(T(IP).LT.PCTL)) THEN
            FLR(IP)=(C(IP)+CL(IP)*80/ABS(PCTH-PCTL))*(T(IP)-PCTL)
            T(IP)=PCTL
          END IF
          IF ((BT(IP).LT.PCTL).AND.(T(IP).GT.PCTL)) THEN
            FLR(IP)=(C(IP)+CL(IP)*.5)*(T(IP)-PCTL)
            T(IP)=PCTL
            IL(IP)=1
          END IF
          IF ((BT(IP).LT.PCTH).AND.(T(IP).GT.PCTH)) THEN
            FLR(IP)=(C(IP)+CL(IP)*80/ABS(PCTH-PCTL))*(T(IP)-PCTH)
            T(IP)=PCTH
          END IF
        END IF        
        IF (NSD(IP).EQ.1) THEN
          IF ((LEVEL.EQ.1).AND.(T(IP).LT.0)) THEN
            L1=1
          END IF
          IF ((LEVEL.EQ.2).AND.(T(IP).LT.0)) THEN
            L2=1
          END IF
          IF (SS(IP).LT.0) THEN
            SS(IP)=0
          END IF
          IF (WW(IP).LT.0) THEN
            WW(IP)=0
          END IF
          IF (DPH(IP).GT.FD) THEN
            FD=DPH(IP)
          END IF
          SCA=SCA+SS(IP)*AS(IP)/AT
          WA=WA+WW(IP)*AS(IP)/AT
          EVA=EVA+EV(IP)*DTM*AS(IP)/AT
          ssmf=ssmf+smf(ip)*AS(IP)/AT
          IF ((DPH(IP).GT.0).AND.(SS(IP).GT.0)) THEN
            GM=SS(IP)/DPH(IP)
            EE=16*EXP(21E-3*GM)
            GM=GM*EXP(SS(IP)/2/EE*DTM/24)
            IF (GM.GT.916) THEN
              GM=916
            END IF  
            DPH(IP)=SS(IP)/GM
          END IF  
        END IF
   90 CONTINUE
C
      IF (FD.GT.0) THEN
        IFD=FD*100
        DO 180 K=1,IFD+1
          NFD(K)=NFD(K)+1
  180   CONTINUE
        IF (IFD.GT.IFDX) THEN
          IFDX=IFD
          write(*,*)ifdx,ifd,fd
        END IF
      END IF
      IF (SCA.GT.0) THEN
        ISCA=SCA
        DO 185 K=1,ISCA+1
          NSCA(K)=NSCA(K)+1
  185   CONTINUE
        IF (ISCA.GT.ISCAX) THEN
          ISCAX=ISCA
        END IF
      END IF
      TL1=TL1+L1
      TL2=TL2+L2
      SE=E(IH)+Q(IH)+EROT
      IF (SE.GT.0) THEN
        RT=RT+1
        SSE=SSE+SE
        RSE=SSE/RT
      END IF  
      ssmft=ssmft+ssmf
C
  200 IF (IDS.GE.0) THEN
       IF (LPR.EQ.LPR0) THEN
        ITN=0
        IQN=0
        IEN=0
        IDP=0
        DO 170 IP=1,IPX
          IF (IOPT(IP).EQ.1) THEN
            ITN=ITN+1
            TD(ITN)=T(IP)
          END IF
          IF (IOPQ(IP).EQ.1) THEN
            IQN=IQN+1
            QD(IQN)=Q(IP)
          END IF
          IF (IOPE(IP).EQ.1) THEN
            IEN=IEN+1
            ED(IEN)=E(IP)
          END IF
          IF (NSD(IP).EQ.1) THEN
            IDP=IDP+1
            DPD(IDP)=DPH(IP)
            SSD(IDP)=SS(IP)
            TSD(IDP)=T(IP)
          END IF          
  170   CONTINUE
        WRITE(7,610)MT,MD,TM,T0,X0,SUN/.86,RNT/.86,
     &  VEL,PRE,FD,SCA,WA,EVA,ssmf,SE/.86/at,flow,
     &   (TD(I),I=1,ITN),(QD(I)/.86,I=1,IQN),(ED(I)/.86,I=1,IEN)
        WRITE(9,610)MT,MD,TM,(DPD(I),I=1,IDP),(SSD(I),I=1,IDP),
     &  (TSD(I),I=1,IDP)
  610   FORMAT(2(I2,','),200(E10.4,','))
        LPR=0
       END IF
       LPR=LPR+1
      END IF 
C
          IF ((PRE.GT.0).AND.(T(I2).GT.0)) THEN
            TS=CF*T(I2)+(1-CF)*BT(I2)
            QF2=QF2+(BF(I2)*200*TS+(1-BF(I2))*HTR(I2)*(TS-SAT(I2)))
          END IF
      DO 100 IP=1,IPX
        IF (NSD(IP).EQ.1) THEN
          IF (((SCA+PRE).EQ.0).AND.(WW(IP).GT.FLO)) THEN
            WW(IP)=FLO
          END IF
          BW(IP)=WW(IP)
          BS(IP)=SS(IP)
          BDPH(IP)=DPH(IP)
          IF ((PRE.GT.0).AND.(T(IP).GT.0)) THEN
            TS=CF*T(IP)+(1-CF)*BT(IP)
            QF=QF+(BF(IP)*200*TS+(1-BF(IP))*HTR(IP)*(TS-SAT(IP)))
     &         *AS(IP)
          END IF
        END IF
        BT(IP)=T(IP)
  100 CONTINUE        
   96 CONTINUE
   95 CONTINUE
C 
      IF ((LD.NE.LEND).OR.(LY.NE.LYE)) GOTO 40
C
      CLOSE(7)
      OPEN(7,FILE='hist.rst')
      RT=RT*DTM
      TL1=TL1*DTM
      TL2=TL2*DTM
      SSE=SSE*DTM
      QF=QF*DTM
      qf2=qf2*dtm
      DO 150 I=1,IFDX+1
        CFD(I)=NFD(I)*DTM
  150 CONTINUE
      DO 155 I=1,ISCAX+1
        CSCA(I)=NSCA(I)*DTM
  155 CONTINUE
      DO 157 I=1,3
        DO 158 J=1,2
          CTIME(I,J)=NTIME(I,J)*DTM
  158   CONTINUE
  157 CONTINUE   
      WRITE(7,650)RT,SSE/at*4.186,QF/at*4.186,qf2*4.186
     &            ,ssmft,at*1000
      write(7,650)((CTIME(I,J),I=1,3),J=1,2),TL1,TL2
      WRITE(7,650)(CFD(I),I=1,IFDX+1)
      WRITE(7,650)(CSCA(I),I=1,ISCAX+1)
      WRITE(7,655)DTM,LPR0,LST,LYE,LEND,MODE
      WRITE(7,650)FLO,HIX,AR0
      WRITE(7,650)QB/at/.86,TEMP,I1,I2,tmoff
      WRITE(7,650)A1,B1,C1,A2,B2,C2,A3,B3,C3,TAOFF
  650 FORMAT(100F10.1)        
  655 format(f10.2,5i10)
C
      OPEN(8,FILE='file2.net')
      WRITE(8,*)IPX 
      DO 5 IP=1,IPX
        WRITE(8,*)IP,BT(IP),ID(IP),C(IP),CL(IP)
        WRITE(8,*)NSD(IP),AS(IP),NQD(IP),QQ(IP)
        WRITE(8,*)IOPT(IP),IOPQ(IP),IOPE(IP)
        WRITE(8,*)NPN(IP)
        DO 6 J=1,NPN(IP)
          WRITE(8,*)NP(IP,J),NKR(IP,J),HR(IP,J),ndmy
    6   CONTINUE
    5 CONTINUE
      WRITE(*,*)'LPMX=',LPMX
      STOP
      END
C
C
      FUNCTION FUNX(T)
      P=51.7063*10**(6.21147-2886.37/(1.8*T+491.69)-337269.46/(1.8
     &  *T+491.69)**2)
      FUNX=.622*P/(760-P)
      RETURN
      END   
C
      FUNCTION FUNA(V)
      FUNA=5+3.4*V
      RETURN
      END
C      
      SUBROUTINE SUNPO(TM,SH,CH,SA,CA)
C     *SUN POSITION
      COMMON /POSITION/RD,LD,PIDO,PKED
      W=2*3.1416*LD/366
      DEL=.362213-23.2476*COS(W+.153231)-.336891*COS(2*W+.207099)
     \     -.185265*COS(3*W+.620129)
      AE=-2.78641E-04+.122772*COS(W+1.49831)-.165458*COS(2*W-1.26155)
     \   -5.35383E-03*COS(3*W-1.1571)
      SJIK=15*(TM-12+AE)+PKED-135
      SH=SIN(PIDO*RD)*SIN(DEL*RD)+COS(PIDO*RD)*COS(DEL*RD)*COS(SJIK*RD)
      IF (SH.LE.0) THEN
          SH=0
          CH=0
          SA=0
          CA=0
        ELSE
          CH=SQRT(1-SH**2)
          SA=COS(DEL*RD)*SIN(SJIK*RD)/CH
          CA=(SH*SIN(PIDO*RD)-SIN(DEL*RD))/CH/COS(PIDO*RD)
      END IF
      RETURN
      END
C
      SUBROUTINE QDS(SH,CH,SA,CA,RDN,RSH,RTV)
      COMMON /POSITION/RD,LD,PIDO,PKED
      WH=0
      WK=0
      SFE=(1-COS(WK*RD))/2
      SFA=1-SFE
      CG=CA*COS(WH*RD)+SA*SIN(WH*RD)
      IF (SH.GT.0) THEN 
        TG=SQRT(1/CG**2-1)
      END IF
      CT=SH*COS(WK*RD)+CH*SIN(WK*RD)*CG
      IF (CT.GT.0) THEN 
          RDV=RDN*CT 
        ELSE 
          RDV=0
      END IF
      RFV=SFE*.1*(RDN*SH+RSH)+SFA*RSH
      RTV=RDV+RFV
      RETURN
      END
C
