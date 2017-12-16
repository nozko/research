C      MD-2b.FOR
C
      DIMENSION IP(0:100,0:100),NP(5000,8),COND(5000,8),C(5000),
     +	      ID(5000),
     +          BT(5000),NSD(5000),AS(5000),NQD(5000),QQ(5000),
     +          NPN(5000),CL(5000),CK(6),KD(50),KDM(5000,8),IOPT(5000)
      DIMENSION IMN(1000),IMX(1000),XL(1000),XD(0:1000),NPP(10000,8),
     +          JMN(1000),JMX(1000),YL(1000),YD(0:1000),TI(1000)
      DIMENSION HC(1000),CP(1000),CBL(1000),RX(1000),RY(1000),nc(1000)
      DIMENSION IMNb(100),IMXb(100),XLb(100)
	dimension JMNb(100),JMXb(100),YLb(100)
      DIMENSION RXb(100),RYb(100),nb(100,100)
      CHARACTER FN1*20,FN2*20,fn*20,chr*5
C
      WRITE(*,*)'file name to input'
      READ(*,*)FN1
      OPEN(2,FILE=FN1)
      write(*,*)'file name(.num)'
      read(*,*)fn
      OPEN(8,FILE=fn)
      WRITE(*,*)'file name for output(+.mdl)'
      READ(*,*)FN2
      OPEN(7,FILE=FN2)
C
      READ(2,*)ZL
      read(2,*)ibx
      do 10 i=1,ibx
        read(2,*)imnb(i),imxb(i),xlb(i),rxb(i)
   10 continue
      read(2,*)jbx
      do 21 j=1,jbx
        read(2,*)jmnb(j),jmxb(j),ylb(j),ryb(j)
   21 continue
      read(2,*)chr
      do 22 j=1,jbx
        read(2,*)(nb(i,j),i=1,ibx)
   22 continue
      do 24 j=1,jbx
        do 26 i=1,ibx
        if (nb(i,j).gt.0) then
          ismx=ismx+1
          nc(ismx)=nb(i,j)
          imn(ismx)=imnb(i)
          imx(ismx)=imxb(i)
          xl(ismx)=xlb(i)
          rx(ismx)=rxb(i)
          jmn(ismx)=jmnb(j)
          jmx(ismx)=jmxb(j)
          yl(ismx)=ylb(j)
          ry(ismx)=ryb(j)
        end if
        if (nb(i,j).gt.ncx) ncx=nb(i,j)
   26   continue
   24 continue
      read(2,*)ncxx
      write(*,*)ncxx
      if (ncxx.lt.ncx) then
        write(*,*)'error  ctrl+c'
        read(*,*)a
      end if
      DO 30 I=1,ncxx
        READ(2,*)NN,KD(nn),HC(nn),CP(nn),CBL(nn)
   30 CONTINUE
      IX=IMXb(Ibx)
      JX=JMXb(jbx)
C
C
      DO 40 IS=1,ISMX
        DO 44 J=JMN(IS),JMX(IS)
          DO 46 I=IMN(IS),IMX(IS)
            IF (IP(I,J).EQ.0) THEN 
              NX=NX+1
              IP(I,J)=NX
            END IF  
   46     CONTINUE
   44   CONTINUE
   40 CONTINUE  
C
      DO 200 I=2,IX
        II=1
  205   IF ((I.GT.IMN(II)).AND.(I.LE.IMX(II))) THEN
            IF (RX(II).GT.0) THEN
                XD(I)=XD(I-1)+XL(II)*RX(II)**(I-IMN(II)-1)
              ELSE 
                XD(I)=XD(I-1)+XL(II)/(IMX(II)-IMN(II))
            END IF    
            GOTO 199
          ELSE
             II=II+1
             GOTO 205
        END IF
  199   WRITE(*,*)XD(I)
  200 CONTINUE
      DO 210 J=2,JX
        II=1
  215   IF ((J.GT.JMN(II)).AND.(J.LE.JMX(II))) THEN
            IF (RY(II).GT.0) THEN
                YD(J)=YD(J-1)+YL(II)*RY(II)**(J-JMN(II)-1)
              ELSE   
                YD(J)=YD(J-1)+YL(II)/(JMX(II)-JMN(II))
            END IF  
            GOTO 209
          ELSE
             II=II+1
             GOTO 215
        END IF
  209   WRITE(*,*)YD(J)      
  210 CONTINUE
C
C
      DO 50 IS=1,ISMX
          DO 54 J=JMN(IS),JMX(IS)
            DO 56 I=IMN(IS),IMX(IS)
              NP(IP(I,J),1)=IP(I,J-1)
              NP(IP(I,J),2)=IP(I-1,J)
              NP(IP(I,J),3)=IP(I,J+1)
              NP(IP(I,J),4)=IP(I+1,J)
              COND(IP(I,J),1)=0
              COND(IP(I,J),2)=0
              COND(IP(I,J),3)=0
              COND(IP(I,J),4)=0
              NPN(IP(I,J))=0
              IF (I.EQ.1) NP(IP(I,J),2)=0
              IF (I.EQ.IX) NP(IP(I,J),4)=0
              IF (J.EQ.1) NP(IP(I,J),1)=0
              IF (J.EQ.JX) NP(IP(I,J),3)=0
   56       CONTINUE
   54     CONTINUE
   50 CONTINUE
      DO 60 IS=1,ISMX
          DO 64 J=JMN(IS),JMX(IS)
            DO 66 I=IMN(IS),IMX(IS)
              DX=(XD(I+1)-XD(I-1))/2
              DY=(YD(J+1)-YD(J-1))/2
              IF (I.EQ.IMN(IS)) DX=(XD(I+1)-XD(I))/2
              IF (I.EQ.IMX(IS)) DX=(XD(I)-XD(I-1))/2
              IF (J.EQ.JMN(IS)) DY=(YD(J+1)-YD(J))/2
              IF (J.EQ.JMX(IS)) DY=(YD(J)-YD(J-1))/2
              V=DX*DY
              IF (J.NE.1) CK(1)=DX/(YD(J)-YD(J-1))
              IF (I.NE.1) CK(2)=DY/(XD(I)-XD(I-1))
              CK(3)=DX/(YD(J+1)-YD(J))
              CK(4)=DY/(XD(I+1)-XD(I))
              IF (I.EQ.IMN(IS)) CK(2)=0
              IF (I.EQ.IMX(IS)) CK(4)=0
              IF (J.EQ.JMN(IS)) CK(1)=0
              IF (J.EQ.JMX(IS)) CK(3)=0
              C(IP(I,J))=C(IP(I,J))+CP(nc(IS))*V
              CL(IP(I,J))=CL(IP(I,J))+CBL(nc(IS))*V
              DO 68 L=1,4
                IF (CK(L).GT.0) THEN
                  NPN(IP(I,J))=NPN(IP(I,J))+1
                  LL=NPN(IP(I,J))
                  COND(IP(I,J),LL)=CK(L)*HC(nc(IS))
                  KDM(IP(I,J),LL)=KD(nc(IS))
                  NPP(IP(I,J),LL)=NP(IP(I,J),L)
                END IF
   68         CONTINUE  
   66       CONTINUE
   64     CONTINUE
   60 CONTINUE
C
      WRITE(8,*)IX,JX 
      WRITE(8,*)(XD(I),I=1,IX)
      WRITE(8,*)(YD(J),J=1,JX)
        DO 72 J=1,JX
          WRITE(8,600)(IP(I,J),I=1,IX)
  600     FORMAT(100I5)
   72   CONTINUE
C
C      WRITE(*,*)'input boundary'
      READ(2,*)NA
      DO 20 IB=1,NA
      READ(2,*)I1,I2,J1,J2,IDD,NSDD
      IF (J1.EQ.J2) THEN
          J=J1
          DO 92 I=I1,I2
            DX=(XD(I+1)-XD(I-1))/2
            IF (I.EQ.I1) DX=(XD(I+1)-XD(I))/2
            IF (I.EQ.I2) DX=(XD(I)-XD(I-1))/2
            AS(IP(I,J))=AS(IP(I,J))+DX
            ID(IP(I,J))=IDD
            NSD(IP(I,J))=NSDD
   92     CONTINUE
      END IF               
      IF (I1.EQ.I2) THEN
          I=I1
          DO 96 J=J1,J2
            DY=(YD(J+1)-YD(J-1))/2
            IF (J.EQ.J1) DY=(YD(J+1)-YD(J))/2
            IF (J.EQ.J2) DY=(YD(J)-YD(J-1))/2
            AS(IP(I,J))=AS(IP(I,J))+DY
            ID(IP(I,J))=IDD
            NSD(IP(I,J))=NSDD
   96     CONTINUE
      END IF        
   20 CONTINUE
C
      READ(2,*)chr
    9 WRITE(*,*)'input opt'
      WRITE(*,*)'I1,I2,J1,J2'
      READ(2,*)I1,I2,J1,J2
      IF (I1.GT.0) THEN
        DO J=J1,J2
        DO I=I1,I2
          IOPT(IP(I,J))=1
        END DO
        END DO
      GOTO 9
      END IF
C
    7 READ(2,*)fn2
      IF (fn2.ne.'end') then
        WRITE(*,*)'input initial temp'
        DO J=1,JX
          READ(2,*)Ti(J)
        END DO
        DO 160 I=1,IX
          DO 160 J=1,JX
            IF (IP(I,J).GT.0) THEN
              BT(IP(I,J))=Ti(J)
            END IF
  160     CONTINUE
        else
         write(*,*)'initial temp.'
         read(*,*)t0
         BT(1:NX)=t0
      end if
C
   14 WRITE(7,*)NX
      DO 100 I=1,NX
        C(I)=C(I)*ZL
        CL(I)=CL(I)*ZL
        AS(I)=AS(I)*ZL
        QQ(I)=QQ(I)*ZL
        WRITE(7,*)I,BT(I),ID(I),C(I),CL(I)
        WRITE(7,*)NSD(I),AS(I),NQD(I),QQ(I)
        WRITE(7,*)Iopt(i),IDM,IDM
        WRITE(7,*)NPN(I)
          DO 102 J=1,NPN(I)
            COND(I,J)=COND(I,J)*ZL
            WRITE(7,*)NPP(I,J),KDM(I,J),COND(I,J),idm
  102     CONTINUE
  100 CONTINUE
      STOP
      END
