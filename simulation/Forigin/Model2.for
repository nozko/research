C     MODEL2.FOR
C    
      PARAMETER (ipMX=30000,npmx=200)
      COMMON BT(ipmx),ID(ipmx),C(ipmx),NSD(ipmx),NQD(ipmx),
     +   NPN(ipmx),CL(ipmx),NSDO(50),NSDN(50),btt(ipmx),nsd1(ipmx)
      COMMON AS(ipmx),QQ(ipmx),IOPT(ipmx),IOPQ(ipmx),IOPE(ipmx)
      COMMON NP(ipmx,npmx),NKR(ipmx,npmx),HR(ipmx,npmx),NK(ipmx,npmx)
      DIMENSION NFN(npmx)
      CHARACTER FNAME*20,FNAMP(ipmx)*20,Fip(ipmx)*12,
     +          FNL(npmx)*20,NFP(ipmx,npmx)*20
C     
      WRITE(*,*)'file name for input'
      READ(*,*)FNAME
      OPEN(7,FILE='file.net')
      OPEN(8,FILE='model.dat')
      OPEN(3,FILE=FNAME)
C
      IS=1
      LF=0
    5 READ(3,*)FNAME,A
      IF (FNAME.NE.'*') THEN
        read(3,*)ns
        DO 7 I=1,NS
          READ(3,*)NSDO(I),NSDN(I)
    7   CONTINUE
        WRITE(*,*)FNAME
        LF=LF+1
        FNL(LF)=FNAME
        NFN(LF)=IS
        WRITE(8,*)FNAME,IS
        OPEN(2,FILE=FNAME)
        READ(2,*)NX
        DO 30 I=IS,IS+NX-1
          WRITE(*,*)I
          FNAMP(I)=FNAME
          READ(2,*)FIP(I),BT(I),ID(I),C(I),CL(I)
          READ(2,*)NSD(I),AS(I),NQD(I),QQ(I)
          C(I)=A*C(I)
          CL(I)=A*CL(I)
          AS(I)=A*AS(I)
          QQ(I)=A*QQ(I)
          READ(2,*)IOPT(I),IOPQ(I),IOPE(I)
          READ(2,*)NPN(I)
          DO 40 J=1,NPN(I)
            READ(2,*)NP(I,J),NKR(I,J),HR(I,J)
            NFP(I,J)='NUL'
            HR(I,J)=A*HR(I,J)
            NK(I,J)=1
            IF (NP(I,J).GT.0) THEN
              NP(I,J)=NP(I,J)+IS-1
            END IF
   40     CONTINUE
          nsd1(i)=nsd(i)
          DO 9 J=1,NS
            if (NSD(I).EQ.NSDO(J)) then
              write(*,*)nsd(i),NSDN(J)
              nsd1(i)=NSDN(J)
            end if
    9     CONTINUE
          nsd(i)=nsd1(i)
   30   CONTINUE
        IS=IS+NX
        CLOSE(2)
        GOTO 5
      END IF
      IPX=IS
C
   60 READ(3,*)FNAME
      IF (FNAME.EQ.'*') GOTO 50
      WRITE(*,*)FNAME
      IP=IPX
      WRITE(*,*)IP
      LF=LF+1
      FNL(LF)=FNAME
      NFN(LF)=IP
      WRITE(8,*)FNAME,IP
      FNAMP(IP)=FNAME
      READ(3,*)ID(IP)
      IF (ID(IP).EQ.1) GOTO 25
      READ(3,*)BT(IP),NQD(IP)
      READ(3,*)A,D,C1
      AS(IP)=A*D
      C(IP)=A*D*C1
      CL(IP)=0
      READ(3,*)NPN(IP)
      DO 15 I=1,NPN(IP)
      READ(3,*)NFP(IP,I),A,B,IDX,KD
      write(8,*)'     ',nfp(ip,i),idx,kd
C      if (IPO.lt.0) ipo=ip+ipo
C      IF (A.LT.0) THEN
C        A=AS(IPO)
C        NPN(IPO)=NPN(IPO)+1
C        NP(IPO,NPN(IPO))=IP
C        NKR(IPO,NPN(IPO))=1
C        HR(IPO,NPN(IPO))=A*B
C      END IF
C      NP(IP,I)=IPO
      NKR(IP,I)=IDX         
      HR(IP,I)=A*B
      NK(IP,I)=KD
   15 CONTINUE
C
   75 READ(3,*)HC,NNSD,KOD
      IF (HC.EQ.0) GOTO 25
        DO 22 I=1,IPX
          IF (NSD(I).EQ.NNSD) THEN
             WRITE(8,*)NNSD,'  ',FNAMP(I)
             NPN(IP)=NPN(IP)+1
             NP(IP,NPN(IP))=I
             NKR(IP,NPN(IP))=KOD
             HR(IP,NPN(IP))=AS(I)*HC
             NPN(I)=NPN(I)+1
             NP(I,NPN(I))=IP
             NKR(I,NPN(I))=1
             HR(I,NPN(I))=AS(I)*HC
             NFP(IP,NPN(IP))='NUL'
             NFP(I,NPN(I))='NUL'
          END IF
   22   CONTINUE
      GOTO 75
C
   25 IPX=IPX+1
      GOTO 60
C
   50 ipx=ipx-1
   51 READ(3,*)FNAME
      IF (FNAME.EQ.'*') GOTO 70
      OPEN(4,FILE=FNAME)
      READ(4,*)NX,NXJ
      DO 80 I=1,NX
        READ(4,*)IPP
        J0=NPN(IPP)
        DO 81 J=1,NXJ
          READ(4,*)NPP,HH
          J0=J0+1
          NP(IPP,J0)=NPP
          NKR(IPP,J0)=3
          HR(IPP,J0)=HH
          NFP(IPP,J0)='NUL'
   81   CONTINUE
        NPN(IPP)=J0
   80 CONTINUE
      CLOSE(4)
      GOTO 51
c
   70 WRITE(*,*)'point no. to output temp'
   85 IP=0
      READ(3,*)FNAME
	DO 86 K=1,LF
	  IF (FNL(K).EQ.FNAME) IP=NFN(K)
   86 CONTINUE
      WRITE(*,*)FNAME,IP
      IF (IP.GT.0) THEN 
        IOPT(IP)=1
	  GOTO 85
      END IF
      WRITE(*,*)'point no. to output cal'
   90 IP=0
      READ(3,*)FNAME
	DO 91 K=1,LF
	  IF (FNL(K).EQ.FNAME) IP=NFN(K)
   91 CONTINUE
      WRITE(*,*)FNAME,IP
      IF (IP.GT.0) THEN
        IOPQ(IP)=1
        GOTO 90
      END IF
      WRITE(*,*)'point no. to output load'
   95 IP=0
      READ(3,*)FNAME
	DO 96 K=1,LF
	  IF (FNL(K).EQ.FNAME) IP=NFN(K)
   96 CONTINUE
      WRITE(*,*)FNAME,IP
      IF (IP.GT.0) THEN 
        IOPE(IP)=1
	  GOTO 95
      END IF
      write(*,*)'NSD to output'
   97 READ(3,*)NSDP
      write(*,*)nsdp
      IF (NSDP.GT.0) THEN
        DO 98 IP=1,IPX
          IF (NSDP.EQ.NSD(IP)) THEN
            IOPT(IP)=1
            IOPQ(IP)=1
          END IF
   98   CONTINUE
        GOTO 97
      END IF
C
      read(3,*)FNAME
      write(*,*)fname
      if (FNAME.ne.'*') then
        OPEN(2,FILE=FNAME)
        READ(2,*)NXX
        DO 106 IP=1,NXX
          READ(2,*)iii,BT(IP),DM1,DM2,dm3
          READ(2,*)DM3,DM4,DM5,DM6
          read(2,*)ii1,ii2,ii3
          READ(2,*)NPND
          DO 107 J=1,NPND
  	        READ(2,*)n1,n2,dm,nd
  107     CONTINUE
  106   CONTINUE
      end if
C
  115 WRITE(7,*)IPX
      DO 100 IP=1,IPX
        WRITE(7,*)IP,BT(IP),ID(IP),C(IP),CL(IP)
        WRITE(7,*)NSD(IP),AS(IP),NQD(IP),QQ(IP)
    	WRITE(7,*)IOPT(IP),IOPQ(IP),IOPE(IP)
        WRITE(7,*)NPN(IP)
        DO 110 J=1,NPN(IP)
          IF (NFP(IP,J).NE.'NUL') THEN
            IF (NFP(IP,J).EQ.'B') THEN
                NP(IP,J)=IP-1
              ELSE
                IC=0
                DO 300 K=1,LF
                  IF (NFP(IP,J).EQ.FNL(K)) THEN
                     NP(IP,J)=NFN(K)
                     IC=1
                  END IF
  300           CONTINUE
                IF (IC.EQ.0) THEN
                  do 303 k=1,ipx
                    IF (NFP(IP,J).EQ.FIP(K)) THEN
                     NP(IP,J)=K
                     IC=1
                    END IF
  303             CONTINUE
                END IF
                IF (IC.EQ.0) WRITE(*,*)'ERR',ip,j,NFP(IP,J)
             END IF
          END IF
          IF (NP(IP,J).LT.0) THEN
            WRITE(*,*)'IP<0',IP,J,NP(IP,J)
  302       READ(*,*)FNAME
            IC=0
            DO 301 K=1,LF
                  IF (FNAME.EQ.FNL(K)) THEN
                     NP(IP,J)=NFN(K)
                     IC=1
                  END IF
  301       CONTINUE
            IF (IC.EQ.0) THEN
               WRITE(*,*)'ERR'
               GOTO 302
            END IF
          END IF
	  WRITE(7,*)NP(IP,J),NKR(IP,J),HR(IP,J),NK(IP,J)
  110   CONTINUE
  100 CONTINUE
      OPEN(9,FILE='NAMELST')
	WRITE(9,*)LF
	DO 120 K=1,LF
	  WRITE(9,*)FNL(K),NFN(K)
  120 CONTINUE
      STOP
      END
 
