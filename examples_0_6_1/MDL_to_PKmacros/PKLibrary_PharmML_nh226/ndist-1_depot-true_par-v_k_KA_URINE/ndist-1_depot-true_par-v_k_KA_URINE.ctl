$PROB WARFARIN PK WITH URINE COLLECTION

$INPUT ID 
TIME ; units="h"
WT ; units="kg"
AGE ; units="y"
SEX
AMT ; units="mg"
CMT
DVID
DV ; units="mg/L"
MDV
$DATA ../Data/warfarin_urine.csv IGNORE=#

$EST METHOD=COND INTER 
MAX=9990 SIG=3 NOABORT ;PRINT=1
$COV

$THETA
(0.001,0.05) ; POP_CLNR  units="L/h/70kg"
(0.001,0.05) ; POP_CLR  units="L/h/70kg"
(0.001,8)   ; POP_V  units="L/70kg"
(0.001,2) ; POP_KA units="1/h"
(0.001,1) ; POP_TLAG units="h"
$OMEGA BLOCK(2)
0.1 ; PPV_CL
0.01 0.1 ; PPV_V
$OMEGA
0.1 ; PPV_KA
0.1 ; PPV_TLAG

$SIGMA
0.01 ; RUV_PROP
0.05 ; RUV_ADD units="mg/L"
1 ; RUV_ADDU units="mg"

$SUBR ADVAN2 TRANS2

$PK

   ; Covariate model

   GRPCL=(THETA(1) + THETA(2))*(WT/70)**0.75
   GRPV=THETA(3)*WT/70
   GRPKA=THETA(4)
   GRPLG=THETA(5)

   ; Individual parameters
   CL=GRPCL*EXP(ETA(1))
   V=GRPV*EXP(ETA(2))
   KA=GRPKA*EXP(ETA(3))
   TLAG=GRPLG*EXP(ETA(4))
   FURINE=THETA(2)/(THETA(1)+THETA(2))

   ;Translate to PK Library standard names
   S2=V
   FO=FURINE
   ALAG1=TLAG
   K=CL/V

$ERROR
   CONC=A(2)/V
   AMTU=A(3)
   IF (DVID.EQ.1) THEN
      Y=CONC*(1+ERR(1))+ERR(2)
   ELSE
      Y=AMTU+ERR(3)
   ENDIF

$TABLE ID TIME WT SEX AGE ; covariates
CL V KA TLAG ; EBE estimates
CMT Y ; predictions
ONEHEADER NOPRINT FILE=warf.fit
