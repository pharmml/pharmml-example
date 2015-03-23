
$PROB WARFARIN PK
;O'Reilly RA, Aggeler PM, Leong LS. Studies of the coumarin anticoagulant
;drugs: The pharmacodynamics of warfarin in man.
;Journal of Clinical Investigation 1963;42(10):1542-1551

;O'Reilly RA, Aggeler PM. Studies on coumarin anticoagulant drugs
;Initiation of warfarin therapy without a loading dose.
;Circulation 1968;38:169-177

$INPUT ID 
TIME ; units="h"
WT ; units="kg"
AGE ; units="y"
SEX
AMT ; units="mg"
CMT
DV ; units="mg/L"
MDV
$DATA ../Data/warfarin_ivoral.csv IGNORE=#

$EST METHOD=COND INTER 
MAX=9990 SIG=3 NOABORT ;PRINT=1
$COV

$THETA
(0.001,0.1) ; POP_CL  units="L/h/70kg"
(0.001,8)   ; POP_V  units="L/70kg"
(0.001,2) ; POP_KA units="1/h"
(0.001,1) ; POP_TLAG units="h"
(0.001,0.5,1) ; POP_FORAL
$OMEGA BLOCK(2)
0.1 ; PPV_CL
0.01 0.1 ; PPV_V
$OMEGA
0.1 ; PPV_KA
0.1 ; PPV_TLAG
0.1 ; PPV_FORAL

$SIGMA
0.01 ; RUV_PROP
0.05 ; RUV_ADD units="mg/L"

$SUBR ADVAN2 TRANS1

$PK

   ; Covariate model
   GRPCL=THETA(1)*(WT/70)**0.75
   GRPV=THETA(2)*WT/70
   GRPKA=THETA(3)
   GRPLG=THETA(4)
   GRPFO=THETA(5)

   ; Individual parameters
   CL=GRPCL*EXP(ETA(1))
   V=GRPV*EXP(ETA(2))
   KA=GRPKA*EXP(ETA(3))
   TLAG=GRPLG*EXP(ETA(4))
   FORAL=GRPFO*EXP(ETA(5))

   ;Translate to PK Library standard names
   S2=V
   F1=FORAL
   ALAG1=TLAG
   K=CL/V

$ERROR
   CONC=F
   Y=CONC*(1+ERR(1))+ERR(2)

$TABLE ID TIME WT SEX AGE ; covariates
CL V KA TLAG F1 ; EBE estimates
CMT Y ; predictions
ONEHEADER NOPRINT FILE=warf.fit
