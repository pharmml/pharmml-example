; Script generated by the pharmML2Nmtran Converter v.0.3.0
; Source	: PharmML 0.6.1
; Target	: NMTRAN 7.3.0
; Model 	: UseCase17_1
; Dated 	: Wed Dec 09 19:37:42 GMT 2015

$PROBLEM "generated by MDL2PharmML v.6.0"

$INPUT  ID TIME=TIME WT=DROP AGE=DROP SEX=DROP AMT SS ADDL II DVID DV MDV LOGTWT
$DATA "warfarin_conc_SSADDL.csv" IGNORE=@
$SUBS ADVAN13 TOL=6

$MODEL 
COMP (COMP1) 	;GUT
COMP (COMP2) 	;CENTRAL



$PK 
POP_CL = THETA(1)
POP_V = THETA(2)
POP_KA = THETA(3)
RUV_PROP = THETA(4)
RUV_ADD = THETA(5)
BETA_CL_WT = THETA(6)
BETA_V_WT = THETA(7)

ETA_CL =  ETA(1)
ETA_V =  ETA(2)
ETA_KA =  ETA(3)


MU_1 = LOG(POP_CL) + BETA_CL_WT * LOGTWT
CL =  EXP(MU_1 +  ETA(1)) ;

MU_2 = LOG(POP_V) + BETA_V_WT * LOGTWT
V =  EXP(MU_2 +  ETA(2)) ;

MU_3 = LOG(POP_KA)
KA =  EXP(MU_3 +  ETA(3)) ;

A_0(1) = 0
A_0(2) = 0

$DES 
GUT_DES = A(1)
CENTRAL_DES = A(2)
RATEIN_DES = (GUT_DES*KA)
CC_DES = (CENTRAL_DES/V)
DADT(1) = -(RATEIN_DES)
DADT(2) = (RATEIN_DES-((CL*CENTRAL_DES)/V))

$ERROR 
GUT = A(1)
CENTRAL = A(2)
RATEIN = (GUT*KA)
CC = (CENTRAL/V)
IPRED = CC
W = RUV_ADD+RUV_PROP*IPRED
Y = IPRED+W*EPS(1)
IRES = DV - IPRED
IWRES = IRES/W

$THETA 
( 0.001 , 0.1 )	;POP_CL
( 0.001 , 8.0 )	;POP_V
( 0.001 , 0.362 )	;POP_KA
( 0.0 , 0.1 )	;RUV_PROP
( 0.0 , 0.1 )	;RUV_ADD
(0.75  FIX )	;BETA_CL_WT
(1.0  FIX )	;BETA_V_WT

$OMEGA BLOCK(2) CORRELATION SD
(0.1 )	;PPV_CL
(0.01 )	;ETA_CL_ETA_V
(0.1 )	;PPV_V

$OMEGA 
(0.1 SD )	;PPV_KA

$SIGMA 
1.0 FIX


$EST METHOD=SAEM AUTO=1 PRINT=100 CINTERVAL=30 ATOL=6 SIGL=6


$TABLE  ID TIME AMT SS ADDL II DVID MDV LOGTWT PRED IPRED RES IRES WRES IWRES Y DV NOAPPEND NOPRINT FILE=sdtab

$TABLE  ID CL V KA ETA_CL ETA_V ETA_KA NOAPPEND NOPRINT FILE=patab

$TABLE  ID LOGTWT NOAPPEND NOPRINT FILE=cotab


