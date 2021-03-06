; Script generated by the pharmML2Nmtran Converter v.0.3.0
; Source	: PharmML 0.6.1
; Target	: NMTRAN 7.3.0
; Model 	: UseCase7
; Dated 	: Wed Dec 09 19:38:00 GMT 2015

$PROBLEM "generated by MDL2PharmML v.6.0"

$INPUT  ID TIME=TIME WT=DROP AMT CMT DVID DV MDV LOGTWT
$DATA "warfarin_conc_cmt.csv" IGNORE=@
$SUBS ADVAN13 TOL=6

$MODEL 
COMP (COMP1) 	;AD1
COMP (COMP2) 	;CENTRAL



$PK 
POP_CL = THETA(1)
POP_V = THETA(2)
POP_KA = THETA(3)
POP_TLAG = THETA(4)
RUV_PROP = THETA(5)
RUV_ADD = THETA(6)
BETA_CL_WT = THETA(7)
BETA_V_WT = THETA(8)
POP_F1 = THETA(9)

ETA_CL =  ETA(1)
ETA_V =  ETA(2)
ETA_KA =  ETA(3)
ETA_TLAG =  ETA(4)


MU_1 = LOG(POP_CL) + BETA_CL_WT * LOGTWT
CL =  EXP(MU_1 +  ETA(1)) ;

MU_2 = LOG(POP_V) + BETA_V_WT * LOGTWT
V =  EXP(MU_2 +  ETA(2)) ;

MU_3 = LOG(POP_KA)
KA =  EXP(MU_3 +  ETA(3)) ;

MU_4 = LOG(POP_TLAG)
TLAG =  EXP(MU_4 +  ETA(4)) ;

NM_F1 = POP_F1

ALAG1 = TLAG 
F1 = NM_F1 

A_0(1) = 0.0
A_0(2) = 0.0

$DES 
AD1_DES = A(1)
CENTRAL_DES = A(2)
CONC_DES = (CENTRAL_DES/V)
DADT(1) = -((KA*AD1_DES))
DADT(2) = ((KA*AD1_DES)-((CL/V)*CENTRAL_DES))

$ERROR 
AD1 = A(1)
CENTRAL = A(2)
CONC = (CENTRAL/V)
IPRED = CONC
W = RUV_ADD+RUV_PROP*IPRED
Y = IPRED+W*EPS(1)
IRES = DV - IPRED
IWRES = IRES/W

$THETA 
( 0.001 , 0.1 )	;POP_CL
( 0.001 , 8.0 )	;POP_V
( 0.001 , 0.362 )	;POP_KA
( 0.001 , 1.0 , 10.0 )	;POP_TLAG
( 0.001 , 0.1 )	;RUV_PROP
( 0.001 , 0.1 )	;RUV_ADD
(0.75  FIX )	;BETA_CL_WT
(1.0  FIX )	;BETA_V_WT
(1.0  FIX )	;POP_F1

$OMEGA BLOCK(2) CORRELATION SD
(0.1 )	;PPV_CL
(0.01 )	;ETA_CL_ETA_V
(0.1 )	;PPV_V

$OMEGA 
(0.1 SD )	;PPV_KA
(0.1  FIX SD )	;PPV_TLAG

$SIGMA 
1.0 FIX


$EST METHOD=SAEM AUTO=1 PRINT=100 CINTERVAL=30 ATOL=6 SIGL=6


$TABLE  ID TIME AMT CMT DVID MDV LOGTWT PRED IPRED RES IRES WRES IWRES Y DV NOAPPEND NOPRINT FILE=sdtab

$TABLE  ID CL V KA TLAG NM_F1 ETA_CL ETA_V ETA_KA ETA_TLAG NOAPPEND NOPRINT FILE=patab

$TABLE  ID LOGTWT NOAPPEND NOPRINT FILE=cotab


