; Script generated by the pharmML2Nmtran Converter v.0.3.0
; Source	: PharmML 0.6.1
; Target	: NMTRAN 7.3.0
; Model 	: UseCase14
; Dated 	: Wed Dec 09 19:37:30 GMT 2015

$PROBLEM "generated by MDL2PharmML v.6.0"

$INPUT  ID TIME=TIME TRT AMT=DROP RATE=DROP ADDL=DROP II=DROP WT=DROP DVID DV MDV REP=DROP CPX=DROP PCAX=DROP INRX=DROP ICL=DROP IV=DROP ITABS=DROP ITLAG=DROP IF1=DROP IPCA0=DROP IEMAX=DROP IC50=DROP ITEQ=DROP NEVT=DROP
$DATA "warfarin_TTE_exact.csv" IGNORE=#
$SUBS ADVAN13 TOL=9

$MODEL
	; ONE HARDCODED COMPARTMENT 
	COMP=COMP1

$PK 
POP_HBASE = THETA(1)
POP_BTATRT = THETA(2)

DUMMY = ETA(1)

BTATRT = POP_BTATRT

H_BASE = POP_HBASE

$DES 
HBASE_DES = (H_BASE/365)
HAZTRT_DES = (BTATRT*TRT)
HAZ_DES = (HBASE_DES*(1+HAZTRT_DES))
HAZARD_FUNC_DES = HAZ_DES
DADT(1) = HAZARD_FUNC_DES

$ERROR 
HBASE = (H_BASE/365)
HAZTRT = (BTATRT*TRT)
HAZ = (HBASE*(1+HAZTRT))
CUMHAZ=A(1)        ; CUMHAZ SINCE LAST EVENT
HAZARD_FUNC = HAZ

IF (DV.EQ.0) THEN
		 Y=EXP(-CUMHAZ) ; LIKELIHOOD OF CENSORED EVENT
ENDIF
IF (DV.EQ.1) THEN
		 Y=HAZARD_FUNC*EXP(-CUMHAZ)   ; LIKELIHOOD OF EVENT AT EXACT TIME
ENDIF

$THETA 
( 0.0 , 0.1 )	;POP_HBASE
( 0.0 , 0.4 )	;POP_BTATRT

$OMEGA 0 FIX

;Sim_start

$EST METHOD=SAEM AUTO=1 PRINT=100 CINTERVAL=30 ATOL=6 SIGL=6
 LIKELIHOOD LAPLACE NUMERICAL NOINTERACTION
;$SIM (12345) (12345 UNIFORM) ONLYSIM NOPREDICTION
;Sim_end

$TABLE  ID TIME TRT DVID MDV Y DV NOAPPEND NOPRINT FILE=sdtab

$TABLE  ID BTATRT H_BASE NOAPPEND NOPRINT FILE=patab

$TABLE  ID TRT NOAPPEND NOPRINT FILE=cotab


