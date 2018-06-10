DEFINE_CALL GLA(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
                INPUT,INPUT,                                     OUTPUT)

TITLE LINTUL1
*----------------------------------------------------------------------*
*  LINTUL, Light INTerception and UtiLization simulator                *
*          A simple general crop growth model, which simulates dry     *
*          matter production as the result of light interception and   *
*          utilization with a constant light use efficiency.           *
*  LINTUL1 is the version of LINTUL for optimal growing conditions.    *
*                                                                      *
*          Example for spring wheat                                    *
*                                                                      *
*  DLO-Research Institute for Agrobiology and Soil Fertility (AB-DLO)  *
*  Dept of Theor. Prod. Ecology, Wageningen Agric. Univ. (TPE-WAU)     *
*                                                                      *
*  Reference: Spitters, C.J.T. & A.H.C.M. Schapendonk, 1990.           *
*  Evaluation of breeding strategies for drought tolerance in potato   *
*  by means of crop growth simulation. Plant and Soil 123: 193-203.    *
*----------------------------------------------------------------------*

***   1. Initial conditions and run control

INITIAL

*     Initial conditions
INCON ZERO = 0.
      WLVI = LAII / SLA

*     Run control
FINISH TSUM > 2080.
TIMER STTIME = 58.; FINTIM = 300.; DELT = 1.; PRDEL = 1.
TRANSLATION_GENERAL DRIVER='EUDRIV'
PRINT LAI, WSOTHA, WSO, WST, WLV, WRT, TSUM, DAVTMP, DTR ,TRDD, TRAIN,...
FRT,FLV,FST,FSO,TSUM

***   2. Environmental data and temperature sum

DYNAMIC

WEATHER WTRDIR=' '; CNTR='NLD'; ISTN=1; IYEAR=1971
*     Reading weather data from weather file:
*     RDD    Daily global radiation        J/(m2*d)
*     TMMN   Daily minimum temperature     degree C
*     TMMX   Daily maximum temperature     degree C

      DTR    = RDD/1.E+6
      DAVTMP = 0.5 * (TMMN + TMMX)
      DTEFF  = MAX ( 0., DAVTMP-TBASE )
      EMERG  = INSW(TIME-DOYEM, 0., 1.)
      TSUM   = INTGRL(ZERO, RTSUM)
      RTSUM  = DTEFF*EMERG

*      Extra for exercises 21 January 2008
TRDD = INTGRL(ZERO, RDD)
TRAIN= INTGRL(ZERO, RAIN)
***   3. Leaf growth and senescence

      CALL GLA(TIME,DOYEM,DTEFF,TSUM,LAII,RGRL,DELT,SLA,LAI,GLV,...
               GLAI)
      GLV   = FLV * GTOTAL

      DLAI  = LAI * RDR 
      RDR   = MAX(RDRDV, RDRSH)
      RDRDV = INSW(TSUM-TSUMAN, 0., AFGEN(RDRT, DAVTMP))
      RDRSH = LIMIT(0., RDRSHM, RDRSHM * (LAI-LAICR) / LAICR)
      DLV   = WLVG * RDR

      RLAI  = GLAI - DLAI
      LAI   = INTGRL(ZERO, RLAI)

***   4. Light interception and total crop growth rate

      PARINT = 0.5 * DTR    * (1. - EXP(-K*LAI))
      GTOTAL = LUE * PARINT
                 
***   5. Growth rates and dry matter production of plant organs

      FRT    = AFGEN( FRTTB, TSUM )
      FLV    = AFGEN( FLVTB, TSUM )
      FST    = AFGEN( FSTTB, TSUM )
      FSO    = AFGEN( FSOTB, TSUM )

      WLVG   = INTGRL( WLVI, RWLVG)
      WLVD   = INTGRL( ZERO, DLV  )
      WST    = INTGRL( ZERO, RWST )
      WSO    = INTGRL( ZERO, RWSO )
         WSOTHA = WSO / 100.
      WRT    = INTGRL( ZERO, RWRT )
      WLV    = WLVG + WLVD
      RWLVG  = GTOTAL * FLV - DLV
      RWST   = GTOTAL * FST
      RWSO   = GTOTAL * FSO
      RWRT   = GTOTAL * FRT

***   6. Functions and parameters for spring wheat

*     Section 1
PARAM LAII  = 0.012; SLA = 0.022

*     Section 2
PARAM TBASE = 0.

*     Section 3
PARAM DOYEM = 60.
PARAM RGRL  = 0.009; TSUMAN = 1110.; LAICR = 4.; RDRSHM = 0.03
FUNCTION RDRT = -10.,0.03, 10.,0.03, 15.,0.04, 30.,0.09, 50.,0.09

*     Section 4
PARAM LUE = 3.0; K = 0.6

*     Section 5
*     Partitioning tables for leaves (LV), stems (ST), 
*     storage organs (SO) and roots (RT):
FUNCTION FRTTB =     0.,0.50,   110.,0.50,   275.,0.34,   555.,0.12, ...
      780.,0.07,  1055.,0.03,  1160.,0.02,  1305.,0.  ,  2500.,0. 
FUNCTION FLVTB =     0.,0.33,   110.,0.33,   275.,0.46,   555.,0.44, ...
      780.,0.14,  1055.,0.  ,                            2500.,0.
FUNCTION FSTTB =     0.,0.17,   110.,0.17,   275.,0.20,   555.,0.44, ...
      780.,0.79,  1055.,0.97,  1160.,0.  ,               2500.,0.
FUNCTION FSOTB =     0.,0.  ,                                        ...
                  1055.,0.  ,  1160.,0.98,  1305.,1.  ,  2500.,1.

************************************************************************
END
STOP

* ---------------------------------------------------------------------*
*  SUBROUTINE GLA                                                      *
*  Purpose: This subroutine computes daily increase of leaf area index *
*           (ha leaf/ ha ground/ d)                                    *
* ---------------------------------------------------------------------*

      SUBROUTINE GLA(TIME,DOYEM,DTEFF,TSUM,LAII,RGRL,DELT,SLA,LAI,GLV,
     $               GLAI)
      IMPLICIT REAL (A-Z)

*---- Growth during maturation stage:
      GLAI = SLA * GLV

*---- Growth during juvenile stage:
      IF ((TSUM.LT.330.).AND.(LAI.LT.0.75))
     $   GLAI = LAI * (EXP(RGRL * DTEFF * DELT) - 1.) / DELT

*---- Growth at day of seedling emergence:
      IF ((TIME.GE.DOYEM).AND.(LAI.EQ.0.))
     $   GLAI = LAII / DELT

*---- Growth before seedling emergence:
      IF (TIME.LT.DOYEM) GLAI = 0.

      RETURN
      END
