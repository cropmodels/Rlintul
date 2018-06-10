*-------------------------------------------------------------------------*
* Copyright 2013. Wageningen University, Plant Production Systems group,  *
* P.O. Box 430, 6700 AK Wageningen, The Netherlands.                      *
* You may not use this work except in compliance with the Licence.        *
* You may obtain a copy of the Licence at:                                *
*                                                                         *
* http://models.pps.wur.nl/content/licence-agreement                      *
*                                                                         *
* Unless required by applicable law or agreed to in writing, software     *
* distributed under the Licence is distributed on an "AS IS" basis,       *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
*-------------------------------------------------------------------------*

DEFINE_CALL GLA(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
                INPUT,INPUT,INPUT,INPUT,INPUT,                   OUTPUT)
DEFINE_CALL PENMAN(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,   OUTPUT,OUTPUT)
DEFINE_CALL EVAPTR(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
                   INPUT,INPUT,INPUT,                     OUTPUT,OUTPUT)
DEFINE_CALL DRUNIR(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
                   INPUT,INPUT,INPUT,              OUTPUT,OUTPUT,OUTPUT)

TITLE LINTUL2
*----------------------------------------------------------------------*
*  LINTUL, Light INTerception and UtiLization simulator                *
*          A simple general crop growth model, which simulates dry     *
*          matter production as the result of light interception and   *
*          utilization with a constant light use efficiency.           *
*  LINTUL2 is an extended version of LINTUL1 (the version of LINTUL    *
*          for optimal growing conditions). LINTUL2 includes a simple  *
*          water balance for studying effects of drought. The water    *
*          balance can be found in section 6 of the program, and the   *
*          effect of drought on light use efficiency in section 4.     *
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
INCON ZERO = 0.; ROOTDI = 0.1
      WLVI = LAII / SLA
      WAI  = 1000. * ROOTDI * WCI

*     Run control
FINISH TSUM > 2080.
TIMER STTIME = 58.; FINTIM = 300.; DELT = 1.; PRDEL = 5.
TRANSLATION_GENERAL DRIVER='EUDRIV'
PRINT LAI, WSOTHA, WSO, WST, WLV, WRT, TSUM, DAVTMP, DTR, ...
      ROOTD, TRAN, EVAP, TRANRF, WA, WC, YEAR , RAIN , TRAIN

***   2. Environmental data and temperature sum

DYNAMIC

WEATHER WTRDIR='E:\homeFSTWin_4.12\FSTWin_Examples'; CNTR='NLD'; ISTN=1;...
                                                       IYEAR=1971
*     Reading weather data from weather file:
*     RDD    Daily global radiation        J/(m2*d)
*     TMMN   Daily minimum temperature     degree C
*     TMMX   Daily maximum temperature     degree C
*     VP     Vapour pressure               kPa
*     WN     Wind speed                    m/s
*     RAIN   Precipitation                 mm

TRAIN = INTGRL(ZERO , RAIN)

      DTR    = RDD/1.E+6
      DAVTMP = 0.5 * (TMMN + TMMX)
      DTEFF  = MAX ( 0., DAVTMP-TBASE )
      EMERG  = MAX ( REAAND(TIME-DOYEM+1.,WC-WCWP), INSW(-LAI,1.,0.) )
      TSUM   = INTGRL(ZERO, RTSUM)
      RTSUM  = DTEFF*EMERG

***   3. Leaf growth and senescence

      CALL GLA( TIME,DOYEM,DTEFF,TSUM,LAII,RGRL,DELT,SLA,LAI,GLV,...
                TRANRF,WC,WCWP,...
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
      GTOTAL = LUE * PARINT * TRANRF

***   5. Growth rates and dry matter production of plant organs

      FRTWET = AFGEN( FRTTB, TSUM )
         FRTMOD = MAX( 1., 1./(TRANRF+0.5) )
      FRT    = FRTWET * FRTMOD
         FSHMOD = (1.-FRT) / (1.-FRT/FRTMOD)
      FLV    = AFGEN( FLVTB, TSUM ) * FSHMOD
      FST    = AFGEN( FSTTB, TSUM ) * FSHMOD
      FSO    = AFGEN( FSOTB, TSUM ) * FSHMOD

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

      RROOTD = RRDMAX * INSW( WC-WCWP, 0., 1. ) * ...
               REAAND( ROOTDM-ROOTD, TSUMAN-TSUM ) * EMERG
      ROOTD  = INTGRL( ROOTDI, RROOTD)

***   6. Soil moisture balance

      EXPLOR = 1000. * RROOTD * WCFC
      RNINTC = MIN( RAIN, 0.25*LAI )

      CALL PENMAN( DAVTMP,VP,DTR,LAI,WN,RNINTC, ...
                   PEVAP,PTRAN)
      CALL EVAPTR( PEVAP,PTRAN,ROOTD,WA,WCAD,WCWP,WCFC,WCWET,WCST,...
                   TRANCO,DELT,...
                   EVAP,TRAN)
          TRANRF = TRAN / NOTNUL(PTRAN)
      CALL DRUNIR( RAIN,RNINTC,EVAP,TRAN,IRRIGF,...
                   DRATE,DELT,WA,ROOTD,WCFC,WCST,...
                   DRAIN,RUNOFF,IRRIG)
      RWA = (RAIN+EXPLOR+IRRIG) - (RNINTC+RUNOFF+TRAN+EVAP+DRAIN)
      WA  = INTGRL( WAI,RWA)
      WC  = 0.001 * WA/ROOTD

***   7. Functions and parameters for spring wheat

*     Section 1
PARAM WCI   = 0.36
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
PARAM ROOTDM = 1.2; RRDMAX = 0.012
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

*     Section 6
PARAM WCAD = 0.08; WCWP = 0.23; WCFC = 0.36; WCWET = 0.48; WCST = 0.55
PARAM TRANCO = 8.; DRATE = 50.; IRRIGF = 1.

***********************************************************************

END

STOP

* ---------------------------------------------------------------------*
*  SUBROUTINE GLA                                                      *
*  Purpose: This subroutine computes daily increase of leaf area index *
*           (ha leaf/ ha ground/ d)                                    *
* ---------------------------------------------------------------------*

      SUBROUTINE GLA(TIME,DOYEM,DTEFF,TSUM,LAII,RGRL,DELT,SLA,LAI,GLV,
     $               TRANRF,WC,WCWP,
     $               GLAI)
      IMPLICIT REAL (A-Z)

*---- Growth during maturation stage:
      GLAI = SLA * GLV

*---- Growth during juvenile stage:
      IF ((TSUM.LT.330.).AND.(LAI.LT.0.75))
     $   GLAI = LAI * (EXP(RGRL * DTEFF * DELT) - 1.) / DELT * TRANRF

*---- Growth at day of seedling emergence:
      IF ((TIME.GE.DOYEM).AND.(LAI.EQ.0.).AND.(WC.GT.WCWP))
     $   GLAI = LAII / DELT

*---- Growth before seedling emergence:
      IF (TIME.LT.DOYEM) GLAI = 0.

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE PENMAN                                                   *
*  Purpose: Computation of the PENMAN EQUATION                         *
* ---------------------------------------------------------------------*

      SUBROUTINE PENMAN(DAVTMP,VP,DTR,LAI,WN,RNINTC,
     $                  PEVAP,PTRAN)
      IMPLICIT REAL (A-Z)

      DTRJM2 = DTR * 1.E6
      BOLTZM = 5.668E-8
      LHVAP  = 2.4E6
      PSYCH  = 0.067

      BBRAD  = BOLTZM * (DAVTMP+273.)**4 * 86400.
      SVP    = 0.611 * EXP(17.4 * DAVTMP / (DAVTMP + 239.))
      SLOPE  = 4158.6 * SVP / (DAVTMP + 239.)**2
      RLWN   = BBRAD * MAX(0.,0.55*(1.-VP/SVP))
      NRADS  = DTRJM2 * (1.-0.15) - RLWN
      NRADC  = DTRJM2 * (1.-0.25) - RLWN
      PENMRS = NRADS * SLOPE/(SLOPE+PSYCH)
      PENMRC = NRADC * SLOPE/(SLOPE+PSYCH)

      WDF    = 2.63 * (1.0 + 0.54 * WN)
      PENMD  = LHVAP * WDF * (SVP-VP) * PSYCH/(SLOPE+PSYCH)

      PEVAP  =     EXP(-0.5*LAI)  * (PENMRS + PENMD) / LHVAP
      PTRAN  = (1.-EXP(-0.5*LAI)) * (PENMRC + PENMD) / LHVAP
      PTRAN  = MAX( 0., PTRAN-0.5*RNINTC )

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE EVAPTR                                                   *
*  Purpose: To compute actual rates of evaporation and transpiration   *
* ---------------------------------------------------------------------*

      SUBROUTINE EVAPTR(PEVAP,PTRAN,ROOTD,WA,WCAD,WCWP,WCFC,WCWET,WCST,
     $                  TRANCO,DELT,
     $                  EVAP,TRAN)
      IMPLICIT REAL (A-Z)

      WC   = 0.001 * WA   / ROOTD
      WAAD = 1000. * WCAD * ROOTD
      WAFC = 1000. * WCFC * ROOTD

      EVAP  = PEVAP * LIMIT( 0., 1., (WC-WCAD)/(WCFC-WCAD) )

         WCCR = WCWP + MAX( 0.01, PTRAN/(PTRAN+TRANCO) * (WCFC-WCWP) )
         IF (WC.GT.WCCR) THEN 
             FR = LIMIT( 0., 1., (WCST-WC)/(WCST-WCWET) )
         ELSE
             FR = LIMIT( 0., 1., (WC-WCWP)/(WCCR-WCWP)  )
         ENDIF
      TRAN = PTRAN * FR

         AVAILF = MIN( 1., ((WA-WAAD)/DELT)/NOTNUL(EVAP+TRAN) )
      EVAP = EVAP * AVAILF
      TRAN = TRAN * AVAILF

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE DRUNIR                                                   *
*  Purpose: To compute rates of drainage, runoff and irrigation        *
* ---------------------------------------------------------------------*

      SUBROUTINE DRUNIR(RAIN,RNINTC,EVAP,TRAN,IRRIGF,
     $                  DRATE,DELT,WA,ROOTD,WCFC,WCST,
     $                  DRAIN,RUNOFF,IRRIG)
      IMPLICIT REAL (A-Z)

      WC   = 0.001 * WA   / ROOTD
      WAFC = 1000. * WCFC * ROOTD
      WAST = 1000. * WCST * ROOTD

      DRAIN  = LIMIT( 0., DRATE, (WA-WAFC)/DELT +
     $               (RAIN - RNINTC - EVAP - TRAN)                  )

      RUNOFF =          MAX( 0., (WA-WAST)/DELT +
     $               (RAIN - RNINTC - EVAP - TRAN - DRAIN)          )

      IRRIG  = IRRIGF * MAX( 0., (WAFC-WA)/DELT -
     $               (RAIN - RNINTC - EVAP - TRAN - DRAIN - RUNOFF) )

      RETURN
      END
