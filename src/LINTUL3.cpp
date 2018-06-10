/*
Author Robert Hijmans
Date: May 2016
License: GNU General Public License (GNU GPL) v. 2 

Based on LINTUL5.FOR by Joost Wolf. Date of last revision:  August 2011
This is the same as the LINTUL-3 FST model
*/

using namespace std;
#include <vector>
#include <math.h>

#include "SimUtil.h"
#include "LINTUL3.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif


/* ---------------------------------------------------------------------------*
*  SUBROUTINE PENMAN                                                        *
*  Author : Daniel van Kraalingen                                           *
*           based on an earlier version written by: Kees van Diepen         *
*           and Gert-Jan Nooy                                               *
*  Date   : 9-JAN-1987; CO2 update in 1995                                                     *
*  Purpose: This subroutine calculates atmospheric transmission,            *
*           potential evaporation from a free water surface (E0),           *
*           a bare soil surface (ES0) and the potential                     *
*           transpiration of a closed crop canopy (ET0); CO2 effect on      *
*           crop transpiration is included.
*                                                                           *
*  FORMAL PARAMETERS:  (I = input,O = output,C = control,IN = init,T = time)          *
*  name    meaning                                     units  class         *
*  ----    -------                                     -----  -----         *
*  DOY    Julian day number                             -      I,T         *
*  DAYL    Astronomical day length                       h      I           *
*  sinLD   Intermediate variable                         -      I           *
*  cosLD   Intermediate variable                         -      I           *
*  ALTI    Altitude                                      m      I           *
*  TMIN    Minimum temperature during day                C      I           *
*  TMAX    Maximum temperature during day                C      I           *
*  DTR     Daily total irradiation                     MJ m-2 d-1  I         *
*  VAP     Vapour pressure                             mbar     I           *
*  WIND    Average windspeed                           m s-1    I           *
*  CO      Atmospheric CO2 concentration               mbar     I           *
*  E0      Potential evaporation of open water         cm d-1   O           *
*  ES0     Potential evaporation of soil               cm d-1   O           *
*  ETC     Potential evapotranspiration of crop (CO2 corr.)  cm d-1   O     *
*  AVRAD   Daily total irradiation                     J m-2 d-1  O         *                                                                           *
*---------------------------------------------------------------------------*/


std::vector<double> PENMAN(int DOY, double DAYL, double sinLD, double cosLD, double ALTI, double TMIN, double TMAX, double DTR, double WIND, double VAPR, double CO2) {

/* correction factor for Penman evapo-transpiration rate as a function of atmospheric CO2 concentration
   (increase in LAI and light interception  with CO2 increase taken into account in function values). */

	std::vector<std::vector<double> > FPENMTB (5, std::vector<double>(2));
	FPENMTB[0][0] = 40;
	FPENMTB[0][1] = 1.05;
	FPENMTB[1][0] = 360;
	FPENMTB[1][1] = 1;
	FPENMTB[2][0] = 720;
	FPENMTB[2][1] = 0.95;
	FPENMTB[3][0] = 1000;
	FPENMTB[3][1] = 0.92;
	FPENMTB[4][0] = 2000;
	FPENMTB[4][1] = 0.92;

//*---- initial section
      double A = 0.20;
      double B = 0.56;

//*---- Albedo for water surface, soil surface and canopy
      double REFCFW = 0.05;
      double REFCFS = 0.15;
      double REFCFC = 0.25;

//*---- Latent heat of evaporation of water (J kg-1 = J mm-1 m-2) and
//     Stefan Boltzmann constant (J m-2 d-1 K-1)) Psychrometric  instrument constant (K-1)
      double LHVAP = 2.45E6;
      double STBC  = 4.9E-3;
      double PSYCON = 0.000662;

//*---- average daily temperature
      double TMPA = (TMIN+TMAX)/2.;
/*---- error checks on some input variable ranges
      if (TMIN.GT.TMAX) STOP 'ERR in PENMAN:TMIN > TMAX'
      if (WIND<0.)   STOP 'ERR in PENMAN:WIND < 0'
      if (DTR<0.)  STOP 'ERR in PENMAN:AVRAD < 0'
      if (VAP<0.)  STOP 'ERR in PENMAN:VAP < 0'
*/
//*---- Temperature difference (Celsius)
      double TDIF = TMAX-TMIN;

//*---- Coefficient Bu in wind function, dependent on temperature difference
      double BU = 0.54+0.35 * clamp(0, 1, (TDIF - 12) / 4);

//*---- Barometric pressure (mbar), Psychrometric constant (mbar K-1)
      double PBAR = 1013. * exp(-0.034*ALTI/(TMPA+273.));
      double GAMMA = PSYCON*PBAR;

//*---- Saturated vapour pressure according to equation of Goudriaan (1977)
      double SVAP = 6.11 * exp(17.4*TMPA/(TMPA+239.));

//*---- Measured vapour pressure should not be greater than saturated vapour pressure
      double VAP = std::min(SVAP, VAPR);

//*---- Derivative of SVAP with respect to temperature, i.e. slope of the
//*     SVAP-temperature curve (mbar K-1)
      double DELTA = 239 * 17.4 * SVAP/ pow((TMPA+239), 2);

//*---- integral of sinB (DsinB) with sinB equal to sine of solar elevation
      double AOB = clamp(-1, 1, sinLD/cosLD);
      double DsinB = 3600 * (DAYL * sinLD + 24 * cosLD * sqrt(1-AOB*AOB) / M_PI);

//*---- solar constant (SC) and daily extraterrestrial radiation (ANGOT)
      double SC = 1370.*(1.+0.033*cos(2.* M_PI * DOY/365));
      double ANGOT = std::max(0.0001, SC * DsinB) ;

//*---- atmospheric transmission (ATMTR)
      double AVRAD = std::min(0.80 * ANGOT, DTR * 1.0E6);
      double ATMTR = AVRAD/ANGOT;

//*---- The expression n/N (RELSSD) from the Penman formula is estimated
//*     from the Angstrom formula: RI = RA(A+B.n/N) -> n/N = (RI/RA-A)/B,
//*     where AVRAD = RI and DSO = RA, the Angot radiation, obtained by a CALL
//*     to ASTRO and three statements from RADIAT:
      double RELSSD = clamp (0.,1.,(ATMTR-A)/B);

//*---- Terms of the Penman formula, for water surface, soil
//*     surface and canopy
//*     Net outgoing long-wave radiation (J m-2 d-1) according to Brunt (1932)
      double RB = STBC * pow((TMPA+273.), 4) * (0.56-0.079 * sqrt(VAP)) * (0.1 + 0.9 * RELSSD) ;

//*---- Net absorbed radiation
      double RNW = AVRAD * (1.-REFCFW)-RB;
      double RNS = AVRAD * (1.-REFCFS)-RB;
      double RNC = AVRAD * (1.-REFCFC)-RB;

//*---- Evaporative demand of the atmosphere (mm d-1)
      double EA = 0.26*(SVAP-VAP)*(0.5+BU*WIND);
      double EAC = 0.26*(SVAP-VAP)*(1.0+BU*WIND);

//*---- Penman formula (1948), and conversion to cm d-1
      double E0 = 0.1*(DELTA*(RNW/LHVAP)+GAMMA*EA)/(DELTA+GAMMA);
      double ES0 = 0.1*(DELTA*(RNS/LHVAP)+GAMMA*EA)/(DELTA+GAMMA);
      double ET0 = 0.1*(DELTA*(RNC/LHVAP)+GAMMA*EAC)/(DELTA+GAMMA);

//*---- Evapo(transpi)ration is limited in order to prevent negative values
       E0 = std::max(0.,E0);
      ES0 = std::max(0.,ES0);
      ET0 = std::max(0.,ET0);

//*---- correction of potential evapo-transpiration for atmospheric CO2 concentration
      double ETC = ET0 * approx(FPENMTB, CO2);

	return std::vector<double> {E0, ES0, ETC, AVRAD};
}


/*---------------------------------------------------------------------------
  SUBROUTINE ASTRO
  Author : Daniel van Kraalingen
           documented by Rob Groot
  Date   : 15 march 1987
  Purpose: This subroutine calculates astronomical daylength (DAYL)
           and photoperiodically active daylength (DAYLP), which
           are used in the calculation of canopy assimilation and in the
           calculation of the development rate.

  FORMAL PARAMETERS:  (I = input,O = output,C = control,IN = init,T = time)
  name     meaning                                    units  class
  ----     -------                                    -----  -----
  DOY     Julian daynumber                                    I,T
  LAT      latitude                                  degrees   I
  DAYL     astronomical daylength                       h      O
  DAYLP    photoperiodically active daylength           h      O
  sinLD    intermediate variable                       =     O
  cosLD    intermediate variable                       =     O
---------------------------------------------------------------------------*/

std::vector<double> ASTRO(int DOY, double LAT) {

//*---- conversion factor from degrees to radians
      double RAD = M_PI / 180.;

//*---- declination of the sun as function of daynumber (DOY)
      double DEC = -asin(sin(23.45*RAD)*cos(2 * M_PI * (DOY+10) / 365));

//*---- sinLD, cosLD are intermediate variables
      double sinLD = sin(RAD * LAT) * sin(DEC);
      double cosLD = cos(RAD * LAT) * cos(DEC);

//*---- daylength (DAYL) and photoperiodically active daylength (DAYLP)
    double DAYL = 12.0*(1.+2. * asin(clamp(-1.,+1.,sinLD/cosLD))/M_PI);
    double DAYLP = 12.0*(1.+2. * asin(clamp(-1.,+1.,(-sin(-4 * RAD) + sinLD ) / cosLD)) / M_PI);

	return std::vector<double> {DAYL, DAYLP, sinLD, cosLD};
}



double SWEAF(double ET0, double CGNR) {
/* Chapter 20 in documentation WOFOST Version 4.1 (1988)
  The fraction of easily available soil water between field capacity and wilting point is a function of the
  potential evapotranspiration rate (for a closed canopy) in cm/day, ET0, and the crop group number, CGNR (from
  1 (=drought-sensitive) to 5 (=drought-resistent)). The function SWEAF describes this relationship given in tabular
  form by Doorenbos & Kassam (1979) and by Van Keulen & Wolf (1986; p.108, table 20).
  Authors: D.M. Jansen and C.A. van Diepen, October 1986. */

      double A = 0.76;
      double B = 1.5;
//*     curve for CGNR 5, and other curves at fixed distance below it
      double sweaf = 1./(A + B * ET0) - (5. - CGNR) * 0.10;
//*     correction for lower curves (CGNR less than 3)
      if (CGNR <  3) sweaf = sweaf + (ET0-0.6)/(CGNR * (CGNR + 3));

      sweaf = clamp (0.10, 0.95, sweaf);
	  return(sweaf);
}





void Lintul3Model::weather_step() {
	//if (time > wth.tmax.size()) { 	//stop 	}

// Average daily temperature (C)
    TMPA = 0.5 * (wth.tmin[time] + wth.tmax[time]);
	wth.CO2 = 360;

// calculate daylength DAYL,DAYLP,SINLD,COSLD
	DOY = doy_from_days(control.start + step);

	std::vector<double> astro = ASTRO(DOY, wth.latitude);
	DAYLP = astro[1];

// calculate potential soil evaporation and crop transpiration E0,ES0,ETC,AVRAD
	std::vector<double> penman = PENMAN (DOY, astro[0] ,astro[2], astro[3], wth.elevation, wth.tmin[time], wth.tmax[time], wth.srad[time], wth.wind[time], wth.vapr[time], wth.CO2);
	E0 = penman[0];
	ES0 = penman[1];
	ETC = penman[2];
//	double AVRAD = penman[4];

// Daily photosynthetically active radiation (PAR, MJ/m2)
 //  PAR = AVRAD * 0.50;

	PAR = wth.srad[time] * 0.5;

//---- Total photosynthetically active radiation (MJ/m2)
//    TPAR = TPAR + PAR;

}


void Lintul3Model::crop_initialize() {
	crop.TDW = crop.TDWI;
	crop.DVS = crop.DVSI;
	crop.DVR = 0;
	crop.WRTI =  approx (crop.FRTB, crop.DVSI) * crop.TDW;
	crop.WRT = crop.WRTI;
	crop.TAGB = crop.TDW - crop.WRT;
	crop.WLVGI = approx (crop.FLTB, crop.DVSI) * crop.TAGB;
	crop.WLVG = crop.WLVGI;
	crop.LAII = crop.WLVGI * approx(crop.SLATB, crop.DVSI);
	crop.LAI = crop.LAII;
	crop.RLAI = 0;
	crop.WLVD = 0;
	crop.WSTI =  approx (crop.FSTB, crop.DVSI) * crop.TAGB;
	crop.WST = crop.WSTI;
	crop.WSOI =  approx (crop.FOTB, crop.DVSI) * crop.TAGB;
	crop.WSO = crop.WSOI;
	crop.WSTD = 0;
	crop.WRTD = 0;
	crop.RWLVG = 0;
	crop.RWST = 0;
	crop.RWRT = 0;
	crop.RWSO = 0;
	crop.DLV = 0;
	crop.DRST = 0;
	crop.DRRT = 0;
	crop.GRT = 0;
	crop.TSULP = 0;
	crop.TSUM = 0;
	crop.TSUML = 0;
	crop.DTSULP = 0;
	crop.DTSUM = 0;
	crop.DTSUML = 0;
	crop.DVRED = 1.;
	crop.TPAR = 0;
	crop.TPARINT = 0;
	crop.PARINT = 0;
	crop.ATN = 0 ;
	crop.ATP = 0;
	crop.ATK = 0;
	crop.GTSUM = 0;
	crop.RD = crop.RDI;
	crop.NMAXLVI = approx (crop.NMXLV, crop.DVSI);
	crop.NMAXSTI = crop.LSNR * crop.NMAXLVI;
	crop.NMAXRTI = crop.LRNR * crop.NMAXLVI;
	crop.ANLVI = crop.NMAXLVI * crop.WLVGI;
	crop.ANSTI = crop.NMAXSTI * crop.WSTI;
	crop.ANRTI = crop.NMAXRTI * crop.WRTI;
	crop.ANLV = crop.ANLVI;
	crop.ANST = crop.ANSTI;
	crop.ANRT = crop.ANRTI;
	crop.ANSOI = 0;
	crop.ANSO = crop.ANSOI;
	crop.NLOSSL = 0;
	crop.NLOSSR = 0;
	crop.NLOSSS = 0;
	crop.NUPTT = 0;
	crop.NFIXTT = 0;
	crop.PMAXLVI = approx (crop.PMXLV, crop.DVSI);
	crop.PMAXSTI = crop.LSPR * crop.PMAXLVI;
	crop.PMAXRTI = crop.LRPR * crop.PMAXLVI;
	crop.APLVI = crop.PMAXLVI * crop.WLVGI;
	crop.APSTI = crop.PMAXSTI * crop.WSTI;
	crop.APRTI = crop.PMAXRTI * crop.WRTI;
	crop.APLV = crop.APLVI;
	crop.APST = crop.APSTI;
	crop.APRT = crop.APRTI;
	crop.APSOI = 0;
	crop.APSO = crop.APSOI;
	crop.PLOSSL = 0;
	crop.PLOSSR = 0;
	crop.PLOSSS = 0;
	crop.PUPTT = 0;

	crop.KMAXLVI = approx (crop.KMXLV, crop.DVSI);
	crop.KMAXSTI = crop.LSKR * crop.KMAXLVI;
	crop.KMAXRTI = crop.LRKR * crop.KMAXLVI;
	crop.AKLVI = crop.KMAXLVI * crop.WLVGI;
	crop.AKSTI = crop.KMAXSTI * crop.WSTI;
	crop.AKRTI = crop.KMAXRTI * crop.WRTI;
	crop.AKLV = crop.AKLVI;
	crop.AKST = crop.AKSTI;
	crop.AKRT = crop.AKRTI;
	crop.AKSOI = 0;
	crop.AKSO = crop.AKSOI;
	crop.KLOSSL = 0;
	crop.KLOSSR = 0;
	crop.KLOSSS = 0;
	crop.KUPTT = 0;
	crop.CTRAN = 0;
	crop.CNNI = 0;
	crop.CPNI = 0;
	crop.CKNI = 0;
	crop.CNPKI = 0;
	crop.TRANRF = 1;
	crop.NNI = 1;
	crop.PNI = 1;
	crop.KNI = 1;
	crop.NPKI = 1;
	crop.FINTT = 0;

	crop.RD = crop.RDI;


	crop.alive = true;
	crop.FLOW = false;
	crop.EMERG = false;

}


void Lintul3Model::soil_initialize() {
// initialization of water balance
	soil.TDRAIN = 0;
	soil.PERC3 = 0;
	soil.TRAIN = 0;
	soil.RAIN0= 0;
	soil.TTRANS = 0;
	soil.TESOIL = 0;
	soil.TRA = 0;
	soil.EVA = 0;
	soil.TIRR = 0;
	soil.RIRR = 0;
	soil.TRUNOF = 0;
	soil.RIRR = 0;
	soil.RUNOF = 0;

/*----- Initialization of soil moisture content (SMACT, cm3/cm3), available amounts of water
(above wilting point) in maximum effective rooted zone, both actual (WAVT)
and total amount of water in maximum rooted zone (WTOT), all in cm. */
	soil.SMACT = std::max(soil.SMW, std::min(soil.SMI, soil.SMFC));
	soil.SMACTL = std::max(soil.SMW, std::min(soil.SMLOWI, soil.SMFC));
	soil.RDM = std::min(soil.RDMSO, crop.RDMCR);
	soil.WAVTI = crop.RDI * (soil.SMACT - soil.SMW);
	soil.WAVTLI = (soil.RDM - crop.RDI) * (soil.SMACTL - soil.SMW);
	soil.WAVT = soil.WAVTI;
	soil.WAVTL = soil.WAVTLI;
	soil.WTOT = crop.RDI * soil.SMACT;
	soil.WTOTL = (soil.RDM - crop.RDI) * soil.SMACTL;
	soil.WTOTN = soil.WTOT;
	soil.WTOTLN = soil.WTOTL;
	soil.DWAT = 0;
	soil.DWATL = 0;
	soil.DWOT = 0;
	soil.DWOTL = 0;
	soil.WDR = 0;
	soil.TWDR = 0;
	soil.DSLR = 3;


	soil.NMINI = soil.NMINS;
	soil.NMIN = soil.NMINI;
	soil.PMINI = soil.PMINS;
	soil.PMIN = soil.PMINI;
	soil.KMINI = soil.KMINS;
	soil.KMIN = soil.KMINI;
	soil.NMINT = 0;
	soil.PMINT = 0;
	soil.KMINT = 0;
	soil.RNMINS = 0;
	soil.RNMINT = 0;
	soil.RPMINS = 0;
	soil.RPMINT = 0;
	soil.RKMINS = 0;
	soil.RKMINT = 0;

}


void Lintul3Model::model_initialize() {
	step = 0;
	crop_initialize();
	soil_initialize();
	
	time = control.emergence - wth.date[0];
	emergence = time + control.emergence - control.start;

	DOY = doy_from_days(control.start);

//	control.DAYPL = emergence;
//	control.DAYPL = control.planting[run];
//	control.PL = false;

	DIRR = 0;
	DIRRO = 0;
}


void Lintul3Model::soil_rates() {
/*  SUBROUTINE WATBALS
  Author: Joost Wolf
  Date  : May 2011
  Purpose: This subroutine, as adapted to the LINTUL model,  simulates
  the soil water status in the
  maximum effective rooted zone, of which the depth can be limited by
  both crop and soil characteristics. The following processes are taken
  into account:
    - soil surface evaporation
    - crop transpiration
    - precipitation
    - surface runoff
    - irrigation
    - drainage to the sub-soil

  FORMAL PARAMETERS:(I= input, O= output, C= control, IN= init., T-time)
  name     meaning                                  units       class
  ----     -------                                  -----       -----
   ISOIL   number of soil data set                   -            I
   INITI   indicates initialization of run           -            I,C
   IOPT    indicates optimal (=1), water limited (=2),
           water and N limited (=3) and
           water and N, P and K limited run (=4)     -            I,C
   IRRI    automatic irrigation (=1), actual irrigation
           from table (=2) or non irrigated(=0)      -            I,C
   TERMIN  indicates terminal section                -            I,C
   EMERG   indicates crop emergence                  -            I,C
   DOY    Julian day number                         -            I
   ES0     potential evaporation from soil surface   cm d-1       I
   ETC     potential transpiration of crop           cm d-1       I
   RAIN    precipitation                             cm d-1       I
   FINTT   fractional light interception              -           I
   DEPNR   crop group number for soil water depletion -           I
   RD      actual rooting depth                      cm           I
   RDMCR   crop specific maximum rooting depth       cm           I
   RR      root growth rate                          cm d-1       I
   RDM     crop/soil related maximal rooting depth   cm           O
   CFET    crop specific correction for transpiration -           I
   IAIRDU air ducts in roots present (=1) or not(=0)  -           I,C
   SMACT   actual soil mosture content               cm3 cm-3     O
   TTRANS  cumulative crop transpiration             cm           O
   TDRAIN  cumulative drainage to the sub-soil       cm           O
   TRAIN   cumulative precipitation                  cm           O
   TESOIL  cumulative soil evaporation               cm           O
   TRUNOF  cumulative surface runoff                 cm           O
   TIRR    cumulative irrigation                     cm           O
   TRANRF  reduction as a result of drought/wetness  -            O
   RUNFR   fraction of precipitation lost by runoff  -            O
   WTOT    total water in rooted zone                cm           O
   WTOTL   total water in lower zone                 cm           O
   WAVT    total available water in rooted zone      cm           O
   WAVTL   total available water in lower zone       cm           O

       int ISOIL,IRRI,IRRI1, DOY, IAIRDU, ILIRRT,ICROP,IOPT
       bool INITI, EMERG
	   REAL IRRTAB (30)
       COMMON /WATBOUT/ TTRANS1, TDRAIN1, TRAIN1, TESOIL1, DWAVT,DWAVTL
       COMMON /WATBOUT/ TRUNOF1, TIRR1, RUNFR1, IRRI1 */


	double RWET, RDRY;

//----- Infiltration of precipitation (RAIN) and irrigation (RIRR) (cm)
	double RIRR = DIRR;
	double PERC = (1 - soil.RUNFR) * RAIN + RIRR;

//----- Water loss by surface runoff, preliminary
	double RUNOFP = soil.RUNFR * RAIN;

//------water added to root zone by root growth (cm/d), resp. total and available water
	if (crop.EMERG) {
		soil.WDR = crop.RR * soil.SMACTL;
		soil.WDRA = crop.RR * (soil.SMACTL - soil.SMW);
	}

//----- Maximum transpiration rate (cm/d) as function of light interception (FINTT) and
//*       crop-specific correction factor of potential transpiration (CFET)

	double TRMAX = std::max(0.0001, crop.CFET * ETC * crop.FINTT);

//----- Actual transpiration (cm/d) as function soil moisture content (SMACT)

	double SWDEP = SWEAF(ETC, crop.DEPNR);
	double SMCR = (1 - SWDEP) * (soil.SMFC - soil.SMW) + soil.SMW;
	RDRY = clamp(0, 1, (soil.SMACT - soil.SMW)/(SMCR - soil.SMW));


	if (! crop.IAIRDU) {
//        reduction in transpiration in case of oxygen shortage for non-rice crops (IAIRDU not 1)
//        critical soil moisture content for aeration
		double SMAIR = soil.SM0 - soil.CRAIRC;
//        count days since start oxygen shortage (up to 4 days)
		if (soil.SMACT > SMAIR) {
			soil.DSOS = std::min ((soil.DSOS + 1), 4);
		} else {
			soil.DSOS = 0.;
		}
//        maximum reduction reached after 4 days
		double RWETMX = clamp (0., 1.,(soil.SM0 - soil.SMACT)/(soil.SM0 - SMAIR));
		RWET   = RWETMX + (1 - soil.DSOS / 4) * (1 - RWETMX);
	} else {
//       no reduction for rice crops
		RWET   = 1.;
	}

// actual crop transpiration rate reduced for drought and for oxygen shortage; not done for optimal conditions
	if (control.IOPT == 1) {
		crop.TRA = TRMAX;
	} else {
		crop.TRA = std::max(0.0, std::min(soil.WAVT, RDRY * RWET * TRMAX));
//       growth reduction function for drought or oxygen shortage
	}
	crop.TRANRF = crop.TRA / TRMAX;


//------ Maximum soil evaporation rate (cm/d) as function of light interception (FINTT)
	double EVMAX = ES0 * (1 - crop.FINTT);

//------ Actual soil evaporation rate (cm/d) in dependence of days since last rain, last rainfall and topsoil moisture
	if (PERC > 0.5) {
		soil.EVA = EVMAX;
		soil.DSLR = 1;
	} else {
		soil.DSLR = soil.DSLR + 1.;
		double EVMAXT = EVMAX * clamp(0., 1., (sqrt (soil.DSLR) - sqrt(soil.DSLR - 1)) * soil.CFEV);
		soil.EVA = std::max(0., minvalue(std::vector<double> {EVMAX, EVMAXT + PERC, 10 * (soil.SMACT-soil.SMDRY)}));
	}

// Water capacity of rooted and lower zone, at field capacity and at soil saturation (cm)
	double CAP = (soil.SMFC - soil.SMACT) * crop.RD;
	double CAPL= (soil.SMFC - soil.SMACTL) * (soil.RDM - crop.RD);
	double CAP0 = (soil.SM0 - soil.SMACT) * crop.RD;
	double CAPL0 = (soil.SM0 - soil.SMACTL) * (soil.RDM - crop.RD);

// Effective percolation incl. ET losses, potential and next dependent on soil water holding capacity (cm/d)
	double PERC1P = PERC - soil.EVA - crop.TRA;
	double PERC1 = std::min(soil.KSUB + CAP0, PERC1P);

// water loss by surface runoff, final
	soil.RUNOF = RUNOFP + std::max(0.,PERC1P -PERC1);

// Rooted zone is filled up to field capacity and the excess amount is lost by percolation to lower zone

	if (CAP <= PERC1) {
		soil.PERC2 = std::min(soil.KSUB + CAPL0, PERC1 - CAP);
	} else {
		soil.PERC2 = 0;
	}

// Lower zone is filled up to field capacity and the excess amount is lost by drainage to subsoil

	if (CAPL <= soil.PERC2) {
		soil.PERC3 = std::min(soil.KSUB, soil.PERC2 - CAPL);
	} else {
		soil.PERC3 = 0;
	}

// Effective irrigation water application on previous day
	DIRROLD = DIRR ;

// Demand for irrigation water
// Automatic irrigation (cm)
	double WAVFC = crop.RD * (soil.SMFC - soil.SMW);

	if ((control.IRRI == 1) && (soil.SMACT <= (soil.SMCR + 0.02)) && (RAIN < 1)) {
		DIRR= clamp(0., 3., 0.7 * (WAVFC - soil.WAVT));
		DIRRO = 0.;
	} else if (control.IRRI == 2) {
// Actual effective irrigation from table
		DIRR1= approx (control.IRRTAB, DOY) + DIRRO;
		DIRR= std::min(3., DIRR1);
		DIRRN = DIRR1 - DIRR;
// Irrigation left for the next day
		DIRRO = DIRRN;
	} else {
		DIRR= 0.;
		DIRRO= 0.;
	}


//---- precipitation during last day
	RAIN0 = RAIN;

//----- Change in total water and available water (DWAT) in rooted and lower zones
	soil.DWOT = soil.PERC1 - soil.PERC2 + soil.WDR;
	soil.DWOTL = soil.PERC2 - soil.PERC3 - soil.WDR;
	soil.DWAT = soil.PERC1 - soil.PERC2 + soil.WDRA;
	soil.DWATL = soil.PERC2 - soil.PERC3 - soil.WDRA;

}


//void Lintul3Model::crop_rates() {}




void Lintul3Model::crop_states() {

//----- Temperature sums (C.d) from sowing/planting (P) and from emergence with and without daylength effect
    crop.TSULP = crop.TSULP + crop.DTSULP;
    crop.TSUM = crop.TSUM + crop.DTSUM;
    crop.TSUML = crop.TSUML + crop.DTSUML;

//----- Start of flowering
    if ((!crop.FLOW) & (crop.TSUML >= crop.TSUM1)) {
        crop.IDFLOW = DOY;
        crop.FLOW = true;
    }

//---- Total intercepted radiation (MJ/m2)
    crop.TPARINT = crop.TPARINT + crop.PARINT;

//---- Dry weights of total biomass, living crop organs, and total above-ground living biomass (kg DM/ha)
    crop.GTSUM = crop.GTSUM + crop.GRT;
    crop.WLVG = crop.WLVG + crop.RWLVG;

	crop.WST = crop.WST + crop.RWST;
	crop.WRT = crop.WRT + crop.RWRT;
	crop.WSO = crop.WSO + crop.RWSO;
	crop.TAGBG = crop.WLVG + crop.WST + crop.WSO;

//---- Development stage (-)
	crop.DVS = crop.DVS + crop.DVR;

//---- Dry weights of dead crop organs and total above-ground biomass incl. dead crop organs (kg DM/ha)
    crop.WLVD = crop.WLVD + crop.DLV;
	crop.WSTD = crop.WSTD + crop.DRST;
	crop.WRTD = crop.WRTD + crop.DRRT;
	crop.TAGB = crop.TAGBG + crop.WLVD + crop.WSTD;

// Total leaf weight, both green and dead (kg DM ha-1)
//	crop.WLV = crop.WLVG + crop.WLVD;

// ---  Carbon balance check
//	double CBALAN = abs(crop.GTSUM + (crop.WRTI + crop.WLVGI + crop.WSTI + crop.WSOI) - (crop.WLVG + crop.WST + crop.WSO + crop.WRT + crop.WLVD + crop.WRTD + crop.WSTD));
//	if (CBALAN >= 1) cout << "carbon balance CBALAN not 0" << endl;



//----- Rooting depth and Leaf area index
    crop.RD = crop.RD + crop.RR;
    crop.LAI = crop.LAI + crop.RLAI;


//----- Total N/P/K uptake by crop over time (kg N/P/K ha-1) from soil and by biological fixation
    crop.NUPTT = crop.NUPTT + crop.NUPTR;
    crop.PUPTT = crop.PUPTT + crop.PUPTR;
    crop.KUPTT = crop.KUPTT + crop.KUPTR;
    crop.NFIXTT = crop.NFIXTT + crop.NFIXTR;

//-----Actual N/P/K amount in various living organs and total living N/P/K amount(kg N/P/K ha-1)
    crop.ANLV = crop.ANLV + crop.RNLV;
    crop.ANST = crop.ANST + crop.RNST;
    crop.ANRT = crop.ANRT + crop.RNRT;
    crop.ANSO = crop.ANSO + crop.RNSO;
	crop.NLIVT = crop.ANLV + crop.ANST + crop.ANRT + crop.ANSO;

	crop.APLV = crop.APLV + crop.RPLV;
    crop.APST = crop.APST + crop.RPST;
    crop.APRT = crop.APRT + crop.RPRT;
    crop.APSO = crop.APSO + crop.RPSO;
	crop.PLIVT = crop.APLV + crop.APST + crop.APRT + crop.APSO;

	crop.AKLV = crop.AKLV + crop.RKLV;
    crop.AKST = crop.AKST + crop.RKST;
    crop.AKRT = crop.AKRT + crop.RKRT;
    crop.AKSO = crop.AKSO + crop.RKSO;
	crop.KLIVT = crop.AKLV + crop.AKST + crop.AKRT + crop.AKSO;

//-----N/P/K losses from leaves, roots and stems due to senescence and total N/P/K loss (kg N/P/K ha-1)
    crop.NLOSSL = crop.NLOSSL + crop.RNLDLV;
	crop.NLOSSR = crop.NLOSSR + crop.RNLDRT;
	crop.NLOSSS = crop.NLOSSS + crop.RNLDST;
    crop.NLOSST = crop.NLOSSL + crop.NLOSSR + crop.NLOSSS;

	crop.PLOSSL = crop.PLOSSL + crop.RPLDLV;
	crop.PLOSSR = crop.PLOSSR + crop.RPLDRT;
	crop.PLOSSS = crop.PLOSSS + crop.RPLDST;
    crop.PLOSST = crop.PLOSSL + crop.PLOSSR + crop.PLOSSS;

	crop.KLOSSL = crop.KLOSSL + crop.RKLDLV;
	crop.KLOSSR = crop.KLOSSR + crop.RKLDRT;
	crop.KLOSSS = crop.KLOSSS + crop.RKLDST;
    crop.KLOSST = crop.KLOSSL + crop.KLOSSR + crop.KLOSSS;

//---- total N/P/K in living and dead roots
    crop.NROOT = crop.ANRT + crop.NLOSSR;
	crop.PROOT = crop.APRT + crop.PLOSSR;
	crop.KROOT = crop.AKRT + crop.KLOSSR;


// cumulative values for TRANRF, NNI etc. and NPKI over growth period
	if (crop.EMERG) {
		crop.CTRAN = crop.CTRAN + crop.TRANRF;
		crop.CNNI = crop.CNNI + crop.NNI;
		crop.CPNI = crop.CPNI + crop.PNI;
		crop.CKNI = crop.CKNI + crop.KNI;
		crop.CNPKI = crop.CNPKI + crop.NPKI;
	}




}


void Lintul3Model::soil_states() {
/* Integrals for drainage, rainfall, runoff, irrigation, transpiration,
   evaporation, total and available water in root and lower (between maximal
   and actual root) zones (cm), and water added to root zone by root growth */
	soil.TDRAIN = soil.TDRAIN + soil.PERC3;
    soil.TRAIN = soil.TRAIN + soil.RAIN0;
    soil.TRUNOF = soil.TRUNOF + soil.RUNOF;
    soil.TIRR = soil.TIRR + soil.RIRR;
    soil.TESOIL = soil.TESOIL + soil.EVA;
    soil.TTRANS = soil.TTRANS + soil.TRA;
    soil.WTOT = soil.WTOT + soil.DWOT;
    soil.WTOTL = soil.WTOTL + soil.DWOTL;
    soil.WAVT = soil.WAVT + soil.DWAT;
	soil.WAVTL = soil.WAVTL + soil.DWATL;
	soil.TWDR = soil.TWDR + soil.WDR;
// Actual soil moisture content in root and lower zones(cm3/cm3)
    soil.SMACT = soil.WTOT / std::max(0.0001, crop.RD);
	soil.SMACTL = soil.WTOTL / std::max(0.0001, soil.RDM - crop.RD);

//----- New amount of water in rooted and lower zones (after rate calculations)
	soil.WTOTN = soil.WTOTN + soil.DWOT;
	soil.WTOTLN = soil.WTOTLN + soil.DWOTL;


//----- Check on soil water balance: WBAL and WBALL should be zero.
//	double WBAL= abs(soil.WTOT - crop.TRA - soil.EVA + RAIN + soil.WDR + control.DIRROLD - soil.RUNOF - soil.PERC2 - soil.WTOTN);
//	if (WBAL > 1.E-4) cout << "soil moisture balance WBAL not 0" << endl;
//	double WBALL= abs(soil.WTOTL + soil.PERC2 - soil.PERC3 - soil.WDR - soil.WTOTLN);
//	if (WBALL > 1.E-4) cout <<  "soil moisture balance WBALL not 0" << endl;


//----- Soil mineral N/P/K  and Total mineral N/P/K available from both fertiliser and soil (kg N/P/K ha-1)
    soil.NMIN = soil.NMIN + soil.RNMINS;
    soil.NMINT = soil.NMINT + soil.RNMINT;
    soil.PMIN = soil.PMIN + soil.RPMINS;
    soil.PMINT = soil.PMINT + soil.RPMINT;
    soil.KMIN = soil.KMIN + soil.RKMINS;
    soil.KMINT = soil.KMINT + soil.RKMINT;


}


void Lintul3Model::output_initialize() {
	out.step.resize(0);
	out.TSUM.resize(0);
	out.DVS.resize(0);
	out.LAI.resize(0);
	out.WLVD.resize(0);
	out.WLVG.resize(0);
	out.WLV.resize(0);
	out.WST.resize(0);
	out.WRT.resize(0);
	out.WSO.resize(0);
	out.ES0.resize(0);
	out.ETC.resize(0);
	out.TRANRF.resize(0);
	out.GLAI.resize(0);
	out.NNI.resize(0);
	out.NPKI.resize(0);
	out.NMINT.resize(0);
	out.NMIN.resize(0);
	out.NUPTT.resize(0);
	out.NFIXTT.resize(0);
	out.NLIVT.resize(0);
	out.NLOSST.resize(0);
	out.PMINT.resize(0);
	out.PMIN.resize(0);
	out.PUPTT.resize(0);
	out.PLIVT.resize(0);
	out.PLOSST.resize(0);
	out.KMINT.resize(0);
	out.KMIN.resize(0);
	out.KUPTT.resize(0);
	out.KLIVT.resize(0);
	out.KLOSST.resize(0);	
}


void Lintul3Model::model_output(){
	out.step.push_back(step);
	out.TSUM.push_back(crop.TSUM);
	out.DVS.push_back(crop.DVS);
	out.LAI.push_back(crop.LAI);
	out.WLVD.push_back(crop.WLVD);
	out.WLVG.push_back(crop.WLVG);
	out.WST.push_back(crop.WST);
	out.WRT.push_back(crop.WRT);
	out.WSO.push_back(crop.WSO);
	out.ES0.push_back(ES0);
	out.ETC.push_back(ETC);
	out.TRANRF.push_back(crop.TRANRF);
	out.GLAI.push_back(crop.GLAI);
	out.NNI.push_back(crop.NNI);
	out.NPKI.push_back(crop.NPKI);
	out.NMINT.push_back(soil.NMINT);
	out.NMIN.push_back(soil.NMIN);
	out.NUPTT.push_back(crop.NUPTT);
	out.NFIXTT.push_back(crop.NFIXTT);
	out.NLIVT.push_back(crop.NLIVT);
	out.NLOSST.push_back(crop.NLOSST);

	out.PMINT.push_back(soil.PMINT);
	out.PMIN.push_back(soil.PMIN);
	out.PUPTT.push_back(crop.PUPTT);
	out.PLIVT.push_back(crop.PLIVT);
	out.PLOSST.push_back(crop.PLOSST);

	out.KMINT.push_back(soil.KMINT);
	out.KMIN.push_back(soil.KMIN);
	out.KUPTT.push_back(crop.KUPTT);
	out.KLIVT.push_back(crop.KLIVT);
	out.KLOSST.push_back(crop.KLOSST);
}


void Lintul3Model::model_run() {

	crop.alive =true;
	model_initialize();

	while ((crop.alive) & (step < control.maxdur)) {
		weather_step();
		crop_rates();

		if (control.IOPT == 3) {
			crop_ratesNPK();
			crop.PNI = 1;
			crop.KNI = 1;
		} else if (control.IOPT == 4) {
			crop_ratesNPK();
		}

		if (control.IOPT > 1) {
			soil_rates();
		}
		model_output();
		crop_states();
		if (control.IOPT > 1) {
			soil_states();
		}
		time++;
		step++;
	}
}

