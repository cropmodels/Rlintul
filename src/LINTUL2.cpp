/*
Author Robert Hijmans
Date: May 2016
License: GNU General Public License (GNU GPL) v. 2 
*/


using namespace std;
#include <vector>
#include <algorithm>
#include <cmath>

#include "SimUtil.h"
#include "LINTUL2.h"


double Penman(double tavg, double vapr, double srad, double wind, double albedo) {
  static double BOLTZM = 5.668E-8;
  static double LHVAP  = 2.4E6;
  static double PSYCH  = 0.067;
  
  double DTRJM2 = srad * 1.E6;
  double BBRAD  = BOLTZM * std::pow((tavg + 273.15), 4) * 86400;
  double SVP    = 0.611 * exp(17.4 * tavg / (tavg + 239.));
  double SLOPE  = 4158.6 * SVP / std::pow((tavg + 239.), 2);
  double RLWN   = BBRAD * std::max(0., 0.55*(1. - vapr/SVP));
  double NRAD  = DTRJM2 * (1 - albedo) - RLWN;
  double PENMR = NRAD * SLOPE/(SLOPE+PSYCH);
  double WDF    = 2.63 * (1 + 0.54 * wind);
  double PENMD  = LHVAP * WDF * (SVP-vapr) * PSYCH/(SLOPE+PSYCH);
  double ET  = (PENMR + PENMD) / LHVAP;
  return(ET);
}



void Lintul2Model::weather_step() {
	//if (time > wth.tmax.size()) { 	//stop 	}
	Tavg = (wth.tmin[time] + wth.tmax[time] ) / 2;
	Teff = std::max(0.0, Tavg - crop.Tbase);

	// compute soil evaporation and crop transpiration
	// potential
	soil.PEVAP = exp(-0.5 * crop.LAI) * Penman(Tavg, wth.vapr[time], wth.srad[time], wth.wind[time], 0.15);
	crop.PTRAN = (1 - exp(-0.5 * crop.LAI)) * Penman(Tavg, wth.vapr[time], wth.srad[time], wth.wind[time], 0.25);
	RainIntercepted = std::min( wth.prec[time], 0.25 * crop.LAI );			
	crop.PTRAN = std::max(0., crop.PTRAN-0.5 * RainIntercepted);
	
	// actual
	double FR, AVAILF;
	double WAAD = 1000 * soil.WCAD * crop.ROOTD;
	double EVAP = soil.PEVAP * clamp( 0, 1, (soil.WC - soil.WCAD)/(soil.WCFC - soil.WCAD) );
	double WCCR = soil.WCWP + std::max( 0.01, crop.PTRAN/(crop.PTRAN + crop.TRANCO) * (soil.WCFC - soil.WCWP) );
	if (soil.WC > WCCR) {
	  	FR = clamp(0, 1, (soil.WCST - soil.WC)/(soil.WCST - soil.WCWET));
	} else {
	    FR = clamp(0, 1, (soil.WC - soil.WCWP)/(WCCR - soil.WCWP));
	}
	double TRAN = crop.PTRAN * FR;
	double ET = EVAP + TRAN;
	if (ET > 0) {
		AVAILF = (soil.WA - WAAD) / ET;
	} else {
		AVAILF = soil.WA - WAAD;		
	}  
	AVAILF = min(AVAILF, 1.);
	soil.EVAP = std::max(0., EVAP * AVAILF);
	crop.TRAN = std::max(0., TRAN * AVAILF);
	
	if (crop.PTRAN > 0.01) {
		crop.TRANRF = crop.TRAN / crop.PTRAN;
	} else {
		crop.TRANRF = 1;
	}	
	
}


void Lintul2Model::crop_initialize() {
	crop.GLV = 0; // growth rate of leaves
	crop.rLAI = 0;
	crop.RWLVG = 0;
	crop.DLV = 0;
	crop.RWST = 0;
	crop.RWSO = 0;
	crop.RWRT = 0;
	crop.RROOTD = 0;
	crop.LAI = 0;
	crop.WLVD = 0;
	crop.WST = 0;
	crop.WSO = 0;
	crop.WRT = 0;
	crop.WLVG = 0;
	crop.WLV = 0;
	crop.WLVi = crop.LAIi / crop.SLA;
	crop.ROOTD = crop.ROOTDi; 
	crop.TRANRF = 1;	
	crop.alive = true;
	crop.emerged = false;
	crop.emergday = true;
	crop.TRAN = 0;
}


void Lintul2Model::soil_initialize() {
	soil.WA = 1000 * crop.ROOTDi * soil.WCi;	
	soil.WC   = 0.001 * soil.WA / crop.ROOTD;
	soil.EXPLOR = 0;
	soil.RWA = 0;
	soil.EVAP = 0;
}



void Lintul2Model::model_initialize() {
	step = 0;
	Tsum = 0;
	Teff = 0;
	crop_initialize();
	soil_initialize();
	time = control.start - wth.date[0];
	emergence = time + control.emergence - control.start;

}


void Lintul2Model::soil_rates() {

	double WAFC = 1000. * soil.WCFC * crop.ROOTD;
	double WAST = 1000. * soil.WCST * crop.ROOTD;
	soil.DRAIN  = clamp(0., soil.DRATE, (soil.WA - WAFC) + (wth.prec[time] - RainIntercepted - soil.EVAP - crop.TRAN));
	soil.RUNOFF = std::max(0., (soil.WA - WAST) + (wth.prec[time] - RainIntercepted - soil.EVAP - crop.TRAN - soil.DRAIN));
	soil.IRRIG  = soil.IRRIGF * std::max(0., (WAFC - soil.WA - wth.prec[time] - RainIntercepted - soil.EVAP - crop.TRAN - soil.DRAIN - soil.RUNOFF) );
	soil.RWA = (wth.prec[time] + soil.EXPLOR + soil.IRRIG) - (RainIntercepted + soil.RUNOFF + crop.TRAN + soil.EVAP + soil.DRAIN);
}

		
void Lintul2Model::crop_rates() {

	if (!crop.emerged) { 
		crop.emerged = (time >= emergence) & (soil.WC > soil.WCWP);
		if (!crop.emerged) { 
			return; 
		} else {
			crop.emergday = true;
		}
	}

	double RDR, RDRDV, RDRSH, dLAI, gLAI; 
	double PARint, FRT, FST, FSO, GTotal;
	double FRTWET, FRTMOD, FSHMOD, FLV;
	
	// LAI growth rate	
	if (crop.emergday) {
		gLAI = crop.LAIi;
		crop.WLVG = crop.WLVi;
		crop.emergday = false;
	} else if ((Tsum < 330.) & (crop.LAI < 0.75)) { 
		gLAI = crop.LAI * (exp(crop.RGRL * Teff) - 1) * crop.TRANRF;
	} else {
		gLAI = crop.SLA * crop.GLV;
	}

		// radiation interception
	PARint = 0.5 * wth.srad[time] * (1 - exp(-crop.K * crop.LAI));
		
		// Total crop growth rate				
	GTotal = crop.RUE * PARint * crop.TRANRF;
		
		// Development stage dependent partitioning of dry matter to organs
	FRTWET = approx(crop.FRTTB, Tsum );
	FRTMOD = std::max( 1., 1/(crop.TRANRF + 0.5) );
	FRT    = FRTWET * FRTMOD;
	FSHMOD = (1 - FRT) / (1 - FRT/FRTMOD);
	FLV = approx(crop.FLVTB, Tsum) * FSHMOD;
	FST = approx(crop.FSTTB, Tsum) * FSHMOD;
	FSO = approx(crop.FSOTB, Tsum) * FSHMOD;

			
	// LAI senescence
	// relative death rate from development (ageing) 
	RDRDV = ((Tsum - crop.Tsum1) < 0) ? 0 : approx(crop.RDRT, Tavg);
	// relative death rate from self-shading
	RDRSH = clamp(0, crop.RDRSHM, crop.RDRSHM * (crop.LAI - crop.LAIcr) / crop.LAIcr);
	// relative death rate
	RDR = std::max(RDRDV, RDRSH);
	// loss of LAI				
	dLAI = crop.LAI * RDR * crop.TRANRF;
	// rate of LAI change		
	crop.rLAI = gLAI - dLAI;
	// death rate (by weight)
	crop.DLV = crop.WLVG * RDR;
		
	crop.RWLVG = GTotal * FLV - crop.DLV;
	crop.RWST = GTotal * FST;
	crop.RWSO = GTotal * FSO;
	crop.RWRT = GTotal * FRT;
	crop.GLV = GTotal * FLV;

		
	// increase in root depth
	if ((crop.ROOTD < crop.ROOTDM) & (crop.Tsum1 >  Tsum ) & (soil.WC > soil.WCWP)) {
		crop.RROOTD = crop.RRDMAX;
	} else {
		crop.RROOTD = 0;
	}
		
	// Exploration rate of new soil water layers by root depth growth
	soil.EXPLOR = 1000. * crop.RROOTD * soil.WCFC;

	if (wth.tmax[time] <  -10) {
		crop.alive = false;	
	}
	if (Tsum > crop.Tsum2) {
		crop.alive = false;
	}			
}		



void Lintul2Model::crop_states() {
	if (crop.emerged) {
		crop.LAI  = crop.LAI + crop.rLAI;
		crop.WLVG = crop.WLVG + crop.RWLVG;
		crop.WLVD = crop.WLVD + crop.DLV;
		crop.WLV  = crop.WLVG + crop.WLVD;
		crop.WST  = crop.WST + crop.RWST;
		crop.WSO  = crop.WSO + crop.RWSO;
		crop.WRT  = crop.WRT + crop.RWRT;
		crop.ROOTD = crop.ROOTD + crop.RROOTD;
		Tsum = Tsum + Teff;
	}
}

void Lintul2Model::soil_states() {	
	soil.WA = soil.WA + soil.RWA;
	soil.WC = 0.001 * soil.WA / crop.ROOTD;
}


void Lintul2Model::output_initialize() {
	out.step.resize(0);
	out.TSUM.resize(0);
	out.LAI.resize(0);
	out.WLVD.resize(0);
	out.WLV.resize(0);
	out.WLVG.resize(0);
	out.WST.resize(0);
	out.WRT.resize(0);
	out.WSO.resize(0);
	out.EVAP.resize(0);
	out.TRAN.resize(0);
	out.TRANRF.resize(0);
	out.WA.resize(0);
	out.WC.resize(0);
	out.RWA.resize(0);
}


void Lintul2Model::model_output(){
	out.step.push_back(step);
	out.TSUM.push_back(Tsum);
	out.LAI.push_back(crop.LAI);
	out.WLVD.push_back(crop.WLVD);
	out.WLV.push_back(crop.WLV);
	out.WLVG.push_back(crop.WLVG);
	out.WST.push_back(crop.WST);
	out.WRT.push_back(crop.WRT);
	out.WSO.push_back(crop.WSO);
	out.EVAP.push_back(soil.EVAP);
	out.TRAN.push_back(crop.TRAN);
	out.TRANRF.push_back(crop.TRANRF);
	out.WA.push_back(soil.WA);
	out.WC.push_back(soil.WC);
	out.RWA.push_back(soil.RWA);
}



void Lintul2Model::model_run() {

	output_initialize();
	model_initialize(); 

	// --- START SIMULATION --- 
	while ((crop.alive) & (step < control.maxdur)) {

		weather_step();
		crop_rates();
		soil_rates();	
		model_output();
		crop_states();
		soil_states();
		time++;
		step++;
						
	} // while alive		
}	  

