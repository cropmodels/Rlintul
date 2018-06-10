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
#include "LINTUL1.h"

void Lintul1Model::output_initialize() {
	out.step.resize(0);
	out.TSUM.resize(0);
	out.DLV.resize(0);
	out.LAI.resize(0);
	out.WLVD.resize(0);
	out.WLV.resize(0);
	out.WLVG.resize(0);
	out.WST.resize(0);
	out.WRT.resize(0);
	out.WSO.resize(0);
}


void Lintul1Model::model_output(){
	out.step.push_back(step);
	out.TSUM.push_back(crop.TSUM);
	out.DLV.push_back(crop.DLV);
	out.LAI.push_back(crop.LAI);
	out.WLVD.push_back(crop.WLVD);
	out.WLV.push_back(crop.WLV);
	out.WLVG.push_back(crop.WLVG);
	out.WST.push_back(crop.WST);
	out.WRT.push_back(crop.WRT);
	out.WSO.push_back(crop.WSO);
}


void Lintul1Model::weather_step() {
	Tavg = (wth.tmin[time] + wth.tmax[time] ) / 2;
	Teff = std::max(0.0, Tavg - crop.Tbase);
}


void Lintul1Model::crop_initialize() {
	crop.TSUM = 0;
	crop.LAI = 0;
	crop.WLVD = 0;
	crop.WST = 0;
	crop.WSO = 0;
	crop.WLV = 0;
	crop.WRT = 0;
	crop.WLVi = crop.LAIi / crop.SLA;
	crop.WLVG = crop.WLVi;
	crop.alive = true;
}



void Lintul1Model::model_initialize() {
	step = 0;
	crop_initialize();
	output_initialize();
	crop.emergday = true;	
	time = control.emergence - wth.date[0];
	// need to check for out of bounds times (before of after start)
	// stop if time < 0 || time > wth.date.size();
}



void Lintul1Model::crop_rates() {
	  // rates  
	double RDR, RDRDV, RDRSH; // leaf relative death rate / from ageing / from self-shading
	double dLAI, gLAI; // LAI rates
	double PARint; // intercepted radiation
	double FRT, FLV, FST, FSO; // partitioning
	double GTotal;
  
    // radiation interception
	PARint = 0.5 * wth.srad[time] * (1 - exp(-crop.K * crop.LAI));
    
    // Total crop growth rate				
	GTotal = crop.RUE * PARint;
    
    // Development stage dependent partitioning of dry matter to organs
	FRT = approx(crop.FRTTB, crop.TSUM);
	FLV = approx(crop.FLVTB, crop.TSUM);
	FST = approx(crop.FSTTB, crop.TSUM);
	FSO = approx(crop.FSOTB, crop.TSUM);
    
    
    // LAI growth rate	
	if (crop.emergday) {
		gLAI = crop.LAIi;
		crop.emergday = false;
	} else if ((crop.TSUM < 330) & (crop.LAI < 0.75)) { 
		gLAI = crop.LAI * (exp(crop.RGRL * Teff) - 1);
	} else {
		gLAI = crop.SLA * crop.GLV;
	}
    
    // LAI senescence
    // relative death rate from development (ageing) 
	RDRDV = (crop.TSUM - crop.Tsum1) < 0 ? 0 : approx(crop.RDRT, Tavg);
    // relative death rate from self-shading
	RDRSH = clamp(0, crop.RDRSHM, crop.RDRSHM * (crop.LAI - crop.LAIcr) / crop.LAIcr);
    // relative death rate
	RDR = std::max(RDRDV, RDRSH);
    // loss of LAI				
	dLAI = crop.LAI * RDR;
    // rate of LAI change		
	crop.rLAI = gLAI - dLAI;
    // death rate (by weight)
	crop.DLV = crop.WLVG * RDR;
    
	crop.RWLVG = GTotal * FLV - crop.DLV;
	crop.RWST = GTotal * FST;
	crop.RWSO = GTotal * FSO;
	crop.RWRT = GTotal * FRT;
	crop.GLV = GTotal * FLV;
	
	if (crop.TSUM > crop.Tsum2) {
		crop.alive = false;
	} else if (wth.tmax[time] <  -10) {
		crop.alive = false;	
	}
	
}

	

void Lintul1Model::crop_states() {
  
	crop.LAI  = crop.LAI + crop.rLAI;
	crop.WLVG = crop.WLVG + crop.RWLVG;
	crop.WLVD = crop.WLVD + crop.DLV;
	crop.WLV  = crop.WLVG + crop.WLVD;
	crop.WST  = crop.WST + crop.RWST;
	crop.WSO  = crop.WSO + crop.RWSO;
	crop.WRT  = crop.WRT + crop.RWRT;

	crop.TSUM = crop.TSUM + Teff;  
}

	
	
void Lintul1Model::model_run() {
  
	model_initialize(); 
				
	while ((crop.alive) & (step < control.maxdur)) {  
		weather_step();
		crop_rates();
		model_output();
		crop_states();
		time++;
		step++;		
	}
}

