/*
Author Robert Hijmans
Date: May 2016
License: GNU General Public License (GNU GPL) v. 2 
*/

#include "SimUtil.h"

struct Lintul2Output {
	std::vector<unsigned> step;
	std::vector<double> TSUM, LAI, WLVD, WLV, WLVG, WST, WRT, WSO, EVAP, TRAN, TRANRF, WA, WC, RWA; 
};


struct Lintul2Soil {
	Lintul2Soil(){ EXPLOR=0; PEVAP=0; EVAP=0; RUNOFF=0; DRAIN=0; IRRIG=0; RWA=0; WA=0; WC=0; };
// parameters
	double  WCi, WCAD, WCWP, WCFC, WCWET, WCST, DRATE, IRRIGF;
// variables 
	double EXPLOR, PEVAP, EVAP, RUNOFF, DRAIN, IRRIG;
// rates 
	double RWA;
// states
	double WA, WC;
}; 


struct Lintul2Crop {
// LINTUL1
// parameters
	double WLVi, LAIi, LAIcr, SLA, Tbase, RGRL, Tsum1, Tsum2, RDRSHM, RUE, K;
	//  Partitioning tables for leaves (LV), stems (ST), storage organs (SO) and roots (RT):
	std::vector<double> RDRT, FLVTB, FSTTB, FSOTB, FRTTB;
// RATES
	double rLAI, RWLVG, DLV, RWST, RWSO, RWRT;
// STATES  
	double LAI, WLVD, WST, WSO, WRT, WLVG, WLV;

// VARIABLES  
    bool emerged, emergday, alive;
	double GLV;
	double ROOTDi, TRANCO, ROOTDM, RRDMAX;
	double PTRAN, TRAN, TRANRF;
// rates?
	double RROOTD; 
// states?   
	double ROOTD;
};


struct Lintul2Control {
  long start, emergence;
  unsigned maxdur = 365;
}; 


struct Lintul2Model {

	DailyWeather wth;

	Lintul2Control control;
	Lintul2Crop crop;
	Lintul2Soil soil;
	Lintul2Output out;

	double Tavg, Teff, Tsum, RainIntercepted;
	int time, emergence; 
	unsigned step;

//	Lintul2Model(Lintul2Crop c, Lintul2Soil s, LintulControl t, DailyWeather w) : crop(c), soil(s), control(t), wth(w) { };
	
	void weather_step();
	void output_initialize();
	void crop_initialize();
	void crop_rates();
	void crop_states();
	void soil_initialize();
	void soil_rates();
	void soil_states();
	void model_initialize();	
	void model_output();
	void model_run();
	
};

