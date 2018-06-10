#include "SimUtil.h"

struct Lintul1Output {
	std::vector<unsigned> step;
	std::vector<double> TSUM, DLV, LAI, WLVD, WLV, WLVG, WST, WRT, WSO; 
};


struct Lintul1Crop {
// LINTUL1
	Lintul1Crop(){
		LAIi = 0.012;
		SLA = 0.022;
		Tbase = 0;
		RGRL  = 0.009;
		Tsum1 = 1110;
		Tsum2 = 2180;
		LAIcr = 4;
		RDRSHM = 0.03;
		RUE = 3.0;
		K = 0.6;
	}

// PARAMETERS
	double LAIi, LAIcr, SLA, Tbase, RGRL, Tsum1, Tsum2, RDRSHM, RUE, K;
//  Partitioning tables for leaves (LV), stems (ST), storage organs (SO) and roots (RT):
	std::vector<double> RDRT, FLVTB, FSTTB, FSOTB, FRTTB;

// VARIABLES
    bool emerged, emergday, alive;
	double WLVi, GLV;
	
// RATES
	double rLAI, RWLVG, DLV, RWST, RWSO, RWRT;
	
// STATES
	double TSUM, LAI, WLVD, WST, WSO, WRT, WLVG, WLV;
};

/*
struct Lintul1Weather {
//	double longitude, latitude, elevation;
	double CO2 = 400; 
	std::vector<long> date;
	std::vector<double> srad, tmin, tmax;
};
*/

struct Lintul1Control {
  long start, emergence;
  unsigned maxdur = 365;
}; 

struct Lintul1Model {
	//Lintul1Model(Lintul1Crop c, LintulControl t, Lintul1Weather w) : crop(c), control(t), wth(w) { };

	DailyWeather wth;

	Lintul1Control control;
	Lintul1Crop crop;
	Lintul1Output out;

	double Tavg, Teff;
	int time;
	unsigned step;

	void weather_step();
	void output_initialize();
	
	void crop_initialize();
	void crop_rates();
	void crop_states();

	void model_initialize();
	void model_run();
	void model_output();
};

