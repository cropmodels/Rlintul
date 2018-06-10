#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#include <vector>
#include "SimUtil.h"
#include "LINTUL2.h"
#include "R_interface_util.h"



// [[Rcpp::export(name = ".lintul2")]]
NumericMatrix lintul2(List crop, DataFrame weather, List soil, List control) {
 
// crop parameters
	struct Lintul2Crop crp;
	crp.LAIi = doubleFromList(crop, "LAIi"); 
	crp.SLA = doubleFromList(crop, "SLA"); 
	crp.Tbase = doubleFromList(crop, "Tbase");
	crp.RGRL  = doubleFromList(crop, "RGRL");  
	crp.Tsum1 = doubleFromList(crop, "Tsum1"); 
	crp.Tsum2 = doubleFromList(crop, "Tsum2");  
	crp.LAIcr = doubleFromList(crop, "LAIcr");  
	crp.RDRSHM = doubleFromList(crop, "RDRSHM");
	crp.RUE = doubleFromList(crop, "RUE");  
	crp.K = doubleFromList(crop, "K");
	 
	crp.RDRT = TBFromList(crop, "RDRT");
	crp.FRTTB = TBFromList(crop, "FRTTB");
	crp.FLVTB = TBFromList(crop, "FLVTB");
	crp.FSTTB = TBFromList(crop, "FSTTB");
	crp.FSOTB = TBFromList(crop, "FSOTB");
       
	crp.TRANCO = doubleFromList(crop, "TRANCO");
	crp.ROOTDi = doubleFromList(crop, "ROOTDi");
	crp.ROOTDM = doubleFromList(crop, "ROOTDM");
	crp.RRDMAX = doubleFromList(crop, "RRDMAX");
	  
	struct DailyWeather wth;
	wth.tmin = doubleFromDF(weather, "tmin");
	wth.tmax = doubleFromDF(weather, "tmax");
	wth.srad = doubleFromDF(weather, "srad");	
	wth.prec = doubleFromDF(weather, "prec");	
	wth.vapr = doubleFromDF(weather, "vapr");
	wth.wind = doubleFromDF(weather, "wind");
	wth.date = longFromDF(weather, "date");	
//	DateVector wdate = dateFromDF(weather, "SimDate");
//	wth.startdate = SimDate(wdate[0].getDay(), wdate[0].getMonth(), wdate[0].getYear());


	struct Lintul2Soil sol;
	sol.WCi = doubleFromList(soil, "WCi");
	sol.WCAD = doubleFromList(soil, "WCAD");
	sol.WCWP = doubleFromList(soil, "WCWP");
	sol.WCFC = doubleFromList(soil, "WCFC");
	sol.WCWET = doubleFromList(soil, "WCWET");
	sol.WCST = doubleFromList(soil, "WCST");
	sol.DRATE = doubleFromList(soil, "DRATE");
	sol.IRRIGF = doubleFromList(soil, "IRRIGF");
	
	
	Lintul2Control ctr;
	ctr.maxdur = intFromList(control, "maxdur");
	ctr.start = intFromList(control, "start");
	ctr.emergence = intFromList(control, "emergence");
//	ctr.long_output = boolFromList(control, "long_output"); 

	
	//for (int s=0; s < nsim; s++) {
//	int s = 0;
/*	if (emergence[s] < wdate[0]) {
		stop("emergence requested before the beginning of the weather data");
	} else if (emergence[s] > wdate[nwth-1]) {
		stop("emergence requested after the end of the weather data");
	} else if (emergence[s] < start[s]) {
		stop("emergence requested before the start of simulation");	
	}	
*/	
	Lintul2Model m;
	m.crop = crp;
	m.soil = sol;
	m.control = ctr;
	m.wth = wth;
	
	m.model_run();

	size_t nr = m.out.step.size();
	NumericMatrix out(nr, 15) ;
	for( size_t i=1; i<nr; i++){
		out(i,0) = m.out.step[i];
		out(i,1) = m.out.TSUM[i];
		out(i,2)= m.out.LAI[i];
		out(i,3) = m.out.WLVG[i];
		out(i,4) = m.out.WLVD[i];
		out(i,5) = m.out.WLV[i];
		out(i,6) = m.out.WST[i];
		out(i,7) = m.out.WRT[i];
		out(i,8) = m.out.WSO[i];	
		out(i,9) = m.out.EVAP[i];	
		out(i,10) = m.out.TRAN[i];	
		out(i,11) = m.out.TRANRF[i];	
		out(i,12) = m.out.WA[i];	
		out(i,13) = m.out.WC[i];	
		out(i,14) = m.out.RWA[i];	
	}
	colnames(out) = CharacterVector::create("step", "Tsum", "LAI", "WLVG", "WLVD", "WLV", "WST", "WRT", "WSO", "EVAP", "TRAN", "TRANRF", "WA", "WC", "RWA");		
	return out;
	
}
