#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#include <vector>
#include "SimUtil.h"
#include "LINTUL3.h"
#include "R_interface_util.h"

// [[Rcpp::export(name = ".lintul3")]]
NumericMatrix lintul3(List crop, DataFrame weather, List soil, List control) {
 
// crop parameters
	struct Lintul3Crop crp;
      
	crp.IDSL = intFromList(crop, "IDSL"); 
	crp.TSUM1 = doubleFromList(crop, "TSUM1"); 
	crp.TSUM2 = doubleFromList(crop, "TSUM2"); 
	crp.DVSI = doubleFromList(crop, "DVSI"); 
	crp.DVSEND = doubleFromList(crop, "DVSEND"); 
	crp.TDWI = doubleFromList(crop, "TDWI"); 
	crp.RGRLAI = doubleFromList(crop, "RGRLAI"); 
	crp.SPA = doubleFromList(crop, "SPA"); 
	crp.DTSMTB = TBFromList(crop, "DTSMTB");
	crp.SLATB = TBFromList(crop, "SLATB");
	crp.SSATB = TBFromList(crop, "SSATB");
	crp.TBASE = doubleFromList(crop, "TBASE"); 
	crp.KDIFTB = TBFromList(crop, "KDIFTB");
	crp.RUETB = TBFromList(crop, "RUETB");
	crp.TMPFTB = TBFromList(crop, "TMPFTB");
	crp.KDIFTB = TBFromList(crop, "KDIFTB");
	crp.RUETB = TBFromList(crop, "RUETB");
	crp.TMPFTB = TBFromList(crop, "TMPFTB");
	crp.TMNFTB = TBFromList(crop, "TMNFTB");
	crp.COTB = TBFromList(crop, "COTB");
	crp.FRTB = TBFromList(crop, "FRTB");
	crp.FLTB = TBFromList(crop, "FLTB");
	crp.FSTB = TBFromList(crop, "FSTB");
	crp.FOTB = TBFromList(crop, "FOTB");
	crp.RDRL = doubleFromList(crop, "RDRL"); 
	crp.RDRLTB = TBFromList(crop, "RDRLTB");
	crp.RDRSHM = doubleFromList(crop, "RDRSHM"); 
	crp.RDRNS = doubleFromList(crop, "RDRNS"); 
	crp.RDRRTB = TBFromList(crop, "RDRRTB");
	crp.RDRSTB = TBFromList(crop, "RDRSTB");
	crp.CFET = doubleFromList(crop, "CFET");
	crp.DEPNR = doubleFromList(crop, "DEPNR");
	crp.IAIRDU = intFromList(crop, "IAIRDU");
	crp.RDI = doubleFromList(crop, "RDI");
	crp.RRI = doubleFromList(crop, "RRI");
	crp.RDMCR = doubleFromList(crop, "RDMCR");
	crp.DVSDR = doubleFromList(crop, "DVSDR");
	crp.DVSDLT = doubleFromList(crop, "DVSDLT");
	crp.DVSNLT = doubleFromList(crop, "DVSNLT");
	crp.DVSNT = doubleFromList(crop, "DVSNT");
	crp.TBASEM = doubleFromList(crop, "TBASEM");
	crp.TEFFMX = doubleFromList(crop, "TEFFMX");
	crp.TSUMEM = doubleFromList(crop, "TSUMEM");
	crp.FNTRT = doubleFromList(crop, "FNTRT");
	crp.FRNX = doubleFromList(crop, "FRNX");
	crp.FRPX = doubleFromList(crop, "FRPX");
	crp.FRKX = doubleFromList(crop, "FRKX");
	crp.LAICR = doubleFromList(crop, "LAICR");
	crp.LRNR = doubleFromList(crop, "LRNR");
	crp.LSNR = doubleFromList(crop, "LSNR");
	crp.LRPR = doubleFromList(crop, "LRPR");
	crp.LSPR = doubleFromList(crop, "LSPR");
	crp.LRKR = doubleFromList(crop, "LRKR");
	crp.LSKR = doubleFromList(crop, "LSKR");
	crp.NLAI = doubleFromList(crop, "NLAI");
	crp.NLUE = doubleFromList(crop, "NLUE");
	crp.NMAXSO = doubleFromList(crop, "NMAXSO");
	crp.PMAXSO = doubleFromList(crop, "PMAXSO");
	crp.KMAXSO = doubleFromList(crop, "KMAXSO");
	crp.NPART = doubleFromList(crop, "NPART");
	crp.NFIXF = doubleFromList(crop, "NFIXF");
	crp.NSLA = doubleFromList(crop, "NSLA");
	crp.RNFLV = doubleFromList(crop, "RNFLV");
	crp.RNFRT = doubleFromList(crop, "RNFRT");
	crp.RNFST = doubleFromList(crop, "RNFST");
	crp.TCNT = doubleFromList(crop, "TCNT"); 
	crp.NMXLV = TBFromList(crop, "NMXLV");
	crp.RPFLV = doubleFromList(crop, "RPFLV");
	crp.RPFRT = doubleFromList(crop, "RPFRT");
	crp.RPFST = doubleFromList(crop, "RPFST");
	crp.TCPT = doubleFromList(crop, "TCPT"); 
	crp.PMXLV = TBFromList(crop, "PMXLV");
	crp.RKFLV = doubleFromList(crop, "RKFLV");
	crp.RKFRT = doubleFromList(crop, "RKFRT");
	crp.RKFST = doubleFromList(crop, "RKFST");
	crp.TCKT = doubleFromList(crop, "TCKT"); 
	crp.KMXLV = TBFromList(crop, "KMXLV");
	crp.PHOTTB = TBFromList(crop, "PHOTTB");

    crp.RDI = doubleFromList(crop, "RDI");

	DailyWeather wth;
	wth.tmin = doubleFromDF(weather, "tmin");
	wth.tmax = doubleFromDF(weather, "tmax");
	wth.srad = doubleFromDF(weather, "srad");	
	wth.prec = doubleFromDF(weather, "prec");	
	wth.vapr = doubleFromDF(weather, "vapr");
	wth.wind = doubleFromDF(weather, "wind");	
	wth.date = longFromDF(weather, "date");	
	
	struct Lintul3Soil sol;
    sol.SMDRY =  doubleFromList(soil, "SMDRY");
    sol.SMW = doubleFromList(soil, "SMW");
    sol.SMFC = doubleFromList(soil, "SMFC");
    sol.SM0 = doubleFromList(soil, "SM0");
    sol.SMI = doubleFromList(soil, "SMI");
    sol.SMLOWI = doubleFromList(soil, "SMLOWI");
    sol.RDMSO = doubleFromList(soil, "RDMSO");
    sol.RUNFR = doubleFromList(soil, "RUNFR");
    sol.CFEV = doubleFromList(soil, "CFEV");
	sol.KSUB = doubleFromList(soil, "KSUB");
    sol.CRAIRC = doubleFromList(soil, "CRAIRC");
    sol.NRFTAB = TBFromList(soil, "NRFTAB");	  
    sol.PRFTAB = TBFromList(soil, "PRFTAB");	  
    sol.KRFTAB = TBFromList(soil, "KRFTAB");	  
    sol.NMINS = doubleFromList(soil, "NMINS");
	sol.RTNMINS = doubleFromList(soil, "RTNMINS");
	sol.PMINS = doubleFromList(soil, "PMINS");
	sol.RTPMINS = doubleFromList(soil, "RTPMINS");
	sol.KMINS = doubleFromList(soil, "KMINS");
	sol.RTKMINS = doubleFromList(soil, "RTKMINS");

	struct Lintul3Control ctr;
//* actual irrigation data
    ctr.IRRTAB  = TBFromList(control, "IRRTAB");
    ctr.FERNTAB = TBFromList(control, "FERNTAB");
    ctr.FERPTAB = TBFromList(control, "FERPTAB"); 
    ctr.FERKTAB = TBFromList(control, "FERKTAB"); 
	ctr.IOPT = intFromList(control, "IOPT");
	ctr.PL = boolFromList(control, "PL");
	ctr.IRRI = intFromList(control, "IRRI");
 
	ctr.start = intFromList(control, "start"); 
	ctr.emergence = intFromList(control, "emergence"); 
//	ctr.long_output = boolFromList(control, "long_output"); 


/*
	if (emergence[s] < wdate[0]) {
		stop("emergence requested before the beginning of the weather data");
	} else if (emergence[s] > wdate[nwth-1]) {
		stop("emergence requested after the end of the weather data");
	} else if (emergence[s] < start[s]) {
		stop("emergence requested before the start of simulation");	
	}
*/

	Lintul3Model m;
	m.crop = crp;
	m.soil = sol;
	m.control = ctr;
	m.wth = wth;
	
	m.model_run();

	size_t nr = m.out.step.size();
	NumericMatrix out(nr, 32) ;
	for( size_t i=1; i<nr; i++){
		out(i,0) = m.out.step[i];
		out(i,1) = m.out.TSUM[i];
		out(i,2) = m.out.DVS[i];
		out(i,3) = m.out.LAI[i];
		out(i,4) = m.out.WLVG[i];
		out(i,5) = m.out.WLVD[i];
		out(i,6) = m.out.WLV[i];
		out(i,7) = m.out.WST[i];
		out(i,8) = m.out.WRT[i];
		out(i,9) = m.out.WSO[i];	
		out(i,10) = m.out.ES0[i];	
		out(i,11) = m.out.ETC[i];	
		out(i,12) = m.out.TRANRF[i];	
		out(i,13) = m.out.GLAI[i];	
		out(i,14) = m.out.NNI[i];	
		out(i,15) = m.out.NPKI[i];	
		out(i,16) = m.out.NMINT[i];	
		out(i,17) = m.out.NMIN[i];	
		out(i,18) = m.out.NUPTT[i];	
		out(i,19) = m.out.NFIXTT[i];	
		out(i,20) = m.out.NLIVT[i];	
		out(i,21) = m.out.NLOSST[i];	
		
		out(i,22) = m.out.PMINT[i];	
		out(i,23) = m.out.PMIN[i];	
		out(i,24) = m.out.PUPTT[i];	
		out(i,25) = m.out.PLIVT[i];	
		out(i,26) = m.out.PLOSST[i];	
		
		out(i,27) = m.out.KMINT[i];	
		out(i,28) = m.out.KMIN[i];	
		out(i,29) = m.out.KUPTT[i];	
		out(i,30) = m.out.KLIVT[i];	
		out(i,31) = m.out.KLOSST[i];	
	}
	
	CharacterVector nms = {"step", "TSUM", "DVS", "LAI", "WLVG", "WLVD", "WLV", "WST", " WRT", "WSO", " ES0", " ETC", "TRANRF", "GLAI", "NNI", "NPKI", "NMINT", "NMIN", "NUPT", "NFIX", "NLIV", "NLOSS", "PMINT", "PMIN", "PUPT", "PLIV", "PLOSS", "KMINT", "KMIN", "KUPT", "KLIV", "KLOSS"};
	colnames(out) = nms;
	return out;	
}


