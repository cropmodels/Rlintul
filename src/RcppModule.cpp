/*
Author Robert Hijmans
Date: May 2016
License: GNU General Public License (GNU GPL) v. 2 
*/

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#include "R_interface_util.h"
#include "LINTUL1.h"
#include "LINTUL2.h"
#include "LINTUL3.h"


void setWeather1(Lintul1Model* m, NumericVector date, NumericVector tmin, NumericVector tmax, NumericVector srad) {
	DailyWeather wth;
	wth.tmin = Rcpp::as<std::vector<double>>(tmin);
	wth.tmax = Rcpp::as<std::vector<double>>(tmax);
	wth.srad = Rcpp::as<std::vector<double>>(srad);
	wth.date = Rcpp::as<std::vector<long>>(date);
	m->wth = wth;
}


void setWeather2(Lintul2Model* m, NumericVector date, NumericVector tmin, NumericVector tmax, NumericVector srad, NumericVector prec, NumericVector wind, NumericVector vapr) {
	DailyWeather wth;
	wth.tmin = Rcpp::as<std::vector<double>>(tmin);
	wth.tmax = Rcpp::as<std::vector<double>>(tmax);
	wth.srad = Rcpp::as<std::vector<double>>(srad);
	wth.wind = Rcpp::as<std::vector<double>>(wind);
	wth.vapr = Rcpp::as<std::vector<double>>(vapr);
	wth.prec = Rcpp::as<std::vector<double>>(prec);
	wth.date = Rcpp::as<std::vector<long>>(date);
//	wth.latitude  = location[1];
//	wth.elevation = location[2];
	m->wth = wth;
}

void setWeather3(Lintul3Model* m, NumericVector date, NumericVector tmin, NumericVector tmax, NumericVector srad, NumericVector prec, NumericVector wind, NumericVector vapr) {
	DailyWeather wth;
	wth.tmin = Rcpp::as<std::vector<double>>(tmin);
	wth.tmax = Rcpp::as<std::vector<double>>(tmax);
	wth.srad = Rcpp::as<std::vector<double>>(srad);
	wth.wind = Rcpp::as<std::vector<double>>(wind);
	wth.vapr = Rcpp::as<std::vector<double>>(vapr);
	wth.prec = Rcpp::as<std::vector<double>>(prec);
	wth.date = Rcpp::as<std::vector<long>>(date);
//	wth.latitude  = location[1];
//	wth.elevation = location[2];
	m->wth = wth;
}

RCPP_EXPOSED_CLASS(DailyWeather)

RCPP_EXPOSED_CLASS(Lintul1Control)
RCPP_EXPOSED_CLASS(Lintul1Crop)
RCPP_EXPOSED_CLASS(Lintul1Model)
RCPP_EXPOSED_CLASS(Lintul1Output)

RCPP_EXPOSED_CLASS(Lintul2Control)
RCPP_EXPOSED_CLASS(Lintul2Crop)
RCPP_EXPOSED_CLASS(Lintul2Soil)
RCPP_EXPOSED_CLASS(Lintul2Model)
RCPP_EXPOSED_CLASS(Lintul2Output)

RCPP_EXPOSED_CLASS(Lintul3Crop)
RCPP_EXPOSED_CLASS(Lintul3Soil)
RCPP_EXPOSED_CLASS(Lintul3Control)
RCPP_EXPOSED_CLASS(Lintul3Model)
RCPP_EXPOSED_CLASS(Lintul3Output)


RCPP_MODULE(LINTUL3){
    using namespace Rcpp;

    class_<Lintul3Control>("Lintul3Control")
		.field("start", &Lintul3Control::start) 
		.field("emergence", &Lintul3Control::emergence) 
		.field("maxdur", &Lintul3Control::maxdur) 
		.field("IOPT", &Lintul3Control::IOPT) 
		.field("IDPL", &Lintul3Control::IDPL) 
		.field("DAYPL", &Lintul3Control::DAYPL) 
		.field("PL", &Lintul3Control::PL) 
		.field("IRRI", &Lintul3Control::IRRI) 
//		.field("RDMCR", &Lintul3Control::RDMCR) 
//		.field("DIRROLD", &Lintul3Control::DIRROLD) 
//		.field("DIRRO", &Lintul3Control::DIRRO) 
//		.field("DIRRN", &Lintul3Control::DIRRN) 
//		.field("DIRR", &Lintul3Control::DIRR) 
//		.field("DIRR1", &Lintul3Control::DIRR1) 
		.field("IRRTAB", &Lintul3Control::IRRTAB) 
		.field("FERNTAB", &Lintul3Control::FERNTAB) 
		.field("FERPTAB", &Lintul3Control::FERPTAB) 
		.field("FERKTAB", &Lintul3Control::FERKTAB) 		
	;

    class_<DailyWeather>("DailyWeather")
//		.field("longitude", &DailyWeather::longitude) 
//		.field("latitude", &DailyWeather::latitude) 
//		.field("elevation", &DailyWeather::elevation) 
		.field("CO2",  &DailyWeather::CO2) 
		.field("date", &DailyWeather::date) 
		.field("srad", &DailyWeather::srad) 
		.field("tmin", &DailyWeather::tmin) 
		.field("tmax", &DailyWeather::tmax) 
		.field("prec", &DailyWeather::prec) 
		.field("wind", &DailyWeather::wind) 
		.field("vapr", &DailyWeather::vapr) 
	;
	
    class_<Lintul3Crop>("Lintul3Crop")
		.field("IDSL", &Lintul3Crop::IDSL)
		.field("TSUM1", &Lintul3Crop::TSUM1)
		.field("TSUM2", &Lintul3Crop::TSUM2)
		.field("DVSI", &Lintul3Crop::DVSI)
		.field("DVSEND", &Lintul3Crop::DVSEND)
		.field("TDWI", &Lintul3Crop::TDWI)
		.field("RGRLAI", &Lintul3Crop::RGRLAI)
		.field("SPA", &Lintul3Crop::SPA)
		.field("DTSMTB", &Lintul3Crop::DTSMTB)
		.field("SLATB", &Lintul3Crop::SLATB)
		.field("SSATB", &Lintul3Crop::SSATB)
		.field("TBASE", &Lintul3Crop::TBASE)
		.field("KDIFTB", &Lintul3Crop::KDIFTB)
		.field("RUETB", &Lintul3Crop::RUETB)
		.field("TMPFTB", &Lintul3Crop::TMPFTB)
		.field("KDIFTB", &Lintul3Crop::KDIFTB)
		.field("RUETB", &Lintul3Crop::RUETB)
		.field("TMPFTB", &Lintul3Crop::TMPFTB)
		.field("TMNFTB", &Lintul3Crop::TMNFTB)
		.field("COTB", &Lintul3Crop::COTB)
		.field("FRTB", &Lintul3Crop::FRTB)
		.field("FLTB", &Lintul3Crop::FLTB)
		.field("FSTB", &Lintul3Crop::FSTB)
		.field("FOTB", &Lintul3Crop::FOTB)
		.field("RDRL", &Lintul3Crop::RDRL)
		.field("RDRLTB", &Lintul3Crop::RDRLTB)
		.field("RDRSHM", &Lintul3Crop::RDRSHM)
		.field("RDRNS", &Lintul3Crop::RDRNS)
		.field("RDRRTB", &Lintul3Crop::RDRRTB)
		.field("RDRSTB", &Lintul3Crop::RDRSTB)
		.field("CFET", &Lintul3Crop::CFET)
		.field("DEPNR", &Lintul3Crop::DEPNR)
		.field("IAIRDU", &Lintul3Crop::IAIRDU )
		.field("RDI", &Lintul3Crop::RDI)
		.field("RRI", &Lintul3Crop::RRI)
		.field("RDMCR", &Lintul3Crop::RDMCR)
		.field("DVSDR", &Lintul3Crop::DVSDR)
		.field("DVSDLT", &Lintul3Crop::DVSDLT)
		.field("DVSNLT", &Lintul3Crop::DVSNLT)
		.field("DVSNT", &Lintul3Crop::DVSNT)
		.field("TBASEM", &Lintul3Crop::TBASEM)
		.field("TEFFMX", &Lintul3Crop::TEFFMX)
		.field("TSUMEM", &Lintul3Crop::TSUMEM)
		.field("FNTRT", &Lintul3Crop::FNTRT)
		.field("FRNX", &Lintul3Crop::FRNX)
		.field("FRPX", &Lintul3Crop::FRPX)
		.field("FRKX", &Lintul3Crop::FRKX)
		.field("LAICR", &Lintul3Crop::LAICR)
		.field("LRNR", &Lintul3Crop::LRNR)
		.field("LSNR", &Lintul3Crop::LSNR)
		.field("LRPR", &Lintul3Crop::LRPR)
		.field("LSPR", &Lintul3Crop::LSPR)
		.field("LRKR", &Lintul3Crop::LRKR)
		.field("LSKR", &Lintul3Crop::LSKR)
		.field("NLAI", &Lintul3Crop::NLAI)
		.field("NLUE", &Lintul3Crop::NLUE)
		.field("NMAXSO", &Lintul3Crop::NMAXSO)
		.field("PMAXSO", &Lintul3Crop::PMAXSO)
		.field("KMAXSO", &Lintul3Crop::KMAXSO)
		.field("NPART", &Lintul3Crop::NPART)
		.field("NFIXF", &Lintul3Crop::NFIXF)
		.field("NSLA", &Lintul3Crop::NSLA)
		.field("RNFLV", &Lintul3Crop::RNFLV)
		.field("RNFRT", &Lintul3Crop::RNFRT)
		.field("RNFST", &Lintul3Crop::RNFST)
		.field("TCNT", &Lintul3Crop::TCNT)
		.field("NMXLV", &Lintul3Crop::NMXLV)
		.field("RPFLV", &Lintul3Crop::RPFLV)
		.field("RPFRT", &Lintul3Crop::RPFRT)
		.field("RPFST", &Lintul3Crop::RPFST)
		.field("TCPT", &Lintul3Crop::TCPT)
		.field("PMXLV", &Lintul3Crop::PMXLV)
		.field("RKFLV", &Lintul3Crop::RKFLV)
		.field("RKFRT", &Lintul3Crop::RKFRT)
		.field("RKFST", &Lintul3Crop::RKFST)
		.field("TCKT", &Lintul3Crop::TCKT)
		.field("KMXLV", &Lintul3Crop::KMXLV)
		.field("PHOTTB", &Lintul3Crop::PHOTTB)
		.field("RDI", &Lintul3Crop::RDI)
	;

    class_<Lintul3Soil>("Lintul3Soil")
		.field("SMDRY", &Lintul3Soil::SMDRY)
		.field("SMW", &Lintul3Soil::SMW)
		.field("SMFC", &Lintul3Soil::SMFC)
		.field("SM0", &Lintul3Soil::SM0)
		.field("SMI", &Lintul3Soil::SMI)
		.field("SMLOWI", &Lintul3Soil::SMLOWI)
		.field("RDMSO", &Lintul3Soil::RDMSO)
		.field("RUNFR", &Lintul3Soil::RUNFR)
		.field("CFEV", &Lintul3Soil::CFEV)
		.field("KSUB", &Lintul3Soil::KSUB)
		.field("CRAIRC", &Lintul3Soil::CRAIRC)
		.field("PMINS", &Lintul3Soil::PMINS)
		.field("NMINS", &Lintul3Soil::NMINS)
		.field("KMINS", &Lintul3Soil::KMINS)
		.field("RTPMINS", &Lintul3Soil::RTPMINS)
		.field("RTNMINS", &Lintul3Soil::RTNMINS)
		.field("RTKMINS", &Lintul3Soil::RTKMINS)
		.field("PRFTAB", &Lintul3Soil::PRFTAB)
		.field("NRFTAB", &Lintul3Soil::NRFTAB)
		.field("KRFTAB", &Lintul3Soil::KRFTAB)	
	;

	
    class_<Lintul3Model>("Lintul3Model")
		.constructor()
		.method("run", &Lintul3Model::model_run, "run the model")		
		.method("setWeather", &setWeather3)
		.field("crop", &Lintul3Model::crop, "crop")
		.field("soil", &Lintul3Model::soil, "soil")
		.field("control", &Lintul3Model::control, "control")
		.field("out", &Lintul3Model::out, "out")
		.field("weather", &Lintul3Model::wth, "weather")
		
	;			

    class_<Lintul3Output>("Lintul3Output")
		.field_readonly("step", &Lintul3Output::step)
		.field_readonly("TSUM", &Lintul3Output::TSUM)
		.field_readonly("LAI", &Lintul3Output::LAI)
		.field_readonly("WLV", &Lintul3Output::WLV)
		.field_readonly("WLVD", &Lintul3Output::WLVD)
		.field_readonly("WLVG", &Lintul3Output::WLVG)
		.field_readonly("WST", &Lintul3Output::WST)
		.field_readonly("WRT", &Lintul3Output::WRT)
		.field_readonly("WSO", &Lintul3Output::WSO)
		.field_readonly("TRANRF", &Lintul3Output::TRANRF)
	;
	
}



RCPP_MODULE(LINTUL2){
    using namespace Rcpp;

    class_<Lintul2Control>("Lintul2Control")
		.field("start", &Lintul2Control::start) 
		.field("emergence", &Lintul2Control::emergence) 
		.field("maxdur", &Lintul2Control::maxdur) 
	;

    class_<DailyWeather>("DailyWeather")
//		.field("longitude", &DailyWeather::longitude) 
//		.field("latitude", &DailyWeather::latitude) 
//		.field("elevation", &DailyWeather::elevation) 
		.field("CO2",  &DailyWeather::CO2) 
		.field("date", &DailyWeather::date) 
		.field("srad", &DailyWeather::srad) 
		.field("tmin", &DailyWeather::tmin) 
		.field("tmax", &DailyWeather::tmax) 
		.field("prec", &DailyWeather::prec) 
		.field("wind", &DailyWeather::wind) 
		.field("vapr", &DailyWeather::vapr) 
	;
	
    class_<Lintul2Crop>("Lintul2Crop")
		.field("LAIi",   &Lintul2Crop::LAIi) 
		.field("SLA",    &Lintul2Crop::SLA) 
		.field("Tbase",  &Lintul2Crop::Tbase) 
		.field("RGRL",   &Lintul2Crop::RGRL) 
		.field("Tsum1",  &Lintul2Crop::Tsum1) 
		.field("Tsum2",  &Lintul2Crop::Tsum2) 
		.field("LAIcr",  &Lintul2Crop::LAIcr) 
		.field("RDRSHM", &Lintul2Crop::RDRSHM) 
		.field("RUE",    &Lintul2Crop::RUE) 
		.field("K",      &Lintul2Crop::K) 
		.field("RDRT",   &Lintul2Crop::RDRT)   
		.field("FRTTB",  &Lintul2Crop::FRTTB)  
		.field("FLVTB",  &Lintul2Crop::FLVTB)  
		.field("FSTTB",  &Lintul2Crop::FSTTB)  
		.field("FSOTB",  &Lintul2Crop::FSOTB)  
		.field("ROOTDi", &Lintul2Crop::ROOTDi)  
		.field("ROOTDM", &Lintul2Crop::ROOTDM)  
		.field("RRDMAX", &Lintul2Crop::RRDMAX)  
		.field("TRANCO", &Lintul2Crop::TRANCO)  
	;

    class_<Lintul2Soil>("Lintul2Soil")
		.field("WCi", &Lintul2Soil::WCi)  
		.field("WCAD", &Lintul2Soil::WCAD)  
		.field("WCWP", &Lintul2Soil::WCWP)  
		.field("WCFC", &Lintul2Soil::WCFC)  
		.field("WCWET", &Lintul2Soil::WCWET)  
		.field("WCST", &Lintul2Soil::WCST)  
		.field("DRATE", &Lintul2Soil::DRATE)  
		.field("IRRIGF", &Lintul2Soil::IRRIGF)  
	;
		
    class_<Lintul2Model>("Lintul2Model")
		.constructor()
		.method("run", &Lintul2Model::model_run, "run the model")		
		.method("setWeather", &setWeather2)

		.field("crop", &Lintul2Model::crop, "crop")
		.field("soil", &Lintul2Model::soil, "soil")
		.field("control", &Lintul2Model::control, "control")
		.field("out", &Lintul2Model::out, "out")
		.field("weather", &Lintul2Model::wth, "weather")
		
	;			

    class_<Lintul2Output>("Lintul2Output")
		.field_readonly("step", &Lintul2Output::step)
		.field_readonly("TSUM", &Lintul2Output::TSUM)
		.field_readonly("LAI", &Lintul2Output::LAI)
		.field_readonly("WLV", &Lintul2Output::WLV)
		.field_readonly("WLVD", &Lintul2Output::WLVD)
		.field_readonly("WLVG", &Lintul2Output::WLVG)
		.field_readonly("WST", &Lintul2Output::WST)
		.field_readonly("WRT", &Lintul2Output::WRT)
		.field_readonly("WSO", &Lintul2Output::WSO)
		.field_readonly("EVAP", &Lintul2Output::EVAP)
		.field_readonly("TRAN", &Lintul2Output::TRAN)
		.field_readonly("TRANRF", &Lintul2Output::TRANRF)
		.field_readonly("WA", &Lintul2Output::WA)
		.field_readonly("WC", &Lintul2Output::WC)
		.field_readonly("RWA", &Lintul2Output::RWA)
	;
	
}

	

RCPP_MODULE(LINTUL1){
    using namespace Rcpp;

    class_<Lintul1Control>("Lintul1Control")
		.field("emergence", &Lintul1Control::emergence) 
		.field("maxdur", &Lintul1Control::maxdur) 
	;

    class_<DailyWeather>("DailyWeather")
		.field("CO2",  &DailyWeather::CO2) 
		.field("date", &DailyWeather::date) 
		.field("srad", &DailyWeather::srad) 
		.field("tmin", &DailyWeather::tmin) 
		.field("tmax", &DailyWeather::tmax) 
	;
	
    class_<Lintul1Crop>("Lintul1Crop")
		.field("LAIi", &Lintul1Crop::LAIi) 
		.field("SLA", &Lintul1Crop::SLA) 
		.field("Tbase", &Lintul1Crop::Tbase) 
		.field("RGRL", &Lintul1Crop::RGRL) 
		.field("Tsum1", &Lintul1Crop::Tsum1) 
		.field("Tsum2", &Lintul1Crop::Tsum2) 
		.field("LAIcr", &Lintul1Crop::LAIcr) 
		.field("RDRSHM", &Lintul1Crop::RDRSHM) 
		.field("RUE", &Lintul1Crop::RUE) 
		.field("K", &Lintul1Crop::K) 
		.field("RDRT", &Lintul1Crop::RDRT)   
		.field("FRTTB", &Lintul1Crop::FRTTB)  
		.field("FLVTB", &Lintul1Crop::FLVTB)  
		.field("FSTTB", &Lintul1Crop::FSTTB)  
		.field("FSOTB", &Lintul1Crop::FSOTB)  
	;
	
    class_<Lintul1Model>("Lintul1Model")
		.constructor()
		.method("run", &Lintul1Model::model_run, "run the model")		
		.method("setWeather", &setWeather1)
		.field("crop", &Lintul1Model::crop, "crop")
		.field("control", &Lintul1Model::control, "control")
		.field("out", &Lintul1Model::out, "out")
		.field("weather", &Lintul1Model::wth, "weather")
	;			

    class_<Lintul1Output>("Lintul1Output")
		.field_readonly("step", &Lintul1Output::step)
		.field_readonly("TSUM", &Lintul1Output::TSUM)
		.field_readonly("DLV", &Lintul1Output::DLV)
		.field_readonly("LAI", &Lintul1Output::LAI)
		.field_readonly("WLV", &Lintul1Output::WLV)
		.field_readonly("WLVD", &Lintul1Output::WLVD)
		.field_readonly("WLVG", &Lintul1Output::WLVG)
		.field_readonly("WST", &Lintul1Output::WST)
		.field_readonly("WRT", &Lintul1Output::WRT)
		.field_readonly("WSO", &Lintul1Output::WSO)
	;
	
};

	


