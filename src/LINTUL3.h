//#include "SimDate.h"
#include "SimUtil.h"
#include <string>


struct Lintul3Output {
	std::vector<unsigned> step;
	std::vector<double> TSUM, DVS, LAI, WLVD, WLVG, WLV, WST, WRT, WSO, ES0, ETC, TRANRF, GLAI, 
						NNI, NPKI, NMINT, NMIN, NUPTT, NLIVT, NLOSST, NFIXTT,  
									PMINT, PMIN, PUPTT, PLIVT, PLOSST, 
									KMINT, KMIN, KUPTT, KLIVT, KLOSST;
};


struct Lintul3Control {
	long start, emergence;
	unsigned maxdur;
	int IOPT;
	int IDPL, DAYPL;
	bool PL;
	std::vector<double> IRRTAB, FERNTAB, FERPTAB, FERKTAB;
	int IRRI;
}; 


struct Lintul3Soil {
//PARAMETERS
	double SMDRY, SMW, SMFC, SM0, SMI, SMLOWI, RDMSO, RUNFR, CFEV, KSUB, CRAIRC;
	// Soil N, P and K mineralization
		// total mineral soil Npk available at start of growth period [kg/ha]
	double PMINS, NMINS, KMINS;
		//fraction of soil mineral N coming available per day [day-1]
	double RTNMINS, RTPMINS, RTKMINS;
	// Recovery fraction of applied fertiliser N, P and K
	std::vector<double> NRFTAB, KRFTAB, PRFTAB; 
	
// states	
	int DSOS;
	double TDRAIN, TRAIN, TRUNOF, TIRR, TESOIL, TTRANS, WTOT, WTOTL, WTOTN, WTOTLN, WAVT, WAVTL, TWDR, SMACT, SMACTL;
    double DSLR, WAVTLI, WAVTI;
        
// rates
	double PERC1, PERC2, PERC3, RAIN0, RUNOF, RIRR, EVA, TRA, DWOT, DWOTL, DWAT, DWATL, WDR;
	
// ?	
	double RDM, WDRA, SMCR;
	double RKMINT,RNMINT, RPMINT;
	double NMINT, PMINT, KMINT;

	double RKMINS,  RNMINS, RPMINS;
	double NMINI, NMIN, PMINI, PMIN, KMINI, KMIN; 

// params ?	

}; 



struct Lintul3Crop {
	bool alive, IAIRDU, EMERG, FLOW;
	
	double CFET, DEPNR, DVSDLT, DVSDR, DVSEND, DVSI, DVSNLT, DVSNT, FNTRT, FRKX, FRNX, FRPX,  LAICR, LRKR, LRNR, LRPR, LSKR, LSNR, LSPR;
	double NFIXF, NLAI, NLUE, NPART, NSLA, RDI, RDMCR, RDRL, RDRNS, RDRSHM, RGRLAI, RKFLV, RKFRT, RKFST, RNFLV, RNFRT, RNFST, RPFLV; 
	double RPFRT, RPFST, RRI, SPA, TBASE, TBASEM, TCKT, TCNT, TCPT, TDWI, TEFFMX, TSUM1, TSUM2, TSUMEM;
	double TRA, TRANRF;
	int RDSINT, IDSL;
	std::vector<double> COTB, DTSMTB, FLTB, FOTB, FRTB,  FSTB, KDIFTB, KMXLV, NMXLV; 
	std::vector<double> PHOTTB, PMXLV, RDRLTB, RDRRTB, RDRSTB, RUETB, SLATB, SSATB, TMNFTB, TMPFTB;

	double PMAXSO, KMAXSO, NMAXSO;

//
	double FINTT, RR, RD;
	double NPKI;
	double GLAI, SLA;
	
	double NNI, KNI, PNI;
	double NPKREF,RUE;
// RNW
   double TSULP, TSUM, TSUML;


// ?
	double TDW, DVS, DVR, WRTI, WRT, TAGB, WLVGI, WLVG, LAII, LAI, RLAI;
	double WLVD, WSTI, WST, WSOI, WSO, WSTD, WRTD, RWLVG, RWST, RWRT, RWSO, DLV, DRST, DRRT, GRT;
	double DTSULP, DTSUM, DTSUML, DVRED, TPAR, TPARINT, PARINT, ATN, ATP, ATK, GTSUM;
	double NMAXLVI, NMAXSTI, NMAXRTI, ANLVI, ANSTI, ANRTI, ANLV, ANST, ANRT, ANSOI, ANSO, NLOSSL, NLOSSR, NLOSSS, NUPTT;
	double NFIXTT, PMAXLVI, PMAXSTI, PMAXRTI, APLVI, APSTI, APRTI, APLV, APST, APRT, APSOI, APSO, PLOSSL;
	double PLOSSR, PLOSSS, PLOSST, PUPTT, KMAXLVI, KMAXSTI, KMAXRTI, AKLVI, AKSTI, AKRTI, AKLV, AKST, AKRT, AKSOI, AKSO;
	double KLOSSL, KLOSSR, KLOSSS, KLOSST, KUPTT, CTRAN, CNNI, CPNI, CKNI, CNPKI;


	double KROOT, PROOT, NROOT, RKLDST, RKLDRT, RKLDLV, RPLDST, RPLDLV, RKLV, RKST, RKRT, RKSO, KLIVT, RNLDLV, RNLDRT, RNLDST, NLOSST, RPLDRT;
	double TAGBG, NUPTR, PUPTR, KUPTR, NFIXTR, RNLV, RNST, RNRT, RNSO, NLIVT, RPLV, RPST, RPRT, RPSO, PLIVT;

	int IDFLOW;
	int DAYEM, IDEM;
//, IDEMERG;
};



struct Lintul3Model {


	Lintul3Crop crop;
	Lintul3Soil soil;
	Lintul3Control control;
	DailyWeather wth;
	Lintul3Output out;
	
//	std::vector<std::vector<double> > out;
//	std::vector<std::string> out_names;
	
//	double Teff, Tsum, RainIntercepted;
	double RAIN, RAIN0, E0, ES0, ETC, TMPA, PAR, DAYLP; 
	
	int time;
	unsigned step, DOY;
	long emergence; //, today;
	double DIRROLD, DIRRO, DIRRN, DIRR, DIRR1;  

//	Lintul3Model(Lintul3Crop c, Lintul3Soil s, Lintul3Control t, DailyWeather w) : crop(c), soil(s), control(t), wth(w) { };
	
	void weather_step();
	void output_initialize();
	
	void crop_initialize();
	void crop_rates();
	void crop_ratesNPK();
	void crop_states();

	void soil_initialize();
	void soil_rates();
	void soil_states();
	
	void model_initialize();	
	void model_run();
	void model_output();
};

