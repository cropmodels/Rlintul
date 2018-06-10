using namespace std;
#include <vector>
#include <algorithm>
#include <iostream>
#include "SimUtil.h"
#include "LINTUL3.h"
#include <cmath>

/*  Based on FORTRAN code by Joost Wolf                                                      
  Date: Crop model developed on the basis of LINTUL3.fst in August 2011    
                                                                          
  FORMAL PARAMETERS:(I = input, O = output, C = control, IN = init., T-time)  
  name     meaning                                  units       class     
  ----     -------                                  -----       -----     
  INITI    indicates initialization of run           -           I,C      
  IOPT     indicates optimal ( = 1), water limited ( = 2),                    
           water and N limited ( = 3) and				    				 
           water and N, P and K limited run ( = 4)     -           I,C      
  DOY      julian day number                         -           I        
  IDEM     date of emergence                         -           I        
  IDPL     date of planting/sowing                   -           I        
  PL       indicates planting as start of simulation -           I,C      
  TERMIN   indicates terminal section                -           I,C      
  EMERG    indicates crop emergence                  -           I,C      
  TMIN     minimum air temperature                   C           I        
  TMAX     maximum air temperature                   C           I        
  AVRAD    daily total irradiation                   J m-2 d-1   I        
  CO       atmospheric CO2 concentration             ppmv        I        
  TRANRF   reduction factor due to drought/wetness   -           I          
  RDMSO    soil related maxiumum rooting depth       cm          I        
  DAYLP    photoperiodically active daylength        h           I        


  IDEMERG  date of emergence                         -           O        
  IDFLOW   date of flowering                                     O        
  IDHALT   end date of crop growth                               O        

  TAGB     total above-ground biomass                kg DM ha-1  O        
  WLVG     weight of living leaves                   kg DM ha-1  O        
  WLVD     weight of dead leaves                     kg DM ha-1  O        
  WST      weight of stems                           kg DM ha-1  O        
  WRT      weight of roots                           kg DM ha-1  O        
  WSO      weight of storage organs                  kg DM ha-1  O        
  RD       actual rooting depth                      cm          O        
  RDMCR    crop specific maximum rooting depth       cm          O        
  RR       root growth rate                          cm d-1      O        
  RDM      soil/crop related maximal rooting depth   cm          I        
  LAI      leaf area index                           m2 m-2      O        
  DEPNR   crop group number for soil water depletion -           O        
  CFET    crop specific correction for transpiration -           O         
  IAIRDU  air ducts in roots present ( = 1) or not( = 0)  -      O,C      
  TSUM     temperature sum from emergence            C d         O        
  DVS      development stage                         -           O        
  DVSEND   development stage at end of growth period -           O        
  TSULP    temperature sum from sowing/planting      C d         O        
  FINT     fractional light interception for PAR     -           O        
  FINTT    fractional light interception for total radiation -   O        
  TPARINT  total intercepted radiation (PAR)         MJ m-2      O        
  TPAR     total photosynthetically active radiation MJ m-2      O        
  TSUML    temperature sum from emergence incl. dayl.effect C.d  O        
  NNI      nitrogen nutrition index                  -           O        
  NPKI     NPK nutrition index ( = minimum of N/P/K-index) -     O        
  NMINT    total mineral N from soil and fertiliser  kg N ha-1   O        
  NMIN     mineral N available from soil for crop    kg N ha-1   O        
  NUPTT    total N uptake by crop from soil          kg N ha-1   O        
  NFIXTT   total N uptake by crop from biol.fixation kg N ha-1   O        
  NLIVT    Amount of N in living crop organs         kg N ha-1   O        
  NLOSST   Amount of N in dead crop organs           kg N ha-1   O        
  PMINT    total mineral P from soil and fertiliser  kg P ha-1   O        
  PMIN     mineral P available from soil for crop    kg P ha-1   O        
  PUPTT    total P uptake by crop from soil          kg P ha-1   O        
  PLIVT    Amount of P in living crop organs         kg P ha-1   O        
  PLOSST   Amount of P in dead crop organs           kg P ha-1   O        
  KMINT    total mineral K from soil and fertiliser  kg K ha-1   O        
  KMIN     mineral K available from soil for crop    kg K ha-1   O        
  KUPTT    total K uptake by crop from soil          kg K ha-1   O        
  KLIVT    Amount of K in living crop organs         kg K ha-1   O        
  KLOSST   Amount of K in dead crop organs           kg K ha-1   O        
  YCH      indicates year change                     -           I         
                                                                          
------------------------------------------------------------------------
*/ 



void Lintul3Model::crop_rates() {


	if (!crop.EMERG) {
		if (control.PL) {
// start at planting, compute emergence date
			if (crop.TSULP >= crop.TSUMEM) {
				crop.DAYEM = DOY;
				crop.EMERG = true;
				crop.DTSULP = 0;
			} else {
				crop.DTSULP = clamp(0, crop.TEFFMX - crop.TBASEM, TMPA - crop.TBASEM);
			}
		} else if (time >= emergence) {
// start at given emergence date
			crop.DAYEM = DOY;
			crop.EMERG = true;
		}
	}

	if (crop.EMERG) {
	
//// Development rate /////
		crop.DTSUM = std::max(0., approx(crop.DTSMTB, TMPA));
		double DTEFF = max(0., TMPA - crop.TBASE);


// Reduction of development rate until flowering by  day length
		double RDAYL;
		if ((crop.IDSL == 1) & (! crop.FLOW)) {
			RDAYL = approx(crop.PHOTTB, DAYLP);
		} else {
			RDAYL = 1;
		}
// Change in temperature sums from emergence with and without day length effect
		crop.DTSUML = crop.DTSUM * RDAYL;
	
// Development stage rate       
		if (crop.DVS < 1.0) { 
			// vegetative phase
			crop.DVR = crop.DTSUML / crop.TSUM1; 
		} else { 
			// generative phase
			crop.DVR = crop.DTSUML / crop.TSUM2;
		}
		   
		   
//// Radiation use efficiency ////

// Radiation use efficiency as dependent on development stage (g DM MJ-1)
		crop.RUE = approx(crop.RUETB, crop.DVS);

		
// Correction of radiation use efficiency for change in atmospheric CO2 concentration (-)
		double RCO = approx(crop.COTB, wth.CO2);
		
// Reduction of radiation use efficiency for day-time temperature and for low minimum temperature
		double DTEMP = wth.tmax[time] - 0.25 * (wth.tmax[time] - wth.tmin[time]);
		double RTMP = approx(crop.TMPFTB, DTEMP) * approx(crop.TMNFTB, wth.tmin[time]);
// Correction of RUE for both non-optimal temperatures and atmospheric CO2
		crop.RUE = crop.RUE * RTMP * RCO;

		

//// GROWTH /// 
		double KDIF = approx(crop.KDIFTB, crop.DVS);
		crop.PARINT = PAR * (1- exp(-KDIF * crop.LAI)); 

// Fractional light interception [-] for resp. PAR and total radiation
		double FINT =  crop.PARINT / PAR;
		double KGLOB = 0.75 * KDIF;
		crop.FINTT =  1 - exp(-KGLOB * crop.LAI);  
		
// Nutrient reduction factor       
// NREF =  exp(-NLUE*(1.0-NNI))
		crop.NPKREF =  clamp(0, 1, 1 - crop.NLUE * pow(1 - crop.NPKI, 2));
// Growth rate GRT in kg DM ha-1 d-1 calculated from PARINT in MJ/m2 and RUE in g/MJ --> multiply with 10
		crop.GRT = 10 * crop.RUE * crop.PARINT * min(crop.TRANRF, crop.NPKREF);
		
//		 cout << crop.PARINT << " - " << crop.GRT << endl;

		
// Water-Nutrient stress factor
// is not used
//		crop.RNW = min(crop.TRANRF, crop.NPKI);
	

// Specific Leaf area(ha/kg), as dependent on nutrient stress.
		double SLA = approx(crop.SLATB, crop.DVS) * exp(-crop.NSLA * (1. - crop.NPKI));
	// CALL GLA(DAY,EMERG,DTEFF,LAII,RGRLAI,DELT,SLA,LAI,GLV,NLAI,DVS,TRANRF,NPKI,GLAI);
	// daily increase of leaf area index (ha leaf area/ ha ground area/ d)                          *


//// Death rates ////	
		
// Relative death rate of roots (d-1)
		double RDRRT = approx(crop.RDRRTB, crop.DVS);
		double RDRST = approx(crop.RDRSTB, crop.DVS);
		
		if (crop.DVS < crop.DVSDR){
			crop.DRRT =  0.;
			crop.DRST =  0.;
		} else {
			crop.DRRT = crop.WRT * RDRRT;
			crop.DRST = crop.WST * RDRST;
		}
				

//// Partitioning ////		
// Biomass partitioning functions under non-stressed situations (-)
		double FRT = approx(crop.FRTB, crop.DVS );
		double FLV   = approx(crop.FLTB, crop.DVS );
		double FST   = approx(crop.FSTB, crop.DVS);
		double FSO   = approx(crop.FOTB, crop.DVS );
	
// Root growth (cm d-1)
		crop.RR = min(crop.RRI * INSW(crop.TRANRF - 0.01, 0., 1.),  soil.RDM - crop.RD);
			
// CALL SUBPAR (NPART,TRANRF,NNI,FRTWET,FLVT,FSTT,FSOT,FSHMOD,FLVMOD,FRT,FLV,FST,FSO);
// Modification of dry matter partitioning to leaves, stems, roots and storage organs in dependence of water and N stress
	
		double FRTMOD, FLVMOD;
		if(crop.TRANRF <= crop.NNI) {
	// Water stress is more severe than nitrogen stress. 
			FRTMOD = std::max(1., 1/(crop.TRANRF + 0.5));
			FRT    = std::min(0.6, FRT * FRTMOD);
			FLV    = FLV * (1-FRT);
			FST    = FST * (1-FRT); 
			FSO    = FSO * (1-FRT);
		} else {
	// Nitrogen stress is more severe as compared to water stress resulting in less partitioning to leaves and more to stems
			double FLVT = FLV;
			FLVMOD = exp(- crop.NPART * (1 - crop.NNI));
			FLV    = FLVT * FLVMOD * (1-FRT);
			FST    = (FST + FLVT - FLV) * (1-FRT) ;
			FSO    = FSO  * (1-FRT);
		}
		
//		double FCHECK = FRT + FLV + FST + FSO;
//		if (FCHECK >= 1.01) cout << "assimilate allocation check over crop organs FCHECK not 0" << endl;
	


//// Death rates ////
	
// Relative death rate of leaves due to senescence/ageing as dependent on mean daily temperature (d-1)
		double RDRTMP = approx(crop.RDRLTB, TMPA);	

// CALL DEATHL(DAY,EMERG,DVS,DVSDLT,RDRTMP,RDRSHM,RDRL,TRANRF,LAI, LAICR,WLVG,RDRNS,NPKI,SLA,RDRDV,RDRSH,RDR,DLV,DLVS,DLVNS,DLAIS,DLAINS,DLAI)
// The relative death rate (d-1) of leaves due to age,shading and drought and due to nutrient (NPK) stress and the death rate of leaves in total (kg ha-1 d-1)                   *
		double RDRDV;
		
		if (crop.DVS < crop.DVSDLT){
			RDRDV = 0.;
		} else {
			RDRDV = RDRTMP;
		}
		
		double RDRSH = std::max(0., crop.RDRSHM * (crop.LAI - crop.LAICR) / crop.LAICR);
		double RDRDRY = (1 - crop.TRANRF) * crop.RDRL;
		double RDR  = maxvalue(std::vector<double> {RDRDV, RDRSH, RDRDRY});
		
		double DLAINS, DLVNS;
		if (crop.NPKI < 1) {
			DLVNS  = crop.WLVG * crop.RDRNS * (1 - crop.NPKI);
			DLAINS = DLVNS * SLA;
		} else {
			DLVNS  = 0;
			DLAINS = 0;
		}
		
		double DLVS, DLAIS;
		DLVS = crop.WLVG * RDR;
		DLAIS = crop.LAI * RDR;
		
		crop.DLV = DLVS + DLVNS;
		double DLAI = DLAIS + DLAINS;

	
// growth rate of roots, leaves, stem and storage organs (kg ha-1 d-1)                        *
		crop.RWRT = crop.GRT * FRT - crop.DRRT;
		crop.RWLVG = crop.GRT * FLV - crop.DLV;
		crop.RWST = crop.GRT * FST - crop.DRST;
		crop.RWSO = crop.GRT * FSO;

	
		if (crop.LAI == 0) {
// Growth at day of seedling emergence:
			crop.GLAI = crop.LAII;
		} else if ((crop.DVS < 0.2) & (crop.LAI < 0.75)) {
// Growth during juvenile stage:
			crop.GLAI = (crop.LAI * (exp(crop.RGRLAI * DTEFF) - 1.))* crop.TRANRF * exp(- crop.NLAI* (1.0 - crop.NPKI));
		} else {
// Growth during maturation stage:
	// Leaf growth
			crop.GLAI = SLA * crop.GRT * FLV;
		}

	// Net rate of change of Leaf area (m2 leaf area m-2 d-1)
		crop.RLAI = crop.GLAI - DLAI;


// Finish conditions
		if ((FINT < 0.05) & (crop.TAGB > 200)) {
			crop.alive = false;
		}	
		
	}	
}


	

void Lintul3Model::crop_ratesNPK() {


// translocatable N/P/K in leaves, stem, roots and storage organs (kg N/P/K ha-1)
// CALL NTRLOC(ANLV,ANST,ANRT,WLVG,WST,WRT,RNFLV,RNFST,RNFRT,FNTRT,ATNLV,ATNST,ATNRT,ATN, APLV,APST,APRT,AKLV,AKST,AKRT, RPFLV,RPFST,RPFRT,RKFLV,RKFST,RKFRT, ATPLV,ATPST,ATPRT,ATP,ATKLV,ATKST,ATKRT,ATK);
		
	double ATNLV = std::max (0., crop.ANLV - crop.WLVG * crop.RNFLV);
	double ATNST = std::max (0., crop.ANST - crop.WST * crop.RNFST);
	double ATNRT = std::min((ATNLV + ATNST) * crop.FNTRT, crop.ANRT - crop.WRT * crop.RNFRT);
	crop.ATN = ATNLV + ATNST + ATNRT;
	
	double ATPLV = std::max (0., crop.APLV - crop.WLVG * crop.RPFLV);
	double ATPST = std::max (0., crop.APST - crop.WST * crop.RPFST);
	double ATPRT = std::min((ATPLV + ATPST) * crop.FNTRT, crop.APRT - crop.WRT * crop.RPFRT);
	crop.ATP = ATPLV +  ATPST + ATPRT;
	
	double ATKLV = std::max (0., crop.AKLV - crop.WLVG * crop.RKFLV);
	double ATKST = std::max (0., crop.AKST - crop.WST * crop.RKFST);
	double ATKRT = std::min((ATKLV + ATKST) * crop.FNTRT, crop.AKRT - crop.WRT * crop.RKFRT);
	crop.ATK = ATKLV + ATKST + ATKRT;  
		
	
// * Total vegetative living above-ground biomass (kg DM ha-1)
	double TBGMR = crop.WLVG + crop.WST;
	
// N/P/K concentrations (kg N/P/K kg-1 DM) in the living leaves, stem, roots and storage organs
/*	double NFLV = crop.ANLV / NOTNUL(crop.WLVG);
	double NFST = crop.ANST / NOTNUL(crop.WST);
	double NFRT = crop.ANRT / NOTNUL(crop.WRT);
	double NFSO = crop.ANSO / NOTNUL(crop.WSO);
	
	double PFLV = crop.APLV / NOTNUL(crop.WLVG);
	double PFST = crop.APST / NOTNUL(crop.WST);
	double PFRT = crop.APRT / NOTNUL(crop.WRT);
	double PFSO = crop.APSO / NOTNUL(crop.WSO);
	
	double KFLV = crop.AKLV / NOTNUL(crop.WLVG);
	double KFST = crop.AKST / NOTNUL(crop.WST);
	double KFRT = crop.AKRT / NOTNUL(crop.WRT);
	double KFSO = crop.AKSO / NOTNUL(crop.WSO);
*/
// Total N/P/K in vegetative living above-ground biomass (kg N/P/K ha-1)
	double NUPGMR = crop.ANLV + crop.ANST;
	double PUPGMR = crop.APLV + crop.APST;
	double KUPGMR = crop.AKLV + crop.AKST;
	
// Fertilizer N/P/K application (kg N/P/K ha-1 d-1) and its recovery fraction (-)--------------*
	double FERTN  = approx(control.FERNTAB, DOY);
	double NRF    = approx(soil.NRFTAB, DOY);
	double FERTNS = FERTN * NRF;
	double FERTP  = approx(control.FERPTAB, DOY);
	double PRF    = approx(soil.PRFTAB, DOY);
	double FERTPS = FERTP * PRF;
	double FERTK  = approx(control.FERKTAB, DOY);
	double KRF    = approx(soil.KRFTAB, DOY);
	double FERTKS = FERTK * KRF;
	
// Check on N balance
//	double NBALAN = abs(crop.NUPTT + crop.NFIXTT + (crop.ANLVI + crop.ANSTI + crop.ANRTI + crop.ANSOI)-(crop.ANLV + crop.ANST + crop.ANRT + crop.ANSO + crop.NLOSSL + crop.NLOSSR + crop.NLOSSS)) ; 
//	if (NBALAN >= 1) cout << "nitrogen balance NBALAN not 0" << endl;
	
// Check on P balance
//	double PBALAN = abs(crop.PUPTT + (crop.APLVI + crop.APSTI + crop.APRTI + crop.APSOI)-(crop.APLV + crop.APST + crop.APRT + crop.APSO + crop.PLOSSL + crop.PLOSSR + crop.PLOSSS));  
//	if (PBALAN >= 1) cout << "phosphorus balance PBALAN not 0" << endl;
	
// Check on K balance
//	double KBALAN = abs(crop.KUPTT + (crop.AKLVI + crop.AKSTI + crop.AKRTI + crop.AKSOI)-(crop.AKLV + crop.AKST + crop.AKRT + crop.AKSO + crop.KLOSSL + crop.KLOSSR + crop.KLOSSS)); 
//	if (KBALAN >= 1) cout << "potassium balance KBALAN not 0" << endl;


// Total N/P/K in living above-ground crop organs (kg N/P/K ha-1)   
//	double NTAG = crop.ANLV + crop.ANST + crop.ANSO;
//	double PTAG = crop.APLV + crop.APST + crop.APSO;
//	double KTAG = crop.AKLV + crop.AKST + crop.AKSO;


// Maximum N/P/K concentration in the leaves, from which the N/P/K conc. in the
// stem and roots are derived, as a function of development stage (kg N/P/K kg-1 DM)  
	double NMAXLV = approx(crop.NMXLV, crop.DVS);
	double PMAXLV = approx(crop.PMXLV, crop.DVS);
	double KMAXLV = approx(crop.KMXLV, crop.DVS);
	
// N/P/K concentration in above-ground living biomass (kg N/P/K kg-1 DM)
//	double NTAC = NTAG/ crop.TAGBG;
//	double PTAC = PTAG/ crop.TAGBG;
//	double KTAC = KTAG/ crop.TAGBG;
	
// N/P/K supply to the storage organs (kg N/P/K ha-1 d-1)
	double NSUPSO = INSW(crop.DVS - crop.DVSNT, 0., crop.ATN / crop.TCNT);
	double PSUPSO = INSW(crop.DVS - crop.DVSNT, 0., crop.ATP / crop.TCPT);
	double KSUPSO = INSW(crop.DVS - crop.DVSNT, 0., crop.ATK / crop.TCKT);
	
// N/P/K concentrations in total vegetative living above-ground biomass  (kg N/P/K kg-1 DM) 
	double NFGMR = NUPGMR / NOTNUL(TBGMR);
	double PFGMR = PUPGMR / NOTNUL(TBGMR);
	double KFGMR = KUPGMR / NOTNUL(TBGMR);
	
// * Residual N/P/K concentrations in total vegetative living above-ground biomass  (kg N/P/K kg-1 DM) 
	double NRMR = (crop.WLVG * crop.RNFLV + crop.WST * crop.RNFST) / NOTNUL(TBGMR);
	double PRMR = (crop.WLVG * crop.RPFLV + crop.WST * crop.RPFST) / NOTNUL(TBGMR);
	double KRMR = (crop.WLVG * crop.RKFLV + crop.WST * crop.RKFST) / NOTNUL(TBGMR);


// ---  Maximum N/P/K concentrations in stems and roots (kg N/P/K kg-1 DM)
	double NMAXST = crop.LSNR * NMAXLV;
	double NMAXRT = crop.LRNR * NMAXLV;
	double PMAXST = crop.LSPR * PMAXLV;
	double PMAXRT = crop.LRPR * PMAXLV;
	double KMAXST = crop.LSKR * KMAXLV;
	double KMAXRT = crop.LRKR * KMAXLV;



// Calling the subroutine for calculating optimal N/P/K concentrations in leaves and stems
// CALL NOPTM(FRNX,NMAXLV,NMAXST,NOPTLV,NOPTST,FRPX,PMAXLV, PMAXST,POPTLV,POPTST,FRKX,KMAXLV,KMAXST,KOPTLV,KOPTST);
// Purpose: To compute the optimal N/P/K concentrations of crop organs (kg N/P/K kg-1 DM)*
	double NOPTLV = crop.FRNX * NMAXLV;
	double NOPTST = crop.FRNX * NMAXST;
	double POPTLV = crop.FRPX * PMAXLV;
	double POPTST = crop.FRPX * PMAXST;
	double KOPTLV = crop.FRKX * KMAXLV;
	double KOPTST = crop.FRKX * KMAXST;
	
// Optimal amount of N/P/K i/n vegetative above-ground living biomass and its N/P/K concentration
	double NOPTS = NOPTST * crop.WST;
	double NOPTL = NOPTLV * crop.WLVG; 
	double NOPTMR = (NOPTL + NOPTS) / NOTNUL(TBGMR);
	
	double POPTS = POPTST * crop.WST;
	double POPTL = POPTLV * crop.WLVG; 
	double POPTMR = (POPTL + POPTS) / NOTNUL(TBGMR);
	
	double KOPTS = KOPTST * crop.WST;
	double KOPTL = KOPTLV * crop.WLVG; 
	double KOPTMR = (KOPTL + KOPTS) / NOTNUL(TBGMR);
	
// Calling the subroutine for calculating the NPK and N/P/K Nutrition Indices (NPKI & NNI etc.)
// CALL NNINDX(DAY,DAYEM,EMERG,NFGMR,NRMR,NOPTMR,NNI,PFGMR,PRMR,POPTMR,PNI,KFGMR,KRMR,KOPTMR,KNI,NPKI);
// Purpose: To compute NPK and N/P/K Nutrition Indices (-)
// double NNINDX(DAY,DAYEM,EMERG,NFGMR,NRMR,NOPTMR,NNI,PFGMR,PRMR,POPTMR,PNI,KFGMR,KRMR,KOPTMR,KNI,NPKI) {
	
	double TINY = 0.001;

	if (crop.EMERG) {
		crop.NNI = clamp(TINY, 1.0, ((NFGMR-NRMR)/NOTNUL(NOPTMR-NRMR)));
		crop.PNI = clamp(TINY, 1.0, ((PFGMR-PRMR)/NOTNUL(POPTMR-PRMR)));
		crop.KNI = clamp(TINY, 1.0, ((KFGMR-KRMR)/NOTNUL(KOPTMR-KRMR)));
		crop.NPKI = minvalue(std::vector<double> {crop.NNI, crop.PNI, crop.KNI} );
	} else {
		crop.NNI = 0.0;
		crop.PNI = 0.0;
		crop.KNI = 0.0;
		crop.NPKI = 0.0;
	}


// Calling the subroutine for N/P/K demand of leaves, roots and stem storage organs (kg N/P/K ha-1 d-1)
// CALL NDEMND(NMAXLV,NMAXST,NMAXRT,NMAXSO,WLVG,WST,WRT,WSO,PMAXLV,PMAXST,PMAXRT,PMAXSO,KMAXLV,KMAXST,KMAXRT,KMAXSO, ANLV,ANST,ANRT,ANSO,TCNT, APLV,APST,APRT,APSO,TCPT, AKLV,AKST,AKRT,AKSO,TCKT,NDEML,NDEMS,NDEMR,NDEMSO,PDEML,PDEMS,PDEMR,PDEMSO,KDEML,KDEMS,KDEMR,KDEMSO);
// Purpose: To compute the N/P/K demands of crop organs (kg N/P/K ha-1) *
	double NDEML = std::max(NMAXLV * crop.WLVG - crop.ANLV, 0.);
	double NDEMS = std::max(NMAXST * crop.WST - crop.ANST, 0.);
	double NDEMR = std::max(NMAXRT * crop.WRT - crop.ANRT, 0.);
	double NDEMSO = std::max(crop.NMAXSO * crop.WSO - crop.ANSO, 0.) /  crop.TCNT;
	
	double PDEML = std::max(PMAXLV * crop.WLVG - crop.APLV, 0.); 
	double PDEMS = std::max(PMAXST * crop.WST - crop.APST, 0.);
	double PDEMR = std::max(PMAXRT * crop.WRT - crop.APRT, 0.);
	double PDEMSO = std::max(crop.PMAXSO * crop.WSO - crop.APSO, 0.) /  crop.TCPT;
	
	double KDEML = std::max(KMAXLV * crop.WLVG - crop.AKLV, 0.);
	double KDEMS = std::max(KMAXST * crop.WST - crop.AKST, 0.);
	double KDEMR = std::max(KMAXRT * crop.WRT - crop.AKRT, 0.);
	double KDEMSO = std::max(crop.KMAXSO * crop.WSO - crop.AKSO, 0.) /  crop.TCKT;
	
// Total N/P/K demand (kg N/P/K ha-1)
	double NDEMTO = max(0.0, (NDEML + NDEMS + NDEMR));
	double PDEMTO = max(0.0, (PDEML + PDEMS + PDEMR));
	double KDEMTO = max(0.0, (KDEML + KDEMS + KDEMR));
	
// Rate of N/P/K uptake in grains (kg N/P/K ha-1 d-1)
	crop.RNSO = std::min(NDEMSO, NSUPSO);
	crop.RPSO = std::min(PDEMSO, PSUPSO);
	crop.RKSO = std::min(KDEMSO, KSUPSO);


	double NUPTR, NFIXTR, PUPTR, KUPTR;

// Nutrient uptake limiting factor (-) at low moisture conditions in the
// rooted soil layer before anthesis. After DVSNLT, there is no nutrient uptake from the soil
	double NLIMIT = INSW(crop.DVS - crop.DVSNLT, INSW(crop.TRANRF - 0.01, 0, 1) , 0);

	
	if (crop.EMERG) {
// Total N/P/K uptake (kg N/P/K ha-1 d-1) from soil and by biological N fixation
		NUPTR = (std::max (0., std::min((1 - crop.NFIXF) * NDEMTO, soil.NMINT)) * NLIMIT);
		NFIXTR = std::max (0., NUPTR * crop.NFIXF / std::max (0.02, 1- crop.NFIXF) );
// 	Total P/K uptake (kg P/K ha-1 d-1) from soil 
		PUPTR = (std::max (0., std::min(PDEMTO, soil.PMINT)) * NLIMIT);
		KUPTR = (std::max (0., std::min(KDEMTO, soil.KMINT)) * NLIMIT);
	
// No N/P/K limitation for optimal and water limited production
		if ((control.IOPT == 1) | (control.IOPT == 2)) {
			NUPTR = (max (0., (1 - crop.NFIXF) * NDEMTO) * NLIMIT ); 
			NFIXTR = (max (0., crop.NFIXF * NDEMTO) * NLIMIT );
			PUPTR = (max (0., PDEMTO) * NLIMIT);
			KUPTR = (max (0., KDEMTO) * NLIMIT);
		} else if (control.IOPT == 3) {
// No P/K limitation for nitrogen limited production
			PUPTR = (max (0., PDEMTO) * NLIMIT);
			KUPTR = (max (0., KDEMTO) * NLIMIT);
		}
	} else {
		NUPTR = 0.;
		NFIXTR = 0.;
		PUPTR = 0.;
		KUPTR = 0.;
	}    
	
// Calling the subroutine for calculating N/P/K translocated from leaves, stem, and roots (kg N/P/K ha-1 d-1)
// CALL NTRANS(RNSO,ATNLV,ATNST,ATNRT,ATN,RNTLV,RNTST,RNTRT,RPSO,ATPLV,ATPST,ATPRT,ATP,RPTLV,RPTST,RPTRT,RKSO,ATKLV,ATKST,ATKRT,ATK,RKTLV,RKTST,RKTRT);
// compute the N/P/K translocated from different organs (kg N/P/K ha-1 d-1) *
	double RNTLV =  crop.RNSO * ATNLV/ NOTNUL(crop.ATN);
	double RNTST =  crop.RNSO * ATNST/ NOTNUL(crop.ATN);
	double RNTRT =  crop.RNSO * ATNRT/ NOTNUL(crop.ATN);
	double RPTLV =  crop.RPSO * ATPLV/ NOTNUL(crop.ATP);
	double RPTST =  crop.RPSO * ATPST/ NOTNUL(crop.ATP);
	double RPTRT =  crop.RPSO * ATPRT/ NOTNUL(crop.ATP);
	double RKTLV =  crop.RKSO * ATKLV/ NOTNUL(crop.ATK);
	double RKTST =  crop.RKSO * ATKST/ NOTNUL(crop.ATK);
	double RKTRT =  crop.RKSO * ATKRT/ NOTNUL(crop.ATK);
	
	
	
// Calling the subroutine to compute the partitioning of the total
// N/P/K uptake rates (NUPTR,PUPTR,KUPTR) over the leaves, stem and roots (kg N/P/K ha-1 d-1)
// CALL RNUSUB(DAY,DAYEM,EMERG,NDEML,NDEMS,NDEMR,NUPTR,PDEML,PDEMS,PDEMR,PUPTR,KDEML,KDEMS,KDEMR,KUPTR,NFIXTR,NDEMTO, RNULV,RNUST,RNURT, PDEMTO, RPULV,RPUST,RPURT,KDEMTO,RKULV,RKUST,RKURT);
// Purpose: To compute the partitioning of the total N/P/K uptake rates (NUPTR,PUPTR,KUPTR) over leaves, stem, and roots (kg N/P/K ha-1 d-1)*
	double RNULV, RNUST, RNURT, RPULV, RPUST, RPURT, RKULV, RKUST, RKURT;

	if (crop.EMERG) {
		RNULV = (NDEML / NOTNUL(NDEMTO)) * (NUPTR + NFIXTR);
		RNUST = (NDEMS / NOTNUL(NDEMTO)) * (NUPTR + NFIXTR);
		RNURT = (NDEMR / NOTNUL(NDEMTO)) * (NUPTR + NFIXTR);
		
		RPULV = (PDEML / NOTNUL(PDEMTO)) * PUPTR;
		RPUST = (PDEMS / NOTNUL(PDEMTO)) * PUPTR;
		RPURT = (PDEMR / NOTNUL(PDEMTO)) * PUPTR;
		
		RKULV = (KDEML / NOTNUL(KDEMTO)) * KUPTR;
		RKUST = (KDEMS / NOTNUL(KDEMTO)) * KUPTR;
		RKURT = (KDEMR / NOTNUL(KDEMTO)) * KUPTR;
	} else {
		RNULV = 0.;
		RNUST = 0.;
		RNURT = 0.;
		RPULV = 0.;
		RPUST = 0.;
		RPURT = 0.;
		RKULV = 0.;
		RKUST = 0.;
		RKURT = 0.;
	}
	
// Soil N/P/K supply (g N m-2 d-1) through mineralization during crop growth
	if (crop.EMERG) {
		soil.RNMINS = - std::max(0., std::min(soil.RTNMINS * soil.NMINI * NLIMIT, soil.NMIN));
		soil.RPMINS = - std::max(0., std::min(soil.RTPMINS * soil.PMINI * NLIMIT, soil.PMIN));
		soil.RKMINS = - std::max(0., std::min(soil.RTKMINS * soil.KMINI * NLIMIT, soil.KMIN));
	}
	
// Change in total inorganic N/P/K in soil as function of fertilizer input, soil N/P/K mineralization and crop uptake.
	soil.RNMINT = FERTNS - NUPTR - soil.RNMINS;
	soil.RPMINT = FERTPS - PUPTR - soil.RPMINS;
	soil.RKMINT = FERTKS - KUPTR - soil.RKMINS;
	
// Rate of change of N/P/K in crop organs   
	crop.RNST = RNUST - RNTST - crop.RNLDST;
	crop.RNRT = RNURT - RNTRT - crop.RNLDRT;
	crop.RNLV = RNULV - RNTLV - crop.RNLDLV;  
	crop.RPST = RPUST - RPTST - crop.RPLDST;
	crop.RPRT = RPURT - RPTRT - crop.RPLDRT;
	crop.RPLV = RPULV - RPTLV - crop.RPLDLV;  
	crop.RKST = RKUST - RKTST - crop.RKLDST;
	crop.RKRT = RKURT - RKTRT - crop.RKLDRT;
	crop.RKLV = RKULV - RKTLV - crop.RKLDLV;   

// CALL RNLD(DVS,WRT,WST,RDRRT,RDRST,RNFLV,DLV,RNFRT,RPFLV,RPFRT,RKFLV,RKFRT,RNFST,RPFST,RKFST,DVSDR,DRRT,DRST,RNLDLV,RNLDRT,RNLDST,RPLDLV,RPLDRT,RPLDST,RKLDLV,RKLDRT,RKLDST);
// N/P/K losses due to death of leaves, stems and roots (kg N/P/K ha-1 d-1)

	crop.RNLDLV =  crop.RNFLV * crop.DLV;
	crop.RNLDRT =  crop.RNFRT * crop.DRRT;
	crop.RNLDST =  crop.RNFST * crop.DRST;
	crop.RPLDLV =  crop.RPFLV * crop.DLV;
	crop.RPLDRT =  crop.RPFRT * crop.DRRT;
	crop.RPLDST =  crop.RPFST * crop.DRST;
	crop.RKLDLV =  crop.RKFLV * crop.DLV;
	crop.RKLDRT =  crop.RKFRT * crop.DRRT;
	crop.RKLDST =  crop.RKFST * crop.DRST;

}	
	
