#include "abc_000_warning.h"

#include "abc_001_config.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "abc_datatype.h"
#include "abc_ide_util.h"
#include "abc_ts_func.h"
#include "abc_date.h"

#include "beastv2_io.h"
#include "globalvars.h"

void BEAST2_print_options(A(OPTIONS_PTR)  opt)
{	
	if (opt->extra.printOptions== 0 || GLOBAL_QUIET_MODE) return;
	A(METADATA_PTR) meta = &(opt->io.meta);
	BEAST2_IO_PTR   io    = &opt->io;
	TimeVecInfo     *tvec  = &io->T;
	I32PTR          dims  = io->dims;
	
	char* yesno[] = { "Yes", "No" };
 
#if M_INTERFACE==1
	char* dateNumOriginStr = "0000-01-00";
	char comment	= '%';
	char filler		= '.';
	char dashdot    = '_';
	char* emptyList = "[]";
	char* logicals[] = { "false", "true" };
	char* nanstr = "NaN";
	char* em1 = ""; "<strong>";
	char* em2 = ""; "</strong>";
#elif R_INTERFACE==1
	char* dateNumOriginStr = "1970-01-01";
	char comment = '#';
	char filler  = '$';
	char dashdot = '.';
	char* emptyList = "list()";
	char* logicals[] = {"FALSE", "TRUE" };
	char* nanstr = "NaN";
	char* em1 = "";
	char* em2 = "";
#elif P_INTERFACE==1
	char* dateNumOriginStr = "0001-01-01";
	char comment = '#';
	char filler = '.';
	char dashdot = '_';
	char* emptyList  = " rb.args() ### or 'lambda: None': just get an empty object###";
	char* logicals[] = { "False", "True" };
	char* nanstr     ="float('nan')";
	char* em1 = "";
	char* em2 = "";
#endif

	I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
	I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
	I08 hasTrendCmpnt   = 1;
	I08 hasAny          = 1;

	// ##__VA_ARGS__ to eat off the trainling ,
	// __VA_OPT__(,)
	//http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2016/p0306r2.html
	//https://stackoverflow.com/questions/52891546/what-does-va-args-mean
	//https://stackoverflow.com/questions/5891221/variadic-macros-with-zero-arguments
	//https://web.archive.org/web/20070616054343/http://developers.sun.com/sunstudio/support/CCcompare.html

	#define Print1(fmtstr,hasComponent,...) if(hasComponent)  r_printf(fmtstr, ##__VA_ARGS__);
	#define Print2(fmtstr,hasComponent,...) if(hasComponent)  r_printf(fmtstr __VA_OPT__(,)  __VA_ARGS__);
	#define Print(fmtstr,hasComponent,...) {if(hasComponent) r_printf(fmtstr,  __VA_ARGS__);}

	Print("%s", hasAny, "\n");	
	Print("%sINFO%s: To supress printing the parameers in beast(),      set print.options = 0 \n", hasAny, em1, em2);
	Print("%sINFO%s: To supress printing the parameers in beast%cirreg(),set print.options = 0 \n", hasAny, em1, em2, dashdot);
	Print("%sINFO%s: To supress printing the parameers in beast123(),   set extra%cprintOptions = 0  \n", hasAny, em1, em2, filler);
		
	Print("%sINFO%s: To supress warning messages in beast(),            set quiet = 1 \n", hasAny, em1, em2);
	Print("%sINFO%s: To supress warning messages in beast%cirreg(),      set quiet = 1 \n", hasAny, em1, em2, dashdot);
	Print("%sINFO%s: To supress warning messages in beast123(),         set extra%cquiet = 1  \n", hasAny, em1, em2, filler);

	Print("%s", hasAny, "\n");
	Print("%c--------------------------------------------------%c\n", hasAny, comment, comment);
	Print("%c       Brief summary of Input Data                %c\n", hasAny, comment, comment);
	Print("%c--------------------------------------------------%c\n", hasAny, comment, comment);	
	Print("%sData Dimension%s: One signal of length %d\n",                    hasAny && io->ndim == 1, em1, em2,    dims[0]); 
	Print("%sData Dimension%s: [%dx%d] - %d signals of length %d each\n",    hasAny && io->ndim == 2, em1, em2, dims[0], dims[1], io->numOfPixels, dims[io->meta.whichDimIsTime-1]);
	Print("%sData Dimension%s: [%dx%dx%d] - %d signals of length %d each\n", hasAny && io->ndim == 3, em1, em2, dims[0], dims[1], dims[2], io->numOfPixels, dims[io->meta.whichDimIsTime - 1]);
	Print("%sIsOrdered%s     : %s\n", hasAny, em1, em2, tvec->IsOrdered ? "Yes, ordered in time": "No, unordered in time, to be sorted/ordered before running BEAST" );

	int unitFormat = 3;
	if (tvec->isDate &&  tvec->out.asDailyTS)  unitFormat = 1;  // analyazed as daily time series
	if (tvec->isDate && !tvec->out.asDailyTS)  unitFormat = 2;  // unit is fractional year
	
	F32  dT;

	char* yesMsg = "Yes, evenly spaced at interval of ";
	char* noMsg  = "No, unevenly spaced at avg interval of ";
	char msg[100];

	dT = io->T.data_dt;
	if (unitFormat == 1) snprintf(msg, 99, "%g days", dT);
	if (unitFormat == 2) snprintf(msg, 99, "%g year = %g months = %g days", dT, dT*12, dT *365);
	if (unitFormat == 3) snprintf(msg, 99, "%g (unknown unit)", dT);
	Print("%sIsRegular%s     : %s %s\n", hasAny, em1, em2, io->T.isRegular  || io->T.isConvertedFromStartDeltaOnly ? yesMsg : noMsg, msg);
	 	
	dT = io->meta.deltaTime;
	if (unitFormat == 1) snprintf(msg, 99, "%g days", dT);
	if (unitFormat == 2) snprintf(msg, 99, "%g year = %g months = %g days", dT, dT * 12, dT * 365);
	if (unitFormat == 3) snprintf(msg, 99, "%g (unknown unit)", dT);
	if (!tvec->out.needReOrder && tvec->out.needAggregate) {
		if (tvec->isRegular==0 && !io->T.isConvertedFromStartDeltaOnly) {
		Print("%sPreprocessing%s : Aggregate irregular data into a regular interval of %s\n", 1L, em1, em2, msg );
		} else if (tvec->isRegular|| io->T.isConvertedFromStartDeltaOnly) {
		Print("%sPreprocessing%s : Resample regular data to a regular interval of %s\n", 1L, em1, em2, msg);
	}
	}
	if (tvec->isDate && tvec->out.asDailyTS == 0 && fabs(1. / 365 - dT) < dT * 1e-3) {
		Print("%sLeap Years   %s : 366 days in leap years reduced to 365 days\n", 1L, em1, em2);
	}
	
		
	dT = meta->deltaTime;
	if (meta->hasSeasonCmpnt) {
		F32 period = meta->period * dT;
		if (unitFormat == 1) snprintf(msg, 99, "%g days", period);
		if (unitFormat == 2) snprintf(msg, 99, "%g year = %g months = %g days", period, period * 12, period * 365);
		if (unitFormat == 3) snprintf(msg, 99, "%g (unknown unit)", period);

	Print("%sHasSeasonCmpnt%s: %-5s | period = %s. The model 'Y=Trend+Season+Error' is fitted.\n", hasAny, em1, em2, logicals[!!meta->hasSeasonCmpnt], msg);
	Print("              : Num_of_DataPoints_per_Period = period/deltaTime = %g/%g = %g\n", hasAny, meta->period* dT,dT, meta->period);
	Print("%sHasOutlierCmpt%s: %-5s | If true, Y=Trend+Season+Outlier+Error fitted instead of Y=Trend+Season+Error\n", hasAny, em1, em2, logicals[!!meta->hasOutlierCmpnt]);
	} else{
	Print("%shasSeasonCmpnt%s: %-5s | no periodic or seasonal component. The model Y=Trend+Error is fitted.\n", hasAny, em1, em2, logicals[!!meta->hasSeasonCmpnt]);
	Print("%sHasOutlierCmpt%s: %-5s | If true, Y=Trend+Outlier+Error (experimental) is fitted instead of Y=Trend+Error \n", hasAny, em1, em2, logicals[!!meta->hasOutlierCmpnt]);		
	}
	Print("%sDeseasonalize%s : %-5s | If true, remove a global seasonal  cmpnt before running BEAST & add it back after BEAST\n", hasSeasonCmpnt, em1, em2, logicals[!!meta->deseasonalize]);
	Print("%sDetrend%s       : %-5s | If true, remove a global trend component before running BEAST & add it back after BEAST\n", hasAny, em1, em2, logicals[!!meta->detrend]);

	if (IsNaN(meta->missingValue)) msg[0] = 0;
	else                           snprintf(msg, 99, " or %g ", meta->missingValue);
	Print("%sMissingValue%s  : NaN %s flagged as missing values \n", 1L, em1, em2, msg);
	Print("%sMaxMissingRate%s: if more than %g%% of data is missing, BEAST will skip it.\n", 1L, em1, em2, meta->maxMissingRate * 100);

	Print("%s", hasAny, "\n");

	//////////////////////////////////////////////////////////////////////////////////////////////////
	//  Below is for Metadata, Prior, MCMC, and Extra
	//////////////////////////////////////////////////////////////////////////////////////////////////

	Print("%s", hasAny,"\n");
	Print("%c--------------------------------------------------%c\n", hasAny,comment,comment);
	Print("%c      OPTIONS used in the MCMC inference          %c\n", hasAny, comment, comment);
	//Print("%c--------------------------------------------------%c\n", hasAny, comment, comment);
	Print("%c--------------------------------------------------%c\n", hasAny, comment, comment);
	Print("%s", hasAny, "\n");	
	Print("%c......Start of displaying 'MetaData' ......\n", hasAny, comment); 
	Print("metadata                = %-10s %c metadata is used to interpret the input data Y\n",	hasAny, emptyList,comment);
	//Print("  metadata%cisRegular        = %s\n", hasAny, filler, logicals[!!meta->isRegular]);
	//Print("  metadata%cisOrdered        = %s\n", hasAny, filler, logicals[!!meta->isOrdered]);
	//Print("  metadata%cneedAggregate    = %s\n",	hasAny,		filler, logicals[!!meta->needAggregate]);
	Print("metadata%cseason         = 'none'     %c trend-only data with no periodic variation\n",     !hasSeasonCmpnt, filler, comment)                                //no seasonal cmpnt
	Print("metadata%cseason         = 'harmonic' %c fit a harmonic model to the periodic component\n", hasSeasonCmpnt&&meta->seasonForm =='S', filler, comment)  //haromic season basis
	Print("metadata%cseason         = 'dummy'    %c fit a dummy model to the periodic component\n",    hasSeasonCmpnt&&meta->seasonForm =='D', filler, comment)
	Print("metadata%cseason         = 'svd'      %c fit the periodic component with singlular-vector-decomposition bases \n",      hasSeasonCmpnt&&meta->seasonForm =='V', filler, comment)
	

	dT = meta->deltaTime;
	char* extraPeriodMsg = meta->IsPeriodEstimated ? "[Guessed by autocorrelation]" : "";

	if (tvec->isDate) {
		if (tvec->out.asDailyTS) {			
			F64 start  = JDN_to_DateNum( (int) meta->startTime);
			F64 fyear  = FracYear_from_DateNum(round(meta->startTime));

			Print("metadata%cstartTime      = %-10g %c unit: datenum,  %g days lapsed since %s (decimal year of %g)\n", hasAny, filler, start, comment, start, dateNumOriginStr, fyear);
			Print("metadata%cdeltaTime      = %-10g %c unit: days\n", hasAny, filler,  dT, comment);
			Print("metadata%cperiod         = %-10g %c unit: days %s\n", hasSeasonCmpnt, filler, meta->period * dT, comment , extraPeriodMsg);
		}
		else {
			int yr, mon, day;
			FracYear_to_YMD(meta->startTime, &yr, &mon, &day);			
			Print("metadata%cstartTime      = %-10g %c %04d-%02d-%02d\n", hasAny, filler, meta->startTime, comment, yr, mon, day);
			Print("metadata%cdeltaTime      = %-10g %c %g year(s) = %g month(s) = %g day(s)\n",    hasAny, filler, dT, comment, dT, dT * 12, dT * 365);
			Print("metadata%cperiod         = %-10g %c %g year(s) = %g month(s) = %g day(s) %s\n", hasSeasonCmpnt, filler, meta->period * dT, comment, meta->period * dT, meta->period * dT * 12, meta->period * dT * 365, extraPeriodMsg);
		}
	} else {			
			Print("metadata%cstartTime      = %-10g %c unknown unit\n",    hasAny, filler, meta->startTime, comment);
			Print("metadata%cdeltaTime      = %-10g %c unknown unit\n",    hasAny, filler, dT, comment);
			Print("metadata%cperiod         = %-10g %c unknown unit %s\n", hasSeasonCmpnt, filler, meta->period * dT, comment, extraPeriodMsg);
	}

	Print("metadata%cwhichDimIsTime = %d\n",   hasAny && io->ndim > 1, filler, meta->whichDimIsTime); 
	Print("metadata%cmissingValue   = %g\n",   meta->missingValue==meta->missingValue, filler, meta->missingValue);	
	//Print("  metadata%cmissingValue   = %s\n",   meta->missingValue!=meta->missingValue, filler, nanstr);
	Print("metadata%cmaxMissingRate = %-10g %c if more than %g%% of data is missing, BEAST will skip it.\n", hasAny,	filler, meta->maxMissingRate, comment, meta->maxMissingRate * 100);

	Print("metadata%cdeseasonalize  = %-10s %c if true,remove a global seasonal cmpnt before running BEAST & add it back later\n",   hasSeasonCmpnt, filler, logicals[!!meta->deseasonalize],comment);
	Print("metadata%cdetrend        = %-10s %c if true,remove a global trend  cmpnt before running BEAST & add it back later\n",   hasAny,         filler, logicals[!!meta->detrend], comment);
	Print("%c........End of displaying MetaData ........\n\n",	hasAny,  comment);

	A(PRIOR_PTR) prior = &(opt->prior);

	Print("%c......Start of displaying 'prior' ......\n", hasAny, comment);
	Print("prior                   = %-10s %c prior is the true model parameters of BEAST\n", hasAny, emptyList, comment);
	Print("prior%cseasonMinOrder    = %-10d %c sorder.minmax[1]: min harmonic order alllowed\n", hasSeasonCmpnt, filler, prior->seasonMinOrder,comment );
	Print("prior%cseasonMaxOrder    = %-10d %c sorder.minmax[2]: max harmonic order alllowed\n", hasSeasonCmpnt, filler, prior->seasonMaxOrder,comment);
	Print("prior%cseasonMinKnotNum  = %-10d %c scp.minmax[1]   : min num of seasonal chngpts allowed\n", hasSeasonCmpnt, filler, prior->seasonMinKnotNum, comment);
	Print("prior%cseasonMaxKnotNum  = %-10d %c scp.minmax[2]   : max num of seasonal chngpts allowed\n", hasSeasonCmpnt, filler, prior->seasonMaxKnotNum, comment);
	Print("prior%cseasonMinSepDist  = %-10d %c sseg.min        : min seasonal segment length in terms of datapoints\n", hasSeasonCmpnt, filler, prior->seasonMinSepDist,comment);
	Print("prior%cseasonLeftMargin  = %-10d %c sseg.leftmargin : no season chngpts in the first %d datapoints\n", hasSeasonCmpnt, filler, prior->seasonLeftMargin, comment,prior->seasonLeftMargin);
	Print("prior%cseasonRightMargin = %-10d %c sseg.rightmargin: no seoson chngpts in the last %d datapoints\n", hasSeasonCmpnt, filler, prior->seasonRightMargin, comment, prior->seasonRightMargin);

	Print("prior%ctrendMinOrder     = %-10d %c torder.minmax[1]: min trend polynomial order alllowed\n", hasTrendCmpnt,  filler, prior->trendMinOrder, comment);
	Print("prior%ctrendMaxOrder     = %-10d %c torder.minmax[2]: max trend polynomial order alllowed\n", hasTrendCmpnt,  filler, prior->trendMaxOrder, comment);
	Print("prior%ctrendMinKnotNum   = %-10d %c tcp.minmax[1]   : min num of chngpts in trend allowed\n", hasTrendCmpnt,  filler, prior->trendMinKnotNum, comment);
	Print("prior%ctrendMaxKnotNum   = %-10d %c tcp.minmax[2]   : max num of chngpts in trend allowed\n", hasTrendCmpnt,  filler, prior->trendMaxKnotNum, comment);
	Print("prior%ctrendMinSepDist   = %-10d %c tseg.min        : min trend segment length in terms of datapoints\n", hasTrendCmpnt,  filler, prior->trendMinSepDist, comment);
	Print("prior%ctrendLeftMargin   = %-10d %c tseg.leftmargin : no trend chngpts in the first %d datapoints\n", hasTrendCmpnt, filler, prior->trendLeftMargin, comment, prior->trendLeftMargin);
	Print("prior%ctrendRightMargin  = %-10d %c tseg.rightmargin: no trend chngpts in the last %d datapoints\n", hasTrendCmpnt, filler, prior->trendRightMargin, comment, prior->trendRightMargin);


	Print("prior%coutlierMaxKnotNum = %-10d %c ocp             : max num of datapoints treated as outliers\n", hasOutlierCmpnt, filler, prior->outlierMaxKnotNum, comment);
	Print("prior%coutlierSigFactor  = %-10g %c above which datapoints considered outliers\n", hasOutlierCmpnt, filler, prior->outlierSigFactor, comment);
	Print("prior%cK_MAX             = %-10d %c max number of terms in general linear model (relevant only at small values)\n", hasAny,      filler, prior->K_MAX, comment);
	Print("prior%cprecValue         = %-10g %c useful mainly when precPriorType='constant'\n", hasAny,      filler, prior->precValue, comment);
	Print("prior%cmodelPriorType    = %-10d\n", hasAny, filler, prior->modelPriorType);
	//Print("   prior%calpha1           = %f\n", hasAny,      filler, prior->alpha1 );
	//Print("   prior%calpha2           = %f\n", hasAny,      filler, prior->alpha2);
	//Print("   prior%cdelta1           = %f\n", hasAny,      filler, prior->delta1);
	//Print("   prior%cdelta2           = %f\n", hasAny,      filler, prior->delta2);	
	
	if (prior->precPriorType == ConstPrec)
	Print("prior%cprecPriorType     = '%s'\n", hasAny, filler, "constant")
	else if (prior->precPriorType == UniformPrec)
	Print("prior%cprecPriorType     = '%s'\n", hasAny, filler, "uniform")
	else if (prior->precPriorType == ComponentWise)
	Print("prior%cprecPriorType     = '%s'\n", hasAny, filler, "componentwise")
	else if (prior->precPriorType == OrderWise)
	Print("prior%cprecPriorType     = '%s'\n", hasAny, filler, "orderwise");

	Print("%c......End of displaying prior ......\n\n", hasAny, comment);

 
	A(MCMC_PTR) mcmc = &(opt->mcmc);

	Print("%c......Start of displaying 'mcmc' ......\n", hasAny, comment);
	Print("mcmc                           = %-10s %c mcmc is not BEAST parameters but MCMC sampler options\n",	 hasAny, emptyList,comment);
	Print("mcmc%cseed                      = %-10" PRIu64 " %c A nonzero seed to replicate among runs\n", hasAny, filler, mcmc->seed, comment);
	Print("mcmc%csamples                   = %-10d %c Number of samples saved per chain: the larger, the better\n", hasAny, filler, mcmc->samples,comment);
	Print("mcmc%cthinningFactor            = %-10d %c Thinning the chain: the larger, the better \n", hasAny, filler, mcmc->thinningFactor, comment);
	Print("mcmc%cburnin                    = %-10d %c Number of initial samples discarded: the larger, the better\n", hasAny, filler, mcmc->burnin, comment);
	Print("mcmc%cchainNumber               = %-10d %c Number of chains: the larger, the better\n", hasAny, filler, mcmc->chainNumber, comment);
	Print("mcmc%cmaxMoveStepSize           = %-10d %c Max step of jumping from current changepoint: No need to change\n", hasAny, filler, mcmc->maxMoveStepSize,comment);
	Print("mcmc%ctrendResamplingOrderProb  = %-10g %c Proposal probability of sampling trend polynominal order \n", hasTrendCmpnt, filler,  mcmc->trendResamplingOrderProb, comment);
	Print("mcmc%cseasonResamplingOrderProb = %-10g %c Proposal probability of sampling seasoanl order \n", hasSeasonCmpnt,filler, mcmc->seasonResamplingOrderProb,comment);
	Print("mcmc%ccredIntervalAlphaLevel    = %-10g %c The alphal level for Credible Intervals\n", hasAny, filler, mcmc->credIntervalAlphaLevel,comment);
	Print("%c Total number of models randomly visited in BEAST is (burnin+sampples*thinFactor)*chainNumber=%d\n", hasAny, comment, (mcmc->burnin+mcmc->thinningFactor*mcmc->samples)*mcmc->chainNumber);
	Print("%c......End of displaying mcmc ......\n\n", hasAny, comment);

 
	A(EXTRA_PTR) extra = &(opt->extra);

	Print("%c......Start of displaying 'extra' ......\n", hasAny, comment); 
	Print("extra                      = %-5s %c extra is used to configure output/computing options\n", hasAny, emptyList,comment);
	Print("extra%cdumpInputData        = %-5s %c if true, dump a copy of the input data as o%cdata \n", hasAny,      filler, logicals[!!extra->dumpInputData], comment, filler);
	Print("extra%cwhichOutputDimIsTime = %-5d %c 1,2 or 3; which dim of the result is time; used for a 2D/3D input Y\n", hasAny,	   filler, extra->whichOutputDimIsTime,comment );
	Print("extra%ccomputeCredible      = %-5s %c if true, compute  credibiel interval of estimated Y (e.g., o%ctrend%cCI)\n", hasAny,      filler, logicals[!!extra->computeCredible], comment,filler, filler);
	Print("extra%cfastCIComputation    = %-5s %c if true, do not sort but approximiate CI \n", hasAny,      filler, logicals[!!extra->fastCIComputation], comment);
	Print("extra%ccomputeSeasonOrder   = %-5s %c if true, dump the estimated time-varying seasonal order: o.season.order \n", hasSeasonCmpnt, filler, logicals[!!extra->computeSeasonOrder], comment);
	Print("extra%ccomputeTrendOrder    = %-5s %c if true, dump the estimated trend polynomial order \n", hasTrendCmpnt,  filler, logicals[!!extra->computeTrendOrder], comment);
	Print("extra%ccomputeSeasonChngpt  = %-5s %c if true, dump the seasoanl changepoints (scp) in the output \n", hasSeasonCmpnt, filler, logicals[!!extra->computeSeasonChngpt], comment);
	Print("extra%ccomputeTrendChngpt   = %-5s %c if true, dump the trend changepoints (tcp) in the output \n", hasTrendCmpnt,  filler, logicals[!!extra->computeTrendChngpt], comment);
	Print("extra%ccomputeOutlierChngpt = %-5s %c if true, dump the outliers cangepoints (ocp) in the output \n", hasOutlierCmpnt,filler, logicals[!!extra->computeOutlierChngpt], comment);
	Print("extra%ccomputeSeasonAmp     = %-5s %c  compute time-varying seasonal mangitude if season=harmonic  \n", hasSeasonCmpnt, filler, logicals[!!extra->computeSeasonAmp], comment);
	Print("extra%ccomputeTrendSlope    = %-5s %c if true, dump the time-varying slope in trend\n", hasTrendCmpnt,  filler, logicals[!!extra->computeTrendSlope], comment);
	
	Print("extra%ctallyPosNegSeasonJump= %-5s %c differentiate postive/negative jumps at scp\n", hasSeasonCmpnt, filler, logicals[!!extra->tallyPosNegSeasonJump], comment);
	Print("extra%ctallyPosNegTrendJump = %-5s %c differentiate postive/negative jumps at tcp\n", hasTrendCmpnt,  filler, logicals[!!extra->tallyPosNegTrendJump], comment);
	Print("extra%ctallyIncDecTrendJump = %-5s %c differentiate increased/decreased slopes at tcp\n", hasTrendCmpnt,  filler, logicals[!!extra->tallyIncDecTrendJump], comment);
	Print("extra%ctallyPosNegOutliers  = %-5s %c differentiate postive/negative outliers\n", hasOutlierCmpnt,filler, logicals[!!extra->tallyPosNegOutliers], comment);

	Print("extra%cprintProgressBar     = %-5s %c if true, show an ascii progressbar\n", hasAny, filler, logicals[!!extra->printProgressBar], comment);
	Print("extra%cprintOptions         = %-5s %c if true, print the option of the BEAST run\n", hasAny, filler, logicals[!!extra->printOptions], comment);
	Print("extra%cconsoleWidth         = %-5d %c an integer specifying the console width for printing\n", hasAny, filler, extra->consoleWidth, comment);
	//Print("  extra%cnumCPUCoresToUse     = %d\n", filler, extra->numCPUCoresToUse, hasAny);
	Print("extra%cnumThreadsPerCPU     = %-5d %c each cpu core spawns %d concurrent threads (for beast123())\n", hasAny, filler, extra->numThreadsPerCPU, comment, extra->numThreadsPerCPU);
	Print("extra%cnumParThreads        = %-5d %c total number of threads (for beast123() only)\n", hasAny, filler, extra->numParThreads, comment);
	Print("%c......End of displaying extra ......\n\n", hasAny, comment);
 
#if M_INTERFACE==1
	extern void matlab_IOflush(void);
	matlab_IOflush();
#endif
}

#include "abc_000_warning.h"