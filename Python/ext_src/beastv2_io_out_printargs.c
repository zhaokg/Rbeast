#include "abc_000_warning.h"

#include "abc_001_config.h"


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "abc_datatype.h"
#include "abc_ide_util.h"
#include "abc_ts_func.h"

#include "beastv2_io.h"


void BEAST2_print_options(A(OPTIONS_PTR)  opt)
{	
	if (opt->extra.printOptions== 0) return;

#if M_INTERFACE==1
	char comment	= '%';
	char filler		= '.';
	char* emptyList = "[]";
	char* logicals[] = { "false", "true" };
	char* nanstr = "NaN";
#elif R_INTERFACE==1
	char comment = '#';
	char filler  = '$';
	char* emptyList = "list()";
	char* logicals[] = {"FALSE", "TRUE" };
	char* nanstr = "NaN";
#elif P_INTERFACE==1
	char comment = '#';
	char filler = '.';
	char* emptyList  = " rb.args() ### or 'lambda: None': just get an empty object###";
	char* logicals[] = { "False", "True" };
	char* nanstr     ="float('nan')";
#endif

	I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
	I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
	I08 hasTrendCmpnt   = 1;
	I08 hasAlways       = 1;

	// ##__VA_ARGS__ to eat off the trainling ,
	// __VA_OPT__(,)
	//http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2016/p0306r2.html
	//https://stackoverflow.com/questions/52891546/what-does-va-args-mean
	//https://stackoverflow.com/questions/5891221/variadic-macros-with-zero-arguments
	//https://web.archive.org/web/20070616054343/http://developers.sun.com/sunstudio/support/CCcompare.html

	#define Print1(fmtstr,hasComponent,...) if(hasComponent)  r_printf(fmtstr, ##__VA_ARGS__);
	#define Print2(fmtstr,hasComponent,...) if(hasComponent)  r_printf(fmtstr __VA_OPT__(,)  __VA_ARGS__);
	#define Print(fmtstr,hasComponent,...) {if(hasComponent) r_printf(fmtstr,  __VA_ARGS__);}

	Print("%s", hasAlways,"\n");
	Print("%c--------------------------------------------------%c\n", hasAlways,comment,comment);
	Print("%c      OPTIONS used in the MCMC inference          %c\n", hasAlways, comment, comment);
	Print("%c--------------------------------------------------%c\n", hasAlways, comment, comment);
	Print("%c  Set extra%cprintOptions=0 to suppress printing   %c\n", hasAlways, comment,filler,comment);
	Print("%c--------------------------------------------------%c\n", hasAlways, comment, comment);
	Print("%s", hasAlways, "\n");

 	A(METADATA_PTR) meta = &(opt->io.meta);
 
	
	Print("%c......Start of displaying 'MetaData' ......\n", hasAlways, comment);
	Print("  metadata = %s\n",								hasAlways, emptyList);
	Print("  metadata%cisRegularOrdered = %s\n",	hasAlways,		filler, logicals[!!meta->isRegularOrdered]);
	Print("  metadata%cseason           = '%s'\n", !hasSeasonCmpnt, filler, "none")                                //no seasonal cmpnt
	Print("  metadata%cseason           = '%s'\n", hasSeasonCmpnt&&meta->seasonForm =='S', filler, "harmonic")  //haromic season basis
	Print("  metadata%cseason           = '%s'\n", hasSeasonCmpnt&&meta->seasonForm =='D', filler, "dummy")
	Print("  metadata%cperiod           = %f\n",   hasSeasonCmpnt, filler,  meta->period*meta->deltaTime);
	Print("  metadata%cstartTime        = %.5f\n", hasAlways,		filler,  meta->startTime);
	Print("  metadata%cdeltaTime        = %.5f\n", hasAlways,		filler,  meta->deltaTime);
	Print("  metadata%cwhichDimIsTime   = %d\n",   hasAlways,      filler, meta->whichDimIsTime);	
	Print("  metadata%cmissingValue     = %.5f\n",	meta->missingValue==meta->missingValue, filler, meta->missingValue);	
	Print("  metadata%cmissingValue     = %s\n",   meta->missingValue!=meta->missingValue, filler, nanstr);
	Print("  metadata%cmaxMissingRate   = %.4f\n",	hasAlways,      filler, meta->maxMissingRate);
	Print("  metadata%cdeseasonalize    = %s\n",   hasSeasonCmpnt, filler, logicals[!!meta->deseasonalize]);
	Print("  metadata%cdetrend          = %s\n",  hasAlways,      filler, logicals[!!meta->detrend]);
	Print("%c........End of displaying MetaData ........\n\n",	hasAlways,  comment);

	A(PRIOR_PTR) prior = &(opt->prior);

	Print("%c......Start of displaying 'prior' ......\n", hasAlways, comment);
	Print("  prior = %s\n",							  hasAlways, emptyList);
	Print("  prior%cmodelPriorType	  = %d\n", hasAlways,      filler, prior->modelPriorType );
	Print("  prior%cseasonMinOrder   = %d\n", hasSeasonCmpnt, filler, prior->seasonMinOrder );
	Print("  prior%cseasonMaxOrder   = %d\n", hasSeasonCmpnt, filler, prior->seasonMaxOrder);
	Print("  prior%cseasonMinKnotNum = %d\n", hasSeasonCmpnt, filler, prior->seasonMinKnotNum);
	Print("  prior%cseasonMaxKnotNum = %d\n", hasSeasonCmpnt, filler, prior->seasonMaxKnotNum);
	Print("  prior%cseasonMinSepDist = %d\n", hasSeasonCmpnt, filler, prior->seasonMinSepDist);
	Print("  prior%ctrendMinOrder	  = %d\n", hasTrendCmpnt,  filler, prior->trendMinOrder);
	Print("  prior%ctrendMaxOrder	  = %d\n", hasTrendCmpnt,  filler, prior->trendMaxOrder );
	Print("  prior%ctrendMinKnotNum  = %d\n", hasTrendCmpnt,  filler, prior->trendMinKnotNum);
	Print("  prior%ctrendMaxKnotNum  = %d\n", hasTrendCmpnt,  filler, prior->trendMaxKnotNum);
	Print("  prior%ctrendMinSepDist  = %d\n", hasTrendCmpnt,  filler, prior->trendMinSepDist);
	Print("  prior%coutlierMaxKnotNum= %d\n", hasOutlierCmpnt, filler, prior->outlierMaxKnotNum);
	Print("  prior%coutlierSigFactor = %.3f\n", hasOutlierCmpnt, filler, prior->outlierSigFactor);
	Print("  prior%cK_MAX            = %d\n", hasAlways,      filler, prior->K_MAX);
	Print("  prior%cprecValue        = %f\n", hasAlways,      filler, prior->precValue );
	//Print("   prior%calpha1           = %f\n", hasAlways,      filler, prior->alpha1 );
	//Print("   prior%calpha2           = %f\n", hasAlways,      filler, prior->alpha2);
	//Print("   prior%cdelta1           = %f\n", hasAlways,      filler, prior->delta1);
	//Print("   prior%cdelta2           = %f\n", hasAlways,      filler, prior->delta2);	
	
	if (prior->precPriorType == ConstPrec)
	Print("  prior%cprecPriorType    = '%s'\n", hasAlways, filler, "constant")
	else if (prior->precPriorType == UniformPrec)
	Print("  prior%cprecPriorType    = '%s'\n", hasAlways, filler, "uniform")
	else if (prior->precPriorType == ComponentWise)
	Print("  prior%cprecPriorType    = '%s'\n", hasAlways, filler, "componentwise")
	else if (prior->precPriorType == OrderWise)
	Print("  prior%cprecPriorType    = '%s'\n", hasAlways, filler, "orderwise");

	Print("%c......End of displaying pripr ......\n\n", hasAlways, comment);

 
	A(MCMC_PTR) mcmc = &(opt->mcmc);

	Print("%c......Start of displaying 'mcmc' ......\n", hasAlways, comment);
	Print("  mcmc = %s\n",								 hasAlways, emptyList);
	Print("  mcmc%cseed                      = %d\n", hasAlways, filler, mcmc->seed);
	Print("  mcmc%csamples                   = %d\n", hasAlways, filler, mcmc->samples);
	Print("  mcmc%cthinningFactor            = %d\n", hasAlways, filler, mcmc->thinningFactor);
	Print("  mcmc%cburnin                    = %d\n", hasAlways, filler, mcmc->burnin);
	Print("  mcmc%cchainNumber               = %d\n", hasAlways, filler, mcmc->chainNumber);
	Print("  mcmc%cmaxMoveStepSize           = %d\n", hasAlways, filler, mcmc->maxMoveStepSize);
	Print("  mcmc%ctrendResamplingOrderProb  = %.4f\n", hasTrendCmpnt, filler,  mcmc->trendResamplingOrderProb);
	Print("  mcmc%cseasonResamplingOrderProb = %.4f\n", hasSeasonCmpnt,filler, mcmc->seasonResamplingOrderProb );
	Print("  mcmc%ccredIntervalAlphaLevel    = %.3f\n", hasAlways, filler, mcmc->credIntervalAlphaLevel);
	Print("%c......End of displaying mcmc ......\n\n", hasAlways, comment);

 
	A(EXTRA_PTR) extra = &(opt->extra);

	Print("%c......Start of displaying 'extra' ......\n", hasAlways, comment);
	Print("  extra = %s\n", hasAlways, emptyList);
	Print("  extra%cdumpInputData        = %s\n", hasAlways,      filler, logicals[!!extra->dumpInputData]);
	Print("  extra%cwhichOutputDimIsTime = %d\n", hasAlways,	   filler, extra->whichOutputDimIsTime );
	Print("  extra%ccomputeCredible      = %s\n", hasAlways,      filler, logicals[!!extra->computeCredible]);
	Print("  extra%cfastCIComputation    = %s\n", hasAlways,      filler, logicals[!!extra->fastCIComputation] );
	Print("  extra%ccomputeSeasonOrder   = %s\n", hasSeasonCmpnt, filler, logicals[!!extra->computeSeasonOrder]);
	Print("  extra%ccomputeTrendOrder    = %s\n", hasTrendCmpnt,  filler, logicals[!!extra->computeTrendOrder]  );
	Print("  extra%ccomputeSeasonChngpt  = %s\n", hasSeasonCmpnt, filler, logicals[!!extra->computeSeasonChngpt]);
	Print("  extra%ccomputeTrendChngpt   = %s\n", hasTrendCmpnt,  filler, logicals[!!extra->computeTrendChngpt] );
	Print("  extra%ccomputeOutlierChngpt = %s\n", hasOutlierCmpnt,filler, logicals[!!extra->computeOutlierChngpt]  );
	Print("  extra%ccomputeSeasonAmp     = %s\n", hasSeasonCmpnt, filler, logicals[!!extra->computeSeasonAmp]);
	Print("  extra%ccomputeTrendSlope    = %s\n", hasTrendCmpnt,  filler, logicals[!!extra->computeTrendSlope]);
	
	Print("  extra%ctallyPosNegSeasonJump= %s\n", hasSeasonCmpnt, filler, logicals[!!extra->tallyPosNegSeasonJump]);
	Print("  extra%ctallyPosNegTrendJump = %s\n", hasTrendCmpnt,  filler, logicals[!!extra->tallyPosNegTrendJump]);
	Print("  extra%ctallyIncDecTrendJump = %s\n", hasTrendCmpnt,  filler, logicals[!!extra->tallyIncDecTrendJump]);
	Print("  extra%ctallyPosNegOutliers  = %s\n", hasOutlierCmpnt,filler, logicals[!!extra->tallyPosNegOutliers] );

	Print("  extra%cprintProgressBar     = %s\n", hasAlways, filler, logicals[!!extra->printProgressBar]);
	Print("  extra%cprintOptions         = %s\n", hasAlways, filler, logicals[!!extra->printOptions]);
	Print("  extra%cconsoleWidth         = %d\n", hasAlways, filler, extra->consoleWidth);
	//Print("  extra%cnumCPUCoresToUse     = %d\n", filler, extra->numCPUCoresToUse, hasAlways);
	Print("  extra%cnumThreadsPerCPU     = %d\n", hasAlways, filler, extra->numThreadsPerCPU);
	Print("  extra%cnumParThreads        = %d\n", hasAlways, filler, extra->numParThreads);
	Print("%c......End of displaying extra ......\n\n", hasAlways, comment);
 
#if M_INTERFACE==1
	extern void matlab_IOflush(void);
	matlab_IOflush();
#endif
}

#include "abc_000_warning.h"