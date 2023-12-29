#include <math.h>
#include <string.h>
#include "abc_000_warning.h"

#include "abc_common.h"
#include "abc_001_config.h"
#include "abc_sort.h"
#include "abc_vec.h" //f32_sumfilter
#include "abc_blas_lapack_lib.h"
#include "abc_rand_pcg_local.h"


/*Discarded code*/
//clock_t tic, toc; tic = clock(); toc = clock();   
/**********************************************/
static INLINE int isSpace(char c) {return (c == ' ' || c == '\t' || c == '\n');}

int get_word(char* str) {
	int i = 0, wordlen = 0;
	while (isSpace(str[i])) {
		str[i] = ' ';
		wordlen++;
		i++;
	}
	while (str[i] && !isSpace(str[i++])) { wordlen++;}

	return wordlen;}

// http://www.shedai.net/c/new/#section1group25
char* word_wrap(char* str, int LINE_LENGTH) {

	int wordLen;                     /* length of current word */
	int current_lineLen = 0;             /* current length of line */
	int start_curLen = 0;              /* index of beginning if line */

	while ((wordLen = get_word(str + start_curLen+ current_lineLen)) > 0) {
		if (current_lineLen + wordLen < LINE_LENGTH) {
			current_lineLen += wordLen;
		} else {
			str[start_curLen + current_lineLen] = '\n';
			start_curLen += current_lineLen + 1L;
			current_lineLen = 0L;
		}

	}
	return str;
}

char* word_wrap_indented(char* str, int LINE_LENGTH, int nspace) {
	char* space = "                      ";
	char newstr[2000];
	int  newlen = 0;
	int  isfirstline = 1;

	int wordLen;                     /* length of current word */
	int current_lineLen = 0;             /* current length of line */
	int start_curLen = 0;              /* index of beginning if line */

	while ((wordLen = get_word(str + start_curLen + current_lineLen)) > 0) {

		if (current_lineLen + wordLen < LINE_LENGTH) {
			current_lineLen += wordLen;
		}
		else {

			str[start_curLen + current_lineLen] = '\n';

			if (isfirstline) {
				isfirstline = 0;
				memcpy(newstr + newlen, str + start_curLen, current_lineLen);
				newlen += current_lineLen;
				newstr[newlen++] = '\n';
				LINE_LENGTH -= nspace;
			}
			else {
				memcpy(newstr + newlen, space, nspace);
				newlen += nspace;
				memcpy(newstr + newlen, str + start_curLen, current_lineLen);
				newlen += current_lineLen;
				newstr[newlen++] = '\n';
			}

			start_curLen += current_lineLen + 1L;
			current_lineLen = 0;

		}
	}
	
	if (isfirstline) {
		memcpy(newstr + newlen, str + start_curLen, current_lineLen);
	} else {
		memcpy(newstr + newlen, space, nspace);
		newlen += nspace;
		memcpy(newstr + newlen, str + start_curLen, current_lineLen);
	}
	newlen += current_lineLen;
	if (newstr[newlen] != '\n') {
		newstr[newlen++] = '\n';
	}
	newstr[newlen++] = 0;
	strcpy(str, newstr);
	return str;
}


void ToUpper(char* s) { for (int i = 0; s[i] != '\0'; i++) 	s[i] = s[i] >= 'a' && s[i] <= 'z' ? s[i] - 32 : s[i]; }

I32 strcicmp(char const * _restrict a, char const * _restrict b) {
	for (;; a++, b++) {
		I32 d = ((*a) | (U08)32) - ((*b) | (U08)32);
		if (d != 0 || !*a)	
			return d;
	}
}

I32 strcicmp_nfirst(char const* _restrict a, char const* _restrict b, int nfirst) {
	// only compare the strings up to the first 'nfirst' bytes
	if (nfirst == 0) {
		nfirst = (int) strlen(a) + 1;
	}

	int i = 0;
	for (;; a++, b++) {
		I32 d = ((*a) | (U08)32) - ((*b) | (U08)32);
		i++;
		if (d != 0 || !*a || i >= nfirst) {
			return d;
		}		
	}
}
I32 strmatch(char const* _restrict full, char const* _restrict part) {
	for (;; full++, part++) {
		I32 d = ((*full) | (U08)32) - ((*part) | (U08)32);
		if (d != 0 || !*part) {
			return (!*part)? 0: d;
		}
		
	}
}

F32 DeterminePeriod(F32PTR Y, I32 N)
{
	F32PTR TMP   = (F32PTR)malloc(sizeof(F32)*N * 6);

	F32PTR X    =		  TMP;        //the first 4*N elems are used for the 4 rows of x
	F32PTR Yfit =		  TMP+N*4;    //the fifth row is used for yIft
	U08PTR isNA = (U08PTR)(Yfit+N);   //the sixth row is used for isNA

	F32    delta    = 1.f / N;
	f32_fill_val(1.0, X, N);        //const term: 1,1,1,
	f32_seq(X + N,   0.0, delta, N);  //the linear term
	memcpy( X + 2*N, X + N, sizeof(F32) * N); f32_pow_vec_inplace(X + 2 * N, 2, N);//3rd-order
	memcpy( X + 3*N, X + N, sizeof(F32) * N); f32_pow_vec_inplace(X + 3 * N, 3, N);//4th-order


	memset(isNA, 0, sizeof(char)*N); 
	for ( I32 i = 0; i < N; i++) 	{
		if (Y[i] != Y[i]  ) 		{
			isNA[i]      = 1L;
			X[i]         = 0.f;
			X[N     + i] = 0.f;
			X[N+N   + i] = 0.f;
			X[N+N+N + i] = 0.f;
			Y[i] = 0.f;
		}
	}

	{
		F32 XtX_tmp[16];
		F32 B[4];		
		linear_regression(Y, X, N, N, 4L, B, Yfit, Y, XtX_tmp);
		//Up to now, Y is the residual: the detrended Y
	}
	
	//YFIT and X are not used after this point
	//But isNA is still used. It is safe to use the start of TMP as buffers of ans and INDEX

	/***********************************************/
	// The start of calculating auto-correlation function
	/***********************************************/
 
	I32    M     = (int)   ceil(N / 2);
	F32PTR ans   =		   TMP;       // ans is an array of lenght M	

	for (I32 i = 1; i <= M; i++){
		I32 nLength   = N - i;
		I32 start = i + 1;
		F32 sumXY = 0, sumXX = 0, sumYY = 0, sumX= 0, sumY = 0;
		I32 nValidNum = 0;
		for (I32 j = 1; j <= nLength; j++)		{
			I32 Ix = j - 1;
			I32 Iy = start + (j - 1) - 1;
			if ( (isNA[Ix] + isNA[Iy]) == 0)	{
				nValidNum++;
				F32 x = Y[Ix];
				F32 y = Y[Iy];
				sumX  += x;	
				sumY  += y;
				sumXY += x*y; 				
				sumXX += x*x;
				sumYY += y*y;
			}
		}
		F32 MX = sumX / (F32)nValidNum;
		F32 MY = sumY / (F32)nValidNum;
		ans[i - 1] = (sumXY/nValidNum - MX*MY) / sqrtf((sumXX / N - MX*MX)*(sumYY / N - MY*MY));
	}


	//isNA is reused as a tmp buff
	U08PTR isPeak = isNA;
	I32PTR INDEX  = (I32PTR)(TMP + M); // INDEX is an array of lenght M

	memset(isNA, 0, (size_t) M);
	I32  numPeaks = 0;
	for ( I32 i = 2; i <= (M - 1); i++)	{
		if (ans[(i)-1] > ans[(i - 1) - 1] && ans[(i)-1] > ans[(i + 1) - 1]) {
			isPeak[i - 1] = 1;
			INDEX[numPeaks++] = i;
		}
	}
	 
	//Up to this point, isPeak(TMP+5*N) and INDEX (TMP+M) are still used
 
	I32 period = -1L;	 
	if (numPeaks != 0) 	{
		for (I32 pID = 1; pID <= max(1, (int)floorf(numPeaks / 3.f)); pID++)	{

			period = INDEX[pID - 1];			
			I32  NumOfPeriod = (int)floorf((F32)(M - 1) / period);
			I32  goodTimes   = 0;
			for (I32 i = 1; i <= NumOfPeriod; i++)	{			 
				I32 segLength = period * i;
				goodTimes += isPeak[segLength - 1] || isPeak[segLength + 1 - 1] || isPeak[segLength-1 - 1];
			}

			if (goodTimes >= min(3, NumOfPeriod))			{
				break;
			}
			period = -1L;
		}

	}
	 
	free(TMP);
	return (F32)period;
}

 static F32 confidenceInterval(F32PTR half, I32 n, char leftOrRight)
{
	//r_cblas_scopy(w + 1, prob + cpt[i] - w, 1, half, 1);
	
	F32 sum      = f32_sum(half, n);
	//r_ippsSum_32f( half, n,   &sum, ippAlgHintAccurate);
	//r_ippsMulC_32f_I(1.f/sum, half, n);
		
	if (leftOrRight == 'R') {		
		F32 cumSum = 0;		
		I32 j      = 0;
		for (; j < n; j++)	{
			cumSum += half[j];
			if ( cumSum/sum >= 0.95f) break;
		}
		F32 J  = j + 1.f;
		return J - (cumSum - 0.95f*sum) / half[j];
		 
	} 

	if (leftOrRight == 'L')
	{		
		F32 cumSum = 0;		
		I32 j      =n-1;
		for (; j >=0; j--) 	{
			cumSum += half[j];
			if (cumSum / sum >= 0.95f) break;
		}
		F32 J = (F32)(n-j);
		F32 delta = J - (cumSum - 0.95f*sum) / half[j];
		return delta;
	}

	return -999;
	 	
}
static I32 find_changepoint_v0(F32PTR prob, F32PTR mem, F32 threshold, I32PTR cpt, F32PTR cptCI, I32 N, I32 minSepDist, I32 maxCptNumber)
{
	if (maxCptNumber == 0)	{return maxCptNumber;}

	//half the inter-changepoint distacne 	
	I32 w = (I32) round((minSepDist - 1) / 2);
	w = w >= 0 ? w : 0;
	I32 w2 = w * 2 + 1;
	//Set zeros to the mem array--used to hold the sum
	r_ippsSet_32f(0, mem, N);			
	
	I32PTR cpfromSum_Pos = (I32PTR)mem + N;
	F32PTR cpfromSum_Val = (F32PTR)mem + N * 2;

	I32PTR cpfromProb_Pos = (I32PTR) mem + N * 3;
	F32PTR cpfromProb_Val = (F32PTR)mem + N * 4;

	//A fast way to apply a w-sum filter to the prob arrary
	for (I32 i = -w; i <= w; i++)
	{
		I32 len = i > 0 ? i : -i;
		I32 startIdx_mem  = i <= 0 ? 0 : i;
		I32 startIdx_prob = i <= 0 ? -i : 0;
		r_ippsAdd_32f_I(prob + startIdx_prob, mem + startIdx_mem, N - len);
	}
	// Now mem holds the sum-filterd probability profile: the value at a point is 
	// the sume of the prob values within a window

	I32  UPPERIDX  = N - w;
	I32  numCpt = 0;
	for (I32 i = w; i < UPPERIDX; i++)
	{
		if (mem[i] < threshold) continue;
		
		Bool isLargeThanNearestNeighor		= (mem[i] >= mem[i - 1]) && (mem[i] >= mem[i + 1]);
		Bool isLargeThanNearestTwoNeighors	= (mem[i] * 4.0) > (mem[i + 1] + mem[i + 2] + mem[i - 1] + mem[i - 2]);
		if (!(isLargeThanNearestNeighor && isLargeThanNearestTwoNeighors)) continue;

		// If the current point is a local maximum in the filtered prob curve
		// , then find the maximum in the prob curve within a (w2+1) window at the current point
		I32		upperIdx_1	= i + w;
		I32		maxIdx		= -999;
		F32		maxVal		= -999;
		for (I32 j = i - w; j <= upperIdx_1; j++)
		{
			if ((prob[j] > prob[j - 1] && prob[j] >= prob[j + 1]) || (prob[j] >= prob[j - 1] && prob[j] > prob[j + 1]))
			{
				if (prob[j] > maxVal) 	maxIdx = j, maxVal = prob[j];
			
			}			
		}		
		if ( maxVal < 0.f	)	continue;
			

		//If the current identified peak is within a distance less (2*w+1)
		I32 diff_btw_twoNeighbors	=	maxIdx-cpfromProb_Pos[numCpt - 1]; 
		if ((numCpt == 0) || diff_btw_twoNeighbors >= w2 || diff_btw_twoNeighbors <= -w2)
		{
			cpfromSum_Pos[numCpt] = i;
			cpfromSum_Val[numCpt] = mem[i];
			cpfromProb_Pos[numCpt] = maxIdx;
			cpfromProb_Val[numCpt] = maxVal;
			numCpt++;
			continue;
		}
		else
		{
			if (maxVal >= cpfromProb_Val[numCpt - 1])
			{
				cpfromSum_Pos[numCpt - 1] = i;
				cpfromSum_Val[numCpt - 1] = mem[i];
				cpfromProb_Pos[numCpt - 1] = maxIdx;
				cpfromProb_Val[numCpt - 1] = maxVal;
				continue;
			}
		}		
		 
	}//for (I32 i = w; i < N - w; i++)
	if (numCpt == 0) { return numCpt; }
	 

	f32_QuickSortD(cpfromProb_Val, cpfromProb_Pos, 0, numCpt - 1);

	numCpt = min(numCpt, maxCptNumber);
	r_cblas_scopy(numCpt, (F32PTR)cpfromProb_Pos, 1, (F32PTR) cpt, 1);

	I32PTR INDEX = (I32PTR) mem;
	F32PTR CPT_float = mem + N;
	for (I32 i = 0; i < numCpt; i++)
	{
		*INDEX++ = i;
		*CPT_float++ = (F32)cpt[i];

	}
	INDEX		= INDEX - numCpt;
	CPT_float	= CPT_float - numCpt;	
	f32_QuickSortA(CPT_float, INDEX, 0, numCpt - 1);

	//Compute confidence intervals for indentified changepoints
	for (I32 i = 0; i < numCpt; i++)
	{
		cptCI[i] = -9999.f;
		cptCI[numCpt+i] = -9999.f;
	}

	F32 delta;
	delta = confidenceInterval(prob, ((I32) CPT_float[0]-0+1), 'L');
	cptCI[0]						 =   delta;
	delta = confidenceInterval(prob + (I32)CPT_float[numCpt - 1], (N - (I32)CPT_float[numCpt - 1] + 1), 'R');
	cptCI[numCpt + numCpt - 1] =   delta;
	if (numCpt == 1) {
		cptCI[0] = CPT_float[0] - cptCI[0];
		cptCI[1] = CPT_float[0] + cptCI[1];
		return numCpt; 
	}

	for (I32 i = 0; i < (numCpt-1); i++)
	{ 
		F32 del1,del2,del;
		del1 = cptCI[numCpt + i] > 0 ? cptCI[numCpt + i] : cptCI[i];
		del2 = cptCI[i + 1] > 0 ? cptCI[i + 1] : ((cptCI[numCpt + i + 1] > 0) ? cptCI[numCpt + i + 1] : -9999.f);

    	del = CPT_float[i + 1] - CPT_float[i];
		if (del2 <= 0)
		{
				del1 =  del1 * 2.f;
				del = (del1 > del) ? del/2 : del1;
		}else
		{
			del = del*del1 / (del1 + del2);
		}
			
		delta = confidenceInterval(prob + (I32)CPT_float[i], (I32) ceil(del), 'R');
		cptCI[numCpt + i] = delta;


		del = CPT_float[i + 1] - CPT_float[i];
		if (del2 <= 0)
		{
			delta = del - delta * 2;
			del = delta <= 0 ? del / 2 : delta;
		}
		else
		{
		
			del2 = del2 * 2.f;
			del = (del2 >= del) ? del/2 : del2;
		}
		I32 len = (I32)ceil(del);
		delta = confidenceInterval(prob + (I32)CPT_float[i+1]-(len-1), len, 'L');
		cptCI[i+1] = delta;

	}

	F32PTR temp = mem + 2 * N;
	r_cblas_scopy(2*numCpt, cptCI, 1, temp, 1);

	for (I32 i = 0; i < numCpt; i++)
	{
		I32 idx  ;
		idx = INDEX[i];
		cptCI[idx] = CPT_float[i] - temp[i];
		cptCI[numCpt + idx] = CPT_float[i] + temp[numCpt+i];
	}

	return numCpt;
}
 I32 FindChangepointv1(F32PTR prob, F32PTR mem, F32 threshold, I32PTR cpt, F32PTR cptCI, I32 N, I32 minSepDist, I32 maxCptNumber)
{
	if (maxCptNumber == 0)	{ return maxCptNumber; }
	// Set zeros to the mem array--used to hold the sum-filtered prob
	r_ippsSet_32f(0, mem, N);

	I32PTR cpfromSumP_Pos  =  (I32PTR) mem + N;
	F32PTR cpfromSumP_Val  =  (F32PTR) mem + N * 2;
	I32PTR cpfromProb_Pos =  (I32PTR) mem + N * 3;
	F32PTR cpfromProb_Val =  (F32PTR) mem + N * 4;

	// Half the inter-changepoint distacne
	I32 w0 = minSepDist / 2;   // leftside half
	I32 w1 = minSepDist - w0;  // rightside half
	/*
	// A fast way to apply a w-sum filter to the prob arrary
	for (I32 i = -w1; i <= w0; i++)
	{
		I32 len = i > 0 ? i : -i;
		I32 startIdx_mem = i <= 0 ? 0 : i;
		I32 startIdx_prob = i <= 0 ? -i : 0;
		r_ippsAdd_32f_I(prob + startIdx_prob, mem + startIdx_mem, N - len);
	}	
	// Now mem holds the sum-filterd probability profile: the value at a point 
	// is the sume of the prob values within a window [i-w0, i+w1]
	*/
   
	// Now mem holds the sum-filterd probability profile: the value at a point 
	// is the sume of the prob values within a window [i-w0,i+w1]
	//In this funciton, the two parts of the window are computed differently from the above commented routine
	//	I32 w0 = minSepDist / 2;   // leftside half
	// I32 w1 = (minSepDist - w0)-1;  // rightside half
	f32_sumfilter(prob, mem, N, minSepDist);
 

	I32  LOWERIDX  = (minSepDist + 1);
	I32  UPPERIDX  = N - (minSepDist + 1);
	I32  numCpt    = 0;
	for (I32 i = LOWERIDX; i < UPPERIDX; i++)
	{
		if (mem[i] < threshold) continue;

		Bool isLargeThanNearestNeighor     = (mem[i] >= mem[i - 1]) && (mem[i] >= mem[i + 1]);
		Bool isLargeThanNearestTwoNeighors = (mem[i] * 4.0) > (mem[i + 1] + mem[i + 2] + mem[i - 1] + mem[i - 2]);
		if ( isLargeThanNearestNeighor==0 || isLargeThanNearestTwoNeighors==0 ) continue;

		// If the current point is a local maximum in the filtered prob curve
		// , then find the maximum in the PROB curve within a (w2+1) window at the current point
		I32		UPPERIDX_1 = i + w1;
		I32		maxIdx = -9999;
		F32		maxVal = -9999.f;
		for (I32 j = i - w0; j <= UPPERIDX_1; j++) 	{
			if ((prob[j] > prob[j - 1] && prob[j] >= prob[j + 1]) || (prob[j] >= prob[j - 1] && prob[j] > prob[j + 1]))			{
				if (prob[j] > maxVal) {
					maxIdx = j;
					maxVal = prob[j];
				}				
			}
		}
		if (maxVal <= 0.f)	continue;


		//If the current identified peak is within a distance less (2*w+1)
		I32 dist_to_prevCpt = maxIdx - cpfromProb_Pos[numCpt - 1];
		if ((numCpt == 0) || dist_to_prevCpt > minSepDist || dist_to_prevCpt < -minSepDist)	{
			cpfromSumP_Pos[numCpt]  = i;
			cpfromSumP_Val[numCpt]  = mem[i];
			cpfromProb_Pos[numCpt] = maxIdx;
			cpfromProb_Val[numCpt] = maxVal;
			numCpt++;
			continue;
		} else	{  
			// If the current peak has a larer magnitude than the previously identified cp plus 
			// it is within a distance less than minSepDist , then replace the old changepoint
			// with the new peak
			/*
			if (maxVal >= cpfromProb_Val[numCpt - 1]){
				cpfromSum_Pos[numCpt - 1]  = i;
				cpfromSum_Val[numCpt - 1]  = mem[i];
				cpfromProb_Pos[numCpt - 1] = maxIdx;
				cpfromProb_Val[numCpt - 1] = maxVal;
				continue;
			}
			*/
			if (mem[i] >= cpfromSumP_Val[numCpt - 1]) {
				cpfromSumP_Pos[numCpt - 1] = i;
				cpfromSumP_Val[numCpt - 1] = mem[i];
				cpfromProb_Pos[numCpt - 1] = maxIdx;
				cpfromProb_Val[numCpt - 1] = maxVal;
				continue;
			}
		}

	}//for (I32 i = w; i < N - w; i++)

	if (numCpt == 0) { return numCpt; }
	
	f32_QuickSortD(cpfromSumP_Val, cpfromProb_Pos, 0, numCpt - 1);
	numCpt  = min(numCpt, maxCptNumber);

	f32_copy( (F32PTR)cpfromProb_Pos, (F32PTR)cpt, numCpt);	
	// At this point, cpfromSum_val has the prob values sorted from largest to least
	// cpfromProb_Pos and cpt (integers) are the same, containing the cpt locations

	// Now cpfromSum_Val(mem+2*N) should not be touched and will be used later
	I32PTR INDEX_timeToProbAmp = (I32 *)mem ;
	F32PTR cpt_f32             = (F32 *)mem+N;	
	for (I32 i = 0; i < numCpt; i++) {		
		cpt_f32[i]             = (F32)cpt[i];
		INDEX_timeToProbAmp[i] = i;
	}
	f32_QuickSortA(cpt_f32, INDEX_timeToProbAmp, 0, numCpt - 1);
	// Now, cpt_f32 contains the cpts ranked in the order of time
	// INDEX saves the indices mapping to the list ranked in magnitdue

	//Compute confidence intervals for indentified changepoints
	// Fill the ctpCI with negative values
	f32_fill_val(-9999.f, cptCI, 2*numCpt);
		
	F32PTR tmpSeg  = (F32*) mem + 3 * N;
	I32PTR nullSeg = (I32*) mem + 4 * N;  
	for (I32 i = 0; i < numCpt; i++)
	{
		I32 startIdx, endIdx,len;

		endIdx	 = (I32) cpt_f32[i];
		startIdx = i==0 ? 0 : (I32) cpt_f32[i-1];
		startIdx = (startIdx + endIdx) / 2;
		len = endIdx-startIdx + 1;		

		f32_copy(prob+startIdx, tmpSeg, len);		
		f32_QuickSortA(tmpSeg, nullSeg, 0, len - 1); // nullSeg is just provided as an input but the result is not used
		cptCI[i] = confidenceInterval(tmpSeg, len, 'L');

    	//-------------------------------------------------------------------
		//-------------------------------------------------------------------

		startIdx = (I32)cpt_f32[i];
		endIdx   = i == (numCpt - 1) ? (N - 1) : (I32)cpt_f32[i + 1];
		endIdx   = (startIdx + endIdx) / 2;
	    len      = endIdx - startIdx + 1;

		f32_copy(prob + startIdx, tmpSeg, len);		
		f32_QuickSortD(tmpSeg, nullSeg, 0, len - 1); // nullSeg is just provided as an input but the result is not used
		cptCI[numCpt + i] = confidenceInterval(tmpSeg, len, 'R');
 	}
	
	//Now, cpt contains knots sorted in magnitude
	// cptCI contains CI for knots sorted in time

	F32PTR cptCI_backup = mem + 3*N;
	f32_copy(cptCI, cptCI_backup, 2 * numCpt);
 
	F32PTR cpt_summedProb = mem;
	/*
	mem   : Index/cpt_summedProb
	mem+N : CPT_float
	mem+2*N:cpfromSum_Val
	mem+3*N:cpt_CI_backup
	*/
	for (I32 i = 0; i < numCpt; i++)	{
		//Duel use of mem: one for INDEX and another for cpt_summrProb
		//The index is feteched first and the same elment is overwritten with the summProb
		I32 idx             = INDEX_timeToProbAmp[i];
		cptCI[idx]          = cpt_f32[i] - cptCI_backup[i];
		cptCI[numCpt + idx] = cpt_f32[i] + cptCI_backup[numCpt + i];

		cpt_summedProb[i] = cpfromSumP_Val[i]>1 ? 1.f : cpfromSumP_Val[i];
	}

	return numCpt;
}

static I32 FindChangepoint_minseg_is_1(F32PTR prob, F32PTR mem, F32 threshold, I32PTR cpt, F32PTR cptCI, I32 N, I32 minSepDist, I32 maxCptNumber) {

	// minSepDist may be zero.
	
	F32PTR sump           = mem;
	I32PTR cpfromSumP_Pos = (I32PTR)mem + N;
	F32PTR cpfromSumP_Val = (F32PTR)mem + N * 2;
	I32PTR cpfromProb_Pos = (I32PTR)mem + N * 3;
	F32PTR cpfromProb_Val = (F32PTR)mem + N * 4;

	// Set zeros to the mem array--used to hold the sum-filtered prob
	r_ippsSet_32f(0, sump, N);        

    // Now mem holds the sum-filterd probability profile: the value at a point  is the sume 
	// of the prob values within a window [i-w0,i+w1].In this funciton, the two parts of 
	// the window are computed differently from the above commented routine [i-w0, i+w1]
	I32 minSegLen = minSepDist + 1;
	I32 w0 = minSegLen / 2;       // leftside half
	I32 w1 = (minSegLen - w0)-1;  // rightside half	  and we have  w0 >= w1		
	f32_sumfilter(prob, sump, N, minSegLen);

	//Special treatment when minSepDist=0; this is for the outlier component only
	I32  LOWERIDX = 0;       // LOWERIDX=(minSepDist + 1);  
	I32  UPPERIDX = N - 1;   /// N - (minSepDist + 1);
	I32  numCpt = 0;

	for (I32 i = LOWERIDX; i < UPPERIDX; i++) 	{
		if (sump[i] < threshold) continue;
		cpfromSumP_Pos[numCpt] = i;
		cpfromSumP_Val[numCpt] = sump[i];
		numCpt++;	 
	}

	// The cpts identified above are based on peak in the summed probability, and now for each of them
	// find the best peak in the original prob.
	for (I32 i = 0; i < numCpt; i++) {
		I32     cpt         = cpfromSumP_Pos[i];		
		I32     LOWERIDX    = cpt-w0;
		I32		UPPERIDX    = cpt+w1;

		cpfromProb_Pos[i] = cpt;
		cpfromProb_Val[i] = prob[cpt];

	}
	//cpfromProb_Pos may contian changepoinss that are within a distance of minSepDist from each other
 
	if (numCpt == 0) { return numCpt; }

	f32_QuickSortD(cpfromSumP_Val, cpfromProb_Pos, 0, numCpt - 1);
	numCpt = min(numCpt, maxCptNumber);

	f32_copy((F32PTR)cpfromProb_Pos, (F32PTR)cpt, numCpt);
	// At this point, cpfromSumP_val has the prob values sorted from largest to least
	// cpfromProb_Pos and cpt (integers) are the same, containing the cpt locations

	// Now cpfromSumP_Val(mem+2*N) should not be touched and will be used later
	I32PTR INDEX_timeToProbAmp = (I32*)mem;
	F32PTR cpt_f32              = (F32*)mem + N;
	for (I32 i = 0; i < numCpt; i++) {
		cpt_f32[i]             = (F32)cpt[i];
		INDEX_timeToProbAmp[i] = i;
	}
	f32_QuickSortA(cpt_f32, INDEX_timeToProbAmp, 0, numCpt - 1);
	// Now, cpt_f32 contains the cpts ranked in the order of time
	// INDEX saves the indices mapping to the list ranked in magnitdue

	//Compute confidence intervals for indentified changepoints
	// Fill the ctpCI with negative values
	

	//Now, cpt contains knots sorted in magnitude
	// cptCI contains CI for knots sorted in time

	F32PTR cpt_summedProb = mem;
	/*
	mem   : Index/cpt_summedProb
	mem+N : CPT_float
	mem+2*N:cpfromSum_Val
	mem+3*N:cpt_CI_backup
	*/
	for (I32 i = 0; i < numCpt; i++) {
		//Duel use of mem: one for INDEX and another for cpt_summrProb
		//The index is feteched first and the same elment is overwritten with the summProb
		I32 idx = INDEX_timeToProbAmp[i];
		cptCI[idx]          = cpt_f32[i] - 0;
		cptCI[numCpt + idx] = cpt_f32[i] + 0;

		cpt_summedProb[i] = cpfromSumP_Val[i] > 1 ? 1.f : cpfromSumP_Val[i];
	}

	return numCpt;
}

 
I32 FindChangepoint(F32PTR prob, F32PTR mem, F32 threshold, I32PTR cpt, F32PTR cptCI, I32 N, I32 minSepDist, I32 maxCptNumber)
{
	// minSepDist may be zero.
	// The seg pattern is  "c"     when minseg=1
	// The seg pattern is  "xc"    when minseg=2
	// The seg pattern is  "xcx"   when minseg=3
	// The seg pattern is  "xxcx"   when minseg=4
	// The seg pattern is  "xxcxx"  when minseg=5

	if (maxCptNumber == 0) { return maxCptNumber; }

	if (minSepDist == 0) {
	// this is only for the outlier component
		return FindChangepoint_minseg_is_1(prob, mem, threshold, cpt, cptCI, N, minSepDist, maxCptNumber);
	}

	F32PTR sump           = mem;
	I32PTR cpfromSumP_Pos = (I32PTR)mem + N;
	F32PTR cpfromSumP_Val = (F32PTR)mem + N * 2;
	I32PTR cpfromProb_Pos = (I32PTR)mem + N * 3;
	F32PTR cpfromProb_Val = (F32PTR)mem + N * 4;

	/* // A fast way to apply a w-sum filter to the prob arrary
	for (I32 i = -w1; i <= w0; i++)	{
		I32 len          = i > 0 ? i : -i;
		I32 startIdx_mem = i <= 0 ? 0 : i;
		I32 startIdx_prob = i <= 0 ? -i : 0;
		r_ippsAdd_32f_I(prob + startIdx_prob, mem + startIdx_mem, N - len);
	} // Now mem holds the sum-filterd probability profile: the value at a point is the sume of the prob values within a window [i-w0, i+w1]
	*/

	// Set zeros to the mem array--used to hold the sum-filtered prob
	r_ippsSet_32f(0, sump, N);        

    // Now mem holds the sum-filterd probability profile: the value at a point  is the sume 
	// of the prob values within a window [i-w0,i+w1].In this funciton, the two parts of 
	// the window are computed differently from the above commented routine [i-w0, i+w1]
	I32 minSegLen = minSepDist + 1;
	I32 w0        = minSegLen / 2;       // leftside half
	I32 w1        = (minSegLen - w0)-1;  // rightside half	  and we have  w0 >= w1		
	f32_sumfilter(prob, sump, N, minSegLen);

	I32  LOWERIDX = (minSepDist + 1);
	I32  UPPERIDX = N - (minSepDist + 1);
	I32  numCpt = 0;
	for (I32 i = LOWERIDX; i < UPPERIDX; i++) 	{

		if (sump[i] < threshold) continue;

		// Imporant: stricltly larger than the keft, larger or equal to the right: oick up the rising side of the plateau peak
		Bool isLargeThanNearestNeighor     = (sump[i] > sump[i - 1] && sump[i] >= sump[i + 1]) || (sump[i] >= sump[i - 1] && sump[i] >  sump[i + 1]);
		Bool isLargeThanNearestTwoNeighors = (sump[i] * 4.0) >= (sump[i + 1] + sump[i + 2] + sump[i - 1] + sump[i - 2]);
		// if minSegLen = 1,2, 3, or 4, we only check the two neihgoring points, not the four neighors
		if (isLargeThanNearestNeighor == 0 || (minSegLen>4 && isLargeThanNearestTwoNeighors == 0)  ) continue;

		/*********************************************************************/
		// If the current point is a local maximum in the filtered prob curve
		/*********************************************************************/
		I32 dist_to_prevCpt = i - cpfromSumP_Pos[numCpt - 1]; //invalid if numCpt==0
		if ((numCpt == 0) || dist_to_prevCpt > minSepDist || dist_to_prevCpt < -minSepDist) {
			//If the current identified peak is NOT within a distance less (2*w+1)
			cpfromSumP_Pos[numCpt] = i;
			cpfromSumP_Val[numCpt] = sump[i];
			numCpt++;
			continue;
		} else {
			//if it is within a distance less than minSepDist PLUS the current peak has 
			// a larer magnitude than the previously identified cp then replace the old changepoint
			// with the new peak.			
			if (sump[i] > cpfromSumP_Val[numCpt - 1]) { // Possible buggy if using sump[i] ">=" cpfromSumP_Val[numCpt - 1]
				cpfromSumP_Pos[numCpt - 1] = i;
				cpfromSumP_Val[numCpt - 1] = sump[i];
				continue;
			}
		}

	}

	/************************************************************************************************/
	// The cpts identified above are based on peak in the summed probability, and now for each of them
	// find the best peak in the original prob.
	/************************************************************************************************/
	if (minSepDist <= 4) {
		for (I32 i = 0; i < numCpt; i++) {
			I32     cpt = cpfromSumP_Pos[i];
			I32     LOWERIDX = cpt - w0;
			I32		UPPERIDX = cpt + w1;

			// Pre-assign the best cpt to the sumP-based cpt: this is needed when minSepdist=0 (i.e., outlier componet)
			// such that the following loop only runs once and no local min is found to update maxIdx and maxVAl
			I32		maxIdx = cpt;
			F32		maxVal = prob[cpt];
			for (I32 j = LOWERIDX; j <= UPPERIDX; j++) {
				F32 curProb = prob[j];
				if ((curProb > prob[j - 1] && curProb >= prob[j + 1]) || (curProb >= prob[j - 1] && curProb > prob[j + 1])) {
					if (curProb > maxVal) {
						maxIdx = j;
						maxVal = curProb;
					}
				}
			}
			cpfromProb_Pos[i] = maxIdx;
			cpfromProb_Val[i] = maxVal;

		}
	} else {
		// if minSeptDist >4, we will use a more complex rule to pinpoint the best peak bcz  the simple alg for the minsepdist,4 case
		// may not find the location with the highest amassed probability
		for (I32 i = 0; i < numCpt; i++) {
			I32     cpt = cpfromSumP_Pos[i];
			I32     LOWERIDX = cpt - w0;
			I32		UPPERIDX = cpt + w1;

			//First, find all the local peaks up to a number of 10
			I32 topTenPeaksLoc[10];
			F32 topTenPeaksPrb[10];
			I32 numTopTenCpt = 0;

			for (I32 j = LOWERIDX; j <= UPPERIDX; j++) {
				F32 curProb = prob[j];
				if ((curProb > prob[j - 1] && curProb >= prob[j + 1]) || (curProb >= prob[j - 1] && curProb > prob[j + 1])) {
					if (numTopTenCpt < 10) {
						topTenPeaksLoc[numTopTenCpt] = j;
						topTenPeaksPrb[numTopTenCpt] = curProb;
						numTopTenCpt++;
					}	else {
						// If the buf is full, find the worst peak and replace if the curprob is higher.
						F32 minProbInList    = 100.f;
						int minProbListIndex = 0;
						for (int k = 0; k < 10; k++) {
							if (minProbInList > topTenPeaksPrb[k]) {
								minProbListIndex = k;
								minProbInList    = topTenPeaksPrb[k];
							}
						}

						if (curProb > minProbInList) {
							topTenPeaksLoc[minProbListIndex] = j;
							topTenPeaksPrb[minProbListIndex] = curProb;
						}
					}					
				} //(curProb > prob[j - 1] && curProb >= prob[j + 1]) || (curProb >= prob[j - 1] && curProb > prob[j + 1]))
			} //for (I32 j = LOWERIDX; j <= UPPERIDX; j++) 

			if (numTopTenCpt == 0) {
				// If no local peask were found (this is unlikely the case),just use the best cpt to the sumP-based cpt: 
				cpfromProb_Pos[i] = cpt;
				cpfromProb_Val[i] = prob[cpt];
			}
			else if (numTopTenCpt == 1) {
				cpfromProb_Pos[i] = topTenPeaksLoc[0];
				cpfromProb_Val[i] = topTenPeaksPrb[0];
			}
			else {
				// there are at least 2 cpts in the top list
				f32_QuickSortD(topTenPeaksPrb, topTenPeaksLoc, 0, numTopTenCpt - 1);
				if (topTenPeaksPrb[0] > 1.5 * topTenPeaksPrb[1]) {
					// if the best peak is far beter than the the seond best peak
					cpfromProb_Pos[i] = topTenPeaksLoc[0];
					cpfromProb_Val[i] = topTenPeaksPrb[0];
				}
				else {
					int bestK = 0;
					F32 bestProb_3Neighors = -9999;
					for (int k = 0; k < numTopTenCpt; k++) {
						int cpLoc = topTenPeaksLoc[k];
						F32 cpProb_3Neighors = prob[cpLoc] + prob[cpLoc + 1] + prob[cpLoc - 1];
						if (cpProb_3Neighors > bestProb_3Neighors) {
							bestProb_3Neighors = cpProb_3Neighors;
							bestK              = k;
						}
					}
					cpfromProb_Pos[i] = topTenPeaksLoc[bestK];
					cpfromProb_Val[i] = topTenPeaksPrb[bestK];
				}
			} //if (numTopTenCpt == 0) 

		}//for (I32 i = 0; i < numCpt; i++)
	}

	//cpfromProb_Pos may contian changepoinss that are within a distance of minSepDist from each other
 


	if (numCpt == 0) { return numCpt; }

	f32_QuickSortD(cpfromSumP_Val, cpfromProb_Pos, 0, numCpt - 1);
	numCpt = min(numCpt, maxCptNumber);

	f32_copy((F32PTR)cpfromProb_Pos, (F32PTR)cpt, numCpt);
	// At this point, cpfromSumP_val has the prob values sorted from largest to least
	// cpfromProb_Pos and cpt (integers) are the same, containing the cpt locations

	// Now cpfromSumP_Val(mem+2*N) should not be touched and will be used later
	I32PTR INDEX_timeToProbAmp = (I32*)mem;
	F32PTR cpt_f32              = (F32*)mem + N;
	for (I32 i = 0; i < numCpt; i++) {
		cpt_f32[i]             = (F32)cpt[i];
		INDEX_timeToProbAmp[i] = i;
	}
	f32_QuickSortA(cpt_f32, INDEX_timeToProbAmp, 0, numCpt - 1);
	// Now, cpt_f32 contains the cpts ranked in the order of time
	// INDEX saves the indices mapping to the list ranked in magnitdue

	//Compute confidence intervals for indentified changepoints
	// Fill the ctpCI with negative values
	f32_fill_val(-9999.f, cptCI, 2 * numCpt);

	F32PTR tmpSeg = (F32*)mem + 3 * N;
	I32PTR nullSeg = (I32*)mem + 4 * N;
	for (I32 i = 0; i < numCpt; i++) 	{

		I32 startIdx, endIdx, len;

		//-------------------------------------------------------------------
		endIdx = (I32)cpt_f32[i];
		startIdx = i == 0 ? 0 : (I32)cpt_f32[i - 1];
		startIdx = (startIdx + endIdx) / 2;
		len = endIdx - startIdx + 1;

		f32_copy(prob + startIdx, tmpSeg, len);
		f32_QuickSortA(tmpSeg, nullSeg, 0, len - 1); // nullSeg is just provided as an input but the result is not used
		cptCI[i] = confidenceInterval(tmpSeg, len, 'L');		
		//-------------------------------------------------------------------

		//-------------------------------------------------------------------
		startIdx = (I32)cpt_f32[i];
		endIdx = i == (numCpt - 1) ? (N - 1) : (I32)cpt_f32[i + 1];
		endIdx = (startIdx + endIdx) / 2;
		len = endIdx - startIdx + 1;

		f32_copy(prob + startIdx, tmpSeg, len);
		f32_QuickSortD(tmpSeg, nullSeg, 0, len - 1); // nullSeg is just provided as an input but the result is not used
		cptCI[numCpt + i] = confidenceInterval(tmpSeg, len, 'R');
		//-------------------------------------------------------------------
	}

	//Now, cpt contains knots sorted in magnitude
	// cptCI contains CI for knots sorted in time

	F32PTR cptCI_backup = mem + 3 * N;
	f32_copy(cptCI, cptCI_backup, 2 * numCpt);

	F32PTR cpt_summedProb = mem;
	/*
	mem   : Index/cpt_summedProb
	mem+N : CPT_float
	mem+2*N:cpfromSum_Val
	mem+3*N:cpt_CI_backup
	*/
	for (I32 i = 0; i < numCpt; i++) {
		//Duel use of mem: one for INDEX and another for cpt_summrProb
		//The index is feteched first and the same elment is overwritten with the summProb
		I32 idx = INDEX_timeToProbAmp[i];
		cptCI[idx]          = cpt_f32[i] - cptCI_backup[i];
		cptCI[numCpt + idx] = cpt_f32[i] + cptCI_backup[numCpt + i];

		cpt_summedProb[i] = cpfromSumP_Val[i] > 1 ? 1.f : cpfromSumP_Val[i];
	}

	return numCpt;
}

I32 FindChangepoint_LeftRightMargins(F32PTR prob, F32PTR mem, F32 threshold, I32PTR cpt, F32PTR cptCI, I32 N, I32 minSepDist, I32 maxCptNumber, I32 leftMargin, I32 rightMargin)
{
	if (maxCptNumber == 0) { return maxCptNumber; }

	I32 minSegLen = minSepDist + 1;

	// minSepDist may be zero.
	// The seg pattern is  "c"      when minseg=1 (minSeptDist=0)
	// The seg pattern is  "xc"     when minseg=2
	// The seg pattern is  "xcx"    when minseg=3
	// The seg pattern is  "xxcx"   when minseg=4
	// The seg pattern is  "xxcxx"  when minseg=5
		
	// this is only for the outlier component
	if (minSepDist  == 0) { return FindChangepoint_minseg_is_1(prob, mem, threshold, cpt, cptCI, N, minSepDist, maxCptNumber);}

	F32PTR sump           = mem;
	I32PTR cpfromSumP_Pos = (I32PTR)mem + N;
	F32PTR cpfromSumP_Val = (F32PTR)mem + N * 2;
	I32PTR cpfromProb_Pos = (I32PTR)mem + N * 3;
	F32PTR cpfromProb_Val = (F32PTR)mem + N * 4;

	/* // A fast way to apply a w-sum filter to the prob arrary
	for (I32 i = -w1; i <= w0; i++)	{
		I32 len           = i > 0 ? i : -i;
		I32 startIdx_mem  = i <= 0 ? 0 : i,  startIdx_prob = i <= 0 ? -i : 0;
		r_ippsAdd_32f_I(prob + startIdx_prob, mem + startIdx_mem, N - len);
	} // Now mem holds the sum-filterd probability profile: the value at a point is the sume of the prob values within a window [i-w0, i+w1]
	*/
		 
	// Now mem holds the sum-filterd probability profile: the value at a point  is the sume 
	// of the prob values within a window [i-w0,i+w1].In this funciton, the two parts of 
	// the window are computed differently from the above commented routine [i-w0, i+w1]
	I32 w0        = minSegLen / 2;         // leftside half
	I32 w1        = (minSegLen - w0) - 1;  // rightside half	  and we have  w0 >= w1		
	f32_sumfilter(prob, sump, N, minSegLen);  //  hold the sum-filtered prob
	sump[0] = 0; // force the first elemtn to be zero bcz the first pos can never be a changepoint

   #define IsLocalPeak(j, p)      (  ( p[j]>p[j-1] && p[j]>=p[j+1] )  || ( p[j]>=p[j-1]&&p[j] > p[j+1]   ) )
   #define IsPeak_4Neigbor(j, p)  (  p[j]>=p[j-1] && p[j]>=p[j+1] && p[j]>=p[j-2] && p[j] >= p[j+2]   )
   
	I32  LOWERIDX = 0     + leftMargin + 1; //leftMargin>=0; lowerindx >=1
	I32  UPPERIDX = (N-1) - rightMargin;    //rightMargin>=0; upperix <= N-1
	I32  numCpt   = 0;
	for (I32 i = LOWERIDX; i <= UPPERIDX; i++) {

		if (sump[i] < threshold) continue;

		// i alsways has a lefthand neighobor, but may not have a righthand neighor
		Bool hasRightNeighor = (i < N-1);
		Bool hasTwoNeighors  = i > 1 && (i < N- 2);

		// Imporant: stricltly larger than the left, larger or equal to the right: pick up the rising side of the plateau peak
		Bool isLargeThanNearestNeighor     = IsLocalPeak(i, sump);
		Bool isLargeThanNearestTwoNeighors = IsPeak_4Neigbor(i, sump);
		Bool isLargeThanLeftNeighor        = (sump[i] > sump[i - 1]);  // used only when  rightMargin=0

		// if minSegLen = 1,2, 3, or 4, we only check the two neihgoring points, not the four neighors
		if ( ( hasRightNeighor  && isLargeThanNearestNeighor==0 )                      ||
			 ( minSegLen > 4    && hasTwoNeighors && isLargeThanNearestTwoNeighors==0) ||
			 ( !hasRightNeighor && isLargeThanLeftNeighor==0) )
			continue;

		/*********************************************************************/
		// If the current point is a local maximum in the filtered prob curve
		/*********************************************************************/
		I32 dist_to_prevCpt = i - cpfromSumP_Pos[numCpt - 1]; //invalid if numCpt==0
		if ( (numCpt == 0) || dist_to_prevCpt > minSepDist || dist_to_prevCpt < -minSepDist) {
			//If the current identified peak is NOT within a distance less (2*w+1)
			cpfromSumP_Pos[numCpt] = i;
			cpfromSumP_Val[numCpt] = sump[i];
			numCpt++;			
		} else {
			// if it is within a distance less than minSepDist PLUS the current peak has a larer magnitude 
			// than the previously identified cp then replace the old changepoint with the new peak.			
			if (sump[i] > cpfromSumP_Val[numCpt - 1]) { // Possible buggy if using sump[i] ">=" cpfromSumP_Val[numCpt - 1]
				cpfromSumP_Pos[numCpt - 1] = i;
				cpfromSumP_Val[numCpt - 1] = sump[i];				
			}
		}

	}


	/************************************************************************************************/
	// The cpts given above are peaks in the summed prob, and now for each, find the best peak in the original prob.
	/************************************************************************************************/
	if (minSepDist <= 4) {
		for (I32 i = 0; i < numCpt; i++) {
			I32     cpt      = cpfromSumP_Pos[i];
			I32     LOWERIDX = max(cpt - w0, 1  );  // the pos 0 is neverr a changepont
			I32		UPPERIDX = min(cpt + w1, N-1);  // N-1 may be a changepoint when rightMargin=0

			// Pre-assign the best cpt to the sumP-based cpt: this is needed when minSepdist=0 (i.e., outlier componet)
			// such that the following loop only runs once and no local min is found to update maxIdx and maxVAl
			I32		maxIdx = cpt;
			F32		maxVal = prob[cpt];
			for (I32 j = LOWERIDX; j <= UPPERIDX; j++) {		 
				if (prob[j] > maxVal) {
						maxIdx = j;
						maxVal = prob[j];
				}
			}
			cpfromProb_Pos[i] = maxIdx;
			cpfromProb_Val[i] = maxVal;
		}
	}
	else {
		// if minSeptDist >4, we will use a more complex rule to pinpoint the best peak bcz  the simple alg 
		// for the minsepdist,4 case may not find the location with the highest amassed probability
		for (I32 i = 0; i < numCpt; i++) {
			I32     cpt      = cpfromSumP_Pos[i];
			I32     LOWERIDX = max(cpt - w0, 1);     // the pos 0 is neverr a changepont
			I32		UPPERIDX = min(cpt + w1, N - 1); // N-1 may be a changepoint when rightMargin=0

			// (1) First, find all the local peaks up to a number of 25
			I32 topTenPeaksLoc[15];
			F32 topTenPeaksPrb[15];
			I32 numTopTenCpt  = 0;

			for (I32 j = LOWERIDX; j <= UPPERIDX; j++) {
				F32 curProb           = prob[j];
				Bool hasRightNeighbor = (j < N - 1);
				if (  (hasRightNeighbor&&IsLocalPeak(j, prob) ) ||  ( !hasRightNeighbor && curProb>prob[j-1])    )
				{
					if (numTopTenCpt < 15) {
						topTenPeaksLoc[numTopTenCpt] = j;
						topTenPeaksPrb[numTopTenCpt] = curProb;
						numTopTenCpt++;
					} else {
						// If the buf is full, find the worst peak and replace if the curprob is higher.
						F32 minProbInList    = 100.f;
						int minProbListIndex = 0;
						for (int k = 0; k < 15; k++) {
							if (minProbInList > topTenPeaksPrb[k]) {
								minProbListIndex = k;
								minProbInList = topTenPeaksPrb[k];
							}
						}

						if (curProb > minProbInList) {
							topTenPeaksLoc[minProbListIndex] = j;
							topTenPeaksPrb[minProbListIndex] = curProb;
						}
					}
				} //(curProb > prob[j - 1] && curProb >= prob[j + 1]) || (curProb >= prob[j - 1] && curProb > prob[j + 1]))
			} //for (I32 j = LOWERIDX; j <= UPPERIDX; j++) 

			// (2) Second, pick up the best peak
			if (numTopTenCpt == 0) {
				// If no local peask were found (this is unlikely the case),just use the best cpt to the sumP-based cpt: 
				cpfromProb_Pos[i] = cpt;
				cpfromProb_Val[i] = prob[cpt];
			} else if (numTopTenCpt == 1) {
				cpfromProb_Pos[i] = topTenPeaksLoc[0];
				cpfromProb_Val[i] = topTenPeaksPrb[0];
			} else {
				// there are at least 2 cpts in the top list
				f32_QuickSortD(topTenPeaksPrb, topTenPeaksLoc, 0, numTopTenCpt - 1);
				if (topTenPeaksPrb[0] > 1.5 * topTenPeaksPrb[1]) {
					// if the best peak is far beter than the the seond best peak
					cpfromProb_Pos[i] = topTenPeaksLoc[0];
					cpfromProb_Val[i] = topTenPeaksPrb[0];
				} else {
					// determine the best peak based on the summed prob over the ttwo neighboring points
					int bestK              = 0;
					F32 bestProb_3Neighors = -9999;
					for (int k = 0; k < numTopTenCpt; k++) {
						int cpLoc             = topTenPeaksLoc[k];
						F32 cpProb_3Neighors = cpLoc==N-1? (prob[cpLoc] + prob[cpLoc - 1]): (prob[cpLoc] + prob[cpLoc + 1] + prob[cpLoc - 1]);
						if (cpProb_3Neighors > bestProb_3Neighors) {
							bestProb_3Neighors = cpProb_3Neighors;
							bestK = k;
						}
					}
					cpfromProb_Pos[i] = topTenPeaksLoc[bestK];
					cpfromProb_Val[i] = topTenPeaksPrb[bestK];
				}
			} //if (numTopTenCpt == 0) 

		}//for (I32 i = 0; i < numCpt; i++)
	}

	/*************************************************************************************************/
	//cpfromProb_Pos may contian changepoinss that are within a distance of minSepDist from each other
	/*************************************************************************************************/


	if (numCpt == 0) { return numCpt; }

	f32_QuickSortD(cpfromSumP_Val, cpfromProb_Pos, 0, numCpt - 1);
	numCpt = min(numCpt, maxCptNumber);

	f32_copy((F32PTR)cpfromProb_Pos, (F32PTR)cpt, numCpt);
	// At this point, cpfromSumP_val has the prob values sorted from largest to least
	// cpfromProb_Pos and cpt (integers) are the same, containing the cpt locations

	/*************************************************************************************************/
	// Now cpfromSumP_Val(mem+2*N) should not be touched and will be used later
	/*************************************************************************************************/
	I32PTR INDEX_timeToProbAmp = (I32*)mem;
	F32PTR cpt_f32             = (F32*)mem + N;
	for (I32 i = 0; i < numCpt; i++) {
		cpt_f32[i]             = (F32)cpt[i];
		INDEX_timeToProbAmp[i] = i;
	}
	f32_QuickSortA(cpt_f32, INDEX_timeToProbAmp, 0, numCpt - 1);
	// Now, cpt_f32 contains the cpts ranked in the order of time
	// INDEX saves the indices mapping to the list ranked in magnitdue

	//Compute confidence intervals for indentified changepoints
	// Fill the ctpCI with negative values
	f32_fill_val(-9999.f, cptCI, 2 * numCpt);

	F32PTR tmpSeg = (F32*)mem + 3 * N;
	I32PTR nullSeg = (I32*)mem + 4 * N;
	for (I32 i = 0; i < numCpt; i++) {

		I32 startIdx, endIdx, len;

		// Get the CI lenght on the lefthand side
		endIdx = (I32)cpt_f32[i];
		startIdx = i == 0 ? 0 : (I32)cpt_f32[i - 1];
		startIdx = (startIdx + endIdx) / 2;
		len = endIdx - startIdx + 1;

		f32_copy(prob + startIdx, tmpSeg, len);
		f32_QuickSortA(tmpSeg, nullSeg, 0, len - 1); // nullSeg is just provided as an input but the result is not used
		cptCI[i] = confidenceInterval(tmpSeg, len, 'L');
		//-------------------------------------------------------------------

		// Get the CI lenght on the righthand side
		startIdx = (I32)cpt_f32[i];
		endIdx   = i == (numCpt - 1) ? (N - 1) : (I32)cpt_f32[i + 1];
		endIdx   = (startIdx + endIdx) / 2;
		len      = endIdx - startIdx + 1;

		f32_copy(prob + startIdx, tmpSeg, len);
		f32_QuickSortD(tmpSeg, nullSeg, 0, len - 1); // nullSeg is just provided as an input but the result is not used
		cptCI[numCpt + i] = confidenceInterval(tmpSeg, len, 'R');
		//-------------------------------------------------------------------
	}

	// Now, cpt contains knots sorted in magnitude
	// cptCI contains CI for knots sorted in time

	F32PTR cptCI_backup = mem + 3 * N;
	f32_copy(cptCI, cptCI_backup, 2 * numCpt);

	F32PTR cpt_summedProb = mem;
	/*  	mem   :  Index/cpt_summedProb
	        mem+N :  CPT_float
	        mem+2*N: cpfromSum_Val
	        mem+3*N: cpt_CI_backup 	*/

	for (I32 i = 0; i < numCpt; i++) {
		// Dual use of mem: one for INDEX and another for cpt_summrProb
		// The index is feteched first and the same elment is overwritten with the cpt_summedProb
		I32 idx             = INDEX_timeToProbAmp[i];
		cptCI[idx]          = cpt_f32[i] - cptCI_backup[i];
		cptCI[numCpt + idx] = cpt_f32[i] + cptCI_backup[numCpt + i];

		cpt_summedProb[i] = min(1.0,  cpfromSumP_Val[i]);
	}

	return numCpt;
}

/*
int NormalizeF32ArrayWithNaNOmitted(F32PTR Y, U32PTR rowsMissing, I32 N,F32PTR mean, F32PTR sd, F32PTR yty) {

	I32     nMissing = 0;	
	F64     sum      = 0;
	F64     sqsum    = 0;
	for (rI32 i = 0; i < N; i++) { 
		F64 y = *Y++;
		if (y != y) {
			rowsMissing[nMissing++] = i;
		} else {
			sum   += y;
			sqsum += y*y;
		}
	}
	Y -= N;

	I32 n    = N - nMissing;
	F64 MEAN64 = sum / n;
	F64 SD64   = (sqsum / n - MEAN64*MEAN64); //https://www.mun.ca/biology/scarr/Simplified_calculation_of_variance.html
	SD64	   = SD64 > 0 ? sqrt(SD64): 1.f;

	F32 MEAN32    = MEAN64;
	F32 INV_STD32 = 1.f / (F32)SD64;
	F32 YtY = 0;
	{
		
		I32   jOmit   = 0;
		
		for (I32 i = 0; i< N; ++i) {
			if (jOmit < nMissing && i == rowsMissing[jOmit]) 	{
				Y[i] = 0.f; // fill bad points with zeros
				jOmit++;
			} else {
				Y[i]  = (Y[i] - MEAN32)*INV_STD32;
				YtY   += Y[i] * Y[i];
			}
		} 

		mean[0] = MEAN32;
		sd[0]   = (F32) SD64;
		yty[0]  = YtY;

	}

	return nMissing;
}
*/

#if defined(WIN64_OS) || defined(WIN32_OS) 
	#include "float.h"
	#if defined(MSVC_COMPILER)
	void EnableFloatExcetion(void) {
		unsigned int _oldState;
		errno_t err = _controlfp_s(&_oldState, EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL | EM_INVALID, MCW_EM);
		//errno_t err = _controlfp_s(&_oldState, EM_OVERFLOW | EM_ZERODIVIDE | EM_DENORMAL | EM_INVALID, MCW_EM);
		//errno_t err = _controlfp_s(&_oldState, EM_OVERFLOW | EM_ZERODIVIDE | EM_INVALID, MCW_EM); 
	}
	#elif defined(GCC_COMPILER)
void EnableFloatExcetion(void) {
		unsigned int _oldState;
		errno_t err = _controlfp_s(&_oldState, _EM_OVERFLOW | _EM_UNDERFLOW | _EM_ZERODIVIDE | _EM_DENORMAL | _EM_INVALID, _MCW_EM);
		//errno_t err = _controlfp_s(&_oldState, EM_OVERFLOW | EM_ZERODIVIDE | EM_DENORMAL | EM_INVALID, MCW_EM);
		//errno_t err = _controlfp_s(&_oldState, EM_OVERFLOW | EM_ZERODIVIDE | EM_INVALID, MCW_EM); 
	}
	#endif
#else
	#include "fenv.h" 
void EnableFloatExcetion(void) {
     // https://lists.gnu.org/archive/html/bug-gsl/2009-01/msg00001.html 
     // for unknow reasons, feenableexcept is not in fenv.h for Rtools.
	 //https://stackoverflow.com/questions/33191530/how-to-solve-undefined-reference-to-functions-in-fenv-h-when-using-uclibc
	 //lm is need for the linker 

	//stackoverflow.com/questions/37819235/how-do-you-enable-floating-point-exceptions-for-clang-in-os-x
	//https://gitlab.ikp.kit.edu/AirShowerPhysics/corsika/-/issues/415
	#if defined(LINUX_OS) 
	   // Not availalbe in C23 for clang16
	   //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW); 
	#endif

	
    // MacOS: feenableexcept not defined on MacOS  //https://philbull.wordpress.com/2012/07/27/portability-issues-with-floating-point-exception-handling/
	// Rtools: feenableexcept not defined for Rtools
	// Linux:  seems to be only defined for C99    https://linux.die.net/man/3/feenableexcept
	// Check C99: https://stackoverflow.com/questions/19674401/is-there-a-preprocessor-macro-to-detect-c99-across-platforms
}
#endif


/*
static INLINE void zeroOut_Xmars_zero(F32PTR Xt_mars, F32PTR Xt_zeroBackup, U32PTR rowsMissing,
	U32 nMissing, I32 N, I32 Npad, int k)
{
	//https: //stackoverflow.com/questions/3174850/what-is-the-correct-type-for-array-indexes-in-c#
	//register  ptrdiff_t tmpidx;

 
	for (j = 1; j <= k; j++) 	{
		for (i = 1; i <= nMissing; i++)	{
			I64 tmpidx = rowsMissing[i - 1] - 1;
			*Xt_zeroBackup++ = Xt_mars[tmpidx];
			Xt_mars[tmpidx] = 0.f;
			}
			Xt_mars = Xt_mars + Npad;
	}
 

	for (I32 j = k; j > 0; j--) {
		for (I32 i = nMissing; i > 0; i--) {
			I64 tmpidx = (*rowsMissing++);
			tmpidx--;
			*Xt_zeroBackup++ = Xt_mars[tmpidx];
			Xt_mars[tmpidx] = 0.f;
		}
		rowsMissing -= nMissing;
		Xt_mars += Npad;
	}

}
*/
 
 

#include "abc_cpu.h"
#include "abc_ide_util.h"
void SetupRoutines_AutoByCPU(Bool quite) {

	static int IS_CPU_INSTRUCTON_CHECKED = 0;
	if (IS_CPU_INSTRUCTON_CHECKED) {
		return;
	}
	if (!quite)	r_printf("\nOn the first run, check the CPU instruction set ... \n");
	 
	 struct cpu_x86 cpuinfo;
	 detect_host(&cpuinfo);
	 if (!quite) print_cpuinfo(&cpuinfo);
	 i386_cpuid_caches(quite);

#if !defined(SOLARIS_COMPILER) && defined(TARGET_64) && !defined(ARM64_OS)
	if (cpuinfo.HW_AVX512_F /*Foundation*/ && cpuinfo.HW_AVX512_BW && cpuinfo.HW_AVX512_DQ && cpuinfo.HW_AVX512_VL) {
		 SetupVectorFunction_AVX512();
		 SetupPCG_AVX512();
		 if (!quite)	  r_printf("CPU checking result: The AVX512-enabled library is used ... \n\n");
	}
	else if (cpuinfo.HW_AVX &&cpuinfo.HW_AVX2 && cpuinfo.HW_FMA3) {
		 SetupVectorFunction_AVX2();
		 SetupPCG_AVX2();
		 if (!quite)	 r_printf("CPU checking result: The AVX2-enabled library is used ...\n\n");
	}
	else {
		 SetupVectorFunction_Generic();
		 SetupPCG_GENERIC();
		 if (!quite)	 r_printf("CPU checking result: No AVX2/AVX512 is supported and the default library is used ...\n\n");
	}
#else
	 SetupVectorFunction_Generic();
	 SetupPCG_GENERIC();
	 if (!quite) r_printf("CPU checking result: No AVX2 is supported and the default library is used ...");
#endif

	 IS_CPU_INSTRUCTON_CHECKED = 1;

}

void SetupRoutines_UserChoice(int avxOption) {

#if !defined(SOLARIS_COMPILER) && defined(TARGET_64) && !defined(ARM64_OS)

	if (avxOption == 0) {
		SetupVectorFunction_Generic();
		SetupPCG_GENERIC();
	}
	else if (avxOption == 2) {
		SetupVectorFunction_AVX2();
		SetupPCG_AVX2();
	}
	else if (avxOption == 512) {
		SetupVectorFunction_AVX512();
		SetupPCG_AVX512();
	}
	else {
		SetupRoutines_AutoByCPU(1L); // Do not print out
	} 
#else
	//for SOaris system or Win32
	SetupVectorFunction_Generic();
	SetupPCG_GENERIC();
#endif

}


//https:// stackoverflow.com/questions/152016/detecting-cpu-architecture-compile-time
 
  const char *getBuild( void ) {
	  //Get current architecture, detectx nearly every architecture. Coded by Freak
        #if defined(__x86_64__) || defined(_M_X64)
        return "x86_64";
        #elif defined(i386) || defined(__i386__) || defined(__i386) || defined(_M_IX86)
        return "x86_32";
        #elif defined(__ARM_ARCH_2__)
        return "ARM2";
        #elif defined(__ARM_ARCH_3__) || defined(__ARM_ARCH_3M__)
        return "ARM3";
        #elif defined(__ARM_ARCH_4T__) || defined(__TARGET_ARM_4T)
        return "ARM4T";
        #elif defined(__ARM_ARCH_5_) || defined(__ARM_ARCH_5E_)
        return "ARM5"
        #elif defined(__ARM_ARCH_6T2_) || defined(__ARM_ARCH_6T2_)
        return "ARM6T2";
        #elif defined(__ARM_ARCH_6__) || defined(__ARM_ARCH_6J__) || defined(__ARM_ARCH_6K__) || defined(__ARM_ARCH_6Z__) || defined(__ARM_ARCH_6ZK__)
        return "ARM6";
        #elif defined(__ARM_ARCH_7__) || defined(__ARM_ARCH_7A__) || defined(__ARM_ARCH_7R__) || defined(__ARM_ARCH_7M__) || defined(__ARM_ARCH_7S__)
        return "ARM7";
        #elif defined(__ARM_ARCH_7A__) || defined(__ARM_ARCH_7R__) || defined(__ARM_ARCH_7M__) || defined(__ARM_ARCH_7S__)
        return "ARM7A";
        #elif defined(__ARM_ARCH_7R__) || defined(__ARM_ARCH_7M__) || defined(__ARM_ARCH_7S__)
        return "ARM7R";
        #elif defined(__ARM_ARCH_7M__)
        return "ARM7M";
        #elif defined(__ARM_ARCH_7S__)
        return "ARM7S";
        #elif defined(__aarch64__) || defined(_M_ARM64)
        return "ARM64";
        #elif defined(mips) || defined(__mips__) || defined(__mips)
        return "MIPS";
        #elif defined(__sh__)
        return "SUPERH";
        #elif defined(__powerpc) || defined(__powerpc__) || defined(__powerpc64__) || defined(__POWERPC__) || defined(__ppc__) || defined(__PPC__) || defined(_ARCH_PPC)
        return "POWERPC";
        #elif defined(__PPC64__) || defined(__ppc64__) || defined(_ARCH_PPC64)
        return "POWERPC64";
        #elif defined(__sparc__) || defined(__sparc)
        return "SPARC";
        #elif defined(__m68k__)
        return "M68K";
        #else
        return "UNKNOWN";
        #endif
 }
 



#include "abc_000_warning.h"
