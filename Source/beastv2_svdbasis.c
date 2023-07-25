#include "abc_000_warning.h"
#include "abc_mat.h"

#include "abc_ide_util.h"
#include "abc_blas_lapack_lib.h"
#include "abc_mem.h"

static F32 __compute_a__(F32PTR x, F32PTR m, int P) {
	F32 xm = 0, mm = 0;
	for (int i = 0; i < P; i++) {
	// Bad values have been set to zeros
		xm += x[i] * m[i];
		mm += m[i] * m[i];
	}
	return  (mm == 0) ? 0: xm / mm;
}

static void __compute_intial_mean_(F32PTR mean, F32PTR Y, I32PTR nPtsPertime, int Ncycle, int P) {
	
	int p = 0;
	memset(mean, 0, sizeof(F32) * P);
	for (int i = 0; i < Ncycle * P; i++) {
		mean[p++] += Y[i];
		p = (p == P) ? 0 : p;
	}
	F32 sum = 0;
	for (int i = 0; i < P; i++) {
		mean[i] /= nPtsPertime[i];
		sum += mean[i] * mean[i];
	}
	sum = sqrtf(sum); 
	f32_mul_val_inplace(1 / sum, mean, P);
}

static  void commpute_residual(F32PTR Yerror, F32PTR Y, F32PTR X, I32PTR Igood, int Ngood, int N, int  K,F32PTR Xt, F32PTR B, F32PTR XtX ) {

	if (Ngood <=K) {
		memset(Yerror, 0, sizeof(F32) * N);
		return;
	}

	int p = 0;
	for (int i = 0; i < K; ++i) {
		for (int j = 0; j < N; ++j) {
			Xt[p] = Igood[j] ? X[p] : 0;
			p++;
		}
	}

	// Get XtY	
	r_cblas_sgemv(CblasColMajor, CblasTrans, N, K, 1.f, X, N, Y, 1L, 0.f, B, 1L);

	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, K, K, N, 1.0f, X, N, X, N, 0., XtX, K);

	r_LAPACKE_spotrf(LAPACK_COL_MAJOR, 'U', K, XtX, K);
	r_LAPACKE_spotrs(LAPACK_COL_MAJOR, 'U', K, 1L, XtX, K, B, K);

	F32PTR Yfit = Yerror;
	r_cblas_sgemv(CblasColMajor, CblasNoTrans, N, K, 1.f, X, N, B, 1L, 0.f, Yfit, 1L);

	F32 tmp = Yfit[0] * 0;
	if (tmp == tmp) {
		// not NAN or INF
		r_ippsSub_32f(Yfit, Y, Yerror, N); //Yerror=Y-Yfit 
	}
	else {
		f32_fill_val(0, Yerror, N);
	}


}

typedef struct SVDBasisMEM {
	int N, P, Ncycle;
	I32PTR nPtsPerTime; //P
	I32PTR goodTimeIndices; //P
	F32PTR Ytrue;// Ncycley*P;
	F32PTR Ycur;  //Ncyle(*P
	I32PTR Igood;// Ncycley*P;
	I32PTR NgoodPerPeriod;// Ncycley;
	F32PTR A, B; //P


	F32PTR M; // PXP
	F32PTR Mcopy; // PXP
	F32PTR XtX; // P*P;
	F32PTR Bcoeff; // P;

} SVDBasisMEM;



I64 Get_Alloc_SVDBasisMEM(int N, int P, SVDBasisMEM* s, VOID_PTR bufBase) {

	s = (s == NULL) ? (void *)(intptr_t)1 : s; // Make s NON-NULL, so memnodes_calc_offset won't stop at the first node

	int Ncycle = (N + P - 1) / P;
	
	MemNode nodes[] = {
		{.addr=&s->nPtsPerTime, .size=sizeof(I32)*P,      .align=4},
		{.addr = &s->goodTimeIndices     , .size = sizeof(I32) * P, .align = 4},
		{.addr = &s->Ytrue     , .size = sizeof(I32) * Ncycle*P, .align = 4},
		{.addr = &s->Ycur     , .size = sizeof(I32) * Ncycle * P, .align = 4},
		{.addr = &s->Igood     , .size = sizeof(I32) * Ncycle * P, .align = 4},
		{.addr = &s->NgoodPerPeriod     , .size = sizeof(I32) * Ncycle , .align = 4},
		{.addr = &s->A     , .size = sizeof(I32) * P , .align = 4},
		{.addr = &s->B     , .size = sizeof(I32) * P , .align = 4},
		{.addr = &s->M     , .size = sizeof(I32) * P*P , .align = 4},
		{.addr = &s->Mcopy     , .size = sizeof(I32) * P * P , .align = 4},
		{.addr = &s->XtX     , .size = sizeof(I32) * P * P , .align = 4},
		{.addr = &s->Bcoeff     , .size = sizeof(I32) * P , .align = 4},
		{.addr =NULL     , .size =0 , .align = 4},
	};
 	I64 totalSize= memnodes_calc_offsets(nodes, NULL);

	if (s> (SVDBasisMEM*)(intptr_t)1) {
		s->N = N;
		s->P = P;
		s->Ncycle = Ncycle;
		memnodes_assign_from_alignedbase(nodes, bufBase);
	}

	return totalSize;

}

void get_detrended_seasonalcmpt(F32PTR y, SVDBasisMEM* mem) {
	int N = mem->N;
	int P = mem->P;
	//F32 mean = f32_nanmean(y, N, NULL);
    //	f32_add_val_inplace(-mean, y, N);

	F32PTR X       = NULL;
	F32PTR Yerror  = y;
	simple_linear_regression_nan(y, X, N, NULL, Yerror);
	// Now Yerror is the detrended seasonal residual signal	
}



void compute_seasonal_svdbasis(F32PTR y, F32PTR Yout, int  Kmax, SVDBasisMEM * mem) {

	int N = mem->N;
	int P = mem->P; 


	//  FInd out how many valid points at each time of the period
	f32_compute_seasonal_avg(y, N, P, NULL, mem->nPtsPerTime);  

	int  Ptrue = 0;
	// Save the indices of the times that have vald values
	for (int i = 0; i < P; i++) {
		mem->goodTimeIndices[Ptrue] = i;
		Ptrue += mem->nPtsPerTime[i] > 0;
	}
	

	F32PTR Ytrue = mem->Ytrue;
	int    Ntrue;		
	if (Ptrue < P) {
		int p = 0;
		Ntrue = 0;
		for (int i = 0; i < N; i++) {
			Ytrue[Ntrue] = y[i];
			Ntrue        += mem->nPtsPerTime[p] > 0;
			p++;
			p = (p == P) ? 0 : p;
		}
				
		int pidx = 0;  // Remove all the zeros from nPtsPertime
		for (int i = 0; i < P; i++) {
			mem->nPtsPerTime[pidx] = mem->nPtsPerTime[i];
			pidx                   += mem->nPtsPerTime[i] > 0;
		}
	} else {
		Ntrue = N;
		f32_copy(y, Ytrue, N);
	}
	int Ncycle = mem->Ncycle;
	for (int i = Ntrue; i < Ptrue * Ncycle; i++) {
		// Add the extra NAN if there is a remainder period
		Ytrue[i] = getNaN();
	}

	I32PTR Igood          = mem->Igood;
	I32PTR NgoodPerPeriod = mem->NgoodPerPeriod;
	int idx = 0;
	for (int i = 0; i < Ncycle; i++) {
		int Ngood = 0;
		for (int p = 0; p < Ptrue; p++) {
			Ngood      += (Ytrue[idx] == Ytrue[idx]);
			Igood[idx]  = (Ytrue[idx] == Ytrue[idx]);
			idx++;
		}
		NgoodPerPeriod[i] = Ngood;
	}
	// Change NANs to zeros
	int Ngood=0;
	for (int i = 0; i < Ncycle * Ptrue; i++) {
		Ngood   += Ytrue[i] == Ytrue[i];
		Ytrue[i] = (Ytrue[i] == Ytrue[i]) ? Ytrue[i] : 0;
	}

	F32PTR Ycur = mem->Ycur;
	f32_copy(Ytrue, Ycur, Ncycle * Ptrue);
	
	//assert(Ptrue > 0);

	int KmaxTrue   = (Ngood+Ptrue-1)/Ptrue *Ptrue;
	KmaxTrue = min(KmaxTrue, Ptrue);
	KmaxTrue = min(KmaxTrue, Kmax);
 
	int kBasis = 0;
	while ( kBasis < KmaxTrue) {

		F32PTR mean = mem->M + kBasis * Ptrue;
		// Bad values in Ycur should have been replaced with zeros
		__compute_intial_mean_(mean, Ycur, mem->nPtsPerTime, Ncycle, Ptrue);

	   int iter_outer = 0;
	   F32 cLarange   = 0;
		while (iter_outer < 15) {

			F32PTR A = mem->A;  memset(A, 0, sizeof(F32) * Ptrue);
			F32PTR B = mem->B;  memset(B, 0, sizeof(F32) * Ptrue);
			
			for (int i = 0; i < Ncycle; i++) {
				I32PTR igood = Igood + i * Ptrue;
				F32PTR Ycol  = Ycur  + i * Ptrue;
				F32 a   = __compute_a__(Ycol, mean, Ptrue); // bad values in y have been replaced with zeros
				F32 a2  = a * a;
				for (int j = 0; j < Ptrue; j++) {
					A[j] += igood[j] ? a * Ycol[j] : 0;
					B[j] += igood[j] ? a2 : 0;
				}
			}
						
			int iter = 0;
			while (iter < 10) {
				F32 F0 = 0;
				F32 F1 = 0;
				for (int i = 0; i < Ptrue; i++) {
					F32 A2 = A[i] * A[i];
					F32 BpC2 = (B[i] + cLarange) * (B[i] + cLarange);
					F32 BpC3 = BpC2 * (B[i] + cLarange);
					F0 += A2 / BpC2;
					F1 += A2 / BpC3;
				}
				F0 = F0 - 1;
				F1 = -2 * F1;
				F32 cold = cLarange;
				cLarange = cLarange - F0 / F1;
				F32 diff = fabsf(cold - cLarange);
				if (diff < 1e-6) {
					break;
				}
				iter++;
			}

			F32 maxError = 0;
			for (int i = 0; i < Ptrue; ++i) {
				A[i] = A[i] / (B[i] + cLarange);
				maxError = max(maxError,  fabsf(mean[i] - A[i]) );
			}
		 
			if (maxError < 1e-5) {
				break;
			}

			//Save the temp result to mean
			f32_copy(A, mean, Ptrue);
			iter_outer++;
		}  //while (iter_outer < 5)
		kBasis++;
		
		for (int i = 0; i < Ncycle; i++) {
			commpute_residual(Ycur + i * Ptrue, Ytrue + i * Ptrue, mem->M, Igood + i * Ptrue, NgoodPerPeriod[i],
				             Ptrue, kBasis,	mem->Mcopy, mem->B, mem->XtX);
		}

		F32 maxAbsValaue=f32_absmax(Ycur, Ncycle*Ptrue);
		F32 meanSqr     = 0;
		for (int i = 0; i < Ncycle * Ptrue; i++) {
			meanSqr += Ycur[i] * Ycur[i];
		}

		//r_printf("maxAbsValaue: %f | MSE %f\n", maxAbsValaue, sqrtf(meanSqr/(Ncycle*Ptrue)) );
		if (maxAbsValaue < 1e-4) {
			break;
		}
		
	}

	for (int i = 0; i < kBasis; i++) {
		F32PTR mean_ext = Yout   + i * N;
		F32PTR mean     = mem->M + i * P;
		// mean_ext[goodIndices]=mean
		f32_scatter_vec_byindex(mean_ext, mem->goodTimeIndices, mean, Ptrue);
		
		// There are bad values to be filled
		if (Ptrue < P) {	        
			f32_interp1dvec_cycled_inplace(mean_ext, P, mem->goodTimeIndices, Ptrue);
		}
		F32 avg = f32_sum(mean_ext, P) / P;
		f32_add_val_inplace(-avg, mean_ext, P);
	}

	for (int i = kBasis; i < Kmax; i++) {
		F32PTR mean_ext = Yout + i * N;
		memset(mean_ext, 0, sizeof(F32)* P);
	}
	//f32_copy(mem->M, Yout, Ptrue* kBasis);
	for (int i = 0; i < Kmax; i++) {
		F32PTR mean_ext = Yout + i * N;
		f32_rep_vec1d_upto_inplace(mean_ext, P, N);
	}
}

void compute_seasonal_svdbasis_from_originalY(F32PTR y, int N, int P, F32PTR Yout, int  Kmax, VOID_PTR BUF) {

	SVDBasisMEM MEM;
	Get_Alloc_SVDBasisMEM(N, P, &MEM, BUF);
	F32PTR Y = MEM.Ycur;
	f32_copy(y, Y, N);
	get_detrended_seasonalcmpt(Y, &MEM);
	compute_seasonal_svdbasis(Y, Yout, Kmax, &MEM);
}

void compute_seasonal_svdbasis_from_seasonalY(F32PTR y, int N, int P, F32PTR Yout, int  Kmax, VOID_PTR BUF) {

	SVDBasisMEM MEM;
	Get_Alloc_SVDBasisMEM(N, P, &MEM, BUF);
	compute_seasonal_svdbasis(y, Yout, Kmax, &MEM);
}
#include "abc_000_warning.h"
 

