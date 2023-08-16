#include "math.h"

#include "abc_000_warning.h"
#include "abc_math.h"
 

//http: //www.machinedlearnings.com/2011/06/fast-approximate-logarithm-exponential.html
F32 fastlog1(F32 x)
{
	register union { F32 f; U32 i; } vx = { x };
	register union{ U32 i; F32 f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
	vx.f = (F32)vx.i* 1.1920928955078125e-7f*0.69314718f;
	vx.f = vx.f - 124.22551499f*0.69314718f	- 1.498030302f*0.69314718f * mx.f
			- 1.72587999f*0.69314718f / (0.3520887068f + mx.f);
	return vx.f;
}
F32 fastlog2(F32 x) {
	register union { F32 f; U32 i; } vx = { x };
	register union { U32 i; F32 f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
	vx.f = (F32) ( (double)vx.i* (double)(1.1920928955078125e-7f*0.69314718f) );
	vx.f = vx.f - 124.2098722217282f*0.69314718f - 1.502704726439149f*0.69314718f * mx.f
				- 1.746553042329640f*0.69314718f / (0.356745518976678f + mx.f);
	return vx.f;
}
F32 fastlog(F32 x)
{
	register union { F32 f; U32 i; } vx = { x };
	register union { U32 i; F32 f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
	vx.f = (F32)vx.i* (1.1920928955078125e-7f*0.69314718f);
	//vx.f = vx.f - 124.9998025915288f*0.69314718f+mx.f*(-0.967552023489078*0.69314718f - 0.151061326562556*0.69314718f*mx.f) +
	//	- 1.091678328444786f*0.69314718f / (0.238301594719865f + mx.f);
	vx.f = vx.f - 125.5233166734556f*0.69314718f + mx.f*(-0.413356886671142 + mx.f*(-0.472721975352920f + 0.078018528401178f*mx.f))*0.69314718f +
		-0.787757784962750f*0.69314718f / (0.1781810261970705f + mx.f);
	return vx.f;
}
F32 sum_log_diag(F32PTR p, I32 K)
{
	F32 x = 0;
	for (I32 i = 0; i < K; i++)
	{
		register union { F32 f; U32 i; } vx = { *p };
		register union{ U32 i; F32 f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
		vx.f = (F32)vx.i* (1.1920928955078125e-7f*0.69314718f);
		//vx.f = vx.f - 124.9998025915288f*0.69314718f+
		//mx.f*(-0.967552023489078*0.69314718f - 0.151061326562556*0.69314718f*mx.f) +
		//	- 1.091678328444786f*0.69314718f / (0.238301594719865f + mx.f);
		vx.f = vx.f - 125.5233166734556f*0.69314718f + mx.f*(-0.413356886671142 + mx.f*(-0.472721975352920 + 0.078018528401178*mx.f))*0.69314718f +
			-0.787757784962750f*0.69314718f / (0.1781810261970705f + mx.f);
		x += vx.f;
		p = p + K + 1;
	}
	return x;
}
F32 sum_log_diagv2(F32PTR p, I32 K)
{
	F64 sumlog  = 0;

	F64 cumprod = 1.;
	for (I32 i = 0; i < K; i++)	{
		F64  x           = *p;
		F64  cumprod_new = cumprod*x;
		if (cumprod_new > 1e-305 && cumprod_new < 1e305){
			cumprod = cumprod_new;
		} else	{
			sumlog  += log(cumprod);
			cumprod  = x;
		}
		p = p + K + 1;
	}

	sumlog += log(cumprod);
	return (F32) sumlog;
}
F32 sumlog(F32PTR p, I32 K)
{
	F64 sumlog = 0;

	F64 cumprod = 1.;
	for (I32 i = 0; i < K; i++) {
		F64  x = *p++;
		F64  cumprod_new = cumprod * x;
		if (cumprod_new > 1e-305 && cumprod_new < 1e305) {
			cumprod = cumprod_new;
		} else	{
			sumlog += log(cumprod);
			cumprod = x;
		}	 
	}

	sumlog += log(cumprod);
	return (F32) sumlog;
}
F32 fastexp(F32 x){

	x = (x > 88.7f) ? 88.7f : x;
	x = (x < -51.f) ? -51.f : x;

	x = 1.442695040f*x;	
	register float z = x - (F32)((int)x) + ((x < 0) ? 1.0f : 0.0f);
	register  union { U32 i; F32 f; } v;
	//v.i=  (1 << 23) * (x + 121.2740575f + 27.7280233f / (4.84252568f - z) - 1.49012907f * z);
	v.i = (U32)     (    8388608.f * (x + 121.2740575f + 27.7280233f / (4.84252568f - z) - 1.49012907f * z)          );
	return v.f;
}


/*
//https://bduvenhage.me/performance/machine_learning/2019/06/04/fast-exp.html
double fast_exp(double x) {
	union { double d; int32_t i[2]; } u;
	u.i[0] = 0;
	u.i[1] = int32_t(double((1 << 20) / log(2.0)) * x + double((1 << 20) * 1023 - 0));
	return u.d;
}
*/
/*
F32 mFast_Log2(F32 val) {
union { F32 val; I32 x; } u = { val };
register float log_2 = (F32)(((u.x >> 23) & 255) - 128);
u.x = (u.x & 0x007FFFFF) | 0x3f800000;
log_2 += ((-0.3358287811f) * u.val + 2.0f) * u.val - 0.65871759316667f;
return (log_2*0.69314718f);
}
*/
F32 fastsqrt (F32 x)
{//https ://en.wikipedia.org/wiki/Methods_of_computing_square_roots
	register union {
		U32 i;
		F32   f;
	} v;
	v.f = x; /* Same bits, but as an int */
	v.i -= 1 << 23; /* Subtract 2^m. */
	v.i >>= 1; /* Divide by 2. */
	v.i += 1 << 29; /* Add ((b + 1) / 2) * 2^m. */
	v.f = (v.f + x / v.f);
	v.f = (v.f*0.25f + x / v.f);
	return v.f;
}









 
double r8_normal_01_cdf_inverse(double p);
 

static double r8poly_value(int n, double a[], double x) {
	//  Purpose:   R8POLY_VALUE evaluates a double precision polynomial.
	//   Discussion:
	//
	//    For sanity's sake, the value of N indicates the NUMBER of 
	//    coefficients, or more precisely, the ORDER of the polynomial,
	//    rather than the DEGREE of the polynomial.  The two quantities
	//    differ by 1, but cause a great deal of confusion.
	//
	//    Given N and A, the form of the polynomial is:
	//      p(x) = a[0] + a[1] * x + ... + a[n-2] * x^(n-2) + a[n-1] * x^(n-1)
	//  Parameters:
	//
	//    Input, int N, the order of the polynomial.
	//    Input, double A[N], the coefficients of the polynomial.
	//    A[0] is the constant term.
	//    Input, double X, the point at which the polynomial is to be evaluated.
	//    Output, double R8POLY_VALUE, the value of the polynomial at X.

	double value = 0;
	for (int i = n - 1; 0 <= i; i--) {
		value = value * x + a[i];
	}
	return value;
}
 
// people.math.sc.edu/Burkardt/cpp_src/asa241/asa241.html
double normal_cdf_inverse(double p) {
//  Purpose:   R8_NORMAL_01_CDF_INVERSE inverts the standard normal CDF.
//  Discussion:   The result is accurate to about 1 part in 10**16.
//
//  Reference:
//    Michael Wichura,
//    The Percentage Points of the Normal Distribution,
//    Algorithm AS 241,
//    Applied Statistics,
//    Volume 37, Number 3, pages 477-484, 1988.
//
//  Parameters:
//
//    Input, double P, the value of the cumulative probability 
//    densitity function.  0 < P < 1.  If P is outside this range, an "infinite"
//    value is returned.

//
    static double a[8] = {
      3.3871328727963666080,     1.3314166789178437745e+2,      1.9715909503065514427e+3,  1.3731693765509461125e+4,
      4.5921953931549871457e+4,  6.7265770927008700853e+4,      3.3430575583588128105e+4,  2.5090809287301226727e+3 };

    static double b[8] = {
      1.0,                       4.2313330701600911252e+1,      6.8718700749205790830e+2,  5.3941960214247511077e+3,
      2.1213794301586595867e+4,  3.9307895800092710610e+4,       2.8729085735721942674e+4,  5.2264952788528545610e+3 };

    static double c[8] = {
      1.42343711074968357734,     4.63033784615654529590,      5.76949722146069140550,     3.64784832476320460504,
      1.27045825245236838258,     2.41780725177450611770e-1,     2.27238449892691845833e-2,  7.74545014278341407640e-4 };

    static double d[8]   = {
      1.0,                        2.05319162663775882187,       1.67638483018380384940,     6.89767334985100004550e-1,
      1.48103976427480074590e-1,  1.51986665636164571966e-2,      5.47593808499534494600e-4,  1.05075007164441684324e-9 };

    static double e[8] = {
      6.65790464350110377720,     5.46378491116411436990,      1.78482653991729133580,     2.96560571828504891230e-1,
      2.65321895265761230930e-2,  1.24266094738807843860e-3,   2.71155556874348757815e-5,  2.01033439929228813265e-7 };

    static double f[8] = {
      1.0,                        5.99832206555887937690e-1,       1.36929880922735805310e-1,  1.48753612908506148525e-2,
      7.86869131145613259100e-4,  1.84631831751005468180e-5,      1.42151175831644588870e-7,  2.04426310338993978564e-15 };


    if (p <= 0.0)  return -INFINITY;
    if (p >= 1.)   return INFINITY;
 

 

    double split1 = 0.425;    
    double q      = p - 0.5;
    if (fabs(q) <= split1)    {
		double const1 = 0.180625;
        double r      = const1 - q * q;
		return q * r8poly_value(8, a, r) / r8poly_value(8, b, r);
    }  else   {
		
		double r;

		r = (q<0) ? p: 1.0-p; // r >0;
        r = sqrt(-log(r));

		double split2 = 5.0;
		double result;
        if (r <= split2)      {
			double const2 = 1.6;
            r = r - const2;
			result=r8poly_value(8, c, r) / r8poly_value(8, d, r);
        }  else  {
            r = r - split2;
			result = r8poly_value(8, e, r) / r8poly_value(8, f, r);
        }

        if (q < 0.0)     {
			result = -result;
        }

		return result;
    }

 
}
 


 


#include "abc_000_warning.h"