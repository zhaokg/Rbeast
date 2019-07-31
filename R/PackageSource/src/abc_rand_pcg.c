
#include "abc_rand_pcg.h" 
#include "abc_datatype.h" 
#include <math.h>  
typedef struct pcg32_random_struct {
	uint64_t state;
	uint64_t inc;
} pcg32_random_t;
static pcg32_random_t globalStream={ 34354354ULL,56567456ULL };
void pcg_set_seed(uint64_t initstate,uint64_t initseq)
{
	uint32_t rnd;
	globalStream.state=0U;
	globalStream.inc=(initseq << 1u)|1u;
	pcg_random(&rnd,1);
	globalStream.state+=initstate;
	pcg_random(&rnd,1);
}
void pcg_random(uint32_t * _restrict rnd,int32_t N)
{
	uint64_t oldstate=globalStream.state;
	uint64_t shift=globalStream.inc;
	for (int i=0; i < N; i++)
	{
		uint32_t xorshifted=((oldstate >> 18u) ^ oldstate) >> 27u   ;
		uint32_t rot=oldstate >> 59u;
		*rnd++=(xorshifted >> rot)|(xorshifted << ((-rot) & 31));
		oldstate=oldstate * 6364136223846793005ULL+shift;
	}
	globalStream.state=oldstate;
}
typedef struct GAUSS_CONSTANT {
	float    x[64];
	uint32_t yRatio[63];
	float    yInverse[63];
	float    PARAM_R,INV_PARAM_R;
} GAUSS_CONSTANT;
static GAUSS_CONSTANT GAUSS={ .x={ 0,0.019584285230127,0.039176085503098,0.058782936068943,0.078412412733112, 
0.098072152488661,0.117769874579095,0.137513402144336,0.157310684610171,0.177169820991740,0.197099084294312,
0.217106947210130,0.237202109328788,0.257393526100938,0.277690439821577,0.298102412930487,0.318639363964375,
0.339311606538817,0.360129891789569,0.381105454763557,0.402250065321725,0.423576084201200,0.445096524985516,
0.466825122852590,0.488776411114670,0.510965806738247,0.533409706241281,0.556125593618691,0.579132162255556,
0.602449453164424,0.626099012346421,0.650104070647995,0.674489750196082,0.699283302383220,0.724514383492366,
0.750215375467941,0.776421761147928,0.803172565597918,0.830510878205399,0.858484474141832,0.887146559018876,
0.916556667533113,0.946781756301046,0.977897543940542,1.009990169249582,1.043158263318454,1.077515567040280,
1.113194277160929,1.150349380376008,1.189164350199337,1.229858759216589,1.272698641190536,1.318010897303537, 
1.366203816372098,1.417797137996268,1.473467577947101,1.534120544352547,1.601008664886076,1.675939722773444,
1.761670410363067,1.862731867421652,1.987427885929896,2.153874694061456,2.41755901623650 },
.yRatio={ 16773999,16767562,16761112,16754640,16748136,16741589,16734989,16728325,16721587,16714763, 
16707840,16700807,16693651,16686358,16678913,16671302,16663507,16655512,16647297,16638843,16630128,16621128,
16611818,16602170,16592154,16581736,16570879,16559544,16547684,16535250,16522186,16508431,16493913,16478554, 
16462265,16444942,16426470,16406715,16385522,16362712,16338074,16311363,16282287,16250497,16215576,16177016, 
16134195,16086343,16032497,15971429,15901559,15820813,15726414,15614561,15479906,15314686,15107183,14838856, 
14478538,13969524,13196768,11886086,9182634 },
.yInverse={ 0.167772160000000e8,0.167804337107029e8,0.167900954887396e8,0.168062273321792e8,0.168288727728840e8, 
0.168580931819422e8,0.168939682029590e8,0.169365963190389e8,0.169860955611753e8,0.170426043678181e8,0.171062826076636e8, 
0.171773127802660e8,0.172559014119663e8,0.173422806679528e8,0.174367102051039e8,0.175394792947299e8,0.176509092495599e8, 
0.177713561954893e8,0.179012142359128e8,0.180409190651708e8,0.181909520980666e8,0.183518451949447e8,0.185241860769783e8, 
0.187086245447135e8,0.189058796353774e8,0.191167478819898e8,0.193421128712637e8,0.195829563393308e8,0.198403710967247e8,
0.201155761397016e8,0.204099343877364e8,0.207249735920063e8,0.210624110937242e8,0.214241832835665e8,0.218124808367550e8,
0.222297910899452e8,0.226789493099905e8,0.231632011146443e8,0.236862789891804e8,0.242524967693829e8,0.248668672301306e8, 
0.255352496766927e8,0.262645369022485e8,0.270628943829048e8,0.279400696443366e8,0.289077971595910e8,0.299803352204569e8,
0.311751880874552e8,0.325140929293103e8,0.340243927554388e8,0.357409846346236e8,0.377091470037256e8,0.399887489452612e8, 
0.426607037407803e8,0.458372068395134e8,0.496786432318855e8,0.544228824001546e8,0.604390940468508e8,0.683340892537681e8,
0.791831158870118e8,0.950978898378280e8,1.208991317887296e8,1.706491807080511 },
.PARAM_R=2.41755901623650,
.INV_PARAM_R=0.413640367529367
};
void pcg_gauss(F32PTR RND,int N)
{
	for (int i=1; i <=N; i++)
	{
		uint32_t rnd[2];
		pcg_random(rnd,2);
		uint32_t U24=rnd[0];
		uint8_t  u8=U24 & 0xff;
		U24=U24 >> 8;
		int sign=(u8 & 0x80) ?+1 : -1;
		u8=u8 & 0x3f;
		float x;
		int64_t IDX=u8;
		if (u8 < 63)
		{
			float delta;
			delta=(GAUSS.x[IDX+1] - GAUSS.x[IDX]) * 2.328306436538696e-10f;
			while (1)
			{
				x=GAUSS.x[IDX]+delta*(float)rnd[1];
				if (U24 <=GAUSS.yRatio[IDX])
					break;
				if (U24 <=expf(-0.5f*x*x)* GAUSS.yInverse[IDX])
					break;
				pcg_random(rnd,2);
				U24=rnd[0] >> 8;
			}
		}
		else
		{
			float U1,U2;
			U2=(float)U24* 5.960464477539063e-08f; 
			while (1)
			{
				U1=rnd[1] * 2.328306436538696e-10f; 
				U1=1 - U1;
				x=GAUSS.PARAM_R - logf(U1) *GAUSS.INV_PARAM_R;
				if (exp(-GAUSS.PARAM_R * (x - 0.5f * GAUSS.PARAM_R)) * U2  < expf(-0.5f * x * x))
					break;
				pcg_random(rnd,2);
				U2=rnd[0] * 2.328306436538696e-10f;
			}
		}
		*RND++=x*sign;
	}
}
DISABLE_MANY_WARNINGS
void pcg_gamma(F32PTR rnd,float a,int N)
{
	float b;
	float c;
	float INV_MAX=2.328306436538696e-10f;
	float gam;
	float u[2];
	if (a >=1)
	{
		b=a - 1.0f;
		c=a+a+a - 0.75f;
		for (int i=0; i < N; i++)
		{
			while (1)
			{
				pcg_random((uint32_t*)u,2);
				u[0]=(float)(*(uint32_t *)&u[0]) * INV_MAX;
				u[1]=(float)(*(uint32_t *)&u[1]) * INV_MAX;
				float w=u[0] * (1 - u[0]);
				float y=sqrtf(c/w) * (u[0] - 0.5f);
				gam=b+y;
				if (gam >=0)
				{
					float z=64.0f * (w*w*w) * (u[1] * u[1]);
					if (z <=(1 - 2 * (y*y)/gam))
						break;
					float logZ=logf(z);
					if (b==0)
					{
						if (logZ <=-2 * y) break;
					}
					else
					{
						if (logZ <=2 * (b * logf(gam/b) - y)) break;
					}
				} 
			}
			*rnd++=gam;
		}
		return;
	}
	if (a > 0) 
	{
		b=1+.3678794f * a;
		for (int i=0; i < N; i++)
		{
			while (1)
			{
				pcg_random((uint32_t*)u,2);
				u[0]=(float)(*(uint32_t *)&u[0]) * INV_MAX;
				u[1]=(float)(*(uint32_t *)&u[1]) * INV_MAX;
				c=b * u[0];
				if (c < 1)
				{
					gam=expf(logf(c)/a);
					if (-logf(1 - u[1]) >=gam) break;
				}
				else
				{
					gam=-logf((b - c)/a);
					if (-logf(1 - u[1]) >=(1 - a) * logf(gam))
						break;
				}
			}
			*rnd++=gam;
		}
		return;
	}
	if (a < 0.f)
	{
		float nan=(float) 1e300;
		nan=nan*nan*0.f;
		for (int i=0; i < N; i++)
		{
			*rnd++=nan;
		}
		return;
	}
	if (a==0.f)
	{
		for (int i=0; i < N; i++)
		{
			*rnd++=0.f;
		}
		return;
	}
}
ENABLE_MANY_WARNINGS
void pcg_wishart_unit(F32PTR rnd,F32PTR tmp,int32_t m,int32_t df,char upperOrLower,char zerofill)
{
	int32_t N=(m - 1)*m/2;
	pcg_gauss(tmp,N);
	F32PTR tpt;
	tpt=rnd+1;
	for (int i=1; i<m  ; i++)
	{			
		int32_t  movElem=m - i;
		for (int j=0; j < movElem; j++)
		{
			*tpt++=*tmp++;
		}
		tpt=tpt+i+1;
	}
	tmp=rnd;
	for (int i=1; i <=m; i++)
	{
		pcg_gamma(tmp,(float)(df-i+1.)/2.f,1);
		*tmp=sqrtf((*tmp) * 2.f);
		tmp=tmp+(m+1);
	}
}
