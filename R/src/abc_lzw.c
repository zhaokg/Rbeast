
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>  //malloc
#include <string.h> //memset

#include "abc_000_warning.h"

#include "abc_ide_util.h" //r_printf



#undef LZW_COMPAT              /* include backwards compatibility code */

#define LZW_CHECKEOS            /* include checks for strips w/o EOI code */

#define MAXCODE(n)	((1L<<(n))-1)
/*
* The TIFF spec specifies that encoded bit
* strings range from 9 to 12 bits.
*/
#define BITS_MIN        9               /* start with 9 bits */
#define BITS_MAX        12              /* max of 12 bit strings */
/* predefined codes */
#define CODE_CLEAR      256             /* code to clear string table */
#define CODE_EOI        257             /* end-of-information code */
#define CODE_FIRST      258             /* first free code entry */
#define CODE_MAX        MAXCODE(BITS_MAX)
#define HSIZE           9001L           /* 91% occupancy */
#define HSHIFT          (13-8)
#ifdef LZW_COMPAT
/* NB: +1024 is for compatibility with old files */
#define CSIZE           (MAXCODE(BITS_MAX)+1024L)
#else
#define CSIZE           (MAXCODE(BITS_MAX)+1L)
#endif

/*
* State block for each open TIFF file using LZW
* compression/decompression.  Note that the predictor
* state block must be first in this data structure.
*/
typedef struct {
	//TIFFPredictorState predict;     /* predictor super class */

	unsigned short  nbits;          /* # of bits/code */
	unsigned short  maxcode;        /* maximum code for lzw_nbits */
	unsigned short  free_ent;       /* next free entry in hash table */
	unsigned long   nextdata;       /* next bits of i/o */
	long            nextbits;       /* # of valid bits in lzw_nextdata */

	int             rw_mode;        /* preserve rw_mode from init */
} LZWBaseState;

#define lzw_nbits       base.nbits
#define lzw_maxcode     base.maxcode
#define lzw_free_ent    base.free_ent
#define lzw_nextdata    base.nextdata
#define lzw_nextbits    base.nextbits

/*
* Encoding-specific state.
*/
typedef uint16_t hcode_t;			/* codes fit in 16 bits */
typedef struct {
	long	hash;
	hcode_t	code;
} hash_t;

/*
* Decoding-specific state.
*/
typedef struct code_ent {
	struct code_ent *next;
	unsigned short	length;		/* string len, including this token */
	unsigned char	value;		/* data value */
	unsigned char	firstchar;	/* first token of string */
} code_t;


typedef struct {
	LZWBaseState base;

	/* Decoding specific data */
	long    dec_nbitsmask;		/* lzw_nbits 1 bits, right adjusted */
	long    dec_restart;		/* restart count */
#ifdef LZW_CHECKEOS
	uint64_t  dec_bitsleft;		/* available bits in raw data */
	int64_t old_tif_rawcc;         /* value of tif_rawcc at the end of the previous TIFLZWDecode() call */
#endif
	//decodeFunc dec_decode;		/* regular or backwards compatible */
	code_t* dec_codep;		/* current recognized code */
	code_t* dec_oldcodep;		/* previously recognized code */
	code_t* dec_free_entp;		/* next free entry */
	code_t* dec_maxcodep;		/* max available entry */
	code_t* dec_codetab;		/* kept separate for small machines */

	/* Encoding specific data */
	int     enc_oldcode;		/* last code encountered */
	long    enc_checkpoint;		/* point at which to clear table */
#define CHECK_GAP	10000		/* enc_ratio check interval */
	long    enc_ratio;		/* current compression ratio */
	long    enc_incount;		/* (input) data bytes encoded */
	long    enc_outcount;		/* encoded (output) bytes */
	uint8_t*  enc_rawlimit;		/* bound on tif_rawdata buffer */
	hash_t* enc_hashtab;		/* kept separate for small machines */
} LZWCodecState;

#define LZWState(tif)		((LZWBaseState*) (tif)->tif_data)
#define DecoderState(tif)	((LZWCodecState*) LZWState(tif))
#define EncoderState(tif)	((LZWCodecState*) LZWState(tif))

 

/*
* LZW Decoder.
*/

#ifdef LZW_CHECKEOS
/*
* This check shouldn't be necessary because each
* strip is suppose to be terminated with CODE_EOI.
*/
#define	NextCode(_tif, _sp, _bp, _code, _get) {				\
if ((_sp)->dec_bitsleft < (uint64_t)nbits) { \
	r_printf("LZWDecode: Strip not terminated with EOI code"); \
	_code = CODE_EOI;					\
} \
else {	\
	_get(_sp, _bp, _code);					\
	(_sp)->dec_bitsleft -= nbits;				\
}								\
}
#else
#define	NextCode(tif, sp, bp, code, get) get(sp, bp, code)
#endif

#define	GetNextCode(sp, bp, code) {				\
	nextdata = (nextdata<<8) | *(bp)++;			\
	nextbits += 8;						\
	if (nextbits < nbits) {					\
		nextdata = (nextdata<<8) | *(bp)++;		\
		nextbits += 8;					\
	}							\
	code = (hcode_t)((nextdata >> (nextbits-nbits)) & nbitsmask);	\
	nextbits -= nbits;					\
}


static LZWCodecState* LZWSetupDecode(void)
{
	LZWCodecState* sp = NULL;
	int code;


		/*
		* Allocate state block so tag methods have storage to record
		* values.
		*/
	sp = (LZWCodecState*)malloc(sizeof(LZWCodecState));

	sp->dec_codetab = NULL;
	//sp->dec_decode = NULL;

		/*
		* Setup predictor setup.
		*/
		//(void)TIFFPredictorInit(tif);


	//assert(sp != NULL);

	if (sp->dec_codetab == NULL) {
		sp->dec_codetab = (code_t*)malloc(CSIZE*sizeof (code_t));
		if (sp->dec_codetab == NULL) {
			r_printf("No space for LZW code table");
			return (0);
		}
		/*
		* Pre-load the table.
		*/
		code = 255;
		do {
			sp->dec_codetab[code].value = (unsigned char)code;
			sp->dec_codetab[code].firstchar = (unsigned char)code;
			sp->dec_codetab[code].length = 1;
			sp->dec_codetab[code].next = NULL;
		} while (code--);
		/*
		* Zero-out the unused entries
		*/
		/* Silence false positive */
		/* coverity[overrun-buffer-arg] */
		memset(&sp->dec_codetab[CODE_CLEAR], 0,
			(CODE_FIRST - CODE_CLEAR) * sizeof (code_t));
	}
	return (sp);
}

/*
* Setup state for decoding a strip.
*/
static int LZWPreDecode(LZWCodecState* sp )
{
 
 
	//assert(sp != NULL);
	if (sp->dec_codetab == NULL)
	{
		//tif->tif_setupdecode(tif);
		if (sp->dec_codetab == NULL)
			return (0);
	}

	/*
	* Check for old bit-reversed codes.
	*/

{
		sp->lzw_maxcode = MAXCODE(BITS_MIN) - 1;
		//sp->dec_decode = LZWDecode;
	}
	sp->lzw_nbits = BITS_MIN;
	sp->lzw_nextbits = 0;
	sp->lzw_nextdata = 0;

	sp->dec_restart = 0;
	sp->dec_nbitsmask = MAXCODE(BITS_MIN);
#ifdef LZW_CHECKEOS
	sp->dec_bitsleft = 0;
	sp->old_tif_rawcc = 0;
#endif
	sp->dec_free_entp = sp->dec_codetab + CODE_FIRST;
	/*
	* Zero entries that are not yet filled in.  We do
	* this to guard against bogus input data that causes
	* us to index into undefined entries.  If you can
	* come up with a way to safely bounds-check input codes
	* while decoding then you can remove this operation.
	*/
	memset(sp->dec_free_entp, 0, (CSIZE - CODE_FIRST)*sizeof (code_t));
	sp->dec_oldcodep = &sp->dec_codetab[-1];
	sp->dec_maxcodep = &sp->dec_codetab[sp->dec_nbitsmask - 1];
	return (1);
}



int  LZWDecode_TIFF(uint8_t* tif_rawcp, uint8_t* op0,
int64_t tif_rawcc, int64_t occ0 )
{
	LZWCodecState * sp = LZWSetupDecode();
	LZWPreDecode(sp);

	char *op = (char*)op0;
	long occ = (long)occ0;
	char *tp;
	unsigned char *bp;
	hcode_t code;
	int len;
	long nbits, nextbits, nbitsmask;
	unsigned long nextdata;
	code_t *codep, *free_entp, *maxcodep, *oldcodep;

	//(void)s;
	//assert(sp != NULL);
	//assert(sp->dec_codetab != NULL);

	/*
	Fail if value does not fit in long.
	*/
	if ((int64_t)occ != occ0)
		return (0);
	/*
	* Restart interrupted output operation.
	*/
	if (sp->dec_restart) {
		long residue;

		codep = sp->dec_codep;
		residue = codep->length - sp->dec_restart;
		if (residue > occ) {
			/*
			* Residue from previous decode is sufficient
			* to satisfy decode request.  Skip to the
			* start of the decoded string, place decoded
			* values in the output buffer, and return.
			*/
			sp->dec_restart += occ;
			do {
				codep = codep->next;
			} while (--residue > occ && codep);
			if (codep) {
				tp = op + occ;
				do {
					*--tp = codep->value;
					codep = codep->next;
				} while (--occ && codep);
			}
			return (1);
		}
		/*
		* Residue satisfies only part of the decode request.
		*/
		op += residue;
		occ -= residue;
		tp = op;
		do {
			int t;
			--tp;
			t = codep->value;
			codep = codep->next;
			*tp = (char)t;
		} while (--residue && codep);
		sp->dec_restart = 0;
	}

	bp = (unsigned char *) tif_rawcp;
#ifdef LZW_CHECKEOS
	sp->dec_bitsleft += (((uint64_t) tif_rawcc - sp->old_tif_rawcc) << 3);
#endif
	nbits = sp->lzw_nbits;
	nextdata = sp->lzw_nextdata;
	nextbits = sp->lzw_nextbits;
	nbitsmask = sp->dec_nbitsmask;
	oldcodep = sp->dec_oldcodep;
	free_entp = sp->dec_free_entp;
	maxcodep = sp->dec_maxcodep;

	while (occ > 0) {
		NextCode(tif, sp, bp, code, GetNextCode);
		if (code == CODE_EOI)
			break;
		if (code == CODE_CLEAR) {
			do {
				free_entp = sp->dec_codetab + CODE_FIRST;
				//_TIFFmemset
				memset(free_entp, 0,
					(CSIZE - CODE_FIRST) * sizeof (code_t));
				nbits = BITS_MIN;
				nbitsmask = MAXCODE(BITS_MIN);
				maxcodep = sp->dec_codetab + nbitsmask - 1;
				NextCode(tif, sp, bp, code, GetNextCode);
			} while (code == CODE_CLEAR);	/* consecutive CODE_CLEAR codes */
			if (code == CODE_EOI)
				break;
			if (code > CODE_CLEAR) {
				r_printf("LZWDecode: Corrupted LZW table at scanline  d");
				return (0);
			}
			*op++ = (char)code;
			occ--;
			oldcodep = sp->dec_codetab + code;
			continue;
		}
		codep = sp->dec_codetab + code;

		/*
		* Add the new entry to the code table.
		*/
		if (free_entp < &sp->dec_codetab[0] ||
			free_entp >= &sp->dec_codetab[CSIZE]) {
			r_printf(			"Corrupted LZW table at scanline");
			return (0);
		}

		free_entp->next = oldcodep;
		if (free_entp->next < &sp->dec_codetab[0] ||
			free_entp->next >= &sp->dec_codetab[CSIZE]) {
			r_printf("Corrupted LZW table at scanline ");
			return (0);
		}
		free_entp->firstchar = free_entp->next->firstchar;
		free_entp->length = free_entp->next->length + 1;
		free_entp->value = (codep < free_entp) ?
			codep->firstchar : free_entp->firstchar;
		if (++free_entp > maxcodep) {
			if (++nbits > BITS_MAX)		/* should not happen */
				nbits = BITS_MAX;
			nbitsmask = MAXCODE(nbits);
			maxcodep = sp->dec_codetab + nbitsmask - 1;
		}
		oldcodep = codep;
		if (code >= 256) {
			/*
			* Code maps to a string, copy string
			* value to output (written in reverse).
			*/
			if (codep->length == 0) {
				r_printf(
					"Wrong length of decoded string: "
					"data probably corrupted at scanlin "
					 );
				return (0);
			}
			if (codep->length > occ) {
				/*
				* String is too long for decode buffer,
				* locate portion that will fit, copy to
				* the decode buffer, and setup restart
				* logic for the next decoding call.
				*/
				sp->dec_codep = codep;
				do {
					codep = codep->next;
				} while (codep && codep->length > occ);
				if (codep) {
					sp->dec_restart = (long)occ;
					tp = op + occ;
					do  {
						*--tp = codep->value;
						codep = codep->next;
					} while (--occ && codep);
					if (codep)
						r_printf(	"Bogus encoding, loop in the code table; scanline ");
				}
				break;
			}
			len = codep->length;
			tp = op + len;
			do {
				int t;
				--tp;
				t = codep->value;
				codep = codep->next;
				*tp = (char)t;
			} while (codep && tp > op);
			if (codep) {
				r_printf("Bogus encoding, loop in the code table; scanline ");
				break;
			}
			//assert(occ >= len);
			op += len;
			occ -= len;
		}
		else {
			*op++ = (char)code;
			occ--;
		}
	}

	//tif->tif_rawcc -= (tmsize_t)((uint8_t*)bp - tif->tif_rawcp);
	//tif->tif_rawcp = (uint8_t*)bp;
#ifdef LZW_CHECKEOS
	sp->old_tif_rawcc = tif_rawcc;
#endif
	sp->lzw_nbits = (unsigned short)nbits;
	sp->lzw_nextdata = nextdata;
	sp->lzw_nextbits = nextbits;
	sp->dec_nbitsmask = nbitsmask;
	sp->dec_oldcodep = oldcodep;
	sp->dec_free_entp = free_entp;
	sp->dec_maxcodep = maxcodep;

	if (occ > 0) {

		r_printf("Not enough data at scanline  (shortbytes)");
		return (0);
	}
	free(sp);
	return (1);
}


/*
* LZW Encoding.
*/







/*

static void
LZWCleanup()
{
	(void)TIFFPredictorCleanup(tif);

	assert(tif->tif_data != 0);

	if (DecoderState(tif)->dec_codetab)
		_TIFFfree(DecoderState(tif)->dec_codetab);

	if (EncoderState(tif)->enc_hashtab)
		_TIFFfree(EncoderState(tif)->enc_hashtab);

	_TIFFfree(tif->tif_data);
	tif->tif_data = NULL;

	_TIFFSetDefaultCompressionState(tif);
}
*/



#include "abc_000_warning.h"