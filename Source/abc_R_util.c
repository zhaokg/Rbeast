#include "abc_R_util.h"
DISABLE_MANY_WARNINGS
static char fileID UNUSED_DECORATOR='a';
#if R_INTERFACE==1
SEXP getListElement(SEXP list,const char *str)
{
	SEXP elmt=NULL; 
	SEXP names=getAttrib(list,R_NamesSymbol);
	for (int i=0; i < length(list); i++)
	if (strcmp(CHAR(STRING_ELT(names,i)),str)==0) {
		elmt=VECTOR_ELT(list,i);
		break;
	}
	return elmt;
}
#endif
ENABLE_MANY_WARNINGS
