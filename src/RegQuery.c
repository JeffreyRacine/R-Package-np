#include <Rdefines.h>

#ifdef WIN32
#define NOGDI
#define _WINGDI_H
#include <windows.h>

SEXP RegQuery(SEXP sexp_hkey, SEXP sexp_key) { 
    DWORD   cValues;              // number of values for key 
    DWORD   cchMaxValue;          // longest value name 
    DWORD   cbMaxValueData;       // longest value data 
    DWORD 	i, retCode= ERROR_SUCCESS; 
    DWORD 	cchValue; 
	DWORD 	cchData;

    HKEY hKey;
	SEXP sexp_valuedata;
	switch(INTEGER(sexp_hkey)[0]){
		case 1:
			retCode=RegOpenKeyEx(HKEY_CLASSES_ROOT, CHAR(STRING_ELT(sexp_key,0)), 0,  KEY_READ, &hKey);
			break;
		case 2:
  			retCode=RegOpenKeyEx(HKEY_CURRENT_USER, CHAR(STRING_ELT(sexp_key,0)), 0,  KEY_READ, &hKey);
			break;
		case 3:
  			retCode=RegOpenKeyEx(HKEY_LOCAL_MACHINE, CHAR(STRING_ELT(sexp_key,0)), 0,  KEY_READ, &hKey);
			break;
		case 4:
  			retCode=RegOpenKeyEx(HKEY_USERS, CHAR(STRING_ELT(sexp_key,0)), 0,  KEY_READ, &hKey);
			break;
	}
	if(retCode != ERROR_SUCCESS)
		return R_NilValue;
		
    retCode = RegQueryInfoKey(hKey,NULL,NULL,NULL,NULL,NULL,NULL, 
        &cValues,                // number of values for this key 
        &cchMaxValue,            // longest value name 
        &cbMaxValueData,         // longest value data 
        NULL, NULL);      
 
    if (cValues){
     	CHAR	cDATA[cbMaxValueData+1];
    	TCHAR  	achValue[cchMaxValue+1]; 
	 		
		PROTECT(sexp_valuedata=allocVector(STRSXP,2*cValues));
       	for (i=0, retCode=ERROR_SUCCESS; i<cValues; i++) {  
        	cchValue = cchMaxValue+1;
			cchData = cbMaxValueData+1; 
         	achValue[0] = '\0'; 
            retCode = RegEnumValue(hKey, i, achValue, &cchValue, NULL, NULL, cDATA,&cchData); 
            if (retCode == ERROR_SUCCESS ){ 
			  	SET_STRING_ELT(sexp_valuedata, 2*i, COPY_TO_USER_STRING(achValue));
			  	SET_STRING_ELT(sexp_valuedata, 2*i+1, COPY_TO_USER_STRING(cDATA));
            } 
        }
		UNPROTECT(1);
		return sexp_valuedata;    
	}
	return R_NilValue;
}
#endif

