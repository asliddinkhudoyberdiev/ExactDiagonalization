#include "einlesen.h"	

using namespace std;


double get_zahl1(char *str){
    return strtod(str, NULL);
}


double get_zahl2(char *str){
	while (*str && (isdigit(*str) || *str == '-' || *str == '+' || *str == '.' || *str == 'e'))
			str++;
		
	while (*str && !isdigit(*str) && *str != '-' && *str != '+' && *str != '.' && *str != 'e')
			str++;
			
    return strtod(str, NULL);
}


double get_zahl3(char *str){
	while (*str && (isdigit(*str) || *str == '-' || *str == '+' || *str == '.' || *str == 'e'))
			str++;
		
	while (*str && !isdigit(*str) && *str != '-' && *str != '+' && *str != '.' && *str != 'e')
			str++;
		

	while (*str && (isdigit(*str) || *str == '-' || *str == '+' || *str == '.' || *str == 'e'))
			str++;
		
	while (*str && !isdigit(*str) && *str != '-' && *str != '+' && *str != '.' && *str != 'e')
			str++;
			
    return strtod(str, NULL);
}


double get_zahl4(char *str){
	while (*str && (isdigit(*str) || *str == '-' || *str == '+' || *str == '.' || *str == 'e'))
			str++;
		
	while (*str && !isdigit(*str) && *str != '-' && *str != '+' && *str != '.' && *str != 'e')
			str++;
		

	while (*str && (isdigit(*str) || *str == '-' || *str == '+' || *str == '.' || *str == 'e'))
			str++;
		
	while (*str && !isdigit(*str) && *str != '-' && *str != '+' && *str != '.' && *str != 'e')
			str++;	
		
	while (*str && (isdigit(*str) || *str == '-' || *str == '+' || *str == '.' || *str == 'e'))
			str++;
		
	while (*str && !isdigit(*str) && *str != '-' && *str != '+' && *str != '.' && *str != 'e')
			str++;
			
    return strtod(str, NULL);
}