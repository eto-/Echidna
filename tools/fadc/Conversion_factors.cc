/* BOREXINO Conversion units program
 *
 * Author: Gromov Maxim <gromov@physics.msu.ru>
 * Maintainer: Gromov Maxim
 * 
 */

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#undef DEBUG
#define FIRST_RUN_TIME 1259918554 //It is the beginning time of the first run

#define CF_ASUM 1
#define CF_DSUM 2

using namespace std;

double GetConversionFactor(int Current_TimeStart, int Choose_Sum)
{
	ifstream ASum_DSum_to_MeV;
    char  Run_ConversionFactor[256]; 
    memset(Run_ConversionFactor,0,256);
    const size_t len =50;
    int TimeStart_int = 0;
    double ConversionFactor  = 0.;
    
    //Check the possibility that the run exist
    if(Current_TimeStart < FIRST_RUN_TIME)
    {
        #ifdef DEBUG
        cout << "Error: This run doesn't exist\n" << endl;
        #endif
        return -1.0;
    }
    
    //Read strings one by one from data file and 
    //compare Current_TimeStart with current value of the right time edge
    
    ASum_DSum_to_MeV.open("ASum_DSum_to_MeV.txt");
    string tmp;
    while(!ASum_DSum_to_MeV.eof())
    {
        getline(ASum_DSum_to_MeV,tmp);
        memcpy(Run_ConversionFactor,tmp.c_str(),tmp.size());
         char TimeStart_c[len];
         char ConversionFactor_c[len];
         char *start_copy  = strchr(Run_ConversionFactor,'\t');
         char *finish_copy = strrchr(Run_ConversionFactor,'\t');
         int i = 0;
         
         memset(TimeStart_c,0,len);
         memset(ConversionFactor_c,0,len);

          start_copy++;
          while(*start_copy != '\t')
          {
              strncpy(&TimeStart_c[i],start_copy,1);
			  i++;
              start_copy++;
           }
           TimeStart_int = atoi(TimeStart_c);

           #ifdef DEBUG
				cout << "TimeStart:  " << TimeStart_int << endl;
           #endif
           
           //The conversion factor is chosen for ASum or DSum
           if(Choose_Sum<=1) //ASum
           {
				strncpy(ConversionFactor_c,start_copy+2,6);
			}
			else //DSum
			{
				 strncpy(ConversionFactor_c,finish_copy+2,6);
			}
			ConversionFactor = atof(ConversionFactor_c); 
            
            if(Current_TimeStart <= TimeStart_int)
            {
				#ifdef DEBUG
					cout << "Conversion factor:  " << ConversionFactor << endl;
                 #endif
                  break;
            }
    }
    ASum_DSum_to_MeV.close();
    
    //Return the value or interpolation beyond the known values
    #ifdef DEBUG
		cout << "Conversion factor: " << ConversionFactor << endl;
    #endif
    return ConversionFactor;
}
