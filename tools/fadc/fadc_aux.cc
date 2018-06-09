#ifndef ECHIDNA_CYCLE
#define ECHIDNA_CYCLE 16
#endif

#ifndef FADC_ROOT_FILES_PATH
//#define FADC_ROOT_FILES_PATH "/bxstorage/fadc/rootfiles/cycle_16"
#define FADC_ROOT_FILES_PATH "/storage/gpfs_data/borexino/fadc/rootfiles/cycle_16"
#endif

#ifndef FADC_VALIDATED_LIST
//#define FADC_VALIDATED_LIST "/home/production/FADC_run_validation_out.txt"
#define FADC_VALIDATED_LIST "/storage/gpfs_data/borexino/users/litvinov/FADC_run_validation_out.txt"
#endif

//#define ECHIDNA_ROOT_FILES_PATH "/bxstorage/rootfiles/cycle_16"
#define ECHIDNA_ROOT_FILES_PATH "/storage/gpfs_data/borexino/rootfiles/cycle_16"

//#include "TSystem.h"

// Execute program in shell and get output in string
int exec_in_shell(string cmd, string &output) 
{
    FILE* pipe = gSystem->OpenPipe(cmd.c_str(), "r");
    
    if (!pipe) 
        return -1;
        
    char buffer[128];
    
    while(!feof(pipe)) 
    {
    	if(fgets(buffer, 128, pipe) != NULL)
    		output += buffer;
    }
    
    gSystem->ClosePipe(pipe);
    return 0;
}


// Check FADC run in validation list and search rootfile path using 'find'
string get_fadc_run_path(int runnum)
{
    string result, line;
    char runname[1000];
    
    sprintf(runname, "Run%06d", runnum);
    
    // search in FADC validated file
    ifstream validated_file(FADC_VALIDATED_LIST);
    
    if(!validated_file.is_open())
        return "";

    bool found = false;
    while(!validated_file.eof())
    {
        getline(validated_file, line);
        
        if (line.find(runname) != string::npos) 
        {
            found = true;
            break;
        }
    }
    
    validated_file.close();
    
    if (!found)
        return "";  // run is not validated
    
    sprintf(runname, "Run%06d_fadc_c16.root", runnum);
    
    if (!exec_in_shell(string("find ") + string(FADC_ROOT_FILES_PATH) + " -name " + string(runname), result))
    {
        // remove newline from the end
        return result.substr(0, result.length()-1); // will be empty if not found
    }
    else
        return "";  // error
}


// Find Echidna rootfile path using 'find'
string get_echidna_run_path(int runnum)
{
    string result, line;
    char runname[1000];
    
    sprintf(runname, "Run%06d", runnum);
    
    // search in Echidna validated file
    //ifstream validated_file("/home/production/run_validation_out.txt");
    ifstream validated_file("/storage/gpfs_data/borexino/users/litvinov/run_validation_out.txt");
    
    if(!validated_file.is_open())
        return "";

    bool found = false;
    while(!validated_file.eof())
    {
        getline(validated_file, line);
        
        if (line.find(runname) != string::npos) 
        {
            found = true;
            break;
        }
    }
    
    validated_file.close();
    
    if (!found)
        return "";  // run is not validated
    
    sprintf(runname, "Run%06d_c16.root", runnum);
    
    if (!exec_in_shell(string("find ") + string(ECHIDNA_ROOT_FILES_PATH) + " -name " + string(runname), result))
    {
        // remove newline from the end
        return result.substr(0, result.length()-1); // will be empty if not found
    }
    else
        return "";  // error
}
