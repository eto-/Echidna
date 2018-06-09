//
// root macro to check for echidna root file integrity
//
	
int CheckBarn(TFile* f);

int Check(const char* fname, int last_ev) {

	TFile *f=TFile::Open(fname);
	if ( f->IsZombie() ) {
		std::cerr << " File " << fname << " not found\n";
		f->Close();
		return 1;
	}
	f->cd();
	if (last_ev) {
	  TTree *bxtree=(TTree*)f->Get("bxtree");
	  int nentries = bxtree->GetEntries();
	  float diff = last_ev - nentries;
	  float a = diff / last_ev;
	  if ( a > 0.002 ) {
	    std::cerr << " Too few events in file " << fname << "Expected: " << last_ev << "  Found: " << nentries << "\n";
	    return 2;
	  }
	}
	int isThereBarn = CheckBarn(f);
	if (!isThereBarn){
		std::cerr<<"barn not found in file "<<fname<<"\n";
		return 3;
	}
	f->Close(); 
	//std::cout << " File " << fname << " is fine with " << nentries << " out of expected " << last_ev << std::endl;
	return 0;
}

int CheckBarn(TFile* f){
	bx_barn	*barn = (bx_barn*)f->FindObjectAny("bxbarn");
	if (!barn) return 0;
	return 1;
}

int CheckFiles(int start=0, int end=999999, bool root=false, const char* path_log=".", const char* path_root="/shed/cycle9/", const char* fname="valid_runs.list") {

	gROOT->Reset();

	int filenum[100000];
	int last_event[100000];
        int bad_list[5000];
        int bad_counter=0;

	FILE *fp = fopen(fname,"rt");
	if ( !fp ) {
		std::cerr << " File not found\n";
		return 0;
	}
	
	int nfiles=0;
	int a,b;
	while (fscanf(fp,"%d %d\n",&a,&b) == 2) {
		filenum[nfiles] = a;
		last_event[nfiles] = b;
		nfiles++;
	}
	fclose(fp);

	int err_log = 0;
	int err_critic = 0;
	int err_root = 0;
	int err_size = 0;
	int err_end = 0;
	int err_barn = 0;

	for(int i=0; i<nfiles; i++) {
		if ( filenum[i] >= start && filenum[i] <= end) {
			char s[500];

			sprintf(s,"%s/Run%06d_c9.log",path_log,filenum[i]);
			FILE *fp=fopen(s,"rt");
			if ( !fp ) {
				std::cerr << " Log file " << s << " not found.\n";
				err_log++;
				bad_list[bad_counter++] = filenum[i];
			} else {
				int error=0;
				int ended=0;
				while(fgets(s,450,fp)) {
					if ( strstr(s,"critic"))
						error = 1;
					if ( strstr(s,"Echidna ENDED"))
						ended=1;
				}
				fclose(fp);
				if ( error ) {
					err_critic++;
					if (!root) std::cerr << "Run " << filenum[i] << " has critic errors \n";
				        bad_list[bad_counter++] = filenum[i];
				}
				else {
					if (ended) {
						if ( root ) {
							sprintf(s,"%s/Run%06d_c9.root",path_root,filenum[i]);
							int ret = Check( s, last_event[i] );
							if ( ret == 1 )
								err_root++;
							if ( ret == 2 )
								err_size++;
							if ( ret != 0 )
								bad_list[bad_counter++] = filenum[i];
						}
					} else {
						err_end++;
						//if (!root) 
						std::cerr << "Run " << filenum[i] << " did not END \n";
						bad_list[bad_counter++] = filenum[i];
					}
				}
			}
		}
	}

    	std::cerr << " Errors because log file not found: " << err_log << std::endl;
    	std::cerr << " Errors because log file contains critic errors: " << err_critic << std::endl;
    	std::cerr << " Errors because job was not ended: " << err_end << std::endl;
    	std::cerr << " Errors because root file not found: " << err_root << std::endl;
    	std::cerr << " Errors because root file size is wrong: " << err_size << std::endl;
	std::cerr << std::endl;

	FILE *fp=fopen("wrong_runs.list","a+");
	for(int i=0; i<bad_counter; i++)
		fprintf(fp,"%d\n",bad_list[i]);
	fclose(fp);

	std::cerr << "Return value:="<< (err_root + err_size + err_log + err_critic + err_end) << std::endl;
	return  (err_root + err_size + err_log + err_critic + err_end);

}

