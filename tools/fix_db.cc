/* function to fix db data from raw data file     */
/* read data from a file (compressed) and         */
/* update run table in daq_config accordingly     */
/* M. Pallavicini - 31-12-2002 created            */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>
#include <string.h>
#include <libpq-fe.h>

/* event data structure (header only) */
struct Event_t {
  long bytes;
  long run;
  long id;
  long crateflag;
  long pointers[4];
  long time;
  long nsec;
  short trg1[4];
  long trg2[8];
};

/* helper structure */
struct RunData_t {
  int run;
  int files;
  time_t start;
  time_t end;
  int duration;
  int events;
  char start_s[200];
  char end_s[200];
};
static struct RunData_t data;

/* functions */
int GetRunData(int run, int nFiles, struct RunData_t *d,int);
int CheckDb(struct RunData_t *d); 
void DumpData(struct RunData_t *d);
int CountFiles(int);
int GetYesNo(const char*);

/* main program */
int main(int argc, char **argv) {
  
  /* check format */
  if ( argc < 2 ) {
    printf("Usage: FixRunTable run [runs]\n"); 
  }
  
  int nRuns = argc-1;
  
  for( int i=0; i<nRuns; i++) {
     
    /* count number of files for this run */
    int run = atol(argv[i+1]);
    int nFiles = CountFiles(run); 
  
    if ( nFiles >= 0 ) {
	    
      //printf("Processing run number %d\n",run);
      //printf("Number of files found: %d\n",nFiles); 
      
      /* get data from raw data files and compare with DB */
      int ret = GetRunData(run,nFiles,&data,9999999);
      if (ret) // try again with less events 
      	ret = GetRunData(run,nFiles,&data,ret-5);
      if ( ret == 0) {
        if ( CheckDb( &data ) ) {
	  printf("Run %d is up to date\n",data.run);
	}
      }
      else {
      	printf("Run %d corrupted\n",run);
      }
    } else {
      printf("No files found for run %d\n",run);
    }
  }
  exit(0);
}

/* count number of files for a given run      */
/* it assumes all files are available on disk */
int CountFiles(int run) {
  int nFiles=0;
  while(nFiles<100) {
    char fname[30];
    sprintf(fname,"Run%06d_%02d.out.gz",run,nFiles+1);
    FILE *fp = fopen(fname,"r");
    if (fp) {
      fclose(fp);
      //printf("Found %s\n",fname);
      nFiles++;
      continue;
    } 
    else {
        printf("File %s not found\n",fname);
        return nFiles;
    }
  }
  printf("Too many files for run %d\n",run);
  return nFiles;
}

static unsigned int evbuff[500000];

int GetRunData(int run, int nfiles, struct RunData_t *d, int max_events) {
  if (nfiles == 0) {
    d->run = run;
    d->start = d->end =  d->duration = 0;
    d->files = 0;
    d->events=0;
    DumpData(d);
    return 0;
  }
  //printf("Analyzing run %d (%d files)\n",run,nfiles);
  printf("Please be patient.");
  char file1[30];
  char file2[30];

  sprintf(file1,"Run%06d_01.out.gz",run);
  sprintf(file2,"Run%06d_%02d.out.gz",run,nfiles);
   
  /* get firsr event from first file */
  FILE* fp = (FILE*)gzopen(file1,"rb");
  if ( fp == NULL ) {
    printf("Error opening %s. Exiting. \n",file1);
    exit(1);
  }
  
  int len;
  //  int ret = gzread((void*)fp,&len,4);
  int ret = gzread((gzFile)fp,&len,4);
  if ( ret==0 ) {
    printf("Unexpected end of file %s (1)\n",file1);
    return -1;
  } 
  if ( ret==1 ) {
    printf("Error in file %s (1)\n",file1);
    return -1;
  } 
  evbuff[0]=len;
  len -= 4;
  //  ret = gzread((void*)fp,&(evbuff[1]),len);
  ret = gzread((gzFile)fp,&(evbuff[1]),len);
  if ( ret==0 ) {
    printf("Unexpected end of file %s (2)\n",file1);
    return -1;
  } 
  if ( ret==1 ) {
    printf("Error in file %s (2)\n",file1);
    return -1;
  } 
  
  /* get start time from first event */
  struct Event_t *ev = (struct Event_t*)&evbuff[0];
  time_t t_start = ev->trg2[2];
  /* do some checks */
  if ( ev->run != run || ev->id != 0 ) {
    printf("Error in first event!! Exiting \n");
    exit(1);
  } 

  /* if more than 1 file, go to last file */
  if ( nfiles>1) {
    //    gzclose((void*)fp); 
    gzclose((gzFile)fp); 
    fp = (FILE*)gzopen(file2,"rb");
    if ( fp == NULL ) {
      printf("Error opening %s. Exiting. \n",file2);
      exit(1);
    }
  }
  int count=0;
  /* read the file up to the end to get last event */ 
  while( 1 ) {
    if ((count++%4000)==0) {printf(".");fflush(stdout);}
    if (count > max_events) break;

    //    ret = gzread((void*)fp,&len,4);
    ret = gzread((gzFile)fp,&len,4);
    if ( ret == 0 ) {
      printf("End of file in %s\n",file2);
      break;
    }
    if ( ret==1 ) {
      printf("Error in file %s (3)\n",file2);
      return -1;
    } 
    evbuff[0]=len;
    len -=4;
    //    ret = gzread((void*)fp,&(evbuff[1]),len); 
    ret = gzread((gzFile)fp,&(evbuff[1]),len); 
    if ( ret==0 ) {
      printf("Unexpected end of file %s (3)\n",file2);
      return -1;
    } 
    if ( ret==1 ) {
      printf("Error in file %s (4)\n",file2);
      return -1;
    } 
  }
  printf("\nGetting data\n");
  ev = (struct Event_t*)&evbuff[0];
  time_t t_end = ev->trg2[2];
  /* do some checks */
  if ( ev->run != run || ev->id != ev->trg2[3] ) {
    printf("\n\n****WARNING!!! Last event is corrupted (%ld,%ld,%ld,%ld)\n",
		    run,ev->run,ev->id,ev->trg2[3]);
    return count;
  }
  
  /* file structure */
  d->run = run;
  d->start = t_start;
  d->end = t_end;
  d->duration = t_end - t_start;
  d->files = nfiles;
  d->events=ev->id+1;

  DumpData(d);

  return 0;
}

void DumpData(struct RunData_t *d) {
  printf("\t\tRun:\t%ld\n",d->run);
  printf("\t\tEvents:\t%ld\n",d->events);
  if ( d->start )
    printf("\t\tStart:\t%s",asctime(localtime(&(d->start))));
  else
    printf("\t\tStart:\t%s\n",d->start_s);
  if ( d->end )
    printf("\t\tEnd:\t%s",asctime(localtime(&(d->end))));
  else
    printf("\t\tEnd:\t%s\n",d->end_s);
  printf("\t\tDuration:\t%d\n",d->duration);
  printf("\t\tFiles:\t%ld\n",d->files);
}

int CheckDb(struct RunData_t *d) {

  PGconn *Connection;
  PGresult *Result;

  printf("Check Data Base (Table Run in daq_config) for run %d\n",d->run);

  char host[100];
  strcpy(host,"bxdb.lngs.infn.it");
  char myhost[200];
  if ( getenv("HOSTNAME") )
      strcpy(myhost,getenv("HOSTNAME"));
  else
      strcpy(myhost,"bxdb");
  //printf("My host %s\n",myhost);
  if ( strcmp(myhost,"bxweb")==0 ||
       strcmp(myhost,"bxdb")==0 ||
       strcmp(myhost,"bxoff1")==0 ||
       strcmp(myhost,"bxterm1")==0 ||
       strcmp(myhost,"bxbuild")==0 ) {
     strcpy(host,"bxdb");
  }
  const char* port = "5432";
  //printf("Connecting to %s at port %d\n",host,port);
  
  Connection = PQsetdbLogin( host, port, 0, 0, "daq_config", "storage", 0);

  if ( PQstatus(Connection) == CONNECTION_BAD ) {
    printf("Cannot connect to host %s, port %d, data base daq_config\n",
		    host, port);
    PQfinish( Connection );
    return 0; 
  }

  /* get data for this run */
  char query[1000];
  sprintf(query,"SELECT * FROM \"Run\" WHERE \"RunNumber\"=%d",d->run);
  
  Result = PQexec( Connection, query );
  
  if ( PQresultStatus( Result ) != PGRES_TUPLES_OK ) {
    printf("Error getting data from table Run for run %d\n",d->run);
    PQfinish( Connection );
    return 0;
  }

  int nRecs = PQntuples( Result );

  if ( nRecs < 1 ) {
    printf("No record found for run %d\n");
    return 0;
  }

  struct RunData_t db;
  db.start = 0;
  db.end = 0;
  db.duration=0;
  const char *val = PQgetvalue( Result, 0, PQfnumber(Result,"\"RunNumber\"") );
  sscanf(val,"%d",&(db.run));
  val = PQgetvalue( Result, 0, PQfnumber(Result,"\"Events\"") );
  sscanf(val,"%d",&(db.events));
  val = PQgetvalue( Result, 0, PQfnumber(Result,"\"NumberOfFiles\"") );
  sscanf(val,"%d",&(db.files));
  val = PQgetvalue( Result, 0, PQfnumber(Result,"\"Duration\"") );
  sscanf(val,"%d",&(db.duration));
  val = PQgetvalue( Result, 0, PQfnumber(Result,"\"StartTime\"")); 
  strcpy(db.start_s,val);
  val = PQgetvalue( Result, 0, PQfnumber(Result,"\"StopTime\""));
  strcpy(db.end_s,val);

  printf("The data base content is now:\n");
  DumpData( &db );

  if ( d->run != db.run ) {
    printf("Fatal error. Inconsistent run number!!!!! Exiting. \n");
    exit(1);
  }

  /* check number of events */
  if ( d->events != db.events ) {
    printf("The number of events does not match.\n");
    if ( GetYesNo("Should I update the data base") ) {
      sprintf(query,"UPDATE \"Run\" SET \"Events\"=%d WHERE \"RunNumber\"=%d",
		      d->events,d->run);
      Result = PQexec( Connection, query );
      if ( PQresultStatus( Result ) != PGRES_COMMAND_OK ) {
        printf("Error updating table Run for run %d\n",d->run);
        PQfinish( Connection );
        return 0;
      }
    }
    else
      printf("Number of events in Data base not changed\n");
  }
  /* check number of files */
  if ( d->files != db.files ) {
    printf("The number of files does not match.\n");
    if ( GetYesNo("Should I update the data base") ) {
      sprintf(query,"UPDATE \"Run\" SET \"NumberOfFiles\"=%d WHERE \"RunNumber\"=%d",
		      d->files,d->run);
      Result = PQexec( Connection, query );
      if ( PQresultStatus( Result ) != PGRES_COMMAND_OK ) {
        printf("Error updating table Run for run %d\n",d->run);
        PQfinish( Connection );
        return 0;
      }
    }
    else
      printf("Number of files in Data base not changed\n");
  }
  /* check time and duration */
  if ( d->duration != db.duration && (d->events!=db.events) && d->duration) {
    printf("The time stamps do not match.\n");
    if ( GetYesNo("Should I update the data base") ) {
      sprintf(query,"UPDATE \"Run\" SET \"Duration\"=%d WHERE \"RunNumber\"=%d",
		      d->duration,d->run);
      Result = PQexec( Connection, query );
      if ( PQresultStatus( Result ) != PGRES_COMMAND_OK ) {
        printf("Error 1 updating table Run for run %d\n",d->run);
        PQfinish( Connection );
        return 0;
      }
      sprintf(query,"UPDATE \"Run\" SET \"StopTime\"='%s' WHERE \"RunNumber\"=%d",
		      asctime(localtime(&(d->end))),d->run);
      Result = PQexec( Connection, query );
      if ( PQresultStatus( Result ) != PGRES_COMMAND_OK ) {
        printf("Error 2 updating table Run for run %d\n",d->run);
        PQfinish( Connection );
        return 0;
      }  
    }
    else
      printf("Number of files in Data base not changed\n");
  }
  

  PQfinish( Connection );
    
  return 1;
}

int GetYesNo(const char* question) {
  char c;
  printf("%s? [y/n] : ",question);
  while(1) {
    c=getchar();
    if (c=='y' || c=='Y' ) return 1;
    if (c=='n' || c=='N' ) return 0;
    if (c!='\n') printf("\nInvalid answer!\n%s? [y/n] : ",question);
  }
  return 0;
}
