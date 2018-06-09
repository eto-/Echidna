#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>

void usage () {
  std::cerr << "Usage fake_omon [-t time(s)] [-j jump] [-e events] [input_file=stdin]" << std::endl;
}

int main (int argc, char * argv[]) {
  struct timespec sleep = { tv_sec:1, tv_nsec:0 };
  int jump = 0;
  int nevents = -1;
  FILE *pfile = 0;
  gzFile gzfile;


  while (1) {
    const char valid_options[] = "t:j:e:";
    int c = ::getopt (argc, argv, valid_options);
    if (c == -1) break;

    char *error_ptr = 0;
    float t;
    switch (c) {
      case 't':
	t = ::strtod (optarg, &error_ptr);
	if (t <= 0 || (error_ptr && *error_ptr)) usage ();
	sleep.tv_sec = int (t);
	sleep.tv_nsec = int ((t - int (t)) * 1e9);
	break;
      case 'j':
	jump = ::strtol (optarg, &error_ptr, 0);
	if (jump <= 0 || (error_ptr && *error_ptr)) usage ();
	break;
      case 'e':
	nevents = ::strtol (optarg, &error_ptr, 0);
	if (nevents <= 0 || (error_ptr && *error_ptr)) usage ();
	break;
      default:
	usage ();
    }
  }
  

  if (optind == argc) gzfile = ::gzdopen (0, "r");;
  if ((optind + 1) == argc) {
    std::string fname = argv[optind++];
    if (fname.substr (0, 7) == "http://") {
      std::string command = std::string ("wget -nv --timeout=60 -O - ") + fname;
      pfile = ::popen(command.c_str () , "r");
      gzfile = ::gzdopen (::fileno (pfile), "r");
    } else gzfile = ::gzopen (fname.c_str (), "r");
  }
  if (optind < argc) usage ();

  while (1) {
    char *buffer = 0;
    int buffer_size = 0;
    if (nevents > 0 && !--nevents) break;
    unsigned long size;
    if (::gzread (gzfile, &size, sizeof (unsigned long)) !=  sizeof (unsigned long)) break;
    if (size > buffer_size) {
      buffer_size = int(size * 1.5);
      buffer = new char[buffer_size];
    }
    *(unsigned long*)buffer = size;
    if (::gzread (gzfile, buffer + sizeof (unsigned long), size - sizeof (unsigned long)) != size - sizeof (unsigned long)) break;
    if (jump > 0 && --jump) continue;
    nanosleep (&sleep, NULL);
    std::cout.write (buffer, size);
  }

  return 0;
}


      

      






