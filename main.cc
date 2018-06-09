#include "bx_reco_framework.hh"
#include "bx_options.hh"
#include "messenger.hh"
#include "barn_interface.hh"

#include <iostream>
#include <stdexcept>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>


int main(int argc, char *argv[]) {
  time_t start_time = ::time (0);
  try {	
    barn_interface::disable_root_signals ();

    bx_options::initialize_echidna (argc, argv);

    bx_reco_framework frame;
    frame.run();
    bx_message msg(bx_message::info);
    if (messenger::get ()->get_error_count () || messenger::get ()->get_warn_count ())
      msg << "Echidna ENDED with " << messenger::get ()->get_error_count () << " errors and " << messenger::get ()->get_warn_count () << " warns" << dispatch;
    else
      msg << "Echidna NORMAL END" << dispatch;
  } catch (std::exception& ex) {
    std::cout << ">>> exception: " << ex.what() << std::endl;
  } catch (...) {
    std::cout << ">>> exception: unknown type" << std::endl;
  }
  time_t end_time = ::time (0);
  struct rusage usage_data;
  if (!::getrusage (RUSAGE_SELF, &usage_data)) {
    bx_message msg(bx_message::info);
    msg << "Timing statistics: " << newline;
    msg << "time real " << ::difftime (end_time, start_time) << "s" << newline;
    msg << "time user " << usage_data.ru_utime.tv_sec << "s" << newline;
    msg << "time sys  " << usage_data.ru_stime.tv_sec << "s" << dispatch;
  }
  messenger::get ()->flush ();
  bx_options::finish ();
}
