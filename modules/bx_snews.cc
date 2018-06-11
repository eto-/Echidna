/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * 	   Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Davide D'Angelo <Davide.Dangelo@lngs.infn.it>
 *
 * $Id: bx_snews.cc,v 1.2 2009/06/26 22:15:34 razeto Exp $
 *
 * Implementation of bx_snews
 * If you cut and paste from this file, 
 * remember to remove comments
 * 
 * A few lines of sample code are commented out
 * to allow clean compilation.
 * C-style comments are used in this cases, 
 * while C++style are used for explanations
 */
#include "bx_snews.hh"
#include "bx_echidna_event.hh"
#include <string>
#include <iostream>

///////////////METHODS TO TALK TO THE BRAIN//////////////////////
#ifndef NOXML
#include <XmlRpcCpp.h> 
namespace {
  //global variables for XLM-RPC.
  std::string server_address = "http://snews.daqnet.borex:8080/RPC2";
  std::string my_name = "Echidna";
  int32_t key = 0;

  void BrainInit () {
    XmlRpcClient::Initialize("online echidna snews client", "1");
  }
  void BrainEnd () {
    XmlRpcClient::Terminate();
  }

  int32_t BrainUnregister () {
    // Build our parameter array.
    XmlRpcValue param_array = XmlRpcValue::makeArray();
    param_array.arrayAppendItem( XmlRpcValue::makeString( my_name ) );
    param_array.arrayAppendItem( XmlRpcValue::makeInt( key ) );

    // Create an object to resprent the server, and make our call.
    XmlRpcClient server ( server_address );
    XmlRpcValue result = server.call("snews.unregister", param_array);

    XmlRpcValue::int32 res = result.getInt();

    return res;
  }

  int32_t BrainRegister () {
    // Build our parameter array.
    XmlRpcValue param_array = XmlRpcValue::makeArray();
    param_array.arrayAppendItem( XmlRpcValue::makeString( my_name ) );

    // Create an object to resprent the server, and make our call.
    XmlRpcClient server ( server_address );
    XmlRpcValue result = server.call("snews.register", param_array);

    XmlRpcValue::int32 keys = result.getInt();
    key = keys;

    return keys;
  }

  int32_t BrainTrigger (double alarm_time, double level, double duration, double signal, double probability) {
    // Build our parameter array.
    XmlRpcValue param_array = XmlRpcValue::makeArray();
    param_array.arrayAppendItem( XmlRpcValue::makeString( my_name ) );
    param_array.arrayAppendItem( XmlRpcValue::makeInt( key ) );
    param_array.arrayAppendItem( XmlRpcValue::makeDouble( alarm_time ) );
    param_array.arrayAppendItem( XmlRpcValue::makeDouble( level ) );
    param_array.arrayAppendItem( XmlRpcValue::makeDouble( duration ) );
    param_array.arrayAppendItem( XmlRpcValue::makeDouble( signal ) );
    param_array.arrayAppendItem( XmlRpcValue::makeDouble( probability ) );	

    // Create an object to resprent the server, and make our call.
    XmlRpcClient server ( server_address );
    XmlRpcValue result = server.call("snews.trigger", param_array);

    XmlRpcValue::int32 status = result.getInt();

    return status;
  }

  int32_t BrainLinger () {
    // Build our parameter array.
    XmlRpcValue param_array = XmlRpcValue::makeArray();
    param_array.arrayAppendItem( XmlRpcValue::makeString( my_name ) );
    param_array.arrayAppendItem( XmlRpcValue::makeInt( key ) );

    // Create an object to resprent the server, and make our call.
    XmlRpcClient server ( server_address );
    XmlRpcValue result = server.call("snews.linger", param_array);

    XmlRpcValue::int32 status = result.getInt();

    return status;
  }
};
#else
namespace {
  void BrainInit () { }
  void BrainEnd () { }
  int32_t BrainUnregister () { return true; }
  int32_t BrainRegister () { return true; }
  int32_t BrainTrigger (double alarm_time, double level, double duration, double signal, double probability) { return true; }
  int32_t BrainLinger () { return true; }
};
#endif
/////////////////////////////////////////////////////////////////////


// ctor
bx_snews::bx_snews (): bx_base_module("bx_snews", bx_base_module::main_loop) {
  // Option: you can require to run only on events which have reached a given stage (for a given detector segment)
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
  // Option: you can require to run only on events with a given trigger type
  require_trigger_type (bx_trigger_event::neutrino);
}


// module interface
void bx_snews::begin () {
  registered=false;
  last_linger=0;
  connect_brain=get_parameter("connect_brain").get_bool();
  // Sample message
  get_message (bx_message::debug) << "test begin" << dispatch;

// Option: histogram to test how the module work, to be written to file
//  my_histo_check = new TH1F ("my_histo_check", "My useful plot of this and that", 255, 0, 255);
//  barn_interface::get ()->store (barn_interface::file, my_histo_check, this);

  if (connect_brain) BrainInit ();
}


bx_echidna_event* bx_snews::doit (bx_echidna_event *ev) {
  try {
    if (connect_brain) {
      if(!registered) registered = BrainRegister();
      time_t curr_time = time (0);
      if (curr_time - last_linger > 30) { 
	linger=BrainLinger();
	if(linger==-1) registered = false; 
	last_linger=time (0);
      }
    }
    if (ev->get_laben().get_nclusters() != 1) return ev;
    if (ev->get_laben().get_cluster(0).get_clustered_nhits() <= 250) return ev;
    if (ev->get_trigger ().get_trgtype () > 2) return ev;

    if ((ev->get_trigger ().get_btb_inputs () == 4 
	|| ev->get_muon ().has_cluster ()) 
	&&(ev->get_laben ().get_cluster(0).get_clustered_nhits() >= 2100 
	&& ev->get_laben().get_cluster(0).get_mean_time () > 100) 
	&& (ev->get_laben ().get_cluster(0).get_clustered_nhits() < 2100 
	&& ev->get_laben ().get_cluster(0).get_clustered_nhits() >= 900 
	&& ev->get_laben().get_cluster(0).get_split_npeaks ()> 0 
	&& ev->get_laben().get_cluster(0).get_split_peak(0).get_start_time () > 30) 
	&& (ev->get_laben ().get_cluster(0).get_clustered_nhits() < 900 
	&& ev->get_laben().get_cluster(0).get_split_npeaks () > 0 
	&& ev->get_laben().get_cluster(0).get_split_peak(0).get_start_time () > 40)) {
      EventData CurEv;
      uint32_t Sec,NSec;
      ev->get_trigger().get_gps_time(Sec,NSec);
      CurEv.Time = Sec + double(NSec) * 10e-9;
      CurEv.NHits = ev->get_laben().get_cluster(0).get_clustered_nhits();
      CurEv.ID = ev->get_event_number();
      Train.push_back(CurEv);
      while(Train.size () > 1) {
	double dt = CurEv.Time - Train.front().Time;
	if (dt<=30) break;
	Train.pop_front();
      }
      if (Train.size () >= 15) {
	double allarm_time = Train.front().Time;
	double probability = 0;
	double level=0;
	double duration=CurEv.Time-Train.front().Time;
	double signal=Train.size();     
	if (connect_brain) allarmreturn=BrainTrigger(allarm_time,level,duration,
	    signal,probability);
	//get_message (bx_message::warn) << "lunghezza catena = " << Train.size () << "ID primo evento = " <<
	//Train.front().ID << dispatch;
      }
    }
#ifndef NOXML
  } catch (XmlRpcFault &f) {
    std::cerr << "XMLRPC excpetion: " << f.getFaultString () << std::endl;
    throw;
  }
#else
  } catch (...) { throw; }
#endif
  return ev;
}

void bx_snews::end () {
  if (connect_brain) {
    BrainUnregister();
    BrainEnd ();
  }
  get_message (bx_message::debug) << "test end" << dispatch;
}


/*
 * $Log: bx_snews.cc,v $
 * Revision 1.2  2009/06/26 22:15:34  razeto
 * May disable XML-rpc with environment variable NOXML
 *
 * Revision 1.1  2008-10-23 09:13:46  dicienzo
 * Added snews monitor
 *
 *
 */
