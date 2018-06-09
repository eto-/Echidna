/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_laben_decoder.cc,v 1.87 2017/05/12 14:15:45 misiaszek Exp $
 *
 * Implementation of bx_laben_decoder
 *
 */
#include "bx_laben_decoder.hh"
#include "messenger.hh"
#include "laben_time_hit.hh"
#include "laben_charge_hit.hh"
#include "constants.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "db_channel.hh"
#include "cmap.hh"
#include "barn_interface.hh"

#include <map>
#include <algorithm>


// ctor
bx_laben_decoder::bx_laben_decoder (): bx_base_module("bx_laben_decoder", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::raw);
}

// module interface
void bx_laben_decoder::begin () {
  get_message(bx_message::debug) << "begin" << dispatch;

  b_discard_out_of_gate_hits = get_parameter ("discard_out_of_gate_hits").get_bool ();
  b_discard_reference_hits = get_parameter ("discard_reference_hits").get_bool ();
  b_discard_retrigger_hits = get_parameter ("discard_retrigger_hits").get_bool ();
  b_discard_reflection_hits = get_parameter ("discard_reflection_hits").get_bool ();
  b_discard_disabled_lg = get_parameter ("discard_disabled_lg").get_bool ();
  b_discard_calib_data = get_parameter ("discard_calib_data").get_bool ();
  b_warn_empy_channel = get_parameter ("warn_empy_channel").get_bool ();
  b_mc_use_charge_calib_data = get_parameter ("mc_use_charge_calib_data").get_bool ();
  i4_nhits_threshold = get_parameter ("nhits_threshold").get_float ();
  f4_shift_min = get_parameter ("shift_min").get_float ();
  f4_shift_max = get_parameter ("shift_max").get_float ();

  trigger_charge = new TH1F ("trigger_charge", "Trigger hit charge values", 255, 0, 255);
  laser_charge = new TH1F ("laser_charge", "Laser hit charge values", 255, 0, 255);
  gray_cross_h = new TH2F ("gray_cross", "Gray cross population", 50, 0, 1, 350, 0, 3500000);
  barn_interface::get ()->store (barn_interface::file, trigger_charge, this);
  barn_interface::get ()->store (barn_interface::file, laser_charge, this);
  barn_interface::get ()->store (barn_interface::file, gray_cross_h, this);

    // channels properties (disabled, ch_info, empty)
  p_disabled_lg = new unsigned char[constants::laben::channels + 1];
  memset (p_disabled_lg, 0, constants::laben::channels + 1); // Init vector at 0
  i4_ordinary_pmt = i4_n_disabled_channels = i4_n_disabled_charge = i4_n_disabled_pmts = i4_n_disabled_pmts_charge = 0;
  ch_info_v = new const db_channel_laben*[constants::laben::channels + 1];
  p_empty_lg = new bool[constants::laben::channels + 1];
  for (unsigned i = 1; i <= constants::laben::channels; i++) {
    ch_info_v[i] = &dynamic_cast<const db_channel_laben&>(bx_dbi::get ()->get_channel (i));
    p_empty_lg[i] = ch_info_v[i]->is_empty () || ch_info_v[i]->is_pmt_disconnected ();
    if (ch_info_v[i]->is_ordinary ()) i4_ordinary_pmt++;
  }
}

namespace {
  struct in_gate_checker {
    in_gate_checker (double low_bound, double high_bound): f8_low_bound(low_bound), f8_high_bound(high_bound) {}
    bool operator() (const bx_laben_decoded_hit& hit) { 
      return ((hit.get_raw_time () > f8_low_bound) && (hit.get_raw_time () < f8_high_bound)); 
    }
    double f8_low_bound, f8_high_bound;
  };
};

bx_echidna_event* bx_laben_decoder::doit (bx_echidna_event *ev) {
  detector_interface::get ()->read_disabled_channels (ev->get_event_number ());
  const std::vector<int>& v = detector_interface::get ()->get_disabled_channels ();
  const std::vector<int>& vc = detector_interface::get ()->get_disabled_charge ();
  if (v.size () != i4_n_disabled_channels || vc.size () != i4_n_disabled_charge) {
    i4_n_disabled_pmts = 0;
    i4_n_disabled_pmts_charge = 0;
    memset (p_disabled_lg, 0, constants::laben::channels + 1); // Init vector at 0
    for (unsigned i = 0; i < v.size (); i++) if (v[i] <= constants::laben::channels) {
      if (p_disabled_lg[v[i]] == db_run::timing) continue;
      p_disabled_lg[v[i]] = db_run::timing;
      if (ch_info_v[v[i]]->is_ordinary ()) i4_n_disabled_pmts++;
    }
    for (unsigned i = 0; i < vc.size (); i++) if (vc[i] <= constants::laben::channels) {
      if (p_disabled_lg[vc[i]] == db_run::timing || p_disabled_lg[vc[i]] == db_run::charge) continue;
      p_disabled_lg[vc[i]] = db_run::charge;
      if (ch_info_v[vc[i]]->is_ordinary ()) i4_n_disabled_pmts_charge++;
    }
    i4_n_disabled_charge = vc.size ();
    i4_n_disabled_channels = v.size ();
  }

  bx_laben_event& er = ev->get_laben ();
  bx_laben_decoded_event& e = dynamic_cast<bx_laben_decoded_event&>(er);

    // Set number of alive pmt
  e.i4_n_live_pmts = i4_ordinary_pmt - i4_n_disabled_pmts;
  if (e.i4_n_live_pmts < 0) get_message (bx_message::critic) << "e.i4_n_live_pmts can never be negative " << e.i4_n_live_pmts << dispatch;
  e.i4_n_live_charge = i4_ordinary_pmt - i4_n_disabled_pmts - i4_n_disabled_pmts_charge;
  e.i4_nhits_on_empty = 0;

    // Check if to apply calib data
  bool use_calib_data = true;
  if (b_discard_calib_data || ev->get_trigger ().is_pulser () || ev->is_mctruth_enabled ()) use_calib_data = false;
  bool use_charge_calib_data = (ev->is_mctruth_enabled ()) ? (b_mc_use_charge_calib_data) : use_calib_data;

    // Create an empty laben time hit and an empty laben charge hit
  laben_time_hit t_hit(use_calib_data);
  laben_charge_hit c_hit(use_charge_calib_data);

    // Check for gray counter cross in the trigger gate.
  float gray_cross_ratio = m_check_gray_cross (er);
  bool gray_cross = gray_cross_ratio > 0.7;

    // 1) Create the matrix and an iterator
  std::map <int, bx_laben_decoded_hit_list> channel_hits_map;
  std::map <int, bx_laben_decoded_hit_list>::iterator it;

    // 2) Fill the matrix (even fill the time since the laben_time_hit is create to check if the hit is valid)
  int invalid_on_good = 0, invalid_on_charge = 0;
  for (int i = 0; i < er.get_raw_nhits (); i++) {
    const bx_laben_raw_hit &hit = er.get_raw_hit (i);
    int lg = hit.get_logical_channel ();
    const db_channel_laben *ch_info = ch_info_v[lg];

      // Ignore laben invalid hits (0xffff)
    if (hit.check_flag (bx_laben_raw_hit::invalid)) {
      if (ch_info->is_ordinary ()) {
        unsigned j = 0;
        if (p_disabled_lg[lg] != db_run::charge) {
	if(v.size()){
          while (j <= i4_n_disabled_channels){
        if (lg == v[j]) break;
        j++;}
        if (j > i4_n_disabled_channels){
          invalid_on_good ++;
          invalid_on_charge ++;}
        }
	else{
	          invalid_on_good ++;
          invalid_on_charge ++;}}
        if (p_disabled_lg[lg] == db_run::charge){
          if(vc.size()){
	while (j <= i4_n_disabled_charge){
        if (lg == vc[j]) break;
        j++;}
        if (j > i4_n_disabled_charge)
          invalid_on_charge ++;}
	else
	 invalid_on_charge ++;
        }
      }
      continue;
    } 

      // Ignore counter hits (new fw)
    if (hit.check_flag (bx_laben_raw_hit::counter)) continue;

      // Counts hits on empty
    if (p_empty_lg[lg]) e.i4_nhits_on_empty ++;

      // Ignore disabled channels
    if (b_discard_disabled_lg && p_disabled_lg[lg] == db_run::timing) continue;

      // Init the laben time hit
    t_hit.init (hit, gray_cross);

      // Skip invalid hits (laben_time_hit can be invalid even if hit is valid)
    if (!t_hit.is_valid ()) continue;

      // Create a decoded hit
    bx_laben_decoded_hit decoded_hit(hit, i);

      // Associate the db_channel
    decoded_hit.p_db_channel = ch_info;

      // Fill the time data (charge will be filled later)
    decoded_hit = t_hit;

      // add state
    if (p_disabled_lg[lg] == db_run::timing) decoded_hit.u1_flag |= bx_laben_decoded_hit::disabled;

      // Push the decoded hit in the column
    channel_hits_map[lg].push_back (decoded_hit);
  }


    // 3) Process several time the matrix looping on every column
    // 3a) Fill the charge data
  for (it = channel_hits_map.begin (); it != channel_hits_map.end (); it++) {
    if (p_disabled_lg[it->first] == db_run::charge) {
      for (bx_laben_decoded_hit_list::iterator jt = it->second.begin (); jt != it->second.end (); jt++) 
      { jt->f4_charge_bin = jt->f4_uncorrected_charge_bin = jt->i4_charge_npe = 0; jt->f4_charge_pe = jt->f4_uncorrected_charge_pe = 0; }
      continue;
    }

    bx_laben_decoded_hit_list &c_list = it->second;

      // Sort in time the column vector (all hits have the same channel)
    c_list.sort ();

      // Fill the charge data: the first is special since ho previus hits 
      // can be overlapped, so use the simple laben_charge_hit ctor.
    bx_laben_decoded_hit& hit = *(c_list.begin ());
    c_hit.init (hit);
    hit = c_hit;
      // If other elements are present use the laben_charge_hit with 2 raw events
    bx_laben_decoded_hit_list::iterator it_prev = c_list.begin ();
    bx_laben_decoded_hit_list::iterator it_curr = ++(c_list.begin ());
    while (it_curr != c_list.end ()) {
      c_hit.init (*it_curr, *it_prev);
      *it_curr = c_hit;
      it_curr++; it_prev++;
    }
  }

    // 3ab) Remove retriggering hits
  for (it = channel_hits_map.begin (); it != channel_hits_map.end (); it++) {
    remove_retrigger_hits (it->second);
    remove_cable_reflection_hits (it->second);
  }

    // 3b) Fill internal vector used by m_get_moda_mean
  std::vector<double> trigger_times, laser_times;
  for (it = channel_hits_map.begin (); it != channel_hits_map.end (); it++) {
    bx_laben_decoded_hit_list &c_list = it->second;
    const db_channel* ch_info = c_list.begin ()->get_db_channel ();
    if (ch_info->is_trigger ()) {
      const bx_laben_decoded_hit &hit = c_list.back ();
      if (ev->get_run_number () > 2070) {
	if (hit.get_uncorrected_charge_bin () > 70 || ev->get_trigger ().get_btb_inputs () == 4) {
	  trigger_times.push_back (hit.get_raw_time ());
	  trigger_charge->Fill (hit.get_uncorrected_charge_bin ());
	}
      } else {
	trigger_times.push_back (hit.get_raw_time ());
	trigger_charge->Fill (hit.get_uncorrected_charge_bin ());
      }
    } else if (ch_info->is_laser ()) {
      for (bx_laben_decoded_hit_list::iterator it_curr = c_list.begin (); it_curr != c_list.end (); it_curr++) {
	laser_charge->Fill (it_curr->get_uncorrected_charge_bin ());
	if (it_curr->get_uncorrected_charge_bin () > 65) laser_times.push_back (it_curr->get_raw_time ());
      }
    }
  }
    // And calculate the times
  if (!trigger_times.size ()) {
    get_message (bx_message::error) << "no trigger hits in event " << ev->get_event_number () << " of type " << int(ev->get_trigger ().get_trgtype ()) << dispatch;
    return ev;
  }
  if (ev->get_trigger ().is_laser () && !laser_times.size ()) {
    get_message (bx_message::error) << "no laser hits in laser event " << ev->get_event_number () << dispatch;
    return ev;
  }
  e.f8_trigger_rawt = m_get_moda_mean (trigger_times);
  e.f8_laser_rawt = m_get_moda_mean (laser_times);

    // 3bb Handle gray counter cyclicity
  static const double gray_range = (1UL << 16) * 50;
  static const double end_of_gate_limit = 3e3;
  if (e.f8_trigger_rawt < (gray_range - end_of_gate_limit)) { // rearrange hits time for events in the lower gray counter part 1.8 over 3.2 ms
    double bound = e.f8_trigger_rawt + end_of_gate_limit; // Every hit later than 2us is from previus gray cycle
    for (it = channel_hits_map.begin (); it != channel_hits_map.end (); it++) {
      bx_laben_decoded_hit_list &c_list = it->second;
      for (bx_laben_decoded_hit_list::iterator it_curr = c_list.begin (); it_curr != c_list.end (); it_curr++) {
	if (it_curr->f8_raw_time < bound) it_curr->f8_raw_time += gray_range;
      }
    }
    e.f8_trigger_rawt += gray_range;
    e.f8_laser_rawt += gray_range;
  }

    // Fill histogram of gray_cross probability vs trigger_time
  gray_cross_h->Fill (gray_cross_ratio, e.f8_trigger_rawt);

  if (ev->get_trigger ().is_laser () && laser_times.size () > 0 && ::fabs(e.f8_laser_rawt - e.f8_trigger_rawt) > 2e4) {  // the laser will be out of gate
    get_message (bx_message::warn) << "event " << ev->get_event_number () << " laser_time - trigger_time too distant (dt " << e.f8_laser_rawt-e.f8_trigger_rawt << " ns) laser_time is " << e.f8_laser_rawt << " ns and trg_time is " << e.f8_trigger_rawt << " ns with gray_cross_ratio " << gray_cross_ratio << dispatch;
    return 0;
  }

    // 3c) Remove out of gate hits
    // and assign order_in_channel
  double low_bound = e.get_trigger_rawt () - 30000;
  double high_bound = e.get_trigger_rawt () + end_of_gate_limit;
  if (ev->get_trigger ().is_neutron ()) low_bound = e.get_trigger_rawt () - 1.6e6;
  if (low_bound < 0) {
    get_message (bx_message::warn) << "internal error low_bound negative (" << low_bound << ") or event " << ev->get_event_number ()<<  dispatch;
    return 0;
  }
  for (it = channel_hits_map.begin (); it != channel_hits_map.end (); it++) {
    bx_laben_decoded_hit_list &c_list = it->second;
    bx_laben_decoded_hit_list::iterator middle;
    middle = std::stable_partition (c_list.begin (), c_list.end (), in_gate_checker (low_bound, high_bound));
    if (b_discard_out_of_gate_hits) c_list.erase (middle, c_list.end ());
    else for (; middle != c_list.end (); middle++) middle->u1_flag |= bx_laben_decoded_hit::out_of_gate;
    int i = 0;
    for (bx_laben_decoded_hit_list::iterator it_curr = c_list.begin (); it_curr != c_list.end (); it_curr++) it_curr->u1_order_in_channel = ++i;
  }

    // 3d) Calculate the vector lenghts
  int std_decoded_nhits = 0, trigger_decoded_nhits = 0, laser_decoded_nhits = 0;
  for (it = channel_hits_map.begin (); it != channel_hits_map.end (); it++) {
    bx_laben_decoded_hit_list &c_list = it->second;
    if (c_list.size ()) {
      const db_channel* ch_info = c_list.begin ()->get_db_channel ();
      if (ch_info->is_ordinary ()) std_decoded_nhits += c_list.size ();
      else if (ch_info->is_trigger ()) trigger_decoded_nhits += c_list.size ();
      else if (ch_info->is_laser ()) laser_decoded_nhits += c_list.size ();
    }
  }
    // Ask to the decoded hit vectors to reserve the required space (avoiding too many relocations)
  e.decoded_hits.reserve (std_decoded_nhits);
  e.trigger_reference_decoded_hits.reserve (trigger_decoded_nhits);
  e.laser_reference_decoded_hits.reserve (laser_decoded_nhits);

    // 3e) Finally fill the vectors
    // Now process the matrix looping on every column
  e.i4_npmts = 0;
  for (it = channel_hits_map.begin (); it != channel_hits_map.end (); it++) {
    bx_laben_decoded_hit_list &c_list = it->second;
    if (c_list.size ()) {
      const db_channel* ch_info = c_list.begin ()->get_db_channel ();

        // Add the column to the event decoded_hit list.
      if (ch_info->is_ordinary () || (!b_discard_reference_hits && !ch_info->is_empty ()))
	e.decoded_hits.insert (e.decoded_hits.end (), c_list.begin (), c_list.end ());

        // Fill npmts
      if (ch_info->is_ordinary ()) e.i4_npmts ++;
 
        // Add the rerference hits
      if (ch_info->is_laser ()) 
	e.laser_reference_decoded_hits.insert (e.laser_reference_decoded_hits.end (), c_list.begin (), c_list.end ());
      else if (ch_info->is_trigger ())
	e.trigger_reference_decoded_hits.insert (e.trigger_reference_decoded_hits.end (), c_list.begin (), c_list.end ());
    }
  }
 
  //pup
  if(get_parameter ("enable_pile_up").get_bool ()) {
    for (bx_laben_decoded_hit::bx_laben_decoded_hit_vector::iterator it = e.decoded_hits.begin (); it != e.decoded_hits.end (); it++) {
      double d_hit_time = it->get_raw_time ();
      double ref_time = e.get_trigger_rawt () + f4_trigger_start_;
      double subt_time = d_hit_time - ref_time;
      if(subt_time >= f4_shift_min && subt_time < f4_shift_max)
        it->f8_raw_time = d_hit_time - f4_shift_min;
    }
  }  


    // Sort the list in time  
  std::sort (e.decoded_hits.begin (), e.decoded_hits.end ());

    // calculate charge and npe
  e.f4_npe = e.f4_charge = 0.;
  for (bx_laben_decoded_hit::bx_laben_decoded_hit_vector::const_iterator it = e.decoded_hits.begin (); it != e.decoded_hits.end (); it++) {
    e.f4_npe += it->get_charge_npe ();
    e.f4_charge += it->get_charge_pe ();
  }
  e.i4_invalid_pmts = invalid_on_good;
  e.i4_invalid_charge = invalid_on_charge;
    
    // Put an energy threshold
  if (int(e.decoded_hits.size ()) < i4_nhits_threshold) return ev;

  er.mark_stage (bx_base_event::decoded);
  return ev;
}

void bx_laben_decoder::end () {
  delete [] p_disabled_lg;
  delete [] ch_info_v;

  get_message(bx_message::debug) << "end" << dispatch;
}

  // Collect statistic to sse if we are in a gray cross boundary
float bx_laben_decoder::m_check_gray_cross (const bx_laben_event& er) {
    // Create an empty laben time hit
  laben_time_hit t_hit;
  
  int n_cross = 0;
  int n_valid = 0;
  
  for (int i = 0; i < er.get_raw_nhits (); i++) {
    const bx_laben_raw_hit &hit = er.get_raw_hit (i);
    
       // only check logical channels
    if (bx_dbi::get ()->get_channel (hit.get_logical_channel ()).is_trigger ()) {

	// Init the laben time hit
      t_hit.init (hit);
      
	// Skip invalid hits
      if (!t_hit.is_valid ()) continue;
      
	// Check if the hit is close to a crossing reagion, the increment ncross
      if (t_hit.is_gray_crossing_window ()) n_cross ++;

      n_valid ++;
    }
  }

    // a gray cross window is present if 70% of hits are in the window
  float ratio = float (n_cross) / float (n_valid);
  return ratio;
}

double bx_laben_decoder::m_get_moda_mean (const std::vector<double>& v) {
  if (!v.size ()) return 0;

  double *times = new double[v.size ()];
  for (unsigned i = 0; i < v.size (); i++) times[i] = v[i];
  std::sort (times, times + v.size ());
  
  double moda = m_search_moda (times, v.size ());
  double sum = 0;
  double count = 0;
  for (unsigned i = 0; i < v.size (); i++) {
    if (times[i] < (moda - 5) || times[i] > (moda + 5)) continue;
    sum += times[i];
    count ++;
  }
  
  delete [] times;
  return count ? sum / count : 0.;
}

double bx_laben_decoder::m_search_moda (const double *times, unsigned short size) {
  if (!size) return 0;
  if (size == 1) return times[0];
  
  double min = times[0];
  double max = times[size - 1];
  double dt = max - min;
  if (dt < 10) return (min + max) / 2; // If the are too close use simple interpolation

  unsigned char *histo = new unsigned char[int(dt / 5) + 1];
  std::fill (histo, histo + int(dt / 5) + 1, 0);
  
  for (unsigned i = 0; i < size; i++) {
    int bin = int((times[i] - min) / 5);
    histo[bin]++;
  }

  int moda_bin = std::max_element (histo, histo + int(dt / 5) + 1) - histo;
    
  delete [] histo;
  return moda_bin * 5 + min;
}

void bx_laben_decoder::remove_retrigger_hits (bx_laben_decoded_hit_list &c_list) {
  bx_laben_decoded_hit_list::iterator it_prev = --(c_list.end ());
  bx_laben_decoded_hit_list::iterator it_curr = it_prev--; // note affix and postfix operators
  while (it_curr != c_list.begin ()) {
    if (it_curr->get_raw_time () - it_prev->get_raw_time () < 180) {
      if (b_discard_retrigger_hits) c_list.erase (it_curr);
      else it_curr->u1_flag |= bx_laben_decoded_hit::retrigger;
    }
    it_curr = it_prev--;
  }
}

void bx_laben_decoder::remove_cable_reflection_hits (bx_laben_decoded_hit_list &c_list) {
  bx_laben_decoded_hit_list::iterator it_curr = --(c_list.end ());
  while (it_curr != c_list.begin ()) {
    bx_laben_decoded_hit_list::iterator it_prev = it_curr;
    bx_laben_decoded_hit_list::iterator it_next = it_curr; it_next--;
    while (it_prev != c_list.begin ()) {
      it_prev--;
      if (it_prev->get_charge_pe () > 0 && it_prev->get_charge_pe () < 3) continue; // ignore reflecion for small signals
      if (::fabs (it_curr->get_raw_time () - it_prev->get_raw_time () - 615) < 10) {
	if (b_discard_reflection_hits) c_list.erase (it_curr);
	else it_curr->u1_flag |= bx_laben_decoded_hit::reflection;
	break;
      }
    }
    it_curr = it_next;
  }
}

/*
 * $Log: bx_laben_decoder.cc,v $
 * Revision 1.87  2017/05/12 14:15:45  misiaszek
 * denis patch
 *
 *
 * Revision 1.85  2015/07/31 21:29:30  koun
 * laben invalid hits bug fixed
 *
 * Revision 1.84  2012/04/02 11:39:30  razeto
 * Avoid warning
 *
 * Revision 1.83  2012-03-22 19:19:57  razeto
 * Added nhits_threshold
 *
 * Revision 1.82  2012-03-22 19:08:37  razeto
 * Added cngs reference channels
 *
 * Revision 1.81  2011-07-18 12:13:16  razeto
 * Bugfix of muon peak
 *
 * Revision 1.80  2011-02-28 14:13:23  razeto
 * Fixed memset usage
 *
 * Revision 1.79  2011-02-28 12:32:09  razeto
 * Zero p_disabled_lg
 *
 * Revision 1.78  2011-02-18 16:07:25  razeto
 * n_live_charge improved
 *
 * Revision 1.77  2011-01-19 10:42:44  davini
 * added the possibility to use charge calibration data for channels in montecarlo
 *
 * Revision 1.76  2010-02-10 09:58:38  razeto
 * disabled_charge bug fixed
 *
 * Revision 1.75  2010-01-20 13:07:00  razeto
 * Keep some more hits in past
 *
 * Revision 1.74  2009-10-26 11:19:37  ddangelo
 * trgtype 'muon_tct' (or 'muon_totc') renamed as 'neutron' after real use
 * trgtype 'muon_mtb' renamed as 'muon' as duplicity with the above one disappeared and 'mtb' is associated to btb_inputs flag
 *
 * Revision 1.73  2009-10-26 11:15:49  razeto
 * Removed a warning
 *
 * Revision 1.72  2009-10-09 10:34:35  ddangelo
 * added computation of n_invalid_on_good
 *
 * Revision 1.71  2009-10-07 12:46:10  razeto
 * Fixed n_live_charge calculation (broken earlier): moreover skip timing disabled channels
 *
 * Revision 1.70  2009-10-06 13:28:44  razeto
 * Added invalid_pmts calculation
 *
 * Revision 1.69  2009-10-05 15:26:28  razeto
 * Calculate npmts for decoded event
 *
 * Revision 1.68  2009-07-16 16:02:38  ddangelo
 * removing a warning
 *
 * Revision 1.67  2008-10-07 14:03:56  razeto
 * Using new flag
 *
 * Revision 1.66  2008-09-25 13:11:32  razeto
 * Zeroing charge for charge disabled events
 *
 * Revision 1.65  2008-09-24 17:15:13  razeto
 * Updated detector map for disabled channels
 *
 * Revision 1.64  2008-02-06 11:39:17  razeto
 * Added livia patch for muons
 *
 * Revision 1.63  2008-02-05 22:15:48  razeto
 * Do not trash events when decoding fails, but do not mark the stage
 *
 * Revision 1.62  2007-12-10 15:45:43  ddangelo
 * debugging
 *
 * Revision 1.61  2007-12-10 11:17:56  razeto
 * use only ordinary channels as base for n_live_pmts
 *
 * Revision 1.60  2007-12-07 14:27:07  ddangelo
 * complying with a renamed variable.
 * added filling of n_live_charge variable.
 * (auth by mainatainer)
 *
 * Revision 1.59  2007-11-20 16:18:15  razeto
 * fixed a bug
 *
 * Revision 1.58  2007-11-20 12:27:39  razeto
 * Do not trash empty channels unless requested
 *
 * Revision 1.57  2007-11-12 12:03:21  razeto
 * Empty now is bx_geometry::empty or disconnected
 *
 * Revision 1.56  2007-11-07 15:29:35  ludhova
 * added charge control for trigger reference hits to be > 70 (agreed with Ale)
 *
 * Revision 1.55  2007-11-05 23:39:53  razeto
 * hits_on_empty filled
 *
 * Revision 1.54  2007-10-30 17:47:35  razeto
 * Using new names for flag
 *
 * Revision 1.53  2007-10-30 17:34:26  razeto
 * Experimental code:
 * - better algorithm for finding tigger time
 * - handle better gray window
 * - handle large gate for neutron events
 *
 * Revision 1.52  2007-05-25 16:11:58  razeto
 * Decoded total charge added
 *
 * Revision 1.51  2007-05-08 23:54:23  razeto
 * Live pmt calculation integrated
 *
 * Revision 1.50  2007-05-02 16:38:57  razeto
 * Updated to new bx_detector: now channels can be disabled runtime (to be tested)
 *
 * Revision 1.49  2007/01/24 15:08:09  razeto
 * Fixed a bug in the charge. Use a low threshold for reference hits
 *
 * Revision 1.48  2006/12/05 13:48:51  razeto
 * Fixed printout
 *
 * Revision 1.47  2006/12/05 13:46:54  razeto
 * Print real trigger type
 *
 * Revision 1.46  2006/12/05 13:35:57  razeto
 * More validity checks
 *
 * Revision 1.45  2006/11/28 13:24:45  razeto
 * Initialize to 0 disabled channels vector
 *
 * Revision 1.44  2006/11/27 10:59:01  razeto
 * Upgraded bx_detector to have more detailed check on bad channels multiplicity:
 * now bad_multiplicity do not exist any more but there are other fields.
 * Upgraded laben_decoder: the fix_bad_channel flag do not exists anymore. To
 * avoid bad channels skipping, set detector_interface.bad_laben_channel_properties
 * to {}.
 * Configurations updated.
 *
 * Revision 1.43  2006/11/20 16:42:16  razeto
 * Fixed a typo
 *
 * Revision 1.42  2006/11/18 14:06:35  razeto
 * Upgraded cable reflection handling
 *
 * Revision 1.41  2006/11/16 11:17:18  razeto
 * Fixed a bug
 *
 * Revision 1.40  2006/11/16 10:08:15  razeto
 * Experimental version of decoder with new peaks suppression (to be debugged)
 *
 * Revision 1.39  2006/09/09 14:08:21  razeto
 * Upgraded to remove retrigger hits and added a fix_bad_channels parameter
 *
 * Revision 1.38  2006/08/21 16:01:51  razeto
 * Now bad laser events are skipped
 *
 * Revision 1.37  2006/08/21 11:21:05  razeto
 * Updated to new barn_interface + added disabled channels from bx_detector
 *
 * Revision 1.36  2006/07/17 18:27:59  razeto
 * Upgraded gray cross detection, added quality check on laser events
 *
 * Revision 1.35  2006/07/13 14:27:30  razeto
 * Upgraded gray cross handling (enlarged again window, lowered probability requirement and add debug histo)
 *
 * Revision 1.34  2006/07/12 14:25:53  ludhova
 * added charge cut for laser refernce lg's, corrected a bug in indexing
 *
 * Revision 1.33  2006/06/29 15:03:25  razeto
 * Added indexes to lower level hits
 *
 * Revision 1.32  2006/05/08 17:31:31  razeto
 * Added db_channel patch for fadc (sent to the mailing list)
 *
 * Revision 1.31  2006/04/01 16:25:41  ludhova
 * changes in laser_time
 *
 * Revision 1.30  2006/03/21 15:03:45  razeto
 * Now time resolution is better calculated in bx_calib_laben_decoding
 *
 * Revision 1.29  2006/01/02 21:23:46  razeto
 * Changed test target to file target for db_barn
 *
 * Revision 1.28  2005/11/20 18:12:36  razeto
 * Updated to new laben_*_hit interface
 *
 * Revision 1.27  2005/06/29 13:31:27  razeto
 * Fixed the charge check on reference hits for older runs
 *
 * Revision 1.26  2005/06/27 16:15:02  razeto
 * Added a new histogram to have the detector time resolution witout recalculating the precalibrations
 *
 * Revision 1.25  2005/06/08 10:23:04  razeto
 * Implemented energy cut for reference hits
 *
 * Revision 1.24  2005/05/09 15:00:07  razeto
 * Reintroduced filling of out_of_gate (since it is still usefull). Fixed a small bug
 *
 * Revision 1.23  2005/05/05 16:54:41  razeto
 * Do not use laser calibration on montecarlo data
 *
 * Revision 1.22  2005/03/14 18:27:05  razeto
 * Added order_in_channel filling
 *
 * Revision 1.21  2005/03/11 15:29:48  razeto
 * Added calibration time data usage in the decoder
 *
 * Revision 1.20  2005/03/03 15:16:24  razeto
 * Upgraded
 *
 * Revision 1.19  2005/03/03 09:33:25  razeto
 * Fixed a typo
 *
 * Revision 1.18  2005/03/03 09:30:35  razeto
 * Added a parameter to even copy reference channels hits and
 * fixed the relevant bug.
 * Changed the window for out of gate hits
 *
 * Revision 1.17  2005/03/02 17:34:59  razeto
 * Updated to new laben_time_hit
 * Removed useless lines
 * Upgraded to really remove gateless hits.
 *
 * Revision 1.16  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.15  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.14  2004/10/19 16:24:26  razeto
 * Updated to use vdt::get_bool
 *
 * Revision 1.13  2004/09/22 13:26:19  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.12  2004/09/22 12:21:03  razeto
 * Updated to follow new sub event getters in bx_base_even
 *
 * Revision 1.11  2004/09/22 10:36:08  razeto
 * Updated to follow sub_detector enum in bx_detector
 *
 * Revision 1.10  2004/09/13 10:11:31  razeto
 * Removed some debug instructions and fixed some typos
 *
 * Revision 1.9  2004/09/09 11:53:18  razeto
 * Doit method mainly rewritten; now several loops on the hit matrix
 * are done now to allow a better (and hopefully cleaner) decoding.
 * New decoded event format introduced (to fully have time reference
 * hits).
 * Decoder is now almost ready to recognise time reference hits from
 * predecoding hits using the charge samples (for runid > 1800).
 *
 * Revision 1.8  2004/07/14 11:49:00  razeto
 * Fixed a bug (array overflow) and upgraded some routines
 *
 * Revision 1.7  2004/06/22 16:04:12  razeto
 * Fixed a bug
 *
 * Revision 1.6  2004/06/08 10:25:19  razeto
 * Updated to new hit name (rawt -> raw_time)
 *
 * Revision 1.5  2004/06/07 09:39:41  razeto
 * Removed the time_sort functor, using decoded hit operator<.
 * Added a global time sorting for the decoded hits.
 *
 * Revision 1.4  2004/06/01 15:44:48  razeto
 * Added event stage marking
 *
 * Revision 1.3  2004/06/01 15:38:25  razeto
 * A lot of development to support time references, out_of_gate hits and much more
 *
 * Revision 1.2  2004/05/21 11:28:23  razeto
 * Updated conforming to the new conventions of event::get_logical_channel
 *
 * Revision 1.1  2004/05/21 08:40:33  razeto
 * Added
 *
 */
