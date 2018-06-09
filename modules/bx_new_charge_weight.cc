/* BOREXINO Reconstruction program
 *
 * Author: Koun Choi <kounchoi@hawaii.edu>
 *
 * $Id: bx_new_charge_weight.cc,v 1.5 2015/08/18 09:37:59 koun Exp $
 *
 */
#include "bx_new_charge_weight.hh"
#include "bx_echidna_event.hh"
#include "db_run.hh"
#include "bx_dbi.hh"
#include "constants.hh"
#include "db_channel.hh"
#include <string>
#include <iostream>
#include <math.h>
#define pi 3.141592653589793238

// ctor
bx_new_charge_weight::bx_new_charge_weight (): bx_base_module("bx_new_charge_weight", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
}

// module interface
void bx_new_charge_weight::begin () {

  nlg = constants::laben::channels;

  // map of disabled channels from bx_detector
  p_disabled_pmts_lg = new bool[nlg + 1];
  p_disabled_charge_lg = new bool[nlg + 1];
  std::fill_n (p_disabled_pmts_lg, nlg + 1, false); // Init vector at 0
  std::fill_n (p_disabled_charge_lg, nlg + 1, false); // Init vector at 0

  p_disabled_lg = new unsigned char[nlg + 1];
  memset (p_disabled_lg, 0, nlg + 1); // Init vector at 0
  ch_info_v = new const db_channel_laben*[nlg + 1];
  for (int i = 1; i <= nlg; i++) 
    ch_info_v[i] = &dynamic_cast<const db_channel_laben&>(bx_dbi::get ()->get_channel (i));
}

bx_echidna_event* bx_new_charge_weight::doit (bx_echidna_event *ev) {

  std::fill_n (p_disabled_pmts_lg, nlg + 1, false); // Init vector at 0
  std::fill_n (p_disabled_charge_lg, nlg + 1, false); // Init vector at 0

  //vector of disabled lg
    const std::vector<int>& v = detector_interface::get ()->get_disabled_channels ();
    for (int i = 0; i < (int) v.size (); i++) if (v[i] <= nlg) {
	if ((p_disabled_lg[v[i]] == db_run::timing)) continue;
	if (ch_info_v[v[i]]->is_ordinary ()){
    int disabled_pmts_lg = v[i];
    p_disabled_pmts_lg[disabled_pmts_lg] = true;
    int disabled_charge_lg = v[i];
    p_disabled_charge_lg[disabled_charge_lg] = true;}
  }
    const std::vector<int>& vc = detector_interface::get ()->get_disabled_charge ();
    for (int i = 0; i < (int) vc.size (); i++) if (vc[i] <= nlg) {
      if (p_disabled_lg[vc[i]] == db_run::timing || p_disabled_lg[vc[i]] == db_run::charge) continue;
	if (ch_info_v[vc[i]]->is_ordinary ()){
    int disabled_charge_lg = vc[i];
    p_disabled_charge_lg[disabled_charge_lg] = true;}
  }
    i4_n_disabled_charge = vc.size ();
    i4_n_disabled_channels = v.size ();


  bx_laben_event& er = ev->get_laben ();
  for (int i = 0; i < er.get_raw_nhits (); i++) {
    const bx_laben_raw_hit &hit = er.get_raw_hit (i);
    int lg = hit.get_logical_channel ();
    const db_channel_laben *ch_info = ch_info_v[lg];
      // Ignore laben invalid hits (0xffff)
         if (hit.check_flag (bx_laben_raw_hit::invalid) && ch_info->is_ordinary ()) {
        unsigned j = 0;
        if (p_disabled_lg[lg] != db_run::charge) {
          while (j <= i4_n_disabled_channels){
        if (lg == v[j]) break;
        j++;}
        if (j > i4_n_disabled_channels){
    p_disabled_pmts_lg[lg] = true;
    p_disabled_charge_lg[lg] = true;}
        }
        if (p_disabled_lg[lg] == db_run::charge) {
          while (j <= i4_n_disabled_charge){
        if (lg == vc[j]) break;
        j++;}
        if (j > i4_n_disabled_charge){
    p_disabled_charge_lg[lg] = true;}
        }
      continue;
    }//if flag
}//for raw_hit

  Int_t size_clusters = ev->get_laben ().get_nclusters ();
  Int_t size_mclusters = ev->get_laben ().get_nclusters_muons ();
  for (int i = 0; i < size_clusters + size_mclusters; i++) {
    bx_laben_cluster& cluster = (i < size_clusters ) ? ev->get_laben().get_cluster(i) : ev->get_laben().get_cluster_muon(i-size_clusters);
    bx_position_lngs& pos = cluster.get_position_lngs ();
        if (!ev->get_laben ().check_stage (bx_base_event::reconstructed_lngs)) {
        get_message (bx_message::error) << "lngs position not available for event " << ev->get_event_number () << " cluster " << i << dispatch;
        break;}

    float x=pos.get_x();
    float y=pos.get_y();
    float z=pos.get_z();
    float r=sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    R_pmt = 0.203;
    d_pmt = 6.85 - 0.3333 + 0.11 - 0.0424;
    angt_pmts = 0;
    QEt_pmts = 0;
    allt_pmts = 0;
    angt_charge = 0;
    QEt_charge = 0;
    allt_charge = 0;

    for (int ilg = 1;  ilg < (1 + nlg); ilg++) {
	QE = 0.;
        const db_channel_laben& channel_info = dynamic_cast<const db_channel_laben&>(bx_dbi::get()->get_channel(ilg));
	xp = channel_info.pmt_x();
	yp = channel_info.pmt_y();
	zp = channel_info.pmt_z();
        float rp=sqrt(pow(xp,2)+pow(yp,2)+pow(zp,2));
        float theta = acos(zp/rp);
        float phi = atan2(xp,yp);
	xp = d_pmt*sin(theta)*cos(phi);
	yp = d_pmt*sin(theta)*sin(phi);
	zp = d_pmt*cos(theta);
        float d = sqrt(pow((xp-x),2) + pow((yp-y),2) + pow((zp-z),2));
	float ang = pi*pow(R_pmt,2)/pow(d,2)*(xp*(xp-x)+yp*(yp-y)+zp*(zp-z))/sqrt(pow(xp,2)+pow(yp,2)+pow(zp,2))/sqrt(pow((x-xp),2)+pow((y-yp),2)+pow((z-zp),2));
        if(r > d_pmt*0.9) ang = 0.00298623;
	//if(channel_info.is_ordinary () && !p_disabled_charge_lg[ilg] && !p_disabled_pmts_lg[ilg] && !bx_dbi::get ()->get_run ().get_laben_qe_nocorr (ilg))
        //get_message (bx_message::error) << "QE not available for channel " << ilg << " event " << ev->get_event_number () << " cluster " << i << dispatch;
        QE = bx_dbi::get ()->get_run ().get_laben_qe_nocorr (ilg);
if (channel_info.is_ordinary () && !p_disabled_pmts_lg[ilg]){
        angt_pmts += ang;
	QEt_pmts += QE;
	allt_pmts += ang*QE;
}
if (channel_info.is_ordinary () && !p_disabled_charge_lg[ilg]){
        angt_charge += ang;
        QEt_charge += QE;
        allt_charge += ang*QE;
}
} // end of lg loop
cluster.npmt_geo_weight = angt_pmts/0.00298623; //solid angle for a pmt seen at the center of detector
cluster.npmt_QE_weight = QEt_pmts;
cluster.npmt_geo_QE_weight = allt_pmts/0.00298623;
cluster.charge_geo_weight = angt_charge/0.00298623;
cluster.charge_QE_weight = QEt_charge;
cluster.charge_geo_QE_weight = allt_charge/0.00298623;
} //end of cluster loop
  return ev;
}

void bx_new_charge_weight::end () {

  delete [] p_disabled_lg;
  delete [] p_disabled_pmts_lg;
  delete [] p_disabled_charge_lg;
  delete [] ch_info_v;
}

/*  
 *  $Log: bx_new_charge_weight.cc,v $
 *  Revision 1.5  2015/08/18 09:37:59  koun
 *  *** empty log message ***
 *
 *  Revision 1.4  2015/08/17 19:49:45  koun
 *  Bug correction
 *
 *  Revision 1.3  2015/08/17 17:21:02  misiaszek
 *  Id & log keywords added
 *
 */
