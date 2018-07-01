/* BOREXINO Reconstruction program
 * 
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: BxEventLinkDef.hh,v 1.14 2011/02/18 17:10:05 ddangelo Exp $
 *
 * CINT LinkDef for event 
 *
 */

#ifndef _BXEVENTLINKDEF_HH
#define _BXEVENTLINKDEF_HH
#if defined(__ROOTCLING__) || defined(__CINT__)
#define _ROOT_CINT_
#endif

#ifdef _ROOT_CINT_

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class BxTrigger;

#pragma link C++ class BxTrack;
#pragma link C++ class BxTrackByPoints;
#pragma link C++ class BxTrackFitted;
#pragma link C++ class BxDistance;

#pragma link C++ class BxLabenRawHit;
#pragma link C++ class std::vector<BxLabenRawHit>;
#pragma link C++ class BxLabenDecodedHit;
#pragma link C++ class std::vector<BxLabenDecodedHit>;
#pragma link C++ class BxLabenClusteredHit;
#pragma link C++ class std::vector<BxLabenClusteredHit>;
#pragma link C++ class BxLabenCluster;
#pragma link C++ class std::vector<BxLabenCluster>;
#pragma link C++ class BxLabenRecHit;
#pragma link C++ class std::vector<BxLabenRecHit>;
#pragma link C++ class BxLabenRecCluster;
#pragma link C++ class std::vector<BxLabenRecCluster>;
#pragma link C++ class BxLaben;
#pragma link C++ class BxPosition;
#pragma link C++ class BxEnergy;

#pragma link C++ class BxMuonRawHit;
#pragma link C++ class std::vector<BxMuonRawHit>;
#pragma link C++ class BxMuonDecodedHit;
#pragma link C++ class std::vector<BxMuonDecodedHit>;
#pragma link C++ class BxMuonClusteredHit;
#pragma link C++ class std::vector<BxMuonClusteredHit>;
#pragma link C++ class BxMuonCluster;
#pragma link C++ class std::vector<BxMuonCluster>;
#pragma link C++ class BxMuon;

#pragma link C++ class BxMcTruthHit;
#pragma link C++ class std::vector<BxMcTruthHit>;
#pragma link C++ class BxMcTruthDaughter;
#pragma link C++ class std::vector<BxMcTruthDaughter>;
#pragma link C++ class BxMcTruthDeposit;
#pragma link C++ class std::vector<BxMcTruthDeposit>;
#pragma link C++ class BxMcTruthUser;
#pragma link C++ class std::vector<BxMcTruthUser>;
#pragma link C++ class BxMcTruthFrame;
#pragma link C++ class std::vector<BxMcTruthFrame>;
#pragma link C++ class BxMcTruth;

#pragma link C++ class BxNeutronPulse;
#pragma link C++ class std::vector<BxNeutronPulse>;
#pragma link C++ class BxNeutron;

#pragma link C++ class BxPhysTags;
#pragma link C++ class std::vector<BxPhysTags>;

#pragma link C++ class BxFwfdCluster;
#pragma link C++ class std::vector<BxFwfdCluster>;
#pragma link C++ class BxFwfd;

#pragma link C++ class BxEvent;

#endif
#endif
/*
 * $Log: BxEventLinkDef.hh,v $
 * Revision 1.14  2011/02/18 17:10:05  ddangelo
 * major code cleanup: removed fadc throughout the program
 *
 * Revision 1.13  2010-06-16 12:54:20  litvinov
 * added FWFD stuff. Commit authorized by ddangelo
 *
 * Revision 1.12  2008-07-11 17:06:27  ddangelo
 * added classes for neutron system (code by S. Davini)
 *
 * Revision 1.11  2008-05-10 23:34:06  ddangelo
 * BxDistance class defined and implemented
 * GetDistance and GetDisatnceError implemented in BxTrackFitted
 * some more work on track and position classes
 *
 * Revision 1.10  2008-05-09 15:46:23  ddangelo
 * muon track implemented as "by_points" and "fitted" classes inheriting from common interface (pure virtual class bx_track)
 * laben and muon host an istance "by_points"
 * echidna event hosts an istance "fitted"
 * structure duplicated at root level and filled
 * GetDistance duplicated in BxPosition class
 *
 * Revision 1.9  2008-02-21 15:58:31  ddangelo
 * added muon track as a separate class.
 * Currently used twice: laben and muon (later also globally fit).
 * both inner and roor event
 * added tracked stage in base event enum
 * bx_muon_rec_event renamed as bx_muon_tracked_event
 *
 * Revision 1.8  2007-12-06 16:48:09  ddangelo
 * upgraded to new muon clustering.
 * muon tracking side debugged and improved.
 *
 * Revision 1.7  2007-10-11 11:23:53  ddangelo
 * new mctruth format (internal AND root event)
 *
 * Revision 1.6  2007-06-22 16:38:22  razeto
 * Tags moved
 *
 * Revision 1.5  2007-05-07 15:06:38  ddangelo
 * BxPhysTags added in BxLabenRecCluster
 *
 * Revision 1.4  2007-03-22 16:09:35  ddangelo
 * added muon clustered level
 *
 * Revision 1.3  2007-03-15 19:17:19  ddangelo
 * pid event removed.
 * laben event upgraded with classes: bx_laben_shaped_cluster and bx_laben_ab_cluster
 * bx_laben_rec_cluster is now a parent class for the 2 new ones.
 * BxEvent modified accordingly: BxLabenRecHit and BxLabenRecCluster added.
 * BxPidEvent removed.
 *
 * Revision 1.2  2006/05/15 09:23:06  ddangelo
 * added class BxPhysTags to be filled by physics tools
 * Removed IgnoreObjectStreamer() calls from all constructors.
 *
 * Revision 1.1  2005/12/12 19:37:57  razeto
 * Moved BxEvent from root to event
 *
 * Revision 1.14  2005/12/12 19:08:09  razeto
 * Split position and energy reco results (auth from maintainers)
 *
 * Revision 1.13  2005/08/22 11:26:47  ddangelo
 * added Pid classes
 *
 * Revision 1.12  2005/07/11 17:16:59  ddangelo
 * removed global event
 * added pid event (partially)
 * untested
 *
 * Revision 1.11  2005/03/01 15:19:10  razeto
 * Merged with cycle_2
 *
 * Revision 1.10  2004/12/15 18:18:18  ddangelo
 * added position reco results. First draft, to be checked.
 *
 * Revision 1.9.2.1  2004/12/06 11:45:17  ddangelo
 * added pragmas for STL classes
 *
 * Revision 1.9  2004/12/01 15:13:41  ddangelo
 * added classes for fadc event.
 * added a few vairiables.
 * Work in progress, more stuff to come.
 *
 * Revision 1.8  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.7  2004/11/26 14:06:23  razeto
 * Added Mantainer field
 *
 * Revision 1.6  2004/07/25 16:29:31  ddangelo
 * some development, not much active at the moment
 *
 * Revision 1.5  2004/07/13 14:50:48  ddangelo
 * added BxClusteredHit. Currently used due to the TClonesArray problem.
 *
 * Revision 1.4  2004/07/13 13:37:10  ddangelo
 * added McTruth and McTruthFrame to root event.
 * McTruthHit commented out for the moment, due to ROOT problems.
 * To be debugged.
 *
 * Revision 1.3  2004/07/07 15:45:26  ddangelo
 * added BxLabenCluster.
 * Some minor debugging.
 *
 * Revision 1.2  2004/06/07 17:14:03  ddangelo
 * conditional compile macros introduced.
 * Laben raw and decoded hit introduced.
 * Muon raw and decoded hit implemented with 2 different classes.
 *
 * Revision 1.1  2004/05/30 11:54:48  ddangelo
 * A first working version of root file classes.
 * Not many physical variables yet;
 * Global Event still commented.
 * names match ROOT standards not echidna ones.
 * Makefile updated (file names).
 *
 *
 */
