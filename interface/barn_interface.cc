/* BOREXINO Reconstruction program
 * 
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: barn_interface.cc,v 1.11 2008/10/10 09:19:19 razeto Exp $
 *
 * Implementation of bx_root_barn
 *
 */

#include "barn_interface.hh"
#include "bx_barn.hh"
#include "messenger.hh"

#include <TKey.h>
#include <TSystem.h>
#include <TSocket.h>
#include <TMessage.h>
#include <TObjString.h>
#include <algorithm>

struct check_name {
  check_name (const std::string& name): s(name) {}
  bool operator() (const root_item& item) { return item.get_name () == s; }
  const std::string& s;
};

struct check_type {
  check_type (barn_interface::item_type type): t(type) {}
  bool operator() (const root_item& item) { return item.get_type () == t; }
  barn_interface::item_type t;
};


barn_interface* barn_interface::me = 0;

barn_interface::barn_interface(): bx_named("barn_interface"), socket(0), retries(10) {
 // items.clear ();
  send = get_parameter ("network_send").get_bool ();
  if (send) {
    socket = new TSocket (get_parameter ("network_address").get_string ().c_str (), get_parameter ("network_port").get_int ());
    if (socket && socket->IsValid ()) socket->Send (10001);
  }
}

barn_interface *barn_interface::get () {
  if (!me) me = new barn_interface;
  return me;
}

void barn_interface::network_send (TObject *obj, const bx_named *owner) {
  if (!obj || !socket) return;

  if (!socket->IsValid ()) {
    delete socket;
    get_message (bx_message::warn) << "network connection to " << get_parameter ("network_address") << ":" << get_parameter ("network_port") <<" broken" << dispatch;
    if (retries > 0 && obj->IsA () != TObjString::Class ()) {
      socket = new TSocket (get_parameter ("network_address").get_string ().c_str (), get_parameter ("network_port").get_int ());;
      retries --;
    } else socket = 0;
    return;
  }

  TMessage message (kMESS_OBJECT);
  TObjString owner_string (owner->get_name ().c_str ());
  message.WriteObject (&owner_string);
  message.WriteObject (obj);
  socket->Send (message);
}

void barn_interface::disable_root_signals () {
  gSystem->ResetSignal(kSigSegmentationViolation);
  gSystem->ResetSignal(kSigIllegalInstruction);
  gSystem->ResetSignal(kSigSystem);
  gSystem->ResetSignal(kSigPipe);
  gSystem->ResetSignal(kSigFloatingException);
}

void barn_interface::store (item_type type, TNamed *obj, const bx_named *owner) {

    // Check if a object is already present with that name
    std::list<root_item>::const_iterator item = std::find_if (items.begin (), items.end (), check_name (obj->GetName ()));
    
  if (item != items.end ()) {
    get_message (bx_message::error) << owner->get_name () << " registering a TNamed with a name already registered by " << item->get_owner_name () << " bx_named" << dispatch;
    return;
  }

  items.push_back (root_item (type, obj, owner));
}
TDirectory* barn_interface::m_get_dir (TDirectory* cwd, const std::string& dir_name) {
  TDirectory *dir;

  TKey *key = cwd->FindKey (dir_name.c_str ());
  if (key) { 
    dir = dynamic_cast<TDirectory*>(cwd->Get (dir_name.c_str ()));
  } else {
    dir = cwd->mkdir (dir_name.c_str ());
  }

  return dir;
}

TObject* barn_interface::m_get_obj (const root_item& item) {
  TObject *obj = gROOT->FindObjectAny (item.get_name ().c_str ());

  if (!obj) {
    get_message (bx_message::warn) << "object " << item.get_name () << " owned by " << item.get_owner_name () << " disappeared" << dispatch; 
  }
  
  return obj;
}

void barn_interface::m_obj_write (TDirectory *parent_dir, const root_item& item) {
  TObject *obj = m_get_obj (item);
  if (!obj) return;

  TDirectory* dir = m_get_dir (parent_dir, item.get_owner_name ());
  
  dir->cd ();
  int w = obj->Write ();
  if (w <= 0) {
    get_message (bx_message::warn) << "object " << item.get_name () << " owned by " << item.get_owner_name () << " was not written" << dispatch; 
  }
  dir->Save ();

  get_message (bx_message::debug) << "object " << obj->GetName () << " written" << dispatch;
}



void barn_interface::m_db_write (TDirectory *dir, bx_barn *obj) {		
  if (!obj) return;
    
  dir->cd ();
  int w = obj->Write ();
  if (w <= 0) {
    get_message (bx_message::warn) << "bx_barn was not written" << dispatch; 
  }
  dir->Save ();

  get_message (bx_message::debug) << "bx_barn written" << dispatch;
}
 
void barn_interface::write (TFile* file) {
  file->cd ();
  TDirectory* barn_dir = 0;

  if (get_parameter ("write_histos").get_bool ()) {
    barn_dir = m_get_dir (file, "barn");
    std::list<root_item>::const_iterator item;
    for (item = items.begin (); item != items.end (); item++) {
      if (item->get_type () != barn_interface::file) {
        get_message(bx_message::debug) << "skipping " << item->get_name () << dispatch; 
        continue;
      }
      barn_dir->cd ();
      m_obj_write (barn_dir, *item);
      file->Flush ();
    }
  }
  if (get_parameter ("write_barn").get_bool ()) {
    if (!barn_dir) barn_dir = m_get_dir (file, "barn");
    bx_barn* barn = new bx_barn("bxbarn", "BxBarn containing echidna histograms and database visitors");
    m_db_write(barn_dir, barn);
  }
  file->Flush ();
}

root_item::root_item (barn_interface::item_type type, TNamed *obj, const bx_named *owner): i_type(type), 
  name(obj->GetName ()), owner_name(owner->get_name ()) {}

/*
 * $Log: barn_interface.cc,v $
 * Revision 1.11  2008/10/10 09:19:19  razeto
 * Do not crash if not able to open socket
 *
 * Revision 1.10  2008-10-01 19:16:31  razeto
 * A little better
 *
 * Revision 1.9  2008-10-01 16:55:32  razeto
 * Fixed a bug
 *
 * Revision 1.8  2008-10-01 16:38:52  razeto
 * Auto-reconnect to viewer
 *
 * Revision 1.7  2008-09-29 21:28:27  razeto
 * Sending reset command for event display
 *
 * Revision 1.6  2007-07-17 08:48:09  razeto
 * barn working again (but root really sucks)
 *
 * Revision 1.5  2007-04-11 15:42:26  razeto
 * Delete sockek in case of failure (hope it fixes network_send)
 *
 * Revision 1.4  2006/09/09 19:00:53  razeto
 * Bug fix
 *
 * Revision 1.3  2006/09/09 14:18:58  razeto
 * Now use the write_histogram and write_barn parameters
 *
 * Revision 1.2  2006/08/28 16:06:57  ludhova
 * Fixed a bug in barn_interface (by Ale)
 *
 * Revision 1.1  2006/08/21 11:05:50  razeto
 * Moved bx_root_barn to barn_interface and db_barn to bx_barn
 *
 * Revision 1.11  2006/07/04 11:03:37  razeto
 * Fixed a bug in root barn
 *
 * Revision 1.10  2006/05/24 15:11:26  razeto
 * Upgraded to send owner name to network server
 *
 * Revision 1.9  2006/05/17 14:59:40  razeto
 * Added experimental interface sending object to network
 *
 * Revision 1.8  2006/01/02 21:21:45  razeto
 * Removed test target, renamed barn folder, moved and renamed db_barn
 *
 * Revision 1.7  2005/11/04 00:51:12  misiaszek
 * Implementation of m_db_write member function added. Now bx_root_barn
 * is storing db_barn object in the Echidna output file.
 *
 * Revision 1.6  2005/02/25 16:46:40  razeto
 * Upgraded to store histos even from bx_named
 *
 * Revision 1.5  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.4  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.3  2004/05/21 11:34:48  razeto
 * Disabled root signal handlers
 *
 * Revision 1.2  2004/05/21 08:34:52  razeto
 * Updated to working
 *
 * Revision 1.1  2004/05/18 15:07:01  razeto
 * Added
 *
 */
