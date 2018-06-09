/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * based on original design by Marco Pallavicini
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: BxEvent.hh,v 1.235 2015/08/26 11:22:04 misiaszek Exp $
 *
 * This is the event stored into root TTree
 * Names match ROOT style and not echidna style.
 * 
 */

#ifndef _BXEVENT_HH
#define _BXEVENT_HH
#ifndef CYCLE_NUMBER
#define CYCLE_NUMBER 18
#endif

#include <TMath.h>
#include <vector>
#include <map>
#include <math.h>
#include <time.h>
#include "Mach4Event.hh"
class db_channel_laben;
class db_channel_muon;
class bx_barn;

// struct of flags for write options
struct bx_write_opts {
  struct {
    bool raw;
    bool decoded;
    int clustered;
    int rec;
  } laben;
  struct {
    bool raw;
    bool decoded;
    bool clustered;
  } muon;
  int mctruth;
};

#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
class bx_echidna_event;
class bx_trigger_event;
class bx_track;
class bx_track_by_points;
class bx_track_fitted;
class bx_laben_event;
class bx_laben_raw_hit;
class bx_laben_decoded_hit;
class bx_laben_clustered_hit;
class bx_laben_cluster;
class bx_laben_rec_hit;
//class bx_laben_ab_mach4_cluster;
class bx_laben_positron_cluster;
class bx_base_position;
class bx_base_energy;
class bx_muon_event;
class bx_muon_raw_hit;
class bx_muon_decoded_hit;
class bx_muon_clustered_hit;
class bx_muon_cluster;
class bx_mctruth_hit;
class bx_mctruth_daughter;
class bx_mctruth_deposit;
class bx_mctruth_user;
class bx_mctruth_frame;
class bx_mctruth_event;
class bx_neutron_event;
class bx_neutron_pulse;
#endif

class BxEvent;
class BxPosition;
class BxTrackFitted;
class BxLaben;

/*********** BxTrigger *************/
class BxTrigger  {
  public:
    BxTrigger() {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    virtual void operator=(const bx_trigger_event&);
#endif
    virtual ~BxTrigger () {};

    enum TrgType {
	neutrino = 1, // Std trg of the ID (comparator of the BTB, Borexino Trigger Board)
	mtb = 2, // Muon Trigger Board (digital OD trigger)
	l355 = 4, // Laser 355nm, radial or oblique
	l394 = 8, // Laser 394nm, timing, radial or oblique
	l266 = 16, // Laser 266nm, radial or oblique
	calib = 32, // Calibration, pulse to Front-End
	random = 64, // "Empty" triggers
	neutron = 128   // thought for analog OD trigger, now neutron.
      };

    enum BtbFlag {
    	flag_mtb = 4,  // Muon Trigger Board (digital OD trigger)
	flag_neutron = 8, // thought for analog OD trigger, now neutron 
	flag_l355 = 16,  // Laser 355nm, radial or oblique
	flag_l266 = 32, // Laser 266nm, radial or oblique
	flag_lcr = 64 // Laser for timing OR Calibration OR Random
      };

    // getters
    UChar_t  GetTrgType               () const { return trgtype          ; }
    Bool_t   IsTrgType     ( TrgType t ) const { return trgtype==t       ; }
    ULong_t  GetTrgTime               () const { return trgtime          ; }
    UChar_t  GetBtbInputs             () const { return btb_inputs       ; }
    Bool_t   HasBtbFlag ( BtbFlag flag ) const { return btb_inputs & flag; }
    UShort_t GetBtbThreshold          () const { return btb_threshold    ; }
    const ULong_t* GetGpsTimes        () const { return gpstimes         ; }
    ULong_t  GetGpsTimeSec            () const { return gpstimes[0]      ; }
    ULong_t  GetGpsTimeNs             () const { return gpstimes[1]      ; }
    time_t GetTimeT 	              () const { return timet		 ; }
    Double_t GetMjd                   () const { return double ( GetTimeT () ) / 86400. + 40587.; }
    static time_t Mjd2TimeT (double mjd) { return time_t((mjd - 40587.) * 86400.); }
    Long_t GetRunDay () const { 
      time_t t = GetTimeT (); 
      struct tm * ts = gmtime (&t); 
      return (ts->tm_year + 1900 - 2007) * 365 + ts->tm_yday - 136; // 136:2007 = May 16 2007 => first run 
    } 
    Bool_t IsDay	     () const { Float_t v; return GetSunAltitude (v) > 0; }
    Bool_t IsNight           () const { return !IsDay(); }
    Float_t GetSunAltitude     (Float_t& azimuth) const { return GetSunAltitude (GetTimeT (), azimuth); }
    time_t GetSunRise	     () const { return GetSunRise (GetTimeT ()); }
    time_t GetSunSet	     () const { return GetSunSet (GetTimeT ()); }
    time_t GetMidday	     () const { return GetMidday (GetTimeT ()); }
    static float GetSunAltitude (time_t t, float& azimuth);
    static time_t GetSunRise (time_t t);
    static time_t GetSunSet  (time_t t);
    static time_t GetMidday  (time_t t);

    virtual Bool_t IsFolder() const { return kTRUE; }
  private:
    UChar_t  trgtype;         // trigger type
    ULong_t  trgtime;         // ppc0 cpu time 
    UShort_t btb_threshold;   // BTB threshold, n. channels 
    UChar_t  btb_inputs;      // BTB active inputs during event. 
			      // bifield (0,1,7 unused;2 mtb; 3 neutron; 4 l355; 5 l266; 6 l394+calib+random) 
    ULong_t gpstimes[2];      // GPS time, [0] seconds since 01.01.00 0:0:0, [1] nanoseconds
    time_t timet;	      // Time_t = number of seconds since 1970

    friend void fmerge_fake_bxevent(BxEvent* ev, Int_t run, Int_t evnum, UChar_t trgtype, Int_t nclu);
    ClassDef (BxTrigger, CYCLE_NUMBER)
};

// BxTrack class (used 3 times: BxEvent, BxMuon, BxLaben)
class BxTrack {
  public:
    virtual ~BxTrack () {};

    // getters
    virtual Float_t GetTheta   () const = 0;
    virtual Float_t GetPhi     () const = 0;
    virtual Float_t GetImpact  () const = 0;
    virtual Float_t GetDistance      (const BxPosition& p) const = 0;
    virtual Float_t GetDistanceError (const BxPosition& p) const = 0;
    virtual Float_t GetDistanceSigma (const BxPosition& p) const = 0;
    virtual BxPosition GetProjection (const BxPosition& p) const = 0;

    Float_t GetPathSSS   () const { Float_t i = GetImpact(); return (i<6.85) ? 2*::sqrtf(::powf(6.85,2)-::powf(i,2)) : 0.; }
    Float_t GetPathIV    () const { Float_t i = GetImpact(); return (i<4.25) ? 2*::sqrtf(::powf(4.25,2)-::powf(i,2)) : 0.; }
    Float_t GetPathBuffer() const { return GetPathSSS()-GetPathIV(); }
    virtual Float_t GetDedx () const = 0;

    virtual Bool_t IsUpward    () const = 0;
    virtual Bool_t IsDownward  () const = 0;

    virtual Bool_t IsFolder() const { return kTRUE; }

    ClassDef (BxTrack, CYCLE_NUMBER)
};

class BxTrackByPoints: public BxTrack {
  public:
    BxTrackByPoints() {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxTrackByPoints(const bx_write_opts&);
    virtual void operator=(const bx_track&);
#endif
    virtual ~BxTrackByPoints () {};

    // getters
    Float_t GetX1 () const { return x1;  }
    Float_t GetY1 () const { return y1;  }
    Float_t GetZ1 () const { return z1;  }
    Float_t GetT1 () const { return t1;  }
    Float_t GetX2 () const { return x2;  }
    Float_t GetY2 () const { return y2;  }
    Float_t GetZ2 () const { return z2;  }
    Float_t GetT2 () const { return t2;  }
    Float_t GetDX1() const { return dx1; }
    Float_t GetDY1() const { return dy1; }
    Float_t GetDZ1() const { return dz1; }
    Float_t GetDX2() const { return dx2; }
    Float_t GetDY2() const { return dy2; }
    Float_t GetDZ2() const { return dz2; }

    Float_t GetPedalT () const { return t1 + (t2-t1) * GetOA() * GetCosAlpha() / GetAB(); }
    Float_t GetPedalX () const { return x1 + (x2-x1) * GetOA() * GetCosAlpha() / GetAB(); }
    Float_t GetPedalY () const { return y1 + (y2-y1) * GetOA() * GetCosAlpha() / GetAB(); }
    Float_t GetPedalZ () const { return z1 + (z2-z1) * GetOA() * GetCosAlpha() / GetAB(); }

    virtual Float_t GetTheta   () const { return theta  ; }
    virtual Float_t GetPhi     () const { return phi    ; }
    virtual Float_t GetDTheta  () const { return dtheta ; }
    virtual Float_t GetDPhi    () const { return dphi   ; }
    virtual Float_t GetImpact  () const { return impact ; }
    virtual Float_t GetDImpact () const { return dimpact; }
    virtual Float_t GetDedx    () const ;
    Int_t GetError   () const { return error  ; }
//    virtual Float_t GetImpact  () const { return ::sqrtf ( ::powf(GetPedalX(), 2) + ::powf(GetPedalY(), 2) + ::powf(GetPedalZ(), 2) ); }
//    virtual Float_t GetTheta   () const { return ::acos ( (z1-z2) / GetAB()); }
//    virtual Float_t GetPhi     () const { Float_t phi = ::atan2(y1-y2,x1-x2); return (phi<0) ? phi + 2*TMath::Pi() : phi; }
    virtual Bool_t IsUpward    () const { return !downward; }
    virtual Bool_t IsDownward  () const { return downward ; }

    virtual Float_t GetDistance      (const BxPosition& p) const;
    virtual Float_t GetDistanceError (const BxPosition& p) const { return 1.; } // AAA To be implemented
    virtual Float_t GetDistanceSigma (const BxPosition& p) const { return GetDistance(p)/GetDistanceError(p); }

    virtual BxPosition GetProjection (const BxPosition& p) const;

    virtual Bool_t IsFolder() const { return kTRUE; }

  private:
    Float_t GetOA       () const { return ::sqrt( x1*x1 + y1*y1 + z1*z1 ); }
    Float_t GetAB       () const { return ::sqrt( pow((x1-x2), 2) + pow((y1-y2), 2) + pow((z1-z2), 2) ); }
    Float_t GetCosAlpha () const { return (-x1*(x2-x1)-y1*(y2-y1)-z1*(z2-z1))/(GetAB()*GetOA()); }
 
    Float_t t1;         // 'entry' point: time of crossing; 0-8500ns
    Float_t x1;         // 'entry' point: x (m);
    Float_t y1;         // 'entry' point: y (m);
    Float_t z1;         // 'entry' point: z (m);
    Float_t t2;         // 'exit' point: time of crossing; 0-8500/ns
    Float_t x2;         // 'exit' point: x (m);
    Float_t y2;         // 'exit' point: y (m);
    Float_t z2;         // 'exit' point: z (m);
    Float_t dx1;        // 'entry' point: error on x (m); 
    Float_t dy1;        // 'entry' point: error on y (m);
    Float_t dz1;        // 'entry' point: error on z (m);
    Float_t dx2;        // 'exit' point: error on x (m);
    Float_t dy2;        // 'exit' point: error on y (m);
    Float_t dz2;        // 'exit' point: error on z (m);
    Float_t theta;      // zenith angle: in radians;
    Float_t phi;        // azimuth angle: in radians;
    Float_t dtheta;     // zenith angle error: in radians;
    Float_t dphi;       // azimuth angle error: in radians;
    Float_t impact;     // impact parameter: in meters;
    Float_t dimpact;    // impact parameter error: in meters;
    Float_t idnhits;    // normalized laben decoded hits (for dedx computation);
    Int_t   error;      // error flag: used by cmt (codes to be defined); 
    Bool_t  downward;   // true for downward going muons

    ClassDef (BxTrackByPoints, CYCLE_NUMBER)
};

// internal use only, same getters are called from BxTrack and from BxPosition
class BxDistance {
  public:
    //BxDistance(const BxTrackFitted&, const BxPosition&);
    virtual ~BxDistance () {};

    Float_t GetDistance      () const { return distance               ; }
    Float_t GetDistanceError () const { return distance_error         ; }
    Float_t GetDistanceSigma () const { return distance/distance_error; } 

  private:
    Float_t distance;
    Float_t distance_error;

  ClassDef (BxDistance, CYCLE_NUMBER)
};

class BxTrackFitted: public BxTrack {
  public:
    BxTrackFitted() {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxTrackFitted(const bx_write_opts&);
    virtual void operator=(const bx_track&);
#endif
    virtual ~BxTrackFitted () {};

    // Class specific getters
    Double_t GetAlpha        () const { return alpha      ; }
    Double_t GetBeta         () const { return beta       ; }
    Double_t GetGamma        () const { return gamma      ; }
    Double_t GetDelta        () const { return delta      ; }
    Double_t GetAlphaError   () const { return alpha_error; }
    Double_t GetBetaError    () const { return beta_error ; }
    Double_t GetGammaError   () const { return gamma_error; }
    Double_t GetDeltaError   () const { return delta_error; }
    Float_t  GetChi2         () const { return chi2       ; }
    Bool_t   HasPoint (Int_t i) const { return points >> (i-1) & 0x01; }
    Int_t    GetNPoints      () const { Int_t n_points = 0; for (int i=0; i < 8; i++) n_points += points >> i & 0x01; return n_points; }

    virtual Float_t GetTheta   () const { return theta  ; }
    virtual Float_t GetPhi     () const { return phi    ; }
    virtual Float_t GetDTheta  () const { return dtheta ; }
    virtual Float_t GetDPhi    () const { return dphi   ; }
    virtual Float_t GetImpact  () const { return impact ; }
    virtual Float_t GetDImpact () const { return dimpact; }
    virtual Float_t GetDedx    () const ;
    
    // Get Entry and Exit Point

     virtual Double_t GetX1 (Double_t rad=6.85) const { 
   	Double_t a = 1+beta*beta+delta*delta; 
	Double_t b = 2*(alpha*beta+gamma*delta); 
	Double_t c = alpha*alpha+gamma*gamma-rad*rad;
	Double_t x1 = (-b+::sqrt(b*b-4*a*c))/(2*a);
	Double_t x2 = (-b-::sqrt(b*b-4*a*c))/(2*a);
	if ( (delta*x1>delta*x2 && downward) || (delta*x2>delta*x1 && !downward) ) return x1;
	else return x2;
    }
    virtual Double_t GetY1 (Double_t rad=6.85) const { return alpha+beta*GetX1(rad); } 
    virtual Double_t GetZ1 (Double_t rad=6.85) const { return gamma+delta*GetX1(rad); } 

    virtual Double_t GetX2 (Double_t rad=6.85) const { 
    	Double_t a = 1+beta*beta+delta*delta; 
	Double_t b = 2*(alpha*beta+gamma*delta); 
	Double_t c = alpha*alpha+gamma*gamma-rad*rad;
	Double_t x1 = (-b+sqrt(b*b-4*a*c))/(2*a);
	Double_t x2 = (-b-sqrt(b*b-4*a*c))/(2*a);
	if ( (delta*x1<delta*x2 && downward) || (delta*x2<delta*x1 && !downward) ) return x1;
	else return x2;
    }
    virtual Double_t GetY2 (Double_t rad=6.85) const { return alpha+beta*GetX2(rad); } 
    virtual Double_t GetZ2 (Double_t rad=6.85) const { return gamma+delta*GetX2(rad); } 

    virtual Float_t GetDistance      (const BxPosition& p) const; 
    virtual Float_t GetDistanceError (const BxPosition& p) const { return 1.; } 
    virtual Float_t GetDistanceSigma (const BxPosition& p) const { return GetDistance(p)/GetDistanceError(p); } 
    virtual BxPosition GetProjection (const BxPosition& p) const; 

    virtual Bool_t  IsUpward   () const { return !downward; }
    virtual Bool_t  IsDownward () const { return downward ; }

    virtual Bool_t IsFolder() const { return kTRUE; }

  private:
    Double_t alpha;       // 3-D fit parameter
    Double_t beta;        // 3-D fit parameter
    Double_t gamma;       // 3-D fit parameter
    Double_t delta;       // 3-D fit parameter
    Double_t alpha_error; // 3-D fit parameter error
    Double_t beta_error;  // 3-D fit parameter error 
    Double_t gamma_error; // 3-D fit parameter error
    Double_t delta_error; // 3-D fit parameter error=
    Float_t  chi2;        // 3-D fit chi2
    Float_t  phi;         // track orientiation
    Float_t  theta;
    Float_t  impact;
    Float_t  dphi;         // track orientiation
    Float_t  dtheta;
    Float_t  dimpact;
    Float_t  idnhits;     // normalized laben decoded hits (for dedx computation);
    Char_t   points;      // points used in fit; bitfield 0=entry OD; 1=exit OD; 
                          // 2=entry ID; 3=exit ID; 4-7=unused
    Bool_t   downward;    // true for downward going muons 

    ClassDef (BxTrackFitted, CYCLE_NUMBER)
};


/*********** BxLaben *************/
class BxLabenRawHit {
  public:
    BxLabenRawHit () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxLabenRawHit (const bx_laben_raw_hit&);
#endif
    virtual ~BxLabenRawHit () {};
 
    // getters
    UShort_t GetLg           () const { return lg;          }
    UChar_t  GetTime1        () const { return time1;       }
    UChar_t  GetTime2        () const { return time2;       }
    UShort_t GetGray         () const { return gray;        }
    UChar_t  GetBase         () const { return base;        }
    UChar_t  GetPeak         () const { return peak;        }
    UChar_t  GetOrder        () const { return order;       }
    UShort_t GetErrors       () const { return errors;      }
    UChar_t  GetFlagsBoard   () const { return flags_board; }
    UChar_t  GetFlagsChannel () const { return flags_ch;    }
    Bool_t   IsGood          () const { return !(flags_board & 0x7 || flags_ch & 0x60); }
    Bool_t   HasFifoFull     () const { return flags_ch    & 0x40; }
    Bool_t   HasFifoEmpty    () const { return flags_ch    & 0x20; }
    Bool_t   HasTrgJump      () const { return flags_board & 0x1;  }
    Bool_t   HasTrgJumpLarge () const { return flags_board & 0x2;  }
    Bool_t   HasTrgInBusy    () const { return flags_board & 0x4;  }
    Bool_t   IsInvalid       () const { return flags_board & 0x80; }

    virtual Bool_t IsFolder() const { return kTRUE; }

  private:
    UShort_t lg;          // logical channel (1 based)
    UChar_t  time1;       // First triangular waveform sampling (ADC bins, 0-255)
    UChar_t  time2;       // Second (~80ns later) triangular waveform sampling (ADC bins, 0-255)
    UShort_t gray;        // Gray counter read out in bins 0->(1<<16-1)
    UChar_t  base;        // Baseline charge sampling(ADC bins 0-255)}
    UChar_t  peak;        // Peak (~80ns later) charge sampling (ADC bins 0-255)}
    UChar_t  order;       // Order in channel as in daq fifo; 1-based
    UChar_t  flags_board; // Flags of board header: 1 = N.TRG not consecutive, 2 = N.TRG > DELTA, 3 = TRG in BUSY, 4 = DPR FULL, 5 = FIFO EMPTY, 6 FIFO FULL
    UChar_t  flags_ch;    // 5 = FIFO EMPTY, 6 = FIFO FULL, 7 = INVALID (...)
    UShort_t errors;      // bitfield. Laben board error flags. A.Razeto refused to provide info (DD)
    ClassDef (BxLabenRawHit, CYCLE_NUMBER)
};

class BxLabenDecodedHit {
  public:
    BxLabenDecodedHit () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxLabenDecodedHit (const bx_laben_decoded_hit&, unsigned short raw_index);
#endif
    virtual ~BxLabenDecodedHit () {};

    enum flag_type {
      out_of_gate = 1,
      reflection  = 2,
      reference   = 4,
      retrigger   = 8,
      disabled    = 16,
    };

    // getters
    UShort_t GetLg         () const { return lg;             }
    Double_t GetRawTime    () const { return raw_time;       }
//    Float_t  GetTimeError () const { return time_error;     }
//    Float_t  GetD80       () const { return d80;            }
    UChar_t  GetFlag       () const { return flag;           }
    Bool_t   IsGood        () const { return !flag;          } 
    Bool_t   IsOutOfGate   () const { return flag & out_of_gate; }
    Bool_t   IsReflection  () const { return flag & reflection;  } 
    Bool_t   IsReference   () const { return flag & reference;   } 
    Bool_t   IsRetrigger   () const { return flag & retrigger;   } 
    Float_t  GetRawCharge  () const { return raw_charge;     }
    Float_t  GetCharge     () const { return charge;         }
    Float_t  GetChargeMean () const { return charge_mean;    }
//    Int_t    GetNpe       () const { return npe;            }
    UChar_t  GetOrder      () const { return order;          }
    UShort_t GetRawIndex   () const { return raw_index;      }
    UShort_t GetNumCluster () const { return num_cluster;    }
    Float_t  GetRecTime    () const { return rec_time;       }
    Bool_t   IsShortCluster() const { return short_cluster;  }
//    BxLabenRawHit& GetRawHit (BxEvent* e) const { return *e->GetLaben().GetRawHits()[raw_index]; }    
//   const db_channel_laben* GetDbChannel (bx_barn* db) const;
    
    virtual Bool_t IsFolder() const { return kTRUE; }

  private:
    UShort_t lg;             // logical channel (1 based)
    Double_t raw_time;       // hit time; Units: ns; range: 0->6.4ms; 0=last gray crossing if < 3.2ms, 0=forelast gray crossing if > 3.2ms
//    Float_t  time_error;     // 0-0.5 hit time validity (good if <0.4)
    UChar_t  flag;           // Flag for out_of_gate=1, reflection=2, reference=4, retrigger=8, disabled=16;
    UChar_t  order;          // Order in channel counting decoded hits only; 1-based
//    Float_t  d80;            // Real 80ns gap in ns (from time2-time1)
    Float_t  raw_charge;     // peak-base in bins (without pileup corr.)
    Float_t  charge;         // photoelectrons (with pileup corr.)
    Float_t  charge_mean;    // photoelectrons normalized with mean (with pileup corr.)
//    Int_t    npe;	     // number of photoelectrons (with a threshold)
    UShort_t raw_index;      // Index in vector of corresponding raw hit
    UShort_t num_cluster;    // Number of cluster the hit belongs to
    Float_t  rec_time;       // Time after TOF subtraction
    Bool_t   short_cluster;  // Flag to say if the hit belonged to the old (c11) clustering 

    ClassDef (BxLabenDecodedHit, CYCLE_NUMBER)
  friend class BxLaben;
};

class BxLabenClusteredHit {
  public:
    BxLabenClusteredHit () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxLabenClusteredHit (const int cluster, const bx_laben_clustered_hit&, unsigned short decoded_index);
#endif
    virtual ~BxLabenClusteredHit () {};

    // getters
    Int_t    GetNumCluster   () const { return num_cluster;   }
    Float_t  GetTime         () const { return time;          }
    Float_t  GetCharge       () const { return charge;        }
  //    BxDecodedHit& GetDecodedHit (BxEvent* e) const { return *(e->GetLaben().GetDecodedHits())[decoded_hit_index]; }
    virtual Bool_t IsFolder() const { return kTRUE; }

//    const db_channel_laben* GetDbChannel (bx_barn* db) const;

  private:
    Int_t    num_cluster;    // Number of cluster (1,2,3) the hit belongs to in cluster array
    Float_t  time;           // Time relative to cluster start time in ns
    Float_t  charge;         // photoelectrons (with pileup corr.)

    ClassDef (BxLabenClusteredHit, CYCLE_NUMBER)
};

class BxPosition {
  public:
    BxPosition () {}
    BxPosition (Float_t xx, Float_t yy, Float_t zz): x(xx), y(yy), z(zz) {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxPosition (const bx_base_position&);
#endif
    virtual ~BxPosition () {};

    // getters
    Float_t GetTime   () const { return time;  }
    Float_t GetX      () const { return x;     }
    Float_t GetY      () const { return y;     }
    Float_t GetZ      () const { return z;     }
    Float_t GetDT     () const { return dt;    }
    Float_t GetDX     () const { return dx;     }
    Float_t GetDY     () const { return dy;     }
    Float_t GetDZ     () const { return dz;     }
    Float_t GetR      () const { return ::sqrt(x*x+y*y+z*z); }
    Float_t GetTheta  () const { return ::sqrt(x*x+y*y+z*z) > 1.e-6 ? ::acos(z/sqrt(x*x+y*y+z*z)) : 100.; }
    Float_t GetPhi    () const { if      ( x >  1.e-6 ) return y > 0. ? ::atan(y/x) : ::atan(y/x) + 2*::acos(-1.);
                                 else if ( x < -1.e-6 ) return ::atan(y/x) + ::acos(-1.);
                                 else                   return y > 0. ? ::acos(-1.)/2 : 3*::acos(-1.)/2; }
    Float_t GetLikelihood () const { return user;  }
    Float_t GetRefIndex   () const { return user;  }
    Float_t GetUser       () const { return user;  }
    Int_t   GetMatrix     () const { return matrix; }
    Bool_t  IsDefPos      () const { return matrix== 1; }
    Bool_t  IsNotDefPos   () const { return matrix==-1; }
    Bool_t  IsApproximate () const { return matrix== 0; }
    Bool_t  IsConverged   () const { return converged; }

    Float_t GetDistance (const BxPosition& p) const { return ::sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y) + (z - p.z) * (z - p.z)); }
    Float_t GetDistance (const BxTrackByPoints& t) const { return t.GetDistance(*this); };
    Float_t GetDistance (const BxTrackFitted&   t) const { return t.GetDistance(*this); };
    Float_t GetDistanceError (const BxTrackFitted& t) const { return t.GetDistanceError(*this); };
    Float_t GetDistanceSigma (const BxTrackFitted& t) const { return t.GetDistanceSigma(*this); };

    virtual Bool_t IsFolder() const { return kTRUE; }

  private:
    Float_t time;       // time of event start from fit, no standard zero available yet
    Float_t x;          // x from fit (m)
    Float_t y;          // y from fit (m)
    Float_t z;          // z from fit (m)
    Float_t dt;         // dt from fit (m)
    Float_t dx;         // dx from fit (m)
    Float_t dy;         // dy from fit (m)
    Float_t dz;         // dz from fit (m)
    Float_t user;       // mi: likelihood, value of the fcn at the minimum (normalized); mach4: index of refraction
    Bool_t  converged;  // Fit converged
    Int_t   matrix;     // Matrix def pos (+1), approximate (0) or not def pos (-1);

    ClassDef (BxPosition, CYCLE_NUMBER)
};

class BxEnergy {
  public:
    BxEnergy () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxEnergy (const bx_base_energy&);
#endif
    virtual ~BxEnergy () {}

    Int_t   GetNHits    () const { return nhits;    }
    Int_t   GetNpe      () const { return npe;      }
    Float_t GetCharge   () const { return charge;   }

    virtual Bool_t IsFolder() const { return kTRUE; }
    
  private:
    Int_t nhits;     // nhits corrected for position 
    Int_t npe;       // npe     "            "
    Float_t charge;  // charge  "            "

    ClassDef (BxEnergy, CYCLE_NUMBER)
};
    

class BxLabenCluster {
  public:
    BxLabenCluster () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxLabenCluster (const bx_laben_cluster&, unsigned short decoded_index);
#endif
    virtual ~BxLabenCluster () {};

    // getters
    Int_t    GetNPmts        () const { return npmts;         }
    Int_t    GetNPmtsConc    () const { return npmts_conc;    }
    Int_t    GetNPmtsShort   () const { return npmts_short;   }
    Int_t    GetNPmtsDt1     () const { return npmts_dt1;     }
    Int_t    GetNPmtsDt2     () const { return npmts_dt2;     }
    Int_t    GetNHits        () const { return nhits;         }
    Int_t    GetNHitsConc    () const { return nhits_conc;    }
    Int_t    GetNHitsShort   () const { return nhits_short;   }
//    Int_t    GetNpe          () const { return npe;           }
//    Int_t    GetNpeConc      () const { return npe_conc;      }
    Float_t  GetCharge       () const { return charge;        }
    Float_t  GetChargeConc   () const { return charge_conc;   }
    Float_t  GetChargeShort  () const { return charge_short;  }
    Float_t  GetChargeDt1    () const { return charge_dt1;    }
    Float_t  GetChargeDt2    () const { return charge_dt2;    }
    Float_t  GetChargeNpmts  () const { return charge_npmts;  }
    Float_t  GetChargeNoavgDt1  () const { return charge_noavg_dt1;  }
    Float_t  GetChargeNoavgDt2  () const { return charge_noavg_dt2;  }
    Float_t  GetChargeNoavg     () const { return charge_noavg;      }
    Float_t  GetChargeNoavgShort() const { return charge_noavg_short;}
    Double_t GetStartTime    () const { return start_time;    }
//    Double_t GetRoughTime    () const { return rough_time;    }
    Float_t  GetMeanTime     () const { return mean_time;     }    
    Float_t  GetMeanTimeShort() const { return mean_time_short;}    
    Float_t  GetRMSTime      () const { return rms_time;      }    
    Float_t  GetRMSTimeShort () const { return rms_time_short;}    
    Float_t  GetDuration     () const { return duration;      }    
    Float_t  GetDurationShort() const { return duration_short;}    
    enum flag_type {
      out_of_gate = 1,
      broad       = 2,
      trigger     = 4,
    };
    UChar_t GetFlag	     () const { return flag;               }
    Bool_t  IsGood	     () const { return !flag;              }
    Bool_t  IsOutOfGate      () const { return flag & out_of_gate; }
    Bool_t  IsBroad          () const { return flag & broad;       }
    Bool_t  IsTrigger        () const { return flag & trigger;     }
    Bool_t  IsNeutron        () const { return is_neutron;    }
    UShort_t GetDecodedIndex () const { return decoded_index; }
    Int_t    GetNPeaks       () const { return npeaks;        }

    const BxPosition& GetBaricenter  () const { return baricenter;   }
    const BxPosition& GetPositionLNGS  () const { return position_lngs;  }
    const BxPosition& GetPositionNoavg  () const { return position_noavg;  }
    Float_t 	      GetPeakTime    (Int_t i) const { return peak_times[i];   }
    Float_t	      GetPeakCharge  (Int_t i) const { return peak_charges[i]; }
    const std::vector<Float_t>& GetPeakTimes  () const { return peak_times;   }
    const std::vector<Float_t>& GetPeakCharges() const { return peak_charges; }

    virtual Bool_t IsFolder() const { return kTRUE; }

    Int_t    Normalize_geo_Pmts   (Int_t arg )   const { return npmt_geo_weight   ? Int_t((Float_t(arg) / Float_t(npmt_geo_weight) * 2000.)) : 0 ; }
    Float_t  Normalize_geo_Pmts   (Float_t arg ) const { return npmt_geo_weight   ? (arg / npmt_geo_weight * 2000) : 0. ; }
    Float_t  Normalize_geo_Charge (Float_t arg ) const { return charge_geo_weight ? (arg / charge_geo_weight * 2000) : 0. ; }
    Int_t    Normalize_QE_Pmts   (Int_t arg )   const { return npmt_QE_weight   ? Int_t((Float_t(arg) / Float_t(npmt_QE_weight) * 2000.)) : 0 ; }
    Float_t  Normalize_QE_Pmts   (Float_t arg ) const { return npmt_QE_weight   ? (arg / npmt_QE_weight * 2000) : 0. ; }
    Float_t  Normalize_QE_Charge (Float_t arg ) const { return charge_QE_weight ? (arg / charge_QE_weight * 2000) : 0. ; }
    Int_t    Normalize_geo_QE_Pmts   (Int_t arg )   const { return npmt_geo_QE_weight   ? Int_t((Float_t(arg) / Float_t(npmt_geo_QE_weight) * 2000.)) : 0 ; }
    Float_t  Normalize_geo_QE_Pmts   (Float_t arg ) const { return npmt_geo_QE_weight   ? (arg / npmt_geo_QE_weight * 2000) : 0. ; }
    Float_t  Normalize_geo_QE_Charge (Float_t arg ) const { return charge_geo_QE_weight ? (arg / charge_geo_QE_weight * 2000) : 0. ; }

  private:
    Int_t    npmts;            // number of hit pmts
    Int_t    npmts_conc;       // number of hit pmts with light concentrators
    Int_t    npmts_short;      // number of hit pmts (old c11 clustering)
    Int_t    npmts_dt1;        // number of hit pmts in 230ns after cluster start time
    Int_t    npmts_dt2;        // number of hit pmts in 400ns after cluster start time
    Int_t    npmts_pos;        // number of hit pmts position corrected
    Int_t    nhits;            // number of hits in cluster
    Int_t    nhits_conc;       // number of hits in cluster with light concentrators
    Int_t    nhits_short;      // number of hits in old (c11) clustering style
    Int_t    nhits_pos;        // number of hits position corrected
//    Int_t    npe;              // number of photoelectrons
//    Int_t    npe_conc;         // number of photoelectrons with light concentrators
    Float_t  charge;           // photoelectron charge
    Float_t  charge_conc;      // photoelectron charge with light concentrators
    Float_t  charge_short;     // photoelectron charge in old (c11) clustering style
    Float_t  charge_dt1;       // photoelectron charge in dt1
    Float_t  charge_dt2;       // photoelectron charge in dt2
    Float_t  charge_pos;       // photoelectron charge position corrected
    Float_t  charge_npmts;     // photoelectron charge for first hit on pmt
    Float_t  charge_noavg_dt1; // photoelectron charge noavg in dt1
    Float_t  charge_noavg_dt2; // photoelectron charge noavg in dt2
    Float_t  charge_noavg;     // photoelectron charge noavg
    Float_t  charge_noavg_short;// photoelectron charge noavg short
    Double_t start_time;       // cluster start time from clustering module; Units: ns; 
                               // range: 0->6.4ms; 0=last gray crossing if < 3.2ms, 0=forelast gray crossing if > 3.2ms
//    Double_t rough_time;       // cluster rough start time from clustering module; same unit as start_time
    Float_t  mean_time;        // cluster mean time in ns relative to start_time 
    Float_t  mean_time_short;  // cluster mean time in old (c11) clustering style
    Float_t  rms_time;         // cluster rms time in ns relative to mean_time
    Float_t  rms_time_short;   // cluster rms time in ns relative to mean_time (old c11 clustering)
    Float_t  duration;         // time of the last hit in ns relative to start time
    Float_t  duration_short;   // time of the last hit in ns relative to start time (old c11 clutering)
    UChar_t  flag;             // Flag for out_of_gate=1, broad=2, trigger=4
    Bool_t   is_neutron;       // shows whether cluster meets neutron definition
    UShort_t decoded_index;    // Index in vector of decoded hit corresponding to first hit
    Int_t    npeaks;           // # of peaks identified by the splitting algorythm
    BxPosition baricenter;     // reconstructed position from baricentrator
    BxPosition position_lngs;  // reconstructed position from LNGS reco
    BxPosition position_noavg;  // reconstructed position from LNGS reco with noavg corrected charge
    Float_t  npmt_geo_weight;     // npmt normalization factor considering geometry
    Float_t  npmt_QE_weight;      // npmt normalization factor considering QE
    Float_t  npmt_geo_QE_weight;  // npmt normalization factor considering geometry and QE
    Float_t  charge_geo_weight;   // charge normalization factor considering geometry
    Float_t  charge_QE_weight;    // charge normalization factor considering QE
    Float_t  charge_geo_QE_weight;// charge normalization factor considering geometry and QE

    std::vector<Float_t> peak_times; //->
    std::vector<Float_t> peak_charges; //->
                                                                                
    ClassDef (BxLabenCluster, CYCLE_NUMBER)
};


class BxLabenRecHit {
  public:
    BxLabenRecHit () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxLabenRecHit (const int cluster, const bx_laben_rec_hit&);
#endif
    virtual ~BxLabenRecHit () {};

    // getters
    Int_t   GetNumCluster () const { return num_cluster; }
    Float_t GetTime       () const { return time;        }
    virtual Bool_t IsFolder() const { return kTRUE; }

//    const db_channel_laben* GetDbChannel (bx_barn* db) const;

  private:
    Int_t   num_cluster; // Number of cluster (1,2,3) the hit belongs to in cluster array
    Float_t time;        // Time relative to energy deposition in ns (i.e. TOF receorrected)
    ClassDef (BxLabenRecHit, CYCLE_NUMBER)
};


class BxLabenRecCluster: public BxPosition/*, public BxEnergy*/ {
  public:
    BxLabenRecCluster () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxLabenRecCluster (const bx_laben_positron_cluster&);
#endif
    virtual ~BxLabenRecCluster () {};

    // getters pid shape
    Float_t        GetNsAsymmetry        () const { return ns_asymmetry;   }
    Float_t        GetSphereChi2         () const { return sphere_chi2;    }
    Float_t        GetSphereLkl          () const { return sphere_lkl;     }
    Float_t        GetSphereRelVar       () const { return sphere_rel_var; }
    Float_t        GetPlaneCos           () const { return plane_cos;      }
    Float_t        GetPlaneChi2          () const { return plane_chi2;     }
    Float_t        GetHPlaneChi2         () const { return h_plane_chi2;   }
    Float_t        GetSHPower (Int_t order) const { return sh_power[order];}
    Char_t         GetQualityFlags       () const { return quality_flags;  }

    // getters pid ab
    const Float_t* GetTailTot            () const { return tailtot;        }
    const Float_t* GetTailTotAbMlp       () const { return tailtot_ab_mlp; }
    const Float_t* GetTailTotC11Mva      () const { return tailtot_c11_mva;}
    Float_t        GetTailTot  (Int_t tail) const; 
    Float_t        GetTailTotAbMlp   (Int_t tail) const; 
    Float_t        GetTailTotC11Mva  (Int_t tail) const; 
    Float_t        GetGatti              () const { return gatti;          }
    Float_t        GetGattiC             () const { return gattic;         }
    Float_t        GetLkl                () const { return lkl;            }
    Float_t        GetLklC               () const { return lklc;           }
    Float_t        GetRiseTime           () const { return rise_time;      }
    Float_t        GetRms           	 () const { return rms;	           }
    Float_t        GetRmsC11           	 () const { return rms_c11;        }
    Float_t        GetKurtosis         	 () const { return kurtosis;	   }
    Float_t        GetKurtosisC11      	 () const { return kurtosis_c11;   }
    Float_t        GetMlpAb         	 () const { return mlp_ab;	   }

    // getters pid ab mach4
/*    const Float_t* GetTailTotMach4         () const { return tailtot_mach4;  }
    Float_t        GetTailTotMach4  (Int_t t) const; 
    Float_t        GetGattiMach4    (Int_t i) const { return gatti_mach4[i]; }
    Float_t        GetMeanMach4            () const { return mean_mach4;     }
    Float_t        GetPeakMach4            () const { return peak_mach4;     }
    Float_t        GetRmsMach4             () const { return rms_mach4;      }
    Float_t        GetSkewMach4            () const { return skew_mach4;     }
    Float_t        GetKurtMach4            () const { return kurt_mach4;     }
*/

    // getters pid positron 
    Float_t        GetGattiOpsBeta         () const { return gatti_ops_beta; }
    Float_t        GetGattiC11Beta         () const { return gatti_c11_beta; }
    Float_t        GetGattiOpsNops         () const { return gatti_ops_nops; }

    virtual Bool_t IsFolder() const { return kTRUE; }

  private:
    // pid shape data members
    Float_t ns_asymmetry;         // north-south charge asymmetry
    Float_t sphere_chi2;          // chi2 to sphericity
    Float_t sphere_lkl;           // likelihood to sphericity
    Float_t sphere_rel_var;       // relative variance = sigma/mean of the hits in the cos(theta) - phi plane centered in the event position and normalized per hit charge. 0 to 2 approx.
    Float_t plane_cos;            // cos of the angle made by the plane fitting the theta-phi parameters space
    Float_t plane_chi2;           // chi2/NDF of the fit of a plane in a cos(theta)-phi distribution
    Float_t h_plane_chi2;         // chi2/NDF of the fit of a horizontal plane in a cos(theta)-phi distribution
    Float_t sh_power[4]; 	  // power of spherical armonics		
    Char_t  quality_flags;        // bitfield, meaning to be assigned

    // alpha-beta data members
    Float_t tailtot[10];          // tailtot for tail value 40-130ns step 10ns
    Float_t tailtot_ab_mlp[10];   // tailtot for a/b discrimination with mlp
    Float_t tailtot_c11_mva[10];  // tailtot for c11/b discrimination with mva
    Float_t gatti;                // gatti optimal filter variable
    Float_t lkl;                  // likelihood ratio
    Float_t gattic;               // gatti optimal filter variable (cumulative)
    Float_t lklc;                 // likelihood ratio (cumulative)
    Float_t rise_time;            // rise time
    Float_t rms;	          // rms
    Float_t rms_c11;	          // rms
    Float_t kurtosis;	          // kurtosis
    Float_t kurtosis_c11;         // kurtosis
    Float_t mlp_ab;		  // mlp value for a/b

    // pid positron gatti members
    Float_t gatti_ops_beta; 	  // gatti parameter for otho-positronium/beta discrimination; ops shapes from bxmc2 (tau = 3.1 ns); beta shapes from DATA 214Bi
    Float_t gatti_c11_beta; 	  // gatti parameter for c11/beta discrimination; c11 shapes from DATA TFC; beta shapes from DATA 214Bi
    Float_t gatti_ops_nops; 	  // gatti parameter for otho-positronium/no-otho-positronium c11 discrimination; ops shapes from bxmc2 (tau = 3.1 ns); nops shapes from bxmc2

    // alpha-beta from mach4 module
/*    Float_t tailtot_mach4[10];    // tailtot for tail value 30-110ns step 5ns
    Float_t gatti_mach4[4];       // gatti optimal filter variable from 4 reference shapes.
    Float_t peak_mach4;           // peak     of tof-recorrected hit time distribution. 
    Float_t mean_mach4;           // mean     of tof-recorrected hit time distribution.
    Float_t rms_mach4;            // rms      of tof-recorrected hit time distribution.
    Float_t skew_mach4;           // skewness of tof-recorrected hit time distribution.
    Float_t kurt_mach4;           // kurtosis of tof-recorrected hit time distribution.
*/
    ClassDef (BxLabenRecCluster, CYCLE_NUMBER)
};


class BxLaben {
  public:
    BxLaben() {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxLaben(const bx_write_opts&);
    virtual void  operator=(const bx_laben_event&);
#endif
    virtual ~BxLaben () {};

    // getters
    Int_t    GetEmptyBoards        () const { return empty_boards;          }
    Double_t GetTriggerTime        () const { return trigger_time;          }
    Double_t GetLaserTime          () const { return laser_time;            }
    Int_t    GetNLivePmts          () const { return n_live_pmts;           }
    Int_t    GetNLiveCharge        () const { return n_live_charge;         }
    Int_t    NormalizePmts   (Int_t arg )   const { return n_live_pmts   ? Int_t((Float_t(arg) / Float_t(n_live_pmts-n_invalid_pmts) * 2000.)) : 0 ; }
    Float_t  NormalizePmts   (Float_t arg ) const { return n_live_pmts   ? (arg / (n_live_pmts-n_invalid_pmts)   * 2000) : 0. ; }
    Float_t  NormalizeCharge (Float_t arg ) const { return n_live_charge ? (arg / (n_live_charge-n_invalid_charge) * 2000) : 0. ; }
//    Float_t  GetNpe                () const { return npe;                   }	
    Int_t    GetNPmts              () const { return npmts;                }
    Float_t  GetCharge             () const { return charge;                }
//    Float_t  GetClusterWindowLimit () const { return cluster_window_limit;  }
    Int_t    GetNHitsOnEmpty	   () const { return n_hits_on_empty;       }

    Int_t    GetNRawHits       () const { return n_raw_hits;                }
    Int_t    GetNRawHitsFW     () const { return n_raw_hits_fw;             }
    Int_t    GetNGood          () const { return n_raw_hits_flags[0];       }
    Int_t    GetNFifoFull      () const { return n_raw_hits_flags[1];       }
    Int_t    GetNFifoEmpty     () const { return n_raw_hits_flags[2];       }
    Int_t    GetNCounter       () const { return n_raw_hits_flags[3];       }
    Int_t    GetNTrgJump       () const { return n_raw_hits_flags[4];       }
    Int_t    GetNTrgJumpLarge  () const { return n_raw_hits_flags[5];       }
    Int_t    GetNTrgInBusy     () const { return n_raw_hits_flags[6];       }
    Int_t    GetNInvalid       () const { return n_raw_hits_flags[7];       }
    Int_t    GetNInvalidPmts   () const { return n_invalid_pmts;            }
    Int_t    GetNInvalidCharge () const { return n_invalid_charge;          }
    Int_t    GetNDecodedHits   () const { return n_decoded_hits;            }
    Int_t    GetNClusters      () const { return n_clusters;                }
    Int_t    GetNClustersMuons () const { return n_clusters_muons;          }
    Int_t    GetNClustersFound () const { return n_clusters_found;          }
    Int_t    GetNClustersOld   () const { return n_clusters_old;            }
    Int_t    GetNClustersNeutron () const { return n_clusters_neutron;      }
    Int_t    GetNClusteredHits () const { return n_clustered_hits;          }
    Bool_t   IsTrackedEnergy   () const { return is_tracked_energy;         }
    Bool_t   IsTrackedTof      () const { return is_tracked_tof;            }
    Bool_t   IsClusterOld      (const BxLabenCluster& c) const { return c.GetNHits() >= 40. * (1.-Float_t(empty_boards)/280.) * Float_t(n_live_pmts-n_invalid_pmts) / 2000.; }

    Bool_t   HasRawHits        () const { return has_raw;                   }
    Bool_t   HasDecodedHits    () const { return has_decoded;               }
    Bool_t   HasClusters       () const { return Bool_t(has_clustered > 0); }
    Bool_t   HasClusteredHits  () const { return Bool_t(has_clustered > 1); }
    Bool_t   HasRecClusters    () const { return Bool_t(has_rec > 0);       }
    Bool_t   HasRecHits        () const { return Bool_t(has_rec > 1);       }

    const std::vector<BxLabenRawHit>&       GetRawHits             () const { return raw_hits;          }
    const BxLabenRawHit&       		    GetRawHit       (Int_t i) const { return raw_hits[i];       }
    const std::vector<BxLabenDecodedHit>&   GetDecodedHits         () const { return decoded_hits;      }
    const BxLabenDecodedHit&   		    GetDecodedHit   (Int_t i) const { return decoded_hits[i];   }
    const std::vector<BxLabenCluster>&      GetClusters            () const { return clusters;          }
    const BxLabenCluster&      		    GetClusterMuon  (Int_t i) const { return clusters_muons[i]; }
    const std::vector<BxLabenCluster>&      GetClustersMuons       () const { return clusters_muons;    }
    const BxLabenCluster&      		    GetCluster      (Int_t i) const { return clusters[i];       }
    const std::vector<BxLabenClusteredHit>& GetClusteredHits       () const { return clustered_hits;    }
    const BxLabenClusteredHit& 		    GetClusteredHit (Int_t i) const { return clustered_hits[i]; }
    const std::vector<BxLabenRecCluster>&   GetRecClusters         () const { return rec_clusters;      }
    const BxLabenRecCluster&      	    GetRecCluster   (Int_t i) const { return rec_clusters[i];   }
    const std::vector<BxLabenRecHit>&       GetRecHits             () const { return rec_hits;          }
    const BxLabenRecHit& 	            GetRecHit       (Int_t i) const { return rec_hits[i];       }
    const std::vector<Int_t>&               GetNPmtsWin1Vec        () const { return npmts_win1;        }
    const std::vector<Int_t>&               GetNPmtsWin2Vec        () const { return npmts_win2;        }
    const std::vector<Float_t>&             GetChargeWin1Vec       () const { return charge_win1;        }
    const std::vector<Float_t>&             GetChargeWin2Vec       () const { return charge_win2;        }
    const BxTrackByPoints&                  GetTrackEnergy         () const { return track_energy;      } 
    const BxTrackByPoints&                  GetTrackTof            () const { return track_tof;         } 

    virtual Bool_t IsFolder() const { return kTRUE; }

  private:
    Int_t    empty_boards;        // debugging; number of empty boards;
    Double_t trigger_time;        // trigger time from (averaged) reference laben channels
                                  // Units: ns; range: 0->6.4ms; 0=last gray crossing if < 3.2ms, 0=forelast gray crossing if > 3.2ms;
    Double_t laser_time;          // laser time from (averaged) reference laben channels
                                  // Units: ns; range: 0->6.4ms; 0=last gray crossing if < 3.2ms, 0=forelast gray crossing if > 3.2ms;
//    Float_t  npe;                 // Sum of all decoded hits charge. 1pe/hit minimum charge assigned. 
    Int_t    npmts;               // Number of Pmts hit
    Float_t  charge;              // Sum of all decoded hits charge in photoelectrons
    Int_t    n_live_pmts;         // Number of ordinary channels alive, computed run time;
    Int_t    n_live_charge;       // Number of ordinary channels with good charge signal, computed run time;
    Int_t    n_hits_on_empty;	  // Number of hits on empty channels (epmty or dead)
    Int_t    n_raw_hits;          // number of raw hits       (even if vector is not written)
    Int_t    n_raw_hits_fw;       // number of raw hits calculated by laben new firmwar (even if vector is not written)
    Int_t    n_raw_hits_flags[8]; // n_good, n_fifo_full, n_fifo_empty, n_trg_jump, n_trg_jump_large, n_trg_in_busy, n_invalid (as from bx_laben_raw_hit::flags)
    Int_t    n_invalid_pmts;      // number of ordinary channels (and not disabled timing) for which FE or FF is set
    Int_t    n_invalid_charge;    // number of ordinary channels (and not disabled charge) for which FE or FF is set
    Int_t    n_decoded_hits;      // number of decoded hits   (even if vector is not written)
    Int_t    n_clusters;          // number of clusters       (even if vector is not written)
    Int_t    n_clusters_muons;    // number of clusters for neutron algorythm in muon gate (even if vector is not written)
    Int_t    n_clusters_found;    // number of clusters found by algorythm (may differ from previous for high multiplicity events)
    Int_t    n_clusters_old;      // tt128 only: number of clusters found by the old (c13) clustering algorythm
    Int_t    n_clusters_neutron;  // tt128 only: number of clusters above the neutron threshold
    Int_t    n_clustered_hits;    // number of clustered hits (even if vector is not written)
    Bool_t   has_raw;             // = true if raw hits are written in file;
    Bool_t   has_decoded;         // = true if decoded hits are written in file
    Int_t    has_clustered;       // = 1 if clusters are written in file; = 2 if also clustered hits are written
    Int_t    has_rec;             // = 1 if clusters are written in file; = 2 if also rec hits are written
//    Float_t  cluster_window_limit;// the maximum window for clustering
    Bool_t   is_tracked_energy;   // Is the track valid
    Bool_t   is_tracked_tof;      // Is the track valid
    std::vector<BxLabenRawHit>       raw_hits;       //->
    std::vector<BxLabenDecodedHit>   decoded_hits;   //->
    std::vector<BxLabenCluster>      clusters;       //->
    std::vector<BxLabenCluster>      clusters_muons; //->
    std::vector<BxLabenClusteredHit> clustered_hits; //->
    std::vector<BxLabenRecCluster>   rec_clusters;   //->
    std::vector<BxLabenRecHit>       rec_hits;       //->
    std::vector<Int_t>               npmts_win1;     //->
    std::vector<Int_t>               npmts_win2;     //->
    std::vector<Float_t>             charge_win1;    //->
    std::vector<Float_t>             charge_win2;    //->
    BxTrackByPoints                  track_energy;     
    BxTrackByPoints                  track_tof;     

    friend class n_live_charge;
    friend void fmerge_fake_bxevent(BxEvent* ev, Int_t run, Int_t evnum, UChar_t trgtype, Int_t nclu);
    ClassDef (BxLaben, CYCLE_NUMBER)
};


/*********** BxMuon *************/
class BxMuonRawHit {
  public:
    BxMuonRawHit () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxMuonRawHit (const bx_muon_raw_hit&);
#endif
    virtual ~BxMuonRawHit () {};
 
    // getters
    UShort_t GetMch       () const { return mch;        }
    Float_t  GetLeadTime  () const { return lead_time;  }
    Float_t  GetTrailTime () const { return trail_time; }

    virtual Bool_t IsFolder () const { return kTRUE; }

  private:
    UShort_t mch;        // muon channel (0-255)
    ULong_t  lead_time;  // leading edge time; Units: TDC ticks (1.0416ns); range: 0->8192
    ULong_t  trail_time; // trailing edge time; Units: TDC ticks (1.0416ns); range: 0->8192

    ClassDef (BxMuonRawHit, CYCLE_NUMBER)
};

class BxMuonDecodedHit {
  public:
    BxMuonDecodedHit () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxMuonDecodedHit (const bx_muon_decoded_hit&);
#endif
    virtual ~BxMuonDecodedHit () {};
 
    // getters
    UShort_t GetMch    () const { return mch;    }
    Float_t  GetTime   () const { return time;   }
    Float_t  GetCharge () const { return charge; }
//    const db_channel_muon* GetDbChannel (bx_barn* db) const;

    virtual Bool_t IsFolder () const { return kTRUE; }

  private:
    UShort_t run;	     // run number
    UShort_t mch;    // muon channel (0-255)
    Float_t  time;   // time after gate start (neutrino) or after led ref (laser); ns; 0->8.5$us
    Float_t  charge; // Charge (pe)

    ClassDef (BxMuonDecodedHit, CYCLE_NUMBER)
};

class BxMuonClusteredHit {
  public:
    BxMuonClusteredHit () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxMuonClusteredHit (const bx_muon_clustered_hit&);
#endif
    virtual ~BxMuonClusteredHit () {};
 
    // getters
    UShort_t GetMch    () const { return mch;    }
    Float_t  GetTime   () const { return time;   }
    Float_t  GetCharge () const { return charge; }
    //const db_channel_muon* GetDbChannel (bx_barn* db) const;

    virtual Bool_t IsFolder () const { return kTRUE; }

  private:
    UShort_t run;     // run number
    UShort_t mch;    // muon channel (0-255)
    Float_t  time;   // time after gate start (neutrino) or after led ref (laser); ns; 0->8.5$us
    Float_t  charge; // Charge (pe)

    ClassDef (BxMuonClusteredHit, CYCLE_NUMBER)
};

class BxMuonCluster {
  public:
    BxMuonCluster () {} 
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxMuonCluster (const bx_muon_cluster&);
#endif
    virtual ~BxMuonCluster () {};

    Int_t   GetId        () const { return id; }
    Float_t GetCharge    () const { return charge; }
    Float_t GetX         () const { return x; }
    Float_t GetY         () const { return y; }
    Float_t GetZ         () const { return z; }
    Float_t GetRadius    () const { return sqrtf(x*x+y*y+z*z); }
    Float_t GetRc        () const { return sqrtf(x*x+y*y); }
    Float_t GetTheta     () const { return ::acos(z/GetRadius()); }
    Float_t GetPhi       () const { return (y > 0) ? ::acos(x/GetRc()) : 2*TMath::Pi()-::acos(x/GetRc()); }
    Float_t GetStartTime () const { return start_time; }
    Bool_t  IsUp         () const { return z > 0.; }
    Bool_t  IsDown       () const { return (z < 0.) && (z > -5.); }
    Bool_t  IsSss        () const { return z > -5.; }
    Bool_t  IsFloor      () const { return z <= -5.; }

  private:
    Int_t   id;         // Internal cluster Id
    Float_t x;          // x of charge baricenter, exp time weighted
    Float_t y;          // y of "
    Float_t z;          // z of "
    Float_t charge;     // sum of hits charge, pe 
    Float_t start_time; // time of first hit, ns after gate start

    ClassDef (BxMuonCluster, CYCLE_NUMBER)
};



class BxMuon {
  public:
    BxMuon() {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxMuon(const bx_write_opts&);
    virtual void operator=(const bx_muon_event&);
#endif
    virtual ~BxMuon () {};

    // getters
    Bool_t  IsAligned              () const { return is_aligned;                }
    Int_t   GetNRawHits            () const { return n_raw_hits;                }
    Int_t   GetNDecodedHits        () const { return n_decoded_hits;            }
    Int_t   GetNClusters           () const { return n_clusters;                }
    Int_t   GetNClusteredHitsSSS   () const { return n_clustered_hits_sss;      }
    Int_t   GetNClusteredHitsFloor () const { return n_clustered_hits_floor;    }
    Int_t   GetNClusteredHits      () const { return n_clustered_hits_sss + n_clustered_hits_floor; }
    Int_t   GetDecodedNPmts        () const { return decoded_npmts;             }
    Float_t GetDecodedCharge       () const { return decoded_charge;            }
    Bool_t  HasClusterSSS          () const { return has_cluster_sss;           }
    Bool_t  HasClusterFloor        () const { return has_cluster_floor;         }
    Bool_t  HasCluster             () const { return has_cluster_sss || has_cluster_floor; }
    Int_t   GetNPmts               () const { return npmts;                     }
    Float_t GetStartTimeSSS        () const { return start_time_sss;            }
    Float_t GetStartTimeFloor      () const { return start_time_floor;          }
    Float_t GetStartTime           () const { return (start_time_sss < start_time_floor) ? start_time_sss : start_time_floor; }
    Float_t GetChargeSSS           () const { return charge_sss;                }
    Float_t GetChargeFloor         () const { return charge_floor;              }
    Float_t GetCharge              () const { return charge_sss + charge_floor; }
    Bool_t  IsTracked              () const { return is_tracked;                }
    Bool_t  HasRaw                 () const { return has_raw;                   }
    Bool_t  HasDecoded             () const { return has_decoded;               }
    Bool_t  HasClustered           () const { return has_clustered;             }
    const std::vector<BxMuonRawHit      >& GetRawHits       () const { return raw_hits;       }
    const std::vector<BxMuonDecodedHit  >& GetDecodedHits   () const { return decoded_hits;   }
    //const std::vector<BxMuonClusteredHit>& GetClusteredHits () const { return clustered_hits; }
    const std::vector<BxMuonCluster     >& GetClusters      () const { return clusters;       }
    const BxTrackByPoints                & GetTrack         () const { return track;          } 

    virtual Bool_t IsFolder() const { return kTRUE; }

  private:

    Bool_t   is_aligned;             // did not loose sycronicyty
    Int_t    n_raw_hits;             // number of raw hits (even if vector is not written)
    Int_t    n_decoded_hits;         // number of decoded hits (even if vector is not written)
    Int_t    n_clusters;             // number of clusters (even if vector is not written)
    Int_t    n_clustered_hits_sss;   // number of clustered hits on sss(even if vector is not written)
    Int_t    n_clustered_hits_floor; // number of clustered hits on floor (even if vector is not written)
    Int_t    decoded_npmts;          // number of fired pmts after decoding
    Float_t  decoded_charge;         // photoelectron charge after decoding
    Bool_t   has_cluster_sss;        // data present
    Bool_t   has_cluster_floor;      // data present
    Int_t    npmts;                  // number of fired pmts
    Float_t  start_time_sss;         // cluster start time from clustering module; Units: ns; 0=trigger_time+8500;
    Float_t  start_time_floor;       // cluster start time from clustering module; Units: ns; 0=trigger_time+8500;
    Float_t  charge_sss;             // photoelectron charge
    Float_t  charge_floor;           // photoelectron charge
    Bool_t   is_tracked;             // track present
    Bool_t   has_raw;                // = true if raw hits are written in file
    Bool_t   has_decoded;            // = true if decoded hits are written in file
    Int_t    has_clustered;          // = 1 if clusters are written; = 2 if clustered hits are written in file
    std::vector<BxMuonRawHit>        raw_hits;       //->
    std::vector<BxMuonDecodedHit>    decoded_hits;   //->
    std::vector<BxMuonCluster>       clusters;       //->
    //std::vector<BxMuonClusteredHit>  clustered_hits; //->
    BxTrackByPoints                  track;

    ClassDef (BxMuon, CYCLE_NUMBER)
};


/*********** BxMcTruth *************/
class BxMcTruthHit {
  public:
    BxMcTruthHit () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxMcTruthHit (const Int_t frame, const bx_mctruth_hit&);
#endif
    virtual ~BxMcTruthHit () {};
 
    // getters
    Int_t    GetNumFrame () const { return num_frame; }
    UShort_t GetLg       () const { return lg;        }
    Float_t  GetTime     () const { return time;      }

    virtual Bool_t IsFolder () const { return kTRUE; }

  private:
    Int_t    num_frame; // frame the hit belongs to in frames array
    UShort_t lg;        // logical channel (1 based)
    Float_t  time;      // tracking time (ns)

    ClassDef (BxMcTruthHit, CYCLE_NUMBER)
};

class BxMcTruthDaughter {
  public:
    BxMcTruthDaughter () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxMcTruthDaughter ( Int_t frame, const bx_mctruth_daughter&);
#endif
    virtual ~BxMcTruthDaughter () {};
 
    // getters
    Int_t    GetNumFrame           () const { return num_frame;    }
    Int_t    GetId                 () const { return id;           }
    Int_t    GetPdg                () const { return pdg;          }
    Double_t GetTime               () const { return time;         }
    Float_t  GetEnergy             () const { return energy;       }
    Float_t  GetPosition  ( Int_t i ) const { return position[i];  }
    Float_t  GetDirection ( Int_t i ) const { return direction[i]; }

    virtual Bool_t IsFolder () const { return kTRUE; }

  private:
    Int_t    num_frame;  // frame the daughter belongs to in frames array
    Int_t    id; 
    Int_t    pdg;
    Double_t time;
    Float_t  energy;
    Float_t  position[3];
    Float_t  direction[3];

    ClassDef (BxMcTruthDaughter, CYCLE_NUMBER)
};


class BxMcTruthDeposit {
  public:
    BxMcTruthDeposit () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxMcTruthDeposit ( Int_t frame, const bx_mctruth_deposit&);
#endif
    virtual ~BxMcTruthDeposit () {}

    // getters
    Int_t   GetNumFrame         () const { return num_frame;     }
    Int_t   GetPdgParent        () const { return pdg_parent;    }
    Float_t GetEnergy          () const { return energy;        } 
    Float_t GetPosition ( int i ) const { return position[i]; } 

    virtual Bool_t IsFolder () const { return kTRUE; }

  private:
    Int_t   num_frame; // frame the deposit belongs to in frames array
    Int_t   pdg_parent;
    Float_t energy;
    Float_t position[3];

    ClassDef (BxMcTruthDeposit, CYCLE_NUMBER)
};

class BxMcTruthUser {
  public:
    BxMcTruthUser () {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxMcTruthUser ( Int_t frame, const bx_mctruth_user&);
#endif
    virtual ~BxMcTruthUser () {}

    // getters
    Int_t    GetNumFrame () const { return num_frame; }
    Int_t    GetInt1     () const { return int1;      }
    Int_t    GetInt2     () const { return int2;      }
    Float_t  GetFloat1   () const { return float1;    } 
    Float_t  GetFloat2   () const { return float2;    } 
    Double_t GetDouble   () const { return double1;   } 
	
    virtual Bool_t IsFolder () const { return kTRUE; }

  private:
    Int_t    num_frame; // frame it refers to in frames array
    Int_t    int1;      // free for users' usage
    Int_t    int2;      // free for users' usage
    Float_t  float1;    // free for users' usage
    Float_t  float2;    // free for users' usage
    Double_t double1;   // free for users' usage

    ClassDef (BxMcTruthUser, CYCLE_NUMBER)
};


class BxMcTruthFrame {
  public:
    BxMcTruthFrame () { }
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxMcTruthFrame (const bx_mctruth_frame&);
#endif
    virtual ~BxMcTruthFrame () {};
 
    // getters
    UShort_t GetFileId          () const { return file_id;         }
    Double_t GetElecEventTime   () const { return elec_event_time; }
    Int_t    GetEventId         () const { return event_id;        }
    Int_t    GetNSequence       () const { return n_sequence;      }
    Int_t    GetIsotopeCoinc    () const { return isotope_coinc;   }
    Int_t    GetPdg             () const { return pdg;             }
    Double_t GetTime            () const { return time;            }
    Float_t  GetEnergy          () const { return energy;          }
    Float_t  GetVisibleEnergy   () const { return visible_energy;  }
    Float_t  GetPosition   (int i) const { return position  [i];   }
    Float_t  GetBaricenter (int i) const { return baricenter[i];   }
    Float_t  GetDirection  (int i) const { return direction [i];   }
    Int_t    GetIdNpe           () const { return id_npe;          } 
    Int_t    GetOdNpe           () const { return od_npe;          } 
    Int_t    GetNDaughters      () const { return n_daughters;     } 
    Int_t    GetNDeposits       () const { return n_deposits;      } 
    Int_t    GetNUsers          () const { return n_users;         } 
    Int_t    GetNIdPhotons      () const { return n_id_photons;    } 
    Int_t    GetNOdPhotons      () const { return n_od_photons;    } 
    
    virtual Bool_t IsFolder () const { return kTRUE; }

  private:
    UShort_t file_id;
    Double_t elec_event_time;    // relative to trg time
    Int_t    event_id;           // Event ID
    Int_t    n_sequence;         // Sequence number of isotope in the chain
    Int_t    isotope_coinc;      // 1 if it's a chain
    Int_t    pdg;                // Particle Data Group Code
    Double_t time;               // Time
    Float_t  energy;             // Kinetic energy
    Float_t  visible_energy;     // Quenched energy
    Float_t  position  [3];      // Position
    Float_t  baricenter[3];      // Baricenter
    Float_t  direction [3];      // Direction
    Int_t    id_npe;             // Number of photoelectrons in the ID (size of the vector hit_id_v)
    Int_t    od_npe;             // Number of photoelectrons in the OD (size of the vector hit_od_v)
    Int_t    n_daughters;        // Number of daughters (size of the vector daughters_v)
    Int_t    n_deposits;         // Number of deposits  (size of the vector deposits_v)
    Int_t    n_users;            // Size of the vector users_v
    Int_t    n_id_photons;       // Number of generated photons in the ID
    Int_t    n_od_photons;       // Number of generated photons in the OD

    ClassDef (BxMcTruthFrame, CYCLE_NUMBER)
};

class BxMcTruth {
  public:
    BxMcTruth() {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxMcTruth (const bx_write_opts&);
    virtual void operator=(const bx_mctruth_event& e);
#endif
    virtual ~BxMcTruth () {};

    // getters
    Int_t  GetNHitsId    () const {return n_hits_id;   }
    Int_t  GetNHitsOd    () const {return n_hits_od;   }
    Int_t  GetNFrames    () const {return n_frames;    }
    Int_t  GetNDaughters () const {return n_daughters; }
    Int_t  GetNDeposits  () const {return n_deposits;  }
    Int_t  GetNUsers     () const {return n_users;     }

    Bool_t HasFrames  () const { return Bool_t(write_flag > 0); }
    Bool_t HasHits    () const { return Bool_t(write_flag > 1); }

    const BxMcTruthFrame& GetFrame            ( Int_t i ) const { return frames[i]; }
    const std::vector<BxMcTruthFrame>&    GetFrames    () const { return frames;    }
    const std::vector<BxMcTruthDaughter>& GetDaughters () const { return daughters; }
    const std::vector<BxMcTruthDeposit>&  GetDeposits  () const { return deposits;  }
    const std::vector<BxMcTruthUser>&     GetUsers     () const { return users;     }
    const std::vector<BxMcTruthHit>&      GetHitsId    () const { return hits_id;   }
    const std::vector<BxMcTruthHit>&      GetHitsOd    () const { return hits_od;   }

    virtual Bool_t IsFolder() const { return kTRUE; }

  private:
    Int_t n_hits_id;   // number hits       (even if vector is not written)
    Int_t n_hits_od;   // number hits       (even if vector is not written)
    Int_t n_daughters; // number daughters  (even if vector is not written)
    Int_t n_deposits;  // number deposists  (even if vector is not written)
    Int_t n_users;     // number users      (even if vector is not written)
    Int_t n_frames;    // number frames     (even if vector is not written)
    Int_t write_flag;  // =1 if frames are written in file; =2 if also hits, daughters, deposits and users are written
    std::vector<BxMcTruthFrame>    frames;    //->
    std::vector<BxMcTruthDaughter> daughters; //->
    std::vector<BxMcTruthDeposit>  deposits;  //->
    std::vector<BxMcTruthUser>     users;     //->
    std::vector<BxMcTruthHit>      hits_id;   //->
    std::vector<BxMcTruthHit>      hits_od;   //->

    ClassDef (BxMcTruth, CYCLE_NUMBER)
};




/*********** BxNeutron *************/
class BxNeutronPulse  {
  public:
    BxNeutronPulse() {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxNeutronPulse (const bx_neutron_pulse&);
#endif
    virtual ~BxNeutronPulse() {}

    Float_t GetCharge    () const { return charge    ; }
    Float_t GetAmplitude () const { return amplitude ; }
    Float_t GetPeakTime  () const { return peak_time ; }
    Float_t GetRiseTime  () const { return rise_time ; }
    Float_t GetFallTime  () const { return fall_time ; }
    Float_t GetX         () const { return x         ; }
    Float_t GetY         () const { return y         ; }
    Float_t GetZ         () const { return z         ; }
    Float_t GetDX        () const { return dx        ; }
    Float_t GetDY        () const { return dy        ; }
    Float_t GetDZ        () const { return dz        ; }

  private:
    Float_t charge;    // charge as pulse area
    Float_t amplitude; // charge as peak amplitude
    Float_t peak_time; // peak time position in ns after muon time
    Float_t rise_time; // peak standard deviation; ns
    Float_t fall_time; // 1/tau from the exponential fit; ns
    Float_t x;         // position of pulse, x
    Float_t y;         // position of pulse, y
    Float_t z;         // position of pulse, z
    Float_t dx;        // position of pulse, error on x
    Float_t dy;        // position of pulse, error on y
    Float_t dz;        // position of pulse, error on z

  ClassDef(BxNeutronPulse, CYCLE_NUMBER)
};


class BxNeutron {
  public:
    BxNeutron() {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    virtual void operator=(const bx_neutron_event&);
#endif
    virtual ~BxNeutron () {};

    // getters
    Bool_t IsEnabled     () const { return is_enabled;    }
    Bool_t IsAssociated  () const { return is_associated; }
    Int_t  GetNPulses    () const { return pulses.size(); }
    Int_t  GetNNeutrons  () const { return n_neutrons;    }
    const std::vector<BxNeutronPulse>& GetNeutronPulses  () const { return pulses; }
    const BxNeutronPulse& GetNeutronPulse (Int_t i) const { return pulses.at(i);}

  private:
    Bool_t is_enabled;                         // file found
    Bool_t is_associated;                      // trgid/time associated
    Int_t n_neutrons;                          // number of neutron candidates
    std::vector<BxNeutronPulse> pulses;        // ->

  ClassDef (BxNeutron, CYCLE_NUMBER)
};





/*********** BxPhysTags *************/
class BxPhysTags {
   public:
     BxPhysTags ();
     enum tag_id {
         // Firsts 5 are reserved
       good		= 0,
       anytag		= 1,
       debug		= 4,
	 // Fast and other services 5 - 9
       fast		= 5,
       fast_large	= 6,
         // Muon stuff 10 - 24
       muon		= 10,
       muon_daughter	= 11,
       muon_1		= 12,
       muon_daughter_1  = 13,
       muon_2		= 14,
       muon_daughter_2	= 15,
       muon_cngs	= 16,
       muon_daughter_cngs = 17,
       muon_1_strict	= 18,
       muon_1_special   = 19,
         // Neutron and Cosmogenics 25 - 34
       C11  		= 25,
       C11_strict  	= 26,
       C10		= 27,
       LiHe_promt	= 28,
       LiHe_delayed	= 29,
         // Short coincidences 35 - 49
       Bi212		= 35,
       Po212		= 36,
       Bi212_large	= 37,
       Po212_large	= 38,
       Bi214		= 39,
       Po214		= 40,
       Bi214_large	= 41,
       Po214_large	= 42,
       Kr85 		= 43,
       Rb85		= 44,
       Kr85_large	= 45,
       Rb85_large	= 46,
         // Other stuff 54 - 63
       Pb214		= 54,
       Pb214_large	= 55,
       nubar_prompt	= 61,
       nubar_delayed	= 62,
       neutron		= 63,
         // MAX
       max 		= 64 // never >=
     };
     bool get_tag (tag_id id) const { return tags & (1ULL << id); }

     enum user_tag_id {
       user_max		= 32 // never >=
     };
     bool get_tag_user (user_tag_id id) const { return tags_user & (1 << id); }

     enum solar_tag_id {
       solar_max	= 32 // never >=
     };
     bool get_tag_solar (solar_tag_id id) const { return tags_solar & (1 << id); }

     enum cno_tag_id {
       cno_max	= 32 // never >= 32
     };
     bool get_tag_cno(cno_tag_id id) const { return tags_cno & (1 << id); }

     float get_var_user (int i) const {
       switch (i) {
	 case 0: return var_user_0;
	 case 1: return var_user_1;
	 case 2: return var_user_2;
	 case 3: return var_user_3;
	 case 4: return var_user_4;
	 case 5: return var_user_5;
	 case 6: return var_user_6;
	 case 7: return var_user_7;
	 case 8: return var_user_8;
	 case 9: return var_user_9;
	 case 10: return var_user_10;
	 case 11: return var_user_11;
	 case 12: return var_user_12;
	 case 13: return var_user_13;
	 case 14: return var_user_14;
	 case 15: return var_user_15;
	 case 16: return var_user_16;
	 case 17: return var_user_17;
	 case 18: return var_user_18;
	 case 19: return var_user_19;
       }
       return 0;
     }

       // Coincidence values
     int get_coincidence_id () const { return coincidence_id; }
     int get_coincidence_index () const { return coincidence_index; }
     int get_nsiblings () const  { return nsiblings; }
     float get_coincidence_dt () const { return coincidence_dt; }
     float get_coincidence_duration () const { return coincidence_duration; }
     float get_coincidence_distance () const { return coincidence_distance; }
   private:
     void set_tag (tag_id id, bool state) { state ? tags |= (1ULL << id) : tags &= ~(1ULL << id); }
     void set_tag_user (user_tag_id id, bool state) { state ? tags_user |= (1 << id) : tags_user &= ~(1 << id); }
     void set_tag_solar (solar_tag_id id, bool state) { state ? tags_solar |= (1 << id) : tags_solar &= ~(1 << id); }
     void set_tag_cno (cno_tag_id id, bool state) { state ? tags_cno |= (1 << id) : tags_cno &= ~(1 << id); }
     void set_mlp_ab     (float v) { mlp_ab      = v; }
     void set_mlp_ab_m4  (float v) { mlp_ab_m4   = v; }
     void set_bdt_c11b   (float v) { bdt_c11b    = v; }
     void set_bdt_c11b_m4(float v) { bdt_c11b_m4 = v; }
     void set_coincidence_id (int v) { coincidence_id = v; }
     void set_coincidence_index (int v) { coincidence_index = v; }
     void set_nsiblings (int v) { nsiblings = v; }
     void set_coincidence_dt (float v) { coincidence_dt = v; }
     void set_coincidence_duration (float v) { coincidence_duration = v; }
     void set_coincidence_distance (float v) { coincidence_distance = v; }

     void set_var_user (int i, float v) { 
       switch (i) {
	 case 0:  var_user_0 = v; return;
	 case 1:  var_user_1 = v; return;
	 case 2:  var_user_2 = v; return;
	 case 3:  var_user_3 = v; return;
	 case 4:  var_user_4 = v; return;
	 case 5:  var_user_5 = v; return;
	 case 6:  var_user_6 = v; return;
	 case 7:  var_user_7 = v; return;
	 case 8:  var_user_8 = v; return;
	 case 9:  var_user_9 = v; return;
	 case 10: var_user_10 = v; return;
	 case 11: var_user_11 = v; return;
	 case 12: var_user_12 = v; return;
	 case 13: var_user_13 = v; return;
	 case 14: var_user_14 = v; return;
	 case 15: var_user_15 = v; return;
	 case 16: var_user_16 = v; return;
	 case 17: var_user_17 = v; return;
	 case 18: var_user_18 = v; return;
	 case 19: var_user_19 = v; return;
       }
     }

     float var_user_0, var_user_1, var_user_2, var_user_3, var_user_4, var_user_5;
     float var_user_6, var_user_7, var_user_8, var_user_9, var_user_10, var_user_11, var_user_12;
     float var_user_13, var_user_14, var_user_15, var_user_16, var_user_17, var_user_18, var_user_19;
     float mlp_ab, mlp_ab_m4, bdt_c11b, bdt_c11b_m4;
     unsigned long long tags;
     unsigned long tags_user;
     unsigned long tags_solar;
     unsigned long tags_cno;
     short int coincidence_id, coincidence_index, nsiblings;
     float coincidence_dt, coincidence_distance, coincidence_duration;
     friend class bx_candidate;
     ClassDef (BxPhysTags, CYCLE_NUMBER)    
};




/*********** BxFwfd ****************/

class BxFwfdCluster  {
  public:
    BxFwfdCluster();
    BxFwfdCluster(Int_t, Double_t, Float_t, Float_t, Float_t, Float_t, Float_t, Int_t, Float_t, Float_t, Float_t);
    virtual ~BxFwfdCluster() {}

    Int_t    GetPeakPos    () const { return peak_pos; }
    Double_t GetTimePrev   () const { return time_prev; }
    Float_t  GetPeakAmpl   () const { return peak_ampl; }
    Float_t  GetASumCharge () const { return asum_charge; }
    Float_t  GetAttenASum  () const { return atten_asum; }
    Float_t  GetDSumCharge () const { return dsum_charge; }
    Float_t  GetDSumChargeCorr () const { return dsum_charge_corr; }
    Int_t    GetNumEchCluster() const { return num_ech_cluster; }
    Float_t  GetX          () const { return x; }
    Float_t  GetY          () const { return y; }
    Float_t  GetZ          () const { return z; }
    Int_t    GetMuonTag    () const { return muon_tag; }
    Int_t    GetNoiseTag   () const { return noise_tag; }
    void     SetMuonTag    (Int_t tag) { muon_tag = tag; }
    void     SetNoiseTag   (Int_t tag) { noise_tag = tag; }
    Float_t  GetUser1      () const { return user1; }
    Float_t  GetUser2      () const { return user2; }
    void     SetUser1      (Float_t user) { user1 = user; }
    void     SetUser2      (Float_t user) { user2 = user; }

  private:
    Int_t    peak_pos;    // position of the cluster's peak (0 - 511)
    Double_t time_prev;   // time after previous event (ms), related to peak_pos
    Float_t  peak_ampl;   // amplitude of the cluster's peak
    Float_t  asum_charge; // charge of the analog sum
    Float_t  atten_asum;  // charge of the attenuated (factor 3.2) analog sum
    Float_t  dsum_charge; // charge of the digital sum
    Float_t  dsum_charge_corr;  // dsum_charge corrected for saturated channels
    Int_t    num_ech_cluster;   // index number of Echidna cluster corresponding to FADC cluster
    Float_t  x;                 // x coordinate, m.
    Float_t  y;                 // y coordinate, m.
    Float_t  z;                 // z coordinate, m.
    Int_t    muon_tag;          // muon tagging flag
    Int_t    noise_tag;         // noise tagging flag
    Float_t  user1;
    Float_t  user2;

  ClassDef(BxFwfdCluster, CYCLE_NUMBER)
};


class BxFwfd {
  public:
    BxFwfd();
    BxFwfd(Bool_t, Bool_t, Int_t, Int_t, ULong_t, ULong_t, Int_t, Int_t, Int_t, UChar_t, Double_t, Int_t, UInt_t, UChar_t*);

    virtual ~BxFwfd() {}

    Bool_t  IsPresent          () const { return is_present; }
    Bool_t  IsODSumPresent     () const { return is_odsum_present; }
    Int_t   GetNFwfdEvs        () const { return n_fwfd_evs; }
    Int_t   GetUnixTime        () const { return unix_time; }
    ULong_t GetGpsTimeSec      () const { return gpstimes[0]; }
    ULong_t GetGpsTimeNs       () const { return gpstimes[1]; }
    Int_t   GetRun             () const { return run; }
    Int_t   GetEvNum           () const { return evnum; }
    Int_t   GetEvNumBx         () const { return evnum_bx; }
    UChar_t GetError           () const { return error; }
    Double_t  GetRawTime         () const { return raw_time; }
    Int_t   GetTrgType         () const { return trgtype; }
    UInt_t  GetNClusters       () const { return n_clusters; }
    UChar_t GetDCode      (int i) const { return dcode[i]; }
    const Float_t& GetWFormAsum (int i) const { return wform_asum.at(i); }
    const Float_t& GetWFormDsum (int i) const { return wform_dsum.at(i); }
    const Float_t& GetWFormODsum (int i) const { return wform_odsum.at(i); }
    const std::vector<Float_t>& GetWFormA () const { return wform_asum; }
    const std::vector<Float_t>& GetWFormD () const { return wform_dsum; }
    const std::vector<Float_t>& GetWFormOD () const { return wform_odsum; }

    const std::vector<BxFwfdCluster>& GetClusters  () const { return clusters; }
    const BxFwfdCluster& GetCluster (int i) const { return clusters.at(i);}

    void SetStatus(bool v)	{ is_present = v; }
    void SetStatusODsum(bool v)	{ is_odsum_present = v; }
    void SetNFwfdEvs(int v)     { n_fwfd_evs = v; }
    void SetUnixTime(int v)	{ unix_time = v; }
    void SetGpsTimeSec(unsigned long v) { gpstimes[0] = v; }
    void SetGpsTimeNs(unsigned long v)  { gpstimes[1] = v; }
    void SetRun(int v)          { run = v; }
    void SetEvNum(int v)        { evnum = v; }
    void SetEvNumBx(int v)      { evnum_bx = v; }
    void SetError(unsigned char v)	{ error = v; }
    void SetRawTime(double v)	{ raw_time = v; }
    void SetTrgType(int v)      { trgtype = v; }
    void SetNClusters(unsigned int v)  { n_clusters = v; } 
    void SetDCode(int channel, unsigned char v)  { dcode[channel] = v; }
    void SetWFormAsum(std::vector<float> v)  { wform_asum = v; }
    void SetWFormDsum(std::vector<float> v)  { wform_dsum = v; }
    void SetWFormODsum(std::vector<float> v)  { wform_odsum = v; }

    void SetCluster(UInt_t n_clusters, UInt_t num_cluster, Int_t v1, Double_t v2, Float_t v3, Float_t v4, Float_t v5, Float_t v6, Float_t v7, Int_t v8, Float_t v9, Float_t v10, Float_t v11);
    void ClearClusters();
    void ClearWForms();

  private:
    Bool_t  is_present;      // FWFD data present or not
    Bool_t  is_odsum_present;      // OD analog sum present or not in FADC data
    Int_t   n_fwfd_evs;      // number of FWFD events, corresponding to one BX event - to be used in dsts only
    Int_t   unix_time;       // UNIX time, seconds since 01/01/1970 00:00
    ULong_t gpstimes[2];     // GPS time, [0] seconds since 01.01.00 0:0:0, [1] nanoseconds
    Int_t   run;             // FWFD run (normally coincides with main BX run)
    Int_t   evnum;           // FWFD event number
    Int_t   evnum_bx;        // BX event number as received from BX pattern
    UChar_t error;           // =0 - no error, =1 - memory overflow
    Double_t  raw_time;      // hardware time to previous trigger (ns)
    Int_t   trgtype;         // FWFD trigger type:
                             // 1 - neutrino
                             // 3 - 100 Hz pulser in the beginning of the run
                             // 4 - laser
                             // 5 - laser
                             // 8 - muon
                             // 9 - muon
                             // 11 - pulser 0.1 Hz during the run
    UInt_t n_clusters;       // number of clusters
    UChar_t dcode[16];       // flag on discriminator's channel: =0 - channel slept, =1 - channel fired
    std::vector<Float_t> wform_asum; // waveform of the analog sum
    std::vector<Float_t> wform_dsum; // waveform of the digital sum (deconvoluted)
    std::vector<Float_t> wform_odsum; // waveform of the OD analog sum

    std::vector<BxFwfdCluster> clusters;        // ->

  ClassDef(BxFwfd, CYCLE_NUMBER)
};





/*********** BxEvent *************/
class BxEvent : public TObject  {
  public:
    BxEvent() {}
#if !defined(_ECHIDNA_ROOTLIB_) && !defined(__CINT__)
    BxEvent(const bx_write_opts&);
    virtual void operator=(const bx_echidna_event& ev);
#endif
    virtual ~BxEvent () {};

    // event header getters
    Int_t GetRun   () const { return run;   }
    Int_t GetEvNum () const { return evnum; }
    Bool_t IsCrateEnabled ( int i ) const { return (enabled_crates >> (i-1)) & 1;}
    UShort_t GetEnabledCrates () const { return enabled_crates; }
    Bool_t IsMuonEnabled      () const { return (enabled_crates >> 14) & 1; }
    Bool_t IsMuonAligned      () const { return IsMuonEnabled() && muon.IsAligned (); }
    Bool_t IsMcTruthEnabled   () const { return mctruth.GetNFrames() > 0; }
    //Bool_t IsNeutronEnabled   () const { return neutron.IsAssociated(); }
    Bool_t IsTrackedGlobal    () const { return is_tracked_global;}
    //Bool_t IsTrackedCmt       () const { return is_tracked_cmt;}

    // subclass getters
    const BxTrigger&     GetTrigger       () const { return trigger; }
    const BxLaben&       GetLaben         () const { return laben;   }
    const BxMuon&        GetMuon          () const { return muon;    }
    const BxMcTruth&     GetMcTruth       () const { return mctruth; }
    //const BxNeutron&     GetNeutron       () const { return neutron; }
    const BxTrackFitted& GetTrackGlobal   () const { return track_global;   }
    //const BxTrackByPoints& GetTrackCmt      () const { return track_cmt;   }
    const BxPhysTags&    GetPhysTags (int i) const { return tags[i]; }
    int   		 GetNM4Events     () const { return m4s.size (); }
    const Mach4Event&    GetMach4Event (int i) const { return m4s[i]; }
    const BxFwfd&        GetFwfd          () const { return fwfd; }

    // Calculate time difference as Dt = this - prev. Units are us, with a precision
    // that can arrive to ns if |Dt| < 3ms. The precision of ~0.3us is from GPS absolute
    // time and is kept up to 10 years considering the double precision (16 significant
    // digits) -if you trust the atomic time for so long-.
    double GetTimeDifference (const BxEvent* prev) const { return GetTimeDifference (*prev); }
    double GetTimeDifference (const BxEvent& prev) const { return GetTimeDifference (prev.GetTrigger ().GetGpsTimes (), prev.GetLaben ().GetTriggerTime ()); }
    double GetTimeDifference (const ULong_t* prev_gps_times, Double_t prev_laben_trigger_time) const;
    // Calculate the time difference between 2 clusters, as
    // Dt = this_event.current_cluster - prev_event.prev_cluster. 
    // Same timing arguments of previus methods.
    double GetTimeDifference (int current_cluster, const BxEvent* prev_event, int prev_cluster) const { return GetTimeDifference (current_cluster, *prev_event, prev_cluster); }
    double GetTimeDifference (int current_cluster, const BxEvent& prev_event, int prev_cluster) const;
    double GetTimeDifference (int current_cluster, const ULong_t* prev_gps_times, Double_t prev_laben_cluster_time) const;

    virtual Bool_t IsFolder() const { return kTRUE; }

  private:
    Int_t         run;             // run number
    Int_t         evnum;           // event number (trigger ID)
    UShort_t      enabled_crates;  // bitfield: 0=trigger, 1-14=laben, 15=muon
    //Bool_t        is_tracked_cmt;    // the cmt tracker worked fine
    Bool_t        is_tracked_global; // the global tracker worked fine

    BxTrigger     trigger;         // sub-object with trigger info
    BxLaben       laben;           // sub-object with laben info
    BxMuon        muon;            // sub-object with muon info
    BxMcTruth     mctruth;         // sub-object with mc truth info
    //BxNeutron     neutron;         // sub-object with V1731 flash adc boards info used for neutron detection
    //BxTrackByPoints track_cmt;     // muon track from the CMT
    BxTrackFitted   track_global;    // muon track from the global fitter
    std::vector<BxPhysTags> tags;  // written by physics tools
    std::vector<Mach4Event> m4s;   // written by M4_merger
    int n_m4;                      // written by M4_merger
    BxFwfd        fwfd;            // written by physics tools

    friend class bx_candidate;
    friend void merge_events (const BxEvent*, const Mach4Event*);
    friend void fmerge_fake_bxevent(BxEvent* ev, Int_t run, Int_t evnum, UChar_t trgtype, Int_t nclu);
    ClassDef (BxEvent, CYCLE_NUMBER * 100 + M4VERSION)    
};



#endif

/*
 * $Log: BxEvent.hh,v $
 * Revision 1.235  2015/08/26 11:22:04  misiaszek
 * new rms/kurtosis for c11 psd
 *
 * Revision 1.234  2015/08/25 12:28:26  misiaszek
 * mlp_ab variable added
 *
 * Revision 1.233  2015/08/19 03:14:25  koun
 * bug correction in Normalize_geo/QE_Charge
 *
 * Revision 1.232  2015/08/16 21:12:37  misiaszek
 * cvs email test as puck.lngs.infn.it smtp server is not connecting I changed it to gmail smtp
 *
 * Revision 1.231  2015/08/16 18:58:41  misiaszek
 * test of cvs mail
 *
 * Revision 1.230  2015/08/06 12:25:22  koun
 * added Normalize_geo/QE_Pmts/Charge variables
 *
 * Revision 1.229  2015/08/06 09:26:15  misiaszek
 * test of cvs mails
 *
 * Revision 1.228  2015/08/03 15:42:15  misiaszek
 * cno tags added
 *
 * Revision 1.227  2015/07/29 09:23:12  ilia.drachnev
 * added position_reco_noavg variable
 *
 * Revision 1.226  2015/07/28 15:22:26  misiaszek
 * set functions and new mva variables for a/b & c11/b discrimination added
 *
 * Revision 1.225  2015/07/28 09:14:54  misiaszek
 * tailtot for c11/b discrimination added
 *
 * Revision 1.224  2015/07/14 13:55:00  lukyanch
 * Remove inheritance from TObject for BxFwfd and add two user variables to BxFwfdCluster
 *
 * Revision 1.223  2015/07/14 09:54:46  misiaszek
 * charge_noavg_dt1,  charge_noavg_dt2,  charge_noavg , charge_noavg_short added
 *
 * Revision 1.222  2015/07/13 17:36:21  misiaszek
 * var_user variables & solar tags added
 *
 * Revision 1.221  2015/07/13 15:06:44  misiaszek
 * tailtot_ab_mlp for a/b discrimination with MLP algorithm added
 *
 * Revision 1.220  2015/07/02 15:35:25  lukyanch
 * Add muon and noise tags to FADC event
 *
 * Revision 1.219  2015/01/09 15:03:07  misiaszek
 * cycle_18 new unstable
 *
 * Revision 1.218  2015/01/02 02:05:59  misiaszek
 * npmts, nhits & charge postition corrected energy variables added to BxLabenCluster
 *
 * Revision 1.217  2014/12/12 20:23:07  misiaszek
 * charge_win1/win2 added with getters
 *
 * Revision 1.216  2014/12/12 20:11:43  misiaszek
 * track_cmt removed
 *
 * Revision 1.215  2014/12/11 21:27:12  wurm
 * added getters and variables to BxEvent that concern C-11 likelihood tagging
 * * for tracks: new distance-to- and projection-on-track getters, dE/dx calculation
 * * counting of neutron-like clusters in tt128, new flag of neutron clusters
 *
 * Revision 1.214  2014/12/09 16:46:32  misiaszek
 * TObject ineritance restored for BxEvent & BxFwfd
 *
 * Revision 1.213  2014/12/07 16:19:44  misiaszek
 * Many not used variables removed
 *
 * Revision 1.212  2014/12/05 17:24:34  misiaszek
 * position_mi/_dbn/_msk/_mach4/_mach4_fixed and energy_mc/_lik/_msk/_dbn removed
 *
 * Revision 1.211  2014/12/05 16:06:16  misiaszek
 * TObject inheritance removed as it is not needed to save objects to TTree. It was done to remove fUniqueID and fBits from bxtree
 *
 * Revision 1.210  2014/12/03 17:07:53  misiaszek
 * charge_dt1 and charge_dt2 to clusters added
 *
 * Revision 1.209  2014/11/25 17:50:59  lukyanch
 * Add friend function allowing FADC merger to create fake BxEvents
 *
 * Revision 1.208  2014/11/07 12:54:44  litvinov
 * update in FADC stuff: (1) num_ech_cluster and (2) corrected dsum_charge added in clusters, (3) wform of the OD analog sum and (4) whether it is present or not added in base BxFwfd
 *
 * Revision 1.207  2013/06/18 18:56:37  razeto
 * cycle_17 new unstable
 *
 * Revision 1.206  2013-02-02 09:01:49  razeto
 * Incremented to cycle_16 (cycle 15 was lost)
 *
 * Revision 1.205  2013-01-29 12:32:28  litvinov
 * reconstructed x,y,z for FADC (spatial reco is being implemented in fadc reader)
 *
 * Revision 1.204  2013-01-22 09:22:17  mosteiro
 * dt1,2 vars calculated for each cluster
 *
 * Revision 1.203  2012-10-30 16:21:48  ddangelo
 * added vectors for npmts calculation in tt64
 * internal and root event with getters and copy.
 *
 * Revision 1.202  2012-10-22 15:56:04  ddangelo
 * added npmts_dt1, npmts_dt2 to laben (clustered) event.
 * internal and root, with getters and cpy operator.
 *
 * Revision 1.201  2012-10-18 14:15:16  litvinov
 * FADC raw_time is being double, ulong is not enough at 32bits (bxmaster)
 *
 * Revision 1.200  2012-10-17 18:24:57  litvinov
 * FADC raw_time now ULong_t (was UInt_t)
 *
 * Revision 1.199  2011-05-31 15:30:34  razeto
 * muon stuff committed back
 *
 * Revision 1.198  2011-05-31 15:27:27  razeto
 * New LiHe tag
 *
 * Revision 1.197  2011-05-31 15:10:36  razeto
 * Added kr85 large
 *
 * Revision 1.196  2011-05-17 15:08:41  razeto
 * Redefined tags
 *
 * Revision 1.195  2011-04-19 05:54:58  razeto
 * Moved to cycle 15 unstable
 *
 * Revision 1.194  2011-04-17 13:23:27  meindl
 * Updated definition of the Getter "IsClusterOld"
 *
 * Revision 1.193  2011-04-13 15:42:11  ddangelo
 * added a getter to check if neutron cluster is above c13 variable threshold
 *
 * Revision 1.192  2011-04-13 10:37:35  ddangelo
 * added variable nclusters_old
 *
 * Revision 1.191  2011-03-23 16:04:48  ddangelo
 * added charge_clean, zeroing, copy. internal and root event
 *
 * Revision 1.190  2011-02-24 15:51:46  wurm
 * added entry/exit point getters for global tracking
 *
 * Revision 1.189  2011-02-22 08:11:27  razeto
 * GPS start from 2000 UTC + added leap seconds + timet
 *
 * Revision 1.188  2011-02-19 13:29:11  davini
 * bx_pid_positron stuffs
 *
 * Revision 1.187  2011-02-18 18:22:01  ddangelo
 * event support to mach4 a/b discrimination commented out. variable writing with the module commented out.
 * added event support for gatti ops. to be improved.
 *
 * Revision 1.186  2011-02-18 17:10:05  ddangelo
 * major code cleanup: removed fadc throughout the program
 *
 * Revision 1.185  2011-02-18 16:04:23  davini
 * added npmts_400 nhits_400 charge_400 on Oleg's request; added charge_npmts based on Alessandro, Stefano and Livia idea;
 *
 * Revision 1.184  2011-02-11 16:28:46  litvinov
 * FWFW raw_time unsigned instead of signed int
 *
 * Revision 1.183  2011-01-15 20:57:59  razeto
 * Fix indexes
 *
 * Revision 1.182  2010-12-10 10:32:43  litvinov
 * added FADC raw_time. Waveforms are now vectors
 *
 * Revision 1.181  2010-11-30 10:37:00  wurm
 * set getters to phi/theta computed by global_tracker
 *
 * Revision 1.180  2010-11-29 16:03:01  litvinov
 * added run to BxFwfd; n_clusters was UChar_t, now integer
 *
 * Revision 1.179  2010-08-26 08:25:20  wurm
 * fixed theta getter
 *
 * Revision 1.178  2010-08-06 17:20:16  razeto
 * Moving to cycle 14
 *
 * Revision 1.177  2010-08-03 15:58:05  wurm
 * introduced theta, phi, impact variables for fitted tracks
 *
 * Revision 1.176  2010-08-01 15:57:43  razeto
 * Added m4s.size() proxy getter
 *
 * Revision 1.175  2010-07-01 18:21:36  razeto
 * Added counter hit flag for new laben fw and nhits_fw from laben boards
 *
 * Revision 1.174  2010-07-01 13:28:25  ddangelo
 * nhits_bkg changed to float to accomodate info for tt1
 *
 * Revision 1.173  2010-06-30 10:58:01  litvinov
 * added variable n_fwfd_evs to BxFwfd
 *
 * Revision 1.172  2010-06-24 16:38:33  litvinov
 * added constructor with parameters for BxFwfd
 *
 * Revision 1.171  2010-06-17 14:44:22  litvinov
 * more FWFD stuff: default constructor for BxFwfdCluster, dtors, ClearClusters()
 *
 * Revision 1.170  2010-06-16 12:54:19  litvinov
 * added FWFD stuff. Commit authorized by ddangelo
 *
 * Revision 1.169  2010-05-21 15:55:14  ddangelo
 * different things on muon tracks
 *
 * 1.a) old laben_track renamed as laben_track_energy
 * new laben_track_tof added
 *
 * 1.b) (global) track renemed as track_global at base event level
 * track_cmt added at base event level (track by points)
 *
 * 1) all getters updated/integrated
 * is_tracked variable updated/integrated accordingly. inizialization.
 * job ported to root event as well. copy done.
 * friendship with old/new module updated
 *
 * 2) bxtrack_by_points class:
 * - theta, phi and impact added as variables.
 * - errors added on all of the above.
 * - error code variable requested by cmt tracker added
 *
 * Revision 1.168  2010-05-03 19:11:37  razeto
 * Friending to n_live_charge
 *
 * Revision 1.167  2010-04-01 14:06:08  ddangelo
 * re-introduced m4 merger variables, excluded by mistake
 *
 * Revision 1.166  2010-03-29 14:44:10  ddangelo
 * debugging
 *
 * Revision 1.162.2.1  2010-03-29 13:54:16  ddangelo
 * debugging
 *
 * Revision 1.162  2009-11-18 11:43:31  ddangelo
 * added rms_time, duration and npmts short version for back compatibility.
 * npmt variables renamed to npmts also in internal event.
 * mach4_n1700 renamed as mach4_fixed throughout the event classes.
 *
 * Revision 1.161  2009-10-26 19:17:23  ddangelo
 * added bx_position_mach4_n1700 class (internal and root event, copy and getters)
 * in (base)postion class variable likelihood renamed as user. parallel getters to use it as likelihood or refraction index
 *
 * Revision 1.160  2009-10-26 11:19:44  ddangelo
 * trgtype 'muon_tct' (or 'muon_totc') renamed as 'neutron' after real use
 * trgtype 'muon_mtb' renamed as 'muon' as duplicity with the above one disappeared and 'mtb' is associated to btb_inputs flag
 *
 * Revision 1.159  2009-10-23 16:42:22  ddangelo
 * cluster_neutrons renamed as clusters_muons
 *
 * Revision 1.158  2009-10-23 14:00:03  koshio
 * Add the lngs postion reconstruction
 *
 * Revision 1.157  2009-10-23 09:07:15  ddangelo
 * added a laben cluster vector for parallel neutron clustering in the muon gate. empty for non-muon events.
 *
 * Revision 1.156  2009-10-22 16:06:46  ddangelo
 * in bx_laben_cluster added nhits_bkg for evaluate of bkg contribution in case of piled-up events.
 * internal and root event, copy and getters
 *
 * Revision 1.155  2009-10-22 15:10:48  ddangelo
 * fixed a typo
 * totc renamed as neutron
 *
 * Revision 1.154  2009-10-08 15:45:59  ddangelo
 * implemented pid event variables removal and addition.
 * Internal and root, inizialization, copy and getters.
 *
 * Revision 1.153  2009-10-06 17:49:58  ddangelo
 * 1) laben decoded event: added n_invalid_charge (internal and root event, copy and getters)
 * 2) laben raw hit (root): getter IsInvalid now reads dedicated bit 0x80;
 * 3) NormalizePmts() and NormalizeCharge() methods now correct for invalid channels as well.
 *
 * Revision 1.152  2009-10-06 13:36:17  ddangelo
 * laben decoded event: added variable n_invalid_pmts to account for FE||FF condition only on ordinary && !disabled channels
 * internal and root event, getters and copy
 *
 * Revision 1.151  2009-10-05 15:18:16  ddangelo
 * added npmts at laben decoded event level (inner and outer, copy and getters);
 *
 * Revision 1.150  2009-10-05 14:43:07  ddangelo
 * added invalid flag getter in bx laben raw hit
 * corrected for correct array size raw_flags in raw event
 * added 7 getters for individual flag counters
 *
 * Revision 1.149  2009-09-18 22:20:12  ddangelo
 * added sphere_rel_var to laben shaped cluster. Internal and root event, inizialization and copy.
 *
 * Revision 1.148  2009-09-18 10:39:40  wurm
 * once again the phi getter for pointy tracks after looking at CNGS muons; uses atan2() now
 *
 * Revision 1.147  2009-09-17 14:28:51  ddangelo
 * in laben clustered/decoded (internal/root) hit added the flag short_cluster to say if hit belonged to the old (c11) cluster
 *
 * Revision 1.146  2009-09-17 13:57:46  ddangelo
 * in laben cluster added nhits, charge and mean_time '_short' for old clustering values
 *
 * Revision 1.145  2009-09-16 15:49:14  ddangelo
 * added getters for (BxLabenDecodedHits) num_cluster and rec_time
 *
 * Revision 1.144  2009-09-16 15:44:28  ddangelo
 * npe (and npe_conc) removed from BxLaben, BxLabenCluster and BxLabenDecodedHit classes (root event only)
 *
 * Revision 1.143  2009-09-16 15:36:52  ddangelo
 * BxLabenCluster::end_time renamed as duration
 * neutron, tags and raw_index branches disabled
 *
 * Revision 1.142  2009-09-10 11:47:50  razeto
 * Typo fixed (non maintainer trivial upload)
 *
 * Revision 1.141  2009-08-31 10:53:56  ddangelo
 * laben decoded hits: added num cluster and rec time
 *
 * Revision 1.140  2009-07-30 16:20:44  ddangelo
 * - removed IsTracked() getters from BxTrack classes (root event)
 * + added is_tracked variable and getter in BxEvent and BxLaben classes (root event)
 * + fixed the copy of is_tracked variable in BxMuon class (root event)
 * + added is_tracked variable and getter for bx_echidna event class (internal event). to be filled by global tracker
 *
 * Revision 1.139  2009-07-22 10:41:22  ddangelo
 * in position classes, both internal and root event:
 * - n_iterations removed
 * + matrix and converged variables added with getters.
 *
 * Revision 1.138  2009-07-17 17:46:17  ddangelo
 * In root event only:
 * charge (and charge_mean for now) copied also to clustered hits.
 * d80 and time_error commented out from decoded hits
 *
 * Revision 1.137  2009-07-17 15:39:51  ddangelo
 * laben cluster:
 * + added npmts_thresh, nhits_thresh and charge_thresh, to be computed with >0.2pe hits only (internal and root event).
 * - removed cluster rough time (internal and root).
 * - removed npe and npe_conc (root only)
 *
 * Revision 1.136  2009-07-17 14:44:08  ddangelo
 * ppc0 cpu time reintroduced in trg event (for cmp with GPS)
 *
 * Revision 1.135  2009-07-17 13:08:34  ddangelo
 * alpha beta varaibles reorganized, rise_time added
 *
 * Revision 1.134  2009-07-16 15:50:32  ddangelo
 * copyng data implemented
 *
 * Revision 1.133  2009-07-16 15:17:57  ddangelo
 * mach a/b ported to root event
 * debugging tailtot getter for m4
 * other debugging
 *
 * Revision 1.132  2009-07-16 10:17:42  ddangelo
 * infrastructure for m4 position reco (patch by steve&ale)
 *
 * Revision 1.131  2009-04-23 23:58:43  ddangelo
 * comments
 *
 * Revision 1.130  2009-04-20 14:01:14  ddangelo
 * IsTracked() method virtualized and inherited across the 3 track classes.
 *
 * Revision 1.129  2009-04-20 13:50:40  ddangelo
 * added errors on rec coordinates and time in class bx_position and BxPosition. initialization and copy included.
 *
 * Revision 1.128  2009-04-16 17:00:06  ddangelo
 * better comments
 *
 * Revision 1.127  2009-04-16 07:35:03  wurm
 * Changing condition for IsTracked() in fitted track
 *
 * Revision 1.126  2009-04-15 17:12:38  ddangelo
 * n_point changed into a bitfield, internal and root event (BxTrackFitted class)
 *
 * Revision 1.125  2009-04-15 16:54:25  ddangelo
 * GetTrack() returning derived classes for all 3 instances
 *
 * Revision 1.124  2009-04-14 18:04:33  ddangelo
 * minor upgrades
 *
 * Revision 1.123  2009-03-31 17:11:06  ddangelo
 * debugging
 *
 * Revision 1.122  2009-03-27 17:21:31  ddangelo
 * cosmetics
 *
 * Revision 1.121  2009-03-24 09:40:14  razeto
 * Added more muon tags
 *
 * Revision 1.120  2008-12-15 17:13:55  razeto
 * New cycle (12)
 *
 * Revision 1.119  2008-12-15 12:12:44  razeto
 * Added rms_time to cluster (to improve PID)
 *
 * Revision 1.118  2008-12-11 08:03:30  wurm
 * simplified get_theta for track_by_points
 *
 * Revision 1.117  2008-12-10 16:26:12  wurm
 * significant improvement of the fitted track impact getter
 *
 * Revision 1.116  2008-12-10 11:40:24  razeto
 * Added mean charge and peak from tt1
 *
 * Revision 1.115  2008-12-02 13:53:32  ddangelo
 * added n_points to fitted track class and getters. internal and root event
 *
 * Revision 1.114  2008-11-11 17:21:45  razeto
 * Reorganized tags
 *
 * Revision 1.113  2008-10-26 12:24:30  wurm
 * added IsTracked() and GetImpact() for fitted track
 *
 * Revision 1.112  2008-10-17 13:41:12  razeto
 * new development cycle (11)
 *
 * Revision 1.111  2008-10-01 16:44:48  wurm
 * changed getters for theta and phi of fitted track
 *
 * Revision 1.110  2008-10-01 16:17:04  ddangelo
 * removed is_pointlike variables and related stuff (bx_filters will perform this task)
 *
 * Revision 1.109  2008-09-24 12:15:58  razeto
 * Few more variables for phys tags
 *
 * Revision 1.108  2008-08-26 15:30:30  ddangelo
 * added neutron enabled and association flag, muon aligned flags
 *
 * Revision 1.107  2008-08-25 14:46:21  razeto
 * Fixed a warning
 *
 * Revision 1.106  2008-08-21 21:47:27  razeto
 * added midday
 *
 * Revision 1.105  2008-08-21 16:00:29  razeto
 * Added Mjd2TimeT
 *
 * Revision 1.104  2008-08-21 15:56:54  razeto
 * Reorganization and upgrades in sun code
 *
 * Revision 1.103  2008-08-08 12:48:03  razeto
 * New tag format redesign completed
 *
 * Revision 1.102  2008-08-05 16:30:52  ddangelo
 * added fall time, some variables renamed
 *
 * Revision 1.101  2008-08-05 15:30:50  ddangelo
 * added chi2
 *
 * Revision 1.100  2008-07-31 20:09:34  razeto
 * Position ctor for BxPosition
 *
 * Revision 1.99  2008-07-31 20:07:48  razeto
 * Reorganization of tags structure
 *
 * Revision 1.98  2008-07-30 09:50:42  razeto
 * Fixed Mjd calculation
 *
 * Revision 1.97  2008-07-29 18:28:19  razeto
 * Azimut and latitude experimental calculation added
 *
 * Revision 1.96  2008-07-17 15:55:00  ddangelo
 * added track direction, removed unused getter
 *
 * Revision 1.95  2008-07-11 17:06:27  ddangelo
 * added classes for neutron system (code by S. Davini)
 *
 * Revision 1.94  2008-07-11 14:28:14  ddangelo
 * fixed a typo
 *
 * Revision 1.93  2008-07-11 14:21:36  ddangelo
 * debugging
 *
 * Revision 1.92  2008-05-11 15:56:45  ddangelo
 * added theta and phi calculation for BxTrackFitted
 *
 * Revision 1.91  2008-05-10 23:34:05  ddangelo
 * BxDistance class defined and implemented
 * GetDistance and GetDisatnceError implemented in BxTrackFitted
 * some more work on track and position classes
 *
 * Revision 1.90  2008-05-09 15:46:23  ddangelo
 * muon track implemented as "by_points" and "fitted" classes inheriting from common interface (pure virtual class bx_track)
 * laben and muon host an istance "by_points"
 * echidna event hosts an istance "fitted"
 * structure duplicated at root level and filled
 * GetDistance duplicated in BxPosition class
 *
 * Revision 1.89  2008-05-08 19:06:36  razeto
 * Added get_time_t and GetTimeT from gps time
 *
 * Revision 1.88  2008-04-27 20:46:40  ddangelo
 * added 3 getters to BxTrack class
 *
 * Revision 1.87  2008-04-14 14:14:18  ddangelo
 * added BxLaben::GetTrack(), forgotten before
 *
 * Revision 1.86  2008-03-03 13:51:32  razeto
 * Added muon_2 experimental tag
 *
 * Revision 1.85  2008-02-27 20:46:13  razeto
 * new development cycle (10)
 *
 * Revision 1.84  2008-02-27 20:26:30  razeto
 * New clasdef(9) version and new cycle version
 *
 * Revision 1.83  2008-02-26 18:29:25  ddangelo
 * added is_pointlike/tracklike variable and getters in rec (shaped) cluster
 * both inner and root event
 *
 * Revision 1.82  2008-02-26 17:27:44  ddangelo
 * added
 * BxTrack::GetDistance(BxPosition&)
 * BxPosition::GetDistance(BxTrack&)
 * P.S.: thanks to yura
 *
 * Revision 1.81  2008-02-21 15:58:31  ddangelo
 * added muon track as a separate class.
 * Currently used twice: laben and muon (later also globally fit).
 * both inner and roor event
 * added tracked stage in base event enum
 * bx_muon_rec_event renamed as bx_muon_tracked_event
 *
 * Revision 1.80  2008-01-14 16:28:28  ddangelo
 * overloading NormalizePmts for Int_t type
 *
 * Revision 1.79  2007-12-20 18:41:29  ddangelo
 * handling of individual nhits for sss and floor
 * filling some varaibles left over before.
 * some more debugging
 *
 * Revision 1.78  2007-12-07 14:09:47  ddangelo
 * added n_live_charge (internal and root), getters and normalize method
 * renamed a variable
 *
 * Revision 1.77  2007-12-06 16:48:09  ddangelo
 * upgraded to new muon clustering.
 * muon tracking side debugged and improved.
 *
 * Revision 1.76  2007-11-29 15:42:19  litvinov
 * added muon_cngs tag
 *
 * Revision 1.75  2007-11-16 13:33:48  ddangelo
 * added getters for muon track
 * plain getters and smart getters (angles, pedal, impact)
 * fixed a typo on mctruth
 *
 * Revision 1.74  2007-11-14 19:02:11  ddangelo
 * added tracking variables to root event
 *
 * Revision 1.73  2007-11-14 17:07:55  ddangelo
 * new muon clustering variables.
 * indipendent up/floor clustered hits vectors
 * (internal and root event)
 * filling and inizialization. tested.
 *
 * Revision 1.72  2007-11-12 15:31:13  razeto
 * Added cluster flags and removed broad variable (davide auth)
 *
 * Revision 1.71  2007-11-07 14:48:17  ddangelo
 * debugging
 *
 * Revision 1.70  2007-11-07 13:32:47  ddangelo
 * debugging
 *
 * Revision 1.69  2007-11-06 14:45:10  ddangelo
 * debugging mctruth
 *
 * Revision 1.68  2007-11-05 23:39:03  razeto
 * Added new hits_on_empty variable to laben data
 *
 * Revision 1.67  2007-10-31 17:12:53  razeto
 * Added a new variable for clustering
 *
 * Revision 1.66  2007-10-31 15:42:14  ddangelo
 * added quenched energy (mctruth)
 * disk format, internal and root event
 *
 * Revision 1.65  2007-10-30 18:40:42  ddangelo
 * added multi condition flag to laben decoded hit
 * getter, smart getters, copy.
 *
 * Revision 1.64  2007-10-30 18:12:27  ddangelo
 * added empty_boards to root event too.
 * getters
 *
 * Revision 1.63  2007-10-30 15:45:02  ddangelo
 * added # of laben cluster found by algorythm (different from saved one for high multiplicity events)
 * added end hit time of a cluster
 * internal and root event. getters, copy, initialization, etc...
 *
 * Revision 1.62  2007-10-29 17:20:44  ddangelo
 * adding nphotons for msk energy reco
 *
 * Revision 1.60  2007-10-25 15:46:03  ddangelo
 * added a bit field to shaped event, meaning to be assigned. smart getters to be added at that time. internal and external event.
 * all variables zeroed in ctor for shaped event.
 *
 * Revision 1.59  2007-10-25 14:43:22  ddangelo
 * added a variable to rec clsuter
 *
 * Revision 1.58  2007-10-11 11:23:52  ddangelo
 * new mctruth format (internal AND root event)
 *
 * Revision 1.57  2007-10-11 10:49:54  razeto
 * Cycle 8 deployed
 *
 * Revision 1.56  2007-10-01 14:53:06  dfranco
 * added a method to evaluate the modified julian date for each event
 *
 * Revision 1.55  2007-09-28 16:03:20  dfranco
 * added C11 and C10 tags
 *
 * Revision 1.54  2007-09-28 15:27:04  razeto
 * Finally normalize at 2000 (as decided at sw meeting)
 *
 * Revision 1.53  2007-09-20 14:53:46  razeto
 * Renamed krypton tags to right names
 *
 * Revision 1.52  2007-09-19 11:07:01  razeto
 * Added anytag
 *
 * Revision 1.51  2007-09-13 06:01:02  razeto
 * const attribute for GetTimeDiff family
 *
 * Revision 1.50  2007-09-06 13:01:51  razeto
 * Added const to distance getter
 *
 * Revision 1.49  2007-08-27 14:58:38  dfranco
 * changed tag names
 *
 * Revision 1.48  2007-08-26 09:35:23  dfranco
 * Added new tags
 *
 * Revision 1.47  2007-08-18 22:39:52  razeto
 * Fixed a bug, added debug tag
 *
 * Revision 1.46  2007-08-08 21:04:17  razeto
 * GetTimeDiff behaves good even if there are no clusters
 *
 * Revision 1.45  2007-07-20 15:03:35  ddangelo
 * added a const in fadc getter
 *
 * Revision 1.44  2007-06-22 16:38:22  razeto
 * Tags moved
 *
 * Revision 1.43  2007-06-22 15:15:26  razeto
 * Moved to cycle 7
 *
 * Revision 1.42  2007-06-09 18:13:46  razeto
 * Added few tags (do not break event format)
 *
 * Revision 1.41  2007-06-06 10:51:16  razeto
 * Modified the bx_phys tags and added BxPosition::GetDistance (without changing the event format)
 *
 * Revision 1.40  2007-06-05 15:07:39  razeto
 * Modified the bx_phys tags (without changing the event format)
 *
 * Revision 1.39  2007-06-01 15:57:00  ddangelo
 * added the enabled crates word to root event.
 * One plain getter in echidna event
 *
 * Revision 1.38  2007-05-30 16:02:48  ddangelo
 * muon has_cluster variable converted to int. -1 for mcr disabled condition.
 * runtime enabled condition of crates recovered from builder logic
 *
 * Revision 1.36  2007-05-25 14:57:54  ddangelo
 * cosmetics
 *
 * Revision 1.35  2007-05-25 14:37:34  ddangelo
 * added management of btb flags.
 * Redesigned internal low level access.
 *
 * Revision 1.34  2007-05-07 16:41:42  ddangelo
 * n_live_channels renamed as n_live_pmts as requested by ale.
 * getters and normalize improved
 *
 * Revision 1.33  2007-05-07 15:06:38  ddangelo
 * BxPhysTags added in BxLabenRecCluster
 *
 * Revision 1.32  2007-05-07 14:40:47  ddangelo
 * cleaned up.
 * BxPhysTags reintroduced and upgraded
 *
 * Revision 1.31  2007-05-07 13:40:26  ddangelo
 * applying patch to flag TObjects with cycle numbers
 *
 * Revision 1.30  2007-05-07 12:49:33  ddangelo
 * cleaned up useless methods
 *
 * Revision 1.29  2007-05-04 16:33:24  pallas
 * Changed variables for alpha beta
 * Now tailtot is an array of 10 values
 *
 * Revision 1.28  2007-05-03 17:25:39  ddangelo
 * added n_live_channels for run-time information. added relative getter.
 * added functions Normalize() for common data types.
 * Both internal and root event
 *
 * Revision 1.27  2007-04-27 14:23:29  pallas
 * Small changes to alpha beta variables
 *
 * Revision 1.26  2007-03-28 17:40:24  ddangelo
 * introduced dbn "energy" class.
 * introduced muon decoded charge/npmts.
 * completed the job on Laben Rec Cluster.
 *
 * Revision 1.25  2007-03-27 15:18:20  ddangelo
 * variables npe_conc, charge_conc, nhits_conc added to laben cluster
 * f4_pe renamed as f4_charge in bx_laben_event.hh
 * decoded_charge and decoded_npmts added tu muon event
 *
 * Revision 1.24  2007-03-23 19:49:05  ddangelo
 * debugging
 *
 * Revision 1.23  2007-03-22 19:36:57  ddangelo
 * implemented muon db channel
 *
 * Revision 1.22  2007-03-22 16:09:35  ddangelo
 * added muon clustered level
 *
 * Revision 1.21  2007-03-15 19:17:19  ddangelo
 * pid event removed.
 * laben event upgraded with classes: bx_laben_shaped_cluster and bx_laben_ab_cluster
 * bx_laben_rec_cluster is now a parent class for the 2 new ones.
 * BxEvent modified accordingly: BxLabenRecHit and BxLabenRecCluster added.
 * BxPidEvent removed.
 *
 * Revision 1.20  2007-03-15 10:43:58  razeto
 * Calculate time difference between clusters
 *
 * Revision 1.19  2007-03-08 16:10:58  razeto
 * Add include to allow compiling with newest root (v5.15)
 *
 * Revision 1.18  2007-02-22 16:10:43  pallas
 * 	Avoiding warnings
 *
 * Revision 1.17  2006-12-14 15:01:02  ddangelo
 * added BTB trhreshold in root file
 *
 * Revision 1.16  2006/11/06 17:46:56  razeto
 * Reorganized few names
 *
 * Revision 1.15  2006/10/23 15:34:34  ddangelo
 * applied ale's patch to inlcude laben integrity flags
 *
 * Revision 1.14  2006/10/13 15:31:18  razeto
 * Added variables from raw data
 *
 * Revision 1.13  2006-08-21 16:52:53  ddangelo
 * added decoded index to clustered hit and plain getter.
 * Computation from map in BxLaben::operator=.
 * More fancy getters to come.
 *
 * Revision 1.12  2006/08/21 16:17:40  ddangelo
 * fixed a typo
 *
 * Revision 1.11  2006/08/21 11:12:56  razeto
 * Updated to use new bx_barn
 *
 * Revision 1.10  2006/05/16 10:43:39  ddangelo
 * added description of new fadc variables
 *
 * Revision 1.9  2006/05/15 09:50:03  ddangelo
 * added 5 variables for improved fadc daq
 *
 * Revision 1.8  2006/05/15 09:23:06  ddangelo
 * added class BxPhysTags to be filled by physics tools
 * Removed IgnoreObjectStreamer() calls from all constructors.
 *
 * Revision 1.7  2006/05/08 17:31:33  razeto
 * Added db_channel patch for fadc (sent to the mailing list)
 *
 * Revision 1.6  2006/05/05 12:28:16  razeto
 * Added getter for lg in clutstered hit (auth from Davide)
 *
 * Revision 1.5  2006/01/09 16:18:49  razeto
 * Moved db_barn getters to .cc, to avoid db_barn.hh inclusion in BxEven.hh
 *
 * Revision 1.4  2005/12/30 14:08:22  razeto
 * Added direct object getter (other than vector of objects) for clusters, peaks and frames (auth from Davide)
 *
 * Revision 1.3  2005/12/30 12:17:12  razeto
 * Updated the comment to GetTimeDifference
 *
 * Revision 1.2  2005/12/18 14:35:54  razeto
 * Implememted GetTimeDifferece method as asked in the echidna meeting
 *
 * Revision 1.1  2005/12/12 19:37:57  razeto
 * Moved BxEvent from root to event
 *
 * Revision 1.47  2005/12/12 19:08:09  razeto
 * Split position and energy reco results (auth from maintainers)
 *
 * Revision 1.46  2005/12/06 16:15:04  ddangelo
 * added GPS time
 *
 * Revision 1.45  2005/12/02 18:11:39  misiaszek
 *
 * GetDbChannel( db_barn* db ) function added to decoded/clustered laben hits
 * lg and run number added to decoded/clustered laben hits (auth from Davide D'Angelo)
 *
 * Revision 1.44  2005/11/18 16:42:55  razeto
 * Added new integer charge variable to decoded hit and cluster (auth from davide)
 *
 * Revision 1.43  2005/10/13 13:48:52  razeto
 * Added radius getter (auth from davide)
 *
 * Revision 1.42  2005/10/04 20:20:36  razeto
 * Added iterations value during minimization
 *
 * Revision 1.41  2005/09/26 15:56:20  razeto
 * Added rough time from clustering (auth from davide)
 *
 * Revision 1.40  2005/09/26 12:37:57  razeto
 * Fixed a bug and introduced likelihood in BxPosition
 *
 * Revision 1.39  2005/09/20 17:17:12  razeto
 * Fixed pid shape event (auth from mantainer)
 *
 * Revision 1.38  2005/09/19 13:03:45  ddangelo
 * added peaks from splitting.
 * A tmp size is also added to account for ROOT mishandling of vector<float>.
 *
 * Revision 1.37  2005/08/22 11:26:47  ddangelo
 * added Pid classes
 *
 * Revision 1.36  2005/07/27 16:58:57  ddangelo
 * added flag for "broad" (i.e. spaparanzeted) clusters
 *
 * Revision 1.35  2005/07/11 17:16:59  ddangelo
 * removed global event
 * added pid event (partially)
 * untested
 *
 * Revision 1.34  2005/07/06 16:36:42  ddangelo
 * 'Bool_t is_out_of_gate' re-introduced in class BxLabenDecodedHit
 *
 * Revision 1.33  2005/07/06 16:14:18  ddangelo
 * in BxLabenClusteredHit: 'Double_t time' moved to 'Float_t time'
 * in BxLabenCluster: 'Float_t start_time' moved to 'Double_t start_time'
 *
 * Revision 1.32  2005/06/20 16:46:46  ddangelo
 * added bx_position_reco_dbn
 *
 * Revision 1.31  2005/04/21 16:02:28  ddangelo
 * - removed out_of_gate variable from laben decode hit since no longer used.
 * - fixed a small bug in filling the n_decoded_hits and the the n_clusters variables in laben event.
 *
 * Revision 1.30  2005/03/14 19:34:54  ddangelo
 * added indexes for backtracing laben hits
 *
 * Revision 1.29  2005/03/14 13:45:22  ddangelo
 * added "order" in laben raw, decoded and clustered stages
 *
 * Revision 1.28  2005/03/01 15:19:10  razeto
 * Merged with cycle_2
 *
 * Revision 1.27  2005/02/10 20:05:04  ddangelo
 * some comments added
 *
 * Revision 1.26  2005/02/10 16:59:36  ddangelo
 * added a variable with size for every vector in every class and relative getters.
 * added some comments
 * removed compilation warnings
 *
 * Revision 1.25  2004/12/22 17:02:48  ddangelo
 * minor improvements
 *
 * Revision 1.24  2004/12/15 18:18:18  ddangelo
 * added position reco results. First draft, to be checked.
 *
 * Revision 1.23.2.1  2005/02/04 10:06:16  ddangelo
 * added elec event time to mctruth
 *
 * Revision 1.23  2004/12/03 11:59:47  ddangelo
 * added dsum waveform
 * added raw channel waveform (very haevy, disabled by default)
 * Both temporarily use fixed length C-style arrays. To be improved sometime.
 *
 * Revision 1.22  2004/12/02 16:44:09  ddangelo
 * added digital sum ampl+charge+peak
 * waveform still missing
 *
 * Revision 1.21  2004/12/01 15:13:41  ddangelo
 * added classes for fadc event.
 * added a few vairiables.
 * Work in progress, more stuff to come.
 *
 * Revision 1.20  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.19  2004/11/26 14:06:23  razeto
 * Added Mantainer field
 *
 * Revision 1.18  2004/11/24 21:08:21  ddangelo
 * added lists of hits outside clusters/fragments
 * Use of std::vector introduced instead of TClonesArray
 * Required ROOT v4.00/??
 *
 * Revision 1.17  2004/10/06 13:49:35  ddangelo
 * fixed a memory leak in BxLabenCluster
 *
 * Revision 1.16  2004/09/23 16:20:39  ddangelo
 * added TClonesArray with hits in cluster class.
 * Splitting doesn't work yet, but file is written correctly.
 *
 * Revision 1.15  2004/09/22 13:25:39  ddangelo
 * updated bx_reco_event to bx_echidna_event
 *
 * Revision 1.14  2004/09/22 11:00:55  ddangelo
 * added support for bx_baricentrator (by Ale)
 *
 * Revision 1.13  2004/08/31 13:28:41  ddangelo
 * added charge to decoded hit
 *
 * Revision 1.12  2004/07/26 17:36:33  ddangelo
 * fBits and fUniqueID suppressed in ROOT file.
 *
 * Revision 1.11  2004/07/26 16:08:35  ddangelo
 * ClassDef macros moved at the end of class definitions.
 *
 * Revision 1.10  2004/07/25 16:41:43  ddangelo
 * added logical channel to laben decoded hits.
 *
 * Revision 1.9  2004/07/25 16:29:31  ddangelo
 * some development, not much active at the moment
 *
 * Revision 1.8  2004/07/13 14:50:48  ddangelo
 * added BxClusteredHit. Currently used due to the TClonesArray problem.
 *
 * Revision 1.7  2004/07/13 13:37:10  ddangelo
 * added McTruth and McTruthFrame to root event.
 * McTruthHit commented out for the moment, due to ROOT problems.
 * To be debugged.
 *
 * Revision 1.6  2004/07/07 15:45:26  ddangelo
 * added BxLabenCluster.
 * Some minor debugging.
 *
 * Revision 1.5  2004/06/25 14:46:39  ddangelo
 * fixed a bug.
 * removed a comment.
 *
 * Revision 1.4  2004/06/22 13:02:52  ddangelo
 * added laben laser and trigger time in BxLaben obj
 *
 * Revision 1.3  2004/06/07 17:14:03  ddangelo
 * conditional compile macros introduced.
 * Laben raw and decoded hit introduced.
 * Muon raw and decoded hit implemented with 2 different classes.
 *
 * Revision 1.2  2004/06/03 15:00:46  ddangelo
 * (Copy) constructor replaced by operator=() for every class.
 * Still in a tmp state.
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
