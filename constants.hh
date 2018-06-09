/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: constants.hh,v 1.22 2012/01/16 09:48:02 ddangelo Exp $
 *
 * Container class for all program constants,
 * to be implemented as public static const members for easy access.
 *
 */

#ifndef _CONSTANTS_HH
#define _CONSTANTS_HH

class constants {
  public:
    class number {
      public:
	static const float pi;
	static float rad_to_deg (float deg) { return deg * f4_rad_to_deg_factor; }
	static float deg_to_rad (float rad) { return rad / f4_rad_to_deg_factor; }
      private:
	static const float f4_rad_to_deg_factor;
    };
    class physic {
      public:
	static const float c;
	static const float c_m_ns;
    };
    class fiber {
      public:
	static const unsigned char number_of_bundles;
	static const unsigned char number_of_fibers_in_bundle[];
    };
    class trigger {
      public:
        static const unsigned char crate_number;
    };
    class laben {
      public:
	static const unsigned short nboards;
        static const unsigned char ncrates;
        static const unsigned short pmts;
        static const unsigned short channels;
	static const unsigned char board_per_rack;
	static const unsigned char channels_per_rack;
	static const unsigned char channels_per_board;
	static int board_in_rack (int lg) { return ((lg - 1) % channels_per_rack) / channels_per_board + 1; }
	static int channel_in_board (int lg) { return ((lg - 1) % channels_per_rack) % channels_per_board + 1; }
        static bool is_laben (int lg) { return lg > 0 && lg <= channels; }
	class frontend {
	  public:
	    static const unsigned char board_per_rack;
	    static int board_in_rack (int lg) {
	      int fe_couple = (laben::board_in_rack (lg) - 1) / 3 + 1;
	      if (laben::channel_in_board (lg) & 0x1) return fe_couple * 2 - 1;
	      else return fe_couple * 2;
	    }
	};
    };
    class muon {
      public:
        class tdc {
          public:
            static const unsigned char max_boards;
            static const unsigned char chips_per_board;
            static const unsigned char channels_per_chip;
            static const unsigned char channels_per_board;
            static const float ns_per_clock;
            static const float new_ns_per_clock;
        };
        static const unsigned char crate_number;
        static const unsigned short channels;
        static const unsigned short channel_offset;
        static const unsigned long max_hits;
        static const float pe_per_clock; //tmp, to be removed
        static bool is_muon (int lg) { return lg > channel_offset && lg <= (channels + channel_offset); }
    };
    static const unsigned short precalibration_nevents;
    static int crate (int lg) { return (lg - 1) / 160 + 1; }
};


#endif
/*
 * $Log: constants.hh,v $
 * Revision 1.22  2012/01/16 09:48:02  ddangelo
 * added gate width for new tdc
 *
 * Revision 1.21  2011-10-13 17:40:01  ddangelo
 * added patch to handle new OD TDCs
 *
 * Revision 1.20  2011-02-18 17:10:04  ddangelo
 * major code cleanup: removed fadc throughout the program
 *
 * Revision 1.19  2007-12-07 17:33:32  ddangelo
 * removed fixed muon tdc clock range
 *
 * Revision 1.18  2007-10-30 17:37:12  razeto
 * Added nboards for laben namespace
 *
 * Revision 1.17  2007-05-04 09:09:12  ddangelo
 * added laben::pmts
 *
 * Revision 1.16  2007-03-28 17:41:07  ddangelo
 * added a shortcut
 *
 * Revision 1.15  2006-12-15 13:42:49  pelczar
 * *** empty log message ***
 *
 * Revision 1.14  2006-09-08 10:30:43  ddangelo
 * added getters to check channel nature
 *
 * Revision 1.13  2005/09/23 14:01:19  zavatare
 * Fixed some formulae in laben namespace
 *
 * Revision 1.12  2005/09/20 16:01:09  razeto
 * Added some constants
 *
 * Revision 1.11  2005/03/18 16:35:57  razeto
 * Added some values
 *
 * Revision 1.10  2004/11/29 16:38:24  ddangelo
 * some fadc constants added (by Sergei)
 *
 * Revision 1.7  2004/10/01 10:24:06  razeto
 * Added some laser constants
 *
 * Revision 1.6  2004/06/01 11:35:12  ddangelo
 * debugging
 *
 * Revision 1.5  2004/06/01 11:21:55  ddangelo
 * reorganized with nested classes
 *
 * Revision 1.4  2004/04/27 13:35:23  ddangelo
 * trigger and fadc sub-obj added. crate numbers added.
 *
 * Revision 1.3  2004/04/12 16:11:38  razeto
 * Added a laben constant
 *
 * Revision 1.2  2004/04/02 11:24:26  ddangelo
 * some constants added
 *
 * Revision 1.1  2004/04/01 10:50:33  ddangelo
 * added static class for wide spread constants.
 *
 *
 */
