#include "bx_position_reco_mach4.hh"
#include "bx_echidna_event.hh"
#include "bx_mach4_minuit_wrapper.hh"
#include "db_channel.hh"
#include "db_calib.hh"
#include "bx_dbi.hh"
#include "constants.hh"
#include "TH1F.h"
#include "barn_interface.hh"

#include <cmath>

#define SPEED_OF_LIGHT 0.299792458 /* m/ns */
#define IDX_REF_PC 1.53
#define RAW_TIME_HIGH_CUT 54 /* ns */

static double n_c = IDX_REF_PC / SPEED_OF_LIGHT; ;

const static bx_position_reco_mach4 *pm;


//Wrapper functions to put member functions into Minuit
static void minuit_likelihood_fcn(int32_t &npar, double * grad, double & fcnval,
								  double *xval, int32_t iflag)
{
	//                         x        y       z        t
	fcnval = pm->Likelihood(xval[0], xval[1], xval[2], xval[3]);
}

//define index n as a function of npe to compensate the energy dependence of
//position reconstruction --Jingke July 26, 2009

double Get_Index(double npe){
	double n_idx = IDX_REF_PC;
	if(npe <= 2000) n_idx = 1.88 - 0.18*exp(-npe/2466.3);
	else n_idx = 1.96 - 0.28444*exp(-npe/3476.06);
	return n_idx;
}

////////

//ctor
bx_position_reco_mach4::bx_position_reco_mach4() : bx_base_module("bx_position_reco_mach4", bx_base_module::main_loop), mn(0)
{
	require_event_stage(bx_detector::laben, bx_base_event::clustered);
}


void bx_position_reco_mach4::begin()
{

	//Load in parameters from config file
	refidx = get_parameter("refidx").get_float ();

	//Cache local copies of PMT positions and timing data usability
	
	for(size_t channel = 1; channel <= constants::laben::channels; channel++)
	{
		const db_channel_laben& channel_info = dynamic_cast<const db_channel_laben&>(bx_dbi::get()->get_channel(channel));

		xpmt[channel] = channel_info.pmt_x();
		ypmt[channel] = channel_info.pmt_y();
		zpmt[channel] = channel_info.pmt_z();
	}
	
	// Set up scintillator PDF
	
	//One photon
	scint_pdf[0].have_reflection = false;
	scint_pdf[0].SetQ0 (0.5509);
	scint_pdf[0].SetSigma (3.6679);
	scint_pdf[0].SetTau (0.7415, 13.8752);
	scint_pdf[0].SetReflectionAmplitude (0.005628092);
	scint_pdf[0].SetReflectionCenter (45.2413);
	scint_pdf[0].SetReflectionSigma (25.4981);
	scint_pdf[0].CalculateCoefficients();
	
	//Two photons
	scint_pdf[1].have_reflection = false;
	scint_pdf[1].SetQ0 (0.2066);
	scint_pdf[1].SetSigma (3.635);
	scint_pdf[1].SetTau (6.8285, 0.6126);
	scint_pdf[1].SetReflectionAmplitude (0.005628092);
	scint_pdf[1].SetReflectionCenter (45.2413);
	scint_pdf[1].SetReflectionSigma (25.4981);
	scint_pdf[1].CalculateCoefficients();
	
	//Three photons
	scint_pdf[2].have_reflection = false;
	scint_pdf[2].SetQ0 (0.9006);
	scint_pdf[2].SetSigma (3.631);
	scint_pdf[2].SetTau (0.4979, 4.3609);
	scint_pdf[2].SetReflectionAmplitude (0.005628092);
	scint_pdf[2].SetReflectionCenter (45.2413);
	scint_pdf[2].SetReflectionSigma (25.4981);
	scint_pdf[2].CalculateCoefficients();
	
	
	//Set up MINUIT
	pm = this;
	if(!mn)
	{
		try 
		{
			mn = new MinuitWrapper(minuit_likelihood_fcn);
		}
		catch (...) 
		{
			mn = NULL;
		}
	}
}






bx_echidna_event* bx_position_reco_mach4::doit(bx_echidna_event *ev)
{
	
	size_t n_clusters = ev->get_laben().get_nclusters();
	
	//Loop over all clusters in the event
	
	for(size_t i = 0; i < n_clusters; i++)
		
	{
		
		//First hit in this cluster
		cluster = &ev->get_laben().get_cluster(i);
		end_hit = cluster->get_clustered_nhits() - 1;
		first_hit = 0;
		double t_1st = cluster->get_clustered_hit(first_hit).get_time();
		
		
		// We want to avoid the situation where a PMT has triggered several ns
		// earlier than all the other ones ("pre-pulse"), which messes up the
		// position recon.  Ignore all hits more than 10 ns before the peak time
		// (this should also exclude light leaks from fiber optic cables in
		// diffuser ball calibration runs).
		// XXX: make the time before peak (10 ns) configurable
		
		while( first_hit + 1 < end_hit &&
			  cluster->get_clustered_hit(first_hit + 1).get_time() > t_1st + 3)
		{
			first_hit ++;
			t_1st = cluster->get_clustered_hit(first_hit).get_time();
		}
		
		if(first_hit + 1 >= end_hit)
		{
			//got to end of cluster without finding any close-enough hits
			//something wrong -- get out of here
			//	  std::cout << "Got to end of cluster without finding close enough hits" << std::endl;
			return ev;
		}
		//Get reference to mach4 position -- we will fill it later
		bx_position_mach4& pos = ev->get_laben().get_cluster(i).get_position_mach4();
		bx_position_mach4& pos_fixed = ev->get_laben().get_cluster(i).get_position_mach4_fixed();		
		
		// Tweak Minuit a bit
		mn->SetParameterName(0, "x");
		mn->SetParameterName(1, "y");
		mn->SetParameterName(2, "z");
		mn->SetParameterName(3, "t");
		mn->SetParLimits(3, t_1st - 100, t_1st + 50);
		// Allow reconstructed times after ^^^^^^^^ the first hit, as it's
		// possible (albeit statistically unlikely) that the first few hits in
		// the cluster were dark noise.
		
		
		double pos_info[8] = {0.0}; //x,y,z,t,dx,dy,dz,dt
		bool is_valid = false;
		
		// With the PMT photocathodes ~ 6.5 m from the vessel center, and a
		// vessel of radius 4.25 m, the first hit should occur between 11.3 and
		// 32.6 ns after the event.  It is most important that events near the
		// Fiducial Volume boundary (r = 3 m, distance to PMTs = 3.5 m) are
		// correctly reconstructed, hence assume a time offset of 17.6 ns.
		pos_info[3] = t_1st - 17.6;
		
		if(cluster->get_charge() >= 0) /* In Mach4 set this to ~130 if generating shapes file */
		{
			// You can add more calls to Reconstruct() below, but note it's the *last*
			// one that "counts" for the purposes of calculating dt and setting the
			// emission times for the alpha/beta module.  Reconstruct takes the values
			// of x,y,z,t by reference, so they can be re-used to start from a better
			// initial condition.  (If Reconstruct returns false, it doesn't change
			// them.)  
			
			scint_pdf[0].SetPDFType(ScintPDF::STANDARD);
			scint_pdf[1].SetPDFType(ScintPDF::STANDARD);
			scint_pdf[2].SetPDFType(ScintPDF::STANDARD);
			
			//npe-dependent index of refraction position reconstruction.
			
			double n_idx = Get_Index( ev->get_laben().normalize_charge(cluster->get_charge()) );
			pos.f4_user = n_idx;
			
			if( Reconstruct( n_idx, ev, i, pos_info) ){
				
				pos.b_converged = true;
				is_valid = true;
			}
			
			else
				pos.b_converged = false;
			
			pos.f4_x = pos_info[0];
			pos.f4_y = pos_info[1];
			pos.f4_z = pos_info[2];
			pos.f4_t = pos_info[3];
			
			pos.f4_dx = pos_info[4];
			pos.f4_dy = pos_info[5];
			pos.f4_dz = pos_info[6];
			pos.f4_dt = pos_info[7];
			
			//stupid hack with effective index of refraction to make it look
			//like the collaboration actually understands light propagation
			//
			  
			//initialize variables
			for(int32_t k =0; k<8; k++) pos_info[k] = 0; 
			pos_info[3] = t_1st - 17.6;
			
			pos_fixed.f4_user = refidx;
			
			if( Reconstruct(refidx, ev, i, pos_info) )
				pos_fixed.b_converged = true;
			
			else
				pos_fixed.b_converged = false;
			
			pos_fixed.f4_x = pos_info[0];
			pos_fixed.f4_y = pos_info[1];
			pos_fixed.f4_z = pos_info[2];
			pos_fixed.f4_t = pos_info[3];
			
			pos_fixed.f4_dx = pos_info[4];
			pos_fixed.f4_dy = pos_info[5];
			pos_fixed.f4_dz = pos_info[6];
			pos_fixed.f4_dt = pos_info[7];
			
		}
		
		
		if( !is_valid )
		{
			// The reconstructed value of t may be completely meaningless, but we
			// still need a value, however rough, to calculate dt.  So estimate
			// via the time of the 1st hit minus an offset (let's say 20 ns).
			pos.f4_t = t_1st - 20;
		}
		
		
		
		
		/* XXX  TO BE FIGURED OUT XXX */
		/* -------------------------- 
		 // Determine the time of the physical event relative to start of run.
		 cluster->rel_time = ev->Trigger()->GetRelTime() // run start to trigger
		 + (-ladata->GetPulseReference()    // trigger to start of trigger window
		 + t) * 1e-9;                    // start of window to physical event
		 
		 -------------------------- */
		 
	}
	ev->get_laben().mark_stage (bx_base_event::reconstructed_mach4);
	
	return ev;
}


void bx_position_reco_mach4::end()
{
	//  delete xpmt;
	//  delete ypmt;
	//  delete zpmt;
	//  delete pm;
}



bool bx_position_reco_mach4::Reconstruct(double n_idx, 
										 bx_echidna_event *ev, size_t cl,
										 double *pos_info)
{
	
	int32_t errflag, convergence, is_valid = 0;
	//size_t n = end_hit - first_hit;
	//  double r = std::sqrt(x*x + y*y + z*z);
	
	//set refractive index to desired value
	n_c = n_idx / SPEED_OF_LIGHT;
	
	
	//set maximum |x,y,z| value to permit for posn recon
	const double limit = 8;
	for(size_t i = 0; i < 3; i++)
	{
		mn->SetParLimits(i, -limit, limit);
	}
	
	const size_t n_iterations = 1;
	for(iteration = 0; iteration < n_iterations; iteration++)
	{
		mn->SetParameter(0, pos_info[0]);
		mn->SetParameter(1, pos_info[1]);
		mn->SetParameter(2, pos_info[2]);
		mn->SetParameter(3, pos_info[3]);
		mn->Execute(errflag, convergence);
		
		//x,y,z,t
		pos_info[0] = mn->GetParameter(0), pos_info[1] = mn->GetParameter(1);
		pos_info[2] = mn->GetParameter(2), pos_info[3] = mn->GetParameter(3);
		
		double r = std::sqrt(pos_info[0]*pos_info[0] + pos_info[1]*pos_info[1] + pos_info[2]*pos_info[2]);
		is_valid = (errflag == 0 && convergence == 3 && r < limit) ? 1 : 0;
		
		//dx, dy, dz, dt
		pos_info[4] = mn->GetParError(0); pos_info[5] = mn->GetParError(1);
		pos_info[6] = mn->GetParError(2); pos_info[7] = mn->GetParError(3);
	}
	return is_valid;
}


double time_of_flight(double x , double y , double z ,
					  double x0, double y0, double z0)
{
	return n_c * sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0));
}

///////

double bx_position_reco_mach4::Likelihood(double x, double y, double z, double t) const
{
	double result = 0;
	//  double tpmt0 = cluster->get_time();
	for (size_t i = 0; i < (size_t)cluster->get_clustered_nhits(); i++)
	{
		size_t channel = cluster->get_clustered_hit(i).get_decoded_hit().get_raw_hit().get_logical_channel();
		
		double tpmt = cluster->get_clustered_hit(i).get_time();
		if(tpmt < RAW_TIME_HIGH_CUT)
		{
			
			double tof = time_of_flight(x, y, z,
										xpmt[channel], ypmt[channel], zpmt[channel]);
			
			int32_t photons = int32_t(cluster->get_clustered_hit(i).get_decoded_hit().get_charge_pe() + 0.5);
			//Correct PDF for multiple occupancies
			
			if(photons <= 1)
				result -= scint_pdf[0].Log(tpmt - t - tof);
			
			if(photons == 2)
				result -= scint_pdf[1].Log(tpmt - t - tof);
			
			if(photons >= 3)
				result -= scint_pdf[2].Log(tpmt - t - tof);
		}
	}
	return result;
}

//dtor
bx_position_reco_mach4 :: ~bx_position_reco_mach4() {}
