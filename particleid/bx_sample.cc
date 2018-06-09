/* BOREXINO Reconstruction program
 * 
 * implementation for class bx_sample
 * 02-08-02 M. Pallavicini - created
 * 20-01-04 D. Manuzio - modified
 * 20-06-05 MP into echidna
 * Maintainer: Marco Pallavicini <marco.pallavicini@ge.infn.it>
 *
 * $Id: bx_sample.cc,v 1.9 2007/07/06 09:18:13 pallas Exp $
 *
 *
 */

#include "bx_sample.hh"

//ClassImp(bx_sample)

// empty ctor
bx_sample::bx_sample() : name("UNKNOWN") {
}

// dtor
bx_sample::~bx_sample() {
}

//
// compute mean trucated around peak position
//
float bx_sample::truncated_mean(float& moda) {
	TH1F *h = new TH1F("temp_xxxxxxx","temp",50,0.,50.);
	int n = size();
	for ( int i=0; i<n; i++) 
		h->Fill( (*this)[i] );

	int NumHits=0;
	float MeanTime=0.;
	moda = h->GetBinCenter(h->GetMaximumBin());
	float low  = moda - 5.;
	if (low<0.) low=0.;
	float high = moda + 20.;

	for ( int i=0; i<n; i++ ) {
		float t = (*this)[i];
		if ( t>low && t<high ) {
			MeanTime += t;
			NumHits++;
		}
	}

	if ( NumHits ) 
		MeanTime /= (float)NumHits;
	else
		MeanTime = 0.;

	delete h;
	return MeanTime;
}

// operator <<
std::ostream& operator<< (std::ostream& os, const bx_sample& v) {
  os << "Name=" << v.sample_name() << std::endl;
  os << "Size=" << v.size() << std::endl;
  os << "Values=";
  for(unsigned int i=0; i<v.size(); i++)
    os << v[i] << ",";
  os << std::endl;
  return os;
}

// operator >>
std::istream& operator>> (std::istream& is, bx_sample& v) {
  std::string s;
  char s2[100] = "";

  // get name
  getline(is,s);
  sscanf(s.c_str(),"Name=%s",s2);
  v.set_name( std::string(s2) );
  
  // get size
  int size;
  getline(is,s);
  sscanf(s.c_str(),"Size=%d",&size);
  
  // clear sample
  v.clear();

  // get values
  getline(is,s);
  const char* buff = s.c_str()+7;
  for(int i=0; i<size; i++) {
    float ft;
    char s3[100] = "";
    int j=0;
    while(buff[j]!=',') {
      s3[j]=buff[j];
      j++;
    }
    buff += (j+1);
    sscanf(s3,"%f",&ft);
    v.push_back(ft);
  }
  return is;
}


