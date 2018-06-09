// $Id: bx_sample.hh,v 1.21 2015/01/09 15:03:08 misiaszek Exp $
// class bx_sample
// description of an event for a generic particle id algorythm
// a bx_sample is basically a vector of numbers
// it inherits from <vector> and it is not much more than that
// MP 22-02-07 re-design from old code
//

#ifndef __BX_SAMPLE_HH_
#define __BX_SAMPLE_HH_
#ifndef CYCLE_NUMBER
#define CYCLE_NUMBER 18
#endif

#include <iostream>
#include <string>
#include <vector>
#include "TH1.h"

class bx_sample: public std::vector<float> {
public:
  
  // empty ctor
  bx_sample();
  
  // dtor
  virtual ~bx_sample();

  // set and get sample name
  const std::string &sample_name() const {return name;}
  void set_name(const std::string& n) {name=n;}

  // compute truncated mean; returns also moda 
  float truncated_mean(float& moda);
  
  //ClassDef(bx_sample,CYCLE_NUMBER)
  
private:
  std::string name;
  friend std::ostream& operator<<(std::ostream& os,const bx_sample& v);
  friend std::istream& operator>>(std::istream& is, bx_sample& v);
};

// stream output operator
std::ostream& operator << (std::ostream& os,const bx_sample& v);

// stream input operator
std::istream& operator >> (std::istream& is, bx_sample& v);


#endif
