/****************************************************************
 SwiftSample.H
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SwiftSample_H
#define INCL_SwiftSample_H
#include <iostream>
#include "BOOM/Comparator.H"
using namespace std;
using namespace BOOM;

class SwiftSample {
  float p, q, theta;
  void computeTheta();
public:
  SwiftSample();
  SwiftSample(float p,float q);
  float getP() const;
  float getQ() const;
  float getTheta() const;
};

class SwiftSampleComparator : public Comparator<SwiftSample>
{
public:
  bool equal(SwiftSample &a,SwiftSample &b)
  { return a.getTheta()==b.getTheta(); }
  bool greater(SwiftSample &a,SwiftSample &b)
  { return a.getTheta()>b.getTheta(); }
  bool less(SwiftSample &a,SwiftSample &b)
  { return a.getTheta()<b.getTheta(); }
};

#endif

