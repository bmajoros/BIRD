/****************************************************************
 SwiftSample.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SwiftSample.H"
using namespace std;


SwiftSample::SwiftSample()
{
  // ctor
}



SwiftSample::SwiftSample(float p,float q)
  : p(p), q(q)
{
  computeTheta();
  // ctor
}



float SwiftSample::getP() const
{
  return p;
}



float SwiftSample::getQ() const
{
  return q;
}



void SwiftSample::computeTheta()
{
  theta=q/(1-q)/(p/(1-p));
}



float SwiftSample::getTheta() const
{
  return theta;
}




