// req_feature2dissim.cpp

/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* Author: Vyacheslav Brover
*
* File Description:
*   Print dissimilarities computed from features for requested pairs of objects
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



namespace 
{
  
  
struct Feature : Named
{
  bool optional {false};
  
  explicit Feature (string &&line)
    : Named (move (line))
    { trim (name);
      if (trimSuffix (name, " 0"))
        ;
      else if (trimSuffix (name, " 1"))
        optional = true;
      trim (name);
    }
    
  bool operator< (const Feature &other) const
    { return name < other. name; }
  bool operator== (const Feature &other) const
    { return name == other. name; }
};



Vector<Feature> readFeatures (const string &fName)
{
  Vector<Feature> vec;  vec. reserve (16);  // PAR

  LineInput f (fName);
  while (f. nextLine ())
    vec << Feature (move (f. line));
    
  vec. sort ();
  QC_ASSERT (vec. isUniq ());

  return vec;
}



double optional_weight = 0.0;



double vecs2dissim_half (const Vector<Feature> &vec1,
                         const Vector<Feature> &vec2)
{
  double n = 0.0;
  for (const Feature& f : vec2)
  {
  #if 0
    if (f. optional)
      continue;
    if (! vec1. containsFast (f))
      n++;
  #else
    const size_t index = vec1. binSearch (f);
    if (f. optional)
    {
      if (index != NO_INDEX && ! vec1 [index]. optional)
        n += optional_weight;
    }
    else
      if (index == NO_INDEX)
        n += 1.0;
      else
        if (vec1 [index]. optional)
          n += optional_weight;
  #endif
  }
  return n;
}



double vecs2dissim (const Vector<Feature> &vec1,
                    const Vector<Feature> &vec2)
{
  return   vecs2dissim_half (vec1, vec2)
         + vecs2dissim_half (vec2, vec1);
}



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print dissimilarities computed from features for requested pairs of objects")
    {
      version = VERSION;
  	  addPositional ("req", "File with pairs of objects");
  	  addPositional ("dir", "Directory with <object> files containing feature-optional(0/1) pairs");
  	  addPositional ("optional_weight", "Weight of optional-nonoptional match (0..1)");
  	}



	void body () const
	{
		const string reqFName        =       getArg ("req");
		const string dir             =       getArg ("dir");
		             optional_weight = stod (getArg ("optional_weight"));
		QC_ASSERT (optional_weight >= 0.0);
		QC_ASSERT (optional_weight <= 1.0);
		
		
    LineInput reqF (reqFName);
    string obj1, obj2, rest;
    Istringstream iss;
    Vector<Feature> vec1;
    Vector<Feature> vec2;
    while (reqF. nextLine ())
    {
      trim (reqF. line);
      iss. reset (reqF. line);
      rest. clear ();
      iss >> obj1 >> obj2 >> rest;
      QC_ASSERT (! obj1. empty ());
      QC_ASSERT (! obj2. empty ());
      QC_ASSERT (rest. empty ());

      vec1 = readFeatures (dir + "/" + obj1);
      vec2 = readFeatures (dir + "/" + obj2);
      
      cout << obj1 << '\t' << obj2 << '\t' << vecs2dissim (vec1, vec2) << endl;
    }
	}
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



