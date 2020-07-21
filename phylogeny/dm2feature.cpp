// dm2feature.cpp

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
*   Create files in makeFeatureTree format forma Data Master file
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace 
{
  
  
VectorPtr<Attr1> attrs;



void savePhen (size_t from, 
               size_t to, 
               Notype /*&res*/,
               const string &featureDirName,
               const Dataset &ds)
{
  Progress prog (to - from);  
  FOR_START (size_t, i, from, to)
  {
    prog ();
    OFStream f (featureDirName, ds. objs [i] -> name, ""); 
    for (const Attr1* attr : attrs)
      if (const BoolAttr1* boolAttr = attr->asBoolAttr1 ()) 
      {       
        const ebool value = boolAttr->getBool (i);
      #if 1
        if (value != EFALSE)
           f << attr->name << ' ' << (value == ETRUE ? '0' : '1') << endl;
      #else
        if (value != UBOOL)
          f << attr->name << ':' << (value == ETRUE ? '1' : '0') << endl;
      #endif
      }
      else if (const NominAttr1* nominAttr = attr->asNominAttr1 ())
      {
        if (! attr->isMissing (i))
          f << attr->name << ':' << nominAttr->value2str (i) << endl;
      }
      else
        ERROR;
  }
}




struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Create files in makeFeatureTree format forma Data Master file")
    {
  	  version = VERSION;
  	  addPositional ("data", dmSuff + "-file"); 
  	  addPositional ("attrs", "List of Boolean or nominal attributes in <data>");
  	  addPositional ("feature_dir", "Directory where the files will be created");
    }


	
	void body () const final
  {
	  const string dsFName     = getArg ("data");
	  const string attrsFName  = getArg ("attrs");
    const string featureDirName = getArg ("feature_dir");


    const Dataset ds (dsFName);
    ds. qc ();
            
    {
      LineInput f (attrsFName);
      string name;
      Istringstream iss;
      while (f. nextLine ())
      {
        iss. reset (f. line);
        name. clear ();
        iss >> name;
        QC_ASSERT (! name. empty ());
        const Attr* attr = ds. name2attr (name);
        QC_ASSERT (attr);
        if (const BoolAttr1* boolAttr = attr->asBoolAttr1 ())
          attrs << boolAttr;
        else if (const NominAttr1* nominAttr = attr->asNominAttr1 ())
          attrs << nominAttr;
        else
          throw runtime_error ("Attribute " + strQuote (name) + " is neither Boolean nor nominal");
      }
    }

    vector<Notype> notypes;
	  arrayThreads (savePhen, ds. objs. size (), notypes, cref (featureDirName), cref (ds));
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



