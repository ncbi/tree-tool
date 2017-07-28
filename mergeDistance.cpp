// mergeDistance.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "matrix.hpp"
#include "dataset.hpp"
using namespace DM_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print a " + dmSuff + "-file with merged <Distance> attrributes of a list of " + dmSuff + "-files")
    {
      // Input
  	  addPositional ("file_list", "List of " + dmSuff + "-files without \"" + dmSuff + "\" with a binary <Distance> attribute");
  	  addPositional ("attr", "Distance attribute");
  	  // Output
  	  addPositional ("output", dmSuff + "-file without \"" + dmSuff);
  	}



	void body () const
	{
		const string file_list = getArg ("file_list");
		const string attrName  = getArg ("attr");
		const string output    = getArg ("output");
		ASSERT (! output. empty ());
		
		
		Set<string> objNames;
		size_t fileNum = 0;
    {
      // Pass 1
      LineInput li (file_list);
      while (li. nextLine ())
      {
        fileNum++;
        const Dataset ds (li. line);
        ds. qc ();
        ASSERT (ds. getUnitMult ());
        for (const Obj* obj : ds. objs)
          objNames << obj->name;
      }
    }
    ASSERT (! objNames. empty ());
    cout << "# Objects: " << objNames. size () << endl;
    cout << "# Files: " << fileNum << endl;


    Dataset ds_new;
    for (const string& objName : objNames)
      ds_new. appendObj (objName);
    ds_new. setName2objNum ();
    ASSERT (ds_new. objs. size () == objNames. size ());
    

		Vector<PositiveAttr2*> dists;  dists. reserve (fileNum);
		Vector<Real> means;            means. reserve (fileNum);
    {
      // Pass 2
      LineInput li (file_list);
      Progress prog;
      while (li. nextLine ())
      {
        prog (li. line);
        const Dataset ds (li. line);
        const Attr* attr = ds. name2attr (attrName);
        ASSERT (attr);
        const PositiveAttr2* dist = attr->asPositiveAttr2 ();
        ASSERT (dist);
        auto dist1 = new PositiveAttr2 (attrName + "_" + toString (li. lineNum), ds_new, dist->decimals);
        dists << dist1;
        MeanVar mv;        
        FOR (size_t, row, ds. objs. size ())
        {
          const size_t row_new = ds_new. getName2objNum (ds. objs [row] -> name);
          ASSERT (row_new != NO_INDEX);
          FOR (size_t, col, ds. objs. size ())
          {
            const size_t col_new = ds_new. getName2objNum (ds. objs [col] -> name);
            ASSERT (col_new != NO_INDEX);
            const Real val = dist->get (row, col);
            if (! isNan (val))
            {
              dist1->put (row_new, col_new, val);  
              if (row_new < col_new)
                mv << val;
            }
          }
        }
        const Real mean = mv. getMean (); 
        means << mean;
        cout << dist1->name << ": " << "Average = " << mean << endl;
        ASSERT (positive (mean));
        FOR (size_t, row, ds_new. objs. size ())
          FOR (size_t, col, ds_new. objs. size ())
            if (! dist1->isMissing2 (row, col))
              dist1->put (row, col, dist1->get (row, col) / mean);
      }
    }
    ASSERT (dists. size () == fileNum);
    ASSERT (means. size () == fileNum);    
    
    
    auto dist_new = new PositiveAttr2 (attrName, ds_new, 6);  // PAR
    Vector<Real> vars (fileNum, 1);    
    for (;;)
    {
      FOR (size_t, row, ds_new. objs. size ())
        FOR (size_t, col, ds_new. objs. size ())
        {
          Real val = 0;
          Real wSum = 0;
          FOR (size_t, i, dists. size ())
          {
            const PositiveAttr2* dist1 = dists [i];
            if (! dist1->isMissing2 (row, col))
            {
              const Real weight = 1 / vars [i];
              val  += weight * dist1->get (row, col);
              wSum += weight;
            }
          }
          if (positive (wSum))
            dist_new->put (row, col, val / wSum);
          else
            dist_new->setMissing (row, col);
        }
        
      Real varDiff = 0;
      FOR (size_t, i, dists. size ())
      {
        const PositiveAttr2* dist1 = dists [i];
        Real resid2 = 0;
        size_t n = 0;
        FOR (size_t, row, ds_new. objs. size ())
          FOR (size_t, col, row)
            if (! dist1->isMissing2 (row, col))
            {
              resid2 += sqr (dist1->get (row, col) - dist_new->get (row, col));
              n++;
            }
        ASSERT (n);        
        const Real var = resid2 / (Real) n;
        varDiff += abs (vars [i] - var);
        vars [i] = var;
      }
      cerr << "varDiff = " << varDiff << endl;  // ??
      if (varDiff / (Real) vars. size () < 1e-6)  // PAR
        break;
    }
    FOR (size_t, i, dists. size ())
      cout << dists [i] -> name << " variance = " << vars [i] << endl;
            
    
    Real mean = 0;
    Real wSum = 0;
    FOR (size_t, i, fileNum)
    {
      const Real weight = 1 / vars [i];
      mean += weight * means [i];
      wSum += weight;
    }
    mean /= wSum;  // geometric average ??
    FOR (size_t, row, ds_new. objs. size ())
      FOR (size_t, col, ds_new. objs. size ())
        if (! dist_new->isMissing2 (row, col))
          dist_new->put (row, col, dist_new->get (row, col) * mean);
      
      
    {
      OFStream f ("", output, dmExt);
      const Sample sm (ds_new);
      VectorPtr<Attr> attrs (1, dist_new);
      sm. save (attrs, f);
    }
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



