// convex.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "dataset.hpp"
#include "prediction.hpp"
using namespace DM_sp;



namespace
{



// --> dataset ??
struct Convex : MultiVariate<PositiveAttr1>
// Convex components
// Model:
//   attr_i = k_i * combo + epsilon_i
//   E (epsilon_i) = 0
//   cov (epsilon_i, combo) = 0
//   => cov (attr_i, combo) = k_i var (combo)
{
  typedef  MultiVariate<PositiveAttr1>  P;
  
  Dataset& ds;
  
  // Output
  PositiveAttr1* combo {nullptr};
  MVector beta;
    // >= 0; sum() = 1
  

  Convex (const Sample &sample_arg,
          const Space1<PositiveAttr1> &space_arg,
          const string &comboAttrName)
    : P (sample_arg, space_arg)
    , ds (* const_cast <Dataset*> (sample. ds))
    , combo (new PositiveAttr1 (comboAttrName, ds, decimals_def))
	  , beta (space. size (), 1 / (Real) space. size ())
    {
      Matrix attrSim (space. size ());
      {
      	Verbose vb;
			  setAttrSim (attrSim, true);
			}
		  if (attrSim. inverse (). isZero ())
		    throw runtime_error ("Cannot invert VC");

		  Space1<PositiveAttr1> comboSp (ds, false); 
		  comboSp << combo;

      Progress prog;
      for (;;)
      {
	    	const MVector beta_prev (beta);	
	
	    	setCombo ();
	
			  MVector xy (beta. size ());  
			  FFOR (size_t, i, space. size ())
			  {
				  L2LinearNumPrediction lr (sample, comboSp, * space [i]);
				  lr. solveUnconstrained ();
				  lr. qc ();  
				  if (isNan (lr. absCriterion))
				    throw runtime_error ("No data in " + space [i] -> name);
				  xy [i] = lr. beta [0];
			  }
			  
			  beta. multiply (false, attrSim, false, xy, false);
			  
			  FFOR (size_t, i, beta. size ())
			    maximize (beta [i], 0.0);
			  
			  const Real beta_sum = beta. sum ();
			  ASSERT (beta_sum > 0);  // prove ??
			  beta. putProdAll (1 / beta_sum);

		 		const Real diff = beta. maxAbsDifferenceVec (beta_prev);
		    prog (toString (diff));
		 		if (diff < 1e-6)  // PAR
		 			break;
			}
    }
  Convex* copy () const final
    { return new Convex (*this); }
  void qc () const override
    {
    	if (! qc_on)
    		return;
    	P::qc ();
    	ASSERT (combo);
    	ASSERT (& combo->ds == sample. ds);
    	ASSERT (beta. min () >= 0);
    	ASSERT (eqReal (beta. sum (), 1));
    }
  void saveText (ostream &os) const override
    { FFOR (size_t, j, space. size ())
    	  os << space [j] -> name << '\t' << beta [j] << endl;
    }
    

  void setCombo ()
    {
    	FFOR (size_t, i, sample. size ())
    	{
    		Real sum = 0;
    		Real mult_sum = 0;
    		FFOR (size_t, j, space. size ())
    		{
    			const Real x = (*space [j]) [i];
    		  if (isNan (x))
    		  	continue;
    		//ASSERT (beta [j] >= 0);
    		  sum      += beta [j] * x;
    		  mult_sum += beta [j];
    		}
    		(*combo) [i] = sum / mult_sum;
    	}
    }
};




struct ThisApplication : Application
{	
  ThisApplication ()
	  : Application ("Convex components. Infinities are replaced by NANs.")
		{
		  addPositional ("file", dmSuff + "-file without the extension with Positive attributes");	  
		  // Output
		  addKey ("output", "Output combined attribute");
		}
	
	
	
	void body () const final
	{
		const string inFName     = getArg ("file");
		const string outputFName = getArg ("output");
		

    Dataset ds (inFName);
    
    const Sample sm (ds);
    if (! positive (sm. mult_sum))
      throw runtime_error ("Too small data size");

    Space1<PositiveAttr1> sp (ds, true);
    {
	    size_t n = 0;
	    for (const PositiveAttr1* attr : sp)
	    	n += const_cast <PositiveAttr1*> (attr) -> inf2missing ();
	    const Real size = (Real) ds. objs. size () * (Real) ds. attrs. size ();
	    cerr << "# Infinities: " << n << " (" << ((Real) n * 100.0 / size) << "%)" << endl;    
	  }    
    ds. qc ();
    
    cerr << "Convex ..." << endl;
    const Convex conv (sm, sp, "combo");  // PAR
    conv. qc ();  
    conv. print (cout);
    
    if (! outputFName. empty ())
    {
    	OFStream f (outputFName);
    	FFOR (size_t, i, ds. objs. size ())
    	  f << ds. objs [i] -> name << '\t' << (* conv. combo) [i] << endl;
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


