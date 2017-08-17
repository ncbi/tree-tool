// uniKernel.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "dataset.hpp"
using namespace DM_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
  : Application ("Kernel p.d.f. estimation. Print distribution with the step = window")
  { 
	  addPositional ("file", dmSuff + "-file without the extension");
	  addPositional ("attrName", "name of a real attribute");
	  addKey ("window", "Window of UniKernel; 0 <=> default", "0");  
  }


	void body () const
	{
		const string fName    = getArg ("file");  
		const string attrName = getArg ("attrName");
		const Real window     = str2real (getArg ("window"));
		
		
    Dataset ds (fName);
    ds. qc ();
    
    const Sample sm (ds);
    
    const NumAttr1* attr = nullptr;
    {
      const Attr* attr_ = ds. name2attr (attrName);
      ASSERT (attr_);
      attr = attr_->asNumAttr1 ();
    }
    ASSERT (attr);
    
    const UniVariate<NumAttr1> analysis (sm, *attr);
    
    UniKernel uniKernel (analysis);
    if (window)
      uniKernel. setParam (window);
    else
      uniKernel. estimate ();
    uniKernel. qc ();
    cout << uniKernel. nameParam () << endl;    
    cout << "Mode: " << (*attr) [uniKernel. getModeObjNum ()] << endl;

    if (uniKernel. halfWindow)
    {
      const Real step = 2 * uniKernel. halfWindow;  // PAR
      Real x = uniKernel. attr_min - step;
      while (leReal (x, uniKernel. attr_max + step))
      {
        cout << x << "\t" << uniKernel. pdf (x) << endl;
        x += step;
      }
    }
	}
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
//return ThisApplication().AppMain(argc, argv); 
  ThisApplication app;
  return app. run (argc, argv);
}



