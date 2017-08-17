// attr2_2phylip.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/matrix.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print an Attr2 in Phylip format")
    {
  	  addPositional ("file", dmSuff + "-file");
  	  addPositional ("attr2Name", "Attribute name of an object-object table in the " + dmSuff + "-file");
  	  addPositional ("map", "Output map file for nw_rename");
  	}



	void body () const
	{
		const string fName     = getArg ("file");
		const string attr2Name = getArg ("attr2Name");
		const string mapFName  = getArg ("map");		
		
		
    Dataset ds (fName);
    
    // dist
    const RealAttr2* dist = nullptr;
    {
      const Attr* attr = ds. name2attr (attr2Name);
      ASSERT (attr);
      dist = attr->asRealAttr2 ();
    }
    ASSERT (dist);

    Matrix& matr = const_cast <RealAttr2*> (dist) -> matr;
    
    OFStream f ("", mapFName, "");    

    cout << ds. objs. size () << endl;
    ONumber on (cout, 5, false);
    FOR (size_t, row, ds. objs. size ())
    {
      const string name (("X" + toString (row + 1) + "          "). substr (0, 10));
      cout << name;
      f << name << ' ' << ds. objs [row] -> name << endl;
      FOR (size_t, col, ds. objs. size ())
        cout << ' ' << matr. get (false, row, col);
      cout << endl;
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



