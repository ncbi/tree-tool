// evolution.cpp

#undef NDEBUG
#include "../common.inc"

#include "evolution.hpp"



namespace DistTree_sp
{


Real intersection2dissim (size_t size1,
                          size_t size2,
                          size_t intersection,
                          size_t intersection_min,
	                        Prob sizes_ratio_min,
	                        bool ave_arithP)
{
  ASSERT (intersection <= size1);
  ASSERT (intersection <= size2);
  ASSERT (isProb (sizes_ratio_min));
  
  
  if (verbose ())
    cout << '\t' << intersection 
         << '\t' << size1 
         << '\t' << size2 
         << endl;
         
  if (  (Real) min (size1, size2) 
  	  / (Real) max (size1, size2)
  	    < sizes_ratio_min
  	 )
  	return NAN;
         
  if (intersection < intersection_min)
  	return NAN;
  	
  const Real dissim1 = log ((Real) size1 / (Real) intersection);
  const Real dissim2 = log ((Real) size2 / (Real) intersection);
  ASSERT (dissim1 >= 0);
  ASSERT (dissim2 >= 0);
  
  if (ave_arithP)
    return ave_arith (dissim1, dissim2);
  //return ave_geom  (dissim1, dissim2);
  else
    return ave_harm  (dissim1, dissim2);
}




// Hashes

Hashes::Hashes (const string &fName)
{
  reserve (10000);  // PAR

  LineInput hf (fName);
  size_t prev = 0;
  while (hf. nextLine ())
  {
    const size_t hash = str2<size_t> (hf. line);
    ASSERT (hash);
    if (hash <= prev)
      throw runtime_error ("Hash " + hf. line + " is not greater than the previous hash " + toString (prev));
    *this << hash;
    prev = hash;
  } 
  searchSorted = true;
/*
  if (empty ())
    throw runtime_error ("No hashes for " + fName);
*/
}



// Read Hashes from a binary file ??




// DissimAverage::DissimAttr

DissimAverage::DissimAttr::DissimAttr (const PositiveAttr1* attr_arg,
		                                   Real center_arg)
: Named (attr_arg->name)
, attr (attr_arg)
, center (center_arg)
{ 
	ASSERT (center > 0);
	ASSERT (attr);
}



DissimAverage::DissimAttr::DissimAttr (const string &line,
		                                   bool loadStat)
{
	// Cf. saveText()
	istringstream iss (line);
	iss >> name; 
	if (loadStat)
	{
		iss >> center >> var;
		if (! (center > 0))
			throw runtime_error ("Center should be > 0");
    if (! (center < INF))
			throw runtime_error ("Center should be finite");
		if (! (var >= 0))
			throw runtime_error ("Variance should be >= 0");
	}
}



void DissimAverage::DissimAttr::qc () const
{
	if (! qc_on)
		return;

	Named::qc ();
	
	ASSERT (goodName (name));
	IMPLY (! isNan (center), center > 0);
	ASSERT (var >= 0);
}



void DissimAverage::DissimAttr::setOutlier (Real value_target) const
{ 
	outlier =    /*value > 1   // PAR
	          &&*/ /*abs*/ (value - value_target) / getSD () > 1;   // PAR
}



void DissimAverage::DissimAttr::setValue (size_t objNum)
{
	ASSERT (attr);
	value = (*attr) [objNum] / center;
}



void DissimAverage::DissimAttr::setVar (const PositiveAttr1& averageAttr)
{
	ASSERT (attr);
	ASSERT (& averageAttr. ds == & attr->ds);
		
	var = 0;
	Real mult_sum = 0;
	FFOR (size_t, i, attr->ds. objs. size ())
	{
		const Real x = (*attr)     [i];
		const Real y = averageAttr [i];
		if (! goodValue (x))
			continue;
		if (! goodValue (y))
			continue;
		const Real weight = 1 /*1 / value2var (y)*/;
		ASSERT (weight >= 0);
		if (weight != INF)
		{
			var      += weight * sqr (x - y);  
			mult_sum += weight;  
		}
	}
	if (mult_sum == 0)
		throw runtime_error ("Attribute " + name + " has no data");
	ASSERT (mult_sum > 0);
	var /= mult_sum;
	ASSERT (var >= 0);
}




// DissimAverage

DissimAverage::DissimAverage (const string &fName,
                              bool loadStat)
{ 
	LineInput f (fName);
	while (f. nextLine ())
	  if (! f. line. empty ())
		{ 
			const DissimAttr attr (f. line, loadStat);
			attr. qc ();
		  attrs << attr;
		}
}    


void DissimAverage::qc () const
{
	if (! qc_on)
		return;
	
	const Dataset* ds = nullptr;
	for (const DissimAttr& dissimAttr : attrs)
	{
		dissimAttr. qc ();
	  if (ds)
	  	{ ASSERT (ds == & dissimAttr. attr->ds); }
	  else
	  	ds = & dissimAttr. attr->ds;
  }
}



Real DissimAverage::get () const
{
	ASSERT (! attrs. empty ());

  for (const DissimAttr& attr : attrs)
   	attr. outlier = false;

	Real ave = NAN;
	for (;;)
	{
		const Real ave_prev = ave;
		
    Real sum = 0;
    Real weights = 0;
    for (const DissimAttr& attr : attrs)
    	if (   DissimAttr::goodValue (attr. value)
    		  && ! attr. outlier
    		  && ! attr. bad ()
    		 )
    	{
    		ASSERT (attr. value >= 0);
    		const Real weight = attr. getWeight ();
    		ASSERT (weight >= 0);
    		if (weight != INF)
    		{
	    		sum     += weight * attr. value;
	    		weights += weight;
	    	}
    	}  		
    ave = sum / weights;
  //cout << ave << endl;
    
    if (isNan (ave) || abs (ave - ave_prev) < 1e-6)  // PAR
    	break;
    
    for (const DissimAttr& attr : attrs)
    {
    	attr. setOutlier (ave);
  	//cout << attr. ave << '\t' << attr. getSD () << '\t' << attr. outlier << endl;
    }
  }
  
  return ave;
}



void DissimAverage::calibrate (PositiveAttr1& averageAttr)
{
  for (DissimAttr& attr : attrs)
  {
  	attr. var = 0.1;  
  	ASSERT (! attr. bad ());
  	attr. mv. clear ();
  }

  Progress prog;
  for (;;)
  {
    // averageAttr
  	FFOR (size_t, i, averageAttr. ds. objs. size ())
  	{
  		setValues (i);
  		averageAttr [i] = get ();
		  for (DissimAttr& attr : attrs)
		  	if (   ! attr. bad ()
		  		  && DissimAttr::goodValue (attr. value)
		  		 )
		  		attr. mv. add (attr. outlier);
  	}  	

    const MVector vars_old (getVars ());

    // DissimAttr::var
    const Real sd = setVars (averageAttr);
 		
 		const Real diff = getVars (). maxAbsDifferenceVec (vars_old);
    prog (toString (diff) + " " + toString (sd));
 		if (diff < 1e-6)  // PAR
 			break;
 	}
}



MVector DissimAverage::getVars () const
{ 
	MVector vars (attrs. size ());
	FFOR (size_t, i, attrs. size ())
	  vars [i] = attrs [i]. var;
	return vars;
}




}


