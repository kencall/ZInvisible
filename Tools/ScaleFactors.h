#ifndef SCALEFACTORS_H
#define SCALEFACTORS_H

class ScaleFactors
{
public:
    static double sf_norm0b() 
    {
      return 0.783;//0.788
    }
    static double sfunc_norm0b() 
    {
	return 0.152;
    }
    static double sf_norm0bElec() 
    {
	return 1.059;
    }
    static double sfunc_norm0bElec() 
    {
	return 0.191;
    }

};

#endif
