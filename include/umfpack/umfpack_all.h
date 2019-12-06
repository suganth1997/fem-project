#include <iostream>
extern "C"
{
  #include "umfpack.h"
}
void UMFPACK_return(int ret)
{
#ifdef DEBUG
#define _DEBUG_PRINT(x) std::cout << x;
    if (ret == 0) 
    {
        _DEBUG_PRINT("solved sucessfully (" << ret << ")" << std::endl);
    }
    else if (ret ==   1) 
    {
        _DEBUG_PRINT("solved failed (" << ret << "): Matrix singular" << std::endl);
    }
    else if (ret ==   2)
    {
        _DEBUG_PRINT("solved failed (" << ret << "): det(A) != 0 but < eps" << std::endl);
    }
    else if (ret ==   3)
    {
        _DEBUG_PRINT("solved failed (" << ret << "): det(A) != 0 but > inf" << std::endl);
    }
    else if (ret ==  -1)
    {
        _DEBUG_PRINT("solved failed (" << ret << "): Not enough memory!" << std::endl);
    }
    else if (ret ==  -3)
    {
        _DEBUG_PRINT("solved failed (" << ret << "): Used Numeric object is invalided!" << std::endl);
    }
    else if (ret ==  -4)
    {
        _DEBUG_PRINT("solved failed (" << ret << "): Used Symbolic object is invalided!" << std::endl);
    }
    else if (ret ==  -5)
    {
        _DEBUG_PRINT("solved failed (" << ret << "): Argument missing!" << std::endl);
    }
    else if (ret ==  -6)
    {
        _DEBUG_PRINT("solved failed (" << ret << "): Number of rows and columns must be greater 0!" << std::endl);
    }
    else if (ret ==  -8)
    {
        _DEBUG_PRINT("solved failed (" << ret << "): Invalidid matrix!" << std::endl);
    }
    else if (ret == -11)
    {
        _DEBUG_PRINT("solved failed (" << ret << "): Different pattern" << std::endl);
    }
    else if (ret == -13)
    {
        _DEBUG_PRINT("solved failed (" << ret << "): Invalidid system!" << std::endl);
    }
    else if (ret == -15)
    {
        _DEBUG_PRINT("solved failed (" << ret << "): Invalidid permutation!" << std::endl);
    }
    else if (ret == -17)
    {
        _DEBUG_PRINT("solved failed (" << ret << "): Error due to I/O!!!" << std::endl);
    }
    else if (ret == -911)
    {
        _DEBUG_PRINT("solved failed (" << ret << "): Internal error!!!" << std::endl);
    }
#endif
//   if (ret != 0) exit(4711);
}
