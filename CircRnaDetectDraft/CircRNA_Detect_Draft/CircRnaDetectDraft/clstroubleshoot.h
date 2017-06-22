#ifndef CLSTROUBLESHOOT_H
#define CLSTROUBLESHOOT_H
#include <string>
#include "../../ShareLibrary/clsbasealgorithm.h"
using namespace std;


class ClsTroubleShoot
{
public:
    ClsTroubleShoot();
    ~ClsTroubleShoot();

public:
    //Find out valid reads by mapping
#ifdef USEBWA
    void CheckBWAMappingResult(string strBamFilePath);
#endif
};

#endif // CLSTROUBLESHOOT_H
