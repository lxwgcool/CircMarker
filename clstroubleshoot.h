#ifndef CLSTROUBLESHOOT_H
#define CLSTROUBLESHOOT_H
#include <string>
using namespace std;

class ClsTroubleShoot
{
public:
    ClsTroubleShoot();
    ~ClsTroubleShoot();

public:
    //Find out valid reads by mapping
    void CheckBWAMappingResult(string strBamFilePath);
};

#endif // CLSTROUBLESHOOT_H