#ifndef CLSJELLYFISH_H
#define CLSJELLYFISH_H
#include <string>
using namespace std;

class ClsJellyfish
{
public:
    static ClsJellyfish& GetInstance()
    {
        static ClsJellyfish instance;
        return instance;
    }

private:
    ClsJellyfish(){}
    ClsJellyfish(const ClsJellyfish&);
    ClsJellyfish& operator = (const ClsJellyfish&);

public:
    string GetSolidKmer(string strReadsFaPath, int iKmerLen, bool bForceNew = false,
                        string strAnchorKmerFileName = "AnchorKmer.fa",
                        bool bUseFixThreshould = false, int iCoverageThreshold=10);

    string GetAllKmer(string strReadsFaPath, int iKmerLen);

};

#endif // CLSJELLYFISH_H
