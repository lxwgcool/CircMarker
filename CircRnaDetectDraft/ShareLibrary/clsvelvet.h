#ifndef CLSVELVET_H
#define CLSVELVET_H
#include <string>
using namespace std;

class ClsVelvet
{
public:
    static ClsVelvet& GetInstance()
    {
        static ClsVelvet instance;
        return instance;
    }
private:
    ClsVelvet(){}
    ClsVelvet(const ClsVelvet&);
    ClsVelvet& operator = (const ClsVelvet&);

public:
    string LocalAssembly(string strFilePath, int iMinOutputLen = 100, int iKmerLen = 31,
                         bool bAppendOrgReads = false);
private:
    string RelaxMergeByContigsMerger(string strFilePath);
};

#endif // CLSVELVET_H
