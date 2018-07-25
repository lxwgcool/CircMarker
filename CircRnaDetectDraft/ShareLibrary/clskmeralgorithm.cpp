#include "clskmeralgorithm.h"
#include <iostream>
#include <string.h>

unsigned int ConvertKmerToNum32(string strKmer) // The maximun length is 16
{
    unsigned int iValue = 0;
    for(unsigned int i = 0; i < strKmer.length(); i++)
    {
        if(strKmer[i] == 'A' ||  strKmer[i] == 'a') // 0
        {
            iValue |= 0;
        }
        else if(strKmer[i] == 'T' ||  strKmer[i] == 't') // 1
        {
            iValue |= 1;
        }
        else if(strKmer[i] == 'G' ||  strKmer[i] == 'g') // 2
        {
            iValue |= 2;
        }
        else if(strKmer[i] == 'C' ||  strKmer[i] == 'c') // 3
        {
            iValue |= 3;
        }
        else
        {
            cout << "Abnormal Kmer String" << endl;
            cout << strKmer << endl;
            return -1;
        }

        if(i < strKmer.length() - 1) // not the end
            iValue <<= 2;
        else
            break;
    }
    return iValue;
}

//uint64_t
uint64_t ConvertKmerToNum64(string strKmer)
{
    uint64_t iValue = 0;
    for(unsigned int i = 0; i < strKmer.length(); i++)
    {
        if(strKmer[i] == 'A' ||  strKmer[i] == 'a') // 0
        {
            iValue |= 0;
        }
        else if(strKmer[i] == 'T' ||  strKmer[i] == 't') // 1
        {
            iValue |= 1;
        }
        else if(strKmer[i] == 'G' ||  strKmer[i] == 'g') // 2
        {
            iValue |= 2;
        }
        else if(strKmer[i] == 'C' ||  strKmer[i] == 'c') // 3
        {
            iValue |= 3;
        }
        else
        {
            cout << "Abnormal Kmer String" << endl;
            return -1;
        }

        if(i < strKmer.length() - 1) // not the end
            iValue <<= 2;
        else
            break;
    }
    return iValue;
}

string ConvertNum64ToKmer(uint64_t i64Value, int iKmerLen)
{
    uint64_t iTemp = 0;
    string strKmer;
    for(int i=0; i<iKmerLen; i++)
    {
        iTemp = i64Value;
        iTemp >>= 2;
        iTemp <<= 2;

        if( (iTemp | 0) == i64Value ) // this is 'A'
            strKmer = "A" + strKmer;
        else if( (iTemp | 1) == i64Value ) // this is 'T'
            strKmer = "T" + strKmer;
        else if( (iTemp | 2) == i64Value ) // this is 'G'
            strKmer = "G" + strKmer;
        else if( (iTemp | 3) == i64Value ) // this is 'C'
            strKmer = "C" + strKmer;

        i64Value >>= 2;
    }
    return strKmer;
}

string ConvertNum32ToKmer(unsigned int i32Value, int iKmerLen)
{
    unsigned int iTemp = 0;
    string strKmer;
    for(int i=0; i<iKmerLen; i++)
    {
        iTemp = i32Value;
        iTemp >>= 2;
        iTemp <<= 2;

        if( (iTemp | 0) == i32Value ) // this is 'A'
            strKmer = "A" + strKmer;
        else if( (iTemp | 1) == i32Value ) // this is 'T'
            strKmer = "T" + strKmer;
        else if( (iTemp | 2) == i32Value ) // this is 'G'
            strKmer = "G" + strKmer;
        else if( (iTemp | 3) == i32Value ) // this is 'C'
            strKmer = "C" + strKmer;

        i32Value >>= 2;
    }
    return strKmer;
}
