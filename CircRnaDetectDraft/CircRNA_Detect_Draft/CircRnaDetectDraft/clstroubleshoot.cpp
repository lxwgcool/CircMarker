#include "clstroubleshoot.h"
#include "../../ShareLibrary/clsbasealgorithm.h"

ClsTroubleShoot::ClsTroubleShoot()
{
}

ClsTroubleShoot::~ClsTroubleShoot()
{

}

void ClsTroubleShoot::CheckBWAMappingResult(string strBamFilePath)
{    
    ofstream ofsStart;
    ofsStart.open("./Start_Support_Reads.txt");
    ofstream ofsEnd;
    ofsEnd.open("./End_Support_Reads.txt");
    int offset = 5 * 2;

    //Go to check sorted bam file
    //1: parse bam file & output the expected reads
    BamReader* pBamReader = new BamReader();
    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;
    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        //Case 1: do not map  ------->
        if(!al.IsMapped() || al.QueryBases == "")
            continue;

        //Case 2: If the reads maps fine
        ///1: Check if there is soft clip
        int iStart = al.Position;
        bool bFindSoftClip = false;
        for(std::vector<CigarOp>::iterator itr = al.CigarData.begin();
            itr != al.CigarData.end(); itr++)
        {
            switch(itr->Type)
            {
                case 'M': // alignment match (can be a sequence match or mismatch)
                    iStart += itr->Length;
                    break;
                case 'I': // insertion to the reference
                    iStart += itr->Length;
                    break;
                case 'D': // deletion from the reference
                case 'N':  // skipped region from the reference
                    break;
                case 'S':  // soft clipping (clipped sequences present in SEQ)
                case 'H':  // hard clipping (clipped sequences NOT present in SEQ)
                    bFindSoftClip = true;
                    break;
                case 'P': // padding (silent deletion from padded reference)
                case '=': // sequence match
                case 'X': // sequence mismatch
                    iStart += itr->Length;
                    break;
            }
            if(bFindSoftClip)
                break;
        }

        if(bFindSoftClip)
        {
            //Check start 9991949:
            if(abs(iStart - 1191425) < offset)
            {
                ofsStart << "Query Base: " << al.QueryBases << endl;
                ofsStart << "Align Base: " << al.AlignedBases << endl;
                ofsStart << (al.IsReverseStrand() ? "-" : "+") << endl << endl;                
                ///ofsStart << IntToStr(al.AlignedBases.length()) << " " << IntToStr(iStart - al.Position) << endl;
                ///ofsStart << IntToStr(al.Position) << " " << IntToStr(al.Position + al.AlignedBases.length())
                ///         << endl << ">>>>>>>>>>>" << endl;
            }
            //Check end 9992956:
            if(abs(iStart - 1203372) < offset)
            {
                ofsEnd << "Query Base: " << al.QueryBases << endl;
                ofsEnd << "Align Base: " << al.AlignedBases << endl;
                ofsEnd << (al.IsReverseStrand() ? "-" : "+") << endl << endl;                
                ///ofsEnd << IntToStr(al.AlignedBases.length()) << " " << IntToStr(iStart - al.Position) << endl;
                ///ofsEnd << IntToStr(al.Position) << " " << IntToStr(al.Position + al.AlignedBases.length())
                ///       << endl << ">>>>>>>>>>>" << enrdl;
            }
        }
    }
    delete pBamReader;
    pBamReader = NULL;

    ofsStart.close();
    ofsEnd.close();
}

