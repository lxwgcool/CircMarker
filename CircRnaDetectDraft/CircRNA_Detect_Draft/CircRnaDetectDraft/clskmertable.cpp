#include "clskmertable.h"
#include <iostream>
#include "../../ShareLibrary/clskmeralgorithm.h"
#include <algorithm>

ClsKmerTable::ClsKmerTable()
{
}

ClsKmerTable::~ClsKmerTable()
{
    m_mpKT.clear();
}

void ClsKmerTable::CreateKmerTable(string strDNARef, int iReadsLen, float fKmerRatio,
                                   vector<St_Row_Chrom>& vChrom)
{
    //Get paramter
    int iBoundLen = iReadsLen * fKmerRatio;

    //Read Fasta:
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    vector<St_Fasta> vFasta;
    pFastaReader->ReadFastaRegular(strDNARef, vFasta);
    delete pFastaReader;
    pFastaReader = NULL;

    //Create Kmer Table:
    //Get tag value for each exon
    St_PosInfo stPosInfo;
    m_mpKT.clear();

    stPosInfo.ucChromIndex = 0;
    //I just remember: (1)the chrmosone should be consecutive,
    //                 (2) only the chromosone which could be found back in reference will be used for Kmer Table
    for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin(); itrChrom != vChrom.end(); itrChrom++)
    {
        St_Fasta* pCurRefFa = NULL;

        for(vector<St_Fasta>::iterator itrRef = vFasta.begin(); itrRef != vFasta.end(); itrRef++)
        {
            //cout << itrRef->strName << endl;
            string strRefName = "";
            if(itrRef->strName.find(' ') == -1)
                strRefName = itrRef->strName;
            else
                strRefName = itrRef->strName.substr(0, itrRef->strName.find(' '));
            //cout << "strRefName: " << strRefName << endl;

            if(strRefName == itrChrom->strName)
            {
                pCurRefFa = &(*itrRef);
                break;
            }            
        }
        if(pCurRefFa == NULL)
        {
            cout << "Do not find related reference chromosone: " << itrChrom->strName << endl;
            stPosInfo.ucChromIndex++;
            continue;
        }

        stPosInfo.uiGeneIndex = 0;
        for(vector<St_Raw_Gene>::iterator itrGene = itrChrom->vRG.begin();
            itrGene != itrChrom->vRG.end(); itrGene++)
        {

            //Get corresponding reference sequence
            stPosInfo.ucTranscriptIndex = 0;
            for(vector<St_Raw_Transcript>::iterator itrRT = itrGene->vRT.begin();
                itrRT != itrGene->vRT.end(); itrRT++)
            {
                stPosInfo.ucExonIndex = 0;
                for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin();
                    itrExon != itrRT->vRExon.end(); itrExon++)
                {
                    if(itrExon->iStart > itrExon->iEnd)
                    {
                        cout << "Exon Start > End, Strange!" << endl;
                    }
                    else
                    {
                        CollectKmerInfo(pCurRefFa, &(*itrExon), stPosInfo, iBoundLen);
                    }
                    stPosInfo.ucExonIndex++;
                }
                stPosInfo.ucTranscriptIndex++;
            }
            stPosInfo.uiGeneIndex++;
        }
        stPosInfo.ucChromIndex++;
    }
}

//我们在建立Kmer表的时候并没有特别去考虑相应的RC或者regular的direction的区别问题，也就是一样的对待
void ClsKmerTable::CollectKmerInfo(St_Fasta* pCurRefFa, St_Raw_Exon* pExon,
                                   St_PosInfo& stPosInfo, int iBoundLen)
{

    //看看这里的reverse complememntory 我怎么去处理！！！ --> 这里我们不需要处理，我们直接根据 正向序列去进行建表
    //1: check 这个 exon值不值得建表 -->
    if( !pExon->GetIsSupportCircRNA() ) // 如果不support 我们就不建表
        return;
    //2: 这里我们应用我们的新逻辑
    //(1) we do not use every kmer to build the kmer table
    //(2) 我们仅仅考虑前后40个bp的kmer：前后reads长度的40%
    //(3) 并且将阈值定到能够完全大于一遍的kmer的数量:目前的策略 len(reads) - (2*len(kmer) +  1).
    //    就我们手上的case而言，是42 bps，刚好大于一点点单边的最大kmer数目
    // 我们要保证每一边有“iBoundLen”个kmer ！！！

    if(pExon->GetLength() <= 2 * iBoundLen + KMERLEN*2 - 1)
    {
        //stPosInfo.cPart = 'U';
        int iSplitPoint = (pExon->GetLength() - KMERLEN)/2;
        for(int i = pExon->iStart - 1; i<= pExon->iEnd - KMERLEN; i++)  // Notice: 在这里我们的start 和 end 都属于 真实的position，也就是从0开始的
        {            
            if(i < (pExon->iStart - 1) + iSplitPoint)
                stPosInfo.cPart = 'S';
            else
                stPosInfo.cPart = 'E';
            string strCurKmer = pCurRefFa->strSeq.substr(i, KMERLEN);
            InsertCurKmer(strCurKmer, stPosInfo);
        }
    }
    else // Exon的长度足够去进行左右两边的采样 --> 那我们开始吧
    {
        //对于当前exon的start -->
        stPosInfo.cPart = 'S';
        for(int i = pExon->iStart - 1; i<= pExon->iStart + iBoundLen - 1; i++)
        {
            string strCurKmer = pCurRefFa->strSeq.substr(i, KMERLEN);
            InsertCurKmer(strCurKmer, stPosInfo);
        }
        //对于当前eoxn的end <--
        stPosInfo.cPart = 'E';
        int iStart = pExon->iEnd - KMERLEN - iBoundLen;
        for(int i = iStart; i<= iStart + iBoundLen; i++)
        {
            string strCurKmer = pCurRefFa->strSeq.substr(i, KMERLEN);
            InsertCurKmer(strCurKmer, stPosInfo);
        }
    }
}

void ClsKmerTable::InsertCurKmer(string strCurKmer, St_PosInfo& stPosInfo) //只是不把exon中的还有N的kmer去建立表
{
    if(strCurKmer.find('n') != string::npos || strCurKmer.find('N') != string::npos )
        return;

    unsigned int uiKmer = ConvertKmerToNum32(strCurKmer);

    if(m_mpKT.find(uiKmer) == m_mpKT.end()) // 找不到
    {
        m_mpKT[uiKmer].push_back(stPosInfo);
    }
    else
    {
        vector<St_PosInfo>& vPosInfo = m_mpKT[uiKmer];
        if( find(vPosInfo.begin(), vPosInfo.end(), stPosInfo) == vPosInfo.end() ) // 找不到
        {
            vPosInfo.push_back(stPosInfo);
        }
    }
}

map<unsigned int, vector<St_PosInfo> >& ClsKmerTable::GetKT()
{
    return m_mpKT;
}
