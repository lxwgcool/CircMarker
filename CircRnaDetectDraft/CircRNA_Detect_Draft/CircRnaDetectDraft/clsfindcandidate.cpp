#include "clsfindcandidate.h"
#include "../../ShareLibrary/clskmeralgorithm.h"
#include "../../ShareLibrary/clsbasealgorithm.h"
#include "clsresultcomparison.h"
#include "unistd.h"
#include <algorithm>

const int SAMPLINGNUM = 8;
const float ARRYPOSRATIO[SAMPLINGNUM] = {.1, .2, .3, .4, .5, .6, .7, .8};
const int CSAMPLINGHITMIN = 2;

//#define DEBUG

int g_EE = 0;
int g_ES = 0;
int g_SS = 0;
int g_SE = 0;
int g_else = 0;
int g_large1 = 0;

ClsFindCandidate::ClsFindCandidate():m_iTotalNum(0), m_iRegNum(0)
{
    /*m_mpRegCircCandi[stCurCandi1] = 1;

    bool bFind = false;
    for(map<St_Candidate, int>::iterator itr = m_mpRegCircCandi.begin(); itr != m_mpRegCircCandi.end(); itr++)
    {
        if(itr->first.iStartPos == 21221994 && itr->first.iEndPos == 21220010)
        {
            bFind = true;
            break;
        }
    }
    if(bFind)
        cout << "1" << endl;
    else
        cout << "Failed: " << endl;

    St_Candidate stCurCandi2(1, 173744981, 173735285);
    m_mpRegCircCandi[stCurCandi2] = 1;
    bFind = false;
    for(map<St_Candidate, int>::iterator itr = m_mpRegCircCandi.begin(); itr != m_mpRegCircCandi.end(); itr++)
    {
        if(itr->first.iStartPos == 21221994 && itr->first.iEndPos == 21220010)
        {
            bFind = true;
            break;
        }
    }
    if(bFind)
        cout << "2" << endl;
    else
        cout << "Failed: " << endl;


    St_Candidate stCurCandi3(1, 21268823, 21220010);
    m_mpRegCircCandi[stCurCandi3] = 1;
    bFind = false;
    for(map<St_Candidate, int>::iterator itr = m_mpRegCircCandi.begin(); itr != m_mpRegCircCandi.end(); itr++)
    {
        if(itr->first.iStartPos == 21221994 && itr->first.iEndPos == 21220010)
        {
            bFind = true;
            break;
        }
    }
    if(bFind)
        cout << "3" << endl;
    else
        cout << "Failed: " << endl;
    */
}

// <PosInfo Pair, support count>
/* 我们的策略如下：
 * 1: 滤掉所有看有N的reads
 * 2: Evenly 取 10 个 点 进行匹配
 *  (1) 只有至少两个匹配上了我们才进行全部的匹配
 *  (2) 看匹配上的大多数kmer的方向
 *     a) 整体方向如果是reverse complememtory 那么我们将read进行翻转
 *     b) 如果整体方向是正的方向，那么我们保持原始reads序列不改变
 *  (3) 我们考虑几种情况
 *     a)  ------
 *         ===---=====
 *     b)  ------
 *         ===---====---====
 *     c)  -----------
 *            ====
 *
 * Notice:
 * (1) 我们不需要将reads强行转成大写，因为去进行kmer table 查找的时候总是需要转换成int32去查找的
 */
void ClsFindCandidate::CheckHitting(int iMinSupportReads, float fKmerRatio, int iReadsLen,
                                    map<unsigned int, vector<St_PosInfo> >& mpKT,
                                    vector<St_Row_Chrom>& vChrom,
                                    vector<St_Fastq>& vFastq, string strChromName)
{           
    //cout << endl << "------------------CheckHitting------------------" << endl << endl;

    vector<St_Candidate> vSelfCircCandi;
    vector<St_Candidate> vRegCircCandi;

    //2: Iterator each reads
    unsigned int arrySamplingKmer[SAMPLINGNUM];

    for(vector<St_Fastq>::iterator itr = vFastq.begin(); itr != vFastq.end(); itr++)
    {
        //(1) delete all of reads contain letter 'N'
        if(itr->strSeq.find('N') != string::npos) // 这个表示找到了
            continue;

        if(itr->strSeq.length() < iReadsLen - 10)
            continue;

        //如果没有N
        //(2) evenly sampling reads
        string strCurSeq = ""; // If this is valid reads could be used for further detection

        // a) normal direction
        //   i) get kmer ID
        bool bSeqRC = false;
        for(int i = 0; i < SAMPLINGNUM; i++)
        {
            string strTmp = itr->strSeq.substr(itr->strSeq.length()*ARRYPOSRATIO[i], KMERLEN);
            arrySamplingKmer[i] = ConvertKmerToNum32(strTmp);
        }
        //   ii) Check how may could be found back in vChrom
        if(CheckSampling(arrySamplingKmer, mpKT))
        {
            strCurSeq = itr->strSeq;
        }
        else // b) RC direction
        {
            for(int i = 0; i < SAMPLINGNUM; i++)
            {
                string strTmp = GetReverseCompelement(itr->strSeq.substr(itr->strSeq.length()*ARRYPOSRATIO[i], KMERLEN));
                arrySamplingKmer[i] = ConvertKmerToNum32(strTmp);
            }
            if(CheckSampling(arrySamplingKmer, mpKT))
            {
                strCurSeq = GetReverseCompelement(itr->strSeq);
                bSeqRC = true;
            }
        }

        //(3)Go Further Detection
        if(strCurSeq != "")
        {
            // Go Further Deteciton ->Go
            //cout << strCurSeq << endl;
            CheckHitForCurReads(strCurSeq, bSeqRC, mpKT, vChrom, fKmerRatio,
                                vSelfCircCandi, vRegCircCandi);
        }
    }

    /*
    ofstream ofs;
    ofs.open("./Vector.txt");
    ofs << IntToStr(m_vRegCandi.size()) << endl;
    for(vector<St_Candidate>::iterator itr = m_vRegCandi.begin(); itr != m_vRegCandi.end(); itr++)
    {
        ofs << IntToStr(itr->iStartPos) << ", " << IntToStr(itr->iEndPos) << endl;
    }
    ofs.close();*/

    /*
    cout << "we back to outside" << endl;
    for(map<St_Candidate, int>::iterator itr = m_mpRegCircCandi.begin(); itr != m_mpRegCircCandi.end(); itr++)
    {
        if(itr->first.iStartPos == 21221994 && itr->first.iEndPos == 21220010)
        {
            cout << "Got 21221994 & 21220010" << endl;
            cout << IntToStr(itr->second) << endl;
            break;
        }
    }*/

//    cout << "g_EE     --> " << IntToStr(g_EE) << endl;
//    cout << "g_ES     --> " << IntToStr(g_ES) << endl;
//    cout << "g_SS     --> " << IntToStr(g_SS) << endl;
//    cout << "g_SE     --> " << IntToStr(g_SE) << endl;
//    cout << "g_else   --> " << IntToStr(g_else) << endl;
//    cout << "g_large1 --> " << IntToStr(g_large1) << endl;

    RefineCandiate(iMinSupportReads, vSelfCircCandi, vRegCircCandi, strChromName);
}

void ClsFindCandidate::AssembleReads(string strReads1Path, string strReads2Path, vector<St_Fastq>& vFastq)
{
    //1: Reads1 是必须有值的
    ClsFastqReader* pFqReader = new ClsFastqReader();
    pFqReader->ReadFastqFile(strReads1Path, vFastq, false); // 在这里我们不强制转换成大写！
    pFqReader->ReadFastqFile(strReads2Path, vFastq, false); // 在这里我们不强制转换成大写！

//    if(vFastq.empty())
//    {
//        delete pFqReader;
//        pFqReader = NULL;
//        return;
//    }

//    //2: 如果Reads2 也有值：我们读出来之后去跟reads1取一部分拼接到一起，然后将组成的新的结果塞给vFastq中
//    if(strReads2Path != "")
//    {
//        int iReadsLen = vFastq[0].strSeq.length();

//        vector<St_Fastq> vAdditionalFastq;
//        pFqReader->ReadFastqFile(strReads2Path, vAdditionalFastq, true);
//        if(vAdditionalFastq.size() == vFastq.size()) // 确定他们都是一对一的关系
//        {
//            vector<St_Fastq> vTempFastq;
//            vTempFastq.resize(vAdditionalFastq.size());
//            for(unsigned int i = 0; i < vAdditionalFastq.size(); i++)
//            {
//                vTempFastq.at(i).strName = IntToStr(i);
//                vTempFastq.at(i).strComments = "+";
//                vTempFastq.at(i).strQuality = "...";
//                //第一个的头和第二个的头合并到一起
//                vTempFastq.at(i).strSeq = vFastq.at(i).strSeq.substr(0, iReadsLen/2) +
//                                          GetReverseCompelement(vAdditionalFastq.at(i).strSeq).substr(0, iReadsLen/2);
//            }
//            //Merge into vFastq
//            vFastq.insert(vFastq.end(), vTempFastq.begin(), vTempFastq.end());
//            vTempFastq.clear();
//        }
//        vAdditionalFastq.clear();
//    }

    delete pFqReader;
    pFqReader = NULL;
}

//这里有一个需要注意：reads反向过来进行比对，依然符合我们最简单的判断的归类方式
bool ClsFindCandidate::CheckSampling( unsigned int* arrySamplingKmer,
                                      map<unsigned int, vector<St_PosInfo> >& mpKT)
{
    int iHitNumber = 0;
    for(int i = 0; i< SAMPLINGNUM; i++)
    {
        if(mpKT.find(arrySamplingKmer[i]) == mpKT.end()) // do not find
            continue;
        iHitNumber++;
    }

    if(iHitNumber < CSAMPLINGHITMIN)
        return false;
    else
        return true;
}

void ClsFindCandidate::CheckHitForCurReads(string& strSeq, bool bSeqRC,
                                           map<unsigned int, vector<St_PosInfo> >& mpKT,
                                           vector<St_Row_Chrom>& vChrom, float fKmerRatio,
                                           vector<St_Candidate>& vSelfCircCandi,
                                           vector<St_Candidate>& vRegCircCandi)
{
    //不能用只考虑最高两组的想法，因为者们以来，我们会忽略掉一些跨了两个以上exon的reads
    vector<St_HitCase> vCases;

//    //1: Collect Hit Time and Info Order
//    cout << "Collect Hit Time and Info Order" << endl;
//    cout << "mpKT size: " << mpKT.size() << endl;
//    cout << "vChrom Size: " << vChrom.size() << endl;
//    cout << "We got size!" << endl;


//    int iTemp = 0;
//    for(map<unsigned int, vector<St_PosInfo> >::iterator itr = mpKT.begin(); itr != mpKT.end(); itr++)
//    {
//        if(iTemp > 10)
//            break;

//        cout << "---" << endl;
//        cout << itr->first << endl;
//        cout << IntToStr(itr->second.begin()->ucChromIndex) << " --- "
//             << IntToStr(itr->second.begin()->ucTranscriptIndex) << endl;
//        cout << "---" << endl;

//        iTemp++;
//    }

//    cout << "Collect Cases" << endl;
    for(unsigned int i = 0 ; i <= (strSeq.length() - KMERLEN + 1); i++)
    {
        string strCurKmer = strSeq.substr(i, KMERLEN);
        unsigned int uiKmer = ConvertKmerToNum32(strCurKmer);
        if(mpKT.find(uiKmer) == mpKT.end())
            continue;
        //Can find it
 //       cout << strCurKmer << endl;
        for(vector<St_PosInfo>::iterator itr = mpKT[uiKmer].begin(); itr != mpKT[uiKmer].end(); itr++)
        {
            if(vCases.empty()) // 第一次搞起
            {
                int iCutHitExonLen = vChrom.at(itr->ucChromIndex).vRG.at(itr->uiGeneIndex).vRT.at(itr->ucTranscriptIndex).vRExon.at(itr->ucExonIndex).GetLength();
                vCases.push_back(St_HitCase(*itr, 1, iCutHitExonLen));
            }
            else
            {
                bool bFind = false;
                for(vector<St_HitCase>::iterator itrCase = vCases.begin();
                    itrCase != vCases.end(); itrCase++)
                {
                    for(vector<St_HitExon>::iterator itrHE = itrCase->vHitExons.begin();
                        itrHE != itrCase->vHitExons.end(); itrHE++)
                    {
                        if( *itr == itrHE->stPI )
                        {
                            itrHE->iCount++;
                            itrHE->UpdateVPart(itr->cPart);
                            bFind = true;
                            break;
                        }
                    }
                    if(!bFind)
                    {
                        //如果是属于一个transcript,那么我们还是把你丢进去
                        if(itr->CheckSameTranscript(itrCase->vHitExons.at(0).stPI))
                        {
                            int iCutHitExonLen = vChrom.at(itr->ucChromIndex).vRG.at(itr->uiGeneIndex).vRT.at(itr->ucTranscriptIndex).vRExon.at(itr->ucExonIndex).GetLength();
                            itrCase->vHitExons.push_back(St_HitExon(*itr, 1, iCutHitExonLen));
                            bFind = true;
                        }
                    }

                    if(bFind)
                        break;
                }
                if(!bFind) //找不到的话，我们就插入
                {
                    int iCutHitExonLen = vChrom.at(itr->ucChromIndex).vRG.at(itr->uiGeneIndex).vRT.at(itr->ucTranscriptIndex).vRExon.at(itr->ucExonIndex).GetLength();
                    vCases.push_back(St_HitCase(*itr, 1, iCutHitExonLen));
               //     cout << "Get New Case: " << itr->ucChromIndex << " " << itr->ucTranscriptIndex << endl;
                }
            }
        }
    }

//    cout << "Step 1 finished" << endl;

#ifdef DEBUG
    ///这里我们debug一下--> 将包含两个目标exon的这样的reads输出 -->
    //我们在这里沿用如下逻辑
    //1: 寻找相应的index
    St_PosInfo stStartExonInfo;
    St_PosInfo stEndExonInfo;
    bool bFindBoth = false;
    int iChromIndex = 0;
    for(vector<St_Row_Chrom>::iterator itr = vChrom.begin(); itr != vChrom.end(); itr++)
    {
        int iGeneIndex = 0;
        for(vector<St_Raw_Gene>::iterator itrRG = itr->vRG.begin(); itrRG != itr->vRG.end(); itrRG++)
        {
            int iTranscriptIndex = 0;
            for(vector<St_Raw_Transcript>::iterator itrRT = itrRG->vRT.begin();
                itrRT != itrRG->vRT.end(); itrRT++)
            {
                bool bFindStartExon = false;
                bool bFindEndExon = false;
                int iExonIndex = 0;
                for(vector<St_Raw_Exon>::iterator itrRE = itrRT->vRExon.begin();
                    itrRE != itrRT->vRExon.end(); itrRE++)
                {
                    if(itrRE->iStart == 233270759 && itrRE->iEnd == 233270936)
                    {
                        stStartExonInfo.ucChromIndex = iChromIndex;
                        stStartExonInfo.uiGeneIndex = iGeneIndex;
                        stStartExonInfo.ucTranscriptIndex = iTranscriptIndex;
                        stStartExonInfo.ucExonIndex = iExonIndex;
                        bFindStartExon = true;
                    }

                    if(itrRE->iStart == 233275460 && itrRE->iEnd == 233275601)
                    {
                        stEndExonInfo.ucChromIndex = iChromIndex;
                        stEndExonInfo.uiGeneIndex = iGeneIndex;
                        stEndExonInfo.ucTranscriptIndex = iTranscriptIndex;
                        stEndExonInfo.ucExonIndex = iExonIndex;
                        bFindEndExon = true;
                    }

                    if(bFindStartExon && bFindEndExon)
                    {
                        bFindBoth = true;
                        break;
                    }
                    iExonIndex++;
                }
                if(bFindBoth)
                    break;
                iTranscriptIndex++;
            }
            if(bFindBoth)
                break;

            iGeneIndex++;
        }
        if(bFindBoth)
            break;
        iChromIndex++;
    }

    if(bFindBoth)
    {
        bFindBoth = false;
        int iStartExonKmerHitNum = 0;
        int iEndExonKmerHitNum = 0;
        //2: 看看是否存在一个被 hit 的transcript能够包含这两个exon，也就是说这两个exon被同时hit了
        for(vector<St_HitCase>::iterator itrCase = vCases.begin(); itrCase != vCases.end(); itrCase++)
        {
            bool bFindStartExon = false;
            bool bFindEndExon = false;
            for(vector<St_HitExon>::iterator itrHitE = itrCase->vHitExons.begin();
                itrHitE != itrCase->vHitExons.end(); itrHitE++)
            {
                if(itrHitE->stPI == stStartExonInfo)
                {
                    bFindStartExon= true;
                    iStartExonKmerHitNum = itrHitE->iCount;
                }

                if(itrHitE->stPI == stEndExonInfo)
                {
                    bFindEndExon= true;
                    iEndExonKmerHitNum = itrHitE->iCount;
                }

                if(bFindStartExon && bFindEndExon)
                {
                    bFindBoth = true;
                    //我们需要在这里输出hit的情况
                    for(vector<St_HitExon>::iterator itrTmp = itrCase->vHitExons.begin();
                        itrTmp != itrCase->vHitExons.end(); itrTmp++)
                    {
                        cout << IntToStr(itrTmp->stPI.ucExonIndex) << "("
                             << IntToStr(itrTmp->iCount) << "), ";
                    }
                    cout << endl;
                    break;
                }
            }
            if(bFindBoth)
                break;
        }
        if(bFindBoth)
        {
            cout << strSeq << endl;
            cout << IntToStr(iStartExonKmerHitNum) << "  " << IntToStr(iEndExonKmerHitNum) << endl;
            cout << endl;
        }
    }
    ///<-------
    ///
#endif
    /* 2:查看是否属于有效的hit，这是一个复杂的过程, 我们大致可以分为以下几步
     * (1) 提取最高的两组被hit的exon
     *
     */
    //2: 选取最高得票的那一个 -->找到hit kmer 最多的那个transcript    
    St_HitCase* pHitCase = NULL;
#ifdef DEBUG
    cout << "How many Hit --> start" << endl;
#endif
    for(vector<St_HitCase>::iterator itr = vCases.begin(); itr != vCases.end(); itr++)
    {

#ifdef DEBUG
        cout << IntToStr(itr->GetHitSum()) << " --> " ;
        //output the related exons
        for(vector<St_HitExon>::iterator itrExon = itr->vHitExons.begin();
            itrExon != itr->vHitExons.end(); itrExon++)
        {
            St_PosInfo& stPI = itrExon->stPI;
            St_Raw_Exon& stExon = vChrom.at(stPI.ucChromIndex).vRG.at(stPI.uiGeneIndex).vRT.at(stPI.ucTranscriptIndex).vRExon.at(stPI.ucExonIndex);
            cout << "<" << IntToStr(stExon.iStart) << ", " << IntToStr(stExon.iEnd) << "> - "
                 << IntToStr(stExon.GetLength()) << ", ";
        }
        cout << endl;
#endif

        if(pHitCase == NULL)
            pHitCase = &(*itr);
        else
        {            
            //这里的经验是需要同时满足两个条件
            //(1) Maximum number of hit
            //(2) shortest length

            //if((float)pHitCase->GetHitSum() / pHitCase->GetHitLen() < (float)itr->GetHitSum() / itr->GetHitLen())
            //    pHitCase = &(*itr);
            if(pHitCase->GetHitSum() <= (itr->GetHitSum() - 5)) // 满足条件1
                pHitCase = &(*itr);
            else if(abs(pHitCase->GetHitSum() - itr->GetHitSum()) < 5) //满足条件2
            {
                if(pHitCase->GetHitLen() > itr->GetHitLen()) // we need to find the shortest one
                {
                    //for the case only 1 exon been hit
                    if(itr->vHitExons.size() == 1)
                    {
                        St_PosInfo& stPI = itr->vHitExons[0].stPI;
                        St_Raw_Exon& stExon = vChrom.at(stPI.ucChromIndex).vRG.at(stPI.uiGeneIndex).vRT.at(stPI.ucTranscriptIndex).vRExon.at(stPI.ucExonIndex);
                        if(stExon.GetIsCRNAHeadExon() && stExon.GetIsCRNATailExon())
                        {
                            pHitCase = &(*itr);
                        }
                    }
                    else // for the case of two exons
                        pHitCase = &(*itr);
                }
            }
        }
    }

#ifdef DEBUG
    cout << "<-- End" << endl;
#endif

    if(pHitCase == NULL)
    {
#ifdef DEBUG
        cout << "Do not find pHitCase" << endl;
#endif
        return;
    }

#ifdef DEBUG
    cout << "Hit Sum: " << IntToStr(pHitCase->GetHitSum()) << endl;
#endif

    //1: filter the too small number of hitting
    //cout << "Erase Small" << endl;
    int iExtaSmall = 5;  //this threshold is OK! (only 1,2,3,4 will be erased)
    for( vector<St_HitExon>::iterator itr = pHitCase->vHitExons.end() - 1;
         itr >= pHitCase->vHitExons.begin(); itr--)
    {
        if(itr->iCount < iExtaSmall)
            pHitCase->vHitExons.erase(itr);
    }

#ifdef DEBUG
    //Cout Best Hit case
    cout << endl << "Best Hit Case --> " << endl;
    cout << IntToStr(pHitCase->GetHitSum()) << " --> " ;
    for(vector<St_HitExon>::iterator itrExon = pHitCase->vHitExons.begin();
        itrExon != pHitCase->vHitExons.end(); itrExon++)
    {
        St_PosInfo& stPI = itrExon->stPI;
        St_Raw_Exon& stExon = vChrom.at(stPI.ucChromIndex).vRG.at(stPI.uiGeneIndex).vRT.at(stPI.ucTranscriptIndex).vRExon.at(stPI.ucExonIndex);
        cout << "<" << IntToStr(stExon.iStart) << ", " << IntToStr(stExon.iEnd) << "> - "
             << IntToStr(stExon.GetLength()) << ", ";
    }
    cout << "<----" << endl;
#endif

    //2: Make filter based on the TWO threshoulds
    //Threashold 1: 针对hit多少去进行抉择: Get threshold --> fKmerRatio
    bool bTheasholdHitNum = CheckHitNum(strSeq, fKmerRatio, pHitCase);
    if(!bTheasholdHitNum)
        return;

    //Threashold 2: 针对是否好似前后
    bool bTheasholdHitPart = CheckHitPart(vChrom, pHitCase,
                                          strSeq, fKmerRatio, bSeqRC);
    if(!bTheasholdHitPart)
        return;

    /* 到了这一步可以认为，基本上是hit上了，但是是否符合条件还需要看下面两个非常重要的因素
     * (1) 过于小的hit给过滤掉
     * (2) 看看剩下的case是属于 self back 还是 regular back
     *   a) 如果是self back
     *     1) 首先 case里面的 exon只有一个
     *     2）看看这个exon的首位是不是符合我们的signal
     *     3）看看这个exon的长度是不是小于我们的reads
     *     4) 符合条件1，2，3那么我们认为符合我们的预期属于相应的candidate
     *
     *   b)看是否是regular circular RNA
     *     1) 寻找翻转junction处的两个不同的exon
     *     2）前面的是否符合tail signal， 后面的是否符合head siginal
     *     3) 考虑这是一个符合预期的candiate 如果同时符合1，2条件
     */

    //3: Check How many exons in each Case
    //cout << "Check How many exons in each Case" << endl;
    if(pHitCase->vHitExons.size() == 1)
    {
#ifdef DEBUG
        cout << "CheckSelfCircRNA>>>>>>>" << endl;
#endif
        // 3.1: Potential Self-circular RNA
        bool bFind = CheckSelfCircRNA(pHitCase, strSeq.length(), vChrom, vSelfCircCandi);

        //For debug        
        if(bFind)
        {
            /*
            St_PosInfo& stPI = pHitCase->vHitExons[0].stPI;
            St_Raw_Exon& stCurExon = vChrom.at(stPI.ucChromIndex).vRG.at(stPI.uiGeneIndex).vRT.at(stPI.ucTranscriptIndex).vRExon.at(stPI.ucExonIndex);
            if(stCurExon.iStart == 26156109 && stCurExon.iEnd == 26156321)
            {
                cout << "<26156109, 26156321>: " << strSeq << endl;
            }
            */
        }
    }
    else if(pHitCase->vHitExons.size() > 1)
    {
#ifdef DEBUG
        cout << "CheckRegularCircRNA<<<<<<" << endl;
#endif
        // 3.2: Potential Regular-Circular RNA
        St_Candidate stCandiRecord;
        bool bFind = CheckRegularCircRNA(pHitCase, vChrom, stCandiRecord, vRegCircCandi);

        //For Debug

        if(bFind)
        {            
            //m_iTotalNum++;
            //if(!bSeqRC)
            //    m_iRegNum++;

            St_PosInfo& stPI1 = pHitCase->vHitExons[0].stPI;
            St_Raw_Exon& stCurExon1 = vChrom.at(stPI1.ucChromIndex).vRG.at(stPI1.uiGeneIndex).vRT.at(stPI1.ucTranscriptIndex).vRExon.at(stPI1.ucExonIndex);

            St_PosInfo& stPI2 = pHitCase->vHitExons[1].stPI;
            St_Raw_Exon& stCurExon2 = vChrom.at(stPI2.ucChromIndex).vRG.at(stPI2.uiGeneIndex).vRT.at(stPI2.ucTranscriptIndex).vRExon.at(stPI2.ucExonIndex);


//            cout << "Donor(Start, End) & Acceptor(Start, End): "
//                 << "(" << IntToStr(stCurExon1.iStart) << ", " << IntToStr(stCurExon1.iEnd) << ")"
//                 << " -- "
//                 << "(" << IntToStr(stCurExon2.iStart) << ", " << IntToStr(stCurExon2.iEnd) << ")" << endl;

//            cout << "Candi Start & End: "<< IntToStr(stCandiRecord.iStartPos) << " : "
//                 << IntToStr(stCandiRecord.iEndPos) << endl;

//            cout << "RC: " << (stCandiRecord.bRC ? "Yes" : "No") << endl;

//            cout << "Current Reads: " << (bSeqRC ? "RC" : "Reg") << endl
//                 << strSeq << endl << endl;

//            if( stCurExon1.iStart == 46659946 && stCurExon1.iEnd == 46660073 &&  //这里的是小的
//                stCurExon2.iStart == 46660225 && stCurExon2.iEnd == 46660323 )   //这里的是大的
//            {
//                cout << IntToStr(stPI1.ucChromIndex) << ", " << IntToStr(stPI1.uiGeneIndex) << ", "
//                     << IntToStr(stPI1.ucTranscriptIndex) << ", " << IntToStr(stPI1.ucExonIndex) << endl;
//                cout << IntToStr(stPI2.ucChromIndex) << ", " << IntToStr(stPI2.uiGeneIndex) << ", "
//                     << IntToStr(stPI2.ucTranscriptIndex) << ", " << IntToStr(stPI2.ucExonIndex) << endl;
//                cout << IntToStr(abs(46660225 - 46660073)) << endl;
//                cout << "<46659946, 46660073>: " << strSeq << endl;
//                cout << strSeq << endl;
//            }
        }
    }
    else
    {}

    //release resources
    vCases.clear();
}

bool ClsFindCandidate::CheckHitNum(string& strSeq, float fKmerRatio, St_HitCase* pHitCase)
{
    if(pHitCase == NULL)
        return false;

    int iBoundLen = strSeq.length() * fKmerRatio; //This could be viewed as the number of kmer
    int iLimit = 2 * iBoundLen + KMERLEN*2 - 1; //This is the sequence length which used to generate kmer
    int iMinHit = 0;
    if(pHitCase->GetHitLen() <= iLimit)
        iMinHit = pHitCase->GetHitLen() - KMERLEN + 5; // just a little bit larger，从而保证这个exon能够被覆盖一次以上
                                                       // 这个重要
    else
    {
        iMinHit = iBoundLen * 1.2;
    }

    //int iMinHit = (strSeq.length() - (2*KMERLEN + 1)) * .6; //.7; //之前是0.5
    if(pHitCase->GetHitSum() < iMinHit)
        return false;
    else
        return true;
}

bool ClsFindCandidate::CheckHitPart(vector<St_Row_Chrom>& vChrom, St_HitCase* pHitCase,
                                    string& strSeq, float fKmerRatio, bool bSeqRC)
{
    if(pHitCase == NULL || pHitCase->vHitExons.empty())
        return false;

    /*
    cout << "pHitCase->vExons size: " << IntToStr(pHitCase->vHitExons.size()) <<  endl;
    for(vector<St_HitExon>::iterator itr = pHitCase->vHitExons.begin();
        itr != pHitCase->vHitExons.end(); itr++)
    {
        St_PosInfo& stPI1 = itr->stPI;
        St_Raw_Exon& stCurExon1 = vChrom.at(stPI1.ucChromIndex).vRG.at(stPI1.uiGeneIndex).vRT.at(stPI1.ucTranscriptIndex).vRExon.at(stPI1.ucExonIndex);
        cout << stCurExon1.GetLength() << " -- Hit Num: " << itr->iCount << endl;
    }*/

    //-->Just for erase the warnning
    if(bSeqRC)
    {}
    //<--

    int iBoundLen = strSeq.length() * fKmerRatio;
    int iLimit = 2 * iBoundLen + KMERLEN*2 - 1;

    //如果是1个exon
    if(pHitCase->vHitExons.size() == 1)
    {
        if(pHitCase->GetHitLen() <= iLimit)  //短的我们就直接发返回了
            return true;
        else
        {
            //1:Get this exon
            St_PosInfo& stPI = pHitCase->vHitExons[0].stPI;
            St_Raw_Exon& stCurExon = vChrom.at(stPI.ucChromIndex).vRG.at(stPI.uiGeneIndex).vRT.at(stPI.ucTranscriptIndex).vRExon.at(stPI.ucExonIndex);
            if(!stCurExon.bRC) // 如果是正向序列 --> 在此我们正向反向都这么去做，看看最后结果如何
            {
                // 我觉得还是要首先看长度 --> 在长度确定之后，我们再去观察相应的hitting part
                // 因为现在，无论exon的长短，我们都会有相应的 hitting part信息


                //对于正向序列，我们的vPart就需要是 EEEE --> SSSS
                vector<char>& vPart = pHitCase->vHitExons[0].vPart;
                vector<char> vTag;
                vector<int> vNum;
                for(vector<char>::iterator itr = vPart.begin(); itr != vPart.end(); itr++)
                {
                    if(vTag.empty())
                    {
                        vTag.push_back(*itr);
                        vNum.push_back(1);
                    }
                    else
                    {
                        if( (*(vTag.end() - 1)) == *itr )
                            (*(vNum.end() -1))++;
                        else
                        {
                            vTag.push_back(*itr);
                            vNum.push_back(1);
                        }
                    }
                }

                for(int i = vNum.size() - 1; i >= 0; i--) //erase the random hit
                {
                    if(vNum.at(i) <= 3)  // must > 3
                    {
                        vNum.erase(vNum.begin() + i);
                        vTag.erase(vTag.begin() + i);
                    }
                }

                //我们现在得到了浓缩后的了，那么我们就看看是不是make sense
                if(vTag.size() == 2)
                {
                    //if(!bSeqRC) // normal direction 注意此处通过测试的结果显示，跟direction没有关系
                    if(vTag[0] == 'E' && vTag[1] == 'S')
                        return true;
                    else
                        return false;
                }
                //else if(vTag.size() == 1)
                //{
                //    if(vTag[0] == 'U')
                //        return true;
                //    else
                //        return false;
                //}
                else
                {
                    return false;
                }
            }
            else if(stCurExon.bRC) //如果是负向序列
            {                
                // 我觉得还是要首先看长度 --> 在长度确定之后，我们再去观察相应的hitting part
                // 因为现在，无论exon的长短，我们都会有相应的 hitting part信息
                //对于正向序列，我们的vPart就需要是 EEEE --> SSSS
                vector<char>& vPart = pHitCase->vHitExons[0].vPart;
                vector<char> vTag;
                vector<int> vNum;
                for(vector<char>::iterator itr = vPart.begin(); itr != vPart.end(); itr++)
                {
                    if(vTag.empty())
                    {
                        vTag.push_back(*itr);
                        vNum.push_back(1);
                    }
                    else
                    {
                        if( (*(vTag.end() - 1)) == *itr )
                            (*(vNum.end() -1))++;
                        else
                        {
                            vTag.push_back(*itr);
                            vNum.push_back(1);
                        }
                    }
                }

                for(int i = vNum.size() - 1; i >= 0; i--) //erase the random hit
                {
                    if(vNum.at(i) <= 3)  // must > 3
                    {
                        vNum.erase(vNum.begin() + i);
                        vTag.erase(vTag.begin() + i);
                    }
                }

                //我们现在得到了浓缩后的了，那么我们就看看是不是make sense
                if(vTag.size() == 2)
                {
                    //if(!bSeqRC) // normal direction 注意此处通过测试的结果显示，跟direction没有关系
                    if(vTag[0] == 'E' && vTag[1] == 'S')
                        return true;
                    else
                        return false;
                }

                //else if(vTag.size() == 1)
                //{
                //    if(vTag[0] == 'U')
                //        return true;
                //    else
                //        return false;
                //}
                else
                {
                    return false;
                }
            }
            else
                return false;
        }
    }    
    else if(pHitCase->vHitExons.size() == 2)
    {                
        St_PosInfo& stPI1 = pHitCase->vHitExons[0].stPI;
        St_Raw_Exon& stCurExon1 = vChrom.at(stPI1.ucChromIndex).vRG.at(stPI1.uiGeneIndex).vRT.at(stPI1.ucTranscriptIndex).vRExon.at(stPI1.ucExonIndex);

        St_PosInfo& stPI2 = pHitCase->vHitExons[1].stPI;
        St_Raw_Exon& stCurExon2 = vChrom.at(stPI2.ucChromIndex).vRG.at(stPI2.uiGeneIndex).vRT.at(stPI2.ucTranscriptIndex).vRExon.at(stPI2.ucExonIndex);

        // 我们在这里同样进行一个check，限制reads的比对，需要是先匹配上后面的再匹配上前面的 --> 搞起
        // 这个策略跟 self circle的策略是一样的
        if(!stCurExon1.bRC && !stCurExon2.bRC)  // 如果是正向序列            
        {            
            //对于正向序列，我们的vPart就需要是 EEEE --> SSSS
            //For the first one:
            vector<char>& vPart1 = pHitCase->vHitExons[0].vPart;
            vector<char> vTag1;
            vector<int> vNum1;
            for(vector<char>::iterator itr = vPart1.begin(); itr != vPart1.end(); itr++)
            {
                if(vTag1.empty())
                {
                    vTag1.push_back(*itr);
                    vNum1.push_back(1);
                }
                else
                {
                    if( (*(vTag1.end() - 1)) == *itr )
                        (*(vNum1.end() -1))++;
                    else
                    {
                        vTag1.push_back(*itr);
                        vNum1.push_back(1);
                    }
                }
            }
            /// make filter for this raw tag vector
            /// 1: the tag number (avoid the ramdon tag)
            for(int i = vNum1.size() - 1; i >= 0; i--)
            {
                if(vNum1.at(i) <= 3)  // must > 3
                {
                    vNum1.erase(vNum1.begin() + i);
                    vTag1.erase(vTag1.begin() + i);
                }
            }
            ///2: transfer "S+E" --> "E"
            if(vTag1.size() == 2)
            {
                if(vTag1[0] == 'S' && vTag1[1] == 'E')
                {
                    vTag1.erase(vTag1.begin());
                    vNum1.erase(vNum1.begin());
                }
            }

            //For the second one:
            vector<char>& vPart2 = pHitCase->vHitExons[1].vPart;
            vector<char> vTag2;
            vector<int> vNum2;
            for(vector<char>::iterator itr = vPart2.begin(); itr != vPart2.end(); itr++)
            {
                if(vTag2.empty())
                {
                    vTag2.push_back(*itr);
                    vNum2.push_back(1);
                }
                else
                {
                    if( (*(vTag2.end() - 1)) == *itr )
                        (*(vNum2.end() -1))++;
                    else
                    {
                        vTag2.push_back(*itr);
                        vNum2.push_back(1);
                    }
                }
            }

            /// make filter for this raw tag vector
            /// 1: the tag number (avoid the ramdon tag)
            for(int i = vNum2.size() - 1; i >= 0; i--)
            {
                if(vNum2.at(i) <= 3)  // must > 3
                {
                    vNum2.erase(vNum2.begin() + i);
                    vTag2.erase(vTag2.begin() + i);
                }
            }
            ///2: transfer "S+E" --> "S"
            if(vTag2.size() == 2)  // we have combined  --> //check length
            {
                if(vTag2[0] == 'S' && vTag2[1] == 'E')
                {
                    vTag2.erase(vTag2.end()-1);
                    vNum2.erase(vNum2.end()-1);
                }
            }

            if(vTag1.size() == 1 && vTag2.size() == 1)
            {
                //return true;

                if(vTag1[0] == 'E' && vTag2[0] == 'E')
                    return false;
                else if(vTag1[0] == 'E' && vTag2[0] == 'S')
                    return true;
                else if(vTag1[0] == 'S' && vTag2[0] == 'S')
                    return false;
                else if(vTag1[0] == 'S' && vTag2[0] == 'E')
                    return false;
                else                
                    return false;
            }
            else
                return false; //new tag should satisfy all previous requirement!!!
        }
        else if (stCurExon1.bRC && stCurExon2.bRC) //负向序列看看 -->
        {            
            //对于正向序列，我们的vPart就需要是 EEEE --> SSSS
            //For the first one:
            vector<char>& vPart1 = pHitCase->vHitExons[0].vPart;
            vector<char> vTag1;
            vector<int> vNum1;
            for(vector<char>::iterator itr = vPart1.begin(); itr != vPart1.end(); itr++)
            {
                if(vTag1.empty())
                {
                    vTag1.push_back(*itr);
                    vNum1.push_back(1);
                }
                else
                {
                    if( (*(vTag1.end() - 1)) == *itr )
                        (*(vNum1.end() -1))++;
                    else
                    {
                        vTag1.push_back(*itr);
                        vNum1.push_back(1);
                    }
                }
            }

            /// make filter for this raw tag vector
            /// 1: the tag number (avoid the ramdon tag)
            for(int i = vNum1.size() - 1; i >= 0; i--)
            {
                if(vNum1.at(i) <= 3)  // must > 3
                {
                    vNum1.erase(vNum1.begin() + i);
                    vTag1.erase(vTag1.begin() + i);
                }
            }
            //2: transfer "S+E" --> "E"
            if(vTag1.size() == 2)
            {
                if(vTag1[0] == 'S' && vTag1[1] == 'E')
                {
                    vTag1.erase(vTag1.begin());
                    vNum1.erase(vNum1.begin());
                }
            }

            //For the second one:
            vector<char>& vPart2 = pHitCase->vHitExons[1].vPart;
            vector<char> vTag2;
            vector<int> vNum2;
            for(vector<char>::iterator itr = vPart2.begin(); itr != vPart2.end(); itr++)
            {
                if(vTag2.empty())
                {
                    vTag2.push_back(*itr);
                    vNum2.push_back(1);
                }
                else
                {
                    if( (*(vTag2.end() - 1)) == *itr )
                        (*(vNum2.end() -1))++;
                    else
                    {
                        vTag2.push_back(*itr);
                        vNum2.push_back(1);
                    }
                }
            }

            /// make filter for this raw tag vector
            /// 1: the tag number (avoid the ramdon tag)
            for(int i = vNum2.size() - 1; i >= 0; i--)
            {
                if(vNum2.at(i) <= 3)  // must > 3
                {
                    vNum2.erase(vNum2.begin() + i);
                    vTag2.erase(vTag2.begin() + i);
                }
            }
            //2: transfer "S+E" --> "S"
            if(vTag2.size() == 2)  // we have combined  --> //check length
            {
               if(vTag2[0] == 'S' && vTag2[1] == 'E')
                {
                    vTag2.erase(vTag2.end()-1);
                    vNum2.erase(vNum2.end()-1);
                }
            }

            if(vTag1.size() == 1 && vTag2.size() == 1)
            {
                //return true;

                if(vTag1[0] == 'E' && vTag2[0] == 'E')
                    return false;
                else if(vTag1[0] == 'E' && vTag2[0] == 'S')
                    return true;
                else if(vTag1[0] == 'S' && vTag2[0] == 'S')
                    return false;
                else if(vTag1[0] == 'S' && vTag2[0] == 'E')
                    return false;
                else
                    return false;
            }
            else
                return false; //false; //new tag should satisfy all previous requirement!!!
        }
        else
            return false; // 这个case是一个负一个正，这个case的出现本身就是荒谬的
    }
    else
    {
        return false;
    }
}

bool ClsFindCandidate::CheckSelfCircRNA(St_HitCase* pHitCase, int iReadLen,
                                        vector<St_Row_Chrom>& vChrom, vector<St_Candidate>& vSelfCircCandi)
{    
    St_PosInfo& stPI = pHitCase->vHitExons[0].stPI;
    St_Raw_Exon& stCurExon = vChrom.at(stPI.ucChromIndex).vRG.at(stPI.uiGeneIndex).vRT.at(stPI.ucTranscriptIndex).vRExon.at(stPI.ucExonIndex);

    //我们看看如果相应的exon短于我们的reads的长度，那么我们需要更多的hit才算ok
    if(stCurExon.GetLength() < iReadLen) //至少我们的Hitting 要大于exon本身的长度
    {
        if(pHitCase->GetHitSum() < (stCurExon.GetLength() - KMERLEN + 1))
            return false;
    }

    //1: check the signal
    bool bSupportHead = stCurExon.GetIsCRNAHeadExon();
    bool bSupportTail = stCurExon.GetIsCRNATailExon();

#ifdef DEBUG
    if(!bSupportHead)
        cout << "Head Do not support"  << endl;

    if(!bSupportTail)
        cout << "Tail Do not support"  << endl;
#endif

    if(bSupportHead && bSupportTail)
    {
#ifdef DEBUG
        cout << "Both Support" << endl;
        cout << "Boundary: <" << IntToStr(stCurExon.iStart) << ", "
             << IntToStr(stCurExon.iEnd) << ">" << endl;
        cout << "Both Hit" << endl;
#endif

        //if(stCurExon.GetLength() < iReadLen)
        //{
            //我们认为是一个候选者
            //注意这个地方: 因为我们是找circular rna的jucntion，因此circular rna的start应该是原始exon的tail，反之毅然
            St_Candidate stCurCandi(stPI.ucChromIndex, stCurExon.iEnd, stCurExon.iStart, 1, stCurExon.bRC);
            stCurCandi.SetCircType(ctSelf);

            //--->For vector
            bool bFind = false;
            for(vector<St_Candidate>::iterator itr = vSelfCircCandi.begin();
                itr != vSelfCircCandi.end(); itr++)
            {
                if(*itr == stCurCandi)
                {
                    itr->iSupportNum++;
                    bFind = true;
                    break;
                }
            }
            if(!bFind) //找不到
            {
                vSelfCircCandi.push_back(stCurCandi);

                /*
                //outout the tag result
                St_PosInfo& stPI = pHitCase->vHitExons[0].stPI;
                St_Raw_Exon& stCurExon = vChrom.at(stPI.ucChromIndex).vRG.at(stPI.uiGeneIndex).vRT.at(stPI.ucTranscriptIndex).vRExon.at(stPI.ucExonIndex);
                if(!stCurExon.bRC) // 如果是正向序列
                {
                    //对于正向序列，我们的vPart就需要是 EEEE --> SSSS
                    vector<char>& vPart = pHitCase->vHitExons[0].vPart;
                    vector<char> vTag;
                    vector<int> vNum;
                    for(vector<char>::iterator itr = vPart.begin(); itr != vPart.end(); itr++)
                    {
                        if(vTag.empty())
                        {
                            vTag.push_back(*itr);
                            vNum.push_back(1);
                        }
                        else
                        {
                            if( (*(vTag.end() - 1)) == *itr )
                                (*(vNum.end() -1))++;
                            else
                            {
                                vTag.push_back(*itr);
                                vNum.push_back(1);
                            }
                        }
                    }

                    //我们现在得到了浓缩后的了，那么我们就看看是不是make sense
                    if(vTag.size() == 2)
                    {
                    }
                    else if(vTag.size() == 1)
                    {
                        if(vTag[0] == 'U')
                        {
                            cout << IntToStr(stCurExon.GetLength())
                                 << "----------------------> self_Circ" << endl;
                        }
                    }
                }*/
            }
        //}
        //<---
            /*
            //我们认为是一个候选者
            //注意这个地方: 因为我们是找circular rna的jucntion，因此circular rna的start应该是原始exon的tail，反之毅然
            St_Candidate stCurCandi(stPI.ucChromIndex, stCurExon.iEnd, stCurExon.iStart);
            if(m_mpSelfCircCandi.find(stCurCandi) == m_mpSelfCircCandi.end())
                m_mpSelfCircCandi[stCurCandi] = 1;
            else
                m_mpSelfCircCandi[stCurCandi]++;

            if(m_mpSelfCircCandi[stCurCandi] == 0)
                cout << "Why is 0 !!!" << endl;
        //}*/
        return true;
    }
    else
        return false;
}

bool ClsFindCandidate::CheckRegularCircRNA(St_HitCase* pHitCase,
                                           vector<St_Row_Chrom>& vChrom,
                                           St_Candidate& stCandiRecord, vector<St_Candidate>& vRegCircCandi)
{
    //Get the direction of gene    
    int iChromIndex = pHitCase->vHitExons.begin()->stPI.ucChromIndex;
    int iGeneIndex = pHitCase->vHitExons.begin()->stPI.uiGeneIndex;
    bool bRC = vChrom[iChromIndex].vRG[iGeneIndex].bRC;

    //1: find out the junction point (exon index decreased)
    St_PosInfo* pHeadPI = NULL;
    St_Raw_Exon* pHeadExon = NULL;

    St_PosInfo* pTailPI = NULL;
    St_Raw_Exon* pTailExon = NULL;

    if(!bRC) //正向序列
    {
        unsigned int uiCurEI = 0; //EI: Exon Index

#ifdef DEBUG
        cout << "!bRC" << endl;
#endif
        for(vector<St_HitExon>::iterator itr = pHitCase->vHitExons.begin();
            itr != pHitCase->vHitExons.end(); itr++)
        {

#ifdef DEBUG
            cout << IntToStr(itr->stPI.ucExonIndex) << "("
                 << IntToStr(itr->iCount) << ")" << endl;
#endif

            if(itr->stPI.ucExonIndex >= uiCurEI)
            {
                uiCurEI = itr->stPI.ucExonIndex;
                continue;
            }
            else //Find the junction point
            {
                //For Head Exon
                pHeadPI = &(itr-1)->stPI;
                pHeadExon = &vChrom.at(pHeadPI->ucChromIndex).vRG.at(pHeadPI->uiGeneIndex).vRT.at(pHeadPI->ucTranscriptIndex).vRExon.at(pHeadPI->ucExonIndex);
#ifdef DEBUG
                cout << "Head (Head_tag, Tail_tag): "
                     << pHeadExon->strHead2Bp << ", " << pHeadExon->strTail2Bp
                     << " " << IntToStr(pHeadExon->iStart) << " " << IntToStr(pHeadExon->iEnd) << endl;
#endif

                //For Tail Exon
                pTailPI = &itr->stPI;
                pTailExon = &vChrom.at(pTailPI->ucChromIndex).vRG.at(pTailPI->uiGeneIndex).vRT.at(pTailPI->ucTranscriptIndex).vRExon.at(pTailPI->ucExonIndex);
#ifdef DEBUG
                cout << "Tail (Head_tag, Tail_tag): "
                     << pTailExon->strHead2Bp << ", " << pTailExon->strTail2Bp
                     << " " << IntToStr(pTailExon->iStart) << " " << IntToStr(pTailExon->iEnd) << endl;
#endif
                break;
            }
        }
    }
    else // 反向序列 在这个里面 exon是从 “大的位置” 到 “小的位置”排下来的！！！
    {
        unsigned int uiCurEI = 999; //EI: Exon Index
        /*
        pHeadPI = &pHitCase->vHitExons[0].stPI;
        pHeadExon = &vChrom.at(pHeadPI->ucChromIndex).vRG.at(pHeadPI->uiGeneIndex).vRT.at(pHeadPI->ucTranscriptIndex).vRExon.at(pHeadPI->ucExonIndex);

        pTailPI = &pHitCase->vHitExons[1].stPI;
        pTailExon = &vChrom.at(pTailPI->ucChromIndex).vRG.at(pTailPI->uiGeneIndex).vRT.at(pTailPI->ucTranscriptIndex).vRExon.at(pTailPI->ucExonIndex);
        */

#ifdef DEBUG
        cout << "bRC: " << endl;
#endif           
        for(vector<St_HitExon>::iterator itr = pHitCase->vHitExons.begin();
            itr != pHitCase->vHitExons.end(); itr++)
        {
            if(itr->stPI.ucExonIndex <= uiCurEI) // 这样就一路从大到小了 -->
            {
                uiCurEI = itr->stPI.ucExonIndex;
                continue;
            }
            else //Find the junction point  这里我们发现了一个从小到大的 -->
            {
                //For Head Exon
                pHeadPI = &(itr-1)->stPI;
                pHeadExon = &vChrom.at(pHeadPI->ucChromIndex).vRG.at(pHeadPI->uiGeneIndex).vRT.at(pHeadPI->ucTranscriptIndex).vRExon.at(pHeadPI->ucExonIndex);
#ifdef DEBUG
        cout << "Head (Head_tag, Tail_tag): "
             << pHeadExon->strHead2Bp << ", " << pHeadExon->strTail2Bp << endl;
#endif
                //For Tail Exon
                pTailPI = &itr->stPI;
                pTailExon = &vChrom.at(pTailPI->ucChromIndex).vRG.at(pTailPI->uiGeneIndex).vRT.at(pTailPI->ucTranscriptIndex).vRExon.at(pTailPI->ucExonIndex);
#ifdef DEBUG
        cout << "Tail (Head_tag, Tail_tag): "
             << pTailExon->strHead2Bp << ", " << pTailExon->strTail2Bp << endl;
#endif
                //if(abs(pHeadPI->ucExonIndex - pTailPI->ucExonIndex) == 1)
                //    return false;
                //else
                    break;
            }
        }
    }

    if(pTailPI == NULL || pHeadPI == NULL)
    {
#ifdef DEBUG
        if(pTailPI == NULL)
            cout << "Do not find Tail" << endl;
        if(pHeadPI == NULL)
            cout << "Do not find Head" << endl;
#endif
        return false;
    }

    //2: 我们的预期应该是尾巴连到头上面，于是如下
    //Check Signal
//    bool bSupportHead = pHeadExon->GetIsCRNAHeadExon();
//    bool bSupportTail = pTailExon->GetIsCRNATailExon();
    bool bSupportHead = pTailExon->GetIsCRNAHeadExon();
    bool bSupportTail = pHeadExon->GetIsCRNATailExon();

    if(bRC)
    {
        bSupportHead = pTailExon->GetIsCRNAHeadExon();
        bSupportTail = pHeadExon->GetIsCRNATailExon();
    }

#ifdef DEBUG
    if(!bSupportHead)
        cout << "Head Do not support"  << endl;

    if(!bSupportTail)
        cout << "Tail Do not support"  << endl;
#endif

    if(bSupportHead && bSupportTail)
    {
#ifdef DEBUG
        cout << "Both Support" << endl;
        cout << "Start --> <" << IntToStr(pHeadExon->iStart) << ", " << IntToStr(pHeadExon->iEnd) << ">" << endl;
        cout << "End   --> <" << IntToStr(pTailExon->iStart) << ", " << IntToStr(pTailExon->iEnd) << ">" << endl;
        //我们认为是一个候选者
        //注意这个地方: 因为我们是找circular rna的jucntion，因此circular rna的start应该是原始exon的tail，反之毅然
        if(pTailExon->iEnd == 17742)
        {
            cout << "AAAAAA" << endl;
        }
#endif       

        St_Candidate stCurCandi;
        stCurCandi.bRC = bRC;
        stCurCandi.SetCircType(ctRegular);
        if(!bRC) //如果是正向
        {                        
            stCurCandi.SetChromStartEnd(pHeadPI->ucChromIndex, pHeadExon->iEnd, pTailExon->iStart);
        }
        else //反向
        {
            stCurCandi.SetChromStartEnd(pHeadPI->ucChromIndex, pHeadExon->iEnd, pTailExon->iStart);
        }

        //--->For vector
        bool bFind = false;
        for(vector<St_Candidate>::iterator itr = vRegCircCandi.begin();
            itr != vRegCircCandi.end(); itr++)
        {
            if(*itr == stCurCandi)
            {
                itr->iSupportNum++;
                bFind = true;
                break;
            }
        }
        if(!bFind) //找不到
        {                       
            //if(!bRC) //正向的Exon
            {
                //UpdateCurCandiForRegularDirection(pHitCase, vChrom, stCurCandi);
            }
            vRegCircCandi.push_back(stCurCandi);
        }

        stCandiRecord = stCurCandi;
        //<---

        /*
        if(m_mpRegCircCandi.find(stCurCandi) == m_mpRegCircCandi.end())  // 这里都是原来的老的map的结构体，现在摒弃了
        {
            m_mpRegCircCandi[stCurCandi] = 1;

            if(stCurCandi.iStartPos == 148328392 && stCurCandi.iEndPos == 148253553)
                cout << "New 148328392, 148253553" << endl;

            bFind = false;
            for(map<St_Candidate, int>::iterator itr = m_mpRegCircCandi.begin(); itr != m_mpRegCircCandi.end(); itr++)
            {
                if(itr->first.iStartPos == 148328392 && itr->first.iEndPos == 148253553)
                {
                    bFind = true;
                    break;
                }
            }
            if(bFind)
                cout << "Find itr->first.iStartPos == 148328392 && itr->first.iEndPos == 148253553" << " --> " << IntToStr(m_mpRegCircCandi.size()) << " --> "
                     << IntToStr(stCurCandi.iStartPos) << ", " << IntToStr(stCurCandi.iEndPos) << endl;
            else
                cout << "Failed: " << IntToStr(stCurCandi.iStartPos) << ", " << IntToStr(stCurCandi.iEndPos) << " --> " << IntToStr(m_mpRegCircCandi.size()) << endl;*/

            /*
            if(m_mpRegCircCandi.size() == 4034)
            {
                ofstream ofs;
                ofs.open("./4034.txt");
                for(map<St_Candidate, int>::iterator itr = m_mpRegCircCandi.begin(); itr != m_mpRegCircCandi.end(); itr++)
                {
                    ofs << IntToStr(itr->first.iStartPos) << ", " << IntToStr(itr->first.iEndPos) << endl;
                }
                ofs.close();
            }

            if(m_mpRegCircCandi.size() == 4035)
            {
                ofstream ofs;
                ofs.open("./4035.txt");
                for(map<St_Candidate, int>::iterator itr = m_mpRegCircCandi.begin(); itr != m_mpRegCircCandi.end(); itr++)
                {
                    ofs << IntToStr(itr->first.iStartPos) << ", " << IntToStr(itr->first.iEndPos) << endl;
                }
                ofs.close();
            }*/
        /*}
        else
        {
            m_mpRegCircCandi[stCurCandi]++;*/

            /*
            bool bFind = false;
            for(map<St_Candidate, int>::iterator itr = m_mpRegCircCandi.begin(); itr != m_mpRegCircCandi.end(); itr++)
            {
                if(itr->first.iStartPos == 21221994 && itr->first.iEndPos == 21220010)
                {
                    bFind = true;
                    break;
                }
            }

            if(stCurCandi.iStartPos == 21221994 && stCurCandi.iEndPos == 21220010)
                cout << "21221994++" << " " << IntToStr(m_mpRegCircCandi.size());
            else
                cout << "++" << " " << IntToStr(m_mpRegCircCandi.size());

            if(bFind)
                cout << " Y" << endl;
            else
                cout << " NN" << endl;*/

        //}
        return true;
    }
    else
        return false;
}

void ClsFindCandidate::UpdateCurCandiForRegularDirection(St_HitCase* pHitCase,
                                                         vector<St_Row_Chrom>& vChrom,
                                                         St_Candidate& stCurCandi)
{
    //记录相应的情况
    if(pHitCase->vHitExons.size() == 2)
    {
        St_PosInfo& stPI1 = pHitCase->vHitExons[0].stPI;
        St_Raw_Exon& stCurExon1 = vChrom.at(stPI1.ucChromIndex).vRG.at(stPI1.uiGeneIndex).vRT.at(stPI1.ucTranscriptIndex).vRExon.at(stPI1.ucExonIndex);

        St_PosInfo& stPI2 = pHitCase->vHitExons[1].stPI;
        St_Raw_Exon& stCurExon2 = vChrom.at(stPI2.ucChromIndex).vRG.at(stPI2.uiGeneIndex).vRT.at(stPI2.ucTranscriptIndex).vRExon.at(stPI2.ucExonIndex);

        // 我们在这里同样进行一个check，限制reads的比对，需要是先匹配上后面的再匹配上前面的 --> 搞起
        // 这个策略跟 self circle的策略是一样的
        if( (!stCurExon1.bRC && !stCurExon2.bRC) ||
            (stCurExon1.bRC && stCurExon2.bRC) ) // 如果是正向序列
        {
             //对于正向序列，我们的vPart就需要是 EEEE --> SSSS
             //For the first one:
             vector<char>& vPart1 = pHitCase->vHitExons[0].vPart;
             vector<char> vTag1;
             vector<int> vNum1;
             for(vector<char>::iterator itr = vPart1.begin(); itr != vPart1.end(); itr++)
             {
                 if(vTag1.empty())
                 {
                     vTag1.push_back(*itr);
                     vNum1.push_back(1);
                 }
                 else
                 {
                     if( (*(vTag1.end() - 1)) == *itr )
                         (*(vNum1.end() -1))++;
                     else
                     {
                         vTag1.push_back(*itr);
                         vNum1.push_back(1);
                     }
                 }
             }

             /*
             /// make filter for this raw tag vector
             /// 1: the tag number (avoid the ramdon tag)
             for(int i = vNum1.size() - 1; i >= 0; i--)
             {
                 if(vNum1.at(i) <= 3)  // must > 3
                 {
                     vNum1.erase(vNum1.begin() + i);
                     vTag1.erase(vTag1.begin() + i);
                 }
             }
             ///2: transfer "S+E" --> "E"
             if(vTag1.size() == 2)
             {
                 if(vTag1[0] == 'S' && vTag1[1] == 'E')
                 {
                     vTag1.erase(vTag1.begin());
                     vNum1.erase(vNum1.begin());
                 }
             }*/

             //For the second one:
             vector<char>& vPart2 = pHitCase->vHitExons[1].vPart;
             vector<char> vTag2;
             vector<int> vNum2;
             for(vector<char>::iterator itr = vPart2.begin(); itr != vPart2.end(); itr++)
             {
                 if(vTag2.empty())
                 {
                     vTag2.push_back(*itr);
                     vNum2.push_back(1);
                 }
                 else
                 {
                     if( (*(vTag2.end() - 1)) == *itr )
                         (*(vNum2.end() -1))++;
                     else
                     {
                         vTag2.push_back(*itr);
                         vNum2.push_back(1);
                     }
                 }
             }

             /*
             /// make filter for this raw tag vector
             /// 1: the tag number (avoid the ramdon tag)
             for(int i = vNum2.size() - 1; i >= 0; i--)
             {
                 if(vNum2.at(i) <= 3)  // must > 3
                 {
                     vNum2.erase(vNum2.begin() + i);
                     vTag2.erase(vTag2.begin() + i);
                 }
             }
             ///2: transfer "S+E" --> "S"
             if(vTag2.size() == 2)
             {
                 if(vTag2[0] == 'S' && vTag2[1] == 'E')
                 {
                     vTag2.erase(vTag2.end()-1);
                     vNum2.erase(vNum2.end()-1);
                 }
             }*/

             if(vTag1.size() == 1 && vTag2.size() == 1)
             {
                 if(vTag1[0] == 'E' && vTag2[0] == 'E')
                 {
                     stCurCandi.strTag = "EE";
                     cout << IntToStr(vNum1.at(0)) << ", " << IntToStr(vNum2.at(0)) << " --> "
                          << stCurCandi.iStartPos << ", " << stCurCandi.iEndPos
                          << "   ------------   " << IntToStr(abs(stPI2.ucExonIndex - stPI1.ucExonIndex)) << endl;
                     g_EE++;
                 }
                 else if(vTag1[0] == 'E' && vTag2[0] == 'S')
                 {
                     stCurCandi.strTag = "ES";
                     g_ES++;
                 }
                 else if(vTag1[0] == 'S' && vTag2[0] == 'S')
                 {
                     stCurCandi.strTag = "SS";
                     g_SS++;
                 }
                 else if(vTag1[0] == 'S' && vTag2[0] == 'E')
                 {
                     stCurCandi.strTag = "SE";
                     g_SE++;
                 }
                 else
                 {
                     cout << vTag1[0] << ", " << vTag2[0] << " --> "
                          << IntToStr(stCurExon1.GetLength()) << ", "
                          << IntToStr(stCurExon2.GetLength()) << endl;
                     g_else++;
                 }
             }
             else
             {

                 for(vector<char>::iterator itr = vTag1.begin(); itr != vTag1.end(); itr++)
                 {
                     cout << *itr << " ";
                 }
                 cout << " | ";
                 for(vector<char>::iterator itr = vTag2.begin(); itr != vTag2.end(); itr++)
                 {
                     cout << *itr << " ";
                 }
                 cout << " ---> ";
                 for(vector<int>::iterator itr = vNum1.begin(); itr != vNum1.end(); itr++)
                 {
                     cout << IntToStr(*itr) << ", ";
                 }
                 cout << " | ";
                 for(vector<int>::iterator itr = vNum2.begin(); itr != vNum2.end(); itr++)
                 {
                     cout << IntToStr(*itr) << ", ";
                 }
                 //cout the length of each exon
                 cout << " --> "
                      << IntToStr(stCurExon1.GetLength()) << ", "
                      << IntToStr(stCurExon2.GetLength()) << endl;
                 //cout << endl;
                 stCurCandi.strTag = "LG";
                 g_large1++;
             }
         }
    }
}

void ClsFindCandidate::RefineCandiate(int iMinSupportReads,
                                      vector<St_Candidate>& vSelfCircCandi,
                                      vector<St_Candidate>& vRegCircCandi, string strChromName)
{   
    //make the filter for self circ Candi and Ref Circ Candi
    ///For self circle
    if(!vSelfCircCandi.empty())
    {
        for(vector<St_Candidate>::iterator itr = vSelfCircCandi.end() - 1;
            itr >= vSelfCircCandi.begin(); itr--)
        {
            if(itr->iSupportNum < iMinSupportReads)
                vSelfCircCandi.erase(itr);
        }
    }

    ///For regular circle
    if(!vRegCircCandi.empty())
    {
        for(vector<St_Candidate>::iterator itr = vRegCircCandi.end() - 1;
            itr >= vRegCircCandi.begin(); itr--)
        {
            if(itr->iSupportNum < iMinSupportReads)
                vRegCircCandi.erase(itr);
        }
    }

    /*
    for(map<St_Candidate, int>::iterator itr = m_mpSelfCircCandi.begin(); itr != m_mpSelfCircCandi.end();)
    {
        if(itr->second < iMinSupportReads)
        {
            m_mpSelfCircCandi.erase(itr++);
        }
        else
            itr++;
    }

    for(map<St_Candidate, int>::iterator itr = m_mpRegCircCandi.begin(); itr != m_mpRegCircCandi.end();)
    {
        if(itr->second < iMinSupportReads)
        {
            m_mpRegCircCandi.erase(itr++);
        }
        else
            itr++;
    }*/

    ofstream ofs;
    string strName = "./Detection_Result/Candidate_" + strChromName + ".txt";
    ofs.open(strName.c_str());

    ofstream ofsBrief; // 这个是为了解析起来比较方便的文件
    strName = "./Detection_Result/Brief_" + strChromName + ".txt";
    ofsBrief.open(strName.c_str());

    //1: Statistic Result
    ofs << "=============Statistic Number=============" <<endl;
    ofs << "Self Circular RNA   : " << IntToStr(vSelfCircCandi.size()) << endl;
    ofs << "Regular Circular RNA: " << IntToStr(vRegCircCandi.size()) << endl << endl;

    //2: For the full set of Self Circular RNA
    ofs << "=============Self Circular RNA=============" <<endl;
    for(vector<St_Candidate>::iterator itr = vSelfCircCandi.begin(); itr != vSelfCircCandi.end(); itr++)
    {
        ofs << strChromName << "\t"
             << "<" << IntToStr(itr->iStartPos) << ", " << IntToStr(itr->iEndPos) << ">"
             << "\t" << IntToStr(itr->iSupportNum) << endl;
        ofsBrief << strChromName << " "
                 << IntToStr(itr->iStartPos) << " " << IntToStr(itr->iEndPos) << " "
                 << IntToStr(itr->iSupportNum) << " "
                 //<< itr->strTag << " "
                 << itr->GetStrand() << " "
                 << "S" << endl;
    }
    /*for(map<St_Candidate, int>::iterator itr = m_mpSelfCircCandi.begin(); itr != m_mpSelfCircCandi.end(); itr++)
    {
        ofs << IntToStr(itr->first.ucChromIndex) << "\t"
             << "<" << IntToStr(itr->first.iStartPos) << ", " << IntToStr(itr->first.iEndPos) << ">"
             << "\t" << IntToStr(itr->second) << endl;
        ofsBrief << IntToStr(itr->first.ucChromIndex) << " "
                 << IntToStr(itr->first.iStartPos) << " " << IntToStr(itr->first.iEndPos) << " "
                 << IntToStr(itr->second) << " "
                 << "S" << endl;
    }*/

    //3: For the full set of Regular Circular RNA
    ofs << endl << "=============Regular Circular RNA=============" <<endl;
    for(vector<St_Candidate>::iterator itr = vRegCircCandi.begin(); itr != vRegCircCandi.end(); itr++)
    {
        ofs << strChromName << "\t"
             << "<" << IntToStr(itr->iStartPos) << ", " << IntToStr(itr->iEndPos) << ">"
             << "\t" << IntToStr(itr->iSupportNum) << endl;

        ofsBrief << strChromName << " "
                 << IntToStr(itr->iStartPos) << " " << IntToStr(itr->iEndPos) << " "
                 << IntToStr(itr->iSupportNum) << " "
                 //<< itr->strTag << " "
                 << itr->GetStrand() << " "
                 << "R" << endl;
    }
    /*
    for(map<St_Candidate, int>::iterator itr = m_mpRegCircCandi.begin(); itr != m_mpRegCircCandi.end(); itr++)
    {
        ofs << IntToStr(itr->first.ucChromIndex) << "\t"
             << "<" << IntToStr(itr->first.iStartPos) << ", " << IntToStr(itr->first.iEndPos) << ">"
             << "\t" << IntToStr(itr->second) << endl;

        ofsBrief << IntToStr(itr->first.ucChromIndex) << " "
                 << IntToStr(itr->first.iStartPos) << " " << IntToStr(itr->first.iEndPos) << " "
                 << IntToStr(itr->second) << " "
                 << "R" << endl;
    }*/

    ofs.close();
    ofsBrief.close();
}
