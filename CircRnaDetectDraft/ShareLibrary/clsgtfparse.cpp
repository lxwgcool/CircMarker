#include "clsgtfparse.h"
//#include <stdio.h>
#include "unistd.h" // for "access"
//#include <iostream>
//#include <fstream>
//#include <string.h>
//#include <stdlib.h>
#include "clsbasealgorithm.h"
#include "clsfastareader.h"


//#define DEBUG
//#define PRINGTFSEQ

ClsGTFParse::ClsGTFParse(): m_strGtfPath(""), m_strDNARefPath(""), m_iKmerLen(-1), m_iReadLen(-1),
                            m_fKmerRatio(-1)
{
}

void ClsGTFParse::Init(string strGtfPath, string strDNARefPath, int iKmerLen,
                       int iReadLen, float fKmerRatio)
{
    m_strGtfPath = strGtfPath;
    m_strDNARefPath = strDNARefPath;
    m_iKmerLen = iKmerLen;
    m_iReadLen = iReadLen;
    m_fKmerRatio = fKmerRatio;
}

bool ClsGTFParse::ReadGTF(vector<St_Row_Chrom>& vChrom, string strGtfPath)
{
    if(strGtfPath != "")
    {
        m_strGtfPath = strGtfPath;
    }

    if(::access(m_strGtfPath.c_str(), 0) != 0)
    {
        cout << "Error: GTF File fail to be detected!" << endl;
        return false;
    }

    vChrom.clear();
    //Parse the contig file
    ifstream infile;
    infile.open(m_strGtfPath.c_str(), ios::in);
    string strLine = "";
    string strCurChromName = "";

    St_Row_Chrom stRChrom;
    St_Raw_Gene stRGene;
    St_Raw_Transcript stRT;
    St_Raw_Exon stRExon;

    while(!infile.eof()) // check if reached the end of file
    {
        getline(infile, strLine);
        string::size_type iPos = strLine.find('\t');
        if(iPos == string::npos)
            continue;

        int iStart = 0;
        //1: chromoson name
        string strChromoson = strLine.substr(iStart, iPos-iStart);
        ///Check if need to be saved to vChrom
        if(strCurChromName == "")
        {
            strCurChromName = strChromoson;
            stRChrom.strName = strChromoson;
        }
        else
        {
            if(strCurChromName == strChromoson)
            {}
            else //New one if "!="
            {
                vChrom.push_back(stRChrom);
                stRChrom.Refresh();
                stRGene.Refresh();
                stRT.Refresh();
                stRExon.Refresh();

                stRChrom.strName = strChromoson;
                strCurChromName = strChromoson;
            }
        }

        //2: null
        iPos++;
        iPos = strLine.find('\t', iPos);

        //3: Type
        iStart = iPos + 1;
        iPos = strLine.find('\t', iStart);
        string strType = strLine.substr(iStart, iPos-iStart);

        //4: Start Pos
        iStart = iPos + 1;
        iPos = strLine.find('\t', iStart);
        int iDataStart = atoi(strLine.substr(iStart, iPos-iStart).c_str());

        //5: End Pos
        iStart = iPos + 1;
        iPos = strLine.find('\t', iStart);
        int iDataEnd = atoi(strLine.substr(iStart, iPos-iStart).c_str());

        //6: Null
        iPos++;
        iPos = strLine.find('\t', iPos);

        //7: Direction
        iStart = iPos + 1;
        iPos = strLine.find('\t', iStart);
        string strDirection = strLine.substr(iStart, iPos-iStart);
        bool bRC = false;
        if(strDirection == "-")
            bRC = true;

        //8: NULL
        iPos++;
        iPos = strLine.find('\t', iPos);

        //10: Gene ID and Transcript ID, and Gene Name
        iStart = iPos + 1;
        string strValue = strLine.substr(iStart, strLine.length()-iStart);

        //10.1 Gene ID
        iPos = strValue.find("gene_id");
        iStart = strValue.find("\"", iPos) + 1;
        iPos = strValue.find("\"", iStart);
        string strGeneID = strValue.substr(iStart, iPos - iStart);

        //10.2 Transcript ID
        iPos = strValue.find("transcript_id");
        iStart = strValue.find("\"", iPos) + 1;
        iPos = strValue.find("\"", iStart);
        string strTranscriptID = strValue.substr(iStart, iPos - iStart);

        //10.3 Gene Name
        iPos = strValue.find("gene_name");
        iStart = strValue.find("\"", iPos) + 1;
        iPos = strValue.find("\"", iStart);
        string strGeneName = strValue.substr(iStart, iPos - iStart);

        //10.4 Gene Biotype --> this is important to distinguish the coding gene and non-coding gene
        iPos = strValue.find("gene_biotype");
        iStart = strValue.find("\"", iPos) + 1;
        iPos = strValue.find("\"", iStart);
        string strGeneBioType = strValue.substr(iStart, iPos - iStart);

        if(strType == "gene")
        {
            if(!stRT.vRExon.empty()) // save the last transcript
            {
                stRGene.vRT.push_back(stRT);
                stRT.Refresh();
            }

            if(stRGene.vRT.empty()) // it is the first one
            {}
            else
            {
                stRChrom.vRG.push_back(stRGene);
                stRGene.Refresh();
            }

            stRGene.iStart = iDataStart;
            stRGene.iEnd = iDataEnd;
            stRGene.bRC = bRC;
            stRGene.strChromoson = strChromoson;
            stRGene.strID = strGeneID;
            stRGene.strBioType = strGeneBioType;
            stRGene.strName = strGeneName;
        }
        else if(strType == "transcript")
        {
            if(stRT.vRExon.empty())
            {}
            else
            {
                stRGene.vRT.push_back(stRT);
                stRT.Refresh();
            }

            stRT.bRC = bRC;
            stRT.iStart = iDataStart;
            stRT.iEnd = iDataEnd;
            stRT.strID = strTranscriptID;
        }
        else if(strType == "exon")
        {
            stRExon.bRC = bRC;
            stRExon.iStart = iDataStart;
            stRExon.iEnd = iDataEnd;
            stRT.vRExon.push_back(stRExon);
            stRExon.Refresh();
        }
    }

    if(!stRT.vRExon.empty())
    {
        stRGene.vRT.push_back(stRT);
        stRT.Refresh();
    }

    if(!stRGene.vRT.empty())
    {
        stRChrom.vRG.push_back(stRGene);

        vChrom.push_back(stRChrom);

        stRGene.Refresh();
    }

    ///--> Try to Fill the missing part
    if(!vChrom.empty())
    {
        int iStartChrom = atoi(vChrom.begin()->strName.c_str()) - 1;
        St_Row_Chrom stRChrom;
        for(int i=0 ; i<iStartChrom; i++)
        {
            vChrom.insert(vChrom.begin(), stRChrom);
        }
    }

    ///---------> we try to erase the gene with "-" direction           
    /*
    for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin(); itrChrom != vChrom.end(); itrChrom++)
    {
        if(itrChrom->vRG.empty())
            continue;

        for(vector<St_Raw_Gene>::iterator itrGene = itrChrom->vRG.end()-1; itrGene >= itrChrom->vRG.begin();
            itrGene--)
        {
            if(!itrGene->bRC) //我们现在只保存负向序列
                itrChrom->vRG.erase(itrGene);
        }
    }
    */
    ///<---------


//    ofstream ofs;
//    ofs.open("./gtf_brif.txt");

//    for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin(); itrChrom != vChrom.end(); itrChrom++)
//    {
//        ofs << "Chromoson Name: " << itrChrom->strName << endl;
//        for(vector<St_Raw_Gene>::iterator itrGene = itrChrom->vRG.begin(); itrGene != itrChrom->vRG.end(); itrGene++)
//        {
//            ofs << "\t" << "Gene ID: " << itrGene->strID << " | " << "Name: " << itrGene->strName << " " << (itrGene->bRC ? "-" : "+") << endl;

//            for(vector<St_Raw_Transcript>::iterator itrRT = itrGene->vRT.begin(); itrRT != itrGene->vRT.end(); itrRT++)
//            {
//                ofs << "\t\t" << "Transcript ID: " << itrRT->strID << " --- " << IntToStr(itrRT->iStart)
//                    << " | " << IntToStr(itrRT->iEnd) << " " << (itrRT->bRC ? "-" : "+") << endl;

//                int iIndex = 0;
//                for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin(); itrExon != itrRT->vRExon.end(); itrExon++)
//                {
//                    ofs << "\t\t\t" << "Exon" << IntToStr(iIndex) << ": "<< IntToStr(itrExon->iStart)
//                        << " | " << IntToStr(itrExon->iEnd) << " " << (itrExon->bRC ? "-" : "+") << endl;
//                    iIndex++;
//                }
//            }
//        }
//    }

//    //Cout the statistic result
//    for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin(); itrChrom != vChrom.end(); itrChrom++)
//    {
//         ofs << "Chromoson Name: " << itrChrom->strName << endl;
//         ofs << "\t" << "Gene Num: " << IntToStr(itrChrom->vRG.size()) << endl;
//         for(vector<St_Raw_Gene>::iterator itrGene = itrChrom->vRG.begin(); itrGene != itrChrom->vRG.end(); itrGene++)
//         {
//            ofs << IntToStr(itrGene->vRT.size()) << "  ";
//         }
//         ofs << endl << endl;
//    }

//    ofs.close();
    return true;
}

//-->The format of this "RefRNA" is identified by us!!!
void ClsGTFParse::GetRNARef(string strRNARefPath, string strChromName,
                            vector<St_Raw_Gene>& vGenes, bool bExportSeq)
{
    //Read DNA Ref
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    vector<St_Fasta> vFasta;
    pFastaReader->ReadFastaRegular(m_strDNARefPath, vFasta);
    delete pFastaReader;
    pFastaReader = NULL;

    ofstream ofs;
    ofs.open(strRNARefPath.c_str());

    St_Fasta* pFasta = NULL;

    for(vector<St_Fasta>::iterator itrRef = vFasta.begin(); itrRef != vFasta.end(); itrRef++)
    {
        //get the chromson ID:
        string strRefID = itrRef->strName.substr(0, itrRef->strName.find(" "));
        if(itrRef->strName.find(" ") == string::npos)
            strRefID = itrRef->strName;

        if(strRefID == strChromName)
        {
            pFasta = &(*itrRef);
            break;
        }
    }

    if(pFasta == NULL)
        return;

    for(vector<St_Raw_Gene>::iterator itr = vGenes.begin(); itr != vGenes.end(); itr++)
    {        

        for(vector<St_Raw_Transcript>::iterator itrRT = itr->vRT.begin(); itrRT != itr->vRT.end(); itrRT++)
        {
            string strSeq = "";
            string strName = "";
            string strIndex = "";
            string strExonBundary = "";
            int iExonIndex = 0;
            for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin(); itrExon != itrRT->vRExon.end(); itrExon++)
            {
                int iStart = itrExon->iStart < itrExon->iEnd ? itrExon->iStart : itrExon->iEnd;
                int iEnd = itrExon->iStart > itrExon->iEnd ? itrExon->iStart : itrExon->iEnd;
                int iLen = iEnd - iStart + 1;

                strSeq += pFasta->strSeq.substr(iStart, iLen);

                strIndex += IntToStr(iExonIndex) + ", ";
                strExonBundary += "(" + IntToStr(itrExon->iStart) + ", " + IntToStr(itrExon->iEnd) + "), ";
                iExonIndex++;
            }
            strName = itr->strChromoson + " | " + itr->strName + " | " + IntToStr(itr->iStart) + ":" + IntToStr(itr->iEnd) + " | " +
                      itrRT->strID + " | " + IntToStr(itrRT->iStart) + ":" + IntToStr(itrRT->iEnd) +
                      " | " + strIndex + " | " + strExonBundary + " | " + IntToStr(strSeq.length());
            ofs << ">" << strName << endl;
            if(bExportSeq)
                ofs << strSeq << endl;
            else
                ofs << "" << endl;
        }

        pFasta = NULL;
    }

    ofs.close();
    vFasta.clear();
}

// Default: do NOT load sequence info
void ClsGTFParse::LoadRNARef(string strRNARefPath, vector<St_Raw_Gene>& vGenes, bool bLoadSeq)
{
    vGenes.clear();

    //Load RNA Ref File
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    vector<St_Fasta> vFasta;
    pFastaReader->ReadFastaRegular(strRNARefPath, vFasta, bLoadSeq);
    delete pFastaReader;
    pFastaReader = NULL;

    St_Raw_Gene stRG; // raw Gene   ----------------> Need to change the code to fit the data stucture !!!!!
    St_Raw_Transcript stRT; //raw Transcript
    St_Raw_Exon stRE; //raw exon
    for(vector<St_Fasta>::iterator itr = vFasta.begin(); itr != vFasta.end(); itr++)
    {
        int iStart = 0;
        int iEnd = itr->strName.find('|') - 1;

        //Chromosome Name
        string strChromoson = itr->strName.substr(iStart, iEnd - iStart);

        //Gene Name
        iStart = iEnd + 3;
        iEnd = itr->strName.find('|', iStart) - 1;
        string strGeneName = itr->strName.substr(iStart, iEnd - iStart);

        //Gene Start/End
        iStart = iEnd + 3;
        iEnd = itr->strName.find('|', iStart) - 1;
        string strTemp = itr->strName.substr(iStart, iEnd - iStart);
        int iGeneStart = atoi(strTemp.substr(0, strTemp.find(':')).c_str());
        int iGeneEnd = atoi(strTemp.substr(strTemp.find(':') + 1, strTemp.length() - strTemp.find(':') - 1).c_str());

        //---->Check if we need to start a new gene or we should save this transcript in current gene
        if(stRG.strChromoson != strChromoson || stRG.strName != strGeneName ||
           stRG.iStart != iGeneStart || stRG.iEnd != iGeneEnd)
        {
            if(stRG.strChromoson != "") // Do not the initial statement
            {
                vGenes.push_back(stRG);
                stRG.Refresh();
                //Set the value for stRG
                stRG.strChromoson = strChromoson;
                stRG.strName = strGeneName;
                stRG.iStart = iGeneStart;
                stRG.iEnd = iGeneEnd;
            }
            else
            {
                stRG.strChromoson = strChromoson;
                stRG.strName = strGeneName;
                stRG.iStart = iGeneStart;
                stRG.iEnd = iGeneEnd;
            }
        }

        //Transcription Name
        iStart = iEnd + 3;
        iEnd = itr->strName.find('|', iStart) - 1;
        stRT.strID = itr->strName.substr(iStart, iEnd - iStart);

        //Gene Start & End
        iStart = iEnd + 3;
        iEnd = itr->strName.find('|', iStart) - 1;
        strTemp = itr->strName.substr(iStart, iEnd - iStart);
        stRT.iStart = atoi(strTemp.substr(0, strTemp.find(':')).c_str());
        stRT.iEnd = atoi(strTemp.substr(strTemp.find(':') + 1, strTemp.length() - strTemp.find(':') - 1).c_str());

        //Exons
        iStart = iEnd + 3;
        iEnd = itr->strName.find('|', iStart) - 1;
        iStart = iEnd + 3;
        iEnd = itr->strName.find('|', iStart) - 1;
        strTemp = itr->strName.substr(iStart, iEnd - iStart);

        iStart = 0;

        while(strTemp.find('(', iStart) != string::npos)
        {
            //For exon start
            iStart = strTemp.find('(', iStart) + 1;
            iEnd = strTemp.find(',', iStart);
            stRE.iStart = atoi(strTemp.substr(iStart, iEnd - iStart).c_str());

            //For exon end
            iStart = iEnd + 2;
            iEnd = strTemp.find(')', iStart);
            stRE.iEnd = atoi(strTemp.substr(iStart, iEnd - iStart).c_str());

            stRT.vRExon.push_back(stRE);
        }

        //push transcript to current gene
        stRG.vRT.push_back(stRT);
        stRT.Refresh();

        if(itr == vFasta.end() - 1) // The last one
        {
            vGenes.push_back(stRG);
        }
    }

    vFasta.clear();
}
void ClsGTFParse::GetTagValue(vector<St_Row_Chrom>& vChrom) //This Tag means: the cicular splicing tag
{    
    //cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>GetTagValue<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    //Read DNN Reference
    //Read DNA Ref
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    vector<St_Fasta> vFasta;
    pFastaReader->ReadFastaRegular(m_strDNARefPath, vFasta);
    delete pFastaReader;
    pFastaReader = NULL;

#ifdef DEBUG
    bool bShowStartExonInfo = true;
    bool bShowEndExonInfo = true;
#endif

    //Get tag value for each exon
    for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin(); itrChrom != vChrom.end(); itrChrom++)
    {
        St_Fasta* pCurRefFa = NULL;
        for(vector<St_Fasta>::iterator itrRef = vFasta.begin(); itrRef != vFasta.end(); itrRef++)
        {
            //cout << itrRef->strName << endl;
            string strRefName = "";
            if(itrRef->strName.find(' ') == string::npos)
                strRefName = itrRef->strName;
            else
                strRefName = itrRef->strName.substr(0, itrRef->strName.find(' '));

            //cout << "strRefName: " << strRefName << endl;
            //cout << "itrChrom->strName" << itrChrom->strName << endl;

            if(strRefName == itrChrom->strName)
            {
                pCurRefFa = &(*itrRef);

//                ///For Debug -->
//                cout << "strRefName       : " << itrRef->strName << endl;
//                cout << "itrChrom->strName: " << itrChrom->strName << endl;
//                ///<---

                break;
            }
        }
        if(pCurRefFa == NULL)
        {
            //cout << "Do not find related reference chromosone: " << itrChrom->strName << endl;
            continue;
        }

        for(vector<St_Raw_Gene>::iterator itrGene = itrChrom->vRG.begin();
            itrGene != itrChrom->vRG.end(); itrGene++)
        {
            //Get corresponding reference sequence
            for(vector<St_Raw_Transcript>::iterator itrRT = itrGene->vRT.begin();
                itrRT != itrGene->vRT.end(); itrRT++)
            {
                for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin();
                    itrExon != itrRT->vRExon.end(); itrExon++)
                {
                    if(itrExon->iStart > itrExon->iEnd)
                    {
                        cout << "Exon Start > End, Strange!" << endl;
                    }
                    else
                    {
                        //cout << IntToStr(itrExon->iStart-1-2) << ": " << IntToStr(itrExon->iEnd-1+1) << endl;
                        //cout << IntToStr(pCurRefFa->strSeq.length()) << endl;

                        //这里我们需要存储真实值，也就是需要考虑RC值
                        /*
                        if(itrExon->bRC)
                        {
                            //Beginning tag
                            itrExon->strHead2Bp = GetReverseCompelement(pCurRefFa->strSeq.substr(itrExon->iEnd-1+1, 2));
                            ToUpper(itrExon->strHead2Bp);
                            //Ending tag
                            itrExon->strTail2Bp = GetReverseCompelement(pCurRefFa->strSeq.substr(itrExon->iStart-1-2, 2));
                            ToUpper(itrExon->strTail2Bp);
                        }
                        else*/

                        {
                            //Beginning tag
                            itrExon->strHead2Bp = pCurRefFa->strSeq.substr(itrExon->iStart-1-2, 2);
                            ToUpper(itrExon->strHead2Bp);
                            //Ending tag
                            itrExon->strTail2Bp = pCurRefFa->strSeq.substr(itrExon->iEnd-1+1, 2);
                            ToUpper(itrExon->strTail2Bp);
                        }

                        //cout << itrExon->strHead2Bp << endl;
                        //cout << itrExon->strTail2Bp << endl;
                        //cout << "Done--->" << endl;
#ifdef DEBUG
                        //Get Boundary Length
                        int iExtractLen = m_iReadLen * m_fKmerRatio + m_iKmerLen - 1;
                        if(itrExon->iStart == 1191425 && itrExon->iEnd == 1191505)
                        {
                            //看看本身的长度--> 将相应的形成kmer的那一段拿出来
                            if(itrExon->GetLength() <= 2 * m_iReadLen * m_fKmerRatio + m_iKmerLen*2 - 1)
                            {
                                //stPosInfo.cPart = 'U';
                                iExtractLen = (itrExon->GetLength() - m_iKmerLen)/2 + m_iKmerLen;
                            }

                            if(bShowEndExonInfo)
                            {                                
                                cout << "<" << IntToStr(itrExon->iStart) << ", " << IntToStr(itrExon->iEnd)
                                     << ">" << endl;
                                cout << "Head Tag: " << itrExon->strHead2Bp << endl;
                                cout << "Tail Tag: " << itrExon->strTail2Bp << endl;
                                if(itrExon->bRC)
                                    cout << "-" << endl;
                                else
                                    cout << "+" << endl;

                                if(itrExon->GetIsSupportCircRNA())
                                    cout << "Support Tag" << endl;
                                else
                                    cout << "Do NOT Support Tag" << endl;
                                cout << IntToStr(itrExon->GetLength()) << endl;

                                //cout the exon sequence
                                cout << pCurRefFa->strSeq.substr(itrExon->iStart-1-2,
                                                                 itrExon->iEnd - itrExon->iStart + 5) << endl;
                                cout << "Left Extract: "
                                     << pCurRefFa->strSeq.substr(itrExon->iStart - 1, iExtractLen) << endl;
                                     //<< pCurRefFa->strSeq.substr(itrExon->iStart-50, 50) << endl;
                                cout << "Right Extract: "
                                     << pCurRefFa->strSeq.substr(itrExon->iEnd - iExtractLen,
                                                                 iExtractLen) << endl;
                                     //<< pCurRefFa->strSeq.substr(itrExon->iEnd, 50) << endl;
                                cout << endl;
                                bShowEndExonInfo = false;
                            }
                        }

                        if(itrExon->iStart == 1203242 && itrExon->iEnd == 1203372)
                        {
                            //看看本身的长度--> 将相应的形成kmer的那一段拿出来
                            if(itrExon->GetLength() <= 2 * m_iReadLen * m_fKmerRatio + m_iKmerLen*2 - 1)
                            {
                                //stPosInfo.cPart = 'U';
                                iExtractLen = (itrExon->GetLength() - m_iKmerLen)/2 + m_iKmerLen;
                            }

                            if(bShowStartExonInfo)
                            {
                                cout << "<" << IntToStr(itrExon->iStart) << ", " << IntToStr(itrExon->iEnd)
                                     << ">" << endl;
                                cout << "Head Tag: " << itrExon->strHead2Bp << endl;
                                cout << "Tail Tag: " << itrExon->strTail2Bp << endl;
                                if(itrExon->bRC)
                                    cout << "-" << endl;
                                else
                                    cout << "+" << endl;

                                if(itrExon->GetIsSupportCircRNA())
                                    cout << "Support Tag" << endl;
                                else
                                    cout << "Do NOT Support Tag" << endl;
                                cout << IntToStr(itrExon->GetLength()) << endl;

                                //cout the exon sequence
                                cout << pCurRefFa->strSeq.substr(itrExon->iStart-1-2,
                                                                 itrExon->iEnd - itrExon->iStart + 5) << endl;
                                cout << "Left Extract: "
                                     << pCurRefFa->strSeq.substr(itrExon->iStart - 1, iExtractLen) << endl;
                                     //<< pCurRefFa->strSeq.substr(itrExon->iStart-50, 50) << endl;
                                cout << "Right Extract: "
                                     << pCurRefFa->strSeq.substr(itrExon->iEnd - iExtractLen,
                                                                 iExtractLen) << endl;
                                     //<< pCurRefFa->strSeq.substr(itrExon->iEnd, 50) << endl;
                                cout << endl;
                                bShowStartExonInfo = false;
                            }
                        }                        
#endif
                    }
                }
            }
        }
    }

#ifdef PRINGTFSEQ
    //Print Out the Collection Result
    cout << "Output gtf_brif.txt with Seq Info >>>>>>>>>>>>>>>>>>>" << endl;
    ofstream ofs;
    ofs.open("./gtf_brif.txt");

    int iChromIndex = 0;
    for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin(); itrChrom != vChrom.end();
        itrChrom++, iChromIndex++)
    {
        ///Find related reference file --->
        St_Fasta* pCurRefFa = NULL;
        for(vector<St_Fasta>::iterator itrRef = vFasta.begin(); itrRef != vFasta.end(); itrRef++)
        {
            cout << itrRef->strName << endl;
            string strRefName = "";
            if(itrRef->strName.find(' ') == string::npos)
                strRefName = itrRef->strName;
            else
                strRefName = itrRef->strName.substr(0, itrRef->strName.find(' '));

            cout << "strRefName: " << strRefName << endl;
            cout << "itrChrom->strName" << itrChrom->strName << endl;

            if(strRefName == itrChrom->strName)
            {
                pCurRefFa = &(*itrRef);
                break;
            }
        }
        if(pCurRefFa == NULL)
        {
            cout << "Do not find related reference chromosone: " << itrChrom->strName << endl;
            continue;
        }
        ///<----

        ofs << iChromIndex << " " << "Chromoson Name: " << itrChrom->strName << endl;
        int iGeneIndex = 0;
        for(vector<St_Raw_Gene>::iterator itrGene = itrChrom->vRG.begin(); itrGene != itrChrom->vRG.end();
            itrGene++, iGeneIndex++)
        {
            ofs << "\t" << iGeneIndex << " " <<"Gene ID: " << itrGene->strID << " | " << "Name: " << itrGene->strName << " " << (itrGene->bRC ? "-" : "+") << endl;

            int iTranscriptIndex= 0;
            for(vector<St_Raw_Transcript>::iterator itrRT = itrGene->vRT.begin(); itrRT != itrGene->vRT.end();
                itrRT++, iTranscriptIndex++)
            {
                ofs << "\t\t" << iTranscriptIndex << " " << "Transcript ID: " << itrRT->strID << " --- " << IntToStr(itrRT->iStart)
                    << " | " << IntToStr(itrRT->iEnd) << " " << (itrRT->bRC ? "-" : "+") << endl;

                int iIndex = 0;
                for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin(); itrExon != itrRT->vRExon.end(); itrExon++)
                {
                    ofs << "\t\t\t" << "Exon " << IntToStr(iIndex) << ": "<< itrExon->strHead2Bp
                        << " | " << itrExon->strTail2Bp << " " << (itrExon->bRC ? "-" : "+")
                        << "   <" << IntToStr(itrExon->iStart) << ", " << IntToStr(itrExon->iEnd) << ">"
                        << " Len: " << IntToStr(abs(itrExon->iEnd - itrExon->iStart)) << endl;

                    //-->Display Exon Sequence
                    ofs << "Seq:" << endl;

                    string strExonSeq = "Nil";
                    if(pCurRefFa != NULL)
                    {
                        strExonSeq = pCurRefFa->strSeq.substr(itrExon->iStart-1, itrExon->GetLength());
                    }
                    DisplayString(ofs, strExonSeq);
                    ofs << endl;
                    //<--

                    iIndex++;
                }
            }
        }
    }

    ofs.close();
#endif

    //release reference sequence
    vFasta.clear();
}

//void ClsGTFParse::ColletcPossibleCRNA(vector<St_Raw_Gene>& vGenes)
//{
//    /*
//    vector<St_Raw_Exon> vHeadExon;
//    vector<int> vHeadIndex;
//    vector<St_Raw_Exon> vTailExon;
//    vector<int> vTailIndex;
//    St_CandiAtom stAtom;
//    int iIndex = 0;
//    for(vector<St_Raw_Gene>::iterator itr = vGenes.begin(); itr != vGenes.end(); itr++)
//    {
//        for(vector<St_Raw_Transcript>::iterator itrRT = itr->vRT.begin(); itrRT != itr->vRT.end(); itrRT++)
//        {
//            vHeadExon.clear();
//            vHeadIndex.clear();
//            vTailExon.clear();
//            vTailIndex.clear();
//            iIndex = 0;
//            for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin(); itrExon != itrRT->vRExon.end(); itrExon++)
//            {
//                //For Head Exon
//                if(itrExon->GetIsCRNAHeadExon())
//                {
//                    vHeadExon.push_back(*itrExon);
//                    vHeadIndex.push_back(iIndex);
//                }
//                //For Tail Exon
//                if(itrExon->GetIsCRNATailExon())
//                {
//                    vTailExon.push_back(*itrExon);
//                    vTailIndex.push_back(iIndex);
//                }
//                iIndex++;
//            }
//            //Create Possible Transcription Exon Pair
//            if(vHeadExon.empty() || vTailExon.empty())
//                continue;  //Next Transcript
//            //begin detection possible CRNA Exon Pair --->
//            for(int i=0; i<vHeadExon.size(); i++)
//            {
//                for(int j=0; j<vTailExon.size(); j++)
//                {
//                    if(vHeadIndex[i] <= vTailIndex[j]) // 头在尾巴的前面(重合也是允许的)
//                    {
//                        int iStart = vHeadExon[i].iStart;
//                        int iEnd = vTailExon[j].iEnd;
//                        if(iStart <= iEnd)
//                        {
//                            stAtom.iStart =iStart;
//                            stAtom.iEnd =iEnd;
//                        }
//                        else //Switch the value
//                        {
//                            stAtom.iStart =iEnd;
//                            stAtom.iEnd =iStart;
//                        }
//                        itrRT->vPossibleCRNA.push_back(stAtom);
//                    }
//                }
//            }
//            //<--
//        }
//    }*/
//}

