#ifndef H_Genome
#define H_Genome

#include "IncludeDefine.h"
#include "Parameters.h"
#include "PackedArray.h"
#include "SharedMemory.h"
#include "Variation.h"
#include "SuperTranscriptome.h"

class GTF;

typedef struct mappability {
    vector<string> chrNames;
    vector<uint64> start;
    vector<uint64> end;
    vector<float> ratings;
    int32 kmerSize;

    uint findIntervalIdx(uint startingPos) {
        for (uint i = 0; i < start.size(); i++) {
            if (startingPos >= start[i] && startingPos < end[i]) return i;
        }

        return -1;
    }
} mapRatings;

class Genome {
private:
    key_t shmKey;
    char *shmStart;
    uint OpenStream(string name, ifstream & stream, uint size);
    void HandleSharedMemoryException(const SharedMemoryException & exc, uint64 shmSize);
    void mapInfoLoad(string mapFilePath, mapRatings& mapInfo);
    void processMapLine(ifstream& mapStream, string mapFilePath, char* mapCharLine, string& line, mapRatings& mapInfo);
    void parseMapLine(string& line, mapRatings& mapInfo);
public:
    Parameters &P;
    ParametersGenome &pGe;
    SharedMemory *sharedMemory;

    enum {exT,exS,exE,exG,exL}; //indexes in the exonLoci array from GTF
     
    char *G, *G1;
    char *GMap; // contains the corresponding mappability rating for the base of the genome
    uint64 nGenome, nG1alloc;
    uint64 nGenomeMap; // size of GMap
    PackedArray SA,SAinsert,SApass1,SApass2;
    PackedArray SAi;
    Variation *Var;

    uint nGenomeInsert, nGenomePass1, nGenomePass2, nSAinsert, nSApass1, nSApass2;

    //chr parameters
    vector <uint64> chrStart, chrLength, chrLengthAll; // keeps track of chromosome start positions and lengths
    uint genomeChrBinNbases, chrBinN, *chrBin;
    vector <string> chrName, chrNameAll;
    map <string,uint64> chrNameIndex;

    uint *genomeSAindexStart;//starts of the L-mer indices in the SAindex, 1<=L<=pGe.gSAindexNbases

    uint nSA, nSAbyte, nChrReal;//genome length, SA length, # of chromosomes, vector of chromosome start loci
    uint nGenome2, nSA2, nSAbyte2, nChrReal2; //same for the 2nd pass
    uint nSAi; //size of the SAindex
    unsigned char GstrandBit, SAiMarkNbit, SAiMarkAbsentBit; //SA index bit for strand information
    uint GstrandMask, SAiMarkAbsentMask, SAiMarkAbsentMaskC, SAiMarkNmask, SAiMarkNmaskC;//maske to remove strand bit from SA index, to remove mark from SAi index

    //SJ database parameters
    uint sjdbOverhang, sjdbLength; //length of the donor/acceptor, length of the sj "chromosome" =2*pGe.sjdbOverhang+1 including spacer
    uint sjChrStart,sjdbN; //first sj-db chr
    uint sjGstart; //start of the sj-db genome sequence
    uint *sjDstart,*sjAstart,*sjStr, *sjdbStart, *sjdbEnd; //sjdb loci
    uint8 *sjdbMotif; //motifs of annotated junctions
    uint8 *sjdbShiftLeft, *sjdbShiftRight; //shifts of junctions
    uint8 *sjdbStrand; //junctions strand, not used yet

   //sequence insert parameters
    uint genomeInsertL; //total length of the sequence to be inserted on the fly
    uint genomeInsertChrIndFirst; //index of the first inserted chromosome

    // mappability rating parameters
    // we invlude mappability scores for large and small kmers in order to compute the in-between interpolated scores
    mapRatings ratingsSmall {};
    mapRatings ratingsLarge {};
    uint64 nMapRatings;

    //SuperTranscriptome genome
    SuperTranscriptome *superTr;

    Genome (Parameters &P, ParametersGenome &pGe);
    //~Genome();

    void freeMemory();
    void genomeLoad();
    void genomeOutLoad();
    void chrBinFill();
    void chrInfoLoad();
    void loadAllMappability();
    void genomeSequenceAllocate(uint64 nGenomeIn, uint64 &nG1allocOut, char*& Gout, char*& G1out);
    void loadSJDB(string &genDir);

    void insertSequences();

    //void consensusSequence(); DEPRECATED
    
    void genomeGenerate();
    void writeChrInfo(const string dirOut);
    void concatenateChromosomes(const vector<vector<uint8>> &vecSeq, const vector<string> &vecName, const uint64 padBin);
    void writeGenomeSequence(const string dirOut);
    
    //transform genome coordinates
    struct {
        bool convYes;
        bool gapsAreJunctions;
        Genome *g;
        string convFile;
        vector<array<uint64,3>> convBlocks;
        uint64 nMinusStrandOffset;//offset for the (-) strand, typically=nGenomeReal
    } genomeOut;
    
    typedef struct {
        uint64 pos;
        int32 len;//0: SNV, <0: deletion; >0: insertion
        array<string,2> seq;//sequence for SNV and insertions, empty for deletions
    } VariantInfo;
    
    void transformGenome(GTF *gtf) ;
    void transformChrLenStart(map<string,vector<VariantInfo>> &vcfVariants, vector<uint64> &chrStart1, vector<uint64> &chrLength1);
    void transformGandBlocks(map<string,vector<VariantInfo>> &vcfVariants, vector<uint64> &chrStart1, vector<uint64> &chrLength1, vector<array<uint64,3>> &transformBlocks, char *Gnew);
    void transformBlocksWrite(vector<array<uint64,3>> &transformBlocks);
    void transformExonLoci(vector<array<uint64,exL>> &exonLoci, vector<array<uint64,3>> &transformBlocks);
};
#endif
