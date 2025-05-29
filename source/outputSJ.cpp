#include "ReadAlignChunk.h"
#include "Parameters.h"
#include "OutSJ.h"
#include <limits.h>
#include "ErrorWarning.h"
#include <stdlib.h>

#define REDEMPTION_THRESHOLD 10 // defines the maximum difference in lengths between the overhang and the min allowable overhang length for which we consider redemption
#define MAPPABILITY_THRESHOLD 0.5 // SJs overhangs that fall within the redemption threshold have to have a mappability >= MAPPABILITY_THRESHOLD to be output

float computeFinalScore(uint16 ohSize, uint16 smallKmerSize, uint16 largeKmerSize, float smallMap, float largeMap) {
    // differences between the length of the overhang and the smaller and larger kmers
    uint32 distanceSmall = abs(ohSize - smallKmerSize);
    uint32 distanceLarge = abs(ohSize - largeKmerSize);

    // computing weights based on differences
    float smallWeight = 1 - (distanceSmall / (distanceSmall + distanceLarge));
    float largeWeight = 1 - (distanceLarge / (distanceSmall + distanceLarge));

    // computing final mappability scores
    float mappability = (smallWeight * smallMap) + (largeWeight * largeMap);
    return mappability;
}

bool isSJRedeemed(Junction& oneSJ, Parameters& P, Genome& G) {
    if (G.pGe.mappabilityFileSmall == "" || G.pGe.mappabilityFileLarge == "") return false;

    // start positions of the right and left overhang
    uint32 startOverhangRight = *oneSJ.start + *oneSJ.gap;
    uint32 startOverhangLeft = *oneSJ.start - *oneSJ.overhangLeft;

    // chromosome names of the right and left overhang
    uint chrRightIdx = G.chrBin[startOverhangRight >> G.pGe.gChrBinNbits];
    uint chrLeftIdx = G.chrBin[startOverhangLeft >> G.pGe.gChrBinNbits];
    string chrRight = G.chrName.at(chrRightIdx);
    string chrLeft = G.chrName.at(chrLeftIdx);

    // chromosome starting positions of the right and left overhang
    uint32 ohChrRightStart = startOverhangRight + 1 - G.chrStart[chrRightIdx];
    uint32 ohChrLeftStart = startOverhangLeft + 1 - G.chrStart[chrLeftIdx];

    // mappability scores of the right and left overhang
    float rightOhScoreSmall = G.ratingsSmall.getScore(chrRight, ohChrRightStart);
    float rightOhScoreLarge = G.ratingsLarge.getScore(chrRight, ohChrRightStart);
    float leftOhScoreSmall = G.ratingsSmall.getScore(chrLeft, ohChrLeftStart);
    float leftOhScoreLarge = G.ratingsLarge.getScore(chrLeft, ohChrLeftStart);

    // cannot find the associated mapability score so just return false
    if (rightOhScoreLarge == -1 || rightOhScoreSmall == -1 || leftOhScoreLarge == - 1 || leftOhScoreSmall== -1) {
        ostringstream errOut;
        errOut <<"EXITING Could not find mapability score for chr: " << chrLeft << " " << ohChrRightStart << " and " << chrRight << " " << ohChrLeftStart << "\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    } 

    // computing the mapabilities
    float mappabilityLeft = computeFinalScore(*oneSJ.overhangLeft, G.ratingsSmall.kmerSize, G.ratingsLarge.kmerSize, leftOhScoreSmall, leftOhScoreLarge);
    float mappabilityRight = computeFinalScore(*oneSJ.overhangLeft, G.ratingsSmall.kmerSize, G.ratingsLarge.kmerSize, leftOhScoreSmall, leftOhScoreLarge);

    // computing the difference in length between the overhang and the minimum allowed overhang
    int32 leftDiff = (uint) P.outSJfilterOverhangMin[(*oneSJ.motif+1)/2] - *oneSJ.overhangLeft;
    int32 rightDiff = (uint) P.outSJfilterOverhangMin[(*oneSJ.motif+1)/2] - *oneSJ.overhangRight;

    bool redemptionFlagLeft = (leftDiff > 0) && (leftDiff <= REDEMPTION_THRESHOLD);
    bool redemptionFlagRight = (rightDiff > 0) && (rightDiff <= REDEMPTION_THRESHOLD);

    bool isRedeemedLeft = redemptionFlagLeft && (mappabilityLeft >= MAPPABILITY_THRESHOLD);
    bool isRedeemedRight = redemptionFlagRight && (mappabilityRight >= MAPPABILITY_THRESHOLD);

    return isRedeemedLeft && isRedeemedRight;
}

int compareUint(const void* i1, const void* i2) {//compare uint arrays
    uint s1=*( (uint*)i1 );
    uint s2=*( (uint*)i2 );

    if (s1>s2) {
        return 1;
    } else if (s1<s2) {
        return -1;
    } else {
        return 0;
    };
};

void outputSJ(ReadAlignChunk** RAchunk, Parameters& P) {//collapses junctions from all therads/chunks; outputs junctions to file

    Junction oneSJ(RAchunk[0]->RA->genOut);
    char** sjChunks = new char* [P.runThreadN+1];
    #define OUTSJ_limitScale 2
    OutSJ allSJ (P.limitOutSJcollapsed*OUTSJ_limitScale, P, RAchunk[0]->RA->genOut);

    if (P.outFilterBySJoutStage!=1) {//chunkOutSJ
        for (int ic=0;ic<P.runThreadN;ic++) {//populate sjChunks with links to data
            sjChunks[ic]=RAchunk[ic]->chunkOutSJ->data;
            memset(sjChunks[ic]+RAchunk[ic]->chunkOutSJ->N*oneSJ.dataSize,255,oneSJ.dataSize);//mark the junction after last with big number
        };
    } else {//chunkOutSJ1
        for (int ic=0;ic<P.runThreadN;ic++) {//populate sjChunks with links to data
            sjChunks[ic]=RAchunk[ic]->chunkOutSJ1->data;
            memset(sjChunks[ic]+RAchunk[ic]->chunkOutSJ1->N*oneSJ.dataSize,255,oneSJ.dataSize);//mark the junction after last with big number
        };
    };

    while (true) {
        int icOut=-1;//chunk from which the junction is output
        for (int ic=0;ic<P.runThreadN;ic++) {//scan through all chunks, find the "smallest" junction
            if ( *(uint*)(sjChunks[ic]) < ULONG_MAX && (icOut==-1 || compareSJ((void*) sjChunks[ic], (void*) sjChunks[icOut]) < 0 ) ) {
                    icOut=ic;
                };
        };

        if (icOut<0) break; // no more junctions to output

        for (int ic=0; ic<P.runThreadN; ic++) {//scan through all chunks, find the junctions equal to icOut-junction
            if (ic!=icOut && compareSJ((void*) sjChunks[ic], (void*) sjChunks[icOut])==0) {
                oneSJ.collapseOneSJ(sjChunks[icOut],sjChunks[ic],P);//collapse ic-junction into icOut
                sjChunks[ic] += oneSJ.dataSize;//shift ic-chunk by one junction
            };
        };

        // write out the junction - saves some info about the junction
        oneSJ.junctionPointer(sjChunks[icOut], 0); //point to the icOut junction


        //filter the junction
        bool isRedeemed = isSJRedeemed(oneSJ, P, RAchunk[0]->mapGen); 
        *oneSJ.strand = isRedeemed ? '+' : '-'; // replacing the strand with the redemption flag because I cannot get my own output to work...
        // bool isRedeemed = false;
        // *oneSJ.isRedeemed = isRedeemed;
        bool sjFilter;
        sjFilter = *oneSJ.annot > 0 \
                || ( ( *oneSJ.countUnique >= (uint) P.outSJfilterCountUniqueMin[(*oneSJ.motif+1)/2] \
                    || (*oneSJ.countMultiple + *oneSJ.countUnique) >= (uint) P.outSJfilterCountTotalMin[(*oneSJ.motif+1)/2] )\
                    /*
                    If the overhang is too small we disregard the splice junction
                    However what we want is that if the overhang is smaller than the parameter however its mappability is high (1 or within a certain threshold) we output the splice junction
                    */
                && ((*oneSJ.overhangLeft >= (uint) P.outSJfilterOverhangMin[(*oneSJ.motif+1)/2] \
                && *oneSJ.overhangRight >= (uint) P.outSJfilterOverhangMin[(*oneSJ.motif+1)/2]) || isRedeemed) \
                && ( (*oneSJ.countMultiple + *oneSJ.countUnique) > P.outSJfilterIntronMaxVsReadN.size() || *oneSJ.gap <= (uint) P.outSJfilterIntronMaxVsReadN[*oneSJ.countMultiple+*oneSJ.countUnique-1]) );

        if (sjFilter) { //record the junction in all SJ
            memcpy(allSJ.data + allSJ.N * oneSJ.dataSize, sjChunks[icOut], oneSJ.dataSize);
            allSJ.N++;
            if (allSJ.N == allSJ.Nstore-1 ) {
                /*
                ostringstream errOut;
                errOut <<"EXITING because of fatal error: buffer size for SJ output is too small\n";
                errOut <<"Solution: increase input parameter --limitOutSJcollapsed\n";
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                */
                allSJ.dataSizeIncrease();
                P.inOut->logMain << "Increased the size of chunkOutSJ to " << allSJ.Nstore <<'\n';
            };
        };

        sjChunks[icOut] += oneSJ.dataSize;//shift icOut-chunk by one junction
    };

    bool* sjFilter=new bool[allSJ.N];
    if (P.outFilterBySJoutStage != 2) {
        //filter non-canonical junctions that are close to canonical
        uint* sjA = new uint [allSJ.N*3];
        for (uint ii=0;ii<allSJ.N;ii++) {//scan through all junctions, filter by the donor ditance to a nearest donor, fill acceptor array
            oneSJ.junctionPointer(allSJ.data,ii);

            sjFilter[ii]=false;
            uint x1=0, x2=-1;
            if (ii>0)         x1=*( (uint*)(allSJ.data+(ii-1)*oneSJ.dataSize) ); //previous junction donor
            if (ii+1<allSJ.N) x2=*( (uint*)(allSJ.data+(ii+1)*oneSJ.dataSize) ); //next junction donor
            uint minDist=min(*oneSJ.start-x1, x2-*oneSJ.start);
            sjFilter[ii] = minDist >= (uint) P.outSJfilterDistToOtherSJmin[(*oneSJ.motif+1)/2];
            sjA[ii*3] = *oneSJ.start+(uint)*oneSJ.gap;//acceptor
            sjA[ii*3+1] = ii;

            if (*oneSJ.annot==0) {
                sjA[ii*3+2]=*oneSJ.motif;
            } else {
                sjA[ii*3+2]=SJ_MOTIF_SIZE+1;
            };

        };
        qsort((void*) sjA, allSJ.N, sizeof(uint)*3, compareUint);
        for (uint ii=0;ii<allSJ.N;ii++) {//
            if (sjA[ii*3+2]==SJ_MOTIF_SIZE+1) {//no filtering for annotated junctions
                sjFilter[sjA[ii*3+1]]=true;
            } else {
                uint x1=0, x2=-1;
                if (ii>0)         x1=sjA[ii*3-3]; //previous junction donor
                if (ii+1<allSJ.N) x2=sjA[ii*3+3]; //next junction donor
                uint minDist=min(sjA[ii*3]-x1, x2-sjA[ii*3]);
                sjFilter[sjA[ii*3+1]] = sjFilter[sjA[ii*3+1]] && ( minDist >= (uint) P.outSJfilterDistToOtherSJmin[(sjA[ii*3+2]+1)/2] );
            };
        };
    };

    //output junctions
    P.sjAll[0].reserve(allSJ.N);
    P.sjAll[1].reserve(allSJ.N);

    if (P.outFilterBySJoutStage != 1) {//output file
        ofstream outSJfileStream((P.outFileNamePrefix+"SJ.out.tab").c_str());
        ofstream outSJtmpStream((P.outFileTmp+"SJ.start_gap.tsv").c_str());
        for (uint ii=0; ii<allSJ.N; ii++) {//write to file
            if ( P.outFilterBySJoutStage == 2 || sjFilter[ii]  ) {
                oneSJ.junctionPointer(allSJ.data,ii);
                oneSJ.outputStream(outSJfileStream);//write to file
                outSJtmpStream << *oneSJ.start <<'\t'<< *oneSJ.gap <<'\n';
                P.sjAll[0].push_back(*oneSJ.start);
                P.sjAll[1].push_back(*oneSJ.gap);
            };
        };
        outSJfileStream.close();
    } else {//make sjNovel array in P
        P.sjNovelN=0;
        for (uint ii=0;ii<allSJ.N;ii++) { //count novel junctions
            if (sjFilter[ii]) { //only those passing filter
                oneSJ.junctionPointer(allSJ.data,ii);
                if (*oneSJ.annot==0) P.sjNovelN++;
            };
        };
        P.sjNovelStart = new uint [P.sjNovelN];
        P.sjNovelEnd = new uint [P.sjNovelN];
        P.inOut->logMain <<"Detected " << P.sjNovelN <<" novel junctions that passed filtering, will proceed to filter reads that contained unannotated junctions"<<endl;

        uint isj=0;
        for (uint ii=0;ii<allSJ.N;ii++) {//write to file
            if (sjFilter[ii]) {
                oneSJ.junctionPointer(allSJ.data,ii);
                if (*oneSJ.annot == 0) {//unnnotated only
                    P.sjNovelStart[isj]=*oneSJ.start;
                    P.sjNovelEnd[isj]=*oneSJ.start+(uint)(*oneSJ.gap)-1;
                    isj++;
                };
            };
        };
    };
};
