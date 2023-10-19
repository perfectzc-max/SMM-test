import sys
import datetime

from params import Params
from utils import Utils

class Read(object):

    def __init__(self, aFam, read, sUMI, isPositive, isOverlapped=False):
        '''
        Constructor
        '''
        self.fam = aFam  # My Family
        self.origRead = read
        self.sUMI = sUMI
        self.isPositive = isPositive
        self.isOverlapped = isOverlapped  # True if is second in the pair

        self.myBrother = None  # Overlapped read index

        if isPositive:
            self.sUMIMod = sUMI
        else:
            self.sUMIMod = Utils.umiSwap(sUMI)

        self.baseList = None

    def sameOrientation(self, read):
        res = (self.origRead.is_reverse == read.is_reverse) and (self.origRead.is_read1 == read.is_read1)
        res = res or (self.origRead.is_reverse != read.is_reverse) and (self.origRead.is_read1 != read.is_read1)
        return res

    def calcBaseList(self):
        resList = []  # List of string
        read = self.origRead

        nReadPos = self.fam.nStart - read.reference_start
        if nReadPos < 0:
            resList = [None] * (-nReadPos)
            nReadPos = 0

        flagM = 0  # Indicates that M section was not processed yet
        for i in range(len(read.cigartuples)):
            nSect, nCount = read.cigartuples[i][0], read.cigartuples[i][1]

            if nSect == Params.SectM:  # Put bases (M section)
                for j in range(nReadPos, nReadPos + nCount):
                    if read.query_qualities[j] < Params.minBaseQuality:
                        resList += Params.chNoQual
                    else:
                        resList += read.query_sequence[j]
                nReadPos += nCount
                flagM = 1
            elif nSect == Params.SectI:  # Put inserted bases (I section)
                if flagM:  # if I is after M then
                    resList[-1] += read.query_sequence[nReadPos: nReadPos + nCount]  # Add insertion to last base
                # Otherwise treat it as soft-clipped
                nReadPos += nCount

            elif nSect == Params.SectD:  # Put 'D' if deleted (D section)
                resList += [Params.chDel] * nCount

            elif nSect == Params.SectS:  # Skip soft-clipped part (S section)
                nReadPos += nCount

            elif nSect == Params.SectH:  # Ignore hard-clipped part (H section)
                pass
            else:
                sys.exit('Unrecognized section code in Cigar tuple:\t%s' % nSect)

        nRefLen = self.fam.nEnd - self.fam.nStart  # 0-based, nEnd - behind the
        nMyLen = len(resList)

        # Set length
        if nMyLen < nRefLen:
            resList += [None] * (nRefLen - nMyLen)
        else:
            del resList[nRefLen:]

        self.baseList = resList

    def printOut(self):
        read = self.origRead
        myStr = "{}\t{}\t{}".format(read.query_name, read.reference_start, read.reference_end)
        print(myStr)

    def logOut(self):
        read = self.origRead
        myStr = "{}\t{}\t{}".format(read.query_name, read.reference_start, read.reference_end)
        Utils.logOut(myStr, True)

    def getBaseAt(self, nRefStart):
        sBase = None
        if nRefStart in range(self.fam.nStart, self.fam.nEnd):
            if not self.baseList:
                self.calcBaseList()
            sBase = self.baseList[nRefStart - self.fam.nStart]
        return sBase

    def backTrackPos(self, btResData, nRefStart):
        chReadBase = self.getBaseAt(nRefStart)
        rd = self.origRead
        # Read info string
        outString = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chReadBase, int(self.isPositive), rd.query_name, \
                                                             rd.reference_start, rd.reference_end, rd.next_reference_start,
                                                             rd.cigarstring, rd.mapping_quality)
        btResData.write(outString)


class Family(object):

    def __init__(self, aChrom, sUMI, nStart, nEnd):
        '''
        Constructor
        '''
        self.currChrom = aChrom
        self.sUMI = sUMI
        self.nStart = nStart  # Left position of the family (0-based)
        self.nEnd = nEnd  # Right position of the family (0-based)
        self.readList = []  # Set of Read objects

    def addRead(self, read):
        self.readList.append(read)

    def inFamily(self, read):
        return self.nStart <= read.reference_start and self.nEnd >= read.reference_end

    def printOut(self):
        myStr = "Chrom: {}, UMI: {}, Start: {}, End: {}, Num reads: {}".format(self.currChrom, self.sUMI, self.nStart,
                                                                              self.nEnd, len(self.readList))
        print(myStr)

    def logOut(self):
        myStr = "Chrom: {}, UMI: {}, Start: {}, End: {}, Num reads: {}".format(self.currChrom, self.sUMI, self.nStart,
                                                                              self.nEnd, len(self.readList))
        Utils.logOut(myStr, True)


class Group(object):
    # Overlap part calculation in progress
    def __init__(self, aChrom, sUMI, nStart, nEnd):
        self.familyList = []  # List of Families
        self.nOverlapStart = nStart  # (0-based)
        self.nOverlapEnd = nEnd  # (0-based)
        self.chrom = aChrom
        self.sUMI = sUMI
        self.readList = []  # Set of Read objects

    def addRead(self, read):
        self.readList.append(read)

    def printOut(self):
        for fam in self.familyList:
            fam.printOut()

    def logOut(self):
        for fam in self.familyList:
            fam.logOut()

    def addFamily(self, fam):
        self.familyList.append(fam)
        if len(self.familyList) == 1:
            self.nOverlapStart = fam.nStart
            self.nOverlapEnd = fam.nEnd
        else:
            if fam.nStart < self.nOverlapStart:
                self.nOverlapStart = fam.nStart
            if fam.nEnd > self.nOverlapEnd:
                self.nOverlapEnd = fam.nEnd

    def calcOverlapped(self, refSeq):
        for fam in self.familyList:
            for read in fam.readList:
                read.calcBaseList()

        bCorrect = False  # Found good base
        for nRefPos in range(self.nOverlapStart, self.nOverlapEnd):
            chGoodBase = refSeq[nRefPos]  # Expected
            chCurBase = None
            chBestBase = None
            nBestCount = 0
            bOnlyBad = True  # All bases are bad
            bSkip = False  # Skip position

            for fam in self.familyList:
                for read in fam.readList:
                    chReadBase = read.getBaseAt(nRefPos)
                    if chReadBase and chReadBase != Params.chDel:
                        bOnlyBad = False
                        if chCurBase is None:  # Current base is not set yet
                            chCurBase = chReadBase
                            nBestCount = 1
                            chBestBase = chCurBase
                        elif chReadBase == chCurBase:
                            nBestCount += 1
                        else:
                            if nBestCount == 1:
                                break
                            else:
                                bSkip = True
                                break
                if bSkip:
                    break

            if not bOnlyBad and chBestBase and nBestCount == 1:
                bCorrect = True
                refSeq[nRefPos] = chBestBase  # Set corrected base
                for fam in self.familyList:
                    for read in fam.readList:
                        chReadBase = read.getBaseAt(nRefPos)
                        if chReadBase == chBestBase:
                            read.backTrackPos(btResData, nRefPos)
                            chReadBase = None
                        elif chReadBase and chReadBase != Params.chDel:
                            read.backTrackPos(btErrData, nRefPos)
                            chReadBase = None
        if not bCorrect:
            return False

        return True


class ClusterData(object):
    # One Chromosome + UMI
    def __init__(self, aChrom, sUMI):
        '''
        Constructor
        '''
        self.myGroups = []  # List of Group objects
        self.currChrom = aChrom
        self.sUMI = sUMI
        self.familyList = []  # List of Families

    def addRead(self, read):
        sUMI = read.sUMIMod
        inGroup = False
        inFamily = False
        inFamilyObject = None

        # 1. Check for the same group
        for gr in self.myGroups:
            if gr.sUMI == sUMI:
                inGroup = True
                if read.sameOrientation(gr.readList[0]):  # If the same orientation as in group
                    gr.addRead(read)
                else:
                    bAdded = False
                    for fam in gr.familyList:
                        if fam.inFamily(read):
                            inFamily = True
                            inFamilyObject = fam
                            inFamilyObject.addRead(read)
                            bAdded = True
                            break
                    if not inFamily:
                        newFam = Family(self.currChrom, sUMI, read.reference_start, read.reference_end)
                        newFam.addRead(read)
                        gr.addFamily(newFam)
                        self.familyList.append(newFam)

        # 2. Otherwise add new group
        if not inGroup:
            newGr = Group(self.currChrom, sUMI, read.reference_start, read.reference_end)
            newGr.addRead(read)
            self.myGroups.append(newGr)
            newFam = Family(self.currChrom, sUMI, read.reference_start, read.reference_end)
            newFam.addRead(read)
            newGr.addFamily(newFam)
            self.familyList.append(newFam)

    def calcOverlapped(self, refSeq):
        for gr in self.myGroups:
            gr.calcOverlapped(refSeq)

    def printOut(self):
        for fam in self.familyList:
            fam.printOut()

    def logOut(self):
        for fam in self.familyList:
            fam.logOut()

    def createResultStrings(self):
        resultStrings = []
        for gr in self.myGroups:
            for fam in gr.familyList:
                # Add header
                headerString = "{}\t{}\t{}\t{}".format(gr.chrom, fam.sUMI, fam.nStart, fam.nEnd)
                resultStrings.append(headerString)

                # Prepare corrected sequence
                seqList = [Params.chDel] * (fam.nEnd - fam.nStart)
                for read in fam.readList:
                    for i in range(fam.nStart, fam.nEnd):
                        chReadBase = read.getBaseAt(i)
                        if chReadBase != Params.chDel:
                            seqList[i - fam.nStart] = chReadBase
                seqString = "".join(seqList)
                resultStrings.append(seqString)

        return resultStrings


class CorrectionData(object):
    '''
    classdocs
    '''

    def __init__(self, refSeq, BTLength):
        '''
        Constructor
        '''
        self.refSeq = refSeq
        self.BTLength = BTLength
        self.clusterData = {}
        self.bFirstUMI = True
        self.nCorrCount = 0

    def addRead(self, read):
        if read is not None and read.reference_start is not None and read.reference_end is not None and read.reference_start < read.reference_end:
            if self.bFirstUMI:  # We are interested only in the first UMI
                sUMI = read.sUMIMod
                self.bFirstUMI = False
            else:
                sUMI = "0"  # Create dummy UMI

            # Create if not exist
            if not (sUMI in self.clusterData):
                self.clusterData[sUMI] = ClusterData(read.reference_name, sUMI)
            self.clusterData[sUMI].addRead(read)

    def correct(self):
        for cl in self.clusterData:
            clData = self.clusterData[cl]
            clData.calcOverlapped(self.refSeq)
            self.nCorrCount += 1
        return self.nCorrCount

    def getResult(self):
        resultStrings = []
        for cl in self.clusterData:
            clData = self.clusterData[cl]
            resultStrings += clData.createResultStrings()
        return resultStrings
