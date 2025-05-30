#!/usr/bin/env python

import argparse
import sys
import gzip
import numpy as np
import itertools

import genomics

from time import sleep

from threading import Thread

from multiprocessing import Process

if sys.version_info>=(3,0):
    from multiprocessing import SimpleQueue
else:
    from multiprocessing.queues import SimpleQueue

import egglib


####################################################################################################################################


'''A function that reads from the window queue, calls some other function and writes to the results queue
This function needs to be tailored to the particular analysis funcion(s) you're using. This is the function that will run on each of the N cores.'''
def stats_wrapper(windowQueue, resultQueue, windType, genoFormat, sampleData, outgroupSamples, minSites, stats, doPops, skipPairs, indHet, eggW, eggB, maxMissing, analysis, minData, addWindowID=False, roundTo=4):
    while True:
        windowNumber,window = windowQueue.get() # retrieve window

        if windowNumber == -1:
            resultQueue.put((-1,None,False)) # this is the way of telling everything we're done
            break

        if windType == "coordinate" or windType == "predefined":
            scaf,start,end,mid,sites = (window.scaffold, window.limits[0], window.limits[1], window.midPos(),window.seqLen())
        else: scaf,start,end,mid,sites = (window.scaffold, window.firstPos(), window.lastPos(),window.midPos(),window.seqLen())
        if sites >= minSites:

            isGood = True
            #make alignment object
            Aln = genomics.genoToAlignment(window.seqDict(), sampleData, genoFormat = genoFormat)
            statsDict = {}
            if doPops: statsDict.update(genomics.popDiv(Aln, doPairs = not skipPairs))
            if indHet:
                hetDict = Aln.sampleHet()
                for key in hetDict.keys(): hetDict["het_" + key] = hetDict.pop(key)
                statsDict.update(hetDict)

            eggList = []    
            if eggW or eggB or args.calcMissing:             
                Alnsequences = []
                AlnsequencesOutgroup = []
                for i in range(len(Aln.groups)):
                    if outgroupSamples is not None:
                        if re.sub("_[AB]","",Aln.names[i]) not in outgroupSamples:
                            Alnsequences.append((Aln.names[i], "".join([str(x) for x in Aln.array[i]]), list(str(popNr[Aln.groups[i]]))))
                        else:
                            AlnsequencesOutgroup.append((Aln.names[i], "".join([str(x) for x in Aln.array[i]]), list(str(popNr[Aln.groups[i]]))))
                    else:
                        Alnsequences.append((Aln.names[i], "".join([str(x) for x in Aln.array[i]]), list(str(popNr[Aln.groups[i]]))))
       
                if eggW:
                    eggListW.append('lseff')
                    for pop in popNames:
                        if pop != 'OUT':
                            AlnsequencesPop = []
                            sampleNrPop = 0
                            for s in Alnsequences:
                                if s[2] == list(str(popNr[pop])):
                                    sampleNrPop += 1
                                    AlnsequencesPop.append(s)

                            AlnEgg = egglib.Align.create(AlnsequencesPop, alphabet=egglib.alphabets.DNA)

                            if outgroupSamples is not None:
    
                                for AlnSeqO in AlnsequencesOutgroup:
                                    AlnEgg.add_outgroup(AlnSeqO[0], AlnSeqO[1], int('999'))

                            eggStatistics = egglib.stats.ComputeStats() 
                            eggStatistics.add_stats(*eggListW)

                            eggStats = eggStatistics.process_align(AlnEgg, max_missing = maxMissing)

                            eggStats_copy = {**eggStats}
                            for key in eggStats_copy.keys():
                                eggStats[key+'_'+pop+'_egg'] = eggStats.pop(key)
                            statsDict.update(eggStats)

                
                if eggB:
                    eggListB.append('lseff')
                    for x in range(len(popNames)-1):
                        for y in range(x+1, len(popNames)):

                            popN1 = popNames[x]
                            popN2 = popNames[y]
                            
                            AlnsequencesPair = []
                            sampleNrPair = 0
                            sampleNrPop = 0
                            for s in Alnsequences:
                                
                                if s[2] == list(str(popNr[popN1])) or s[2] == list(str(popNr[popN2])):
                                    sampleNrPop += 1
                                    AlnsequencesPair.append(s)
                            
                            AlnEgg = egglib.Align.create(AlnsequencesPair, alphabet=egglib.alphabets.DNA)
                            struct = egglib.struct_from_labels(AlnEgg, lvl_pop=0)
                            
                            eggStatistics = egglib.stats.ComputeStats()
                            eggStatistics.set_structure(struct)
                            eggStatistics.add_stats(*eggListB)
                            
                            eggStats = eggStatistics.process_align(AlnEgg)
                            
                            eggStats_copy = {**eggStats}
                            for key in eggStats_copy.keys():
                                eggStats[key+'_'+popN1+'_'+popN2+'_egg'] = eggStats.pop(key)
                        
                            statsDict.update(eggStats)
                            
            if args.calcMissing:
                missing = {}
                for pop in popNames:

                    AlnsequencesPop = []
                    sampleNrPop = 0
                    for s in Alnsequences:
                        if s[2] == popNr[pop]:
                            sampleNrPop += 1

                            AlnsequencesPop.append(list(s[1]))

                    AlnNp = np.array(AlnsequencesPop)
                    
                    Ncolumn = (AlnNp != 'N').sum(0)
                    
                    genoPS = sum(Ncolumn/sampleNrPop)/len(Ncolumn)
                    
                    missing['genoPS_'+pop] = genoPS
                    
                    statsDict.update(missing)
                    
            values = [round(statsDict[stat], 4) if isinstance(statsDict[stat], float) else statsDict[stat] for stat in stats]
            #print(values)
                
        else:
            isGood = False
            values = [np.NaN]*len(stats)
        results = [] if not addWindowID else [window.ID]
        results += [scaf,start,end,mid,sites] + values
        resultString = ",".join([str(x) for x in results])
        resultQueue.put((windowNumber, resultString, isGood))


'''a function that watches the result queue and sorts results. This should be a generic funcion regardless of the result, as long as the first object is the result number, and this increases consecutively.'''
#def sorter(resultQueue, writeQueue, verbose):
  #global resultsReceived
  #sortBuffer = {}
  #expect = 0
    #threadsComplete = 0 #this will keep track of the worker threads and once they're all done this thread will break
  #while True:
    #resNumber,result,isGood = resultQueue.get()
    #resultsReceived += 1
    #if verbose:
      #sys.stderr.write("Sorter received result " + str(resNumber))
    #if resNumber == expect:
      #writeQueue.put((resNumber,result,isGood))
      #if verbose:
        #sys.stderr.write("Result {} sent to writer".format(resNumber))
      #expect +=1
      ##now check buffer for further results
      #while True:
        #try:
          #result,isGood = sortBuffer.pop(str(expect))
          #writeQueue.put((expect,result,isGood))
          #if verbose:
            #sys.stderr.write("Result {} sent to writer".format(expect))
          #expect +=1
        #except:
          #break
    #else:
      ##otherwise this line is ahead of us, so add to buffer dictionary
      #sortBuffer[str(resNumber)] = (result,isGood)

def sorter(resultQueue, writeQueue, verbose, nWorkerThreads):
    global resultsReceived
    sortBuffer = {}
    expect = 0
    threadsComplete = 0 #this will keep track of the worker threads and once they're all done this thread will break
    while True:
        windowNumber,result,isGood = resultQueue.get()
        #check if we're done
        if windowNumber == -1: threadsComplete += 1
        if threadsComplete == nWorkerThreads:
            writeQueue.put((-1,None,False))
            break #this is the way of telling everything we're done
        resultsReceived += 1
        if verbose:
            sys.stderr.write("Sorter received window {}\n".format(windowNumber))
        if windowNumber == expect:
            writeQueue.put((windowNumber,result,isGood))
            if verbose:
                sys.stderr.write("Slice {} sent to writer\n".format(windowNumber))
            expect +=1
            #now check buffer for further results
            while True:
                try:
                    result,isGood = sortBuffer.pop(str(expect))
                    writeQueue.put((expect,result,isGood))
                    if verbose:
                        sys.stderr.write("Slice {} sent to writer\n".format(expect))
                    expect +=1
                except:
                    break
        else:
            #otherwise this line is ahead of us, so add to buffer dictionary
            sortBuffer[str(windowNumber)] = (result,isGood)



'''a writer function that writes the sorted result. This is also generic'''
def writer(writeQueue, out, writeFailedWindows=False):
    global resultsWritten
    global resultsHandled
    while True:
        windowNumber,result,isGood = writeQueue.get()
        #check if we're done
        if windowNumber == -1: break
        if verbose:
            sys.stderr.write("Writer received result {}\n".format(windowNumber))
        if isGood or writeFailedWindows:
            out.write(result + "\n")
            resultsWritten += 1
        resultsHandled += 1


'''loop that checks stats'''
def checkStats():
  while True:
    sleep(10)
    sys.stderr.write("\n{} windows queued, {} results received, {} results written.\n".format(windowsQueued, resultsReceived, resultsWritten))


####################################################################################################################


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument("--windType", help="Type of windows to make", action = "store", choices = ("sites","coordinate","predefined"), default = "coordinate")
    parser.add_argument("-w", "--windSize", help="Window size in bases", type=int, action = "store", required = False, metavar="sites")
    parser.add_argument("-s", "--stepSize", help="Step size for sliding window", type=int, action = "store", required = False, metavar="sites")
    parser.add_argument("-m", "--minSites", help="Minumum good sites per window", type=int, action = "store", required = False, metavar="sites", default = 1)
    parser.add_argument("-O", "--overlap", help="Overlap for sites sliding window", type=int, action = "store", required = False, metavar="sites")
    parser.add_argument("-D", "--maxDist", help="Maximum span distance for sites window", type=int, action = "store", required = False)
    parser.add_argument("--windCoords", help="Window coordinates file (scaffold start end)", required = False)
    parser.add_argument("--maxMissing", help="Minumum sites that have to be genotyped at a position to be considered by egglib", type=float, action = "store", required = False, metavar="maxMissing", default = 0.25)

    parser.add_argument("--minData", help="Minumum proportion of individuals (or pairs) with >=minSites data", type=float, action = "store", required = False, metavar="prop", default = 0.01)
    
    parser.add_argument("-p", "--population", help="Pop name and optionally sample names (separated by commas)",
                        required = False, action='append', nargs="+", metavar=("popName","[samples]"))
    parser.add_argument("--popsFile", help="Optional file of sample names and populations", action = "store", required = False)
    parser.add_argument("--samples", help="Samples to include for individual analysis", action = "store", metavar = "sample names")
    
    parser.add_argument("--ploidy", help="Ploidy for each sample", action = "store", type=int, nargs="+")
    parser.add_argument("--ploidyFile", help="File with samples names and ploidy as columns", action = "store")
    parser.add_argument("--haploid", help="Alternatively just name samples that are haploid (comma separated)", action = "store", metavar = "sample names")
    parser.add_argument("--inferPloidy", help="Ploidy will be inferred in each window (NOT RECOMMENED)", action = "store_true")
    parser.add_argument("--skipPairs", help="Do not do pairwise statistics", action = "store_true")
    parser.add_argument("--indHet", help="Calculate individual heterozygosity", action = "store_true")
    parser.add_argument("--outgroup", help="list of outgroup samples (comma separated)", action = "store", metavar = "outgroup names")

    parser.add_argument("--analysis", help="Type of statistics to get (you can add multiple, separated by spaces)", action = "store", nargs = "+",
                        choices = ("popFreq","popDist", "popPairDist", "indPairDist","indHet"),
                        default = ("popDist", "popPairDist",))
    
    parser.add_argument("--roundTo", help="Round stats to X decimal places", type=int, default=4)
    
    parser.add_argument("-g", "--genoFile", help="Input genotypes file", required = False)
    parser.add_argument("-o", "--outFile", help="Results file", required = False)
    parser.add_argument("--exclude", help="File of scaffolds to exclude", required = False)
    parser.add_argument("--include", help="File of scaffolds to analyse", required = False)
    parser.add_argument("-f", "--genoFormat", help="Format of genotypes in genotypes file", action='store', choices = ("phased","pairs","haplo","diplo"), required = True)
    parser.add_argument("--header", help="Header text if no header in input", action = "store")
    parser.add_argument("-eggW","--egglibStatsWithin", help="Calculate egglib stats within populations. Give list of stats seperated by ,", action="store", required = False)
    parser.add_argument("-eggB","--egglibStatsBetween", help="Calculate egglib stats between populations. Give list of stats seperated by ,", action="store", required = False)
    parser.add_argument("--calcMissing", help="Calculate missing data for each populations as the average genotype coverage per site", action = "store_true")

    parser.add_argument("-T", "--threads", help="Number of worker threads for parallel processing", type=int, default=1, required = False, metavar="threads")
    parser.add_argument("--verbose", help="Verbose output", action="store_true")
    parser.add_argument("--addWindowID", help="Add window name or number as first column", action="store_true")
    parser.add_argument("--writeFailedWindows", help="Write output even for windows with too few sites.", action="store_true")
    
    
    args = parser.parse_args()
    
    #window parameters
    windType = args.windType
    
    if args.windType == "coordinate":
        assert args.windSize, "Window size must be provided."
        windSize = args.windSize
        stepSize = args.stepSize
        if not stepSize: stepSize = windSize
        assert not args.overlap, "Overlap does not apply to coordinate windows. Use --stepSize instead."
        assert not args.maxDist, "Maximum distance only applies to sites windows."
    
    elif args.windType == "sites":
        assert args.windSize, "Window size (number of sites) must be provided."
        windSize = args.windSize
        overlap = args.overlap
        if not overlap: overlap = 0
        maxDist = args.maxDist
        if not maxDist: maxDist = np.inf
        assert not args.stepSize, "Step size only applies to coordinate windows. Use --overlap instead."
    else:
        assert args.windCoords, "Please provide a file of window coordinates."
        assert not args.overlap, "Overlap does not apply for predefined windows."
        assert not args.maxDist, "Maximum does not apply for predefined windows."
        assert not args.stepSize,"Step size does not apply for predefined windows."
        assert not args.include,"You cannot only include specific scaffolds if using predefined windows."
        assert not args.exclude,"You cannot exclude specific scaffolds if using predefined windows."
        with open(args.windCoords,"rt") as wc: windCoords = tuple([(x,int(y),int(z),) for x,y,z in [line.split()[:3] for line in wc]])
    
    minSites = args.minSites
    if not minSites: minSites = windSize
    
    #file info
    genoFormat = args.genoFormat
    
    outFileName = args.outFile
    
    exclude = args.exclude
    include = args.include
    
    #other
    verbose = args.verbose
    
    
    ############## parse samples and populations
    popNames = []
    popInds = []
    allInds = []
    popNr = {}
    if args.population is not None:
        doPops = True
        n = 0
        for p in args.population:
            n+=1
            popNames.append(p[0])
            popNr[p[0]] = n
            if len(p) > 1: popInds.append(p[1].split(","))
            else: popInds.append([])

        if args.popsFile:
            with open(args.popsFile, "r") as pf: popDict = dict([ln.split() for ln in pf])
            for ind in popDict.keys():
                try: popInds[popNames.index(popDict[ind])].append(ind)
                except: pass

        for p in popInds: assert len(p) >= 1, "All populations must be represented by at least one sample."

        allInds += list(set([i for p in popInds for i in p]))
    else: doPops = False

    if args.samples is not None:
        allInds = list(set(allInds + args.samples.split(",")))
        args.indHet = True

    assert doPops or args.indHet, "Populations not specified, and individual het not requested. Nothing to do."

    outgroupSamples = None
    if args.outgroup is not None:
        outgroupSamples = args.outgroup.split(",")
        allInds = list(set(allInds + outgroupSamples))
    
    #if populations and samples not specified, just get all sample names from file
    if len(allInds) == 0:
        with gzip.open(args.genoFile, "r") if args.genoFile.endswith(".gz") else open(args.genoFile, "r") as gf:
            allInds = gf.readline().split()[2:]

    ploidyDict = dict(zip(allInds,[2]*len(allInds)))

    if args.haploid:
        for sample in args.haploid.split(","):
            ploidyDict[sample] = 1

    if args.ploidy is not None:
        ploidy = args.ploidy if len(args.ploidy) != 1 else args.ploidy*len(allInds)
        assert len(ploidy) == len(allInds), "Incorrect number of ploidy values supplied."
        ploidyDict = dict(zip(allInds,ploidy))
    elif args.ploidyFile is not None:
        with open(args.ploidyFile, "r") as pf: ploidyDict = dict([[s[0],int(s[1])] for s in [l.split() for l in pf]])
    elif args.inferPloidy:
        ploidyDict = dict(zip(allInds,[None]*len(allInds)))
    else:
        if args.genoFormat == "haplo": ploidyDict = dict(zip(allInds,[1]*len(allInds)))
        else: ploidyDict = dict(zip(allInds,[2]*len(allInds)))
        if args.haploid:
            for sample in args.haploid.split(","): ploidyDict[sample] = 1

    sampleData = genomics.SampleData(indNames = allInds, popNames = popNames, popInds = popInds, ploidyDict = ploidyDict)

    
    ############################################################################################################################################
    
    #open files
    
    if args.genoFile: genoFile = gzip.open(args.genoFile, "rt") if args.genoFile.endswith(".gz") else open(args.genoFile, "rt")
    else: genoFile = sys.stdin
    
    if args.outFile: outFile = gzip.open(args.outFile, "wt") if args.outFile.endswith(".gz") else open(args.outFile, "wt")
    else: outFile = sys.stdout
    
    if not args.addWindowID: outFile.write("scaffold,start,end,mid,sites,")
    else: outFile.write("windowID,scaffold,start,end,mid,sites,")
    
    ############################################################################################################################################
    
    #stats to output

    stats = []

    if args.indHet: stats += ["het_" + n for n in allInds]

    if doPops:
        stats += ["pi_" + n for n in popNames]

        if not args.skipPairs:
            stats += ["dxy_" + x + "_" + y for x,y in itertools.combinations(popNames, 2)]
            stats += ["Fst_" + x + "_" + y for x,y in itertools.combinations(popNames, 2)]

    popNames2 = popNames[:]
    if "OUT" in popNames2:
        popNames2.remove("OUT")

    if args.egglibStatsWithin:
        eggListW = args.egglibStatsWithin.split(',')
        for statistic in eggListW:
            stats += [statistic + "_" + n + "_egg" for n in popNames2]
        stats += ["lseff" + "_" + n + "_egg" for n in popNames2]
        
    if args.egglibStatsBetween:
        eggListB = args.egglibStatsBetween.split(',')
        for statistic in eggListB:
            stats += [statistic +  "_" + x + "_" + y + "_egg" for x,y in itertools.combinations(popNames2, 2)]
        stats += ["lseff" +  "_" + x + "_" + y + "_egg" for x,y in itertools.combinations(popNames2, 2)]

    if args.calcMissing:
        stats += ['genoPS'+ "_" + n for n in popNames2]
    
    # print stats
    outFile.write(",".join(stats) + "\n")
    
    ##############################################################
    
    #scafs to exclude
    
    if exclude:
        scafsFile = open(exclude, "rU")
        scafsToExclude = [line.rstrip() for line in scafsFile.readlines()]
        sys.stderr.write("{} scaffolds will be excluded.".format(len(scafsToExclude)))
        scafsFile.close()
    else:
        scafsToExclude = None
    
    if include:
        scafsFile = open(include, "rU")
        scafsToInclude = [line.rstrip() for line in scafsFile.readlines()]
        sys.stderr.write("{} scaffolds will be analysed.".format(len(scafsToInclude)))
        scafsFile.close()
    else:
        scafsToInclude = None
    
    
    ##########################################################################################################
    
    #counting stat that will let keep track of how far we are
    windowsQueued = 0
    resultsReceived = 0
    resultsWritten = 0
    resultsHandled = 0
    
    '''Create queues to hold the data one will hold the line info to be passed to the analysis'''
    windowQueue = SimpleQueue()
    #one will hold the results (in the order they come)
    resultQueue = SimpleQueue()
    #one will hold the sorted results to be written
    writeQueue = SimpleQueue()


    '''start worker Processes for analysis. The comand should be tailored for the analysis wrapper function
    of course these will only start doing anything after we put data into the line queue
    the function we call is actually a wrapper for another function.(s) This one reads from the line queue, passes to some analysis function(s), gets the results and sends to the result queue'''
    workerThreads = []
    sys.stderr.write("\nStarting {} worker threads\n".format(args.threads))
    for x in range(args.threads):
        workerThread = Process(target=stats_wrapper, args = (windowQueue, resultQueue, windType, genoFormat, sampleData, outgroupSamples, minSites,
                                                 stats, doPops, args.skipPairs, args.indHet, args.egglibStatsWithin, args.egglibStatsBetween, args.maxMissing, args.analysis, args.minData, args.addWindowID, args.roundTo))
        workerThread.daemon = True
        workerThread.start()
        workerThreads.append(workerThread)


    '''thread for sorting results'''
    sorterThread = Thread(target=sorter, args=(resultQueue,writeQueue,verbose,args.threads,))
    sorterThread.daemon = True
    sorterThread.start()
    
    '''start thread for writing the results'''
    writerThread = Thread(target=writer, args=(writeQueue, outFile, args.writeFailedWindows,))
    writerThread.daemon = True
    writerThread.start()
    
    
    '''start background Thread that will run a loop to check run statistics and print
    We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
    checkerThread = Thread(target=checkStats)
    checkerThread.daemon = True
    checkerThread.start()    
    
    
    
    ##########################################################
    
    #get windows and analyse
    if windType == "coordinate": windowGenerator = genomics.slidingCoordWindows(genoFile, windSize, stepSize,
                                                                                headerLine = args.header,
                                                                                names = sampleData.indNames,
                                                                                include = scafsToInclude,
                                                                                exclude = scafsToExclude)
    elif windType == "sites": windowGenerator = genomics.slidingSitesWindows(genoFile, windSize, overlap,
                                                                            maxDist, minSites,
                                                                            headerLine = args.header,
                                                                            names = sampleData.indNames,
                                                                            include = scafsToInclude,
                                                                            exclude = scafsToExclude)
    else: windowGenerator = genomics.predefinedCoordWindows(genoFile, windCoords,
                                                            headerLine = args.header,
                                                            names = sampleData.indNames)
    
    
    for window in windowGenerator:
        windowQueue.put((windowsQueued,window))
        windowsQueued += 1
    
    ############################################################################################################################################
    
    #Now we send completion signals to all worker threads
    for x in range(args.threads):
        windowQueue.put((-1,None,)) # -1 tells the threads to break
    
    sys.stderr.write("\nWaiting for all threads to finish\n".format(args.threads))
    for x in range(len(workerThreads)):
        workerThreads[x].join()
    
    sorterThread.join()
    writerThread.join()
    
    genoFile.close()
    outFile.close()
    
    sys.stderr.write(str(windowsQueued) + " windows were tested.\n")
    sys.stderr.write(str(resultsWritten) + " results were written.\n")
    
    sys.stderr.write("\nDone.\n")
    
    sys.exit()

