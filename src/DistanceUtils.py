# All the distance calculation functions

import numpy as np
import parasail
import zlib
import concurrent.futures
from multiprocessing import Pool
from Levenshtein import distance

# Calculate Smith-Waterman distance between two sequences.
def calculateSwDistance(args):
    seq1, seq2, matchScore, mismatchPenalty, gapOpen, gapExtend = args

    matrix = parasail.matrix_create("ACGT", matchScore, mismatchPenalty)
    result = parasail.sw_trace(seq1, seq2, gapOpen, gapExtend, matrix)

    score = result.score
    maxPossibleScore = min(len(seq1), len(seq2)) * matchScore

    if maxPossibleScore == 0:
        return 0.0 if seq1 == seq2 else 1.0

    return 1.0 - (score / maxPossibleScore)

# Create a distance matrix parallel
def createDistanceMatrixParallel(originalGenes, inheritedGenes,
                                 matchScore=2, mismatchPenalty=-1,
                                 gapOpen=1, gapExtend=1, processes=None):
    pairs = [
        (originalGenes[i], inheritedGenes[j],
         matchScore, mismatchPenalty, gapOpen, gapExtend)
        for i in range(len(originalGenes))
        for j in range(len(inheritedGenes))
    ]

    with Pool(processes=processes) as pool:
        distances = pool.map(calculateSwDistance, pairs)

    distanceMatrix = np.array(distances).reshape(len(originalGenes),
                                                 len(inheritedGenes))
    return distanceMatrix

# Pairwise Smith-Waterman distances between two genomes.
def smithWatermanDistance(genome1, genome2):
    epsilon = 1e-6
    smithWatermanDists = np.zeros((genome1.genesNumber, genome2.genesNumber))
    go, ge = 1, 1
    matchScore = 5
    maxPossibleScore = 5000  # adjustable

    for i in range(genome1.genesNumber):
        for j in range(genome2.genesNumber):
            g1 = genome1.genes[i].seq
            g2 = genome2.genes[j].seq
            res = parasail.sw_trace(g1, g2, go, ge, parasail.dnafull)
            distance = 1.0 - (res.score / maxPossibleScore)
            smithWatermanDists[i][j] = distance + epsilon

    return smithWatermanDists


def klDivergence(p, q):
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

# Levenshtein
def distLev2(genome1, genome2):
    epsilon = 1e-6
    levensteinDistances = np.zeros((genome1.genesNumber, genome2.genesNumber))
    for i in range(genome1.genesNumber):
        for j in range(genome2.genesNumber):
            g1 = genome1.genes[i].seq
            g2 = genome2.genes[j].seq
            levensteinDistances[i][j] = distance(g1, g2) + epsilon
    return levensteinDistances

# Create an NCD (Normalized Compression Distance) given 2 genomes
def createNCDMatrix(genome1, genome2):
    epsilon = 1e-6
    ncdDistances = np.zeros((len(genome1.genes), len(genome2.genes)))
    ncdMeans = np.zeros((len(genome1.genes), len(genome2.genes)))

    for i in range(len(genome1.genes)):
        for j in range(len(genome2.genes)):
            seq1 = str.encode(str(genome1.genes[i].seq))
            seq2 = str.encode(str(genome2.genes[j].seq))

            z1 = len(zlib.compress(seq1))
            z2 = len(zlib.compress(seq2))
            z11 = len(zlib.compress(seq1 + seq1))
            z22 = len(zlib.compress(seq2 + seq2))
            z12 = len(zlib.compress(seq1 + seq2))
            z21 = len(zlib.compress(seq2 + seq1))

            ncdDistances[i][j] = epsilon + ((z12 + z21 - (z11 + z22)) / (z12 + z21))
            ncdMeans[i][j] = (z12 + z21) / 2.0

    return ncdDistances, ncdMeans


def initiateTripletDictionary(base):
    dictionary = {}
    triplet = ""
    for i in range(len(base)):
        triplet += base[i]
        for j in range(len(base)):
            triplet += base[j]
            for k in range(len(base)):
                triplet = triplet[0:2]
                triplet += base[k]
                dictionary[str(triplet)] = 0
            triplet = triplet[0:1]
        triplet = ""
    return dictionary


def createTripletDictionary(dictionary, gene):
    """
   Count how often each triplet (3-character sequence) occurs in the gene and normalize to get frequencies.

   Args:
       Dictionary mapping triplets to counts (initially 0)
       Gene object with a 'seq' attribute containing the sequence

   Returns:
       dict: updated dictionary mapping triplets to their normalized frequencies
    """
    for i in range(len(gene.seq) - 2):
        str1 = str(gene.seq[i: i + 3])
        dictionary[str1] += 1
    for key in dictionary.keys():
        dictionary[key] = dictionary[key] / (len(gene.data) - 2)
    return dictionary

def hellinger(genome1, genome2, base):
    """
    Computes the Hellinger distances between all genes of two genomes.

    Args:
        2 genomes
        baseTriplets: list of all possible triplets used to initialize dictionaries

    Returns:
        numpy.ndarray: matrix of Hellinger distances between genes of genome1 and genome2
    """
    epsilon = 1e-6
    baseTripletDictionary = initiateTripletDictionary(base)
    hellingerDistances = np.zeros((genome1.genesNumber, genome2.genesNumber))

    try:
        for i in range(len(genome1.genes)):
            tripletDict1 = createTripletDictionary(baseTripletDictionary.copy(), genome1.genes[i])
            for j in range(genome2.genesNumber):
                tripletDict2 = createTripletDictionary(baseTripletDictionary.copy(), genome2.genes[j])
                distanceSum = epsilon
                for triplet in tripletDict1.keys():
                    diffSquared = (np.sqrt(tripletDict1[triplet]) - np.sqrt(tripletDict2[triplet])) ** 2
                    distanceSum += diffSquared
                hellingerDistances[i][j] = np.sqrt(distanceSum) / np.sqrt(2.0)
        return hellingerDistances
    except Exception as e:
        print(f"Hellinger distance exception at gene pair ({i}, {j}):",e)
        exit()

# Run a function with timeout, return None if it fails.
def runWithTimeout(func, timeoutSeconds, *args, **kwargs):
    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(func, *args, **kwargs)
        try:
            return future.result(timeout=timeoutSeconds)
        except concurrent.futures.TimeoutError:
            print(f"Function {func.__name__} timed out after {timeoutSeconds} seconds")
            return None
        except Exception as e:
            print("runWithTimeout exception:", e)
            return None
