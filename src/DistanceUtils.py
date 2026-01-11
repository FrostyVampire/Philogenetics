# All the distance calculation functions
import itertools

import numpy as np
import parasail
import zlib
import concurrent.futures
from multiprocessing import Pool
from Levenshtein import distance

k = 3   # The k used by k-mers

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
def createDistanceMatrixParallel(originalGenes, inheritedGenes,matchScore=2, mismatchPenalty=-1, gapOpen=1, gapExtend=1, processes=None):
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

# Levenshtein distance between 2 given genomes
def levenshteinDistance(genome1, genome2):
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

            # ncdDistances[i][j] = epsilon + ((z12 + z21 - (z11 + z22)) / (z12 + z21))
            ncdDistances[i][j] = epsilon + ((z12+z21-2*min(z1,z2)) / 2*max(z1,z2))
            ncdMeans[i][j] = (z12 + z21) / 2.0

    return ncdDistances, ncdMeans

# Initiate the dictionary for k-mers of a given length (e.g. AGCT with length 3 will generate all possibilities AAA, AAC,..., TTT)
def initiateKMerDictionary(base, length):
    return {''.join(p): 0 for p in itertools.product(base, repeat=length)}

# Count how often each sequence of length k occurs in the gene and normalize to get frequencies
def computeKMerFrequencies(dictionary, gene, length):
    """
   Args:
        dictionary (dict): mapping k-mers to counts (initially 0)
        gene: object with a 'seq' attribute containing the sequence
        length (int): size of the k-mers (e.g., 3 for triplets)
   Returns:
       dict: updated dictionary mapping k-mers to their normalized frequencies
    """
    seq = gene.seq
    for i in range(len(seq) - length + 1):
        kmer = seq[i: i + length]
        dictionary[kmer] += 1

    total = len(seq) - length + 1
    for key in dictionary:
        dictionary[key] /= total
    return dictionary

# Computes the Hellinger distances between all genes of two genomes using k-mers
def hellinger(genome1, genome2, base):
    """
    Args:
        genome1, genome2: genome objects, each with .genes (list) and .genesNumber
        base (list[str]): alphabet of valid bases (e.g. ["A","C","G","T"])
    Returns:
        np.ndarray: matrix of Hellinger distances between genes of genome1 and genome2
    """
    epsilon = 1e-6
    baseKMerDict = initiateKMerDictionary(base, k)
    hellingerDistances = np.zeros((genome1.genesNumber, genome2.genesNumber))

    try:
        for i in range(genome1.genesNumber):
            kmerDict1 = computeKMerFrequencies(baseKMerDict.copy(), genome1.genes[i], k)
            for j in range(genome2.genesNumber):
                kmerDict2 = computeKMerFrequencies(baseKMerDict.copy(), genome2.genes[j], k)
                distanceSum = epsilon
                for kmer in kmerDict1:
                    diffSquared = (np.sqrt(kmerDict1[kmer]) - np.sqrt(kmerDict2[kmer])) ** 2
                    distanceSum += diffSquared
                hellingerDistances[i][j] = np.sqrt(distanceSum) / np.sqrt(2.0)
        return hellingerDistances
    except Exception as e:
        print(f"Hellinger distance exception at gene pair ({i}, {j}): {e}")
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

def calculateLempelZivDistance(seq1, seq2):
    """
    Compute the Lempel–Ziv distance between two sequences.
    This version is adapted from the prototype, cleaned up.
    """
    s1, s2 = "".join(map(str, seq1)), "".join(map(str, seq2))
    concat = s1 + s2

    c1 = lempelZivComplexity(s1)
    c2 = lempelZivComplexity(s2)
    c12 = lempelZivComplexity(concat)

    return (c12 - min(c1, c2)) / max(c1, c2)

# Compute Lempel-Ziv complexity (number of distinct substrings found sequentially).
def lempelZivComplexity(seq: str) -> int:
    i, k, l, n = 0, 1, 1, len(seq)
    c = 1
    while True:
        if i + k == n or seq[i + k] != seq[l + k]:
            if k > l:
                l = k
            i += 1
            if i == l:
                c += 1
                l = 1
                i = 0
                if l + 1 > n:
                    break
            k = 1
        else:
            k += 1
            if l + k > n:
                c += 1
                break
    return c

# Create distance matrix based on Lempel-Ziv complexity difference.
def createLZMatrix(genome1, genome2):
    epsilon = 1e-6
    lzDistances = np.zeros((genome1.genesNumber, genome2.genesNumber))
    for i in range(genome1.genesNumber):
        for j in range(genome2.genesNumber):
            g1 = str(genome1.genes[i].seq)
            g2 = str(genome2.genes[j].seq)
            c1, c2 = lempelZivComplexity(g1), lempelZivComplexity(g2)
            c12 = lempelZivComplexity(g1 + g2)
            lzDistances[i][j] = epsilon + (c12 - min(c1, c2)) / max(c1, c2)
    return lzDistances

# Compute the synteny index between two genomes.
def syntenyIndex(genome1, genome2):
    # SI = (2 * common genes) / (len(genome1) + len(genome2))
    set1, set2 = set(genome1), set(genome2)
    common = len(set1 & set2)
    return (2.0 * common) / (len(genome1) + len(genome2))

# Convert Synteny Index (SI) into a distance measure.
def exactDistanceFromSyntenyIndex(genome1, genome2):
    si = syntenyIndex(genome1, genome2)
    return 1.0 - si