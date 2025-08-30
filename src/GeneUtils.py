import random
import numpy as np

def changeData(geneRoot, letter, base):
    if geneRoot.children != []:
        length = np.random.exponential(geneRoot.children[0].edgeLength)
        if length < 1:
            geneRoot.children[0].data = geneRoot.children[0].data + letter
            changeData(geneRoot.children[0], letter, base)
        else:
            letter1 = random.choice(base)
            indle = random.uniform(0, 1)
            if indle > 0.125:
                geneRoot.children[0].data = geneRoot.children[0].data + letter1
                changeData(geneRoot.children[0], letter1, base)

        length = np.random.exponential(geneRoot.children[1].edgeLength)
        if length < 1:
            geneRoot.children[1].data = geneRoot.children[1].data + letter
            changeData(geneRoot.children[1], letter, base)
        else:
            letter2 = random.choice(base)
            indle = random.uniform(0, 1)
            if indle > 0.125:
                geneRoot.children[1].data = geneRoot.children[1].data + letter2
                changeData(geneRoot.children[1], letter2, base)
    return geneRoot


def createGene(gene, geneRoot):
    if geneRoot.children != []:
        gene = createGene(gene, geneRoot.children[0])
        gene = createGene(gene, geneRoot.children[1])
    else:
        gene.append(geneRoot.data)
    return gene


def makeMatrixFromTree(geneRoot, myArray):
    if geneRoot.children != []:
        myString1 = makeMatrixFromTree(geneRoot.children[0], myArray)
        myString2 = makeMatrixFromTree(geneRoot.children[1], myArray)
    else:
        myArray.append(geneRoot.data)
    return myArray
