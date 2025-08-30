import random
import numpy as np
from Genome import Genome

class Tree:
    def __init__(self, leaves, geneNum, geneLength, exp, insRate=0.0, delRate=0.0):
        """
        Build a random binary tree until we have `leaves` leaf genomes.
        Each new edge length is drawn from Exponential(mean=exp).
        """
        self.root = Genome(geneNum, geneLength, "root")
        self.leaves = [self.root]
        self.exp = float(exp)
        self.insRate = float(insRate)
        self.delRate = float(delRate)

        while len(self.leaves) < leaves:
            parent = random.choice(self.leaves)
            self.leaves.remove(parent)
            e1, e2 = np.random.exponential(self.exp, 2)
            c1 = parent.evolve(e1, f"{parent.name}_1", self.insRate, self.delRate)
            c2 = parent.evolve(e2, f"{parent.name}_2", self.insRate, self.delRate)
            self.leaves.extend([c1, c2])

    def printSubtree(self, node=None, level=0, showGenes=False):
        if node is None:
            node = self.root
        indent = "  " * level
        print(f"{indent}{node.name} (edge={node.edgeLength:.3f})")
        if showGenes:
            for i, g in enumerate(node.genes, start=1):
                print(f"{indent}  gene{i}: {g}")
        for child in node.children:
            self.printSubtree(child, level + 1, showGenes)
