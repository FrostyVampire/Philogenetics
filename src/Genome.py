from Gene import Gene

class Genome:
    def __init__(self, genesNumber, geneLength, name="root"):
        self.name = name
        self.genes = [Gene(geneLength) for _ in range(genesNumber)]
        self.parent = None
        self.edgeLength = 0.0
        self.children = []

    # Create a child
    def evolve(self, edgeLength, name, insRate=0.0, delRate=0.0):
        child = Genome(0, 0, name)
        child.genes = [g.evolve(edgeLength, insRate, delRate) for g in self.genes]
        child.parent = self
        child.edgeLength = edgeLength
        self.children.append(child)
        return child

    def __str__(self):
        return f"{self.name}: {[str(g) for g in self.genes]}"
