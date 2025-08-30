from TreeUtils import chooseLeaf

class Node:
    def __init__(self, data="", parent=None):
        self.data = data
        self.children = []
        self.parent = parent
        self.edgeLength = 0.0

    def addChildren(self, numberOfLeaves, node1):
        numberOfLeaves = chooseLeaf(self, node1, numberOfLeaves)
        return numberOfLeaves
