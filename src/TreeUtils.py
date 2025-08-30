import random
import numpy as np
from Node import Node   # assuming Node is defined in node.py

def chooseLeaf(root, nodeToAppend, numberOfLeaves):
    if root.children == []:
        root.children.append(nodeToAppend)
        return numberOfLeaves
    else:
        appendOrNot = random.randint(0, 1)
        if appendOrNot == 0:
            root.children.append(nodeToAppend)
            return numberOfLeaves + 1
        else:
            numOfSons = len(root.children)
            node = random.randint(1, numOfSons)
            node1 = root.children[node - 1]
            numberOfLeaves = chooseLeaf(node1, nodeToAppend, numberOfLeaves)
            return numberOfLeaves


def newickFormat(root, str1):
    if root.children != []:
        str1 = str1 + "("
        str1 = newickFormat(root.children[0], str1) + ","
        str1 = newickFormat(root.children[1], str1) + "):" + str(root.edgeLength)
        return str1
    else:
        str1 = str1 + root.data + ":" + str(root.edgeLength)
        return str1


def BuildTreeByText(newick, i, root, count):
    while i != count:
        if newick[i] == "(":
            root.children.append(Node("", root))
            root.children.append(Node("", root))
            root.children[0].parent = root
            root.children[1].parent = root
            root.children[0], i = BuildTreeByText(newick[0:], i+1, root.children[0], count)

        if newick[i] == ",":
            root.children[1], i = BuildTreeByText(newick[0:], i+1, root.children[1], count)

        if newick[i] == ":":
            if newick[i-1] == "(":
                ind = newick.find(",", i)
                edgeLen1 = newick[i+1:ind]
                root.edgeLength = edgeLen1
                return root, ind
            else:
                ind = newick.find(")", i)
                edgeLen1 = newick[i+1:ind]
                root.edgeLength = edgeLen1
                return root, ind

        if newick[i+1] == ":":
            ind1 = newick.find(",", i+1)
            ind2 = newick.find(")", i+1)
            ind = min(ind1, ind2)
            edgeLen1 = newick[i+2:ind]
            root.edgeLength = edgeLen1
            return root, ind

        i = i + 1

    return root, i
