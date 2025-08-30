# Main file
import numpy as np
import os
from Tree import Tree
from DistanceUtils import createNCDMatrix
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk, messagebox
from DistanceUtils import createNCDMatrix
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment

NumTries = 10

def main():
    # --- Input parameters ---
    numberOfLeaves = 30 # int(input("Number of leaves: "))
    genesPerGenome = 50 # int(input("Genes per genome: "))
    geneLength = 10 # int(input("Gene length: "))
    exponentialMean = 0.1 #float(input("Exponential mean (scale): "))
    insertionRate = 0.05 #float(input("Insertion rate (0-1): "))
    deletionRate = 0.05 #float(input("Deletion rate (0-1): "))

    # Create the tree
    tree = Tree(
        leaves=numberOfLeaves,
        geneNum=genesPerGenome,
        geneLength=geneLength,
        exp=exponentialMean,
        insRate=insertionRate,
        delRate=deletionRate
    )

    # # Print the tree
    # print("\nGenerated Tree (structure only):")
    # tree.printSubtree(showGenes=False)

    # # Print all genomes at leaves
    # print("\nLeaf genomes:")
    # for leaf in tree.leaves:
    #     print(leaf)

    # Get all leaves
    leaves = tree.leaves
    numLeaves = len(leaves)

    # Initialize matrices
    ncdDistances = np.zeros((numLeaves, numLeaves))
    ncdMeans = np.zeros((numLeaves, numLeaves))

    # Compute NCD for every pair of leaves
    for i in range(numLeaves):
        for j in range(i, numLeaves):  # only compute upper triangle, symmetric
            dist, mean = createNCDMatrix(leaves[i], leaves[j])
            # Since createNCDMatrix returns matrices per gene, we can average or sum
            ncdDistances[i][j] = ncdDistances[j][i] = np.mean(dist)
            ncdMeans[i][j] = ncdMeans[j][i] = np.mean(mean)

    # Save results to file
    outputFilePath = os.path.join("input", "resultsNCD2.txt")
    with open(outputFilePath, "w") as file:
        file.write("NCD Distance Matrix:\n")
        for row in ncdDistances:
            file.write("\t".join(f"{dist:.4f}" for dist in row) + "\n")

    print(f"\nNCD matrix saved to {outputFilePath}")

    plotTree()

# Create a graph
def plotTree(ncdDistances, canvas, ax):
    x, y = [], []
    for i in range(min(15, len(ncdDistances))):
        x.append(0.05 * i)
        if i < len(ncdDistances) and len(ncdDistances[i]) > 2:
            try:
                ncdOr1 = ncdDistances[i][1] / (ncdDistances[i][2] + 1e-6)
                y.append(ncdOr1 * 100)
            except ZeroDivisionError:
                y.append(0)

    ax.clear()
    ax.plot(x[:len(y)], y, "m", linestyle="solid", label="mean NCD % ratio")
    ax.set_xlabel("Edge length")
    ax.set_ylabel("Percentage")
    ax.legend(loc="best")
    canvas.draw()

class NCDApp:
    def __init__(self, root):
        self.root = root
        self.root.title("NCD Tree Generator")

        # Default parameters
        self.params = {
            "Number of leaves": 30,
            "Genes per genome": 50,
            "Gene length": 10,
            "Exponential mean": 0.1,
            "Insertion rate": 0.05,
            "Deletion rate": 0.05,
        }

        # Input form
        self.entries = {}
        row = 0
        for label, default in self.params.items():
            ttk.Label(root, text=label).grid(row=row, column=0, sticky="w", padx=5, pady=5)
            entry = ttk.Entry(root)
            entry.insert(0, str(default))
            entry.grid(row=row, column=1, padx=5, pady=5)
            self.entries[label] = entry
            row += 1

        # Buttons
        run_button = ttk.Button(root, text="Run", command=self.runAnalysis)
        run_button.grid(row=row, column=0, columnspan=2, pady=10)

        save_button = ttk.Button(root, text="Save Plot", command=self.savePlot)
        save_button.grid(row=row, column=2, padx=5, pady=5)


        # Matplotlib figure
        self.fig, self.ax = plt.subplots(figsize=(8, 5))
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().grid(row=row+1, column=0, columnspan=2, padx=5, pady=5)

    def runAnalysis(self):
        try:
            # Parse input values
            numberOfLeaves = int(self.entries["Number of leaves"].get())
            genesPerGenome = int(self.entries["Genes per genome"].get())
            geneLength = int(self.entries["Gene length"].get())
            exponentialMean = float(self.entries["Exponential mean"].get())
            insertionRate = float(self.entries["Insertion rate"].get())
            deletionRate = float(self.entries["Deletion rate"].get())
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numbers")
            return

        # Create the tree
        tree = Tree(
            leaves=numberOfLeaves,
            geneNum=genesPerGenome,
            geneLength=geneLength,
            exp=exponentialMean,
            insRate=insertionRate,
            delRate=deletionRate
        )

        # Compute NCD matrices
        # Open file
        os.makedirs("output", exist_ok=True)
        orthologFilePath = os.path.join("output", "ortholog pairs.txt")
        orthologFile = open(orthologFilePath, "w")
        orthologFile.write("Pair\tGeneA\tGeneB\n")  # header

        leaves = tree.leaves
        numLeaves = len(leaves)
        ncdDistances = np.zeros((numLeaves, numLeaves))
        ncdMeans = np.zeros((numLeaves, numLeaves))

        # Initialize distance matrices
        NcdDistances = np.zeros((numLeaves, numLeaves))
        NcdMeans = np.zeros((numLeaves, numLeaves))

        for i in range(numLeaves):
            for j in range(i + 1, numLeaves):
                AllDistances = []
                AllMeans = []

                # Monte Carlo / repeated trials
                for t in range(NumTries):
                    Dist, Mean = createNCDMatrix(leaves[i], leaves[j])
                    AllDistances.append(Dist)
                    AllMeans.append(Mean)

                # Average over trials
                AvgDistance = np.mean(AllDistances, axis=0)
                AvgMean = np.mean(AllMeans, axis=0)

                # Store mean values in matrices
                NcdDistances[i][j] = NcdDistances[j][i] = np.mean(AvgDistance)
                NcdMeans[i][j] = NcdMeans[j][i] = np.mean(AvgMean)

                # Compute optimal gene matching
                RowIndices, ColIndices = linear_sum_assignment(AvgDistance)
                OrthologPairs = list(zip(RowIndices, ColIndices))

                # Write ortholog pairs in readable format
                orthologFile.write(f"Leaf pair {i}-{j}\n")
                for a, b in OrthologPairs:
                    orthologFile.write(f"{a}\t{b}\n")
                orthologFile.write("\n")


        orthologFile.close()
        print(f"Ortholog pairs saved to {orthologFilePath}")

        # Save results
        os.makedirs("input", exist_ok=True)
        outputFilePath = os.path.join("input", "resultsNCD2.txt")
        with open(outputFilePath, "w") as file:
            file.write("NCD Distance Matrix:\n")
            for row in ncdDistances:
                file.write("\t".join(f"{dist:.4f}" for dist in row) + "\n")

        # Update plot
        plotTree(ncdDistances, self.canvas, self.ax)

        messagebox.showinfo("Done", f"NCD matrix saved to {outputFilePath}")

    # Save graph
    def savePlot(self):
        os.makedirs("output", exist_ok=True)
        output_file = os.path.join("output", "treePlot.png")
        self.fig.savefig(output_file, dpi=300)
        messagebox.showinfo("Saved", f"Plot saved to {output_file}")

if __name__ == "__main__":
    root = tk.Tk()
    app = NCDApp(root)
    root.mainloop()