# Main file
import numpy as np
import os
from Tree import Tree
import tkinter as tk
from tkinter import ttk, messagebox
import DistanceUtils
import matplotlib

from src import PlotUtils

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment

NumTries = 10

class NCDApp:
    def __init__(self, root):
        self.root = root
        self.root.title("NCD Tree Generator")

        # Default parameters
        self.params = {
            "Number of leaves": 30,
            "Genes per COG": 1,
            "COGs per genome": 50,
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
            genesPerCOG = int(self.entries["Genes per COG"].get())
            numCOGs = int(self.entries["COGs per genome"].get())
            geneLength = int(self.entries["Gene length"].get())
            exponentialMean = float(self.entries["Exponential mean"].get())
            insertionRate = float(self.entries["Insertion rate"].get())
            deletionRate = float(self.entries["Deletion rate"].get())
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numbers")
            return

        # --- Generate tree and genomes with multiple COGs ---
        final_leaves = None

        for cog_index in range(numCOGs):
            # Create a tree for this COG
            tree = Tree(
                leaves=numberOfLeaves,
                geneNum=genesPerCOG,
                geneLength=geneLength,
                exp=exponentialMean,
                insRate=insertionRate,
                delRate=deletionRate
            )

            if cog_index == 0:
                final_leaves = tree.leaves
            else:
                # Append genes to existing leaves
                for leaf_old, leaf_new in zip(final_leaves, tree.leaves):
                    leaf_old.genes.extend(leaf_new.genes)

        leaves = final_leaves
        numLeaves = len(leaves)

        # --- Compute NCD matrices and orthologs ---
        os.makedirs("output", exist_ok=True)
        orthologFilePath = os.path.join("output", "ortholog pairs.txt")
        orthologFile = open(orthologFilePath, "w")
        orthologFile.write("Pair\tGeneA\tGeneB\n")  # header

        ncdDistances = np.zeros((numLeaves, numLeaves))
        ncdMeans = np.zeros((numLeaves, numLeaves))

        for i in range(numLeaves):
            for j in range(i + 1, numLeaves):
                AllDistances = []
                AllMeans = []

                for t in range(NumTries):
                    Dist, Mean = DistanceUtils.createNCDMatrix(leaves[i], leaves[j])
                    AllDistances.append(Dist)
                    AllMeans.append(Mean)

                AvgDistance = np.mean(AllDistances, axis=0)
                AvgMean = np.mean(AllMeans, axis=0)

                ncdDistances[i][j] = ncdDistances[j][i] = np.mean(AvgDistance)
                ncdMeans[i][j] = ncdMeans[j][i] = np.mean(AvgMean)

                RowIndices, ColIndices = linear_sum_assignment(AvgDistance)
                OrthologPairs = list(zip(RowIndices, ColIndices))

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
        PlotUtils.plotNCD(ncdDistances, self.canvas, self.ax)

        messagebox.showinfo("Done", f"NCD matrix saved to {outputFilePath}")

    def savePlot(self):
        os.makedirs("output", exist_ok=True)
        output_file = os.path.join("output", "treePlot.png")
        self.fig.savefig(output_file, dpi=300)
        messagebox.showinfo("Saved", f"Plot saved to {output_file}")

if __name__ == "__main__":
    root = tk.Tk()
    app = NCDApp(root)
    root.mainloop()