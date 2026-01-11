# PlotUtils.py
import pandas as pd
import matplotlib.pyplot as plt

# Plot results from a CSV file with filtering and formatting.
def plotCsvResults(csv_file, output_image=None):
    data = pd.read_csv(csv_file, skiprows=1, header=None)
    selected = data.iloc[:, [0, 1, 2, 3, 6, 8]]
    selected.columns = ['Field0', 'Field1', 'Field2', 'Field3', 'Field6', 'Field8']

    filtered = selected[selected['Field1'] == 2].copy()
    filtered['YValue'] = filtered['Field6'] / filtered['Field3']

    colors = {v: c for v, c in zip(filtered['Field3'].unique(), ['blue', 'red', 'green', 'purple'])}
    line_styles = {v: s for v, s in zip(filtered['Field2'].unique(), ['-', '--', ':'])}

    plt.figure(figsize=(10, 6))
    for (field2, field3), subset in filtered.groupby(['Field2', 'Field3']):
        subset = subset.sort_values('Field0')
        plt.plot(subset['Field0'], subset['YValue'],
                 color=colors[field3], linestyle=line_styles[field2],
                 marker='o', label=f'alpha={field2}, gene_len={field3}')

    plt.xlabel('edge length')
    plt.ylabel('LZ-Time / geneLength')
    plt.title('LZ Times (|A-B|=02)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    if output_image:
        plt.savefig(output_image, dpi=300, bbox_inches='tight')
    plt.show()


# Create a graph based on NCD distances
def plotNCD(ncdDistances, canvas, ax):
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

def plotTrials(filePath):
    """
    Plot trial results from a results file.
    Expects a tab- or comma-delimited file with columns:
    seq, a, ei, fa, sc
    """

    # Read file into DataFrame
    try:
        df = pd.read_csv(filePath, sep=None, engine="python")  # auto-detect comma/tab
    except Exception as e:
        print(f"Error reading {filePath}: {e}")
        return

    # Validate required columns
    required_cols = {"seq", "a", "ei", "fa", "sc"}
    if not required_cols.issubset(df.columns):
        print(f"File {filePath} is missing required columns: {required_cols - set(df.columns)}")
        return

    # Derived values
    df["abLen"] = df["seq"].astype(str).str.len()
    df["edgeLen"] = df["ei"] * 0.05
    df["succ"] = 100 - df["fa"]
    df["meanScore"] = (50 - df["sc"]) * 2

    # Set up plot
    plt.figure(figsize=(8, 6))
    colors = {50: "r", 75: "g", 100: "b", 125: "m", 150: "c"}
    linestyles = {0: "solid", 1: "dashed", 2: "dotted", 3: "dashdot"}

    # Plot by (abLen, alpha)
    for (abLen, alpha), group in df.groupby(["abLen", "a"]):
        color = colors.get(abLen, "k")  # default black
        linestyle = linestyles.get(alpha, "solid")
        plt.plot(
            group["edgeLen"], group["meanScore"],
            label=f"{abLen} (α={alpha})",
            color=color,
            linestyle=linestyle
        )

    plt.xlabel("Edge Length")
    plt.ylabel("Mean Score")
    plt.title("Trial Results")
    plt.legend()
    plt.grid(True)
    plt.show()
