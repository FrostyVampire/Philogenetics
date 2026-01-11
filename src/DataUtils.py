# DataUtils.py
def reduceData(input_file, output_file):
    """
    Clean raw results file into a CSV-like format.
    """
    with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
        next(f_in)  # skip header
        for line in f_in:
            cleaned = (
                line.replace("timeNCD and timeEvoluate and lngthRelTot for edge length", "")
                .replace("for A-B length ", ",")
                .replace("for alfa ", ",")
                .replace(" for gene length ", ",")
                .replace(" genesInGenome ", ",")
                .replace(",,", ",")
                .replace(" loops ", "")
                .replace(" times", "")
            )
            f_out.write(cleaned)