import random
import numpy as np

BASES = ["A", "C", "G", "T"]

class Gene:
    def __init__(self, length, seq=None):
        self.seq = seq if seq is not None else "".join(random.choice(BASES) for _ in range(length))
        self.length = len(self.seq)

    # Evolve a gene through mutations
    def evolve(self, edgeLength, insRate=0.0, delRate=0.0):
        """
        Per-letter mutation:
        For each base, draw X ~ Poisson(edgeLength). If X >= 1, perform one event:
          - change (prob = 1 - insRate - delRate)
          - insert (prob = insRate)   -> insert one random base before current base
          - delete (prob = delRate)   -> drop current base
        """
        # Verify validity of insertions and deletions
        if insRate < 0 or delRate < 0 or insRate + delRate > 1:
            raise ValueError("Rates must satisfy: insRate >= 0, delRate >= 0, insRate + delRate <= 1")

        newSeq = []
        for base in self.seq:
            mutate = np.random.poisson(edgeLength) >= 1
            if mutate:
                action = random.choices(
                    ["change", "insert", "delete"],
                    weights=[1 - insRate - delRate, insRate, delRate],
                    k=1
                )[0]

                if action == "change":
                    base = random.choice([b for b in BASES if b != base])
                elif action == "insert":
                    newSeq.append(random.choice(BASES))
                elif action == "delete":
                    # skip adding this base
                    continue

            newSeq.append(base)

        return Gene(len(newSeq), "".join(newSeq))

    # To string
    def __str__(self):
        return self.seq
