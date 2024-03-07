from collections import Counter

import itertools
import json
import sys
from sys import stderr
import minorminer
import dimod

import numpy as np
import matplotlib
matplotlib.use("agg")    # must select backend before importing pyplot
import matplotlib.pyplot as plt
from dimod import BinaryQuadraticModel
from dwave.system import LeapHybridSampler
from neal import SimulatedAnnealingSampler
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import FixedEmbeddingComposite, EmbeddingComposite
from pulp import *

from exact_cover import exact_cover_bqm

def build_subsets(n):
    """Optimized and corrected version to generate subsets of constraints for N-Queens problem."""
    subsets = []
    for x in range(n):
        for y in range(n):
            col = x
            row = y + n

            subset = {col, row}

            # Diagonal
            diag = x + y
            subset.add(diag + 2*n)

            # Anti-diagonal
            anti_diag = (n - 1 - x + y)
            subset.add(anti_diag + 4*n - 1) 

            subsets.append(subset)

    return subsets

def add_constraint(bqm, S, x1, x2):
    lamda = 3

    bqm.add_variable(S, lamda * (1/2))
    bqm.add_variable(x1, lamda * (1/2))
    bqm.add_variable(x2, lamda * (1/2))
    
    bqm.add_interaction(S, x1, lamda * (-1))
    bqm.add_interaction(S, x2, lamda * (-1))
    bqm.add_interaction(x1, x2, lamda * 2)

def handle_diag_constraints_S(bqm, subsets, diag_constraints):
    """Update bqm with diagonal (and anti-diagonal) constraints.
    Duplicates are penalized.
    """
    for constraint in diag_constraints:
        S = bqm.num_variables
        subset_indices = [i for i, subset in enumerate(subsets) if constraint in subset]

        if len(subset_indices) >= 2:
            add_constraint(bqm, S, subset_indices[0], subset_indices[1])
            S += 1
            for i in range(2, len(subset_indices)):
                add_constraint(bqm, S, S-1, subset_indices[i])

    return bqm

def handle_diag_constraints(bqm, subsets, diag_constraints):
    """Update bqm with diagonal (and anti-diagonal) constraints.
    Duplicates are penalized.
    """
    for constraint in diag_constraints:
        subset_indices = [i for i, subset in enumerate(subsets) if constraint in subset]

        for i in range(len(subset_indices)):
            for j in range(i):
                if subset_indices[i] != subset_indices[j]:
                    bqm.add_interaction(subset_indices[i], subset_indices[j], 2)
    return bqm

def n_queens(n, sampler=None, num_reads=2000):
    """Returns a potential solution to the n-queens problem in a list of sets,
    each containing constraint IDs representing a queen's location.

    Args:
        n: Number of queens to place.

        sampler: A binary quadratic model sampler. Defaults to dwave-system's LeapHybridSampler.
    """
    num_row_col_constraints = 2 * n
    row_col_constraint_ids = set(range(num_row_col_constraints))

    num_diag_constraints = 4 * n - 6   # includes anti-diag constraints
    diag_constraint_ids = set(range(num_row_col_constraints, num_row_col_constraints + num_diag_constraints))

    # Build subsets of constraint IDs. Each subset will become a variable in our BQM.
    subsets = build_subsets(n)

    # Build BQM with only row/col constraints
    bqm = exact_cover_bqm(row_col_constraint_ids, subsets)

    out_path = "embedding_output.json" 

    # Add diag/anti-diag constraints - duplicates are penalized.
    bqm = handle_diag_constraints_S(bqm, subsets, diag_constraint_ids)

    if sampler is None:
        # sampler = LeapHybridSampler()
        # sampler = SimulatedAnnealingSampler()
        sampler = EmbeddingComposite(DWaveSampler(token='DEV-e0ac368d04813d5d0a2019a61e39c30c446c6397'))
        
    sampleset = sampler.sample(bqm, num_reads=num_reads)

    sample = sampleset.first.sample
    # for run in sample:
    #     print(run.energy)
    print(sample)

    return [subsets[i] for i in range(len(subsets)) if sample[i]]

def is_valid_solution(n, solution):
    """Check that solution is valid by making sure all constraints were met.

    Args:
        n: Number of queens in the problem.

        solution: A list of sets, each containing constraint IDs that represent
                  a queen's location.
    """
    count = Counter()

    for queen in solution:
        count = count + Counter(queen)

    # Check row/col constraints
    for i in range(2*n):
        if count[i] != 1:
            if i < n:
                col = i
                print("Column {} has {} queens.".format(col, count[i]))
            else:
                row = np.abs(i - (2*n - 1)) # Convert constraint id to row index
                print("Row {} has {} queens.".format(row, count[i]))

            return False

    # Check diag/anti-diag constraints
    for i in range(2*n, 4*n - 6):
        if count[i] > 1:
            if i <= 4*n - 4:
                print("Top-left to bottom-right diagonal {} has {} queens.".format(i, count[i]))
            else:
                print("Bottom-left to top-right diagonal {} has {} queens.".format(i, count[i]))

            return False

    return True

def plot_chessboard(n, queens):
    """Create a chessboard with queens using matplotlib. Image is saved
    in the root directory. Returns the image file name.
    """
    chessboard = np.zeros((n,n))
    chessboard[1::2,0::2] = 1
    chessboard[0::2,1::2] = 1

    # Adjust fontsize for readability
    if n <= 10:
        fontsize = 30
    elif n <= 20:
        fontsize = 10
    else:
        fontsize = 5

    plt.xticks(np.arange(n))
    plt.yticks(np.arange(n))

    plt.imshow(chessboard, cmap='binary')

    # Place queens
    for subset in solution:
        x = y = -1
        for constraint in subset:
            if constraint < n:
                x = constraint
            elif constraint >= n and constraint < 2*n:
                y = np.abs(constraint - (2*n - 1)) # Convert constraint ID to row index

        if x != -1 and y != -1:
            plt.text(x, y, u"\u2655", fontsize=fontsize, ha='center',
                     va='center', color='black' if (x - y) % 2 == 0 else 'white')

    # Save file in root directory
    file_name = "{}-queens-solution.png".format(n)
    plt.savefig(file_name)

    return file_name

def get_sanitized_input():
    while True:
        print("Enter the number of queens to place (n > 0):")
        n = input()

        try:
            n = int(n)
            if n <= 0:
                print("Input must be greater than 0.")
                continue
            if n >= 200:
                # Run but give a warning
                print("Problems with large n will run very slowly.")

        except ValueError:
            print("Input type must be int.")
            continue

        return n

if __name__ == "__main__":
    n = get_sanitized_input()

    if n > 20:
        print("Solution image is large and may be difficult to view.")
        print("Plot settings in plot_chessboard() may need adjusting.")

    print("Trying to place {n} queens on a {n}*{n} chessboard.".format(n=n))
    solution = n_queens(n)

    if is_valid_solution(n, solution):
        print("Solution is valid.")
    else:
        print("Solution is invalid.")

    file_name = plot_chessboard(n, solution)
    print("Chessboard created. See: {}".format(file_name))
