import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, RegularPolygon
from numpy import radians as rad
from collections import Counter
from itertools import chain
from HMM_Sequence import HMM_Sequence
from HMM_transitions import HMM_transitions
from HMM_draw import HMM_draw
import sys


print(f"Arguments count: {len(sys.argv)}")
for i, arg in enumerate(sys.argv):
    print(f"Argument {i:>6}: {arg}")

filepath = sys.argv[1]#'C:\\Users\\hpatterton\\PycharmProjects\\HMMGraphic\\sequence3.txt'
hmm_transitions = HMM_transitions(filepath)
transition_probabilities = hmm_transitions.get_transition_probabilities()
print(transition_probabilities)
number_of_nodes = len(transition_probabilities)-1
hmm_draw = HMM_draw(number_of_nodes, transition_probabilities)
hmm_draw.draw()
hmm_draw.show_drawing()

