import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, RegularPolygon
from numpy import radians as rad
from collections import Counter
from itertools import chain
from HMM_draw import HMM_draw

class HMM:

    def __init__(self, number_nodes, molecule_type):
        self.dictionary_of_transitions = {'M0_M1':0,'M0_I0':1,'M0_D1':2,'I0_M1':3,'I0_I0':4,'I0_D1':5,'MX_MY':6,'MX_IX':7,'MX_DY':8,'IX_MY':9,'IX_IX':10,'IX_DY':11,'DX_MY':12,'DX_IX':13,'DX_DY':14,'ML_Mend':15,'ML_IL':16,'IL_Mend':17,'IL_IL':18,'DL_Mend':19,'DL_IL':20}
        self.dictionary_of_emmission_nodes = {'M':0,'I':1}
        self.dictionary_of_amino_acids = {'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19}
        self.dictionary_of_nucleotides = {'G':0,'A':1,'T':2,'C':3}
        # 0 is nucleotide/aa, 1 is node number (0 is start, N is end), 2 is M or I
        self.molecule_type = molecule_type
        if self.molecule_type == 'P': # make sure it is not a 'shallow' list
            self.array_of_emmisions = [[[0 for i in range(2)] for j in range(number_nodes+2)] for k in range(len(self.dictionary_of_amino_acids))]
        else:
            self.array_of_emmisions = [[[0 for i in range(2)] for j in range(number_nodes+2)] for k in range(len(self.dictionary_of_nucleotides))]
        print(self.array_of_emmisions)
        self.number_of_nodes = number_nodes
        self.array_of_transitions = [[0]*21]*(number_nodes+2) # 21 is the sum of all possible transitions listed in self.dictionary_of_transitions
        self.array_of_nodes = [[0],[0],[0]]
        self.node_type_indexes = {'M':0,'I':1,'D':2}
        self.make_end_node()
        for node_number in range(number_nodes):
            for node_type in range(len(self.node_type_indexes)):
                self.make_node([], node_type)
        self.make_start_node()
        return

    def set_transition_value(self,value, node_index, sub_node_index):
        self.array_of_transitions[node_index][sub_node_index] = value
        return

    def get_transition_value(self, node_index, sub_node_index):
        return(self.array_of_transitions[node_index][sub_node_index])

    def make_end_node(self):
        node = HMM_node([],'M',[-1,-1,-1])
        self.array_of_nodes[self.node_type_indexes['M']][0] = node
        #self.array_of_nodes[0][self.node_type_indexes['I']] = 0
        return(node)

    def make_node(self, emmission_array, node_type_index):
        # number of M nodes currently
        index_of_this_node = len(self.array_of_nodes[0])
        # # point to earlier nodes
        if node_type_index == self.node_type_indexes['I']:
            child_pointers = [index_of_this_node-1,-1,index_of_this_node-1]
        elif node_type_index == self.node_type_indexes['M'] or node_type_index == self.node_type_indexes['D']:
            child_pointers = [index_of_this_node-1,index_of_this_node-1,index_of_this_node-1]
        node = HMM_node([], node_type_index, child_pointers)
        self.array_of_nodes[node_type_index].append(node)
        return

    def make_start_node(self):
        index_of_this_node = len(self.array_of_nodes[0])
        child_pointers = [index_of_this_node-1, index_of_this_node-1, index_of_this_node-1]
        node = HMM_node([], self.node_type_indexes['M'], child_pointers)
        self.array_of_nodes[self.node_type_indexes['M']].append(node)
        return

    def list_nodes(self):
        #print(self.array_of_nodes)
        number_of_nodes = self.array_of_nodes[0]
        print(number_of_nodes)
        # print(self.array_of_nodes[0][0].get_node_info())
        # print(self.array_of_nodes[1][0].get_node_info())
        # print(self.array_of_nodes[1][1].get_node_info())
        # print(self.array_of_nodes[1][2].get_node_info())
        # for node in range(0,number_of_nodes):
        #     for i in range(0,3):
        #         if self.array_of_nodes[node][i] > 0:
        #             print(node, self.array_of_nodes[node][i].get_node_info())
        return

    def Read_FastA(path):
        f = open(path, 'r')
        name = []
        sequence = []
        current_sequence = ''
        for line in f:
            if line[0] == '>':
                name.append(line[1:].rstrip('\n'))
                if len(current_sequence) > 0:
                    sequence.append(current_sequence)
                current_sequence = ''
            else:
                current_sequence += line.rstrip('\n')
        sequence.append(current_sequence)
        return (name, sequence)







class HMM_node:

    def __init__(self, list_of_symbols, type_of_node, indexes_of_child_nodes):
        #self.emmission_symbols = list_of_symbols
        self.node_type = type_of_node
        self.child_indexes = indexes_of_child_nodes
        return

    def get_node_type(self):
        print(self.node_type)
        return

    def get_node_info(self):
        print('type of mode = ',self.node_type, "child pointers=",self.child_indexes)
        return

class HMM_sequence:

    def __init__(self):
        self.number_of_sequences = 0
        self.start_block_is_conserved = True
        self.name = []
        self.sequence = []
        self.sequence_array = np.array([['A','-','-','-','-','T','C'],['A','C','-','-','A','A','A'],['G','A','-','-','C','A','A'],['G','A','T','-','C','A','A'],['G','T','G','A','C','-','A']])
        self.sequence_blocks = [[],[]]
        self.order_of_blocks = []
        self.indexes_of_block_order = []
        self.column_categories = []
        self.current_M_column = []
        self.next_M_column = []
        self.D_column = []
        self.current_I_array = []
        self.current_hmm_pattern_node = 0
        self.conservation_cutoff = 0.75
        self.dictionary_of_amino_acids = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8,
                                          'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16,
                                          'V': 17, 'W': 18, 'Y': 19}
        self.dictionary_of_nucleotides = {'G': 0, 'A': 1, 'T': 2, 'C': 3}
        self.dictionary_of_transition_paths = {'MM':0,'MI':1,'MD':2,'IM':3,'II':4,'ID':5,'DM':6,'DI':7,'DD':8}
        self.number_of_paths = len(self.dictionary_of_transition_paths)
        self.transition_paths = []

        return

    def Read_FastA(self, path):
        try:
            with open(path, 'r') as f:
                pass
        except FileNotFoundError:
            print('File not found')
            return False
        current_sequence = ''
        for line in f:
            if line[0] == '>':
                self.name.append(line[1:].rstrip('\n'))
                if len(current_sequence) > 0:
                    self.sequence.append(list(current_sequence))
                current_sequence = ''
            else:
                current_sequence += line.rstrip('\n')
        self.sequence.append(list(current_sequence))
        self.sequence_array = np.array(self.sequence)
        if len(self.sequence) > 0:
            number_of_sequences = len(self.sequence)
        else:
            return False
        sequence_length = self.sequence[0]
        result = True
        for i in range(1,number_of_sequences):
            if sequence_length == len(self.sequence[i]):
                result = False
        if result == False:
            return False

        return True

    def get_pseudocounts(self, numerator, denominator, molecule_type):
        if molecule_type == 'P':
            return((numerator+1)/(denominator+20))
        else:
            return((numerator+1)/(denominator+4))

    def get_emmission_list(self,list_of_symbols, molecule_type):
        if molecule_type == 'P':
            emmision_array = [0 for i in range(len(self.dictionary_of_amino_acids))]
            symbol_count = Counter(list_of_symbols)
            total_symbols = sum(symbol_count.values())
            for character in self.dictionary_of_amino_acids:
                numerator = symbol_count[character]
                denominator = total_symbols
                # You will need to call pseudocounts here if you want to use it
                emmision_array[self.dictionary_of_amino_acids[character]] = numerator/denominator
        elif molecule_type == 'D':
            emmision_array = emmision_array = [0 for i in range(len(self.dictionary_of_nucleotides))]
            symbol_count = Counter(list_of_symbols)
            total_symbols = sum(symbol_count.values())
            for character in self.dictionary_of_nucleotides:
                numerator = symbol_count[character]
                denominator = total_symbols
                # You will need to call pseudocounts here if you want to use it
                entry = self.dictionary_of_nucleotides[character]
                emmision_array[entry] = numerator#/denominator
        else:
            emmision_array = []

        return emmision_array

    def convert_list_to_1D(self, array_block):
        list_array = array_block.tolist()
        list_1D = list(chain.from_iterable(list_array))
        return list_1D

    def remove_dashes_from_list(self, list_1D):
        list_dashes_stripped = [symbol for symbol in list_1D if symbol != '-']
        return list_dashes_stripped

    def remove_dashes_from_np_array(self, array):
        new_array = []
        list_array = array.tolist()
        for i in range(len(list_array)):
            row = [symbol for symbol in list_array[i][:] if symbol != '-']
            new_array.append(row)
        return new_array

    def remove_dashes_from_array(self, array):
        new_array = []
        for i in range(len(array)):
            row = [symbol for symbol in array[i][:] if symbol != '-']
            new_array.append(row)
        return new_array

    def is_conserved(self, row):
        row_as_list = row.tolist()
        length_of_block = len(row_as_list)
        number_of_hyphens = row_as_list.count('-')
        if (length_of_block-number_of_hyphens)/length_of_block >= self.conservation_cutoff:
            return True
        else:
            return False

    def is_block_conserved(self, index):
        if (index < len(self.column_categories)):
            if (self.column_categories[index] == 0):
                return True
            else:
                return False
        else:
            return (-1)

    def make_column_category_list(self):
        number_of_columns = len(self.sequence_array[0,:].tolist())
        for i in range(number_of_columns):
            if(self.is_conserved(self.sequence_array[:,i]) == True):
                self.column_categories.append(0) # 0 is for a conserved columns
            else:
                self.column_categories.append(1) # 1 is for an insert column
        return number_of_columns

    def get_column_or_block(self, index):
        number_of_columns = len(self.sequence_array[0, :].tolist())
        if index >= number_of_columns:
            return
        if len(self.column_categories[index]) == 0:  # conserved, return only column
            return (self.sequence_array[:,index].tolist(), index+1)
        else:
            # how many non_conserved columns do we have
            end_of_block = 1
            while (index + end_of_block < number_of_columns) and (self.column_categories[index + end_of_block] == 1):
                end_of_block += 1
            return (self.sequence_array[:,index:index+end_of_block].tolist(), index + end_of_block)

    def is_single_column(self, array):
        if len(array[0] == 1):
            return True
        else:
            return False

    def identify_blocks(self):
        start = 0
        current_column = start
        length = len(self.sequence_array[0,:].tolist())
        while current_column < length:
            while (current_column < length) and (self.is_conserved(self.sequence_array[:,current_column])):
                current_column += 1
            if (start < length) and (current_column <= length) and start != current_column:
                self.sequence_blocks[0].append((start,current_column))
                #print(start,current_column)
                self.order_of_blocks.append(0)
                start = current_column
            while (current_column < length) and not(self.is_conserved(self.sequence_array[:,current_column])):
                current_column += 1
            if (start < length) and (current_column <= length):
                self.sequence_blocks[1].append((start,current_column))
                #print(start, current_column)
                self.order_of_blocks.append(1)
                start = current_column
        # this gives the index into self.sequence_blocks for each entry (0 or 1) in order_of_blocks
        null_index = 0
        one_index = 0
        for i in range(len(self.order_of_blocks)):
            if self.order_of_blocks[i] == 0:
                self.indexes_of_block_order.append(null_index)
                null_index += 1
            else:
                self.indexes_of_block_order.append(one_index)
                one_index += 1

        return self.order_of_blocks, self.indexes_of_block_order


    def number_of_conserved_nodes(self):
        sum=0
        length = len(self.sequence_blocks[0])
        for i in range(length):
            sum += self.sequence_blocks[0][i][1]-self.sequence_blocks[0][i][0]
        return(sum)

    def map_symbols_to_MID(self, row):
        self.current_M_column.clear()
        row_as_list = row.tolist()
        for i in range(len(row_as_list)):
            if row_as_list[i] == '-':
                self.current_M_column[i] = 'D'
            else:
                self.current_M_column[i] = 'M'
        return self.current_M_column

    # calculate the number of self jumps and jumps to next conserved block for insert block
    def number_of_jumps(self, array):
        pass
        return

    def GetTotalNumberOfConservedBocks(self):
        number_of_conserved_blocks = len(self.order_of_blocks) - self.order_of_blocks.count(1)
        return number_of_conserved_blocks

    def IsNextBlockConserved(self, current_index):
        if (current_index + 1) < self.GetTotalNumberOfConservedBocks():
            if(self.order_of_blocks[current_index+1] == 0):
                return True
            else:
                return False
        return

    def insert_block_interpret(self, array):
        number_of_pass_throughs = 0
        for i in range(len(self.next_M_column)):
            if len(self.next_M_column[i]) == 0:
                number_of_pass_throughs += 1
        number_of_inserts = len(self.next_M_column) - number_of_pass_throughs
        number_of_return_loops = 0
        for i in range(len(self.next_M_column)):
            if len(self.next_M_column[i]) > 0:
                number_of_return_loops += (len(self.next_M_column[i]) - 1)

        return (number_of_pass_throughs, number_of_inserts, number_of_return_loops)

    def number_of_deletions_in_M_block(self, array):
        number_of_deletions = array.count('-')
        return number_of_deletions

    def get_indexes_of_deletions_in_M_block(self, array):
        return [index for index in range(len(array)) if array[index] == '-']

    def number_of_rows(self):
        return(len(self.sequence_array.tolist()))

    def list_of_insert_to_deletion(self, current_M_block, new_M_block):
        pass
        return

    # make iD list where any row with insertion becomes 'I', all else 'P' for pass-through
    def collapse_insert_block_to_list(self, array):
        insert_list = []
        number_of_columns = len(a[0])
        for i in range(len(a)):
            if a[i].count('-') < number_of_columns:
                insert_list.append('I')
            else:
                insert_list.append('P')
        return insert_list



    def step_through(self):
        self.D_column = [0]*self.number_of_rows()
        self.transition_paths = [[0] * self.number_of_paths for i in range(len(self.sequence_array))]
        current_index = 0
        self.current_M_column = ['M']*self.number_of_rows()
        while current_index < self.column_categories:
            self.next_M_column, next_index = self.get_next_column_or_block(current_index)
            if (self.is_block_conserved(current_index)):

            else:
                number_of_pass_throughs, number_of_inserts, number_of_return_loops = self.insert_block_interpret(self.next_M_column)
                # dictionary_of_transition_paths = {'MM':0,'MI':1,'MD':2,'IM':3,'II':4,'ID':5,'DM':6,'DI':7,'DD':8}

                self.current_M_column = self.next_M_column
                current_index = next_index
                self.next_M_column, next_index = self.get_column_or_block(current_index)

                number_of_rows = self.number_of_rows()
                insert_list = self.collapse_insert_block_to_list(self.current_M_column)

                # any returns from deletions to insertions:
                number_of_deletions_to_insertions = 0
                number_of_insertions_to_deletions = 0
                number_of_passthrough_to_deletions = 0
                number_of_current_active_deletions = len(self.D_column)
                if len(self.D_column > 0):
                    for i in range(len(insert_list)):
                        if insert_list[i] == 'I' and self.D_column[i] == 'D':
                            number_of_deletions_to_insertions += 1
                            self.D_column[i] = ''

                # any insertions to deletions?
                if self.number_of_deletions_in_M_block(self.next_M_column) > 0:
                    for i in range(len(insert_list)):
                        if insert_list[i] == 'I' and self.next_M_column[i] == '-':
                            number_of_insertions_to_deletions += 1
                            self.D_column[i] = 'D'
                        if insert_list[i] == 'P' and self.next_M_column[i] == '-':
                            number_of_passthrough_to_deletions += 1
                            self.D_column[i] = 'D'

                # Calculate fractions
                number_of_pass_throughs, number_of_inserts, number_of_return_loops

                if current_index == 0: # we are at the start node. No previous deletions

                    self.transition_paths[current_index][
                        self.dictionary_of_transition_paths('MI')] = number_of_inserts / len(self.next_M_column)
                    self.transition_paths[current_index][
                        self.dictionary_of_transition_paths('II')] = number_of_return_loops / (
                                number_of_inserts + number_of_return_loops)
                    self.transition_paths[current_index][
                        self.dictionary_of_transition_paths('IM')] = (number_of_inserts-number_of_insertions_to_deletions) / (
                                number_of_inserts + number_of_return_loops)
                    self.transition_paths[current_index][
                        self.dictionary_of_transition_paths('ID')] = number_of_insertions_to_deletions / (
                                number_of_inserts + number_of_return_loops)
                    self.transition_paths[current_index][
                        self.dictionary_of_transition_paths('MI')] = number_of_inserts / len(self.next_M_column)
                    self.transition_paths[current_index][
                        self.dictionary_of_transition_paths('MD')] = number_of_passthrough_to_deletions / len(self.next_M_column)
                    self.transition_paths[current_index][
                        self.dictionary_of_transition_paths('MM')] = (len(self.next_M_column) - number_of_passthrough_to_deletions - number_of_inserts) / len(self.next_M_column)














    # make a list of tuples [(0),(1),(2,4),(4),...,(7)]

    # def GetConservedBlock(self, index):
    #     if index < self.GetTotalNumberOfConservedBocks():
    #         bl
    #         column =
    #         self.sequence_array



    def step_through_sequence_blocks(self):
        total_number_of_conserved_blocks = self.GetTotalNumberOfConservedBocks()

        # current_block = 0
        # current_M_node = 0
        # self.path_probabilities = [[]*(self.number_of_conserved_nodes()+2)]
        # path = [[0]*self.number_of_paths]
        # while current_block < len(self.order_of_blocks):
        #     if self.order_of_blocks[current_block] == 1: # is this an insertion block
        #         start, end = self.sequence_blocks[self.indexes_of_block_order[current_block]]
        #         sequence_row = self.sequence_array[:, start]
        #         current_state = self.map_symbols_to_MID(sequence_row)
        #         # what is the symbols in the next conserved column if there is one
        #         M_to_D = 0
        #         if current_block < len(self.order_of_blocks-1):
        #             start, end = self.sequence_blocks[self.indexes_of_block_order[current_block+1]]
        #             sequence_row = self.sequence_array[:, start]
        #             new_state = self.map_symbols_to_MID(sequence_row)
        #             for i in range(len(sequence_row)):
        #                 if current_state[i] == [] and new_state[i] == 'D':
        #                     M_to_D += 1
        #         fraction_M_to_D = M_to_D/len(new_state)
        #
        #         number_of_deletions =0
        #         # count the number of deletions that occurs immediately, without insertion
        #         sequence_row = self.sequence_array[:, start]
        #
        #
        #     else: # this is a conserved block
        #         start,end = self.sequence_blocks[self.indexes_of_block_order[current_block]]
        #         sequence_row = self.sequence_array[:,start]
        #         current_state = self.map_symbols_to_MID(sequence_row)
        #         number_of_Ms = current_state.count('M')
        #         number_of_Ds = current_state.count('D')
        #         total = len(current_state)
        #         #call pseudocounts here if you want to use them
        #         path[self.dictionary_of_transition_paths('MM')] = number_of_Ms/total
        #         path[self.dictionary_of_transition_paths('MD')] = number_of_Ms / total
        #         self.path_probabilities[current_M_node].append(path)
        #         current_block += 1
        #         current_M_node += 1
        #
        #
        #
        # for i in range(len(self.order_of_blocks)):
        #     block_start_and_end = self.sequence_blocks[self.order_of_blocks[i]][self.indexes_of_block_order[i]]
        #     print(block_start_and_end)
        return



#HMM = HMM(3,'N')
#HMM.list_nodes()
# hmm_sequences = HMM_sequences()
# print(hmm_sequences.sequence_array)
# hmm_sequences.identify_blocks()
# print(hmm_sequences.sequence_blocks)
# print(hmm_sequences.number_of_conserved_nodes())
# print(HMM.number_of_nodes)
# print(HMM.array_of_nodes)
#HMM.make_end_node()
# hmm_draw = HMM_draw()
# hmm_draw.draw()
# hmm_draw.draw_figure_title('My HMM')
# #hmm_draw.draw_arrow(1,1,4,4)
# hmm_draw.show_drawing()
hmm_sequence = HMM_sequence()
print(hmm_sequence.sequence_array)
hmm_sequence.make_column_category_list()
current_index = 0
while current_index < len(hmm_sequence.column_categories):
    block, new_index = hmm_sequence.get_next_column_or_block(current_index)
    print(current_index, block, hmm_sequence.remove_dashes_from_array(block))
    current_index = new_index
# order_of_blocks, indexes_of_block_order = hmm_sequence.identify_blocks()
# print(order_of_blocks, hmm_sequence.sequence_blocks, indexes_of_block_order)
# #hmm_sequence.step_through_sequence_blocks()
# my_array = np.array([['-','-','-'],['G','-','T'],['G','G','G']])
#
# new_array = hmm_sequence.remove_dashes_from_array(my_array)
# print(new_array)


# result = hmm_sequence.get_emmission_list(['G','G','A','T','C','C','C',], 'D')
