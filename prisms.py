#################################################################################
# 
# Code accompanying Chapter 4: Constructing Infinite Series, of my essay 
# 
# Crystallisations and Combinatorics for Many Triangulated 3-Manifolds
#
# which was written in partial fulfillment of the requirements for the degree of
# Bachelor of Science (Advanced Mathematics) (Honours)
# at the University of Sydney
#
# Max Tobin, 2020
# 
#################################################################################

################################################################################## 
#
# Notes: 
#
# Colours are always given as -1 0, and 1 and colour 0 always represents the 
# colour fixed under 180 degree rotation.
# Hence 180 degree rotation performs the colour swap c -> -c or (-1 1).
# 
# Both degrees and colours are given as sequences around the boundary of the 
# square, of length half the number of vertices on the boundary and using the 
# same direction and starting point when given in pairs.
#
#################################################################################

import math
import copy
import itertools
import random

def is_valid_prism(seq_pair):
    """
    Check if a pair of degree and colour sequences for a fundamental square boundary satisfy the conditions for a (1^1,3^k)-type 2-sphere.
    Returns as a boolean.
    
    seq_pair: pair (tuple) of lists of integers, the first giving degrees (positive integers) and the second giving colours (-1, 0 or 1) 
    """
    deg = seq_pair[0]
    colours = seq_pair[1]
    l = len(deg)
    for i in range(-1,l-1):
        if deg[i]+deg[i+1] == 5:
            valid = True
            for j in range(1,int(math.ceil(l/2-1))+1):
                if deg[(i-j)%l]+deg[(i+1+j)%l] != 6:
                    valid = False
                    break
            if valid:
                realigned_colours = colours.copy()
                if i != -1:
                    for j in range(l-i-1):
                        c = -1*realigned_colours.pop()
                        realigned_colours.insert(0,c)
                other_side_colours = [-1*c for c in realigned_colours]
                other_side_colours.reverse()
                if other_side_colours == [c%3-1 for c in realigned_colours] or other_side_colours == [(c-1)%3-1 for c in realigned_colours]:
                    deg[:] = deg[i+1:]+deg[:i+1]
                    colours[:] = realigned_colours
                    return True
    return False

def add_triangles(seq_pair,times):
    """
    Compute all sequence pairs obtainable from a given one by adding opposite triangles a given number of times, while only creating
        new interior vertices of degree 6.
    Returns as a list of sequence pairs (tuples of lists of integers).

    seq_pair: pair (tuple) of lists of integers, the first giving degrees (positive integers) and the second giving colours (-1, 0 or 1)
    times: positive integer
    """
    deg = seq_pair[0]
    l = len(deg)
    colours = seq_pair[1]
    sequence_pairs = []
    seen = set()
    if times == 1:
        v = deg[-2]
        after = deg[-1]
        for i in range(len(deg)):
            before = v
            v = after
            after = deg[i]
            if after != 5:
                if v == 5:
                    valid_add = False
                    if before != 5:
                        if i <= 1:
                            if colours[i-2] != (-1)*colours[i]:
                                valid_add = True
                                new_deg = [after+1]+deg[i+1:i-2]+[before+1]
                                new_col = colours[i:l-1+i]
                        else:
                            if colours[i-2] != colours[i]:
                                valid_add = True
                                new_deg = deg[:i-2]+[before+1,after+1]+deg[i+1:]
                                new_col = colours[:i-1]+colours[i:]
                    if valid_add:
                        new_seq_tup = (tuple(new_deg),tuple(new_col))
                        if new_seq_tup not in seen:
                            seen.add(new_seq_tup)
                            sequence_pairs.append((new_deg,new_col))
                else:
                    if i == 0:
                        new_deg = [after+1]+deg[1:-1]+[v+1,1]
                        new_col = colours+[colours[0]+(-1)*colours[-1]] 
                            # new vertex is different colour to both its neighbours, noting 1 <-> -1 when wrapping around end of list
                    else:
                        new_deg = deg[:i-1]+[v+1,1,after+1]+deg[i+1:]
                        new_col = colours[:i]+[(-1)*(colours[i-1]+colours[i])]+colours[i:]
                            # new vertex is different colour to both its neighbours, no wrapping needed
                    new_seq_tup = (tuple(new_deg),tuple(new_col))
                    if new_seq_tup not in seen:
                            seen.add(new_seq_tup)
                            sequence_pairs.append((new_deg,new_col))
        return sequence_pairs
    elif times > 1:
        for prev_seq_pair in add_triangles(seq_pair,times-1):
            for new_seq_pair in add_triangles(prev_seq_pair,1):
                new_seq_tup = (tuple(new_seq_pair[0]),tuple(new_seq_pair[1]))
                if new_seq_tup not in seen:
                    seen.add(new_seq_tup)
                    sequence_pairs.append(new_seq_pair)
        return sequence_pairs

def extend_prisms(seq_pairs,times):
    """
    Compute all degree sequences of fundamental squares defining (1^1,3^k)-type 2-spheres obtainable from a given list of sequence 
        pairs by adding opposite triangles a given number of times.
    Assumes the original sequence is itself from a (1^1,3^k)-type 2-sphere fundamental square.
    Returns as a list of sequence pairs (tuples of lists of integers).

    seq_pairs: list of pairs (tuples) of lists of integers, the first giving degrees (positive integers) and the second giving colours (-1, 0 or 1)
    times: positive integer
    """
    valids = []
    seen = set()
    extensions = []
    for seq_pair in seq_pairs:
        extensions.extend(add_triangles(seq_pair,times))
    for ex in extensions:
        if is_valid_prism(ex):
            ex_tup = (tuple(ex[0]),tuple(ex[1]))
            if ex_tup not in seen:
                seen.add(ex_tup)
                valids.append(ex)
    return valids

def prism_search(seq_pair,n,times):
    """
    Attempt to extend a given (1^1,3^k)-type 2-sphere fundamental square (with n=3k+1) to successively larger ones by applying random 
        valid extensions.
    Extension will be attempted the given number of times.
    An attempt involves computing extend_prisms() and if it contains at least one element (success) choosing one at random.
        If the previous attempt succeeded, attempt to extend the result by adding triangles once.
        If the previous attempt failed, attempt to extend the same base by adding triangles one additional time.

    seq_pair: pair (tuple) of lists of integers, the first giving degrees (positive integers) and the second giving colours (-1, 0 or 1)
    n: positive integer, giving permutation length or half number of tetrahedra for original sphere
    times: positive integer
    """
    i = 0
    j = 1
    print("n = {} : degrees {}, colours: {}".format(n,seq_pair[0],seq_pair[1]))
    while i < times:
        extensions = extend_prisms([seq_pair],j)
        if extensions:
            seq_pair = random.choice(extensions)
            n = n+3*j
            print("n = {} : degrees {}, colours: {}".format(n,seq_pair[0],seq_pair[1]))
            i = i+1
            j = 1
        else:
            print("    extending {} time(s) failed".format(j))
            j = j+1
            i = i+1