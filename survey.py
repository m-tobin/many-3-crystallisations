#################################################################################
# 
# Code accompanying Chapter 3: Survey of $(1^1,3^k)$-Type Gems, of my essay 
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

import math
import copy
import itertools
import random

##################################################
# Basic operations on permutations
##################################################

def oneline_to_cycle(oneline):
    """
    Convert a permutation in one-line notation to cycle notation.
    Returns as a list of lists (representing cycles) sorted by increasing length.

    oneline: list of integers
    """
    cycles = []
    remaining = list(range(len(oneline)))
    while len(remaining) != 0:
        first_ind = remaining.pop(0)
        entry = oneline[first_ind]
        if first_ind == entry:
            cycles.append([first_ind])
        else:
            cycle = [first_ind,entry]
            ind = entry
            remaining.remove(ind)
            entry = oneline[ind]
            while entry != first_ind:
                cycle.append(entry)
                ind = entry
                remaining.remove(ind)
                entry = oneline[ind]
            cycles.append(cycle)
    cycles.sort(key=len)
    return cycles

def cycle_to_oneline(cycles,n):
    """
    Convert a permutation in cycle notation to one-line notation.
    Returns as a list of integers.

    cycles: list of lists of integers
    n: integer (length of permutation)
    """
    oneline = list(range(n))
    for cycle in cycles:
        k = len(cycle)
        for i in range(k-1):
            oneline[cycle[i]] = cycle[i+1]
        oneline[cycle[k-1]] = cycle[0]
    return oneline

def inverse(oneline):
    """
    Find the inverse of a permutation.
    Returns as a permutation in one-line notation (list of integers).

    oneline: permutation in one-line notation (list of integers)
    """
    inv = list(range(len(oneline)))
    for i in range(len(oneline)):
        inv[oneline[i]] = i
    return inv

def compose(oneline1,oneline2):
    """
    Calculate the composition (using "right-to-left" function composition order) of two permutations.
    Returns as a permutation in one-line notation (list of integers).

    oneline1, oneline2: permutations in one-line notation (list of integers)
    """
    new = []
    for entry in oneline2:
        new.append(oneline1[entry])
    return new

def relabel(oneline,new_labels):
    """
    Find the image of a permutation after relabelling letters {0,1,...,n-1} according to another given permutation.
    Returns as a permutation in one-line notation (list of integers).

    oneline, new_labels: permutations in one-line notation (list of integers)    
    """
    new = []
    for i in range(len(oneline)):
        old_label_i = new_labels.index(i)
        old_image = oneline[old_label_i]
        image = new_labels[old_image]
        new.append(image)
    return new


##################################################
# Generating permutations
##################################################

def recurse_cycle_type_perms(remaining):
    """
    Find all partitions of a list into (ordered) 3-cycles.
    Recursive step of given_structure_perms().
    Returns as a list containing lists of length 3.

    remaining: list with length a multiple of 3.
    """
    first = remaining[0]
    if len(remaining) == 3:
        yield [[first,remaining[1],remaining[2]]]
        yield [[first,remaining[2],remaining[1]]]
    else:
        remaining.remove(first)
        for i in remaining:
            for j in remaining:
                if i==j:
                    break
                else:
                    new_remaining = remaining.copy()
                    new_remaining.remove(i)
                    new_remaining.remove(j)
                    cycle_1 = [[first,j,i]]
                    cycle_2 = [[first,i,j]]
                    for tail in recurse_cycle_type_perms(new_remaining):
                        yield cycle_1 + tail
                        yield cycle_2 + tail

def cycle_type_perms(n):
    """
    Construct all permutations of length n with cycle type (1^1,3^{(n-1)/3}).
    Returns as a generator object for a list of permutations in cycle notation (lists of lists of integers).

    n: positive integer of the form 3k+1.
    """
    for fixed_point in range(n):
        remaining = list(range(n))
        remaining.remove(fixed_point)
        for tail in recurse_cycle_type_perms(remaining):
            yield [[fixed_point]]+tail

def cycle_type_perms_restricted(n):
    """
    Construct all permutations of length n with cycle type (1^1,3^{(n-1)/3}) containing cycles (1) and (0 4 7).
    Returns as a generator object for a list of permutations in cycle notation (lists of lists of integers).

    n: integer greater than or equal to 10 of the form 3k+1
    """
    remaining = [2,3,5,6]+list(range(8,n))
    for tail in recurse_cycle_type_perms(remaining):
        yield [[1],[0,4,7]]+tail


##################################################
# Searching for 2-spheres
##################################################

def mu_oneline(n):
    """
    Construct the permutation mu_n := (0)(1 2 3)(4 5 6)...(n-3 n-2 n-1) in oneline notation.
    Returns as a list of integers (permutation in one-line notation).

    n: positive integer of the form n=3k+1
    """
    mu = [0]
    for i in range(int((n-1)/3)):
        mu.append(3*i+2)
        mu.append(3*i+3)
        mu.append(3*i+1)
    return mu

def perm_mu_inv(oneline):
    """
    Calculate sigma * mu^-1, where sigma is a given permutation and mu = (0)(1 2 3)(4 5 6)...(n-3 n-2 n-1).
    Returns as a list of integers (permutation in one-line notation).

    Note: in priciple we will want the cycle type of mu^-1 sigma. This has the same cycle type but is
        quicker to compute.

    oneline: list of integers
    """
    new = [oneline[0]]
    n_three_cycles = int((len(oneline)-1)/3)
    for i in range(n_three_cycles):
        k = 3*i + 1
        new.append(oneline[k+2])
        new.append(oneline[k])
        new.append(oneline[k+1])
    return new

def has_correct_cycle_structure(cycles):
    """
    Check if a permutation has the cycle type (1^1,3^{(n-1)/3}).
    Returns as a boolean.

    cycles: permutation in cycle notation (list of lists of integers)
    """
    if len(cycles[0]) != 1:
        return False
    else:
        for cycle in cycles[1:]:
            if len(cycle) != 3:
                return False
    return True

def component_graph(oneline):
    """
    Construct a graph from a potential gem of the specified form indicating how 2-coloured edges join the (0,1)-connected components.
    Vertices are connected components of the graph specified by id and (0)(1 2 3)...(n-3 n-2 n-1), with an edge if an edge specified by
        the given permutation joins them.
    Returns as a dictinary, with keys 0,1,...,(n-1)/3+1 representing connected components {0,0'}, {1,1',2,2',3,3'}, ...

    oneline: permutation (specifying third colour) in one-line notation (list of integers)
    """
    n = len(oneline)
    num_components = int((n-1)/3)+1
    graph = {comp:[] for comp in range(num_components)}
    for i in range(n):
        i_comp = math.ceil(i/3)
        image_comp = math.ceil(oneline[i]/3)
        if i_comp != image_comp:
            if image_comp not in graph[i_comp]:
                graph[i_comp].append(image_comp)
            if i_comp not in graph[image_comp]:
                graph[image_comp].append(i_comp)
    return graph

def is_connected(graph):
    """
    Check if a graph in dictionary format is connected.
    Returns as a Boolean.

    graph: dictionary with keys representing vertices and values which are lists of all adjacent vertices
    """
    found = [0]
    unchecked = graph[0].copy()
    found.extend(unchecked)
    while len(found) < len(graph):
        new = []
        for v in unchecked:
            for image in graph[v]:
                if image not in found and image not in new:
                    new.append(image)
        if new == []:
            return False
        else:
            found.extend(new)
            unchecked = new
    return True

def family_two_spheres(n):
    """
    Find all 2-gems (2-spheres) specified by permutations id, mu = (0)(1 2 3)...(n-3 n-2 n-1) and a third permutation sigma
        such that sigma and mu^-1 * sigma both have cycle type (1^1,3^{(n-1)/3}).
    Returns as a list of permutations representing sigma in cycle notation (lists of lists of integers).

    n: positive integer of the form 3k+1
    """
    family = []
    for perm in cycle_type_perms(n):
        perm_oneline = cycle_to_oneline(perm,n)
        muinvp_cycles = oneline_to_cycle(perm_mu_inv(perm_oneline))
        if has_correct_cycle_structure(muinvp_cycles):
            if is_connected(component_graph(perm_oneline)):
                family.append(perm)
    return family

def family_two_spheres_restricted(n):
    """
    Find all 2-gems (2-spheres) specified by permutations id, mu = (0)(1 2 3)...(n-3 n-2 n-1) and a third permutation sigma
        containing the cycles (1) and (0 4 7) such that sigma and mu^-1 * sigma both have cycle type (1^1,3^{(n-1)/3}).
    Returns as a list of permutations representing sigma in cycle notation (lists of lists of integers).

    n: positive integer of the form 3k+1 greater than or equal to 10
    """
    family = []
    for perm in cycle_type_perms_restricted(n):
        perm_oneline = cycle_to_oneline(perm,n)
        muinvp_cycles = oneline_to_cycle(perm_mu_inv(perm_oneline))
        if has_correct_cycle_structure(muinvp_cycles):
            if is_connected(component_graph(perm_oneline)):
                family.append(perm)
    return family


##################################################
# Isomorphisms of 2-spheres
##################################################

def reflect(oneline):
    """
    Find the image of a permutation under the standard reflection isomorphism.
    Returns as a permutation in one-line notation (list of integers).

    oneline: permutation in one-line notation (list of integers)  
    """
    new_labels = [0]
    for i in range(int((len(oneline)-1)/3)):
        new_labels.append(3*i+1)
        new_labels.append(3*i+3)
        new_labels.append(3*i+2)
    return inverse(relabel(oneline,new_labels))

def rotations(num_blocks):
    """
    List all relabellings of {0,1...,3*num_blocks} corresponding to cycling elements of blocks (1,2,3), (4,5,6), (7,8,9), ... etc.
    Returns as a generator object.

    num_blocks: positive integer
    """
    if num_blocks == 0:
        yield [0]
    else:
        z = 3*num_blocks
        y = z-1
        x = y-1
        for l in rotations(num_blocks-1):
            yield l+[x,y,z]
            yield l+[y,z,x]
            yield l+[z,x,y]

def block_permutations(num_blocks):
    """
    List all relabellings of {0,1...,3*num_blocks} corresponding to permuting the blocks (1,2,3), (4,5,6), (7,8,9), ... etc.
    Returns as a generator object.

    num_blocks: positive integer 
    """
    blocks = [[3*i+1,3*i+2,3*i+3] for i in range(num_blocks)]
    for perm in itertools.permutations(range(num_blocks)):
        relabelling = [0]
        for i in perm:
            relabelling.extend(blocks[i])
        yield relabelling

def isomorphism_class(sphere,n,oneline=False,reduce=False,index=False):
    """
    Find all spheres isomorphic to a given sphere by vertex-relabelling isomorphisms.
    Returns as a list of the permutations which specify the spheres in cycle notation (lists of lists of integers).

    sphere: permutation in cycle notation specifying a sphere (list of lists of integers)
    """
    num_blocks = int((n-1)/3)
    if oneline:
        sphere_oneline = sphere
    else:
        sphere_oneline = cycle_to_oneline(sphere,n)
    reflections = [sphere_oneline,reflect(sphere_oneline)]
    rots = list(rotations(num_blocks))
    block_perms = list(block_permutations(num_blocks))
    iso_class = []
    for i in range(len(reflections)):
        for j in range(len(rots)):
            rotated = relabel(reflections[i],rots[j])
            for k in range(len(block_perms)):
                new_sphere = relabel(rotated,block_perms[k])
                if not reduce or new_sphere not in iso_class:
                    if index:
                        iso_class.append((new_sphere,(i,j,k)))
                    else:
                        iso_class.append(new_sphere)
    if index:
        if oneline:
            return iso_class
        else:
            return [(oneline_to_cycle(pair[0]),pair[1]) for pair in iso_class]
    else:
        if oneline:
            return iso_class
        else:
            return [oneline_to_cycle(sphere) for sphere in iso_class]


##################################################
# Searching for 3-crystallisations
##################################################

def three_colour_graphs(oneline1,oneline2,n):
    """
    Construct two graphs describing connectivity after removing colour 1 or colour 0 from a potential 3-crystallisation.
    First graph is produced by removing 1-coloured edges and joining pairs (i,i').
    Second graph is produced by removing 0-coloured edges.
    Returns as a pair (tuple) of dictionaries representing graphs, with keys representing vertices and values lists of vertices joined to it by an edge.
        First has keys 0,1,...,n-1.
        Second has keys 0,1,...,2n-1 reperesenting 0,1,...,n-1,0',1',...,(n-1)' respectively.

    oneline1, oneline2: permutations in one-line notation which together define a 3-crystallisation (lists of integers)
    n: length of permutations oneline1, oneline2 (positive integer)
    """
    mu = mu_oneline(n)
    graph_023 = {i:[] for i in range(n)}
    graph_123 = {i:[] for i in range(2*n)}
    for i in range(n):
        image_mu = mu[i]
        image1 = oneline1[i]
        image2 = oneline2[i]
        graph_023[i].append(image1)
        graph_023[i].append(image2)
        graph_023[image1].append(i)
        graph_023[image2].append(i)
        graph_123[i].append(image_mu+n)
        graph_123[i].append(image1+n)
        graph_123[i].append(image2+n)
        graph_123[image_mu+n].append(i)
        graph_123[image1+n].append(i)
        graph_123[image2+n].append(i)
    return (graph_023,graph_123)

def is_contracted(oneline1,oneline2,n):
    """
    Check if the two remaining subgraphs to be checked for a potential 3-crystallisation are a connected.
    Assumes it has already been confirmed that both defining permutations individually define (connected) 2-sphere gems.

    Note the terminology ''contracted'' means that removing any colour gives a connected graph.

    oneline1, oneline2: permutations in one-line notation which together define a 3-crystallisation (lists of integers)
    n: length of permutations oneline1, oneline2 (positive integer)
    """
    (graph_023,graph_123) = three_colour_graphs(oneline1,oneline2,n)
    return is_connected(graph_023) and is_connected(graph_123)

def family_three_crystallisations(n, sphere1=[], spheres=[], return_product=False):
    """
    Find all 3-crystallisations specified by permutations id and three permutations of cycle type (1^1,3^{(n-1)/3}): 
        mu = (0)(1 2 3)...(n-3 n-2 n-1)
        a fixed permutation sphere1 defining a 2-sphere gem
        a permutation from a given list of spheres (permutations defining 2-sphere gems). 
    If sphere1 is not specified, the first element of spheres is used.
    If spheres is not specified, the full list of 2-spheres family_two_spheres(n) is used.
    Returns as a list of pairs (tuples) of permutations in cycle notation (lists of lists of integers).
    
    n: positive integer of the form 3k+1
    sphere1: permutation in cycle notation
    spheres: list of permutations in cycle notation
    return_product: if true, returns instead as a list of triples (sigma1, sigma2, sigma1^-1 * sigma2)
    """
    if spheres == []:
        spheres = family_two_spheres(n)
    if sphere1 == []:
        sphere1 = spheres[0]
    gems = []
    oneline1 = cycle_to_oneline(sphere1,n)
    for sphere2 in spheres:
        oneline2 = cycle_to_oneline(sphere2,n)
        s1_inv_s2_cycles = oneline_to_cycle(compose(inverse(oneline1),oneline2))
        if has_correct_cycle_structure(s1_inv_s2_cycles) and is_contracted(oneline1,oneline2,n):
            if return_product:
                gems.append((sphere1,sphere2,s1_inv_s2_cycles))
            else:
                gems.append((sphere1,sphere2))
    return gems


##################################################
# Colour swaps for 3-crystallisations
##################################################

def colour_swap(gem,swap):
    """
    Apply a colour swap to a 3-crystallisation specified by id, mu = (0)(1 2 3)(4 5 6)... and two given permutations.
    The swap is specified by a permutation on {0,1,2,3} given in oneline notation.
    Returns as a pair (tuple) of permutations in cycle notation.

    Procedure is to relabel the colours as specified, then relabel the vertices to such that the first 3 colours,
        or 2 if this is not possible, are the same permutations as for the original gem.
    The vertex relabelling is achieved by:
        First, relabeling only bottom-row vertices such that colour 0 becomes id.
        Second, relabelling colour-0 pairs together such that the cycles of colour 1 (ordered lexicographically in cycle notation)
            become the corresponding cycles of mu.
        Third, finding the first 2-sphere isomorphism (as ordered by isomorphism_class()) mapping colour 2 to the given
            original colour 2 permutation (if one exists) and applying it to the gem.

    gem: pair of permutations in cycle notation representing colours 2 and 3 respectively (tuple of lists of lists)
    swap: permutation on {0,1,2,3} in oneline notation (list) 
    """
    n = len(gem[0])*3-2
    mu = mu_oneline(n)
    sigma = cycle_to_oneline(gem[0],n)
    new_perms = [[],[],[],[]]
    new_perms[swap[0]] = list(range(n))
    new_perms[swap[1]] = mu.copy()
    new_perms[swap[2]] = sigma.copy()
    new_perms[swap[3]] = cycle_to_oneline(gem[1],n)
    # realign first colour to id
    if swap[0] != 0:
        inv = inverse(new_perms[0])
        new_perms[0] = list(range(n))
        for i in range(1,4):
            new_perms[i] = compose(new_perms[i],inv)
    # realign second colour to mu
    if new_perms[1] != mu:
        cycles = oneline_to_cycle(new_perms[1])
        new_labels = list(range(n))
        new_labels[cycles[0][0]] = 0
        for i in range(1,len(cycles)):
            new_labels[cycles[i][0]] = 3*i-2
            new_labels[cycles[i][1]] = 3*i-1
            new_labels[cycles[i][2]] = 3*i
        new_perms[1] = mu.copy()
        for i in range(2,4):
            new_perms[i] = relabel(new_perms[i],new_labels)
    # realign third colour to original sigma
    if new_perms[2] != sigma:
        iso_class = isomorphism_class(new_perms[2],n,oneline=True,index=True)
        for pair in iso_class:
            if pair[0] == sigma:
                iso_index = pair[1]
                num_blocks = int((n-1)/3)
                rotations_gen = itertools.islice(rotations(num_blocks),iso_index[1],None)
                rot = rotations_gen.__next__()
                block_perms_gen = itertools.islice(block_permutations(num_blocks),iso_index[2],None)
                block_perm = block_perms_gen.__next__()
                new_perms[2] = sigma
                if iso_index[0] == 0:
                    new_perms[3] = relabel(relabel(new_perms[3],rot),block_perm)
                else:
                    new_perms[3] = relabel(relabel(reflect(new_perms[3]),rot),block_perm)
                break
    return (oneline_to_cycle(new_perms[2]),oneline_to_cycle(new_perms[3]))

def colour_swap_class(gem):
    """
    Compute the 3-crystallisations obtained from the given one for each of the 24 colour swaps and list them with all corresponding swaps.
    Returns as a list of pairs (tuples) where:
        the first is a 3-crystallisation given as a pair (tuple) of permutations in cycle notation.
        the second is a list of swaps given as 4-tuples representing a permutation of {0,1,2,3}.

    gem: pair of permutations in cycle notation representing colours 2 and 3 respectively (tuple of lists of lists)
    """
    gem_class = []
    swaps = []
    for swap in itertools.permutations(range(4)):
        new_gem = colour_swap(gem,swap)
        try:
            i = gem_class.index(new_gem)
            swaps[i].append(swap)
        except ValueError:
            gem_class.append(new_gem)
            swaps.append([swap])
    return list(zip(gem_class,swaps))

def all_colour_swap_classes(gems):
    """
    Sort a list of 3-crystallisations into colour swap classes.
    Returns as a list where each element is the output of colour_swap_class() for one member of each colour swap class.

    gems: list of pairs of permutations in cycle notation representing colours 2 and 3 respectively (list of tuples of lists of lists)
    """
    remaining = copy.deepcopy(gems)
    classes = []
    while len(remaining) != 0:
        gem = remaining[0]
        gem_class = []
        swaps = []
        for swap in itertools.permutations(range(4)):
            new_gem = colour_swap(gem,swap)
            try:
                i = gem_class.index(new_gem)
                swaps[i].append(swap)
            except ValueError:
                gem_class.append(new_gem)
                swaps.append([swap])
                remaining.remove(new_gem)
        classes.append(list(zip(gem_class,swaps)))
    return classes


##################################################
# Fundamental groups
##################################################

def fundamental_group(triple,simplify=False):
    """
    Computes a finite presentation of the fundamental group of a 3-crystallisation specified by permutations
        id, mu = (0)(1 2 3)...(n-3 n-2 n-1) and two permutations sigma1, sigma2 with cycle type (1^1,3^{(n-1)/3}).
    Input format is the output format of family_three_crystallisations(n, return_product=True).
    Returns as a pair of lists (generators, relations).
        Generators are strings of the form "x0", "x1", "x2", ...
        Relations are lists of generators "xi" and their inverses "xi^-1".

    Uses the standard procedure for crystallisation fundamental groups (as in Ferri et. al. 1986) restricted to this special case.
    The generators represent (0,1)-coloured cycles, ordered by the corresponding cycles in mu and excluding the last.
        In particular "x0" is the cycle (0,0') of length 2, and all other cycles are of length 6.
    The relations are obtained from (2,3)-coloured cycles, ordered by the corresponding cycles in sigma1^-1 * sigma2 and excluding the last.
        In particular the first is from the cycle of length 2, so it has length at most 2 and all others have length at most 6.

    triple: 3-tuple (sigma1, sigma2, sigma1^-1 * sigma2) of permutations in cycle notation (lists of lists of integers).
    simplify: if true, cancel any occurences of "xi xi^-1" or "xi^-1 xi" from relations.
    """
    n = len(triple[0])*3-2
    num_three_cycles = int((n-1)/3)
    p2invp1_cycles = triple[2]
    p2_oneline = cycle_to_oneline(triple[1],n)
    generators = ["x{}".format(i) for i in range(num_three_cycles)]
    short_cycle_1 = p2invp1_cycles[0][0]
    short_cycle_1_gen = math.ceil(short_cycle_1/3)
    short_cycle_2 = p2_oneline[short_cycle_1]
    short_cycle_2_gen = math.ceil(short_cycle_2/3)
    short_relation = []
    if simplify and short_cycle_1_gen == short_cycle_2_gen:
        pass
    else:
        if short_cycle_1_gen < num_three_cycles:
            short_relation.append("x{}".format(short_cycle_1_gen))
        if short_cycle_2_gen < num_three_cycles:
            short_relation.append("x{}^-1".format(short_cycle_2_gen))
    relations = [short_relation]
    for i in range(1,num_three_cycles):
        relation = []
        for j in p2invp1_cycles[i]:
            colour2_gen = math.ceil(j/3)
            if colour2_gen < num_three_cycles:
                colour2_gen_string = "x{}".format(colour2_gen)
                if simplify:
                    if relation != [] and relation[-1] == colour2_gen_string + "^-1":
                        relation.pop()
                    else:
                        relation.append(colour2_gen_string)
                else:
                    relation.append(colour2_gen_string)
            colour3_gen = math.ceil(p2_oneline[j]/3)
            if colour3_gen < num_three_cycles:
                colour3_gen_string = "x{}".format(colour3_gen)
                if simplify:
                    if relation != [] and relation[-1] == colour3_gen_string:
                        relation.pop()
                    else:
                        relation.append(colour3_gen_string+"^-1")
                else:
                    relation.append(colour3_gen_string+"^-1")
        relations.append(relation)
    return (generators,relations)

def all_fundamental_groups(n, gems=[], simplify=False, reduce=False):
    """
    Compute a finite presentation of the fundamental group for all 3-crystallisations specified by permutations id, 
        mu = (0)(1 2 3)...(n-3 n-2 n-1) and two permutations with cycle type (1^1,3^{(n-1)/3}).
    Returns as a list of pairs (generators,relations) as in the output format of fundamental_group().

    n: positive integer of the form 3k+1
    gems: optionally provide already computed output of family_three_crystallisations(n, return_product=True)
    simplify: if True, cancel any occurences of "xi xi^-1" or "xi^-1 xi" from relations.
    reduce: if True, remove duplicates of the same presentation (after simplification if used with simplify)
    """
    if gems == []:
        gems = family_three_crystallisations(n, return_product=True)
    presentations = []
    if reduce:
        for gem in gems:
            fg = fundamental_group(gem,simplify)
            if fg not in presentations:
                presentations.append(fg)
    else:
        for gem in gems:
            presentations.append(fundamental_group(gem,simplify))
    return presentations

def gap_readable_relations(group):
    """
    Convert the relations of a fundamental group to a form readable by GAP.
    Input format is the output format of fundamental_group().
    Returns as a string.

    group: pair of lists containing generators (strings) and relations (lists of strings)
    """
    g = copy.deepcopy(group)
    for i in range(len(g[1])):
        for j in range(len(g[1][i])):
            g[1][i][j] = 'g.'+str(int(g[1][i][j][1])+1)+g[1][i][j][2:]
    rel_string_list = ['*'.join(rel) for rel in g[1]]
    rel_single_string = '['+', '.join(rel_string_list)+']'
    return rel_single_string