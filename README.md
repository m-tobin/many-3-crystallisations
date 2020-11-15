# Many 3-Manifold Crystallisations

## About

This code accompanies my essay:

<em>Crystallisations and Combinatorics for Many Triangulated 3-Manifolds</em> by Max Tobin, 2020

which I wrote as part of a Bachelor of Science (Advanced Mathematics) (Honours) degree at the University of Sydney.

It is hosted here primarily for the purpose of reproducibility for the enumerative search results presented in the essay.

## Files and Usage

The aforementioned essay, in particular Sections 3.2 and 4.2, describes the usage and purpose of this code.

### survey.py

Implementation of all procedures described in Section 3.2, by which all data in Table 3.1 and 
(in conjunction with GAP https://www.gap-system.org/ and Regina http://regina-normal.github.io/ as described) Table 3.2 and Appendix A was obtained.

Note in particular the functions:

```python
family_two_spheres()
```
for (1^1,3^k)-type 2-spheres with 2 colours fixed
```python
family_two_spheres_restricted()
```
for the restricted set as in Proposition 3.2.1
```python
isomorphism_class()
```
for isomorphism classes of (1^1,3^k)-type 2-spheres (possibly with duplicates)
```python
family_three_crystallisations()
```
for (1^1,3^k)-type 3-crystallisations with 3 colours fixed
```python
family_three_crystallisations_isos()
```
for 3-crystallisations from isomorphisms of a single base sphere, in particular for n = 25
```python
find_three_crystallisation()
```
for finding larger examples from a given base sphere, in particular for n = 28
```python
all_colour_swap_classes()
```
for sorting 3-crystallisations into colour swap classes, and
```python
fundamental_group()
```
for computing a finite presentation of the fundamental group.

### prisms.py

For extending fundamental squares by adding opposite triangles as in Section 4.2, 
in terms of degree and colour sequences. This is not required for reproducing any data,
but it is an effective search method which can be used to look for patterns, such as the 
beginnings of the first infinite series presented in the essay.

For an example of this, try:
```python
n = 13
square = ([3,2],[0,-1])
deg = [3,2]
print('n = 13: {}'.format(square))
for i in range(1,4):
    deg.insert(0,3)
    for pair in prisms.extend_prisms([square],2*i+3):
        if pair[0] == deg:
            square = pair
            break
    n = n+3*(2*i+3)
    print('n = {}: {}'.format(n,square))
```

### regina_converter.py

An short example script showing how 3-crystallisations as obtained in survey.py can be easily 
converted into combinatorial triangulations in Regina.
