#################################################################################
# 
# A scipt for converting 3-crystallisations in the format of survey.py into
# triangulations in Regina
# 
#################################################################################

import survey as surv
import regina

# These are input parameters
# All gems should have the same first permutation (third colour)
# Example: all (1^1,3^k)-type 3-crystallisation for n=7
n = 7
gems = surv.family_three_crystallisations(7,sphere1=[[1], [0, 4, 2], [3, 5, 6]])

idn = list(range(n))
mu = surv.mu_online(n)
sigmas = surv.cycle_to_oneline(gems[0][0],n)
taus = [surv.cycle_to_oneline(g[1]) for g in gems]
triangulations = []

for tau in taus:
	tau = []
	perms = [idn,mu,sigma,tau]
	T = regina.Triangulation3()
	simplices = [T.newTetrahedron() for i in range(2*n)]
	id4 = regina.Perm4()
	for c in range(4):
		for i in range(n):
			simplices[i].join(c,simplices[perms[c][i]+n],id4)
	triangulations.append(T)

# Example: computing each unique isomorphism signature for triangulations

unique_sigs = []
for T in triangulations:
	s = T.isoSig()
	if s not in unique_sigs:
		unique_sigs.append(s)

print "Distinct Isomorphism Signatures:\n"
for s in unique_sigs:
	print s