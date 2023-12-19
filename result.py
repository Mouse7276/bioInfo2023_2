from Bio.Seq import Seq
import numpy as np
import random
import matplotlib.pyplot as plt

import variables
import criteria
import mutation
import crossover

# Create variable
seq = variables.createDNASeq(500)
Pt = variables.createPt(seq)
PtPrime = variables.createPtPrime(Pt)
P1 = variables.getP1(seq, Pt)

# Print result
print('sequence: ', seq, '\n')
print('Pt: ', Pt, '\n')
print('Pt\': ', PtPrime, '\n')
print('P1: ', P1, '\n')
print('length: ', variables.getDNALength(P1), '\n')
print('Tm: ', variables.getDNAMeltingTemp(P1), '\n')
print('GC content: ', variables.getDNAGCContent(P1), '\n')
print('annealing position appearance: ', variables.checkDNAAnnealingPositionAppearance(P1, seq), '\n')
print('3\' end primer: ', variables.checkDNA3EndPrimer(P1), '\n')

max_gen = 100
pop_size = 100
M = 50
Pe = 0.1
Pm = 0.1
max_itr = 50
ptlist, scorelist = [], []
for _ in range(max_gen):
    pt = variables.createPt(seq)
    ptlist.append(tuple(pt))
    scorelist.append(criteria.getFitnessValuePt(seq, pt))
gen = dict(zip(ptlist, scorelist))

itr = 1
while True:
    itr += 1

    if 0 in list(gen.keys())[0] or itr == max_itr:
        print('hello')
        break
    prob = [1/fitness for seq, fitness in gen.items()]
    cur = list(gen.keys())
    mating_pool_indices = np.random.choice(len(cur), M, prob)
    mating_pool = np.array(cur)[mating_pool_indices.astype(int)]
    
    # crossover
    for i in range(M):
        for j in range(M):
            rand_gen = np.random.rand()
            if rand_gen <= Pe:
                new_gen = crossover.crossoverPt(mating_pool[i], mating_pool[j])
                gen[tuple(variables.PtprimetoPt(new_gen[0]))] = criteria.getFitnessValuePt(seq, variables.PtprimetoPt(new_gen[0]))
                gen[tuple(variables.PtprimetoPt(new_gen[1]))] = criteria.getFitnessValuePt(seq, variables.PtprimetoPt(new_gen[1]))
    
    # mutation
    mutation_list = []
    for i in range(M):
        rand_gen = np.random.rand()
        if rand_gen <= Pm:
            mutation_list.append(variables.createPtPrime(mating_pool[i]))
    if mutation_list == []:
        pass
    elif len(mutation_list) == 1:
        tmp = mutation_list[0]
        mutation_list.append(tmp)
    else:
        new_gen = mutation.makeMutation(mutation_list)
        for e in new_gen:
            pt = variables.PtprimetoPt(e)
            score = criteria.getFitnessValuePt(seq, pt)
            gen[tuple(pt)] = score
    
    # evaluating
    gen = dict(sorted(gen.items(), key=lambda x:x[1]))
    seq_list = list(gen.keys())
    fit_list = list(gen.values())
    seq_list = seq_list[:pop_size]
    fit_list = fit_list[:pop_size]
    gen = dict(zip(seq_list, fit_list))
    
    print(f'{itr}th iteration. Best primer : {seq_list[0]}, fitness : {fit_list[0]}')
    
print(gen)
        
