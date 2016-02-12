from itertools import chain, combinations, permutations
import numpy as np
from pulp import *
from scipy.linalg import eig 

N = 6
M = (N*(N-1))/2
K = 3
EDGES = []

for x in xrange(N):
    for y in xrange(x+1, N):
        EDGES.append((x,y))


def all_subsets(x):
    return chain.from_iterable(combinations(x, r) for r in range(1,len(x)+1))

def hb(x, i):
    return (x&(1<<i))>0

def winner(x):
    scores = [0]*N
    for y in xrange(M):
        if hb(x,y):
            scores[EDGES[y][0]] += 1
        else:
            scores[EDGES[y][1]] += 1
            
    for y in xrange(N):
        if scores[y] == N-1:
            return y
    
    return -1

def edge_ind(i, j):
    if i>j:
        return edge_ind(j,i)
    else:
        return (i*(2*N-(i+1)))/2 + (j-i-1)

def match_winner(x, i, j):
    if i==j:
        return -1
    else:
        b = edge_ind(i,j)
        if hb(x, b):
            return EDGES[b][0]
        else:
            return EDGES[b][1]
        
def bracket_winner(x, bracket):
    if len(bracket)==1:
        return bracket[0]
    elif len(bracket)==2:
        return match_winner(x, bracket[0], bracket[1])
    else:
        lhalf = bracket[:len(bracket)/2]
        rhalf = bracket[len(bracket)/2:]
        lwinner = bracket_winner(x, lhalf)
        rwinner = bracket_winner(x, rhalf)
        return match_winner(x, lwinner, rwinner)


def solveLP():
    LP = LpProblem("lp", LpMinimize)

    # Variables
    probs = [[LpVariable("v{}_{}".format(str(x), str(y)), 0, 1) \
        for y in xrange(N)] for x in xrange(1<<M)]
    uppers = [[LpVariable("upper{}_{}".format(str(x), str(ind)), 0, 1) \
        for ind, y in enumerate(combinations(range(N), K))] for x in xrange(1<<M)]
    alpha = LpVariable("alpha", 0, 2)

    # Objective
    LP += alpha

    # Constraints
    for x in xrange(1<<M):
        LP += (lpSum(probs[x]) == 1)
        if winner(x)>=0:
            LP += (probs[x][winner(x)] == 1)
            
        for ind, y in enumerate(combinations(range(N), K)):
            flip_edges = [1<<edge_ind(z[0], z[1]) for z in combinations(y, 2)]
            flip_bms = [sum(flips) for flips in all_subsets(flip_edges)]
            
            cur = lpSum([probs[x][z] for z in y])
            for bm in flip_bms:
                nx = x^bm
                ncur = lpSum([probs[nx][z] for z in y])
                LP += (uppers[x][ind] >= ncur)
    
            LP += (alpha >= uppers[x][ind]-cur)
    
    GLPK().solve(LP)

    for v in LP.variables():
        print v.name, "=", v.varValue
        
    print "objective=", value(LP.objective)


def rivest(x):
    prob = LpProblem("rivest", LpMinimize)
    
    lpvars = [LpVariable("p_{}".format(i)) for i in xrange(N)]
    
    prob += lpSum(lpvars)
    
    def coeff(i, j):
        if i==j:
            return 1
        elif i>j:
            return 2-coeff(j,i)
        else:
            return 0 if hb(x, edge_ind(i,j)) else 2
    
    for i in xrange(N):
        prob += (lpvars[i] >= 0)
    
    for i in xrange(N):
        prob += (lpSum([coeff(i,j)*lpvars[j] for j in xrange(N)]) >= 1)
    
    GLPK(msg=0).solve(prob)
    
    #for v in prob.variables():
        #print v.name, "=", v.varValue
        
    #print "objective=", value(prob.objective)
    
    return [v.varValue for v in prob.variables()]

BRACKETS = list(permutations(range(N)))
def random_sing_elim(x):
    p = [0 for y in xrange(N)]
    
    for bracket in BRACKETS:
        b_winner = bracket_winner(x, bracket)
        p[b_winner] += 1.0/len(BRACKETS)
    
    return p

def page_rank(x):
    def coeff(i, j):
        frac = 1./(N)
        if i==j:
            deg = len([k for k in xrange(N) if match_winner(x, i, k)==i])
            return (deg+1)*frac
        elif i>j:
            return frac - coeff(j,i)
        else:
            return 0 if hb(x, edge_ind(i,j)) else frac
            
    transition_mat = np.matrix([[coeff(i,j) for j in xrange(N)] for i in xrange(N)])
    #print transition_mat

    S, U = eig(transition_mat.T)
    stationary = np.array(U[:, np.where(np.abs(S - 1.) < 1e-8)[0][0]].flat)
    stationary = stationary / np.sum(stationary)
    
    return [x.real for x in stationary]


def check_rule(rule):
    ps = [rule(x) for x in xrange(1<<M)]
    print 'Done computing ps'

    alpha = 0
    for x in xrange(1<<M):
        #print x, ps[x]
        for y in combinations(range(N), K):
            flip_edges = [1<<edge_ind(z[0], z[1]) for z in combinations(y, 2)]
            flip_bms = [sum(flips) for flips in all_subsets(flip_edges)]
            
            cur = sum([ps[x][z] for z in y])
            for bm in flip_bms:
                nx = x^bm
                new = sum([ps[nx][z] for z in y])
                alpha = max(alpha, new - cur)
                
                #if abs(new-cur - 0.5454545454) < 1e-5:
                    #print x, y, bm, nx, ps[x], ps[nx]
                
        
        #for y in xrange(M):
            #nx = x^(1<<y)
            #alpha = max(alpha, ps[x][EDGES[y][0]] + ps[x][EDGES[y][1]] - ps[nx][EDGES[y][0]] - ps[nx][EDGES[y][1]])
            
            #if not hb(x, y) and ps[nx][EDGES[y][0]] < ps[x][EDGES[y][0]]:
                #print 'not monotone!:', x, nx, y, ps[x], ps[nx]
            
            ## find equality cases
            #best = (N-2)/(1.*(2*N-3))  # best = 1./3
            #if abs(ps[x][EDGES[y][0]] + ps[x][EDGES[y][1]] - ps[nx][EDGES[y][0]] - ps[nx][EDGES[y][1]] - best) < 1e-5:
                #print x, EDGES[y][0], EDGES[y][1], ps[x], ps[nx]

    print 'alpha:', alpha
                                                                                                          
    
if __name__ == '__main__':
    solveLP()
    #check_rule(random_sing_elim)
    