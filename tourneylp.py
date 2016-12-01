from itertools import chain, combinations, permutations
import numpy as np
from pulp import *
from scipy.linalg import eig

N = 5
M = (N*(N-1))/2
K = 2
R = M - (K*(K-1))/2
EDGES = []

for x in xrange(N):
    for y in xrange(x+1, N):
        EDGES.append((x,y))

# miscellaneous functions

def all_subsets(x):
    return chain.from_iterable(combinations(x, r) for r in range(len(x)+1))

def hb(x, i):
    return (x&(1<<i))>0

def pad(s, l):
    return '0'*(l-len(s)) + s

def fact(n):
    ans = 1.
    for i in xrange(1,n+1):
        ans*=i
    return ans

# operations on tournament graphs

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

def permute(perm, x):
    nx = 0
    for i in xrange(N):
        for j in xrange(i+1, N):
            w = i if perm[i] < perm[j] else j
            if match_winner(x, i, j)==w:
                nx += 1<<edge_ind(perm[i], perm[j])

    return nx


# lp solving

def solveLP():
    LP = LpProblem("lp", LpMinimize)

    # Variables
    padlen = len(str(1<<M))
    probs = [[LpVariable("v{}_{}".format(pad(str(x), padlen), str(y)), 0, 1) \
        for y in xrange(N)] for x in xrange(1<<M)]

    uppers = [[LpVariable("upper{}_{}".format(str(w), str(ind)), 0, 1) \
        for ind, y in enumerate(combinations(range(N), K))] for w in xrange(1<<R)]
    lowers = [[LpVariable("lower{}_{}".format(str(w), str(ind)), 0, 1) \
        for ind, y in enumerate(combinations(range(N), K))] for w in xrange(1<<R)]
    alpha = LpVariable("alpha", 0, 1)

    # Objective
    LP += alpha

    # Constraints
    for x in xrange(1<<M):
        LP += (lpSum(probs[x]) == 1)
        if winner(x)>=0:
            LP += (probs[x][winner(x)] == 1)

    def gen_x(y, w):
        x = 0
        curpow = 1
        for i in xrange(N):
            for j in xrange(i+1, N):
                if (i not in y) or (j not in y):
                    x+=(w%2)*curpow
                    w/=2
                curpow*=2
        return x

    for ind, y in enumerate(combinations(range(N), K)):
        flip_edges = [1<<edge_ind(z[0], z[1]) for z in combinations(y, 2)]
        flip_bms = [sum(flips) for flips in all_subsets(flip_edges)]
        for w in xrange(1<<R):
            x = gen_x(y, w)

            #cur = lpSum([probs[x][z] for z in y])
            for bm in flip_bms:
                nx = x^bm
                ncur = lpSum([probs[nx][z] for z in y])
                LP += (uppers[w][ind] >= ncur)
                LP += (lowers[w][ind] <= ncur)

            LP += (alpha >= uppers[w][ind]-lowers[w][ind])

    GLPK(msg=0).solve(LP)

    #for v in LP.variables():
        #print v.name, "=", v.varValue

    #print "objective=", value(LP.objective)

    return [[probs[x][y].varValue for y in xrange(N)] for x in xrange(1<<M)]

# rules

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

def half_snm(x):
    prob = LpProblem("rivest", LpMinimize)

    lpvars = [LpVariable("p_{}".format(i), 0, 1) for i in xrange(N)]
    alpha = LpVariable("alpha", 0, 1)

    prob += alpha

    def losers(x, i):
        return [j for j in xrange(N) if match_winner(x, i, j) == i]

    prob += (lpSum([lpvars[i] for i in xrange(N)]) == 1)
    for i in xrange(N):
        prob += (lpSum([lpvars[j] for j in losers(x, i)]) <= alpha)

    GLPK(msg=0).solve(prob)

    #print bin(x)
    #for v in prob.variables():
        #print v.name, "=", v.varValue
    if alpha.varValue > 0.4:
        print alpha.varValue

    print x, [v.varValue for v in lpvars], alpha.varValue
    return alpha.varValue

    #TODO: return below instead of alpha.varValue (to make this a rule again)
    #return [v.varValue for v in lpvars]

def half_snm_matt(x):
    prob = LpProblem("rivest", LpMinimize)

    lpvars = [LpVariable("p_{}".format(i), 0, 1) for i in xrange(N)]
    alpha = LpVariable("alpha", 0, 1)

    prob += alpha

    # extend losers function to sets
    def losers(x, i):
        return [j for j in xrange(N) if match_winner(x, i, j) == i]

    def losers_set(x, S):
        ret = set()
        for i in S:
            ret.update(losers(x,i))
        return [x for x in ret if x not in S]

    prob += (lpSum([lpvars[i] for i in xrange(N)]) == 1)
    for S in all_subsets(range(N)):
        prob += (lpSum([lpvars[j] for j in losers_set(x, S)]) <= alpha)

    GLPK(msg=0).solve(prob)

    #print bin(x)
    #for v in prob.variables():
        #print v.name, "=", v.varValue
    if alpha.varValue > 0.5:
        print x, alpha.varValue
        print bin(x)
        for v in prob.variables():
            print v.name, "=", v.varValue

    # print x, [v.varValue for v in lpvars], alpha.varValue
    return alpha.varValue

    #TODO: return below instead of alpha.varValue (to make this a rule again)
    #return [v.varValue for v in lpvars]

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

def ariel(x):
    freq = [0]*N

    w = winner(x)
    if w!=-1:
        freq[w] = 1.
        return freq

    for comb in combinations(range(N), K):
        if any([all([match_winner(x, i, j)==1 for j in xrange(N) if j not in comb]) for i in comb]):
            for i in comb:
                freq[i]+=1

    m = min(freq)
    for i in xrange(N):
        freq[i]-=m

    s = sum(freq)
    if s>0:
        for i in xrange(N):
            freq[i]/=(1.0*s)
    else:
        for i in xrange(N):
            freq[i] = 1.0/N

    return freq


# operations on rules

def mult(rule1, rule2):
    def mult_rule(x):
        p = rule1(x)
        q = rule2(x)

        ans = [0. for _ in xrange(N)]

        for i in xrange(N):
            for j in xrange(N):
                if j==i:
                    ans[i] += p[i]*q[i]
                else:
                    if match_winner(x, i, j)==i:
                        ans[i] += p[i]*q[j] + q[i]*p[j]

        return ans

    return mult_rule

def square(rule):
    return mult(rule, rule)

def symmetrize(table):
    ntable = [[0. for y in xrange(N)] for x in xrange(1<<M)]
    nfac = fact(N)
    print 'factorial is', nfac
    for x in xrange(1<<M):
        for perm in permutations(range(N)):
            px = permute(perm, x)
            for y in xrange(N):
                ntable[px][perm[y]] += table[x][y]/nfac

    return ntable

def table_rule(table):
    def rule(x):
        return table[x]
    return rule

def check_rule(rule):
    ps = [rule(x) for x in xrange(1<<M)]
    print 'Done computing ps'

    alpha = 0
    for x in xrange(1<<M):
        print x, ps[x]
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
    #table = solveLP()
    #table = symmetrize(table)
    #check_rule(table_rule(table))
    #check_rule(ariel)
    maxalpha = 0.
    for x in xrange(1<<M):
        cur = half_snm_matt(x)
        maxalpha = max(maxalpha, cur)

    print 'Largest alpha:', maxalpha
