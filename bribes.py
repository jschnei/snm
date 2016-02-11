import random

VERBOSE = False

def verbose(s):
    if VERBOSE:
        print s

class Node:
    def __init__(self,left=None,right=None,player=None):
        self.left = left
        self.right = right
        self.parent = None
        
        self.player = player
        self.prize = 0.
        
    def __str__(self):
        if not self.left and not self.right:
            return str(self.player)
        else:
            return "{}[{}] ({} | {})".format(self.player, self.prize, self.left, self.right)
        
    def total_prize(self):
        ans = self.prize
        if self.left:
            ans += self.left.total_prize()
        if self.right:
            ans += self.right.total_prize()
            
        return ans

def win(p1, p2, matchups):
    return p1 if (matchups[p1][p2]==1) else p2

def build_SE_bracket(matchups, st, end):
    if end-st==1:
        return Node(player=st)
    
    mid = (st+end)/2
    lhalf = build_SE_bracket(matchups, st, mid)
    rhalf = build_SE_bracket(matchups, mid, end)
    winner = win(lhalf.player, rhalf.player, matchups)
    
    ret = Node(left=lhalf,right=rhalf,player=winner)
    lhalf.parent = ret
    rhalf.parent = ret
    
    return ret


def potential_earnings(match, player, matchups):
    ans = 0.
    while match and matchups[player][match.player]>=0:
        ans += match.prize
        match = match.parent
    
    return ans


def check_bribes(match, matchups):
    if not match.left or not match.right:
        return False
    
    p1 = match.left.player
    p2 = match.right.player
    
    winner = match.player
    loser = (p1+p2)-winner
    
    winner_earnings = potential_earnings(match.parent, winner, matchups)
    loser_earnings = potential_earnings(match.parent, loser, matchups)
    
    if loser_earnings > winner_earnings:
        verbose('Incentive to bribe in {} vs {}'.format(p1, p2))
        return True
    
    lans = check_bribes(match.left, matchups)
    rans = check_bribes(match.right, matchups)
    return (lans or rans)
    
def assign_prizes(bracket, matchups):
    if not bracket.left or not bracket.right:
        return
    
    if not bracket.parent:
        bracket.prize = 1.
    else:
        p1 = bracket.left.player
        p2 = bracket.right.player
        
        winner = bracket.player
        loser = (p1+p2)-winner
        
        winner_earnings = potential_earnings(bracket.parent, winner, matchups)
        loser_earnings = potential_earnings(bracket.parent, loser, matchups)
        
        bracket.prize = max(loser_earnings - winner_earnings, 0)
    
    assign_prizes(bracket.left, matchups)
    assign_prizes(bracket.right, matchups)

def random_matchups(N):
    matchups = [[0 for _ in xrange(N)] for _ in xrange(N)]
    for i in xrange(N):
        for j in xrange(i+1, N):
            matchups[i][j] = 2*random.randint(0,1) - 1
            matchups[j][i] = -matchups[i][j]
    return matchups

def trial(N):
    matchups = random_matchups(N)

    bracket = build_SE_bracket(matchups, 0, len(matchups))
    assign_prizes(bracket, matchups)
    
    return bracket.total_prize()
    
    
if __name__ == '__main__':
    samples = 10
    avg = sum([trial(1024) for _ in xrange(samples)])/samples
    print avg
    