# SageMath certificate accompanying
# "A Detailed Finite-Carry Proof of a Conjecture on
#  Permutation Binomials in Characteristic Two"
# Author: Yi Li

# Part I. Every one of the 24 gcds in Section 3.
R.<Y> = PolynomialRing(GF(2))
P = [1+Y+Y^8, 1+Y+Y^5, 1+Y^4+Y^5, 1+Y^7+Y^8]
H = [1+Y+Y^26, 1+Y+Y^23, 1+Y^22+Y^23,
     1+Y^2+Y^25, 1+Y^25+Y^26, 1+Y^23+Y^25]
Q = Y^2+Y+1
for i in range(4):
    for j in range(6):
        expected = P[1] if (i,j)==(1,0) else (
                   P[2] if (i,j)==(2,4) else Q)
        assert gcd(P[i],H[j]) == expected
assert P[1] == Q*(Y^3+Y^2+1)
assert P[2] == Q*(Y^3+Y+1)
print('PASS: all 24 polynomial gcds')

# Part II. All type +/- Fano conditions for the seven indices.
from fractions import Fraction
C = (1,3,5,7,9,11,13)
BASE = {'+': (0,1,3), '-': (0,2,3)}
LINES = {typ: {frozenset((h+a)%7 for a in base) for h in range(7)}
         for typ,base in BASE.items()}

def vanishes(c,b,typ,terms=3):
    alpha=3*c-1
    e=[(t-3*(alpha & ((b+43*t)%64)))%7 for t in range(terms)]
    if terms==2:
        return e[0]==e[1]
    return len(set(e))==3 and frozenset(e) in LINES[typ]

# The exact residue sets in the paper, stated modulo d.
TABLE = {
  1: ( 4,{1},              {0}),
  3: (16,{10,11,12},       {13,14,15}),
  5: (16,{10,12,14},       {11,13,15}),
  7: (32,{9,16,18,19},     {0,22,23,25}),
  9: (32,{1,7,9,11,13,16,19},{0,2,8,22,25,28,30}),
 11: (64,set(range(32,42)),set(range(0,10))),
 13: (64,{0,8,13,15,38,49,58,60},{3,26,28,33,41,45,47,56})
}
for c in C:
    d,plus,minus = TABLE[c]
    for b in range(64):
        assert vanishes(c,b,'+') == ((b%d) in plus)
        assert vanishes(c,b,'-') == ((b%d) in minus)
print('PASS: all seven rows of the trinomial residue table')

# Part III. Verify the four regions after c=1,3,5.
EXPECTED = {
 ('+',Fraction(1,3),Fraction(2,5)): [2,18,34,50],
 ('+',Fraction(2,3),Fraction(4,5)): [6,22,38,54],
 ('-',Fraction(1,5),Fraction(1,3)): [5,21,37,53],
 ('-',Fraction(3,5),Fraction(2,3)): [9,25,41,57]
}
result={}
for typ in ('+','-'):
    for k3 in range(3):
        for k5 in range(5):
            lo=max(Fraction(k3,3),Fraction(k5,5))
            hi=min(Fraction(k3+1,3),Fraction(k5+1,5))
            if lo>=hi: continue
            found=[]
            for rho in range(64):
                if all(vanishes(c,(c*rho-43*k-1)%64,typ)
                       for c,k in ((1,0),(3,k3),(5,k5))):
                    found.append(rho)
            if found:
                result[(typ,lo,hi)]=found
assert result==EXPECTED,(result,EXPECTED)
print('PASS: Table 2 (four nonboundary regions)')

# Part IV. Exact check of every carry cell and allowable boundary.
D=sorted({Fraction(j,c) for c in C for j in range(c+1)})
assert len(D)==42

def survivor(theta,rho,typ):
    for c in C:
        v=c*theta
        if v.denominator==1:             # two-term boundary
            b,terms=42,2
        else:
            k=v.numerator//v.denominator
            b,terms=(c*rho-43*k-1)%64,3
        if not vanishes(c,b,typ,terms):
            return False
    return True

for typ in ('+','-'):
    for lo,hi in zip(D,D[1:]):
        mid=(lo+hi)/2
        assert not any(survivor(mid,rho,typ) for rho in range(64))
    for theta in D[:-1]:
        j,d=theta.numerator,theta.denominator
        if d not in (1,3,9,11):continue
        rho=(43*j*pow(d,-1,64))%64
        assert not survivor(theta,rho,typ)
print('PASS: every open and arithmetically admissible boundary state')
