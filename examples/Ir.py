import sys
sys.path.insert(0, '..')
from cluster.crystal import RectangularCrystal


α = 4
geom = [
    ['Ir', [0,     0,   0]],
    ['Ir', [α/2, α/2,   0]],
    ['Ir', [α/2,   0, α/2]],
    ['Ir', [0,   α/2, α/2]],
]
Ir = RectangularCrystal(α, α, α, geom, 'cF')
ttt = Ir.tile(0, 2, 0, 2, 0, 2)
x, y, z = ttt.xyz.sum(axis=0)/len(ttt)
print(len(ttt))
print()
print(ttt.shifted(-x, -y, -z))
