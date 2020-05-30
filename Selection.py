import random as _random
from math import exp as _exp

'''
На вход - отсортированный набор пар генотип-пригодность
(первая пара с лучшей пригодностью, последняя - с худшей)

На выход - набор из N таких же пар
'''


#Отбор усечением
def truncationSelection(cells, N, t=0.5):
    return [
        cells[_random.randint(0, int(len(cells)*t-1) )]
        for i in range(N)
    ]


#Элитарный отбор
def eliteSelection(cells, N):
    return cells[:N]





Selection = {
    "funcs": {
        "truncationSelection": truncationSelection,
        "eliteSelection": eliteSelection
    }
}