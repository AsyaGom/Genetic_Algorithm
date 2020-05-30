import random as _random
from math import pi as _pi

'''
На вход - генотип
меняем в нём же
На выход - ничего
'''

#Набор излучений, изменяющих определённый ген в хромосоме на delta (изменения не выходят за рамки borders)
def _mutator_0(hrom, delta = 1e-2, borders=[0, _pi]):
    gen = hrom[0]
    
    l = max(gen-delta, borders[0])
    r = min(borders[1], gen+delta)
    
    hrom[0] = _random.uniform(l,r)
    
    
def _mutator_1(hrom, delta = 1e-2, borders=[0.1, 1.5]):
    gen = hrom[1]
    
    l = max(gen-delta, borders[0])
    r = min(borders[1], gen+delta)
    
    hrom[1] = _random.uniform(l,r)
    

def _mutator_2(hrom, delta = 1e-2, borders=[0., 0.]):
    gen = hrom[2]
    
    l = gen-delta
    r = gen+delta
    
    hrom[2] = _random.uniform(l,r)
    
    

def _mutator_3(hrom, delta = 1e-2, borders=[0., 0.]):
    gen = hrom[3]
    
    l = gen-delta
    r = gen+delta
    
    hrom[3] = _random.uniform(l,r)
    
    
    
#Массив для более быстрого обращения 
# к мутаторам
_delta = 1e-1
_genMutator = [
    lambda hrom: _mutator_0(hrom, delta = _delta, borders=[-10  , 10]),
    lambda hrom: _mutator_1(hrom, delta = _delta, borders=[0.001, 10]),
    lambda hrom: _mutator_2(hrom, delta = _delta, borders=[0.   , 0.]),
    lambda hrom: _mutator_3(hrom, delta = _delta, borders=[0.   , 0.]),
]

    
    
    
##########################################
##########################################
#Одноточечная мутация
def oncePointMutation(genotip):
    i = _random.randint(0, len(genotip   )-1) #номер хромосомы, в которой изменим
    j = _random.randint(0, len(genotip[0])-1) #номер гена, который изменим
    
    _genMutator[j](genotip[i])
    
    
#Мутация для вещественных особей
# Не очень понятно описание,
#  поэтому сведено к банальному - 
#   "если случайное число меньше 1/m - меняем"

def realValuedMutation(genotip, m=5):
    for hrom in genotip:
        for j in range(len(hrom)):
            if _random.random()<1/m:
                _genMutator[j](hrom)
                


#Функция для работы "первого" варианта ГА
def fullMutation(genotips):
    #Ф-я многократного копиорвания исходного генотипа
    def duplicator(genotip, N):
        return [ [hrom.copy() for hrom in genotip] for i in range(N)]
    
    #В зависимости от этого параметра будет изменяться
    #3 первых, 2 последних или все 4 гена в хромосоме
    var = [3,2,4]; t=2
    
    #Копируем каждый генотип необходимое число раз
    _newGenotips = []
    for genotip in genotips:
        _newGenotips += duplicator(genotip, 1+len(genotips[0])*var[t])
        
    #В каждой копии изменяем свой ген
    for j in range(len(genotips[0])):
        for i in [[0,1,2],[2,3],range(4)][t]:
            for genotip in _newGenotips[1+var[t]*j+i::1+len(genotips)*var[t]]: #код, благодоря которому изменяется i,j ген для копий всех поступивших генотипов
                _genMutator[i](genotip[j])
                
    return _newGenotips



######################
#########
Mutation = {
    "funcs": {
        "oncePointMutation": oncePointMutation,
        "realValuedMutation": realValuedMutation
    }
}