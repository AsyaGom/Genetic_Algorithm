import random as _random
from Base import genotipHemmingDistance

'''
Всем даным алгоритмам на вход подаётся массив пар [genotip, приспособл],
желательно отсортированных по уменьшению пригодности

т.е. если целевую функцию МИНИМИЗИРУЕМ - сортировка от меньшего к большему
если МАКСИМИЗИРУЕМ - от больщего к меньшему

на возврате - набор массивов с такими же парами [genotip, приспособл]
'''

'''kol - сколько наборов родителей требуется'''


#Панмиксия
def panmixia(cells):
    #Случайно пронумеровали все элементы популяции
    lst = [_random.randint(0,len(cells)-1) for i in cells]
    
    ret = [
        [cell for cell,n in zip(cells,lst) if n==i] 
        for i in range(len(cells)-1)
    ]
    
    return [r for r in ret if len(r)>1]


#############################
#Генотипный ин- и аутбридинг
def _breeding(cells, t, kol = 1): 
    ret = [None for i in range(kol)]
    
    for i in range(kol):
        #Выбираем случайного первого радителя
        base = cells[_random.randint(0,len(cells)-1)]

        #Создаём список с расстояниями от всех клеток до выбранной
        lst = [
            [cell, genotipHemmingDistance(base[0], cell[0])]
            for cell in cells if not cell is base
        ]

        #сортируем по "удалению"
        lst.sort(key = lambda i: i[1])

        #Сохраняем пару родителей
        ret[i] = [
            min( base, lst[t][0],  key = lambda i: i[1]),
            max( base, lst[t][0],  key = lambda i: i[1])
        ]
        
    #Возвращаем родителей
    return ret

outbreeding = lambda cells, kol=1: _breeding(cells, -1, kol=kol) #t=-1 - аутобридинг



#############################################
#Словарь доступных функций
ParentsSelection = {
    "funcs": {
        "panmixia": panmixia,
        "outbreeding": outbreeding
    }
}