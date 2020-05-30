import random as _random
'''
На вход - набор массивов с парами генотип-пригодность
На выход набор массивов с генотипами
'''


#Рекомбинация для рулеточного отбора
# (настроена для МИНИМИЗИРУЕМОГО функционала)
def rouletteWheelRecombination(parents):
    ret = [None for i in parents]
    
    for k, _parents in enumerate(parents):
        #Массив с вероятностями отбора
        f = [_parents[-1][1] + 1e-5 - parent[1] for parent in _parents] #перевод от МИН к МАКС
        s = sum(f)
        intervals = [
            _f/s
            for _f in f
        ]
        for i in range(1, len(intervals)):
            intervals[i] += intervals[i-1]

        #Функция для поиска, какого родителя выбрать
        def _getParent(num):
            for i in range(len(intervals)):
                if num < intervals[i]: return i
            return -1 #если попали сюда - num==1, что соответствует последнему родителю

        ret[k] = [
            [
                _parents[ _getParent(_random.random()) ][0][i][j]
                for j in range(len(_parents[0][0][0])) #j - номер гена в хромосоме
            ]#создаём хромосому потомка
            for i in range(len(_parents[0][0])) #i - номер хромосомы
        ]#создаём генотип потомка
    
    return ret





'''
Кроссинговер реализуется на уровне хромосом, а не генов,
как остальные мутаторы
'''

#Одноточечный кроссинговер
def singlePointCrossover(parents):
    ret = [None for i in parents]
    
    for k, _parents in enumerate(parents):
        #Копируем родителей
        children = [
                [hrom.copy() for hrom in cell[0]]
                for cell in _parents 
        ]

        #Выбираем точку разреза
        sc = _random.randint(0, len(children[0])-1)

        #Перетасовываем 
        children[0][:sc], children[1][:sc] = children[1][:sc], children[0][:sc]
        
        ret[k] = children
        
    return ret



#Однородных кроссовер
def HUX(parents):
    ret = [None for i in parents]
    
    for k, _parents in enumerate(parents):
        #Копируем родителей
        children = [
                [hrom.copy() for hrom in cell[0]]
                for cell in _parents 
        ]

        #Выбираем точку разреза
        sc = int(len(children[0])/2)

        #Перетасовываем 
        children[0][:sc], children[1][:sc] = children[1][:sc], children[0][:sc]
        
        ret[k] = children
        
    return ret

Recombination = {
    "funcs": {
        "rouletteWheelRecombination": rouletteWheelRecombination,
        "singlePointCrossover": singlePointCrossover,
        "HUX": HUX
        
    }
}