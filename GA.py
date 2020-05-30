from matplotlib import pyplot as plt


#Модуль для прогрессбара
from tqdm import tqdm

from Base import Base
from ParentsSelection import ParentsSelection
from Recombination import Recombination
from Mutation import Mutation, fullMutation
from Selection import Selection



'''
Для всех вариантов ГА критерий останова - прошло _step_ итераций
'''

#################################################
#Вариант генетического алгоритма из первой версии
#<точный, но долгий>
def ga_0(baseFuncs, data, steps=1000):
    #Более краткая запись генератора функций
    getCells = lambda genotips: Base["funcs"]["getCells"](genotips, data, baseFuncs)
    
    genotips = [Base["funcs"]["randGenotip"](baseFuncs)]
    cells = getCells(genotips)
    
    for i in tqdm(range(steps)):
        #Берём 2 лучшие клетки
        cells = Selection["funcs"]["eliteSelection"](cells, 2)
        
        genotips = [cell[0] for cell in cells] #извлекаем из них генотипы
        genotips = fullMutation(genotips)       # новая популяция
        
        cells = getCells(genotips) #оченка пригодности
        
    return cells[0]


#################
#канонический ГА
#<Быстрый но плохо работает>
def ga_1(baseFuncs, data, steps=1000, N = 64):
    #Короткая запись функции перерасчёта пригодности
    estimation = lambda genotip: Base["funcs"]["estimation"](genotip, data, baseFuncs)
    
    genotips = [Base["funcs"]["randGenotip"](baseFuncs) for i in range(N)]
    cells = Base["funcs"]["getCells"](genotips, data, baseFuncs)
    
    for i in tqdm(range(steps)):
        #пропорциональный отбор заменён отбором усечением
        lst = Selection["funcs"]["truncationSelection"](cells, N)     # промежуточный массив
        parents = [ParentsSelection["funcs"]["panmixia"](lst)[0][:2]] # из которого выбираются 2 родителя
        
        #одноточечным кроссинговером скрещиваем
        children = Recombination["funcs"]["singlePointCrossover"](parents)[0]
        
        #Мутация
        for child in children:
            Mutation["funcs"]["oncePointMutation"](child)
            
        #Заменяем родительские хромосомы на хромосомы потомков; перерасч. пригодность
        parents[0][0][0][:] = children[0]; parents[0][0][1] = estimation(children[0])
        parents[0][1][0][:] = children[1]; parents[0][1][1] = estimation(children[1])
        
        #Сортируем
        cells.sort(key = lambda i: i[1])
                
    return cells[0]



########
#Генитор
#<Чуть медленнее канонического, чуть хуже первого>
def ga_2(baseFuncs, data, steps=1000, N = 64):
    #Более короткая запись функции расчёта прригодности
    estimation = lambda genotip: Base["funcs"]["estimation"](genotip, data, baseFuncs)
    
    #Создаём популяцию и расчитываем пригодность особей
    genotips = [Base["funcs"]["randGenotip"](baseFuncs) for i in range(N)]
    cells = Base["funcs"]["getCells"](genotips, data, baseFuncs)
    
    flag = False
    
    for i in tqdm(range(steps)):
        parents = ParentsSelection["funcs"]["panmixia"](cells)
        
        
        #Две особи, из которых получим одну - больше похоже на описание рулеточной 
        # рекомбинации
        children = Recombination["funcs"]["rouletteWheelRecombination"](parents)
        
        #Мутация - случайно выбрана одноточечная
        for child in children:
            Mutation["funcs"]["oncePointMutation"](child)
        
        
        flag ^= True #меняем состояние флага
        if flag:#если флаг==Истина - заменяем генотип худшего родителя
            for _parents, child in zip(parents, children):
                l = max(_parents, key = lambda i: i[1])
                l[0][:] = child
                l[1] = estimation(child) #перерасчёт пригодности
        else:#иначе - худшего в популяции
            cells[-1][0][:] = children[0]
            cells[-1][1] = estimation(children[0]) #перерасчёт пригодности
        
        #Сортируем
        cells.sort(key = lambda i: i[1])
        
    return cells[0]



####
#CHC
def ga_3(baseFuncs, data, steps=1000, N = 64):
    #Короткая запись функции перерасчёта пригодности
    estimation = lambda genotip: Base["funcs"]["estimation"](genotip, data, baseFuncs)
    
    genotips = [Base["funcs"]["randGenotip"](baseFuncs) for i in range(N)]
    cells = Base["funcs"]["getCells"](genotips, data, baseFuncs)
    
    for i in tqdm(range(1,steps+1)):
        cells = Selection["funcs"]["eliteSelection"](cells, N)
        
        #Раз в 500 итераций происходит cataclysmic mutation
        if i%500==0:
            for j in range(1, N):
                Mutation["funcs"]["realValuedMutation"](cells[j][0], m=3)
                cells[j][1] = estimation(cells[j][0])
        
        
        parents = ParentsSelection['funcs']["outbreeding"](cells, kol=5)
        children = Recombination["funcs"]["HUX"](parents)
        
        #Добавляем детей в виде клеток
        for _children in children:
            cells += [
                [genotip, estimation(genotip)]
                for genotip in _children
            ]
        
        #Сортируем
        cells.sort(key = lambda i: i[1])
                
    return cells[0]


##############################
#С нефиксир размером популяции
#<Оооооооооооооочень долгий, точность не определена>
def ga_4(baseFuncs, data, steps=1000, N = 64, maxAge = 3):
    estimation = lambda genotip: Base["funcs"]["estimation"](genotip, data, baseFuncs)
    
    
    genotips = [Base["funcs"]["randGenotip"](baseFuncs) for i in range(N)]
    cells = [
        [genotip, estimation(genotip), maxAge]
        for genotip in genotips
    ]
    
    minE = maxE = cells[0][1]
    
    
    for i in tqdm(range(steps)):
        for cell in cells:
            if maxE < cell[1]: maxE = cell[1]
            if minE > cell[1]: minE = cell[1]
                
        
        #print(i, '\t', len(cells), minE)
        
        #Выбираем родителей, скрещиваем
        parents = ParentsSelection["funcs"]["panmixia"](cells)        # выбрали родителей
        children = Recombination["funcs"]["rouletteWheelRecombination"](parents) #перекомбинировали
        
        
        #Отсчитываем каждой клетке год жизни
        # и если пришёл срок - убиваем
        cells = [
            [cell[0], cell[1], cell[2]-1]
            for cell in cells if cell[2]>0
        ]
        
        cells_c = [None for child in children]
        for i, child in enumerate(children):
            Mutation["funcs"]["oncePointMutation"](child)
            e = estimation(child)
            cells_c[i] = [
                child, 
                e, 
                int(maxAge - (e-minE)/(maxE-minE)*maxAge)
            ]
        cells += cells_c
        
    cells.sort(key = lambda i: i[1])
    return cells[0][0:2]






# Инструкция по выбору:
# какой вариант алгоритма будет использоваться, может принимать значения:
# - 0 - первая версия (выбираем 2 лучших представителя, делаем кучу копий, в каждой копии заменяем по одному нуклеотиду)
# - 1 - Канонический ГА
# - 2 - Генитор
# - 3 - СНС
# - 4 - С нефиксированным размером популяции (оооооооооочень медленный)

GA = [
    ga_0, # Первый
    ga_1, # Канонический
    ga_2, # Генитор
    ga_3, # CHC
    ga_4, # С нефиксир размером популяции
]

