import random as _random
import PARAMS as P
#import data_generator as P



###########################
#Создание генотипов

#Случайное заполнение генотипов
def randGenotip(baseFuncs, dx_range = [0, 1.5], ax_range = [0.1, 1.5], a_range = [0, 4], dy_range = [-3,3]):
    return [
        [
            _random.uniform(dx_range[0],dx_range[1]),
            _random.uniform(ax_range[0],ax_range[1]),
            _random.uniform( a_range[0], a_range[1]),
            _random.uniform(dy_range[0],dy_range[1])
        ]
        for i in baseFuncs
    ]


############################
#Целевая функия по-умолчанию
def _E(Y1, Y2):
    return  sum(((y1-y2))**2 for y1,y2 in zip(Y1,Y2))


###########################
#Функция запуска алгоритма
def go(ga, func_list, data, steps = 10000):
    genotip, prig = ga(func_list, data, steps)
    h = Base["funcs"]["assembly"](genotip, func_list)

    #Магия для генерации строки
    h_str = ''
    dy = 0
    for hrom, func in zip(genotip, func_list):
        if hrom[2]==0: continue

        dy += hrom[-1]
        h_str += ' {:+3.10f} × {} \n'.format(
            hrom[2], 
            P.func_name[func].format('{:-3.6f}x {:+3.6f}'.format(hrom[1],hrom[0]))
        )
    h_str += ' {:+3.6f}'.format(dy)
        
    return (h, prig, h_str)


######################
#Сбор функции из генотипа
def assembly(genotip, baseFuncs):
    return lambda x: sum([
        bF(x, dx=hrom[0], ax=hrom[1], a=hrom[2]) + hrom[3] 
        for hrom,bF in zip(genotip, baseFuncs)
    ])


#######################
#Вычисление пригодности для заданного генотипа
# при заданных данных
def estimation(genotip, data, baseFuncs, E = _E):
    h = assembly(genotip, baseFuncs)
    return E(data[1], [h(x) for x in data[0]])




#######################
#Сборка и сортировка пар [genotip, пригодность]
def getCells(genotips, data, baseFuncs, E = _E):
    cells = [
            [genotip, estimation(genotip, data, baseFuncs)]
            for genotip in genotips
        ]
    cells.sort(key = lambda i: i[1])
    return cells

    
    
########################################################
# Расстояние Хемминга между двумя хромосомами, двумя генотипами
def hromHemmingDistance(hrom1, hrom2, eps = 1e-4):
    #число совпавших генов
    return sum([1 for g1,g2 in zip(hrom1,hrom2) if abs(g1-g2)>eps])
    

def genotipHemmingDistance(genotip1, genotip2, eps = 1e-4):
    return sum([
        hromHemmingDistance(hrom1,hrom2,eps) 
        for hrom1, hrom2 in zip(genotip1,genotip2)
    ])



##############################
#Словарь с основными функциями
Base = {
    "funcs":{
        "E": _E,
        "Go": go,
        "assembly": assembly,
        "getCells": getCells,
        "randGenotip": randGenotip,
        "estimation": estimation,
        "hemmingDistance": genotipHemmingDistance
    }
}