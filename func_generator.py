import functions as F 
import PARAMS as P
import numpy as np
import math as m

#Создание массива случайных значений размером shape
def get_random_in_interval(interval = [-1, 1], shape = 1):
    return np.random.rand(shape) * (interval[1] - interval[0]) + interval[0]

#Создаё случайную функцию (на основе одной базовой)
def generate_single_function(function):
    dx = get_random_in_interval(P.dx_interval) # создаём массив dx
    ax = get_random_in_interval(P.ax_interval) # создаём массив ax 
    a = get_random_in_interval(P.a_interval)   # создаём массив a
    
    func = lambda x: function(x, dx, ax, a)
    
    tune_coeff = max(abs(func(P.x_interval[0])),abs(func(P.x_interval[1])))
    
    #случайно выбираем 1 или -1
    is_neg = (np.random.randint(0,2,1) - 0.5) * 2
    
    func_str = ' {:+3.10f} × {}'.format((a/tune_coeff*is_neg)[0],P.func_name[function].format('{:+3.6f}x {:+3.6f}'.format(ax[0],dx[0])))
    
    return (lambda x: a/tune_coeff*is_neg*func(x), func_str)

#Создаём функцию из complexity базовых
def generate_complex_function(func_list, complexity = 3, chances_count = 1, dy_range = [-2, 2]):
    if len(func_list)*chances_count < complexity:
        chances_count = m.ceil(complexity/len(func_list))
    
    chances = np.full(len(func_list),chances_count)
    res_func_list = np.empty(complexity,dtype = type(func_list[0]))
    
    func_str = ''
    
    for i in range(complexity):
        probs = chances / np.sum(chances) 
        random_func = np.random.rand(1) 
        fits_list = np.logical_and(np.cumsum(probs) >= random_func, chances > 0) 
        chosen_func = func_list[fits_list][0]
        new_func, new_func_str = generate_single_function(chosen_func)
        res_func_list[i] = new_func
        func_str += new_func_str
        index = np.where(func_list == chosen_func)[0][0]
        chances[index] = chances[index] - 1
    
    dy = get_random_in_interval(dy_range)
    func_str += ' {:+3.6f}'.format(dy[0])  
    res_func = lambda x: np.sum([i(x) for i in res_func_list]) + dy[0]
    
    #Возврат: ф-я от x и запись этой функции
    return (res_func, func_str)

#Создаём датасет для функции из "число" точек на интервале. 
def generate_dataset(func, dots_count = 100, x_range = P.x_interval, noise_ratio = 0.25):
    #случайно определяем точки, в которых будем считать функцию
    x = np.random.rand(dots_count) * (x_range[1] - x_range[0]) + x_range[0]
    
    #на этих точках - считаем функцию с маленькой ошибкой
    y = np.array([func(i) + np.random.randn(1) * noise_ratio for i in x])
    return (x,y)
