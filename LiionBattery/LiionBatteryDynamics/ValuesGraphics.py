from matplotlib import pyplot as plt

#Функция построения графика одной величины
def OneTimeValueGraphic(times,#Моменты времени
                        values,#Величины в моменты времени
                        graphicsName,#Имя полотна
                        axeName,#Имя оси
                        needGrid = True#Нужна ли сетка
                        ):#Построение графика одной величины во времени
    #Строим полотно
    plt.figure()
    
    #Строим график
    plt.plot(times,values)
    
    #Добавляем надписи
    plt.title(graphicsName)#Имя полотна
    plt.xlabel("Время, с")#Имя оси времени
    plt.ylabel(axeName)#Имя оси величины
    
    #Добавляем сетку
    plt.grid(needGrid)
def TimesValuesGraphics(times,#Моменты времени
                        listValues,#Список величин в моменты времени
                        listValuesNames,#Список имен величин
                        graphicsName,#Имя полотна
                        axeName,#Имя оси
                        needGrid = True#Нужна ли сетка
                        ):#Построение графиков нескольких величин во времени
    #Строим полотно
    plt.figure()
    ax = plt.subplot(111)#Область графика
    
    #Строим графики
    nValues = len(listValuesNames)#Число величин
    inds = range(nValues)
    for ind in inds:
        plt.plot(times,listValues[ind],label=listValuesNames[ind])
    
    #Добавляем надписи
    plt.title(graphicsName)#Имя полотна
    plt.xlabel("Время, с")#Имя оси времени
    plt.ylabel(axeName)#Имя оси величины
    
    #Добавляем сетку
    plt.grid(needGrid)
    
    #Выделяем место для легенды
    if nValues == 2:
        posLegend = 0.843
    elif nValues == 3:
        posLegend = 0.783
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height*(1 - posLegend), box.width, box.height*posLegend])
    
    #Добавляем легенду
    plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.165))
    