from MathProtEnergyProcSynDatas import RandomGenerateAttributesAndDynamicParametersBase


def OptimizeModeling(attributesBorder,  # Границы аттрибутов
                     dynamicParametersBorder,  # Границы начального состояния

                     nOptimizeModes,  # Число аттрибутов режима
                     attributesNPoints,  # Число точек аттрибутов
                     dynamicParametersNDyblicates,  # Число состояний, определяющих конкретную динамику

                     PathResultOptimize,  # Путь к результату

                     modelDynamicsFun  # Функция моделирования динамик
                     ):
    # Генерируем случайные знаяения аттрибутов и начальных состояний
    (attributes,
     dynamicParameters) = RandomGenerateAttributesAndDynamicParametersBase(attributesNPoints,  # Число точек аттрибутов
                                                                           nOptimizeModes,  # Число режимов работы
                                                                           dynamicParametersNDyblicates,  # Число состояний, определяющих конкретную динамику

                                                                           attributesBorder,  # Границы аттрибутов
                                                                           dynamicParametersBorder  # Границы динамических параметров
                                                                           )

    # Моделируем динамики
    modelDynamicsFun(dynamicParameters,  # Параметры динамики
                     attributes  # Аттрибуты
                     )

    # Сохраняем значения аттрибутов аккумулятора в файл
    attributes.to_csv(PathResultOptimize[0],
                      sep=PathResultOptimize[2],
                      decimal=PathResultOptimize[3],
                      index=False)

    # Сохраняем значения начальных состояний аккумулятора в файл
    dynamicParameters.to_csv(PathResultOptimize[1],
                             sep=PathResultOptimize[2],
                             decimal=PathResultOptimize[3],
                             index=False)

    # Возвращаем отсутствие ошибки
    return (dynamicParameters,
            attributes)
