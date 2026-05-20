import os


def PostModeling(allPars,  # Параметры моделирования с индексами
                 saveDynamicFun,  # Функтор сохранения динамики
                 PathResult  # Путь к результатам
                 ):
    # Формируем имя файла параметров
    ParametersFileName = os.path.join(PathResult[0], "Parameters.csv")

    # Сохраняем параметры
    allPars.to_csv(ParametersFileName,
                   sep=PathResult[1],
                   decimal=PathResult[2],
                   index=False)

    # Возвращаем отсутствие ошибки
    return 0
