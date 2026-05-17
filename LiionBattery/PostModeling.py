import os


def PostModeling(allPars,  # Параметры моделирования с индексами
                 saveDynamicFun,  # Функтор сохранения динамики
                 PathResult,  # Путь к результатам

                 # Файл CSV
                 sep,  # Сепаратор CSV
                 dec  # Десятичный разделитель
                 ):
    # Формируем имя файла параметров
    ParametersFileName = os.path.join(PathResult, "Parameters.csv")

    # Сохраняем параметры
    allPars.to_csv(ParametersFileName,
                   sep=sep, decimal=dec,
                   index=False)

    # Возвращаем отсутствие ошибки
    return 0
