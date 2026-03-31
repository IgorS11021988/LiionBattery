import numpy as np

from MathProtEnergyProcSynDatas.TimesMoments import LinearTimesMoments
from MathProtEnergyProcSynDatas.Indicate import PlotGraphicIndicate, SaveDynamicToFileIndicate
from MathProtEnergyProcSynDatas.File import DynamicSaveAndSaveGraphics


# Функция расчета динамики
def InputArrayCreate(Pars,  # Параметры

                     integrateAttributes  # Аттрибуты интегрирования
                     ):  # Формирование массивов входных параметров
    # Корректируем токи во внешней цепи их аттрибуты
    Pars[["I",  # Ток во внешней цепи
          "IStep",  # Ток разряда после ступеньки, Cnom
          "AI"  # Амплитуда колебаний тока во внешней цепи
          ]] *= Pars[["Cnom"]].to_numpy()
    integrateAttributes[["bI0DCh",  # Граница разрядного нулевого тока
                         "bI0Ch"  # Граница зарядного нулевого тока
                         ]] *= Pars[["Cnom"]].to_numpy()

    # Добавляем границу по току в параметры
    Pars["bI0Ch"] = integrateAttributes["bI0Ch"].to_numpy()

    # Корректируем аттрибуты качающейся частоты
    Pars["AosI"] *= Pars["fI"] / 100  # Амплитуда качающейся частоты

    # Относительная емкость
    rCnom = Pars["Cnom"] / 3.003

    # Корректируем частотные характеристики
    Pars[["AosI", "fI", "fosI"]] *= 2 * np.pi

    # Емкости двойных слоев и мембраны
    Pars["Cbin0p"] *= rCnom  # Емкость положительного двойного слоя, Ф
    Pars["Cm"] *= rCnom  # Емкость мембраны, Ф
    Pars["Cbin0n"] *= rCnom  # Емкость отрицательного двойного слоя, Ф

    # Сопротивления двойных слоев и мембраны
    Pars["Rbin0p"] /= rCnom  # Сопротивление положительного двойного слоя, Ом
    Pars["Rm0"] /= rCnom  # Сопротивление мембраны, Ом
    Pars["Rbin0n"] /= rCnom  # Сопротивление отрицательного двойного слоя, Ом

    # Тепловые свойства литий-ионных аккумуляторов
    Pars["KInAkk"] *= rCnom  # Коэффициент теплопередачи содержимого литий-ионного элемента, Вт/К
    Pars["CInAkk"] *= rCnom  # Теплоемкость содержимого литий-ионного элемента, Дж/К
    Pars["KBAkk"] *= rCnom  # Коэффициент теплопередачи корпуса литий-ионного элемента, Вт/К
    Pars["CBAkk"] *= rCnom  # Теплоемкость корпуса литий-ионного элемента, Дж/К

    # Сопротивление клемм, Ом
    Pars["Rkl"] /= rCnom

    # Переводим номинальную емкость в кулоны
    Pars["Cnom"] *= 3600

    # Корректируем доли неактивных ячеек
    nonCeilNames = ["nuNonCeilp", "nuNonCeiln"]  # Имена долей разрешенных ячеек
    bIndNonCeilMin = (Pars[nonCeilNames] < 0)
    if np.any(bIndNonCeilMin):
        Pars.loc[bIndNonCeilMin["nuNonCeilp"].to_numpy(), "nuNonCeilp"] = 0
        Pars.loc[bIndNonCeilMin["nuNonCeiln"].to_numpy(), "nuNonCeiln"] = 0
    bIndNonCeilMax = (Pars[nonCeilNames] > 100)
    if np.any(bIndNonCeilMax):
        Pars.loc[bIndNonCeilMax["nuNonCeilp"].to_numpy(), "nuNonCeilp"] = 100
        Pars.loc[bIndNonCeilMax["nuNonCeiln"].to_numpy(), "nuNonCeiln"] = 100
    Pars[nonCeilNames] /= 100

    # Корректируем долю разрушенного положительного электрода
    bIndNuDegPosElMin = (Pars["nuDegPosEl"] < 0).to_numpy()
    if np.any(bIndNuDegPosElMin):
        Pars.loc[bIndNuDegPosElMin, "nuDegPosEl"] = 0
    bIndNuDegPosElMax = (Pars["nuDegPosEl"] > 100).to_numpy()
    if np.any(bIndNuDegPosElMax):
        Pars.loc[bIndNuDegPosElMax, "nuDegPosEl"] = 100
    Pars["nuDegPosEl"] /= 100

    # Полные зарядовые числа молей электродов
    Pars[["qMatAllp", "qMatAlln"]] = Pars[["cMatAllp", "cMatAlln"]].to_numpy() * Pars[["Cnom"]].to_numpy()

    # Число молей неразрушенного положительного электрода
    Pars["qNoDegPosEl"] = Pars["qMatAllp"].to_numpy() * (1 - Pars["nuDegPosEl"].to_numpy())

    # Начальные числа молей деградированных и недеградированных материалов электродов
    Pars["qDegPosEl"] = Pars["qMatAllp"].to_numpy() * Pars["nuDegPosEl"].to_numpy()
    Pars[["qMatDegElp", "qMatDegEln"]] = Pars[nonCeilNames].to_numpy() * Pars[["qNoDegPosEl", "qMatAlln"]].to_numpy()
    Pars[["qMatElp", "qMatEln"]] = Pars[["qNoDegPosEl", "qMatAlln"]].to_numpy() - Pars[["qMatDegElp", "qMatDegEln"]].to_numpy()

    # Удаляем зарядовое число молей неразрушенного положительного электрода
    Pars.drop(columns=["qNoDegPosEl"], axis=1, inplace=True)

    # Начальное состояние
    Pars["qbinp"] *= (Pars["EbinpC"] * (1 - Pars["q"]) + Pars["EbinpD"] * Pars["q"]) * Pars["Cbin0p"]  # Заряд на положительном двойном слое, Кл
    Pars["qm"] *= rCnom  # Заряд мембраны, Кл
    Pars["qbinn"] *= (Pars["EbinnC"] * (1 - Pars["q"]) + Pars["EbinnD"] * Pars["q"]) * Pars["Cbin0n"]  # Заряд на отрицательном двойном слое, Кл
    Pars["q"] *= Pars["Cnom"]  # Перенесенный через внешнюю цепь заряд, Кл
    Pars["TInAkk"] += Pars["Tokr"]  # Начальная температура сождержимого литий-ионного элемента, град С
    Pars["TBAkk"] += Pars["Tokr"]  # Начальная температура корпуса литий-ионного элемента, град С

    # Корректируем температуры
    Pars[["TInAkk", "TBAkk", "Tokr", "bRTp", "bRTm", "bRTn"]] += 273.15

    # Время интегрирования
    Tints = np.array(integrateAttributes["TintI0"], dtype=np.double)  # Времена интегрирования при нулевых токах или при заряде
    bIDisCharge = (Pars["I"].to_numpy() > integrateAttributes["bI0DCh"].to_numpy())
    if np.any(bIDisCharge):
        # Рассчитываем время интегрирования для ненулевого тока
        IZeros = 1.011 * Pars["I"][bIDisCharge]  # Ненулевой ток (по модулю)
        Tints[bIDisCharge] = (Pars["Cnom"][bIDisCharge] - Pars["q"][bIDisCharge]) / IZeros  # Время интегрирования

    # Время начала второй ступени
    Pars["cTIStep"] = np.array(Pars["cTIStep"], dtype=np.double).reshape(-1,)  # Приводим коэффициент времени ступенчатого перехода к массиву
    Pars["cTIStep"] *= Tints  # Время начала конца разряда

    # Корректируем время разряда
    if np.any(bIDisCharge):
        # Оставшееся время
        tost = Tints[bIDisCharge] - Pars.loc[bIDisCharge, "cTIStep"]
        bIDisCharge[bIDisCharge] = (tost > 0)  # Индексы ненулевого оставшегося времени

        if np.any(bIDisCharge):
            # Оставшийся заряд
            cqost = Pars["Cnom"][bIDisCharge] - Pars["q"][bIDisCharge] - Pars["I"][bIDisCharge] * Pars["cTIStep"][bIDisCharge]

            # Времена интегрирования
            IStepZeros = 1.011 * np.abs(Pars["IStep"][bIDisCharge])  # Ненулевой ток (по модулю)
            Tints[bIDisCharge] = Pars["cTIStep"][bIDisCharge] + cqost / IStepZeros  # Обновленные времена интегрирования

    #  Моменты времени
    NPoints = np.array(integrateAttributes["NPoints"], dtype=np.int32)  # Числа точек интегрирования
    ts = LinearTimesMoments(Tints,  # Времена интегрирования
                            NPoints  # Числа точек интегрирования
                            )

    # Возвращаем исходные данные динамики системы
    return (Pars, Tints, ts)


# Обработка результатов моделирования динамик
def OutputValues(dyns, fileName,
                 sep, dec, index,
                 plotGraphics=False  # Необходимость построения графиков
                 ):
    # Получаем величины из кортежа
    (t, Ukl, Ubinp, Ubinn, Um,
     TInAkk, TBAkk, q, Ibinp, Im, Ibinn, Icur, Tokr,
     qMatElp, qMatEln,
     qMatDegElp, qMatDegEln, qDegPosEl) = dyns

    # Заголовки и динамики
    dynamicsHeaders = {"Time": t,
                       "Ukl": Ukl,
                       "Ubinp": Ubinp,
                       "Ubinn": Ubinn,
                       "Um": Um,
                       "TInAkk": TInAkk,
                       "TBAkk": TBAkk,
                       "q": q,
                       "Icur": Icur,
                       "Ibinp": Ibinp,
                       "Ibinn": Ibinn,
                       "Im": Im,
                       "Tokr": Tokr,
                       "qMatElp": qMatElp,
                       "qMatEln": qMatEln,
                       "qMatDegElp": qMatDegElp,
                       "qMatDegEln": qMatDegEln,
                       "qDegPosEl": qDegPosEl
                       }

    # Одиночные графики на полотне
    oneTimeValueGraphics = [{"values": Ukl,  # Величины в моменты времени
                             "graphName": "Напряжение на клеммах",  # Имя полотна
                             "yAxesName": "Напряжение, В",  # Имя оси ординат
                             "graphFileBaseName": "AkkVoltage"  # Имя файла графика
                             },

                            {"values": qDegPosEl,  # Величины в моменты времени
                             "graphName": "Зарядовое число молей разрушенного положительного электрода",  # Имя полотна
                             "yAxesName": "Зарядовое число молей",  # Имя оси ординат
                             "graphFileBaseName": "AkkDegPosEl"  # Имя файла графика
                             }]

    # Группы графиков на полотне
    timesValuesGraphics = [{"listValues": [TInAkk, TBAkk],  # Список величин в моменты времени
                            "listValuesNames": ["Содержимое", "Корпус"],  # Список имен величин (в моменты времени)
                            "graphName": "Температура в литий-ионном аккумуляторе",  # Имя полотна
                            "yAxesName": "Температура, град С",  # Имя оси
                            "graphFileBaseName": "AkkTemperatures"  # Имя файла графика
                            },

                           {"listValues": [Ubinn, Ubinp, Um],  # Список величин в моменты времени
                            "listValuesNames": ["Отрицательный двойной слой",
                                                "Положительный двойной слой",
                                                "Мембрана"],  # Список имен величин (в моменты времени)
                            "graphName": "Напряжения в литий-ионном аккумуляторе",  # Имя полотна
                            "yAxesName": "Напряжение, В",  # Имя оси
                            "graphFileBaseName": "InAkkVoltage"  # Имя файла графика
                            },

                           {"listValues": [Ibinn, Ibinp, Im, Icur],  # Список величин в моменты времени
                            "listValuesNames": ["Отрицательный двойной слой",
                                                "Положительный двойной слой",
                                                "Мембрана",
                                                "Ток во внешней цепи"],  # Список имен величин (в моменты времени)
                            "graphName": "Токи в литийионном аккумуляторе",  # Имя полотна
                            "yAxesName": "Ток, Cnom",  # Имя оси
                            "graphFileBaseName": "AkkCurrents"  # Имя файла графика
                            },

                           {"listValues": [qMatElp, qMatDegElp],  # Список величин в моменты времени
                            "listValuesNames": ["Недеградированный материал",
                                                "Деградированный материал"],  # Список имен величин (в моменты времени)
                            "graphName": "Материалы положительного электрода",  # Имя полотна
                            "yAxesName": "Зарядовое число молей",  # Имя оси
                            "graphFileBaseName": "AkkDegPositiveElectrode"  # Имя файла графика
                            },

                           {"listValues": [qMatEln, qMatDegEln],  # Список величин в моменты времени
                            "listValuesNames": ["Недеградированный материал",
                                                "Деградированный материал"],  # Список имен величин (в моменты времени)
                            "graphName": "Материалы отрицательного электрода",  # Имя полотна
                            "yAxesName": "Зарядовое число молей",  # Имя оси
                            "graphFileBaseName": "AkkDegNegativeElectrode"  # Имя файла графика
                            }]

    # Сохраняем динамику в .csv файл и отображаем графики
    DynamicSaveAndSaveGraphics(dynamicsHeaders,  # Словарь динамик с заголовками
                               fileName,  # Имя файла динамик

                               t,  # Моменты времени
                               oneTimeValueGraphics,  # Один график на одном полотне
                               timesValuesGraphics,  # Несколько графиков на одном полотне

                               plotGraphics,  # Необходимость построения графиков

                               sep, dec,   # Разделители (csv и десятичный соответственно)

                               saveDynamicIndicator=SaveDynamicToFileIndicate,  # Индикатор сохранения динамики
                               saveGraphicIndicator=PlotGraphicIndicate,  # Индикатор отображения графиков
                               index=index  # Индекс динамики
                               )
