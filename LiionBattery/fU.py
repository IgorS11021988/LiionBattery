import numpy as np


# Функция условий протекания процессов
def fU(t,  # Моменты времени
       UParametersSystemParameters  # U-параметры системы
       ):
    # Выделяем параметры динамики
    [I,  # Ток во внешней цепи
     IStepDischarge,  # Ток разряда после ступеньки
     fI,  # Частота колебаний внешнего тока
     AI,  # Амплитуда колебаний тока во внешней цепи
     fosI,  # Частота качания частоты колебаний внешнего тока
     AosI,  # Амплитуда качания частоты колебаний тока во внешней цепи
     tIStep,  # Время начала конца заряда
     tauEndCharge,  # Постоянная времени конца заряда
     alphaEndCharge,  # Экспоненциальный коэффициент постоянной времени разряда
     powEndCharge  # Степень постоянной времени разряда
     ] = UParametersSystemParameters  # Постоянная составляющая тока

    # Рассчитываем ток во внешней цепи
    Icur = np.full_like(t, I, dtype=np.double)  # Массив токов во внешней цепи
    bIcur = (t > tIStep)  # Моменты времени после ступеньки
    if np.any(bIcur):  # Есть моменты времени после ступеньки
        if I < 0:
            # Относительные моменты времени конца разряда
            rt = (t[bIcur] - tIStep) / tauEndCharge

            # Экспоненциальный коэффициент
            rLambda = 1 + alphaEndCharge * np.power(rt, powEndCharge)

            # Ток во внешней цепи
            Icur[bIcur] *= np.exp(-rLambda * rt)
        else:
            # Токи после ступеньки
            Icur[bIcur] = IStepDischarge
    fI += AosI * np.sin(2 * np.pi * fosI * t)  # Колебания частоты
    Icur += AI * np.sin(2 * np.pi * fI * t)  # Учитываем колебания

    # Выводим результат
    return Icur
