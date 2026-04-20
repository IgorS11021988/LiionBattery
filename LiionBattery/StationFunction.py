import numpy as np

from .StationFunctions import funEbin, funCbin, funRbin, funNuEMat, funMuMat, funADNu, funCf0, funKDDegPosEl, funCKbAkk
from MathProtEnergyProc import NonEqSystemQBase


# Функция состояния для литий-ионного аккумулятора
class IndepStateFunction(object):
    # Тело функции
    def __call__(self,

                 stateCoordinates,
                 reducedTemp,
                 systemParameters):
        # получаем электрические заряды
        [qbinp,  # Электрический заряд положительного двойного слоя
         qm,  # Электрический заряд мембраны
         qbinn,  # Электрический заряд отрицательного двойного слоя
         q,  # Перенесенный через внешнюю цепь электрический заряд
         qMatElp,  # Зарядовое сило молей недеградированного положительного электрода
         qMatEln,  # Зарядовое сило молей недеградированного отрицательного электрода
         qMatDegElp,  # Зарядовое сило молей деградированного положительного электрода
         qMatDegEln,  # Зарядовое сило молей деградированного отрицательного электрода
         qDegPosEl  # Зарядовое сило молей разрушенного положительного электрода
         ] = stateCoordinates

        # Получаем температуру
        [TInAkk,  # Температура содержимого аккумулятора
         TBAkk  # Температура корпуса аккумулятора
         ] = reducedTemp

        # Получаем параметры
        [I,  # Ток во внешней цепи
         UChMax,  # Граничное напряжение заряда
         bI0Ch,  # Граница зарядного нулевого тока
         Tokr,  # Температура окружающей среды
         EbinpC,  # ЭДС положительного двойного слоя в заряженном состоянии
         EbinnC,  # ЭДС отрицательного двойного слоя в заряженном состоянии
         EbinpD,  # ЭДС положительного двойного слоя в разряженном состоянии
         EbinnD,  # ЭДС отрицательного двойного слоя в разряженном состоянии
         Cbin0p,  # Емкость положительного двойного слоя
         Cm,  # Емкость мембраны
         Cbin0n,  # Емкость отрицательного двойного слоя
         Rbin0p,  # Сопротивление положительного двойного слоя
         Rm0,  # Сопротивление мембраны
         Rbin0n,  # Сопротивление отрицательного двойного слоя
         KInAkk,  # Коэффициент теплопередачи литий-ионного элемента
         CInAkk,  # Теплоемкость литий-ионного элемента
         KBAkk,  # Коэффициент теплоотдачи корпуса литий-ионного элемента
         CBAkk,  # Теплоемкость корпуса литий-ионного элемента
         Cnom,  # Номинальная емкость элемента
         rLiEpE,  # Приведенное зарядовое число положительного электрода
         rLiEnE,  # Приведенное зарядовое число отрицательного электрода
         alphaRIp,  # Коэффициент сопротивления по току положительного двойного слоя
         alphaRIn,  # Коэффициент сопротивления по току отрицательного двойного слоя
         alphaRQp,  # Коэффициент сопротивления по перенесенному через внешнюю цепь заряду положительного двойного слоя
         alphaRQn,  # Коэффициент сопротивления по перенесенному через внешнюю цепь заряду отрицательного двойного слоя
         nRQp,  # Степенной коэффициент сопротивления по перенесенному через внешнюю цепь заряду положительного двойного слоя
         nRQn,  # Степенной коэффициент сопротивления по перенесенному через внешнюю цепь заряду отрицательного двойного слоя
         alphaRTp,  # Экспоненциальный коэффициент сопротивления по температуре положительного электрода
         alphaRTm,  # Экспоненциальный коэффициент сопротивления по температуре мембраны
         alphaRTn,  # Экспоненциальный коэффициент сопротивления по температуре отрицательного электрода
         bRTp,  # Граничная температура по сопротивлению положительного электрода
         bRTm,  # Граничная температура по сопротивлению мембраны
         bRTn,  # Граничная температура по сопротивлению отрицательного электрода
         rCRTp,  # Постоянный коэффициент температурной зависимости положительного электрода
         rCRTm,  # Постоянный коэффициент температурной зависимости мембраны
         rCRTn,  # Постоянный коэффициент температурной зависимости отрицательного электрода
         alphaCQp,  # Зарядовый коэффициент емкости положительного электрода
         alphaCQn,  # Зарядовый коэффициент емкости отрицательного электрода
         qMatAllp,  # Общее зарядовое число молей материала положительного электрода
         qMatAlln,  # Общее зарядовое число молей материала отрицательного электрода
         muActsp,  # Зарядовый потенциал деградационной активации положительного электрода
         muActsn,  # Зарядовый потенциал деградационной активации отрицательного электрода
         bMuDegp,  # Зарядовый потенциал деградационного порога положительного электрода
         bMuDegn,  # Зарядовый потенциал деградационного порога отрицательного электрода
         kActElp,  # Коэффициент активации положительного электрода
         kActEln,  # Коэффициент активации отрицательного электрода
         kDegElp,  # Коэффициент деградации положительного электрода
         kDegEln,  # Коэффициент деградации отрицательного электрода
         aActElsp,  # Коэффиицент токовой активации положительного электрода
         aActElsn,  # Коэффиицент токовой активации отрицательного электрода
         bMuDegPosEl,  # Барьерный потенциал начала разрушения положительного электрода при переразрядке
         kDDegps,  # Коэффициент разрушения положительного электрода при перезарядке
         bDeltaTtoOkr,  # Граничная разность температур, при которой ухудшается теплообмен
         cDeltaTtoOkr,  # Коэффициент разности температур, при которой ухудшается теплообмен
         rNuLi0p,  # Относительное число молей лития в положительном электроде

         betaRI2p,
         betaRI2n,
         betaRI3p,
         betaRI3n,
         betaRQ2p,
         betaRQ2n,
         betaRQ3p,
         betaRQ3n,
         betaRT2p,
         betaRT2m,
         betaRT2n,
         betaRT3p,
         betaRT3m,
         betaRT3n,
         betaCQ2p,
         betaCQ2n,
         betaCQ3p,
         betaCQ3n,
         betaEQ2p,
         betaEQ2n,
         betaEQ3p,
         betaEQ3n,
         betaMuAct2p,
         betaMuAct2n,
         betaMuAct3p,
         betaMuAct3n,
         betaNonCeilp1,
         betaNonCeiln1,
         betaNonCeilp2,
         betaNonCeiln2,
         betaNonCeilp3,
         betaNonCeiln3,
         betaNonCeilQp1,
         betaNonCeilQn1,
         betaNonCeilQp2,
         betaNonCeilQn2,
         betaNonCeilQp3,
         betaNonCeilQn3,
         betaADNuMat1p,
         betaADNuMat1n,
         betaADNuMat2p,
         betaADNuMat2n,
         betaADNuMat3p,
         betaADNuMat3n,
         betaADNuMatDeg1p,
         betaADNuMatDeg1n,
         betaADNuMatDeg2p,
         betaADNuMatDeg2n,
         betaADNuMatDeg3p,
         betaADNuMatDeg3n,
         betaMuDegPos1,
         betaMuDegPos2,
         betaMuDegPos3,
         betaNuLiDegPos1,
         betaNuLiDegPos2,
         betaNuLiDegPos3,
         betaDeltaTtoOkr2,
         betaDeltaTtoOkr3,

         Rkl  # Сопротивление клемм
         ] = systemParameters

        # Рассчитываем приведенные (в единицах зарядовой емкости) числа молей интеркалированных в электроды ионов лития
        nuLip = qbinp + q  # Положительный электрод
        nuLin = qbinn + q  # Отрицательный электрод

        # Определяем ЭДС двойных слоев
        (Ebinp, Ebinn) = funEbin(EbinpC, EbinnC, EbinpD, EbinnD, nuLip,
                                 nuLin, rLiEpE, rLiEnE, Cnom,
                                 betaEQ2p, betaEQ2n, betaEQ3p, betaEQ3n)

        # Определяем емкости двойных слоев
        (Cbinp, Cbinn) = funCbin(qbinp, qbinn, alphaCQp, alphaCQn, Cbin0p, Cbin0n,
                                 betaCQ2p, betaCQ2n, betaCQ3p, betaCQ3n)

        # Рассчитываем напряжения двойных слоев
        Ubinp = qbinp / Cbinp  # Положительный двойной слой
        Um = qm / Cm  # Мембрана
        Ubinn = qbinn / Cbinn  # Отрицательный двойной слой
        self.__Ubinp = Ubinp
        self.__Um = Um
        self.__Ubinn = Ubinn

        # Напряжение на клеммах
        Uin = Ubinp + Um + Ubinn  # Напряжение внутри аккумулятора
        self.__Uin = Uin
        Ukl = Uin - I * Rkl

        # Определяем ток
        if (I < -bI0Ch) and (Ukl > UChMax):
            I = (Uin - UChMax) / Rkl
        self.__Icur = I

        # Внешний поток теплоты на корпус аккумулятора
        heatStreambEnPow = Rkl * np.power(I, 2)

        # Определяем падения напряжения на двойных слоях
        dissUbinp = Ebinp - Ubinp  # Положительный двойной слой
        dissUbinn = Ebinn - Ubinn  # Отрицательный двойной слой

        # Определяем числа и приведенные числа молей активированных материалов электродов
        (nuActElp, rNuActElp) = funNuEMat(qMatElp, qDegPosEl + qMatDegElp, qMatAllp)  # Положительный электрод
        (nuActEln, rNuActEln) = funNuEMat(qMatEln,             qMatDegEln, qMatAlln)  # Отрицательный электрод

        # Определем активационные и деградационные потенцилы электродов
        (muActp, muDegp) = funMuMat(nuActElp, rNuActElp,
                                    muActsp, bMuDegp,
                                    betaMuAct2p, betaMuAct3p)  # Положительный электрод
        (muActn, muDegn) = funMuMat(nuActEln, rNuActEln,
                                    muActsn, bMuDegn,
                                    betaMuAct2n, betaMuAct3n)  # Отрицательный электрод

        # Матрица Якоби приведенной энтропии по электрическим зарядам
        JSq = np.array([dissUbinp, -Um, dissUbinn, Ebinp + Ebinn,
                        muActp, muActn, muDegp, muDegn,
                        muActp - bMuDegPosEl], dtype=np.double) / TInAkk

        # Матрица Гесса приведенной энтропии по температуре и электрическим зарядам
        HSqT = np.vstack([-JSq / TInAkk,
                          np.zeros_like(JSq)])

        # Приведенные первые и вторые производные приведенной энтропии по температуре
        JST = np.array([CInAkk, CBAkk], dtype=np.double) / reducedTemp
        HSTT = -JST / reducedTemp

        # Вычисляем корректировочные коэффициенты по базовому сопротивлению двойных слоев
        corrRbinp = funCf0(qDegPosEl + qMatDegElp, qMatAllp,
                           betaNonCeilp1, betaNonCeilp2, betaNonCeilp3)  # Положительный электрод
        corrRbinn = funCf0(qMatDegEln, qMatAlln,
                           betaNonCeiln1, betaNonCeiln2, betaNonCeiln3)  # Отрицательный электрод

        # Вычисляем корректировочные коэффициенты по зарядовой емкости электродов
        corrQp = funCf0(qDegPosEl + qMatDegElp, qMatAllp,
                        betaNonCeilQp1, betaNonCeilQp2, betaNonCeilQp3)  # Положительный электрод
        corrQn = funCf0(qMatDegEln, qMatAlln,
                        betaNonCeilQn1, betaNonCeilQn2, betaNonCeilQn3)  # Отрицательный электрод

        # Определяем сопротивления двойных слоев и мембраны
        rAkk = funRbin(alphaRIp, alphaRIn, dissUbinp, dissUbinn,
                       nRQp, nRQn, alphaRQp / corrQp, alphaRQn / corrQn,
                       nuLip, nuLin, alphaRTp, alphaRTn, bRTp, bRTn,
                       rCRTp, rCRTn, alphaRTm, bRTm, rCRTm, TInAkk, Cnom,
                       Rbin0p * corrRbinp, Rbin0n * corrRbinn, Rm0, rNuLi0p,
                       betaRI2p, betaRI2n, betaRI3p, betaRI3n,
                       betaRQ2p, betaRQ2n, betaRQ3p, betaRQ3n,
                       betaRT2p, betaRT2m, betaRT2n, betaRT3p,
                       betaRT3m, betaRT3n) / (TInAkk / NonEqSystemQBase.GetTbase())
        [rbinp, rbinn, rm] = rAkk

        # Определям кинетические коэффициенты активации и деградации
        (kActp, kDegp) = funADNu(rNuActElp, TInAkk,
                                 muActp, muDegp,
                                 kActElp, kDegElp,
                                 betaADNuMat1p, betaADNuMat2p, betaADNuMat3p,
                                 betaADNuMatDeg1p, betaADNuMatDeg2p, betaADNuMatDeg3p)  # Положительный электрод
        (kActn, kDegn) = funADNu(rNuActEln, TInAkk,
                                 muActn, muDegn,
                                 kActEln, kDegEln,
                                 betaADNuMat1n, betaADNuMat2n, betaADNuMat3n,
                                 betaADNuMatDeg1n, betaADNuMatDeg2n, betaADNuMatDeg3n)  # Отрицательный электрод

        # Определяем коэффициент разрушения положительного электрода
        rMuDegPos = dissUbinp - bMuDegPosEl
        kDDegp = funKDDegPosEl(nuLip, rMuDegPos, TInAkk,
                               Cnom, bMuDegPosEl, kDDegps,
                               betaMuDegPos1, betaMuDegPos2, betaMuDegPos3,
                               betaNuLiDegPos1, betaNuLiDegPos2, betaNuLiDegPos3)

        # Коэффициент деградации электродов
        KDegEl = np.array([kDegp, kDegn, kDDegp], dtype=np.double)

        # Коэффициент активации электродов
        aActp = np.array([-aActElsp * np.sign(dissUbinp)], dtype=np.double)
        aActn = np.array([-aActElsn * np.sign(dissUbinn)], dtype=np.double)

        # Коэффициенты теплообмена
        cKbAkk = funCKbAkk(TBAkk, Tokr,
                           cDeltaTtoOkr, bDeltaTtoOkr,
                           betaDeltaTtoOkr2, betaDeltaTtoOkr3)
        KQAkk = np.array([KInAkk * TInAkk, KBAkk * cKbAkk * Tokr], dtype=np.double) * TBAkk / NonEqSystemQBase.GetTbase()

        # Выводим результат
        return (I, Tokr,
                heatStreambEnPow,
                JSq, JST, HSqT, HSTT,
                rbinp, rbinn, rm,
                aActp, aActn,
                kActp, kActn,
                KDegEl, KQAkk)

    # Выводим внутреннее напряжение
    def GetUin(self):
        return self.__Uin

    # Выводим ток
    def GetIcur(self):
        return self.__Icur

    # Выводим напряжение положительного двойного слоя
    def GetUbinp(self):
        return self.__Ubinp

    # Выводим напряжение мембраны
    def GetUm(self):
        return self.__Um

    # Выводим напряжение отрицательного двойного слоя
    def GetUbinn(self):
        return self.__Ubinn
