from MathProtEnergyProc import standartIntegrateDyn


# Интегратор динамики в прямых задачах
integDynamic = standartIntegrateDyn(method="LSODA")


# Интегратор динамики в задачах оптимизации
integDynamicOptimize = standartIntegrateDyn(method="LSODA")
