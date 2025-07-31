#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include "Constants.h"

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ДАННЫЕ ПО СЛОЯМ ГРУНТОВ ВВОДИМ В ФУНКЦИИ GetSoilLayers() !!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// Структура для представления слоя грунта
struct SoilLayer
{
    double thickness; // Толщина слоя [м]
    double unitWeight; // Объёмный вес сухого грунта [кН/м3]
    double saturatedWeight; // Объемный вес водонасыщенного грунта [кН/м3]
    double deformationModulusFirst; // Модуль деформации для первичной кривой [МПа]
    double deformationModulusSecond; // Модуль деформации для вторичной кривой [МПа]
    double depth; // Глубина слоя относительно подошвы фундамента [м]

    SoilLayer(double _thickness, double _unitWeight, double _depth, double _deformationModulusFirst, double _deformationModulusSecond = 0)
        : thickness(_thickness),
        unitWeight(_unitWeight),
        saturatedWeight(_unitWeight - WATER_DENSITY),
        deformationModulusFirst(_deformationModulusFirst),
        deformationModulusSecond(_deformationModulusSecond == 0 ? _deformationModulusFirst * DEFORMATION_MODULUS_SCALE_FACTOR : _deformationModulusSecond),
        depth(_depth)
    {}
};


struct CalculationLayer
{
    const SoilLayer& soil;
    static int s_index;
    static double s_depthSummary;
    static double s_soilStressSummary;

    int index; // Порядковый номер слоя
    double thickness; // Толщина расчетного слоя [м]
    double middleElevation; // Отметка середины расчетного слоя относительно подошвы фундамета[м]
    double bottomElevation; // Отметка низа расчетного слоя относительно подошвы фундамента [м]
    double soilStress; // Среднее значение вертикального напряжения от собственного веса грунта [кПа]
    double soilStressReductionFactor; // Коэффициент перехода напряжений. 0.2 для E < 7 [МПа], 0.5 при E >= 7 [МПа].
    double soilStressReducted; // Среднее значение вертикльного напряжения от собственного веса грунта с учетом коэффициента 0.2 или 0.5 [кПа]
    double foundationPitStress; // Среднее значение вертикального напряжения от собственного веса выбранного при отрывке котлована грунта [кПа]
    double structureStress; // Среднее значение вертикального нормального напряжения от внешней нагрузки в грунте [кПа]
    double alpha; // Коэффициент "a" считаем по формуле (не по таблице)
    double settlement1; // Осадка слоя от левого слагаемого [м]
    double settlement2; // Осадка слоя от правого слагаемого [м]
    double settlementTotal; // Осадка суммарная [м]

    CalculationLayer(const SoilLayer& _soilLayer, double _thickness, const CalculationLayer* _prevLayer = nullptr);
    void CalculateSettlement();
    void CalculateStresses(const CalculationLayer* _prevLayer);
    void CalculateElevations();
    void CalculateAlpha();
    friend std::ostream& operator<<(std::ostream& s, const CalculationLayer& layer);
};

void ProcessSoilLayersCalculation(const std::vector<SoilLayer>& soilLayers, std::vector<CalculationLayer>& calculationLayers);
double GetTotalSettlement(const std::vector<CalculationLayer>& layers);
