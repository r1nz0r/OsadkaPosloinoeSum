#include "GEO.h"

double CalculationLayer::s_depthSummary = 0.0;
double CalculationLayer::s_soilStressSummary = SOIL_INITIAL_STRESS;
int CalculationLayer::s_index = 1;


CalculationLayer::CalculationLayer(const SoilLayer& _soilLayer, double _thickness, const CalculationLayer* _prevLayer)
    : soil(_soilLayer),
    index(s_index++),
    thickness(_thickness),
    middleElevation(0.0),
    bottomElevation(0.0),
    soilStress(0.0),
    foundationPitStress(0.0),
    structureStress(0.0),
    alpha(0.0),
    settlement1(0.0),
    settlement2(0.0),
    settlementTotal(0.0)
{
    CalculateElevations();
    CalculateAlpha();
    CalculateStresses(_prevLayer);
    CalculateSettlement();
}


void CalculationLayer::CalculateSettlement()
{
    settlement1 = BETTA * (structureStress - foundationPitStress) * thickness / soil.deformationModulusFirst;
    settlement2 = BETTA * foundationPitStress * thickness / soil.deformationModulusSecond;
    settlementTotal = FOUNDATION_DEPTH >= 5.0 ? settlement1 + settlement2 : settlement1; // Если глубина заложения ФП более 5 [м] то учитываем второе слагаемое.
}

void CalculationLayer::CalculateStresses(const CalculationLayer* _prevLayer)
{
    double waterHeight = std::max((middleElevation + FOUNDATION_DEPTH - GROUND_WATER_DEPTH), -0.001);
    soilStress = s_soilStressSummary + thickness / 2 * (waterHeight > 0.0 ? soil.saturatedWeight : soil.unitWeight); // Если слой ниже УГВ то учитываем это в расчете водонасыщенным весом грунта.

    if (_prevLayer)
    {
        waterHeight = std::max((_prevLayer->middleElevation + FOUNDATION_DEPTH - GROUND_WATER_DEPTH), -0.001);
        soilStress += _prevLayer->thickness / 2 * (waterHeight > 0.0 ? _prevLayer->soil.saturatedWeight : _prevLayer->soil.unitWeight);
    }

    s_soilStressSummary = soilStress;
    soilStressReductionFactor = soil.deformationModulusFirst >= 7000 ? 0.5 : 0.2;
    soilStressReducted = soilStress * soilStressReductionFactor;
    foundationPitStress = SOIL_INITIAL_STRESS * alpha;

    structureStress = (Q - GROUND_WATER_PLATE_PRESSURE) * alpha;

    if (soilStress < 0)
    {
        std::cerr << "Ошибка схождения данных. Отрицательные напряжения в расчетном слое грунта.\n";
        std::cout << *this;
        std::terminate();
    }
}

void CalculationLayer::CalculateElevations()
{
    middleElevation = s_depthSummary + thickness / 2.0;
    bottomElevation = s_depthSummary + thickness;
    s_depthSummary += thickness;
}

// Формула А. Лява для величины сжимающих напряжений. Альтернативный вариант - интерполяция таблицы из СП 22.13330.
void CalculationLayer::CalculateAlpha()
{
    double D = sqrt(L_HALF * L_HALF + B_HALF * B_HALF + middleElevation * middleElevation);
    double s1 = (L_HALF * B_HALF * middleElevation / D) * (L_HALF * L_HALF + B_HALF * B_HALF + 2 * middleElevation * middleElevation) / (D * D * middleElevation * middleElevation + L_HALF * L_HALF * B_HALF * B_HALF);
    double s2 = asin((L_HALF * B_HALF) / (sqrt(L_HALF * L_HALF + middleElevation * middleElevation) * sqrt(B_HALF * B_HALF + middleElevation * middleElevation)));
    alpha = (2 / PI) * (s1 + s2);
}

std::ostream& operator<<(std::ostream& s, const CalculationLayer& layer)
{
    //s << "Layer index: " << layer.index << std::endl;
    //s << "Layer thickness: " << layer.thickness << " m." << std::endl;
    //s << "Layer bottom elevation: " << layer.bottomElevation << " m." << std::endl;
    //s << "Layer soil stress: " << layer.soilStress << " kPa." << std::endl;
    //s << "Layer soil stress reduction factor: " << layer.soilStressReductionFactor << std::endl;
    //s << "Layer soil stress reducted: " << layer.soilStressReducted << " kPa." << std::endl;
    //s << "Layer foundation pit stress: " << layer.foundationPitStress << " kPa." << std::endl;
    //s << "Layer sctructure stress: " << layer.structureStress << " kPa." << std::endl;
    //s << "Layer alpha factor: " << layer.alpha << " kPa." << std::endl;
    //s << "Layer settlement 1: " << layer.settlement1 << " m." << std::endl;
    //s << "Layer settlement 2: " << layer.settlement2 << " m." << std::endl;
    //s << "Layer settlement summary: " << layer.settlementTotal << " m." << std::endl;
    //s << "---------------------------------------------" << std::endl;

    s << layer.index << ";" << layer.thickness << ";" << layer.bottomElevation << ";" << layer.soilStress << ";" << layer.soilStressReductionFactor << ";"
        << layer.soilStressReducted << ";" << layer.alpha << ";" << layer.foundationPitStress << ";" << layer.structureStress << ";"
        << layer.settlement1 << ";" << layer.settlement2 << ";" << layer.settlementTotal << ";" << std::endl;
    return s;
}



void ProcessSoilLayersCalculation(const std::vector<SoilLayer>& soilLayers, std::vector<CalculationLayer>& calculationLayers)
{
    for (const auto& soil : soilLayers)
    {
        double remainingThickness = soil.thickness;

        while (remainingThickness > 0)
        {
            if (calculationLayers.size() != 0 && calculationLayers.back().soilStressReducted >= calculationLayers.back().structureStress)
                break;

            double thickness = std::min(SUB_LAYER_THICKNESS, remainingThickness);
            CalculationLayer* previousLayer = calculationLayers.size() < 1 ? nullptr : &calculationLayers[calculationLayers.size() - 1];
            calculationLayers.push_back(CalculationLayer(soil, thickness, previousLayer));
            remainingThickness -= thickness;
        }
    }
}

double GetTotalSettlement(const std::vector<CalculationLayer>& layers)
{
    double settlement = 0.0;

    for (const auto& layer : layers)
        settlement += layer.settlementTotal;

    return settlement;
}
