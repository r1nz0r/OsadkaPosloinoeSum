#pragma once
#include <vector>
struct SoilLayer;

std::vector<SoilLayer> GetSoilLayers();

#include "SoilParameters.h"
#include "GEO.h"

std::vector<SoilLayer> GetSoilLayers()
{
    static std::vector<SoilLayer> soilLayers
    {
        {9.8,18.0, 9.8, 22900, 114500},
        {6,19.5, 15.8, 20000, 100000},
        {100,19.3, 115.8, 36200, 181000},
        //{4.3, 21.14, 4.3, 30000, 109000}, // Толщина 4.3 м, γ=21.14 [кН/м3], E1=30000 [кПа], E2=109000 [кПа]
        //{100.0, 21.22, 104.3, 43000, 88000}, // Толщина 100 м, γ=21.22 [кН/м3], E1=43000 [кПа], E2=88000 [кПа]
    };

    return soilLayers;
}


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


#pragma once
#include "BuildingParameters.h"
#include <string>

// Вспомогательный константы
constexpr double BETTA = 0.8; // 5.6.31 SP 22.13330.2016
constexpr double DEFORMATION_MODULUS_SCALE_FACTOR = 5; // Примечание 1, раздел 5.6.31 SP 22.13330.2016
constexpr double DEFORMATION_MODULUS_BOUNDARY = 7000; // [kN/m2], 5.6.41 SP 22.13330.2016
constexpr double PI = 3.14159265;
constexpr double WATER_DENSITY = 10; // Плотность воды [кН/м3]

// Данные для расчетных слоев
constexpr double SUB_LAYER_THICKNESS_MANUAL = 0.01; // Желаемая толщина подслоя hi [м]
constexpr double SUB_LAYER_THICKNESS = std::min(SUB_LAYER_THICKNESS_MANUAL, 0.2 * B); // Толщина подслоя hi [м]. Принята равной SUB_LAYER_THICKNESS_MANUAL [м], если hi < 0.4b, иначе hi = 0.2b
constexpr double BACKFILL_DENSITY = 20; // Плотность грунта котлована [кН/м3]
constexpr double GROUND_WATER_PLATE_PRESSURE = std::max((FOUNDATION_DEPTH - GROUND_WATER_DEPTH) * WATER_DENSITY, 0.0); // Давление под подошвой фундамента с учетом грунтовой воды [кПа]

// Предварительный расчет напряжений в основании фундамента
constexpr double SOIL_INITIAL_STRESS = BACKFILL_DENSITY * FOUNDATION_DEPTH;

// Путь к файлу с результатами:
const std::string RESULT_PATH = R"(E:\WORK\PROJECTS\Praktika\TROSHIN\test.csv)";


#pragma once
#include <iostream>

// Ввод данных о фундаменте
// 
// SECTION 3 30UJZ
constexpr double Q = 357; // Нагрузка на подошву фундамента без учета грунтовой воды [кПа]
constexpr double B = 14.4; // Ширина подошвы фундамента [м]
constexpr double L = 66.7; // Длина фундамента [м]
constexpr double FOUNDATION_DEPTH = 7.8; // Глубина заложения фундамента [м]
constexpr double B_HALF = B / 2; // Половина ширины фундамента [м]
constexpr double L_HALF = L / 2; // Половина длины фундамента [м]
constexpr double GROUND_WATER_DEPTH = 1000; // Отметка верха водоносного слоя относительно земли [м]. Если ГВ нет - указать большое значение.


#include "GEO.h"
#include "Constants.h"
#include "BuildingParameters.h"
#include "SoilParameters.h"
#include <fstream>

void PrintResults(double settlement, const std::vector<CalculationLayer>& calculationLayers);

int main()
{
	setlocale(LC_ALL, "rus");

	std::vector<CalculationLayer> calculationLayers;
	ProcessSoilLayersCalculation(GetSoilLayers(), calculationLayers);
	double settlement = GetTotalSettlement(calculationLayers); // Расчетная осадка [м]
	PrintResults(settlement, calculationLayers);

	return 0;
}

void PrintResults(double settlement, const std::vector<CalculationLayer>& calculationLayers)
{
	std::cout << "\nTotal settlement: " << settlement * 1000 << " mm.\n\n-----------------------------------------\n";
	std::cout << "For more details see file: " << RESULT_PATH << std::endl;

	std::ofstream fileStream;
	fileStream.open(RESULT_PATH);

	if (fileStream.is_open())
	{
		fileStream << "INITIAL CALCULATION PARAMETERS:\n--------------------------------------\n";
		// Характеристики сооружения
		fileStream << "q = " << Q << " kPa\n";
		fileStream << "B = " << B << " m\n";
		fileStream << "L = " << L << " m\n";
		fileStream << "FOUNDATION_DEPTH = " << FOUNDATION_DEPTH << " m\n";
		fileStream << "GROUND_WATER_DEPTH = " << GROUND_WATER_DEPTH << " m\n";
		// Вспомогательные константы
		fileStream << "BETTA = " << BETTA << std::endl;
		fileStream << "DEFORMATION_MODULUS_SCALE_FACTOR = " << DEFORMATION_MODULUS_SCALE_FACTOR << std::endl;
		fileStream << "DEFORMATION_MODULUS_BOUNDARY = " << DEFORMATION_MODULUS_BOUNDARY << " kPa\n";
		fileStream << "WATER_DENSITY = " << WATER_DENSITY << " kPa\n";
		fileStream << "SUB_LAYER_THICKNESS_MANUAL = " << SUB_LAYER_THICKNESS_MANUAL << " m\n";
		fileStream << "SUB_LAYER_THICKNESS = " << SUB_LAYER_THICKNESS << " m\n";
		fileStream << "BACKFILL_DENSITY = " << BACKFILL_DENSITY << " kN/m3\n";
		fileStream << "GROUND_WATER_PLATE_PRESSURE = " << GROUND_WATER_PLATE_PRESSURE << " kPa\n";
		fileStream << "SOIL_INITIAL_STRESS = " << SOIL_INITIAL_STRESS << " kPa\n";
		fileStream << "-----------------------------------------" << std::endl;
		//Грунты
		fileStream << "INITIAL SOIL PARAMETERS:\n--------------------------------------\n";
		int soilIndex = 1;
		fileStream << "Index;Thickness, m;yf, kN/m3;yf_sat, kN/m3;E1, kPa;E2, kPa" << std::endl;

		for (const auto& soil : GetSoilLayers())
		{
			fileStream << soilIndex << ";" << soil.thickness << ";" << soil.unitWeight << ";" << soil.saturatedWeight
				<< ";" << soil.deformationModulusFirst << ";" << soil.deformationModulusSecond << std::endl;
			++soilIndex;
		}
		
		//Итоговая осадка основания
		fileStream << "\n-----------------------------------------\n" << std::endl;
		fileStream << "TOTAL DEFORMATION: " << settlement * 1000 << " mm.\n";
		fileStream << "\n-----------------------------------------\n" << std::endl;
		//Вывод всех параметров расчетных слоев
		fileStream << "Index;Layer Thickness, m;Bottom elevation, m;sigma_zg, kPa;sigma_red_factor;sigma_zg_red, kPa;alpha;sigma_zy, kPa;sigma_zp, kPa;s1, m;s2, m;s_tot, m" << std::endl;

		for (const auto& layer : calculationLayers)
		{
			fileStream << layer.index << ";" << layer.thickness << ";" << layer.bottomElevation << ";" << layer.soilStress
				<< ";" << layer.soilStressReductionFactor << ";" << layer.soilStressReducted << ";" << layer.alpha
				<< ";" << layer.foundationPitStress << ";" << layer.structureStress << ";" << layer.settlement1
				<< ";" << layer.settlement2 << ";" << layer.settlementTotal << std::endl;
		}
	}
}
