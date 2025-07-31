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
