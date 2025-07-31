#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <numeric>

// Для работы с JSON используется библиотека nlohmann/json.
// Скачайте файл json.hpp и поместите его в папку с проектом.
#include "json.hpp"

using json = nlohmann::json;

// --- Структуры для хранения данных ---

// Входные данные по слою грунта (из JSON)
struct SoilLayerData {
    double thickness_m;
    double unitWeight_kN_m3;
    double saturatedUnitWeight_kN_m3;
    double deformationModulusFirst_kPa;
    double deformationModulusSecond_kPa;
};

// Все входные параметры, загружаемые из JSON
struct InputData {
    // Параметры фундамента
    double load_kPa;
    double width_m;
    double length_m;
    double foundationDepth_m;
    double groundWaterDepth_m;

    // Параметры расчета
    double betta_coefficient;
    double backfillDensity_kN_m3;
    double sublayerManualThickness_m;
    std::string resultFilePath;

    // Слои грунта
    std::vector<SoilLayerData> soilLayers;
};

// Результаты расчета для одного элементарного слоя (для вывода в CSV)
struct CalculationLayerResult {
    int index;
    const SoilLayerData& soil; // Ссылка на исходные свойства грунта
    double thickness;
    double bottomElevation;
    double soilStress;
    double soilStressReductionFactor;
    double soilStressReducted;
    double foundationPitStress;
    double structureStress;
    double alpha;
    double settlement1;
    double settlement2;
    double settlementTotal;
};

// --- Основной класс для выполнения расчетов ---

class SettlementCalculator {
public:
    explicit SettlementCalculator(InputData data) : m_input(std::move(data)) {}

    // Главный метод, запускающий вычисления
    void run() {
        // --- Инициализация состояния расчета ---
        std::vector<CalculationLayerResult> calculationResults;
        double depthSummary = 0.0;
        double soilStressSummary = m_input.backfillDensity_kN_m3 * m_input.foundationDepth_m;
        double totalSettlement = 0.0;
        double compressibleDepth = 0.0;

        // --- Основной цикл расчета по слоям (логика из ProcessSoilLayersCalculation) ---
        for (const auto& soil : m_input.soilLayers) {
            double remainingThickness = soil.thickness_m;
            bool stopConditionMet = false;

            while (remainingThickness > 0) {
                // Проверка условия остановки (нижняя граница сжимаемой толщи)
                if (!calculationResults.empty()) {
                    const auto& lastLayer = calculationResults.back();
                    if (lastLayer.soilStressReducted >= lastLayer.structureStress) {
                        compressibleDepth = depthSummary;
                        stopConditionMet = true;
                        break;
                    }
                }

                const double sublayerThickness = getSublayerThickness();
                const double thickness = std::min(sublayerThickness, remainingThickness);

                // Создаем и заполняем расчетный слой, полностью повторяя оригинальную логику
                CalculationLayerResult currentLayer = createCalculationLayer(
                    soil,
                    thickness,
                    depthSummary,
                    soilStressSummary,
                    calculationResults.empty() ? nullptr : &calculationResults.back()
                );

                // Обновляем состояние для следующего шага
                depthSummary += thickness;
                soilStressSummary = currentLayer.soilStress;
                totalSettlement += currentLayer.settlementTotal;

                calculationResults.push_back(currentLayer);
                remainingThickness -= thickness;
            }
            if (stopConditionMet) {
                break;
            }
        }

        // Если условие не сработало, сжимаемая толща равна всей пройденной глубине
        if (compressibleDepth == 0.0) {
            compressibleDepth = depthSummary;
        }

        // --- Вывод результатов ---
        printConsoleSummary(totalSettlement, compressibleDepth);
        writeCsvReport(calculationResults, totalSettlement, compressibleDepth);
    }

private:
    // Создание и расчет одного элементарного слоя
    CalculationLayerResult createCalculationLayer(
        const SoilLayerData& soil,
        double thickness,
        double currentDepthSummary,
        double currentSoilStressSummary,
        const CalculationLayerResult* prevLayer)
    {
        CalculationLayerResult layer{
            .index = static_cast<int>(prevLayer ? prevLayer->index + 1 : 1),
            .soil = soil,
            .thickness = thickness
        };

        // --- Расчет отметок (эквивалент CalculateElevations) ---
        layer.bottomElevation = currentDepthSummary + thickness;
        double middleElevation = currentDepthSummary + thickness / 2.0;

        // --- Расчет коэффициента alpha (эквивалент CalculateAlpha) ---
        layer.alpha = calculateAlpha(middleElevation);

        // --- Расчет напряжений (эквивалент CalculateStresses) ---
        double waterHeight = std::max((middleElevation + m_input.foundationDepth_m - m_input.groundWaterDepth_m), 0.0);
        double saturatedWeight = soil.unitWeight_kN_m3 - 10.0; // Как в оригинальном коде
        
        layer.soilStress = currentSoilStressSummary + thickness / 2.0 * (waterHeight > 0.0 ? saturatedWeight : soil.unitWeight_kN_m3);
        
        if (prevLayer) {
            double prevMiddleElevation = prevLayer->bottomElevation - prevLayer->thickness / 2.0;
            double prevWaterHeight = std::max((prevMiddleElevation + m_input.foundationDepth_m - m_input.groundWaterDepth_m), 0.0);
            double prevSaturatedWeight = prevLayer->soil.unitWeight_kN_m3 - 10.0;
            layer.soilStress += prevLayer->thickness / 2.0 * (prevWaterHeight > 0.0 ? prevSaturatedWeight : prevLayer->soil.unitWeight_kN_m3);
        }

        layer.soilStressReductionFactor = soil.deformationModulusFirst_kPa >= 7000 ? 0.5 : 0.2;
        layer.soilStressReducted = layer.soilStress * layer.soilStressReductionFactor;
        
        double initialStress = m_input.backfillDensity_kN_m3 * m_input.foundationDepth_m;
        layer.foundationPitStress = initialStress * layer.alpha;
        
        double waterPressure = std::max((m_input.foundationDepth_m - m_input.groundWaterDepth_m) * 10.0, 0.0);
        layer.structureStress = (m_input.load_kPa - waterPressure) * layer.alpha;

        // --- Расчет осадки (эквивалент CalculateSettlement) ---
        layer.settlement1 = m_input.betta_coefficient * (layer.structureStress - layer.foundationPitStress) * thickness / soil.deformationModulusFirst_kPa;
        layer.settlement2 = m_input.betta_coefficient * layer.foundationPitStress * thickness / soil.deformationModulusSecond_kPa;
        layer.settlementTotal = m_input.foundationDepth_m >= 5.0 ? layer.settlement1 + layer.settlement2 : layer.settlement1;

        return layer;
    }

    // Формула А. Лява для коэффициента alpha
    double calculateAlpha(double middleElevation) const {
        constexpr double PI = 3.14159265;
        double l_half = m_input.length_m / 2.0;
        double b_half = m_input.width_m / 2.0;
        double D = sqrt(l_half * l_half + b_half * b_half + middleElevation * middleElevation);
        if (D == 0.0) return 1.0;
        double denominator = (D * D * middleElevation * middleElevation + l_half * l_half * b_half * b_half);
        if (denominator == 0.0) return 1.0;

        double s1 = (l_half * b_half * middleElevation / D) * (l_half * l_half + b_half * b_half + 2 * middleElevation * middleElevation) / denominator;
        double s2_arg_den = (sqrt(l_half * l_half + middleElevation * middleElevation) * sqrt(b_half * b_half + middleElevation * middleElevation));
        if (s2_arg_den == 0.0) return (2 / PI) * s1;
        double s2_arg = (l_half * b_half) / s2_arg_den;
        if (s2_arg > 1.0) s2_arg = 1.0;
        if (s2_arg < -1.0) s2_arg = -1.0;

        double s2 = asin(s2_arg);
        return (2 / PI) * (s1 + s2);
    }

    // Определение толщины расчетного подслоя hi
    double getSublayerThickness() const {
        return std::min(m_input.sublayerManualThickness_m, 0.2 * m_input.width_m);
    }

    // Вывод итогов в консоль
    void printConsoleSummary(double totalSettlement, double compressibleDepth) const {
        std::cout << "\n--- РЕЗУЛЬТАТЫ РАСЧЕТА ---" << std::endl;
        std::cout << std::fixed;
        std::cout << "Итоговая осадка: " << std::setprecision(2) << totalSettlement * 1000 << " мм." << std::endl;
        std::cout << "Глубина сжимаемой толщи: " << std::setprecision(2) << compressibleDepth << " м." << std::endl;
        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "Подробный отчет сохранен в файл: " << m_input.resultFilePath << std::endl;
    }

    // Запись подробного отчета в CSV
    void writeCsvReport(const std::vector<CalculationLayerResult>& results, double totalSettlement, double compressibleDepth) const {
        std::ofstream file(m_input.resultFilePath);
        if (!file.is_open()) {
            std::cerr << "Ошибка: не удалось создать файл для записи результатов: " << m_input.resultFilePath << std::endl;
            return;
        }

        file << std::fixed << std::setprecision(4);
        
        file << "ИТОГОВЫЕ РЕЗУЛЬТАТЫ\n";
        file << "Итоговая осадка, мм;" << totalSettlement * 1000 << "\n";
        file << "Глубина сжимаемой толщи, м;" << compressibleDepth << "\n\n";

        file << "ИСХОДНЫЕ ДАННЫЕ\n";
        file << "Параметр;Значение;Ед. изм.\n";
        file << "q;" << m_input.load_kPa << ";кПа\n";
        file << "B;" << m_input.width_m << ";м\n";
        file << "L;" << m_input.length_m << ";м\n";
        file << "d;" << m_input.foundationDepth_m << ";м\n";
        file << "УГВ;" << m_input.groundWaterDepth_m << ";м\n\n";

        file << "ДЕТАЛЬНЫЙ РАСЧЕТ ПО СЛОЯМ\n";
        file << "Index;h_i, м;z_bot, м;sigma_zg, кПа;k;sigma_zg_red, кПа;alpha;sigma_zy, кПа;sigma_zp, кПа;s1, м;s2, м;s_tot, м\n";
        for (const auto& layer : results) {
            file << layer.index << ";" << layer.thickness << ";" << layer.bottomElevation << ";" << layer.soilStress << ";"
                 << layer.soilStressReductionFactor << ";" << layer.soilStressReducted << ";" << layer.alpha << ";"
                 << layer.foundationPitStress << ";" << layer.structureStress << ";" << layer.settlement1 << ";"
                 << layer.settlement2 << ";" << layer.settlementTotal << "\n";
        }
    }

    InputData m_input;
};

// Загрузка параметров из JSON файла
InputData loadInputFromJSON(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл с входными данными: " + filename);
    }
    json data = json::parse(file);
    
    InputData input;
    input.load_kPa = data["foundation"]["load_kPa"];
    input.width_m = data["foundation"]["width_m"];
    input.length_m = data["foundation"]["length_m"];
    input.foundationDepth_m = data["foundation"]["foundation_depth_m"];
    input.groundWaterDepth_m = data["foundation"]["ground_water_depth_m"];

    input.betta_coefficient = data["calculation_params"]["betta_coefficient"];
    input.backfillDensity_kN_m3 = data["calculation_params"]["backfill_density_kN_m3"];
    input.sublayerManualThickness_m = data["calculation_params"]["sublayer_manual_thickness_m"];
    input.resultFilePath = data["calculation_params"]["result_file_path"];
    
    for (const auto& layer_json : data["soil_layers"]) {
        // В оригинальном коде E2 вычислялся, если был 0. Повторяем эту логику.
        double e2 = layer_json.contains("E2_kPa") ? layer_json["E2_kPa"].get<double>() : 0.0;
        if (e2 == 0.0) {
            e2 = layer_json["E1_kPa"].get<double>() * 5.0; // DEFORMATION_MODULUS_SCALE_FACTOR
        }
        input.soilLayers.push_back({
            layer_json["thickness_m"],
            layer_json["unit_weight_kN_m3"],
            layer_json["saturated_unit_weight_kN_m3"],
            layer_json["E1_kPa"],
            e2
        });
    }
    return input;
}


int main() {
    setlocale(LC_ALL, "Russian");
    try {
        InputData input = loadInputFromJSON("input.json");
        SettlementCalculator calculator(input);
        calculator.run();
    } catch (const std::exception& e) {
        std::cerr << "Критическая ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}


{
  "foundation": {
    "load_kPa": 357,
    "width_m": 14.4,
    "length_m": 66.7,
    "foundation_depth_m": 7.8,
    "ground_water_depth_m": 1000
  },
  "calculation_params": {
    "betta_coefficient": 0.8,
    "backfill_density_kN_m3": 20.0,
    "sublayer_manual_thickness_m": 0.4,
    "result_file_path": "calculation_results.csv"
  },
  "soil_layers": [
    {
      "thickness_m": 9.8,
      "unit_weight_kN_m3": 18.0,
      "saturated_unit_weight_kN_m3": 18.0,
      "E1_kPa": 22900,
      "E2_kPa": 114500
    },
    {
      "thickness_m": 6.0,
      "unit_weight_kN_m3": 19.5,
      "saturated_unit_weight_kN_m3": 19.5,
      "E1_kPa": 20000,
      "E2_kPa": 100000
    },
    {
      "thickness_m": 100.0,
      "unit_weight_kN_m3": 19.3,
      "saturated_unit_weight_kN_m3": 19.3,
      "E1_kPa": 36200,
      "E2_kPa": 181000
    }
  ]
}

