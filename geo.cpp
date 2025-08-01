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
struct SoilLayerData
{
    double thickness_m;
    double unitWeight_kN_m3;
    double saturatedUnitWeight_kN_m3;
    double deformationModulusFirst_kPa;
    double deformationModulusSecond_kPa;
};

// Все входные параметры, загружаемые из JSON
struct InputData
{
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
struct CalculationLayerResult
{
    int index; // Порядковый номер слоя
    const SoilLayerData& soil; // Ссылка на исходные свойства грунта
    double thickness; // Толщина расчетного слоя [м]
    double bottomElevation; // Отметка низа расчетного слоя относительно подошвы фундамента [м]
    double soilStress;  // Среднее значение вертикального напряжения от собственного веса грунта [кПа]
    double soilStressReductionFactor; // Коэффициент перехода напряжений. 0.2 для E < 7 [МПа], 0.5 при E >= 7 [МПа].
    double soilStressReducted; // Среднее значение вертикльного напряжения от собственного веса грунта с учетом коэффициента 0.2 или 0.5 [кПа]
    double foundationPitStress; // Среднее значение вертикального напряжения от собственного веса выбранного при отрывке котлована грунта [кПа]
    double structureStress; // Среднее значение вертикального нормального напряжения от внешней нагрузки в грунте [кПа]
    double alpha; // Коэффициент "a" считаем по формуле (не по таблице)
    double settlement1; // Осадка слоя от левого слагаемого [м]
    double settlement2; // Осадка слоя от правого слагаемого [м]
    double settlementTotal; // Осадка суммарная [м]
};

// --- Основной класс для выполнения расчетов ---

class SettlementCalculator
{
public:
    explicit SettlementCalculator(InputData data) : m_input(std::move(data)) {}
    void Run();

private:
    // Создание и расчет одного элементарного слоя
    CalculationLayerResult CreateCalculationLayer(
        const SoilLayerData& soil,
        double thickness,
        double currentDepthSummary,
        double currentSoilStressSummary,
        const CalculationLayerResult* prevLayer);    

    // Формула А. Лява для коэффициента alpha
    double CalculateAlpha(double middleElevation) const;

    // Вывод итогов в консоль
    void PrintConsoleSummary(double totalSettlement, double compressibleDepth) const;

    // Запись подробного отчета в CSV
    void WriteCsvReport(const std::vector<CalculationLayerResult>& results, double totalSettlement, double compressibleDepth) const;
    
    // Определение толщины расчетного подслоя hi
    inline double GetSublayerThickness() const
    {
        return std::min(m_input.sublayerManualThickness_m, 0.2 * m_input.width_m);
    }

private:
    InputData m_input;
};

// Загрузка параметров из JSON файла
InputData LoadInputFromJSON(const std::string& filename)
{
    std::ifstream file(filename);

    if (!file.is_open())
    {
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

    for (const auto& layer_json : data["soil_layers"])
    {        
        double e2 = layer_json.contains("E2_kPa") ? layer_json["E2_kPa"].get<double>() : 0.0;
       
        if (e2 == 0.0)
        {
            e2 = layer_json["E1_kPa"].get<double>() * 5.0; // По СП 22.13330 если нет данных то E2 = 5*E1
        
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

std::string GetInputFilePath()
{
    std::string path;
    std::cout << "Введите путь к файлу данных и нажмите Enter." << std::endl;
    std::cout << "(Если файл (input.json) в той же папке, что и программа, просто нажмите Enter): ";
    std::getline(std::cin, path);

    // Если пользователь ничего не ввел, используем имя файла по умолчанию
    if (path.empty())
    {
        return "input.json";
    }

    // Простая проверка, чтобы убрать случайные кавычки, если пользователь скопировал путь
    if (path.front() == '"' && path.back() == '"')
    {
        path = path.substr(1, path.length() - 2);
    }

    return path;
}

void SettlementCalculator::Run()
{
    // --- Инициализация состояния расчета ---
    std::vector<CalculationLayerResult> calculationResults;
    double depthSummary = 0.0;
    double soilStressSummary = m_input.backfillDensity_kN_m3 * m_input.foundationDepth_m;
    double totalSettlement = 0.0;
    double compressibleDepth = 0.0;

    // --- Основной цикл расчета по слоям ---
    for (const auto& soil : m_input.soilLayers)
    {
        double remainingThickness = soil.thickness_m;
        bool stopConditionMet = false;

        while (remainingThickness > 0)
        {
            // Проверка условия остановки (нижняя граница сжимаемой толщи)
            if (!calculationResults.empty())
            {
                const auto& lastLayer = calculationResults.back();
                if (lastLayer.soilStressReducted >= lastLayer.structureStress)
                {
                    compressibleDepth = depthSummary;
                    stopConditionMet = true;
                    break;
                }
            }

            const double sublayerThickness = GetSublayerThickness();
            const double thickness = std::min(sublayerThickness, remainingThickness);

            // Создаем и заполняем расчетный слой
            CalculationLayerResult currentLayer = CreateCalculationLayer(
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
        if (stopConditionMet)
        {
            break;
        }
    }

    // Если условие не сработало, сжимаемая толща равна всей пройденной глубине
    if (compressibleDepth == 0.0)
    {
        compressibleDepth = depthSummary;
    }

    // --- Вывод результатов ---
    PrintConsoleSummary(totalSettlement, compressibleDepth);
    WriteCsvReport(calculationResults, totalSettlement, compressibleDepth);

    // --- Дать пользователю возможность увидеть результаты ---
    std::cout << "Нажмите Enter для выхода...";
    std::cin.get();
}

CalculationLayerResult SettlementCalculator::CreateCalculationLayer(
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

    // --- Расчет отметок
    layer.bottomElevation = currentDepthSummary + thickness;
    double middleElevation = currentDepthSummary + thickness / 2.0;

    // --- Расчет коэффициента alpha--
    layer.alpha = CalculateAlpha(middleElevation);

    // --- Расчет напряжений
    double waterHeight = std::max((middleElevation + m_input.foundationDepth_m - m_input.groundWaterDepth_m), 0.0);
    double saturatedWeight = soil.unitWeight_kN_m3 - 10.0;

    layer.soilStress = currentSoilStressSummary + thickness / 2.0 * (waterHeight > 0.0 ? saturatedWeight : soil.unitWeight_kN_m3);

    if (prevLayer)
    {
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

    // --- Расчет осадки
    layer.settlement1 = m_input.betta_coefficient * (layer.structureStress - layer.foundationPitStress) * thickness / soil.deformationModulusFirst_kPa;
    layer.settlement2 = m_input.betta_coefficient * layer.foundationPitStress * thickness / soil.deformationModulusSecond_kPa;
    layer.settlementTotal = m_input.foundationDepth_m >= 5.0 ? layer.settlement1 + layer.settlement2 : layer.settlement1;

    return layer;
}

double SettlementCalculator::CalculateAlpha(double middleElevation) const
{
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
    s2_arg = std::clamp(s2_arg, -1.0, 1.0);

    double s2 = asin(s2_arg);
    return (2 / PI) * (s1 + s2);
}

void SettlementCalculator::PrintConsoleSummary(double totalSettlement, double compressibleDepth) const
{
    std::cout << "\n--- РЕЗУЛЬТАТЫ РАСЧЕТА ---" << std::endl;
    std::cout << std::fixed;
    std::cout << "Осадка: " << std::setprecision(2) << totalSettlement * 1000 << " мм." << std::endl;
    std::cout << "Глубина сжимаемой толщи: " << std::setprecision(2) << compressibleDepth << " м." << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "Подробный отчет сохранен в файл: " << m_input.resultFilePath << std::endl;
}

void SettlementCalculator::WriteCsvReport(const std::vector<CalculationLayerResult>& results, double totalSettlement, double compressibleDepth) const
{
    std::ofstream file(m_input.resultFilePath);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось создать файл для записи результатов: " << m_input.resultFilePath << std::endl;
        return;
    }

    file << std::fixed << std::setprecision(4);

    file << "ИТОГОВЫЕ РЕЗУЛЬТАТЫ\n";
    file << "Осадка, мм;" << totalSettlement * 1000 << "\n";
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

    for (const auto& layer : results)
    {
        file << layer.index << ";" << layer.thickness << ";" << layer.bottomElevation << ";" << layer.soilStress << ";"
            << layer.soilStressReductionFactor << ";" << layer.soilStressReducted << ";" << layer.alpha << ";"
            << layer.foundationPitStress << ";" << layer.structureStress << ";" << layer.settlement1 << ";"
            << layer.settlement2 << ";" << layer.settlementTotal << "\n";
    }
}

int main()
{
    setlocale(LC_ALL, "Russian");
    try
    {
        std::string config_path = GetInputFilePath();
        InputData input = LoadInputFromJSON(config_path);
        SettlementCalculator calculator(input);
        calculator.Run();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nКритическая ошибка: " << e.what() << std::endl;
        std::cout << "Нажмите Enter для выхода...";
        std::cin.get();
        return 1;
    }

    return 0;
}
