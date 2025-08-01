#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <numeric>
#include <algorithm> // для std::min, std::max
#include <utility>   // для std::move

// Для работы с JSON используется библиотека nlohmann/json.
// Скачайте файл json.hpp и поместите его в папку с проектом.
#include "json.hpp"

using json = nlohmann::json;

// --- Структуры для хранения данных ---

struct SoilLayerData
{
    double thickness_m;
    double unitWeight_kN_m3;
    double saturatedUnitWeight_kN_m3;
    double deformationModulusFirst_kPa;
    double deformationModulusSecond_kPa;
};

struct InputData
{
    double load_kPa;
    double width_m;
    double length_m;
    double foundationDepth_m;
    double groundWaterDepth_m;
    double betta_coefficient;
    double backfillDensity_kN_m3;
    double sublayerManualThickness_m;
    std::string resultFilePath;
    std::vector<SoilLayerData> soilLayers;
};

struct CalculationLayerResult
{
    int index;
    const SoilLayerData& soil;
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

class SettlementCalculator
{
public:
    explicit SettlementCalculator(InputData data) : m_input(std::move(data)) {}
    void Run();

private:
    CalculationLayerResult CreateCalculationLayer(
        const SoilLayerData& soil,
        double thickness,
        double currentDepthSummary,
        double currentSoilStressSummary,
        const CalculationLayerResult* prevLayer);

    double CalculateAlpha(double middleElevation) const;
    void PrintConsoleSummary(double totalSettlement, double compressibleDepth) const;
    void WriteCsvReport(const std::vector<CalculationLayerResult>& results, double totalSettlement, double compressibleDepth) const;
    void WriteHtmlReport(const std::vector<CalculationLayerResult>& results, double compressibleDepth) const;

    inline double GetSublayerThickness() const
    {
        return std::min(m_input.sublayerManualThickness_m, 0.4 * m_input.width_m);
    }

private:
    InputData m_input;
};

// --- Реализация методов класса ---

void SettlementCalculator::Run()
{
    std::vector<CalculationLayerResult> calculationResults;
    double depthSummary = 0.0;
    double soilStressSummary = m_input.backfillDensity_kN_m3 * m_input.foundationDepth_m;
    double totalSettlement = 0.0;
    double compressibleDepth = 0.0;

    for (const auto& soil : m_input.soilLayers)
    {
        double remainingThickness = soil.thickness_m;
        bool stopConditionMet = false;

        while (remainingThickness > 0)
        {
            if (!calculationResults.empty())
            {
                const auto& lastLayer = calculationResults.back();
                if (lastLayer.structureStress <= lastLayer.soilStressReducted)
                {
                    compressibleDepth = depthSummary;
                    stopConditionMet = true;
                    break;
                }
            }

            const double sublayerThickness = GetSublayerThickness();
            const double thickness = std::min(sublayerThickness, remainingThickness);

            CalculationLayerResult currentLayer = CreateCalculationLayer(
                soil,
                thickness,
                depthSummary,
                soilStressSummary,
                calculationResults.empty() ? nullptr : &calculationResults.back());

            depthSummary += thickness;
            double middleOfLayer = currentLayer.bottomElevation - currentLayer.thickness / 2.0;
            double waterHeight = std::max((middleOfLayer + m_input.foundationDepth_m - m_input.groundWaterDepth_m), 0.0);
            double currentUnitWeight = (waterHeight > 0.0) ? soil.saturatedUnitWeight_kN_m3 : soil.unitWeight_kN_m3;
            soilStressSummary += (currentUnitWeight - (waterHeight > 0 ? 10.0 : 0.0)) * thickness;
            
            totalSettlement += currentLayer.settlementTotal;

            calculationResults.push_back(currentLayer);
            remainingThickness -= thickness;
        }
        if (stopConditionMet)
        {
            break;
        }
    }

    if (compressibleDepth == 0.0)
    {
        compressibleDepth = depthSummary;
    }

    PrintConsoleSummary(totalSettlement, compressibleDepth);
    WriteCsvReport(calculationResults, totalSettlement, compressibleDepth);
    WriteHtmlReport(calculationResults, compressibleDepth);

    std::cout << u8"\nРабота программы успешно завершена. Нажмите Enter для выхода...";
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
        .thickness = thickness};

    layer.bottomElevation = currentDepthSummary + thickness;
    double middleElevation = currentDepthSummary + thickness / 2.0;

    layer.alpha = CalculateAlpha(middleElevation);

    double waterHeight = std::max((middleElevation + m_input.foundationDepth_m - m_input.groundWaterDepth_m), 0.0);
    double currentUnitWeight = (waterHeight > 0.0) ? soil.saturatedUnitWeight_kN_m3 : soil.unitWeight_kN_m3;
    layer.soilStress = currentSoilStressSummary + (currentUnitWeight - (waterHeight > 0 ? 10.0 : 0.0)) * thickness / 2.0;

    layer.soilStressReductionFactor = soil.deformationModulusFirst_kPa >= 7000 ? 0.5 : 0.2;
    layer.soilStressReducted = layer.soilStress * layer.soilStressReductionFactor;

    double initialStress = m_input.backfillDensity_kN_m3 * m_input.foundationDepth_m;
    layer.foundationPitStress = initialStress * layer.alpha;

    double p0 = m_input.load_kPa - initialStress;
    layer.structureStress = p0 * layer.alpha;

    layer.settlement1 = m_input.betta_coefficient * (layer.structureStress - layer.foundationPitStress) * thickness / soil.deformationModulusFirst_kPa;
    layer.settlement2 = m_input.betta_coefficient * layer.foundationPitStress * thickness / soil.deformationModulusSecond_kPa;
    layer.settlementTotal = m_input.foundationDepth_m >= 5.0 ? layer.settlement1 + layer.settlement2 : layer.settlement1;

    return layer;
}

double SettlementCalculator::CalculateAlpha(double middleElevation) const
{
    if (middleElevation < 1e-6) middleElevation = 1e-6;

    constexpr double PI = 3.141592653589793;
    double l_half = m_input.length_m / 2.0;
    double b_half = m_input.width_m / 2.0;
    double z = middleElevation;

    double m = l_half / b_half;
    double n = z / b_half;

    double term1_num = m * n;
    double term1_den = sqrt(1 + m*m + n*n);
    double term1_mult_num = 1;
    double term1_mult_den1 = m*m + n*n;
    double term1_mult_den2 = 1 + n*n;

    if (term1_den == 0 || term1_mult_den1 == 0 || term1_mult_den2 == 0) return 0;

    double term1 = (term1_num / term1_den) * ( (term1_mult_num/term1_mult_den1) + (term1_mult_num/term1_mult_den2) );
    
    double term2_arg_den = sqrt(n*n+1) * sqrt(m*m+n*n);
    if (term2_arg_den == 0) return 0;
    double term2_arg = m / term2_arg_den;
    
    term2_arg = std::clamp(term2_arg, -1.0, 1.0);
    double term2 = atan(term2_arg);

    return (term1 + term2) / PI;
}

void SettlementCalculator::PrintConsoleSummary(double totalSettlement, double compressibleDepth) const
{
    std::cout << u8"\n--- РЕЗУЛЬТАТЫ РАСЧЕТА ---" << std::endl;
    std::cout << std::fixed;
    std::cout << u8"Итоговая осадка: " << std::setprecision(2) << totalSettlement * 1000 << u8" мм." << std::endl;
    std::cout << u8"Глубина сжимаемой толщи: " << std::setprecision(2) << compressibleDepth << u8" м." << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << u8"Подробный CSV отчет сохранен в файл: " << m_input.resultFilePath << std::endl;
    std::cout << u8"Интерактивный HTML отчет сохранен в файл: calculation_report.html" << std::endl;
}

void SettlementCalculator::WriteCsvReport(const std::vector<CalculationLayerResult>& results, double totalSettlement, double compressibleDepth) const
{
    std::ofstream file(m_input.resultFilePath);
    if (!file.is_open())
    {
        std::cerr << u8"Ошибка: не удалось создать CSV файл для записи результатов." << std::endl;
        return;
    }
    
    // Пишем BOM для UTF-8, чтобы Excel корректно распознал кириллицу
    file << (char)0xEF << (char)0xBB << (char)0xBF;

    file << std::fixed << std::setprecision(4);
    file << u8"ИТОГОВЫЕ РЕЗУЛЬТАТЫ\n";
    file << u8"Осадка, мм;" << totalSettlement * 1000 << "\n";
    file << u8"Глубина сжимаемой толщи, м;" << compressibleDepth << "\n\n";
    file << u8"ДЕТАЛЬНЫЙ РАСЧЕТ ПО СЛОЯМ\n";
    file << u8"Index;h_i, м;z_bot, м;sigma_zg, кПа;k;sigma_zg_red, кПа;alpha;sigma_zy, кПа;sigma_zp, кПа;s1, м;s2, м;s_tot, м\n";
    for (const auto& layer : results)
    {
        file << layer.index << ";" << layer.thickness << ";" << layer.bottomElevation << ";" << layer.soilStress << ";"
             << layer.soilStressReductionFactor << ";" << layer.soilStressReducted << ";" << layer.alpha << ";"
             << layer.foundationPitStress << ";" << layer.structureStress << ";" << layer.settlement1 << ";"
             << layer.settlement2 << ";" << layer.settlementTotal << "\n";
    }
}

void SettlementCalculator::WriteHtmlReport(const std::vector<CalculationLayerResult>& results, double compressibleDepth) const
{
    std::ofstream file("calculation_report.html");
    if (!file.is_open())
    {
        std::cerr << u8"Ошибка: не удалось создать HTML-отчет." << std::endl;
        return;
    }

    // --- ИСПРАВЛЕНИЕ: Пишем UTF-8 BOM в начало файла ---
    file << (char)0xEF << (char)0xBB << (char)0xBF;

    // Используем u8"" для всех строк с кириллицей
    file << u8R"(<!DOCTYPE html>
<html lang="ru">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Отчет по расчету осадки</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-annotation@3.0.1/dist/chartjs-plugin-annotation.min.js"></script>
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif; margin: 0; background-color: #f4f7fa; color: #333; }
        .container { max-width: 1200px; margin: 20px auto; padding: 20px; background-color: #fff; box-shadow: 0 4px 8px rgba(0,0,0,0.1); border-radius: 8px; }
        h1, h2 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        .chart-container { position: relative; margin: 20px 0; padding: 20px; border: 1px solid #ddd; border-radius: 8px; }
        .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-top: 20px; }
        @media (max-width: 900px) { .grid { grid-template-columns: 1fr; } }
    </style>
</head>
<body>
<div class="container">
    <h1>Результаты расчета осадки</h1>
    <div class="grid">
        <div class="chart-container">
            <h2>Эпюры напряжений</h2>
            <canvas id="stressChart"></canvas>
        </div>
        <div class="chart-container">
            <h2>График осадки</h2>
            <canvas id="settlementChart"></canvas>
        </div>
    </div>
    <script>
)";
    file << std::fixed << std::setprecision(4);

    double p0 = m_input.load_kPa - (m_input.backfillDensity_kN_m3 * m_input.foundationDepth_m);
    double initial_sigma_zp = p0 * CalculateAlpha(0.001);
    double initial_sigma_zg = m_input.backfillDensity_kN_m3 * m_input.foundationDepth_m;
    double initial_sigma_zg_red = initial_sigma_zg * (results.empty() ? 0.2 : results.front().soil.deformationModulusFirst_kPa >= 7000 ? 0.5 : 0.2);

    file << "const sigma_zp_data = [{x: " << initial_sigma_zp << ", y: 0}";
    for (const auto& r : results) { file << ", {x: " << r.structureStress << ", y: " << r.bottomElevation << "}"; }
    file << "];\n";

    file << "const sigma_zg_red_data = [{x: " << initial_sigma_zg_red << ", y: 0}";
    for (const auto& r : results) { file << ", {x: " << r.soilStressReducted << ", y: " << r.bottomElevation << "}"; }
    file << "];\n";

    file << "const settlement_data = [{x: 0, y: 0}";
    double cumulativeSettlement = 0.0;
    for (const auto& r : results) {
        cumulativeSettlement += r.settlementTotal;
        file << ", {x: " << cumulativeSettlement * 1000 << ", y: " << r.bottomElevation << "}";
    }
    file << "];\n";
    
    file << "const compressibleDepth = " << compressibleDepth << ";\n";

    file << u8R"(
        Chart.register(ChartAnnotation);
        const stressCtx = document.getElementById('stressChart').getContext('2d');
        new Chart(stressCtx, {
            type: 'scatter',
            data: {
                datasets: [{
                    label: 'σ_zp (доп. от нагрузки)',
                    data: sigma_zp_data,
                    borderColor: 'rgb(255, 99, 132)',
                    backgroundColor: 'rgba(255, 99, 132, 0.5)',
                    showLine: true,
                    tension: 0.1
                }, {
                    label: 'σ_zg,red (бытовое с коэф.)',
                    data: sigma_zg_red_data,
                    borderColor: 'rgb(54, 162, 235)',
                    backgroundColor: 'rgba(54, 162, 235, 0.5)',
                    showLine: true,
                    tension: 0.1
                }]
            },
            options: {
                scales: {
                    y: { reverse: true, title: { display: true, text: 'Глубина, м' } },
                    x: { position: 'top', title: { display: true, text: 'Напряжение, кПа' }, min: 0 }
                },
                plugins: {
                    annotation: {
                        annotations: {
                            line1: {
                                type: 'line',
                                yMin: compressibleDepth,
                                yMax: compressibleDepth,
                                borderColor: 'rgb(75, 192, 192)',
                                borderWidth: 2,
                                borderDash: [6, 6],
                                label: {
                                    content: `Hс = ${compressibleDepth.toFixed(2)} м`,
                                    display: true,
                                    position: 'start',
                                    backgroundColor: 'rgba(75, 192, 192, 0.8)'
                                }
                            }
                        }
                    }
                }
            }
        });

        const settlementCtx = document.getElementById('settlementChart').getContext('2d');
        new Chart(settlementCtx, {
            type: 'scatter',
            data: {
                datasets: [{
                    label: 'Накопленная осадка',
                    data: settlement_data,
                    borderColor: 'rgb(255, 159, 64)',
                    backgroundColor: 'rgba(255, 159, 64, 0.5)',
                    showLine: true,
                    tension: 0.1
                }]
            },
            options: {
                scales: {
                     y: { reverse: true, title: { display: true, text: 'Глубина, м' } },
                     x: { position: 'top', title: { display: true, text: 'Осадка, мм' }, min: 0 }
                }
            }
        });
    </script>
</div>
</body>
</html>
)";
}

// --- Глобальные функции и main ---

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
            e2 = layer_json["E1_kPa"].get<double>() * 5.0;
        }
        input.soilLayers.push_back({
            layer_json["thickness_m"],
            layer_json["unit_weight_kN_m3"],
            layer_json["saturated_unit_weight_kN_m3"],
            layer_json["E1_kPa"],
            e2});
    }
    return input;
}

std::string GetInputFilePath()
{
    std::string path;
    std::cout << u8"Введите путь к файлу данных и нажмите Enter." << std::endl;
    std::cout << u8"(Если файл (input.json) в той же папке, что и программа, просто нажмите Enter): ";
    std::getline(std::cin, path);

    if (path.empty())
    {
        return "input.json";
    }

    if (path.front() == '"' && path.back() == '"')
    {
        path = path.substr(1, path.length() - 2);
    }
    return path;
}

int main()
{
    // setlocale больше не нужен для корректной записи файла, но полезен для вывода в консоль
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
        std::cout << u8"Нажмите Enter для выхода...";
        std::cin.get();
        return 1;
    }
    return 0;
}

