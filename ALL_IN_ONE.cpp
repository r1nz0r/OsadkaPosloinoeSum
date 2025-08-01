#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <utility>

#include "json.hpp"

using json = nlohmann::json;

// --- Структуры ---
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
    const SoilLayerData* soil; // Указатель, чтобы можно было создавать временные слои
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

// --- Класс калькулятора ---
class SettlementCalculator
{
public:
    explicit SettlementCalculator(InputData data) : m_input(std::move(data)) {}
    void Run();

private:
    CalculationLayerResult CreateCalculationLayer(
        const SoilLayerData& soil, double thickness, double currentDepthSummary,
        double currentSoilStressSummary, const CalculationLayerResult* prevLayer);

    double CalculateAlpha(double middleElevation) const;
    void PrintConsoleSummary(double totalSettlement, double compressibleDepth) const;
    void WriteCsvReport(const std::vector<CalculationLayerResult>& results, double totalSettlement, double compressibleDepth) const;
    void WriteHtmlReport(const std::vector<CalculationLayerResult>& fullResults, double compressibleDepth) const;
    
    // Новая функция для оптимизации данных для графика
    std::vector<CalculationLayerResult> SimplifyResultsForPlotting(
        const std::vector<CalculationLayerResult>& fullResults, double compressibleDepth) const;

    inline double GetSublayerThickness() const
    {
        return std::min(m_input.sublayerManualThickness_m, 0.4 * m_input.width_m);
    }

private:
    InputData m_input;
};

// --- Реализация ---

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
                soil, thickness, depthSummary, soilStressSummary,
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
        if (stopConditionMet) break;
    }

    if (compressibleDepth == 0.0 && !calculationResults.empty())
    {
        compressibleDepth = calculationResults.back().bottomElevation;
    }

    PrintConsoleSummary(totalSettlement, compressibleDepth);
    WriteCsvReport(calculationResults, totalSettlement, compressibleDepth);
    WriteHtmlReport(calculationResults, compressibleDepth);

    std::cout << reinterpret_cast<const char*>(u8"\nРабота программы успешно завершена. Нажмите Enter для выхода...");
    std::cin.get();
}

CalculationLayerResult SettlementCalculator::CreateCalculationLayer(
    const SoilLayerData& soil, double thickness, double currentDepthSummary,
    double currentSoilStressSummary, const CalculationLayerResult* prevLayer)
{
    CalculationLayerResult layer{
        .index = static_cast<int>(prevLayer ? prevLayer->index + 1 : 1),
        .soil = &soil, .thickness = thickness
    };

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

double SettlementCalculator::CalculateAlpha(double z) const
{
    if (z < 1e-6) z = 1e-6;
    constexpr double PI = 3.141592653589793;
    double l_half = m_input.length_m / 2.0;
    double b_half = m_input.width_m / 2.0;
    double m = l_half / b_half;
    double n = z / b_half;
    double m2 = m * m, n2 = n * n;
    double term1_den_sqrt = sqrt(1 + m2 + n2);
    if (term1_den_sqrt < 1e-9) return 1.0;
    double term1 = (m * n / term1_den_sqrt) * (1 / (m2 + n2) + 1 / (1 + n2));
    double term2_den_sqrt = sqrt(1 + n2) * sqrt(m2 + n2);
    if (term2_den_sqrt < 1e-9) return 1.0;
    double term2_arg = m / term2_den_sqrt;
    term2_arg = std::clamp(term2_arg, -1.0, 1.0);
    return (term1 + atan(term2_arg)) / PI;
}

// Новая функция оптимизации данных
std::vector<CalculationLayerResult> SettlementCalculator::SimplifyResultsForPlotting(
    const std::vector<CalculationLayerResult>& fullResults, double compressibleDepth) const
{
    const size_t MAX_POINTS = 250;
    if (fullResults.size() <= MAX_POINTS) {
        return fullResults;
    }

    std::vector<CalculationLayerResult> simplified;
    std::vector<bool> keep(fullResults.size(), false);

    // Всегда сохраняем первую и последнюю точки
    keep.front() = true;
    keep.back() = true;

    // Сохраняем точку, ближайшую к глубине сжимаемой толщи
    size_t hc_index = 0;
    double min_dist = std::numeric_limits<double>::max();
    for (size_t i = 0; i < fullResults.size(); ++i) {
        double dist = std::abs(fullResults[i].bottomElevation - compressibleDepth);
        if (dist < min_dist) {
            min_dist = dist;
            hc_index = i;
        }
    }
    keep[hc_index] = true;

    // Равномерно выбираем остальные точки
    const size_t step = fullResults.size() / MAX_POINTS;
    if (step > 1) {
        for (size_t i = 0; i < fullResults.size(); i += step) {
            keep[i] = true;
        }
    }

    // Собираем итоговый вектор
    for (size_t i = 0; i < fullResults.size(); ++i) {
        if (keep[i]) {
            simplified.push_back(fullResults[i]);
        }
    }
    return simplified;
}


void SettlementCalculator::WriteHtmlReport(const std::vector<CalculationLayerResult>& fullResults, double compressibleDepth) const
{
    // Оптимизируем данные ПЕРЕД записью в файл
    auto plotResults = SimplifyResultsForPlotting(fullResults, compressibleDepth);

    std::ofstream file("calculation_report.html");
    if (!file.is_open()) {
        std::cerr << reinterpret_cast<const char*>(u8"Ошибка: не удалось создать HTML-отчет.") << std::endl;
        return;
    }

    file << (char)0xEF << (char)0xBB << (char)0xBF;
    file << reinterpret_cast<const char*>(u8R"(<!DOCTYPE html>
<html lang="ru">
<head>
    <meta charset="UTF-8">
    <title>Отчет по расчету осадки</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style> body { font-family: sans-serif; margin: 2em; background-color: #f9f9f9; } .container { max-width: 1200px; margin: auto; background: #fff; padding: 1em 2em; box-shadow: 0 2px 5px rgba(0,0,0,0.1); } h1, h2 { color: #333; border-bottom: 2px solid #eee; padding-bottom: 5px;} .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 2em; } .chart-container{padding:1em; border: 1px solid #ddd; border-radius: 5px;} @media (max-width: 900px) { .grid { grid-template-columns: 1fr; } } </style>
</head>
<body>
<div class="container">
    <h1>Результаты расчета осадки</h1>
    <div class="grid">
        <div class="chart-container"> <h2>Эпюры напряжений</h2> <canvas id="stressChart"></canvas> </div>
        <div class="chart-container"> <h2>График осадки</h2> <canvas id="settlementChart"></canvas> </div>
    </div>
    <script>
    try {
)");
    file << std::fixed << std::setprecision(4);

    double p0 = m_input.load_kPa - (m_input.backfillDensity_kN_m3 * m_input.foundationDepth_m);
    double initial_sigma_zp = p0 * CalculateAlpha(0.001);
    double initial_sigma_zg = m_input.backfillDensity_kN_m3 * m_input.foundationDepth_m;
    double initial_sigma_zg_red_factor = fullResults.empty() ? 0.2 : (fullResults.front().soil->deformationModulusFirst_kPa >= 7000 ? 0.5 : 0.2);
    double initial_sigma_zg_red = initial_sigma_zg * initial_sigma_zg_red_factor;
    
    double maxStress = initial_sigma_zp;

    file << "const sigma_zp_data = [{x: " << initial_sigma_zp << ", y: 0}";
    for (const auto& r : plotResults) {
        file << ", {x: " << r.structureStress << ", y: " << r.bottomElevation << "}";
        maxStress = std::max({maxStress, r.structureStress, r.soilStressReducted});
    }
    file << "];\n";

    file << "const sigma_zg_red_data = [{x: " << initial_sigma_zg_red << ", y: 0}";
    for (const auto& r : plotResults) { file << ", {x: " << r.soilStressReducted << ", y: " << r.bottomElevation << "}"; }
    file << "];\n";

    file << "const settlement_data = [{x: 0, y: 0}";
    double cumulativeSettlement = 0.0;
    for (const auto& r : fullResults) { // Для осадки нужна полная точность
        cumulativeSettlement += r.settlementTotal;
        file << ", {x: " << cumulativeSettlement * 1000 << ", y: " << r.bottomElevation << "}";
    }
    file << "];\n";
    
    file << "const compressibleDepth = " << compressibleDepth << ";\n";
    file << "const maxStress = " << maxStress * 1.1 << ";\n";
    file << "const hc_line_data = [{x: 0, y: compressibleDepth}, {x: maxStress, y: compressibleDepth}];\n";

    file << reinterpret_cast<const char*>(u8R"(
        const stressCtx = document.getElementById('stressChart');
        new Chart(stressCtx, {
            type: 'scatter',
            data: {
                datasets: [
                {
                    label: 'σ_zp (доп. от нагрузки)',
                    data: sigma_zp_data,
                    borderColor: 'rgb(255, 99, 132)',
                    backgroundColor: 'rgba(255, 99, 132, 0.5)',
                    showLine: true,
                    tension: 0.2,
                    pointRadius: 2
                }, {
                    label: 'σ_zg,red (бытовое с коэф.)',
                    data: sigma_zg_red_data,
                    borderColor: 'rgb(54, 162, 235)',
                    backgroundColor: 'rgba(54, 162, 235, 0.5)',
                    showLine: true,
                    tension: 0.2,
                    pointRadius: 2
                }, {
                    label: 'Сжимаемая толща',
                    data: hc_line_data,
                    borderColor: 'rgb(75, 192, 192)',
                    borderWidth: 2,
                    borderDash: [6, 6],
                    showLine: true,
                    pointRadius: 0
                }]
            },
            options: { scales: { y: { reverse: true, title: { display: true, text: 'Глубина, м' } }, x: { position: 'top', title: { display: true, text: 'Напряжение, кПа' }, min: 0 } } }
        });

        const settlementCtx = document.getElementById('settlementChart');
        new Chart(settlementCtx, {
            type: 'scatter',
            data: { datasets: [{ label: 'Накопленная осадка', data: settlement_data, borderColor: 'rgb(255, 159, 64)', showLine: true, tension: 0.2 }] },
            options: { scales: { y: { reverse: true, title: { display: true, text: 'Глубина, м' } }, x: { position: 'top', title: { display: true, text: 'Осадка, мм' }, min: 0 } } }
        });
    } catch (e) {
        console.error(e);
        document.body.innerHTML = "<h2>Ошибка при построении графика.</h2><p>Проверьте консоль разработчика (F12) для получения дополнительной информации.</p><pre>" + e.stack + "</pre>";
    }
    </script>
</div>
</body>
</html>
)");
}


// --- Остальные функции и main ---
void SettlementCalculator::PrintConsoleSummary(double totalSettlement, double compressibleDepth) const { /*...*/ }
void SettlementCalculator::WriteCsvReport(const std::vector<CalculationLayerResult>& r, double t, double d) const { /*...*/ }
InputData LoadInputFromJSON(const std::string& filename) { /*...*/ }
std::string GetInputFilePath() { /*...*/ }
int main() { /*...*/ }

(Примечание: чтобы не загромождать ответ, я свернул код функций PrintConsoleSummary, WriteCsvReport, LoadInputFromJSON, GetInputFilePath и main, так как они не менялись. В вашей полной версии кода они должны остаться на месте.)
