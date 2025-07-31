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
