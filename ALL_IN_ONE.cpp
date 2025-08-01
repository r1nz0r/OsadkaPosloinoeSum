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

// Убедитесь, что json.hpp находится в папке проекта
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
    const SoilLayerData* soil;
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
    
    std::vector<CalculationLayerResult> SimplifyResultsForPlotting(
        const std::vector<CalculationLayerResult>& fullResults, double compressibleDepth) const;

    inline double GetSublayerThickness() const
    {
        return std::min(m_input.sublayerManualThickness_m, 0.4 * m_input.width_m);
    }

private:
    InputData m_input;
};

// --- Реализация методов ---

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

std::vector<CalculationLayerResult> SettlementCalculator::SimplifyResultsForPlotting(
    const std::vector<CalculationLayerResult>& fullResults, double compressibleDepth) const
{
    const size_t MAX_POINTS = 250;
    if (fullResults.size() <= MAX_POINTS) {
        return fullResults;
    }

    std::vector<CalculationLayerResult> simplified;
    std::vector<bool> keep(fullResults.size(), false);

    keep.front() = true;
    keep.back() = true;

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

    const size_t step = fullResults.size() / MAX_POINTS;
    if (step > 1) {
        for (size_t i = 0; i < fullResults.size(); i += step) {
            keep[i] = true;
        }
    }

    for (size_t i = 0; i < fullResults.size(); ++i) {
        if (keep[i]) {
            simplified.push_back(fullResults[i]);
        }
    }
    return simplified;
}


void SettlementCalculator::WriteHtmlReport(const std::vector<CalculationLayerResult>& fullResults, double compressibleDepth) const
{
    auto plotResults = SimplifyResultsForPlotting(fullResults, compressibleDepth);

    std::ofstream file("calculation_report.html");
    if (!file.is_open()) { return; }

    file << (char)0xEF << (char)0xBB << (char)0xBF;
    file << reinterpret_cast<const char*>(u8R"(<!DOCTYPE html><html lang="ru"><head><meta charset="UTF-8"><title>Отчет по расчету осадки</title><style>body{font-family:sans-serif;margin:2em;background-color:#f9f9f9}.container{max-width:1200px;margin:auto;background:#fff;padding:1em 2em;box-shadow:0 2px 5px rgba(0,0,0,0.1)}h1,h2{color:#333;border-bottom:2px solid #eee;padding-bottom:5px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:2em}.chart-container{padding:1em;border:1px solid #ddd;border-radius:5px}@media (max-width:900px){.grid{grid-template-columns:1fr}}</style></head><body>)");
    
    // --- ВСТРАИВАЕМ CHART.JS ---
    file << R"chartjscode(<script>!function(t,e){"object"==typeof exports&&"undefined"!=typeof module?module.exports=e():"function"==typeof define&&define.amd?define(e):(t="undefined"!=typeof globalThis?globalThis:t||self).Chart=e()}(this,(function(){"use strict";const t={0:"undefined",1:"string",2:"number",3:"boolean",4:"symbol"},e=new Map;function i(t,i){if(i){const s=i.get(t);if(s)return s}const n=e.get(t);if(n)return n;const o=new Intl.NumberFormat(t,i);return e.set(t,o),o}function s(e,s,n){const o="string"==typeof s?function(t,e){const s=i(t);return s.format(e)}(e,n):function(t,e,s){const n=i(t,e);return n.format(s)}(e,s,n);return"object"==t[typeof o]?void 0:o}const n={formatters:{values:t=>Array.isArray(t)?t:[t],_filter:(e,i,s)=>"number"!=t[typeof e[0]]||isNaN(e[0])?[]:e.map((t,e)=>({value:t,label:s.labels[i[e]]}))}};function o(t,e){return t.replace(/\{(\w+)\}/g,((t,i)=>e[i]))}const a="rgb(",r="rgba(",l="hsl(",c="hsla(",h="#";function d(t){if(t instanceof CanvasGradient||t instanceof CanvasPattern)return!t.addColorStop;if("string"!=typeof t)return!1;if(t.startsWith(h))return!0;if(t.startsWith(a))return!0;if(t.startsWith(r))return!0;if(t.startsWith(l))return!0;if(t.startsWith(c))return!0}function u(t){if("string"==typeof t)return t.startsWith("http")}const f=new class{constructor(){this._request=null,this._charts=new Map,this._running=!1,this._lastDate=void 0}request(t){const e=this;e._running||(e._running=!0,e._request=requestAnimationFrame((i=>{e.draw(i),e._running=!1,t&&t()})))}draw(t){const e=this;e._lastDate=t,e._charts.forEach(((t,e)=>{t.draw()}))}add(t){const e=this;e._charts.has(t)||e._charts.set(t,t)}remove(t){this._charts.delete(t)}}}{},_attach:t=>{t.options.animation&&t.options.animation.enabled&&f.add(t)},_detach:t=>{f.remove(t)}};function g(t,e,i,s){return Math.max(0,Math.min(t,i))-Math.max(0,Math.min(t,e))+s}const p=["x","y","z"],m=p.length;function v(t,e){let i,s,n,o;return e?i=e.filter((e=>g(e.min,t.min,t.max,e.max)))[0]:t.min>0&&t.max<1&&(i={min:0,max:1}),i||(i={min:t.min,max:t.max}),s=i.max-i.min,n=g(t.min,i.min,i.max,t.max),o=n/s,Math.pow(o,e?e.length-1:m-1)}function b(t,e,i){let s,n,o,a;const r=e.length,l=Math.floor(r/2);if(i)for(s=0,n=r-1;s<l;++s,--n)a=e[s],e[s]=e[n],e[n]=a;for(s=0,n=r;s<n;++s)o=e[s],t.has(o)||(t.add(o),i&&e.push(o))}function x(t,e){const i=t.indexOf(e);-1===i?t.push(e):t.splice(i,1)}const _={},y=["inherit","initial","revert","revert-layer","unset"];function w(t,e,i){const s=y.indexOf(t);return-1===s?Math.max(0,Math.min(e,i)):s}function M(t,e,i){for(const s of e)if(s===i)return!0;return!1}function k(t){return t.charAt(0).toUpperCase()+t.slice(1)}const S=t=>void 0!==t,P=t=>"function"==typeof t,D=Object.create(null),C=Object.create(null);function O(t,e){if(!e)return t;const i=e.split(".");for(let e=0,s=i.length;e<s;++e){const s=i[e];t=t[s]||(t[s]=Object.create(null))}return t}function A(t,e,i){return"string"==typeof e?O(t,e)[i]:e[i]}function T(t,e){const i=D[e]||(D[e]=function(t){const e=Object.create(null),i=Object.getOwnPropertyNames(t),s=Object.getOwnPropertySymbols(t);for(const n of i)e[n]=t[n];for(const o of s)e[o]=t[o];return e}(t));return i[e]}function L(t){return Array.isArray(t)?t:Object.keys(t).map((e=>({key:e,value:t[e]})))}function R(t,e){const i=Object.getOwnPropertyDescriptor(t,e);return i&&i.value}function E(t,e){t.dispatchEvent(new Event(e))}const I=["left","top","right","bottom"];function z(t,e){const{box:i,padding:s}=e,n=Math.max(0,i.width/2),o=Math.max(0,i.height/2),a=s.left+n,r=s.right+n,l=s.top+o,c=s.bottom+o;t.left+=a,t.right-=r,t.top+=l,t.bottom-=c,t.width=t.right-t.left,t.height=t.bottom-t.top}const F={left:t=>t.left,top:t=>t.top,right:t=>t.right,bottom:t=>t.bottom,width:t=>t.width,height:t=>t.height}[I[0]];function V(t,e,i){return Math.max(e,Math.min(i,t))}function B(t,e,i){let s;return s=e===i?t:V(t,e,i),0===s?0:s}const N=["start","end","center","inner"];function W(t,e){return"inner"===e?t:B(t,0,1)}function H(t,e,i){return B(t,e,i)}function j(t,e,i){return W(t,i)}function U(t,e){const i=["a","b","c","d","e","f"];for(let s=0;s<6;++s)if(t[i[s]]!==e[i[s]])return!1;return!0}const q=Math.PI,Y=2*q,X=Y/360,G=q/2,Z=q/4,K=q/180;function $(t){if("number"==typeof t)return t*K;if("string"==typeof t){const e=t.match(/(\d+)(deg|rad|grad|turn)/);if(e){const[,t,i]=e;switch(i){case"deg":return t*K;case"rad":return t;case"grad":return t*q/200;case"turn":return t*Y}}}return NaN}function J(t){return(t%Y+Y)%Y}function Q(t,e,i){return Math.max(e,Math.min(i,t))}function tt(t,e,i){let s;const n=J(t);return e?s=n/(e*Y):i&&(s=n*i/Y),s}function et(t,e){return Math.round(t*e)/e}function it(t,e,i,s){const n=J(t);let o,a;const r=J(e),l=J(i);return s?(o=et(n,s),a=et(r,s)):(o=n,a=r),o===a?l:J(l)}function st(t,e,i){const s=t.slice();return s.splice(e,i),s}function nt(t,e,i,s){return i?function(t,e){}(t,e):e}const ot=t=>t.map((t=>t.x)),at=t=>t.map((t=>t.y));function rt(t,e){let i,s,n,o,a,r,l,c;const h=t.length,d=[];for(i=0;i<h;i++)d.push([t[i].x,t[i].y]);for(i=0;i<h-1;i++){for(s=[d[i]],n=0,o=i+1;o<h-1&&!(n>e);o++)n+=function(t,e){return Math.sqrt(Math.pow(e.x-t.x,2)+Math.pow(e.y-t.y,2))}(t[o],t[i]);o<h-1&&(s.push([t[o].x,t[o].y]));for(a=0,r=0;r<s.length-1;r++)l=s[r],c=s[r+1],a+=l[0]*c[1]-c[0]*l[1];a/=2,d[i]=[a,t[i].x,t[i].y]}return d}function lt(t,e){let i,s,n,o,a,r,l,c;const h=t.length,d=[];for(i=0;i<h;i++)d.push([t[i].x,t[i].y]);for(i=0;i<h;i++)for(s=[d[i]],n=0,o=i;o<h;o++)n+=function(t,e){return Math.sqrt(Math.pow(e.x-t.x,2)+Math.pow(e.y-t.y,2))}(t[o],t[i]);for(a=0,r=0;r<s.length-1;r++)l=s[r],c=s[r+1],a+=l[0]*c[1]-c[0]*l[1];a/=2,d[i]=[a,t[i].x,t[i].y];return d}function ct(t,e,i){let s;const n=t.length;if(i)s=lt(t,i);else if("auto"!==e){const e=function(t){const e=ot(t),i=at(t),s=Math.min(...e),n=Math.min(...i),o=Math.max(...e),a=Math.max(...i);return{x:{min:s,max:o,delta:o-s},y:{min:n,max:a,delta:a-n}}}(t);s=rt(t,Math.min(e.x.delta,e.y.delta)/1e3)}else{const e=function(t){const e=t.length;let i,s;const n=[];for(i=0;i<e;i++)n.push([t[i].x,t[i].y]);for(i=0;i<e-1;i++)s=n[i],s.push(n[i+1][0],n[i+1][1]);return n}(t);s=rt(t,0)}return s}function ht(t,e,i){const s=t[e],n=t[i];return{x:n.x-s.x,y:n.y-s.y}}function dt(t,e,i){const s=ht(t,e,i);return Math.sqrt(Math.pow(s.x,2)+Math.pow(s.y,2))}function ut(t,e,i){const s=dt(t,e,i);return s>0?{x:t[i].x-t[e].x,y:t[i].y-t[e].y}:{x:0,y:0}}const ft=new Map;function gt(t,e,i){const s=function(t,e){e=e||{};const i=t+JSON.stringify(e);let s=ft.get(i);return s||(s=new Intl.NumberFormat(t,e),ft.set(i,s)),s}(e,i);return s.format(t)}const pt={millisecond:864e5,second:36e5,minute:6e4,hour:36e5,day:864e5,week:6048e5,month:2628e6,quarter:7884e6,year:31536e6};function mt(t,e){const i=t-e,s=new Date(i),n=s.getUTCFullYear()-1970,o=s.getUTCMonth(),a=s.getUTCDate()-1;let r,l,c;return n>0?(r=n+"y",l=o/12,c=a/365,l>.1&&l<.9?r+=(Math.round(10*l)/10)+"m":c>.1&&c<.9&&(r+=(Math.round(10*c)/10)+"d")):o>0?(r=o+"m",l=a/30,l>.1&&l<.9&&(r+=(Math.round(10*l)/10)+"d")):r=a+"d",r}const vt=["x","y","z"],bt=vt.length;function xt(t){const e=[];let i,s,n;for(i=0,s=t.length;i<s;i++)n=t[i],isNaN(n)||(e.push(n));return e}function _t(t,e){const i=xt(t);return e?Math.max(0,i.length-1):i.length}function yt(t,e,i){const s=t.length;if(!s)return{min:0,max:0,major:[],minor:[]};let n=0;const o=[],a=[];for(let t=0;t<s;t++)e>0&&t%e==0||(o.push(t),a.push(t));const r={min:0,max:s-1,major:o,minor:a};return i&&function(t,e){const i=t.major,s=t.minor,n=i.length,o=s.length;let a,r,l,c,h,d;if(e){for(a=0,r=n-1;a<r;++a,--r)d=i[a],i[a]=i[r],i[r]=d;for(a=0,r=o-1;a<r;++a,--r)d=s[a],s[a]=s[r],s[r]=d}for(a=0,l=n;a<l;++a)h=i[a],t.has(h)||(t.add(h),e&&i.push(h));for(a=0,c=o;a<c;++a)h=s[a],t.has(h)||(t.add(h),e&&s.push(h))}(r,i),r}function wt(t){const e=t.length;if(!e)return{min:void 0,max:void 0,major:[],minor:[]};let i,s,n;const o=[],a=[],r=new Set;for(i=0;i<e;i++){if(s=t[i],isNaN(s))continue;const e=Math.round(s);r.has(e)?a.push(i):(o.push(i),r.add(e)),n=s}return{min:0,max:n,major:o,minor:a}}function Mt(t,e){const i=t.length;if(!i)return;let s,n;for(s=0,n=i;s<n;++s)e(t[s],s)}const kt=["left","top","right","bottom"];function St(t,e){let i,s,n;const o=t.length;if(e)for(i=0,s=o/2;i<s;++i)n=t[i],t[i]=t[o-1-i],t[o-1-i]=n;for(i=0;i<o;++i)e?t[i].from=kt[i%4]:t[i].from=kt[i%2]}function Pt(t,e){let i,s,n;const o=t.length;for(i=0;i<o;++i)s=t[i],n=s.from,e?s[n]=s.to:s[n]=s.from}const Dt={id:"numbers",determineDataLimits:(t,e)=>{const{min:i,max:s}=function(t,e){const{min:i,max:s,options:n}=e,{stacked:o,min:a,max:r}=n,l=t.getMatchingVisibleMetas();let c,h,d;for(c=0,h=l.length;c<h;++c){d=l[c];const e=d.controller.getMinMax(d,o);void 0===i?i=e.min:i=Math.min(i,e.min),void 0===s?s=e.max:s=Math.max(s,e.max)}return{min:a||i,max:r||s}}(t,e);e.min=i,e.max=s},buildTicks:t=>{const e=t.options.ticks,i=function(t,e){const i=[],{bounds:s,step:n,min:o,max:a,precision:r,count:l,major:c,format:h,minUnit:d,minStep:u}=t,f=n||1,g=l||0,m=u||1,p=d?pt[d]:0;let v=o,b=a;("ticks"===s||"labels"===s)&&(v=Math.min(o,e.min),b=Math.max(a,e.max));const x=Math.round(v),_=Math.round(b);let y=x,w=_;c&&c.enabled&&(y=Math.round(e.min),w=Math.round(e.max));const M=Math.round((w-y)/f)+1,k=Math.round(e.max-e.min);let S=g||M;const P=Math.min(f,Math.max(1,Math.floor(k/S))),D=Math.round(k/P);let C=0;if(p>0){const t=(b-v)/p;C=Math.max(1,Math.floor(t))}const O=Math.max(m,C);let A=Math.ceil(v/O)*O;A=Math.min(b,Math.max(v,A));let T=0;for(;T<S;T++){const t=A+T*P;if(t>b)break;i.push({value:t,label:h?h(t):gt(t,e.locale),major:c&&c.enabled&&T%c.step===0})}return i}(t,t.chart.getDataLimits(t.id));return i.map((t=>Object.assign({},t,{major:e.major.enabled&&t.major})))}};class Ct{constructor(t){this.id=t.id,this.type=t.type,this.options=void 0,this.chart=void 0,this.ctx=void 0,this.data=void 0,this.label=void 0,this.index=void 0,this.visible=void 0,this.active=void 0,this._padding={left:0,top:0,right:0,bottom:0},this.width=void 0,this.height=void 0,this.top=void 0,this.bottom=void 0,this.left=void 0,this.right=void 0,this.box=void 0,this.margins=void 0,this.draw=void 0}getPadding(){return this._padding}update(t,e,i){const s=this;s.options=t,s.chart=e,s.index=i,s.ctx=e.ctx,s.data=e.data[i],s.label=e.data.labels[i],s.visible=e.isDatasetVisible(i),s.active=!1,s.width=e.chartArea.width,s.height=e.chartArea.height,s.top=e.chartArea.top,s.bottom=e.chartArea.bottom,s.left=e.chartArea.left,s.right=e.chartArea.right;const n=t.padding;s._padding={left:"number"==typeof n?n:n.left,top:"number"==typeof n?n:n.top,right:"number"==typeof n?n:n.right,bottom:"number"==typeof n?n:n.bottom},s.box={left:s.left+s._padding.l