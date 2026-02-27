#ifndef TGRAPH_EXTENSIONS_HH
#define TGRAPH_EXTENSIONS_HH
#include "TGraphErrors.h"
void AdjustNpart(TGraph *gData, TGraph *gModel)
{
    if (gData->GetN() != gModel->GetN())
    {
        printf("Graphs have different number of points\n");
        return;
    }
    for (Int_t i = 0; i < gData->GetN(); i++)
    {
        Double_t x, y;
        gData->GetPoint(i, x, y);
        gModel->SetPointX(i, x);
    }
}
std::pair<Double_t, Double_t> FindRange(const TGraph *g, Bool_t IsXaxis = false)
{
    Double_t min = 1e9;
    Double_t max = -1e4;
    for (size_t ip = 0; ip < g->GetN(); ip++)
    {
        Double_t val = IsXaxis ? g->GetX()[ip] : g->GetY()[ip];
        Double_t err = IsXaxis ? g->GetErrorX(ip) : g->GetErrorY(ip);
        min = val - err < min ? val - err : min;
        max = val + err > max ? val + err : max;
    }
    return std::make_pair(min, max);
}
std::pair<Double_t, Double_t> FindRange(std::vector<TGraph *> &graphVec, Bool_t IsXaxis = false)
{
    Double_t min = 1e8;
    Double_t max = -1e8;
    Int_t ig = 0;
    for (auto g : graphVec)
    {
        if (!g)
        {
            std::cerr << "Graph " << ig << " is null" << std::endl;
            continue;
        }
        std::pair<Double_t, Double_t> range = FindRange(g, IsXaxis);
        min = range.first < min ? range.first : min;
        max = range.second > max ? range.second : max;
        ig++;
    }
    return std::make_pair(min, max);
}
std::vector<double> GetRangeFrame(const std::vector<TGraph *> &graph)
{
    std::vector<double> Range{1e9, 1e9, 0, 0};
    for (auto g : graph)
    {
        std::pair<Double_t, Double_t> range = FindRange(g, true);
        Range[0] = range.first < Range[0] ? range.first : Range[0];
        Range[2] = range.second > Range[2] ? range.second : Range[2];
        range = FindRange(g, false);
        Range[1] = range.first < Range[1] ? range.first : Range[1];
        Range[3] = range.second > Range[3] ? range.second : Range[3];
    }

    return Range;
}

TGraphErrors *CombineError(TGraphErrors *g1, TGraphErrors *g2)
{
    if (g1->GetN() != g2->GetN())
    {
        Info("CombineError", "Graphs have different number of points");
        return nullptr;
    }
    TGraphErrors *g = new TGraphErrors();
    for (size_t ip = 0; ip < g1->GetN(); ip++)
    {
        Double_t x1, y1, x2, y2;
        g1->GetPoint(ip, x1, y1);
        g2->GetPoint(ip, x2, y2);
        g->SetPoint(ip, x1, y1);
        g->SetPointError(ip, 0, sqrt(pow(g1->GetErrorY(ip), 2) + pow(g2->GetErrorY(ip), 2)));
    }
    return g;
}

TGraphErrors *GetRatio(TGraphErrors *dividend, TGraphErrors *divisor)
{
    TGraphErrors *Ratio = new TGraphErrors();
    for (int i = 0; i < dividend->GetN(); i++)
    {
        Double_t x, y, x1, y1;
        dividend->GetPoint(i, x, y);
        divisor->GetPoint(i, x1, y1);
        Ratio->SetPoint(i, x, y / y1);
        Ratio->SetPointError(i, 0, sqrt(pow(dividend->GetErrorY(i) / y1, 2) + pow(divisor->GetErrorY(i) * y / pow(y1, 2), 2)));
    }
    return Ratio;
}
TGraphErrors *GetDifference(TGraphErrors *minuend, TGraphErrors *subtrahend)
{
    TGraphErrors *gDifference = new TGraphErrors();
    for (int i = 0; i < minuend->GetN(); i++)
    {
        Double_t x, y, x1, y1;
        minuend->GetPoint(i, x, y);
        subtrahend->GetPoint(i, x1, y1);
        gDifference->SetPoint(i, x, y - y1);
        gDifference->SetPointError(i, 0, sqrt(pow(minuend->GetErrorY(i), 2) + pow(subtrahend->GetErrorY(i), 2)));
    }
    return gDifference;
}
TGraphErrors *GetDeviation(TGraphErrors *gData, TGraphErrors *gBenchmark)
{
    TGraphErrors *gRatio = new TGraphErrors();
    for (Int_t ip = 0; ip < gData->GetN(); ip++)
    {
        Double_t x, y, x1, y1;
        gData->GetPoint(ip, x, y);
        double min{1e5};
        int rightPoint{-1};
        for (Int_t jp = 0; jp < gBenchmark->GetN(); jp++)
        {
            gBenchmark->GetPoint(jp, x1, y1);
            if (fabs(x - x1) < min)
            {
                min = fabs(x - x1);
                rightPoint = jp;
            }
        }
        Double_t yMatched = gBenchmark->GetY()[rightPoint];
        Double_t sigma = sqrt(pow(gData->GetErrorY(ip), 2) + pow(gBenchmark->GetErrorY(rightPoint), 2));
        gRatio->SetPoint(gRatio->GetN(), x, (y - yMatched) / sigma);
        gRatio->SetPointError(gRatio->GetN() - 1, 0, sigma);
        printf("Data: (%.1f, %f), Benchmark: %f, sigma: %f, Deviation: %f\n",x, y, yMatched, sigma, (y - yMatched) / sigma);
    }
    return gRatio;
}

TGraph *GetRatio(TGraph *graph, TGraph *oldgraph)
{
    TGraphErrors *Ratio = new TGraphErrors();
    for (int i = 0; i < graph->GetN(); i++)
    {
        Double_t x, y, x1, y1;
        graph->GetPoint(i, x, y);
        double err = graph->GetEY()[i];
        oldgraph->GetPoint(i, x1, y1);
        double err1 = oldgraph->GetEY()[i];
        Ratio->SetPoint(i, x, y / y1);
        Ratio->SetPointError(i, 0, sqrt(pow(err / y1, 2) + pow(err1 * y / pow(y1, 2), 2)));
    }
    return Ratio;
}

TGraph *GetRatio(TGraphAsymmErrors *graph, TGraphAsymmErrors *oldgraph)
{
    TGraphErrors *Ratio = new TGraphErrors();
    Int_t ipoint = 0;
    for (int i = 0; i < graph->GetN(); i++)
    {
        Double_t x, y, x1, y1;
        graph->GetPoint(i, x, y);
        double err = graph->GetErrorY(i);
        oldgraph->GetPoint(i, x1, y1);
        double err1 = oldgraph->GetErrorY(i);
        Ratio->SetPoint(ipoint, x, y / y1);
        Ratio->SetPointError(ipoint, 0, sqrt(pow(err / y1, 2) + pow(err1 * y / pow(y1, 2), 2)));
        ipoint++;
    }
    return Ratio;
}

Double_t GetMaximum(TGraphErrors *graph, Bool_t IsX = 0)
{
    Double_t *Values = IsX ? graph->GetX() : graph->GetY();
    Double_t *Errors = IsX ? graph->GetEX() : graph->GetEY();
    Double_t max = Values[0] + Errors[0];
    for (int i = 1; i < graph->GetN(); i++)
    {
        Double_t val = Values[i] + Errors[i];
        max = val > max ? val : max;
    }
    return max;
}

Double_t GetMaximum(TGraph *graph, Bool_t IsX = 0)
{
    Double_t *Values = IsX ? graph->GetX() : graph->GetY();
    Double_t max = Values[0];
    for (int i = 1; i < graph->GetN(); i++)
    {
        if (Values[i] > max)
        {
            max = Values[i];
        }
    }
    return max;
}

Double_t GetMinimum(TGraphErrors *graph, Bool_t IsX = 0)
{
    Double_t *Values = IsX ? graph->GetX() : graph->GetY();
    Double_t *Errors = IsX ? graph->GetEX() : graph->GetEY();
    Double_t min = Values[0] - Errors[0];
    for (int i = 1; i < graph->GetN(); i++)
    {
        Double_t val = Values[i] - Errors[i];
        min = val < min ? val : min;
    }
    return min;
}

Double_t GetMinimum(TGraph *graph, Bool_t IsX = 0)
{
    Double_t *Values = IsX ? graph->GetX() : graph->GetY();
    Double_t min = Values[0];
    for (int i = 1; i < graph->GetN(); i++)
    {
        if (Values[i] < min)
        {
            min = Values[i];
        }
    }
    return min;
}

TGraphErrors *Convert2Errors(TGraph *graph)
{
    TGraphErrors *g = new TGraphErrors(graph->GetN());
    for (Int_t i = 0; i < graph->GetN(); i++)
    {
        Double_t x, y;
        graph->GetPoint(i, x, y);
        g->SetPoint(i, x, y);
        g->SetPointError(i, 0, 0);
    }
    return g;
}

TGraphErrors *Convert2Errors(TGraphAsymmErrors *graph)
{
    TGraphErrors *g = new TGraphErrors(graph->GetN());
    for (Int_t i = 0; i < graph->GetN(); i++)
    {
        Double_t x, y;
        graph->GetPoint(i, x, y);
        g->SetPoint(i, x, y);
        g->SetPointError(i, 0, graph->GetErrorY(i));
    }
    return g;
}

TGraphErrors *ShiftXaxis(TGraphErrors *graph, Double_t shift)
{
    TGraphErrors *g = new TGraphErrors(graph->GetN());
    g->SetMarkerStyle(graph->GetMarkerStyle());
    for (Int_t i = 0; i < graph->GetN(); i++)
    {
        Double_t x, y;
        graph->GetPoint(i, x, y);
        g->SetPoint(i, x + shift, y);
        g->SetPointError(i, 0, graph->GetErrorY(i));
    }
    return g;
}

TGraphErrors *ShiftXaxis(TGraphErrors *graph)
{
    TGraphErrors *g = new TGraphErrors(graph->GetN());
    g->SetMarkerStyle(graph->GetMarkerStyle());
    for (Int_t i = 0; i < graph->GetN(); i++)
    {
        Double_t x, y;
        graph->GetPoint(i, x, y);
        g->SetPoint(i, x * 1.1, y);
        g->SetPointError(i, 0, graph->GetErrorY(i));
    }
    return g;
}

Int_t SearchPoint(TGraphErrors *graph, Double_t x)
{
    Int_t point = 0;
    Double_t Minimum = 1e5;
    for (Int_t i = 0; i < graph->GetN(); i++)
    {
        Double_t difference = fabs(graph->GetX()[i] - x);
        if (difference < Minimum)
        {
            Minimum = difference;
            point = i;
        }
    }
    return point;
}

#endif