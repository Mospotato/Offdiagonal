#ifndef CLONE_WITH_EDGE_H
#define CLONE_WITH_EDGE_H
#include "GraphTools.h"
std::map<Style_t, Style_t> EdgeMap{{20, 24}, {21, 25}, {22, 26}, {23, 32}, {29, 30}, {33, 27}, {34, 28}};
std::map<Style_t, Style_t> FillMap{{24, 20}, {25, 21}, {22, 26}, {23, 32}, {29, 30}, {33, 27}, {34, 28}};
void SetBoldEdge(Bool_t Bold = true)
{
    if (Bold)
    {
        EdgeMap[20] = 53;
        EdgeMap[33] = 56;
    }
    return;
}

Style_t GetEdgeStyle(Style_t style)
{
    if (EdgeMap.count(style))
    {
        // std::cout << "Converted: " << style << " to " << EdgeMap[style] << '\n';
        return EdgeMap[style];
    }
    return style;
}

Style_t GetFillStyle(Style_t style)
{
    for (auto pair : EdgeMap)
    {
        if (pair.second == style)
        {
            // std::cout << "Converted: " << style << " to " << pair.first << '\n';
            return pair.first;
        }
    }
    // Info("GetEdgeStyle", "Style %d does not have corresponding fill style.", style);
    return style;
}

void DrawEdge(TGraph *graph)
{
    auto *objClone = Convert2Errors(graph);
    objClone->SetMarkerColor(1);
    objClone->SetMarkerSize(graph->GetMarkerSize());
    objClone->SetMarkerStyle(GetEdgeStyle(graph->GetMarkerStyle()));
    objClone->SetLineColor(graph->GetLineColor());
    objClone->SetLineStyle(graph->GetLineStyle());
    objClone->SetLineWidth(graph->GetLineWidth());
    objClone->Draw(graph->GetDrawOption());
}

void DrawEdge(TGraph *graph, Option_t *option)
{
    auto *objClone = Convert2Errors(graph);
    objClone->SetMarkerColor(1);
    objClone->SetMarkerSize(graph->GetMarkerSize());
    objClone->SetMarkerStyle(GetEdgeStyle(graph->GetMarkerStyle()));
    objClone->SetLineColor(graph->GetLineColor());
    objClone->SetLineStyle(graph->GetLineStyle());
    objClone->SetLineWidth(graph->GetLineWidth());
    objClone->Draw(option);
}

void ClearError(TGraph *graph)
{
    if (graph->InheritsFrom(TGraphErrors::Class()))
    {
        TGraphErrors *gClone = dynamic_cast<TGraphErrors *>(graph);
        for (Int_t i = 0; i < gClone->GetN(); i++)
        {
            gClone->SetPointError(i, 0, 0);
        }
    }
    else if (graph->InheritsFrom(TGraphAsymmErrors::Class()))
    {
        TGraphAsymmErrors *gClone = dynamic_cast<TGraphAsymmErrors *>(graph);
        for (Int_t i = 0; i < gClone->GetN(); i++)
        {
            gClone->SetPointError(i, 0, 0, 0, 0);
        }
    }
}

void GetEdgeFromLegend(TLegend *leg, Bool_t DrawObject = true)
{
    TLegend *legClone = (TLegend *)leg->Clone();
    legClone->SetFillStyle(0);
    legClone->Clear();
    auto primitives = leg->GetListOfPrimitives();
    for (auto primitiveObj : *primitives)
    {
        auto primitive = (TLegendEntry *)primitiveObj;
        TObject *object = primitive->GetObject();
        if (!object)
        {
            legClone->AddEntry(object, primitive->GetLabel(), primitive->GetOption());
            continue;
        }
        Bool_t IsGraph = object->InheritsFrom(TGraph::Class());
        Bool_t IsHisto = object->InheritsFrom(TH1::Class());
        TString DrawOption = object->GetDrawOption();
        if (!IsGraph && !IsHisto)
        {
            legClone->AddEntry(object, primitive->GetLabel(), primitive->GetOption());
            continue;
        }
        if (DrawObject && !(IsGraph && DrawOption.Contains("p")) && !(IsHisto && DrawOption.Contains("HIST")))
        {
            legClone->AddEntry(object, primitive->GetLabel(), primitive->GetOption());
            continue;
        }
        auto *objClone = object->Clone();
        if (IsGraph)
        {
            TGraph *gObject = dynamic_cast<TGraph *>(object);
            TGraph *gClone = dynamic_cast<TGraph*>(objClone);
            ClearError(gClone);
            Style_t MarkerStyle = gObject->GetMarkerStyle();
            Style_t EdgeStyle = GetEdgeStyle(MarkerStyle);
            gClone->SetMarkerColor(1);
            gClone->SetMarkerStyle(EdgeStyle);
        }
        else
        {
            TH1 *hClone = dynamic_cast<TH1*>(objClone);
            hClone->SetFillStyle(1);
            printf("Adding %s with option %s\n", object->GetName(), object->GetDrawOption());
            if (!DrawOption.Contains("same"))
            {
                DrawOption += " same";
            }
        }
        if (DrawObject)
            objClone->Draw(DrawOption);
        legClone->AddEntry(objClone, primitive->GetLabel(), primitive->GetOption());
    }
    if (leg->GetHeader())
    {
        legClone->SetHeader(leg->GetHeader(), "C");
    }
    legClone->Draw(leg->GetDrawOption());
}
#endif