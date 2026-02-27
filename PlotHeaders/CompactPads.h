#ifndef COMPACT_PADS_H
#define COMPACT_PADS_H
#include "include.h"
struct CompactPad
{
public:
    //========================================================
    // This override only aligns the pads.
    // Usage :
    // Suppose you want the pads to be aligned in 2 Rows and 5 Columns, with all the margins to be 0.10 :
    // std::vector<std::vector<TPad *>> pad; // std::vector is essential, the first index is row and the second one is column.
    // CompactPads(pad, 2, 5, 0.1, 0.1, 0.1, 0.1);
    // Suggestions :
    // xLatex->SetTextFont(val); val = 10 * fontnumber + 3, in which case the text size is given in pixels.
    //========================================================
    std::vector<std::vector<TPad *>> getCompactPads(Int_t Rows, Int_t Columns, Double_t LeftMargin, Double_t RightMargin, Double_t TopMargin, Double_t BottomMargin)
    {
        Refresh();
        std::vector<std::vector<TPad *>> pad(Rows);
        Double_t vStep = (1. - BottomMargin - TopMargin - (Rows - 1) * vSpacing) / Rows;
        Double_t hStep = (1. - LeftMargin - RightMargin - (Columns - 1) * hSpacing) / Columns;
        Double_t vposd, vposu, vmard, vmaru, vfactor;
        Double_t hposl, hposr, hmarl, hmarr, hfactor;
        for (Int_t irow = 0; irow < Rows; irow++)
        {
            if (Rows != 1)
            {
                if (irow == 0)
                {
                    vposd = 0.0;
                    vposu = BottomMargin + vStep;
                    vfactor = vposu - vposd;
                    vmard = BottomMargin / vfactor;
                    vmaru = 0.0;
                }
                else if (irow == Rows - 1)
                {
                    vposd = vposu + vSpacing;
                    vposu = vposd + vStep + TopMargin;
                    vfactor = vposu - vposd;
                    vmard = 0.0;
                    vmaru = TopMargin / (vposu - vposd);
                }
                else
                {
                    vposd = vposu + vSpacing;
                    vposu = vposd + vStep;
                    vfactor = vposu - vposd;
                    vmard = 0.0;
                    vmaru = 0.0;
                }
            }
            else
            {
                vposd = 0.0;
                vposu = BottomMargin + vStep;
                vfactor = vposu - vposd;
                vmard = BottomMargin / vfactor;
                vmaru = TopMargin / (vposu);
            }
            for (Int_t jcolumn = 0; jcolumn < Columns; jcolumn++)
            {
                if (Columns != 1)
                {
                    if (jcolumn == 0)
                    {
                        hposl = 0.0;
                        hposr = LeftMargin + hStep;
                        hfactor = hposr - hposl;
                        hmarl = LeftMargin / hfactor;
                        hmarr = 0.0;
                    }
                    else if (jcolumn == Columns - 1)
                    {
                        hposl = hposr + vSpacing;
                        hposr = hposl + hStep + RightMargin;
                        hfactor = hposr - hposl;
                        hmarl = 0.0;
                        hmarr = RightMargin / (hposr - hposl);
                    }
                    else
                    {
                        hposl = hposr + vSpacing;
                        hposr = hposl + hStep;
                        hfactor = hposr - hposl;
                        hmarl = 0.0;
                        hmarr = 0.0;
                    }
                }
                else
                {
                    hposl = 0.0;
                    hposr = LeftMargin + hStep;
                    hfactor = hposr - hposl;
                    hmarl = LeftMargin;
                    hmarr = RightMargin;
                }
                pad[irow].emplace_back(new TPad(Form("pad_%d_%d", irow, jcolumn), "", hposl, vposd, hposr, vposu));
                Configure(pad[irow][jcolumn], hmarl, hmarr, vmard, vmaru);
            }
        }
        return pad;
    }
    //========================================================
    // This override is prepared for 1-dimensional scenario
    // Usage :
    // Suppose you want the pads to be aligned in 5 Rows, with all the margins to be 0.10 :
    // std::vector<TPad *> pad;
    // CompactPads(pad, 5, 0.1, 0.1, 0.1, 0.1, true);
    //========================================================
    std::vector<TPad *> getCompactPads(Int_t nPads, Double_t LeftMargin, Double_t RightMargin, Double_t TopMargin, Double_t BottomMargin, Bool_t IsHorizontal)
    {
        Int_t Rows = IsHorizontal ? 1 : nPads;
        Int_t Columns = IsHorizontal ? nPads : 1;
        std::vector<std::vector<TPad *>> Pad2D = getCompactPads(Rows, Columns, LeftMargin, RightMargin, TopMargin, BottomMargin);
        std::vector<TPad *> pad;
        if (IsHorizontal)
        {
            pad = Pad2D[0];
        }
        else
        {
            pad.resize(nPads);
            for (Int_t i = 0; i < nPads; i++)
            {
                pad[i] = Pad2D[i][0];
            }
        }
        return pad;
    }
    std::vector<TPad *> getPropotionPads(const std::vector<double> &proportions, Bool_t IsHorizontal)
    {
        double sum = std::accumulate(proportions.begin(), proportions.end(), 0);
        std::vector<TPad *> pads(proportions.size());
        double smaller{1.0}, larger{1.0};
        for (size_t i = 0; i < proportions.size(); ++i)
        {
            larger = smaller;
            smaller -= proportions[i] / sum;
            if (IsHorizontal)
            {
                pads[i] = new TPad(Form("pad_%zu", i), "", smaller, 0, larger, 1);
            }
            else
            {
                pads[i] = new TPad(Form("pad_%zu", i), "", 0, smaller, 1, larger);
                printf("pad_%zu: smaller = %f, larger = %f\n", i, smaller, larger);
            }
            Configure(pads[i], 0., 0., 0., 0.);
        }
        return pads;
    }
    //========================================================
    // This method calculated the corrected position of the object in the pad, since margins of different pads are different.
    // Usage :
    // xLatex->DrawLatexNDC(0.75, GetPosition(0.6, false), hello);
    // xLatex->DrawLatexNDC(GetPositionX(0.3, true), 0.75, world);
    //========================================================
    Double_t GetPosition(Double_t Y, Bool_t IsVertical)
    {
        if (IsVertical)
        {
            return gPad->GetBottomMargin() + Y * (1. - gPad->GetBottomMargin() - gPad->GetTopMargin());
        }
        else
        {
            return gPad->GetLeftMargin() + Y * (1. - gPad->GetLeftMargin() - gPad->GetRightMargin());
        }
    }
    //========================================================
    // This method calculated the corrected length of ticks for pads with different sizes.
    // Usage :
    // hFrame->GetXaxis()->SetTickLength(GetTicks(0.04, false));
    // hFrame->GetYaxis()->SetTickLength(GetTicks(0.04, true));
    //========================================================
    Double_t GetTicks(Double_t Length, Bool_t IsVertical)
    {
        double padW = gPad->GetWw() * gPad->GetAbsWNDC();
        double padH = gPad->GetWh() * gPad->GetAbsHNDC();
        double mL = gPad->GetLeftMargin() * padW;
        double mR = gPad->GetRightMargin() * padW;
        double mB = gPad->GetBottomMargin() * padH;
        double mT = gPad->GetTopMargin() * padH;
        double tickScale = IsVertical ? (padH - mB - mT)/padH * padW : (padW - mL - mR)/padW * padH;
        return Length / tickScale;
    }
    std::pair<TPad *, TPad *> SetTwoPads(Bool_t IsVertical, Double_t Fraction)
    {
        TPad *fUpper, *fLower;
        if (IsVertical)
        {
            fUpper = new TPad("fUpper", "", 0, Fraction, 1, 1);
            fLower = new TPad("fLower", "", 0, 0, 1, Fraction);
        }
        else
        {
            fUpper = new TPad("fUpper", "", 0, 0, Fraction, 1);
            fLower = new TPad("fLower", "", Fraction, 0, 1, 1);
        }
        fUpper->SetMargin(0., 0., 0., 0.);
        fLower->SetMargin(0., 0., 0., 0.);
        return std::make_pair(fUpper, fLower);
    }
    std::pair<TPad *, TPad *> SetTwoPads(Bool_t IsVertical, Double_t Fraction, Double_t left, Double_t right, Double_t bottom, Double_t top)
    {
        std::pair<TPad *, TPad *> PadPair = SetTwoPads(IsVertical, Fraction);
        if (IsVertical)
        {
            PadPair.first->SetMargin(left, right, 0., top);
            PadPair.second->SetMargin(left, right, bottom, 0);
        }
        else
        {
            PadPair.first->SetMargin(left, 0., bottom, top);
            PadPair.second->SetMargin(0., right, bottom, top);
        }
        PadPair.first->cd();
        return PadPair;
    }
    static CompactPad &getInstance()
    {
        static CompactPad instance;
        return instance;
    }
    void SetVSpacing(Double_t vSpacing) { this->vSpacing = vSpacing; }
    void SetHSpacing(Double_t hSpacing) { this->hSpacing = hSpacing; }
    Double_t GetVSpacing() const { return vSpacing; }
    Double_t GetHSpacing() const { return hSpacing; }
    void Refresh()
    {
        Count = 0;
    }
    void SetSpacing(Double_t vSpacing, Double_t hSpacing)
    {
        this->vSpacing = vSpacing;
        this->hSpacing = hSpacing;
    }
    std::vector<TPad *> ClonePads(const std::vector<TPad *> &pad)
    {
        std::vector<TPad *> clone;
        for (auto &i : pad)
        {
            clone.emplace_back((TPad *)i->Clone());
        }
        return clone;
    }
    void Configure(TVirtualPad *pad, Double_t left, Double_t right, Double_t bottom, Double_t top)
    {
        pad->SetMargin(left, right, bottom, top);
        pad->SetFrameBorderMode(0);
        pad->SetTicks(1, 1);
        pad->SetBorderMode(0);
        pad->SetBorderSize(0);
        pad->SetFillColor(kWhite);
        return;
    }

private:
    CompactPad() {};
    Int_t Count{0};
    // Double_t AbsWidth{1.0};
    // Double_t AbsHeight{1.0};
    Double_t vSpacing{0.0};
    Double_t hSpacing{0.0};
};

std::vector<TVirtualPad *> getSparsePads(Int_t nColumn, Double_t Shift, Double_t left, Double_t right, Double_t bottom, Double_t top)
{
    std::vector<TVirtualPad *> pads(nColumn);
    Double_t Fraction = (1 + Shift * nColumn - Shift) / nColumn;
    std::pair<TPad *, TPad *> PadPair = CompactPad::getInstance().SetTwoPads(false, Fraction);
    pads[0] = PadPair.first;
    PadPair.first->Draw();
    PadPair.second->Draw();
    PadPair.second->Divide(nColumn - 1, 1, 0., 0);
    for (Int_t iColumn = 1; iColumn < nColumn; iColumn++)
    {
        pads[iColumn] = PadPair.second->cd(iColumn);
        CompactPad::getInstance().Configure(pads[iColumn], left, right, bottom, top);
    }
    CompactPad::getInstance().Configure(PadPair.first, 1 - (1 - Fraction) * (1 - left) / (Fraction * (nColumn - 1)), right * (1 - Fraction) / (nColumn - 1) / Fraction, bottom, top);
    return pads;
}
#endif