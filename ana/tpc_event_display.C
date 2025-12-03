// macro for TPC event display
// original source is from: https://github.com/hyptpc/tpc_ana_example
// check usage in README here
// orginaly created by R.J.Saito

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <utility>
#include <vector>

#include "../include/padHelper.hh"
#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2Poly.h"
#include "TEllipse.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TParticle.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TVector3.h"

static TString g_rootfile_path;
static TCanvas* c_canvas_tpc = nullptr;

// Global geometry parameters
TString dcgeom_path = "../param/DCGEO/DCGeomParam_E90_20251202";
TString dsize_path = "../param/DSIZE/DetSize_E90_20251202";
Double_t gTargetZ_Global = 0.0;
Double_t gTpcPadCenterZ = 0.0;
Double_t gTargetRadius = 0.0;
Double_t gTargetHolderRadius = 0.0;

inline TVector3 LocalToWorld(Double_t x_local, Double_t y_local, Double_t z_local)
{
    // HypTPC placement: rotateX(90 deg) then rotateZ(180 deg).
    // Geant4 applies rotations in sequence, so the lab mapping becomes
    //   X_lab = -X_local,
    //   Y_lab =  Z_local,
    //   Z_lab = -Y_local.
    return TVector3(-x_local, z_local, -y_local);
}

inline std::pair<Double_t, Double_t> ProjectToDisplay(const TVector3& local)
{
    const auto lab = LocalToWorld(local.X(), local.Y(), local.Z());
    const Double_t z_disp = -(lab.Z() + gTpcPadCenterZ);
    const Double_t x_disp = lab.X();
    return {z_disp, x_disp};
}

void read_geometry()
{
    std::string line;

    std::ifstream dc_file(dcgeom_path);
    while (std::getline(dc_file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        int id;
        std::string name;
        double x, y, z;
        ss >> id >> name;
        if (name == "Target") {
            ss >> x >> y >> z;
            gTargetZ_Global = z;
            // std::cout << "Info: Read Target Global Z = " << gTargetZ_Global << " mm" << std::endl;
        }
    }

    std::ifstream size_file(dsize_path);
    while (std::getline(size_file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string name;
        ss >> name;
        if (name == "TpcPadCenterZ") {
            ss >> gTpcPadCenterZ;
            // std::cout << "Info: Read TpcPadCenterZ = " << gTpcPadCenterZ << " mm" << std::endl;
        } else if (name == "Target") {
            double trin, trout, tlength;
            ss >> trin >> trout >> tlength;
            gTargetRadius = trout / 2;
            // std::cout << "Info: Read Target Radius = " << gTargetRadius << " mm" << std::endl;
        } else if (name == "TargetHolder") {
            double thrin, throut, thlength;
            ss >> thrin >> throut >> thlength;
            gTargetHolderRadius = throut / 2;
            // std::cout << "Info: Read Target Holder Radius = " << gTargetHolderRadius << " mm" << std::endl;
        }
    }
}

void TPC_pad_template(TH2Poly* h)
{
    Double_t X[5];
    Double_t Y[5];
    for (Int_t l = 0; l < padHelper::NumOfLayersTPC; ++l) {
        Double_t pLength = padHelper::padParameter[l][5];
        Double_t st = (180. - (360. / padHelper::padParameter[l][3]) * padHelper::padParameter[l][1] / 2.);
        Double_t sTheta = (-1 + st / 180.) * TMath::Pi();
        Double_t dTheta = (360. / padHelper::padParameter[l][3]) / 180. * TMath::Pi();
        Double_t cRad = padHelper::padParameter[l][2];
        Int_t nPad = padHelper::padParameter[l][1];
        for (Int_t j = 0; j < nPad; ++j) {
            Double_t x_std[5], y_std[5];

            x_std[1] = (cRad + (pLength / 2.)) * TMath::Cos(j * dTheta + sTheta);
            x_std[2] = (cRad + (pLength / 2.)) * TMath::Cos((j + 1) * dTheta + sTheta);
            x_std[3] = (cRad - (pLength / 2.)) * TMath::Cos((j + 1) * dTheta + sTheta);
            x_std[4] = (cRad - (pLength / 2.)) * TMath::Cos(j * dTheta + sTheta);
            x_std[0] = x_std[4];

            y_std[1] = (cRad + (pLength / 2.)) * TMath::Sin(j * dTheta + sTheta);
            y_std[2] = (cRad + (pLength / 2.)) * TMath::Sin((j + 1) * dTheta + sTheta);
            y_std[3] = (cRad - (pLength / 2.)) * TMath::Sin((j + 1) * dTheta + sTheta);
            y_std[4] = (cRad - (pLength / 2.)) * TMath::Sin(j * dTheta + sTheta);
            y_std[0] = y_std[4];

            // Rotation 180 deg around center: (z, x) -> (-z, +x) after axis flip
            for (Int_t k = 0; k < 5; ++k) {
                Double_t z_temp = x_std[k] + gTpcPadCenterZ;
                Double_t x_temp = y_std[k];
                X[k] = -1.0 * z_temp;
                Y[k] = x_temp;
            }

            h->AddBin(5, X, Y);
        }
    }
    h->SetMinimum(0.1);
    h->SetMaximum(1000.0);
}

void draw_track(Int_t n_rand)
{
    gStyle->SetOptStat(0);
    Int_t font = 132; // 132: Times New Roman, 42: Helvetica
    gStyle->SetTextFont(font);
    gStyle->SetLabelFont(font, "XYZ");
    gStyle->SetTitleFont(font, "XYZ");
    gStyle->SetTitleFont(font, "");
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetTitleSize(0.04, "XYZ");
    gStyle->SetTitleSize(0.06, "t");
    

    if (!c_canvas_tpc) { c_canvas_tpc = new TCanvas("Event Display", "TPC Event Display", 800, 800); }

    TFile* f = new TFile(g_rootfile_path.Data());
    TTreeReader reader("g4s2s", f);
    TTreeReaderValue<std::vector<TParticle>> tpcHits(reader, "TPC");
    TTreeReaderValue<Int_t> evnum(reader, "evnum");

    Int_t tot_num = reader.GetEntries();
    Int_t n = (n_rand == -1) ? 0 : n_rand % tot_num;
    reader.SetEntry(n);

    std::cout << "Drawing Event: " << *evnum << " (Entry: " << n << ")" << std::endl;

    auto* h_tpc_2d =
        new TH2Poly("h_tpc_2d", Form("TPC Hits (Event %d);z [mm];x [mm]", *evnum), -280.0, 280.0, -280.0, 280.0);
    TPC_pad_template(h_tpc_2d);

    int hitCount = 0;
    for (const auto& hit : *tpcHits) {
        const auto coords = ProjectToDisplay(TVector3(hit.Vx(), hit.Vy(), hit.Vz()));
        const Double_t edep = hit.Energy();
        h_tpc_2d->Fill(coords.first, coords.second, edep > 0 ? edep : 1.0);
        hitCount++;
    }
    std::cout << "TPC Hits: " << hitCount << std::endl;

    c_canvas_tpc->cd(1);
    c_canvas_tpc->SetRightMargin(0.13);
    c_canvas_tpc->SetLeftMargin(0.13);
    gPad->SetLogz();
    h_tpc_2d->GetXaxis()->CenterTitle(kTRUE);
    h_tpc_2d->GetYaxis()->CenterTitle(kTRUE);
    h_tpc_2d->Draw("colz");

    const auto targetCoords = ProjectToDisplay(TVector3(0., 0., gTpcPadCenterZ));
    TEllipse* targetCircle = new TEllipse(targetCoords.first, targetCoords.second, gTargetRadius);
    targetCircle->SetFillStyle(0);
    targetCircle->SetLineColor(kRed);
    targetCircle->SetLineWidth(3);
    targetCircle->Draw("same");

    const auto targetHolderCoords = ProjectToDisplay(TVector3(0., 0., gTpcPadCenterZ));
    TEllipse* targetHolderCircle = new TEllipse(targetHolderCoords.first, targetHolderCoords.second, gTargetHolderRadius);
    targetHolderCircle->SetFillStyle(0);
    targetHolderCircle->SetLineColor(kMagenta);
    targetHolderCircle->SetLineWidth(3);
    targetHolderCircle->Draw("same");

    gPad->Modified();
    gPad->Update();

    f->Close();
    delete f;
}

void event(Int_t n = -1)
{
    if (gTargetZ_Global == 0.0 && gTpcPadCenterZ == 0.0) { read_geometry(); }

    if (n == -1) {
        std::random_device rd;
        std::mt19937 gen(rd());
        Int_t n_rand = gen();
        while (n_rand < 0) n_rand = gen();
        draw_track(n_rand);
    } else {
        draw_track(n);
    }
}

void set_path(TString rootfile_path)
{
    g_rootfile_path = rootfile_path;
}
