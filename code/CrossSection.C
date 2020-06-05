/////Amilkar Quintero
////Temple University
////August 2019
////
//// Unfold the decorrelation angle
///

#define REBIN_DET 1
#define REBIN_HAD 4
#define NSCAN 100

// Change path to the root-files at lines 276+


#define DO_SCALE 1

#if DO_SCALE
# define SCALE_2_7 1
# define SCALE_7_12 10
# define SCALE_12_30 100
#else
# define SCALE_2_7 1
# define SCALE_7_12 1
# define SCALE_12_30 1
#endif

#include <string>
#include <getopt.h>

#include "Riostream.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TString.h"
#include "TDatime.h"
#include "TMath.h"
#include "TStyle.h"
#include <TVector.h>
#include <TLatex.h>
#include <TMatrixDSymEigen.h>
#include <TApplication.h>
#include "TUnfoldDensity.h"
#include <TSpline.h>
#include <libconfig.h++>

#undef assert
#define assert(expr) \
  if (!(expr)) { \
    std::cout << "ASSERTION FALIED: " << __func__ << ":" << __LINE__ << ": " \
              << #expr << std::endl; \
    abort(); \
  }

static double MinY = 0.0001;
static double MaxY = 70;

static int g_counter = 0;

static int markercnt = 0;
static int markers[] = { 24, 25, 26 };

static
bool g_show_stat_unc = true;

static
bool g_fix_tau = false;

static
bool g_limit_x = true;

static
bool g_cut_low_dphi = false;

#define IDX_PT 0
#define IDX_Q2 1
static int g_var_idx = 0;

static
bool g_use_project_z = false;

static
bool g_plot_unc_hists = false;

static struct {
  const char *name, *texname, *legnames[3];
} g_var[2] = {
  { "Pt",
    "P_{T}",
    {
      "2.5 < p_{T} <  7 GeV",
      "  7 < p_{T} < 12 GeV",
      " 12 < p_{T} < 30 GeV",
    },
  },
  {
    "Q2",
    "Q^{2}",
    {
      " 10 < Q^{2} <  50 GeV",
      " 50 < Q^{2} < 100 GeV",
      "100 < Q^{2} < 350 GeV",
    },
  },
};

#define VAR (g_var[g_var_idx])

TH2D* project_z(TH3D *h3, int zmin_min, int zbin_max)
{
  int nbinsx = h3->GetNbinsX();
  double xmin = h3->GetXaxis()->GetBinLowEdge(1);
  double xmax = h3->GetXaxis()->GetBinUpEdge(nbinsx);

  int nbinsy = h3->GetNbinsY();
  double ymin = h3->GetYaxis()->GetBinLowEdge(1);
  double ymax = h3->GetYaxis()->GetBinUpEdge(nbinsy);

  TH2D *h2 = new TH2D(
    Form("projection-z-%d", g_counter++), "",
    nbinsx, xmin, xmax,
    nbinsy, ymin, ymax
  );

  for (int ix = 1; ix <= nbinsx; ++ix) {
    for (int iy = 1; iy <= nbinsy; ++iy) {
      double bin = 0;
      double err2 = 0;
      for (int iz = zmin_min; iz <= zbin_max; ++iz) {
        bin += h3->GetBinContent(ix, iy, iz);
        err2 += pow(h3->GetBinError(ix, iy, iz), 2);
      }
      h2->SetBinContent(ix, iy, bin);
      h2->SetBinError(ix, iy, sqrt(err2));
    }
  }

  return h2;
}

TH2D*
get_migrations(TH3D *h3, int binmin, int binmax)
{
  if (g_use_project_z) {
    return project_z(h3, binmin, binmax);
  } else {
    h3->GetZaxis()->SetRange(binmin, binmax);
    return (TH2D*)h3->Project3D("yxe");
  }
}

int
next_marker()
{
  if (markercnt == 3)
    markercnt = 0;
  return markers[markercnt++];
}

TH1D*
cut_low_dphi_1d(TH1 *hin, double at)
{
  int nbins = hin->GetNbinsX();
  int ibin = hin->GetXaxis()->FindBin(at);

  double from = hin->GetXaxis()->GetBinLowEdge(ibin);
  double to = hin->GetXaxis()->GetBinUpEdge(nbins);

  TH1D *hout = new TH1D {
    Form("%s-cutdphi-%d", hin->GetName(), g_counter++), // name
      "", // title
      nbins - ibin + 1, // N bins
      from, to // range
  };

  for (int i = 1; i <= nbins; ++i) {
    double bin = hin->GetBinContent(i);
    double err = hin->GetBinError(i);
    if (i < ibin) {
      double outbin = hout->GetBinContent(0);
      double outerr = hout->GetBinError(0);
      hout->SetBinContent(0, outbin + bin);
      hout->SetBinError(0, sqrt(outerr*outerr + err*err));
    } else {
      hout->SetBinContent(i - ibin + 1, bin);
      hout->SetBinError(i - ibin + 1, err);
    }
  }

  std::cout << "=== underflow: " << hout->GetBinContent(0) << " +/- " << hout->GetBinError(0) << std::endl;
  return hout;
}

TH2D*
cut_low_dphi_2d(TH2 *hin, double at)
{
  // x axis
  int nbinsx = hin->GetNbinsX();
  int ibinx = hin->GetXaxis()->FindBin(at);
  // y axis
  int nbinsy = hin->GetNbinsY();
  int ibiny = hin->GetYaxis()->FindBin(at);

  // x axis
  double fromx = hin->GetXaxis()->GetBinLowEdge(ibinx);
  double tox = hin->GetXaxis()->GetBinUpEdge(nbinsx);
  // y axis
  double fromy = hin->GetYaxis()->GetBinLowEdge(ibiny);
  double toy = hin->GetYaxis()->GetBinUpEdge(nbinsy);

  TH2D *hout = new TH2D {
    Form("%s-cutdphi-%d", hin->GetName(), g_counter++), // name
      "", // title
      nbinsx - ibinx + 1, // N x bins
      fromx, tox, // x axis range
      nbinsy - ibiny + 1, // N y bins
      fromy, toy // y axis range
  };

  for (int ix = 1; ix <= nbinsx; ++ix) {
    for (int iy = 1; iy  <= nbinsy; ++iy) {
      double bin = hin->GetBinContent(ix, iy);
      double err = hin->GetBinError(ix, iy);
      if (ix < ibinx && iy < ibiny) {
        double outbin = hout->GetBinContent(0);
        double outerr = hout->GetBinError(0);
        hout->SetBinContent(0, outbin + bin);
        hout->SetBinError(0, sqrt(outerr*outerr + err*err));
      } else if (ix < ibinx) {
        int binno = hout->GetBin(0, iy - ibiny + 1);
        double outbin = hout->GetBinContent(binno);
        double outerr = hout->GetBinError(binno);
        hout->SetBinContent(binno, outbin + bin);
        hout->SetBinError(binno, sqrt(outerr*outerr + err*err));
      } else if (iy < ibiny) {
        int binno = hout->GetBin(ix - ibinx + 1, 0);
        double outbin = hout->GetBinContent(binno);
        double outerr = hout->GetBinError(binno);
        hout->SetBinContent(binno, outbin + bin);
        hout->SetBinError(binno, sqrt(outerr*outerr + err*err));
      } else {
        hout->SetBinContent(ix - ibinx + 1, iy - ibiny + 1, bin);
        hout->SetBinError(ix - ibinx + 1, iy - ibiny + 1, err);
      }
    }
  }

  std::cout << "=== underflow: " << hout->GetBinContent(0) << " +/- " << hout->GetBinError(0) << std::endl;
  return hout;
}

TH1*
get_true(double scale, double zmin, double zmax, TH2D *h2tr)
{
  static int hist_id = 0;

  TH1D *htr = h2tr->ProjectionX(
      Form("htr_%s_%d_%d-%d", VAR.name, int(zmin), int(zmax), hist_id++),
      zmin < 0 ? 1                 : h2tr->GetYaxis()->FindBin(zmin),
      zmax < 0 ? h2tr->GetNbinsY() : h2tr->GetYaxis()->FindBin(zmax),
      "e" // compute errors
  );
  htr->SetLineColor(2);
  htr->SetMarkerColor(2);
  //htr->Sumw2();
  htr->Rebin(REBIN_HAD);
  //std::cout << "=== N true bins (plotted) = " << htr->GetNbinsX() << std::endl;
  htr->Scale(scale / htr->Integral());
  htr->SetLineWidth(1);
  htr->GetYaxis()->SetRangeUser(MinY,MaxY);
  //htr->SetTitle("#bf{#font[22]{#scale[1.5]{ZEUS Preliminary}}}");
  htr->SetTitle("");
  htr->GetXaxis()->SetTitle("#scale[1.3]{#Delta#phi (rad)}     ");
  htr->GetYaxis()->SetTitle("#scale[1.3]{#frac{1}{#sigma} #frac{d#sigma}{d#Delta#phi} (rad^{-1})}    ");
  htr->GetYaxis()->SetLabelSize(0.04);
  htr->GetYaxis()->SetTitleOffset(1.2);

  TAxis* a = htr->GetXaxis();
  if (g_limit_x) {
    //std::cout << "!!!!!!!!!!!!!!!!!!!!!!!\n";
    a->SetRangeUser(TMath::Pi() / 2,TMath::Pi()); //In the plot at detector level this is 0.055
    //a->SetNdivisions(4, kFALSE);
    //a->SetLabelSize(0.044);
    //a->SetBit(TAxis::kLabelsHori);
    //a->ChangeLabel(1,-1,-1,-1,-1,-1,"\"#pi/2\"");
    //a->ChangeLabel(2,-1,-1,-1,-1,-1,"2#pi/3");
    //a->ChangeLabel(3,-1,-1,-1,-1,-1,"5#pi/6");
    //a->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
    //a->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
    //a->SetLabelSize(0.13);
  } else {
    a->SetNdivisions(-303);
    a->SetBit(TAxis::kLabelsHori);
    a->SetLabelSize(0.044);
  }

  return htr;
  //return cut_low_dphi_1d(htr, M_PI_2);
}

TH1*
unfold_primary(
    const char *prefix,
    TUnfoldDensity &unfold,
    int *ibest = NULL,
    TGraph **lcurve = NULL,
    TSpline **logtaux = NULL,
    TSpline **logtauy = NULL)
{
  // Scan L-curve with default settings for tau:
  int best = unfold.ScanLcurve(NSCAN, 0, 0, lcurve, logtaux, logtauy);
  //int best = unfold.ScanTau(NSCAN, 0, 0, 0, TUnfoldDensity::kEScanTauRhoAvg, 0, 0, lcurve, logtaux, logtauy);
  //unfold.ScanTau(NSCAN, 0, 0, NULL);
  if (ibest)
    *ibest = best;

  std::cout << "=== get output" << std::endl;
  TH1 *unfolded = unfold.GetOutput(Form("%s-unfolded-primary-%d", prefix, g_counter++));
  return unfolded;
}

template <typename Matrix> void
TH2_to_Matrix(const TH2 *h2, Matrix &m)
{
  size_t nrow = h2->GetNbinsY();
  size_t ncol = h2->GetNbinsX();

  // copy histogram to matrix
  m.ResizeTo(nrow, ncol);
  for (size_t i = 0; i < nrow; ++i) {
    for (size_t j = 0; j < ncol; ++j)
      m(i, j) = h2->GetBinContent(i+1, j+1);
  }
}

std::pair<double, double>
ematrix_eigenvalues(const TH2 *ematrix_hist)
{
  // conver hystogram to matrix
  TMatrixTSym<double> ematrix;
  TH2_to_Matrix(ematrix_hist, ematrix);

  // get eigen values
  TVectorD lambda = TMatrixDSymEigen(ematrix).GetEigenValues();

  // find max and min eigen values
  double lambda_min = DBL_MAX, lambda_max = DBL_MIN;
  for (size_t i = 0; i < lambda.GetNoElements(); ++i) {
    lambda_min = std::min(lambda_min, lambda[i]);
    lambda_max = std::max(lambda_max, lambda[i]);
  }
  assert(lambda_max > 0);
  assert(lambda_min > 0);

  return { lambda_min, lambda_max };
}

TH1*
unfold_secondary(
    const char *prefix,
    TUnfoldDensity &unfold,
    int *ibest = NULL,
    TGraph **lcurve = NULL,
    TSpline **logtaux = NULL,
    TSpline **logtauy = NULL)
{
  // Determine tau_min and tau_max:
  TH2 *ematrix = unfold.GetEmatrixInput(Form("%s-ematrix-%d", prefix, g_counter++));

  auto [lambda_min, lambda_max] = ematrix_eigenvalues(ematrix); // get min and max eigenvalues of error-matrix
  double tau_max = 1. / sqrt(lambda_min);
  double tau_min = 1. / sqrt(lambda_max);

  std::cout << "lambda_min = " << lambda_min << ", " << "tau_max = " << tau_max << std::endl;
  std::cout << "lambda_max = " << lambda_max << ", " << "tau_min = " << tau_min << std::endl;

  // Scan L-curve:
  int best = unfold.ScanLcurve(NSCAN, tau_min, tau_max, lcurve, logtaux, logtauy);
  if (ibest)
    *ibest = best;
  //unfold.ScanTau(NSCAN, tau_min, tau_max, NULL);

  // Plot unfolded data:
  TH1 *unfolded = unfold.GetOutput(Form("%s-unfolded-%d", prefix, g_counter++));
  return unfolded;
}

void
add_deltas(TH1 *h, TH1 *src)
{
  assert(h->GetNbinsX() == src->GetNbinsX());
  int n = src->GetNbinsX();
  for (int i = 1; i <= n; ++i) {
    double sys_err = std::abs(src->GetBinContent(i));
    double sys_err_rel = std::abs(sys_err) / h->GetBinContent(i) * 100.;
    std::cout << "sys. err. = " << sys_err << " (abs), "
                                << sys_err_rel << " (rel %)" << std::endl;
    double old_err = h->GetBinError(i);
    double new_err = sqrt(old_err*old_err + sys_err*sys_err);
    h->SetBinError(i, new_err);
  }
}

struct LCurveData {
  TGraph *lcurve;
  TSpline *logtaux, *logtauy;
  double x, y;
  int ibest;

  void draw()
  {
    lcurve->DrawClone("APLE");
    double t[1], px[1], py[1];
    logtaux->GetKnot(ibest, t[0], px[0]);
    logtauy->GetKnot(ibest, t[0], py[0]);
    auto lcurve_point = new TGraph { 1, px, py };
    lcurve_point->SetMarkerStyle(3);
    lcurve_point->SetMarkerColor(kRed);
    lcurve_point->DrawClone("p same");
  }
};

struct UnfoldData {
  TH2 *h2mig;
  TH2 *h2cov;
  TH2 *h2cor;
};

TH1*
unfold_data(
    double scale,
    double zmin, double zmax,
    TH2 *h2da, TH3D *h3mc, const std::vector<TH3D*> &hsys,
    UnfoldData *unfold_data = NULL,
    LCurveData *lcurve_data = NULL
)
{
  static int sys_cnt = 0;

  TH1D *hda = h2da->ProjectionX(
      Form("hda_%s_%d_%d-%d", VAR.name, int(zmin), int(zmax), g_counter++), // title
      zmin < 0 ? 1                 : h2da->GetYaxis()->FindBin(zmin), // low bin
      zmax < 0 ? h2da->GetNbinsY() : h2da->GetYaxis()->FindBin(zmax),  // high bin
      "e" // compute errors
  );
  //hda->Sumw2();
  //hda->Scale(1 / hda->Integral());
  hda->Rebin(REBIN_DET);

  int zmin_bin = zmin < 0 ? 1                 : h3mc->GetZaxis()->FindBin(zmin);
  int zmax_bin = zmax < 0 ? h3mc->GetNbinsZ() : h3mc->GetZaxis()->FindBin(zmax);
  std::cout << "Z-vals: " <<     zmin << " .. " <<     zmax << std::endl;
  std::cout << "Z-bins: " << zmin_bin << " .. " << zmax_bin << std::endl;

  TH3D *h3mc_clone = (TH3D*)h3mc->Clone();
  //h3mc_clone->GetZaxis()->SetRange(zmin_bin, zmax_bin);
  h3mc_clone->RebinX(REBIN_DET);
  h3mc_clone->RebinY(REBIN_HAD);
  //TH2D *h2mc = (TH2D*)h3mc_clone->Project3D("yxe"); // Correlation Matrix Hadron vs detector (+ compute errors)
  //TH2D *h2mc = project_z(h3mc_clone, zmin_bin, zmax_bin);
  TH2D *h2mc = get_migrations(h3mc_clone, zmin_bin, zmax_bin);
  std::cout << "N entries in cov. matrix: " << h2mc->GetEntries() << std::endl;

  double hda_integral = hda->Integral();
  double h2mc_integral = h2mc->Integral();

  if (g_cut_low_dphi) {
    // rebin histograms so that all bins below Pi/2 are moved into underflow bin
    hda  = cut_low_dphi_1d(hda, M_PI_2);
    h2mc = cut_low_dphi_2d(h2mc, M_PI_2);
  }

  TUnfoldDensity unfold {
      // # Migration matrix
      h2mc,

      // # Axis layout
      TUnfold::kHistMapOutputVert,

      // # Regularization mode
      TUnfold::kRegModeCurvature,
      //TUnfold::kRegModeDerivative,
      //TUnfold::kRegModeSize,

      // # Constraint
      TUnfold::kEConstraintArea,
      //TUnfold::kEConstraintNone,

      // # Density mode
      TUnfoldDensity::kDensityModeBinWidthAndUser,
      //TUnfoldDensity::kDensityModeBinWidth,
      //TUnfoldDensity::kDensityModeNone,

      // # Output binning
      NULL, // default

      // # Input binning
      NULL, // default

      // # Selection of regularized distribution
      NULL, // default

      // # Axis steering
      // Options:
      //   u : exclude underflow bin from derivatives along this axis
      //   o : exclude overflow bin from derivatives along this axis
      //   U : exclude underflow bin
      //   O : exclude overflow bin
      //   b : use bin width for derivative calculation
      //   B : same as 'b', in addition normalize to average bin width
      //   N : completely exclude derivatives along this axis
      //   p : axis is periodic (e.g. azimuthal angle), so include derivatives
      //       built from combinations involving bins at both ends of the axis
      //       "wrap around"
      //"*[UOB]", // default
      "*[uOB]" // use underflow
  };
  unfold.SetInput(hda);

  {
    int isys = 1;
    for (auto h3 : hsys) {
      assert(h3);
      std::cout << "adding source of systematics: " << h3->GetName() << std::endl;
      TH3D *h3_clone = (TH3D*)h3->Clone();
      //h3_clone->GetZaxis()->SetRange(zmin_bin, zmax_bin);
      h3_clone->RebinX(REBIN_DET);
      h3_clone->RebinY(REBIN_HAD);
      //TH2D *h2 = (TH2D*)h3_clone->Project3D("yxe"); // Correlation Matrix Hadron vs detector (+ compute errors)
      TH2D *h2 = get_migrations(h3_clone, zmin_bin, zmax_bin);
      if (g_cut_low_dphi)
        h2 = cut_low_dphi_2d(h2, M_PI_2);

      unfold.AddSysError(h2,
          Form("sys-%d", isys++),
          TUnfold::kHistMapOutputVert,
          TUnfoldSys::kSysErrModeMatrix
      );
    }
  }

  int ibest;
  TGraph *lcurve;
  TSpline *logtaux, *logtauy;

  auto hda_unfold = unfold_primary(Form("%d_%d-%d", int(zmin), int(zmax), g_counter++), unfold, &ibest, &lcurve, &logtaux, &logtauy);
  if (g_fix_tau) {
    std::cout << "run secondary unfolding" << std::endl;
    hda_unfold = unfold_secondary(Form("%d_%d-%d-sec", int(zmin), int(zmax), g_counter++), unfold, &ibest, &lcurve, &logtaux, &logtauy);
  }
  std::cout << "\e[3mChi2 (A) / Ndf\e[0m: " << unfold.GetChi2A()/unfold.GetNdf() << std::endl;
  std::cout << "\e[3mChi2 (L) / Ndf\e[0m: " << unfold.GetChi2L()/unfold.GetNdf() << std::endl;
  std::cout << "\e[3mChi2 (full) / Ndf\e[0m: " << (unfold.GetChi2A() + unfold.GetChi2L()) / unfold.GetNdf() << std::endl;
  hda_unfold->SetMarkerColor(1);
  hda_unfold->SetMarkerStyle(next_marker());
  hda_unfold->SetMarkerSize(0.8);
  hda_unfold->SetFillColor(1);

  if (lcurve_data) {
    lcurve_data->ibest = ibest;
    lcurve_data->lcurve = lcurve;
    lcurve_data->logtaux = logtaux;
    lcurve_data->logtauy = logtauy;
    lcurve_data->x = unfold.GetLcurveX();
    lcurve_data->y = unfold.GetLcurveY();
  }

  if (unfold_data) {
    h2mc->SetName(Form("h2mig-%d", g_counter++));
    unfold_data->h2mig = h2mc;
    unfold_data->h2cov = unfold.GetEmatrixTotal(Form("h2cov-%d", g_counter++));
    unfold_data->h2cor = unfold.GetRhoIJtotal(Form("h2cor-%d", g_counter++));
  }

  // remove statistical uncertainties
  if (not g_show_stat_unc) {
    std::cout << "removing statistical uncertainties" << std::endl;
    int n = hda_unfold->GetNbinsX();
    for (int i = 1; i <= n; ++i)
      hda_unfold->SetBinError(i, 0);
  }

  // add deltas from sys-sources
  {
    std::cout << "adding deltas from sys. sources" << std::endl;
    int isys = 1;
    for (auto _ : hsys) {
      std::cout << "- sys source #" << isys << std::endl;
      std::string src_name = Form("sys-%d", isys);
      std::string h_name = Form("sys-%d_%d", isys, sys_cnt++);
      isys++;
      auto h1sys = unfold.GetDeltaSysSource(src_name.c_str(), h_name.c_str());
      add_deltas(hda_unfold, h1sys);
    }
  }

  // add deltas from uncertanty on tau
  {
    auto h1tau = unfold.GetDeltaSysTau(Form("sys-tau-%d", sys_cnt++));
    if (h1tau) {
      std::cout << "adding deltas from tau uncertainty" << std::endl;
      add_deltas(hda_unfold, h1tau);
    }
  }

  hda_unfold->Scale(scale / hda_integral);

  return hda_unfold;
}

TH1D* get_unc_hist(const char *name, TH1* hin)
{
  auto hout = (TH1D*)hin->Clone(name);
  int n = hout->GetNbinsX();
  for (int i = 1; i <= n; ++i) {
    double binval = hout->GetBinContent(i);
    double binerr = hout->GetBinError(i);
    hout->SetBinContent(i, binerr / binval * 100.);
    hout->SetBinError(i, 0);
  }

  hout->GetYaxis()->SetTitle("unc. (%)");
  hout->GetYaxis()->SetTitleSize(0.1);
  hout->GetYaxis()->SetTitleOffset(0.2);
  hout->GetYaxis()->SetLabelSize(0.075);

  if (g_limit_x)
    hout->GetXaxis()->SetRangeUser(TMath::Pi() / 2, TMath::Pi());
  hout->GetXaxis()->SetTitle("#Delta#phi (rad)");
  hout->GetXaxis()->SetTitleSize(0.1);
  hout->GetXaxis()->SetTitleOffset(0.4);
  hout->GetXaxis()->SetLabelSize(0.075);

  return hout;
}

void
unc_from_mc_to_data(TH1 *hmc, TH1 *hda)
{
  int nbins = hmc->GetNbinsX();
  assert(nbins == hda->GetNbinsX());
  for (int i = 1; i <= nbins; ++i) {
    double emc = hmc->GetBinError(i);
    double eda = hda->GetBinError(i);
    double etot = sqrt(emc*emc + eda*eda);
    hda->SetBinError(i, etot);
  }
}

//--->START MAIN PROGRAM
//_____________________________________________________________________________
void CrossSection(const char *dpath, const char *mpath, const char *tpath,
    const char *out_path, const std::vector<std::string> &syspath)
{
  //  gROOT->cd();
  TDatime now;                                          //Set time in Root
  now.Print();
  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);


  /////Read files --1-->
  TFile *fda = new TFile(dpath);
  TFile *fmc = new TFile(mpath);
  TFile *ftr = new TFile(tpath);
  /////END Read files <--1--

  if (fda->IsZombie() || fmc->IsZombie() || ftr->IsZombie()) {
    std::cerr << "Error: failed to open some data files" << std::endl;
    exit(EXIT_FAILURE);
  }


  ////Define Histograms --2-->
  TH2D* h2daPtDecorrPhi[5];
  TH3D* h3mcPtDecorrPhi[5];
  TH2D* h2trPtDecorrPhi[5];
  TH2D* h2daQ2DecorrPhi[5];
  TH3D* h3mcQ2DecorrPhi[5];
  TH2D* h2trQ2DecorrPhi[5];

  TH2D* h2mcPtDecorrPhi[5];
  TH2D* h2mcQ2DecorrPhi[5];

  for(Int_t ijet = 0; ijet < 5; ijet++) {
    h2daPtDecorrPhi[ijet] = (TH2D*)fda->Get(Form("h2PtDecorrPhi_%d",ijet));
    h3mcPtDecorrPhi[ijet] = (TH3D*)fmc->Get(Form("h3DecorrPhiHadDecPt_%d",ijet));
    h2trPtDecorrPhi[ijet] = (TH2D*)ftr->Get(Form("h2PtDecorrPhi_%d",ijet));
    h2daQ2DecorrPhi[ijet] = (TH2D*)fda->Get(Form("h2Q2DecorrPhi_%d",ijet));
    h3mcQ2DecorrPhi[ijet] = (TH3D*)fmc->Get(Form("h3DecorrPhiHadDecQ2_%d",ijet));
    h2trQ2DecorrPhi[ijet] = (TH2D*)ftr->Get(Form("h2Q2DecorrPhi_%d",ijet));

    h2mcPtDecorrPhi[ijet] = (TH2D*)fmc->Get(Form("h2PtDecorrPhi_%d",ijet));
    h2mcQ2DecorrPhi[ijet] = (TH2D*)fmc->Get(Form("h2Q2DecorrPhi_%d",ijet));
  }

  std::vector<TH3D*> h3mcPtDecorrPhi_sys[5];
  std::vector<TH3D*> h3mcQ2DecorrPhi_sys[5];
  for (auto &path : syspath) {
    std::cout << "loading systematics: " << path << std::endl;
    TFile *fsys = new TFile(path.c_str());
    assert(!fsys->IsZombie());
    for(Int_t ijet = 0; ijet < 5; ijet++) {
      h3mcPtDecorrPhi_sys[ijet].push_back(
          (TH3D*)fsys->Get(Form("h3DecorrPhiHadDecPt_%d",ijet))
      );
      h3mcQ2DecorrPhi_sys[ijet].push_back(
          (TH3D*)fsys->Get(Form("h3DecorrPhiHadDecQ2_%d",ijet))
      );
    }
  }
  ////End Define Histograms <--2--


  /////////////
  TH2D **h2tr, **h2da, **h2mc;
  TH3D **h3mc;
  std::vector<TH3D*> *h3mc_sys;
  double bins[6];
  if (g_var_idx == IDX_PT) {
    h2tr = h2trPtDecorrPhi;
    h2da = h2daPtDecorrPhi;
    h2mc = h2mcPtDecorrPhi;
    h3mc = h3mcPtDecorrPhi;
    h3mc_sys = h3mcPtDecorrPhi_sys;
    bins[0] =  2.5; bins[1] =  7;
    bins[2] =  7.5; bins[3] = 12;
    bins[4] = 12.5; bins[5] = 30;
  } else /* g_var_idx == IDX_Q2 */ {
    h2tr = h2trQ2DecorrPhi;
    h2da = h2daQ2DecorrPhi;
    h2mc = h2mcQ2DecorrPhi;
    h3mc = h3mcQ2DecorrPhi;
    h3mc_sys = h3mcQ2DecorrPhi_sys;
    // 34 bins spanning [10, 350) GeV^2
    // => width of a bin is 10 GeV^2
    bins[0] =  11; bins[1] =  41; // [10, 50)
    bins[2] =  51; bins[3] = 101; // [50, 100)
    bins[4] = 101; bins[5] = 341; // [100, 350)
  }

  // ========================================================================
  //                               One+ jet integrated
  //
  TCanvas *cint = new TCanvas;

  cint->cd();
  TPad *cint_main = new TPad { "cint_main", "", 1.0, 1.0, 0.0, 0.3 };
  cint_main->SetBottomMargin(0);
  cint_main->SetLogy();
  cint_main->Draw();
  cint_main->cd();

  TLegend leg_int { 0.7, 0.13, 0.9, 0.39 };
  leg_int.SetBorderSize(0);

  // true
  TH1 *htr_int = get_true(1, -1, -1, h2tr[0]);
  htr_int->GetYaxis()->SetRangeUser(std::pow(10, -2.5), 1);
  htr_int->SetTitle("");
  htr_int->DrawClone("hist");
  leg_int.AddEntry(htr_int, "MC, hadron level", "lp");

  // data
  TH1 *hda_int = get_true(1, -1, -1, h2da[0]);
  hda_int->SetLineColor(kBlue);
  hda_int->SetMarkerColor(kBlue);
  hda_int->Rebin(3);
  hda_int->DrawClone("same hist e1");
  leg_int.AddEntry(hda_int, "Data", "lp");

  // MC
  TH1 *hmc_int = get_true(1, -1, -1, h2mc[0]);
  hmc_int->SetLineColor(kGreen);
  hmc_int->SetMarkerColor(kGreen);
  hmc_int->Rebin(3);
  hmc_int->DrawClone("same hist e1");
  leg_int.AddEntry(hmc_int, "MC, detector level", "lp");

  // unfolded data
  std::cout << "=== unfolding data (integrated)" << std::endl;
  UnfoldData unfold_data_int;
  LCurveData lcurve_data_int;
  auto hda_int_unfold = unfold_data(1, -1, -1, h2da[0], h3mc[0], h3mc_sys[0], &unfold_data_int, &lcurve_data_int);
  hda_int_unfold->SetMarkerStyle(24);
  hda_int_unfold->DrawClone("same e1");
  leg_int.AddEntry(hda_int_unfold, "Data, unfolded", "p");


  // unfolded MC
  std::cout << "=== unfolding MC (integrated)" << std::endl;
  auto hmc_int_unfold = unfold_data(1, -1, -1, h2mc[0], h3mc[0], h3mc_sys[0]);
  //hmc_int_unfold->SetMarkerStyle(26);
  //hmc_int_unfold->DrawClone("same e1");
  //leg_int.AddEntry(hmc_int_unfold, "MC, unfolded", "p");

  leg_int.DrawClone("same");

  cint->cd();
  TPad *cint_rdamc = new TPad { "cint_rdamc", "", 0.0, 0.3, 1.0, 0.2 };
  cint_rdamc->SetGridy();
  cint_rdamc->SetTopMargin(0);
  cint_rdamc->SetBottomMargin(0);
  cint_rdamc->Draw();
  cint_rdamc->cd();

  auto h_int_rdamc = (TH1D*)hda_int->Clone();
  h_int_rdamc->SetTitle("");
  h_int_rdamc->SetMarkerStyle(20);
  h_int_rdamc->SetMarkerSize(1);
  h_int_rdamc->SetMarkerColor(kBlack);
  h_int_rdamc->SetLineColor(kBlack);
  h_int_rdamc->Divide(hmc_int);
  h_int_rdamc->GetYaxis()->SetTitle("#scale[3]{Data/MC_{det}}");
  h_int_rdamc->GetXaxis()->SetTitle("#Delta#phi (rad)      ");
  h_int_rdamc->GetXaxis()->SetTitleSize(0.13);
  h_int_rdamc->GetYaxis()->SetTitleSize(0.05);
  //h_int_rdamc->GetXaxis()->SetTitleOffset(0.5);
  h_int_rdamc->GetYaxis()->SetTitleOffset(0.75);
  h_int_rdamc->GetYaxis()->SetRangeUser(0.75, 1.25);
  h_int_rdamc->GetXaxis()->SetRangeUser(M_PI_2, M_PI);
  h_int_rdamc->GetYaxis()->SetNdivisions(7);
  h_int_rdamc->GetYaxis()->SetLabelSize(0.12);
  {
    TAxis* a = h_int_rdamc->GetXaxis();
    a->SetNdivisions(-303);
    a->SetBit(TAxis::kLabelsHori);
    a->ChangeLabel(1,-1,-1,-1,-1,-1,"#pi/2");
    a->ChangeLabel(2,-1,-1,-1,-1,-1,"2#pi/3");
    a->ChangeLabel(3,-1,-1,-1,-1,-1,"5#pi/6");
    a->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
    a->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
    a->SetLabelSize(0.13);
  }
  h_int_rdamc->DrawClone("pe");

  cint->cd();
  TPad *cint_rdatr = new TPad { "cint_rdatr", "", 0.0, 0.2, 1.0, 0.1 };
  cint_rdatr->SetGridy();
  cint_rdatr->SetTopMargin(0);
  cint_rdatr->Draw();
  cint_rdatr->cd();

  auto h_int_rdatr = (TH1D*)hda_int_unfold->Clone();
  h_int_rdatr->SetTitle("");
  h_int_rdatr->SetMarkerStyle(20);
  h_int_rdatr->SetMarkerSize(1);
  h_int_rdatr->SetMarkerColor(kBlack);
  h_int_rdatr->SetLineColor(kBlack);
  h_int_rdatr->Divide(htr_int);
  h_int_rdatr->GetYaxis()->SetTitle("#scale[3]{Data_{had}/MC_{had}}");
  //h_int_rdatr->GetYaxis()->LabelsOption("v");
  h_int_rdatr->GetXaxis()->SetTitle("#Delta#phi (rad)      ");
  h_int_rdatr->GetXaxis()->SetTitleSize(0.13);
  h_int_rdatr->GetYaxis()->SetTitleSize(0.05);
  h_int_rdatr->GetXaxis()->SetTitleOffset(0.5);
  //h_int_rdatr->GetYaxis()->SetTitleOffset(0.99);
  h_int_rdatr->GetYaxis()->SetRangeUser(0.75, 1.25);
  h_int_rdatr->GetXaxis()->SetRangeUser(M_PI_2, M_PI);
  h_int_rdatr->GetYaxis()->SetNdivisions(7);
  h_int_rdatr->GetYaxis()->SetLabelSize(0.12);
  {
    TAxis* a = h_int_rdatr->GetXaxis();
    a->SetNdivisions(-303);
    a->SetBit(TAxis::kLabelsHori);
    a->ChangeLabel(1,-1,-1,-1,-1,-1,"#pi/2");
    a->ChangeLabel(2,-1,-1,-1,-1,-1,"2#pi/3");
    a->ChangeLabel(3,-1,-1,-1,-1,-1,"5#pi/6");
    a->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
    a->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
    a->SetLabelSize(0.13);
  }
  h_int_rdatr->DrawClone("pe");


  // Covariance matrix
  gStyle->SetPaintTextFormat("5.2f");
  unfold_data_int.h2cor->SetTitle("");
  unfold_data_int.h2cor->SetMaximum(1);
  unfold_data_int.h2cor->SetMinimum(-1);

  auto c_int_e = new TCanvas;
  unfold_data_int.h2cor->DrawClone("colz");
  unfold_data_int.h2cor->DrawClone("text same");

  // L-curve
  auto c_lcurve_int = new TCanvas;
  lcurve_data_int.draw();


  TCanvas *cre = new TCanvas("Fig7","Results",0,0,1120,800);//1400,1000);
  cre->cd();
  cre->Divide(2,2);


  double scale[3] = { SCALE_2_7, SCALE_7_12, SCALE_12_30 };
  TH1* htr_mult_bin[4][3];
  TH1* hda_mult_bin[4][3];
  UnfoldData uinfo_mult_bin[4][3];

  auto run_get_true_njets = [&](int njets) {
    for (int ibin = 0; ibin < 3; ++ibin) {
      TH1 *htr = get_true(
          scale[ibin],
          bins[2*ibin], bins[2*ibin + 1],
          h2tr[njets]
      );
      htr_mult_bin[njets][ibin] = htr;
    }
  };

  auto run_unfold_njets = [&](int njets) {
    for (int ibin = 0; ibin < 3; ++ibin) {
      TH1 *hda = unfold_data(
          scale[ibin], // scale
          bins[2*ibin], bins[2*ibin + 1], // Z-axis range
          h2da[njets], h3mc[njets], // data and migration matrix
          h3mc_sys[njets], // systematics
          &uinfo_mult_bin[njets][ibin] // unfolding-info
      );
      unc_from_mc_to_data(htr_mult_bin[njets][ibin], hda);
      hda_mult_bin[njets][ibin] = hda;
    }
  };

  auto draw_njets = [&](int njets) {
    // Draw true distributions:
    for (int ibin = 0; ibin < 3; ++ibin) {
      TH1 *htr = htr_mult_bin[njets][ibin];
      htr->DrawClone(ibin == 0 ? "hist" : "hist same");
    }
    // Draw data distributions:
    for (int ibin = 0; ibin < 3; ++ibin) {
      TH1 *hda = hda_mult_bin[njets][ibin];
      hda->DrawClone("same e1");
    }
    // Labels:
    TLatex* latex = new TLatex();
    latex->SetTextFont(22);
    latex->SetTextSize(0.05);
    latex->DrawLatexNDC(0.25, 0.80, Form("Jets #geq %d", njets + 1));
    // Legend:
    TLegend *legend = new TLegend(0.42, 0.13, 0.62, 0.39);
    legend->SetTextSize(0.04);
    legend->SetHeader("NC DIS ZEUS 330 pb^{-1}");
    legend->AddEntry(hda_mult_bin[njets][2], Form("%s (x 100)", VAR.legnames[2]), "p");
    legend->AddEntry(hda_mult_bin[njets][1], Form("%s (x 10)", VAR.legnames[1]), "p");
    legend->AddEntry(hda_mult_bin[njets][0], Form("%s", VAR.legnames[0]), "p");
    legend->AddEntry(htr_mult_bin[njets][0], "Ariadne 4.12", "l");
    legend->SetBorderSize(0);
    legend->Draw();
  };


  // ========================================================================
  //                               One+ jet
  Int_t njets = 0;
  cre->cd(njets+1);
  cre->cd(njets+1)->SetLogy();
  cre->SetLogy();
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  run_get_true_njets(njets);
  run_unfold_njets(njets);
  draw_njets(njets);


  // ========================================================================
  //                               Two+ jet
  njets = 1;
  cre->cd(njets+1);
  cre->cd(njets+1)->SetLogy();
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  run_get_true_njets(njets);
  run_unfold_njets(njets);
  draw_njets(njets);


  // ========================================================================
  //                               Three+ jet
  njets = 2;
  cre->cd(njets+1);
  cre->cd(njets+1)->SetLogy();
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  run_get_true_njets(njets);
  run_unfold_njets(njets);
  draw_njets(njets);


  // ========================================================================
  //                               Four+ jet
  njets = 3;
  cre->cd(njets+1);
  cre->cd(njets+1)->SetLogy();
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  run_get_true_njets(njets);
  run_unfold_njets(njets);
  draw_njets(njets);

  TH1 *hunc_mult_bin[4][3];

  auto run_get_unc_njets = [&](int njets) {
    for (int ibin = 0; ibin < 3; ++ibin) {
      TH1 *hda = hda_mult_bin[njets][ibin];
      TH1 *hunc = get_unc_hist(Form("hunc_mult_bin[%d][%d]", njets, ibin), hda);
      hunc->SetTitle(Form("#scale[1.25]{%s}", VAR.legnames[ibin], njets + 1));
      hunc_mult_bin[njets][ibin] = hunc;
    }
  };

  run_get_unc_njets(0);
  run_get_unc_njets(1);
  run_get_unc_njets(2);
  run_get_unc_njets(3);

  auto draw_unc_njets = [&](int njets) {
    auto pad = gPad;
    pad->Divide(1, 3);
    for (int ibin = 0; ibin < 3; ++ibin) {
      pad->cd(ibin + 1);
      hunc_mult_bin[njets][ibin]->DrawClone("p");
    }
  };

  if (g_plot_unc_hists) {
    auto c_unc = new TCanvas;
    c_unc->Divide(2, 2);
    c_unc->cd(1); draw_unc_njets(0);
    c_unc->cd(2); draw_unc_njets(1);
    c_unc->cd(3); draw_unc_njets(2);
    c_unc->cd(4); draw_unc_njets(3);
  }


  TDatime now1;
  now1.Print();

  if (out_path) {
    TFile *fout = new TFile(out_path, "RECREATE");

    htr_int->SetName("htr_int");
    hda_int->SetName("hda_int");
    hda_int_unfold->SetName("hda_int_unfold");
    hmc_int->SetName("hmc_int");
    unfold_data_int.h2mig->SetName("h2mig_int");
    unfold_data_int.h2cor->SetName("h2cor_int");
    unfold_data_int.h2cov->SetName("h2cov_int");
    h_int_rdamc->SetName("hrdamc_int");
    h_int_rdatr->SetName("hrdatr_int");

    htr_int->Write();
    hda_int->Write();
    hda_int_unfold->Write();
    hmc_int->Write();
    unfold_data_int.h2mig->Write();
    unfold_data_int.h2cor->Write();
    unfold_data_int.h2cov->Write();
    h_int_rdamc->Write();
    h_int_rdatr->Write();

    for (int njets = 0; njets < 4; ++njets) {
      for (int ibin = 0; ibin < 3; ++ibin) {
        htr_mult_bin[njets][ibin]->SetName(Form("htr_mult_bin[%d][%d]", njets, ibin));
        hda_mult_bin[njets][ibin]->SetName(Form("hda_mult_bin[%d][%d]", njets, ibin));
        uinfo_mult_bin[njets][ibin].h2mig->SetName(Form("h2mig_mult_bin[%d][%d]", njets, ibin));
        uinfo_mult_bin[njets][ibin].h2cor->SetName(Form("h2cor_mult_bin[%d][%d]", njets, ibin));
        uinfo_mult_bin[njets][ibin].h2cov->SetName(Form("h2cov_mult_bin[%d][%d]", njets, ibin));
        hunc_mult_bin[njets][ibin]->SetName(Form("hunc_mult_bin[%d][%d]", njets, ibin));

        htr_mult_bin[njets][ibin]->Write();
        hda_mult_bin[njets][ibin]->Write();
        uinfo_mult_bin[njets][ibin].h2mig->Write();
        uinfo_mult_bin[njets][ibin].h2cor->Write();
        uinfo_mult_bin[njets][ibin].h2cov->Write();
        hunc_mult_bin[njets][ibin]->Write();
      }
    }
    fout->Write();
    delete fout;
  }
}

int main(int argc, char *argv[])
{
  std::string dpath, mpath, tpath, opath, cpath;
  bool do_sys = false;
  std::vector<std::string> sys;

  struct option longopts[] = {
    { "help"     , false, NULL,    'h' },
    { "output"   , true , NULL,    'o' },
    { "data"     , true , NULL,    'd' },
    { "mc"       , true , NULL,    'm' },
    { "true"     , true , NULL,    't' },
    { "cfg"      , true , NULL,    'c' },
    { "sys"      , false, NULL,    's' },
    { "no-stat"  , false, NULL,    'S' },
    { "q2"       , false, NULL,    'Q' },
    { "pt"       , false, NULL,    'P' },
    { "fix-tau"  , false, NULL, 0x01FF },
    { "full-dphi", false, NULL, 0x02FF },
    { "cut-dphi" , false, NULL, 0x03FF },
    { "no-info"  , false, NULL, 0x04FF },
    { "unc-hists", false, NULL, 0x05FF },
    { 0, 0, 0, 0 },
  };

  {
    int opt;
    while ((opt = getopt_long(argc, argv, "ho:d:m:t:c:s", longopts, NULL)) > 0) {
      switch (opt) {
        case 'h':
          std::cout
            << "usage: " << argv[0] << " [OPTIONS]" << std::endl
            << std::endl
            << "OPTIONS:" << std::endl
            << "  --help    -h                           Show this message and exit." << std::endl
            << "  --output  -o  <path>                   Save normalized xsections to root-file." << std::endl
            << "  --data    -d  <path>                   Specify path to data." << std::endl
            << "  --mc      -m  <path>                   Specify path to MC (detector level)." << std::endl
            << "  --true    -t  <path>                   Specify path to MC-true (hadron level)." << std::endl
            << "  --cfg     -c  <path>                   Load data, MC, MC-true, and sys-files (all optional) from specified config-file." << std::endl
            << "  --sys     -s  [<path-1> <path-2> ...]  Enable estimation of systematics and specify MC files to be used." << std::endl
            << "  --no-stat -S                           Remove statistical uncertainties." << std::endl
            << "  --q2      -Q                           Use Q2 binning (enabled by default)." << std::endl
            << "  --pt      -P                           Use Pt binning." << std::endl
            << "  --fix-tau                              Perform secondary unfolding with fixed regularization parameter." << std::endl
            << "  --full-dphi                            Show whole Dphi-range on X-axis." << std::endl
            << "  --cut-dphi                             Move all bins below Pi/2 into underflow bin." << std::endl
            << "  --no-info                              Disable Info-level messages from ROOT." << std::endl
            << "  --unc-hists                            Plot resulting total uncertainties." << std::endl
            ;
          exit(EXIT_SUCCESS);
          break;

        case 'o':
          opath = optarg;
          break;

        case 'd':
          dpath = optarg;
          break;

        case 'm':
          mpath = optarg;
          break;

        case 't':
          tpath = optarg;
          break;

        case 'c':
          cpath = optarg;
          break;

        case 's':
          do_sys = true;
          break;

        case 'S':
          g_show_stat_unc = false;
          break;

        case 'P':
          g_var_idx = IDX_PT;
          break;

        case 'Q':
          g_var_idx = IDX_Q2;
          break;

        case 0x01FF:
          g_fix_tau = true;
          break;

        case 0x02FF:
          g_limit_x = false;
          break;

        case 0x03FF:
          g_cut_low_dphi = true;
          break;

        case 0x04FF:
          gErrorIgnoreLevel = kWarning;
          break;

        case 0x05FF:
          g_plot_unc_hists = true;
          break;

        default:
          assert(!"undefine command line option");
      }
    }
  }

  // load config-file
  if (!cpath.empty()) {
    libconfig::Config cfg;
    cfg.readFile(cpath.c_str());

    if (dpath.empty() && cfg.lookupValue("data", dpath))
      std::cout << "load data path: \"" << dpath << '"' << std::endl;

    if (mpath.empty() && cfg.lookupValue("MC", mpath) || cfg.lookupValue("MC-det", mpath))
      std::cout << "load MC path: \"" << mpath << '"' << std::endl;

    if (tpath.empty() && (cfg.lookupValue("true", tpath) || cfg.lookupValue("MC-had", tpath)))
      std::cout << "load MC-had path: \"" << tpath << '"' << std::endl;

    if (do_sys) {
      const char *sysstr =
        cfg.exists("sys") ? "sys" :
        cfg.exists("systematic") ? "systematic" :
        cfg.exists("systematics") ? "systematics" :
        NULL;
      if (sysstr) {
        for (auto &spath : cfg.lookup(sysstr))
          sys.emplace_back((const char*)spath);
      }
    }
  }

  if (do_sys) {
    if (sys.empty() && optind == argc) {
      std::cerr << "Error: no data for systematics specified" << std::endl;
      exit(EXIT_FAILURE);
    }
    while (optind < argc)
      sys.emplace_back(argv[optind++]);
  } else {
    if (optind != argc) {
      std::cerr << "Error: unexpected extra arguments are given" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  TApplication app("CrossSection", &argc, argv);

  if (dpath.empty() || mpath.empty() || tpath.empty()) {
    std::cerr << "Error: missing some input" << std::endl;
    exit(EXIT_FAILURE);
  }

  CrossSection(dpath.c_str(), mpath.c_str(), tpath.c_str(),
      opath.empty() ? NULL : opath.c_str(), sys);

  if (opath.empty())
    app.Run();
}

