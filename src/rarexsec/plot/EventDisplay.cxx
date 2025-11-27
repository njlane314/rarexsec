#include "rarexsec/plot/EventDisplay.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>

#include <ROOT/RConfig.h>
#include <TStyle.h>
#include <nlohmann/json.hpp>

#include "rarexsec/plot/Plotter.h"
//____________________________________________________________________________
rarexsec::plot::EventDisplay::EventDisplay(Spec spec, Options opt, DetectorData data)
    : spec_(std::move(spec)), opt_(std::move(opt)), data_(std::move(data)), plot_name_(rarexsec::plot::Plotter::sanitise(spec_.id)), output_directory_(opt_.out_dir) {}
//____________________________________________________________________________
rarexsec::plot::EventDisplay::EventDisplay(Spec spec, Options opt, SemanticData data)
    : spec_(std::move(spec)), opt_(std::move(opt)), data_(std::move(data)), plot_name_(rarexsec::plot::Plotter::sanitise(spec_.id)), output_directory_(opt_.out_dir) {}
//____________________________________________________________________________
void rarexsec::plot::EventDisplay::draw(TCanvas& canvas) 
{
    setup_canvas(canvas);
    build_histogram();

    switch (spec_.mode) {
    case Mode::Detector:
        draw_detector(canvas);
        break;
    case Mode::Semantic:
        draw_semantic(canvas);
        if (opt_.show_legend)
            draw_semantic_legend();
        break;
    }

    canvas.Update();
}
//____________________________________________________________________________
void rarexsec::plot::EventDisplay::draw_and_save(const std::string& image_format) 
{
    std::filesystem::create_directories(output_directory_);
    TCanvas canvas(plot_name_.c_str(), spec_.title.c_str(),
                   opt_.canvas_size, opt_.canvas_size);
    draw(canvas);

    const std::string fmt = image_format.empty() ? "png" : image_format;
    canvas.SaveAs((output_directory_ + "/" + plot_name_ + "." + fmt).c_str());
}
//____________________________________________________________________________
void rarexsec::plot::EventDisplay::draw_and_save(const std::string& image_format,
                                                 const std::string& file_override) 
{
    std::filesystem::create_directories(output_directory_);
    TCanvas canvas(plot_name_.c_str(), spec_.title.c_str(),
                   opt_.canvas_size, opt_.canvas_size);
    draw(canvas);

    if (!file_override.empty()) {
        canvas.SaveAs(file_override.c_str());
    } else {
        const std::string fmt = image_format.empty() ? "png" : image_format;
        canvas.SaveAs((output_directory_ + "/" + plot_name_ + "." + fmt).c_str());
    }
}
//____________________________________________________________________________
void rarexsec::plot::EventDisplay::setup_canvas(TCanvas& c) const 
{
    c.SetCanvasSize(opt_.canvas_size, opt_.canvas_size);
    c.SetBorderMode(0);
    c.SetFrameBorderMode(0);
    c.SetFrameLineColor(0);
    c.SetFrameLineWidth(0);
    c.SetFixedAspectRatio();

    const double m = std::clamp(opt_.margin, 0.02, 0.25);
    c.SetTopMargin(m);
    c.SetBottomMargin(m);
    c.SetLeftMargin(m);
    c.SetRightMargin(m);

    gStyle->SetTitleAlign(23);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleY(1 - m / 3.0);
}
//____________________________________________________________________________
std::pair<int, int> rarexsec::plot::EventDisplay::deduce_grid(int requested_w,
                                                              int requested_h,
                                                              std::size_t flat_size) 
{
    if (requested_w > 0 && requested_h > 0)
        return {requested_w, requested_h};
    if (requested_w > 0 && requested_h == 0) {
        int h = static_cast<int>(flat_size / requested_w);
        return {requested_w, std::max(1, h)};
    }
    if (requested_h > 0 && requested_w == 0) {
        int w = static_cast<int>(flat_size / requested_h);
        return {std::max(1, w), requested_h};
    }
    int dim = static_cast<int>(std::sqrt(static_cast<double>(flat_size)));
    dim = std::max(1, dim);
    return {dim, dim};
}
//____________________________________________________________________________
void rarexsec::plot::EventDisplay::build_histogram() 
{
    const int bin_offset = 1;
    const auto [W, H] = std::visit([&](auto const& vec) {
        return deduce_grid(spec_.grid_w, spec_.grid_h, vec.size());
    },
                                   data_);

    hist_.reset(new TH2F(spec_.id.c_str(), spec_.title.c_str(),
                         W, 0, W,
                         H, 0, H));
    hist_->SetDirectory(nullptr);

    if (std::holds_alternative<DetectorData>(data_)) {
        const auto& v = std::get<DetectorData>(data_);
        const int n = static_cast<int>(v.size());
        for (int r = 0; r < H; ++r) {
            for (int c = 0; c < W; ++c) {
                const int idx = r * W + c;
                if (idx >= n)
                    break;
                float x = v[idx];
                if (opt_.use_log_z && x <= opt_.det_min)
                    x = static_cast<float>(opt_.det_min);
                hist_->SetBinContent(c + bin_offset, r + bin_offset, x);
            }
        }
    } else {
        const auto& v = std::get<SemanticData>(data_);
        const int n = static_cast<int>(v.size());
        for (int r = 0; r < H; ++r) {
            for (int c = 0; c < W; ++c) {
                const int idx = r * W + c;
                if (idx >= n)
                    break;
                hist_->SetBinContent(c + bin_offset, r + bin_offset, v[idx]);
            }
        }
    }
}
//____________________________________________________________________________
void rarexsec::plot::EventDisplay::draw_detector(TCanvas& c) 
{
    c.SetFillColor(kWhite);
    c.SetTicks(0, 0);

    hist_->SetStats(false);
    hist_->SetMinimum(opt_.det_min);
    hist_->SetMaximum(opt_.det_max);

    hist_->GetXaxis()->SetTitle("Local Wire Coordinate");
    hist_->GetYaxis()->SetTitle("Local Drift Coordinate");
    hist_->GetXaxis()->CenterTitle(true);
    hist_->GetYaxis()->CenterTitle(true);

    constexpr double axis_offset = 0.80;
    hist_->GetXaxis()->SetTitleOffset(axis_offset);
    hist_->GetYaxis()->SetTitleOffset(axis_offset);

    hist_->GetXaxis()->SetTickLength(0);
    hist_->GetYaxis()->SetTickLength(0);
    hist_->GetXaxis()->SetLabelSize(0);
    hist_->GetYaxis()->SetLabelSize(0);
    hist_->GetXaxis()->SetAxisColor(0);
    hist_->GetYaxis()->SetAxisColor(0);

    if (opt_.use_log_z)
        c.SetLogz();

    hist_->Draw("COL");
}
//____________________________________________________________________________
void rarexsec::plot::EventDisplay::draw_semantic(TCanvas& c)
{
    constexpr int palette_size = 15;
    const int background = TColor::GetColor(230, 230, 230);

    std::array<int, palette_size> palette = {
        background,
        TColor::GetColor("#666666"),
        TColor::GetColor("#e41a1c"),
        TColor::GetColor("#377eb8"),
        TColor::GetColor("#4daf4a"),
        TColor::GetColor("#ff7f00"),
        TColor::GetColor("#984ea3"),
        TColor::GetColor("#ffff33"),
        TColor::GetColor("#1b9e77"),
        TColor::GetColor("#f781bf"),
        TColor::GetColor("#a65628"),
        TColor::GetColor("#66a61e"),
        TColor::GetColor("#e6ab02"),
        TColor::GetColor("#a6cee3"),
        TColor::GetColor("#b15928")};
    gStyle->SetPalette(palette_size, palette.data());

    c.SetFillColor(kWhite);
    c.SetFrameFillColor(background);
    c.SetTicks(0, 0);

    hist_->SetStats(false);
    hist_->GetZaxis()->SetRangeUser(-0.5, palette_size - 0.5);

    hist_->GetXaxis()->SetTitle("Local Wire Coordinate");
    hist_->GetYaxis()->SetTitle("Local Drift Coordinate");
    hist_->GetXaxis()->CenterTitle(true);
    hist_->GetYaxis()->CenterTitle(true);

    constexpr double axis_offset = 0.80;
    hist_->GetXaxis()->SetTitleOffset(axis_offset);
    hist_->GetYaxis()->SetTitleOffset(axis_offset);

    hist_->GetXaxis()->SetTickLength(0);
    hist_->GetYaxis()->SetTickLength(0);
    hist_->GetXaxis()->SetLabelSize(0);
    hist_->GetYaxis()->SetLabelSize(0);
    hist_->GetXaxis()->SetAxisColor(0);
    hist_->GetYaxis()->SetAxisColor(0);

    hist_->Draw("COL");
}
//____________________________________________________________________________
void rarexsec::plot::EventDisplay::draw_semantic_legend()
{
    constexpr int palette_size = 15;
    const int background = TColor::GetColor(230, 230, 230);

    std::array<int, palette_size> counts{};
    if (std::holds_alternative<SemanticData>(data_)) {
        for (int v : std::get<SemanticData>(data_)) {
            if (v >= 0 && v < palette_size)
                counts[static_cast<std::size_t>(v)]++;
        }
    }

    std::vector<int> order(palette_size - 1);
    std::iota(order.begin(), order.end(), 1);
    std::stable_sort(order.begin(), order.end(),
                     [&](int a, int b) { return counts[a] > counts[b]; });

    legend_.reset(new TLegend(0.12, 0.86, 0.95, 0.975, "", "brNDC"));
    legend_->SetNColumns(std::max(1, opt_.legend_cols));
    legend_->SetFillColor(background);
    legend_->SetFillStyle(1001);
    legend_->SetBorderSize(0);
    legend_->SetTextFont(42);
    legend_->SetTextSize(0.025);

    const std::array<const char*, palette_size> labels = {
        "#emptyset",
        "Cosmic", "#mu", "e^{-}", "#gamma",
        "#pi^{#pm}", "#pi^{0}", "n", "p",
        "K^{#pm}", "K^{0}", "#Lambda", "#Sigma^{#pm}",
        "#Sigma^{0}", "Other"};

    std::array<int, palette_size> palette = {
        background,
        TColor::GetColor("#666666"),
        TColor::GetColor("#e41a1c"),
        TColor::GetColor("#377eb8"),
        TColor::GetColor("#4daf4a"),
        TColor::GetColor("#ff7f00"),
        TColor::GetColor("#984ea3"),
        TColor::GetColor("#ffff33"),
        TColor::GetColor("#1b9e77"),
        TColor::GetColor("#f781bf"),
        TColor::GetColor("#a65628"),
        TColor::GetColor("#66a61e"),
        TColor::GetColor("#e6ab02"),
        TColor::GetColor("#a6cee3"),
        TColor::GetColor("#b15928")};

    legend_entries_.clear();
    for (int idx : order) {
        auto h = std::make_unique<TH1F>((spec_.id + "_leg_" + std::to_string(idx)).c_str(), "", 1, 0, 1);
        if (counts[idx] > 0) {
            h->SetFillColor(palette[idx]);
            h->SetLineColor(palette[idx]);
            h->SetLineWidth(1);
            h->SetFillStyle(1001);
            std::ostringstream lab;
            lab << labels[static_cast<std::size_t>(idx)] << " (" << counts[idx] << ")";
            legend_->AddEntry(h.get(), lab.str().c_str(), "f");
        } else {
            h->SetFillColor(background);
            h->SetLineColor(background);
            h->SetLineWidth(0);
            h->SetFillStyle(1001);
            legend_->AddEntry(h.get(), "", "f");
        }
        legend_entries_.push_back(std::move(h));
    }

    legend_->Draw();
}
//____________________________________________________________________________

static std::string replace_all(std::string s, const std::string& from, const std::string& to) {
    if (from.empty())
        return s;
    std::size_t pos = 0;
    while ((pos = s.find(from, pos)) != std::string::npos) {
        s.replace(pos, from.size(), to);
        pos += to.size();
    }
    return s;
}
//____________________________________________________________________________
static std::string format_tag(std::string pattern, const std::string& plane, int run, int sub, int evt) 
{
    pattern = replace_all(std::move(pattern), "{plane}", plane);
    pattern = replace_all(std::move(pattern), "{run}", std::to_string(run));
    pattern = replace_all(std::move(pattern), "{sub}", std::to_string(sub));
    pattern = replace_all(std::move(pattern), "{evt}", std::to_string(evt));
    return pattern;
}
//____________________________________________________________________________
void rarexsec::plot::EventDisplay::render_from_rdf(ROOT::RDF::RNode df, const BatchOptions& opt) 
{
    std::error_code ec;
    std::filesystem::create_directories(opt.out_dir, ec);
    if (ec) {
        std::cerr << "[EventDisplay] Failed to create output directory '"
                  << opt.out_dir << "': " << ec.message() << '\n';
    }

    auto filtered = df;
    if (!opt.selection_expr.empty())
        filtered = filtered.Filter(opt.selection_expr);

    const auto n_rows =
        static_cast<std::size_t>(filtered.Count().GetValue());
    if (n_rows == 0) {
        std::cerr << "[EventDisplay] No rows matched selection; nothing to render."
                  << '\n';
        return;
    }

    const std::size_t n_to_render =
        (opt.n_events > 0)
            ? std::min<std::size_t>(n_rows,
                                    static_cast<std::size_t>(opt.n_events))
            : n_rows;

    std::clog << "[EventDisplay] Selection matched " << n_rows
              << " rows; rendering up to " << n_to_render << " events."
              << '\n';

    auto limited = filtered.Range(static_cast<ULong64_t>(n_to_render));

    const bool use_combined_pdf = (!opt.combined_pdf.empty() && opt.image_format == "pdf");
    std::filesystem::path combined_path;
    std::size_t total_pages = 0;
    if (use_combined_pdf) {
        combined_path = std::filesystem::path(opt.out_dir) / opt.combined_pdf;
        total_pages   = n_to_render * opt.planes.size();

#if defined(R__HAS_IMPLICITMT)
        if (ROOT::IsImplicitMTEnabled()) {
            ROOT::DisableImplicitMT();
            std::clog << "[EventDisplay] Implicit MT disabled for stable combined PDF output." << '\n';
        }
#else
        if (ROOT::IsImplicitMTEnabled()) {
            ROOT::DisableImplicitMT();
            std::clog << "[EventDisplay] ROOT built without R__HAS_IMPLICITMT; disabling MT for combined PDF." << '\n';
        }
#endif
    }

    using nlohmann::json;
    json manifest = json::array();
    std::mutex manifest_mutex;

    auto display_opts = opt.display;
    display_opts.out_dir = opt.out_dir;

    if (opt.mode == Mode::Detector) {
        const std::vector<std::string> cols{
            opt.cols.run, opt.cols.sub, opt.cols.evt,
            opt.cols.det_u, opt.cols.det_v, opt.cols.det_w};

        std::atomic<std::size_t> page_idx{0};
        limited.Foreach(
            [&](int run, int sub, int evt,
                const std::vector<float>& det_u,
                const std::vector<float>& det_v,
                const std::vector<float>& det_w) {
                std::clog << "[EventDisplay] Rendering detector images for "
                          << "run=" << run
                          << " sub=" << sub
                          << " evt=" << evt
                          << '\n';

                auto pick = [&](const std::string& plane) -> const std::vector<float>& {
                    if (plane == "U")
                        return det_u;
                    else if (plane == "V")
                        return det_v;
                    else
                        return det_w;
                };

                for (const auto& plane : opt.planes) {
                    const auto& img = pick(plane);

                    auto plane_opts = display_opts;
                    if (!img.empty()) {
                        std::vector<float> vals;
                        vals.reserve(img.size());
                        for (float v : img)
                            if (v > 0.0f)
                                vals.push_back(v);

                        if (!vals.empty()) {
                            std::sort(vals.begin(), vals.end());
                            auto q = [&](double f) -> float {
                                std::size_t idx = std::min(vals.size() - 1,
                                                           static_cast<std::size_t>(f * vals.size()));
                                return vals[idx];
                            };

                            const float min_pos = q(0.02);
                            const float max_val = q(0.995);

                            plane_opts.det_min = std::max(min_pos, 1e-4f);
                            plane_opts.det_max = max_val;
                        }
                    }
                    const std::string tag = format_tag(opt.file_pattern, plane, run, sub, evt);
                    const std::string title =
                        "Detector Image, Plane " + plane +
                        " - Run " + std::to_string(run) +
                        ", Subrun " + std::to_string(sub) +
                        ", Event " + std::to_string(evt);

                    rarexsec::plot::EventDisplay::Spec spec{tag, title, Mode::Detector};
                    EventDisplay ed(spec, plane_opts, img);

                    if (use_combined_pdf) {
                        const std::size_t idx = page_idx++;
                        std::string target = combined_path.string();
                        if (idx == 0)
                            target += "(";
                        else if (idx + 1 == total_pages)
                            target += ")";
                        ed.draw_and_save("pdf", target);
                        if (!opt.manifest_path.empty()) {
                            std::lock_guard<std::mutex> lock(manifest_mutex);
                            manifest.push_back({{"run", run}, {"sub", sub}, {"evt", evt}, {"plane", plane}, {"file", combined_path.string()}});
                        }
                    } else {
                        ed.draw_and_save(opt.image_format);
                        if (!opt.manifest_path.empty()) {
                            const std::string file = (std::filesystem::path(opt.out_dir) /
                                                      (rarexsec::plot::Plotter::sanitise(tag) + "." + opt.image_format))
                                                         .string();
                            std::lock_guard<std::mutex> lock(manifest_mutex);
                            manifest.push_back({{"run", run}, {"sub", sub}, {"evt", evt}, {"plane", plane}, {"file", file}});
                        }
                    }
                }
            },
            cols);
    } else {
        const std::vector<std::string> cols{
            opt.cols.run, opt.cols.sub, opt.cols.evt,
            opt.cols.sem_u, opt.cols.sem_v, opt.cols.sem_w};

        std::atomic<std::size_t> page_idx{0};
        limited.Foreach(
            [&](int run, int sub, int evt,
                const std::vector<int>& sem_u,
                const std::vector<int>& sem_v,
                const std::vector<int>& sem_w) {
                std::clog << "[EventDisplay] Rendering semantic images for "
                          << "run=" << run
                          << " sub=" << sub
                          << " evt=" << evt
                          << '\n';

                auto pick = [&](const std::string& plane) -> const std::vector<int>& {
                    if (plane == "U")
                        return sem_u;
                    else if (plane == "V")
                        return sem_v;
                    else
                        return sem_w;
                };

                for (const auto& plane : opt.planes) {
                    const auto& img = pick(plane);
                    const std::string tag = format_tag(opt.file_pattern, plane, run, sub, evt);
                    const std::string title =
                        "Semantic Image, Plane " + plane +
                        " - Run " + std::to_string(run) +
                        ", Subrun " + std::to_string(sub) +
                        ", Event " + std::to_string(evt);

                    rarexsec::plot::EventDisplay::Spec spec{tag, title, Mode::Semantic};
                    EventDisplay ed(spec, display_opts, img);

                    if (use_combined_pdf) {
                        const std::size_t idx = page_idx++;
                        std::string target = combined_path.string();
                        if (idx == 0)
                            target += "(";
                        else if (idx + 1 == total_pages)
                            target += ")";
                        ed.draw_and_save("pdf", target);
                        if (!opt.manifest_path.empty()) {
                            std::lock_guard<std::mutex> lock(manifest_mutex);
                            manifest.push_back({{"run", run}, {"sub", sub}, {"evt", evt}, {"plane", plane}, {"file", combined_path.string()}});
                        }
                    } else {
                        ed.draw_and_save(opt.image_format);
                        if (!opt.manifest_path.empty()) {
                            const std::string file = (std::filesystem::path(opt.out_dir) /
                                                      (rarexsec::plot::Plotter::sanitise(tag) + "." + opt.image_format))
                                                         .string();
                            std::lock_guard<std::mutex> lock(manifest_mutex);
                            manifest.push_back({{"run", run}, {"sub", sub}, {"evt", evt}, {"plane", plane}, {"file", file}});
                        }
                    }
                }
            },
            cols);
    }
    if (!opt.manifest_path.empty()) {
        std::ofstream ofs(opt.manifest_path);
        ofs << manifest.dump(2);
        std::clog << "[EventDisplay] Wrote event display manifest: " << opt.manifest_path << '\n';
    }
}
//____________________________________________________________________________
