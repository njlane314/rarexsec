#include "rarexsec/Processor.h"
#include "rarexsec/proc/Selection.h"
#include "rarexsec/proc/Volume.h"

#include <ROOT/RVec.hxx>
#include <cmath>
#include <cstdlib>
#include <string>
#include <cstdint>
#include <functional>
#include <limits>

namespace {
constexpr double kRecognisedPurityMin = 0.5;
constexpr double kRecognisedCompletenessMin = 0.1;

constexpr std::uint64_t kTrainingSeed = 12345u;
constexpr double kTrainingFracSignal = 0.05;
constexpr double kTrainingFracBackground = 0.01;

std::uint64_t make_event_key(int run, int sub, int evt)
{
    const std::uint64_t r = static_cast<std::uint64_t>(static_cast<std::uint32_t>(run));
    const std::uint64_t s = static_cast<std::uint64_t>(static_cast<std::uint32_t>(sub));
    const std::uint64_t e = static_cast<std::uint64_t>(static_cast<std::uint32_t>(evt));
    return (r << 42) ^ (s << 21) ^ e;
}

double stable_uniform(int run, int sub, int evt)
{
    std::uint64_t key = make_event_key(run, sub, evt) ^ kTrainingSeed;
    std::uint64_t h = std::hash<std::uint64_t>{}(key);
    const double inv = 1.0 / static_cast<double>(std::numeric_limits<std::uint64_t>::max());
    return (static_cast<double>(h) + 0.5) * inv;
}
}

//____________________________________________________________________________
ROOT::RDF::RNode rarexsec::Processor::run(ROOT::RDF::RNode node,
                                          const rarexsec::Entry& rec) const
{
    const bool is_data = (rec.source == Source::Data);
    const bool is_ext = (rec.source == Source::Ext);
    const bool is_mc = (rec.source == Source::MC);

    const double scale_mc = (is_mc && rec.pot_nom > 0.0 && rec.pot_eqv > 0.0) ? (rec.pot_nom / rec.pot_eqv) : 1.0;
    const double scale_ext = (is_ext && rec.trig_nom > 0.0 && rec.trig_eqv > 0.0) ? (rec.trig_nom / rec.trig_eqv) : 1.0;

    node = node.Define("w_base", [is_mc, is_ext, scale_mc, scale_ext] {
        const double scale = is_mc ? scale_mc : (is_ext ? scale_ext : 1.0);
        return static_cast<float>(scale);
    });

    if (is_mc) {
        node = node.Define(
            "w_nominal",
            [](float w, float w_spline, float w_tune) {
                const float out = w * w_spline * w_tune;
                if (!std::isfinite(out))
                    return 0.0f;
                if (out < 0.0f)
                    return 0.0f;
                return out;
            },
            {"w_base", "weightSpline", "weightTune"});
    } else {
        node = node.Define("w_nominal", [](float w) { return w; }, {"w_base"});
    }

    if (is_mc) {
        node = node.Define(
            "in_fiducial",
            [](float x, float y, float z) {
                return rarexsec::fiducial::is_in_truth_volume(x, y, z);
            },
            {"nu_vtx_x", "nu_vtx_y", "nu_vtx_z"});

        node = node.Define(
            "count_strange",
            [](int kplus, int kminus, int kzero, int lambda0, int sigplus, int sigzero, int sigminus) {
                return kplus + kminus + kzero + lambda0 + sigplus + sigzero + sigminus;
            },
            {"n_K_plus", "n_K_minus", "n_K0", "n_lambda", "n_sigma_plus", "n_sigma0", "n_sigma_minus"});

        node = node.Define(
            "is_strange",
            [](int strange) { return strange > 0; },
            {"count_strange"});

        node = node.Define(
            "scattering_mode",
            [](int mode) {
                switch (mode) {
                case 0:
                    return 0;
                case 1:
                    return 1;
                case 2:
                    return 2;
                case 3:
                    return 3;
                case 10:
                    return 10;
                default:
                    return -1;
                }
            },
            {"int_mode"});

        node = node.Define(
            "analysis_channels",
            [](bool fv, int nu, int ccnc, int s, int np, int npim, int npip, int npi0, int ngamma) {
                const int npi = npim + npip;
                if (!fv) {
                    if (nu == 0)
                        return static_cast<int>(Channel::OutFV);
                    return static_cast<int>(Channel::External);
                }
                if (ccnc == 1)
                    return static_cast<int>(Channel::NC);
                if (ccnc == 0 && s > 0) {
                    if (s == 1)
                        return static_cast<int>(Channel::CCS1);
                    return static_cast<int>(Channel::CCSgt1);
                }
                if (std::abs(nu) == 12 && ccnc == 0)
                    return static_cast<int>(Channel::ECCC);
                if (std::abs(nu) == 14 && ccnc == 0) {
                    if (npi == 0 && np > 0)
                        return static_cast<int>(Channel::MuCC0pi_ge1p);
                    if (npi == 1 && npi0 == 0)
                        return static_cast<int>(Channel::MuCC1pi);
                    if (npi0 > 0 || ngamma >= 2)
                        return static_cast<int>(Channel::MuCCPi0OrGamma);
                    if (npi > 1)
                        return static_cast<int>(Channel::MuCCNpi);
                    return static_cast<int>(Channel::MuCCOther);
                }
                return static_cast<int>(Channel::Unknown);
            },
            {"in_fiducial", "nu_pdg", "int_ccnc", "count_strange",
             "n_p", "n_pi_minus", "n_pi_plus", "n_pi0", "n_gamma"});

        node = node.Define(
            "is_signal",
            [](bool is_nu_mu_cc, const ROOT::RVec<int>& lambda_decay_in_fid) {
                if (!is_nu_mu_cc) return false;
                for (auto v : lambda_decay_in_fid)
                    if (v) return true;
                return false;
            },
            {"is_nu_mu_cc", "lambda_decay_in_fid"});

        node = node.Define(
            "recognised_signal",
            [](bool is_sig, float purity, float completeness) {
                return is_sig && purity > static_cast<float>(kRecognisedPurityMin) &&
                       completeness > static_cast<float>(kRecognisedCompletenessMin);
            },
            {"is_signal", "neutrino_purity_from_pfp", "neutrino_completeness_from_pfp"});
    } else {
        const int nonmc_channel =
            is_ext ? static_cast<int>(Channel::External) : (is_data ? static_cast<int>(Channel::DataInclusive) : static_cast<int>(Channel::Unknown));

        node = node.Define("in_fiducial", [] { return false; });
        node = node.Define("is_strange", [] { return false; });
        node = node.Define("scattering_mode", [] { return -1; });
        node = node.Define("analysis_channels", [nonmc_channel] { return nonmc_channel; });
        node = node.Define("is_signal", [] { return false; });
        node = node.Define("recognised_signal", [] { return false; });
    }

    node = node.Define(
        "in_reco_fiducial",
        [](float x, float y, float z) {
            return rarexsec::fiducial::is_in_reco_volume(x, y, z);
        },
        {"reco_neutrino_vertex_sce_x", "reco_neutrino_vertex_sce_y", "reco_neutrino_vertex_sce_z"});

    if (is_mc) {
        node = node.Define(
            "is_training",
            [](int run, int sub, int evt, bool is_signal, float w_nominal) {
                if (!std::isfinite(w_nominal) || w_nominal <= 0.0f)
                    return false;
                const double u = stable_uniform(run, sub, evt);
                if (!(u > 0.0 && u < 1.0))
                    return false;
                const double p = is_signal ? kTrainingFracSignal : kTrainingFracBackground;
                return u < p;
            },
            {"run", "sub", "evt", "is_signal", "w_nominal"});

        node = node.Define(
            "w_analysis",
            [](float w_nominal, bool is_signal, bool is_training) {
                if (!std::isfinite(w_nominal) || w_nominal <= 0.0f)
                    return 0.0f;
                if (is_training)
                    return 0.0f;
                const double p = is_signal ? kTrainingFracSignal : kTrainingFracBackground;
                const double scale = 1.0 / (1.0 - p);
                const double w = static_cast<double>(w_nominal) * scale;
                if (!std::isfinite(w) || w <= 0.0)
                    return 0.0f;
                return static_cast<float>(w);
            },
            {"w_nominal", "is_signal", "is_training"});
    } else {
        node = node.Define("is_training", [] { return false; });
        node = node.Define(
            "w_analysis",
            [](float w_nominal) { return w_nominal; },
            {"w_nominal"});
    }

    return node;
}
//____________________________________________________________________________

//____________________________________________________________________________
const rarexsec::Processor& rarexsec::processor()
{
    static const Processor ep{};
    return ep;
}
//____________________________________________________________________________
