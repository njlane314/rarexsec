#include "rarexsec/Processor.h"
#include "rarexsec/proc/Selection.h"
#include "rarexsec/proc/Volume.h"

#include <ROOT/RVec.hxx>
#include <cmath>
#include <cstdlib>
#include <string>

namespace {
constexpr double kRecognisedPurityMin = 0.5;
constexpr double kRecognisedCompletenessMin = 0.1;
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
            {"nu_vx", "nu_vy", "nu_vz"});

        node = node.Define(
            "count_strange",
            [](int kplus, int kminus, int kzero, int lambda0, int sigplus, int sigzero, int sigminus) {
                return kplus + kminus + kzero + lambda0 + sigplus + sigzero + sigminus;
            },
            {"count_kaon_plus", "count_kaon_minus", "count_kaon_zero",
             "count_lambda", "count_sigma_plus", "count_sigma_zero", "count_sigma_minus"});

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
            {"interaction_mode"});

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
            {"in_fiducial", "neutrino_pdg", "interaction_ccnc", "count_strange",
             "count_proton", "count_pi_minus", "count_pi_plus", "count_pi_zero", "count_gamma"});

        node = node.Define(
            "is_signal",
            [](int ch) { return ch == static_cast<int>(Channel::CCS1) || ch == static_cast<int>(Channel::CCSgt1); },
            {"analysis_channels"});

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
