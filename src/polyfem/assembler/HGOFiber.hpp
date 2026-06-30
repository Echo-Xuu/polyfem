#pragma once

#include <polyfem/assembler/GenericFiber.hpp>
#include <polyfem/assembler/GenericElastic.hpp>

namespace polyfem::assembler
{
	class HGOFiber : public GenericFiber<HGOFiber>
	{
	public:
		HGOFiber();

		// sets material params
		void add_multimaterial(const int index, const json &params, const Units &units, const std::string &root_path) override;

		std::string name() const override { return "HGOFiber"; }
		std::map<std::string, ParamFunc> parameters() const override;

		template <typename T>
		T elastic_energy(
			const RowVectorNd &p,
			const double t,
			const int el_id,
			const DefGradMatrix<T> &def_grad) const
		{
			const double k1 = k1_(p, t, el_id);    // fiber stiffness  (manuscript a_f)
			const double k2 = k2_(p, t, el_id);    // exp. stiffening  (manuscript b_f)
			const double kappa = kappa_(p, t, el_id); // dispersion in [0, 1/d]; absent => 0

			// Modified-anisotropy GOH invariant using FULL invariants (manuscript Eq. 4-5):
			//   E4 = kappa * I1 + (1 - d*kappa) * I4 - 1,   d = spatial dimension
			const double d = static_cast<double>(this->size());
			const T i1 = I1(def_grad);
			const T i4 = I4(p, t, el_id, def_grad);
			const T E4 = kappa * i1 + (1.0 - d * kappa) * i4 - 1.0;

			// Smooth distension/compression switch X(E4) (manuscript Eq. 5), replaces
			// the hard tension-only cutoff so the energy is C-infinity. The logistic is
			// centered at E4 = 0, the physical distension/compression boundary (at the
			// unloaded state I1=d, I4=1 => E4=0). The published "-1" inside the switch
			// was a transcription carry-over from the I4-based switch H(I4-1) of the
			// source models (where I4=1 is the unloaded value) and has been confirmed
			// by the authors as an error; the correct, fits-consistent center is 0.
			const T chi = 1.0 / (1.0 + exp(-k_chi_ * E4));

			return (k1 / (2.0 * k2)) * chi * (exp(k2 * E4 * E4) - 1.0);
		}

	private:
		GenericMatParam k1_;     // fiber stiffness   (manuscript a_f)
		GenericMatParam k2_;     // exp. stiffening   (manuscript b_f)
		GenericMatParam kappa_;  // fiber dispersion  kappa  (in [0, 1/d])
		double k_chi_ = 100.0;     // logistic smoothness k_X; manuscript fixes 100
	};
} // namespace polyfem::assembler