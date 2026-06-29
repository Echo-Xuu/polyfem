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
			// the hard tension-only cutoff so the energy is C-infinity.
			// NOTE: the manuscript centers the logistic at E4 = 1 (e4_center default).
			// Physically the distension/compression boundary is E4 = 0; set the JSON
			// key "e4_center": 0.0 to switch there. See accompanying write-up.
			const T chi = 1.0 / (1.0 + exp(-k_chi_ * (E4 - e4_center_)));

			return (k1 / (2.0 * k2)) * chi * (exp(k2 * E4 * E4) - 1.0);
		}

	private:
		GenericMatParam k1_;     // fiber stiffness   (manuscript a_f)
		GenericMatParam k2_;     // exp. stiffening   (manuscript b_f)
		GenericMatParam kappa_;  // fiber dispersion  kappa  (in [0, 1/d])
		double k_chi_ = 100.0;     // logistic smoothness k_X; manuscript fixes 100
		double e4_center_ = 1.0;   // logistic center; manuscript = 1.0, physical boundary = 0.0
	};
} // namespace polyfem::assembler