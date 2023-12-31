// https://arxiv.org/pdf/1310.6271.pdf (Based on this).
// This is a brute-force way of generating optimal sorting networks.
// Requires OR-Tools Library and C++20.

#include <string>
#include <vector>
#include <map>
#include <format>
#include <algorithm>
#include <tuple>

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
#include "ortools/util/time_limit.h"

const int64_t n = 4; // Channel Size
const int64_t depth = 4; // Depth 
const int64_t d = depth + 1;

enum obj_type {
	MINIMIZE_DEPTH = 0,
	MINIMIZE_TOTAL_COMPARATORS = 1
};
obj_type obj = (obj_type) 0;

namespace operations_research {
	namespace sat {
		void create_comparator(std::map<std::tuple<int64_t, int64_t, int64_t>, BoolVar>& L,
			CpModelBuilder& cp_model)
		{
			for (int64_t i = 0; i < n; i++) {
				for (int64_t j = 0; j < n; j++) {
					for (int64_t k = 0; k < d; k++) {
						L[{i, j, k}] = cp_model.NewBoolVar().WithName(std::format("L{}{}{}", i, j, k));
					}
				}
			}
			for (int64_t i = 0; i < n; i++) {
				for (int64_t j = 0; j < n; j++) {
					for (int64_t k = 0; k < d; k++) {
						cp_model.AddEquality(L[{i, j, k}], L[{j, i, k}]);
						if (i == j || k == 0)
							cp_model.AddEquality(L[{i, j, k}], 0);
					}
				}
			}
			for (int64_t k = 0; k < d; k++) {
				for (int64_t i = 0; i < n; i++) {
					std::vector<BoolVar>expr;
					for (int64_t j = 0; j < n; j++) {
						expr.push_back(L[{i, j, k}]);
					}
					cp_model.AddLessOrEqual(LinearExpr::Sum(expr), 1);
				}
			}
			std::vector<BoolVar>total_comparators;
			for (int64_t k = 0; k < d; k++) {
				for (int64_t i = 0; i < n; i++) {
					for (int64_t j = i + 1; j < n; j++) {
						total_comparators.push_back(L[{i, j, k}]);
					}
				}
			}
			if (obj == obj_type::MINIMIZE_TOTAL_COMPARATORS) {
				cp_model.Minimize(LinearExpr::Sum(total_comparators));
			}
			std::vector<BoolVar>bd(d);
			for (int64_t k = 0; k < d; k++) {
				bd[k] = cp_model.NewBoolVar().WithName(std::format("bd{}", k));
			}
			for (int64_t k = 1; k < d; k++) {
				std::vector<BoolVar>expr;
				for (int64_t i = 0; i < n; i++) {
					for (int64_t j = i + 1; j < n; j++) {
						expr.push_back(L[{i, j, k}]);
					}
				}
				cp_model.AddGreaterOrEqual(LinearExpr::Sum(expr), 1).OnlyEnforceIf(bd[k]);
				cp_model.AddEquality(LinearExpr::Sum(expr), 0).OnlyEnforceIf(bd[k].Not());
			}
			if (obj == obj_type::MINIMIZE_DEPTH) {
				cp_model.Minimize(LinearExpr::Sum(bd));
			}
			
		}

		void create_value_constraint(std::map<std::tuple<int64_t, int64_t, int64_t>, BoolVar>& L,
			CpModelBuilder& cp_model, const std::vector<int64_t>& a, int64_t id)
		{
			std::map<std::tuple<int64_t, int64_t, int64_t>, BoolVar> V;
			for (int64_t i = 0; i < n; i++) {
				for (int64_t k = 0; k < d; k++) {
					V[{i, k, id}] = cp_model.NewBoolVar().WithName(std::format("V{}{}{}", i, k, id));
				}
			}
			for (int64_t i = 0; i < n; i++) {
				cp_model.AddEquality(V[{i, 0, id}], a[i]);
			}
			std::map<std::tuple<int64_t, int64_t, int64_t>, BoolVar> b;
			for (int64_t k = 0; k < d; k++) {
				for (int64_t i = 0; i < n; i++) {
					std::vector<BoolVar>expr;
					for (int64_t j = 0; j < n; j++) {
						expr.push_back(L[{i, j, k}]);
					}
					auto sum_expr = LinearExpr::Sum(expr);
					b[{i, k, id}] = cp_model.NewBoolVar().WithName(std::format("b{}{}{}", i, k, id));
					cp_model.AddEquality(sum_expr, 0).OnlyEnforceIf(b[{i, k, id}]);
					cp_model.AddGreaterThan(sum_expr, 0).OnlyEnforceIf(b[{i, k, id}].Not());
				}
			}
			for (int64_t k = 1; k < d; k++) {
				for (int64_t i = 0; i < n; i++) {
					cp_model.AddEquality(V[{i, k, id}], V[{i, k - 1, id}]).OnlyEnforceIf(b[{i, k, id}]);
				}
			}
			for (int64_t k = 1; k < d; k++) {
				for (int64_t i = 0; i < n; i++) {
					for (int64_t j = i + 1; j < n; j++) {
						auto bt = cp_model.NewBoolVar().WithName(std::format("bt1{}{}{}{}", k, i, j, id));
						cp_model.AddLessOrEqual(V[{i, k - 1, id}], V[{j, k - 1, id}]).OnlyEnforceIf(bt);
						cp_model.AddGreaterThan(V[{i, k - 1, id}], V[{j, k - 1, id}]).OnlyEnforceIf(bt.Not());
						cp_model.AddEquality(V[{i, k, id}], V[{i, k - 1, id}]).OnlyEnforceIf({ bt, L[{i, j, k}] });
						cp_model.AddEquality(V[{i, k, id}], V[{j, k - 1, id}]).OnlyEnforceIf({ bt.Not(), L[{i, j, k}] });
					}
				}
				for (int64_t j = 1; j < n; j++) {
					for (int64_t i = j - 1; i >= 0; i--) {
						auto bt = cp_model.NewBoolVar().WithName(std::format("bt2{}{}{}{}", k, i, j, id));
						cp_model.AddGreaterOrEqual(V[{i, k - 1, id}], V[{j, k - 1, id}]).OnlyEnforceIf(bt);
						cp_model.AddLessThan(V[{i, k - 1, id}], V[{j, k - 1, id}]).OnlyEnforceIf(bt.Not());
						cp_model.AddEquality(V[{j, k, id}], V[{i, k - 1, id}]).OnlyEnforceIf({ bt, L[{j, i, k}] });
						cp_model.AddEquality(V[{j, k, id}], V[{j, k - 1, id}]).OnlyEnforceIf({ bt.Not(), L[{j, i, k}] });
					}
				}
			}
			for (int64_t i = 0; i < n - 1; i++) {
				cp_model.AddLessOrEqual(V[{i, d - 1, id}], V[{i + 1, d - 1, id}]);
			}
		}

		void validate(const std::vector<std::vector<std::pair<int64_t, int64_t>>>& fl, const std::vector<std::vector<int64_t>>& sequences)
		{
			for (std::vector<int64_t> seq : sequences) {
				for (int i = 0; i < (int)fl.size(); i++) {
					for (const auto& p : fl[i]) {
						seq[p.first] = std::min(seq[p.first], seq[p.second]);
						seq[p.second] = std::max(seq[p.first], seq[p.second]);
					}
				}
				if (!std::is_sorted(std::begin(seq), std::end(seq))) {
					std::cout << "Error" << '\n';
					return;
				}
			}
		}

		std::vector<std::vector<std::pair<int64_t, int64_t>>>
			solve(std::map<std::tuple<int64_t, int64_t, int64_t>, BoolVar>& L, CpModelBuilder& cp_model)
		{
			Model model;
			SatParameters parameters;
			//parameters.set_max_time_in_seconds(10.0);
			parameters.set_log_search_progress(true);
			parameters.set_num_workers(0);
			parameters.set_cp_model_presolve(false);
			//parameters.set_preferred_variable_order((SatParameters_VariableOrder) 1);
			model.Add(NewSatParameters(parameters));
			const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);
			std::vector<std::string> status_vec = { "Unknown", "Model_Invalid", "Feasible", "Infeasible", "Optimal" };
			std::cout << "Status = " << status_vec[(int)response.status()] << '\n';
			if (response.status() == CpSolverStatus::OPTIMAL ||
				response.status() == CpSolverStatus::FEASIBLE)
			{
				std::cout << "Channel Size = " << n << '\n';
				std::cout << "Original Depth = " << depth << '\n';
				std::string obj_str;
				if (obj == obj_type::MINIMIZE_DEPTH) obj_str = "MINIMIZE_DEPTH";
				else if (obj == obj_type::MINIMIZE_TOTAL_COMPARATORS) obj_str = "MINIMIZE_TOTAL_COMPARATORS";
				std::cout << std::format("Obj_Type = {}, Obj_Val = {}", obj_str, response.objective_value()) << '\n';
				std::vector<std::vector<std::pair<int64_t, int64_t>>>fl;
				for (int64_t k = 0; k < d; k++) {
					std::vector<std::pair<int64_t, int64_t>>l;
					for (int64_t i = 0; i < n; i++) {
						for (int64_t j = i + 1; j < n; j++) {
							if (SolutionBooleanValue(response, L[{i, j, k}]))
								l.push_back({ i, j });
						}
					}
					if (!l.empty())
						fl.push_back(l);
				}
				for (int i = 0; i < (int)fl.size(); i++) {
					std::cout << "Depth: " << i + 1 << "," << ' ' << "Index:" << ' ';
					std::cout << "[";
					int cnt = 0;
					for (const auto& p : fl[i]) {
						if (cnt == (int)fl[i].size() - 1) {
							std::cout << std::format("({}, {})", p.first, p.second);
						}
						else {
							std::cout << std::format("({}, {}), ", p.first, p.second);
						}
						++cnt;
					}
					std::cout << "]";
					std::cout << '\n';
				}
				return fl;
			}

		}

		void optimal_sorting_network(const std::vector<std::vector<int64_t>>& sequences)
		{
			CpModelBuilder cp_model;
			std::map<std::tuple<int64_t, int64_t, int64_t>, BoolVar>L;
			create_comparator(L, cp_model);
			int64_t seq_id = 0;
			for (const auto& seq : sequences) {
				create_value_constraint(L, cp_model, seq, seq_id);
				++seq_id;
			}
			//cp_model.ExportToFile("sort_network_13.txt");
			auto fl = solve(L, cp_model);
			if (!fl.empty()) {
				std::cout << "Validation Started" << '\n';
				validate(fl, sequences);
				std::cout << "Validation Finished" << '\n';
			}
		}
	}
}

std::vector<std::vector<int64_t>>generate_binary_seq(const int64_t N)
{
	std::vector<std::vector<int64_t>> sequences;
	for (int64_t i = 0; i < (1ll << N); i++) {
		std::vector<int64_t>seq;
		for (int64_t j = 0; j < N; j++) {
			if (i & (1ll << j))
				seq.push_back(1);
			else
				seq.push_back(0);
		}
		sequences.push_back(seq);
	}
	return sequences;
}


int main()
{
	auto seq = generate_binary_seq(n);
	operations_research::sat::optimal_sorting_network(seq);

}
