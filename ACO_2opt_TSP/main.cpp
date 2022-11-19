#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <list>
#include <cassert>

using int64 = long long;
using city_t = std::tuple<int, int64, int64>;
using cities_t = std::vector<city_t>;
using ans_t = std::pair<double, cities_t>;

/*
struct city
{
	int n;
	int64 x, y;
	city(int _n, const int64& _x, const int64& _y): n(_n), x(_x), y(_y) {}
	city(int&& _n, int64&& _x, int64&& _y): n(std::move(_n)), x(std::move(_x)), y(std::move(_y)) {}
	city(): n(0), x(0), y(0) {}
};*/

constexpr int MAXN = 64; // city

constexpr int MAX_ANT_N = 50;  // 50
constexpr double alpha = 0.8;  // 0.8 1
constexpr double beta = 4;     // 4   1
constexpr double Q = 100;      // 100
constexpr double P = 0.8;      // 0.8

double d[MAXN + 1][MAXN + 1];
double phero[MAXN + 1][MAXN + 1];
double dphero[MAXN + 1][MAXN + 1];

std::random_device rd; 
std::mt19937 mt(rd());

int choose_city(std::vector<std::pair<int, double>>& prob)
{
	if (prob.size() <= 0)
	{
		std::cerr << "error in choose_city\nprob.size <= 0\n";
		exit(0);
	}
	if (prob.size() == 1)
		return prob[0].first;

	double psum = 0;
	for (size_t i = 0; i < prob.size(); ++i)
		psum += prob[i].second;
	
	//std::cerr << psum << '\n';
	std::uniform_real_distribution<double> unid(0.0, psum);
	double p = unid(mt), choose = 0;
	for (size_t i = 0; i < prob.size(); ++i)
	{
		choose += prob[i].second;
		if (p <= choose)
			return prob[i].first;
	}
	std::cerr << "error in choose_city\nchoose nothing\n";
	exit(0);
}

double dist(city_t& a, city_t& b)
{
	return sqrt( (std::get<1>(a) - std::get<1>(b))*(std::get<1>(a) - std::get<1>(b)) + (std::get<2>(a) - std::get<2>(b))*(std::get<2>(a) - std::get<2>(b)));
}
void generateSol(int n, int ant_n, cities_t& cities, ans_t& best_sol)
{
	std::uniform_int_distribution<int> unid(1, n);
	std::vector<std::pair<int, double>> prob;
	std::pair<double, std::list<int>> ant_sol;
	for (int k = 0; k < ant_n; ++k)
	{
		ant_sol.first = 0;
		ant_sol.second.clear();
		ant_sol.second.emplace_back(unid(mt));
		long long vis = (1LL << ant_sol.second.back());

		while (ant_sol.second.size() < (size_t)n)
		{
			int from = ant_sol.second.back();
			
			prob.clear();

			for (int j = 1; j <= n; ++j)
			{
				if ((1LL<<j) & vis) continue;
				if (from == j) continue;
				prob.emplace_back(j, (pow(phero[from][j], alpha) * pow(1/d[from][j], beta)));
			}
			int to = choose_city(prob);
			ant_sol.second.emplace_back(to);
			ant_sol.first += d[from][to];
			vis |= (1LL<<to);
		}
		ant_sol.first += d[ant_sol.second.back()][ant_sol.second.front()];
		for (auto it = ant_sol.second.begin(); it != ant_sol.second.end(); ++it)
		{
			//std::cout << *it << " ";
		}

		// 2-opt
		std::pair<size_t, size_t> swap_pair{0, 0};
		double opt_len = ant_sol.first;
		for (size_t i = 0; i < ant_sol.second.size(); ++i)
			for (size_t j = i+1; j < ant_sol.second.size(); ++j)
			{
				double new_len = ant_sol.first;
				if (i == 0 && j == ant_sol.second.size() - 1) continue;
				if (i > 0)
				{
					new_len -= d[*std::next(ant_sol.second.begin(), i-1)][*std::next(ant_sol.second.begin(), i)];
					new_len += d[*std::next(ant_sol.second.begin(), i-1)][*std::next(ant_sol.second.begin(), j)];
				}
				else
				{
					assert(i == 0);
					new_len -= d[ant_sol.second.back()][*std::next(ant_sol.second.begin(), i)];
					new_len += d[ant_sol.second.back()][*std::next(ant_sol.second.begin(), j)];
				}

				if (j+1 < ant_sol.second.size())
				{
					new_len -= d[*std::next(ant_sol.second.begin(), j)][*std::next(ant_sol.second.begin(), j+1)];
					new_len += d[*std::next(ant_sol.second.begin(), i)][*std::next(ant_sol.second.begin(), j+1)];
				}
				else
				{
					assert(j == ant_sol.second.size()-1);
					new_len -= d[*std::next(ant_sol.second.begin(), j)][ant_sol.second.front()];
					new_len += d[*std::next(ant_sol.second.begin(), i)][ant_sol.second.front()];
				}
				if (new_len < opt_len)
				{
					swap_pair = std::make_pair(i, j);
					opt_len = new_len;
				}
			}
		ant_sol.first = opt_len;
		
		if (swap_pair.first < swap_pair.second)
			std::reverse(std::next(ant_sol.second.begin(), swap_pair.first), std::next(ant_sol.second.begin(), swap_pair.second+1));
		
		double L = 0;
		for (auto it = ant_sol.second.begin(); std::next(it, 1) != ant_sol.second.end(); ++it)
		{
			L += d[*it][*std::next(it, 1)];
			dphero[*it][*std::next(it, 1)] += Q/L;
		}
		L += d[ant_sol.second.back()][ant_sol.second.front()];
		dphero[ant_sol.second.back()][ant_sol.second.front()] += Q/L;
		
		assert(fabs(L - ant_sol.first) <= 1e-5);

		if (best_sol.first > ant_sol.first)
		{
			best_sol.first = ant_sol.first;
			best_sol.second.resize(n);
			int i = 0;
			for (auto& city : ant_sol.second)
				best_sol.second[i++] = cities[city];
		}
	}
}
void pheroUpdate(int n, double p = P)
{
	for (int i = 1; i <= n; ++i)
	{
		for (int j = 1; j <= n; ++j)
		{
			phero[i][j] = p * phero[i][j] + dphero[i][j];
			dphero[i][j] = 0;
		}
	}
}

void init(int n)
{
	memset(dphero, 0, sizeof(dphero));
	for (int i = 0; i <= n; ++i)
		for (int j = 0; j <= n; ++j)
			phero[i][j] = 1.0;
}

ans_t ACO(cities_t& cities, int n = 30, int t = 1000, int ant_n = MAX_ANT_N)
{
	ans_t best_sol{std::numeric_limits<double>::infinity(), cities_t{}};
	double avg = 0;
	for (int i = 0; i < n; ++i)
	{
		init(cities.size());
		for (int j = 0; j < t; ++j)
		{
			generateSol(cities.size() - 1, ant_n, cities, best_sol);
			pheroUpdate(cities.size() - 1);
			std::cout << '\r' << std::fixed << (double) (j+i*t+1)*100/(t*n) << '%';
		}
		avg += best_sol.first;
	}
	avg /= 30;
	std::cout << '\n';
	std::cout << "avg = " << avg << '\n';
	return best_sol;
}

int main()
{
	size_t n;
	int64 x, y;
 	cities_t cities;
	cities.emplace_back(0, 0, 0);
	while(std::cin >> n >> x >> y)
		cities.emplace_back(n, x, y);
	
	std::sort(cities.begin(), cities.end());

	memset(d, 0, sizeof(d));
	for (size_t i = 1; i <= cities.size()-1; ++i)
	{
		for (size_t j = 1; j <= cities.size()-1; ++j)
		{
			d[i][j] = dist(cities[i], cities[j]);
			//std::cout << d[i][j] << ' ';
		}
		//std::cout << '\n';
	}

	ans_t ans = ACO(cities, 30, 1000);

	std::fstream plottxt;
	plottxt.open("plot.txt", std::ios::out);
	
	std::cout << "the minimal length is " << ans.first << '\n';
	for (auto& cy : ans.second)
	{
		std::cout << std::get<0>(cy) << '\n';
		plottxt << std::get<0>(cy) << " " << std::get<1>(cy) << " " << std::get<2>(cy) << '\n';
	}
	if (ans.second.size())
		plottxt << std::get<0>(ans.second[0]) << " " << std::get<1>(ans.second[0]) << " " << std::get<2>(ans.second[0]) << '\n';
	plottxt.close();
	return 0;
}
