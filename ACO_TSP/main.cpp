#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>

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

constexpr int MAX_ANT_N = 1000;
constexpr int MAXN = 100;
constexpr double alpha = 1;
constexpr double beta = 2;
constexpr double Q = 1;

std::vector<std::vector<double>> d;
std::vector<std::vector<double>> phero;
std::vector<std::vector<double>> dphero;
long long vis[MAX_ANT_N];

std::random_device rd; 
std::mt19937 mt(rd());

int choose_city(vector<pair<int, double>>& prob)
{
	if (prob.size <= 0)
	{
		fprintf(stderr, "error in choose_city\nprob.size <= 0\n");
		exit(0);
	}
	if (prob.size == 1)
		return prob[0].first;
	
	for (int i = 1; i < prob.size(); ++i)
		prob[i].second += prob[i-1].second;

	std::uniform_real_distribution<double> unid(0.0, prob.back->second);
	double p = unid(mt);
	for (int i = 0; i < prob.size(); ++i)
	{
		if (p <= prob[i].second)
			return prob[i].first;
	}
	fprintf(stderr, "error in choose_city\nchoose nothing\n");
	exit(0);
}

void generateSol(int n, int ant_n, ans_t& best_sol)
{
	std::vector<std::pair<int, double>> prob;
	ant_t ant_sol; // std::pair<double, cities_t>
	for (int k = 0; k < ant_n; ++k)
	{
		ant_sol.second.clear();
		ant_sol.second.emplace_back(1);
		vis[k] = 0b10;
		prob.clear();
		while (ant_sol[k].second.size < n)
		{
			double tmp = 0;
			int from = *ant_sol.second.back();
			/*
			for (int j = 2; j <= n; ++j)
			{
				if ((1<<j) & vis[k]) continue;
				tmp += pow(ph[i][j], alpha) + pow(1/d[i][j], beta);
				prob.emplace_back(j, 0);
			}
			*/

			for (int j = 2; j <= n; ++j)
			{
				int& j = to.first;
				if ((1<<j) & vis[k]) continue;
				//to.second = (pow(ph[i][j], alpha) + pow(1/d[i][j], beta)) / tmp;
				prob.emplace_back(j, (pow(ph[from][j], alpha) + pow(1/d[from][j], beta)));
			}
			ant_sol.second.emplace_back(choose_city(prob));
			ant_sol.first += dist(i, *ant_sol.second.back());
			dphero[from][*ant_sol.second.back()] += Q/ant_sol.first;
		}
		best_sol = std::min(ant_sol, best_sol);
	}
}
void pheroUpdate(int n, double p = 0.5)
{
	for (int i = 1; i <= n; ++i)
	{
		for (int j = 1; j <= n; ++i)
			phero[i][j] = p * phero[i][j] + dhero[i][j];
	}
}

ans_t ACO(cities_t& cities, int n = 30, int t = 1000, int ant_n = MAX_ANT_N)
{
	ans_t best_sol{std::numeric_limits<double>::infinity(), cities_t{}};
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < t; ++j)
		{
			generateSol(cities.size(), ant_n, best_sol);
			pheroUpdate(cities.size());
		}
	}
	return best_sol;
}

void init(int n)
{
	phero.resize(n+1);
	dphero.resize(n+1);
	d.resize(n+1);
	for (int i = 0; i <= n; ++i)
	{
		phero[i].resize(n+1);
		std::fill(phero[i].begin(), phero[i].end(), 1);
		dphero[i].resize(n+1);
		d[i].resize(n+1);
	}
}

int main()
{
	int n, x, y;
 	cities_t cities;
	while(std::cin >> n >> x >> y)
		cities.emplace_back(n, x, y);

	init(n);
	
	std::sort(cities.begin(), cities.end());
	for (int i = 0; i < cities.size(); ++i)
	{
		for (int j = 0; j < cities.size(); ++j)
		{
			city_t& a = cities[i];
			city_t& b = cities[j];
			d[i+1][j+1] = sqrt( (std::get<1>(a) - std::get<1>(b))*(std::get<1>(a) - std::get<1>(b)) + (std::get<2>(a) - std::get<2>(b))*(std::get<2>(a) - std::get<2>(b)))
		}
	}
	ans_t ans = ACO(cities, n, t);

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
