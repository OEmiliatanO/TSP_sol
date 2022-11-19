#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <cmath>
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
};
*/

int64 tot = 1;
constexpr int maxn = 26;
std::vector<std::vector<double>> DP;
std::pair<int, int> from[(1 << maxn)][maxn];

ans_t DP_solve(cities_t& cities)
{
	cities_t ans_ct;
	double ans = std::numeric_limits<double>::infinity();
	int n = cities.size();
	tot = (1 << n);
	auto dist = [](const city_t &a, const city_t &b) { int64 x1 = std::get<1>(a), x2 = std::get<1>(b), y1 = std::get<2>(a), y2 = std::get<2>(b); return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)); };
	
	for (int i = 0; i < tot; ++i)
		for (int j = 0; j < n; ++j)
			DP[i][j] = std::numeric_limits<double>::infinity();
	DP[1][0] = 0;
	from[1][0] = std::make_pair(-1, -1);

	std::cerr.precision(4);
	for (int i = 1; i < (1 << n); ++i)
	{
		for (int j = 1; j < n; ++j)
		{
			if (i & (1 << j)) continue;
			if (!(i & 1)) continue;
			for (int k = 0; k < n; ++k)
				if (i & (1 << k))
				{
					double distkj = dist(cities[k], cities[j]);
					if (DP[(1 << j) | i][j] > DP[i][k] + distkj)
					{
						from[(1 << j) | i][j] = std::make_pair(i, k);
						DP[(1 << j) | i][j] = DP[i][k] + distkj;
					}
				}
		}
	}
	std::pair<int, int> idx;
	for (int i = 0; i < n; ++i)
	{
		if (ans > DP[(1 << n) - 1][i] + dist(cities[i], cities[0]))
		{
			idx = std::make_pair((1 << n) - 1, i);
			ans = DP[(1 << n) - 1][i] + dist(cities[i], cities[0]);
		}
	}

	while (idx.second != -1)
	{
		ans_ct.insert(ans_ct.begin(), {std::get<0>(cities[idx.second]), std::get<1>(cities[idx.second]), std::get<2>(cities[idx.second])});
		idx = from[idx.first][idx.second];
	}
	return std::make_pair(ans, ans_ct);
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv)
{
	int n, x, y;
 	cities_t cities;
	while(std::cin >> n >> x >> y)
		cities.emplace_back(n, x, y);
	
	n = cities.size();
	if (n >= 26)
	{
		std::cout << "the data is too big." << '\n';
		return 0;
	}
	DP.resize(1 << n);
	for (auto& it : DP)
		it.resize(n);

	std::sort(cities.begin(), cities.end());
	ans_t ans = DP_solve(cities);

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
