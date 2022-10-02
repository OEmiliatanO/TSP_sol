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
};*/

int64 tot = 1;

ans_t BF(cities_t& cities)
{
	int cnt = 0;
	cities_t ans_ct;
	double ans = std::numeric_limits<double>::infinity();
	std::cerr << "total iteration: " << tot << '\n';
	std::cerr << "00.00  \% completed";
	std::cerr.precision(3);

	do
	{
		double dist = 0;
		city_t prev = *cities.begin();
		for (const auto& ct : cities)
		{
			const auto& [n1,x1,y1] = prev;
			const auto& [n2,x2,y2] = ct;
			dist += sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
			prev = ct;
		}
		const auto& [n1, x1, y1] = prev;
		const auto& [n2, x2, y2] = *cities.begin();
		dist += sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
		if (ans > dist)
		{
			ans = dist;
			ans_ct = cities;
		}
		++cnt;
		std::cerr << '\r' << std::fixed << (double)cnt*100/tot;
	} while(next_permutation(cities.begin(), cities.end()));
	return std::make_pair(ans, ans_ct);
}

int main()
{
	int n, x, y;
 	cities_t cities;
	// use pipe to get the data in
	while(std::cin >> n >> x >> y)
	{
		cities.emplace_back(n, x, y);
		tot *= cities.size();
	}

	std::sort(cities.begin(), cities.end());
	ans_t ans = BF(cities);
	
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
