#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <cmath>

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

ans_t GE(cities_t& cities)
{
	cities_t ans_ct;
	double ans = 0;
	auto n = cities.size();
	auto dist = [](const city_t &a, const city_t &b) { int64 x1 = std::get<1>(a), x2 = std::get<1>(b), y1 = std::get<2>(a), y2 = std::get<2>(b); return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)); };
	ans_ct.emplace_back(cities[0]);
	while (ans_ct.size() < n)
	{
		std::sort(cities.begin() + ans_ct.size(), cities.end(), [&](const city_t& c1, const city_t& c2) { 
				int64 x0 = std::get<1>(*ans_ct.rbegin()), x1 = std::get<1>(c1), x2 = std::get<1>(c2), y0 = std::get<2>(*ans_ct.rbegin()), y1 = std::get<2>(c1), y2 = std::get<2>(c2); 
				return (x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1) < (x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2);
			});

		//std::cerr << "the nearest city to " << std::get<0>(*ans_ct.rbegin()) << " is " << std::get<0>(*(cities.begin() + ans_ct.size())) << '\n';
		//std::cerr << "the distance is " << dist(*ans_ct.rbegin(), *(cities.begin() + ans_ct.size())) << '\n';
		ans += dist(*ans_ct.rbegin(), *(cities.begin() + ans_ct.size()));
		ans_ct.emplace_back(*(cities.begin() + ans_ct.size()));
	}
	ans += dist(*ans_ct.rbegin(), *ans_ct.begin());

	return std::make_pair(ans, ans_ct);
}

int main()
{
	int n, x, y;
 	cities_t cities;
	while(std::cin >> n >> x >> y)
		cities.emplace_back(n, x, y);

	std::sort(cities.begin(), cities.end());
	ans_t ans = GE(cities);
	for (auto& cy : ans.second)
		std::cout << std::get<0>(cy) << " " << std::get<1>(cy) << " " << std::get<2>(cy) << '\n';
	if (ans.second.size())
		std::cout << std::get<0>(ans.second[0]) << " " << std::get<1>(ans.second[0]) << " " << std::get<2>(ans.second[0]) << '\n';

	std::cout << '\n' << "the possible minimal length is " << ans.first << '\n';

	return 0;
}
