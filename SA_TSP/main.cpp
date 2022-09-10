#include <iostream>
#include <cstring>
#include <algorithm>
#include <random>
#include <vector>
#include <cmath>

using int64 = long long;

std::random_device rd;
std::mt19937 mt(rd());
std::uniform_real_distribution<double> unid(0.0, 1.0);

struct city
{
	int n;
	int64 x, y;
	city(int _n, const int64& _x, const int64& _y): n(_n), x(_x), y(_y) {}
	city(int&& _n, int64&& _x, int64&& _y): n(std::move(_n)), x(std::move(_x)), y(std::move(_y)) {}
	city(): n(0), x(0), y(0) {}
};

double dist(const city& c1, const city& c2)
{
	return sqrt((c1.x - c2.x) * (c1.x - c2.x) + (c1.y - c2.y) * (c1.y - c2.y));
}

double E(const std::vector<city>& vec)
{
	double res = 0;
	for (auto it = vec.begin() + 1; it != vec.end(); ++it)
		res += dist(*it, *(it - 1));
	return res + dist(*vec.rbegin(), *vec.begin());
}

std::vector<city> rdChange(const std::vector<city>& vec)
{
	std::vector<city> res(vec);
	iter_swap(res.begin() + res.size() * unid(mt), res.begin() + res.size() * unid(mt));
	return res;
}

constexpr int MAXIter = (int)1e7;
constexpr double Rt = 0.97;
constexpr double EndT = 1e-5;
constexpr double InitT = 1e3;
std::pair<double, std::vector<city>> SA(const std::vector<city>& cities)
{
	int counter = 0;
	std::vector<city> current(cities);
	std::vector<city> nex;
	double minDis = std::numeric_limits<double>::infinity();
	double T = InitT;
	auto condjmp = [&](const double dt) { return unid(mt) < exp(dt / T); };
	while (counter < MAXIter and T > EndT)
	{
		nex = std::move(rdChange(current));
		double Enex = E(nex);
		double Ecurr = E(current);
		double dt = Enex - Ecurr;
		if (dt < 0 or condjmp(dt))
		{
			current.swap(nex);
			minDis = Enex;
		}
		T *= Rt;
		++counter;
	}
	return make_pair(minDis, current);
}

int main()
{
	int n, x, y;
	std::vector<city> cities;
	while(std::cin >> n >> x >> y)
		cities.emplace_back(n, x, y);

	std::pair<double, std::vector<city>> ans;
	ans = std::move(SA(cities));

	std::cout << "the minimal distance is " << ans.first << '\n';
	for (auto& it : ans.second)
		std::cout << it.n << " ";
	std::cout << (*ans.second.begin()).n << '\n';
	
	return 0;
}
