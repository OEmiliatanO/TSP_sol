#include <iostream>
#include <cstdlib>
#include <fstream>
#include <tuple>
#include <vector>
#include <cmath>

using namespace std;

using city = tuple<int,int,int>;

vector<city> cities;

double dist(city& c1, city& c2)
{
	int x1 = get<1>(c1), x2 = get<1>(c2), y1 = get<2>(c1), y2 = get<2>(c2);
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

int main()
{
	fstream ans;
	fstream data;
	ans.open("eil51.ans");
	data.open("eil51.tsp");

	int n, x, y, prev;
	double dis = 0;
	while(data >> n >> x >> y) cities.emplace_back(n, x, y);

	ans >> prev;
	while(ans >> n)
	{
		dis += dist(cities[n - 1], cities[prev - 1]);
		prev = n;
	}
	cout << "the answer is ";
	cout << dis << '\n';

	return 0;
}
