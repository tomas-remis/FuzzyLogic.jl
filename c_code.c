double GaussianMF(double x, double mean, double sigma) {
    return exp(-0.5 * pow((x - mean) / sigma, 2));
}


double TrapezoidalMF(double x, double a, double b, double c, double d) {
    if (x <= a || x >= d) return 0.0;
    if (x >= b && x <= c) return 1.0;
    if (x > a && x < b) return (x - a) / (b - a);
    if (x > c && x < d) return (d - x) / (d - c);
    return 0.0;
    }


double tipper(double service, double food){
	double service_poor = GaussianMF(service, 0.0, 1.5);
	double service_good = GaussianMF(service, 5.0, 1.5);
	double service_excellent = GaussianMF(service, 10.0, 1.5);
	double food_rancid = TrapezoidalMF(food, -2, 0, 1, 3);
	double food_delicious = TrapezoidalMF(food, 7, 9, 10, 12);

	double rule1 = fmax(service_poor, food_rancid);
	double rule2 = service_good;
	double rule3 = fmax(service_excellent, food_delicious);

	double r1_out = 5.002;
	double r2_out = 15;
	double r3_out = 2.0 * service + 0.5 * food + 5.0;

	double numerator = (rule1 * r1_out) + (rule2 * r2_out) + (rule3 * r3_out);
	double denominator = rule1 + rule2 + rule3;

	return numerator / denominator;

}