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


double traffic_light(double queue, double density, double waiting_T){
	double queue_short = GaussianMF(queue, 0.0, 25.0);
	double queue_medium = GaussianMF(queue, 50.0, 30.0);
	double queue_long = GaussianMF(queue, 100.0, 35.0);
	double density_low = TrapezoidalMF(density, -1, 0, 2, 4);
	double density_medium = TrapezoidalMF(density, 3, 4, 6, 7);
	double density_high = TrapezoidalMF(density, 6, 8, 10, 11);
	double waiting_T_short = TrapezoidalMF(waiting_T, -5, 0, 10, 20);
	double waiting_T_medium = TrapezoidalMF(waiting_T, 15, 25, 35, 45);
	double waiting_T_long = TrapezoidalMF(waiting_T, 40, 50, 60, 65);

	double rule1 = fmax(queue_long, density_low);
	double rule2 = fmin(queue_short, density_high);
	double rule3 = fmin(queue_medium, fmin(density_medium, waiting_T_medium));
	double rule4 = fmax(waiting_T_long, queue_short);

	double r1_out = 2.0 * queue + 0.5 * density + 1.5 * waiting_T + 15.0;
	double r2_out = 10;
	double r3_out = 30;
	double r4_out = 10;

	double numerator = (rule1 * r1_out) + (rule2 * r2_out) + (rule3 * r3_out) + (rule4 * r4_out);
	double denominator = rule1 + rule2 + rule3 + rule4;

	return numerator / denominator;

}