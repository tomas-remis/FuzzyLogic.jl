#include <stdio.h>
#include <math.h>

double gauss(double x, double mean, double sigma) {
    return exp(-0.5 * pow((x - mean) / sigma, 2));
}

double trap(double x, double a, double b, double c, double d) {
    if (x <= a || x >= d) return 0.0;
    if (x >= b && x <= c) return 1.0;
    if (x > a && x < b) return (x - a) / (b - a);
    if (x > c && x < d) return (d - x) / (d - c);
    return 0.0;
}

double tipper(double service, double food) {
    // Fuzzification
    double service_poor = gauss(service, 0, 1.5);
    double service_good = gauss(service, 5, 1.5);
    double service_excellent = gauss(service, 10, 1.5);

    double food_rancid = trap(food, -2, 0, 1, 3);
    double food_delicious = trap(food, 7, 9, 10, 12);

    // Rule evaluation
    double rule1 = fmax(service_poor, food_rancid);
    double rule2 = service_good;
    double rule3 = fmax(service_excellent, food_delicious);

    // Sugeno outputs
    double r1_out = 5.002; // Cheap
    double r2_out = 15;    // Average
    double r3_out = 2.0 * service + 0.5 * food + 5.0; // Generous

    // Weighted average calculation
    double numerator = (rule1 * r1_out) + (rule2 * r2_out) + (rule3 * r3_out);
    double denominator = rule1 + rule2 + rule3;

    return numerator / denominator;
}

int main() {
    double service, food;
    printf("Enter the number of %s:", "service");// 做字符串插入
    scanf("%lf", &service);
    printf("Enter the number of %s:", "food");// 做字符串插入
    scanf("%lf", &food);

    double tip = tipper(service, food);
    printf("Suggested tip: %lf\n", tip);

    return 0;
}
