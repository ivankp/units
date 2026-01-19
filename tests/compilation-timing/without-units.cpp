double KineticEnergy(double m, double v) {
    return 0.5 * m * (v * v);
}

int main() {
    return KineticEnergy(3, 7) > 10.;
}
