import java.util.Random;
public class Terrain {
    private double SEED;
    private double mnoz;
    private double h = -(1+SEED*mnoz);

    public Terrain(double seed, double h){
        this.SEED = seed;
        this.mnoz = h;
    }

    public double[][] getMap() {
        final int DATA_SIZE = 9;
        double[][] data = new double[DATA_SIZE][DATA_SIZE];
        data[0][0] = data[0][DATA_SIZE - 1] = data[DATA_SIZE - 1][0] = data[DATA_SIZE - 1][DATA_SIZE - 1] = SEED;
        Random r = new Random();
        for (int sideLength = DATA_SIZE - 1; sideLength >= 2; sideLength /= 2, h /= 2.0) {
            int halfSide = sideLength / 2;
            for (int x = 0; x < DATA_SIZE - 1; x += sideLength) {
                for (int y = 0; y < DATA_SIZE - 1; y += sideLength) {
                    double avg = data[x][y] + data[x + sideLength][y] + data[x][y + sideLength] + data[x + sideLength][y + sideLength];
                    avg /= 4.0;
                    data[x + halfSide][y + halfSide] = avg + (r.nextDouble() * 2 * h) - h;
                }
            }
            for (int x = 0; x < DATA_SIZE - 1; x += halfSide) {
                for (int y = (x + halfSide) % sideLength; y < DATA_SIZE - 1; y += sideLength) {
                    double avg = data[(x - halfSide + DATA_SIZE) % DATA_SIZE][y] + data[(x + halfSide) % DATA_SIZE][y] +
                            data[x][(y + halfSide) % DATA_SIZE] +
                            data[x][(y - halfSide + DATA_SIZE) % DATA_SIZE];
                    avg /= 4.0;

                    avg = avg + (r.nextDouble() * 2 * h) - h;
                    data[x][y] = avg;
                    if (x == 0) data[DATA_SIZE - 1][y] = avg;
                    if (y == 0) data[x][DATA_SIZE - 1] = avg;
                }
            }
        }
    return data;
    }
}
