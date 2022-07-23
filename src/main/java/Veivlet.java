import java.util.concurrent.atomic.AtomicReference;

public class Veivlet {
    WaveletDOG Dog;
    public Veivlet(double[][] image) throws InterruptedException {
        Dog = new WaveletDOG(NormFactor(image));
        Dog.join();
        }

    protected double[][] NormFactor(double[][] pic){
        double min = pic[0][0], max = pic[0][0];
        for (int i = 0; i < pic.length; i++)for (int j = 0; j < pic[i].length; j++){
            if (pic[i][j]<min)min = pic[i][j];
            if (pic[i][j]>max)max = pic[i][j];
        }
        for (int i = 0; i < pic.length; i++) {
            for (int j = 0; j < pic[i].length; j++) {
                double pix = pic[i][j];
                double res = ((pix-min)*254)/(max-min);
                pic[i][j]=  res;
            }
        }
        return pic;
    }

    abstract class Wavelet extends Thread{
        AtomicReference<double[][]> Wavelet = new AtomicReference<>();
        AtomicReference<double[][]> dX = new AtomicReference<>();
        AtomicReference<double[][]> dY = new AtomicReference<>();
        Thread dXThread, dYThread;
        double[][] normalImage;
        int a = 3, Xquantity, Xdecomposition, Yquantity, Ydecomposition;
        int[] kY,mX,nX, kX, mY, nY;
        public Wavelet(double[][] normalImage){
            this.normalImage = normalImage;
            Xquantity = normalImage.length;
            Yquantity = normalImage[0].length;
            Xdecomposition = (int) ((Math.log(Xquantity)/Math.log(2))-1);
            Ydecomposition = (int) ((Math.log(Yquantity)/Math.log(2))-1);
            nX = new int[Xquantity];for (int i = 0; i < nX.length; i++)nX[i] = i;
            mX = new int[Xdecomposition+1];for (int i = 0; i < mX.length; i++)mX[i] = i;
            kX = new int[Xquantity];for (int i = 0; i < kX.length; i++)kX[i] = i;
            mY = new int[Ydecomposition+1];for (int i = 0; i<mY.length;i++)mY[i] = i;
            nY = new int[Yquantity];for (int i = 0; i < nY.length; i++)nY[i] = i;
            kY = new int[Yquantity];for (int i = 0; i < kY.length; i++)kY[i] = i;
            this.start();
        }
        abstract double WaveletF(double x);
        abstract double WaveletFP1(double x);
        protected double diskretWavelet(int x, double m, int n){return Math.pow(a,-m/2)*WaveletF(Math.pow(a,-m)*x-n);}
        protected double diskretWaveletFP1(int x, double m, int n){return Math.pow(a,-m/2)*WaveletFP1(Math.pow(a,-m)*x-n);}
        protected double[][][] DWTx(double[][] pic){ ;
            double[][][] DWTx = new double[kY.length][mX.length][nX.length];
            for (int y : kY) {
                double[][] DWT = new double[mX.length][nX.length];
                for (int m : mX) {
                    for (int n : nX) {
                        for (int x = 0; x < Xquantity-1; x++){
                            DWT[m][n]+= diskretWavelet(x,Math.pow(2,m-1),n)*(pic[x][y]);
                        }
                    }
                }
                DWTx[y] = DWT;
            }
            return DWTx;
        }
        protected double[][] dX(double[][] pic){ ;
            double[][][] DWTWAVEX = DWTx(pic);
            for (int y : kY) {
                for (int x : kX) {
                    int pix = 0;
                    for (int i = 0; i < Xdecomposition; i++) {
                        for (int j = 0; j < Xquantity-1; j++) {
                            pix +=diskretWaveletFP1(x,Math.pow(2,i-1),j)*DWTWAVEX[y][i][j];
                        }
                    }
                    pic[x][y] = Math.abs(pix);
                }
            }
            return pic;
        }
        protected double[][][] DWTy(double[][] pic){
            double[][][]DWTMHY = new double[kX.length][mY.length][nY.length];
            for (int x : kX) {
                double[][] DWT = new double[mY.length][nY.length];
                for (int m : mY) {
                    for (int n : nY) {
                        for (int y = 0; y < Yquantity-1; y++){
                            DWT[m][n] += diskretWavelet(y,Math.pow(2,m-1),n)*(pic[x][y]);
                        }
                    }
                }
                DWTMHY[x] = DWT;
            }
            return DWTMHY;
        }
        protected double[][] dY(double[][] pic){
            double[][][] DWTWAVEY = DWTy(pic);
            for (int x : kX) {
                for (int y : kY) {
                    int summ = 0;
                    for (int i = 0; i < Ydecomposition; i++) {
                        for (int j = 0; j < Yquantity-1; j++) {
                            summ+=(diskretWaveletFP1(y,Math.pow(2,i-1),j)*DWTWAVEY[x][i][j]);
                        }
                    }
                    pic[x][y]  = Math.abs(summ);
                }
            }
            return pic;
        }
        protected double[][] grab(double[][] DifferentX, double[][] DifferentY){
            double[][] res = new double[DifferentX.length][DifferentX[0].length];
            for (int x = 0; x <DifferentX.length-1 ; x++) {
                for (int y = 0; y < DifferentY.length-1; y++) {
                    double pix1 = DifferentX[x][y];
                    double pix2 = DifferentY[x][y];
                    double resInt = Math.sqrt(Math.pow(pix1, 2) + Math.pow(pix2, 2));
                    res[x][y] = resInt;
                }
            }
            return res;
        }
        protected double[][] NormFactor(double[][] pic){
            double min = pic[0][0], max = pic[0][0];
            for (int i = 0; i < pic.length; i++)for (int j = 0; j < pic[i].length; j++){
                if (pic[i][j]<min)min = pic[i][j];
                if (pic[i][j]>max)max = pic[i][j];
            }
            for (int i = 0; i < pic.length; i++) {
                for (int j = 0; j < pic[i].length; j++) {
                    double pix = pic[i][j];
                    double res = ((pix-min)*254)/(max-min);
                    pic[i][j]= (int) res;
                }
            }
            return pic;
        }
        public void run(){
            dXThread = new Thread(()->{
                dX.set(dX(normalImage));
            });
            dYThread = new Thread(()->{
                dY.set(dY(normalImage));
            });
            try {
                dXThread.start();
                dYThread.start();
                dXThread.join();
                dYThread.join();
                Wavelet.set(NormFactor(grab(dX.get(), dY.get())));
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }

    class WaveletDOG extends Wavelet{

        public WaveletDOG(double[][] normalImage) {
            super(normalImage);
        }

        @Override
        protected double WaveletF(double x) {
            return  (Math.pow(Math.E,-Math.pow(x,2)/2) - 0.5*Math.pow(Math.E,-Math.pow(x,2)/8));
        }

        @Override
        protected double WaveletFP1(double x) {
            return (0.125 * x * Math.pow(Math.E,-Math.pow(x,2)/8) - x*Math.pow(Math.E,-Math.pow(x,2)/2));
        }

    }
}