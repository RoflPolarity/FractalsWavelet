import org.math.plot.Plot3DPanel;
import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.UnsupportedEncodingException;
import java.nio.charset.StandardCharsets;

public class projetion {
    static JFrame frame = new JFrame();
    public static void getParamWindow() throws UnsupportedEncodingException {
        String title = "Сгенерировать";
        byte[] titleBytes= title.getBytes("windows-1251");
        String newTitle = new String(titleBytes,StandardCharsets.UTF_8);
        String WindowTitle = "Параметры генерации поверхности";
        byte[] WindowTitleBytes = WindowTitle.getBytes("windows-1251");
        String newWindowTitle = new String(WindowTitleBytes, StandardCharsets.UTF_8);


        JLabel SeedLabel = new JLabel("Seed");
        JLabel HLabel = new JLabel("h");
        JFrame ParamFrame = new JFrame(newWindowTitle);
        JTextField seed = new JTextField(10);
        JTextField h = new JTextField(10);
        final double[] seedValue = new double[1];
        final double[] hValue = new double[1];
        Button generate = new Button();
        generate.setLabel(newTitle);
        generate.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                seedValue[0] = Double.parseDouble(seed.getText());
                hValue[0] = Double.parseDouble(h.getText());
                getTerrain(seedValue[0], hValue[0]);
                frame.repaint();
            }
        });
        JPanel contents = new JPanel(new VerticalLayout());
        contents.add(SeedLabel);
        contents.add(seed);
        contents.add(HLabel);
        contents.add(h);
        contents.add(generate);
        ParamFrame.setContentPane(contents);
        ParamFrame.setBackground(Color.black);
        ParamFrame.setResizable(false);
        ParamFrame.setSize(new Dimension(160,180));
        ParamFrame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        ParamFrame.setVisible(true);
    }
    public static void getTerrain(double seed, double h){
        Terrain terrain = new Terrain(seed, h);
        Veivlet veivlet = new Veivlet(terrain.getMap());
        double[][] map = veivlet.Dog.Wavelet.get();
        double[] x = new double[map.length];
        for (int i = 0; i < x.length; i++) {
            x[i] = i;
        }
        double[] y = new double[map[1].length];
        for (int i = 0; i < y.length; i++) {
            y[i] = i;
        }
        Plot3DPanel projction = new Plot3DPanel("SOUTH");
        projction.addGridPlot("Terrain",x,y,map);
        frame.setSize(800,600);
        frame.setContentPane(projction);
        frame.setLocation(250,250);
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }
    public static void main(String[] args) throws UnsupportedEncodingException {
        getParamWindow();

    }
    static class VerticalLayout implements LayoutManager {
        private Dimension size = new Dimension();

        // Следующие два метода не используются
        public void addLayoutComponent   (String name, Component comp) {}
        public void removeLayoutComponent(Component comp) {}

        // Метод определения минимального размера для контейнера
        public Dimension minimumLayoutSize(Container c) {
            return calculateBestSize(c);
        }
        // Метод определения предпочтительного размера для контейнера
        public Dimension preferredLayoutSize(Container c) {
            return calculateBestSize(c);
        }
        // Метод расположения компонентов в контейнере
        public void layoutContainer(Container container)
        {
            // Список компонентов
            Component list[] = container.getComponents();
            int currentY = 5;
            for (int i = 0; i < list.length; i++) {
                // Определение предпочтительного размера компонента
                Dimension pref = list[i].getPreferredSize();
                // Размещение компонента на экране
                list[i].setBounds(5, currentY, pref.width, pref.height);
                // Учитываем промежуток в 5 пикселов
                currentY += 5;
                // Смещаем вертикальную позицию компонента
                currentY += pref.height;
            }
        }
        // Метод вычисления оптимального размера контейнера
        private Dimension calculateBestSize(Container c)
        {
            // Вычисление длины контейнера
            Component[] list = c.getComponents();
            int maxWidth = 0;
            for (int i = 0; i < list.length; i++) {
                int width = list[i].getWidth();
                // Поиск компонента с максимальной длиной
                if ( width > maxWidth )
                    maxWidth = width;
            }
            // Размер контейнера в длину с учетом левого отступа
            size.width = maxWidth + 5;
            // Вычисление высоты контейнера
            int height = 0;
            for (int i = 0; i < list.length; i++) {
                height += 5;
                height += list[i].getHeight();
            }
            size.height = height;
            return size;
        }
    }
}