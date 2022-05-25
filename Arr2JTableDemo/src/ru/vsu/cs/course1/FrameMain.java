package ru.vsu.cs.course1;

import com.intellij.uiDesigner.core.GridConstraints;
import com.intellij.uiDesigner.core.GridLayoutManager;
import com.intellij.uiDesigner.core.Spacer;
import ru.vsu.cs.util.ArrayUtils;
import ru.vsu.cs.util.JTableUtils;
import ru.vsu.cs.util.SwingUtils;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.Arrays;
import java.util.Scanner;


public class FrameMain extends JFrame {
    private JPanel panelMain;
    private JTable tableInput;
    private JButton buttonLoadInputFromFile;
    private JButton buttonRandomInput;
    private JButton buttonSaveInputInfoFile;
    private JButton buttonReverseRows;
    private JButton buttonSaveOutputIntoFile;
    private JTable tableOutput;
    private JButton buttonReverseColumns;
    private JButton checkButton;
    private JTable outputA;
    private JTable outputB;
    private JTextField матрицаОбратнаяКМатрицеTextField;
    private JTextField матрицаСвободныхЧленовTextField;
    private JButton check1Button;
    private JButton check2Button;
    private JButton buttonEngine;
    private JTextField Values;
    private JTextField Vectors;

    private JFileChooser fileChooserOpen;
    private JFileChooser fileChooserSave;
    private JMenuBar menuBarMain;
    private JMenu menuLookAndFeel;


    public FrameMain() {
        this.setTitle("FrameMain");
        this.setContentPane(panelMain);
        this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        this.pack();

        JTableUtils.initJTableForArray(tableInput, 40, true, true, true, true);
        JTableUtils.initJTableForArray(tableOutput, 40, true, true, true, true);
        //tableOutput.setEnabled(false);
        tableInput.setRowHeight(25);
        tableOutput.setRowHeight(25);

        fileChooserOpen = new JFileChooser();
        fileChooserSave = new JFileChooser();
        fileChooserOpen.setCurrentDirectory(new File("."));
        fileChooserSave.setCurrentDirectory(new File("."));
        FileFilter filter = new FileNameExtensionFilter("Text files", "txt");
        fileChooserOpen.addChoosableFileFilter(filter);
        fileChooserSave.addChoosableFileFilter(filter);

        fileChooserSave.setAcceptAllFileFilterUsed(false);
        fileChooserSave.setDialogType(JFileChooser.SAVE_DIALOG);
        fileChooserSave.setApproveButtonText("Save");

        menuBarMain = new JMenuBar();
        setJMenuBar(menuBarMain);

        menuLookAndFeel = new JMenu();
        menuLookAndFeel.setText("Вид");
        menuBarMain.add(menuLookAndFeel);
        SwingUtils.initLookAndFeelMenu(menuLookAndFeel);

        JTableUtils.writeArrayToJTable(tableInput, new int[][]{
                {1, -2, 1, 1},
                {2, -1, 1, 2},
                {3, 2, 2, -2}
        });

        this.pack();


        buttonLoadInputFromFile.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                try {
                    if (fileChooserOpen.showOpenDialog(panelMain) == JFileChooser.APPROVE_OPTION) {
                        int[][] arr = ArrayUtils.readIntArray2FromFile(fileChooserOpen.getSelectedFile().getPath());
                        JTableUtils.writeArrayToJTable(tableInput, arr);
                    }
                } catch (Exception e) {
                    SwingUtils.showErrorMessageBox(e);
                }
            }
        });
        buttonRandomInput.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                try {
                    int[][] matrix = ArrayUtils.createRandomIntMatrix(
                            tableInput.getRowCount(), tableInput.getColumnCount(), 100);
                    JTableUtils.writeArrayToJTable(tableInput, matrix);
                } catch (Exception e) {
                    SwingUtils.showErrorMessageBox(e);
                }
            }
        });
        buttonSaveInputInfoFile.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                try {
                    if (fileChooserSave.showSaveDialog(panelMain) == JFileChooser.APPROVE_OPTION) {
                        int[][] matrix = JTableUtils.readIntMatrixFromJTable(tableInput);
                        String file = fileChooserSave.getSelectedFile().getPath();
                        if (!file.toLowerCase().endsWith(".txt")) {
                            file += ".txt";
                        }
                        ArrayUtils.writeArrayToFile(file, matrix);
                    }
                } catch (Exception e) {
                    SwingUtils.showErrorMessageBox(e);
                }
            }
        });
        buttonSaveOutputIntoFile.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                try {
                    if (fileChooserSave.showSaveDialog(panelMain) == JFileChooser.APPROVE_OPTION) {
                        int[][] matrix = JTableUtils.readIntMatrixFromJTable(tableOutput);
                        String file = fileChooserSave.getSelectedFile().getPath();
                        if (!file.toLowerCase().endsWith(".txt")) {
                            file += ".txt";
                        }
                        ArrayUtils.writeArrayToFile(file, matrix);
                    }
                } catch (Exception e) {
                    SwingUtils.showErrorMessageBox(e);
                }
            }
        });
//        buttonReverseRows.addActionListener(new ActionListener() {
//            @Override
//            public void actionPerformed(ActionEvent actionEvent) {
//                try {
//                    int[][] matrix = JTableUtils.readIntMatrixFromJTable(tableInput);
//                    Task.reverseRows(matrix);
//                    JTableUtils.writeArrayToJTable(tableOutput, matrix);
//                } catch (Exception e) {
//                    SwingUtils.showErrorMessageBox(e);
//                }
//            }
//        });
        check1Button.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                try {
                    double[][] matrix = JTableUtils.readDoubleMatrixFromJTable(tableInput);
                    double[][] a = Task.defineA(matrix);
                    double[][] b = Task.defineB(matrix);
                    double[][] dets = new double[a.length][1];

                    double det = Task.determinant(a, a.length);
                    //Проверка определителя на 0
                    if (det != 0) {
                        System.out.println("det = " + det + "\n");
                        if (a.length > 1) {
                            dets = Kramer.determs(a, b);
                        }
                    } else {
                        System.out.println("det = 0, => there is no inverse matrix" + "\n");
                    }

                    double[][] rez = Kramer.answersKramer(dets, det);

                    //JTableUtils.writeArrayToJTable(outputA, dets);
                    JTableUtils.writeArrayToJTable(outputB, b);
                    JTableUtils.writeArrayToJTable(outputA, rez);
                } catch (Exception e) {
                    SwingUtils.showErrorMessageBox(e);
                }
            }
        });
        check2Button.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                try {
                    double[][] matrix = JTableUtils.readDoubleMatrixFromJTable(tableInput);
                    double[][] a = Task.defineA(matrix);
                    double[][] b = Task.defineB(matrix);

                    double[][] rez = Gauss.gauss(matrix);

                    //JTableUtils.writeArrayToJTable(outputA, dets);
                    JTableUtils.writeArrayToJTable(outputB, b);
                    JTableUtils.writeArrayToJTable(outputA, rez);
                } catch (Exception e) {
                    SwingUtils.showErrorMessageBox(e);
                }
            }
        });
        checkButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                try {
                    double[][] matrix = JTableUtils.readDoubleMatrixFromJTable(tableInput);
                    double[][] a = Task.defineA(matrix);
                    double[][] b = Task.defineB(matrix);

                    double det = Task.determinant(a, a.length);
                    //Проверка определителя на 0
                    if (det != 0) {
                        System.out.println("det = " + det + "\n");
                        if (a.length > 1) {
                            Task.algDop(a, a.length);
                            Task.transpose(a, a.length);
                            Task.inverceMatrix(a, a.length, det);
                        }
                    } else {
                        System.out.println("det = 0, => there is no inverse matrix" + "\n");
                    }

                    double[][] rez = Task.rezCount(a, b);

                    //JTableUtils.writeArrayToJTable(outputA, a);
                    JTableUtils.writeArrayToJTable(outputB, b);
                    JTableUtils.writeArrayToJTable(outputA, rez);
                } catch (Exception e) {
                    SwingUtils.showErrorMessageBox(e);
                }
            }
        });

        buttonEngine.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                try {
                    double[][] matrix = JTableUtils.readDoubleMatrixFromJTable(tableInput);

                    int n = matrix.length;


                    double[][] V, H, R, EV;
                    double[] d, e, ort;
                    double det;

                    Scanner scanner = new Scanner(System.in);



                    d = new double[n];
                    e = new double[n];
                    ort = new double[n];
                    V = new double[n][n];
                    H = new double[n][n];
                    R = new double[n][n];


                    //print matrix
                    System.out.println("Matrix: ");
                    for (int i = 0; i < matrix.length; i++) {
                        for (int j = 0; j < matrix[0].length; j++) {
                            System.out.print(matrix[i][j] + " ");
                        }
                        System.out.print("\n");
                    }
                    System.out.print("\n");

                    //Prints rref of matrix
                    R = EigenvalueCalculator.rref(matrix);
                    System.out.println("Row reduced matrix: ");
                    for (int i = 0; i < R.length; i++) {
                        for (int j = 0; j < R[0].length; j++) {
                            System.out.print(R[i][j] + " ");
                        }
                        System.out.println("");
                    }
                    System.out.print("\n");

                    // calculate determinant
                    det = EigenvalueCalculator.findDeterminant(matrix);
                    System.out.println("Determinant: " + det);

                    //calculate eigenvalues
                    if (EigenvalueCalculator.isSymmetric(matrix)) {
                        for (int i = 0; i < n; i++) {
                            for (int j = 0; j < n; j++) {
                                V[i][j] = matrix[i][j];
                            }
                        }

                        // Tridiagonalize.
                        EigenvalueCalculator.tred2(n, d, e, V);

                        // Diagonalize.
                        EigenvalueCalculator.tql2(n, d, e, V);
                    } else {
                        H = new double[n][n];
                        ort = new double[n];

                        for (int j = 0; j < n; j++) {
                            for (int i = 0; i < n; i++) {
                                H[i][j] = matrix[i][j];
                            }
                        }

                        // Reduce to Hessenberg form.
                        EigenvalueCalculator.orthes(n, d, e, V, H, ort);

                        // Reduce Hessenberg to real Schur form.
                        EigenvalueCalculator.hqr2(n, d, e, V, H, ort);
                    }


                    //round eigen values
                    for (int i = 0; i < d.length; i++) {
                        d[i] = (double) Math.round((d[i] * 1000000)) / 1000000;
                    }

                    String values = "Real eigenvalues: [";
                    //print eigenvalues
                    System.out.print("Real eigenvalues: [");
                    Arrays.sort(d);
                    for (int i = 0; i < d.length; i++) {
                        if (i != 0) {
                            values+=", ";
                            System.out.print(", ");
                        }
                        values+=d[i];
                        System.out.print(d[i]);
                    }
                    values+="]";
                    System.out.println("]");
                    System.out.print("Complex eigenvalues: [");
                    Arrays.sort(e);
                    for (int i = 0; i < e.length; i++) {
                        if (i != 0) {
                            System.out.print(", ");
                        }
                        System.out.print(e[i]);
                    }
                    System.out.println("]\n");


                    String vectors = "Eigenvectors: ";






                    for(int k = 0; k < d.length; k++) {
                        String ans ="Для л = " + d[k] + "    ";
                        double[][] matrix1 = JTableUtils.readDoubleMatrixFromJTable(tableInput);
                        double[][] matrix2 = new double[n][n + 1];
                        for (int i = 0; i < n; i++) {
                            matrix1[i][i] = matrix1[i][i] - d[k];
                        }

                        for (int i = 0; i < n; i++) {
                            for(int j = 0; j < n; j++) {
                                matrix2[i][j] = matrix1[i][j];
                            }
                        }

                        double[][] rez = Gauss.gauss(matrix2);

                        for(int j = 0; j < n; j++) {
                            ans += (rez[j][0] + " ");
                        }
                        vectors +=(ans + "        ");
                        //JTableUtils.writeArrayToJTable(outputA, rez);
                    }



                    Vectors.setText(vectors);
                    Values.setText(values);









                    //JTableUtils.writeArrayToJTable(outputA, a);
                   // JTableUtils.writeArrayToJTable(outputB, b);
                   // JTableUtils.writeArrayToJTable(outputA, rez);
                } catch (Exception e) {
                    SwingUtils.showErrorMessageBox(e);
                }
            }
        });
    }



    {
// GUI initializer generated by IntelliJ IDEA GUI Designer
// >>> IMPORTANT!! <<<
// DO NOT EDIT OR ADD ANY CODE HERE!
        $$$setupUI$$$();
    }

    /**
     * Method generated by IntelliJ IDEA GUI Designer
     * >>> IMPORTANT!! <<<
     * DO NOT edit this method OR call it in your code!
     *
     * @noinspection ALL
     */
    private void $$$setupUI$$$() {
        panelMain = new JPanel();
        panelMain.setLayout(new GridLayoutManager(5, 2, new Insets(10, 10, 10, 10), 10, 10));
        final JScrollPane scrollPane1 = new JScrollPane();
        scrollPane1.setVerticalScrollBarPolicy(21);
        panelMain.add(scrollPane1, new GridConstraints(0, 0, 1, 2, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_WANT_GROW, null, new Dimension(-1, 200), null, 0, false));
        tableInput = new JTable();
        scrollPane1.setViewportView(tableInput);
        final JPanel panel1 = new JPanel();
        panel1.setLayout(new GridLayoutManager(1, 4, new Insets(0, 0, 0, 0), -1, -1));
        panelMain.add(panel1, new GridConstraints(1, 0, 1, 2, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        buttonLoadInputFromFile = new JButton();
        buttonLoadInputFromFile.setText("Загрузить из файла");
        panel1.add(buttonLoadInputFromFile, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        buttonRandomInput = new JButton();
        buttonRandomInput.setText("Заполнить случайными числами");
        panel1.add(buttonRandomInput, new GridConstraints(0, 1, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        buttonSaveInputInfoFile = new JButton();
        buttonSaveInputInfoFile.setText("Сохранить в файл");
        panel1.add(buttonSaveInputInfoFile, new GridConstraints(0, 2, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        final Spacer spacer1 = new Spacer();
        panel1.add(spacer1, new GridConstraints(0, 3, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, 1, null, new Dimension(100, -1), null, 0, false));
        final JScrollPane scrollPane2 = new JScrollPane();
        scrollPane2.setVerticalScrollBarPolicy(21);
        panelMain.add(scrollPane2, new GridConstraints(3, 0, 1, 2, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_WANT_GROW, null, new Dimension(-1, 200), null, 0, false));
        tableOutput = new JTable();
        scrollPane2.setViewportView(tableOutput);
        final JPanel panel2 = new JPanel();
        panel2.setLayout(new GridLayoutManager(1, 2, new Insets(0, 0, 0, 0), -1, -1));
        panelMain.add(panel2, new GridConstraints(2, 0, 1, 2, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        checkButton = new JButton();
        checkButton.setText("Выполнить");
        panel2.add(checkButton, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        final Spacer spacer2 = new Spacer();
        panel2.add(spacer2, new GridConstraints(0, 1, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, 1, null, null, null, 0, false));
        final JPanel panel3 = new JPanel();
        panel3.setLayout(new GridLayoutManager(1, 2, new Insets(0, 0, 0, 0), -1, -1));
        panelMain.add(panel3, new GridConstraints(4, 0, 1, 2, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        buttonSaveOutputIntoFile = new JButton();
        buttonSaveOutputIntoFile.setText("Сохранить в файл");
        panel3.add(buttonSaveOutputIntoFile, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        final Spacer spacer3 = new Spacer();
        panel3.add(spacer3, new GridConstraints(0, 1, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, 1, null, null, null, 0, false));
    }

    /**
     * @noinspection ALL
     */
    public JComponent $$$getRootComponent$$$() {
        return panelMain;
    }

    private void createUIComponents() {
        // TODO: place custom component creation code here
    }
}
