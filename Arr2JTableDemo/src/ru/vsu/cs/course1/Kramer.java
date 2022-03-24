package ru.vsu.cs.course1;

public class Kramer extends Task{
    public static double[][] changeCol(double[][] a, double[][] b, int numCol){
        for (int i = 0; i < a.length; i++){
            a[i][numCol] = b[i][0];
        }
        return a;
    }

    public static double[][] copy(double[][] a, double[][] a1){
        for (int r = 0; r < a.length; r++) {
            for (int c = 0; c < a.length; c++){
                a1[r][c] = a[r][c];
            }
        }
        return a1;
    }

    //возвращает массив детерминантов для каждой колонки
    public static double[][] determs(double[][] a, double[][] b) {
        double[][] dets = new double[a.length][1];
        double[][] a1 = new double[a.length][a.length];
        copy(a, a1);
        for (int i = 0; i < a1.length; i++) {
            for (int j = 0; j < a1.length; j++){
                a1[j][i] = b[j][0];
            }
            dets[i][0] = determinant(a1, a1.length);
            copy(a, a1);
        }
        return dets;
    }
    //считает массив ответов
    public static double[][] answersKramer(double[][] determinants, double det){
        double[][] answers = new double[determinants.length][1];
        for (int i = 0; i < determinants.length; i++){
            answers[i][0] = determinants[i][0] / det;
        }
        return answers;
    }
}
