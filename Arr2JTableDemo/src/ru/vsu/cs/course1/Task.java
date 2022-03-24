package ru.vsu.cs.course1;


import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.Scanner;

import static java.lang.System.in;
import static java.lang.System.out;

public class Task {
//3 таск: написать программу, решающую слу методом обратной матрицы. ОБЯЗАТЕЛЬНО с графическим интерфейсом
public static Scanner sc = new Scanner(System.in);
    public void input(double [][]array, int size) {
        for (int r = 0; r < size; r++) {
            for (int c = 0; c < size; c++) {
                System.out.println(("\n Enter a[%d][%d]: "+ r + c));
                array[r][c] = sc.nextDouble();
            }
        }
    }

    static void output(double[][] array, int size) {
        for (int r = 0; r < size; r++) {
            for (int c = 0; c < size; c++) {
                System.out.print(array[r][c]);
                System.out.print(" ");
            }
            System.out.println();;
        }
    }

    public static int width(double[][] matr){
        int width = 0;
        for (int i = 0; i < 1; i++) {
            for (int j = 0; j < matr[i].length; j++) {
                width = j + 1;
            }
        }
        return width;
    }

    //разбиение исх. матрицы на матрицу А и В
    public static double[][] defineA(double[][] matr){
        int width = width(matr);
        double[][] a = new double[matr.length][width - 1];
        for (int i = 0; i < matr.length; i++) {
            for (int j = 0; j < width - 1; j++) {
                a[i][j] = matr[i][j];
            }
        }
        return a;
    }

    public static double[][] defineB(double[][] matr){
        double[][] b = new double[matr.length][1];
        for (int i = 0; i < matr.length; i++) {
            b[i][0] = matr[i][matr.length];
        }
        return b;
    }

    //умножение обр. матр к А на В
    public static double[][] rezCount(double[][] a, double[][] b){
        double[][] answer = new double[a.length][1];
        for (int i = 0; i < a.length; i++) {
            answer[i][0] = 0;
            for (int k = 0; k < a.length; k++)
                answer[i][0] += (a[i][k] * b[k][0]);
        }
        return answer;
    }

    // Получение матрицы без i-й строки и j-го столбца
    static void getMatr(double[][] arr, double[][] p, int i, int j, int m) {
        int ki, kj, di, dj;
        di = 0;
        for (ki = 0; ki < m - 1; ki++) { // проверка индекса строки
            if (ki == i) di = 1;
            dj = 0;
            for (kj = 0; kj < m - 1; kj++) { // проверка индекса столбца
                if (kj == j) dj = 1;
                p[ki][kj] = arr[ki + di][kj + dj];
            }
        }
    }

    // Рекурсивное вычисление определителя
    public static double determinant(double[][] arr, int m) {
        double d, k;
        int i, j, n;
        double [][] p = new double[m][m];
        for (i = 0; i<m; i++)
            p[i] = new double [m];
        j = 0; d = 0;
        k = 1; //(-1) в степени i
        n = m - 1;
        if (m < 1) System.out.println("Определитель вычислить невозможно!");
        if (m == 1) {
            return arr[0][0];
        }
        if (m == 2) {
            return arr[0][0] * arr[1][1] - (arr[1][0] * arr[0][1]);
        }
        if (m > 2) {
            for (i = 0; i < m; i++) {
                getMatr(arr, p, i, 0, m);
                d = d + k * arr[i][0] * determinant(p, n);
                k = -k;
            }
        }
        return d;
    }

    public static void algDop(double[][] arr, int m) {
        int i, j, d, k, n;
        double [][] p = new double[m][m];
        double [][] z = new double[m][m];
        for (i = 0; i < m; i++) {
            p[i] = new double [m];
            z[i] = new double [m];
        }
        if (m >= 2) {
            if (m == 2) {
                for (int r = 0; r < m; r++) {
                    for (int c = 0; c < m; c++) {
                        z[r][c] = Math.pow(-1, r + c) * arr[m - r - 1][m - c - 1];
                    }
                }
            }
            if (m > 2) {
                for (int r = 0; r < m; r++) {
                    for (int c = 0; c < m; c++) {
                        getMatr(arr, p, r, c, m);
                        double det = determinant(p, m - 1);
                        z[r][c] = Math.pow(-1, r + c) * det;
                        if (z[r][c] == -0) {
                            z[r][c] =0;
                        }
                    }
                }
            }
            for (int r = 0; r < m; r++) {
                for (int c = 0; c < m; c++) {
                    arr[r][c] = z[r][c];
                }
            }
        }
    }

    public static void transpose(double[][] arr, int size) {
        double t;
        for(int i = 0; i < size; ++i) {
            for(int j = i; j < size; ++j) {
                t = arr[i][j];
                arr[i][j] = arr[j][i];
                arr[j][i] = t;
            }
        }
    }

    public static void inverceMatrix(double[][] arr, int size, double det) {
        for (int r = 0; r < size; r++) {
            for (int c = 0; c < size; c++) {
                arr[r][c] = arr[r][c] / det;
            }
        }
    }


    public static void check(double[][] arr1, double[][] arr2, int[][] answer, int size){
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                answer[i][j] = 0;
                for (int k = 0; k < size; k++)
                    answer[i][j] += Math.round(arr1[i][k] * arr2[k][j]);
                System.out.print(answer[i][j] + " ");
            }
            System.out.println();
        }
    }
}
