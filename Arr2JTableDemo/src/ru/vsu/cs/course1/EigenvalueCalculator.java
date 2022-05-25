package ru.vsu.cs.course1;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class EigenvalueCalculator {
    public static void main(String[] args) {
        double[][] matrix, V, H, R, EV;
        double[] d, e, ort;
        double det;

        int n;

        Scanner scanner = new Scanner(System.in);

        System.out.println("Enter dimension of the square matrix: ");
        n = scanner.nextInt();

        matrix = new double[n][n];
        d = new double[n];
        e = new double[n];
        ort = new double[n];
        V = new double[n][n];
        H = new double[n][n];
        R = new double[n][n];


        // insert values
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                System.out.printf("Values: " + i + " - " + j);
                System.out.printf("\n");
                matrix[i][j] = scanner.nextDouble();
            }
        }

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
        R = rref(matrix);
        System.out.println("Row reduced matrix: ");
        for (int i = 0; i < R.length; i++) {
            for (int j = 0; j < R[0].length; j++) {
                System.out.print(R[i][j] + " ");
            }
            System.out.println("");
        }
        System.out.print("\n");

        // calculate determinant
        det = findDeterminant(matrix);
        System.out.println("Determinant: " + det);

        //calculate eigenvalues
        if (isSymmetric(matrix)) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    V[i][j] = matrix[i][j];
                }
            }

            // Tridiagonalize.
            tred2(n, d, e, V);

            // Diagonalize.
            tql2(n, d, e, V);
        } else {
            H = new double[n][n];
            ort = new double[n];

            for (int j = 0; j < n; j++) {
                for (int i = 0; i < n; i++) {
                    H[i][j] = matrix[i][j];
                }
            }

            // Reduce to Hessenberg form.
            orthes(n, d, e, V, H, ort);

            // Reduce Hessenberg to real Schur form.
            hqr2(n, d, e, V, H, ort);
        }


        //round eigen values
        for (int i = 0; i < d.length; i++) {
            d[i] = (double) Math.round((d[i] * 1000000)) / 1000000;
        }

        //print eigenvalues
        System.out.print("Real eigenvalues: [");
        Arrays.sort(d);
        for (int i = 0; i < d.length; i++) {
            if (i != 0) {
                System.out.print(", ");
            }
            System.out.print(d[i]);
        }
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

        //print Eigenvectors
        EV = calculateEigenvector(matrix, arraySet(d));
        System.out.print("Eigenvectors: \n");
        for (int i = 0; i < matrix.length; i++) {
            if (countNonZero(EV[i]) == 0) {
                continue;
            }
            System.out.print("[");
            for (int j = 0; j < matrix[0].length; j++) {
                if (j != 0) {
                    System.out.print(", ");
                }
                System.out.print(EV[i][j] + "");
            }
            System.out.println("] - " + d[i]);
        }
        System.out.print("\n");
        System.out.println("Program Finished.");
    }

    ////Methods
    //Find Determinant
    public static double findDeterminant(double[][] matrix) {
        if (matrix.length == 2) {
            return (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
        } else {
            int sum = 0;
            for (int i = 0; i < matrix.length; i++) {
                double[][] cofactor = findCofactor(matrix, 0, i);
                sum += (Math.pow(-1, i) * matrix[0][i] * findDeterminant(cofactor));
            }
            return sum;
        }
    }

    //ищем алгебр. дополнения
    public static double[][] findCofactor(double[][] matrix, int i, int j) {
        double[][] coMatrix = new double[matrix.length - 1][matrix.length - 1];
        int coI = 0;
        int coJ = 0;
        for (int x = 0; x < matrix.length; x++) {
            if (x == i) {
                continue;
            }
            for (int y = 0; y < matrix.length; y++) {
                if (y == j) {
                    continue;
                }
                coMatrix[coI][coJ] = matrix[x][y];
                coJ++;
                if (coJ == coMatrix.length) {
                    coI++;
                    coJ = 0;
                }
            }
        }
        return coMatrix;
    }

    //операции
    public static double[] addMatrixRow(double[] row1, double[] row2, double scale) { //добавление row1 к row2
        for (int x = 0; x < row1.length; x++) {
            row2[x] += row1[x] * scale;
        }
        return row2;
    }

    public static double[] subtractMatrixRow(double[] row1, double[] row2, double scale) {
        return addMatrixRow(row1, row2, (scale * -1));
    }

    public static double[] scaleMatrixRow(double[] a, double scale) {
        for (int i = 0; i < a.length; i++) {
            a[i] /= scale;
        }
        return a;
    }

    public static double[][] replaceMatrixRow(double[][] a, int row1Index, int row2Index) {
        double[] temp = a[row1Index];
        a[row1Index] = a[row2Index];
        a[row2Index] = temp;
        return a;
    }

    //приведение матрицы к форме верхней матрицы Хессенберга
    public static double[][] rref(double[][] matrix) {
        double[][] rref = new double[matrix.length][];
        for (int i = 0; i < matrix.length; i++)
            rref[i] = Arrays.copyOf(matrix[i], matrix[i].length);//копируем матрицы

        int r = 0;
        for (int c = 0; c < rref[0].length && r < rref.length; c++) {
            int j = r;
            for (int i = r + 1; i < rref.length; i++) {
                if (Math.abs(rref[i][c]) > Math.abs(rref[j][c])) {
                    j = i;
                }
            }
            if (Math.abs(rref[j][c]) < 0.00001) {
                continue;
            }

            double[] temp = rref[j];
            rref[j] = rref[r];
            rref[r] = temp;

            double s = 1.0 / rref[r][c];
            for (j = 0; j < rref[0].length; j++) {
                rref[r][j] *= s;
            }
            for (int i = 0; i < rref.length; i++) {
                if (i != r) {
                    double t = rref[i][c];
                    for (j = 0; j < rref[0].length; j++) {
                        rref[i][j] -= t * rref[r][j];
                    }
                }
            }
            r++;
        }

        return rref;
    }

    public static double[][] transpose(double[][] a) {
        int m = a.length;
        int n = a[0].length;
        double[][] b = new double[n][m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                b[j][i] = a[i][j];
        return b;
    }


    public static boolean isSymmetric(double[][] a) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                if (i == j) {
                    break;
                }
                if (a[i][j] != a[j][i]) {
                    return false;
                }
            }
        }
        return true;
    }

    //
    // симметричное преобразование Хаусхолдера к тридиагональному виду
    //   на каждом шаге первая колонка  преобразуется через малоразмерное
    //   преобразование подобия Хаусхолдера к первой колонке или , где это полином степени , который определяет стратегию «сдвигов»
    public static void tred2(int n, double[] d, double[] e, double[][] V) {


    // Это получено из процедур Algol tred2 с помощью справочника

        for (int j = 0; j < n; j++) {
            d[j] = V[n - 1][j];
        }

        // Householder reduction to tridiagonal form.

        for (int i = n - 1; i > 0; i--) {

            // Масштабирование, чтобы избежать переполнения.

            double scale = 0.0;
            double h = 0.0;
            for (int k = 0; k < i; k++) {
                scale = scale + Math.abs(d[k]);
            }
            if (scale == 0.0) {
                e[i] = d[i - 1];
                for (int j = 0; j < i; j++) {
                    d[j] = V[i - 1][j];
                    V[i][j] = 0.0;
                    V[j][i] = 0.0;
                }
            } else {

                // плучение вектора хаусхолдера.

                for (int k = 0; k < i; k++) {
                    d[k] /= scale;
                    h += d[k] * d[k];
                }
                double f = d[i - 1];
                double g = Math.sqrt(h);
                if (f > 0) {
                    g = -g;
                }
                e[i] = scale * g;
                h = h - f * g;
                d[i - 1] = f - g;
                for (int j = 0; j < i; j++) {
                    e[j] = 0.0;
                }

                // Примените преобразование подобия к оставшимся столбцам.

                for (int j = 0; j < i; j++) {
                    f = d[j];
                    V[j][i] = f;
                    g = e[j] + V[j][j] * f;
                    for (int k = j + 1; k <= i - 1; k++) {
                        g += V[k][j] * d[k];
                        e[k] += V[k][j] * f;
                    }
                    e[j] = g;
                }
                f = 0.0;
                for (int j = 0; j < i; j++) {
                    e[j] /= h;
                    f += e[j] * d[j];
                }
                double hh = f / (h + h);
                for (int j = 0; j < i; j++) {
                    e[j] -= hh * d[j];
                }
                for (int j = 0; j < i; j++) {
                    f = d[j];
                    g = e[j];
                    for (int k = j; k <= i - 1; k++) {
                        V[k][j] -= (f * e[k] + g * d[k]);
                    }
                    d[j] = V[i - 1][j];
                    V[i][j] = 0.0;
                }
            }
            d[i] = h;
        }

        // Накапливание преобразований

        for (int i = 0; i < n - 1; i++) {
            V[n - 1][i] = V[i][i];
            V[i][i] = 1.0;
            double h = d[i + 1];
            if (h != 0.0) {
                for (int k = 0; k <= i; k++) {
                    d[k] = V[k][i + 1] / h;
                }
                for (int j = 0; j <= i; j++) {
                    double g = 0.0;
                    for (int k = 0; k <= i; k++) {
                        g += V[k][i + 1] * V[k][j];
                    }
                    for (int k = 0; k <= i; k++) {
                        V[k][j] -= g * d[k];
                    }
                }
            }
            for (int k = 0; k <= i; k++) {
                V[k][i + 1] = 0.0;
            }
        }
        for (int j = 0; j < n; j++) {
            d[j] = V[n - 1][j];
            V[n - 1][j] = 0.0;
        }
        V[n - 1][n - 1] = 1.0;
        e[0] = 0.0;
    }

    // симметричный тридиагональный QL алгоритм.

    public static void tql2(int n, double[] d, double[] e, double[][] V) {

        // Это получено из процедур Algol tred2 с помощью справочника

        for (int i = 1; i < n; i++) {
            e[i - 1] = e[i];
        }
        e[n - 1] = 0.0;

        double f = 0.0;
        double tst1 = 0.0;
        double eps = Math.pow(2.0, -52.0);
        for (int l = 0; l < n; l++) {

            // ищем наименьший субдиагональный элемент

            tst1 = Math.max(tst1, Math.abs(d[l]) + Math.abs(e[l]));
            int m = l;
            while (m < n) {
                if (Math.abs(e[m]) <= eps * tst1) {
                    break;
                }
                m++;
            }

            // Если m == l, d[l] - собственное значение,
            // в противном случае выполняем итерацию.

            if (m > l) {
                int iter = 0;
                do {
                    iter = iter + 1;  // проверяем количество итераций здесь

                    // Вычисляем неявный сдвиг

                    double g = d[l];
                    double p = (d[l + 1] - g) / (2.0 * e[l]);
                    double r = Math.hypot(p, 1.0);
                    if (p < 0) {
                        r = -r;
                    }
                    d[l] = e[l] / (p + r);
                    d[l + 1] = e[l] * (p + r);
                    double dl1 = d[l + 1];
                    double h = g - d[l];
                    for (int i = l + 2; i < n; i++) {
                        d[i] -= h;
                    }
                    f = f + h;

                    // Неявное преобразование QL.

                    p = d[m];
                    double c = 1.0;
                    double c2 = c;
                    double c3 = c;
                    double el1 = e[l + 1];
                    double s = 0.0;
                    double s2 = 0.0;
                    for (int i = m - 1; i >= l; i--) {
                        c3 = c2;
                        c2 = c;
                        s2 = s;
                        g = c * e[i];
                        h = c * p;
                        r = Math.hypot(p, e[i]);
                        e[i + 1] = s * r;
                        s = e[i] / r;
                        c = p / r;
                        p = c * d[i] - s * g;
                        d[i + 1] = h + s * (c * g + s * d[i]);

                        // Накапливание преобразований

                        for (int k = 0; k < n; k++) {
                            h = V[k][i + 1];
                            V[k][i + 1] = s * V[k][i] + c * h;
                            V[k][i] = c * V[k][i] - s * h;
                        }
                    }
                    p = -s * s2 * c3 * el1 * e[l] / dl1;
                    e[l] = s * p;
                    d[l] = c * p;

                    // Проверка на сходимость

                } while (Math.abs(e[l]) > eps * tst1);
            }
            d[l] = d[l] + f;
            e[l] = 0.0;
        }

        // Сортировка собственных значений и соответствующих векторов

        for (int i = 0; i < n - 1; i++) {
            int k = i;
            double p = d[i];
            for (int j = i + 1; j < n; j++) {
                if (d[j] < p) {
                    k = j;
                    p = d[j];
                }
            }
            if (k != i) {
                d[k] = d[i];
                d[i] = p;
                for (int j = 0; j < n; j++) {
                    p = V[j][i];
                    V[j][i] = V[j][k];
                    V[j][k] = p;
                }
            }
        }
    }

    // Несимметричное приведение к форме Хессенберга

    public static void orthes(int n, double[] d, double[] e, double[][] V, double[][] H, double[] ort) {



        int low = 0;
        int high = n - 1;

        for (int m = low + 1; m <= high - 1; m++) {

            // Столбец шкалы

            double scale = 0.0;
            for (int i = m; i <= high; i++) {
                scale = scale + Math.abs(H[i][m - 1]);
            }
            if (scale != 0.0) {

                // Вычислить преобразование хаусхолдера

                double h = 0.0;
                for (int i = high; i >= m; i--) {
                    ort[i] = H[i][m - 1] / scale;
                    h += ort[i] * ort[i];
                }
                double g = Math.sqrt(h);
                if (ort[m] > 0) {
                    g = -g;
                }
                h = h - ort[m] * g;
                ort[m] = ort[m] - g;

                // Применить преобразование подобия хаусхолдера
                // H = (I-u*u'/h)*H*(I-u*u')/h)

                for (int j = m; j < n; j++) {
                    double f = 0.0;
                    for (int i = high; i >= m; i--) {
                        f += ort[i] * H[i][j];
                    }
                    f = f / h;
                    for (int i = m; i <= high; i++) {
                        H[i][j] -= f * ort[i];
                    }
                }

                for (int i = 0; i <= high; i++) {
                    double f = 0.0;
                    for (int j = high; j >= m; j--) {
                        f += ort[j] * H[i][j];
                    }
                    f = f / h;
                    for (int j = m; j <= high; j++) {
                        H[i][j] -= f * ort[j];
                    }
                }
                ort[m] = scale * ort[m];
                H[m][m - 1] = scale * g;
            }
        }

        //Накопление преобразований  (Algol's ortran).

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                V[i][j] = (i == j ? 1.0 : 0.0);
            }
        }

        for (int m = high - 1; m >= low + 1; m--) {
            if (H[m][m - 1] != 0.0) {
                for (int i = m + 1; i <= high; i++) {
                    ort[i] = H[i][m - 1];
                }
                for (int j = m; j <= high; j++) {
                    double g = 0.0;
                    for (int i = m; i <= high; i++) {
                        g += ort[i] * V[i][j];
                    }
                    // Двойное разделение позволяет избежать возможного недостаточного потока
                    g = (g / ort[m]) / H[m][m - 1];
                    for (int i = m; i <= high; i++) {
                        V[i][j] += g * ort[i];
                    }
                }
            }
        }
    }


    // Комплексное скалярное деление.

    public static transient double cdivr, cdivi;

    public static void cdiv(double xr, double xi, double yr, double yi) {
        double r, d;
        if (Math.abs(yr) > Math.abs(yi)) {
            r = yi / yr;
            d = yr + r * yi;
            cdivr = (xr + r * xi) / d;
            cdivi = (xi - r * xr) / d;
        } else {
            r = yr / yi;
            d = yi + r * yr;
            cdivr = (r * xr + xi) / d;
            cdivi = (r * xi - xr) / d;
        }
    }


    // Несимметричная редукция от Хессенберга к реальной форме Шура

    public static void hqr2(int nn, double[] d, double[] e, double[][] V, double[][] H, double[] ort) {


        int n = nn - 1;
        int low = 0;
        int high = nn - 1;
        double eps = Math.pow(2.0, -52.0);
        double exshift = 0.0;
        double p = 0, q = 0, r = 0, s = 0, z = 0, t, w, x, y;

        // Храним корни, изолированные балансом, и вычисляем норму матрицы

        double norm = 0.0;
        for (int i = 0; i < nn; i++) {
            if (i < low | i > high) {
                d[i] = H[i][i];
                e[i] = 0.0;
            }
            for (int j = Math.max(i - 1, 0); j < nn; j++) {
                norm = norm + Math.abs(H[i][j]);
            }
        }

        // Внешний цикл по индексу собственного значения

        int iter = 0;
        while (n >= low) {

            // Ищите один меньший субдиагональный элемент

            int l = n;
            while (l > low) {
                s = Math.abs(H[l - 1][l - 1]) + Math.abs(H[l][l]);
                if (s == 0.0) {
                    s = norm;
                }
                if (Math.abs(H[l][l - 1]) < eps * s) {
                    break;
                }
                l--;
            }

            // Проверка на сходимость
            //  Найден один корень

            if (l == n) {
                H[n][n] = H[n][n] + exshift;
                d[n] = H[n][n];
                e[n] = 0.0;
                n--;
                iter = 0;

                // Найдены два корня

            } else if (l == n - 1) {
                w = H[n][n - 1] * H[n - 1][n];
                p = (H[n - 1][n - 1] - H[n][n]) / 2.0;
                q = p * p + w;
                z = Math.sqrt(Math.abs(q));
                H[n][n] = H[n][n] + exshift;
                H[n - 1][n - 1] = H[n - 1][n - 1] + exshift;
                x = H[n][n];

                // реальные значения

                if (q >= 0) {
                    if (p >= 0) {
                        z = p + z;
                    } else {
                        z = p - z;
                    }
                    d[n - 1] = x + z;
                    d[n] = d[n - 1];
                    if (z != 0.0) {
                        d[n] = x - w / z;
                    }
                    e[n - 1] = 0.0;
                    e[n] = 0.0;
                    x = H[n][n - 1];
                    s = Math.abs(x) + Math.abs(z);
                    p = x / s;
                    q = z / s;
                    r = Math.sqrt(p * p + q * q);
                    p = p / r;
                    q = q / r;

                    // Модификация строки

                    for (int j = n - 1; j < nn; j++) {
                        z = H[n - 1][j];
                        H[n - 1][j] = q * z + p * H[n][j];
                        H[n][j] = q * H[n][j] - p * z;
                    }

                    // Модификация столбца

                    for (int i = 0; i <= n; i++) {
                        z = H[i][n - 1];
                        H[i][n - 1] = q * z + p * H[i][n];
                        H[i][n] = q * H[i][n] - p * z;
                    }

                    // накопление трансформаций

                    for (int i = low; i <= high; i++) {
                        z = V[i][n - 1];
                        V[i][n - 1] = q * z + p * V[i][n];
                        V[i][n] = q * V[i][n] - p * z;
                    }

                    // комплексные значения

                } else {
                    d[n - 1] = x + p;
                    d[n] = x + p;
                    e[n - 1] = z;
                    e[n] = -z;
                }
                n = n - 2;
                iter = 0;

                // Конвергенции пока нет

            } else {

                // Сдвиг формы

                x = H[n][n];
                y = 0.0;
                w = 0.0;
                if (l < n) {
                    y = H[n - 1][n - 1];
                    w = H[n][n - 1] * H[n - 1][n];
                }

                // специальная смена Уилкинсона

                if (iter == 10) {
                    exshift += x;
                    for (int i = low; i <= n; i++) {
                        H[i][i] -= x;
                    }
                    s = Math.abs(H[n][n - 1]) + Math.abs(H[n - 1][n - 2]);
                    x = y = 0.75 * s;
                    w = -0.4375 * s * s;
                }

                // MATLAB

                if (iter == 30) {
                    s = (y - x) / 2.0;
                    s = s * s + w;
                    if (s > 0) {
                        s = Math.sqrt(s);
                        if (y < x) {
                            s = -s;
                        }
                        s = x - w / ((y - x) / 2.0 + s);
                        for (int i = low; i <= n; i++) {
                            H[i][i] -= s;
                        }
                        exshift += s;
                        x = y = w = 0.964;
                    }
                }

                iter = iter + 1;

                // Ищите два последовательных небольших субдиагональных элемента

                int m = n - 2;
                while (m >= l) {
                    z = H[m][m];
                    r = x - z;
                    s = y - z;
                    p = (r * s - w) / H[m + 1][m] + H[m][m + 1];
                    q = H[m + 1][m + 1] - z - r - s;
                    r = H[m + 2][m + 1];
                    s = Math.abs(p) + Math.abs(q) + Math.abs(r);
                    p = p / s;
                    q = q / s;
                    r = r / s;
                    if (m == l) {
                        break;
                    }
                    if (Math.abs(H[m][m - 1]) * (Math.abs(q) + Math.abs(r)) <
                            eps * (Math.abs(p) * (Math.abs(H[m - 1][m - 1]) + Math.abs(z) +
                                    Math.abs(H[m + 1][m + 1])))) {
                        break;
                    }
                    m--;
                }

                for (int i = m + 2; i <= n; i++) {
                    H[i][i - 2] = 0.0;
                    if (i > m + 2) {
                        H[i][i - 3] = 0.0;
                    }
                }

                //
                //Двойной шаг QR, включающий строки l:n и столбцы m:n

                for (int k = m; k <= n - 1; k++) {
                    boolean notlast = (k != n - 1);
                    if (k != m) {
                        p = H[k][k - 1];
                        q = H[k + 1][k - 1];
                        r = (notlast ? H[k + 2][k - 1] : 0.0);
                        x = Math.abs(p) + Math.abs(q) + Math.abs(r);
                        if (x != 0.0) {
                            p = p / x;
                            q = q / x;
                            r = r / x;
                        }
                    }
                    if (x == 0.0) {
                        break;
                    }
                    s = Math.sqrt(p * p + q * q + r * r);
                    if (p < 0) {
                        s = -s;
                    }
                    if (s != 0) {
                        if (k != m) {
                            H[k][k - 1] = -s * x;
                        } else if (l != m) {
                            H[k][k - 1] = -H[k][k - 1];
                        }
                        p = p + s;
                        x = p / s;
                        y = q / s;
                        z = r / s;
                        q = q / p;
                        r = r / p;



                        for (int j = k; j < nn; j++) {
                            p = H[k][j] + q * H[k + 1][j];
                            if (notlast) {
                                p = p + r * H[k + 2][j];
                                H[k + 2][j] = H[k + 2][j] - p * z;
                            }
                            H[k][j] = H[k][j] - p * x;
                            H[k + 1][j] = H[k + 1][j] - p * y;
                        }



                        for (int i = 0; i <= Math.min(n, k + 3); i++) {
                            p = x * H[i][k] + y * H[i][k + 1];
                            if (notlast) {
                                p = p + z * H[i][k + 2];
                                H[i][k + 2] = H[i][k + 2] - p * r;
                            }
                            H[i][k] = H[i][k] - p;
                            H[i][k + 1] = H[i][k + 1] - p * q;
                        }



                        for (int i = low; i <= high; i++) {
                            p = x * V[i][k] + y * V[i][k + 1];
                            if (notlast) {
                                p = p + z * V[i][k + 2];
                                V[i][k + 2] = V[i][k + 2] - p * r;
                            }
                            V[i][k] = V[i][k] - p;
                            V[i][k + 1] = V[i][k + 1] - p * q;
                        }
                    }  // (s != 0)
                } // цикл k
            } // проверка сходимости
        } // в то время как (n >= low)

        // Обратная установка для нахождения векторов верхней треугольной формы

        if (norm == 0.0) {
            return;
        }

        for (n = nn - 1; n >= 0; n--) {
            p = d[n];
            q = e[n];

            // реальные векторы

            if (q == 0) {
                int l = n;
                H[n][n] = 1.0;
                for (int i = n - 1; i >= 0; i--) {
                    w = H[i][i] - p;
                    r = 0.0;
                    for (int j = l; j <= n; j++) {
                        r = r + H[i][j] * H[j][n];
                    }
                    if (e[i] < 0.0) {
                        z = w;
                        s = r;
                    } else {
                        l = i;
                        if (e[i] == 0.0) {
                            if (w != 0.0) {
                                H[i][n] = -r / w;
                            } else {
                                H[i][n] = -r / (eps * norm);
                            }

                            //уравнения

                        } else {
                            x = H[i][i + 1];
                            y = H[i + 1][i];
                            q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                            t = (x * s - z * r) / q;
                            H[i][n] = t;
                            if (Math.abs(x) > Math.abs(z)) {
                                H[i + 1][n] = (-r - w * t) / x;
                            } else {
                                H[i + 1][n] = (-s - y * t) / z;
                            }
                        }

                        // контроль переполнения

                        t = Math.abs(H[i][n]);
                        if ((eps * t) * t > 1) {
                            for (int j = i; j <= n; j++) {
                                H[j][n] = H[j][n] / t;
                            }
                        }
                    }
                }

                // комплексные

            } else if (q < 0) {
                int l = n - 1;

                // Последняя векторная составляющая мнимая, поэтому матрица треугольная

                if (Math.abs(H[n][n - 1]) > Math.abs(H[n - 1][n])) {
                    H[n - 1][n - 1] = q / H[n][n - 1];
                    H[n - 1][n] = -(H[n][n] - p) / H[n][n - 1];
                } else {
                    cdiv(0.0, -H[n - 1][n], H[n - 1][n - 1] - p, q);
                    H[n - 1][n - 1] = cdivr;
                    H[n - 1][n] = cdivi;
                }
                H[n][n - 1] = 0.0;
                H[n][n] = 1.0;
                for (int i = n - 2; i >= 0; i--) {
                    double ra, sa, vr, vi;
                    ra = 0.0;
                    sa = 0.0;
                    for (int j = l; j <= n; j++) {
                        ra = ra + H[i][j] * H[j][n - 1];
                        sa = sa + H[i][j] * H[j][n];
                    }
                    w = H[i][i] - p;

                    if (e[i] < 0.0) {
                        z = w;
                        r = ra;
                        s = sa;
                    } else {
                        l = i;
                        if (e[i] == 0) {
                            cdiv(-ra, -sa, w, q);
                            H[i][n - 1] = cdivr;
                            H[i][n] = cdivi;
                        } else {

                            // Solve complex equations

                            x = H[i][i + 1];
                            y = H[i + 1][i];
                            vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                            vi = (d[i] - p) * 2.0 * q;
                            if (vr == 0.0 & vi == 0.0) {
                                vr = eps * norm * (Math.abs(w) + Math.abs(q) +
                                        Math.abs(x) + Math.abs(y) + Math.abs(z));
                            }
                            cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
                            H[i][n - 1] = cdivr;
                            H[i][n] = cdivi;
                            if (Math.abs(x) > (Math.abs(z) + Math.abs(q))) {
                                H[i + 1][n - 1] = (-ra - w * H[i][n - 1] + q * H[i][n]) / x;
                                H[i + 1][n] = (-sa - w * H[i][n] - q * H[i][n - 1]) / x;
                            } else {
                                cdiv(-r - y * H[i][n - 1], -s - y * H[i][n], z, q);
                                H[i + 1][n - 1] = cdivr;
                                H[i + 1][n] = cdivi;
                            }
                        }

                        // Overflow control

                        t = Math.max(Math.abs(H[i][n - 1]), Math.abs(H[i][n]));
                        if ((eps * t) * t > 1) {
                            for (int j = i; j <= n; j++) {
                                H[j][n - 1] = H[j][n - 1] / t;
                                H[j][n] = H[j][n] / t;
                            }
                        }
                    }
                }
            }
        }

        // Vectors of isolated roots

        for (int i = 0; i < nn; i++) {
            if (i < low | i > high) {
                for (int j = i; j < nn; j++) {
                    V[i][j] = H[i][j];
                }
            }
        }

        // Back transformation to get eigenvectors of original matrix

        for (int j = nn - 1; j >= low; j--) {
            for (int i = low; i <= high; i++) {
                z = 0.0;
                for (int k = low; k <= Math.min(j, high); k++) {
                    z = z + V[i][k] * H[k][j];
                }
                V[i][j] = z;
            }
        }
    }

    public static double[] arraySet(double[] array) {
        ArrayList<Double> set = new ArrayList<Double>();
        for (int i = 0; i < array.length; i++) {
            boolean unique = true;
            for (int j = 0; j < set.size(); j++) {
                if (array[i] == set.get(j)) {
                    unique = false;
                }
            }
            if (unique) {
                set.add(array[i]);
            }
        }
        double[] newSet = new double[set.size()];
        for (int i = 0; i < set.size(); i++) {
            newSet[i] = set.get(i);
        }
        return newSet;
    }












    public static double[][] calculateEigenvector(double[][] matrix, double[] d) { //Finds eigenvector from REAL eigenvalues
        double[][] vectors = new double[matrix.length][matrix[0].length];
        int placeIndex = 0;
        for (int det = 0; det < d.length; det++) {//for each eigenvalue
            double[][] temp1 = new double[matrix.length][matrix.length];
            for (int i = 0; i < matrix.length; i++) {
                for (int j = 0; j < matrix.length; j++) {
                    temp1[i][j] = matrix[i][j];
                }
            }
            for (int i = 0; i < matrix.length; i++) {//A-eigen(I)
                temp1[i][i] -= d[det];
            }
            //Check for any variables not in the function
            double[][] temp2 = transpose(rref(temp1));
            for (int i = 0; i < temp1.length; i++) {

                if (countNonZero(temp2[i]) == 0) {
                    double[] temp = new double[temp2.length];
                    temp[i] = 1.0;

                    vectors[placeIndex] = temp;
                    placeIndex++;
                }
            }

            //Check for equations in rows of functions
            double[][] temp3 = rref(temp1);
            double[] temp = new double[temp3[0].length];
            for (int i = 0; i < temp3.length; i++) {
                if (countNonZero(temp3[i]) <= 1) {
                    continue;
                }
                boolean pivot = false;
                for (int j = 0; j < temp1[0].length; j++) {

                    if (pivot) {
                        temp[j] = 0 + temp3[i][j];
                        if ((j == temp1[0].length - 1) && (countNonZero(temp) != 0)) {
                            for (int k = 0; k < temp.length; k++) {
                                vectors[placeIndex][k] = temp[k];
                            }
                            placeIndex++;
                        }
                    }
                    if (temp3[i][j] != 0) {
                        temp[j] -= temp3[i][j];
                        pivot = true;
                    }

                }
            }







        }
        return vectors;
    }

    public static int countNonZero(double[] row) {
        int count = 0;
        for (int i = 0; i < row.length; i++) {
            if (row[i] != 0.0) {
                count++;
            }
        }
        return count;
    }
}