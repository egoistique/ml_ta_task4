package ru.vsu.cs.course1;

public class Gauss{
    public static final double EPS = 1E-5;

    public static void Swap_Lines(int k1, int k2, int n, double[][] A, Boolean[] mark) {
        for (int j = 0; j < n; j++) {
            double tmp;
            tmp = A[k1][j];
            A[k1][j] = A[k2][j];
            A[k2][j] = tmp;
        }
        Boolean tmp;
        tmp = mark[k1];
        mark[k1] = mark[k2];
        mark[k2] = tmp;
    }

    public static double[][] gauss(double[][] A){
        int m = A.length, n = A[0].length - 1;
        Double[] answer = new Double[n];

        //Выполняем прямой ход метода Гаусса
        int min_size = Math.min(m, n);

        for (int k = 0; k < min_size; k++) {
            double maxv = 0; int maxLinePos = k;
            for (int i = k; i < m; i++) {
                if (Math.abs(A[i][k]) > maxv) {
                    maxv = Math.abs(A[i][k]);
                    maxLinePos = i;
                }
            }
            for (int j = 0; j < n + 1; j++) {
                double tmp = A[k][j];
                A[k][j] = A[maxLinePos][j];
                A[maxLinePos][j] = tmp;
            }

            if (Math.abs(maxv) < EPS) {
                continue;
            }
            //Делаем ступенчатый вид
            for (int i = 0; i < m; i++) {
                if (i == k) continue;

                double multiplier = A[i][k]/A[k][k];
                for (int j = k; j < n+1; j++) {
                    A[i][j] -= multiplier * A[k][j];
                }
            }
        }

        //Делим каждый коэффициент i- го уравнения на первый ненулевой коэффициент этого уравнения
        for (int k = 0; k < min_size; k++) {
            if (Math.abs(A[k][k]) > EPS) {
                double multiplier = A[k][k];
                if (Math.abs(multiplier) < EPS) continue;
                for (int j = k; j < n+1; j++) {
                    A[k][j] /= multiplier;
                }
            }
        }

        //Отмечаем одинаковые строки в расширенной матрице A
        Boolean[] mark = new Boolean[m];
        for (int i = 0; i < m; i++) {
            mark[i] = Boolean.FALSE;
        }

        for (int k1 = 0; k1 < m; k1++) {
            if (mark[k1] == Boolean.TRUE) continue;
            for (int k2 = k1+1; k2 < m; k2++) {
                boolean is_equal = true;
                for (int j = 0; j < n+1; j++) {
                    if (Math.abs(A[k1][j] - A[k2][j]) > EPS) {
                        is_equal = false;
                        break;
                    }
                }
                if (is_equal) {
                    mark[k2] = true;
                }
            }
        }

        //Проверяем на совместность
        for (int i = 0; i < m; i++) {
            int zeroes = 0;
            for (int j = 0; j < n+1; j++) {
                if (Math.abs(A[i][j]) < EPS) {
                    zeroes++;
                    A[i][j] = 0.0;
                }
            }
            if (zeroes == n+1) {
                mark[i] = Boolean.TRUE;
            }
            if (zeroes == n && Math.abs(A[i][n]) > EPS) {
                System.out.println("The system of equations is inconsistent");
                return null;
            }
        }

        //Все ненулевые строки переносим вперёд
        for (int i = 0; i < m; i++) {
            for (int j = i+1; j < m; j++) {
                if (mark[i] == Boolean.TRUE && mark[j] == Boolean.FALSE) {
                    Swap_Lines(i, j, n, A, mark);
                }
            }
        }

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n+1; j++) {
                System.out.print(A[i][j] + " ");
            }
            System.out.println();
        }

        //Если ранг совпадает с количеством уравнений,
        // то система имеет единственное решение
        int marks = 0;
        for (int i = 0; i < m; i++) {
            if (mark[i] == Boolean.TRUE) marks++;
        }
        int rank = m-1-marks;

        if (rank == n-1) {
            for (int k = n-1; k >= 0; k--) {
                answer[k] = A[k][n] / A[k][k];
            }

            System.out.println("Anser:");
            for (int k = 0; k < n-1; k++) {
                System.out.print(answer[k] + " ");
            }
            System.out.println(answer[n-1]);
        }//Иначе отмечаем свободные переменные
        else {
            int freeVars = n-(rank+1);

            Boolean[] markedVars = new Boolean[n];
            for (int i = 0; i < n; i++) {
                markedVars[i] = Boolean.FALSE;
            }

            for (int j = 0; j < n; j++) {
                int zeroes = 0;
                for (int i = 0; i < rank; i++) {
                    if (Math.abs(A[i][j]) < EPS) {
                        zeroes++;
                    }
                }
                if (zeroes == rank+1) {
                    if (freeVars > 0) {
                        markedVars[j] = Boolean.TRUE;
                        freeVars--;
                    }
                }
            }
            //Инициализируем свободные переменные
            for (int i = n-1; i >= 0; i--) {
                if (freeVars == 0) break;
                markedVars[i] = Boolean.TRUE;
                freeVars--;
            }
            System.out.println("Initialization of free variables:");
            for (int i = 0; i < n; i++) {
                if (markedVars[i] == Boolean.TRUE) {
                    answer[i] = 1.0;
                    System.out.println("Let: " + i + "-th variable assigned: 1.0" );
                }
            }
            //Выводим в консоль результирующую матрицу
            System.out.println("Answer:");
            for (int i = 0; i < n; i++) {
                if (markedVars[i] == Boolean.TRUE) {
                    System.out.println(i+"-th variable is free");
                }
            }

            for (int i = rank; i >= 0; i--) {
                double cur_sum = 0;

                int curVars = 0;
                for (int j = 0; j < n; j++) {
                    if (markedVars[j] == Boolean.FALSE && Math.abs(A[i][j]) > EPS) {
                        curVars = j;
                        break;
                    }
                }

                System.out.print("X[" + curVars + "] = ");
                for (int j = 0; j < n; j++) {
                    if (markedVars[j] == Boolean.TRUE) {
                        cur_sum += answer[j] * A[i][j];
                        System.out.print("(" + -A[i][j] + "/" + A[i][curVars] + ")" + "*X[" + j + "] ");
                    }
                }
                System.out.println();

                //Выводим в консоль значения зависимых переменных
                cur_sum *= -1;
                cur_sum += A[i][n];


                for (int j = 0; j < n; j++) {
                    if (markedVars[j] == Boolean.FALSE && Math.abs(A[i][j]) > EPS) {
                        answer[j] = cur_sum / A[i][j];
                        markedVars[j] = Boolean.TRUE;
                        break;
                    }
                }

            }

/*
            for (int i = 0; i < n; i++) {
                if (Math.abs(answer[i]) < EPS) answer[i] = 0.0;
            }
*/
            //Вывод данных в консоль
            double[][] answers = new double[A[0].length - 1][1];
            System.out.println("One of the solutions:");
            for (int k = 0; k < n-1; k++) {
                if (answer[k] ==  null){
                    answers[k][0] =  0;
                }
                else answers[k][0] =  answer[k];
                System.out.print(answer[k] + " ");
            }
            answers[n-1][0] = answer[n-1];
            System.out.println(answer[n-1]);
            return answers;
        }

        double[][] answers = new double[A[0].length - 1][1];
        for (int i = 0; i < A.length; i++){
            answers[i][0] =  answer[i];
        }
        return answers;
    }

    public static void main(String[] args) {
//        Scanner sc = new Scanner(System.in);
//        int n, m;
//        System.out.print("Enter number of equation: ");
//        m = sc.nextInt();
//        System.out.print("Enter number of variables: ");
//        n = sc.nextInt();
//
//        Double[][] A = new Double[m][n+1];
//        Double[] answer = new Double[n];
//
//        System.out.println("Enter matrix:");
//        for (int i = 0; i < m; i++) {
//            for (int j = 0; j < n+1; j++) {
//                A[i][j] = sc.nextDouble();
//            }
//        }
//
//        int min_size = Math.min(m, n);
//
//        for (int k = 0; k < min_size; k++) {
//            double maxv = 0; int position_of_line_with_maxv = k;
//            for (int i = k; i < m; i++) {
//                if (Math.abs(A[i][k]) > maxv) {
//                    maxv = Math.abs(A[i][k]);
//                    position_of_line_with_maxv = i;
//                }
//            }
//            for (int j = 0; j < n+1; j++) {
//                double tmp = A[k][j];
//                A[k][j] = A[position_of_line_with_maxv][j];
//                A[position_of_line_with_maxv][j] = tmp;
//            }
//
//            if (Math.abs(maxv) < EPS) {
//                continue;
//            }
//
//            for (int i = 0; i < m; i++) {
//                if (i == k) continue;
//
//                double multiplier = A[i][k]/A[k][k];
//                for (int j = k; j < n+1; j++) {
//                    A[i][j] -= multiplier * A[k][j];
//                }
//            }
//        }
//
//        for (int k = 0; k < min_size; k++) {
//            if (Math.abs(A[k][k]) > EPS) {
//                double multiplier = A[k][k];
//                if (Math.abs(multiplier) < EPS) continue;
//                for (int j = k; j < n+1; j++) {
//                    A[k][j] /= multiplier;
//                }
//            }
//        }
//
//        Boolean[] mark = new Boolean[m];
//        for (int i = 0; i < m; i++) {
//            mark[i] = Boolean.FALSE;
//        }
//
//        for (int k1 = 0; k1 < m; k1++) {
//            if (mark[k1] == Boolean.TRUE) continue;
//            for (int k2 = k1+1; k2 < m; k2++) {
//                boolean is_equal = true;
//                for (int j = 0; j < n+1; j++) {
//                    if (Math.abs(A[k1][j] - A[k2][j]) > EPS) {
//                        is_equal = false;
//                        break;
//                    }
//                }
//                if (is_equal) {
//                    mark[k2] = true;
//                }
//            }
//        }
//        for (int i = 0; i < m; i++) {
//            int cnt_of_zeroes = 0;
//            for (int j = 0; j < n+1; j++) {
//                if (Math.abs(A[i][j]) < EPS) {
//                    cnt_of_zeroes++;
//                    A[i][j] = 0.0;
//                }
//            }
//            if (cnt_of_zeroes == n+1) {
//                mark[i] = Boolean.TRUE;
//            }
//            if (cnt_of_zeroes == n && Math.abs(A[i][n]) > EPS) {
//                System.out.println("The system of equations is inconsistent");
//                return;
//            }
//        }
//
//        for (int i = 0; i < m; i++) {
//            for (int j = i+1; j < m; j++) {
//                if (mark[i] == Boolean.TRUE && mark[j] == Boolean.FALSE) {
//                    Swap_Lines(i, j, n, A, mark);
//                }
//            }
//        }
//
//        for (int i = 0; i < m; i++) {
//            for (int j = 0; j < n+1; j++) {
//                System.out.print(A[i][j] + " ");
//            }
//            System.out.println();
//        }
//
//        int cnt_of_marks = 0;
//        for (int i = 0; i < m; i++) {
//            if (mark[i] == Boolean.TRUE) cnt_of_marks++;
//        }
//        int bottom_border = m-1-cnt_of_marks;
//
//        if (bottom_border == n-1) {
//            for (int k = n-1; k >= 0; k--) {
//                answer[k] = A[k][n] / A[k][k];
//            }
//
//            System.out.println("Anser:");
//            for (int k = 0; k < n-1; k++) {
//                System.out.print(answer[k] + " ");
//            }
//            System.out.println(answer[n-1]);
//        }
//        else {
//            int cnt_of_free_variables = n-(bottom_border+1);
//
//            Boolean[] marked_variables = new Boolean[n];
//            for (int i = 0; i < n; i++) {
//                marked_variables[i] = Boolean.FALSE;
//            }
//
//            for (int j = 0; j < n; j++) {
//                int cnt_of_zeroes = 0;
//                for (int i = 0; i < bottom_border; i++) {
//                    if (Math.abs(A[i][j]) < EPS) {
//                        cnt_of_zeroes++;
//                    }
//                }
//                if (cnt_of_zeroes == bottom_border+1) {
//                    if (cnt_of_free_variables > 0) {
//                        marked_variables[j] = Boolean.TRUE;
//                        cnt_of_free_variables--;
//                    }
//                }
//            }
//            for (int i = n-1; i >= 0; i--) {
//                if (cnt_of_free_variables == 0) break;
//                marked_variables[i] = Boolean.TRUE;
//                cnt_of_free_variables--;
//            }
//            System.out.println("Initialization of free variables:");
//            for (int i = 0; i < n; i++) {
//                if (marked_variables[i] == Boolean.TRUE) {
//                    answer[i] = 1.0;
//                    System.out.println("Let: " + i + "-th variable assigned: 1.0" );
//                }
//            }
//            System.out.println("Answer:");
//            for (int i = 0; i < n; i++) {
//                if (marked_variables[i] == Boolean.TRUE) {
//                    System.out.println(i+"-th variable is free");
//                }
//            }
//
//            for (int i = bottom_border; i >= 0; i--) {
//                double cur_sum = 0;
//
//                int cur_variable = 0;
//                for (int j = 0; j < n; j++) {
//                    if (marked_variables[j] == Boolean.FALSE && Math.abs(A[i][j]) > EPS) {
//                        cur_variable = j;
//                        break;
//                    }
//                }
//
//                System.out.print("X[" + cur_variable + "] = ");
//                for (int j = 0; j < n; j++) {
//                    if (marked_variables[j] == Boolean.TRUE) {
//                        cur_sum += answer[j] * A[i][j];
//                        System.out.print("(" + -A[i][j] + "/" + A[i][cur_variable] + ")" + "*X[" + j + "] ");
//                    }
//                }
//                System.out.println();
//
//                cur_sum *= -1;
//                cur_sum += A[i][n];
//
//
//                for (int j = 0; j < n; j++) {
//                    if (marked_variables[j] == Boolean.FALSE && Math.abs(A[i][j]) > EPS) {
//                        answer[j] = cur_sum / A[i][j];
//                        marked_variables[j] = Boolean.TRUE;
//                        break;
//                    }
//                }
//
//            }
//
///*
//            for (int i = 0; i < n; i++) {
//                if (Math.abs(answer[i]) < EPS) answer[i] = 0.0;
//            }
//*/
//
//            System.out.println("One of the solutions:");
//            for (int k = 0; k < n-1; k++) {
//                System.out.print(answer[k] + " ");
//            }
//            System.out.println(answer[n-1]);
//        }
    }
}
