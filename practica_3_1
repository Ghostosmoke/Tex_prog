#include <iostream>
#include <cmath>
#include "GnuP.h"
#include <fstream>

//#include <mgl2/qt.h>
using namespace std;
double* SLAU(double** a, double* y, int n) {
    double* x, max;
    int k, index;
    const double eps = 0.00001; // точность
    x = new double[n];
    k = 0;
    while (k < n) {
        // Поиск строки с максимальным a[i][k]
        max = abs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)

            if (abs(a[i][k]) > max) {

                max = abs(a[i][k]);
                index = i;
            }

        // Перестановка строк
        for (int j = 0; j < n; j++) {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }

        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;
        // Нормализация уравнений
        for (int i = k; i < n; i++) {
            double temp = a[i][k];
            if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] / temp;
            y[i] = y[i] / temp;
            if (i == k) continue; // уравнение не вычитать само из себя
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] - a[k][j];
            y[i] = y[i] - y[k];
        }
        k++;
    }
    // обратная подстановка
    for (k = n - 1; k >= 0; k--) {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * x[k];
    }
    return x;
}
double** matrix_new(double** A, int N) {
    for (int i = 0; i < N; i++)
        A[i] = new double[N];
    return A;
}
double** zero_matrix(double** A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = 0;
        }

    }
    return A;
}
double* zero_matrix_b(double* b, int n) {
    for (int i = 0; i < n; i++) {
        b[i] = 0;
    }
    return b;
}
double** polinom_k(double** matrix, int col, double* x, double* y, int) {
    double x_sum = 0;
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < col; j++) {

            for (int k = 0; k < col; k++) {
                x_sum += pow(x[k], i + j);
            }
            matrix[i][j] = 0;
            x_sum = 0;
        }


    }
    return matrix;
}
double* polinom_k_b(int col, double* matrix_b, double* x, double* y) {
    double x_sum = 0;
    for (int i = 0; i < col; i++) {
        for (int k = 0; k < col; k++) {
            x_sum += pow(x[k], i) * y[i];
        }
        matrix_b[i] = x_sum;
        x_sum = 0;

    }
    return matrix_b;
}
double absolute_error(double y_none, double y_get) {
    return abs(y_none - y_get);
}
double relative_error(double y_none, double y_absolute) {
    return abs(y_absolute / y_none * 100);
}

int main() {
	GnuP gp;
    setlocale(LC_ALL, "rus");
    double* a, ** A, * b, * y1, * y2,*x,*y;
    const int n = 9;
    int N = 4;
    x = new double[n];
    y = new double[n];
    cout<<"K(s)=Aе^(bs)"<<endl;
	//K(s)=Aе^(bs) n=9
    //double x[n] = { 0,0.5,1,1.5,2,2.5,3,3.5,4 };
    //double y[n] = { 2.31,2.899,3.534,4.412,5.578,6.92,8.699,10.69,13.39};
	ifstream f("points.txt");
	if (f.is_open()) {
	for (int i = 0; i < n; i++) {
				f >> x[i];
				f >> y[i];
	}
	f.close();}
	//Y=Ax^3+Bx^2+Cx+D n=10
    /*double x[n] = { 1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3 };
    double y[n] = { 1.5,2.7,3.9,5.5,7.1,9.1,11.1,12.9,15.5,17.9 };*/
    

    A = new double* [N];
    b = new double[N];
    a = new double[2];
    y1 = new double[n];
    y2 = new double[n];
    A = matrix_new(A, N);
    A = zero_matrix(A, N);
    b = zero_matrix_b(b, N);
    //A[0][0] = n;

    /*for (int i = 0; i < n; i++) {
        A[0][1] += x[i];
        A[1][0] += x[i]; A[1][1] = x[i] * x[i];
        if (y[i] != 0) b[0] += log(y[i]);
        if (x[i] * y[i] != 0) b[1] += log(y[i] * x[i]);
    }
    a = SLAU(A, b, N);*/
    /*for (int i = 0; i < N; i++) {
        cout << a[i] << endl;
    }*/
    double r = 0, Mx = 0, My = 0,My_log=0;
    for (int i = 0; i < n; i++) {
        Mx += x[i]; My += y[i];
        My_log += log(y[i]);
    }
    double yx = 0, xx = 0,yx_no_log=0;
    for (int i = 0; i < n; i++) {
        yx += log(y[i]) * x[i];
        xx += x[i] * x[i];
        yx_no_log+=y[i]* x[i];
    }
    a[1] = (n * yx - My_log * Mx) / (n * xx - Mx * Mx);
    a[0] = (My_log - a[1] * Mx)/n;
    a[0] = exp(a[0]);
    double *a1;
    a1 = new double[2];
    a1[1] = (n * yx_no_log - My * Mx) / (n * xx - Mx * Mx);
    a1[0] = (My - a1[1] * Mx)/n;
    cout << "a:" << a[0] << "*exp^(" << a[1] << "*x)" << endl;
    cout << "\nЛиния регрессии: у=" << a1[0] << "+" << a1[1] << "х" << endl;
    double sum_y_none = 0,sum_y_r=0;
    cout<<"Точки(данные точки, мои точки, линия регрессии)"<<endl;

    for (int i = 0; i < n; i++)
    {
        y1[i] = a[0] * exp(a[1] * x[i]);
        y2[i] = a1[0] + a1[1] * x[i];
        cout<<y[i]<<"\t"<<y1[i]<<"\t"<<y2[i]<<endl;

        sum_y_none += y1[i];
        sum_y_r+=y2[i];
    }

    

    double chisl = 0, znam1 = 0, znam2 = 0;

    for (int i = 0; i < n; i++) {

        chisl += (x[i] - Mx / n) * (y[i] - My / n);
        znam1 += pow((x[i] - Mx / n),2);
        znam2 += pow((y[i] - My / n),2);

    }
    r = chisl / sqrt(znam1 * znam2);
    cout << "Коэффициент корреляции: r=" << r << endl;
    
    double q = absolute_error(sum_y_none, My);
    cout<<"Сумма данных точек и моих точек:"<<My<<"\t"<<sum_y_none<<endl;
	cout<<"Ошибки с функцией"<<endl;
    cout << "Суммарная ошибка: " << q << endl;
    cout << "Средняя ошибка: " << sqrt(q*q)/n << endl;
    cout << "Относительная ошибка: " << sqrt(q*q)/My << endl;
    
    cout<<endl;
    
    cout<<"Ошибки с регресии"<<endl;
    q = absolute_error(sum_y_r, My);
    cout << "Суммарная ошибка: " << q << endl;
    cout << "Средняя ошибка: " << sqrt(q*q)/n << endl;
    cout << "Относительная ошибка: " << sqrt(q*q)/My << endl;
    
	//gp.plotArray(n, x, y, y1, y2);
	gp.plotArrayPar(n,x,y,1,4,11,"Данная функция");
	gp.plotArrayPar(n,x,y1,3,2,12,"Моя функция");
	gp.plotArrayPar(n,x,y2,3,2,1,"Линия регрессии");

	gp.plot();

    return 0;
}
