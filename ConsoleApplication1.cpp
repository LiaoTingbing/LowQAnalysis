// ConsoleApplication1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <cmath>
#include <complex>
#include<mat.h>
#include<string.h>

#define arrayLength( a ) sizeof( a ) / sizeof( a[0] )
#define arrayLengthPointer( a ) _msize( a ) / sizeof( a[0] )
#define PI 3.14159265358979323846

using namespace std; 

/*  函数  */

static complex<double> fourierw1(double * s, double* t, double w, int slength) {
    double dt = *(t + 1) - *t;
    complex<double> sw1(0,0);
    for (int i = 0; i < slength; i++) {
        sw1 = sw1 + exp(1i * w * (*(t + i))) * (*(s + i)) * dt; 
    }
    return sw1; 
}

static double* fftw1(double* t, size_t num) {
    double dt = *(t + 1) - *(t);
    int  pp = int(ceil(log2(num)));
    int Mp = int(pow(2, pp));

    double* w = new double[Mp]();

    // w[0] = Mp; 
    for (int i = 0; i < Mp; i++) {
        *(w+i) = 2 * PI / dt / Mp * i;
    }
    return w;
}

static complex<double>* realczt(double* signal, double* t, int lengthSignal, double* w, int lengthw) {
    complex<double>* fs = new complex <double>[lengthw]();
    complex<double> s;
    for (int j = 0; j < lengthw; j++) {
        s = (0, 0);
        for (int i = 0; i < lengthSignal; i++) {
            s = s + *(signal + i) * exp(1i * *(t + i) * *(w + j));
        }
        *(fs+j) = s;
    }
    return fs;
}

static double maxData(double* p, int lengthData) {
    // size_t  lengthData = _msize(p) / sizeof(p[0]);
    double ma = *p;
    for (int i = 0; i < lengthData; i++) {
        if (*(p+i) > ma) {
            ma = *(p + i);
        }
    }
    return ma;
}

static complex<double> complex_dot(complex<double>* s1, complex<double>* s2, int lengthdata) {
    //size_t num = _msize(s1) / sizeof(s1[0]);
    complex<double> s(0, 0);
    for (int i = 0; i < lengthdata; i++) {
        s = s + *(s1+i)  *   *(s2 + i);
    }
    return s;
}

static complex<double>* sourcenorm(double* time_signal, int r, int c, double* time, double* phase, double* w, int lengthw) {
    /*   r*c  大小  */
    double dt = *(time+1) - *time;
    complex<double>* sw = new complex<double>[lengthw]();
    complex<double> s1(0,0);
    complex<double> s4(c,0);
    double w1 = 0;
    double* ps=NULL;
    double* pt = NULL;
    /*频率点数*/
    for (int k = 0; k < lengthw; k++) {
        w1 = *(w + k);
        s1 = (0, 0);
        // 监视器个数
        for (int j = 0; j < c; j++) {
            ps = time_signal + r * j;
            pt = time + r * j;
            s1 = s1 + fourierw1(ps, pt, w1, r)* exp(-1i*  ( *(phase+j))/180.0*PI *1.0   );
        }
        *(sw+k) = s1 /s4;
    }
    return sw;
}

static double* complexArrayPointAbsSquare(complex <double>* array, int length) {
    double* result = new double[length]();
    for (int i = 0; i < length; i++) {
        *(result + i) = pow(abs(array[i]), 2);
    }
    return result;
}

static double* realArrayPointAdd(double* s1, double* s2, int dataLength) {
    //double* result = new double[dataLength]();

    for (int i = 0; i < dataLength; i++) {
        *(s1 + i) = *(s1 + i) + *(s2 + i);
    }
    return s1;
}

static double* realArrayPointDivision( double* s1, double* s2, int dataLength) {
    double* result = new double[dataLength]();
    for (int i = 0; i < dataLength; i++) {
        *(result + i) = *(s1 + i) / *(s2 + i);
    }
    return result;

}

static int* findPeaksSortIndex(double* fs, int fslength, int number_peaks) {
    /*     fs length(fs)             */
    //int number_peaks = 10;
    int* peakIndex = new int[fslength]();
    double* peakValue = new double[fslength]();

    int p = 0;
    for (int i = 1; i < fslength - 1; i++) {
        if ((*(fs + i - 1) <= *(fs + i)) && (*(fs + i) >= *(fs + i + 1))) {
            *(peakIndex + p) = i;
            *(peakValue + p) = *(fs + i);
            p = p + 1;
        }
    }
    int num = 0;
    for (int i = 0; i < fslength; i++) {
        if (*(peakIndex + i) != 0) {
            num++;
        }
    }
    /*  峰值排序 */
    int* peakIndexSort = new int[num]();
    double* peakValueSort = new double[num]();
    int* tmpIndexSort = new int[num]();
    for (int i = 0; i < num; i++) {
        p = 0;
        for (int j = 0; j < num; j++) {
            if (*(peakValue + i) <= *(peakValue + j)) {
                p++;
            }
        }
        *(tmpIndexSort + i) = p - 1;
    }
    for (int i = 0; i < num; i++) {
        p = *(tmpIndexSort + i);
        *(peakIndexSort + p) = *(peakIndex + i);
        *(peakValueSort + p) = *(peakValue + i);
    }




    int* peakSortindex = new int[number_peaks]();

    for (int i = 0; i < number_peaks; i++) {
        *(peakSortindex + i) = *(peakIndexSort + i);
    }


    //for (int i = 0; i < number_peaks; i++) {
    //    cout << *(peakSortindex + i) << endl;
    //}

   // cout <<num << endl;

    
    delete[] peakIndex, peakValue, peakIndexSort, peakValueSort, tmpIndexSort;

    return peakSortindex;
}



/* ------------------------------------------ */
int main()
{
    clock_t t1, t2;
    t1 = clock();

    /*加载数据*/
    MATFile* pmatFile = matOpen("E3d.mat" , "r");

    // 采样时间
    mxArray* pMxArray = matGetVariable(pmatFile, "t");
    double *t  = (double*)mxGetData(pMxArray);
    size_t tLengthM = mxGetM(pMxArray);  //M  matlab行数 ， Nmatlab列数 ，matlab: M*N

    // 时域信号时间
    pMxArray = matGetVariable(pmatFile, "time");
    double* time = (double*)mxGetData(pMxArray);
    size_t timeLengthM = mxGetM(pMxArray);
    size_t timeLengthN= mxGetN(pMxArray);
    
    // 时域信号
    pMxArray = matGetVariable(pmatFile, "time_signal");
    double* time_signal = (double*)mxGetData(pMxArray);

    //相位
    pMxArray = matGetVariable(pmatFile, "phase");
    double* phase = (double*)mxGetData(pMxArray);
    size_t phaseLengthM = mxGetM(pMxArray);

    // 电场值
    pMxArray = matGetVariable(pmatFile, "E");
    double * E = (double*)mxGetData(pMxArray);
    size_t ELengthM = mxGetM(pMxArray);
    size_t ELengthN = mxGetN(pMxArray);

    /*       开始计算   */

    int frequency_points =int( pow(2,15 )); 

    double* w1 = fftw1(t, tLengthM);
    int  w1length = int(arrayLengthPointer(w1));
    //cout <<  *(w1+ w1length-1)<<endl;
    
    complex<double>* sw1 = sourcenorm(time_signal, int(timeLengthM), int(timeLengthN), time, phase, w1, int(w1length));
    double* sn1 = complexArrayPointAbsSquare(sw1, w1length);
    cout << "sw1 5 " << endl;
    for (int i = 0; i < 5; i++) {
        cout << *(sw1 + i) << endl; 
    }

    double max_sn1 = maxData(sn1, w1length);
    //cout << max_sn1 << endl;
    
    int p1 = 0 , p2=0;
    for (p1 = 0; *(sn1+p1) < 1e-2 * max_sn1; p1++) {  }
    for (p2 = p1; *(sn1 +p2) < 0.9 * max_sn1; p2 ++) {  }
    for (p2 = p2; (*(sn1 + p2) > 1e-2 * max_sn1) && (p2 < w1length); p2 = p2 + 1) { 1; }

    double* w2 = new double[frequency_points]();
    for (int i = 0; i < frequency_points; i++) {
        *(w2+i) = *(w1+p1) + i * (*(w1 + p2) - *(w1 + p1)) / (frequency_points - 1);
    }
    // cout << frequency_points << endl;
    complex<double>* sw2 = sourcenorm(time_signal, int(timeLengthM), int(timeLengthN), time, phase, w2, frequency_points);
    double* sn2 = complexArrayPointAbsSquare(sw2, frequency_points);
    
    cout << "sw2 5 " << endl;
    for (int i = 0; i < 5; i++) {
        cout << *(sw2 + i) << endl; 
    }

    cout <<"czt now " << endl;

    /*  计算 czt 变换 */
    double* fs = new double[arrayLengthPointer(w2)]();
    double* Ep = NULL;
    double extra_factor = 376.73;
    for (int j = 0; j < 6; j++) {
        if (  (j - int(j/6.0)*6) > 2.5 ) {
            extra_factor = 376.73;
        }
        else
        {
            extra_factor = 1.0;
        }


        cout << "j == " << j <<"/54"<< endl;
        Ep = (E + ELengthM * j);
        complex<double>*  fs1  = realczt(Ep,  t,  int (tLengthM ), w2, frequency_points);
        //double* fs2 = complexArrayPointAbsSquare(fs1, arrayLengthPointer(w2));
        for (int i = 0; i < frequency_points; i++) {
            *(fs + i) = *(fs + i) +pow ( abs(*(fs1 + i)* extra_factor),2);
        }
        delete[]fs1  ;
    }

    double* fs_spectrum = realArrayPointDivision(fs, sn2, frequency_points);

    cout << "fs_spectrum 10 " << endl;
    for (int i = 0; i < 10; i++) {
        cout << *(fs_spectrum + i) << endl;
    }

    /* findpeaks    */
    int numberfrequency = 10;
    int* index = findPeaksSortIndex(fs_spectrum, frequency_points, numberfrequency);
    int  p = 0, continue_search = 0;
    double peak_value=0 , FWHM=0 ,dw=0 , FWHM_delta=0 , peak_delta=0;
    double* Q = new double[numberfrequency];
    double* delta_Q = new double[numberfrequency];
    double* f0 = new double[numberfrequency];
    for (int i = 0; i < numberfrequency; i++) {
         p = *(index + i); 
         peak_value = *(fs_spectrum +  p ); 
         continue_search = 1;
        for (p1 = p; (p1 >= 0) && (continue_search); 1) {
            if (*(fs_spectrum + p1) <= peak_value / 2) {
                continue_search = 0;
            }
            else
            {
                p1--;
            } 
        }
        continue_search = 1;
        for ( p2 = p + 1; (p2 <= frequency_points-1) && continue_search; 1) {
            if (*(fs_spectrum + p2) <= peak_value / 2) {
                continue_search = 0;
            }
            else
            {
                p2++;
            }
        }
        if (p1 < 0) { p1 = 0; }
        if (p2 > frequency_points-1 ) { p2 = frequency_points - 1; }

        FWHM = *(w2 + p2) - *(w2 + p1); 
        dw = *(w2 + 1) - *w2;
        FWHM_delta = 2 * dw;
        peak_delta = dw;

        *(Q + i) = *(w2 + p) / FWHM; 
        *(delta_Q+i) = FWHM_delta / FWHM * *(Q + i) + dw / *(w2 + p) * *(Q + i);
        *(f0 + i) = *(w2 + p) / 2 / PI * 1e-12;
    }
    cout << "f0\t" << "Q\t" << "dQ\t" << endl;
    for (int i = 0; i < numberfrequency; i++) {
        cout <<*(f0+i)<< "\t"<<*(Q+i) << "\t"<< *(delta_Q + i) << "\t" << endl;
    }

    delete[]    w1, sw1, sn1, w2, sw2, sn2, fs , fs_spectrum , Q, delta_Q , f0;

    t2 = clock();
    cout << "Time comsume = " << double(t2 - t1) / CLOCKS_PER_SEC <<"s"<< endl;

}


