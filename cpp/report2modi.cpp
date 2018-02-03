#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

#define mx 401
#define i1 95
#define i2 105
#define my 201
#define j1 95
#define j2 105
#define maxitp 10000
#define omegap 1.0
#define errorp 1.0e-4
#define Re 70.0//レイノルズ数

void bcforp(double p[mx][my]){
    for (int j = 0;j < my;j++){
        p[0][j] = 0.0;
        p[mx - 1][j] = 0.0;
    }
    for (int i = 0;i < mx;i++){
        p[i][0] = 0.0;
        p[i][my - 1] = 0.0;
    }
    p[i1][j1] = p[i1 - 1][j1 - 1];
    p[i1][j2] = p[i1 - 1][j2 + 1];
    p[i2][j1] = p[i2 + 1][j1 - 1];
    p[i2][j2] = p[i2 + 1][j2 + 1];
    for (int j = j1 + 1;j < j2;j++){
        p[i1][j] = p[i1 - 1][j];
        p[i2][j] = p[i2 + 1][j];
    }
    for (int i = i1 + 1;i < i2;i++){
        p[i][j1] = p[i][j1 - 1];
        p[i][j2] = p[i][j2 + 1];
    }
}

void bcforv(double u[mx][my], double v[mx][my]){
    for (int j = 0;j < my;j++){
        u[0][j] = 1.0;
        v[0][j] = 0.0;
        u[mx - 1][j] = 2.0 * u[mx - 2][j] - u[mx - 3][j];
        v[mx - 1][j] = 2.0 * v[mx - 2][j] - v[mx - 3][j];
    }
    for (int i = 0;i < mx;i++){
        u[i][0] = 2.0 * u[i][1] - u[i][2];
        v[i][0] = 2.0 * v[i][1] - v[i][2];
        u[i][my - 1] = 2.0 * u[i][my - 2] - u[i][my - 3];
        v[i][my - 1] = 2.0 * v[i][my - 2] - v[i][my - 3];
    }
    for (int i = i1;i < i2 + 1;i++){
        for (int j = j1;j < j2 + 1;j++){
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }
}

void poiseq(double u[mx][my], double v[mx][my], double p[mx][my], double dx, double dy, double dt){
    double ux;
    double uy;
    double vx;
    double vy;
    /*double ** rhs = new double*[mx];
    for( int i=0;i < mx; i++ ) {
        rhs[i] = new double[my];
    }*/
    double rhs[mx][my];
    double dp;
    double resp;
    int itrp;
    for (int i = 1;i < mx - 1;i++){
        for (int j = 1;j < my - 1;j++){
            if ((i > i1 && i < i2) && (j > j1 && j < j2)){
            }else{
                ux = (u[i + 1][j] - u[i - 1][j]) / (2.0 * dx);
                uy = (u[i][j + 1] - u[i][j - 1]) / (2.0 * dy);
                vx = (v[i + 1][j] - v[i - 1][j]) / (2.0 * dx);
                vy = (v[i][j + 1] - v[i][j - 1]) / (2.0 * dy);
                rhs[i][j] = (ux + vy) / dt - (ux * ux + 2.0 * ux * vy + vy * vy);
            }
        }
    }
    for (int itr = 0;itr < maxitp;itr++){
        double res = 0.0;
        for (int i = 1;i < mx - 1;i++){
            for (int j = 1;j < my - 1;j++){
                if ((i > i1 - 1 && i < i2 + 1) && (j > j1 - 1 && j < j2 + 1)){
                }else{
                    dp = (p[i + 1][j] + p[i - 1][j]) / pow(dx,2.0) + (p[i][j + 1] + p[i][j - 1]) / pow(dy,2.0) - rhs[i][j];
                    dp = dp / (2.0 / pow(dx, 2.0) + 2.0 / pow(dy, 2.0)) - p[i][j];
                    res += pow(dp,2.0);
                    p[i][j] += omegap * dp;
                }
            }
        }
        bcforp(p);
        res = sqrt(res / double(mx * my));
        if(res < errorp){
            resp = res;
            itrp = itr;
            break;
        }
    }
    /*for( int i=0; i<mx; i++ ) {
        delete[] rhs[i];
    }
    delete[] rhs;*/
}

void veloeq(double u[mx][my], double v[mx][my], double p[mx][my], double dx, double dy, double dt){
    double urhs[mx][my];
    double vrhs[mx][my];
    for (int i = 1;i < mx - 1;i++){
        for (int j = 1;j < my - 1;j++){
            if ((i > i1 && i < i2) && (j > j1 && j < j2)){
            }else if((i == i1 || i == i2) && (j > j1 && j < j2)){
                urhs[i][j] = 0.0;
                vrhs[i][j] = -(p[i][j + 1] - p[i][j - 1]) / 2.0 / dy;
            }else if((i > i1 && i < i2) && (j == j1 || j == j2)){
                urhs[i][j] = -(p[i + 1][j] - p[i - 1][j]) / 2.0 / dx;
                vrhs[i][j] = 0.0;
            }else if((i == i1 || i == i2) && (j == j1 || j == j2)){
                urhs[i][j] = 0.0;
                vrhs[i][j] = 0.0;
            }else{
                urhs[i][j] = -(p[i + 1][j] - p[i - 1][j]) / 2.0 / dx;
                vrhs[i][j] = -(p[i][j + 1] - p[i][j - 1]) / 2.0 / dy;
            }
        }
    }

    for (int i = 1;i < mx - 1;i++){
        for (int j = 1;j < my - 1;j++){
            if ((i > i1 && i < i2) && (j > j1 && j < j2)){
            }else{
                urhs[i][j] += (u[i + 1][j] - 2.0 * u[i][j] + u[i -1][j]) / (Re * pow(dx,2.0)) + (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) / (Re * pow(dy,2.0));
                vrhs[i][j] += (v[i + 1][j] - 2.0 * v[i][j] + v[i -1][j]) / (Re * pow(dx,2.0)) + (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]) / (Re * pow(dy,2.0));
            }
        }
    }

    for (int j = j1 + 1;j < j2;j++){
        u[i1 + 1][j] = 2.0 * u[i1][j] - u[i1 - 1][j];
        u[i2 - 1][j] = 2.0 * u[i2][j] - u[i2 + 1][j];
        v[i1 + 1][j] = 2.0 * v[i1][j] - v[i1 - 1][j];
        v[i2 - 1][j] = 2.0 * v[i2][j] - v[i2 + 1][j];
    }
    for (int i = 1;i < mx - 1;i++){
        for (int j = 1;j < my - 1;j++){
            if ((i > i1 && i < i2) && (j > j1 && j < j2)){
            }else if(i == 1){
                urhs[i][j] -= u[i][j] * (-u[i+2][j] + 8.0 * (u[i+1][j] - u[i-1][j]) + 1.0) / 12.0 / dx + abs(u[i][j]) * (u[i+2][j] - 4.0 * u[i+1][j] + 6.0 * u[i][j] - 4.0 * u[i-1][j] + 1.0) / 4.0 / dx;
                vrhs[i][j] -= u[i][j] * (-v[i+2][j] + 8.0 * (v[i+1][j] - v[i-1][j]) + 0.0) / 12.0 / dx + abs(u[i][j]) * (v[i+2][j] - 4.0 * v[i+1][j] + 6.0 * v[i][j] - 4.0 * v[i-1][j] + 0.0) / 4.0 / dx;
            }else if(i == mx - 2){
                urhs[i][j] -= u[i][j] * (-2.0 * u[i+1][j] + u[i][j] + 8 * (u[i+1][j] - u[i-1][j]) + u[i-2][j]) / 12.0 / dx + abs(u[i][j]) * (2.0 * u[i+1][j] - u[i][j] - 4.0 * u[i+1][j] + 6.0 * u[i][j] - 4.0 * u[i-1][j] + u[i-2][j]) / 4.0 / dx;
                vrhs[i][j] -= u[i][j] * (-2.0 * v[i+1][j] + v[i][j] + 8 * (v[i+1][j] - v[i-1][j]) + v[i-2][j]) / 12.0 / dx + abs(u[i][j]) * (2.0 * v[i+1][j] - v[i][j] - 4.0 * v[i+1][j] + 6.0 * v[i][j] - 4.0 * v[i-1][j] + v[i-2][j]) / 4.0 / dx;
            }else{
                urhs[i][j] -= u[i][j] * (-u[i+2][j] + 8.0 * (u[i+1][j] - u[i-1][j]) + u[i-2][j]) / 12.0 / dx + abs(u[i][j]) * (u[i+2][j] - 4.0 * u[i+1][j] + 6.0 * u[i][j] - 4.0 * u[i-1][j] + u[i-2][j]) / 4.0 / dx;
                vrhs[i][j] -= u[i][j] * (-v[i+2][j] + 8.0 * (v[i+1][j] - v[i-1][j]) + v[i-2][j]) / 12.0 / dx + abs(u[i][j]) * (v[i+2][j] - 4.0 * v[i+1][j] + 6.0 * v[i][j] - 4.0 * v[i-1][j] + v[i-2][j]) / 4.0 / dx;
            }
        }
    }

    for (int i = i1 + 1;i < i2;i++){
        u[i][j1 + 1] = 2.0 * u[i][j1] - u[i][j1 - 1];
        u[i][j2 - 1] = 2.0 * u[i][j2] - u[i][j2 + 1];
        v[i][j1 + 1] = 2.0 * v[i][j1] - v[i][j1 - 1];
        v[i][j2 - 1] = 2.0 * v[i][j2] - v[i][j2 + 1];
    }
    for (int i = 1;i < mx - 1;i++){
        for (int j = 1;j < my - 1;j++){
            if ((i > i1 && i < i2) && (j > j1 && j < j2)){
            }else if(j == 1){
                urhs[i][j] -= v[i][j] * (-u[i][j+2] + 8.0 * (u[i][j+1] - u[i][j-1]) + 2.0 * u[i][j-1] - u[i][j]) / 12.0 / dy + abs(v[i][j]) * (u[i][j+2] - 4.0 * u[i][j+1] + 6.0 * u[i][j] - 4.0 * u[i][j-1] + 2.0 * u[i][j-1] - u[i][j]) / 4.0 / dy;
                vrhs[i][j] -= v[i][j] * (-v[i][j+2] + 8.0 * (v[i][j+1] - v[i][j-1]) + 2.0 * v[i][j-1] - v[i][j]) / 12.0 / dy + abs(v[i][j]) * (v[i][j+2] - 4.0 * v[i][j+1] + 6.0 * v[i][j] - 4.0 * v[i][j-1] + 2.0 * v[i][j-1] - v[i][j]) / 4.0 / dy;
            }else if(j == my - 2){
                urhs[i][j] -= v[i][j] * (-2.0 * u[i][j+1] + u[i][j] + 8.0 * (u[i][j+1] - u[i][j-1]) + u[i][j-2]) / 12.0 / dy + abs(v[i][j]) * (2.0 * u[i][j+1] - u[i][j] - 4.0 * u[i][j+1] + 6.0 * u[i][j] - 4.0 * u[i][j-1] + u[i][j-2]) / 4.0 / dy;
                vrhs[i][j] -= v[i][j] * (-2.0 * v[i][j+1] + v[i][j] + 8.0 * (v[i][j+1] - v[i][j-1]) + v[i][j-2]) / 12.0 / dy + abs(v[i][j]) * (2.0 * v[i][j+1] - v[i][j] - 4.0 * v[i][j+1] + 6.0 * v[i][j] - 4.0 * v[i][j-1] + v[i][j-2]) / 4.0 / dy;
            }else{
                urhs[i][j] -= v[i][j] * (-u[i][j+2] + 8.0 * (u[i][j+1] - u[i][j-1]) + u[i][j-2]) / 12.0 / dy + abs(v[i][j]) * (u[i][j+2] - 4.0 * u[i][j+1] + 6.0 * u[i][j] - 4.0 * u[i][j-1] + u[i][j-2]) / 4.0 / dy;
                vrhs[i][j] -= v[i][j] * (-v[i][j+2] + 8.0 * (v[i][j+1] - v[i][j-1]) + v[i][j-2]) / 12.0 / dy + abs(v[i][j]) * (v[i][j+2] - 4.0 * v[i][j+1] + 6.0 * v[i][j] - 4.0 * v[i][j-1] + v[i][j-2]) / 4.0 / dy;
            }
        }
    }


    for (int i = 1;i < mx - 1;i++){
        for (int j = 1;j < my - 1;j++){
            if ((i > i1 && i < i2) && (j > j1 && j < j2)){
            }else{
                u[i][j] += dt * urhs[i][j];
                v[i][j] += dt * vrhs[i][j];
            }
        }
    }
}

int main(){
    double CFL = 0.2;//クーラン数
    int nlast = 5000;
    int nlp = 10;

    double dx = 1.0 / double(i2 - i1);
    double dy = 1.0 / double(j2 - j1);
    int icent = (i1 + i2) / 2;
    int jcent = (j1 + j2) / 2;
    double x[mx][my];
    double y[mx][my];
    for (int i = 0;i < mx;i++){
        for (int j = 0;j < my;j++){
            x[i][j] = dx * double(i - icent);
            y[i][j] = dy * double(j - jcent);
        }
    }

    double dt = CFL * fmin(dx,dy);
    int nbegin = 0;
    double time = 0.0;
    int nstep;
    double u[mx][my];
    double v[mx][my];
    double p[mx][my];
    for (int i = 0;i < mx;i++){
        for (int j = 0;j < my;j++){
            u[i][j] = 1.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;
        }
    }

    bcforp(p);
    bcforv(u,v);
    double cpfore;
    double cpback;
    double cpbtm;
    double cptop;
    double cd;
    double cl;
    double CD[nlast];
    double CL[nlast];
    for (int n = 0;n < nlast;n++){
        nstep = n + nbegin;
        //if(nstep%100 == 0){
            cout << nstep << endl;
        //}
        time = time + dt;
        poiseq(u,v,p,dx,dy,dt);
        bcforp(p);
        veloeq(u,v,p,dx,dy,dt);
        bcforv(u,v);
        cd = 0;
        cl = 0;
        for (int j = j1; j < j2;j++){
            cpfore = (2.0 * p[i1][j] + 2.0 * p[i1][j + 1])/2.0;
            cpback = (2.0 * p[i2][j] + 2.0 * p[i2][j + 1])/2.0;
            cd += (cpfore - cpback) * dy;
        }
        for (int i = i1; i < i2;i++){
            cpbtm = (2.0 * p[i][j1] + 2.0 * p[i + 1][j1])/2.0;
            cptop = (2.0 * p[i][j1] + 2.0 * p[i + 1][j2])/2.0;
            cl += (cpbtm - cptop) * dx;
        }
        CD[n] = cd;
        CL[n] = cl;
    }
    ofstream logu;//ファイル書き込み用
    ofstream logv;//ファイル書き込み用
    ofstream logp;//ファイル書き込み用
	logu.open("logu.csv",ios::trunc);
    logv.open("logv.csv",ios::trunc);
    logp.open("logp.csv",ios::trunc);
    /*for (int i = 0;i < mx;i++){
        for (int j = 0;j < my - 1;j++){
            logu << u[i][j] << ",";
            logv << v[i][j] << ",";
            logp << p[i][j] << ",";
        }
        logu << u[i][my] << endl;
        logv << v[i][my] << endl;
        logp << p[i][my] << endl;
    }*/
    for (int j = 0;j < my;j++){
        for (int i = 0;i < mx - 1;i++){
            logu << u[i][j] << ",";
            logv << v[i][j] << ",";
            logp << p[i][j] << ",";
        }
        logu << u[mx - 1][j] << endl;
        logv << v[mx - 1][j] << endl;
        logp << p[mx - 1][j] << endl;
    }
}
