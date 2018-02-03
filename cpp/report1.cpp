#include <iostream>
#include<fstream>
#include <cmath>
#include <vector>
using namespace std;

int main(){
	double rou0 = 1.225;//密度
	double t0 = 293.0;//温度
	double visc = 1.458 * pow(10, -6) * pow(t0, 1.5) / (t0 + 110.4);//粘性係数
	double width = 0.5;//代表長さ
	double v0 = 20.0;//速度
	double re0 = rou0 * v0 * width / visc;//後縁でのレイノルズ数
	double blt = width * 5.2 / sqrt(re0);//ブラジウス解
	int nx = 50000;//xの分割数
	double dx = width / nx;//xの幅
	int ny = 200;//yの分割数
	double ymax = blt * 2.0;//yの長さ
	double dy = ymax / ny;//yの幅
	vector<double> y(ny + 1);//yの値
	vector<double> u(ny + 1, v0);//uの値
	vector<double> v(ny + 1, 0);//vの値
	vector<double> u_pre(ny + 1, v0);//uの値(一つ前)
	vector<double> v_pre(ny + 1, 0);//vの値(一つ前)
	double dudy1;//uのy偏微分(運動量の式)
	double dudy2;//uのy二階偏微分(運動量の式)
	double dudx;//uのx偏微分(運動量の式)
	double dudx1;//uのx偏微分(連続の式)
	double dudx2;//uのx偏微分(連続の式)
	ofstream u_log;//ファイル書き込み用
	u_log.open("u.csv",ios::trunc);
	ofstream v_log;//ファイル書き込み用
	v_log.open("v.csv",ios::trunc);
	ofstream delta_log;//ファイル書き込み用
	delta_log.open("delta.csv",ios::trunc);
	bool flag = 0;//境界層内かどうかの判定
	double delta_star;//排除厚さ
	ofstream delta_star_log;//ファイル書き込み用
	delta_star_log.open("delta_star.csv",ios::trunc);
	double theta;//運動量厚さ
	ofstream theta_log;//ファイル書き込み用
	theta_log.open("theta.csv",ios::trunc);
	double theta_e;//エネルギー厚さ
	ofstream theta_e_log;//ファイル書き込み用
	theta_e_log.open("theta_e.csv",ios::trunc);
	double cf;//壁面摩擦係数
	ofstream cf_log;//ファイル書き込み用
	cf_log.open("cf.csv",ios::trunc);
	for(int j = 0;j < ny + 1;j++){//yの設定
		y[j] = dy * j;
	}
	for(int i = 1;i < nx + 1;i++){
		u[0] = 0.0;
		u[ny] = v0;
		v[0] = 0.0;
		flag = 0;
		delta_star = 0;
		theta = 0;
		theta_e = 0;
		for(int j = 1;j < ny;j++){//uの計算
			dudy1 = (u_pre[j + 1] - u_pre[j - 1]) / (2.0 * dy);
			dudy2 = (u_pre[j + 1] - 2.0 * u_pre[j] + u_pre[j - 1]) / pow(dy, 2.0);
			dudx = (dudy2 * visc / rou0 - v_pre[j] * dudy1) / u_pre[j];
			u[j] = u_pre[j] + dx * dudx;
		}
		for(int j = 1;j < ny + 1;j++){//vの計算
			dudx1 = (u[j] - u_pre[j]) / dx;
			dudx2 = (u[j + 1] - u_pre[j + 1]) / dx;
			v[j] = v[j - 1] - (dudx1 + dudx2) * dy / 2.0;
		}
		for(int j = 0;j < ny + 1;j++){//境界層厚さなどの計算
			u_pre[j] = u[j];
			v_pre[j] = v[j];
			if(i == nx){
				u_log << y[j] << "," << u[j] << endl;//uのログ書き込み
				v_log << y[j] << "," << v[j] << endl;//vのログ書き込み		
			}
			if(u[j] > v0 * 0.995 && flag == 0){//境界層厚さのログ書き込み
				if(j == 0){
					delta_log << width * i / nx << ","  << y[0] << "," << width * i / nx * 5.2 / sqrt(rou0 * v0 * width * i / nx / visc) << endl;
				}else{
					delta_log << width * i / nx << ","  << y[j - 1] << "," << width * i / nx * 5.2 / sqrt(rou0 * v0 * width * i / nx / visc) << endl;
				}
				flag = 1;
			}
			if(j == ny){
			}else{
				delta_star += (1 - (u[j] + u[j + 1]) / (2 * v0)) * dy;//排除厚さ台形近似積分
				theta += (u[j] / v0 * (1 - u[j] / v0) + u[j + 1] / v0 * (1 - u[j + 1] / v0)) / 2 * dy;//運動量厚さ台形近似積分
				theta_e += (u[j] / v0 * (1 - pow(u[j] / v0, 2.0)) + u[j + 1] / v0 * (1 - pow(u[j + 1] / v0, 2.0))) / 2 * dy;//エネルギー厚さ台形近似積分
			}
		}
		delta_star_log << width * i / nx << "," << delta_star << "," << width * i / nx * 1.72 / sqrt(rou0 * v0 * width * i / nx / visc) << endl;//排除厚さログ書き込み
		theta_log << width * i / nx << "," << theta << "," << width * i / nx * 0.664 / sqrt(rou0 * v0 * width * i / nx / visc) << endl;//運動量厚さのログ書き込み
		theta_e_log << width * i / nx << "," << theta_e << endl;//エネルギー厚さのログ書き込み
		cf_log << width * i / nx << "," << visc * u[1] / dy / (0.5 * rou0 * pow(v0, 2.0)) << "," << 0.664 / sqrt(rou0 * v0 * width * i / nx / visc) << endl;//壁面摩擦係数のログ書き込み
	}
	u_log.close();
	v_log.close();
	delta_log.close();
	delta_star_log.close();
	theta_log.close();
	cf_log.close();
}