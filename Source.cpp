#include <sstream>
#include <iostream>
#include<limits.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include<vector>
//#include<unistd.h>
#include <time.h>
#include <algorithm>
//#include <windows.h>
#include <numeric>
using namespace std;
#define type 1
#define Edgeamount 1824
#define L_bits_smallest 96
#define sigma_p 0.3
int const  n = 576, k = 288;
int const H_row = 288, H_col = 576;
int crc_length = 0;
int switch_case = 0;
int Rmax = 0;
int Vmax = 0;
//vector<int> un_saa(H_row, 0);

vector<int> check_location[288];
vector<int> var_location[576];
vector<double> old_R[288];
vector<double> residu[288];
vector<double> new_R[288];
double static sigma;
int Chk_Eq = false;
vector<int> info;
vector<double> errV;

double sgn(double a)
{
	int s;
	if (a > 0)
	{
		s = 1;
	}
	else if (a < 0)
	{
		s = -1;
	}
	else
	{
		s = 1;
	}
	return s;
}
float rand49()
{
	/*rand_max=7FFF (32767) */
	static int Num = 0;
	double number;
	int    i;
	i = rand();
	number = (double)(i) / ((unsigned)(RAND_MAX + 1));
	Num++;
	if (Num >= RAND_MAX)
	{
		time_t t;
		t = time(NULL);
		//		srand((unsigned)(t%RAND_MAX));
		Num = 0;
	}
	return (float)number;
}

double Normal()
{
	static int iset = 0;
	static double qset;
	double vx, vy, r, temp;
	if (iset == 0)//noise=normal*deviate
	{
		do
		{
			vx = 2.0 * rand49() - 1.0;
			vy = 2.0 * rand49() - 1.0;
			r = vx * vx + vy * vy;
		} while (r >= 1.0 || r == 0);
		temp = sqrt(-2.0 * log(r) / r);
		qset = vy * temp;
		iset = 1;
		return (vx * temp);
	}
	else
	{
		iset = 0;
		return qset;
	}
}

double ptanh(double x)
{
	if (x <= -7.0)
		return -0.999998;
	else if (-7.0 < x && x <= -3.68)
		return 0.0004 * x - 0.9972;
	else if (-3.68 < x && x <= -1.82)
		return 0.0268 * x - 0.90;
	else if (-1.82 < x && x <= -1.24)
		return 0.1781 * x - 0.6247;
	else if (-1.24 < x && x <= -0.66)
		return 0.4605 * x - 0.2745;
	else if (-0.66 < x && x <= 0.66)
		return 0.8764 * x;
	else if (0.66 < x && x <= 1.24)
		return 0.4605 * x + 0.2745;
	else if (1.24 < x && x <= 1.82)
		return 0.1781 * x + 0.6247;
	else if (1.82 < x && x <= 3.68)
		return (0.0268 * x + 0.9);
	else if (3.68 < x && x <= 7.0)
		return 0.0004 * x + 0.9972;
	else
		return 0.999998;
}

double patanh(double x)
{
	if (-1 < x && x <= -0.999998)
		return -7.0;
	else if (-0.999998 < x && x <= -0.9987)
		return ((x + 0.9972) / 0.0004);
	else if ((-0.9987 < x && x <= -0.9488))
		return ((x + 0.9) / 0.0268);
	else if ((-0.9488 < x && x <= -0.8455))
		return ((x + 0.6247) / 0.1781);
	else if ((-0.8455 < x && x <= -0.5784))
		return ((x + 0.2745) / 0.4605);
	else if ((-0.5784 < x && x <= 0.5784))
		return (x / 0.8764);
	else if ((0.5784 < x && x <= 0.8455))
		return  ((x - 0.2745) / 0.4605);
	else if ((0.8455 < x && x <= 0.9488))
		return ((x - 0.6247) / 0.1781);
	else if ((0.9488 < x && x <= 0.9987))
		return ((x - 0.9) / 0.0268);
	else if ((0.9987 < x && x <= 0.999998))
		return  ((x - 0.9972) / 0.0004);
	else if ((x > 0.999998))
		return 7.0;
}

double gfunction(double L1, double L2)
{
	{
		return  0.9375 * sgn(L1) * sgn(L2) * min(fabs(L1), fabs(L2));
	}
}
double phi_func(double a)
{
	double result;
	double absolute;
	double tmp;
	absolute = fabs(a);
	if (absolute < pow(10, -15))
		result = -37.42995;
	else if (absolute > 37.43)
		result = -1.110223e-15;
	else
	{
		tmp = (exp(a) - 1) / (exp(a) + 1);
		tmp = fabs(tmp);
		result = log(tmp);
	}
	return result;
}

void read_BP_location()
{
	vector <int> row_reg[H_row];
	vector <int> col_reg[H_col];

	fstream file3, file4;
	file3.open("Check_location.txt", ios::in);
	file4.open("Variable_location.txt", ios::in);
	if (!file3 || !file4)     //檢查檔案是否成功開啟
	{
		cerr << "Can't open file!\n";
		exit(1);     //在不正常情形下，中斷程式的執行
	}
	for (int i = 0; i < H_row; i++)
	{
		row_reg[i].resize(8);
		for (int j = 0; j < 8; j++)
		{
			file3 >> row_reg[i].at(j);
		}
	}
	for (int i = 0; i < H_col; i++)
	{
		col_reg[i].resize(7);
		for (int j = 0; j < 7; j++)
		{
			file4 >> col_reg[i].at(j);
		}
	}
	for (int i = 0; i < H_row; i++)
	{
		for (int j = 0; j < row_reg[i].size(); j++)
		{
			if (row_reg[i].at(j) == -1)
			{
				break;
			}
			else
			{
				check_location[i].push_back(row_reg[i][j]);
			}
		}
	}
	for (int i = 0; i < H_col; i++)
	{
		for (int j = 0; j < col_reg[i].size(); j++)
		{
			if (col_reg[i].at(j) == -1)
			{
				break;
			}
			else
			{
				var_location[i].push_back(col_reg[i][j]);
			}
		}
	}

	file3.close();
	file4.close();
}


int Parity_Check(vector<double> check_node[], vector<double> channel_LLR)
{
	vector<int> early_decode_reg(H_col, 0);
	vector<bool> early_check(H_row, 0);
	vector<double> early_check_var(H_col, 0);

	for (int i = 0; i < H_col; i++)
	{
		double sum_var = 0;
		for (int j = 0; j < var_location[i].size(); j++)
		{
			sum_var = sum_var + check_node[var_location[i].at(j)].at(i);
		}
		early_check_var.at(i) = sum_var + channel_LLR.at(i);
		/*	if (P_sign_V.at(i) != sgn(early_check_var.at(i)) && P_sign_V.at(i) != 0) {
				OC_V.at(i)++;
				stable_Count.at(i) == 0;
			}
			else if (P_sign_V.at(i) != 0) {
				stable_Count.at(i)++;
			}
			P_sign_V.at(i) = sgn(early_check_var.at(i));*/
	}

	for (int i = 0; i < H_col; i++)
	{
		if (early_check_var.at(i) >= 0)
		{
			early_decode_reg.at(i) = 0;
		}
		else
		{
			early_decode_reg.at(i) = 1;
		}
	}

	for (int i = 0; i < H_row; i++)
	{
		bool early = 0;
		for (int j = 0; j < check_location[i].size(); j++)
		{
			early = early ^ early_decode_reg.at(check_location[i].at(j));
		}
		if (early != 0)
		{
			return 0;
		}
	}

	return 1;

}

int Residual1(vector<double> variable_node[])
{
	int i_check = 0;
	int min_check = 0;
	vector<double> rresidu(H_row, 0);
	vector<vector<double> > SEE(H_row, vector<double>(H_col));

	for (int i = 0; i < H_row; i++)
	{
		i_check = i;
		double check_scr = 0;
		double s_scr = -1.0;
		new_R[i].clear();
		old_R[i].clear();
		residu[i].clear();

		for (int j = 0; j < check_location[i_check].size(); j++)
		{
			check_scr = check_scr + phi_func((variable_node[i_check].at(check_location[i_check].at(j))));
			//            check_scr = check_scr * ptanh(0.5*(rvariable_node[i_check].at(check_location[i_check].at(j))));
			//            check_scr = 2 * patanh(check_scr);
			s_scr = s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(j)));
		}

		for (int j = 0; j < check_location[i_check].size(); j++)
		{
			new_R[i_check].push_back((s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(j)))) * phi_func(check_scr - phi_func((variable_node[i_check].at(check_location[i_check].at(j))))));
			//cout << new_R[i_check].back() << " ";
		}
		SEE[i_check] = new_R[i_check];

	}


	for (int i = 0; i < H_row; i++)
	{
		for (int j = 0; j < new_R[i].size(); j++)
		{
			residu[i].push_back(new_R[i].at(j));
		}
	}


	for (int i = 0; i < H_row; i++)
	{
		double res = 0;
		for (int j = 0; j < new_R[i].size(); j++)
		{
			res = max(res, fabs(new_R[i].at(j)));
		}
		res = fabs(res);

		rresidu.at(i) = res;
	}


	for (int i = 0; i < H_row; i++)
	{
		for (int j = 0; j < new_R[i].size(); j++)
		{
			old_R[i].push_back(new_R[i].at(j));
		}
	}

	//    for(int i = 0; i<H_row;i++)
	//    {
	//        cout<<rresidu.at(i)<<endl;
	//    }
	//    system("pause");




	min_check = distance(begin(rresidu), max_element(begin(rresidu), end(rresidu)));
	vector<double> max_index;
	for (int i = 0; i < residu[min_check].size(); i++)
	{
		max_index.push_back(fabs(residu[min_check].at(i)));
	}

	Rmax = distance(begin(max_index), max_element(begin(max_index), end(max_index)));
	Vmax = check_location[min_check].at(Rmax);

	max_index.clear();

	return min_check;
}

int rand(int x)
{
	return rand() % x;
}

vector<int> FindMaxResdual(vector<double> residu, int n) {
	vector<int> output;
	vector<int> Chk;

	for (int i = 0; i < H_row; i++) {
		Chk.push_back(i);
	}
	for (int i = 0; i < n; i++) {
		//ChkandValue max_index = *max_element(begin(residu), end(residu));
		int current_max = distance(begin(residu), max_element(begin(residu), end(residu)));
		if (residu.at(current_max) == -DBL_MAX) {
			break;
		}
		output.push_back(Chk.at(current_max));
		residu.erase(residu.begin() + current_max);
		Chk.erase(Chk.begin() + current_max);
	}
	return output;
}

int decision(vector<double> final_var) {
	int err = 0;
	vector<int> x_tilde;
	for (int i = 0; i < final_var.size(); i++) {
		x_tilde.push_back(final_var.at(i) > 0 ? 0 : 1);
	}
	for (int i = 0; i < k; i++) {
		err += x_tilde.at(i) == info.at(i) ? 0 : 1;
	}
	return err;
}

vector<double>  BP_decode_flooding(vector<double>receive_bits, int itermax)
{
	int i_check = 0;
	int pre;
	int P_check = 0;
	int init = H_col - n;
	int ko = 0;
	int stop = 0;
	vector<int> decoded_bits(k, 0);
	vector<int> decode_reg(n, 0);

	vector<int> Vk;
	vector<int> Ca;
	vector<double> max_index;
	vector<double> channel_LLR(n, 0);
	vector<double> check_node[H_row];
	vector<double> variable_node[H_row];
	vector<double> check_reg(H_col, 0);
	vector<double> var_reg(H_col, 0);
	vector<double> final_var(n, 0);


	for (int i = 0; i < n; i++)
	{
		channel_LLR.at(i) = (-2.0 * receive_bits.at(i)) / (sigma * sigma);
		//cout << channel_LLR.at(i) << " ";
	}


	for (int i = 0; i < H_row; i++)
	{
		check_node[i].resize(H_col);
		variable_node[i].resize(H_col);
	}


	for (int i = 0; i < H_col; i++)
	{
		for (int j = 0; j < var_location[i].size(); j++)
		{
			variable_node[var_location[i].at(j)].at(i) = channel_LLR.at(i);
		}
	}

	for (int iter = 0; iter < itermax; iter++)
	{
		for (int i_r = 0; i_r < H_row; i_r++)
		{
			double check_scr = 0;
			double s_scr = -1.0;
			i_check = i_r;
			Vk.clear();
			//cout << i_check << " ";
			for (int j = 0; j < check_location[i_check].size(); j++)
			{
				check_scr = check_scr + phi_func((variable_node[i_check].at(check_location[i_check].at(j))));
				s_scr = s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(j)));
				Vk.push_back(check_location[i_check].at(j));
			}

			for (int j = 0; j < check_location[i_check].size(); j++)
			{
				check_node[i_check].at(check_location[i_check].at(j)) = (s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(j)))) * phi_func(check_scr - phi_func((variable_node[i_check].at(check_location[i_check].at(j)))));
				//cout << "CN : " << i_check << " to VN :" << check_location[i_check].at(j) <<" value : "<< check_node[i_check].at(check_location[i_check].at(j)) << endl;
			}


		}
		for (int i_c = 0; i_c < H_col; i_c++) {
			//for (int cc = 0; cc < Vk.size(); cc++)
			{
				double reg = 0.0;
				for (int j = 0; j < var_location[i_c].size(); j++)
				{
					reg = reg + check_node[var_location[i_c].at(j)].at(i_c);
				}

				for (int j = 0; j < var_location[i_c].size(); j++)
				{
					variable_node[var_location[i_c].at(j)].at(i_c) = reg - check_node[var_location[i_c].at(j)].at(i_c) + channel_LLR.at(i_c);
					//cout << "VN : " << Vk.at(cc) << " to CN :" << var_location[Vk.at(cc)].at(j) <<" value : "<< variable_node[var_location[Vk.at(cc)].at(j)].at(Vk.at(cc)) << endl;
				}
			}
		}

		Chk_Eq = Parity_Check(check_node, channel_LLR);
		if (Chk_Eq == true)
		{
			break;
		}

	}


	vector<double> early_check_var(H_col, 0);

	for (int i = 0; i < H_col; i++)
	{
		double sum_var = 0;
		for (int j = 0; j < var_location[i].size(); j++)
		{
			sum_var = sum_var + check_node[var_location[i].at(j)].at(i);
		}
		early_check_var.at(i) = sum_var + channel_LLR.at(i);
		//cout << early_check_var.at(i) << " ";
	}

	return early_check_var;

}

vector<double>  BP_decode_LBP(vector<double>receive_bits, int itermax)
{
	int i_check = 0;
	int pre;
	int P_check = 0;
	int init = H_col - n;
	int ko = 0;
	int stop = 0;
	vector<int> decoded_bits(k, 0);
	vector<int> decode_reg(n, 0);

	vector<int> Vk;
	vector<int> Ca;
	vector<double> max_index;
	vector<double> channel_LLR(n, 0);
	vector<double> check_node[H_row];
	vector<double> variable_node[H_row];
	vector<double> check_reg(H_col, 0);
	vector<double> var_reg(H_col, 0);
	vector<double> final_var(n, 0);


	for (int i = 0; i < n; i++)
	{
		channel_LLR.at(i) = (-2.0 * receive_bits.at(i)) / (sigma * sigma);
		//cout << channel_LLR.at(i) << " ";
	}


	for (int i = 0; i < H_row; i++)
	{
		check_node[i].resize(H_col);
		variable_node[i].resize(H_col);
	}


	for (int i = 0; i < H_col; i++)
	{
		for (int j = 0; j < var_location[i].size(); j++)
		{
			variable_node[var_location[i].at(j)].at(i) = channel_LLR.at(i);
		}
	}

	for (int iter = 0; iter < itermax; iter++)
	{
		for (int i_r = 0; i_r < H_row; i_r++)
		{
			double check_scr = 0;
			double s_scr = -1.0;
			i_check = i_r;
			Vk.clear();
			//cout << i_check << " ";
			for (int j = 0; j < check_location[i_check].size(); j++)
			{
				check_scr = check_scr + phi_func((variable_node[i_check].at(check_location[i_check].at(j))));
				s_scr = s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(j)));
				Vk.push_back(check_location[i_check].at(j));
			}

			for (int j = 0; j < check_location[i_check].size(); j++)
			{
				check_node[i_check].at(check_location[i_check].at(j)) = (s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(j)))) * phi_func(check_scr - phi_func((variable_node[i_check].at(check_location[i_check].at(j)))));
				//cout << "CN : " << i_check << " to VN :" << check_location[i_check].at(j) <<" value : "<< check_node[i_check].at(check_location[i_check].at(j)) << endl;
			}


			for (int cc = 0; cc < Vk.size(); cc++)
			{

				double reg = 0.0;
				for (int j = 0; j < var_location[Vk.at(cc)].size(); j++)
				{
					reg = reg + check_node[var_location[Vk.at(cc)].at(j)].at(Vk.at(cc));
				}

				for (int j = 0; j < var_location[Vk.at(cc)].size(); j++)
				{
					variable_node[var_location[Vk.at(cc)].at(j)].at(Vk.at(cc)) = reg - check_node[var_location[Vk.at(cc)].at(j)].at(Vk.at(cc)) + channel_LLR.at(Vk.at(cc));
					//cout << "VN : " << Vk.at(cc) << " to CN :" << var_location[Vk.at(cc)].at(j) <<" value : "<< variable_node[var_location[Vk.at(cc)].at(j)].at(Vk.at(cc)) << endl;
				}
			}
		}

		Chk_Eq = Parity_Check(check_node, channel_LLR);
		if (Chk_Eq == true)
		{
			break;
		}

	}


	vector<double> early_check_var(H_col, 0);

	for (int i = 0; i < H_col; i++)
	{
		double sum_var = 0;
		for (int j = 0; j < var_location[i].size(); j++)
		{
			sum_var = sum_var + check_node[var_location[i].at(j)].at(i);
		}
		early_check_var.at(i) = sum_var + channel_LLR.at(i);
		//cout << early_check_var.at(i) << " ";
	}

	return early_check_var;

}

vector<double>  BP_decode_LBP_O(vector<double>receive_bits, int itermax)
{
	int i_check = 0;
	int pre;
	int P_check = 0;
	int init = H_col - n;
	int ko = 0;
	int stop = 0;
	vector<int> decoded_bits(k, 0);
	vector<int> decode_reg(n, 0);

	vector<int> Vk;
	vector<int> Ca;
	vector<double> max_index;
	vector<double> channel_LLR(n, 0);
	vector<double> check_node[H_row];
	vector<double> variable_node[H_row];
	vector<double> check_reg(H_col, 0);
	vector<double> var_reg(H_col, 0);
	vector<double> final_var(n, 0);


	for (int i = 0; i < n; i++)
	{
		channel_LLR.at(i) = (-2.0 * receive_bits.at(i)) / (sigma * sigma);
		//cout << channel_LLR.at(i) << " ";
	}


	for (int i = 0; i < H_row; i++)
	{
		check_node[i].resize(H_col);
		variable_node[i].resize(H_col);
	}


	for (int i = 0; i < H_col; i++)
	{
		for (int j = 0; j < var_location[i].size(); j++)
		{
			variable_node[var_location[i].at(j)].at(i) = channel_LLR.at(i);
		}
	}

	for (int iter = 0; iter < itermax; iter++)
	{
		for (int i_r = 0; i_r < H_row; i_r++)
		{
			double check_scr = 0;
			double s_scr = -1.0;
			i_check = i_r;
			Vk.clear();
			//cout << i_check << " ";
			for (int j = 0; j < check_location[i_check].size(); j++)
			{
				check_scr = check_scr + phi_func((variable_node[i_check].at(check_location[i_check].at(j))));
				s_scr = s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(j)));
				Vk.push_back(check_location[i_check].at(j));
			}

			for (int j = 0; j < check_location[i_check].size(); j++)
			{
				check_node[i_check].at(check_location[i_check].at(j)) = (s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(j)))) * phi_func(check_scr - phi_func((variable_node[i_check].at(check_location[i_check].at(j)))));
				//cout << "CN : " << i_check << " to VN :" << check_location[i_check].at(j) << " value : " << check_node[i_check].at(check_location[i_check].at(j)) << endl;
			}


			for (int cc = 0; cc < Vk.size(); cc++)
			{

				double reg = 0.0;
				for (int j = 0; j < var_location[Vk.at(cc)].size(); j++)
				{
					reg = reg + check_node[var_location[Vk.at(cc)].at(j)].at(Vk.at(cc));
				}

				for (int j = 0; j < var_location[Vk.at(cc)].size(); j++)
				{
					if (var_location[Vk.at(cc)].at(j) != i_check) {
						variable_node[var_location[Vk.at(cc)].at(j)].at(Vk.at(cc)) = reg - check_node[var_location[Vk.at(cc)].at(j)].at(Vk.at(cc)) + channel_LLR.at(Vk.at(cc));
						//cout << "VN : " << Vk.at(cc) << " to CN :" << var_location[Vk.at(cc)].at(j) << " value : " << variable_node[var_location[Vk.at(cc)].at(j)].at(Vk.at(cc)) << endl;
					}
				}
			}
		}

		Chk_Eq = Parity_Check(check_node, channel_LLR);
		if (Chk_Eq == true)
		{
			break;
		}

	}


	vector<double> early_check_var(H_col, 0);

	for (int i = 0; i < H_col; i++)
	{
		double sum_var = 0;
		for (int j = 0; j < var_location[i].size(); j++)
		{
			sum_var = sum_var + check_node[var_location[i].at(j)].at(i);
		}
		early_check_var.at(i) = sum_var + channel_LLR.at(i);
		//cout << early_check_var.at(i) << " ";
	}

	return early_check_var;

}


int Q_lim = INT16_MAX;
vector<int> UpdateBuffer;
int globoIterMax = 10;

vector<double>  NWRBP_decode(vector<double>receive_bits, int itermax)
{

	int i_check = 0;
	int ko = 0;
	int pre;
	int P_check = 0;
	const int H_row = k, H_col = n;
	int Vk;
	int iternow = 0;
	double alpha = 0.3;
	vector<int> decoded_bits(k, 0);
	vector<int> decode_reg(n, 0);
	vector<int> check_node_cal(H_row, 0);
	vector<int> index_reg;
	vector<double> res(H_row, 0);
	vector<double> channel_LLR(n, 0);
	vector<double> check_node[H_row];
	vector<double> variable_node[H_row];
	vector<double> final_var(n, 0);
	vector<int> data;


	for (int i = 0; i < n; i++)
	{
		channel_LLR.at(i) = (-2.0 * receive_bits.at(i)) / (sigma * sigma);
	}


	for (int i = 0; i < H_row; i++)
	{
		check_node[i].resize(H_col);
		variable_node[i].resize(H_col);
	}


	int counts = 0;
	for (int i = (H_col - n); i < H_col; i++)
	{
		for (int j = 0; j < var_location[i].size(); j++)
		{
			variable_node[var_location[i].at(j)].at(i) = channel_LLR.at(counts);
		}
		counts++;
	}

	i_check = Residual1(variable_node);


	for (int i = 0; i < H_row; i++)
	{
		int   i_che = i;
		double ress = 0;
		for (int j = 0; j < residu[i_che].size(); j++)
		{
			ress = max(ress, fabs(residu[i_che].at(j)));
		}
		res.at(i_che) = fabs(ress);
	}

	for (int iter = 0; iter < itermax; iter++)
	{
		vector<int> Chk_buffer(H_row, 1);
		vector<int> Chk_history(H_row, 0);
		for (int i_r = 0; i_r < H_row; i_r++)
		{
			//cout<<i_check<<" ";
			check_node_cal.at(i_check)++;

			for (int j = 0; j < check_location[i_check].size(); j++)
			{
				double check_scr = 0;
				double s_scr = -1.0;
				for (int z = 0; z < check_location[i_check].size(); z++)
				{
					check_scr = check_scr + phi_func((variable_node[i_check].at(check_location[i_check].at(z))));
					s_scr = s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(z)));
				}
				check_node[i_check].at(check_location[i_check].at(j)) = (s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(j)))) * phi_func(check_scr - phi_func((variable_node[i_check].at(check_location[i_check].at(j)))));

				pre = i_check;
				res.at(pre) = 0;
				Chk_history.at(pre)++;
				//                residu[pre].at(j) = 0;
				Vk = check_location[i_check].at(j);
				if (Chk_history.at(pre) >= Q_lim) {
					Chk_buffer.at(pre) = 0;
				}
				for (int cc = 0; cc < var_location[Vk].size(); cc++)
				{
					if (var_location[Vk].at(cc) != i_check)
					{
						double reg = 0.0;
						for (int z = 0; z < var_location[Vk].size(); z++)
						{
							reg = reg + check_node[var_location[Vk].at(z)].at(Vk);
						}

						for (int z = 0; z < var_location[Vk].size(); z++)
						{
							variable_node[var_location[Vk].at(z)].at(Vk) = reg - check_node[var_location[Vk].at(z)].at(Vk) + channel_LLR.at(Vk);
						}

						int i_che = 0;
						i_che = var_location[Vk].at(cc);
						double check_scrr = 0.0;
						double  s_scrr = -1.0;

						for (int z = 0; z < check_location[i_che].size(); z++)
						{
							check_scrr = check_scrr + phi_func((variable_node[i_che].at(check_location[i_che].at(z))));
							s_scrr = s_scrr * sgn(variable_node[i_che].at(check_location[i_che].at(z)));
						}


						for (int z = 0; z < check_location[i_che].size(); z++)
						{
							new_R[i_che].at(z) = ((s_scrr * sgn(variable_node[i_che].at(check_location[i_che].at(z)))) * phi_func(check_scrr - phi_func((variable_node[i_che].at(check_location[i_che].at(z))))));

						}

						for (int z = 0; z < residu[i_che].size(); z++)
						{
							residu[i_che].at(z) = (fabs(new_R[i_che].at(z) - check_node[i_che].at(check_location[i_che].at(z))));
						}
						res.at(i_che) = *max_element(begin(residu[i_che]), end(residu[i_che]));
					}
				}
			}
			vector<double> tamp_res = res;
			for (int i = 0; i < tamp_res.size(); i++) {
				if (Chk_buffer[i] == 0) {
					tamp_res[i] = -1;
				}
				//tamp_res[i] = tamp_res[i] * Chk_buffer[i] - 1;
			}
			i_check = distance(begin(res), max_element(begin(res), end(res)));
		}
		//        system("pause");
		P_check = Parity_Check(check_node, channel_LLR);
		if (P_check == true)
		{
			break;
		}
	}

	counts = 0;
	for (int i = 0; i < H_col; i++)
	{
		double sum_var = 0;
		for (int j = 0; j < var_location[i].size(); j++)
		{
			sum_var = sum_var + check_node[var_location[i].at(j)].at(i);
		}

		final_var.at(counts) = sum_var + channel_LLR.at(counts);
		counts++;
	}


	return final_var;
}

vector<double>  RBP_decode(vector<double>receive_bits, int itermax)
{

	int i_check = 0;
	int ko = 0;
	int pre;
	int P_check = 0;
	const int H_row = k, H_col = n;
	int Vk;
	int iternow = 0;
	double alpha = 0.3;
	double beta = 0.9;
	vector<int> decoded_bits(k, 0);
	vector<int> decode_reg(n, 0);
	vector<int> check_node_cal(H_row, 0);
	vector<int> index_reg;
	vector<double> res(H_row, 0);
	vector<double> channel_LLR(n, 0);
	vector<double> check_node[H_row];
	vector<double> variable_node[H_row];
	vector<double> final_var(n, 0);
	vector<int> data;


	for (int i = 0; i < n; i++)
	{
		channel_LLR.at(i) = (-2.0 * receive_bits.at(i)) / (sigma * sigma);
	}


	for (int i = 0; i < H_row; i++)
	{
		check_node[i].resize(H_col);
		variable_node[i].resize(H_col);
	}


	int counts = 0;
	for (int i = (H_col - n); i < H_col; i++)
	{
		for (int j = 0; j < var_location[i].size(); j++)
		{
			variable_node[var_location[i].at(j)].at(i) = channel_LLR.at(counts);
		}
		counts++;
	}

	i_check = Residual1(variable_node);


	for (int i = 0; i < H_row; i++)
	{
		int   i_che = i;
		double ress = 0;
		for (int j = 0; j < residu[i_che].size(); j++)
		{
			ress = max(ress, fabs(residu[i_che].at(j)));
		}
		res.at(i_che) = fabs(ress);
	}

	for (int iter = 0; iter < itermax; iter++)
	{
		vector<int> Chk_buffer(H_row, 1);
		vector<int> Chk_history(H_row, 0);
		vector< vector<int> > edge_update_counter(H_row, vector<int>(H_col));
		for (int i_r = 0; i_r < Edgeamount; i_r++)
		{
			//cout<<i_check<<" ";
			check_node_cal.at(i_check)++;
			int i_Edge = distance(begin(residu[i_check]), max_element(begin(residu[i_check]), end(residu[i_check])));
			//for (int j = 0; j < check_location[i_check].size(); j++)
			{
				double check_scr = 0;
				double s_scr = -1.0;
				for (int z = 0; z < check_location[i_check].size(); z++)
				{
					check_scr = check_scr + phi_func((variable_node[i_check].at(check_location[i_check].at(z))));
					s_scr = s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(z)));
				}
				check_node[i_check].at(check_location[i_check].at(i_Edge)) = (s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(i_Edge)))) * phi_func(check_scr - phi_func((variable_node[i_check].at(check_location[i_check].at(i_Edge)))));
				edge_update_counter[i_check].at(check_location[i_check].at(i_Edge))++;
				pre = i_check;
				//res.at(pre) = 0;
				residu[i_check].at(i_Edge) = 0;
				res.at(i_check) = *max_element(begin(residu[i_check]), end(residu[i_check]));
				Chk_history.at(pre)++;
				//                residu[pre].at(j) = 0;
				Vk = check_location[i_check].at(i_Edge);
				if (Chk_history.at(pre) >= Q_lim) {
					Chk_buffer.at(pre) = 0;
				}
				for (int cc = 0; cc < var_location[Vk].size(); cc++)
				{
					if (var_location[Vk].at(cc) != i_check)
					{
						double reg = 0.0;
						for (int z = 0; z < var_location[Vk].size(); z++)
						{
							reg = reg + check_node[var_location[Vk].at(z)].at(Vk);
						}

						for (int z = 0; z < var_location[Vk].size(); z++)
						{
							variable_node[var_location[Vk].at(z)].at(Vk) = reg - check_node[var_location[Vk].at(z)].at(Vk) + channel_LLR.at(Vk);
						}

						int i_che = 0;
						i_che = var_location[Vk].at(cc);
						double check_scrr = 0.0;
						double  s_scrr = -1.0;

						for (int z = 0; z < check_location[i_che].size(); z++)
						{
							check_scrr = check_scrr + phi_func((variable_node[i_che].at(check_location[i_che].at(z))));
							s_scrr = s_scrr * sgn(variable_node[i_che].at(check_location[i_che].at(z)));
						}


						for (int z = 0; z < check_location[i_che].size(); z++)
						{
							new_R[i_che].at(z) = ((s_scrr * sgn(variable_node[i_che].at(check_location[i_che].at(z)))) * phi_func(check_scrr - phi_func((variable_node[i_che].at(check_location[i_che].at(z))))));

						}

						for (int z = 0; z < residu[i_che].size(); z++)
						{
							residu[i_che].at(z) = (fabs(new_R[i_che].at(z) - check_node[i_che].at(check_location[i_che].at(z))));
						}
						res.at(i_che) = *max_element(begin(residu[i_che]), end(residu[i_che]));
					}
				}
			}
			vector<double> tamp_res = res;
			for (int i = 0; i < tamp_res.size(); i++) {
				if (Chk_buffer[i] == 0) {
					tamp_res[i] = -1;
				}
				//tamp_res[i] = tamp_res[i] * Chk_buffer[i] - 1;
			}
			i_check = distance(begin(res), max_element(begin(res), end(res)));
		}
		//        system("pause");
		P_check = Parity_Check(check_node, channel_LLR);
		if (P_check == true)
		{
			break;
		}
	}

	counts = 0;
	for (int i = 0; i < H_col; i++)
	{
		double sum_var = 0;
		for (int j = 0; j < var_location[i].size(); j++)
		{
			sum_var = sum_var + check_node[var_location[i].at(j)].at(i);
		}

		final_var.at(counts) = sum_var + channel_LLR.at(counts);
		counts++;
	}


	return final_var;
}

vector<double>  RDNWRBP_decode(vector<double>receive_bits, int itermax)
{

	int i_check = 0;
	int ko = 0;
	int pre;
	int P_check = 0;
	const int H_row = k, H_col = n;
	int Vk;
	int iternow = 0;
	double alpha = 0.3;
	double beta = 0.9;
	vector<int> decoded_bits(k, 0);
	vector<int> decode_reg(n, 0);
	vector<int> check_node_cal(H_row, 0);
	vector<int> index_reg;
	vector<double> res(H_row, 0);
	vector<double> channel_LLR(n, 0);
	vector<double> check_node[H_row];
	vector<double> variable_node[H_row];
	vector<double> final_var(n, 0);
	vector<int> data;


	for (int i = 0; i < n; i++)
	{
		channel_LLR.at(i) = (-2.0 * receive_bits.at(i)) / (sigma * sigma);
	}


	for (int i = 0; i < H_row; i++)
	{
		check_node[i].resize(H_col);
		variable_node[i].resize(H_col);
	}


	int counts = 0;
	for (int i = (H_col - n); i < H_col; i++)
	{
		for (int j = 0; j < var_location[i].size(); j++)
		{
			variable_node[var_location[i].at(j)].at(i) = channel_LLR.at(counts);
		}
		counts++;
	}

	i_check = Residual1(variable_node);


	for (int i = 0; i < H_row; i++)
	{
		int   i_che = i;
		double ress = 0;
		for (int j = 0; j < residu[i_che].size(); j++)
		{
			ress = max(ress, fabs(residu[i_che].at(j)));
		}
		res.at(i_che) = fabs(ress);
	}

	for (int iter = 0; iter < itermax; iter++)
	{
		vector<int> Chk_buffer(H_row, 1);
		vector<int> Chk_history(H_row, 0);
		vector< vector<int> > edge_update_counter(H_row, vector<int>(H_col));
		for (int i_r = 0; i_r < H_row; i_r++)
		{
			//cout<<i_check<<" ";
			check_node_cal.at(i_check)++;

			for (int j = 0; j < check_location[i_check].size(); j++)
			{
				double check_scr = 0;
				double s_scr = -1.0;
				for (int z = 0; z < check_location[i_check].size(); z++)
				{
					check_scr = check_scr + phi_func((variable_node[i_check].at(check_location[i_check].at(z))));
					s_scr = s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(z)));
				}
				check_node[i_check].at(check_location[i_check].at(j)) = (s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(j)))) * phi_func(check_scr - phi_func((variable_node[i_check].at(check_location[i_check].at(j)))));
				edge_update_counter[i_check].at(check_location[i_check].at(j))++;
				pre = i_check;
				res.at(pre) = 0;
				Chk_history.at(pre)++;
				//                residu[pre].at(j) = 0;
				Vk = check_location[i_check].at(j);
				if (Chk_history.at(pre) >= Q_lim) {
					Chk_buffer.at(pre) = 0;
				}
				for (int cc = 0; cc < var_location[Vk].size(); cc++)
				{
					if (var_location[Vk].at(cc) != i_check)
					{
						double reg = 0.0;
						for (int z = 0; z < var_location[Vk].size(); z++)
						{
							reg = reg + check_node[var_location[Vk].at(z)].at(Vk);
						}

						for (int z = 0; z < var_location[Vk].size(); z++)
						{
							variable_node[var_location[Vk].at(z)].at(Vk) = reg - check_node[var_location[Vk].at(z)].at(Vk) + channel_LLR.at(Vk);
						}

						int i_che = 0;
						i_che = var_location[Vk].at(cc);
						double check_scrr = 0.0;
						double  s_scrr = -1.0;

						for (int z = 0; z < check_location[i_che].size(); z++)
						{
							check_scrr = check_scrr + phi_func((variable_node[i_che].at(check_location[i_che].at(z))));
							s_scrr = s_scrr * sgn(variable_node[i_che].at(check_location[i_che].at(z)));
						}


						for (int z = 0; z < check_location[i_che].size(); z++)
						{
							new_R[i_che].at(z) = ((s_scrr * sgn(variable_node[i_che].at(check_location[i_che].at(z)))) * phi_func(check_scrr - phi_func((variable_node[i_che].at(check_location[i_che].at(z))))));

						}

						for (int z = 0; z < residu[i_che].size(); z++)
						{
							residu[i_che].at(z) = pow(beta, edge_update_counter[i_che].at(check_location[i_che].at(z))) * (fabs(new_R[i_che].at(z) - check_node[i_che].at(check_location[i_che].at(z))));
						}
						res.at(i_che) = *max_element(begin(residu[i_che]), end(residu[i_che]));
					}
				}
			}
			vector<double> tamp_res = res;
			for (int i = 0; i < tamp_res.size(); i++) {
				if (Chk_buffer[i] == 0) {
					tamp_res[i] = -1;
				}
				//tamp_res[i] = tamp_res[i] * Chk_buffer[i] - 1;
			}
			i_check = distance(begin(res), max_element(begin(res), end(res)));
		}
		//        system("pause");
		P_check = Parity_Check(check_node, channel_LLR);
		if (P_check == true)
		{
			break;
		}
	}

	counts = 0;
	for (int i = 0; i < H_col; i++)
	{
		double sum_var = 0;
		for (int j = 0; j < var_location[i].size(); j++)
		{
			sum_var = sum_var + check_node[var_location[i].at(j)].at(i);
		}

		final_var.at(counts) = sum_var + channel_LLR.at(counts);
		counts++;
	}


	return final_var;
}

vector<double>  RDRBP_decode(vector<double>receive_bits, int itermax)
{

	int i_check = 0;
	int ko = 0;
	int pre;
	int P_check = 0;
	const int H_row = k, H_col = n;
	int Vk;
	int iternow = 0;
	double alpha = 0.3;
	double beta = 0.9;
	vector<int> decoded_bits(k, 0);
	vector<int> decode_reg(n, 0);
	vector<int> check_node_cal(H_row, 0);
	vector<int> index_reg;
	vector<double> res(H_row, 0);
	vector<double> channel_LLR(n, 0);
	vector<double> check_node[H_row];
	vector<double> variable_node[H_row];
	vector<double> final_var(n, 0);
	vector<int> data;


	for (int i = 0; i < n; i++)
	{
		channel_LLR.at(i) = (-2.0 * receive_bits.at(i)) / (sigma * sigma);
	}


	for (int i = 0; i < H_row; i++)
	{
		check_node[i].resize(H_col);
		variable_node[i].resize(H_col);
	}


	int counts = 0;
	for (int i = (H_col - n); i < H_col; i++)
	{
		for (int j = 0; j < var_location[i].size(); j++)
		{
			variable_node[var_location[i].at(j)].at(i) = channel_LLR.at(counts);
		}
		counts++;
	}

	i_check = Residual1(variable_node);


	for (int i = 0; i < H_row; i++)
	{
		int   i_che = i;
		double ress = 0;
		for (int j = 0; j < residu[i_che].size(); j++)
		{
			ress = max(ress, fabs(residu[i_che].at(j)));
		}
		res.at(i_che) = fabs(ress);
	}

	for (int iter = 0; iter < itermax; iter++)
	{
		vector<int> Chk_buffer(H_row, 1);
		vector<int> Chk_history(H_row, 0);
		vector< vector<int> > edge_update_counter(H_row, vector<int>(H_col));
		for (int i_r = 0; i_r < Edgeamount; i_r++)
		{
			//cout<<i_check<<" ";
			check_node_cal.at(i_check)++;
			int i_Edge = distance(begin(residu[i_check]), max_element(begin(residu[i_check]), end(residu[i_check])));
			//for (int j = 0; j < check_location[i_check].size(); j++)
			{
				double check_scr = 0;
				double s_scr = -1.0;
				for (int z = 0; z < check_location[i_check].size(); z++)
				{
					check_scr = check_scr + phi_func((variable_node[i_check].at(check_location[i_check].at(z))));
					s_scr = s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(z)));
				}
				check_node[i_check].at(check_location[i_check].at(i_Edge)) = (s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(i_Edge)))) * phi_func(check_scr - phi_func((variable_node[i_check].at(check_location[i_check].at(i_Edge)))));
				edge_update_counter[i_check].at(check_location[i_check].at(i_Edge))++;
				pre = i_check;
				//res.at(pre) = 0;
				residu[i_check].at(i_Edge) = 0;
				res.at(i_check) = *max_element(begin(residu[i_check]), end(residu[i_check]));
				Chk_history.at(pre)++;
				//                residu[pre].at(j) = 0;
				Vk = check_location[i_check].at(i_Edge);
				if (Chk_history.at(pre) >= Q_lim) {
					Chk_buffer.at(pre) = 0;
				}
				for (int cc = 0; cc < var_location[Vk].size(); cc++)
				{
					if (var_location[Vk].at(cc) != i_check)
					{
						double reg = 0.0;
						for (int z = 0; z < var_location[Vk].size(); z++)
						{
							reg = reg + check_node[var_location[Vk].at(z)].at(Vk);
						}

						for (int z = 0; z < var_location[Vk].size(); z++)
						{
							variable_node[var_location[Vk].at(z)].at(Vk) = reg - check_node[var_location[Vk].at(z)].at(Vk) + channel_LLR.at(Vk);
						}

						int i_che = 0;
						i_che = var_location[Vk].at(cc);
						double check_scrr = 0.0;
						double  s_scrr = -1.0;

						for (int z = 0; z < check_location[i_che].size(); z++)
						{
							check_scrr = check_scrr + phi_func((variable_node[i_che].at(check_location[i_che].at(z))));
							s_scrr = s_scrr * sgn(variable_node[i_che].at(check_location[i_che].at(z)));
						}


						for (int z = 0; z < check_location[i_che].size(); z++)
						{
							new_R[i_che].at(z) = ((s_scrr * sgn(variable_node[i_che].at(check_location[i_che].at(z)))) * phi_func(check_scrr - phi_func((variable_node[i_che].at(check_location[i_che].at(z))))));

						}

						for (int z = 0; z < residu[i_che].size(); z++)
						{
							residu[i_che].at(z) = pow(beta, edge_update_counter[i_che].at(check_location[i_che].at(z))) * (fabs(new_R[i_che].at(z) - check_node[i_che].at(check_location[i_che].at(z))));
						}
						res.at(i_che) = *max_element(begin(residu[i_che]), end(residu[i_che]));
					}
				}
			}
			vector<double> tamp_res = res;
			for (int i = 0; i < tamp_res.size(); i++) {
				if (Chk_buffer[i] == 0) {
					tamp_res[i] = -1;
				}
				//tamp_res[i] = tamp_res[i] * Chk_buffer[i] - 1;
			}
			i_check = distance(begin(res), max_element(begin(res), end(res)));
		}
		//        system("pause");
		P_check = Parity_Check(check_node, channel_LLR);
		if (P_check == true)
		{
			break;
		}
	}

	counts = 0;
	for (int i = 0; i < H_col; i++)
	{
		double sum_var = 0;
		for (int j = 0; j < var_location[i].size(); j++)
		{
			sum_var = sum_var + check_node[var_location[i].at(j)].at(i);
		}

		final_var.at(counts) = sum_var + channel_LLR.at(counts);
		counts++;
	}


	return final_var;
}

template <typename T>
vector<size_t> sort_indexes(const vector<T>& v) {

	// initialize original index locations
	vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	// using std::stable_sort instead of std::sort
	// to avoid unnecessary index re-orderings
	// when v contains elements of equal values 
	stable_sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return fabs(v[i1]) < fabs(v[i2]); });

	return idx;
}

vector<double>Perturbation_add_noise(vector<double>channel_LLR) {
	vector<size_t> Low_index = sort_indexes(channel_LLR);
	for (int i = 0; i < L_bits_smallest; i++)
	{
		channel_LLR.at(Low_index[i]) = (channel_LLR.at(Low_index[i]) + (sigma_p * Normal()));
		//                   receive_bits.at(i) = (codeword.at(i));
	}
	return channel_LLR;
}

vector<double>  RDRBP_Perturbation_decode(vector<double>receive_bits, int itermax, int Perturbation_N)
{

	int i_check = 0;
	int ko = 0;
	int pre;
	int P_check = 0;
	const int H_row = k, H_col = n;
	int Vk;
	int iternow = 0;
	double alpha = 0.3;
	double beta = 0.9;
	vector<int> decoded_bits(k, 0);
	vector<int> decode_reg(n, 0);
	vector<int> check_node_cal(H_row, 0);
	vector<int> index_reg;
	vector<double> res(H_row, 0);
	vector<double> channel_LLR(n, 0);
	vector<double> check_node[H_row];
	vector<double> variable_node[H_row];
	vector<double> final_var(n, 0);
	vector<int> data;


	for (int i = 0; i < n; i++)
	{
		channel_LLR.at(i) = (-2.0 * receive_bits.at(i)) / (sigma * sigma);
	}
	vector<double> Perturbation_channel_LLR = channel_LLR;
	for (int per_N = 0; per_N < Perturbation_N; per_N++) {



		for (int i = 0; i < H_row; i++)
		{
			check_node[i].resize(H_col);
			variable_node[i].resize(H_col);
		}


		int counts = 0;
		for (int i = (H_col - n); i < H_col; i++)
		{
			for (int j = 0; j < var_location[i].size(); j++)
			{
				variable_node[var_location[i].at(j)].at(i) = Perturbation_channel_LLR.at(counts);
			}
			counts++;
		}

		i_check = Residual1(variable_node);


		for (int i = 0; i < H_row; i++)
		{
			int   i_che = i;
			double ress = 0;
			for (int j = 0; j < residu[i_che].size(); j++)
			{
				ress = max(ress, fabs(residu[i_che].at(j)));
			}
			res.at(i_che) = fabs(ress);
		}

		for (int iter = 0; iter < itermax; iter++)
		{
			vector<int> Chk_buffer(H_row, 1);
			vector<int> Chk_history(H_row, 0);
			vector< vector<int> > edge_update_counter(H_row, vector<int>(H_col));
			for (int i_r = 0; i_r < Edgeamount; i_r++)
			{
				//cout<<i_check<<" ";
				check_node_cal.at(i_check)++;
				int i_Edge = distance(begin(residu[i_check]), max_element(begin(residu[i_check]), end(residu[i_check])));
				//for (int j = 0; j < check_location[i_check].size(); j++)
				{
					double check_scr = 0;
					double s_scr = -1.0;
					for (int z = 0; z < check_location[i_check].size(); z++)
					{
						check_scr = check_scr + phi_func((variable_node[i_check].at(check_location[i_check].at(z))));
						s_scr = s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(z)));
					}
					check_node[i_check].at(check_location[i_check].at(i_Edge)) = (s_scr * sgn(variable_node[i_check].at(check_location[i_check].at(i_Edge)))) * phi_func(check_scr - phi_func((variable_node[i_check].at(check_location[i_check].at(i_Edge)))));
					edge_update_counter[i_check].at(check_location[i_check].at(i_Edge))++;
					pre = i_check;
					//res.at(pre) = 0;
					residu[i_check].at(i_Edge) = 0;
					res.at(i_check) = *max_element(begin(residu[i_check]), end(residu[i_check]));
					Chk_history.at(pre)++;
					//                residu[pre].at(j) = 0;
					Vk = check_location[i_check].at(i_Edge);
					if (Chk_history.at(pre) >= Q_lim) {
						Chk_buffer.at(pre) = 0;
					}
					for (int cc = 0; cc < var_location[Vk].size(); cc++)
					{
						if (var_location[Vk].at(cc) != i_check)
						{
							double reg = 0.0;
							for (int z = 0; z < var_location[Vk].size(); z++)
							{
								reg = reg + check_node[var_location[Vk].at(z)].at(Vk);
							}

							for (int z = 0; z < var_location[Vk].size(); z++)
							{
								variable_node[var_location[Vk].at(z)].at(Vk) = reg - check_node[var_location[Vk].at(z)].at(Vk) + Perturbation_channel_LLR.at(Vk);
							}

							int i_che = 0;
							i_che = var_location[Vk].at(cc);
							double check_scrr = 0.0;
							double  s_scrr = -1.0;

							for (int z = 0; z < check_location[i_che].size(); z++)
							{
								check_scrr = check_scrr + phi_func((variable_node[i_che].at(check_location[i_che].at(z))));
								s_scrr = s_scrr * sgn(variable_node[i_che].at(check_location[i_che].at(z)));
							}


							for (int z = 0; z < check_location[i_che].size(); z++)
							{
								new_R[i_che].at(z) = ((s_scrr * sgn(variable_node[i_che].at(check_location[i_che].at(z)))) * phi_func(check_scrr - phi_func((variable_node[i_che].at(check_location[i_che].at(z))))));

							}

							for (int z = 0; z < residu[i_che].size(); z++)
							{
								residu[i_che].at(z) = pow(beta, edge_update_counter[i_che].at(check_location[i_che].at(z))) * (fabs(new_R[i_che].at(z) - check_node[i_che].at(check_location[i_che].at(z))));
							}
							res.at(i_che) = *max_element(begin(residu[i_che]), end(residu[i_che]));
						}
					}
				}
				vector<double> tamp_res = res;
				for (int i = 0; i < tamp_res.size(); i++) {
					if (Chk_buffer[i] == 0) {
						tamp_res[i] = -1;
					}
					//tamp_res[i] = tamp_res[i] * Chk_buffer[i] - 1;
				}
				i_check = distance(begin(res), max_element(begin(res), end(res)));
			}
			//        system("pause");
			P_check = Parity_Check(check_node, Perturbation_channel_LLR);
			if (P_check == true)
			{
				counts = 0;
				for (int i = 0; i < H_col; i++)
				{
					double sum_var = 0;
					for (int j = 0; j < var_location[i].size(); j++)
					{
						sum_var = sum_var + check_node[var_location[i].at(j)].at(i);
					}

					final_var.at(counts) = sum_var + Perturbation_channel_LLR.at(counts);
					counts++;
				}
				return final_var;
			}
		}
		Perturbation_channel_LLR = Perturbation_add_noise(channel_LLR);
	}
	//counts = 0;
	for (int i = 0; i < H_col; i++)
	{
		double sum_var = 0;
		for (int j = 0; j < var_location[i].size(); j++)
		{
			sum_var = sum_var + check_node[var_location[i].at(j)].at(i);
		}

		final_var.at(i) = sum_var + Perturbation_channel_LLR.at(i);
		//counts++;
	}
	return final_var;
}

long int err_amount = 0;

int Decoder_Box(vector<double>receive_bits, int itermax, int decoder_case) {
	switch (decoder_case) {

	case 1:
		return decision(BP_decode_flooding(receive_bits, itermax));
		break;
	case 2:
		return decision(BP_decode_LBP_O(receive_bits, itermax));
		break;
	case 3:
		return decision(NWRBP_decode(receive_bits, itermax));
		break;
	case 4:
		return decision(RBP_decode(receive_bits, itermax));
		break;
	case 5:
		return decision(RDNWRBP_decode(receive_bits, itermax));
		break;
	case 6:
		return decision(RDRBP_decode(receive_bits, itermax));
		break;
	case 7:
		return decision(RDRBP_Perturbation_decode(receive_bits, itermax, 30));
		break;
	default:
		cout << "Non - Decoder";
		system("pause");
		break;
	}
	return 0;
}


void main() {
	srand(time(NULL) + rand());
	read_BP_location();
	double Case_mode = 1;
	double snr_db = 2;
	int decoder_case = 1;
	cout << "1. Flooding  2. LBP  3. NWRBP 4. RBP 5. RD-NWRBP 6. RD-RBP 7. RD-RBP(Perturbation): ";
	cin >> decoder_case;
	cout << "SNR : ";
	cin >> snr_db;
	int itermax = 20;
	cout << "iter max : ";
	cin >> itermax;
	double coderate = 0.5;
	vector<int> codeword(n, -1);
	vector<double> receive_bits(n, 0);
	vector<double> decoder_LLR(H_col, 0);
	vector<int> decode_bits(k, 0);
	int sig;
	double iter_counts = 0;
	double average_iter = 0;
	double kl = 1000;
	double blockerror = 0;
	time_t nStart = time(NULL), nend;
	double times = 0.0;
	for (int i = 0; i < k; i++) {
		info.push_back(0);
	}
	//for (snr_db = 2; snr_db <= 3; snr_db += 0.5) 
	{
		double ErrAmount = 0.0, BlockErr = 0.0, errtamp;
		for (int max_times = 0; max_times < 1000000; max_times++)
		{
			times = max_times + 1;

			sigma = sqrt(1.0 / (2.0 * coderate * (pow(10.0, snr_db / 10.0))));

			for (int i = 0; i < n; i++)
			{
				receive_bits.at(i) = (codeword.at(i) + (sigma * Normal()));
				//                   receive_bits.at(i) = (codeword.at(i));
			}
			errtamp = Decoder_Box(receive_bits, itermax, decoder_case);
			if (errtamp != 0) {
				BlockErr++;
				ErrAmount += errtamp;
				errtamp = 0;
			}
			if (times >= 1) {
				cout << "\r Err : " << ErrAmount << " BER : " << ErrAmount / (k * times) << " Times : " << (int)times << "  ";
			}
			if (ErrAmount >= 3000) {
				break;
			}
		}
		time_t nEnd = time(NULL), nend;
		cout << endl;
		cout << "SNR:" << snr_db << endl;
		cout << "itermax:" << itermax << endl;
		cout << "times:" << (int)times << endl;
		cout << "bit errorcount:" << ErrAmount << endl;
		cout << "bit errorrate:" << ErrAmount / (k * times) << endl;
		cout << "block errorcount:" << BlockErr << endl;
		cout << "block errorrate:" << BlockErr / times << endl;
		cout << "timecost:" << nEnd - nStart << endl;
		cout << endl;
	}

	system("pause");


}