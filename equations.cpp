#include <cmath> 
#include "equations.h"

double deltaT = 100.;
double deltaT12 = 0.0;

// Функции апрокссимации узлов в центр ячейки

double ApproxNode(double vel1, double vel2, double coord1, double coord2)
{
	return (vel2 + vel1) * (coord2 - coord1);
}

//Функции апрокссимации ячеек в узел

double ApproxCell(double stress1, double stress2, const Node& node1, const Node& node2) 
{
    return (stress1 * (node2.z - node1.z) - stress2 * (node2.r - node1.r));
}

//===================================
//Расчет геометрии 
//===================================
void Geometry(const Cell& cell, Cell& pseucell, const std:: vector<Node>& nodes) 
{
	pseucell.perimeter = Perimeter(cell.NID, nodes).first;
	std:: vector<double> lengsid = Perimeter(cell.NID, nodes).second;
	double per = 0.5 * cell.perimeter;
	pseucell.A = sqrt(per * (per - lengsid[0]) * (per - lengsid[1]) * (per - lengsid[2]));
	for (size_t i = 0; i != cell.NID.size(); ++i)
	{
		pseucell.r0 += nodes[cell.NID[i] - 1.].rn;
		pseucell.z0 += nodes[cell.NID[i] - 1.].zn;
	}
	pseucell.r0 /= 3; pseucell.z0 /= 3;
	pseucell.V = 2 * PI * cell.ro0 * pseucell.A * pseucell.r0 / cell.M;
	pseucell.ro = cell.ro0 / pseucell.V;
}

//===================================
//УРС
//===================================
Me::Me()
{
	aa[5] = 3.02e+10; zz[5] = 4.8; gamm[5] = 0; //����
	aa[8] = 1.97e+10; zz[8] = 4.2; gamm[8] = 0; //��������
	aa[11] = 2.15e+10; zz[11] = 5.5; gamm[11] = 0; //�����
	aa[12] = 1.8797e9; zz[12] = 4.668; gamm[12] = 0; //����
	aa[13] = 3.595e8; zz[13] = 1.736; gamm[13] = 0; //���������� 320
	aa[14] = 0.86e+10; zz[14] = 5.3; gamm[14] = 0; //������
	aa[15] = 4.7353e6; zz[15] = 3.46; gamm[15] = 0; //���������� 160
	aa[17] = 1.4506e+11; zz[17] = 1.987; gamm[17] = 0; //������
	aa[19] = 1.955e+10; zz[19] = 3.648; gamm[19] = 0; //LiF
	//aa[19] = 1.472e+9; zz[19] = 8.76; gamm[19] = 0; //������ ��������� ��
	//aa[19] = 1.97e+10; zz[19] = 4.2; gamm[19] = 0;
	aa[21] = 5.71e+10; zz[21] = 4.7; gamm[21] = 0; // ���-90
	aa[22] = 2.3193e8; zz[22] = 3.539; gamm[22] = 0; //���������� 626
	aa[24] = 2.15e+10; zz[24] = 5.5; gamm[24] = 2.6; //�����
};
void Me::rachetcoeff(Cell& cell)
{
	for (int i = 5; i <= 99; i++)
	{
		zz1[i] = zz[i] - 1.0;
		aa1[i] = aa[i] * zz[i] / cell.ro0;
	}
};
void Me::urs(int nmat, const Cell& cell, Cell& pseucell)
{
	double a1 = aa[nmat];
	double z1 = zz[nmat];
	double az1 = aa1[nmat];
	double z11 = zz1[nmat];
	pseucell.p = a1 * (pow((pseucell.ro / cell.ro0), z1) - 1);

	double as2 = az1 * pow((pseucell.ro / cell.ro0), z11);
	if (as2 > 0)
	{
		pseucell.c = sqrt(as2);
	}
	else
	{
		pseucell.c = 100.;
	}
};

//===================================
//Расчет шага по времени
//===================================
void Kurant(const Cell& cell) 
{
    double deltat =  kr * 2.0 * cell.A / (cell.perimeter * cell.c);    
    if(deltat < deltaT) 
	{
        deltaT12 = (deltaT + deltat)/2.0;
        deltaT = deltat;
    }
};

void FreeNodesSC(Node& node, const std:: vector<Node>& nodes, std::vector<Cell>& cells)
{
	double approx_r = 0.0,
			approx_z = 0.0,
			omega = 0.0,
			beta = 0.0,
			gama = 0.0;
	for (size_t i = 0; i != node.idCells.size(); ++i)
	{
		Cell& cell = cells[node.idCells[i] - 1];
		if (cell.NID[0] == node.UID)
		{
			approx_r += ApproxCell(cell.sigma_rr, cell.sigma_rz, nodes[cell.NID[1] - 1.], nodes[cell.NID[2] - 1.]);
			approx_z += ApproxCell(cell.sigma_rz, cell.sigma_zz, nodes[cell.NID[1] - 1.], nodes[cell.NID[2] - 1.]);
		}
		else if (cell.NID[1] == node.UID)
		{
			approx_r += ApproxCell(cell.sigma_rr, cell.sigma_rz, nodes[cell.NID[2] - 1.], nodes[cell.NID[1] - 1.]);
			approx_z += ApproxCell(cell.sigma_rz, cell.sigma_zz, nodes[cell.NID[2] - 1.], nodes[cell.NID[1] - 1.]);
		}
		else
		{
			approx_r += ApproxCell(cell.sigma_rr, cell.sigma_rz, nodes[cell.NID[1] - 1.], nodes[cell.NID[0] - 1.]);
			approx_z += ApproxCell(cell.sigma_rz, cell.sigma_zz, nodes[cell.NID[1] - 1.], nodes[cell.NID[0] - 1.]);
		}
		omega += cell.ro * cell.A;
		beta += (cell.sigma_rr - cell.sigma_00) * cell.A / cell.M;
		gama += cell.sigma_rz * cell.A / cell.M;
	}
	omega /= 3;
	beta /= node.idCells.size();
	gama /= node.idCells.size();
	node.vrn = node.vr + (0.5 * approx_r / omega + beta) * deltaT;
	node.vzn = node.vz + (0.5 * approx_z / omega + gama) * deltaT;
	node.rn = node.r + node.vrn * deltaT12;
	node.zn = node.z + node.vzn * deltaT12;
}

//===================================
//Уравнение движения
//===================================
void SpeedCoord(double time, std:: vector<Node>& nodes, std:: vector<Cell>& cells, std:: vector<Part>& parts, Pressure& press)
{
	for (Node& node : nodes)
	{
		if (node.UID == -1)
		{
			node.r = 0.0;
			node.z = 0.0;
			node.vr = 0.0;
			node.vz = 0.0;
			continue;
		}
	
		switch (node.boundCheck)
		{
		case -1:
		case 1:
		{
			FreeNodesSC(node, nodes, cells);
			break;
		}
		case 0:
		{
			FreeNodesSC(node, nodes, cells);
			node.rn = 0.0;
			node.vrn = 0.0;
			break;
		}
		default:
			break;
		}
		//Контактный алгоритм не доступен в публичной версии. Обращаться по почте lumendany@mail.ru
	}
}

//===================================
//Скорости деформаций и деформации 
//===================================
void Defor(const Cell& cell, Cell& pseucell, std:: vector<Node>& nodes, const double time) 
{
    //УСМ
    pseucell.A12 = (pseucell.A + cell.A) / 2;
    pseucell.V12 = (pseucell.V + cell.V) / 2;
    if (time != 0.0) 
	{
		pseucell.V12tV = (pseucell.V - cell.V) / (deltaT * pseucell.V12);
    } 
	else 
	{
        pseucell.V12tV = 0;
    }
    
	for (size_t i = 0; i + 1 < cell.NID.size(); ++i)
	{
		Node& node1 = nodes[cell.NID[i] - 1.];
		for (size_t j = i + 1; i < cell.NID.size(); ++j)
		{
			Node& node2 = nodes[cell.NID[j] - 1.];
			if(i == 0 && j == 2)
			{
				pseucell.err += ApproxNode(node2.vrn, node1.vrn, node2.zn + node2.z, node1.zn + node1.z);
				pseucell.ezz += ApproxNode(node2.vzn, node1.vzn, node1.rn + node1.r, node2.rn + node2.r);
				pseucell.erz += ApproxNode(node2.vzn, node1.vzn, node2.zn + node2.z, node1.zn + node1.z) + 
									ApproxNode(node2.vrn, node1.vrn, node1.rn + node1.r, node2.rn + node2.r);
			}
			else
			{
				pseucell.err += ApproxNode(node1.vrn, node2.vrn, node1.zn + node1.z, node2.zn + node2.z);
				pseucell.ezz += ApproxNode(node2.vzn, node1.vzn, node2.rn + node2.r, node1.rn + node1.r);
				pseucell.erz += ApproxNode(node2.vzn, node1.vzn, node1.zn + node1.z, node2.zn + node2.z) + 
									ApproxNode(node2.vrn, node1.vrn, node2.rn + node2.r, node1.rn + node1.r);
			}
		}
	}
	pseucell.err /= 4 * pseucell.A12;
    pseucell.ezz /= 4 * pseucell.A12;
    pseucell.erz /= 4 * pseucell.A12;
    pseucell.e00 = pseucell.V12tV - pseucell.err - pseucell.ezz;
    
	pseucell.delta_err = pseucell.err * deltaT;
    pseucell.delta_ezz = pseucell.ezz * deltaT;
    pseucell.delta_e00 = pseucell.e00 * deltaT;
    pseucell.delta_erz = pseucell.erz * deltaT;
};

//===================================
//Расчет псевдовязкости
//===================================
void PseudoVis(const Cell& cell, Cell& pseucell) 
{
    if (pseucell.V12tV < 0) 
	{
        double qk = lamdak * lamdak * cell.ro0 * pseucell.A12 * pseucell.V12tV * pseucell.V12tV / pseucell.V12;
        double ql = pseucell.c * lamdal * cell.ro0 * std:: sqrt(pseucell.A12) * std:: abs(pseucell.V12tV) / pseucell.V12;
        pseucell.qkomb = qk + ql;
        //Вязкая добавка
        double mu = 0.01 * cell.ro0 * sqrt(pseucell.A12) / pseucell.V12;
        pseucell.qrr = 2 * mu * (pseucell.err - pseucell.V12tV / 3);
        pseucell.qzz = 2 * mu * (pseucell.ezz - pseucell.V12tV / 3);
        pseucell.q00 = 2 * mu * (pseucell.e00 - pseucell.V12tV / 3);
        pseucell.qrz = mu * pseucell.erz;
    } 
	else 
	{
        pseucell.qkomb = 0.0;
        pseucell.qrr = 0.0;
        pseucell.qzz = 0.0;
        pseucell.q00 = 0.0;
        pseucell.qrz = 0.0;
    }
};

//===================================
//Расчет девиаторов напряжения
//===================================
void Stress(const Cell& cell, Cell& pseucell, std:: vector<Node>& nodes, std:: vector<Part>& parts, std:: vector<Mat>& mats)
{
	double deltarz = 0.0;
	for (size_t i = 0; i + 1 < cell.NID.size(); ++i)
	{
		Node& node1 = nodes[cell.NID[i] - 1.];
		for (size_t j = i + 1; i < cell.NID.size(); ++j)
		{
			Node& node2 = nodes[cell.NID[j] - 1.];
			if(i == 0 && j == 2)
			{
				deltarz += ApproxNode(node2.vzn + node2.vz, node1.vzn + node1.vz, node1.z, node2.z) + 
							ApproxNode(node2.vrn + node2.vr, node1.vrn + node1.vr, node1.r, node2.r);
			}
			else
			{
				deltarz += ApproxNode(node2.vzn + node2.vz, node1.vzn + node1.vz, node2.z, node1.z) + 
							ApproxNode(node2.vrn + node2.vr, node1.vrn + node1.vr, node2.r, node1.r);
			}
		}
	}
    deltarz /= 4 * cell.A;
    
	double drr = cell.Drr * deltarz * deltaT; 
    double drz = 0.5 * (cell.Dzz - cell.Drr) * deltarz * deltaT;
	double& G = mats[parts[cell.PID - 1.].n_mat - 1.].G;
	double& sigmay = mats[parts[cell.PID - 1.].n_mat - 1.].sigmay;
    
	double intermval = - pseucell.V12tV * deltaT / 3.0;
	pseucell.Drr = cell.Drr + 2.0 * G * (pseucell.delta_err + intermval) + drr;
    pseucell.Dzz = cell.Dzz + 2.0 * G * (pseucell.delta_ezz + intermval) - drr;
    pseucell.D00 = cell.D00 + 2.0 * G * (pseucell.delta_e00 + intermval);
    pseucell.Drz = cell.Drz + G * pseucell.delta_erz + drz;
    
    double f = 2 * (pseucell.Drr * pseucell.Drr + pseucell.Dzz * pseucell.Dzz + pseucell.Drz * pseucell.Drz + pseucell.Drr * pseucell.Dzz);
    if (f > 2 * sigmay * sigmay / 3) {
        double F = sigmay * sqrt(2 / (3 * f));
        pseucell.Drr *= F;
        pseucell.Dzz *= F;
        pseucell.Drz *= F;
        pseucell.D00 *= F;
    }

    //Напряжения 
    pseucell.sigma_rr = (pseucell.Drr + pseucell.qrr) - (pseucell.p + pseucell.qkomb);
    pseucell.sigma_zz = (pseucell.Dzz + pseucell.qzz) - (pseucell.p + pseucell.qkomb);
    pseucell.sigma_00 = (pseucell.D00 + pseucell.q00) - (pseucell.p + pseucell.qkomb);
    pseucell.sigma_rz = (pseucell.Drz + pseucell.qrz); 
	
};

//===================================
//Уравнение энергии
//===================================
void Energy(Cell& cell, Cell& pseucell) 
{
    double D12rr = 0.5 * (pseucell.Drr + cell.Drr); 
    double D12zz = 0.5 * (pseucell.Dzz + cell.Dzz); 
    double D1200 = 0.5 * (pseucell.D00 + cell.D00); 
    double D12rz = 0.5 * (pseucell.Drz + cell.Drz); 
    double Z = pseucell.V12 * (D12rr * pseucell.err + D12zz * pseucell.ezz + D1200 * pseucell.e00 + D12rz * pseucell.erz);
    pseucell.E = cell.E - (cell.p / 2.0 + pseucell.qkomb) * (pseucell.V - cell.V) + Z * deltaT;
}

