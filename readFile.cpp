#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <memory>
#include <numeric>
#include <utility>
#include "readFile.h"
//Данные, полученные из текстового файла
int part_vel = 0; //Парт к которому приложена скорость
double init_vz = 0.0; //Скорость начальная
double kr = 0.0; //Число Куранта
double lamdal = 0.0; //Линейная псевдовязкость
double lamdak = 0.0; //Квадратичная псевдовязкость
double T_lim = 0.0; //Время расчета
double T_zap = 0.0; //Шаг записи
const double PI = 3.1415926535; 

class Cel
{
public: 
	double p; //Давление
	double c; //Скорость звука
	double qkomb, qrr, qzz, q00, qrz; //Комбинированная псевдовязкость
	double err, ezz, e00, erz; //Скорость изменения деформаций
	double delta_err, delta_ezz, delta_e00, delta_erz; //Приращения деформаций
	double Drr, Dzz, D00, Drz; //Девиаторы напряжений
	double sigma_rr, sigma_zz, sigma_00, sigma_rz, intensive_sig; //Компоненты напряжения
	double E; //Энергия
	double perimeter; // Периметр
	double A, A12; //Площадь ячейки
	double r0; // Координата х центра ячейки
	double z0; // Координата z центра ячейки
	double M; //Масса ячейки
	double V, V12, V12tV; //Объем ячейки
	double ro0; //Начальная плотность ячейки
	double ro; //Плотность ячейки
	double roAV0;
private:
	int YID; //ID ячейки
	int PID; //ID парта, к которому принадлежит ячейка
	std:: vector<int> NID;//ID узлов, принадлежащих ячейки
	int boundCheck = -1;
	int boundClass = -1;
public:
	Cel(/* args */);
	~Cel();
};




struct Cell
{
	int YID; //ID ячейки
	int PID; //ID парта, к которому принадлежит ячейка
	std:: vector<int> NID;//ID узлов, принадлежащих ячейки
	int boundCheck = -1;
	int boundClass = -1;
	double p; //Давление
	double c; //Скорость звука
	double qkomb, qrr, qzz, q00, qrz; //Комбинированная псевдовязкость
	double err, ezz, e00, erz; //Скорость изменения деформаций
	double delta_err, delta_ezz, delta_e00, delta_erz; //Приращения деформаций
	double Drr, Dzz, D00, Drz; //Девиаторы напряжений
	double sigma_rr, sigma_zz, sigma_00, sigma_rz, intensive_sig; //Компоненты напряжения
	double E; //Энергия
	double perimeter; // Периметр
	double A, A12; //Площадь ячейки
	double r0; // Координата х центра ячейки
	double z0; // Координата z центра ячейки
	double M; //Масса ячейки
	double V, V12, V12tV; //Объем ячейки
	double ro0; //Начальная плотность ячейки
	double ro; //Плотность ячейки
	double roAV0;
};
struct Node
{
	int UID; //ID узла
	int checkPress = 0;
	int boundCheck = -1;
	double r, z; //Старые координаты узла
	double rn, zn; //Новые координаты узла
	double vr, vz; //Старая скорость узла
	double vrn, vzn; //Новая скорость узла
	std:: vector<std:: shared_ptr<Cell>> idCells; //ID ячеек, в которых числится узел
};
struct Part
{
	int id;
	int n_mat;
	std:: vector<std:: shared_ptr<Cell>> bondaryCells; //Граничные ячейки парта
};
struct PseudoCell
{
	double p1 = 0.0;
	double c1 = 0.0;
	double err_n = 0.0,
			ezz_n = 0.0,
			e00_n = 0.0,
			erz_n = 0.0;
	double delta_err_n = 0.0,
			delta_ezz_n = 0.0,
			delta_e00_n = 0.0,
			delta_erz_n = 0.0;
	double Drr_n = 0.0,
			Dzz_n = 0.0,
			D00_n = 0.0,
			Drz_n = 0.0;
	double E_n = 0.0;
	double sigma_rr_n = 0.0,
			sigma_zz_n = 0.0,
			sigma_00_n = 0.0,
			sigma_rz_n = 0.0;
	double perimeter1 = 0.0;
	double A1 = 0.0;
	double r01 = 0.0;
	double z01 = 0.0;
	double V1 = 0.0;
	double ro1 = 0.0;
	double A12 = 0.0;
	double V12 = 0.0;
	double V12tV = 0.0;
	double qkomb_n = 0.0,
			qrr_n = 0.0,
			qzz_n = 0.0,
			q00_n = 0.0,
			qrz_n = 0.0;
};
class Mat
{
public:
	int id;
	double ro; //Плотность
	double Eu; //Модуль Юнга
	double kpuas; //Коэффициенет Пуассона
	double sigmay; //Предел текучести
	double G;
	double gamm[100];
	double aa[100], zz[100];
	double aa1[100], zz1[100];
};
class Me : public Mat
{
public:
	Me();
	void urs(int nmat, const Cell& cell, Cell& pseucell);
	void rachetcoeff(Cell& cell);
};
struct Pressure
{
	std:: vector<double> time;
	std:: vector<double> press;
};

std:: pair<double, std:: vector<double>> Perimeter(const std:: vector<int>& nodeCell, const std:: vector<Node>& nodes)
{
	std:: vector<double> lengsid;
	lengsid.reserve(3);
	for (size_t i = 0; i + 1 != nodeCell.size(); ++i)
	{
		const Node& node1 = nodes[nodeCell[i] - 1];
		for (size_t j = i + 1; j != nodeCell.size(); ++j)
		{
			const Node& node2 = nodes[nodeCell[j] - 1];
			double delr = node1.r - node2.r;
			double delz = node1.z - node2.z;
			lengsid.push_back(sqrt(delr * delr + delz * delz)); 
		}
	}
	return {std:: accumulate(lengsid.begin(), lengsid.end(), 0.0), lengsid};
}

void InitGeometry(Cell& cell, std:: vector<Part>& parts, std:: vector<Mat>& mats, std:: vector<Node>& nodes) 
	{
		cell.perimeter = Perimeter(cell.NID, nodes).first;
		std:: vector<double> lengsid = Perimeter(cell.NID, nodes).second;
		double per = 0.5 * cell.perimeter;
		cell.ro0 = mats[parts[cell.PID - 1].n_mat - 1].ro;
		cell.A = sqrt(per * (per - lengsid[0]) * (per - lengsid[1]) * (per - lengsid[2]));

		cell.r0 = 0.0; cell.z0 = 0.0;
		for (size_t i = 0; i != cell.NID.size(); ++i)
		{
			cell.r0 += nodes[cell.NID[i] - 1].r / 3;
			cell.z0 += nodes[cell.NID[i] - 1].z / 3;
		}
		cell.M = 2.0 * PI * cell.A * cell.r0 * cell.ro0;
		cell.V = 1.0;
		cell.ro = cell.ro0;
		cell.roAV0 = cell.ro * cell.A; 
	}

void SkipLine(std:: string& s, std::istream& fin, int count)
{	
	for(int i = 0; i != count; ++i)
	{
		std:: getline(fin, s);
	}
}
void SkipValue(std:: string& s, std::istream& fin, int count)
{
	for(int i = 0; i != count; ++i)
	{
		fin >> s;
	}
}

bool CheckBound(const std:: vector<int>& nodeCell, std:: vector<Node>& nodes)
{
	int prevIndex = nodeCell.size() - 2;
	int nextIndex = nodeCell.size() - 1;
	int boundCell = 0;
	for (size_t i = 0; i != nodeCell.size(); ++i)
	{
		int boundNode = 0;	
		Node& node = nodes[nodeCell[i] - 1];
		for (size_t idCell = 0; idCell != node.idCells.size(); ++idCell)
		{
			if (nodeCell[prevIndex] == node.idCells[idCell] -> NID[0] || nodeCell[prevIndex] == node.idCells[idCell] -> NID[1] || 
				nodeCell[prevIndex] == node.idCells[idCell] -> NID[2] || nodeCell[nextIndex] == node.idCells[idCell] -> NID[0] || 
				nodeCell[nextIndex] == node.idCells[idCell] -> NID[1] || nodeCell[nextIndex] == node.idCells[idCell] -> NID[2])
			{
				++boundCell;
				++boundNode;
			}
		}
		if (boundNode < 3 && boundNode != 1)
		{
			node.boundCheck = 1;
		}
		prevIndex = 0;
		nextIndex -= i;
	}
	if (boundCell <= 7)
	{
		return true;
	}
	return false;
}
bool CheckAxis(const std:: vector<int>& nodeCell, std:: vector<Node>& nodes)
{
	for (size_t i = 0; i + 1 != nodeCell.size(); ++i)
	{
		for (size_t j = i + 1; j != nodeCell.size(); ++j)
		{
			Node& node1 = nodes[nodeCell[i] - 1];
			Node& node2 = nodes[nodeCell[j] - 1];
			if (node1.r < 1e-8 && node2.r < 1e-8 && node1.idCells.size() != 1 && node2.idCells.size() != 1 )
			{
				node1.boundCheck = 0;
				node2.boundCheck = 0;
				return true;
			}
		}
	}
	return false;
}
void BoundClass(Cell& cell, std:: vector<Node>& nodes)
{
	Node& node1 = nodes[cell.NID[0] - 1];
	Node& node2 = nodes[cell.NID[1] - 1];
	Node& node3 = nodes[cell.NID[2] - 1];
	if(node1.idCells.size() == 1 || node2.idCells.size() == 1 || node3.idCells.size() == 1)
	{
		cell.boundClass = 1;
	}
	else if (node1.boundCheck == 1 && node2.boundCheck == 1 && node3.boundCheck == 1)
	{
		cell.boundClass = 2;
	}
	else
	{
		cell.boundClass = 0;
	}
}

void ReadFile(std:: vector<Part>& parts, std:: vector<Mat>& mats, std:: vector<Cell>& cells, std:: vector<Node>& nodes, Pressure& press)
{
	std:: string part;
	std:: cin >> part;
	std:: ifstream fin(part);
	std:: string str;
	while (fin >> str) 
	{
		if (str.find("*INITIAL_VEL") != str.npos)
		{
			SkipLine(str, fin, 2);
			fin >> part_vel;
			SkipValue(str, fin, 3);
			fin >> init_vz;
		}
		else if (str.find("*PART") != str.npos)
		{
			while (str.find("pid") == str.npos)
			{
				std:: getline(fin, str);
			}
			Part buf;
			fin >> buf.id >> str >> buf.n_mat;
			parts.push_back(std:: move(buf));
		}
		else if (str.find("*MAT_PL") != str.npos)
		{
			SkipLine(str, fin, 3);
			Mat buf;
			fin >> buf.id >> buf.ro >> buf.Eu >> buf.kpuas >> buf.sigmay;
			buf.G = 0.5 * buf.Eu / (1 + buf.kpuas);
			mats.push_back(std:: move(buf));
		}
		else if (str.find("*DATABASE_B") != str.npos)
		{
			SkipLine(str, fin, 2);
			fin >> T_zap;
		}
		else if (str.find("*CONTROL_BULK") != str.npos)
		{
			SkipLine(str, fin, 2);
			fin >> lamdak >> lamdal;
		}
		else if (str.find("*CONTROL_TERM") != str.npos)
		{
			SkipLine(str, fin, 2);
			fin >> T_lim;
		}
		else if (str.find("*CONTROL_TIME") != str.npos)
		{
			SkipLine(str, fin, 2);
			fin >> str >> kr;
		}
		else if (str.find("*NODE") != str.npos)
		{
			int index = 1;
			SkipLine(str, fin, 2);
			fin >> str;
			do
			{
				if (index == stoi(str))
				{
					Node buf;
					buf.UID = index;
					fin >> buf.r >> buf.z;
					nodes.push_back(buf);
					getline(fin, str);
					fin >> str;
				}
				else
				{
					Node buf;
					buf.UID = -1;
					nodes.push_back(std:: move(buf));
				}
				index++;
			} while (str.find("*") == str.npos);
			fin.seekg(-1000, std::ios::cur);
			getline(fin, str);
		}
		else if (str.find("*ELEMENT") != str.npos)
		{
			SkipLine(str, fin, 2);
			fin >> str;
			do
			{
				Cell buf;
				buf.NID.resize(3);
				buf.YID = stoi(str);
				fin >> buf.PID >> buf.NID[0] >> buf.NID[1] >> buf.NID[2];
				cells.push_back(std:: move(buf));
				nodes[buf.NID[0] - 1].idCells.emplace_back(buf.YID);
				nodes[buf.NID[1] - 1].idCells.emplace_back(buf.YID);
				nodes[buf.NID[2] - 1].idCells.emplace_back(buf.YID);
				getline(fin, str);
				fin >> str;
			} while (str.find("*") == str.npos);
			fin.seekg(-1000, std::ios::cur);
			getline(fin, str);
		}
		else if (str.find("*DEFINE_C") != str.npos)
		{
			SkipLine(str, fin, 4);
			fin >> str;
			do
			{
				press.time.push_back(stod(str));
				fin >> str;
				press.press.push_back(stod(str));
				fin >> str;
			} while (str.find("*") == str.npos);
			fin.seekg(-100, std::ios::cur);
			getline(fin, str);
		}
		else if (str.find("*SET_SEG") != str.npos)
		{
			SkipLine(str, fin, 4);
			fin >> str;
			do
			{
				nodes[stoi(str) - 1].checkPress = 1;
				fin >> str;
				nodes[stoi(str) - 1].checkPress  = 1;
				getline(fin, str);
				fin >> str;
			} while (str.find("*") == str.npos);
			fin.seekg(-1000, std::ios::cur);
			getline(fin, str);
		}
	}
}

void PostRead(std:: vector<Part>& parts, std:: vector<Mat>& mats, std:: vector<Cell>& cells, std:: vector<Node>& nodes)
{
	for (size_t i = 0; i != cells.size(); ++i)
	{
		if (CheckAxis(cells[i].NID, nodes))
		{
			cells[i].boundCheck = 0;
		}
		else if (CheckBound(cells[i].NID, nodes))
		{
			cells[i].boundCheck = 1;
			BoundClass(cells[i], nodes);
			parts[cells[i].PID - 1].bondaryCells.push_back(std:: make_unique<Cell>(cells[i]));
		}
		InitGeometry(cells[i], parts, mats, nodes);
	}

	if (init_vz != 0.0)
	{
		for (size_t i = 0; i != cells.size() && cells[i].PID == part_vel; ++i)
		{
			for (size_t j = 0; j < cells[i].NID.size(); j++)
			{
				nodes[cells[i].NID[j] - 1].vz = init_vz;
			}	
		}
	}
}