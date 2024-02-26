#include <fstream>
#include <iostream>
#include "equations.h"
#include "writeFile.h"

int main() 
{
	std:: vector <Part> parts;
	std:: vector <Mat> mats;
	std:: vector <Cell> cells;
	std:: vector <Node> nodes;
	Pressure press;
	ReadFile(parts, mats, cells, nodes, press);
	PostRead(parts, mats, cells, nodes);
	std:: cout << "Initial data\n" << "Calculation time: " << T_lim << "\n" << "Recording step: " << T_zap << "\n" <<
				"Quadratic pseudoviscosity: " << lamdak << "\n" << "Linear pseudoviscosity: " << lamdal << "\n" << 
				"Kurant: " << kr << "\n==============\n";
	double t = 0.0;
	double step = 0.0;
	int shag = 0;
	while (t < T_lim) 
	{
		std:: cout << "Step:" << shag << "\n";
		
		SpeedCoord(t, nodes, cells, parts, press);

		for (size_t i = 0; i != cells.size(); ++i) 
		{
			Cell pseucell = cells[i];
			Geometry(cells[i], pseucell, nodes);
			Me Met;
			Met.rachetcoeff(cells[i]);
			if (cells[i].PID == 8)
			{
				Met.urs(8, cells[i], pseucell);
			}
			else
			{
				Met.urs(11, cells[i], pseucell);
			}
			Defor(cells[i], pseucell, nodes, t);
			PseudoVis(cells[i], pseucell);
			Stress(cells[i], pseucell, nodes, parts, mats);
			Energy(cells[i], pseucell);
			cells[i] = std:: move(pseucell);
			if ((isnan(cells[i].A) == true) || (isinf(cells[i].A) == true) || (isnan(cells[i].ro) == true) || (isinf(cells[i].ro) == true) ||
				(isnan(cells[i].p) == true) || (isinf(cells[i].p) == true) || (isnan(cells[i].c) == true) || (isinf(cells[i].c) == true) ||
				(isnan(cells[i].err) == true) || (isinf(cells[i].err) == true) || (isnan(cells[i].erz) == true) || (isinf(cells[i].erz) == true) ||
				(isnan(cells[i].ezz) == true) || (isinf(cells[i].ezz) == true) || (isnan(cells[i].e00) == true) || (isinf(cells[i].e00) == true) ||
				(isnan(cells[i].delta_err) == true) || (isinf(cells[i].delta_err) == true) || (isnan(cells[i].delta_erz) == true) || (isinf(cells[i].delta_erz) == true) ||
				(isnan(cells[i].delta_e00) == true) || (isinf(cells[i].delta_e00) == true) || (isnan(cells[i].delta_ezz) == true) || (isinf(cells[i].delta_ezz) == true) ||
				(isnan(cells[i].Drr) == true) || (isinf(cells[i].Drr) == true) || (isnan(cells[i].Drz) == true) || (isinf(cells[i].Drz) == true) ||
				(isnan(cells[i].Dzz) == true) || (isinf(cells[i].Dzz) == true) || (isnan(cells[i].D00) == true) || (isinf(cells[i].D00) == true) ||
				(isnan(cells[i].sigma_rr) == true) || (isinf(cells[i].sigma_rr) == true) || (isnan(cells[i].sigma_rz) == true) || (isinf(cells[i].sigma_rz) == true) ||
				(isnan(cells[i].sigma_zz) == true) || (isinf(cells[i].sigma_zz) == true) || (isnan(cells[i].sigma_00) == true) || (isinf(cells[i].sigma_00) == true) ||
				(isnan(cells[i].sigma_rz) == true) || (isinf(cells[i].sigma_rz) == true) || (isnan(cells[i].sigma_rz) == true) || (isinf(cells[i].sigma_rz) == true) ||
				(isnan(cells[i].qkomb) == true) || (isinf(cells[i].qkomb) == true) || (isnan(cells[i].q00) == true) || (isinf(cells[i].q00) == true) ||
				(isnan(cells[i].qrr) == true) || (isinf(cells[i].qrr) == true) || (isnan(cells[i].qrz) == true) || (isinf(cells[i].qrz) == true) ||
				(isnan(cells[i].qzz) == true) || (isinf(cells[i].qzz) == true) || (isnan(cells[i].E) == true) || (isinf(cells[i].E) == true) ||
				(isnan(cells[i].M) == true) || (isinf(cells[i].M) == true) || (isnan(cells[i].V) == true) || (isinf(cells[i].V) == true) ||
				(isnan(cells[i].r0) == true) || (isinf(cells[i].r0) == true) || (isnan(cells[i].z0) == true) || (isinf(cells[i].z0) == true))
			{
				write(step, nodes, cells, parts);
				std::cout << "NOMBER: " << cells[i].YID << "\n";
				std::cout << "Площадь: " << cells[i].A << "\n";
				std::cout << "Объем: " << cells[i].V << "\n";
				std::cout << "Плотность: " << cells[i].ro << "\n";
				std::cout << "Давление: " << cells[i].p << "\n";
				std::cout << "Скорость звука: " << cells[i].c << "\n";
				std::cout << "Скорость деформации rr: " << cells[i].err << "\n";
				std::cout << "Скорость деформации zz: " << cells[i].ezz << "\n";
				std::cout << "Скорость деформации rz: " << cells[i].erz << "\n";
				std::cout << "Скорость деформации 00: " << cells[i].e00 << "\n";
				std::cout << "Деформации rr: " << cells[i].delta_err << "\n";
				std::cout << "Деформации zz: " << cells[i].delta_ezz << "\n";
				std::cout << "Деформации rz: " << cells[i].delta_erz << "\n";
				std::cout << "Деформации 00: " << cells[i].delta_e00 << "\n";
				std::cout << "Девиатор rr: " << cells[i].Drr << "\n";
				std::cout << "Девиатор zz: " << cells[i].Dzz << "\n";
				std::cout << "Девиатор rz: " << cells[i].Drz << "\n";
				std::cout << "Девиатор 00: " << cells[i].D00 << "\n";
				std::cout << "Напряжение rr: " << cells[i].sigma_rr << "\n";
				std::cout << "Напряжение zz: " << cells[i].sigma_zz << "\n";
				std::cout << "Напряжение rz: " << cells[i].sigma_rz << "\n";
				std::cout << "Напряжение 00: " << cells[i].sigma_00 << "\n";
				std::cout << "Псевдовязкость: " << cells[i].qkomb << "\n";
				std::cout << "Псевдовязкость rr: " << cells[i].qrr << "\n";
				std::cout << "Псевдовязкость zz: " << cells[i].qzz << "\n";
				std::cout << "Псевдовязкость rz: " << cells[i].qrz << "\n";
				std::cout << "Псевдовязкость 00: " << cells[i].q00 << "\n";
				goto q101;
			}
		}
		for (Node& node : nodes) {
			node.vr = node.vrn;
			node.vz = node.vzn;
			node.r = node.rn;
			node.z = node.zn;
		}
		for (Cell& cell : cells) {
			Kurant(cell);
		}
		
		if (t >= step * T_zap)
		{
			write(step, nodes, cells, parts);
			step++;
		}
		shag++;
		t += deltaT;
		std::cout << "Time Step: " << deltaT << "\n";
		std::cout << "New time: " << t << "\n"; 
		std::cout << "==============" << "\n";
	};
	q101:;
	std::cout << "End" << "\n";
	return 0;
}