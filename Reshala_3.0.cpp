#include <vector>
#include <iostream>
//#include "zapis_object.h"
//#include "uravnenia.h"
#include "writeFile.h"
//#include "vivod_xml.cpp"

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
		
			SkorKoord(t, nodes, cells, parts, press);
			
			for (size_t i = 0; i != cells.size(); ++i) 
			{
				Cell pseucell = cells[i];
				Geometry(cells[i], pseucell, nodes);
				Me Met;
				Met.rachetcoeff(cells[i]);
				if (Ych[i].PID == 8)
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
				if ((isnan(Ych[i].A) == true) || (isinf(Ych[i].A) == true) || (isnan(Ych[i].ro) == true) || (isinf(Ych[i].ro) == true) ||
					(isnan(Ych[i].p) == true) || (isinf(Ych[i].p) == true) || (isnan(Ych[i].c) == true) || (isinf(Ych[i].c) == true) ||
					(isnan(Ych[i].err) == true) || (isinf(Ych[i].err) == true) || (isnan(Ych[i].erz) == true) || (isinf(Ych[i].erz) == true) ||
					(isnan(Ych[i].ezz) == true) || (isinf(Ych[i].ezz) == true) || (isnan(Ych[i].e00) == true) || (isinf(Ych[i].e00) == true) ||
					(isnan(Ych[i].delta_err) == true) || (isinf(Ych[i].delta_err) == true) || (isnan(Ych[i].delta_erz) == true) || (isinf(Ych[i].delta_erz) == true) ||
					(isnan(Ych[i].delta_e00) == true) || (isinf(Ych[i].delta_e00) == true) || (isnan(Ych[i].delta_ezz) == true) || (isinf(Ych[i].delta_ezz) == true) ||
					(isnan(Ych[i].Drr) == true) || (isinf(Ych[i].Drr) == true) || (isnan(Ych[i].Drz) == true) || (isinf(Ych[i].Drz) == true) ||
					(isnan(Ych[i].Dzz) == true) || (isinf(Ych[i].Dzz) == true) || (isnan(Ych[i].D00) == true) || (isinf(Ych[i].D00) == true) ||
					(isnan(Ych[i].sigma_rr) == true) || (isinf(Ych[i].sigma_rr) == true) || (isnan(Ych[i].sigma_rz) == true) || (isinf(Ych[i].sigma_rz) == true) ||
					(isnan(Ych[i].sigma_zz) == true) || (isinf(Ych[i].sigma_zz) == true) || (isnan(Ych[i].sigma_00) == true) || (isinf(Ych[i].sigma_00) == true) ||
					(isnan(Ych[i].sigma_rz) == true) || (isinf(Ych[i].sigma_rz) == true) || (isnan(Ych[i].sigma_rz) == true) || (isinf(Ych[i].sigma_rz) == true) ||
					(isnan(Ych[i].qkomb) == true) || (isinf(Ych[i].qkomb) == true) || (isnan(Ych[i].q00) == true) || (isinf(Ych[i].q00) == true) ||
					(isnan(Ych[i].qrr) == true) || (isinf(Ych[i].qrr) == true) || (isnan(Ych[i].qrz) == true) || (isinf(Ych[i].qrz) == true) ||
					(isnan(Ych[i].qzz) == true) || (isinf(Ych[i].qzz) == true) || (isnan(Ych[i].E) == true) || (isinf(Ych[i].E) == true) ||
					(isnan(Ych[i].M) == true) || (isinf(Ych[i].M) == true) || (isnan(Ych[i].V) == true) || (isinf(Ych[i].V) == true) ||
					(isnan(Ych[i].r0) == true) || (isinf(Ych[i].r0) == true) || (isnan(Ych[i].z0) == true) || (isinf(Ych[i].z0) == true))
				{
					vivod(step, Uz, Ych, Pr);
					std::cout << "NOMBER: " << Ych[i].YID << "\n";
					std::cout << "Площадь: " << Ych[i].A << "\n";
					std::cout << "Объем: " << Ych[i].V << "\n";
					std::cout << "Плотность: " << Ych[i].ro << "\n";
					std::cout << "Давление: " << Ych[i].p << "\n";
					std::cout << "Скорость звука: " << Ych[i].c << "\n";
					std::cout << "Площадь на половинном слое: " << A12 << "\n";
					std::cout << "Объем на половинном слое: " << V12 << "\n";
					std::cout << "Удельная скорость изменения объема: " << V12tV << "\n";
					std::cout << "Скорость деформации rr: " << Ych[i].err << "\n";
					std::cout << "Скорость деформации zz: " << Ych[i].ezz << "\n";
					std::cout << "Скорость деформации rz: " << Ych[i].erz << "\n";
					std::cout << "Скорость деформации 00: " << Ych[i].e00 << "\n";
					std::cout << "Деформации rr: " << Ych[i].delta_err << "\n";
					std::cout << "Деформации zz: " << Ych[i].delta_ezz << "\n";
					std::cout << "Деформации rz: " << Ych[i].delta_erz << "\n";
					std::cout << "Деформации 00: " << Ych[i].delta_e00 << "\n";
					std::cout << "Девиатор rr: " << Ych[i].Drr << "\n";
					std::cout << "Девиатор zz: " << Ych[i].Dzz << "\n";
					std::cout << "Девиатор rz: " << Ych[i].Drz << "\n";
					std::cout << "Девиатор 00: " << Ych[i].D00 << "\n";
					std::cout << "Напряжение rr: " << Ych[i].sigma_rr << "\n";
					std::cout << "Напряжение zz: " << Ych[i].sigma_zz << "\n";
					std::cout << "Напряжение rz: " << Ych[i].sigma_rz << "\n";
					std::cout << "Напряжение 00: " << Ych[i].sigma_00 << "\n";
					std::cout << "Псевдовязкость: " << Ych[i].qkomb << "\n";
					std::cout << "Псевдовязкость rr: " << Ych[i].qrr << "\n";
					std::cout << "Псевдовязкость zz: " << Ych[i].qzz << "\n";
					std::cout << "Псевдовязкость rz: " << Ych[i].qrz << "\n";
					std::cout << "Псевдовязкость 00: " << Ych[i].q00 << "\n";
					goto q101;
				}
			}
			for (unsigned int i = 0; i < Uz.size(); i++) {
				Uz[i].vr = Uz[i].vrn;
				Uz[i].vz = Uz[i].vzn;
				Uz[i].r = Uz[i].rn;
				Uz[i].z = Uz[i].zn;
			}
			for (unsigned int i = 0; i < Ych.size(); i++) {
				Kurant(cells[i]);
			}
		
		if (t >= step * T_zap)
		{
			vivod(step, Uz, Ych, Pr);
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